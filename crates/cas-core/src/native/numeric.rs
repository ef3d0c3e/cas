//! `N[e]` and `N[e, d]` (CAS-DOC-5 §9.2, §9.4, CAS-84): the value of `e`
//! with balls for its numbers, certified to `d` digits, or to the working
//! precision.
//!
//! `e` is evaluated first, as any argument. Then, at `p` bits:
//! - its exact numbers, and the constants `Pi` and `E`, where they stand for
//!   values ([`Numeric`]), become balls of `p` bits;
//! - that is evaluated in [`Context::numeric`], at the working precision
//!   `p`, where comparisons of balls decide only what they prove;
//! - the result is done when each of its balls knows the `t` bits of the
//!   target, and no comparison was left undecided by overlapping balls.
//!   Otherwise `p` doubles, from `t + 32` up to `8t + 1024` bits and the
//!   limit on numbers: Ziv's loop, which recomputes the whole expression,
//!   at most twice the cost of its last round.
//!
//! Then each inexact ball of the result is made to show its digits. The
//! printer writes the shortest decimal inside a ball, of at most the digits
//! of the larger of the working precision and the ball's accuracy (§6.1),
//! so:
//! - for `N[e, d]` with `d` within the working precision's digits, the
//!   radius is raised to half a unit of the `d`-th digit, and the printer
//!   writes the nearest decimal of `d` digits. The target is reached when the
//!   radius is that small already, so that `N` of its own result changes
//!   nothing. A complex ball's parts each show their own digits, but are
//!   measured against the larger;
//! - for more digits, and for `N[e]`, the balls are rounded to the `t` bits
//!   of the target, as arithmetic rounds its balls to the working
//!   precision: the printer writes as many digits as the accuracy carries,
//!   `d` for `N[e, d]`, and those of the working precision for `N[e]`.
//!
//! Short of the target at the cap, the result comes with
//! [`Warning::Inaccurate`]. When doubling gains less than 2 bits, the balls
//! of the input limit the result: the loop stops, with the note
//! [`Warning::InputLimited`]. A comparison still undecided then leaves `N`
//! unevaluated, with [`Warning::Undecided`]: either more precision would
//! have told its sides apart, or they are equal and no precision can tell
//! (§9.4).

use cas_num::Integer;
use cas_num::ball::{Ball, ComplexBall, DecimalStyle, Float, Mag};
use foldhash::fast::FixedState;
use hashbrown::HashMap;

use super::ball::shown_decimal;
use crate::attributes::Attributes;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::print::shown_digits;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The bits of a target of `d` digits, `⌈d log₂ 10⌉ + 1` (§9.2): with an
/// upper bound of `log₂ 10`, a bit more at most.
fn digits_target(d: u64) -> u64 {
    let bits = (u128::from(d) * 33_219_280_949).div_ceil(10_000_000_000);
    u64::try_from(bits).unwrap_or(u64::MAX).saturating_add(1)
}

/// The most digits whose target fits in `max` bits.
fn max_digits(max: u64) -> u64 {
    let mut d = shown_digits(max);
    while d > 0 && digits_target(d) > max {
        d -= 1;
    }
    d
}

/// `d` as a positive count of digits.
fn digits(d: &Expr) -> Option<u64> {
    let q = d.to_rational()?;
    if !q.is_integer() || !q.is_positive() {
        return None;
    }
    u64::try_from(q.numer().to_i64()?).ok()
}

/// What `N` aims at: bits (for `N[e]`, and digits past the working
/// precision's), or digits that the working precision shows.
#[derive(Clone, Copy)]
enum Target {
    Bits(u64),
    Digits(u64),
}

impl Target {
    fn bits(self) -> u64 {
        match self {
            Target::Bits(t) => t,
            Target::Digits(d) => digits_target(d),
        }
    }

    /// Whether the balls of `e` know enough of their values.
    fn reached(self, e: &Expr) -> bool {
        match self {
            // The rounding to `t` bits may leave one less.
            Target::Bits(t) => accuracy(e) >= i64::try_from(t).unwrap_or(i64::MAX) - 1,
            Target::Digits(d) => e.distinct().all(|part| match part.view() {
                View::Ball(b) => b.is_exact() || within(b.rad(), half_unit(b.mid(), d).as_ref()),
                View::ComplexBall(z) => {
                    let unit = half_unit(larger(z), d);
                    [z.re(), z.im()]
                        .iter()
                        .all(|b| b.is_exact() || within(b.rad(), unit.as_ref()))
                }
                _ => true,
            }),
        }
    }

    /// `e` with each inexact ball made to show the digits of the target.
    fn finish(self, e: &Expr) -> Expr {
        match self {
            Target::Bits(t) => round(e, t),
            Target::Digits(d) => e.replace(|part| match part.view() {
                View::Ball(b) if !b.is_exact() => Some(Expr::from(&widened(b, d))),
                View::ComplexBall(z) if !z.is_exact() => {
                    let z = ComplexBall::new(widened(z.re(), d), widened(z.im(), d));
                    Some(Expr::from(&z))
                }
                _ => None,
            }),
        }
    }
}

/// Whether `rad` is at most `unit`, when there is one.
fn within(rad: &Mag, unit: Option<&Mag>) -> bool {
    unit.is_some_and(|u| rad <= u)
}

/// The midpoint of the larger part of `z`, which its digits are counted
/// from.
fn larger(z: &ComplexBall) -> &Float {
    let (re, im) = (z.re().mid(), z.im().mid());
    if re.is_zero() || !im.is_zero() && im.exponent() > re.exponent() {
        im
    } else {
        re
    }
}

/// Half a unit of the `d`-th significant digit of `m`, rounded up: the
/// radius at which a ball around `m` shows `m` to `d` digits, the printer
/// writing the shortest decimal inside a ball. `None` for 0, and for a value
/// whose decimal exponent is beyond `i64`.
fn half_unit(m: &Float, d: u64) -> Option<Mag> {
    if m.is_zero() {
        return None;
    }
    let s = Ball::exact(m.clone()).to_decimal(d, DecimalStyle::Digits);
    let (_, digits, exponent) = shown_decimal(&s)?;
    let length = i64::try_from(digits.to_string().len()).ok()?;
    // The leading digit's power of 10, then the `d`-th's.
    let unit = exponent
        .checked_add(length - 1)?
        .checked_sub(i64::try_from(d).ok()? - 1)?;
    let half = Ball::from_scaled(&Integer::from(5), unit.checked_sub(1)?, 64);
    Some(Mag::from_float(half.mid()).add(half.rad()))
}

/// `b` with a radius of half a unit of its `d`-th digit, if it is
/// smaller; as it is when it knows fewer digits, as the inexact 0 of a
/// nearly real value.
fn widened(b: &Ball, d: u64) -> Ball {
    if b.is_exact() {
        return b.clone();
    }
    match half_unit(b.mid(), d) {
        Some(unit) if *b.rad() < unit => Ball::new(b.mid().clone(), unit),
        _ => b.clone(),
    }
}

pub(super) fn n(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match target(ev, e) {
        Some((arg, target)) => approximate(ev, arg, target),
        None => Ok(None),
    }
}

/// The argument of `N[arg, d]` and its target; `None` for digits that are
/// no positive integer within the limit on numbers, with a warning, or
/// still to come.
fn target<'e>(ev: &mut Eval<'_>, e: &'e Expr) -> Option<(&'e Expr, Target)> {
    let max = ev.context().limits().number_bits;
    let working = ev.context().working_precision();
    match e.args() {
        [arg] => Some((arg, Target::Bits(working.min(max)))),
        [arg, d] => match digits(d).filter(|&d| digits_target(d) <= max) {
            Some(d) if d <= shown_digits(working) => Some((arg, Target::Digits(d))),
            Some(d) => Some((arg, Target::Bits(digits_target(d)))),
            None => {
                if d.is_number() {
                    ev.warn(Warning::Native {
                        message: format!(
                            "the number of digits must be a positive integer up to {}",
                            max_digits(max)
                        )
                        .into(),
                        expr: e.clone(),
                    });
                }
                None
            }
        },
        _ => None,
    }
}

/// Ziv's loop on `arg` for `target`.
fn approximate(ev: &mut Eval<'_>, arg: &Expr, target: Target) -> Native {
    let goal = target.bits();
    let cap = goal
        .saturating_mul(8)
        .saturating_add(1024)
        .min(ev.context().limits().number_bits)
        .max(goal);
    let inexact_input = arg.distinct().any(is_inexact);
    let indeterminate = Expr::from(Symbol::INDETERMINATE);
    let had_indeterminate = arg.contains(&indeterminate);
    let mut bits = goal.saturating_add(32).min(cap);
    let mut last: Option<i64> = None;
    loop {
        let numeric = Numeric::new(ev.context(), bits).value(arg);
        let (result, warnings) = ev.eval_in(&ev.context().numeric(bits), &numeric);
        let undecided = warnings
            .iter()
            .rfind(|w| matches!(w, Warning::Undecided { .. }))
            .cloned();
        let result = match result {
            Ok(value) => Some(value),
            // The condition of a procedure's `if` that the balls left open.
            Err(EvalError::UndecidedCondition) if undecided.is_some() => None,
            Err(err) => return Err(err),
        };
        let accuracy = match &result {
            // An indeterminate value that more precision may make finite.
            Some(value) if !had_indeterminate && value.contains(&indeterminate) => i64::MIN,
            Some(value) => accuracy(value),
            None => i64::MIN,
        };
        let done = accuracy != i64::MIN && result.as_ref().is_some_and(|v| target.reached(v));
        if undecided.is_none() && done {
            pass_on(ev, warnings);
            return Ok(result.map(|v| target.finish(&v)));
        }
        let stalled = inexact_input && last.is_some_and(|a| accuracy < a.saturating_add(2));
        if bits >= cap || stalled {
            pass_on(ev, warnings);
            let digits = shown_digits(u64::try_from(accuracy).unwrap_or(0));
            if stalled {
                ev.warn(Warning::InputLimited {
                    expr: arg.clone(),
                    digits,
                });
            }
            if let Some(w) = undecided {
                ev.warn(w);
                return Ok(None);
            }
            if !stalled {
                ev.warn(Warning::Inaccurate {
                    expr: arg.clone(),
                    digits,
                    bits,
                });
            }
            return Ok(result.map(|v| target.finish(&v)));
        }
        last = Some(accuracy);
        bits = bits.saturating_mul(2).min(cap);
    }
}

/// Gives the warnings of the last round, but those of comparisons left
/// undecided, which `N` reports itself.
fn pass_on(ev: &mut Eval<'_>, warnings: Vec<Warning>) {
    for w in warnings {
        if !matches!(w, Warning::Undecided { .. }) {
            ev.warn(w);
        }
    }
}

/// Whether `e` is a ball that is not exact.
fn is_inexact(e: &Expr) -> bool {
    match e.view() {
        View::Ball(b) => !b.is_exact(),
        View::ComplexBall(z) => !z.is_exact(),
        _ => false,
    }
}

/// The least accuracy, in bits, of the balls in `e`; `i64::MAX` without
/// one. A complex ball's is relative to its larger part, so that the inexact
/// 0 of a nearly real value costs nothing.
fn accuracy(e: &Expr) -> i64 {
    e.distinct()
        .filter_map(|part| match part.view() {
            View::Ball(b) => Some(b.rel_accuracy_bits()),
            View::ComplexBall(z) => Some(z.rel_accuracy_bits()),
            _ => None,
        })
        .min()
        .unwrap_or(i64::MAX)
}

/// `e` with its balls that know more than `bits` bits rounded to `bits`.
fn round(e: &Expr, bits: u64) -> Expr {
    let limit = i64::try_from(bits).unwrap_or(i64::MAX);
    e.replace(|part| match part.view() {
        View::Ball(b) if b.rel_accuracy_bits() > limit => {
            Some(Expr::from(&b.round_to_accuracy(bits)))
        }
        View::ComplexBall(z) if z.rel_accuracy_bits() > limit => {
            let z = ComplexBall::new(
                z.re().round_to_accuracy(bits),
                z.im().round_to_accuracy(bits),
            );
            Some(Expr::from(&z))
        }
        _ => None,
    })
}

/// What `N` makes of an expression at `p` bits: balls for its exact numbers
/// and for `Pi` and `E`, where they stand for values. That is in the
/// arguments of numeric functions (the `NumericFunction` heads: arithmetic
/// and the elementary functions), lists, comparisons, logic, `If` and
/// equations, as in Xcas's `evalf`. Elsewhere numbers keep their exact
/// values (`f[1/3]`, `x[2]`), and so do integer exponents (`x^2`), the units
/// `-1` and `±I` of a product (`x - y`), and `E` under an exponent that is
/// no number (`E^x`).
struct Numeric<'c> {
    context: &'c Context,
    bits: u64,
    /// The applications made so far: an expression shares its parts.
    memo: HashMap<Expr, Expr, FixedState>,
}

impl<'c> Numeric<'c> {
    fn new(context: &'c Context, bits: u64) -> Numeric<'c> {
        Numeric {
            context,
            bits,
            memo: HashMap::with_hasher(FixedState::default()),
        }
    }

    fn value(&mut self, e: &Expr) -> Expr {
        let p = self.bits;
        match e.view() {
            View::Int(_) | View::Integer(_) | View::Rational(_) => {
                let q = e.to_rational().expect("a rational");
                Expr::from(&Ball::from_rational(&q, p))
            }
            View::Complex(z) => Expr::from(&ComplexBall::from_complex(z, p)),
            View::Symbol(Symbol::PI) => Expr::from(&Ball::const_pi(p)),
            View::Symbol(Symbol::E) => Expr::from(&Ball::const_e(p)),
            View::Apply { head, args } if self.numeric_head(head) => {
                if let Some(done) = self.memo.get(e) {
                    return done.clone();
                }
                let h = head.as_symbol();
                let parts: Vec<Expr> = args
                    .iter()
                    .enumerate()
                    .map(|(i, a)| {
                        if kept(h, i, a, args) {
                            a.clone()
                        } else {
                            self.value(a)
                        }
                    })
                    .collect();
                let done = if parts == args {
                    e.clone()
                } else {
                    Expr::apply(head, &parts)
                };
                self.memo.insert(e.clone(), done.clone());
                done
            }
            _ => e.clone(),
        }
    }

    /// Whether the arguments of `head` stand for values.
    fn numeric_head(&self, head: &Expr) -> bool {
        head.as_symbol().is_some_and(|s| {
            matches!(
                s,
                Symbol::LIST
                    | Symbol::EQUATION
                    | Symbol::EQUAL
                    | Symbol::UNEQUAL
                    | Symbol::LESS
                    | Symbol::LESS_EQUAL
                    | Symbol::GREATER
                    | Symbol::GREATER_EQUAL
                    | Symbol::AND
                    | Symbol::OR
                    | Symbol::NOT
                    | Symbol::IF
            ) || self
                .context
                .attributes(s)
                .contains(Attributes::NUMERIC_FUNCTION)
        })
    }
}

/// Whether the argument `a`, the `i`-th (from 0) of `args` under `head`,
/// keeps its exact value: an integer exponent, a unit of a product, or `E`
/// under an exponent that is no number.
fn kept(head: Option<Symbol>, i: usize, a: &Expr, args: &[Expr]) -> bool {
    match (head, i) {
        (Some(Symbol::POWER), 1) => a.to_rational().is_some_and(|q| q.is_integer()),
        (Some(Symbol::POWER), 0) => {
            a.as_symbol() == Some(Symbol::E) && !args.get(1).is_some_and(Expr::is_numeric)
        }
        (Some(Symbol::TIMES), _) => is_unit(a),
        _ => false,
    }
}

/// Whether `a` is `-1`, `I` or `-I`.
fn is_unit(a: &Expr) -> bool {
    match a.view() {
        View::Complex(z) => z.re().is_zero() && z.im().abs().is_one(),
        _ => a.to_rational().is_some_and(|q| (-q).is_one()),
    }
}

#[cfg(test)]
mod tests;
