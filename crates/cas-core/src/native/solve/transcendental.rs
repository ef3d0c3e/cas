//! Equations in functions of the unknown (CAS-DOC-11 §4, CAS-213), where
//! the conversion has generators that hold `x` other than `x` itself.
//!
//! - **Substitution and inversion**: with one such generator `g` and no
//!   `x` alone, `f` is a rational function of `g`, solved in a fresh `u`,
//!   then `g = c` inverted for each root `c`, to giac's principal
//!   solutions: `exp(a) = c` is `a = ln c`, none for `c ≤ 0` in the reals
//!   (`csolve(exp(x)=-1,x)` is `[iπ]`); `b^a = c` is `a = ln c / ln b`, an
//!   integer when `c` is a power of `b` (`2^x = 8` is `3`); `ln a = c` is
//!   `a = e^c`; `sin a = c` is `asin c` and `π − asin c`, none for `|c| > 1`
//!   in the reals; `cos a = c` is `±acos c`; `tan a = c` is `atan c`; the
//!   inverse functions the other way, within their ranges; `|a| = c` is
//!   `±c`; `a^(±1/n) = c` is `a = c^(±n)`, none for an even `n` and `c < 0`.
//!   Then `a = v` solved for `x` in turn: `exp(2x) − 3exp(x) + 2` is
//!   `u² − 3u + 2` in `u = eˣ`, so `x = 0` or `ln 2`; `sin(x)² = 1/4` gives
//!   `±π/6`, `5π/6` and `7π/6`.
//! - **Radicals**: with a root `s = r^(1/n)` of something in `x` among the
//!   generators, the numerator as `A + B s` (for a square root, every power
//!   of `s` reduced by `s² = r`) and `(−A)ⁿ − Bⁿ r = 0` solved, repeated
//!   while roots are left; then each candidate checked in the equation it
//!   came from, kept when its value is 0 exactly or its ball holds 0, or
//!   when parameters leave it unknown: `√(x + 1) = x − 1` gives `3`, not
//!   `0`.
//! - Otherwise, as `x eˣ = 1` (giac's `LambertW`) or `sin x = cos x`, the
//!   call stays.

use std::cmp::Ordering;

use cas_num::algebra::{GaussianIntegers, Integers};
use cas_num::{Integer, Rational};
use cas_poly::{Sparse, SparseRing};

use super::super::approximation;
use super::super::polynomial::Numbers;
use super::{BITS, Goal, real_sign, roots_at};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;

/// The most steps of substitution, inversion and raising to a power.
const MAX_DEPTH: u32 = 8;

/// The largest power `b^k = c` is tried for, in `b^a = c`.
const MAX_EXACT_LOG: i64 = 4096;

type Found = Result<Option<Vec<Expr>>, EvalError>;

/// The roots of `f = 0`, converted as `p`, whose generators `others` hold
/// the unknown, which is a generator of its own when `alone`.
pub(super) fn roots(
    ev: &mut Eval<'_>,
    goal: &Goal<'_>,
    f: &Expr,
    p: &Polynomials,
    alone: bool,
    others: &[usize],
    depth: u32,
) -> Found {
    if depth >= MAX_DEPTH {
        return Ok(None);
    }
    if let (false, [j]) = (alone, others)
        && let Some(found) = substituted(ev, goal, p, *j, depth)?
    {
        return Ok(Some(found));
    }
    radicals(ev, goal, f, p, others, depth)
}

/// The fresh unknown of a substitution `depth` steps in.
fn fresh(depth: u32) -> Expr {
    Expr::from(Symbol::new(&format!("Solve$u{depth}")))
}

/// `f`, a rational function of the generator `j` of `p`, solved in it,
/// then that generator inverted at each root.
fn substituted(ev: &mut Eval<'_>, goal: &Goal<'_>, p: &Polynomials, j: usize, depth: u32) -> Found {
    let g = &p.generators()[j];
    let u = fresh(depth);
    let built = p.fraction_expr(&p.fractions()[0]);
    if built.contains(&u) {
        return Ok(None);
    }
    let fu = ev.eval(&built.substitute(&[(g.clone(), u.clone())]))?;
    let inner = Goal {
        e: goal.e,
        x: &u,
        real: goal.real,
    };
    let Some(values) = roots_at(ev, &inner, &fu, depth + 1)? else {
        return Ok(None);
    };
    // An identity in `u`.
    if values.contains(&u) {
        return Ok(None);
    }
    let mut found = Vec::new();
    for c in values {
        let Some(equations) = invert(ev, g, &c, goal)? else {
            return Ok(None);
        };
        for (a, v) in equations {
            let h = ev.eval(&(a - v))?;
            let Some(xs) = roots_at(ev, goal, &h, depth + 1)? else {
                return Ok(None);
            };
            found.extend(xs);
        }
    }
    Ok(Some(found))
}

fn call(head: Symbol, arg: Expr) -> Expr {
    Expr::call(head, [arg])
}

/// The sign of `e`, evaluated, when it is known real.
fn sign(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Ordering>, EvalError> {
    let value = ev.eval(e)?;
    real_sign(ev, &value)
}

/// Whether `c` is known to lie outside `[lo, hi]`, or `(lo, hi)` when
/// `open`.
fn outside(ev: &mut Eval<'_>, c: &Expr, lo: Expr, hi: Expr, open: bool) -> Result<bool, EvalError> {
    let below = sign(ev, &(c.clone() - lo))?;
    let above = sign(ev, &(c.clone() - hi))?;
    Ok(if open {
        matches!(below, Some(Ordering::Less | Ordering::Equal))
            || matches!(above, Some(Ordering::Greater | Ordering::Equal))
    } else {
        below == Some(Ordering::Less) || above == Some(Ordering::Greater)
    })
}

/// The equations `a = v` that `g = c` comes to, `g` a function of the
/// unknown: none when it has no solution, the real ones when only real
/// ones are wanted; `None` when `g` is no function inverted here.
fn invert(
    ev: &mut Eval<'_>,
    g: &Expr,
    c: &Expr,
    goal: &Goal<'_>,
) -> Result<Option<Vec<(Expr, Expr)>>, EvalError> {
    let real = goal.real;
    let pi = || Expr::from(Symbol::PI);
    let half_pi = || Expr::from(Rational::new(Integer::ONE, Integer::from(2))) * pi();
    let Some(head) = g.head().and_then(Expr::as_symbol) else {
        return Ok(None);
    };
    let equations = match (head, g.args()) {
        (Symbol::EXP, [a]) => match log(ev, c, real)? {
            Some(v) => vec![(a.clone(), v)],
            None => Vec::new(),
        },
        (Symbol::LOG, [a]) => vec![(a.clone(), call(Symbol::EXP, c.clone()))],
        (Symbol::SIN, [a]) => {
            if real && outside(ev, c, Expr::from(-1), Expr::from(1), false)? {
                Vec::new()
            } else {
                let asin = call(Symbol::ARC_SIN, c.clone());
                vec![(a.clone(), asin.clone()), (a.clone(), pi() - asin)]
            }
        }
        (Symbol::COS, [a]) => {
            if real && outside(ev, c, Expr::from(-1), Expr::from(1), false)? {
                Vec::new()
            } else {
                let acos = call(Symbol::ARC_COS, c.clone());
                vec![(a.clone(), acos.clone()), (a.clone(), -acos)]
            }
        }
        (Symbol::TAN, [a]) => vec![(a.clone(), call(Symbol::ARC_TAN, c.clone()))],
        (Symbol::ARC_SIN, [a]) => {
            if real && outside(ev, c, -half_pi(), half_pi(), false)? {
                Vec::new()
            } else {
                vec![(a.clone(), call(Symbol::SIN, c.clone()))]
            }
        }
        (Symbol::ARC_COS, [a]) => {
            if real && outside(ev, c, Expr::from(0), pi(), false)? {
                Vec::new()
            } else {
                vec![(a.clone(), call(Symbol::COS, c.clone()))]
            }
        }
        (Symbol::ARC_TAN, [a]) => {
            if real && outside(ev, c, -half_pi(), half_pi(), true)? {
                Vec::new()
            } else {
                vec![(a.clone(), call(Symbol::TAN, c.clone()))]
            }
        }
        (Symbol::ABS, [a]) => {
            if real && sign(ev, c)? == Some(Ordering::Less) {
                Vec::new()
            } else {
                vec![(a.clone(), c.clone()), (a.clone(), -c.clone())]
            }
        }
        (Symbol::POWER, [b, q]) => match power(ev, b, q, c, goal)? {
            Some(equations) => equations,
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    Ok(Some(equations))
}

/// `b^q = c`: an exponential, `b` free of the unknown, or a root of it,
/// `q = ±1/n`; `None` for another power.
fn power(
    ev: &mut Eval<'_>,
    b: &Expr,
    q: &Expr,
    c: &Expr,
    goal: &Goal<'_>,
) -> Result<Option<Vec<(Expr, Expr)>>, EvalError> {
    let real = goal.real;
    if let Some(q) = q.to_rational() {
        // a^(±1/n) = c.
        if !q.numer().unsigned_abs().is_one() || q.is_integer() {
            return Ok(None);
        }
        let n = q.denom();
        let even = !n.unsigned_abs().bit(0);
        let c_sign = sign(ev, c)?;
        if real && even && c_sign == Some(Ordering::Less)
            || q.is_negative() && c_sign == Some(Ordering::Equal)
        {
            return Ok(Some(Vec::new()));
        }
        let inverse = Expr::from(q.checked_recip().expect("not 0"));
        return Ok(Some(vec![(b.clone(), c.clone().pow(inverse))]));
    }
    // b^a = c, b free of the unknown.
    if b.contains(goal.x) {
        return Ok(None);
    }
    let a = q.clone();
    if let (Some(base), Some(value)) = (b.to_rational(), c.to_rational())
        && let Some(k) = exact_log(&base, &value)
    {
        return Ok(Some(vec![(a, Expr::from(k))]));
    }
    Ok(Some(match log(ev, c, real)? {
        Some(v) => vec![(a, v / call(Symbol::LOG, b.clone()))],
        None => Vec::new(),
    }))
}

/// The principal logarithm of `c`, `log|c| + iπ` for a negative `c`;
/// `None` for none: `c = 0`, or `c < 0` when `real`.
fn log(ev: &mut Eval<'_>, c: &Expr, real: bool) -> Result<Option<Expr>, EvalError> {
    Ok(match sign(ev, c)? {
        Some(Ordering::Equal) => None,
        Some(Ordering::Less) if real => None,
        Some(Ordering::Less) => Some(
            call(Symbol::LOG, -c.clone())
                + Expr::call(
                    Symbol::TIMES,
                    [Expr::from(Symbol::I), Expr::from(Symbol::PI)],
                ),
        ),
        _ => Some(call(Symbol::LOG, c.clone())),
    })
}

/// The integer `k` with `b^k = c`, for rationals `b > 0`, `b ≠ 1`, within
/// [`MAX_EXACT_LOG`]: powers of `b`, or of `1/b`, past 1 towards `c` or
/// `1/c`.
fn exact_log(b: &Rational, c: &Rational) -> Option<i64> {
    if !b.is_positive() || *b == Rational::ONE || !c.is_positive() {
        return None;
    }
    let inverse = |q: &Rational| q.checked_recip().expect("not 0");
    let (base, target, sign) = match (*b > Rational::ONE, *c >= Rational::ONE) {
        (true, true) => (b.clone(), c.clone(), 1),
        (true, false) => (b.clone(), inverse(c), -1),
        (false, true) => (inverse(b), c.clone(), -1),
        (false, false) => (inverse(b), inverse(c), 1),
    };
    let mut power = Rational::ONE;
    for k in 0..=MAX_EXACT_LOG {
        match power.cmp(&target) {
            Ordering::Equal => return Some(sign * k),
            Ordering::Greater => return None,
            Ordering::Less => power = &power * &base,
        }
    }
    None
}

/// A root `r^(1/n)` of something in `x`, `n ≥ 2`: `(r, n)`.
fn root_of(g: &Expr, x: &Expr) -> Option<(Expr, u64)> {
    let View::Apply { head, args: [r, q] } = g.view() else {
        return None;
    };
    if head.as_symbol() != Some(Symbol::POWER) || !r.contains(x) {
        return None;
    }
    let q = q.to_rational()?;
    if !q.numer().is_one() {
        return None;
    }
    let n = q.denom().to_i64().and_then(|n| u64::try_from(n).ok())?;
    (n >= 2).then(|| (r.clone(), n))
}

/// A root among the generators `others` isolated and raised to its power,
/// the candidates checked in `f`.
fn radicals(
    ev: &mut Eval<'_>,
    goal: &Goal<'_>,
    f: &Expr,
    p: &Polynomials,
    others: &[usize],
    depth: u32,
) -> Found {
    let Some((j, (r, n))) = others
        .iter()
        .find_map(|&j| root_of(&p.generators()[j], goal.x).map(|root| (j, root)))
    else {
        return Ok(None);
    };
    let k = p.generators().len();
    let parts = match &p.fractions()[0] {
        Fraction::Integers { num, .. } => split(p, &SparseRing::new(Integers, k), num, j, &r, n),
        Fraction::GaussianIntegers { num, .. } => {
            split(p, &SparseRing::new(GaussianIntegers, k), num, j, &r, n)
        }
    };
    let Some((a, b)) = parts else {
        return Ok(None);
    };
    let exponent = i64::try_from(n).expect("a small root");
    let raised = ev.eval(&((-a).pow(exponent) - b.pow(exponent) * r))?;
    let Some(candidates) = roots_at(ev, goal, &raised, depth + 1)? else {
        return Ok(None);
    };
    let ctx = ev.context();
    Ok(Some(
        candidates
            .into_iter()
            .filter(|c| holds(ctx, f, &[(goal.x.clone(), c.clone())]))
            .collect(),
    ))
}

/// The numerator `num` of `p` as `A + B s` in the generator `j`,
/// `s = r^(1/n)`: each power `s^k` is `r^⌊k/n⌋ s^(k mod n)`; `None` when a
/// power of `s` past the first is left.
fn split<R: Numbers>(
    p: &Polynomials,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    j: usize,
    r: &Expr,
    n: u64,
) -> Option<(Expr, Expr)> {
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (k, c) in ring.to_recursive(num, j) {
        let term = p.polynomial_expr(&c) * r.pow(i64::try_from(k / n).ok()?);
        match k % n {
            0 => a.push(term),
            1 => b.push(term),
            _ => return None,
        }
    }
    let sum = |terms: Vec<Expr>| Expr::call(Symbol::PLUS, terms);
    Some((sum(a), sum(b)))
}

/// Whether the values `rules` give the unknowns solve `f = 0`: its value
/// is 0, or a number whose ball holds 0, or unknown for parameters; not a
/// pole.
pub(super) fn holds(ctx: &Context, f: &Expr, rules: &[(Expr, Expr)]) -> bool {
    let Ok(value) = ctx.eval(&f.substitute(rules)) else {
        return true;
    };
    if value == Expr::from(0) {
        return true;
    }
    let infinite = matches!(
        value.as_symbol(),
        Some(Symbol::COMPLEX_INFINITY | Symbol::INDETERMINATE | Symbol::INFINITY)
    ) || value.contains(&Expr::from(Symbol::DIRECTED_INFINITY))
        || value.contains(&Expr::from(Symbol::COMPLEX_INFINITY))
        || value.contains(&Expr::from(Symbol::INDETERMINATE));
    if infinite {
        return false;
    }
    match approximation(ctx, &value, BITS).as_ref().map(Expr::view) {
        Ok(View::Ball(b)) => b.contains_rational(&Rational::ZERO),
        Ok(View::ComplexBall(z)) => z.may_be_zero(),
        _ => true,
    }
}
