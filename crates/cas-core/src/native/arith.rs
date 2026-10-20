//! `Plus`, `Times`, `Power`, `Sqrt` and `DirectedInfinity` (CAS-DOC-4 §6.1).
//!
//! Each receives an application the evaluator has normalised: `Plus` and
//! `Times` are flat and sorted, so like terms and equal bases are adjacent
//! (the canonical order compares products by their factors first, CAS-DOC-3
//! §4.1), and a `Plus` of lists has been threaded. What a native returns is
//! evaluated again, so it may build unevaluated sums and products.

use cas_num::Rational;

use super::linear;
use super::number::{Infinity, Num};
use super::radical::{self, MAX_POWER_BITS};
use crate::assume::{self, Domain, Set, Sign, Truth};
use crate::context::EvalError;
use crate::eval::{Eval, is_list_like, sort_arguments};
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn apply(head: Symbol, args: &[Expr]) -> Expr {
    Expr::apply(&head.into(), args)
}

fn has_head(e: &Expr, s: Symbol) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(s)
}

/// `result`, unless it is `e` itself: no change is no rewrite.
fn changed(e: &Expr, result: Expr) -> Option<Expr> {
    (result != *e).then_some(result)
}

/// `Plus`: numbers added into one, like terms collected (`2 x + x` is
/// `3 x`), zero terms dropped, `Plus[x]` is `x` and `Plus[]` is `0`. An
/// infinite term absorbs the finite ones, `Infinity + x` is `Infinity`;
/// opposite directions, `ComplexInfinity` with another infinity, and
/// `Indeterminate` give `Indeterminate`.
pub(super) fn plus(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if is_normal_sum(e.args()) {
        debug_assert_eq!(sum(e), None, "the shortcut of `Plus` misses a change");
        return Ok(None);
    }
    let result = sum(e);
    warn_indeterminate(ev, e, result.as_ref());
    Ok(result)
}

/// Warns that `e` is an indeterminate form if it gave `Indeterminate` and
/// none of its arguments is.
fn warn_indeterminate(ev: &mut Eval<'_>, e: &Expr, result: Option<&Expr>) {
    let indeterminate = |x: &Expr| x.as_symbol() == Some(Symbol::INDETERMINATE);
    if result.is_some_and(indeterminate) && !e.args().iter().any(indeterminate) {
        ev.warn(Warning::Indeterminate(e.clone()));
    }
}

/// What `Plus` makes of `e`, if anything.
fn sum(e: &Expr) -> Option<Expr> {
    let mut number = Num::zero();
    let mut infinity: Option<Infinity> = None;
    // `(c, m)` for the term `c·m`; like terms are adjacent.
    let mut terms: Vec<(Num, Expr)> = Vec::new();
    for a in e.args() {
        if let Some(n) = Num::of(a) {
            number = number.add(&n);
        } else if let Some(i) = Infinity::of(a) {
            infinity = Some(match infinity {
                Some(j) => j.add(i),
                None => i,
            });
        } else {
            let (c, m) = coefficient(a);
            match terms.last_mut() {
                Some((sum, last)) if *last == m => *sum = sum.add(&c),
                _ => terms.push((c, m)),
            }
        }
    }
    if let Some(i) = infinity {
        return Some(i.to_expr());
    }
    let mut out = Vec::new();
    if !number.is_zero() {
        out.push(number.to_expr());
    }
    out.extend(
        terms
            .into_iter()
            .filter(|(c, _)| !c.is_zero())
            .map(|(c, m)| with_coefficient(&c, m)),
    );
    changed(
        e,
        match out.len() {
            0 => Expr::from(0),
            1 => out.pop().expect("one term"),
            _ => apply(Symbol::PLUS, &out),
        },
    )
}

/// Whether `Plus` leaves these arguments as they are: one nonzero number at
/// most, first, no infinity, and no like terms. Checked without building
/// anything, since most sums are already normal.
fn is_normal_sum(args: &[Expr]) -> bool {
    if args.len() < 2 || args.iter().any(|a| Infinity::of(a).is_some()) {
        return false;
    }
    let numbers = args.iter().filter(|a| a.is_number()).count();
    let first_number =
        numbers == 1 && args[0].is_number() && !Num::of(&args[0]).is_some_and(|n| n.is_zero());
    if numbers > 1 || numbers == 1 && !first_number {
        return false;
    }
    let terms = &args[numbers..];
    terms
        .iter()
        .all(|t| !leading_number(t).is_some_and(|c| c.is_zero()))
        && terms.windows(2).all(|w| monomial(&w[0]) != monomial(&w[1]))
}

/// The number a product starts with, if any.
fn leading_number(term: &Expr) -> Option<Num> {
    if has_head(term, Symbol::TIMES) {
        term.args().first().and_then(Num::of)
    } else {
        None
    }
}

/// The factors of a term after its numeric coefficient: what like terms
/// share. A product without a coefficient is all factors.
fn monomial(term: &Expr) -> &[Expr] {
    if !has_head(term, Symbol::TIMES) {
        return std::slice::from_ref(term);
    }
    match leading_number(term) {
        Some(_) => &term.args()[1..],
        None => term.args(),
    }
}

/// `(c, m)` with `term = c·m`, `c` a number: the first factor of a product
/// that has one, otherwise 1.
fn coefficient(term: &Expr) -> (Num, Expr) {
    if has_head(term, Symbol::TIMES) {
        if let [first, rest @ ..] = term.args() {
            if let Some(c) = Num::of(first) {
                let m = match rest {
                    [m] => m.clone(),
                    _ => apply(Symbol::TIMES, rest),
                };
                return (c, m);
            }
        }
    }
    (Num::one(), term.clone())
}

/// `c·m`, with the factors of a product `m` spliced.
fn with_coefficient(c: &Num, m: Expr) -> Expr {
    if c.is_one() {
        return m;
    }
    let mut factors = vec![c.to_expr()];
    if has_head(&m, Symbol::TIMES) {
        factors.extend_from_slice(m.args());
    } else {
        factors.push(m);
    }
    apply(Symbol::TIMES, &factors)
}

/// `Times`: numbers multiplied into one coefficient, equal bases combined
/// (`x^p x^q` is `x^(p + q)`, so `x x^(-1)` is 1, an accepted unsoundness
/// like O3), powers of numbers with the same exponent multiplied
/// (`2^(1/2) 3^(1/2)` is `6^(1/2)`), `0 x` is 0 (O3), `1 x` is `x`.
/// Lists, in their order, are multiplied as in linear algebra and scaled by
/// the rest (`linear`, CAS-62). Infinite
/// factors: the coefficient's sign turns a direction, `0` times an infinity
/// and `Indeterminate` give `Indeterminate`, and `ComplexInfinity` absorbs
/// the other factors. The other factors of a directed infinity turn it when
/// the assumptions tell their sign: `x Infinity` is `Infinity` for `x > 0`.
pub(super) fn times(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let result = if is_normal_product(e.args()) {
        debug_assert_eq!(
            product(e),
            Ok(None),
            "the shortcut of `Times` misses a change"
        );
        None
    } else {
        product(e)?
    };
    // Most products have no infinity, which their symbols tell in O(1).
    let current = result.as_ref().unwrap_or(e);
    if current.may_contain(Symbol::INFINITY) || current.may_contain(Symbol::DIRECTED_INFINITY) {
        if let Some(turned) = signed_infinity(ev, current) {
            return Ok(Some(turned));
        }
    }
    warn_indeterminate(ev, e, result.as_ref());
    Ok(result)
}

/// A directed infinity times factors of known sign, as the infinity turned
/// by that sign.
fn signed_infinity(ev: &Eval<'_>, e: &Expr) -> Option<Expr> {
    if e.head().and_then(Expr::as_symbol) != Some(Symbol::TIMES) {
        return None;
    }
    let mut direction = None;
    let mut rest = Vec::new();
    for a in e.args() {
        match Infinity::of(a) {
            Some(Infinity::Directed(d)) if direction.is_none() => direction = Some(d),
            Some(_) => return None,
            None if Num::of(a).is_some() => return None,
            None => rest.push(a.clone()),
        }
    }
    let direction = direction?;
    let rest = match rest.as_slice() {
        [] => return None,
        [one] => one.clone(),
        _ => apply(Symbol::TIMES, &rest),
    };
    let turned = match assume::sign_of(ev.context(), &rest) {
        Sign::Positive => direction,
        Sign::Negative => direction.mul(&Num::Real(Rational::from(-1))),
        _ => return None,
    };
    Some(Infinity::Directed(turned).to_expr())
}

/// What `Times` makes of `e`, if anything.
///
/// # Errors
///
/// [`EvalError::Shape`] for lists whose dimensions do not fit.
fn product(e: &Expr) -> Result<Option<Expr>, EvalError> {
    let mut coefficient = Num::one();
    let mut infinity: Option<Infinity> = None;
    // `(base, exponents, the factor itself)`; equal bases are adjacent.
    let mut powers: Vec<(Expr, Vec<Expr>, Expr)> = Vec::new();
    let mut lists = Vec::new();
    for a in e.args() {
        if let Some(n) = Num::of(a) {
            coefficient = coefficient.mul(&n);
        } else if let Some(i) = Infinity::of(a) {
            infinity = Some(match infinity {
                Some(j) => j.mul(i),
                None => i,
            });
        } else if is_list_like(a) {
            lists.push(a.clone());
        } else {
            let (base, exponent) = match a.view() {
                View::Apply { head, args: [b, n] } if head.as_symbol() == Some(Symbol::POWER) => {
                    (b.clone(), n.clone())
                }
                _ => (a.clone(), Expr::from(1)),
            };
            match powers.last_mut() {
                Some((last, exponents, _)) if *last == base => exponents.push(exponent),
                _ => powers.push((base, vec![exponent], a.clone())),
            }
        }
    }
    if let Some(i) = infinity {
        let finite_only = powers.is_empty() && lists.is_empty();
        return Ok(changed(e, times_infinity(e, i, &coefficient, finite_only)));
    }
    if coefficient.is_zero() && lists.is_empty() {
        return Ok(changed(e, Expr::from(0)));
    }
    let mut factors = Vec::new();
    for (base, exponents, factor) in powers {
        factors.push(match exponents.as_slice() {
            [_] => factor,
            _ => apply(Symbol::POWER, &[base, apply(Symbol::PLUS, &exponents)]),
        });
    }
    let factors = merge_number_powers(factors);
    let mut out = Vec::with_capacity(factors.len() + lists.len() + 1);
    if !coefficient.is_one() || factors.is_empty() && lists.is_empty() {
        out.push(coefficient.to_expr());
    }
    out.extend(factors);
    // Vectors and matrices, from the left, scaled by the rest.
    if lists.iter().all(linear::is_list) {
        if let Some(folded) = linear::fold(&lists)? {
            let scalar = match out.len() {
                0 => None,
                1 => out.pop(),
                _ => Some(apply(Symbol::TIMES, &out)),
            };
            let result = match (scalar, linear::shape(&folded)) {
                (None, _) => folded,
                (Some(c), Some(_)) => linear::scale(&c, &folded),
                (Some(c), None) => apply(Symbol::TIMES, &[c, folded]),
            };
            return Ok(changed(e, result));
        }
    }
    out.extend(lists);
    Ok(changed(
        e,
        match out.len() {
            1 => out.pop().expect("one factor"),
            _ => apply(Symbol::TIMES, &out),
        },
    ))
}

/// Whether `Times` leaves these arguments as they are: one number at most,
/// first, neither 0 nor 1, no infinity, no equal bases, and no two powers of
/// numbers with one exponent.
fn is_normal_product(args: &[Expr]) -> bool {
    if args.len() < 2
        || args
            .iter()
            .any(|a| Infinity::of(a).is_some() || linear::is_list(a))
    {
        return false;
    }
    let numbers = args.iter().filter(|a| a.is_number()).count();
    if numbers > 1 || numbers == 1 && !args[0].is_number() {
        return false;
    }
    if numbers == 1 && Num::of(&args[0]).is_some_and(|n| n.is_zero() || n.is_one()) {
        return false;
    }
    let factors: Vec<&Expr> = args[numbers..]
        .iter()
        .filter(|a| !is_list_like(a))
        .collect();
    let base = |f: &Expr| -> Expr {
        match f.view() {
            View::Apply { head, args: [b, _] } if head.as_symbol() == Some(Symbol::POWER) => {
                b.clone()
            }
            _ => f.clone(),
        }
    };
    if factors.windows(2).any(|w| base(w[0]) == base(w[1])) {
        return false;
    }
    let mut exponents = Vec::new();
    for f in factors {
        if let View::Apply { head, args: [b, r] } = f.view() {
            if head.as_symbol() == Some(Symbol::POWER)
                && b.to_rational().is_some_and(|b| b.is_positive())
                && r.to_rational().is_some_and(|r| !r.is_integer())
            {
                if exponents.contains(&r) {
                    return false;
                }
                exponents.push(r);
            }
        }
    }
    true
}

/// A product with an infinite factor `i`, a numeric coefficient, and
/// `finite_only` when nothing else is left.
fn times_infinity(e: &Expr, i: Infinity, coefficient: &Num, finite_only: bool) -> Expr {
    match i {
        Infinity::Indeterminate => i.to_expr(),
        _ if coefficient.is_zero() => Infinity::Indeterminate.to_expr(),
        Infinity::Complex => i.to_expr(),
        Infinity::Directed(d) => {
            // The coefficient turns the direction when its own direction is
            // exact; symbolic factors are kept, their sign being unknown.
            let Some(unit) = coefficient.unit() else {
                return e.clone();
            };
            let infinity = Infinity::Directed(d.mul(&unit)).to_expr();
            if finite_only {
                infinity
            } else {
                let mut out = vec![infinity];
                out.extend(
                    e.args()
                        .iter()
                        .filter(|a| Num::of(a).is_none() && Infinity::of(a).is_none())
                        .cloned(),
                );
                // In canonical order: otherwise the next round sorts it back
                // and this undoes it again, without end.
                let out = sort_arguments(&Symbol::TIMES.into(), &out).unwrap_or(out);
                apply(Symbol::TIMES, &out)
            }
        }
    }
}

/// `Power[a, r] Power[b, r]` as `Power[a b, r]`, for positive rationals `a`
/// and `b` and a rational exponent `r`; the new power is evaluated, which
/// takes out what the product makes a perfect power.
fn merge_number_powers(factors: Vec<Expr>) -> Vec<Expr> {
    let number_power = |f: &Expr| -> Option<(Rational, Rational)> {
        let View::Apply { head, args: [b, r] } = f.view() else {
            return None;
        };
        if head.as_symbol() != Some(Symbol::POWER) {
            return None;
        }
        let (b, r) = (b.to_rational()?, r.to_rational()?);
        (b.is_positive() && !r.is_integer()).then_some((b, r))
    };
    let mut merged: Vec<(Rational, Rational)> = Vec::new();
    let mut seen = 0;
    for f in &factors {
        if let Some((b, r)) = number_power(f) {
            seen += 1;
            match merged.iter_mut().find(|(_, s)| *s == r) {
                Some((product, _)) => *product = &*product * &b,
                None => merged.push((b, r)),
            }
        }
    }
    if merged.len() == seen {
        return factors;
    }
    let mut out: Vec<Expr> = merged
        .into_iter()
        .map(|(b, r)| apply(Symbol::POWER, &[Expr::from(b), Expr::from(r)]))
        .collect();
    out.extend(factors.into_iter().filter(|f| number_power(f).is_none()));
    out
}

/// `Power`:
///
/// - numbers to integer powers, exactly; `0^n` for `n < 0` is
///   `ComplexInfinity`, `0^0` is `Indeterminate`; a result estimated above
///   2²⁴ bits stays unevaluated;
/// - rational powers of rationals in canonical form ([`radical`]);
/// - `x^0` is 1, `x^1` is `x`, `1^x` is 1;
/// - `(b^p)^n` is `b^(p n)` and `(x y)^n` is `x^n y^n` for an integer `n`,
///   the identities that hold for every complex `b`; for a positive
///   rational `b`, `(b^p)^q` is `b^(p q)` for any real `q`;
/// - infinities: a direction to a positive integer power, 0 for a negative
///   power, `Indeterminate` for `Infinity^0`.
pub(super) fn power(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [b, n] = e.args() else {
        return Ok(None);
    };
    if matches!(Infinity::of(b), Some(Infinity::Indeterminate))
        || matches!(Infinity::of(n), Some(Infinity::Indeterminate))
    {
        return Ok(Some(Infinity::Indeterminate.to_expr()));
    }
    // A matrix power, or nothing for a vector: never `x^0` (CAS-62). Apart
    // from the rest, whose frame would stay on the stack while the products
    // evaluate.
    if linear::is_list(b) {
        return linear::power(ev, b, n);
    }
    scalar_power(ev, b, n)
}

/// [`power`] of a base that is not a list, neither side `Indeterminate`.
#[allow(
    clippy::many_single_char_names,
    reason = "b^n, as written in mathematics"
)]
fn scalar_power(ev: &mut Eval<'_>, b: &Expr, n: &Expr) -> Native {
    let real_domain = ev.context().domain() == Domain::Real;
    let indeterminate = || Infinity::Indeterminate.to_expr();
    if let (Some(base), Some(exponent)) = (Num::of(b), Num::of(n)) {
        let result = number_power(&base, &exponent, real_domain);
        // `1/0` and `0^0`.
        if base.is_zero() {
            let power = || apply(Symbol::POWER, &[b.clone(), n.clone()]);
            match result.as_ref().and_then(Expr::as_symbol) {
                Some(Symbol::COMPLEX_INFINITY) => ev.warn(Warning::Infinite(power())),
                Some(Symbol::INDETERMINATE) => ev.warn(Warning::Indeterminate(power())),
                _ => {}
            }
        }
        return Ok(result);
    }
    let exponent = Num::of(n);
    if exponent.as_ref().is_some_and(Num::is_zero) {
        if Infinity::of(b).is_some() {
            ev.warn(Warning::Indeterminate(apply(
                Symbol::POWER,
                &[b.clone(), n.clone()],
            )));
            return Ok(Some(indeterminate()));
        }
        return Ok(Some(Expr::from(1)));
    }
    if exponent.as_ref().is_some_and(Num::is_one) {
        return Ok(Some(b.clone()));
    }
    if Num::of(b).is_some_and(|b| b.is_one()) {
        return Ok(Some(Expr::from(1)));
    }
    let integer = exponent
        .as_ref()
        .and_then(|n| n.real())
        .filter(|n| n.is_integer());
    if let Some(i) = Infinity::of(b) {
        let Some(n) = exponent.as_ref().and_then(Num::real) else {
            return Ok(None);
        };
        return Ok(if n.is_negative() {
            Some(Expr::from(0))
        } else {
            match (i, n.is_integer()) {
                (Infinity::Directed(d), true) => {
                    let k = n.numer().to_i64();
                    k.and_then(|k| d.pow(k))
                        .map(|d| Infinity::Directed(d).to_expr())
                }
                (Infinity::Directed(d), false) if d.is_one() => Some(Symbol::INFINITY.into()),
                (Infinity::Complex, _) => Some(Symbol::COMPLEX_INFINITY.into()),
                _ => None,
            }
        });
    }
    if let View::Apply { head, args } = b.view() {
        match (head.as_symbol(), args) {
            // `(b^p)^n`.
            (Some(Symbol::POWER), [inner, p]) => {
                let positive = inner.to_rational().is_some_and(|q| q.is_positive())
                    && exponent.as_ref().is_some_and(|n| n.real().is_some())
                    && p.to_rational().is_some();
                if integer.is_some() || positive {
                    let product = apply(Symbol::TIMES, &[p.clone(), n.clone()]);
                    return Ok(Some(apply(Symbol::POWER, &[inner.clone(), product])));
                }
                if let Some(rewritten) = real_power_of_power(ev, inner, p, n, real_domain) {
                    return Ok(Some(rewritten));
                }
            }
            // `Abs[x]^n` is `x^n` for a real `x` and an even `n`.
            (Some(Symbol::ABS), [inner])
                if integer.as_ref().is_some_and(|n| {
                    n.numer()
                        .unsigned_abs()
                        .trailing_zeros()
                        .is_some_and(|t| t > 0)
                }) && assume::element(ev.context(), inner, Set::Reals) == Truth::True =>
            {
                return Ok(Some(apply(Symbol::POWER, &[inner.clone(), n.clone()])));
            }
            // `(x y)^n`.
            (Some(Symbol::TIMES), factors) if integer.is_some() => {
                let powers: Vec<Expr> = factors
                    .iter()
                    .map(|f| apply(Symbol::POWER, &[f.clone(), n.clone()]))
                    .collect();
                return Ok(Some(apply(Symbol::TIMES, &powers)));
            }
            _ => {}
        }
    }
    Ok(None)
}

/// `(b^p)^q` for rationals `p` and `q`, when the assumptions allow
/// (CAS-DOC-4 §6.1, §7):
/// - `b^(p q)` for `b >= 0`;
/// - in the Real domain, for a real `b`: `b^(p q)` for an integer `p` and an
///   odd denominator of `q`, real odd roots being multiplicative (O7), and
///   `Abs[b]^(p q)` for an even `p`, since `b^p` is `Abs[b]^p`.
fn real_power_of_power(
    ev: &Eval<'_>,
    b: &Expr,
    p: &Expr,
    q: &Expr,
    real_domain: bool,
) -> Option<Expr> {
    let (pr, qr) = (p.to_rational()?, q.to_rational()?);
    let ctx = ev.context();
    let product = apply(Symbol::TIMES, &[p.clone(), q.clone()]);
    if matches!(
        assume::sign_of(ctx, b),
        Sign::Positive | Sign::NonNegative | Sign::Zero
    ) {
        return Some(apply(Symbol::POWER, &[b.clone(), product]));
    }
    if !real_domain || !pr.is_integer() || assume::element(ctx, b, Set::Reals) != Truth::True {
        return None;
    }
    let even = |k: &cas_num::Integer| k.unsigned_abs().trailing_zeros().is_none_or(|t| t > 0);
    if !even(qr.denom()) {
        return Some(apply(Symbol::POWER, &[b.clone(), product]));
    }
    even(pr.numer()).then(|| {
        let abs = apply(Symbol::ABS, std::slice::from_ref(b));
        apply(Symbol::POWER, &[abs, product])
    })
}

/// `Abs[x]`: of a number, its absolute value, the square root of its norm
/// for a complex one; of an infinity, `Infinity`; `Abs[Abs[x]]` is
/// `Abs[x]`; a numeric coefficient comes out, `Abs[-2 x]` is `2 Abs[x]`; and
/// `x` or `-x` when the assumptions tell the sign of `x`.
pub(super) fn abs(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if let Some(n) = Num::of(x) {
        return Ok(Some(match n {
            Num::Real(q) => Expr::from(q.abs()),
            Num::Complex(z) => {
                let half = Expr::from(Rational::new(1.into(), 2.into()));
                apply(Symbol::POWER, &[Expr::from(z.norm()), half])
            }
        }));
    }
    if matches!(
        Infinity::of(x),
        Some(Infinity::Directed(_) | Infinity::Complex)
    ) {
        return Ok(Some(Symbol::INFINITY.into()));
    }
    if let View::Apply { head, args } = x.view() {
        match (head.as_symbol(), args) {
            (Some(Symbol::ABS), [_]) => return Ok(Some(x.clone())),
            (Some(Symbol::TIMES), [c, rest @ ..]) if Num::of(c).is_some() => {
                let rest = match rest {
                    [one] => one.clone(),
                    _ => apply(Symbol::TIMES, rest),
                };
                let magnitude = apply(Symbol::ABS, std::slice::from_ref(c));
                let inner = apply(Symbol::ABS, &[rest]);
                return Ok(Some(apply(Symbol::TIMES, &[magnitude, inner])));
            }
            _ => {}
        }
    }
    Ok(match assume::sign_of(ev.context(), x) {
        Sign::Positive | Sign::NonNegative | Sign::Zero => Some(x.clone()),
        Sign::Negative | Sign::NonPositive => {
            Some(apply(Symbol::TIMES, &[Expr::from(-1), x.clone()]))
        }
        Sign::Nonzero | Sign::Unknown => None,
    })
}

/// `base^exponent` for two numbers, or `None` to leave it. In the Real
/// domain, an odd root of a negative is the real root (O7).
fn number_power(base: &Num, exponent: &Num, real_domain: bool) -> Option<Expr> {
    let r = exponent.real()?;
    if base.is_zero() {
        return Some(if r.is_zero() {
            Infinity::Indeterminate.to_expr()
        } else if r.is_negative() {
            Infinity::Complex.to_expr()
        } else {
            Expr::from(0)
        });
    }
    if r.is_integer() {
        let n = r.numer().to_i64()?;
        if base.power_bits(n) > MAX_POWER_BITS {
            return None;
        }
        return base.pow(n).map(|v| v.to_expr());
    }
    let result = radical::power(base.real()?, r, real_domain)?;
    let original = Expr::apply(&Symbol::POWER.into(), &[base.to_expr(), exponent.to_expr()]);
    (result != original).then_some(result)
}

/// `Sqrt[x]` is `Power[x, 1/2]`.
pub(super) fn sqrt(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let half = Expr::from(Rational::new(1.into(), 2.into()));
    Ok(Some(apply(Symbol::POWER, &[x.clone(), half])))
}

/// `DirectedInfinity[d]` with its direction of modulus 1: `Infinity` for
/// the direction 1, `ComplexInfinity` for none or 0. A complex direction
/// whose modulus is not rational stays as it is.
pub(super) fn directed_infinity(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [] => Ok(Some(Symbol::COMPLEX_INFINITY.into())),
        [d] => {
            let Some(d) = Num::of(d) else {
                return Ok(None);
            };
            if d.is_zero() {
                return Ok(Some(Symbol::COMPLEX_INFINITY.into()));
            }
            match d.unit() {
                Some(unit) => Ok(changed(e, Infinity::Directed(unit).to_expr())),
                None => Ok(None),
            }
        }
        _ => Ok(None),
    }
}
