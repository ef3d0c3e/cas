//! `PolyGamma` (CAS-DOC-31 §3.2, CAS-355), Xcas's `Psi(x)` and `Psi(x,
//! m)`: the digamma function `ψ = Γ′/Γ` and its `m`-th derivative, in
//! giac's argument order.
//!
//! - `PolyGamma[x, 0]` is `PolyGamma[x]`; an order that is no integer `≥ 0`
//!   leaves the call.
//! - At a rational `x = r + n`, `r` in `(0, 1]` and `n` an integer, the
//!   value at `r` shifted by the recurrence `ψ^(m)(x + 1) = ψ^(m)(x) +
//!   (−1)^m m!/x^(m+1)`, at most [`MAX_SHIFT`] steps and within
//!   `Limits::number_bits`:
//!   - `ψ(1) = −γ` and `ψ^(m)(1) = (−1)^(m+1) m! ζ(m + 1)`, `ζ` through
//!     `Zeta`: powers of `π` for odd `m`, `Zeta[m + 1]` left for even `m`
//!     (giac's `Psi(2,1)` is `pi^2/6 - 1`);
//!   - `ψ(1/2) = −γ − 2 log 2` and `ψ^(m)(1/2) = (−1)^(m+1) m! (2^(m+1) − 1)
//!     ζ(m + 1)`;
//!   - for `m = 0` and a denominator `q ≤` [`MAX_GAUSS`], Gauss's digamma
//!     theorem, `ψ(p/q) = −γ − log(2q) − (π/2) cot(πp/q) + 2 Σ_{k=1}^{⌊(q −
//!     1)/2⌋} cos(2πkp/q) log sin(πk/q)`, the trigonometric values whatever
//!     evaluation makes of them (`ψ(1/3) = −γ − π/(2√3) − (3/2) log 3`).
//! - At the poles, the integers `≤ 0`: `ComplexInfinity` with the warning
//!   of `1/0`.
//! - `ψ(∞) = ∞` and `ψ^(m)(∞) = 0`.
//! - A real ball by `Ball::polygamma`, `Indeterminate` with a warning where
//!   it meets a pole.
//! - Anything else stays: `Psi(n)` and `Psi(x+1)`, as giac's.

#![allow(
    clippy::many_single_char_names,
    reason = "x = r + n and the order m, as written in the comments"
)]

use cas_num::{Integer, Natural, Rational};

use super::number::{Infinity, Num, Precision};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The most steps of the recurrence from `r` in `(0, 1]` to `x`.
const MAX_SHIFT: u64 = 4096;

/// The largest denominator of Gauss's digamma theorem.
const MAX_GAUSS: i64 = 12;

/// The highest order of a ball's polygamma, `Ball::polygamma`'s.
const MAX_ORDER: u64 = 256;

/// `PolyGamma[x]` and `PolyGamma[x, m]`.
pub(super) fn polygamma(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, m) = match e.args() {
        [x] => (x, 0),
        [x, n] => {
            if n.as_i64() == Some(0) {
                return Ok(Some(Expr::call(Symbol::POLY_GAMMA, [x.clone()])));
            }
            let Some(m) = n.as_i64().and_then(|m| u64::try_from(m).ok()) else {
                return Ok(None);
            };
            (x, m)
        }
        _ => return Ok(None),
    };
    if let Some(Infinity::Directed(d)) = Infinity::of(x)
        && d.is_one()
    {
        return Ok(Some(if m == 0 {
            Expr::from(Symbol::INFINITY)
        } else {
            Expr::from(0)
        }));
    }
    if let Some(q) = x.to_rational() {
        return at_rational(ev, e, &q, m);
    }
    if let Some(ref number @ Num::Ball(ref ball)) = Num::of(x) {
        if m > MAX_ORDER {
            return Ok(None);
        }
        let bits = number.precision(Precision::of(ev.context()));
        let value = ball.polygamma(m, bits + 16).round_to_accuracy(bits);
        if !value.is_finite() {
            ev.warn(Warning::Indeterminate(e.clone()));
            return Ok(Some(Symbol::INDETERMINATE.into()));
        }
        return Ok(Some(Num::Ball(Box::new(value)).to_expr()));
    }
    Ok(None)
}

/// `(−1)^k c`.
fn signed(k: u64, c: Expr) -> Expr {
    if k % 2 == 0 { c } else { -c }
}

/// `ψ^(m)` at a rational.
fn at_rational(ev: &mut Eval<'_>, e: &Expr, q: &Rational, m: u64) -> Native {
    if q.is_integer() && !q.is_positive() {
        ev.warn(Warning::Infinite(e.clone()));
        return Ok(Some(Symbol::COMPLEX_INFINITY.into()));
    }
    // x = r + n, r in (0, 1].
    let (r, n) = if q.is_integer() {
        (Rational::ONE, q.numer() - &Integer::ONE)
    } else {
        let floor = q.floor();
        (q - &Rational::from(floor.clone()), floor)
    };
    let Some(n) = n.to_i64() else {
        return Ok(None);
    };
    if n.unsigned_abs() > MAX_SHIFT {
        return Ok(None);
    }
    let Some(base) = at_unit(&r, m) else {
        return Ok(None);
    };
    // Σ 1/(r + j)^(m+1) for j < n, or −Σ 1/(r − j)^(m+1) for 1 ≤ j ≤ −n.
    let power = i64::try_from(m + 1).expect("a small order");
    let max_bits = ev.context().limits().number_bits;
    let mut sum = Rational::ZERO;
    let steps: Vec<i64> = if n >= 0 {
        (0..n).collect()
    } else {
        (1..=-n).map(|j| -j).collect()
    };
    for j in steps {
        let at = &r + &Rational::from(Integer::from(j));
        sum = &sum + &at.pow(-power);
        if sum.numer().bit_length() + sum.denom().bit_length() > max_bits {
            return Ok(None);
        }
    }
    if n < 0 {
        sum = -sum;
    }
    let factorial = Integer::from(Natural::factorial(m));
    let shift = signed(m, Expr::from(&sum * &Rational::from(factorial)));
    let value = ev.eval(&(base + shift))?;
    ev.step("Psi", e, &value);
    Ok(Some(value))
}

/// `ψ^(m)(r)` for `r` in `(0, 1]`, where it is known.
fn at_unit(r: &Rational, m: u64) -> Option<Expr> {
    let gamma = Expr::from(Symbol::EULER_GAMMA);
    let factorial = Expr::from(Integer::from(Natural::factorial(m)));
    let zeta = || {
        Expr::call(
            Symbol::ZETA,
            [Expr::from(i64::try_from(m + 1).expect("small"))],
        )
    };
    if r.is_integer() {
        return Some(if m == 0 {
            -gamma
        } else {
            signed(m + 1, factorial * zeta())
        });
    }
    let p = r.numer().to_i64()?;
    let q = r.denom().to_i64()?;
    if q == 2 {
        let log2 = Expr::call(Symbol::LOG, [Expr::from(2)]);
        return Some(if m == 0 {
            -gamma - Expr::from(2) * log2
        } else {
            let scale = Expr::from(Integer::from(2).pow(m + 1) - Integer::ONE);
            signed(m + 1, factorial * scale * zeta())
        });
    }
    if m != 0 || q > MAX_GAUSS {
        return None;
    }
    Some(gauss(p, q))
}

/// Gauss's digamma theorem at `p/q`, `0 < p < q`.
fn gauss(p: i64, q: i64) -> Expr {
    let pi = || Expr::from(Symbol::PI);
    let fraction = |a: i64, b: i64| Expr::from(Rational::new(Integer::from(a), Integer::from(b)));
    let mut terms = vec![
        -Expr::from(Symbol::EULER_GAMMA),
        -Expr::call(Symbol::LOG, [Expr::from(2 * q)]),
        -(pi() / 2) * Expr::call(Symbol::COT, [fraction(p, q) * pi()]),
    ];
    for k in 1..=(q - 1) / 2 {
        let cosine = Expr::call(Symbol::COS, [fraction(2 * k * p, q) * pi()]);
        let sine = Expr::call(Symbol::SIN, [fraction(k, q) * pi()]);
        terms.push(Expr::from(2) * cosine * Expr::call(Symbol::LOG, [sine]));
    }
    terms.into_iter().sum()
}

#[cfg(test)]
mod tests;
