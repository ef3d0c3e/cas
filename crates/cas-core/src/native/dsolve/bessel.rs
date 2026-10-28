//! Bessel's and Airy's equations (CAS-DOC-34 §6, CAS-385), homogeneous and
//! of order 2, `a₂ y″ + a₁ y′ + a₀ y = 0`, read through `p = a₁/a₂` and `q
//! = a₀/a₂`:
//!
//! - **Bessel's**, `p = 1/x` and `q x² = A x² − n²`, `n` a natural number
//!   and the sign of `A` known: `x² y″ + x y′ + (a² x² − n²) y = 0` up to a
//!   factor, `c_0 J_n(a x) + c_1 Y_n(a x)` with `a = √A` for `A > 0`, and
//!   the modified `c_0 I_n(a x) + c_1 K_n(a x)` with `a = √(−A)` for `A <
//!   0`.
//! - **Airy's**, `p = 0` and `q = −k x`, `k` a number not 0: `y″ = k x y`,
//!   `c_0 Ai(c x) + c_1 Bi(c x)` with `c³ = k`, `c` real.
//!
//! giac solves none of these. Other orders, a right side, or a symbolic `A`
//! or `k` whose sign is not known leave the equation to what follows.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the equation: a₂ y″ + a₁ y′ + a₀ y, p, q, n, k"
)]

use std::cmp::Ordering;

use super::super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The general solution of Bessel's or Airy's equation `Σ aₖ y⁽ᵏ⁾ = f`;
/// `None` when it is neither.
pub(super) fn special(
    ev: &mut Eval<'_>,
    a: &[Expr],
    f: &Expr,
    x: &Expr,
    constants: &[Expr],
) -> Result<Option<Expr>, EvalError> {
    let ([a0, a1, a2], [c0, c1, ..]) = (a, constants) else {
        return Ok(None);
    };
    if *f != Expr::from(0) {
        return Ok(None);
    }
    let together = |ev: &mut Eval<'_>, e: Expr| ev.eval(&Expr::call(Symbol::TOGETHER, [e]));
    let p = together(ev, a1.clone() / a2.clone())?;
    let q = together(ev, a0.clone() / a2.clone())?;
    let general = |ev: &mut Eval<'_>, first: Expr, second: Expr| {
        ev.eval(&(c0.clone() * first + c1.clone() * second))
    };
    // Airy's: p = 0, q = −k x.
    if p == Expr::from(0) {
        let k = together(ev, -(q.clone() / x.clone()))?;
        if k.contains(x) || k == Expr::from(0) {
            return Ok(None);
        }
        let c = match real_sign(ev.context(), &k) {
            Some(Ordering::Greater) => Expr::call(Symbol::POWER, [k, Expr::from(1) / 3]),
            Some(Ordering::Less) => -Expr::call(Symbol::POWER, [-k, Expr::from(1) / 3]),
            _ => return Ok(None),
        };
        let u = ev.eval(&(c * x.clone()))?;
        let ai = Expr::call(Symbol::AIRY_AI, [u.clone()]);
        let bi = Expr::call(Symbol::AIRY_BI, [u]);
        return general(ev, ai, bi).map(Some);
    }
    // Bessel's: p x = 1, q x² = A x² − n².
    if together(ev, p * x.clone())? != Expr::from(1) {
        return Ok(None);
    }
    let w = ev.eval(&Expr::call(Symbol::EXPAND, [q * x.clone().pow(2)]))?;
    let coefficient = |ev: &mut Eval<'_>, k: i64| {
        ev.eval(&Expr::call(
            Symbol::COEFFICIENT,
            [w.clone(), x.clone(), Expr::from(k)],
        ))
    };
    let (big_a, linear, constant) = (
        coefficient(ev, 2)?,
        coefficient(ev, 1)?,
        coefficient(ev, 0)?,
    );
    let rebuilt = ev.eval(&(big_a.clone() * x.clone().pow(2) + constant.clone() - w.clone()))?;
    if linear != Expr::from(0) || rebuilt != Expr::from(0) || big_a.contains(x) {
        return Ok(None);
    }
    // n² = −constant, n a natural number.
    let Some(n_square) = (-constant).to_integer() else {
        return Ok(None);
    };
    if n_square.is_negative() {
        return Ok(None);
    }
    let (n, rest) = n_square.sqrt_rem();
    if !rest.is_zero() {
        return Ok(None);
    }
    let (first, second, scale) = match real_sign(ev.context(), &big_a) {
        Some(Ordering::Greater) => (Symbol::BESSEL_J, Symbol::BESSEL_Y, big_a),
        Some(Ordering::Less) => (Symbol::BESSEL_I, Symbol::BESSEL_K, -big_a),
        _ => return Ok(None),
    };
    let u = ev.eval(&(Expr::call(Symbol::SQRT, [scale]) * x.clone()))?;
    let n = Expr::from(n);
    let y1 = Expr::call(first, [n.clone(), u.clone()]);
    let y2 = Expr::call(second, [n, u]);
    general(ev, y1, y2).map(Some)
}
