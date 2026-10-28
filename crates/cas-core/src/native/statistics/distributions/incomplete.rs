//! The regularized incomplete gamma and beta functions at integer and
//! half-integer parameters (CAS-DOC-23 §4, CAS-304), in closed form, for
//! the cumulative functions of the χ², Student and Fisher distributions.
//!
//! - `P(a, x) = γ(a, x)/Γ(a)`: `1 − e^(−x) Σ_{j<a} xʲ/j!` for a natural
//!   `a`, and from `P(1/2, x) = erf(√x)` by `P(a + 1, x) = P(a, x) − xᵃ
//!   e^(−x)/Γ(a + 1)` for a half-integer.
//! - `I_x(a, b)`: for naturals the binomial sum `Σ_{j=a}^{a+b−1} C(a+b−1, j)
//!   xʲ (1 − x)^(a+b−1−j)`. Otherwise it starts from `I_x(1/2, 1/2) =
//!   (2/π) asin √x`, `I_x(1/2, 1) = √x` or `I_x(1, 1/2) = 1 − √(1 − x)`, by
//!   the parities of `a` and `b`, and steps up by `I_x(a + 1, b) = I_x(a, b) −
//!   xᵃ(1−x)ᵇ/(a B(a, b))` and `I_x(a, b + 1) = I_x(a, b) + xᵃ(1−x)ᵇ/(b
//!   B(a, b))` (DLMF 8.17.20, 8.17.21), with `B` exact through `Gamma`.
//!
//! The parameters come as twice their values, natural numbers; past
//! [`MAX_TERMS`] terms there is no closed form here. The expressions are
//! built unevaluated, for the caller to evaluate.

use cas_num::{Integer, Rational};

use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most terms of a closed form.
const MAX_TERMS: u64 = 1 << 12;

/// `k/2`.
fn half(k: u64) -> Expr {
    Expr::from(Rational::new(Integer::from(k), Integer::from(2)))
}

/// `b^(k/2)`, 1 for `k = 0`.
fn half_power(b: &Expr, k: u64) -> Expr {
    match k {
        0 => Expr::from(1),
        2 => b.clone(),
        _ => b.pow(half(k)),
    }
}

/// `P(a, x)` for `a = twice/2 > 0`; `None` past [`MAX_TERMS`].
pub(super) fn lower_gamma(twice: u64, x: &Expr) -> Option<Expr> {
    if twice == 0 || twice / 2 > MAX_TERMS {
        return None;
    }
    let decay = Expr::call(Symbol::EXP, [-x]);
    if twice % 2 == 0 {
        let terms: Vec<Expr> = (0..twice / 2)
            .map(|j| {
                half_power(x, 2 * j) / Expr::call(Symbol::FACTORIAL, [Expr::from(Integer::from(j))])
            })
            .collect();
        return Some(1 - decay * Expr::call(Symbol::PLUS, terms));
    }
    // a = 1/2 + m.
    let terms: Vec<Expr> = (0..twice / 2)
        .map(|j| half_power(x, 2 * j + 1) / Expr::call(Symbol::GAMMA, [half(2 * j + 3)]))
        .collect();
    let erf = Expr::call(Symbol::ERF, [Expr::call(Symbol::SQRT, [x.clone()])]);
    Some(erf - decay * Expr::call(Symbol::PLUS, terms))
}

/// `I_x(a, b)` for `a = a2/2 > 0` and `b = b2/2 > 0`; `None` past
/// [`MAX_TERMS`].
pub(super) fn beta(a2: u64, b2: u64, x: &Expr) -> Option<Expr> {
    if a2 == 0 || b2 == 0 || a2.saturating_add(b2) > 2 * MAX_TERMS {
        return None;
    }
    let rest = 1 - x;
    if a2 % 2 == 0 && b2 % 2 == 0 {
        let (a, b) = (a2 / 2, b2 / 2);
        let n = a + b - 1;
        let terms: Vec<Expr> = (a..=n)
            .map(|j| {
                let coefficient = Expr::call(
                    Symbol::BINOMIAL,
                    [Expr::from(Integer::from(n)), Expr::from(Integer::from(j))],
                );
                coefficient * half_power(x, 2 * j) * half_power(&rest, 2 * (n - j))
            })
            .collect();
        return Some(Expr::call(Symbol::PLUS, terms));
    }
    // The base case of the parities, then the steps up, a first.
    let (mut a, mut b) = (2 - a2 % 2, 2 - b2 % 2);
    let root = Expr::call(Symbol::SQRT, [x.clone()]);
    let mut value = match (a, b) {
        (1, 1) => 2 * Expr::call(Symbol::ARC_SIN, [root]) / Expr::from(Symbol::PI),
        (1, _) => root,
        (_, 1) => 1 - Expr::call(Symbol::SQRT, [rest.clone()]),
        _ => x.clone(),
    };
    let step = |a: u64, b: u64| {
        let beta = Expr::call(Symbol::BETA, [half(a), half(b)]);
        half_power(x, a) * half_power(&rest, b) / beta
    };
    while a < a2 {
        value -= step(a, b) / half(a);
        a += 2;
    }
    while b < b2 {
        value += step(a, b) / half(b);
        b += 2;
    }
    Some(value)
}

#[cfg(test)]
mod tests;
