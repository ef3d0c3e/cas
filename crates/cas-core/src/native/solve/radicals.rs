//! Roots in radicals of irreducible polynomials of degree 3 and more
//! (CAS-DOC-11 §3, CAS-211), as expressions built, not evaluated.
//!
//! - **Binomials** `a xⁿ + b`, parameters allowed, `c = −b/a`: in the reals
//!   `c^(1/n)`, and `−c^(1/n)` too for an even `n`, the real root of a
//!   negative `c` for an odd one and none for an even one (giac's
//!   `solve(x^4-a,x)` is `[-a^(1/4),a^(1/4)]`); in the complex numbers
//!   `|c|^(1/n) (cos θₖ + i sin θₖ)` for `θₖ = 2kπ/n`, `(2k + 1)π/n` for a
//!   negative `c`.
//! - **Quadratics in `xᵐ`**, `A x²ᵐ + B xᵐ + C` over ℚ: the binomials
//!   `xᵐ = y` for the two roots `y = α ± β √Δ`, the square roots of them
//!   denested when `α² − β² Δ` is a square, as giac's `x^4-10x^2+1` gives
//!   `±√2 ± √3`.
//! - **Cubics** over ℚ, `t = x + b/(3a)` and `t³ + p t + q`, irreducible:
//!   for `D = (q/2)² + (p/3)³ > 0` one real root, Cardano's
//!   `∛(−q/2 + √D) + ∛(−q/2 − √D)` with real cube roots, and the complex
//!   ones `−(u + v)/2 ± i (√3/2)(u − v)`; for `D < 0` three real roots,
//!   `2 √(−p/3) cos(θ/3 − 2πk/3)` with `θ = acos((3q/(2p)) √(−3/p))`. giac
//!   writes `rootof` for both.
//! - **Quartics** over ℚ, `y = x + b/(4a)` and `y⁴ + p y² + q y + r`:
//!   `q = 0` a quadratic in `y²`; otherwise Ferrari's, with a positive
//!   rational root `m` of the resolvent `8m³ + 8pm² + (2p² − 8r)m − q²`, the
//!   quartic as `(y² − σy + p/2 + m + q/(2σ))(y² + σy + p/2 + m − q/(2σ))`,
//!   `σ = √(2m)`. Without such a root, none here: giac gives the roots of
//!   such quartics numerically, and so does the caller.
//!
//! Signs of numbers `α + β√n` are exact, so each root is known real or not.

use std::cmp::Ordering;

use cas_num::algebra::Rationals;
use cas_num::{Integer, Rational};
use cas_poly::{Cancel, Cancelled, DenseRing};

use crate::expr::Expr;
use crate::symbol::Symbol;

fn rational(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn number(q: &Rational) -> Expr {
    Expr::from(q.clone())
}

/// `e^(1/n)`.
fn root(e: Expr, n: u64) -> Expr {
    let exponent = Expr::from(Rational::new(Integer::ONE, Integer::from(n)));
    Expr::call(Symbol::POWER, [e, exponent])
}

fn sqrt(e: Expr) -> Expr {
    root(e, 2)
}

/// `i e`.
fn times_i(e: Expr) -> Expr {
    Expr::call(Symbol::TIMES, [Expr::from(Symbol::I), e])
}

/// The square root of a rational that is a square.
fn rational_sqrt(q: &Rational) -> Option<Rational> {
    if q.is_negative() {
        return None;
    }
    let (n, n_rest) = q.numer().sqrt_rem();
    let (d, d_rest) = q.denom().sqrt_rem();
    (n_rest.is_zero() && d_rest.is_zero()).then(|| Rational::new(n, d))
}

/// The sign of `α + β √n`, `n > 0`, exactly.
pub(in crate::native) fn surd_sign(alpha: &Rational, beta: &Rational, n: &Rational) -> Ordering {
    let a = alpha.cmp(&Rational::ZERO);
    let b = beta.cmp(&Rational::ZERO);
    if b == Ordering::Equal || a == b {
        return if a == Ordering::Equal { b } else { a };
    }
    if a == Ordering::Equal {
        return b;
    }
    // Opposite signs: the larger of α² and β² n wins.
    match (alpha * alpha).cmp(&(&(beta * beta) * n)) {
        Ordering::Greater => a,
        Ordering::Less => b,
        Ordering::Equal => Ordering::Equal,
    }
}

/// `α + β √n`, built.
fn surd(alpha: &Rational, beta: &Rational, n: &Rational) -> Expr {
    number(alpha) + number(beta) * sqrt(number(n))
}

/// `√(α + β √n)` for `α + β √n > 0`, denested when `α² − β² n = s²`:
/// `√((α + s)/2) + sgn(β) √((α − s)/2)`.
pub(in crate::native) fn surd_sqrt(alpha: &Rational, beta: &Rational, n: &Rational) -> Expr {
    if beta.is_zero() {
        return sqrt(number(alpha));
    }
    let d = &(alpha * alpha) - &(&(beta * beta) * n);
    if let Some(s) = rational_sqrt(&d) {
        let two = Rational::from(Integer::from(2));
        let x = sqrt(number(&(&(alpha + &s) / &two)));
        let y = sqrt(number(&(&(alpha - &s) / &two)));
        return if beta.is_negative() { x - y } else { x + y };
    }
    sqrt(surd(alpha, beta, n))
}

/// The roots of `xⁿ = c`, `c ≠ 0`, `sign` its sign when it is known real;
/// the real ones when `real`.
pub(super) fn binomial(c: &Expr, n: u64, sign: Option<Ordering>, real: bool) -> Vec<Expr> {
    let odd = n % 2 == 1;
    if real {
        return match sign {
            Some(Ordering::Less) if odd => vec![-root(-c.clone(), n)],
            Some(Ordering::Less) => Vec::new(),
            _ if odd => vec![root(c.clone(), n)],
            _ => vec![root(c.clone(), n), -root(c.clone(), n)],
        };
    }
    let (modulus, offset) = match sign {
        Some(Ordering::Less) => (root(-c.clone(), n), 1),
        _ => (root(c.clone(), n), 0),
    };
    (0..n)
        .map(|k| {
            let turn = Rational::new(Integer::from(2 * k + offset), Integer::from(n));
            let angle = number(&turn) * Expr::from(Symbol::PI);
            let unit = Expr::call(Symbol::COS, [angle.clone()])
                + times_i(Expr::call(Symbol::SIN, [angle]));
            modulus.clone() * unit
        })
        .collect()
}

/// The roots of `xᵐ = y` for `y = α + β √n` real, `n > 0` not a square:
/// denested square roots for `m = 2`.
fn power_roots(alpha: &Rational, beta: &Rational, n: &Rational, m: u64, real: bool) -> Vec<Expr> {
    let sign = surd_sign(alpha, beta, n);
    if m != 2 {
        return binomial(&surd(alpha, beta, n), m, Some(sign), real);
    }
    match sign {
        Ordering::Less if real => Vec::new(),
        Ordering::Less => {
            let r = times_i(surd_sqrt(&-alpha, &-beta, n));
            vec![r.clone(), -r]
        }
        _ => {
            let r = surd_sqrt(alpha, beta, n);
            vec![r.clone(), -r]
        }
    }
}

/// The roots of `A x²ᵐ + B xᵐ + C`, irreducible over ℚ, plus `shift`.
pub(super) fn quadratic_in_power(
    [c, b, a]: [&Rational; 3],
    m: u64,
    real: bool,
    shift: &Rational,
) -> Vec<Expr> {
    let four = Rational::from(Integer::from(4));
    let delta = &(b * b) - &(&(&four * a) * c);
    let two_a = &Rational::from(Integer::from(2)) * a;
    let alpha = &-b / &two_a;
    let beta = two_a.checked_recip().expect("a quadratic");
    let mut roots = Vec::new();
    if delta.is_positive() {
        for beta in [beta.clone(), -beta] {
            roots.extend(power_roots(&alpha, &beta, &delta, m, real));
        }
    } else if !real {
        // y = α ± β i √(−Δ), complex.
        let width = number(&beta) * sqrt(number(&-delta));
        for y in [
            number(&alpha) + times_i(width.clone()),
            number(&alpha) - times_i(width),
        ] {
            roots.extend(binomial(&y, m, None, false));
        }
    }
    shifted(roots, shift)
}

/// `e + shift` for each `e`.
fn shifted(roots: Vec<Expr>, shift: &Rational) -> Vec<Expr> {
    if shift.is_zero() {
        return roots;
    }
    roots.into_iter().map(|e| e + number(shift)).collect()
}

/// The real cube root of `e`, of the sign given: `e^(1/3)` or
/// `−(−e)^(1/3)`.
fn real_cube_root(e: Expr, positive: bool) -> Expr {
    if positive { root(e, 3) } else { -root(-e, 3) }
}

/// The roots of the irreducible `a x³ + b x² + c x + d` over ℚ, the real
/// ones when `real`.
pub(super) fn cubic([d, c, b, a]: [&Rational; 4], real: bool) -> Vec<Expr> {
    let r = |n: i64| Rational::from(Integer::from(n));
    let shift = &-b / &(&r(3) * a);
    let p = &(&(&(&r(3) * a) * c) - &(b * b)) / &(&(&r(3) * a) * a);
    let q = &(&(&(&(&r(2) * b) * b) * b) - &(&(&(&r(9) * a) * b) * c)
        + &(&(&(&r(27) * a) * a) * d))
        / &(&(&(&r(27) * a) * a) * a);
    if p.is_zero() {
        // t³ = −q.
        let c = -&q;
        let sign = c.cmp(&Rational::ZERO);
        return shifted(binomial(&number(&c), 3, Some(sign), real), &shift);
    }
    let half_q = &q / &r(2);
    let third_p = &p / &r(3);
    let discriminant = &(&half_q * &half_q) + &(&(&third_p * &third_p) * &third_p);
    let roots = if discriminant.is_positive() {
        let s = sqrt(number(&discriminant));
        // −q/2 + √D > 0 unless q > 0 and p < 0; −q/2 − √D > 0 when q < 0 and
        // p < 0.
        let u = real_cube_root(
            number(&-&half_q) + s.clone(),
            !q.is_positive() || p.is_positive(),
        );
        let v = real_cube_root(number(&-&half_q) - s, q.is_negative() && p.is_negative());
        let sum = u.clone() + v.clone();
        if real {
            vec![sum]
        } else {
            let center = number(&rational(-1, 2)) * sum.clone();
            let width = times_i(number(&rational(1, 2)) * sqrt(Expr::from(3)) * (u - v));
            vec![sum, center.clone() + width.clone(), center - width]
        }
    } else {
        // Three real roots: p < 0.
        let modulus = Expr::from(2) * sqrt(number(&-&third_p));
        let cosine = number(&(&(&r(3) * &q) / &(&r(2) * &p))) * sqrt(number(&(&r(-3) / &p)));
        let third = number(&rational(1, 3)) * Expr::call(Symbol::ARC_COS, [cosine]);
        (0..3)
            .map(|k| {
                let turn = number(&rational(-2 * k, 3)) * Expr::from(Symbol::PI);
                modulus.clone() * Expr::call(Symbol::COS, [third.clone() + turn])
            })
            .collect()
    };
    shifted(roots, &shift)
}

/// The roots of the irreducible `a x⁴ + b x³ + c x² + d x + e` over ℚ, the
/// real ones when `real`, by Ferrari's method; `None` when the resolvent
/// cubic has no positive rational root.
///
/// # Errors
///
/// [`Cancelled`] when `cancel` says so, factoring the resolvent.
pub(super) fn quartic(
    [e, d, c, b, a]: [&Rational; 5],
    real: bool,
    cancel: &dyn Cancel,
) -> Result<Option<Vec<Expr>>, Cancelled> {
    let r = |n: i64| Rational::from(Integer::from(n));
    let shift = &-b / &(&r(4) * a);
    let (a2, b2) = (a * a, b * b);
    let p = &(&(&(&r(8) * a) * c) - &(&r(3) * &b2)) / &(&r(8) * &a2);
    let q = &(&(&(&b2 * b) - &(&(&(&r(4) * a) * b) * c)) + &(&(&r(8) * &a2) * d))
        / &(&(&r(8) * &a2) * a);
    let s = &(&(&(&(&r(-3) * &b2) * &b2) + &(&(&(&r(256) * &a2) * a) * e))
        - &(&(&(&r(64) * &a2) * b) * d))
        + &(&(&(&r(16) * a) * &b2) * c);
    let s = &s / &(&(&r(256) * &a2) * &a2);
    if q.is_zero() {
        return Ok(Some(quadratic_in_power(
            [&s, &p, &Rational::ONE],
            2,
            real,
            &shift,
        )));
    }
    // 8m³ + 8pm² + (2p² − 8s)m − q².
    let qx = DenseRing::new(Rationals);
    let resolvent = qx.poly(vec![
        -&(&q * &q),
        &(&r(2) * &(&p * &p)) - &(&r(8) * &s),
        &r(8) * &p,
        r(8),
    ]);
    let factors = qx.factor_with(&resolvent, cancel)?;
    let m = factors.factors().iter().find_map(|(f, _)| {
        let [c0, c1] = f.coeffs() else {
            return None;
        };
        let m = &-c0 / c1;
        m.is_positive().then_some(m)
    });
    let Some(m) = m else {
        return Ok(None);
    };
    // σ = √(2m), irrational as the quartic is irreducible; each quadratic
    // y² − εσy + … has the discriminant α + β σ, α = −2p − 2m, β = −ε q/m.
    let n = &r(2) * &m;
    let sigma = sqrt(number(&n));
    let alpha = &(&r(-2) * &p) - &(&r(2) * &m);
    let mut roots = Vec::new();
    for epsilon in [1, -1] {
        let beta = &(&-&q / &m) * &r(epsilon);
        let center = number(&rational(epsilon, 2)) * sigma.clone();
        let half = number(&rational(1, 2));
        match surd_sign(&alpha, &beta, &n) {
            Ordering::Less if real => {}
            Ordering::Less => {
                let width = half * times_i(surd_sqrt(&-&alpha, &-&beta, &n));
                roots.extend([center.clone() + width.clone(), center - width]);
            }
            _ => {
                let width = half * surd_sqrt(&alpha, &beta, &n);
                roots.extend([center.clone() + width.clone(), center - width]);
            }
        }
    }
    Ok(Some(shifted(roots, &shift)))
}
