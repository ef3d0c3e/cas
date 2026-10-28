//! Integrals that give the exponential, sine and cosine integrals and the
//! error function (CAS-DOC-31 §4, CAS-357), the heuristics after the
//! tabular rule:
//!
//! - **Laurent polynomials**: an exponential, sine, cosine, hyperbolic sine
//!   or cosine `F(u)` of a linear `u = a x + b`, times a rational function
//!   of `x` over ℚ whose denominator is `c (x − x₀)^n`, `n ≥ 1`. In powers of
//!   `y = x − x₀`, with `u = a y + b′`, the powers from 0 go to the
//!   integration of polynomials times `F` (the tabular rule), and the
//!   negative ones come down by parts to `1/y`: `∫ e^u/y^k = −e^u/((k − 1)
//!   y^(k−1)) + a/(k − 1) ∫ e^u/y^(k−1)`, and the sine's and cosine's trade
//!   places, `∫ sin u/y^k = −sin u/((k − 1) y^(k−1)) + a/(k − 1) ∫ cos u/y^(k−1)`
//!   and `∫ cos u/y^k = −cos u/((k − 1) y^(k−1)) − a/(k − 1) ∫ sin u/y^(k−1)`;
//!   then `∫ e^u/y = e^(b′) Ei(a y)`, `∫ sin u/y = cos b′ Si(a y) + sin b′
//!   Ci(a y)` and `∫ cos u/y = cos b′ Ci(a y) − sin b′ Si(a y)`; `sinh` and
//!   `cosh` through their exponentials. `∫ e^x (1 + x)/x²` is `2 Ei(x) −
//!   e^x/x`, and `∫ cos(2x + 1)/x` is `cos 1 Ci(2x) − sin 1 Si(2x)`, where
//!   giac writes half-angle tangents.
//! - **Logarithms**: `x^m/log x` is `Ei((m + 1) log x)` for a rational `m ≠
//!   −1`, giac's `int(x/ln(x),x)`.
//! - **Gaussians**: `P e^q` for a polynomial `P` over ℚ and a quadratic
//!   `q = α x² + β x + γ`, `α` free of `x` and shown negative: `∫ e^q =
//!   √π/(2√(−α)) e^(γ − β²/(4α)) erf(√(−α) (x + β/(2α)))`, and `P = A q′ +
//!   R` by division, so that `∫ P e^q = A e^q − ∫ A′ e^q + R ∫ e^q`, down
//!   the degrees of `A′`. A positive `α` stays: it needs `Erfi`.

use cas_num::algebra::{Rationals, Ring};
use cas_num::{Integer, Rational};
use cas_poly::{Dense, DenseRing};

use super::{call, integral_at, polynomial, rational_function, slope};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

/// The most powers of `1/y` brought down by parts, and the highest degree
/// of a polynomial times a Gaussian.
const MAX_DEGREE: usize = 64;

/// The functions of a linear argument the Laurent rule integrates.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Exp,
    Sin,
    Cos,
}

/// `t` as `Σ w F(u)`: `e^u`, `sin u`, `cos u`, and `sinh u`, `cosh u` as
/// halves of `e^u` and `e^(−u)`.
fn kernels(t: &Expr) -> Option<Vec<(Expr, Kernel, Expr)>> {
    let half = || Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
    match (t.head()?.as_symbol()?, t.args()) {
        (Symbol::EXP, [u]) => Some(vec![(Expr::from(1), Kernel::Exp, u.clone())]),
        (Symbol::POWER, [c, u]) if c.as_symbol() == Some(Symbol::E) => {
            Some(vec![(Expr::from(1), Kernel::Exp, u.clone())])
        }
        (Symbol::SIN, [u]) => Some(vec![(Expr::from(1), Kernel::Sin, u.clone())]),
        (Symbol::COS, [u]) => Some(vec![(Expr::from(1), Kernel::Cos, u.clone())]),
        (Symbol::SINH, [u]) => Some(vec![
            (half(), Kernel::Exp, u.clone()),
            (-half(), Kernel::Exp, -u.clone()),
        ]),
        (Symbol::COSH, [u]) => Some(vec![
            (half(), Kernel::Exp, u.clone()),
            (half(), Kernel::Exp, -u.clone()),
        ]),
        _ => None,
    }
}

fn function(kernel: Kernel, u: Expr) -> Expr {
    match kernel {
        Kernel::Exp => Expr::from(Symbol::E).pow(u),
        Kernel::Sin => call(Symbol::SIN, u),
        Kernel::Cos => call(Symbol::COS, u),
    }
}

/// The heuristics of the module on the factors of `f`, at `depth`.
pub(super) fn integral(
    ev: &mut Eval<'_>,
    e: &Expr,
    factors: &[Expr],
    x: &Expr,
    depth: u8,
) -> Found {
    if let Some(g) = laurent(ev, e, factors, x, depth)? {
        return Ok(Some(g));
    }
    if let Some(g) = logarithmic(ev, factors, x)? {
        return Ok(Some(g));
    }
    gaussian(ev, e, factors, x)
}

/// What is left of `factors` without the one at `i`: their product.
fn others(factors: &[Expr], i: usize) -> Expr {
    factors
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, f)| f.clone())
        .product()
}

/// `(c, x₀, n)` for `p = c (x − x₀)^n`.
fn power_of_linear(p: &Dense<Rational>) -> Option<(Rational, Rational, usize)> {
    let n = p.degree()?;
    let c = p.lead()?.clone();
    if n == 0 {
        return Some((c, Rational::ZERO, 0));
    }
    // x₀ from the next coefficient, c n (−x₀).
    let x0 = -(&p.coeffs()[n - 1] / &(&c * &Rational::from(Integer::from(n))));
    let qx = DenseRing::new(Rationals);
    let linear = qx.poly(vec![-x0.clone(), Rational::ONE]);
    let mut power = qx.poly(vec![c.clone()]);
    for _ in 0..n {
        power = qx.mul(&power, &linear);
    }
    (power == *p).then_some((c, x0, n))
}

/// The coefficients of `p(y + x₀)` in `y`, lowest first.
fn shifted(p: &Dense<Rational>, x0: &Rational) -> Vec<Rational> {
    // Horner's scheme on y + x₀.
    let mut out: Vec<Rational> = Vec::new();
    for c in p.coeffs().iter().rev() {
        // out := out (y + x₀) + c
        let mut next = vec![Rational::ZERO; out.len() + 1];
        for (k, a) in out.iter().enumerate() {
            next[k + 1] = &next[k + 1] + a;
            next[k] = &next[k] + &(a * x0);
        }
        next[0] = &next[0] + c;
        out = next;
    }
    out
}

/// `∫ F(u)/y^k` for `k = 1 … n`, `u = a y + b′` (§ the module comment).
fn reciprocal_powers(
    kernel: Kernel,
    a: &Expr,
    b0: &Expr,
    y: &Expr,
    u: &Expr,
    n: usize,
) -> Vec<Expr> {
    let ay = a.clone() * y.clone();
    let (sin_b, cos_b) = (call(Symbol::SIN, b0.clone()), call(Symbol::COS, b0.clone()));
    let (si, ci) = (
        call(Symbol::SIN_INTEGRAL, ay.clone()),
        call(Symbol::COS_INTEGRAL, ay.clone()),
    );
    let first = |k: Kernel| match k {
        Kernel::Exp => {
            Expr::from(Symbol::E).pow(b0.clone()) * call(Symbol::EXP_INTEGRAL_EI, ay.clone())
        }
        Kernel::Sin => cos_b.clone() * si.clone() + sin_b.clone() * ci.clone(),
        Kernel::Cos => cos_b.clone() * ci.clone() - sin_b.clone() * si.clone(),
    };
    // The sequences of the kernel and, for a sine or cosine, of the other.
    let other = match kernel {
        Kernel::Exp => Kernel::Exp,
        Kernel::Sin => Kernel::Cos,
        Kernel::Cos => Kernel::Sin,
    };
    let mut this = vec![first(kernel)];
    let mut that = vec![first(other)];
    for k in 2..=n {
        let m = Expr::from(i64::try_from(k - 1).expect("small"));
        let boundary = |f: Kernel| -function(f, u.clone()) / (m.clone() * y.pow(m.clone()));
        let step = a.clone() / m.clone();
        let (next_this, next_that) = match kernel {
            Kernel::Exp => {
                let v = boundary(Kernel::Exp) + step.clone() * this[k - 2].clone();
                (v.clone(), v)
            }
            Kernel::Sin => (
                boundary(Kernel::Sin) + step.clone() * that[k - 2].clone(),
                boundary(Kernel::Cos) - step.clone() * this[k - 2].clone(),
            ),
            Kernel::Cos => (
                boundary(Kernel::Cos) - step.clone() * that[k - 2].clone(),
                boundary(Kernel::Sin) + step.clone() * this[k - 2].clone(),
            ),
        };
        this.push(next_this);
        that.push(next_that);
    }
    this
}

/// The Laurent rule (§ the module comment).
fn laurent(ev: &mut Eval<'_>, e: &Expr, factors: &[Expr], x: &Expr, depth: u8) -> Found {
    for (i, t) in factors.iter().enumerate() {
        let Some(parts) = kernels(t) else {
            continue;
        };
        let rest = others(factors, i);
        let Some((num, den)) = rational_function(ev, e, &rest, x)? else {
            continue;
        };
        let Some((c, x0, n)) = power_of_linear(&den) else {
            continue;
        };
        if n == 0 || n > MAX_DEGREE {
            continue;
        }
        let x0e = Expr::from(x0.clone());
        let y = ev.eval(&(x.clone() - x0e.clone()))?;
        let coefficients: Vec<Rational> = shifted(&num, &x0).iter().map(|b| b / &c).collect();
        let mut total = Vec::new();
        for (w, kernel, u) in parts {
            let u = ev.eval(&u)?;
            let Some(a) = slope(&u, x) else {
                return Ok(None);
            };
            let b0 = ev.eval(&u.substitute(&[(x.clone(), x0e.clone())]))?;
            let powers = reciprocal_powers(kernel, &a, &b0, &y, &u, n);
            let mut polynomial_part = Vec::new();
            for (j, cj) in coefficients.iter().enumerate() {
                if cj.is_zero() {
                    continue;
                }
                let cj = Expr::from(cj.clone());
                if j < n {
                    total.push(w.clone() * cj * powers[n - j - 1].clone());
                } else {
                    let k = i64::try_from(j - n).expect("small");
                    polynomial_part.push(cj * y.pow(k));
                }
            }
            if !polynomial_part.is_empty() {
                let integrand = ev.eval(
                    &(polynomial_part.into_iter().sum::<Expr>() * function(kernel, u.clone())),
                )?;
                let Some(g) = integral_at(ev, e, &integrand, x, depth + 1)? else {
                    return Ok(None);
                };
                if g.contains(&Expr::from(Symbol::INTEGRATE)) {
                    return Ok(None);
                }
                total.push(w * g);
            }
        }
        return Ok(Some(total.into_iter().sum()));
    }
    Ok(None)
}

/// `x^m/log x` (§ the module comment).
fn logarithmic(ev: &mut Eval<'_>, factors: &[Expr], x: &Expr) -> Found {
    let reciprocal_log = Expr::call(Symbol::LOG, [x.clone()]).pow(-1);
    let Some(i) = factors.iter().position(|f| *f == reciprocal_log) else {
        return Ok(None);
    };
    let rest = others(factors, i);
    let m = if rest == *x {
        Rational::ONE
    } else {
        match (rest.head().and_then(Expr::as_symbol), rest.args()) {
            (Some(Symbol::POWER), [b, k]) if b == x => match k.to_rational() {
                Some(m) => m,
                None => return Ok(None),
            },
            _ => return Ok(None),
        }
    };
    let next = &m + &Rational::ONE;
    if next.is_zero() {
        return Ok(None);
    }
    let argument = ev.eval(&(Expr::from(next) * Expr::call(Symbol::LOG, [x.clone()])))?;
    Ok(Some(call(Symbol::EXP_INTEGRAL_EI, argument)))
}

/// The quadratic exponent `q = α x² + β x + γ` of a factor `e^q`, and
/// `[γ, β, α]`.
fn quadratic(ev: &mut Eval<'_>, t: &Expr, x: &Expr) -> Result<Option<[Expr; 3]>, EvalError> {
    let q = match (t.head().and_then(Expr::as_symbol), t.args()) {
        (Some(Symbol::EXP), [q]) => q,
        (Some(Symbol::POWER), [c, q]) if c.as_symbol() == Some(Symbol::E) => q,
        _ => return Ok(None),
    };
    if !q.contains(x) {
        return Ok(None);
    }
    let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [q.clone(), x.clone()]))?;
    if degree.as_i64() != Some(2) {
        return Ok(None);
    }
    let mut c = Vec::with_capacity(3);
    for k in 0..3 {
        let ck = ev.eval(&Expr::call(
            Symbol::COEFFICIENT,
            [q.clone(), x.clone(), Expr::from(k)],
        ))?;
        if ck.contains(x) {
            return Ok(None);
        }
        c.push(ck);
    }
    // q must be that polynomial.
    let rebuilt = c[0].clone() + c[1].clone() * x.clone() + c[2].clone() * x.pow(2);
    let rest = ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [q.clone() - rebuilt])],
    ))?;
    if rest != Expr::from(0) {
        return Ok(None);
    }
    Ok(c.try_into().ok())
}

/// Whether `a` is shown negative: a negative number, or by the
/// assumptions.
fn negative(ev: &mut Eval<'_>, a: &Expr) -> Result<bool, EvalError> {
    if let Some(r) = a.to_rational() {
        return Ok(r.is_negative());
    }
    let test = ev.eval(&Expr::call(Symbol::LESS, [a.clone(), Expr::from(0)]))?;
    Ok(test.as_symbol() == Some(Symbol::TRUE))
}

/// The Gaussian rule (§ the module comment).
fn gaussian(ev: &mut Eval<'_>, e: &Expr, factors: &[Expr], x: &Expr) -> Found {
    for (i, t) in factors.iter().enumerate() {
        let Some([gamma, beta, alpha]) = quadratic(ev, t, x)? else {
            continue;
        };
        if !negative(ev, &alpha)? {
            continue;
        }
        let rest = others(factors, i);
        let Some(p) = polynomial(ev, e, &rest, x)? else {
            continue;
        };
        if p.coeffs().len() > MAX_DEGREE {
            continue;
        }
        let root = Expr::call(Symbol::SQRT, [-alpha.clone()]);
        let shift = beta.clone() / (Expr::from(2) * alpha.clone());
        let base = Expr::call(Symbol::SQRT, [Expr::from(Symbol::PI)])
            / (Expr::from(2) * root.clone())
            * Expr::from(Symbol::E)
                .pow(gamma.clone() - beta.pow(2) / (Expr::from(4) * alpha.clone()))
            * call(Symbol::ERF, root * (x.clone() + shift));
        let base = ev.eval(&base)?;
        let exponential = t.clone();
        let derivative = [beta.clone(), Expr::from(2) * alpha.clone()];
        let coefficients: Vec<Expr> = p.coeffs().iter().map(|c| Expr::from(c.clone())).collect();
        let g = reduced(ev, coefficients, &derivative, &exponential, &base, x)?;
        return Ok(Some(g));
    }
    Ok(None)
}

/// `∫ P e^q` for `P` of coefficients `p` (lowest first), `q′ = s + d x`
/// as `[s, d]`, `base = ∫ e^q`.
fn reduced(
    ev: &mut Eval<'_>,
    mut p: Vec<Expr>,
    derivative: &[Expr; 2],
    exponential: &Expr,
    base: &Expr,
    x: &Expr,
) -> Result<Expr, EvalError> {
    while p.last().is_some_and(|c| *c == Expr::from(0)) {
        p.pop();
    }
    if p.len() <= 1 {
        let c = p.pop().unwrap_or_else(|| Expr::from(0));
        return ev.eval(&(c * base.clone()));
    }
    // P = A (d x + s) + R.
    let [s, d] = derivative;
    let degree = p.len() - 1;
    let mut a = vec![Expr::from(0); degree];
    for k in (1..=degree).rev() {
        let ak = ev.eval(&(p[k].clone() / d.clone()))?;
        p[k - 1] = ev.eval(&(p[k - 1].clone() - ak.clone() * s.clone()))?;
        a[k - 1] = ak;
    }
    let r = p[0].clone();
    let a_expr: Expr = a
        .iter()
        .enumerate()
        .map(|(k, c)| c.clone() * x.pow(i64::try_from(k).expect("small")))
        .sum();
    let derivative_of_a: Vec<Expr> = a
        .iter()
        .enumerate()
        .skip(1)
        .map(|(k, c)| Expr::from(i64::try_from(k).expect("small")) * c.clone())
        .collect();
    let inner = reduced(ev, derivative_of_a, derivative, exponential, base, x)?;
    ev.eval(&(a_expr * exponential.clone() - inner + r * base.clone()))
}

#[cfg(test)]
mod tests;
