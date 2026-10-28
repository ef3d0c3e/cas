//! Polynomial tools (CAS-DOC-29, CAS-339 to CAS-341): Xcas's `partfrac`,
//! `symb2poly`, `coeffs` and `poly2symb`, `tcoeff` and `valuation`,
//! `sturmab` and `realroot`, as `Apart`, `CoefficientList`,
//! `FromCoefficients`, `TrailingCoefficient`, `Valuation`, `CountRoots`
//! and `RealRoots`. The variable is `x` when left out, as in giac.
//!
//! - **`Apart[f, x]`** (CAS-339): `f` put together as `p/q`, `q` factored
//!   (`FactorList`) over the rationals or the rational functions of the
//!   other symbols, and an irreducible quadratic with rational
//!   coefficients and a positive discriminant split in radicals, as giac
//!   does. With the quotient `Q` and the remainder `r` of `p` by `q`, a
//!   numerator `A_{g,j}` of degree `< deg g` with unknown coefficients for
//!   each factor `g^e` and `j ≤ e`; `r = Σ A_{g,j} q/g^j` is solved power by
//!   power (`LinSolve`), and `f = Q + Σ A_{g,j}/g^j`. A denominator free of
//!   `x` gives `f`; one of a degree above [`MAX_DEGREE`] leaves the call.
//! - **`CoefficientList[p, x]`** and **`FromCoefficients[l, x]`**
//!   (CAS-340): the coefficients from the leading one down, as giac's
//!   `poly1[1,0,-2,1]`, and back. **`TrailingCoefficient`** and
//!   **`Valuation`**: the lowest term; 0 has valuation `-Infinity`, as
//!   giac's `valuation(0)`. An expression that is no polynomial in `x`
//!   leaves the call.
//! - **`CountRoots[p, x, a, b]`** (CAS-341): the distinct real roots of a
//!   polynomial over ℚ in `[a, b]`, rational bounds or infinities, by
//!   Sturm's sequence. **`RealRoots[p]`** and `RealRoots[p, ε]`: the real
//!   roots in balls of the working precision, or of `−log₂ ε` bits, each
//!   with its multiplicity, `{{root, m}, …}` ascending, from the
//!   square-free factors' roots (`polynomial_roots`, CAS-211).
//!
//! Each call that answers takes one step, its builtin's name.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomials: p over q in x, factors g of multiplicity e, coefficients c"
)]

use cas_num::algebra::Rationals;
use cas_num::ball::{ComplexBall, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::DenseRing;

use super::integrate::polynomial;
use super::solve::linear_system;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The largest degree of a denominator `Apart` takes apart.
const MAX_DEGREE: i64 = 64;

/// The largest degree `CoefficientList` writes out.
const MAX_COEFFICIENTS: i64 = 1 << 16;

fn call(head: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(head, args)
}

/// The expression and the variable of `(p)` and `(p, x)`.
fn with_variable(args: &[Expr]) -> Option<(&Expr, Expr)> {
    match args {
        [p] => Some((p, Expr::symbol("x"))),
        [p, x] if x.as_symbol().is_some() => Some((p, x.clone())),
        _ => None,
    }
}

/// `Exponent[p, x]` as a natural number; `None` for another.
fn degree(ev: &mut Eval<'_>, p: &Expr, x: &Expr) -> Result<Option<i64>, EvalError> {
    let n = ev.eval(&call(Symbol::EXPONENT, vec![p.clone(), x.clone()]))?;
    Ok(n.as_i64().filter(|n| *n >= 0))
}

/// The coefficients of `p` in `x`, from `x^0` up, when `p` is a
/// polynomial in `x`: checked by summing them back.
fn coefficients(ev: &mut Eval<'_>, p: &Expr, x: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    let expanded = ev.eval(&call(Symbol::EXPAND, vec![p.clone()]))?;
    if expanded == Expr::from(0) {
        return Ok(Some(Vec::new()));
    }
    let Some(n) = degree(ev, &expanded, x)?.filter(|n| *n <= MAX_COEFFICIENTS) else {
        return Ok(None);
    };
    let mut cs = Vec::with_capacity(usize::try_from(n + 1).expect("small"));
    for k in 0..=n {
        let c = ev.eval(&call(
            Symbol::COEFFICIENT,
            vec![expanded.clone(), x.clone(), Expr::from(k)],
        ))?;
        if c.contains(x) {
            return Ok(None);
        }
        cs.push(c);
    }
    let back = Expr::call(
        Symbol::PLUS,
        cs.iter()
            .enumerate()
            .map(|(k, c)| c.clone() * x.pow(i64::try_from(k).expect("small")))
            .collect::<Vec<_>>(),
    );
    let difference = ev.eval(&call(Symbol::EXPAND, vec![back - expanded]))?;
    Ok((difference == Expr::from(0)).then_some(cs))
}

/// The roots `(−b ± √d)/(2a)` of `a x² + b x + c` with rational
/// coefficients and a positive discriminant `d` that is no square: the
/// split of an irreducible quadratic, as giac's `partfrac`; `None` for
/// another factor.
fn radical_roots(
    ev: &mut Eval<'_>,
    g: &Expr,
    x: &Expr,
) -> Result<Option<(Expr, [Expr; 2])>, EvalError> {
    let Some(cs) = coefficients(ev, g, x)? else {
        return Ok(None);
    };
    let [c, b, a] = cs.as_slice() else {
        return Ok(None);
    };
    let (Some(a), Some(b), Some(c)) = (a.to_rational(), b.to_rational(), c.to_rational()) else {
        return Ok(None);
    };
    let four = Rational::from(Integer::from(4));
    let d = &(&b * &b) - &(&four * &(&a * &c));
    if !d.is_positive() {
        return Ok(None);
    }
    let root = Expr::from(d).pow(Expr::from(Rational::new(
        Integer::from(1),
        Integer::from(2),
    )));
    let two_a = Expr::from(&Rational::from(Integer::from(2)) * &a);
    let mut roots = Vec::with_capacity(2);
    for s in [-1, 1] {
        roots.push(
            ev.eval(&((Expr::from(-b.clone()) + Expr::from(s) * root.clone()) / two_a.clone()))?,
        );
    }
    let [r1, r2]: [Expr; 2] = roots.try_into().expect("two roots");
    Ok(Some((Expr::from(a), [r1, r2])))
}

pub(super) fn apart(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let together = ev.eval(&call(Symbol::TOGETHER, vec![f.clone()]))?;
    let p = ev.eval(&call(Symbol::NUMERATOR, vec![together.clone()]))?;
    let q = ev.eval(&call(Symbol::DENOMINATOR, vec![together.clone()]))?;
    if !q.contains(&x) {
        return Ok(Some(together));
    }
    let Some(n) = degree(ev, &q, &x)?.filter(|n| *n <= MAX_DEGREE) else {
        return Ok(None);
    };
    // The factors of q with x, and the constant.
    let listed = ev.eval(&call(Symbol::FACTOR_LIST, vec![q.clone()]))?;
    let mut lead = Expr::from(1);
    let mut factors: Vec<(Expr, i64)> = Vec::new();
    for pair in listed.args().chunks(2) {
        let [g, k] = pair else {
            return Ok(None);
        };
        let Some(k) = k.as_i64().filter(|k| *k >= 1) else {
            return Ok(None);
        };
        if !g.contains(&x) {
            lead *= g.pow(k);
            continue;
        }
        if degree(ev, g, &x)? == Some(2)
            && let Some((a, [r1, r2])) = radical_roots(ev, g, &x)?
        {
            lead *= a.pow(k);
            factors.push((x.clone() - r1, k));
            factors.push((x.clone() - r2, k));
        } else {
            factors.push((g.clone(), k));
        }
    }
    let quotient = ev.eval(&call(
        Symbol::POLYNOMIAL_QUOTIENT,
        vec![p.clone(), q.clone(), x.clone()],
    ))?;
    let rest = ev.eval(&call(Symbol::POLYNOMIAL_REMAINDER, vec![p, q, x.clone()]))?;
    // The numerators with unknown coefficients, and their cofactors.
    let mut unknowns = Vec::new();
    let mut parts: Vec<(Expr, Expr, i64)> = Vec::new();
    let mut identity = Vec::new();
    for (i, (g, e)) in factors.iter().enumerate() {
        let Some(d) = degree(ev, g, &x)? else {
            return Ok(None);
        };
        for j in 1..=*e {
            let mut numerator = Vec::new();
            for k in 0..d {
                let c = Expr::symbol(&format!("Apart$c{}", unknowns.len()));
                unknowns.push(c.clone());
                numerator.push(c * x.pow(k));
            }
            let numerator = Expr::call(Symbol::PLUS, numerator);
            let mut cofactor = vec![lead.clone()];
            for (h, (other, k)) in factors.iter().enumerate() {
                let power = if h == i { e - j } else { *k };
                cofactor.push(other.pow(power));
            }
            identity.push(numerator.clone() * Expr::call(Symbol::TIMES, cofactor));
            parts.push((numerator, g.clone(), j));
        }
    }
    let residual = ev.eval(&call(
        Symbol::EXPAND,
        vec![Expr::call(Symbol::PLUS, identity) - rest],
    ))?;
    let mut equations = Vec::with_capacity(usize::try_from(n).expect("small"));
    for k in 0..n {
        let c = ev.eval(&call(
            Symbol::COEFFICIENT,
            vec![residual.clone(), x.clone(), Expr::from(k)],
        ))?;
        equations.push(Expr::call(Symbol::EQUATION, [c, Expr::from(0)]));
    }
    let Some(found) = linear_system(ev, e, &equations, &unknowns)? else {
        return Ok(None);
    };
    let [values] = found.args() else {
        return Ok(None);
    };
    if values.args().len() != unknowns.len() {
        return Ok(None);
    }
    let rules: Vec<(Expr, Expr)> = unknowns
        .into_iter()
        .zip(values.args().iter().cloned())
        .collect();
    let mut terms = vec![quotient];
    for (numerator, g, j) in parts {
        let numerator = ev.eval(&numerator.substitute(&rules))?;
        if numerator != Expr::from(0) {
            terms.push(numerator / g.pow(j));
        }
    }
    let result = ev.eval(&Expr::call(Symbol::PLUS, terms))?;
    ev.step("Apart", e, &result);
    Ok(Some(result))
}

pub(super) fn coefficient_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let Some(mut cs) = coefficients(ev, p, &x)? else {
        return Ok(None);
    };
    if cs.is_empty() {
        cs.push(Expr::from(0));
    }
    cs.reverse();
    Ok(Some(Expr::call(Symbol::LIST, cs)))
}

pub(super) fn from_coefficients(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, x) = match e.args() {
        [l] => (l, Expr::symbol("x")),
        [l, x] if x.as_symbol().is_some() => (l, x.clone()),
        _ => return Ok(None),
    };
    if l.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let n = l.args().len();
    let terms: Vec<Expr> = l
        .args()
        .iter()
        .enumerate()
        .map(|(i, c)| c.clone() * x.pow(i64::try_from(n - 1 - i).expect("small")))
        .collect();
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// The lowest term of a polynomial.
enum Lowest {
    /// The polynomial 0.
    Zero,
    /// The power and its coefficient.
    Term(i64, Expr),
}

/// The lowest term of `p` in `x`; `None` for what is no polynomial in `x`.
fn lowest(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Lowest>, EvalError> {
    let Some((p, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let Some(cs) = coefficients(ev, p, &x)? else {
        return Ok(None);
    };
    let zero = Expr::from(0);
    Ok(Some(
        cs.into_iter()
            .enumerate()
            .find(|(_, c)| *c != zero)
            .map_or(Lowest::Zero, |(k, c)| {
                Lowest::Term(i64::try_from(k).expect("small"), c)
            }),
    ))
}

pub(super) fn trailing_coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(lowest(ev, e)?.map(|low| match low {
        Lowest::Zero => Expr::from(0),
        Lowest::Term(_, c) => c,
    }))
}

pub(super) fn valuation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(lowest(ev, e)?.map(|low| match low {
        Lowest::Zero => Expr::call(
            Symbol::TIMES,
            [Expr::from(-1), Expr::from(Symbol::INFINITY)],
        ),
        Lowest::Term(k, _) => Expr::from(k),
    }))
}

/// A bound of `CountRoots`: a rational, or an infinity of the sign asked.
enum Bound {
    Finite(Rational),
    Infinite,
}

impl Bound {
    fn of(b: &Expr, positive: bool) -> Option<Bound> {
        if let Some(r) = b.to_rational() {
            return Some(Bound::Finite(r));
        }
        // `-Infinity` as `Times[-1, Infinity]` or `DirectedInfinity[-1]`.
        let infinite = b.contains(&Expr::from(Symbol::INFINITY))
            || b.contains(&Expr::from(Symbol::DIRECTED_INFINITY));
        if !infinite || b.contains(&Expr::from(Symbol::COMPLEX_INFINITY)) {
            return None;
        }
        let negative = b.contains(&Expr::from(-1));
        (negative != positive).then_some(Bound::Infinite)
    }

    fn rational(&self) -> Option<&Rational> {
        match self {
            Bound::Finite(r) => Some(r),
            Bound::Infinite => None,
        }
    }
}

pub(super) fn count_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, x, a, b] = e.args() else {
        return Ok(None);
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let (Some(lo), Some(hi)) = (Bound::of(a, false), Bound::of(b, true)) else {
        return Ok(None);
    };
    let Some(dense) = polynomial(ev, e, p, x)? else {
        return Ok(None);
    };
    if dense.is_zero() {
        return Ok(None);
    }
    let qx = DenseRing::new(Rationals);
    let count = qx.count_real_roots(&dense, lo.rational(), hi.rational());
    let result = Expr::from(i64::try_from(count).expect("a count in memory"));
    ev.step("CountRoots", e, &result);
    Ok(Some(result))
}

pub(super) fn real_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, bits) = match e.args() {
        [p] => (p, ev.context().working_precision()),
        [p, eps] => {
            let Some(eps) = eps.to_rational().filter(Rational::is_positive) else {
                return Ok(None);
            };
            // −log₂ ε, from the bit lengths: within a bit or two.
            let bits = eps
                .denom()
                .bit_length()
                .saturating_sub(eps.numer().unsigned_abs().bit_length());
            (p, bits.max(8))
        }
        _ => return Ok(None),
    };
    let x = Expr::symbol("x");
    let Some(dense) = polynomial(ev, e, p, &x)? else {
        return Ok(None);
    };
    if dense.degree().is_none_or(|d| d == 0) {
        return Ok(None);
    }
    let qx = DenseRing::new(Rationals);
    let mut found: Vec<(ComplexBall, i64)> = Vec::new();
    for (factor, m) in qx.squarefree(&dense).factors() {
        if factor.degree().is_none_or(|d| d == 0) {
            continue;
        }
        let cs: Vec<Complex> = factor.coeffs().iter().cloned().map(Complex::from).collect();
        let Some(roots) = polynomial_roots(&cs, bits, 64 * bits.max(64)) else {
            return Ok(None);
        };
        let m = i64::try_from(*m).expect("a multiplicity in memory");
        found.extend(
            roots
                .into_iter()
                .filter(ComplexBall::is_real)
                .map(|z| (z, m)),
        );
    }
    // By midpoint, a total order: the roots' balls are apart.
    found.sort_by(|(a, _), (b, _)| a.re().mid().cmp(b.re().mid()));
    let result = Expr::call(
        Symbol::LIST,
        found
            .into_iter()
            .map(|(z, m)| Expr::call(Symbol::LIST, [Expr::from(z.re()), Expr::from(m)]))
            .collect::<Vec<_>>(),
    );
    ev.step("RealRoots", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
