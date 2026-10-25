//! `Solve` and `CSolve` (CAS-DOC-11 §1, §2, §4, CAS-210): the solutions of
//! an equation in one unknown, as a list, Xcas's `solve` and `csolve`.
//!
//! - `Equation[l, r]` and `Equal[l, r]` are `l − r = 0`; any other
//!   expression `e` is `e = 0`, with giac's warning. `True`, an `Equal` that
//!   evaluation decided, is solved by every `x`, `False` by none.
//! - `l − r` goes to a rational function of `x` whose other generators, free
//!   of `x`, are parameters (CAS-DOC-8 §2), balls among them; anything else
//!   stays, as `Sin[x] = 1/2` (CAS-213).
//! - Its numerator and denominator are factored over ℚ, or ℚ(i) when a
//!   number is complex, with the parameters (CAS-184, CAS-186): the factors
//!   they share cancel, the denominator's are poles, and the numerator's
//!   that hold `x` give the roots, so `(x − 1)/(x² − 1) = 0` has none.
//! - Degree 1: `−c₀/c₁`. Degree 2: `−b/(2a) ± (s/(2a)) √(Δ/s²)` for
//!   `Δ = b² − 4ac` and `s²` the square part of its content, so that
//!   `x² + y² − 1` gives `±√(1 − y²)`; the square root of a complex number
//!   is exact when its norm is a square (`√(−4i)` is `√2 − √2 i`). A higher
//!   degree leaves the call, with a warning (CAS-211).
//! - `Solve` in the Real domain keeps the real roots. With real
//!   coefficients, a quadratic whose `Δ` the assumptions make negative or
//!   not positive, or its ball negative, has none: giac's
//!   `solve(x^2+a^2,x)` is `[]`. With a complex coefficient every root
//!   stays, as giac's `solve(x-i,x)` is `[i]`. `CSolve`, and `Solve` in the
//!   Complex domain, keep them all.
//! - Each root once. When all are numbers, in giac's order: the real ones
//!   ascending, then the others by real part, a larger imaginary part
//!   before a smaller, the positive one of two conjugates first. With
//!   parameters, in the order of the factors, `+√Δ` first.
//! - `l − r` free of `x`: `{x}` when it is 0, every `x` a solution, as
//!   giac's `solve(x=x,x)` is `[x]`; `{}` otherwise.
//!
//! A system, a list of equations or of unknowns, stays (CAS-212). Each call
//! that solves takes one step, `"Solve"` or `"CSolve"`, which the trace
//! leaves unchecked: the roots' order follows the parameters' values.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of equations: a x² + b x + c in x, Δ = s² d"
)]

use std::cmp::Ordering;

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_num::ball::Float;
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Sparse, SparseRing};

use super::approximation;
use super::polynomial::{Dense, Factoring, convert, factored, root_part, within_dense};
use crate::assume::{Domain, Sign, sign_of};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The bits of the balls that order numeric roots and tell signs.
const BITS: u64 = 128;

/// `Solve[eqn, x]`, Xcas's `solve`: the real solutions in the Real domain.
pub(super) fn solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    solving(ev, e, false)
}

/// `CSolve[eqn, x]`, Xcas's `csolve`: the complex solutions.
pub(super) fn complex_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    solving(ev, e, true)
}

/// The sides of an equation, `Equation[l, r]` or `Equal[l, r]`.
fn sides(e: &Expr) -> Option<(&Expr, &Expr)> {
    match (e.head()?.as_symbol()?, e.args()) {
        (Symbol::EQUATION | Symbol::EQUAL, [l, r]) => Some((l, r)),
        _ => None,
    }
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// [`solve`], or [`complex_solve`] when `complex`.
fn solving(ev: &mut Eval<'_>, e: &Expr, complex: bool) -> Native {
    let [eqn, x] = e.args() else {
        return Ok(None);
    };
    if x.as_symbol().is_none() || is_list(eqn) {
        return Ok(None);
    }
    let roots = match eqn.as_symbol() {
        Some(Symbol::TRUE) => vec![x.clone()],
        Some(Symbol::FALSE) => Vec::new(),
        _ => {
            let f = if let Some((l, r)) = sides(eqn) {
                ev.eval(&(l.clone() - r.clone()))?
            } else {
                ev.warn(Warning::Native {
                    message: "argument is not an equation, solving".into(),
                    expr: Expr::call(Symbol::EQUATION, [eqn.clone(), Expr::from(0)]),
                });
                eqn.clone()
            };
            let real = !complex && ev.context().domain() == Domain::Real;
            match roots(ev, e, &f, x, real)? {
                Some(roots) => roots,
                None => return Ok(None),
            }
        }
    };
    let result = Expr::call(Symbol::LIST, roots);
    ev.step(if complex { "CSolve" } else { "Solve" }, e, &result);
    Ok(Some(result))
}

/// The unknown of an equation converted to rational functions.
struct Unknown<'p> {
    /// The conversion.
    p: &'p Polynomials,
    /// The unknown's variable in it.
    k: usize,
    /// Whether only real roots are wanted.
    real: bool,
}

/// The roots of `f = 0` in `x`, each once, in giac's order; `None` when
/// they are not found here.
fn roots(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
    real: bool,
) -> Result<Option<Vec<Expr>>, EvalError> {
    if !f.contains(x) {
        let every = *f == Expr::from(0);
        return Ok(Some(if every { vec![x.clone()] } else { Vec::new() }));
    }
    let Some(p) = convert(ev, e, std::slice::from_ref(f), true)? else {
        return Ok(None);
    };
    let Some(k) = p.variable(x) else {
        return Ok(None);
    };
    if p.generators()
        .iter()
        .enumerate()
        .any(|(j, g)| j != k && g.contains(x))
        || within_dense(ev, e, &p, Dense::All).is_none()
    {
        return Ok(None);
    }
    let unknown = Unknown { p: &p, k, real };
    let n = p.generators().len();
    let found = match &p.fractions()[0] {
        Fraction::Integers { num, den } => {
            factor_roots(ev, e, &unknown, &SparseRing::new(Integers, n), num, den)?
        }
        Fraction::GaussianIntegers { num, den } => factor_roots(
            ev,
            e,
            &unknown,
            &SparseRing::new(GaussianIntegers, n),
            num,
            den,
        )?,
    };
    let Some(found) = found else {
        return Ok(None);
    };
    let mut unique: Vec<Expr> = Vec::with_capacity(found.len());
    for r in found {
        if !unique.contains(&r) {
            unique.push(r);
        }
    }
    Ok(Some(ordered(ev.context(), unique)))
}

/// The roots of the factors of `num / den` that hold the unknown and are
/// not the denominator's, evaluated, in the factors' order; `None` after a
/// warning when one is of degree 3 or more, or factoring stopped.
fn factor_roots<R: Factoring>(
    ev: &mut Eval<'_>,
    e: &Expr,
    u: &Unknown<'_>,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let Some(powers) = factored(ev, e, ring, num, den)? else {
        return Ok(None);
    };
    let mut roots = Vec::new();
    for (f, m) in &powers.factors {
        // A pole.
        if *m < 0 {
            continue;
        }
        let f = &in_unknown(ring, f, u.k);
        match ring.degree(f, u.k).unwrap_or(0) {
            0 => {}
            1 => {
                let [c0, c1] = [0, 1].map(|d| u.p.polynomial_expr(&ring.coeff(f, u.k, d)));
                let root = Expr::call(Symbol::TIMES, [Expr::from(-1), c0, c1.pow(-1)]);
                roots.push(ev.eval(&root)?);
            }
            2 => roots.extend(quadratic_roots(ev, u, ring, f)?),
            d => {
                let factor = ev.eval(&u.p.polynomial_expr(f))?;
                ev.warn(Warning::Native {
                    message: format!("the roots of a factor of degree {d} are not found").into(),
                    expr: factor,
                });
                return Ok(None);
            }
        }
    }
    Ok(Some(roots))
}

/// `f` times the unit that puts the first coefficient of its leading
/// coefficient in `x_k` in normal form, positive over ℤ: the factors are
/// normal in the first generator of canonical order, as giac writes them
/// (`a − x²`), and their roots come from `x² − a`, `+√a` first.
fn in_unknown<R: Factoring>(
    ring: &SparseRing<R>,
    f: &Sparse<R::Elem>,
    k: usize,
) -> Sparse<R::Elem> {
    let base = ring.base();
    let lead = ring.lead_coeff(f, k);
    let Some(c) = lead.coeffs().first() else {
        return f.clone();
    };
    let unit = base.normal_form(c).1;
    if base.is_one(&unit) {
        f.clone()
    } else {
        ring.scale(f, &unit)
    }
}

/// The roots of the irreducible `a x² + b x + c`, `+√Δ` first: none when
/// only real ones are wanted, the coefficients are real, and `Δ` is not
/// positive by the assumptions or negative by its ball. `Δ` is expanded
/// when its products stay within `Limits::terms`, as estimated, and left
/// as `b² − 4ac` otherwise.
fn quadratic_roots<R: Factoring>(
    ev: &mut Eval<'_>,
    u: &Unknown<'_>,
    ring: &SparseRing<R>,
    f: &Sparse<R::Elem>,
) -> Result<Vec<Expr>, EvalError> {
    let [c, b, a] = [0, 1, 2].map(|d| ring.coeff(f, u.k, d));
    let terms = ev.context().limits().terms;
    let within = |p: &Sparse<R::Elem>, q: &Sparse<R::Elem>| {
        p.term_count().saturating_mul(q.term_count()) <= terms
    };
    let delta = (within(&b, &b) && within(&a, &c)).then(|| {
        let four = ring.base().elem_integer(&Integer::from(4));
        ring.sub(&ring.mul(&b, &b), &ring.scale(&ring.mul(&a, &c), &four))
    });
    let [a, b, c] = [a, b, c].map(|p| u.p.polynomial_expr(&p));
    let delta_expr = match &delta {
        Some(delta) => u.p.polynomial_expr(delta),
        None => b.pow(2) - Expr::from(4) * a.clone() * c,
    };
    if u.real && R::REAL && !may_be_positive(ev, &delta_expr)? {
        return Ok(Vec::new());
    }
    let (s, root) = match &delta {
        Some(delta) => {
            let (s, d) = square_part(ring, delta);
            let root = match d.coeffs() {
                // A number, not 0 as the factor is irreducible.
                [n] if d.degrees().iter().all(|&e| e == 0) => complex_root(&R::complex(n)),
                _ => u.p.polynomial_expr(&d).pow(half()),
            };
            (s, root)
        }
        None => (Integer::ONE, delta_expr.pow(half())),
    };
    let over = (Expr::from(2) * a).pow(-1);
    let center = Expr::call(Symbol::TIMES, [Expr::from(-1), b, over.clone()]);
    let radius = Expr::call(Symbol::TIMES, [Expr::from(s), over, root]);
    Ok(vec![
        ev.eval(&(center.clone() + radius.clone()))?,
        ev.eval(&(center - radius))?,
    ])
}

/// Whether `delta`, evaluated, may be positive: not when the assumptions
/// make it negative or not positive (giac's generic answer, as `−4a²`), or
/// when it is a number whose ball is negative.
fn may_be_positive(ev: &mut Eval<'_>, delta: &Expr) -> Result<bool, EvalError> {
    let delta = ev.eval(delta)?;
    match sign_of(ev.context(), &delta) {
        Sign::Negative | Sign::NonPositive => return Ok(false),
        Sign::Positive | Sign::NonNegative | Sign::Zero => return Ok(true),
        Sign::Nonzero | Sign::Unknown => {}
    }
    let Ok(value) = approximation(ev.context(), &delta, BITS) else {
        return Ok(true);
    };
    Ok(match value.view() {
        View::Ball(b) => b.sign() != Some(Ordering::Less),
        _ => true,
    })
}

/// `s` and `p / s²` for the largest integer `s` whose square divides the
/// content of `p`: of the gcd of the real and imaginary parts of the
/// content over ℤ\[i\].
fn square_part<R: Factoring>(
    ring: &SparseRing<R>,
    p: &Sparse<R::Elem>,
) -> (Integer, Sparse<R::Elem>) {
    let (re, im) = R::parts(&ring.content(p));
    let n = re.numer().gcd(im.numer());
    let (s, _) = root_part(&n.unsigned_abs());
    if s.is_one() {
        return (Integer::ONE, p.clone());
    }
    let s = Integer::from(s);
    let base = ring.base();
    let square = base.elem_integer(&(&s * &s));
    let reduced = ring.from_terms(p.terms().map(|(e, c)| (e, base.div_rem(c, &square).0)));
    (s, reduced)
}

/// `1/2`.
fn half() -> Expr {
    Expr::from(Rational::new(Integer::ONE, Integer::from(2)))
}

/// The principal square root of the Gaussian integer `z`: for `z = a + b i`
/// with `b ≠ 0` and `a² + b² = m²`, `√((m + a)/2) + i sgn(b) √((m − a)/2)`;
/// otherwise `z^(1/2)`, which evaluation writes for a real `z`.
fn complex_root(z: &Complex) -> Expr {
    let (a, b) = (z.re(), z.im());
    if b.is_zero() {
        return Expr::from(z).pow(half());
    }
    let norm = &(a * a) + &(b * b);
    let (m, rest) = norm.numer().sqrt_rem();
    if !rest.is_zero() || !norm.denom().is_one() {
        return Expr::from(z).pow(half());
    }
    let m = Rational::from(m);
    let two = Rational::from(Integer::from(2));
    let re = Expr::from(&(&m + a) / &two).pow(half());
    let im = Expr::from(&(&m - a) / &two).pow(half());
    let sign = if b.is_negative() { -1 } else { 1 };
    re + Expr::call(
        Symbol::TIMES,
        [Expr::from(sign), Expr::from(Complex::I), im],
    )
}

/// The roots in giac's order when all are numbers ([`numeric_order`]);
/// otherwise as they are.
fn ordered(ctx: &Context, roots: Vec<Expr>) -> Vec<Expr> {
    let points: Option<Vec<(Float, Float)>> = roots.iter().map(|r| point(ctx, r)).collect();
    let Some(points) = points else {
        return roots;
    };
    let mut pairs: Vec<((Float, Float), Expr)> = points.into_iter().zip(roots).collect();
    pairs.sort_by(|(p, _), (q, _)| numeric_order(p, q));
    pairs.into_iter().map(|(_, r)| r).collect()
}

/// The midpoints of the real and imaginary parts of a number's ball, the
/// imaginary one 0 when its ball holds 0; `None` for no number.
fn point(ctx: &Context, r: &Expr) -> Option<(Float, Float)> {
    let value = approximation(ctx, r, BITS).ok()?;
    match value.view() {
        View::Ball(b) if b.is_finite() => Some((b.mid().clone(), Float::ZERO)),
        View::ComplexBall(z) if z.is_finite() => {
            let im = if z.im().contains_rational(&Rational::ZERO) {
                Float::ZERO
            } else {
                z.im().mid().clone()
            };
            Some((z.re().mid().clone(), im))
        }
        _ => None,
    }
}

/// giac's order of numeric roots: the real ones first, ascending; then by
/// real part, a larger imaginary part in absolute value first, and of two
/// conjugates the positive one first. A total order on the midpoints.
fn numeric_order((a_re, a_im): &(Float, Float), (b_re, b_im): &(Float, Float)) -> Ordering {
    let zero = Float::ZERO;
    match (*a_im == zero, *b_im == zero) {
        (true, true) => a_re.cmp(b_re),
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        (false, false) => a_re
            .cmp(b_re)
            .then_with(|| b_im.cmp_abs(a_im))
            .then_with(|| b_im.cmp(a_im)),
    }
}

#[cfg(test)]
mod tests;
