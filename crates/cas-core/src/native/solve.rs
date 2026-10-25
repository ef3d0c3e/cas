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
//!   is exact when its norm is a square (`√(−4i)` is `√2 − √2 i`).
//! - Higher degrees (CAS-211): binomials `a xⁿ + b` in radicals, parameters
//!   allowed; over ℚ, quadratics in `xᵐ`, cubics and quartics in radicals
//!   where `radicals` finds them; the rest numerically, each root a ball
//!   validated to hold it alone (`cas_num`'s `polynomial_roots`), at the
//!   working precision, as giac's `solve(x^5-x-1,x)` is `[1.16730397826]`.
//!   A factor of a higher degree with parameters, or past
//!   [`MAX_NUMERIC_DEGREE`], leaves the call, with a warning.
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
use cas_num::ball::{ComplexBall, Float, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Cancel, Cancelled, Sparse, SparseRing};

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

/// The highest degree whose roots are sought numerically: each round of
/// the iteration takes the square of the degree in products.
const MAX_NUMERIC_DEGREE: u64 = 64;

/// The bits numeric roots are sought with at most, for the default working
/// precision: 16 times the working precision for a higher one.
const MAX_ROOT_BITS: u64 = 2048;

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
/// warning when one's are not found, or factoring stopped.
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
            d => match higher_roots(ev, u, ring, f, d)? {
                Some(found) => roots.extend(found),
                None => return Ok(None),
            },
        }
    }
    Ok(Some(roots))
}

/// Gives the warning that the roots of the factor `f` are not found, and
/// why.
fn not_found<R: Factoring>(
    ev: &mut Eval<'_>,
    u: &Unknown<'_>,
    f: &Sparse<R::Elem>,
    why: String,
) -> Result<(), EvalError> {
    let factor = ev.eval(&u.p.polynomial_expr(f))?;
    ev.warn(Warning::Native {
        message: why.into(),
        expr: factor,
    });
    Ok(())
}

/// Each of `exprs`, evaluated.
fn evaluated(ev: &mut Eval<'_>, exprs: Vec<Expr>) -> Result<Vec<Expr>, EvalError> {
    exprs.into_iter().map(|e| ev.eval(&e)).collect()
}

/// The roots of the irreducible `f` of degree `d ≥ 3` in the unknown
/// (CAS-DOC-11 §2.4, §3): a binomial in radicals; over ℚ a quadratic in
/// `xᵐ`, a cubic, or a quartic with a rational resolvent root in radicals
/// ([`radicals`]); otherwise numerically ([`numeric_roots`]). `None` after
/// a warning for another factor with parameters, or numeric roots not
/// found.
fn higher_roots<R: Factoring>(
    ev: &mut Eval<'_>,
    u: &Unknown<'_>,
    ring: &SparseRing<R>,
    f: &Sparse<R::Elem>,
    d: u64,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let real = u.real && R::REAL;
    if ring.to_recursive(f, u.k).len() == 2 {
        // a xᵈ + b: xᵈ = −b/a.
        let [b, a] = [0, d].map(|j| u.p.polynomial_expr(&ring.coeff(f, u.k, j)));
        let c = ev.eval(&Expr::call(Symbol::TIMES, [Expr::from(-1), b, a.pow(-1)]))?;
        let sign = if R::REAL { real_sign(ev, &c)? } else { None };
        return evaluated(ev, radicals::binomial(&c, d, sign, real)).map(Some);
    }
    if f.degrees()
        .iter()
        .enumerate()
        .any(|(j, &e)| j != u.k && e > 0)
    {
        not_found::<R>(
            ev,
            u,
            f,
            format!("the roots of a factor of degree {d} with parameters are not found"),
        )?;
        return Ok(None);
    }
    let dense = ring.to_dense(f, u.k).expect("a polynomial in the unknown");
    let coeffs: Vec<Complex> = dense.coeffs().iter().map(R::complex).collect();
    if R::REAL {
        let rationals: Vec<Rational> = coeffs.iter().map(|c| c.re().clone()).collect();
        let cancel: &dyn Cancel = ev.context().interrupt_flag();
        let exact =
            in_radicals(&rationals, real, cancel).map_err(|Cancelled| EvalError::Interrupted)?;
        if let Some(exact) = exact {
            return evaluated(ev, exact).map(Some);
        }
    }
    let found = numeric_roots(ev, &coeffs, real);
    if found.is_none() {
        let why = if d > MAX_NUMERIC_DEGREE {
            format!(
                "the roots of a factor of degree {d} are not sought numerically past degree {MAX_NUMERIC_DEGREE}"
            )
        } else {
            format!("the roots of a factor of degree {d} are not told apart")
        };
        not_found::<R>(ev, u, f, why)?;
    }
    Ok(found)
}

/// The sign of `c`, evaluated, when it is known real: by the assumptions,
/// or by its ball for a number.
fn real_sign(ev: &mut Eval<'_>, c: &Expr) -> Result<Option<Ordering>, EvalError> {
    Ok(match sign_of(ev.context(), c) {
        Sign::Positive => Some(Ordering::Greater),
        Sign::Negative => Some(Ordering::Less),
        Sign::Zero => Some(Ordering::Equal),
        _ => match approximation(ev.context(), c, BITS)?.view() {
            View::Ball(b) => b.sign(),
            _ => None,
        },
    })
}

/// The roots in radicals of the irreducible `Σ cₖ xᵏ` over ℚ of degree 3
/// or more, the real ones when `real`, where [`radicals`] has them: a
/// quadratic in `xᵐ`, a cubic, or a quartic with a rational root of its
/// resolvent.
fn in_radicals(
    c: &[Rational],
    real: bool,
    cancel: &dyn Cancel,
) -> Result<Option<Vec<Expr>>, Cancelled> {
    let d = c.len() - 1;
    // The gcd of the exponents of the terms.
    let m = (1..=d).filter(|&k| !c[k].is_zero()).fold(0, num_gcd);
    if m >= 2 && d == 2 * m {
        let power = u64::try_from(m).expect("a degree in memory");
        return Ok(Some(radicals::quadratic_in_power(
            [&c[0], &c[m], &c[d]],
            power,
            real,
            &Rational::ZERO,
        )));
    }
    match c {
        [d0, d1, d2, d3] => Ok(Some(radicals::cubic([d0, d1, d2, d3], real))),
        [e0, e1, e2, e3, e4] => radicals::quartic([e0, e1, e2, e3, e4], real, cancel),
        _ => Ok(None),
    }
}

/// The gcd of two exponents.
fn num_gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { num_gcd(b, a % b) }
}

/// The roots of the square-free `Σ cₖ xᵏ`, numerically: balls at the
/// working precision, each validated to hold one root; the real ones when
/// `real`, with real coefficients. `None` past [`MAX_NUMERIC_DEGREE`], or
/// when they are not validated within `Limits::number_bits`, or
/// [`MAX_ROOT_BITS`] at the default working precision.
fn numeric_roots(ev: &Eval<'_>, coeffs: &[Complex], real: bool) -> Option<Vec<Expr>> {
    let degree = u64::try_from(coeffs.len() - 1).ok()?;
    if degree > MAX_NUMERIC_DEGREE {
        return None;
    }
    let ctx = ev.context();
    let bits = ctx.working_precision();
    // Clusters of roots need more bits; past a few thousand each round
    // costs more than it is likely to tell.
    let max_bits = ctx
        .limits()
        .number_bits
        .min(bits.saturating_mul(16).max(MAX_ROOT_BITS))
        .max(bits);
    let roots = polynomial_roots(coeffs, bits, max_bits)?;
    Some(
        roots
            .into_iter()
            .filter(|z| !real || z.is_real())
            .map(|z| {
                if z.is_real() {
                    Expr::from(z.re().round_to_accuracy(bits))
                } else {
                    Expr::from(ComplexBall::new(
                        z.re().round_to_accuracy(bits),
                        z.im().round_to_accuracy(bits),
                    ))
                }
            })
            .collect(),
    )
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

mod radicals;
#[cfg(test)]
mod tests;
