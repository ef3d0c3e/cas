//! Roots of polynomials over ℚ (CAS-538, CAS-539): Xcas's `VAS`,
//! `VAS_positive`, `posubLMQ`, `poslbdLMQ`, `sturmseq` and `sturm`,
//! `rationalroot`, `crationalroot` and `complexroot`. The variable is `x`
//! when left out, as in giac; a polynomial of a degree above
//! [`MAX_DEGREE`], or whose coefficients are not rational, leaves the call.
//!
//! - **`RealRootIsolation[p]`** and **`PositiveRootIsolation[p]`**: the real
//!   (positive) roots of the square-free part of `p` isolated by cas-poly's
//!   VAS continued fractions, giac's intervals: each root exact when the
//!   search meets it, otherwise an interval `{a, b}` with exactly one root
//!   inside, ascending: `{{-4, 0}, {1, 3/2}, {3/2, 2}}` for `x³ − 7x + 7`.
//! - **`PositiveRootUpperBound[p]`** and **`PositiveRootLowerBound[p]`**:
//!   the LMQ bounds as giac rounds them, `⌈65 B/64⌉` and a power of two
//!   (`Infinity` without a positive root by Descartes's rule).
//! - **`SturmSequence[p]`**: giac's `sturmseq`. `p` put together as `n/d`
//!   over ℤ, `d`'s leading coefficient positive; for each, a constant, then
//!   the Sturm sequence (`sturm_sequence`, subresultant) of each square-free
//!   factor of odd multiplicity, primitive, the sign of the primitive part
//!   on the factor of highest odd multiplicity (a factor `-1` when there is
//!   none; `x^v` taken out first as giac does, so a factor that is `x`
//!   alone takes it only when `v = 1` and no other is left, CAS-661), the
//!   constant being the leading coefficient over the factors':
//!   `{2, {{1, 0, -2}, {2, 0}, 8}, 1}` for `2x² − 4`. With bounds, giac's
//!   `sturm(p, x, a, b)`: `Σ V(a) − V(b)` over the sequences, the roots and
//!   poles in `(a, b]`; with none, 0 or −1 as the first signs at `a` give
//!   (giac's sign of `p(a)`).
//! - **`RationalRoots[p]`** and **`ComplexRationalRoots[p]`**: the distinct
//!   roots of the linear factors of `p`'s numerator over ℚ (over ℚ(i),
//!   `FactorList` of `i p`, as `ComplexFactor`), giac's `rationalroot` and
//!   `crationalroot`; grouped by multiplicity as giac's square-free
//!   factors, each group's real roots descending, then the others by real
//!   and imaginary part, as giac on small cases (giac's own order follows
//!   its linear factor search).
//! - **`ComplexRoots[p, ε]`** and `ComplexRoots[p, ε, a, b]`: giac's
//!   `complexroot`, every complex root with its multiplicity, `{{z, m},
//!   …}`; a root of a linear factor over ℚ exact, the others balls of `−log₂
//!   ε` bits from the square-free parts (`polynomial_roots`, CAS-211), as
//!   `RealRoots`; ordered by real part, then imaginary part; with the
//!   corners `a` and `b`, only the roots whose balls meet that rectangle.
//!
//! Each call that answers takes one step, its builtin's name.

use std::cmp::Ordering;

use cas_num::algebra::{EuclideanDomain, Rationals};
use cas_num::ball::{Ball, ComplexBall, Float, Round, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::real_roots::{self, Isolated, LowerBound};
use cas_poly::{Dense, DenseRing};

use super::super::integrate::{polynomial, rational_function};
use super::{MAX_DEGREE, call, coefficients, degree, precision, with_variable};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn answer(ev: &mut Eval<'_>, name: &'static str, e: &Expr, result: Expr) -> Native {
    ev.step(name, e, &result);
    Ok(Some(result))
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn small(n: Option<usize>) -> bool {
    n.is_none_or(|n| i64::try_from(n).is_ok_and(|n| n <= MAX_DEGREE))
}

/// `p` over ℚ in `x`: a rational constant, or a polynomial of a degree up
/// to [`MAX_DEGREE`].
fn over_q(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Expr,
    x: &Expr,
) -> Result<Option<Dense<Rational>>, EvalError> {
    if let Some(c) = p.to_rational() {
        return Ok(Some(DenseRing::new(Rationals).constant(c)));
    }
    Ok(polynomial(ev, e, p, x)?.filter(|d| small(d.degree())))
}

/// The positive rational that makes `p / c` integral and primitive.
fn content(coeffs: &[Rational]) -> Rational {
    let (mut numer, mut denom) = (Integer::from(0), Integer::from(1));
    for c in coeffs {
        numer = numer.gcd(c.numer());
        denom = denom.lcm(c.denom());
    }
    Rational::new(numer, denom)
}

fn isolated(entries: Vec<Isolated>) -> Expr {
    list(entries.into_iter().map(|i| match i {
        Isolated::Root(r) => Expr::from(r),
        Isolated::Interval(a, b) => list([Expr::from(a), Expr::from(b)]),
    }))
}

/// `RealRootIsolation[p]` and `PositiveRootIsolation[p]`, and the bounds.
fn vas_family(ev: &mut Eval<'_>, e: &Expr, name: &'static str) -> Native {
    let Some((p, x)) = with_variable(e.args()) else {
        return Ok(None);
    };
    let Some(d) = over_q(ev, e, p, &x)? else {
        return Ok(None);
    };
    let result = match name {
        "RealRootIsolation" => isolated(real_roots::vas(&d)),
        "PositiveRootIsolation" => isolated(real_roots::vas_positive(&d)),
        "PositiveRootUpperBound" => Expr::from(real_roots::posub_lmq(&d)),
        _ => match real_roots::poslbd_lmq(&d) {
            LowerBound::Infinity => Expr::from(Symbol::INFINITY),
            bound => Expr::from(bound.value().expect("finite")),
        },
    };
    answer(ev, name, e, result)
}

/// `RealRootIsolation[p]`, giac's `VAS`.
pub(in crate::native) fn real_root_isolation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vas_family(ev, e, "RealRootIsolation")
}

/// `PositiveRootIsolation[p]`, giac's `VAS_positive`.
pub(in crate::native) fn positive_root_isolation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vas_family(ev, e, "PositiveRootIsolation")
}

/// `PositiveRootUpperBound[p]`, giac's `posubLMQ`.
pub(in crate::native) fn positive_root_upper_bound(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vas_family(ev, e, "PositiveRootUpperBound")
}

/// `PositiveRootLowerBound[p]`, giac's `poslbdLMQ`.
pub(in crate::native) fn positive_root_lower_bound(ev: &mut Eval<'_>, e: &Expr) -> Native {
    vas_family(ev, e, "PositiveRootLowerBound")
}

/// A part of giac's `sturm`: a constant, or a Sturm sequence.
enum Part {
    Constant(Rational),
    Sequence(Vec<Dense<Rational>>),
}

/// A polynomial as giac's list, the leading coefficient first; a constant
/// bare.
fn coefficient_list(p: &Dense<Rational>) -> Expr {
    match p.coeffs() {
        [] => Expr::from(0),
        [c] => Expr::from(c.clone()),
        cs => list(cs.iter().rev().cloned().map(Expr::from)),
    }
}

impl Part {
    fn to_expr(&self) -> Expr {
        match self {
            Part::Constant(c) => Expr::from(c.clone()),
            Part::Sequence(s) => list(s.iter().map(coefficient_list)),
        }
    }
}

/// A numerator and a denominator.
type Fraction = (Dense<Rational>, Dense<Rational>);

/// `p` as `n/d` over ℤ, coprime contents, `d`'s leading coefficient
/// positive, as giac's `e2r`.
fn fraction(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Expr,
    x: &Expr,
) -> Result<Option<Fraction>, EvalError> {
    let qx = DenseRing::new(Rationals);
    let (n, d) = if let Some(c) = p.to_rational() {
        (qx.constant(c), qx.constant(Rational::ONE))
    } else {
        let Some((n, d)) = rational_function(ev, e, p, x)? else {
            return Ok(None);
        };
        if !small(n.degree()) || !small(d.degree()) {
            return Ok(None);
        }
        (n, d)
    };
    let all: Vec<Rational> = n.coeffs().iter().chain(d.coeffs()).cloned().collect();
    let mut s = content(&all).recip();
    if d.lead().is_some_and(Rational::is_negative) {
        s = -s;
    }
    Ok(Some((qx.scale(&n, &s), qx.scale(&d, &s))))
}

/// giac's `sturm` of one side of the fraction, integral.
fn parts(n: &Dense<Rational>) -> Vec<Part> {
    let qx = DenseRing::new(Rationals);
    let Some(lead) = n.lead().filter(|_| n.degree().is_some_and(|k| k > 0)) else {
        return vec![Part::Constant(
            n.coeffs().first().cloned().unwrap_or(Rational::ZERO),
        )];
    };
    let primitive = qx.scale(n, &content(n.coeffs()).recip());
    let mut found: Vec<(Dense<Rational>, u64)> = qx
        .squarefree(&primitive)
        .factors()
        .iter()
        .filter(|(f, m)| m % 2 == 1 && f.degree().is_some_and(|k| k > 0))
        .map(|(f, m)| {
            let mut s = content(f.coeffs()).recip();
            if f.lead().is_some_and(Rational::is_negative) {
                s = -s;
            }
            (qx.scale(f, &s), *m)
        })
        .collect();
    found.sort_by_key(|(_, m)| *m);
    // giac takes `x^v` out before its square-free factorisation (CAS-661):
    // the sign goes on the rest's factor of highest odd multiplicity, so a
    // factor that is `x` alone takes it only when `v = 1` and the rest has
    // none; otherwise it stays apart, `-1` first.
    // A factor that is `x` alone is the multiplicity-`v` one.
    let x_alone = |f: &Dense<Rational>| f.coeffs() == [Rational::ZERO, Rational::ONE];
    let target = found
        .iter()
        .rposition(|(f, _)| !x_alone(f))
        .or_else(|| found.iter().position(|(f, m)| *m == 1 && x_alone(f)));
    let mut odd: Vec<Dense<Rational>> = found.into_iter().map(|(f, _)| f).collect();
    if primitive.lead().is_some_and(Rational::is_negative) {
        match target {
            Some(i) => odd[i] = qx.scale(&odd[i], &-Rational::ONE),
            None => odd.insert(0, qx.constant(-Rational::ONE)),
        }
    }
    let leads = odd
        .iter()
        .fold(Rational::ONE, |a, f| &a * f.lead().expect("not 0"));
    let mut result = vec![Part::Constant(lead / &leads)];
    result.extend(odd.iter().map(|f| Part::Sequence(qx.sturm_sequence(f))));
    result
}

/// The sign changes of a sequence at `a`, zeros left out, and whether its
/// first sign there is negative.
fn changes(sequence: &[Dense<Rational>], a: &Rational) -> (i64, bool) {
    let qx = DenseRing::new(Rationals);
    let (mut count, mut first, mut last) = (0, None, None);
    for q in sequence {
        let value = qx.eval(q, a);
        if value.is_zero() {
            continue;
        }
        let negative = value.is_negative();
        first.get_or_insert(negative);
        if last.is_some_and(|l| l != negative) {
            count += 1;
        }
        last = Some(negative);
    }
    (count, first == Some(true))
}

/// giac's `sturmab`: the roots and poles in `(a, b]` by the sequences; with
/// none, `−1` when the first signs at `a` multiply to a negative, else 0.
fn count(
    n: &Dense<Rational>,
    d: &Dense<Rational>,
    parts: &[Part],
    a: Rational,
    b: Rational,
) -> i64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    let qx = DenseRing::new(Rationals);
    if a == b {
        let (na, da) = (qx.eval(n, &a), qx.eval(d, &a));
        if !na.is_zero() && !da.is_zero() {
            return if na.is_negative() == da.is_negative() {
                0
            } else {
                -1
            };
        }
    }
    let (mut found, mut negative) = (0, false);
    for part in parts {
        match part {
            Part::Sequence(s) => {
                let (at_a, first) = changes(s, &a);
                found += at_a - changes(s, &b).0;
                negative ^= first;
            }
            Part::Constant(c) => negative ^= c.is_negative(),
        }
    }
    if found != 0 {
        found
    } else {
        -i64::from(negative)
    }
}

/// `SturmSequence[p]`, `[p, x]` and `[p, x, a, b]`: giac's `sturmseq` and
/// `sturm`.
pub(in crate::native) fn sturm_sequence(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, x, bounds) = match e.args() {
        [p] => (p, Expr::symbol("x"), None),
        [p, x] if x.as_symbol().is_some() => (p, x.clone(), None),
        [p, x, a, b] if x.as_symbol().is_some() => {
            let (Some(a), Some(b)) = (a.to_rational(), b.to_rational()) else {
                return Ok(None);
            };
            (p, x.clone(), Some((a, b)))
        }
        _ => return Ok(None),
    };
    let Some((n, d)) = fraction(ev, e, p, &x)? else {
        return Ok(None);
    };
    let mut all = parts(&n);
    all.extend(parts(&d));
    let result = match bounds {
        None => list(all.iter().map(Part::to_expr)),
        Some(_) if n.is_zero() => return Ok(None),
        Some((a, b)) => Expr::from(count(&n, &d, &all, a, b)),
    };
    answer(ev, "SturmSequence", e, result)
}

/// The distinct roots of the linear factors of `p`'s numerator, over ℚ(i)
/// with `complex`, each with the multiplicity of its factor; `None` when
/// the numerator is no polynomial in `x` over ℚ (ℚ(i) with `complex`).
fn linear_roots(
    ev: &mut Eval<'_>,
    p: &Expr,
    x: &Expr,
    complex: bool,
) -> Result<Option<Vec<(Complex, i64)>>, EvalError> {
    let together = ev.eval(&call(Symbol::TOGETHER, vec![p.clone()]))?;
    let numerator = ev.eval(&call(Symbol::NUMERATOR, vec![together]))?;
    if numerator.to_complex().is_some() {
        return Ok(Some(Vec::new()));
    }
    let Some(cs) = coefficients(ev, &numerator, x)? else {
        return Ok(None);
    };
    let exact = |c: &Expr| match c.to_complex() {
        Some(z) => complex || z.im().is_zero(),
        None => false,
    };
    if !small(cs.len().checked_sub(1)) || !cs.iter().all(exact) {
        return Ok(None);
    }
    let i = Expr::from(Complex::new(Rational::ZERO, Rational::ONE));
    let target = if complex { i * numerator } else { numerator };
    let listed = ev.eval(&call(Symbol::FACTOR_LIST, vec![target]))?;
    let mut roots = Vec::new();
    for pair in listed.args().chunks(2) {
        let [g, k] = pair else {
            return Ok(None);
        };
        let Some(k) = k.as_i64().filter(|k| *k >= 1) else {
            return Ok(None);
        };
        if !g.contains(x) || degree(ev, g, x)? != Some(1) {
            continue;
        }
        let Some(cs) = coefficients(ev, g, x)? else {
            return Ok(None);
        };
        let [Some(c0), Some(c1)] = [cs[0].to_complex(), cs[1].to_complex()] else {
            return Ok(None);
        };
        let root = &(-&c0) / &c1;
        if complex || root.im().is_zero() {
            roots.push((root, k));
        }
    }
    Ok(Some(roots))
}

/// giac's order of the rational roots, as far as it is cheap: by
/// multiplicity, then the real roots descending, then the others by real
/// and imaginary part.
fn giac_order(a: &(Complex, i64), b: &(Complex, i64)) -> Ordering {
    let real = |z: &Complex| z.im().is_zero();
    a.1.cmp(&b.1)
        .then(real(&b.0).cmp(&real(&a.0)))
        .then_with(|| {
            if real(&a.0) {
                b.0.re().cmp(a.0.re())
            } else {
                (a.0.re(), a.0.im()).cmp(&(b.0.re(), b.0.im()))
            }
        })
}

fn rational_roots_of(ev: &mut Eval<'_>, e: &Expr, complex: bool) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let Some(mut roots) = linear_roots(ev, p, &Expr::symbol("x"), complex)? else {
        return Ok(None);
    };
    roots.sort_by(giac_order);
    let result = list(roots.into_iter().map(|(z, _)| Expr::from(z)));
    let name = if complex {
        "ComplexRationalRoots"
    } else {
        "RationalRoots"
    };
    answer(ev, name, e, result)
}

/// `RationalRoots[p]`, giac's `rationalroot`.
pub(in crate::native) fn rational_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rational_roots_of(ev, e, false)
}

/// `ComplexRationalRoots[p]`, giac's `crationalroot`.
pub(in crate::native) fn complex_rational_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    rational_roots_of(ev, e, true)
}

/// A rational's nearest double, for the order.
fn approximate(q: &Rational) -> f64 {
    Float::from_rational(q, 64, Round::Nearest).0.to_f64()
}

/// Whether the ball `z` meets the rectangle with the corners `a` and `b`.
fn meets(z: &ComplexBall, a: &Complex, b: &Complex) -> bool {
    let within = |x: &Ball, u: &Rational, v: &Rational| {
        let (lo, hi) = if u > v { (v, u) } else { (u, v) };
        x.cmp(&Ball::from_rational(lo, 64)) != Some(Ordering::Less)
            && x.cmp(&Ball::from_rational(hi, 64)) != Some(Ordering::Greater)
    };
    within(z.re(), a.re(), b.re()) && within(z.im(), a.im(), b.im())
}

/// `ComplexRoots[p, ε]` and `ComplexRoots[p, ε, a, b]`, giac's
/// `complexroot`.
pub(in crate::native) fn complex_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, eps, corners) = match e.args() {
        [p, eps] => (p, eps, None),
        [p, eps, a, b] => {
            let (Some(a), Some(b)) = (a.to_complex(), b.to_complex()) else {
                return Ok(None);
            };
            (p, eps, Some((a, b)))
        }
        _ => return Ok(None),
    };
    let Some(bits) = precision(eps) else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let Some(dense) = over_q(ev, e, p, &x)? else {
        return Ok(None);
    };
    if dense.degree().is_none_or(|n| n == 0) {
        return Ok(None);
    }
    let Some(exact) = linear_roots(ev, p, &x, false)? else {
        return Ok(None);
    };
    // (ball, value, multiplicity, order key)
    let qx = DenseRing::new(Rationals);
    let mut rest = dense;
    let mut found: Vec<(ComplexBall, Expr, i64, (f64, f64))> = Vec::new();
    for (z, m) in exact {
        let r = z.re().clone();
        let linear = qx.poly(vec![-r.clone(), Rational::ONE]);
        for _ in 0..m {
            rest = qx.div_rem(&rest, &linear).0;
        }
        let ball = ComplexBall::from_complex(&z, bits);
        found.push((ball, Expr::from(r.clone()), m, (approximate(&r), 0.0)));
    }
    for (factor, m) in qx.squarefree(&rest).factors() {
        if factor.degree().is_none_or(|n| n == 0) {
            continue;
        }
        let cs: Vec<Complex> = factor.coeffs().iter().cloned().map(Complex::from).collect();
        let Some(roots) = polynomial_roots(&cs, bits, 64 * bits.max(64)) else {
            return Ok(None);
        };
        let m = i64::try_from(*m).expect("a multiplicity in memory");
        for z in roots {
            let key = (z.re().mid().to_f64(), z.im().mid().to_f64());
            let value = if z.is_real() {
                Expr::from(z.re().clone())
            } else {
                Expr::from(z.clone())
            };
            found.push((z, value, m, key));
        }
    }
    if let Some((a, b)) = &corners {
        found.retain(|(z, ..)| meets(z, a, b));
    }
    found.sort_by(|u, v| {
        u.3.0
            .partial_cmp(&v.3.0)
            .unwrap_or(Ordering::Equal)
            .then(u.3.1.partial_cmp(&v.3.1).unwrap_or(Ordering::Equal))
    });
    let result = list(
        found
            .into_iter()
            .map(|(_, value, m, _)| list([value, Expr::from(m)])),
    );
    answer(ev, "ComplexRoots", e, result)
}
