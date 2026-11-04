//! Polynomial tools II (CAS-DOC-48 §2, CAS-482): Xcas's `proot`, `pcoeff`,
//! `peval` and `horner`, `ptayl` and `froot`. A polynomial is given as an
//! expression in `x` or as giac's list of coefficients, the leading one
//! first.
//!
//! - **`PolynomialRoots[p]`**: every complex root of `p` over ℚ, each as
//!   often as its multiplicity, by the square-free parts: the roots of a
//!   linear part and of a quadratic part whose discriminant is a square or
//!   minus one exactly (giac's `proot([1,0,1])` is `[-i,i]`), the others
//!   balls by cas-num's certified isolation (CAS-211); ordered by real part,
//!   then imaginary part.
//! - **`FromRoots[l]`**: the coefficients of `Π (x − rᵢ)`, the leading one
//!   first, giac's `poly1[1,-6,11,-6]` for `[1,2,3]`.
//! - **`HornerEval[p, a]`**: `p(a)` by Horner's scheme.
//! - **`TaylorShift[p, a]`**: `p(x + a)`, a list for a list (giac's
//!   `[1,6,14,13]` for `[1,0,2,1]` at 2), expanded for an expression.
//! - **`RootMultiplicities[f]`**: the roots of the numerator and the
//!   denominator of `f` over ℚ, each followed by its multiplicity, negated
//!   for a pole: giac's `froot(x^2*(x-3)^3)` is `[0,2,3,3]`. Each
//!   irreducible factor's roots are `Solve`'s, rational ones first,
//!   ascending; a factor whose roots `Solve` leaves leaves the call.
//!
//! Bézout, Chinese remainders and series (§3, CAS-483), over ℚ:
//!
//! - **`PolynomialExtendedGCD[a, b, x]`**, `egcd` and `gcdex`: `{u, v, d}`
//!   with `u a + v b = d`, by the extended Euclid, `d` the last remainder
//!   as it comes (giac's `[1,-x+1,2]` for `x² + 1` and `x + 1`).
//! - **`PolynomialBezout[a, b, c, x]`**, `abcuv`: `{u, v}` with `u a + v b =
//!   c` and `deg u < deg b`, for `gcd(a, b)` dividing `c`.
//! - **`PolynomialChineseRemainder[{a, p}, {b, q}]`**, `chinrem`: `{c, p
//!   q}`, `c ≡ a mod p` and `c ≡ b mod q`, for coprime `p` and `q`: giac's
//!   `[-x+2,x^2+x]` for `[2,x]` and `[3,x+1]`.
//! - **`SeriesDivision[a, b, n]`**, `divpc`: the power series of `a/b` at
//!   0 up to `xⁿ`, `b(0) ≠ 0`.
//! - **`PadeApproximant[f, x, n, p]`**, `pade`: `P/Q` with `deg P < p`,
//!   `deg Q ≤ n − p + 1`, `Q(0) = 1` and `f Q − P = O(xⁿ⁺¹)`, from `f`'s
//!   series by the linear system of the coefficients of `Q`; a singular
//!   system leaves the call.
//!
//! Factoring and forms (§4, CAS-484):
//!
//! - **`ComplexFactor[p]`**, `cfactor`: `Factor` over ℚ(i), as it factors
//!   when a coefficient is complex: `Factor[i p]/i`, giac's
//!   `(x+i)*(x-i)` for `x² + 1`.
//! - **`CompleteSquare[p, x]`**, `canonical_form`: `a (x + b/(2a))² + c −
//!   b²/(4a)` for `a x² + b x + c`, giac's `(x-3)^2-7`.
//! - **`DigitPolynomial[n, b, x]`**, `genpoly`: `Σ dᵢ xⁱ` with `n = Σ dᵢ bⁱ`,
//!   the balanced digits in `(−b/2, b/2]`, giac's `x^2-1` for 15 in base 4.
//! - **`BerlekampMassey[s]`**, `berlekamp_massey`: the shortest linear
//!   recurrence of a sequence of even length over ℚ, its connection
//!   polynomial `1 + c₁ x + ⋯ + c_L x^L` as a list from the constant up,
//!   giac's `[1,-1,-1]` for the Fibonacci numbers.

use std::cmp::Ordering;

use cas_num::algebra::{EuclideanDomain, Rationals, Ring};
use cas_num::ball::{ComplexBall, Float, Round, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Dense, DenseRing};

use super::super::integrate::polynomial;
use super::{call, coefficients};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// The coefficients of `p` from `x⁰` up: of a list (leading first), or of
/// an expression in `x`.
fn ascending(ev: &mut Eval<'_>, p: &Expr) -> Result<Option<Vec<Expr>>, EvalError> {
    if is_list(p) {
        return Ok(Some(p.args().iter().rev().cloned().collect()));
    }
    coefficients(ev, p, &Expr::symbol("x"))
}

/// `Σ cₖ xᵏ` from ascending coefficients.
fn expression(cs: &[Expr], x: &Expr) -> Expr {
    Expr::call(
        Symbol::PLUS,
        cs.iter()
            .enumerate()
            .map(|(k, c)| c.clone() * x.pow(i64::try_from(k).expect("small")))
            .collect::<Vec<_>>(),
    )
}

/// `p` over ℚ as a dense polynomial.
fn dense(ev: &mut Eval<'_>, e: &Expr, p: &Expr) -> Result<Option<Dense<Rational>>, EvalError> {
    let x = Expr::symbol("x");
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    polynomial(ev, e, &expression(&cs, &x), &x)
}

/// The exact square root of a rational, when it has one.
fn rational_sqrt(q: &Rational) -> Option<Rational> {
    if q.is_negative() {
        return None;
    }
    let (n, d) = (q.numer().isqrt(), q.denom().isqrt());
    let r = Rational::new(n, d);
    (&r * &r == *q).then_some(r)
}

/// A root, exact or a ball, with its real and imaginary parts' midpoints
/// for the order.
struct Root {
    value: Expr,
    key: (f64, f64),
}

/// A rational's nearest double, for the order.
fn approximate(q: &Rational) -> f64 {
    Float::from_rational(q, 64, Round::Nearest).0.to_f64()
}

fn exact(z: Complex) -> Root {
    let key = (approximate(z.re()), approximate(z.im()));
    Root {
        value: Expr::from(z),
        key,
    }
}

/// The roots of a square-free part over ℚ, `None` when they are not
/// isolated.
fn part_roots(factor: &Dense<Rational>, bits: u64) -> Option<Vec<Root>> {
    let cs = factor.coeffs();
    match cs {
        [c, b] => return Some(vec![exact(Complex::from(-(c / b)))]),
        [c, b, a] => {
            // (−b ± √Δ)/(2a), exact when Δ or −Δ is a square.
            let delta = &(b * b) - &(&Rational::from(Integer::from(4)) * &(a * c));
            let two_a = a + a;
            let centre = -(b / &two_a);
            if let Some(r) = rational_sqrt(&delta) {
                let r = &r / &two_a;
                return Some(vec![
                    exact(Complex::from(&centre - &r)),
                    exact(Complex::from(&centre + &r)),
                ]);
            }
            if let Some(r) = rational_sqrt(&-delta) {
                let r = (&r / &two_a).abs();
                return Some(vec![
                    exact(Complex::new(centre.clone(), -r.clone())),
                    exact(Complex::new(centre, r)),
                ]);
            }
        }
        _ => {}
    }
    let complex: Vec<Complex> = cs.iter().cloned().map(Complex::from).collect();
    let roots = polynomial_roots(&complex, bits, 64 * bits.max(64))?;
    Some(
        roots
            .into_iter()
            .map(|z: ComplexBall| {
                let key = (z.re().mid().to_f64(), z.im().mid().to_f64());
                let value = if z.is_real() {
                    Expr::from(z.re().clone())
                } else {
                    Expr::from(z)
                };
                Root { value, key }
            })
            .collect(),
    )
}

/// `PolynomialRoots[p]`.
pub(in crate::native) fn polynomial_roots_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let Some(d) = dense(ev, e, p)? else {
        return Ok(None);
    };
    if d.degree().is_none_or(|n| n == 0) {
        return Ok(None);
    }
    let bits = ev.context().working_precision();
    let qx = DenseRing::new(Rationals);
    let mut roots: Vec<Root> = Vec::new();
    for (factor, m) in qx.squarefree(&d).factors() {
        if factor.degree().is_none_or(|n| n == 0) {
            continue;
        }
        let Some(found) = part_roots(factor, bits) else {
            return Ok(None);
        };
        for r in found {
            for _ in 0..*m {
                roots.push(Root {
                    value: r.value.clone(),
                    key: r.key,
                });
            }
        }
    }
    roots.sort_by(|a, b| {
        a.key
            .0
            .partial_cmp(&b.key.0)
            .unwrap_or(Ordering::Equal)
            .then(a.key.1.partial_cmp(&b.key.1).unwrap_or(Ordering::Equal))
    });
    let result = Expr::call(
        Symbol::LIST,
        roots.into_iter().map(|r| r.value).collect::<Vec<_>>(),
    );
    ev.step("proot", e, &result);
    Ok(Some(result))
}

/// `FromRoots[l]`.
pub(in crate::native) fn from_roots(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    let x = Expr::symbol("x");
    let product = Expr::call(
        Symbol::TIMES,
        l.args()
            .iter()
            .map(|r| x.clone() - r.clone())
            .collect::<Vec<_>>(),
    );
    let Some(cs) = coefficients(ev, &product, &x)? else {
        return Ok(None);
    };
    let result = Expr::call(Symbol::LIST, cs.into_iter().rev().collect::<Vec<_>>());
    ev.step("pcoeff", e, &result);
    Ok(Some(result))
}

/// `HornerEval[p, a]`.
pub(in crate::native) fn horner_eval(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, a] = e.args() else {
        return Ok(None);
    };
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    let value = cs
        .iter()
        .rev()
        .fold(Expr::from(0), |acc, c| acc * a.clone() + c.clone());
    let result = ev.eval(&value)?;
    ev.step("horner", e, &result);
    Ok(Some(result))
}

/// `TaylorShift[p, a]`.
pub(in crate::native) fn taylor_shift(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p, a] = e.args() else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let Some(cs) = ascending(ev, p)? else {
        return Ok(None);
    };
    let shifted = expression(&cs, &x).substitute(&[(x.clone(), x.clone() + a.clone())]);
    let expanded = ev.eval(&call(Symbol::EXPAND, vec![shifted]))?;
    let result = if is_list(p) {
        let Some(new) = coefficients(ev, &expanded, &x)? else {
            return Ok(None);
        };
        Expr::call(Symbol::LIST, new.into_iter().rev().collect::<Vec<_>>())
    } else {
        expanded
    };
    ev.step("ptayl", e, &result);
    Ok(Some(result))
}

/// `RootMultiplicities[f]`.
pub(in crate::native) fn root_multiplicities(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f] = e.args() else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let together = ev.eval(&call(Symbol::TOGETHER, vec![f.clone()]))?;
    let numerator = ev.eval(&call(Symbol::NUMERATOR, vec![together.clone()]))?;
    let denominator = ev.eval(&call(Symbol::DENOMINATOR, vec![together]))?;
    let mut pairs: Vec<(Expr, i64, Option<Rational>)> = Vec::new();
    for (part, sign) in [(numerator, 1), (denominator, -1)] {
        let factors = ev.eval(&call(Symbol::FACTOR_LIST, vec![part]))?;
        if !is_list(&factors) {
            return Ok(None);
        }
        for item in factors.args() {
            let [g, m] = item.args() else {
                return Ok(None);
            };
            if !g.contains(&x) {
                continue;
            }
            let Some(m) = m.as_i64() else {
                return Ok(None);
            };
            let roots = ev.eval(&call(
                Symbol::SOLVE,
                vec![
                    Expr::call(Symbol::EQUATION, [g.clone(), Expr::from(0)]),
                    x.clone(),
                ],
            ))?;
            if !is_list(&roots) {
                return Ok(None);
            }
            for r in roots.args() {
                pairs.push((r.clone(), sign * m, r.to_rational()));
            }
        }
    }
    // Rational roots first, ascending.
    pairs.sort_by(|a, b| match (&a.2, &b.2) {
        (Some(p), Some(q)) => p.cmp(q),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    });
    let result = Expr::call(
        Symbol::LIST,
        pairs
            .into_iter()
            .flat_map(|(r, m, _)| [r, Expr::from(m)])
            .collect::<Vec<_>>(),
    );
    ev.step("froot", e, &result);
    Ok(Some(result))
}

type Poly = Dense<Rational>;

/// A polynomial over ℚ as an expression in `x`.
fn to_expr(d: &Poly, x: &Expr) -> Expr {
    let cs: Vec<Expr> = d.coeffs().iter().cloned().map(Expr::from).collect();
    expression(&cs, x)
}

/// The variable of the trailing argument, `x` when left out.
fn variable(args: &[Expr], n: usize) -> Option<Expr> {
    match args.len() {
        k if k == n => Some(Expr::symbol("x")),
        k if k == n + 1 && args[n].as_symbol().is_some() => Some(args[n].clone()),
        _ => None,
    }
}

/// `p` over ℚ in `x`.
fn over_q(ev: &mut Eval<'_>, e: &Expr, p: &Expr, x: &Expr) -> Result<Option<Poly>, EvalError> {
    polynomial(ev, e, p, x)
}

/// `(u, v, d)` with `u a + v b = d`, the extended Euclid.
fn extended(ring: &DenseRing<Rationals>, a: &Poly, b: &Poly) -> (Poly, Poly, Poly) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1) = (ring.one(), ring.zero());
    let (mut t0, mut t1) = (ring.zero(), ring.one());
    while !r1.is_zero() {
        let (q, r) = ring.div_rem(&r0, &r1);
        (r0, r1) = (r1, r);
        let s = ring.sub(&s0, &ring.mul(&q, &s1));
        (s0, s1) = (s1, s);
        let t = ring.sub(&t0, &ring.mul(&q, &t1));
        (t0, t1) = (t1, t);
    }
    (s0, t0, r0)
}

fn answer(ev: &mut Eval<'_>, name: &'static str, e: &Expr, result: &Expr) -> Native {
    let result = ev.eval(result)?;
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// `PolynomialExtendedGCD[a, b, x]`.
pub(in crate::native) fn polynomial_extended_gcd(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(x) = variable(e.args(), 2) else {
        return Ok(None);
    };
    let (Some(a), Some(b)) = (
        over_q(ev, e, &e.args()[0], &x)?,
        over_q(ev, e, &e.args()[1], &x)?,
    ) else {
        return Ok(None);
    };
    let ring = DenseRing::new(Rationals);
    let (u, v, d) = extended(&ring, &a, &b);
    let result = Expr::call(
        Symbol::LIST,
        [to_expr(&u, &x), to_expr(&v, &x), to_expr(&d, &x)],
    );
    answer(ev, "egcd", e, &result)
}

/// `PolynomialBezout[a, b, c, x]`.
pub(in crate::native) fn polynomial_bezout(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(x) = variable(e.args(), 3) else {
        return Ok(None);
    };
    let args = e.args();
    let (Some(a), Some(b), Some(c)) = (
        over_q(ev, e, &args[0], &x)?,
        over_q(ev, e, &args[1], &x)?,
        over_q(ev, e, &args[2], &x)?,
    ) else {
        return Ok(None);
    };
    if b.is_zero() {
        return Ok(None);
    }
    let ring = DenseRing::new(Rationals);
    let (u0, _, d) = extended(&ring, &a, &b);
    let (k, rest) = ring.div_rem(&c, &d);
    if !rest.is_zero() {
        return Ok(None);
    }
    let u = ring.div_rem(&ring.mul(&u0, &k), &b).1;
    let (v, exact) = ring.div_rem(&ring.sub(&c, &ring.mul(&u, &a)), &b);
    if !exact.is_zero() {
        return Ok(None);
    }
    let result = Expr::call(Symbol::LIST, [to_expr(&u, &x), to_expr(&v, &x)]);
    answer(ev, "abcuv", e, &result)
}

/// `PolynomialChineseRemainder[{a, p}, {b, q}]`.
pub(in crate::native) fn polynomial_chinese_remainder(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [first, second] = e.args() else {
        return Ok(None);
    };
    let ([a, p], [b, q]) = (first.args(), second.args()) else {
        return Ok(None);
    };
    if !is_list(first) || !is_list(second) {
        return Ok(None);
    }
    let x = Expr::symbol("x");
    let (Some(a), Some(p), Some(b), Some(q)) = (
        over_q(ev, e, a, &x)?,
        over_q(ev, e, p, &x)?,
        over_q(ev, e, b, &x)?,
        over_q(ev, e, q, &x)?,
    ) else {
        return Ok(None);
    };
    let ring = DenseRing::new(Rationals);
    let (s, _, d) = extended(&ring, &p, &q);
    if d.degree() != Some(0) {
        return Ok(None);
    }
    // s p = d modulo q: p's inverse is s/d.
    let inverse = ring.div_rem(&s, &d).0;
    let t = ring.div_rem(&ring.mul(&ring.sub(&b, &a), &inverse), &q).1;
    let c = ring.add(&a, &ring.mul(&p, &t));
    let modulus = ring.mul(&p, &q);
    let result = Expr::call(Symbol::LIST, [to_expr(&c, &x), to_expr(&modulus, &x)]);
    answer(ev, "chinrem", e, &result)
}

/// `SeriesDivision[a, b, n]`.
pub(in crate::native) fn series_division(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b, n] = e.args() else {
        return Ok(None);
    };
    let Some(n) = n.as_i64().and_then(|n| usize::try_from(n).ok()) else {
        return Ok(None);
    };
    let x = Expr::symbol("x");
    let (Some(a), Some(b)) = (over_q(ev, e, a, &x)?, over_q(ev, e, b, &x)?) else {
        return Ok(None);
    };
    let at = |d: &Poly, k: usize| d.coeffs().get(k).cloned().unwrap_or(Rational::ZERO);
    let b0 = at(&b, 0);
    if b0.is_zero() {
        return Ok(None);
    }
    let mut q: Vec<Rational> = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let mut v = at(&a, k);
        for j in 1..=k {
            v = &v - &(&at(&b, j) * &q[k - j]);
        }
        q.push(&v / &b0);
    }
    let result = to_expr(&DenseRing::new(Rationals).poly(q), &x);
    answer(ev, "divpc", e, &result)
}

/// `PadeApproximant[f, x, n, p]`.
pub(in crate::native) fn pade_approximant(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, x, n, p] = e.args() else {
        return Ok(None);
    };
    let (Some(n), Some(p)) = (
        n.as_i64().and_then(|n| usize::try_from(n).ok()),
        p.as_i64().and_then(|p| usize::try_from(p).ok()),
    ) else {
        return Ok(None);
    };
    if x.as_symbol().is_none() || p == 0 || p > n + 1 {
        return Ok(None);
    }
    let m = n + 1 - p;
    // The series of f to xⁿ, its order term dropped.
    let series = ev.eval(&Expr::call(
        Symbol::SERIES,
        [
            f.clone(),
            x.clone(),
            Expr::from(0),
            Expr::from(i64::try_from(n).expect("small")),
        ],
    ))?;
    let order = Expr::call(Symbol::ORDER_SIZE, [x.clone()]);
    let truncated = ev.eval(&series.substitute(&[(order, Expr::from(0))]))?;
    let Some(mut c) = coefficients(ev, &truncated, x)? else {
        return Ok(None);
    };
    c.resize(n + 1, Expr::from(0));
    let coefficient = |k: isize| -> Expr {
        usize::try_from(k)
            .ok()
            .and_then(|k| c.get(k).cloned())
            .unwrap_or_else(|| Expr::from(0))
    };
    // Σ_{j=0..m} q_j c_{k−j} = 0 for k = p..n, q_0 = 1: m equations in
    // q_1 … q_m.
    let to_i = |k: usize| isize::try_from(k).expect("small");
    let rows: Vec<Expr> = (p..=n)
        .map(|k| {
            Expr::call(
                Symbol::LIST,
                (1..=m)
                    .map(|j| coefficient(to_i(k) - to_i(j)))
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    let rhs: Vec<Expr> = (p..=n).map(|k| -coefficient(to_i(k))).collect();
    let q: Vec<Expr> = if m == 0 {
        Vec::new()
    } else {
        let solution = ev.eval(&Expr::call(
            Symbol::LIN_SOLVE,
            [
                Expr::call(Symbol::LIST, rows),
                Expr::call(Symbol::LIST, rhs),
            ],
        ))?;
        if !is_list(&solution) || solution.args().len() != m {
            return Ok(None);
        }
        solution.args().to_vec()
    };
    let mut qs = vec![Expr::from(1)];
    qs.extend(q);
    let ps: Vec<Expr> = (0..p)
        .map(|k| {
            Expr::call(
                Symbol::PLUS,
                (0..=m.min(k))
                    .map(|j| qs[j].clone() * coefficient(to_i(k) - to_i(j)))
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    let result = expression(&ps, x) / expression(&qs, x);
    answer(ev, "pade", e, &result)
}

/// `ComplexFactor[p]`.
pub(in crate::native) fn complex_factor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let i = Expr::from(Complex::new(Rational::ZERO, Rational::ONE));
    let factored = ev.eval(&call(Symbol::FACTOR, vec![i.clone() * p.clone()]))?;
    if factored.head() == Some(&Expr::from(Symbol::FACTOR)) {
        return Ok(None);
    }
    answer(ev, "cfactor", e, &(factored / i))
}

/// `CompleteSquare[p, x]`.
pub(in crate::native) fn complete_square(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(x) = variable(e.args(), 1) else {
        return Ok(None);
    };
    let Some(cs) = coefficients(ev, &e.args()[0], &x)? else {
        return Ok(None);
    };
    let [c, b, a] = cs.as_slice() else {
        return Ok(None);
    };
    let shift = ev.eval(&(b.clone() / (2 * a.clone())))?;
    let rest = ev.eval(&(c.clone() - b.clone().pow(2) / (4 * a.clone())))?;
    // The square left as it is: Times and Plus of the parts, unevaluated.
    let square = Expr::call(
        Symbol::POWER,
        [Expr::call(Symbol::PLUS, [x.clone(), shift]), Expr::from(2)],
    );
    let result = Expr::call(
        Symbol::PLUS,
        [Expr::call(Symbol::TIMES, [a.clone(), square]), rest],
    );
    let result = ev.eval(&result)?;
    ev.step("canonical_form", e, &result);
    Ok(Some(result))
}

/// `DigitPolynomial[n, b, x]`.
pub(in crate::native) fn digit_polynomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, b, x] = e.args() else {
        return Ok(None);
    };
    let (Some(mut n), Some(b)) = (
        n.to_integer(),
        b.to_integer().filter(|b| *b >= Integer::from(2)),
    ) else {
        return Ok(None);
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    let mut digits = Vec::new();
    while !n.is_zero() {
        let (_, r) = n.div_rem_euclid(&b);
        let r = if (&r + &r) > b { &r - &b } else { r };
        n = (&n - &r).div_exact(&b);
        digits.push(Expr::from(r));
    }
    answer(ev, "genpoly", e, &expression(&digits, x))
}

/// `BerlekampMassey[s]`.
pub(in crate::native) fn berlekamp_massey(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [s] = e.args() else {
        return Ok(None);
    };
    if !is_list(s) || s.args().len() % 2 == 1 {
        return Ok(None);
    }
    let Some(s) = s
        .args()
        .iter()
        .map(Expr::to_rational)
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let mut c: Vec<Rational> = vec![Rational::ONE];
    let mut prev: Vec<Rational> = vec![Rational::ONE];
    let (mut l, mut m, mut last) = (0usize, 1usize, Rational::ONE);
    for n in 0..s.len() {
        let d = (1..=l).fold(s[n].clone(), |acc, i| &acc + &(&c[i] * &s[n - i]));
        if d.is_zero() {
            m += 1;
            continue;
        }
        let factor = &d / &last;
        let mut next = c.clone();
        if next.len() < prev.len() + m {
            next.resize(prev.len() + m, Rational::ZERO);
        }
        for (j, v) in prev.iter().enumerate() {
            next[j + m] = &next[j + m] - &(&factor * v);
        }
        if 2 * l <= n {
            prev = c;
            l = n + 1 - l;
            last = d;
            m = 1;
        } else {
            m += 1;
        }
        c = next;
    }
    c.resize(l + 1, Rational::ZERO);
    let result = Expr::call(
        Symbol::LIST,
        c.into_iter().map(Expr::from).collect::<Vec<_>>(),
    );
    answer(ev, "berlekamp_massey", e, &result)
}

#[cfg(test)]
mod tests;
