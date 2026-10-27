//! Gosper's algorithm (CAS-DOC-13 §3, CAS-225): the antidifference of a
//! hypergeometric term, when it has one that is hypergeometric \[Gos78\].
//!
//! - The ratio `r(k) = t(k + 1)/t(k)` is taken factor by factor
//!   ([`ratio`]): `c^u` gives `c^(u(k+1) − u(k))` for a `u` linear in `k`;
//!   `u!` for a `u` linear in `k` with an integer slope `α` the product
//!   `(u + 1) ⋯ (u + α)`, or its inverse for a negative slope;
//!   `Binomial[u, v]` as `u!/(v! (u − v)!)`; a rational function of `k`
//!   its quotient at `k + 1` and `k`; integer powers of these. Anything else
//!   is no hypergeometric term here.
//! - `r = (a(k)/b(k)) (c(k + 1)/c(k))` with `gcd(a(k), b(k + h)) = 1` for
//!   every integer `h ≥ 0`: the shifts `h` are the nonnegative integer
//!   roots of `res_k(a(k), b(k + h))`, from its linear factors, and for
//!   each `g = gcd(a(k), b(k + h))` leaves `a` and `g(k − h)` leaves `b`
//!   for `c (g(k − 1) ⋯ g(k − h))`.
//! - A polynomial `x(k)` with `a(k) x(k + 1) − b(k − 1) x(k) = c(k)`, of
//!   Gosper's degree bound, by the linear system of its coefficients
//!   (CAS-212); its free coefficients 0. None: `t` has no hypergeometric
//!   antidifference, as `1/k`.
//! - `T(k) = b(k − 1) x(k)/c(k) t(k)`, and `T(k + 1) − T(k) = t(k)`.

use cas_num::algebra::{Integers, Ring};
use cas_num::{Integer, Rational};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::super::polynomial::{convert, factored};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::{Fraction, Polynomials};
use crate::symbol::Symbol;

/// The Gosper–Petkovšek form `(a, b, c)` of a ratio.
type Form = (Sparse<Integer>, Sparse<Integer>, Sparse<Integer>);

/// The highest degree of `x(k)` tried.
const MAX_DEGREE: i64 = 64;

/// The largest shift tried.
const MAX_SHIFT: i64 = 256;

/// `e` with `k` replaced by `k + 1`.
fn next(e: &Expr, k: &Expr) -> Expr {
    e.substitute(&[(k.clone(), k.clone() + Expr::from(1))])
}

/// `u(k + 1) − u(k)`, expanded: evaluation alone leaves `2 + k − (1 + k)`.
fn step(ev: &mut Eval<'_>, u: &Expr, k: &Expr) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(Symbol::EXPAND, [next(u, k) - u.clone()]))
}

/// Whether `e` is a rational function of `k`: `k` only through sums,
/// products and integer powers.
pub(super) fn rational_in(e: &Expr, k: &Expr) -> bool {
    if e == k || !e.contains(k) {
        return true;
    }
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), args) => args.iter().all(|a| rational_in(a, k)),
        (Some(Symbol::POWER), [b, n]) => n.as_i64().is_some() && rational_in(b, k),
        _ => false,
    }
}

/// The ratio `t(k + 1)/t(k)` of a hypergeometric term, built; `None` when
/// `t` is not one here.
pub(super) fn ratio(ev: &mut Eval<'_>, t: &Expr, k: &Expr) -> Result<Option<Expr>, EvalError> {
    if !t.contains(k) {
        return Ok(Some(Expr::from(1)));
    }
    if rational_in(t, k) {
        return Ok(Some(next(t, k) / t.clone()));
    }
    Ok(match (t.head().and_then(Expr::as_symbol), t.args()) {
        (Some(Symbol::TIMES), args) => {
            let mut ratios = Vec::with_capacity(args.len());
            for a in args {
                let Some(r) = ratio(ev, a, k)? else {
                    return Ok(None);
                };
                ratios.push(r);
            }
            Some(Expr::call(Symbol::TIMES, ratios))
        }
        (Some(Symbol::POWER), [b, u]) if !b.contains(k) => Some(b.clone().pow(step(ev, u, k)?)),
        (Some(Symbol::POWER), [b, n]) if n.as_i64().is_some() => {
            ratio(ev, b, k)?.map(|r| r.pow(n.clone()))
        }
        (Some(Symbol::EXP), [u]) => Some(Expr::call(Symbol::EXP, [step(ev, u, k)?])),
        (Some(Symbol::FACTORIAL), [u]) => factorial_ratio(ev, u, k)?,
        (Some(Symbol::BINOMIAL), [u, v]) => {
            let w = ev.eval(&(u.clone() - v.clone()))?;
            match (
                factorial_ratio(ev, u, k)?,
                factorial_ratio(ev, v, k)?,
                factorial_ratio(ev, &w, k)?,
            ) {
                (Some(ru), Some(rv), Some(rw)) => Some(ru / (rv * rw)),
                _ => None,
            }
        }
        _ => None,
    })
}

/// The ratio of `u!` at `k + 1` and `k`, for `u` linear in `k` with an
/// integer slope `α`: `(u + 1) ⋯ (u + α)`, or `1/(u (u − 1) ⋯ (u + α + 1))`
/// for a negative one.
fn factorial_ratio(ev: &mut Eval<'_>, u: &Expr, k: &Expr) -> Result<Option<Expr>, EvalError> {
    if !u.contains(k) {
        return Ok(Some(Expr::from(1)));
    }
    if !rational_in(u, k) {
        return Ok(None);
    }
    let slope = step(ev, u, k)?;
    let Some(alpha) = slope.as_i64().filter(|a| a.unsigned_abs() <= 64) else {
        return Ok(None);
    };
    let product = |range: std::ops::RangeInclusive<i64>| {
        Expr::call(
            Symbol::TIMES,
            range.map(|i| u.clone() + Expr::from(i)).collect::<Vec<_>>(),
        )
    };
    Ok(Some(match alpha {
        0 => Expr::from(1),
        a if a > 0 => product(1..=a),
        a => product(a + 1..=0).pow(-1),
    }))
}

/// An antidifference `T` of the hypergeometric term `t` in `k`:
/// `T(k + 1) − T(k) = t(k)`, built; `None` when Gosper's algorithm finds
/// none.
pub(super) fn antidifference(
    ev: &mut Eval<'_>,
    e: &Expr,
    t: &Expr,
    k: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some(r) = ratio(ev, t, k)? else {
        return Ok(None);
    };
    let r = ev.eval(&Expr::call(Symbol::TOGETHER, [r]))?;
    let h = Expr::from(Symbol::new("Sum$h"));
    // k with them, a generator even where the ratio is free of it (2^k).
    let Some(p) = convert(ev, e, &[r, h.clone(), k.clone()], false)? else {
        return Ok(None);
    };
    let (Some(kv), Some(hv)) = (p.variable(k), p.variable(&h)) else {
        return Ok(None);
    };
    let inside = p
        .generators()
        .iter()
        .enumerate()
        .any(|(j, g)| j != kv && g.contains(k));
    let Fraction::Integers { num, den } = &p.fractions()[0] else {
        return Ok(None);
    };
    if inside {
        return Ok(None);
    }
    let ring = SparseRing::new(Integers, p.generators().len());
    let cancel: &dyn Cancel = ev.context().interrupt_flag();
    let Some((a, b, c)) = decomposed(ev, e, &ring, num.clone(), den.clone(), kv, hv, cancel)?
    else {
        return Ok(None);
    };
    let one = ring.one();
    let kk = ring.var(kv);
    let b1 = ring.substitute(&b, kv, &ring.sub(&kk, &one));
    let Some(d) = degree_bound(&ring, &a, &b1, &c, kv) else {
        return Ok(None);
    };
    let Some(x) = solved(ev, e, &p, &[&a, &b1, &c], k, d)? else {
        return Ok(None);
    };
    // Together, c against the factors of t: k 2^k has T = (k − 2)/k · k 2^k,
    // which at k = 0 would be 0/0.
    let fraction = p.polynomial_expr(&b1) * x / p.polynomial_expr(&c);
    Ok(Some(ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [fraction * t.clone()],
    ))?))
}

/// The Gosper–Petkovšek form `(a, b, c)` of `num/den` in `kv`, the shift
/// variable `hv` for the resultant; `None` when factoring stops.
#[allow(
    clippy::too_many_arguments,
    reason = "the ring, the fraction and its two variables"
)]
fn decomposed(
    ev: &mut Eval<'_>,
    e: &Expr,
    ring: &SparseRing<Integers>,
    mut a: Sparse<Integer>,
    mut b: Sparse<Integer>,
    kv: usize,
    hv: usize,
    cancel: &dyn Cancel,
) -> Result<Option<Form>, EvalError> {
    let kk = ring.var(kv);
    let shifted = ring.substitute(&b, kv, &ring.add(&kk, &ring.var(hv)));
    let resultant = ring
        .resultant_with(&a, &shifted, kv, cancel)
        .map_err(|_| EvalError::Interrupted)?;
    let mut c = ring.one();
    if resultant.is_zero() {
        return Ok(None);
    }
    let Some(powers) = factored(ev, e, ring, &resultant, &ring.one())? else {
        return Ok(None);
    };
    let mut shifts: Vec<i64> = powers
        .factors
        .iter()
        .filter_map(|(f, _)| integer_root(ring, f, hv))
        .filter(|&h| (0..=MAX_SHIFT).contains(&h))
        .collect();
    shifts.sort_unstable();
    shifts.dedup();
    for h in shifts {
        let constant = |n: i64| ring.constant(Integer::from(n));
        let forward = ring.substitute(&b, kv, &ring.add(&kk, &constant(h)));
        let g = ring
            .gcd_with(&a, &forward, cancel)
            .map_err(|_| EvalError::Interrupted)?;
        if ring.degree(&g, kv).unwrap_or(0) == 0 {
            continue;
        }
        a = ring.div_exact(&a, &g).expect("a divisor");
        let back = ring.substitute(&g, kv, &ring.sub(&kk, &constant(h)));
        b = ring.div_exact(&b, &back).expect("a divisor");
        for i in 1..=h {
            c = ring.mul(&c, &ring.substitute(&g, kv, &ring.sub(&kk, &constant(i))));
        }
    }
    Ok(Some((a, b, c)))
}

/// The root `−β/α` of `f` when it is `α v + β` in the variable `v` alone.
pub(super) fn linear_root(
    ring: &SparseRing<Integers>,
    f: &Sparse<Integer>,
    v: usize,
) -> Option<Rational> {
    if f.degrees()
        .iter()
        .enumerate()
        .any(|(j, &d)| j != v && d > 0)
        || ring.degree(f, v) != Some(1)
    {
        return None;
    }
    let alpha = ring.coeff(f, v, 1).coeffs().first()?.clone();
    let beta = ring
        .coeff(f, v, 0)
        .coeffs()
        .first()
        .cloned()
        .unwrap_or(Integer::ZERO);
    Some(Rational::new(-beta, alpha))
}

/// The integer root of `f` when it is `α h + β` in the variable `hv` alone.
fn integer_root(ring: &SparseRing<Integers>, f: &Sparse<Integer>, hv: usize) -> Option<i64> {
    let root = linear_root(ring, f, hv)?;
    root.is_integer().then(|| root.numer().to_i64()).flatten()
}

/// Gosper's bound on the degree of `x(k)`: `deg c − max(deg a, deg b)`
/// when the degrees or the leading coefficients of `a` and `b(k − 1)`
/// differ; otherwise the larger of `deg c − deg a + 1` and
/// `(B − A)/lead`, `A` and `B` their coefficients of degree one less, when
/// it is an integer. `None` for a negative bound or one past
/// [`MAX_DEGREE`].
fn degree_bound(
    ring: &SparseRing<Integers>,
    a: &Sparse<Integer>,
    b1: &Sparse<Integer>,
    c: &Sparse<Integer>,
    kv: usize,
) -> Option<i64> {
    let degree = |p: &Sparse<Integer>| i64::try_from(ring.degree(p, kv).unwrap_or(0)).ok();
    let (da, db, dc) = (degree(a)?, degree(b1)?, degree(c)?);
    let lead = |p: &Sparse<Integer>, d: i64| ring.coeff(p, kv, u64::try_from(d).unwrap_or(0));
    let d = if da != db || lead(a, da) != lead(b1, db) {
        dc - da.max(db)
    } else {
        let base = dc - da + 1;
        let number = |p: Sparse<Integer>| -> Option<Rational> {
            if p.is_zero() {
                return Some(Rational::ZERO);
            }
            if p.degrees().iter().any(|&d| d > 0) {
                return None;
            }
            p.coeffs().first().map(|c| Rational::from(c.clone()))
        };
        let candidate = if da >= 1 {
            let l = number(lead(a, da));
            let (x, y) = (number(lead(a, da - 1)), number(lead(b1, db - 1)));
            match (l, x, y) {
                (Some(l), Some(x), Some(y)) if !l.is_zero() => {
                    let q = &(&y - &x) / &l;
                    q.is_integer().then(|| q.numer().to_i64()).flatten()
                }
                _ => None,
            }
        } else {
            None
        };
        candidate.map_or(base, |c| c.max(base))
    };
    (0..=MAX_DEGREE).contains(&d).then_some(d)
}

/// The polynomial `x(k)` of degree at most `d` with
/// `a(k) x(k + 1) − b1(k) x(k) = c(k)` for `[a, b1, c]` in `p`'s
/// generators, built; `None` when there is none.
fn solved(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Polynomials,
    [a, b1, c]: &[&Sparse<Integer>; 3],
    k: &Expr,
    d: i64,
) -> Result<Option<Expr>, EvalError> {
    let unknowns: Vec<Expr> = (0..=d)
        .map(|i| Expr::from(Symbol::new(&format!("Sum$x{i}"))))
        .collect();
    // No at^0, which at = 0 would make 0^0.
    let x_at = |at: &Expr| -> Expr {
        Expr::call(
            Symbol::PLUS,
            unknowns
                .iter()
                .enumerate()
                .map(|(i, u)| match i {
                    0 => u.clone(),
                    _ => u.clone() * at.clone().pow(i64::try_from(i).expect("a small degree")),
                })
                .collect::<Vec<_>>(),
        )
    };
    let equation = p.polynomial_expr(a) * x_at(&(k.clone() + Expr::from(1)))
        - p.polynomial_expr(b1) * x_at(k)
        - p.polynomial_expr(c);
    let equation = ev.eval(&equation)?;
    // Its coefficients in k, each 0; itself when k cancels.
    let Some(q) = convert(ev, e, std::slice::from_ref(&equation), false)? else {
        return Ok(None);
    };
    let Fraction::Integers { num, .. } = &q.fractions()[0] else {
        return Ok(None);
    };
    let ring = SparseRing::new(Integers, q.generators().len());
    let zero = |x: Expr| Expr::call(Symbol::EQUATION, [x, Expr::from(0)]);
    let equations: Vec<Expr> = match q.variable(k) {
        Some(kv) => (0..=ring.degree(num, kv).unwrap_or(0))
            .map(|j| zero(q.polynomial_expr(&ring.coeff(num, kv, j))))
            .collect(),
        None => vec![zero(q.polynomial_expr(num))],
    };
    let Some(solutions) = super::super::solve::linear_system(ev, e, &equations, &unknowns)? else {
        return Ok(None);
    };
    let [tuple] = solutions.args() else {
        return Ok(None);
    };
    // Free coefficients 0.
    let zeros: Vec<(Expr, Expr)> = unknowns
        .iter()
        .map(|u| (u.clone(), Expr::from(0)))
        .collect();
    let values: Vec<(Expr, Expr)> = unknowns
        .iter()
        .zip(tuple.args())
        .map(|(u, v)| (u.clone(), v.substitute(&zeros)))
        .collect();
    Ok(Some(x_at(k).substitute(&values)))
}
