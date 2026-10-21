//! The heuristic gcd \[CGG89\]: the gcd of the values of two polynomials
//! at a large integer `ξ`, read back as a polynomial from its digits in
//! base `ξ`, and checked by dividing both. Over ℤ, Parisse's theorem
//! \[Par05\] makes the check enough in any number of variables, once `ξ`
//! is at least twice the smaller of the largest coefficients, plus 2; it
//! holds over ℤ\[i\] too. In several variables the last one that occurs is
//! evaluated, and the gcd of the values comes from the same algorithm in
//! one variable fewer.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the papers: a, b, g, ξ"
)]

use cas_num::algebra::{EuclideanDomain, GaussianIntegers, Integers, Ring};
use cas_num::{Complex, Integer, Natural, Rational};

use super::subres::{Exact, div_coeffs};
use super::{Cofactors, GcdDomain};
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};
use crate::sparse::{Sparse, SparseRing};

/// Points tried before giving up, as SymPy's `dup_zz_heu_gcd`.
const TRIES: usize = 6;

/// The most bits `ξ^(d + 1)` may take, `d` the degree evaluated: past it
/// the values are too large to be worth their gcd, and the heuristic gives
/// up for the modular algorithm (ours, CAS-137).
pub(crate) const MAX_BITS: u64 = 1 << 18;

/// The coefficient domains of the heuristic: ℤ, and ℤ\[i\], whose
/// polynomials are evaluated at integers too.
pub(crate) trait Digits: GcdDomain {
    /// The first point is `BOUND (min h + 1)`, `h` the largest
    /// [`Digits::height`] of each polynomial's coefficients: 2 over ℤ
    /// \[Par05\]. Over ℤ\[i\] the height bounds the modulus only within a
    /// factor `√2`, and 4 keeps Parisse's proof (ours).
    const BOUND: u64;

    /// A bound on the modulus of `c`: `|c|` over ℤ, `|re| + |im|` over
    /// ℤ\[i\].
    fn height(c: &Self::Elem) -> Natural;

    /// `(q, r)` with `c = q ξ + r`, each part of `r` in `(−ξ/2, ξ/2]`: the
    /// lowest digit of `c` in the symmetric representation.
    fn digit(c: &Self::Elem, xi: &Integer) -> (Self::Elem, Self::Elem);
}

impl Digits for Integers {
    const BOUND: u64 = 2;

    fn height(c: &Integer) -> Natural {
        c.unsigned_abs()
    }

    fn digit(c: &Integer, xi: &Integer) -> (Integer, Integer) {
        let (q, r) = c.div_rem_euclid(xi);
        if &r + &r > *xi {
            (q + Integer::ONE, r - xi)
        } else {
            (q, r)
        }
    }
}

impl Digits for GaussianIntegers {
    const BOUND: u64 = 4;

    fn height(c: &Complex) -> Natural {
        &c.re().numer().unsigned_abs() + &c.im().numer().unsigned_abs()
    }

    fn digit(c: &Complex, xi: &Integer) -> (Complex, Complex) {
        let (qr, rr) = Integers::digit(c.re().numer(), xi);
        let (qi, ri) = Integers::digit(c.im().numer(), xi);
        (
            Complex::new(Rational::from(qr), Rational::from(qi)),
            Complex::new(Rational::from(rr), Rational::from(ri)),
        )
    }
}

/// The first point: `BOUND (min(ha, hb) + 1)` for the heights of the two
/// polynomials \[Par05\].
pub(super) fn first_point<R: Digits>(ha: &Natural, hb: &Natural) -> Integer {
    Integer::from(&(ha.min(hb) + 1u64) * R::BOUND)
}

/// The next point after a failure: `73794 ξ ⌊ξ^(1/4)⌋ / 27011`, as SymPy's
/// `dup_zz_heu_gcd`, growing a little faster than `ξ` so that a few tries
/// reach well past a bad one.
pub(super) fn next_point(xi: &Integer) -> Integer {
    (xi * &Integer::from(73_794) * xi.isqrt().isqrt())
        .div_rem_floor(&Integer::from(27_011))
        .0
}

/// Whether the values at `ξ` of polynomials of degree `degree` are past
/// [`MAX_BITS`].
fn too_large(xi: &Integer, degree: u64) -> bool {
    xi.bit_length().saturating_mul(degree.saturating_add(1)) > MAX_BITS
}

/// The gcd `c` of the contents, and `a / c`, `b / c` over their primitive
/// parts: what turns cofactors of the primitive parts into cofactors.
struct Contents<E> {
    gcd: E,
    a: E,
    b: E,
}

impl<E> Contents<E> {
    fn of<R: EuclideanDomain<Elem = E>>(base: &R, ca: &E, cb: &E) -> Contents<E> {
        let gcd = base.gcd(ca, cb);
        let a = base.div_rem(ca, &gcd).0;
        let b = base.div_rem(cb, &gcd).0;
        Contents { gcd, a, b }
    }
}

/// The largest height of the coefficients.
fn height<'a, R: Digits>(coeffs: impl IntoIterator<Item = &'a R::Elem>) -> Natural
where
    R::Elem: 'a,
{
    coeffs
        .into_iter()
        .map(R::height)
        .max()
        .unwrap_or(Natural::ZERO)
}

/// The polynomial whose coefficients are the symmetric digits of `gamma`
/// in base `ξ`, lowest first.
pub(super) fn digits_dense<R: Digits>(
    ring: &DenseRing<R>,
    gamma: &R::Elem,
    xi: &Integer,
) -> Dense<R::Elem> {
    let base = ring.base();
    let mut coeffs = Vec::new();
    let mut v = gamma.clone();
    while !base.is_zero(&v) {
        let (q, r) = R::digit(&v, xi);
        coeffs.push(r);
        v = q;
    }
    ring.poly(coeffs)
}

/// The gcd of non-zero `a` and `b` with its cofactors, or `None` when the
/// heuristic gives up.
pub(crate) fn dense<R: Digits>(
    ring: &DenseRing<R>,
    a: &Dense<R::Elem>,
    b: &Dense<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Dense<R::Elem>>>, Cancelled> {
    let base = ring.base();
    let (ca, cb) = (ring.content(a), ring.content(b));
    let (pa, pb) = (div_coeffs(ring, a, &ca), div_coeffs(ring, b, &cb));
    let contents = Contents::of(base, &ca, &cb);
    let done = |g: &Dense<R::Elem>, qa: &Dense<R::Elem>, qb: &Dense<R::Elem>| {
        (
            ring.scale(g, &contents.gcd),
            ring.scale(qa, &contents.a),
            ring.scale(qb, &contents.b),
        )
    };
    let (Some(m), Some(n)) = (pa.degree(), pb.degree()) else {
        unreachable!("non-zero polynomials");
    };
    if m == 0 || n == 0 {
        // A primitive constant is a unit.
        return Ok(Some(done(&ring.one(), &pa, &pb)));
    }
    let degree = m.max(n) as u64;
    let mut xi = first_point::<R>(&height::<R>(pa.coeffs()), &height::<R>(pb.coeffs()));
    for _ in 0..TRIES {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        if too_large(&xi, degree) {
            return Ok(None);
        }
        let x = base.elem_integer(&xi);
        let gamma = base.gcd(&ring.eval(&pa, &x), &ring.eval(&pb, &x));
        let g = digits_dense(ring, &gamma, &xi);
        let g = div_coeffs(ring, &g, &ring.content(&g));
        if let Some(qa) = ring.div_exact(&pa, &g)
            && let Some(qb) = ring.div_exact(&pb, &g)
        {
            return Ok(Some(done(&g, &qa, &qb)));
        }
        xi = next_point(&xi);
    }
    Ok(None)
}

/// `p` over its content.
fn primitive<R: Digits>(ring: &SparseRing<R>, p: &Sparse<R::Elem>) -> Sparse<R::Elem> {
    let c = ring.content(p);
    if ring.base().is_one(&c) {
        return p.clone();
    }
    ring.map_into(p, ring, |x| ring.base().exact(x, &c))
}

/// The gcd of non-zero `a` and `b` with its cofactors, or `None` when the
/// heuristic gives up.
pub(crate) fn sparse<R: Digits>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Sparse<R::Elem>>>, Cancelled> {
    let base = ring.base();
    let (ca, cb) = (ring.content(a), ring.content(b));
    let divide = |p: &Sparse<R::Elem>, c: &R::Elem| ring.map_into(p, ring, |x| base.exact(x, c));
    let (pa, pb) = (divide(a, &ca), divide(b, &cb));
    let contents = Contents::of(base, &ca, &cb);
    let Some((g, qa, qb)) = primitive_gcd(ring, &pa, &pb, cancel)? else {
        return Ok(None);
    };
    Ok(Some((
        ring.scale(&g, &contents.gcd),
        ring.scale(&qa, &contents.a),
        ring.scale(&qb, &contents.b),
    )))
}

/// The gcd of non-zero primitive `a` and `b`, its cofactors, or `None`:
/// the last variable that occurs evaluated at `ξ`, the gcd of the values
/// from [`sparse`], read back from its digits and checked by division.
fn primitive_gcd<R: Digits>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Sparse<R::Elem>>>, Cancelled> {
    let base = ring.base();
    let (da, db) = (a.degrees(), b.degrees());
    let Some(var) = (0..ring.nvars()).rev().find(|&v| da[v] > 0 || db[v] > 0) else {
        // Primitive constants are units.
        return Ok(Some((ring.one(), a.clone(), b.clone())));
    };
    let degree = da[var].max(db[var]);
    let mut xi = first_point::<R>(&height::<R>(a.coeffs()), &height::<R>(b.coeffs()));
    for _ in 0..TRIES {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        if too_large(&xi, degree) {
            return Ok(None);
        }
        let x = base.elem_integer(&xi);
        let (va, vb) = (ring.eval_var(a, var, &x), ring.eval_var(b, var, &x));
        if !va.is_zero() && !vb.is_zero() {
            let Some((gamma, _, _)) = sparse(ring, &va, &vb, cancel)? else {
                return Ok(None);
            };
            let g = primitive(ring, &digits_sparse(ring, &gamma, var, &xi));
            if let Some(qa) = ring.div_exact_with(a, &g, cancel)?
                && let Some(qb) = ring.div_exact_with(b, &g, cancel)?
            {
                return Ok(Some((g, qa, qb)));
            }
        }
        xi = next_point(&xi);
    }
    Ok(None)
}

/// The polynomial in `x_var` whose coefficients are the symmetric digits
/// in base `ξ` of those of `gamma`, which has no `x_var`.
fn digits_sparse<R: Digits>(
    ring: &SparseRing<R>,
    gamma: &Sparse<R::Elem>,
    var: usize,
    xi: &Integer,
) -> Sparse<R::Elem> {
    let base = ring.base();
    let mut terms = Vec::new();
    for (e, c) in gamma.terms() {
        let mut v = c.clone();
        let mut k = 0;
        while !base.is_zero(&v) {
            let (q, r) = R::digit(&v, xi);
            if !base.is_zero(&r) {
                let mut f = e.clone();
                f[var] = k;
                terms.push((f, r));
            }
            v = q;
            k += 1;
        }
    }
    ring.from_terms(terms)
}
