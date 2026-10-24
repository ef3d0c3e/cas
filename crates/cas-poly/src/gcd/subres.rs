//! Subresultant remainder sequences \[BT71, Bro78\], in Cohen's form
//! \[Coh93, §3.3\]: the gcd and the resultant over any domain with exact
//! division,
//! polynomials in further variables included, so that sparse polynomials
//! take their gcd in one variable over the gcds of their coefficients in
//! the others. The reference the faster algorithms are checked against.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the remainder sequences: a, b, g, h, r"
)]

use cas_num::algebra::{EuclideanDomain, Ring};

use super::{Cofactors, normalized_sparse};
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Coefficients, Dense, DenseRing};
use crate::sparse::{Sparse, SparseRing};

/// What the remainder sequences need of the coefficients: division that is
/// known to be exact.
pub(crate) trait Exact: Coefficients {
    /// `a / b`, for a `b` that divides `a`.
    fn exact(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
}

impl<R: Coefficients + EuclideanDomain> Exact for R {
    fn exact(&self, a: &R::Elem, b: &R::Elem) -> R::Elem {
        let (q, r) = self.div_rem(a, b);
        debug_assert!(self.is_zero(&r), "an exact division");
        q
    }
}

impl<R: EuclideanDomain> Exact for SparseRing<R> {
    fn exact(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> Sparse<R::Elem> {
        self.div_exact(a, b).expect("an exact division")
    }
}

/// `p / c`, coefficient by coefficient, for a `c` that divides each.
pub(crate) fn div_coeffs<C: Exact>(
    ring: &DenseRing<C>,
    p: &Dense<C::Elem>,
    c: &C::Elem,
) -> Dense<C::Elem> {
    let base = ring.base();
    if base.is_one(c) {
        return p.clone();
    }
    ring.poly(p.coeffs().iter().map(|x| base.exact(x, c)).collect())
}

/// The degree of a non-zero polynomial, as a power.
fn degree<E>(p: &Dense<E>) -> u64 {
    p.degree().expect("a non-zero polynomial") as u64
}

/// The last non-zero polynomial of the subresultant sequence of non-zero
/// `a` and `b` with `deg a ≥ deg b` \[Coh93, §3.3.1\], without the contents:
/// an associate of their gcd over the field of fractions of the
/// coefficients, or `None` when that gcd is 1. Each round pseudo-divides
/// and divides the remainder by `g h^δ`, exactly \[BT71\].
fn last_remainder<C: Exact>(
    ring: &DenseRing<C>,
    mut a: Dense<C::Elem>,
    mut b: Dense<C::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Dense<C::Elem>>, Cancelled> {
    let base = ring.base();
    let (mut g, mut h) = (base.one(), base.one());
    loop {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let delta = degree(&a) - degree(&b);
        let r = ring.pseudo_div_rem(&a, &b).1;
        match r.degree() {
            None => return Ok(Some(b)),
            Some(0) => return Ok(None),
            Some(_) => {}
        }
        let divisor = base.mul(&g, &base.pow(&h, delta));
        a = b;
        b = div_coeffs(ring, &r, &divisor);
        g = a.lead().expect("non-zero").clone();
        // h <- h^(1 - δ) g^δ.
        if delta > 0 {
            h = base.exact(&base.pow(&g, delta), &base.pow(&h, delta - 1));
        }
    }
}

/// A gcd of coefficients, which may be cancelled.
type GcdOf<'a, E> = dyn Fn(&E, &E) -> Result<E, Cancelled> + 'a;

/// The gcd of the coefficients of `p`, by `gcd`.
fn content<C: Ring>(
    ring: &C,
    p: &Dense<C::Elem>,
    gcd: &GcdOf<'_, C::Elem>,
) -> Result<C::Elem, Cancelled> {
    let mut c = ring.zero();
    for x in p.coeffs() {
        c = gcd(&c, x)?;
        if ring.is_one(&c) {
            break;
        }
    }
    Ok(c)
}

/// The gcd of non-zero `a` and `b`, up to a unit, over coefficients whose
/// gcd, in normal form, is `gcd` \[Coh93, §3.3.1\]: the gcd of the
/// contents times the primitive part of the last remainder of the
/// primitive parts.
fn gcd_over<C: Exact>(
    ring: &DenseRing<C>,
    a: &Dense<C::Elem>,
    b: &Dense<C::Elem>,
    gcd: &GcdOf<'_, C::Elem>,
    cancel: &dyn Cancel,
) -> Result<Dense<C::Elem>, Cancelled> {
    let base = ring.base();
    let (ca, cb) = (content(base, a, gcd)?, content(base, b, gcd)?);
    let d = gcd(&ca, &cb)?;
    let (pa, pb) = (div_coeffs(ring, a, &ca), div_coeffs(ring, b, &cb));
    let (big, small) = if pa.degree() >= pb.degree() {
        (pa, pb)
    } else {
        (pb, pa)
    };
    // A primitive constant is a unit, so the gcd is the contents': the
    // remainder sequence would pseudo-divide `big` by it in as many rounds
    // as its degree, each through all its coefficients (CAS-166).
    if small.degree() == Some(0) {
        return Ok(ring.constant(d));
    }
    Ok(match last_remainder(ring, big, small, cancel)? {
        None => ring.constant(d),
        Some(last) => {
            let c = content(base, &last, gcd)?;
            ring.scale(&div_coeffs(ring, &last, &c), &d)
        }
    })
}

/// The gcd of non-zero dense `a` and `b`, with its cofactors, by
/// subresultants.
pub(crate) fn dense_gcd<R: Coefficients + EuclideanDomain>(
    ring: &DenseRing<R>,
    a: &Dense<R::Elem>,
    b: &Dense<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Dense<R::Elem>>, Cancelled> {
    let base = ring.base();
    let g = gcd_over(ring, a, b, &|x, y| Ok(base.gcd(x, y)), cancel)?;
    let a_cofactor = ring.div_exact(a, &g).expect("the gcd divides");
    let b_cofactor = ring.div_exact(b, &g).expect("the gcd divides");
    Ok((g, a_cofactor, b_cofactor))
}

/// The gcd of `a` and `b`, normalised: by subresultants in the first
/// variable that occurs, over the gcds of the coefficients in the others,
/// recursively.
pub(crate) fn sparse_gcd_alone<R: EuclideanDomain + Clone>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    if a.is_zero() {
        return Ok(normalized_sparse(ring, b));
    }
    if b.is_zero() {
        return Ok(normalized_sparse(ring, a));
    }
    let (da, db) = (a.degrees(), b.degrees());
    // A constant's gcd is that of the coefficients of both, at once
    // whatever the other's degree: its dense form in a variable has that
    // many coefficients (CAS-166).
    let constant = |d: &[u64]| d.iter().all(|&e| e == 0);
    if constant(&da) || constant(&db) {
        let base = ring.base();
        let mut g = base.zero();
        for c in a.coeffs().iter().chain(b.coeffs()) {
            g = base.gcd(&g, c);
            if base.is_one(&g) {
                break;
            }
        }
        return Ok(ring.constant(g));
    }
    let var = (0..ring.nvars())
        .find(|&v| da[v] > 0 || db[v] > 0)
        .expect("a variable of a polynomial that is not constant");
    let dense = DenseRing::new(ring.clone());
    let gcd = |x: &Sparse<R::Elem>, y: &Sparse<R::Elem>| sparse_gcd_alone(ring, x, y, cancel);
    let g = gcd_over(
        &dense,
        &ring.to_dense_in(a, var),
        &ring.to_dense_in(b, var),
        &gcd,
        cancel,
    )?;
    Ok(normalized_sparse(ring, &ring.from_dense_in(&g, var)))
}

/// The gcd of non-zero sparse `a` and `b`, with its cofactors, by
/// subresultants.
pub(crate) fn sparse_gcd<R: EuclideanDomain + Clone>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Sparse<R::Elem>>, Cancelled> {
    let g = sparse_gcd_alone(ring, a, b, cancel)?;
    let a_cofactor = ring
        .div_exact_with(a, &g, cancel)?
        .expect("the gcd divides");
    let b_cofactor = ring
        .div_exact_with(b, &g, cancel)?
        .expect("the gcd divides");
    Ok((g, a_cofactor, b_cofactor))
}

/// The resultant of non-zero `a` and `b` \[Coh93, §3.3.2\], without the
/// contents: the remainder sequence of [`last_remainder`], with the sign of
/// each exchange of odd degrees, and the last `h` raised to the last
/// degree.
pub(crate) fn resultant<C: Exact>(
    ring: &DenseRing<C>,
    a: &Dense<C::Elem>,
    b: &Dense<C::Elem>,
    cancel: &dyn Cancel,
) -> Result<C::Elem, Cancelled> {
    let base = ring.base();
    let (mut a, mut b) = (a.clone(), b.clone());
    // res(a, b) = (-1)^(deg a deg b) res(b, a).
    let mut negative = false;
    if degree(&a) < degree(&b) {
        negative = degree(&a) % 2 == 1 && degree(&b) % 2 == 1;
        std::mem::swap(&mut a, &mut b);
    }
    let (mut g, mut h) = (base.one(), base.one());
    while degree(&b) > 0 {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let (m, n) = (degree(&a), degree(&b));
        let delta = m - n;
        if m % 2 == 1 && n % 2 == 1 {
            negative = !negative;
        }
        let r = ring.pseudo_div_rem(&a, &b).1;
        if r.is_zero() {
            return Ok(base.zero());
        }
        let divisor = base.mul(&g, &base.pow(&h, delta));
        a = b;
        b = div_coeffs(ring, &r, &divisor);
        g = a.lead().expect("non-zero").clone();
        if delta > 0 {
            h = base.exact(&base.pow(&g, delta), &base.pow(&h, delta - 1));
        }
    }
    // b is a non-zero constant: res = lc(b)^m / h^(m - 1).
    let m = degree(&a);
    let lead = b.lead().expect("non-zero");
    let res = if m == 0 {
        base.one()
    } else {
        base.exact(&base.pow(lead, m), &base.pow(&h, m - 1))
    };
    Ok(if negative { base.neg(&res) } else { res })
}

/// The resultant of non-zero `a` and `b` in `x_var`, by subresultants over
/// the polynomials in the other variables.
pub(crate) fn sparse_resultant<R: EuclideanDomain + Clone>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    var: usize,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    let dense = DenseRing::new(ring.clone());
    resultant(
        &dense,
        &ring.to_dense_in(a, var),
        &ring.to_dense_in(b, var),
        cancel,
    )
}
