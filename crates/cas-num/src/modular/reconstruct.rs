//! Rational reconstruction.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, m, q, r, t)"
)]

use crate::{Integer, Natural, Rational};

#[cfg(test)]
mod tests;

/// The fraction `n/d` with `n ≡ a d (mod m)`, `|n| <= N`, `0 < d <= N` and
/// `gcd(d, m) = 1`, for the largest `N` with `2 N^2 < m`, or `None` if there
/// is none. Within these bounds the fraction is unique.
///
/// This is the usual last step of a multi-modular algorithm whose result is
/// rational: `a` is the result modulo `m`, and `m` has grown past twice the
/// product of the numerator and denominator bounds.
///
/// # Panics
///
/// If `m` is zero.
#[must_use]
pub fn rational_reconstruction(a: &Natural, m: &Natural) -> Option<Rational> {
    assert!(!m.is_zero(), "zero modulus");
    let bound = ((m - 1u64) >> 1).isqrt();
    reconstruct(a, m, &bound, &bound)
}

/// The fraction `n/d` with `n ≡ a d (mod m)`, `|n| <= num_bound`,
/// `0 < d <= den_bound` and `gcd(d, m) = 1`, or `None` if there is none.
/// The bounds must satisfy `2 num_bound den_bound < m`, which makes the
/// fraction unique.
///
/// # Panics
///
/// If `m` is zero, or `2 num_bound den_bound >= m`.
#[must_use]
pub fn rational_reconstruction_with_bounds(
    a: &Natural,
    m: &Natural,
    num_bound: &Natural,
    den_bound: &Natural,
) -> Option<Rational> {
    assert!(!m.is_zero(), "zero modulus");
    assert!(
        (num_bound * den_bound) << 1 < *m,
        "bounds too large: need 2 num_bound den_bound < m"
    );
    reconstruct(a, m, num_bound, den_bound)
}

/// Wang's algorithm: the extended Euclidean algorithm on `(m, a mod m)`,
/// stopped at the first remainder `r <= N`. Every remainder satisfies
/// `r ≡ t a (mod m)` for its cofactor `t`, and when a fraction within the
/// bounds exists it is `r / t` (von zur Gathen and Gerhard, Theorem 5.26).
/// A common factor of `r` and `t` would divide `m` and `d`, so it rules a
/// solution out.
fn reconstruct(
    a: &Natural,
    m: &Natural,
    num_bound: &Natural,
    den_bound: &Natural,
) -> Option<Rational> {
    let (mut r0, mut r1) = (m.clone(), a % m);
    let (mut t0, mut t1) = (Integer::ZERO, Integer::ONE);
    while r1 > *num_bound {
        let (q, r) = r0.div_rem(&r1);
        let t = t0 - &(Integer::from(q) * &t1);
        (r0, r1) = (r1, r);
        (t0, t1) = (t1, t);
    }
    let (negative, d) = t1.into_parts();
    if d > *den_bound || !r1.gcd(&d).is_one() {
        return None;
    }
    Some(Rational::new(
        Integer::from_parts(negative, r1),
        Integer::from(d),
    ))
}
