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
    let (r1, t1) = euclid_until(Integer::from(m.clone()), Integer::from(a % m), num_bound);
    // Floored remainders of non-negative numbers are non-negative.
    let r1 = r1.unsigned_abs();
    let (negative, d) = t1.into_parts();
    if d > *den_bound || !r1.gcd(&d).is_one() {
        return None;
    }
    Some(Rational::new(
        Integer::from_parts(negative, r1),
        Integer::from(d),
    ))
}

/// Xcas's `fracmod` (giac's `in_fracmod`, gen.cc): Wang's loop on `(m, a)`,
/// `a` as given (not reduced modulo `m`) and the quotients floored, stopped
/// at the first remainder `|r| <= ⌊√⌊m/2⌋⌋`; the fraction `r/t`, or `None`
/// when its cofactor `t` is not prime to `m`. No denominator bound, so the
/// fraction need not be [`rational_reconstruction`]'s: `fracmod(13, 121)`
/// is `−4/9`, `fracmod(2, 4)` is `None`, and `m = 0` gives `None` for
/// every `a ≠ 0`.
#[must_use]
pub fn fracmod(a: &Integer, m: &Natural) -> Option<Rational> {
    let bound = (m.clone() >> 1).isqrt();
    let (r, t) = euclid_until(Integer::from(m.clone()), a.clone(), &bound);
    if !t.unsigned_abs().gcd(m).is_one() {
        return None;
    }
    Some(Rational::new(r, t))
}

/// The half-extended Euclidean algorithm on `(r0, r1)`, quotients floored
/// (GMP's `mpz_fdiv_qr`), stopped at the first remainder `|r| <= bound`:
/// that remainder and its cofactor `t`, `r ≡ t r1 (mod r0)`.
fn euclid_until(mut r0: Integer, mut r1: Integer, bound: &Natural) -> (Integer, Integer) {
    let (mut t0, mut t1) = (Integer::ZERO, Integer::ONE);
    while r1.unsigned_abs() > *bound {
        let (q, r) = r0.div_rem_floor(&r1);
        let t = t0 - &(q * &t1);
        (r0, r1) = (r1, r);
        (t0, t1) = (t1, t);
    }
    (r1, t1)
}
