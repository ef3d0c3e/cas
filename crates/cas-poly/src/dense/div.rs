//! Division: with remainder over a field, pseudo-division and exact
//! division over a Euclidean domain such as ℤ \[Knu97, §4.6.1\].

#![allow(
    clippy::many_single_char_names,
    reason = "the names of a = q b + r, with d and c for coefficients"
)]

use cas_num::algebra::{EuclideanDomain, Field};

use super::{Coefficients, Dense, DenseRing};
use crate::cancel::{Cancel, Cancelled};

/// A quotient and a remainder.
pub(crate) type QuoRem<E> = (Dense<E>, Dense<E>);

impl<F: Coefficients + Field> DenseRing<F> {
    /// `(q, r)` with `a = q b + r` and `deg r < deg b`: long division by the
    /// inverse of `b`'s leading coefficient (Knuth's Algorithm D of §4.6.1),
    /// polling `cancel` at each quotient coefficient.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    pub fn div_rem_with(
        &self,
        a: &Dense<F::Elem>,
        b: &Dense<F::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<QuoRem<F::Elem>, Cancelled> {
        let f = &self.base;
        let lead = f
            .inv(b.lead().expect("attempt to divide by zero"))
            .expect("a non-zero leading coefficient");
        let n = b.coeffs.len() - 1;
        if a.coeffs.len() <= n {
            return Ok((Dense::default(), a.clone()));
        }
        let mut r = a.coeffs.clone();
        let mut q = vec![f.zero(); r.len() - n];
        for k in (0..q.len()).rev() {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            let c = f.mul(&r[n + k], &lead);
            if !f.is_zero(&c) {
                for (rj, bj) in r[k..n + k].iter_mut().zip(&b.coeffs) {
                    *rj = f.sub(rj, &f.mul(&c, bj));
                }
            }
            q[k] = c;
        }
        r.truncate(n);
        Ok((self.poly(q), self.poly(r)))
    }
}

impl<R: Coefficients> DenseRing<R> {
    /// Pseudo-division (Knuth's Algorithm R of §4.6.1): `(q, r)` with
    /// `lc(b)^(deg a - deg b + 1) a = q b + r` and `deg r < deg b`, without
    /// dividing in `R`, so over any ring: polynomials in other variables
    /// too. `(0, a)` when `deg a < deg b`.
    ///
    /// Algorithm R multiplies the whole remainder and the quotient so far by
    /// `lc(b)` at each of its `s = deg a - deg b + 1` rounds, `O(s deg a)`
    /// products even when `lc(b)` is 1. Here a round multiplies only the
    /// `deg b` coefficients it subtracts from: a coefficient of `a` takes the
    /// power of `lc(b)` it is owed when the rounds reach it, and the
    /// quotient's coefficients theirs at the end. `O(s deg b)` products, as
    /// for a division, and none by `lc(b)` when it is 1 (CAS-172).
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[must_use]
    pub fn pseudo_div_rem(&self, a: &Dense<R::Elem>, b: &Dense<R::Elem>) -> QuoRem<R::Elem> {
        let r_ = &self.base;
        let d = b.lead().expect("attempt to divide by zero");
        let n = b.coeffs.len() - 1;
        if a.coeffs.len() <= n {
            return (Dense::default(), a.clone());
        }
        let s = a.coeffs.len() - n;
        if n == 0 {
            // d^s a = (d^(s-1) a) d.
            let power = r_.pow(d, (s - 1) as u64);
            return (self.scale(a, &power), Dense::default());
        }
        let monic = r_.is_one(d);
        let mut r = a.coeffs.clone();
        let mut q = vec![r_.zero(); s];
        // Round k, from s - 1 down to 0, is Algorithm R's: q <- d q + c x^k
        // and r <- d r - c x^k b for c = r_(n+k), keeping d^(s-k) a = q b +
        // r. It reads r_(n+k) and changes only r_k to r_(n+k-1), so only
        // those are scaled: r_(k+1) to r_(n+k-1) by d, as the round before
        // brought them up to date, and r_k, still a_k, by the d^(s-k) it is
        // owed. The coefficients below wait for their round; q_k takes its
        // d^k at the end.
        let mut power = r_.one();
        for k in (0..s).rev() {
            let c = r[n + k].clone();
            if !monic {
                for rj in &mut r[k + 1..n + k] {
                    *rj = r_.mul(d, rj);
                }
                power = r_.mul(&power, d);
                r[k] = r_.mul(&power, &r[k]);
            }
            if !r_.is_zero(&c) {
                for (rj, bj) in r[k..n + k].iter_mut().zip(&b.coeffs) {
                    *rj = r_.sub(rj, &r_.mul(&c, bj));
                }
            }
            q[k] = c;
        }
        if !monic {
            let mut power = r_.one();
            for qk in &mut q[1..] {
                power = r_.mul(&power, d);
                *qk = r_.mul(&power, qk);
            }
        }
        r.truncate(n);
        (self.poly(q), self.poly(r))
    }
}

impl<R: Coefficients + EuclideanDomain> DenseRing<R> {
    /// `a / b` when `b` divides `a` in `R[x]`; `None` otherwise. Long
    /// division in which each quotient coefficient must divide exactly.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[must_use]
    pub fn div_exact(&self, a: &Dense<R::Elem>, b: &Dense<R::Elem>) -> Option<Dense<R::Elem>> {
        let r_ = &self.base;
        let d = b.lead().expect("attempt to divide by zero");
        if a.is_zero() {
            return Some(Dense::default());
        }
        let n = b.coeffs.len() - 1;
        if a.coeffs.len() <= n {
            return None;
        }
        let mut r = a.coeffs.clone();
        let mut q = vec![r_.zero(); r.len() - n];
        for k in (0..q.len()).rev() {
            let (c, rest) = r_.div_rem(&r[n + k], d);
            if !r_.is_zero(&rest) {
                return None;
            }
            if !r_.is_zero(&c) {
                for (rj, bj) in r[k..n + k].iter_mut().zip(&b.coeffs) {
                    *rj = r_.sub(rj, &r_.mul(&c, bj));
                }
            }
            q[k] = c;
        }
        r[..n].iter().all(|c| r_.is_zero(c)).then(|| self.poly(q))
    }
}
