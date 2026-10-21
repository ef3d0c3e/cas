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

impl<R: Coefficients + EuclideanDomain> DenseRing<R> {
    /// Pseudo-division (Knuth's Algorithm R of §4.6.1): `(q, r)` with
    /// `lc(b)^(deg a - deg b + 1) a = q b + r` and `deg r < deg b`, without
    /// dividing in `R`. `(0, a)` when `deg a < deg b`.
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
        let mut r = a.coeffs.clone();
        let mut q = vec![r_.zero(); r.len() - n];
        // Each round: q <- d q + c x^k, r <- d r - c x^k b, which keeps
        // d^s a = q b + r and clears r's coefficient of x^(n + k).
        for k in (0..q.len()).rev() {
            let c = r[n + k].clone();
            for qj in &mut q[k + 1..] {
                *qj = r_.mul(d, qj);
            }
            for rj in &mut r[..n + k] {
                *rj = r_.mul(d, rj);
            }
            for (rj, bj) in r[k..n + k].iter_mut().zip(&b.coeffs) {
                *rj = r_.sub(rj, &r_.mul(&c, bj));
            }
            q[k] = c;
        }
        r.truncate(n);
        (self.poly(q), self.poly(r))
    }

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
