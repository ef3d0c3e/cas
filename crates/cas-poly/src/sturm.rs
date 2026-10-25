//! Sturm's sequences (CAS-197, CAS-198): the number of distinct real roots
//! of a polynomial over ℚ, on the line or between rational bounds.
//!
//! For `p` square-free, the sequence is `p₀ = p`, `p₁ = p′`, then the
//! negated remainders `pₖ₊₁ = −(pₖ₋₁ mod pₖ)` down to a constant. With `V(x)`
//! the number of sign changes of the `pₖ(x)`, zeros left out, Sturm's
//! theorem counts the roots in `(a, b]` as `V(a) − V(b)`. At a root `r` of
//! `p`, `p₁(r) ≠ 0` gives the sign `p` takes just after `r`, so `V(r)` is
//! `V` just after `r`: the count holds whether `a` and `b` are roots or
//! not, and a root at `a` itself is counted apart. At `±∞` each sign is
//! that of the leading coefficient, times `(−1)^deg` at `−∞`.

use cas_num::Rational;
use cas_num::algebra::{EuclideanDomain, Rationals};

use crate::dense::{Dense, DenseRing};

/// Sturm's sequence of a square-free `p` of degree at least 1, each
/// remainder over the absolute value of its leading coefficient, which
/// keeps the signs and the coefficients small.
fn sequence(ring: &DenseRing<Rationals>, p: &Dense<Rational>) -> Vec<Dense<Rational>> {
    let mut sequence = vec![p.clone(), ring.derivative(p)];
    loop {
        let n = sequence.len();
        let r = ring.div_rem(&sequence[n - 2], &sequence[n - 1]).1;
        let Some(lead) = r.lead() else {
            return sequence;
        };
        let scale = lead.abs().checked_recip().expect("not 0");
        sequence.push(ring.scale(&r, &-scale));
    }
}

/// The sign changes in `signs`, zeros left out.
fn changes(signs: impl Iterator<Item = std::cmp::Ordering>) -> usize {
    let mut last = None;
    let mut count = 0;
    for s in signs.filter(|s| s.is_ne()) {
        if last.is_some_and(|l| l != s) {
            count += 1;
        }
        last = Some(s);
    }
    count
}

/// Where a sequence's signs are taken: a rational, or `−∞` or `+∞`.
#[derive(Clone, Copy)]
enum At<'a> {
    Point(&'a Rational),
    Infinity { negative: bool },
}

/// `V(at)`.
fn sign_changes(ring: &DenseRing<Rationals>, sequence: &[Dense<Rational>], at: At<'_>) -> usize {
    let zero = Rational::ZERO;
    changes(sequence.iter().map(|q| match at {
        At::Point(x) => ring.eval(q, x).cmp(&zero),
        At::Infinity { negative } => {
            let lead = q.lead().map_or(std::cmp::Ordering::Equal, |c| c.cmp(&zero));
            if negative && q.degree().unwrap_or(0) % 2 == 1 {
                lead.reverse()
            } else {
                lead
            }
        }
    }))
}

impl DenseRing<Rationals> {
    /// The number of distinct real roots of `p` in the closed interval from
    /// `lo` to `hi`, `None` for an infinite end, by Sturm's sequence of the
    /// square-free part of `p`; 0 for a constant, and when `lo > hi`.
    ///
    /// ```
    /// use cas_num::algebra::Rationals;
    /// use cas_num::{Integer, Rational};
    /// use cas_poly::DenseRing;
    ///
    /// let qx = DenseRing::new(Rationals);
    /// // x² − 2: two real roots, ±√2, one of them in [0, 2].
    /// let p = qx.from_i64s(&[-2, 0, 1]);
    /// let (zero, two) = (Rational::from(Integer::from(0)), Rational::from(Integer::from(2)));
    /// assert_eq!(qx.count_real_roots(&p, None, None), 2);
    /// assert_eq!(qx.count_real_roots(&p, Some(&zero), Some(&two)), 1);
    /// ```
    ///
    /// # Panics
    ///
    /// If `p` is zero.
    #[must_use]
    pub fn count_real_roots(
        &self,
        p: &Dense<Rational>,
        lo: Option<&Rational>,
        hi: Option<&Rational>,
    ) -> usize {
        assert!(!p.is_zero(), "the roots of 0");
        if lo.zip(hi).is_some_and(|(lo, hi)| lo > hi) {
            return 0;
        }
        let g = EuclideanDomain::gcd(self, p, &self.derivative(p));
        let q = self.div_rem(p, &g).0;
        if q.degree().unwrap_or(0) == 0 {
            return 0;
        }
        let sequence = sequence(self, &q);
        let at_lo = sign_changes(
            self,
            &sequence,
            lo.map_or(At::Infinity { negative: true }, At::Point),
        );
        let at_hi = sign_changes(
            self,
            &sequence,
            hi.map_or(At::Infinity { negative: false }, At::Point),
        );
        let root_at_lo = lo.is_some_and(|lo| self.eval(&q, lo).is_zero());
        at_lo - at_hi + usize::from(root_at_lo)
    }
}

#[cfg(test)]
mod tests;
