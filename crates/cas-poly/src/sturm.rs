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
//!
//! [`DenseRing::sturm_sequence`] gives giac's sequence instead (CAS-539):
//! subresultant remainders, integral over ℤ, with the signs corrected.

use cas_num::algebra::{EuclideanDomain, Rationals};
use cas_num::{Integer, Rational};

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

    /// giac's Sturm sequence of `p` (its `Tsturm_seq`, CAS-539): `p`, `p′`,
    /// then the subresultant remainders of the primitive part `a` of `p`
    /// and `p′`, down to a constant or to the gcd. With `δ` the degree drop
    /// and `b₀` the divisor's leading coefficient, the pseudo-remainder
    /// `b₀^(δ+1) a mod b` is divided by `−g h^δ` (`g` and `h` the
    /// subresultant's auxiliaries, 1 at first), which leaves a sign error:
    /// `h g` when `δ` is odd and `h ≠ g`, `b₀ g` when `δ` is even and
    /// `b₀ ≠ g`. The errors propagate by twos along the sequence, so two
    /// slots keep them, one for the even steps and one for the odd ones,
    /// and each remainder is multiplied by its slot's. Over the integers
    /// the entries stay integral: `x³ − 7x + 7` gives `x³ − 7x + 7`,
    /// `3x² − 7`, `42x − 63`, `49`. `[p]` for 0, `[p, 0]` for a constant.
    ///
    /// ```
    /// use cas_num::algebra::Rationals;
    /// use cas_poly::DenseRing;
    ///
    /// let qx = DenseRing::new(Rationals);
    /// let p = qx.from_i64s(&[7, -7, 0, 1]);
    /// let expected = [
    ///     p.clone(),
    ///     qx.from_i64s(&[-7, 0, 3]),
    ///     qx.from_i64s(&[-63, 42]),
    ///     qx.from_i64s(&[49]),
    /// ];
    /// assert_eq!(qx.sturm_sequence(&p), expected);
    /// ```
    ///
    /// # Panics
    ///
    /// Not for a degree that fits in memory.
    #[must_use]
    #[allow(
        clippy::many_single_char_names,
        reason = "giac's names: a over b, the auxiliaries g and h, the degrees m and n, the remainder r"
    )]
    pub fn sturm_sequence(&self, p: &Dense<Rational>) -> Vec<Dense<Rational>> {
        let mut sequence = vec![p.clone()];
        if p.is_zero() {
            return sequence;
        }
        let content = rational_content(p);
        let mut a = self.scale(p, &content.recip());
        let mut b = self.derivative(p);
        sequence.push(b.clone());
        let (mut g, mut h) = (Rational::ONE, Rational::ONE);
        let mut sign_error = [Rational::ONE, Rational::ONE];
        for step in 0.. {
            let m = a.degree().unwrap_or(0);
            let n = b.degree().unwrap_or(0);
            if n == 0 {
                break;
            }
            let drop = m - n;
            let exponent = i64::try_from(drop).expect("a degree in memory");
            let b0 = b.lead().expect("degree at least 1").clone();
            let r = self.div_rem(&self.scale(&a, &b0.pow(exponent + 1)), &b).1;
            if r.is_zero() {
                break;
            }
            a = b;
            let temp = h.pow(exponent);
            b = self.scale(&r, &(-(&g * &temp)).recip());
            let slot = &mut sign_error[step % 2];
            if drop % 2 == 1 && h != g {
                *slot = &*slot * &(&h * &g);
            }
            if drop % 2 == 0 && b0 != g {
                *slot = &*slot * &(&b0 * &g);
            }
            sequence.push(self.scale(&b, slot));
            h = if drop == 1 {
                b0.clone()
            } else {
                &(&b0.pow(exponent) * &h) / &temp
            };
            g = b0;
        }
        sequence
    }
}

/// The positive rational `c` that makes `p / c` primitive with integer
/// coefficients: the gcd of the numerators over the lcm of the
/// denominators.
fn rational_content(p: &Dense<Rational>) -> Rational {
    let (mut numer, mut denom) = (Integer::from(0), Integer::from(1));
    for c in p.coeffs() {
        numer = numer.gcd(c.numer());
        denom = denom.lcm(c.denom());
    }
    Rational::new(numer, denom)
}

#[cfg(test)]
mod tests;
