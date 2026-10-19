//! Greatest common divisor, extended GCD, least common multiple and modular
//! inverse.

use std::borrow::Cow;
use std::cmp::Ordering;

use super::{Integer, Repr};
use crate::{Natural, kernel};

#[cfg(test)]
mod tests;

impl Integer {
    /// The greatest common divisor, which is non-negative; `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, other: &Integer) -> Integer {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &other.0) {
            return Integer::from(kernel::gcd_1(a.unsigned_abs(), b.unsigned_abs()));
        }
        Integer::from(self.sign_mag().1.gcd(&other.sign_mag().1))
    }

    /// The least common multiple, which is non-negative; zero if either
    /// operand is zero.
    #[must_use]
    pub fn lcm(&self, other: &Integer) -> Integer {
        Integer::from(self.sign_mag().1.lcm(&other.sign_mag().1))
    }

    /// `(g, s, t)` with `g = gcd(self, other) = s self + t other`.
    ///
    /// The cofactors are the ones GMP's `mpz_gcdext` returns: with `a = self`
    /// and `b = other`, `|s| < |b| / (2g)` and `|t| < |a| / (2g)`, which define
    /// them uniquely, except that
    /// - `s = 0, t = sgn(b)` if `|a| = |b|` or `b` divides `a`;
    /// - `s = sgn(a)` if `b = 0` or `|b| = 2g`;
    /// - `t = sgn(b)` if `a = 0` or `|a| = 2g`;
    /// - `(0, 0, 0)` for `a = b = 0`.
    #[must_use]
    #[allow(
        clippy::many_single_char_names,
        reason = "g, s, t, a, b, m are the documented names"
    )]
    pub fn gcd_ext(&self, other: &Integer) -> (Integer, Integer, Integer) {
        let (a_negative, a) = self.sign_mag();
        let (b_negative, b) = other.sign_mag();
        let (g, s) = a.gcd_cofactor(&b);
        if g.is_zero() {
            return (Integer::ZERO, Integer::ZERO, Integer::ZERO);
        }
        if b.is_zero() {
            return (Integer::from(g), signum(self), Integer::ZERO);
        }
        // The sign of s as the coefficient of `self` rather than of |self|.
        let (s_negative, s) = s.into_parts();
        let s_negative = s_negative != a_negative;
        // Every coefficient of `self` is congruent to s modulo m = |b| / g;
        // pick the one nearest zero, which the tracked one usually is.
        let m = if g.is_one() {
            Cow::Borrowed(&*b)
        } else {
            Cow::Owned(b.div_exact(&g))
        };
        let (s_negative, s) = if m.limbs() == [2] {
            // s is odd, and +1 and -1 are equally near: GMP takes sgn(a).
            (a_negative, Natural::one())
        } else if below_half(&s, &m) {
            (s_negative, s)
        } else {
            // The residues ±r and ∓(m - r); 2r = m only for m = 2.
            let r = s.div_rem(&m).1;
            if below_half(&r, &m) {
                (s_negative, r)
            } else {
                (!s_negative, &*m - &r)
            }
        };
        // t = (g - s self) / other: p = |g - s self|, negative its sign.
        let mut p = &s * &*a;
        let negative = if s_negative == a_negative && !s.is_zero() {
            p -= &g;
            true
        } else {
            p += &g;
            false
        };
        let t = Integer::from_parts(negative != b_negative, p.div_exact(&b));
        (Integer::from(g), Integer::from_parts(s_negative, s), t)
    }

    /// The inverse of `self` modulo `m`, in `0 .. |m|`, or `None` if `m` is
    /// zero or `self` and `m` are not coprime. Modulo `±1` the inverse is 0.
    #[must_use]
    pub fn mod_inverse(&self, m: &Integer) -> Option<Integer> {
        if m.is_zero() {
            return None;
        }
        let (g, s, _) = self.gcd_ext(m);
        g.is_one().then(|| s.div_rem_euclid(m).1)
    }
}

/// Whether `2 x < m`.
fn below_half(x: &Natural, m: &Natural) -> bool {
    let (xb, mb) = (x.bit_length(), m.bit_length());
    if xb + 1 != mb {
        return xb < mb;
    }
    (x + x) < *m
}

/// `-1`, `0` or `1`.
fn signum(x: &Integer) -> Integer {
    Integer::from(match x.sign() {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    })
}
