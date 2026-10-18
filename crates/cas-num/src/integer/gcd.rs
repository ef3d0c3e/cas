//! Greatest common divisor, extended GCD, least common multiple and modular
//! inverse.

use std::cmp::Ordering;

use super::{Integer, Repr};
use crate::kernel;

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
        let b = other.sign_mag().1;
        let (g, s) = a.gcd_cofactor(&b);
        if g.is_zero() {
            return (Integer::ZERO, Integer::ZERO, Integer::ZERO);
        }
        if b.is_zero() {
            return (Integer::from(g), signum(self), Integer::ZERO);
        }
        let s = if a_negative { -s } else { s };
        // Every coefficient of `self` is congruent to s modulo |b| / g; pick
        // the one nearest zero.
        let m = Integer::from(b.div_exact_unchecked(&g));
        let s = if m.to_i64() == Some(2) {
            // s is odd, and +1 and -1 are equally near: GMP takes sgn(a).
            signum(self)
        } else {
            let r = s.div_rem_euclid(&m).1;
            if &r + &r > m { r - &m } else { r }
        };
        let g = Integer::from(g);
        let t = (&g - &(&s * self)).div_exact_unchecked(other);
        (g, s, t)
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

/// `-1`, `0` or `1`.
fn signum(x: &Integer) -> Integer {
    Integer::from(match x.sign() {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    })
}
