//! Division with explicit rounding, exact division and divisibility.

use std::ops::{Div, DivAssign, Rem, RemAssign};

use super::ops::integer_ops;
use super::{Integer, Repr};

const DIV_BY_ZERO: &str = "attempt to divide by zero";

/// How a quotient is rounded; the remainder is then `a - q d`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Rounding {
    /// Towards zero: the remainder has the sign of the dividend (`/` and `%`).
    Trunc,
    /// Towards negative infinity: the remainder has the sign of the divisor.
    Floor,
    /// Towards positive infinity: the remainder has the opposite sign to the
    /// divisor.
    Ceil,
    /// The remainder is non-negative.
    Euclid,
}

/// The adjustment from a truncated quotient to `rounding`: `Some(true)` for
/// `q -= 1, r += d`, `Some(false)` for `q += 1, r -= d`.
fn adjustment(
    rounding: Rounding,
    r_negative: bool,
    r_zero: bool,
    d_negative: bool,
) -> Option<bool> {
    if r_zero {
        return None;
    }
    match rounding {
        Rounding::Trunc => None,
        Rounding::Floor => (r_negative != d_negative).then_some(true),
        Rounding::Ceil => (r_negative == d_negative).then_some(false),
        Rounding::Euclid => r_negative.then_some(!d_negative),
    }
}

impl Integer {
    /// `(q, r)` with `self = q d + r`, `q` rounded as `rounding` says, or `None`
    /// if `d` is zero.
    #[must_use]
    pub fn checked_div_rem(&self, d: &Integer, rounding: Rounding) -> Option<(Integer, Integer)> {
        if d.is_zero() {
            return None;
        }
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &d.0) {
            // Only i64::MIN / -1 overflows; it goes the general way.
            if let (Some(mut q), Some(mut r)) = (a.checked_div(*b), a.checked_rem(*b)) {
                // |r| < |d| and |q| <= |a|, and q is only at an extreme when
                // |d| = 1 (r = 0), so none of these overflow.
                match adjustment(rounding, r < 0, r == 0, *b < 0) {
                    Some(true) => {
                        q -= 1;
                        r += b;
                    }
                    Some(false) => {
                        q += 1;
                        r -= b;
                    }
                    None => {}
                }
                return Some((Integer::from(q), Integer::from(r)));
            }
        }
        let (a_negative, a) = self.sign_mag();
        let (d_negative, dm) = d.sign_mag();
        let (q, r) = a.div_rem(&dm);
        let mut q = Integer::from_parts(a_negative != d_negative, q);
        let mut r = Integer::from_parts(a_negative, r);
        match adjustment(rounding, r.is_negative(), r.is_zero(), d_negative) {
            Some(true) => {
                q -= &Integer::ONE;
                r += d;
            }
            Some(false) => {
                q += &Integer::ONE;
                r -= d;
            }
            None => {}
        }
        Some((q, r))
    }

    fn div_rem_rounded(&self, d: &Integer, rounding: Rounding) -> (Integer, Integer) {
        self.checked_div_rem(d, rounding).expect(DIV_BY_ZERO)
    }

    /// `(q, r)` with `q` rounded towards zero (as `/` and `%`).
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem_trunc(&self, d: &Integer) -> (Integer, Integer) {
        self.div_rem_rounded(d, Rounding::Trunc)
    }

    /// `(q, r)` with `q` rounded towards negative infinity.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem_floor(&self, d: &Integer) -> (Integer, Integer) {
        self.div_rem_rounded(d, Rounding::Floor)
    }

    /// `(q, r)` with `q` rounded towards positive infinity.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem_ceil(&self, d: &Integer) -> (Integer, Integer) {
        self.div_rem_rounded(d, Rounding::Ceil)
    }

    /// `(q, r)` with `0 <= r < |d|`.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem_euclid(&self, d: &Integer) -> (Integer, Integer) {
        self.div_rem_rounded(d, Rounding::Euclid)
    }

    /// `self / d` where `d` is known to divide `self`.
    ///
    /// # Panics
    ///
    /// If `d` is zero or does not divide `self`.
    #[must_use]
    pub fn div_exact(&self, d: &Integer) -> Integer {
        assert!(!d.is_zero(), "{DIV_BY_ZERO}");
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &d.0) {
            assert!(a.wrapping_rem(*b) == 0, "division is not exact");
            if let Some(q) = a.checked_div(*b) {
                return Integer::from(q);
            }
        }
        let (a_negative, a) = self.sign_mag();
        let (d_negative, dm) = d.sign_mag();
        Integer::from_parts(a_negative != d_negative, a.div_exact(&dm))
    }

    /// `self / d` for a non-zero `d` that divides `self` by construction:
    /// [`Integer::div_exact`] without its check (see
    /// [`Natural::div_exact_unchecked`]).
    pub(crate) fn div_exact_unchecked(&self, d: &Integer) -> Integer {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &d.0) {
            if let Some(q) = a.checked_div(*b) {
                debug_assert!(a.wrapping_rem(*b) == 0, "division is not exact");
                return Integer::from(q);
            }
        }
        let (a_negative, a) = self.sign_mag();
        let (d_negative, dm) = d.sign_mag();
        Integer::from_parts(a_negative != d_negative, a.div_exact_unchecked(&dm))
    }

    /// Whether `d` divides `self`. Zero divides only zero.
    #[must_use]
    pub fn is_divisible_by(&self, d: &Integer) -> bool {
        match (&self.0, &d.0) {
            (Repr::Small(a), Repr::Small(0)) => *a == 0,
            // i64::MIN % -1 wraps to 0, which is the right answer.
            (Repr::Small(a), Repr::Small(b)) => a.wrapping_rem(*b) == 0,
            _ => self.sign_mag().1.is_divisible_by(&d.sign_mag().1),
        }
    }
}

impl Div<&Integer> for &Integer {
    type Output = Integer;

    /// Truncating division.
    ///
    /// # Panics
    ///
    /// If `rhs` is zero.
    fn div(self, rhs: &Integer) -> Integer {
        self.div_rem_trunc(rhs).0
    }
}

impl Rem<&Integer> for &Integer {
    type Output = Integer;

    /// Remainder of truncating division (the sign of `self`).
    ///
    /// # Panics
    ///
    /// If `rhs` is zero.
    fn rem(self, rhs: &Integer) -> Integer {
        self.div_rem_trunc(rhs).1
    }
}

impl DivAssign<&Integer> for Integer {
    fn div_assign(&mut self, rhs: &Integer) {
        *self = &*self / rhs;
    }
}

impl RemAssign<&Integer> for Integer {
    fn rem_assign(&mut self, rhs: &Integer) {
        *self = &*self % rhs;
    }
}

integer_ops! {
    Div::div, DivAssign::div_assign, |a, b| a / &b;
    Rem::rem, RemAssign::rem_assign, |a, b| a % &b;
}
