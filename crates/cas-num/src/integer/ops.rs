//! Negation, addition, subtraction and multiplication.
//!
//! Two inline values are combined with `checked_*` arithmetic; on overflow the
//! exact result comes from `i128`. Otherwise the operation runs on sign and
//! magnitude, reusing an owned box where it can, and the result is demoted to
//! the inline form when it fits.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::{Big, Integer, Repr, small_value};
use crate::Natural;

impl Integer {
    /// The canonical form of a boxed value, which may fit an `i64` after an
    /// operation; keeps the box otherwise.
    pub(super) fn from_box(b: Box<Big>) -> Integer {
        match small_value(b.negative, &b.mag) {
            Some(x) => Integer(Repr::Small(x)),
            None => Integer(Repr::Big(b)),
        }
    }

    /// The value as a (possibly non-canonical) box, to operate on in place.
    fn into_box(self) -> Box<Big> {
        match self.0 {
            Repr::Small(x) => Box::new(Big {
                negative: x < 0,
                mag: Natural::from(x.unsigned_abs()),
            }),
            Repr::Big(b) => b,
        }
    }

    /// `self += rhs`, or `self -= rhs` if `subtract`.
    fn add_assign_signed(&mut self, rhs: &Integer, subtract: bool) {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &rhs.0) {
            let (a, b) = (*a, *b);
            let small = if subtract {
                a.checked_sub(b)
            } else {
                a.checked_add(b)
            };
            *self = match small {
                Some(x) => Integer(Repr::Small(x)),
                None if subtract => Integer::from(i128::from(a) - i128::from(b)),
                None => Integer::from(i128::from(a) + i128::from(b)),
            };
            return;
        }
        let mut big = std::mem::take(self).into_box();
        let (rhs_negative, rhs_mag) = rhs.sign_mag();
        let rhs_negative = rhs_negative != subtract;
        if big.negative == rhs_negative {
            big.mag += &*rhs_mag;
        } else if big.mag >= *rhs_mag {
            big.mag -= &*rhs_mag;
        } else {
            big.mag = &*rhs_mag - &big.mag;
            big.negative = rhs_negative;
        }
        *self = Integer::from_box(big);
    }
}

impl Neg for Integer {
    type Output = Integer;

    fn neg(self) -> Integer {
        match self.0 {
            Repr::Small(x) => match x.checked_neg() {
                Some(y) => Integer(Repr::Small(y)),
                None => Integer::from(-i128::from(x)),
            },
            Repr::Big(mut b) => {
                b.negative = !b.negative;
                // -(2^63) is the one boxed value whose negation fits.
                Integer::from_box(b)
            }
        }
    }
}

impl Neg for &Integer {
    type Output = Integer;

    fn neg(self) -> Integer {
        -self.clone()
    }
}

impl AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, rhs: &Integer) {
        self.add_assign_signed(rhs, false);
    }
}

impl SubAssign<&Integer> for Integer {
    fn sub_assign(&mut self, rhs: &Integer) {
        self.add_assign_signed(rhs, true);
    }
}

impl Add<&Integer> for &Integer {
    type Output = Integer;

    fn add(self, rhs: &Integer) -> Integer {
        // Start from the boxed operand, whose buffer the result needs anyway.
        let (mut r, other) = match (&self.0, &rhs.0) {
            (Repr::Small(_), Repr::Big(_)) => (rhs.clone(), self),
            _ => (self.clone(), rhs),
        };
        r += other;
        r
    }
}

impl Sub<&Integer> for &Integer {
    type Output = Integer;

    fn sub(self, rhs: &Integer) -> Integer {
        if let (Repr::Small(_), Repr::Big(_)) = (&self.0, &rhs.0) {
            let mut r = -rhs;
            r += self;
            r
        } else {
            let mut r = self.clone();
            r -= rhs;
            r
        }
    }
}

impl Mul<&Integer> for &Integer {
    type Output = Integer;

    fn mul(self, rhs: &Integer) -> Integer {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &rhs.0) {
            return match a.checked_mul(*b) {
                Some(x) => Integer(Repr::Small(x)),
                None => Integer::from(i128::from(*a) * i128::from(*b)),
            };
        }
        let (a_negative, a) = self.sign_mag();
        let (b_negative, b) = rhs.sign_mag();
        Integer::from_parts(a_negative != b_negative, &*a * &*b)
    }
}

impl MulAssign<&Integer> for Integer {
    /// The product needs a fresh buffer, so this is `*self = &*self * rhs`.
    fn mul_assign(&mut self, rhs: &Integer) {
        *self = &*self * rhs;
    }
}

/// Implements the remaining owned/borrowed operator combinations from
/// `&a op &b` and `a op= &b`; `$by_value_rhs` computes `&a op b` reusing `b`.
macro_rules! integer_ops {
    ($($op:ident::$f:ident, $assign:ident::$fa:ident, $by_value_rhs:expr;)*) => {$(
        impl $assign<Integer> for Integer {
            fn $fa(&mut self, rhs: Integer) {
                $assign::$fa(self, &rhs);
            }
        }

        impl $op<&Integer> for Integer {
            type Output = Integer;

            fn $f(mut self, rhs: &Integer) -> Integer {
                $assign::$fa(&mut self, rhs);
                self
            }
        }

        impl $op<Integer> for Integer {
            type Output = Integer;

            fn $f(mut self, rhs: Integer) -> Integer {
                $assign::$fa(&mut self, &rhs);
                self
            }
        }

        impl $op<Integer> for &Integer {
            type Output = Integer;

            fn $f(self, rhs: Integer) -> Integer {
                let by_value_rhs: fn(&Integer, Integer) -> Integer = $by_value_rhs;
                by_value_rhs(self, rhs)
            }
        }
    )*};
}

pub(super) use integer_ops;

integer_ops! {
    Add::add, AddAssign::add_assign, |a, mut b| { b += a; b };
    // a - b = -(b - a), negated in place.
    Sub::sub, SubAssign::sub_assign, |a, mut b| { b -= a; -b };
    Mul::mul, MulAssign::mul_assign, |a, b| a * &b;
}
