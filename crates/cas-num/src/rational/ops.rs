//! Negation, addition, subtraction, multiplication and division.
//!
//! Sums and products follow Knuth §4.5.1 (Henrici): common factors are
//! cancelled before multiplying, so every result comes out reduced after at
//! most one more GCD of operand-sized values, and intermediate products never
//! grow past the size of the result.
#![allow(
    clippy::many_single_char_names,
    reason = "a/b and c/d are Knuth's names for the operands"
)]

use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::Rational;
use crate::Integer;

impl Rational {
    /// `self / rhs`, or `None` if `rhs` is zero.
    #[must_use]
    pub fn checked_div(&self, rhs: &Rational) -> Option<Rational> {
        // Multiply by the reciprocal, whose sign goes to its numerator.
        match rhs.num.sign() {
            Ordering::Equal => None,
            Ordering::Greater => Some(product(&self.num, &self.den, &rhs.den, &rhs.num)),
            Ordering::Less => Some(product(&self.num, &self.den, &-&rhs.den, &-&rhs.num)),
        }
    }
}

/// `a/b + c/d` for canonical fractions.
fn sum(a: &Integer, b: &Integer, c: &Integer, d: &Integer) -> Rational {
    if b == d {
        // Common for integers (b = d = 1); only the new numerator can share
        // a factor with the denominator.
        return Rational::checked_new(a + c, b.clone()).expect("non-zero denominator");
    }
    let g = b.gcd(d);
    if g.is_one() {
        // gcd(ad + bc, bd) = 1 when b and d are coprime.
        return Rational::from_canonical(a * d + b * c, b * d);
    }
    let (b_g, d_g) = (b.div_exact(&g), d.div_exact(&g));
    let t = a * &d_g + c * &b_g;
    // t = 0 means a/b = -c/d, so b = d, which returned above.
    debug_assert!(!t.is_zero());
    // Any common factor of t and the denominator b_g d is one of t and g.
    let g2 = t.gcd(&g);
    if g2.is_one() {
        Rational::from_canonical(t, b_g * d)
    } else {
        Rational::from_canonical(t.div_exact(&g2), b_g * d.div_exact(&g2))
    }
}

/// `(a/b) (c/d)` for canonical fractions.
fn product(a: &Integer, b: &Integer, c: &Integer, d: &Integer) -> Rational {
    if a.is_zero() || c.is_zero() {
        return Rational::ZERO;
    }
    // a is coprime to b and c to d, so only these cross factors cancel.
    let g1 = a.gcd(d);
    let g2 = c.gcd(b);
    let reduce = |x: &Integer, g: &Integer| {
        if g.is_one() {
            x.clone()
        } else {
            x.div_exact(g)
        }
    };
    Rational::from_canonical(
        reduce(a, &g1) * reduce(c, &g2),
        reduce(b, &g2) * reduce(d, &g1),
    )
}

impl Neg for Rational {
    type Output = Rational;

    fn neg(self) -> Rational {
        Rational::from_canonical(-self.num, self.den)
    }
}

impl Neg for &Rational {
    type Output = Rational;

    fn neg(self) -> Rational {
        -self.clone()
    }
}

impl Add<&Rational> for &Rational {
    type Output = Rational;

    fn add(self, rhs: &Rational) -> Rational {
        sum(&self.num, &self.den, &rhs.num, &rhs.den)
    }
}

impl Sub<&Rational> for &Rational {
    type Output = Rational;

    fn sub(self, rhs: &Rational) -> Rational {
        sum(&self.num, &self.den, &-&rhs.num, &rhs.den)
    }
}

impl Mul<&Rational> for &Rational {
    type Output = Rational;

    fn mul(self, rhs: &Rational) -> Rational {
        product(&self.num, &self.den, &rhs.num, &rhs.den)
    }
}

impl Div<&Rational> for &Rational {
    type Output = Rational;

    /// # Panics
    ///
    /// If `rhs` is zero.
    fn div(self, rhs: &Rational) -> Rational {
        self.checked_div(rhs).expect("attempt to divide by zero")
    }
}

/// Implements the remaining owned/borrowed operator combinations, and both
/// assignment forms, from `&a op &b`.
macro_rules! rational_ops {
    ($($op:ident::$f:ident, $assign:ident::$fa:ident;)*) => {$(
        impl $op<Rational> for Rational {
            type Output = Rational;

            fn $f(self, rhs: Rational) -> Rational {
                $op::$f(&self, &rhs)
            }
        }

        impl $op<&Rational> for Rational {
            type Output = Rational;

            fn $f(self, rhs: &Rational) -> Rational {
                $op::$f(&self, rhs)
            }
        }

        impl $op<Rational> for &Rational {
            type Output = Rational;

            fn $f(self, rhs: Rational) -> Rational {
                $op::$f(self, &rhs)
            }
        }

        impl $assign<&Rational> for Rational {
            fn $fa(&mut self, rhs: &Rational) {
                *self = $op::$f(&*self, rhs);
            }
        }

        impl $assign<Rational> for Rational {
            fn $fa(&mut self, rhs: Rational) {
                *self = $op::$f(&*self, &rhs);
            }
        }
    )*};
}

rational_ops! {
    Add::add, AddAssign::add_assign;
    Sub::sub, SubAssign::sub_assign;
    Mul::mul, MulAssign::mul_assign;
    Div::div, DivAssign::div_assign;
}
