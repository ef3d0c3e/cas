//! Negation, addition, subtraction, multiplication and division.
#![allow(
    clippy::many_single_char_names,
    reason = "a + bi and c + di are the textbook names for the operands"
)]

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::Complex;

impl Complex {
    /// `self / rhs`, or `None` if `rhs` is zero: `self conj(rhs) / norm(rhs)`.
    #[must_use]
    pub fn checked_div(&self, rhs: &Complex) -> Option<Complex> {
        let n = rhs.norm();
        let (a, b, c, d) = (&self.re, &self.im, &rhs.re, &rhs.im);
        let re = (&(a * c) + &(b * d)).checked_div(&n)?;
        let im = (&(b * c) - &(a * d)).checked_div(&n)?;
        Some(Complex { re, im })
    }
}

impl Neg for Complex {
    type Output = Complex;

    fn neg(self) -> Complex {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl Neg for &Complex {
    type Output = Complex;

    fn neg(self) -> Complex {
        -self.clone()
    }
}

impl Add<&Complex> for &Complex {
    type Output = Complex;

    fn add(self, rhs: &Complex) -> Complex {
        Complex {
            re: &self.re + &rhs.re,
            im: &self.im + &rhs.im,
        }
    }
}

impl Sub<&Complex> for &Complex {
    type Output = Complex;

    fn sub(self, rhs: &Complex) -> Complex {
        Complex {
            re: &self.re - &rhs.re,
            im: &self.im - &rhs.im,
        }
    }
}

impl Mul<&Complex> for &Complex {
    type Output = Complex;

    /// `(a + bi)(c + di) = (ac - bd) + (ad + bc)i`, with shortcuts when a
    /// factor is real (the common case of scaling by a rational).
    fn mul(self, rhs: &Complex) -> Complex {
        let (a, b, c, d) = (&self.re, &self.im, &rhs.re, &rhs.im);
        if b.is_zero() {
            return Complex {
                re: a * c,
                im: a * d,
            };
        }
        if d.is_zero() {
            return Complex {
                re: a * c,
                im: b * c,
            };
        }
        Complex {
            re: &(a * c) - &(b * d),
            im: &(a * d) + &(b * c),
        }
    }
}

impl Div<&Complex> for &Complex {
    type Output = Complex;

    /// # Panics
    ///
    /// If `rhs` is zero.
    fn div(self, rhs: &Complex) -> Complex {
        self.checked_div(rhs).expect("attempt to divide by zero")
    }
}

/// Implements the remaining owned/borrowed operator combinations, and both
/// assignment forms, from `&a op &b`.
macro_rules! complex_ops {
    ($($op:ident::$f:ident, $assign:ident::$fa:ident;)*) => {$(
        impl $op<Complex> for Complex {
            type Output = Complex;

            fn $f(self, rhs: Complex) -> Complex {
                $op::$f(&self, &rhs)
            }
        }

        impl $op<&Complex> for Complex {
            type Output = Complex;

            fn $f(self, rhs: &Complex) -> Complex {
                $op::$f(&self, rhs)
            }
        }

        impl $op<Complex> for &Complex {
            type Output = Complex;

            fn $f(self, rhs: Complex) -> Complex {
                $op::$f(self, &rhs)
            }
        }

        impl $assign<&Complex> for Complex {
            fn $fa(&mut self, rhs: &Complex) {
                *self = $op::$f(&*self, rhs);
            }
        }

        impl $assign<Complex> for Complex {
            fn $fa(&mut self, rhs: Complex) {
                *self = $op::$f(&*self, &rhs);
            }
        }
    )*};
}

complex_ops! {
    Add::add, AddAssign::add_assign;
    Sub::sub, SubAssign::sub_assign;
    Mul::mul, MulAssign::mul_assign;
    Div::div, DivAssign::div_assign;
}
