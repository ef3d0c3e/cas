//! Addition and subtraction.

use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Sub, SubAssign};

use super::Natural;
use super::limb_vec::{INLINE, LimbVec};
use crate::{Limb, kernel};

const SUB_OVERFLOW: &str = "attempt to subtract with overflow";

impl Natural {
    /// `self - rhs`, or `None` if `rhs > self`.
    #[must_use]
    pub fn checked_sub(&self, rhs: &Natural) -> Option<Natural> {
        let (a, b) = (self.limbs(), rhs.limbs());
        if kernel::cmp(a, b) == Ordering::Less {
            return None;
        }
        let mut r = LimbVec::zeroed(a.len());
        let borrow = kernel::sub(r.as_mut_slice(), a, b);
        debug_assert_eq!(borrow, 0);
        let mut r = Natural { limbs: r };
        r.normalize();
        r.debug_check();
        Some(r)
    }

    /// `self - rhs` for a single limb, or `None` if `rhs > self`.
    #[must_use]
    pub fn checked_sub_limb(&self, rhs: Limb) -> Option<Natural> {
        let mut r = self.clone();
        r.try_sub_assign_limb(rhs).then_some(r)
    }

    /// `self -= rhs` if `rhs <= self`, leaving `self` unchanged otherwise;
    /// returns whether the subtraction happened.
    pub fn checked_sub_assign(&mut self, rhs: &Natural) -> bool {
        if kernel::cmp(self.limbs(), rhs.limbs()) == Ordering::Less {
            return false;
        }
        let borrow = kernel::sub_assign(self.limbs.as_mut_slice(), rhs.limbs());
        debug_assert_eq!(borrow, 0);
        self.normalize();
        self.debug_check();
        true
    }

    fn try_sub_assign_limb(&mut self, rhs: Limb) -> bool {
        let fits = match self.limbs() {
            [] => rhs == 0,
            [x] => *x >= rhs,
            _ => true,
        };
        if fits {
            let borrow = kernel::sub_1_assign(self.limbs.as_mut_slice(), rhs);
            debug_assert_eq!(borrow, 0);
            self.normalize();
            self.debug_check();
        }
        fits
    }

    fn add_assign_limb(&mut self, rhs: Limb) {
        let carry = kernel::add_1_assign(self.limbs.as_mut_slice(), rhs);
        if carry != 0 {
            self.limbs.push(carry);
        }
        self.debug_check();
    }
}

impl Add<&Natural> for &Natural {
    type Output = Natural;

    fn add(self, rhs: &Natural) -> Natural {
        let (a, b) = if self.limb_len() >= rhs.limb_len() {
            (self.limbs(), rhs.limbs())
        } else {
            (rhs.limbs(), self.limbs())
        };
        let n = a.len();
        if n < INLINE {
            // A carry still fits inline.
            let mut r = LimbVec::zeroed(n);
            let carry = kernel::add(r.as_mut_slice(), a, b);
            if carry != 0 {
                r.push(carry);
            }
            return Natural::from_limb_vec(r);
        }
        if n == INLINE {
            // Through the stack, so that a sum carrying out of the inline
            // capacity allocates once, at its length.
            let mut buf = [0; INLINE + 1];
            buf[n] = kernel::add(&mut buf[..n], a, b);
            let len = n + usize::from(buf[n] != 0);
            return Natural::from_limb_vec(LimbVec::from_slice(&buf[..len]));
        }
        // Room for the carry limb up front: pushing it onto an exactly sized
        // heap buffer would reallocate and copy the whole sum.
        let mut r = LimbVec::zeroed(n + 1);
        let carry = kernel::add(&mut r.as_mut_slice()[..n], a, b);
        if carry != 0 {
            r.as_mut_slice()[n] = carry;
        } else {
            r.truncate(n);
        }
        Natural::from_limb_vec(r)
    }
}

impl AddAssign<&Natural> for Natural {
    fn add_assign(&mut self, rhs: &Natural) {
        if self.limb_len() < rhs.limb_len() {
            self.limbs.resize(rhs.limb_len());
        }
        let carry = kernel::add_assign(self.limbs.as_mut_slice(), rhs.limbs());
        if carry != 0 {
            self.limbs.push(carry);
        }
        self.debug_check();
    }
}

/// Implements the owned/borrowed operator combinations of a commutative
/// operation from `&a op &b` and `a op= &b`, reusing an owned buffer
/// (the longer one when both are owned).
macro_rules! commutative_ops {
    ($($op:ident::$f:ident, $assign:ident::$fa:ident;)*) => {$(
        impl $assign<Natural> for Natural {
            fn $fa(&mut self, rhs: Natural) {
                if self.limb_len() >= rhs.limb_len() {
                    $assign::$fa(self, &rhs);
                } else {
                    let mut rhs = rhs;
                    $assign::$fa(&mut rhs, &*self);
                    *self = rhs;
                }
            }
        }

        impl $op<&Natural> for Natural {
            type Output = Natural;

            fn $f(mut self, rhs: &Natural) -> Natural {
                $assign::$fa(&mut self, rhs);
                self
            }
        }

        impl $op<Natural> for &Natural {
            type Output = Natural;

            fn $f(self, mut rhs: Natural) -> Natural {
                $assign::$fa(&mut rhs, self);
                rhs
            }
        }

        impl $op<Natural> for Natural {
            type Output = Natural;

            fn $f(mut self, rhs: Natural) -> Natural {
                $assign::$fa(&mut self, rhs);
                self
            }
        }
    )*};
}

pub(super) use commutative_ops;

commutative_ops! {
    Add::add, AddAssign::add_assign;
}

impl Sub<&Natural> for &Natural {
    type Output = Natural;

    /// # Panics
    ///
    /// If `rhs > self`.
    fn sub(self, rhs: &Natural) -> Natural {
        self.checked_sub(rhs).expect(SUB_OVERFLOW)
    }
}

impl SubAssign<&Natural> for Natural {
    /// # Panics
    ///
    /// If `rhs > self`.
    fn sub_assign(&mut self, rhs: &Natural) {
        assert!(self.checked_sub_assign(rhs), "{SUB_OVERFLOW}");
    }
}

impl SubAssign<Natural> for Natural {
    fn sub_assign(&mut self, rhs: Natural) {
        *self -= &rhs;
    }
}

impl Sub<&Natural> for Natural {
    type Output = Natural;

    fn sub(mut self, rhs: &Natural) -> Natural {
        self -= rhs;
        self
    }
}

impl Sub<Natural> for Natural {
    type Output = Natural;

    fn sub(mut self, rhs: Natural) -> Natural {
        self -= &rhs;
        self
    }
}

impl Sub<Natural> for &Natural {
    type Output = Natural;

    fn sub(self, rhs: Natural) -> Natural {
        self - &rhs
    }
}

impl AddAssign<Limb> for Natural {
    fn add_assign(&mut self, rhs: Limb) {
        self.add_assign_limb(rhs);
    }
}

impl Add<Limb> for Natural {
    type Output = Natural;

    fn add(mut self, rhs: Limb) -> Natural {
        self += rhs;
        self
    }
}

impl Add<Limb> for &Natural {
    type Output = Natural;

    fn add(self, rhs: Limb) -> Natural {
        let a = self.limbs();
        let mut r = LimbVec::zeroed(a.len());
        let carry = kernel::add_1(r.as_mut_slice(), a, rhs);
        if carry != 0 {
            r.push(carry);
        }
        Natural::from_limb_vec(r)
    }
}

impl SubAssign<Limb> for Natural {
    /// # Panics
    ///
    /// If `rhs > self`.
    fn sub_assign(&mut self, rhs: Limb) {
        assert!(self.try_sub_assign_limb(rhs), "{SUB_OVERFLOW}");
    }
}

impl Sub<Limb> for Natural {
    type Output = Natural;

    fn sub(mut self, rhs: Limb) -> Natural {
        self -= rhs;
        self
    }
}

impl Sub<Limb> for &Natural {
    type Output = Natural;

    fn sub(self, rhs: Limb) -> Natural {
        self.checked_sub_limb(rhs).expect(SUB_OVERFLOW)
    }
}
