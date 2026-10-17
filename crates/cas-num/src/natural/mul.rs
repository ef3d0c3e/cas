//! Multiplication.

use std::ops::{Mul, MulAssign};

use super::Natural;
use super::limb_vec::LimbVec;
use super::ops::commutative_ops;
use crate::{Limb, kernel};

impl Natural {
    /// `self^2`, faster than `self * self`.
    #[must_use]
    pub fn square(&self) -> Natural {
        let a = self.limbs();
        if a.is_empty() {
            return Natural::ZERO;
        }
        let mut r = LimbVec::zeroed(2 * a.len());
        let mut scratch = vec![0; kernel::sqr_scratch_len(a.len())];
        kernel::sqr(r.as_mut_slice(), a, &mut scratch);
        let mut r = Natural { limbs: r };
        r.normalize();
        r.debug_check();
        r
    }
}

impl Mul<&Natural> for &Natural {
    type Output = Natural;

    fn mul(self, rhs: &Natural) -> Natural {
        if std::ptr::eq(self, rhs) {
            return self.square();
        }
        let (a, b) = (self.limbs(), rhs.limbs());
        if a.is_empty() || b.is_empty() {
            return Natural::ZERO;
        }
        let mut r = LimbVec::zeroed(a.len() + b.len());
        let mut scratch = vec![0; kernel::mul_scratch_len(a.len(), b.len())];
        kernel::mul(r.as_mut_slice(), a, b, &mut scratch);
        let mut r = Natural { limbs: r };
        r.normalize();
        r.debug_check();
        r
    }
}

impl MulAssign<&Natural> for Natural {
    /// The product needs a fresh buffer, so this is `*self = &*self * rhs`.
    fn mul_assign(&mut self, rhs: &Natural) {
        *self = &*self * rhs;
    }
}

commutative_ops! {
    Mul::mul, MulAssign::mul_assign;
}

impl MulAssign<Limb> for Natural {
    fn mul_assign(&mut self, rhs: Limb) {
        if rhs == 0 {
            self.limbs.truncate(0);
            return;
        }
        let carry = kernel::mul_1_assign(self.limbs.as_mut_slice(), rhs);
        if carry != 0 {
            self.limbs.push(carry);
        }
        self.debug_check();
    }
}

impl Mul<Limb> for Natural {
    type Output = Natural;

    fn mul(mut self, rhs: Limb) -> Natural {
        self *= rhs;
        self
    }
}

impl Mul<Limb> for &Natural {
    type Output = Natural;

    fn mul(self, rhs: Limb) -> Natural {
        if rhs == 0 {
            return Natural::ZERO;
        }
        let a = self.limbs();
        let mut r = LimbVec::zeroed(a.len());
        let carry = kernel::mul_1(r.as_mut_slice(), a, rhs);
        if carry != 0 {
            r.push(carry);
        }
        Natural::from_limb_vec(r)
    }
}
