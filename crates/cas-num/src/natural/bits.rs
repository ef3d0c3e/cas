//! Shifts, bitwise operations and bit queries.

use std::ops::{BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign};
use std::ops::{Shl, ShlAssign, Shr, ShrAssign};

use super::Natural;
use super::limb_vec::LimbVec;
use super::ops::commutative_ops;
use crate::{Limb, kernel};

/// Splits a bit index into a limb index and a bit within that limb.
///
/// # Panics
///
/// If the limb index does not fit `usize` (impossible on 64-bit targets).
pub(super) fn split_bit_index(k: u64) -> (usize, u32) {
    let bits = u64::from(Limb::BITS);
    let limb = usize::try_from(k / bits).expect("bit index too large");
    let bit = u32::try_from(k % bits).expect("bit within limb");
    (limb, bit)
}

impl Natural {
    /// Number of significant bits (0 for zero).
    #[must_use]
    pub fn bit_length(&self) -> u64 {
        kernel::bit_length(self.limbs())
    }

    /// Number of trailing zero bits, or `None` for zero.
    #[must_use]
    pub fn trailing_zeros(&self) -> Option<u64> {
        kernel::trailing_zeros(self.limbs())
    }

    /// Number of one bits.
    #[must_use]
    pub fn count_ones(&self) -> u64 {
        kernel::count_ones(self.limbs())
    }

    /// Whether bit `k` is set.
    #[must_use]
    pub fn bit(&self, k: u64) -> bool {
        let (limb, bit) = split_bit_index(k);
        self.limbs().get(limb).is_some_and(|x| (x >> bit) & 1 == 1)
    }

    /// Whether the value is `2^k` for some `k`.
    #[must_use]
    pub fn is_power_of_two(&self) -> bool {
        match self.limbs().split_last() {
            None => false,
            Some((top, rest)) => top.is_power_of_two() && rest.iter().all(|&x| x == 0),
        }
    }
}

impl Shl<u64> for &Natural {
    type Output = Natural;

    /// # Panics
    ///
    /// If the result would not fit in memory.
    fn shl(self, k: u64) -> Natural {
        if self.is_zero() {
            return Natural::ZERO;
        }
        let (limbs, bit) = split_bit_index(k);
        let a = self.limbs();
        let len = a.len().checked_add(limbs).expect("shift too large");
        let mut r = LimbVec::zeroed(len);
        let carry = kernel::shl(&mut r.as_mut_slice()[limbs..], a, bit);
        if carry != 0 {
            r.push(carry);
        }
        Natural::from_limb_vec(r)
    }
}

impl ShlAssign<u64> for Natural {
    /// # Panics
    ///
    /// If the result would not fit in memory.
    fn shl_assign(&mut self, k: u64) {
        if self.is_zero() {
            return;
        }
        let (limbs, bit) = split_bit_index(k);
        let n = self.limb_len();
        self.limbs
            .resize(n.checked_add(limbs).expect("shift too large"));
        let r = self.limbs.as_mut_slice();
        r.copy_within(..n, limbs);
        r[..limbs].fill(0);
        let carry = kernel::shl_assign(&mut r[limbs..], bit);
        if carry != 0 {
            self.limbs.push(carry);
        }
        self.debug_check();
    }
}

impl Shl<u64> for Natural {
    type Output = Natural;

    fn shl(mut self, k: u64) -> Natural {
        self <<= k;
        self
    }
}

impl Shr<u64> for &Natural {
    type Output = Natural;

    fn shr(self, k: u64) -> Natural {
        let (limbs, bit) = split_bit_index(k);
        let a = self.limbs();
        if limbs >= a.len() {
            return Natural::ZERO;
        }
        let mut r = Natural {
            limbs: LimbVec::zeroed(a.len() - limbs),
        };
        kernel::shr(r.limbs.as_mut_slice(), &a[limbs..], bit);
        r.normalize();
        r.debug_check();
        r
    }
}

impl ShrAssign<u64> for Natural {
    fn shr_assign(&mut self, k: u64) {
        let (limbs, bit) = split_bit_index(k);
        let n = self.limb_len();
        if limbs >= n {
            self.limbs.truncate(0);
            return;
        }
        self.limbs.as_mut_slice().copy_within(limbs.., 0);
        self.limbs.truncate(n - limbs);
        kernel::shr_assign(self.limbs.as_mut_slice(), bit);
        self.normalize();
        self.debug_check();
    }
}

impl Shr<u64> for Natural {
    type Output = Natural;

    fn shr(mut self, k: u64) -> Natural {
        self >>= k;
        self
    }
}

impl BitAnd<&Natural> for &Natural {
    type Output = Natural;

    fn bitand(self, rhs: &Natural) -> Natural {
        let n = self.limb_len().min(rhs.limb_len());
        let mut r = Natural {
            limbs: LimbVec::zeroed(n),
        };
        kernel::and_n(
            r.limbs.as_mut_slice(),
            &self.limbs()[..n],
            &rhs.limbs()[..n],
        );
        r.normalize();
        r.debug_check();
        r
    }
}

impl BitAndAssign<&Natural> for Natural {
    fn bitand_assign(&mut self, rhs: &Natural) {
        let n = self.limb_len().min(rhs.limb_len());
        self.limbs.truncate(n);
        kernel::and_n_assign(self.limbs.as_mut_slice(), &rhs.limbs()[..n]);
        self.normalize();
        self.debug_check();
    }
}

impl Natural {
    /// `self op= rhs` for or/xor: both keep the longer operand's high limbs.
    fn or_xor_assign(&mut self, rhs: &Natural, op: fn(&mut [Limb], &[Limb])) {
        let n = rhs.limb_len();
        if self.limb_len() < n {
            self.limbs.resize(n);
        }
        op(&mut self.limbs.as_mut_slice()[..n], rhs.limbs());
        // Only xor of equal-length operands can clear the top limb.
        self.normalize();
        self.debug_check();
    }

    fn or_xor(&self, rhs: &Natural, op: fn(&mut [Limb], &[Limb])) -> Natural {
        let (long, short) = if self.limb_len() >= rhs.limb_len() {
            (self, rhs)
        } else {
            (rhs, self)
        };
        let mut r = Natural {
            limbs: LimbVec::from_slice(long.limbs()),
        };
        op(
            &mut r.limbs.as_mut_slice()[..short.limb_len()],
            short.limbs(),
        );
        r.normalize();
        r.debug_check();
        r
    }
}

impl BitOr<&Natural> for &Natural {
    type Output = Natural;

    fn bitor(self, rhs: &Natural) -> Natural {
        self.or_xor(rhs, kernel::or_n_assign)
    }
}

impl BitOrAssign<&Natural> for Natural {
    fn bitor_assign(&mut self, rhs: &Natural) {
        self.or_xor_assign(rhs, kernel::or_n_assign);
    }
}

impl BitXor<&Natural> for &Natural {
    type Output = Natural;

    fn bitxor(self, rhs: &Natural) -> Natural {
        self.or_xor(rhs, kernel::xor_n_assign)
    }
}

impl BitXorAssign<&Natural> for Natural {
    fn bitxor_assign(&mut self, rhs: &Natural) {
        self.or_xor_assign(rhs, kernel::xor_n_assign);
    }
}

commutative_ops! {
    BitAnd::bitand, BitAndAssign::bitand_assign;
    BitOr::bitor, BitOrAssign::bitor_assign;
    BitXor::bitxor, BitXorAssign::bitxor_assign;
}
