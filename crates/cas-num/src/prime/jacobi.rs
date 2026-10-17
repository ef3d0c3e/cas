//! The Jacobi symbol.

use crate::{Integer, Natural};

#[cfg(test)]
mod tests;

impl Integer {
    /// The Jacobi symbol `(self / n)`: -1, 0 or 1. For a prime `n` it is the
    /// Legendre symbol, 1 for a non-zero square modulo `n` and -1 for a
    /// non-square.
    ///
    /// # Panics
    ///
    /// If `n` is not odd and positive.
    #[must_use]
    pub fn jacobi(&self, n: &Integer) -> i32 {
        assert!(
            n.is_positive() && !n.is_divisible_by(&Integer::from(2)),
            "Jacobi symbol needs an odd positive modulus"
        );
        let n = n.unsigned_abs();
        let a = self.div_rem_euclid(&Integer::from(&n)).1.unsigned_abs();
        jacobi(a, n)
    }
}

/// `(a / n)` for odd `n` and `a < n`, by the binary algorithm: factors of
/// two are removed with the second supplement, then reciprocity swaps the
/// arguments (Cohen, Algorithm 1.4.10). Word arithmetic once both fit.
fn jacobi(mut a: Natural, mut n: Natural) -> i32 {
    let mut t = 1;
    loop {
        if let (Ok(a), Ok(n)) = (u64::try_from(&a), u64::try_from(&n)) {
            return t * jacobi_u64(a, n);
        }
        if a.is_zero() {
            return if n.is_one() { t } else { 0 };
        }
        let z = a.trailing_zeros().expect("non-zero");
        a >>= z;
        let n8 = n.limbs()[0] & 7;
        if z % 2 == 1 && (n8 == 3 || n8 == 5) {
            t = -t;
        }
        if a.limbs()[0] & 3 == 3 && n8 & 3 == 3 {
            t = -t;
        }
        (a, n) = (&n % &a, a);
    }
}

fn jacobi_u64(mut a: u64, mut n: u64) -> i32 {
    debug_assert!(n % 2 == 1 && a < n);
    let mut t = 1;
    while a != 0 {
        let z = a.trailing_zeros();
        a >>= z;
        if z % 2 == 1 && (n & 7 == 3 || n & 7 == 5) {
            t = -t;
        }
        if a & 3 == 3 && n & 3 == 3 {
            t = -t;
        }
        (a, n) = (n % a, a);
    }
    if n == 1 { t } else { 0 }
}
