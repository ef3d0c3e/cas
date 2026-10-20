//! Factorials.

use super::Natural;
use crate::Limb;

#[cfg(test)]
mod tests;

/// Up to this, `n!` fits in a limb: `20! < 2^64 < 21!`.
const LIMB_FACTORIAL_MAX: u64 = 20;

impl Natural {
    /// `n!`.
    ///
    /// The odd part as in GMP's method for small `n` \[GMP\]: it is the
    /// product of the odd numbers up to `n` times the odd part of `⌊n/2⌋!`,
    /// and so on down, as in `23! = (23·21·…·3)(11·9·…·3)(5·3)·2^19`. The
    /// power of two, `2^(n - s(n))` with `s(n)` the ones of `n` (Legendre's
    /// formula), is one shift at the end. From the top level down, a running
    /// product takes in the odd numbers each level adds, those in
    /// `(⌊n/2^(i+1)⌋, ⌊n/2^i⌋]`, and multiplies into the odd part. Those
    /// numbers are multiplied in a balanced product tree whose leaves hold
    /// as many of them as fit in a limb. Up to 20, the result is one limb.
    #[must_use]
    pub fn factorial(n: u64) -> Natural {
        if n <= LIMB_FACTORIAL_MAX {
            return Natural::from((2..=n).product::<u64>());
        }
        let mut odd = Natural::one();
        // The product of the odd numbers up to `high`.
        let mut level = Natural::one();
        let mut high = 1;
        for i in (0..=n.ilog2()).rev() {
            let m = n >> i;
            // The largest odd number up to `m`. From `m = 3` on, it grows
            // at every level; below, the level adds nothing.
            let top = m - 1 + m % 2;
            if top > high {
                level *= &odd_product(high + 2, top);
                high = top;
                odd *= &level;
            }
        }
        odd << (n - u64::from(n.count_ones()))
    }
}

/// The product of the odd numbers from `a` to `b`, both odd, `a <= b`.
fn odd_product(a: u64, b: u64) -> Natural {
    let mut leaves: Vec<Limb> = Vec::new();
    let mut leaf: Limb = 1;
    let mut k = a;
    loop {
        if let Some(p) = leaf.checked_mul(k) {
            leaf = p;
        } else {
            leaves.push(leaf);
            leaf = k;
        }
        if k == b {
            break;
        }
        k += 2;
    }
    leaves.push(leaf);
    product(&leaves)
}

/// The product of `leaves`, each half multiplied apart, so that every
/// product has operands of about the same size.
fn product(leaves: &[Limb]) -> Natural {
    match leaves {
        [] => Natural::one(),
        [x] => Natural::from(*x),
        [x, y] => Natural::from(*x) * *y,
        _ => {
            let (low, high) = leaves.split_at(leaves.len() / 2);
            &product(low) * &product(high)
        }
    }
}
