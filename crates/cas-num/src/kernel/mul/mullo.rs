//! Short products: the low half of a product, `a b mod B^n` (GMP's
//! `mpn_mullo_n`), for the quotient of Montgomery's REDC.
//!
//! Schoolbook skips the partial products above `B^n`, half of them. Divide
//! and conquer (Mulders' split, with GMP's proportions) takes the full product
//! of the low `n2` limbs and two short products of `n1 = n - n2` limbs for the
//! cross terms: with `n2` about 0.7n, cheaper than a full product in the
//! Karatsuba and Toom range. In the NTT range a full product is as cheap, as
//! its cost hardly depends on the half dropped.

use super::{Algorithm, addmul_1, addmul_2_carries, algorithm, mul, mul_1, mul_scratch_len};
use crate::Limb;
use crate::kernel::add_n_assign;

/// Smallest length at which [`mullo`] splits by divide and conquer instead of
/// schoolbook. From the tuning configuration.
pub(crate) const MULLO_DC_THRESHOLD: usize = crate::tuning::MULLO_DC_THRESHOLD;
/// Smallest length at which [`mullo`] takes the full product.
pub(crate) const MULLO_MUL_THRESHOLD: usize = crate::tuning::MULLO_MUL_THRESHOLD;

/// The algorithm [`mullo`] uses for `n` limbs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mullo {
    Basecase,
    Dc,
    Full,
}

pub(crate) fn mullo_algorithm(n: usize) -> Mullo {
    if n >= MULLO_MUL_THRESHOLD {
        Mullo::Full
    } else if n >= MULLO_DC_THRESHOLD {
        Mullo::Dc
    } else {
        Mullo::Basecase
    }
}

/// `n1`, the length of the cross terms' short products, for `n >= 2`: GMP's
/// proportions for the algorithm of the full product, in `1..=n / 2`.
fn dc_split(n: usize) -> usize {
    let n1 = match algorithm(n, n) {
        Algorithm::Basecase | Algorithm::Chunked => n / 2,
        Algorithm::Karatsuba => n * 11 / 36,
        Algorithm::Toom3 => n * 9 / 40,
        Algorithm::Toom4 | Algorithm::Ntt => n * 7 / 39,
    };
    n1.clamp(1, n / 2)
}

/// Scratch limbs [`mullo`] needs for `n` limbs.
#[must_use]
pub fn mullo_scratch_len(n: usize) -> usize {
    match mullo_algorithm(n) {
        Mullo::Basecase => 0,
        Mullo::Full => 2 * n + mul_scratch_len(n, n),
        Mullo::Dc => {
            let n1 = dc_split(n);
            let n2 = n - n1;
            (2 * n2 + mul_scratch_len(n2, n2)).max(n1 + mullo_scratch_len(n1))
        }
    }
}

/// `r = a b mod B^n` for `n = r.len() == a.len() == b.len()` and
/// `scratch.len() >= mullo_scratch_len(n)`.
pub fn mullo(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    assert!(a.len() == n && b.len() == n, "mullo: length mismatch");
    assert!(
        scratch.len() >= mullo_scratch_len(n),
        "mullo: scratch too short"
    );
    mullo_unchecked(r, a, b, scratch);
}

fn mullo_unchecked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    match mullo_algorithm(n) {
        Mullo::Basecase => mullo_basecase(r, a, b),
        Mullo::Full => {
            let (t, rest) = scratch.split_at_mut(2 * n);
            mul(t, a, b, rest);
            r.copy_from_slice(&t[..n]);
        }
        Mullo::Dc => {
            let n1 = dc_split(n);
            let n2 = n - n1;
            {
                let (t, rest) = scratch.split_at_mut(2 * n2);
                mul(t, &a[..n2], &b[..n2], rest);
                r.copy_from_slice(&t[..n]);
            }
            // The cross terms a1 b0 + a0 b1 at B^n2, each needed modulo B^n1.
            let (c, rest) = scratch.split_at_mut(n1);
            mullo_unchecked(c, &a[n2..], &b[..n1], rest);
            add_n_assign(&mut r[n2..], c);
            mullo_unchecked(c, &a[..n1], &b[n2..], rest);
            add_n_assign(&mut r[n2..], c);
        }
    }
}

/// [`mullo`] by schoolbook: row `j` adds `a[..n - j] b[j]` at `B^j`, two rows
/// at a time; the carries out of `B^n` are dropped.
pub(crate) fn mullo_basecase(r: &mut [Limb], a: &[Limb], b: &[Limb]) {
    let n = r.len();
    debug_assert!(a.len() == n && b.len() == n);
    if n == 0 {
        return;
    }
    mul_1(r, a, b[0]);
    let mut j = 1;
    while j + 1 < n {
        // Row j + 1 also gets a[n - j - 1] b[j + 1], at B^n: dropped.
        addmul_2_carries(&mut r[j..], &a[..n - j], b[j], b[j + 1], 0);
        j += 2;
    }
    if j < n {
        addmul_1(&mut r[j..], &a[..n - j], b[j]);
    }
}
