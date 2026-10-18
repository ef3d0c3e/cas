//! Short products, high half: a lower bound of `floor(a b / B^n)` with a
//! proven error, for ball midpoints (CAS-DOC-5 §5.2), which need a product
//! near the true one, not the exact top limbs.
//!
//! Schoolbook computes only the columns from `B^(n - 1)` up, about half the
//! partial products; up to 2 limbs the full product is taken instead, as it
//! is faster there. Divide and conquer (Mulders' split, mirrored from
//! [`mullo`](super::mullo)) takes the full product of the high `n2` limbs and
//! two short products of `n1 = n - n2` limbs for the cross terms. In the NTT
//! range the full product is taken, which is exact. The thresholds are
//! [`mullo`](super::mullo)'s: both products drop the same half.
//!
//! Every algorithm returns `r` with `r <= floor(a b / B^n) <= r + e`, where
//! `e = mulhi_error(n)` follows the same recursion as the computation.

use super::mullo::{MULLO_DC_THRESHOLD, MULLO_MUL_THRESHOLD};
use super::{Algorithm, addmul_1, addmul_2_carries, algorithm, mul, mul_1, mul_scratch_len};
use crate::Limb;
use crate::kernel::add_assign;

/// The algorithm [`mulhi`] uses for `n` limbs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mulhi {
    Basecase,
    Dc,
    Full,
}

/// Up to this length the full product is faster than schoolbook's half
/// (measured: 6.5 against 7.1 ns at 2 limbs), and exact.
const MULHI_FULL_MAX_SMALL: usize = 2;

pub(crate) fn mulhi_algorithm(n: usize) -> Mulhi {
    if n >= MULLO_MUL_THRESHOLD || n <= MULHI_FULL_MAX_SMALL {
        Mulhi::Full
    } else if n >= MULLO_DC_THRESHOLD {
        Mulhi::Dc
    } else {
        Mulhi::Basecase
    }
}

/// `n1`, the length of the cross terms' short products, for `n >= 2`, in
/// `1..=n / 2`: the proportions of `mullo`'s split.
fn dc_split(n: usize) -> usize {
    let n1 = match algorithm(n, n) {
        // Toom-32 and Toom-42 take unbalanced shapes only: never chosen here.
        Algorithm::Basecase | Algorithm::Chunked | Algorithm::Toom32 | Algorithm::Toom42 => n / 2,
        Algorithm::Karatsuba => n * 11 / 36,
        Algorithm::Toom3 => n * 9 / 40,
        Algorithm::Toom4 | Algorithm::Ntt => n * 7 / 39,
    };
    n1.clamp(1, n / 2)
}

/// The largest `floor(a b / B^n) - r` over all `n`-limb operands, for the `r`
/// of [`mulhi`]:
///
/// - **Schoolbook**: `n - 1`. The dropped columns `k <= n - 2` hold at most
///   `k + 1` partial products, each below `B^2` at `B^k`, so their sum `D` is
///   below `(B - 1)^2 (n - 1) (B^(n - 1) - 1) / (B - 1) < (n - 1) B^n`. The
///   columns from `n - 1` up are exact (column `n - 2` partly kept, which only
///   lowers the error), and `r` drops the two limbs below `B^n`, so
///   `a b / B^n - r < 1 + (n - 1)`. In particular one limb is exact.
/// - **Divide and conquer**, with `a = a1 B^n1 + a0` (`a1` of `n2` limbs) and
///   likewise `b`: `r` is the top `n` limbs of `a1 b1`, plus `mulhi(a1h, b0)`
///   and `mulhi(a0, b1h)` where `a1h`, `b1h` are the top `n1` limbs of `a1`,
///   `b1`. Against `a b / B^n = a1 b1 / B^(n2 - n1) + (a1 b0 + a0 b1) / B^n2 +
///   a0 b0 / B^n`, the first term loses below 1 to truncation; each cross term
///   loses `a1l b0 / B^n2 < 1` to `a1h`, then below `e(n1) + 1` to the
///   recursive `mulhi`; the last term is below `B^(2 n1) / B^n <= 1`. In all,
///   `a b / B^n - r < 2 e(n1) + 6`, so `e(n) = 2 e(n1) + 5`.
/// - **Full product** (up to 2 limbs, and from `MULLO_MUL_THRESHOLD`): 0.
///
/// `r` never exceeds the true value, so the sums inside never overflow `n`
/// limbs.
#[must_use]
pub fn mulhi_error(n: usize) -> usize {
    match mulhi_algorithm(n) {
        Mulhi::Basecase => n.saturating_sub(1),
        Mulhi::Dc => 2 * mulhi_error(dc_split(n)) + 5,
        Mulhi::Full => 0,
    }
}

/// Scratch limbs [`mulhi`] needs for `n` limbs.
#[must_use]
pub fn mulhi_scratch_len(n: usize) -> usize {
    match mulhi_algorithm(n) {
        Mulhi::Basecase => n + 2,
        Mulhi::Full => 2 * n + mul_scratch_len(n, n),
        Mulhi::Dc => {
            let n1 = dc_split(n);
            let n2 = n - n1;
            (2 * n2 + mul_scratch_len(n2, n2)).max(n1 + mulhi_scratch_len(n1))
        }
    }
}

/// `r` with `r <= floor(a b / B^n) <= r + mulhi_error(n)`, for
/// `n = r.len() == a.len() == b.len()` and
/// `scratch.len() >= mulhi_scratch_len(n)`.
pub fn mulhi(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    assert!(a.len() == n && b.len() == n, "mulhi: length mismatch");
    assert!(
        scratch.len() >= mulhi_scratch_len(n),
        "mulhi: scratch too short"
    );
    mulhi_unchecked(r, a, b, scratch);
}

fn mulhi_unchecked(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    match mulhi_algorithm(n) {
        Mulhi::Basecase => mulhi_basecase(r, a, b, scratch),
        Mulhi::Full => {
            let (t, rest) = scratch.split_at_mut(2 * n);
            mul(t, a, b, rest);
            r.copy_from_slice(&t[n..]);
        }
        Mulhi::Dc => {
            let n1 = dc_split(n);
            let n2 = n - n1;
            {
                let (t, rest) = scratch.split_at_mut(2 * n2);
                mul(t, &a[n1..], &b[n1..], rest);
                r.copy_from_slice(&t[n2 - n1..]);
            }
            // The cross terms a1 b0 + a0 b1, each by the top n1 limbs of its
            // high factor, at B^0 of r.
            let (c, rest) = scratch.split_at_mut(n1);
            mulhi_unchecked(c, &a[n - n1..], &b[..n1], rest);
            let carry = add_assign(r, c);
            debug_assert_eq!(carry, 0, "mulhi: above the true product");
            mulhi_unchecked(c, &a[..n1], &b[n - n1..], rest);
            let carry = add_assign(r, c);
            debug_assert_eq!(carry, 0, "mulhi: above the true product");
        }
    }
}

/// [`mulhi`] by schoolbook, in `t` (`n + 2` limbs, the columns from
/// `B^(n - 2)` up), two rows at a time as in `mul_basecase`.
///
/// Row 0 goes to `B^(n - 1)`. Rows `j` and `j + 1` then take
/// `a[n - 2 - j..] (b[j] + b[j + 1] B)` at `B^(n - 2)`: all their partial
/// products from `B^(n - 1)` up, plus row `j`'s one in column `n - 2`, which
/// only adds part of a dropped column (still a lower bound, within the same
/// error). An odd row left over goes by `addmul_1` at `B^(n - 1)`. Each pass
/// sets the limbs above it from its carries.
pub(crate) fn mulhi_basecase(r: &mut [Limb], a: &[Limb], b: &[Limb], t: &mut [Limb]) {
    let n = r.len();
    debug_assert!(a.len() == n && b.len() == n && t.len() >= n + 2);
    if n == 0 {
        return;
    }
    t[0] = 0;
    t[2] = mul_1(&mut t[1..2], &a[n - 1..], b[0]);
    let mut j = 1;
    while j + 1 < n {
        let (c0, c1) = addmul_2_carries(&mut t[..j + 2], &a[n - 2 - j..], b[j], b[j + 1], 0);
        t[j + 2] = c0;
        t[j + 3] = c1;
        j += 2;
    }
    if j < n {
        t[j + 2] = addmul_1(&mut t[1..j + 2], &a[n - 1 - j..], b[j]);
    }
    r.copy_from_slice(&t[2..n + 2]);
}

#[cfg(test)]
mod tests;
