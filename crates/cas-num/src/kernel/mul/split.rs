//! Products modulo `B^m - 1` below the cyclic transform, split as GMP's
//! `mpn_mulmod_bnm1` does. For even `m = 2h`, `B^m - 1 = (B^h - 1)(B^h + 1)`,
//! and the factors are coprime (the first is odd and their difference is 2).
//! The product is taken modulo each and recombined by the Chinese remainder
//! theorem. Modulo `B^h - 1` it recurses through [`mul_wrap`], which splits
//! again or takes a cyclic transform; modulo `B^h + 1` it is a full product
//! of `h` limbs, reduced. In the Toom range the halves cost about 0.6 of a
//! full product of `m` limbs, which is what reducing it would cost instead.

use super::{
    MUL_WRAP_SPLIT_THRESHOLD, fold, folded, mul, mul_scratch_len, mul_wrap, mul_wrap_scratch_len,
    wraps_cyclically,
};
use crate::Limb;
use crate::kernel::{add_1_assign, shr_assign, sub_1_assign, sub_assign, sub_n_assign};

/// Whether [`mul_wrap`] modulo `B^m - 1` splits it for operands of `an` and
/// `bn` limbs. That takes an even `m` from [`MUL_WRAP_SPLIT_THRESHOLD`] that
/// takes no cyclic transform, and operands long enough that the full
/// product would be well over `m` limbs: the two halves cost about as much
/// as a full product of `3m/4` limbs.
pub(crate) fn splits(m: usize, an: usize, bn: usize) -> bool {
    m % 2 == 0
        && m >= MUL_WRAP_SPLIT_THRESHOLD
        && !wraps_cyclically(m)
        && 2 * (an.min(m) + bn.min(m)) > 3 * m
}

/// The modulus length for a product of two `x`-limb operands known modulo
/// `B^m - 1` by splitting: `x` rounded up so that the length stays even
/// for up to three halvings while it is at least
/// [`MUL_WRAP_SPLIT_THRESHOLD`], as GMP's `mpn_mulmod_bnm1_next_size`.
/// Deeper levels would gain little. The result is below `2x` for `x >= 2`.
pub(crate) fn split_wrap_len(x: usize) -> usize {
    let mut levels = 0;
    while levels < 3 && x >> levels >= MUL_WRAP_SPLIT_THRESHOLD.max(2) {
        levels += 1;
    }
    x.next_multiple_of(1 << levels)
}

/// Scratch limbs [`mul_split`] needs modulo `B^m - 1` for operands of `an`
/// and `bn` limbs.
pub(crate) fn mul_split_scratch_len(m: usize, an: usize, bn: usize) -> usize {
    let h = m / 2;
    let folded = |len: usize| if len > m { m } else { 0 };
    let minus = mul_wrap_scratch_len(h, h, h);
    let plus = 2 * h + mul_scratch_len(h, h);
    folded(an) + folded(bn) + 3 * (h + 1) + minus.max(plus)
}

/// `r = a b mod (B^m - 1)` for even `m = r.len()` and non-empty operands, as
/// a residue in `[0, B^m - 1)`, with `scratch.len() >=
/// mul_split_scratch_len(m, a.len(), b.len())`.
///
/// With `x` the product, `xm = x mod (B^h - 1)` and `xp = x mod (B^h + 1)`,
/// `x = xm + (B^h - 1) y` for `y = (xm - xp) / 2 mod (B^h + 1)`, because
/// `B^h - 1 = -2` modulo `B^h + 1`. With `xm < B^h - 1` and `y <= B^h`,
/// that is below `B^m - 1`.
pub(crate) fn mul_split(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let m = r.len();
    debug_assert!(m % 2 == 0 && m >= 2);
    let h = m / 2;
    let (a, scratch) = folded(a, m, scratch);
    let (b, scratch) = folded(b, m, scratch);
    let (x, rest) = scratch.split_at_mut(3 * (h + 1));
    // xm into r[..h], from the operands modulo B^h - 1.
    {
        let (am, bm) = x[..2 * h].split_at_mut(h);
        fold(am, a);
        fold(bm, b);
        mul_wrap(&mut r[..h], am, bm, rest);
    }
    let (ap, x) = x.split_at_mut(h + 1);
    let (bp, xp) = x.split_at_mut(h + 1);
    fold_plus(ap, a);
    fold_plus(bp, b);
    mul_plus(xp, ap, bp, rest);
    // y = (xm - xp) / 2 mod (B^h + 1), in ap: first xm + (B^h + 1) - xp,
    // which is in [1, 2 B^h), less B^h + 1 if it is at least that.
    let y = ap;
    y[..h].copy_from_slice(&r[..h]);
    y[h] = 1;
    add_1_assign(y, 1);
    let borrow = sub_n_assign(y, xp);
    debug_assert_eq!(borrow, 0);
    if y[h] == 1 && y[..h].iter().any(|&l| l != 0) {
        y[h] = 0;
        sub_1_assign(&mut y[..h], 1);
    }
    // Halving modulo the odd B^h + 1: an odd y gets B^h + 1 first, staying
    // below 2 (B^h + 1).
    if y[0] & 1 == 1 {
        let carry = add_1_assign(&mut y[..h], 1);
        y[h] += 1 + carry;
    }
    shr_assign(y, 1);
    // x = xm + y B^h - y.
    r[h..].copy_from_slice(&y[..h]);
    let borrow = sub_assign(r, y);
    debug_assert_eq!(borrow, y[h]);
}

/// `r = x mod (B^h + 1)` for `h = r.len() - 1` and `x` of at most `2h`
/// limbs, in `[0, B^h]`.
fn fold_plus(r: &mut [Limb], x: &[Limb]) {
    let h = r.len() - 1;
    let (lo, hi) = x.split_at(x.len().min(h));
    r[..lo.len()].copy_from_slice(lo);
    r[lo.len()..].fill(0);
    // B^h = -1: subtract the high part, and on a borrow add B^h + 1 back,
    // which leaves the wrapped difference plus one.
    if sub_assign(&mut r[..h], hi) != 0 {
        r[h] = add_1_assign(&mut r[..h], 1);
    }
}

/// `r = a b mod (B^h + 1)` for `h + 1`-limb `a` and `b` in `[0, B^h]`, in
/// `[0, B^h]`, with `scratch` for the product of their low limbs.
fn mul_plus(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let h = r.len() - 1;
    if a[h] != 0 {
        // a = B^h = -1.
        neg_plus(r, b);
    } else if b[h] != 0 {
        neg_plus(r, a);
    } else {
        let (p, rest) = scratch.split_at_mut(2 * h);
        mul(p, &a[..h], &b[..h], rest);
        fold_plus(r, p);
    }
}

/// `r = -x mod (B^h + 1)` for `x` in `[0, B^h]`: `B^h + 1 - x`, or 0.
fn neg_plus(r: &mut [Limb], x: &[Limb]) {
    let h = r.len() - 1;
    if x[h] != 0 {
        // x = B^h, and B^h + 1 - x = 1.
        r.fill(0);
        r[0] = 1;
    } else if x.iter().all(|&l| l == 0) {
        r.fill(0);
    } else {
        // The complement is B^h - 1 - x; two more make B^h + 1 - x <= B^h.
        for (r, &x) in r.iter_mut().zip(&x[..h]) {
            *r = !x;
        }
        r[h] = add_1_assign(&mut r[..h], 2);
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;
    use rug::Integer;

    use super::{fold_plus, mul_plus, neg_plus};
    use crate::Limb;
    use crate::kernel::mul_scratch_len;
    use crate::testkit::{Mag, ToRug};

    fn val(x: &[Limb]) -> Integer {
        Mag::from_limbs(x.to_vec()).to_rug()
    }

    /// `x` as `h + 1` limbs, for `x <= B^h`.
    fn limbs(x: &Integer, h: usize) -> Vec<Limb> {
        let mut l = Mag::from_rug(x).limbs().to_vec();
        l.resize(h + 1, 0);
        l
    }

    /// Residues modulo `B^h + 1`, `h` from 1 to 6, at the edges: 0, 1,
    /// `B^h - 1` and `B^h` (that is, -1), or random.
    fn residue() -> impl Strategy<Value = (usize, Integer, Integer)> {
        (1usize..=6).prop_flat_map(|h| {
            let one =
                |h: usize| {
                    let p = Integer::from(1) << (64 * u32::try_from(h).unwrap());
                    prop_oneof![
                        Just(Integer::ZERO),
                        Just(Integer::from(1)),
                        Just(Integer::from(&p - 1u32)),
                        Just(p.clone()),
                        proptest::collection::vec(any::<u64>(), h)
                            .prop_map(move |l| Mag::from_limbs(l).to_rug() % (p.clone() + 1u32)),
                    ]
                };
            (Just(h), one(h), one(h))
        })
    }

    proptest! {
        /// The helpers modulo `B^h + 1` give residues in `[0, B^h]`.
        #[test]
        fn plus_helpers((h, x, y) in residue(), wide in proptest::collection::vec(any::<u64>(), 0..=12)) {
            let p = (Integer::from(1) << (64 * u32::try_from(h).unwrap())) + 1u32;
            let (a, b) = (limbs(&x, h), limbs(&y, h));
            let mut r = vec![0xdead_beef; h + 1];
            neg_plus(&mut r, &a);
            prop_assert!(val(&r) < p && val(&r) == Integer::from(-&x).modulo(&p), "neg");
            let mut scratch = vec![0xdead_beef; 2 * h + mul_scratch_len(h, h)];
            mul_plus(&mut r, &a, &b, &mut scratch);
            prop_assert!(val(&r) < p && val(&r) == Integer::from(&x * &y) % &p, "mul");
            let wide = &wide[..wide.len().min(2 * h)];
            fold_plus(&mut r, wide);
            prop_assert!(val(&r) < p && val(&r) == val(wide) % &p, "fold");
        }
    }
}
