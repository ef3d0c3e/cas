//! The fused passes against the same arithmetic on GMP integers, modulo
//! `B^w` where the pass works modulo `B^w`. Lengths reach past a few blocks,
//! so each pass runs both its unrolled blocks and its single-limb tail.

use proptest::prelude::*;
use rug::Integer;

use super::{
    BLOCK, addlsh_assign, addmul_1_rsh1, lsh_add, rsh_add_assign, rsh_rsub_assign,
    rsh_sub_divexact_by3, rsub_divexact_by15, sublsh_assign, sublsh3_divexact_by9,
};
use crate::Limb;
use crate::testkit::strategy;

/// Longest operand: several blocks plus a tail.
const MAX: usize = 5 * BLOCK + 3;

fn val(a: &[Limb]) -> Integer {
    Integer::from_digits(a, rug::integer::Order::Lsf)
}

fn base_pow(w: usize) -> Integer {
    Integer::from(1) << (64 * u32::try_from(w).unwrap())
}

/// `x mod B^w`, nonnegative.
fn wrap(x: Integer, w: usize) -> Integer {
    x.div_rem_euc(base_pow(w)).1
}

/// The `w` limbs of `x mod B^w`.
fn limbs(x: Integer, w: usize) -> Vec<Limb> {
    let mut v = wrap(x, w).to_digits::<Limb>(rug::integer::Order::Lsf);
    v.resize(w, 0);
    v
}

/// A length in `0..=max`.
fn len(max: usize) -> BoxedStrategy<usize> {
    if max == 0 {
        return Just(0).boxed();
    }
    prop_oneof![1 => Just(0), 8 => strategy::len_up_to(max)].boxed()
}

/// Two raw values of `w` limbs.
fn pair(w: usize) -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    (strategy::raw_limbs(w), strategy::raw_limbs(w))
}

/// `x` of `w` limbs and `y` of at most `w`.
fn longer_first() -> impl Strategy<Value = (Vec<Limb>, Vec<Limb>)> {
    len(MAX).prop_flat_map(|w| {
        (
            strategy::raw_limbs(w),
            len(w).prop_flat_map(strategy::raw_limbs),
        )
    })
}

proptest! {
    #[test]
    fn addlsh_assign_matches((mut x, mut y) in longer_first(), k in 1..63u32) {
        // Room for the sum: x < B^w / 2 and y << k < B^w / 2.
        if let Some(xt) = x.last_mut() {
            *xt >>= 1;
        }
        if y.len() == x.len() {
            if let Some(yt) = y.last_mut() {
                *yt >>= k + 1;
            }
        }
        let exact = val(&x) + (val(&y) << k);
        addlsh_assign(&mut x, &y, k);
        prop_assert_eq!(val(&x), exact);
    }

    #[test]
    fn lsh_add_matches(
        (a, b) in len(MAX).prop_flat_map(|n| (strategy::raw_limbs(n), len(n).prop_flat_map(strategy::raw_limbs))),
        swap in any::<bool>(),
        j in 0..64u32,
        k in 0..64u32,
    ) {
        let (a, b) = if swap { (b, a) } else { (a, b) };
        let exact = (val(&a) << j) + (val(&b) << k);
        let mut n = a.len().max(b.len()) + 1;
        while exact >= base_pow(n) {
            n += 1;
        }
        let mut r = vec![0xdead_beef; n];
        lsh_add(&mut r, &a, j, &b, k);
        prop_assert_eq!(val(&r), exact);
    }

    #[test]
    fn sublsh_assign_matches((x, y) in longer_first(), k in 1..64u32) {
        let w = x.len();
        let mut r = x.clone();
        sublsh_assign(&mut r, &y, k);
        prop_assert_eq!(r, limbs(val(&x) - (val(&y) << k), w));
    }

    #[test]
    fn rsh_add_and_rsub_match((x, y) in len(MAX).prop_flat_map(pair), j in 1..64u32) {
        let w = x.len();
        let mut s = x.clone();
        rsh_add_assign(&mut s, &y, j);
        prop_assert_eq!(val(&s), wrap(val(&y) + val(&x), w) >> j);
        let mut d = x.clone();
        rsh_rsub_assign(&mut d, &y, j);
        prop_assert_eq!(val(&d), wrap(val(&y) - val(&x), w) >> j);
    }

    #[test]
    fn addmul_1_rsh1_matches((x, y) in len(MAX).prop_flat_map(pair), m in strategy::limb()) {
        let w = x.len();
        let mut r = x.clone();
        addmul_1_rsh1(&mut r, &y, m);
        prop_assert_eq!(val(&r), wrap(val(&x) + val(&y) * m, w) >> 1);
    }

    /// Each division gets a difference that is an exact multiple modulo
    /// `B^w`, of a quotient that may be "negative" (any `w` limbs).
    #[test]
    fn divexact_passes_match((q, y) in len(MAX).prop_flat_map(pair), j in 1..64u32) {
        let w = q.len();
        let (q, yv) = (val(&q), val(&y));

        // (x >> j) - y = 3 q, with x below B^w: x = (y + 3 q mod B^w) << j
        // plus low bits that the shift drops, when the top j bits are clear.
        let t = wrap(yv.clone() + Integer::from(3) * &q, w);
        if t < base_pow(w) >> j {
            let low = Integer::from(0x5a5a_5a5a_5a5a_5a5a_u64) % (Integer::from(1) << j);
            let mut x = limbs((t << j) + low, w);
            rsh_sub_divexact_by3(&mut x, &y, j);
            prop_assert_eq!(val(&x), q.clone());
        }

        let mut x = limbs(yv.clone() - Integer::from(15) * &q, w);
        rsub_divexact_by15(&mut x, &y);
        prop_assert_eq!(val(&x), q.clone());

        let mut x = limbs((yv << 3u32) + Integer::from(9) * &q, w);
        sublsh3_divexact_by9(&mut x, &y);
        prop_assert_eq!(val(&x), q);
    }
}

/// The shifted division where its input fills every bit: `x >> j` takes the
/// top limb's high bits from zeros above `x`.
#[test]
fn rsh_sub_divexact_by3_full_width() {
    for w in [1, BLOCK - 1, BLOCK, BLOCK + 1, 3 * BLOCK + 2] {
        for j in [1, 2, 63] {
            let x = vec![Limb::MAX; w];
            // (x >> j) - y = 3 q with y = (x >> j) mod 3, so q = (x >> j) div 3.
            let shifted = val(&x) >> j;
            let r = Integer::from(&shifted % 3u32);
            let y = limbs(r, w);
            let mut got = x.clone();
            rsh_sub_divexact_by3(&mut got, &y, j);
            assert_eq!(val(&got), shifted / 3u32, "w = {w}, j = {j}");
        }
    }
}
