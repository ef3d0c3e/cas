//! Proptest strategies producing numbers with adversarial shapes.
//!
//! Uniformly random limbs almost never hit the inputs that break arithmetic
//! (carries rippling through all-ones limbs, values just past a limb boundary,
//! operands of very different sizes), so the generators below weight those
//! shapes heavily. Every generated value is normalised.

use proptest::collection::vec;
use proptest::prelude::*;

use super::ToRug;
use super::model::{Mag, Signed};

/// Default size bound, in limbs, for the unparameterised strategies.
pub(crate) const DEFAULT_MAX_LIMBS: usize = 40;

/// Limb values that sit on carry/borrow and sign boundaries.
const SPECIAL_LIMBS: [u64; 11] = [
    0,
    1,
    2,
    u64::MAX,
    u64::MAX - 1,
    1 << 63,
    (1 << 63) - 1,
    (1 << 63) + 1,
    1 << 32,
    (1 << 32) - 1,
    0x5555_5555_5555_5555,
];

/// A single limb: special values, single-bit and low-mask values, or uniform.
pub(crate) fn limb() -> impl Strategy<Value = u64> + Clone {
    prop_oneof![
        3 => proptest::sample::select(&SPECIAL_LIMBS[..]),
        1 => (0..64u32).prop_map(|k| 1u64 << k),
        1 => (0..64u32).prop_map(|k| (1u64 << k) - 1),
        4 => any::<u64>(),
    ]
}

/// A non-zero limb, for the most significant position.
pub(crate) fn top_limb() -> impl Strategy<Value = u64> + Clone {
    limb().prop_map(|x| if x == 0 { 1 } else { x })
}

/// A length in `1..=max`, weighted towards 1–4 limbs, where most special
/// cases (single limb, first carry into a new limb, inline storage) live.
pub(crate) fn len_up_to(max: usize) -> BoxedStrategy<usize> {
    assert!(max >= 1);
    let small = 1..=max.min(4);
    if max <= 4 {
        small.boxed()
    } else {
        prop_oneof![3 => small, 1 => 5..=max].boxed()
    }
}

/// Lengths `t - 1`, `t` and `t + 1` for each algorithm threshold `t`
/// (never below 1, and clamped by [`super::tested`]), so both sides of every
/// crossover get exercised.
pub(crate) fn len_around(thresholds: &[usize]) -> BoxedStrategy<usize> {
    let mut lens: Vec<usize> = thresholds
        .iter()
        .map(|&t| super::tested(t))
        .flat_map(|t| [t.saturating_sub(1), t, t + 1])
        .filter(|&n| n >= 1)
        .collect();
    lens.sort_unstable();
    lens.dedup();
    proptest::sample::select(lens).boxed()
}

/// `2^k`.
fn pow2(k: usize) -> Mag {
    let mut limbs = vec![0; k / 64 + 1];
    limbs[k / 64] = 1 << (k % 64);
    Mag::from_limbs(limbs)
}

/// `2^k - 1`: `k` one bits.
fn ones(k: usize) -> Mag {
    let mut limbs = vec![u64::MAX; k / 64];
    if k % 64 != 0 {
        limbs.push((1 << (k % 64)) - 1);
    }
    Mag::from_limbs(limbs)
}

/// `2^k + 1` (`k >= 1`).
fn pow2_plus_1(k: usize) -> Mag {
    let mut limbs = pow2(k).limbs().to_vec();
    limbs[0] |= 1;
    Mag::from_limbs(limbs)
}

/// A natural number with exactly `len` limbs (`len >= 1`), filled densely,
/// sparsely (mostly zero limbs), or with all-ones limbs.
pub(crate) fn natural_with_len(len: impl Strategy<Value = usize>) -> impl Strategy<Value = Mag> {
    len.prop_flat_map(|len| {
        assert!(len >= 1);
        let low = len - 1;
        let dense = vec(limb(), low).boxed();
        let sparse = vec(prop_oneof![6 => Just(0u64), 1 => limb()], low).boxed();
        let all_ones = Just(vec![u64::MAX; low]).boxed();
        let fill = prop_oneof![4 => dense, 2 => sparse, 1 => all_ones];
        (fill, top_limb()).prop_map(|(mut limbs, top)| {
            limbs.push(top);
            Mag::from_limbs(limbs)
        })
    })
}

/// A natural number of at most `max_limbs` limbs, mixing zero, one, single
/// limbs, powers of two, `2^k ± 1` around limb boundaries, all-ones values
/// and dense/sparse multi-limb values.
pub(crate) fn natural_up_to(max_limbs: usize) -> BoxedStrategy<Mag> {
    assert!(max_limbs >= 1);
    let max_bits = 64 * max_limbs;
    // Bit counts within one bit of a limb boundary: 64j - 1, 64j, 64j + 1.
    let boundary_bits =
        (1..=max_limbs, 0..3usize).prop_map(move |(j, o)| (64 * j + o - 1).min(max_bits - 1));
    prop_oneof![
        1 => Just(Mag::default()),
        1 => Just(Mag::from_u64(1)),
        3 => top_limb().prop_map(Mag::from_u64),
        1 => (0..max_bits).prop_map(pow2),
        2 => boundary_bits.clone().prop_map(pow2_plus_1),
        2 => boundary_bits.prop_map(ones),
        1 => (1..=max_bits).prop_map(ones),
        6 => natural_with_len(len_up_to(max_limbs)),
    ]
    .boxed()
}

/// [`natural_up_to`] with [`DEFAULT_MAX_LIMBS`].
pub(crate) fn natural() -> BoxedStrategy<Mag> {
    natural_up_to(DEFAULT_MAX_LIMBS)
}

/// A signed integer of at most `max_limbs` limbs.
pub(crate) fn signed_up_to(max_limbs: usize) -> BoxedStrategy<Signed> {
    (any::<bool>(), natural_up_to(max_limbs))
        .prop_map(|(neg, mag)| Signed::new(neg, mag))
        .boxed()
}

/// [`signed_up_to`] with [`DEFAULT_MAX_LIMBS`].
pub(crate) fn signed() -> BoxedStrategy<Signed> {
    signed_up_to(DEFAULT_MAX_LIMBS)
}

/// A pair of naturals of very different sizes: one of `max_limbs / 2 ..= max_limbs`
/// limbs, one of 1–3 limbs, in either order.
pub(crate) fn unbalanced_pair(max_limbs: usize) -> BoxedStrategy<(Mag, Mag)> {
    assert!(max_limbs >= 8);
    let long = natural_with_len(max_limbs / 2..=max_limbs);
    let short = natural_with_len(1..=3usize);
    (long, short, any::<bool>())
        .prop_map(|(long, short, swap)| if swap { (short, long) } else { (long, short) })
        .boxed()
}

/// Raw limbs of exactly `len` limbs for kernel tests: unlike the `natural`
/// strategies, high limbs may be zero (kernels must not assume normalisation).
pub(crate) fn raw_limbs(len: usize) -> BoxedStrategy<Vec<u64>> {
    prop_oneof![
        4 => vec(limb(), len),
        2 => vec(prop_oneof![6 => Just(0u64), 1 => limb()], len),
        1 => Just(vec![u64::MAX; len]),
        1 => Just(vec![0; len]),
    ]
    .boxed()
}

/// Two raw limb vectors of the same length in `0..=max_limbs`.
pub(crate) fn raw_pair_same_len(max_limbs: usize) -> BoxedStrategy<(Vec<u64>, Vec<u64>)> {
    prop_oneof![1 => Just(0), 8 => len_up_to(max_limbs)]
        .prop_flat_map(|n| (raw_limbs(n), raw_limbs(n)))
        .boxed()
}

/// Raw limb vectors `(a, b)` with `b.len() <= a.len() <= max_limbs`.
pub(crate) fn raw_pair_longer_first(max_limbs: usize) -> BoxedStrategy<(Vec<u64>, Vec<u64>)> {
    prop_oneof![1 => Just(0), 8 => len_up_to(max_limbs)]
        .prop_flat_map(|n| (Just(n), 0..=n))
        .prop_flat_map(|(n, m)| (raw_limbs(n), raw_limbs(m)))
        .boxed()
}

/// Signed values within a few units of a promotion edge of the small-integer
/// fast path: `0`, `±2^31`, `±2^32`, `±2^62`, `±2^63` (the `i64` limits),
/// `±2^64`, `±2^65`, `±2^127` and `±2^128` (the `i128` limits).
pub(crate) fn signed_near_i64_bounds() -> BoxedStrategy<Signed> {
    const EXPONENTS: [u32; 9] = [0, 31, 32, 62, 63, 64, 65, 127, 128];
    (
        proptest::sample::select(&EXPONENTS[..]),
        -3i32..=3,
        any::<bool>(),
    )
        .prop_map(|(k, offset, negative)| {
            let base = if k == 0 {
                rug::Integer::new()
            } else {
                rug::Integer::from(1) << k
            };
            let x = base + offset;
            Signed::from_rug(&if negative { -x } else { x })
        })
        .boxed()
}

/// A signed integer for `Integer` tests: arbitrary sizes, promotion edges, or
/// any `i64`.
pub(crate) fn integer() -> BoxedStrategy<Signed> {
    prop_oneof![
        3 => signed(),
        3 => signed_near_i64_bounds(),
        2 => any::<i64>().prop_map(Signed::from_i64),
    ]
    .boxed()
}

/// Pairs of naturals for GCD tests, in either order: independent values
/// (mostly coprime), values sharing a large common factor, consecutive
/// Fibonacci numbers (Euclid's worst case: every quotient is 1), a value and
/// a multiple of it, equal values, and values sharing their leading limbs
/// (where Lehmer's leading-bit estimates are least informative).
pub(crate) fn gcd_pair(max_limbs: usize) -> BoxedStrategy<(Mag, Mag)> {
    assert!(max_limbs >= 3);
    let third = max_limbs / 3;
    let product = |x: &Mag, y: &Mag| Mag::from_rug(&(x.to_rug() * y.to_rug()));
    let fib_max = u32::try_from(64 * max_limbs * 36 / 25).expect("index fits u32");
    let pairs = prop_oneof![
        3 => (natural_up_to(max_limbs), natural_up_to(max_limbs)),
        3 => (natural_up_to(third), natural_up_to(third), natural_up_to(third))
            .prop_map(move |(a, b, c)| (product(&a, &c), product(&b, &c))),
        2 => (0..fib_max).prop_map(|n| {
            let f = |k| Mag::from_rug(&rug::Integer::from(rug::Integer::fibonacci(k)));
            (f(n + 1), f(n))
        }),
        1 => (natural_up_to(2 * third), natural_up_to(third))
            .prop_map(move |(a, k)| (product(&a, &k), a)),
        1 => natural_up_to(max_limbs).prop_map(|a| (a.clone(), a)),
        2 => (2..=max_limbs)
            .prop_flat_map(|n| {
                let m = || natural_with_len(Just(n));
                (m(), m(), 1..n)
            })
            .prop_map(|(a, b, low)| {
                let mut limbs = a.limbs().to_vec();
                limbs[..low].copy_from_slice(&b.limbs()[..low]);
                (a, Mag::from_limbs(limbs))
            }),
    ];
    (pairs, any::<bool>())
        .prop_map(|((a, b), swap)| if swap { (b, a) } else { (a, b) })
        .boxed()
}

/// A canonical fraction `(numerator, denominator)`, reduced by GMP.
fn canonical(num: &rug::Integer, den: &rug::Integer) -> (Signed, Signed) {
    let q = rug::Rational::from((num, den));
    (Signed::from_rug(q.numer()), Signed::from_rug(q.denom()))
}

fn nonzero_integer() -> BoxedStrategy<Signed> {
    integer()
        .prop_filter("non-zero", |d| !d.mag().is_zero())
        .boxed()
}

/// A canonical rational `(numerator, denominator)`: arbitrary fractions,
/// integers, small fractions, and denominators at the promotion edges of
/// `Integer`'s inline representation.
pub(crate) fn rational() -> BoxedStrategy<(Signed, Signed)> {
    let one = || Just(Signed::from_i64(1));
    prop_oneof![
        4 => (integer(), nonzero_integer()),
        1 => (integer(), one()),
        2 => (-30i64..=30, 1i64..=30).prop_map(|(n, d)| (Signed::from_i64(n), Signed::from_i64(d))),
        1 => (integer(), signed_near_i64_bounds().prop_filter("non-zero", |d| !d.mag().is_zero())),
    ]
    .prop_map(|(n, d)| canonical(&n.to_rug(), &d.to_rug()))
    .boxed()
}

/// Pairs of canonical rationals: independent, with denominators sharing a
/// factor (the cancellation in Henrici's sum), or with each numerator sharing
/// a factor with the other's denominator (the cross-cancellation in
/// products), and equal or opposite values.
pub(crate) fn rational_pair() -> BoxedStrategy<((Signed, Signed), (Signed, Signed))> {
    let part = || signed_up_to(DEFAULT_MAX_LIMBS / 3);
    let nonzero = || natural_up_to(DEFAULT_MAX_LIMBS / 3).prop_filter("non-zero", |m| !m.is_zero());
    prop_oneof![
        3 => (rational(), rational()),
        2 => (part(), nonzero(), part(), nonzero(), nonzero()).prop_map(|(a, b, c, d, g)| {
            let g = g.to_rug();
            (canonical(&a.to_rug(), &(b.to_rug() * &g)), canonical(&c.to_rug(), &(d.to_rug() * &g)))
        }),
        2 => (part(), nonzero(), part(), nonzero(), nonzero(), nonzero()).prop_map(|(a, b, c, d, g, h)| {
            let (g, h) = (g.to_rug(), h.to_rug());
            (
                canonical(&(a.to_rug() * &g), &(b.to_rug() * &h)),
                canonical(&(c.to_rug() * &h), &(d.to_rug() * &g)),
            )
        }),
        1 => (rational(), any::<bool>()).prop_map(|((n, d), negate)| {
            let m = if negate { Signed::new(!n.is_negative(), n.mag().clone()) } else { n.clone() };
            ((n, d.clone()), (m, d))
        }),
    ]
    .boxed()
}
