//! The NTT product checked against GMP at every shape, including the
//! largest coefficients (all-ones operands) and transforms long enough to
//! recurse; its modular helpers checked against `u128` remainders.

use proptest::prelude::*;
use rug::Integer;

use super::*;
use crate::testkit::{strategy, tested};

/// [`NTT_BLOCK`] within test sizes.
const BLOCK: usize = tested(NTT_BLOCK);

/// Filler for output and scratch buffers, so results can't rely on zeroed memory.
const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn val(a: &[Limb]) -> Integer {
    Integer::from_digits(a, rug::integer::Order::Lsf)
}

fn run_mul(a: &[Limb], b: &[Limb]) -> Integer {
    let mut r = vec![GARBAGE; a.len() + b.len()];
    let mut scratch = vec![GARBAGE; mul_scratch_len(a.len(), b.len())];
    mul(&mut r, a, b, &mut scratch);
    val(&r)
}

fn run_sqr(a: &[Limb]) -> Integer {
    let mut r = vec![GARBAGE; 2 * a.len()];
    let mut scratch = vec![GARBAGE; sqr_scratch_len(a.len())];
    sqr(&mut r, a, &mut scratch);
    val(&r)
}

#[test]
fn primes_and_roots() {
    let mut product = Integer::from(1);
    for (p, g) in PRIMES {
        assert!(p >> 61 == 1, "{p:#x} not in (2^61, 2^62)");
        assert_eq!((p - 1) % (3 << ROOT_ORDER_LOG2), 0, "{p:#x}");
        assert_ne!(
            Integer::from(p).is_probably_prime(40),
            rug::integer::IsPrime::No
        );
        assert!(crate::prime::is_prime_u64(p));
        let m = Modulus::new(p);
        // Order exactly 3 2^53: g^(3 2^52) = -1 and g^(2^53) != 1.
        assert_eq!(m.pow(g, 3 << (ROOT_ORDER_LOG2 - 1)), p - 1, "{p:#x}");
        assert_ne!(m.pow(g, 1 << ROOT_ORDER_LOG2), 1, "{p:#x}");
        product *= p;
    }
    // Coefficients are below min(an, bn) 2^128, and min(an, bn) is at most
    // half the longest transform, 3 2^52 < 2^54.
    assert!(product > Integer::from(1) << (54 + 128));
}

fn residue(x: u128, p: Limb) -> Limb {
    (x % u128::from(p)) as Limb
}

proptest! {
    #[test]
    fn modular_helpers(
        i in 0..3usize,
        x in prop_oneof![strategy::limb(), Just(Limb::MAX), Just(0)],
        w in strategy::limb(),
        y in strategy::limb(),
    ) {
        let m = Modulus::new(PRIMES[i].0);
        let p = m.p;
        let w = w % p;
        let ws = m.shoup(w);
        prop_assert_eq!(u128::from(ws), (u128::from(w) << 64) / u128::from(p));
        let s = m.mul_shoup(x, w, ws);
        prop_assert!(s < 2 * p);
        prop_assert_eq!(s % p, residue(u128::from(x) * u128::from(w), p));
        let f = m.reduce_limb(x);
        prop_assert!(f < 2 * p);
        prop_assert_eq!(f % p, x % p);
        // Montgomery reduction of a product of lazy residues.
        let (a, b) = (m.reduce_limb(x), m.reduce_limb(y));
        let t = u128::from(a) * u128::from(b);
        let red = m.redc(t);
        prop_assert!(red < 2 * p);
        let r_mod_p = residue(1u128 << 64, p);
        prop_assert_eq!(m.mul(red % p, r_mod_p), residue(t, p));
    }

    #[test]
    fn mul_matches(
        (a, b) in (1..=300usize, 1..=300usize)
            .prop_flat_map(|(n, k)| (strategy::raw_limbs(n), strategy::raw_limbs(k))),
    ) {
        prop_assert_eq!(run_mul(&a, &b), val(&a) * val(&b));
    }

    #[test]
    fn sqr_matches(a in (1..=300usize).prop_flat_map(strategy::raw_limbs)) {
        prop_assert_eq!(run_sqr(&a), val(&a).square());
    }

    /// Transforms past [`NTT_BLOCK`], which split into halves.
    #[test]
    fn mul_matches_recursive(
        (a, b) in (BLOCK / 2..=2 * BLOCK, 1..=2 * BLOCK)
            .prop_flat_map(|(n, k)| (strategy::raw_limbs(n), strategy::raw_limbs(k))),
    ) {
        prop_assert_eq!(run_mul(&a, &b), val(&a) * val(&b));
    }
}

/// `a * b mod (B^n - 1)` from [`mul_cyclic`], as a residue in `[0, B^n - 1)`.
fn run_cyclic(n: usize, a: &[Limb], b: &[Limb]) -> Integer {
    let mut r = vec![GARBAGE; n];
    let mut scratch = vec![GARBAGE; cyclic_scratch_len(n)];
    mul_cyclic(&mut r, a, b, &mut scratch);
    val(&r) % ((Integer::from(1) << (64 * n)) - 1u32)
}

fn modulo_bnm1(x: Integer, n: usize) -> Integer {
    x % ((Integer::from(1) << (64 * n)) - 1u32)
}

proptest! {
    /// Cyclic products of every operand length up to the transform length
    /// (`2^k` and `3 2^k`), through the transform's recursion.
    #[test]
    fn mul_cyclic_matches(
        (n, a, b) in (0..=12u32, any::<bool>())
            .prop_map(|(k, three)| if three { 3usize << k.min(10) } else { 1usize << k })
            .prop_flat_map(|n| (Just(n), 1..=n, 1..=n))
            .prop_flat_map(|(n, an, bn)| (Just(n), strategy::raw_limbs(an), strategy::raw_limbs(bn))),
    ) {
        prop_assert_eq!(run_cyclic(n, &a, &b), modulo_bnm1(val(&a) * val(&b), n));
    }
}

/// All-ones operands of the full transform length: every coefficient is the
/// largest, `n (B - 1)^2`, and the carry out of the top limb is two limbs.
#[test]
fn mul_cyclic_all_ones() {
    for n in [
        1,
        2,
        3,
        4,
        6,
        64,
        96,
        BLOCK,
        3 * BLOCK / 2,
        2 * BLOCK,
        3 * BLOCK,
    ] {
        for (an, bn) in [(n, n), (n, 1), (n / 2 + 1, n)] {
            let (a, b) = (vec![Limb::MAX; an], vec![Limb::MAX; bn]);
            assert_eq!(
                run_cyclic(n, &a, &b),
                modulo_bnm1(val(&a) * val(&b), n),
                "n = {n}, ({an}, {bn})"
            );
        }
    }
}

/// The planner's choices at the shapes these tests use: radix-3 lengths
/// covering the product, radix-3 lengths with wrapped coefficients, and
/// powers of two on either side.
#[test]
fn transform_lengths() {
    // len = an + bn - 1; the cheapest by plan's cost model, which these
    // named lengths assume with the default `ntt_wrap_fraction` of 1.
    if WRAP_FRACTION == 1 {
        assert_eq!(transform_len(2, 2, true), 2);
        assert_eq!(transform_len(5, 4, true), 8);
        assert_eq!(transform_len(300, 300, true), 512);
        assert_eq!(transform_len(384, 384, true), 768);
        assert_eq!(transform_len(1536, 1537, true), 3072);
        assert_eq!(transform_len(1536, 1537, false), 2048);
        assert_eq!(transform_len(1600, 1700, true), 3072);
        assert_eq!(transform_len(1600, 1700, false), 4096);
        assert_eq!(transform_len(2048, 2048, true), 4096);
        assert_eq!(transform_len(3072, 3072, true), 6144);
        assert_eq!(transform_len(3072, 3072, false), 4096);
    }
    for (an, bn) in [(1, 1), (3, 2), (7, 9), (40, 40), (300, 301), (1000, 3)] {
        let n = transform_len(an, bn, true);
        let (half, three) = radix_parts(n);
        assert!(half.is_power_of_two() && n == half * if three { 3 } else { 1 });
        assert!(n >= an.max(bn));
    }
}

/// All-ones operands give the largest coefficients, `min(an, bn) (B - 1)^2`,
/// which the recombination must carry exactly; the lengths straddle the
/// powers of two and the recursion threshold.
#[test]
fn all_ones() {
    let shapes = [
        (1, 1),
        (1, 2),
        (2, 2),
        (3, 2),
        (5, 4),
        (512, 512),
        (513, 512),
        (BLOCK, 1),
        (BLOCK / 2, BLOCK / 2 + 1),
        (BLOCK, BLOCK),
        (BLOCK + 1, BLOCK),
        (3 * BLOCK, 2 * BLOCK + 7),
        (300, 300),
        (1174, 1175),
        (2048, 301),
        (2048, 2),
        // Radix-3 transforms: covering (3, 6, 768, 3072) and wrapped.
        (2, 2),
        (3, 3),
        (384, 384),
        (1536, 1537),
        (1600, 1700),
        (3 * BLOCK / 2, 3 * BLOCK / 2 + 2),
    ];
    for (an, bn) in shapes {
        let (a, b) = (vec![Limb::MAX; an], vec![Limb::MAX; bn]);
        assert_eq!(run_mul(&a, &b), val(&a) * val(&b), "({an}, {bn})");
        assert_eq!(run_mul(&b, &a), val(&a) * val(&b), "({bn}, {an})");
        assert_eq!(run_sqr(&a), val(&a).square(), "{an}");
    }
}

/// Shapes whose top coefficients wrap around the transform, once or (the
/// recomputed top product wrapping again) twice, and the limits of wrapping.
/// The named shapes assume the default `ntt_wrap_fraction` of 1.
#[test]
fn wrap_shapes() {
    let wraps = |an: usize, bn: usize| transform_len(an, bn, false) < an + bn - 1;
    if WRAP_FRACTION == 1 {
        assert!(wraps(300, 300));
        assert!(wraps(2048, 2));
        // Top product of 300 x 300 limbs, which wraps again.
        assert!(wraps(1174, 1175) && wraps(2048, 301));
        assert_eq!(1174 + 1175 - 1 - transform_len(1174, 1175, false), 300);
    }
    // An operand longer than the shorter transform, or a power-of-two length.
    assert!(!wraps(2049, 2));
    assert!(!wraps(512, 512) && !wraps(1, 1));
    // The largest wrapped part e = len - half for balanced operands.
    let half = 1024;
    let e_max = (half / WRAP_FRACTION).div_ceil(2);
    let shape = |e: usize| (half + e + 1 - half.midpoint(e), half.midpoint(e));
    if e_max >= 1 {
        let (an, bn) = shape(e_max);
        assert!(an.max(bn) <= half && wraps(an, bn));
    }
    let (an, bn) = shape(e_max + 1);
    assert!(an.max(bn) <= half && !wraps(an, bn));
}

#[test]
#[should_panic(expected = "empty operand")]
fn empty_operand_panics() {
    let mut r = [0; 1];
    mul(&mut r, &[1], &[], &mut []);
}

/// `a * b mod (B^n - 1)` with `b` prepared once by [`prepare_cyclic`], for
/// each `a` in turn.
fn run_prepared(n: usize, a: &[&[Limb]], b: &[Limb]) -> Vec<Integer> {
    let mut prep = vec![GARBAGE; prepared_len(n)];
    let mut scratch = vec![GARBAGE; cyclic_prepared_scratch_len(n)];
    prepare_cyclic(&mut prep, n, b, &mut scratch);
    a.iter()
        .map(|a| {
            let mut r = vec![GARBAGE; n];
            scratch.fill(GARBAGE);
            mul_cyclic_prepared(&mut r, a, &prep, &mut scratch);
            modulo_bnm1(val(&r), n)
        })
        .collect()
}

proptest! {
    /// A prepared operand gives the same products as [`mul_cyclic`], used
    /// several times.
    #[test]
    fn mul_cyclic_prepared_matches(
        (n, a1, a2, b) in (0..=11u32, any::<bool>())
            .prop_map(|(k, three)| if three { 3usize << k.min(9) } else { 1usize << k })
            .prop_flat_map(|n| (Just(n), 1..=n, 1..=n, 1..=n))
            .prop_flat_map(|(n, x, y, z)| {
                (Just(n), strategy::raw_limbs(x), strategy::raw_limbs(y), strategy::raw_limbs(z))
            }),
    ) {
        let got = run_prepared(n, &[&a1, &a2, &a1], &b);
        let want = |a: &[Limb]| modulo_bnm1(val(a) * val(&b), n);
        prop_assert_eq!(got, vec![want(&a1), want(&a2), want(&a1)]);
    }
}

/// The twiddle cache serves a short transform after a long one (from the
/// long one's prefix), radix-3 lengths next to powers of two, and threads
/// asking for different lengths at once.
#[test]
fn table_cache_lengths_and_threads() {
    let lens = [4096usize, 64, 3 * 512, 2 * BLOCK, 6, 1 << 13, 3 << 11, 32];
    let check = |n: usize, seed: u64| {
        let a: Vec<Limb> = (0..n as u64)
            .map(|i| (i + seed).wrapping_mul(0x9e37_79b9_7f4a_7c15))
            .collect();
        let b: Vec<Limb> = (0..=(n as u64 / 2))
            .map(|i| (i ^ seed).wrapping_mul(0xbf58_476d_1ce4_e5b9))
            .collect();
        assert_eq!(
            run_cyclic(n, &a, &b),
            modulo_bnm1(val(&a) * val(&b), n),
            "n = {n}"
        );
        assert_eq!(run_mul(&a, &b), val(&a) * val(&b), "n = {n}");
    };
    for (i, &n) in lens.iter().enumerate() {
        check(n, i as u64);
    }
    std::thread::scope(|scope| {
        for t in 0..4u64 {
            scope.spawn(move || {
                for (i, &n) in lens.iter().rev().enumerate() {
                    check(n, 100 * t + i as u64);
                }
            });
        }
    });
}
