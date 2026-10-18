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
        assert_eq!((p - 1) % (1 << ROOT_ORDER_LOG2), 0, "{p:#x}");
        assert_ne!(
            Integer::from(p).is_probably_prime(40),
            rug::integer::IsPrime::No
        );
        assert!(crate::prime::is_prime_u64(p));
        let m = Modulus::new(p);
        // Order exactly 2^54: g^(2^53) = -1.
        assert_eq!(m.pow(g, 1 << (ROOT_ORDER_LOG2 - 1)), p - 1, "{p:#x}");
        product *= p;
    }
    // Coefficients are below min(an, bn) 2^128 <= 2^(53 + 128).
    assert!(product > Integer::from(1) << (ROOT_ORDER_LOG2 - 1 + 128));
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
    /// Cyclic products of every operand length up to the transform length,
    /// through the transform's recursion.
    #[test]
    fn mul_cyclic_matches(
        (n, a, b) in (0..=12u32)
            .prop_map(|k| 1usize << k)
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
    for n in [1, 2, 4, 64, BLOCK, 2 * BLOCK] {
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
    let wraps = |an: usize, bn: usize| transform_len(an, bn) < an + bn - 1;
    if WRAP_FRACTION == 1 {
        assert!(wraps(300, 300));
        assert!(wraps(2048, 2));
        // Top product of 300 x 300 limbs, which wraps again.
        assert!(wraps(1174, 1175) && wraps(2048, 301));
        assert_eq!(1174 + 1175 - 1 - transform_len(1174, 1175), 300);
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
