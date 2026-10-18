//! `ModRing` checked against GMP, on both sides of `REDC_N_THRESHOLD` and
//! for odd (Montgomery) and even (plain) moduli.

use proptest::prelude::*;
use rug::Integer;

use super::{Form, ModRing, REDC_N_THRESHOLD, redc_n, redc_n_scratch_len, redc_wraps, window_bits};
use crate::kernel::mul::MUL_WRAP_NTT_THRESHOLD;
use crate::testkit::{Mag, ToRug, differential, strategy, tested};
use crate::{Limb, Natural};

/// A modulus of up to 8 limbs (one limb is `ModWord`'s job, but works here
/// too), mostly odd (Montgomery form), or `2^k ± 1`.
fn small_modulus() -> BoxedStrategy<Mag> {
    modulus(strategy::natural_up_to(8), 64 * 8)
}

/// A modulus of `REDC_N_THRESHOLD - 1` to `REDC_N_THRESHOLD + 1` limbs, so
/// both REDC variants run.
fn large_modulus() -> BoxedStrategy<Mag> {
    modulus(
        strategy::natural_with_len(strategy::len_around(&[REDC_N_THRESHOLD])).boxed(),
        u32::try_from(64 * crate::testkit::tested(REDC_N_THRESHOLD) + 64).unwrap(),
    )
}

fn modulus(any: BoxedStrategy<Mag>, max_bits: u32) -> BoxedStrategy<Mag> {
    let odd = any.clone().prop_map(|m| Mag::from_rug(&(m.to_rug() | 1)));
    let near_power = (max_bits - 128..=max_bits, prop::bool::ANY).prop_map(|(k, plus)| {
        let p = Integer::from(1) << k;
        Mag::from_rug(&if plus { p + 1 } else { p - 1 })
    });
    prop_oneof![3 => odd, 1 => any, 1 => near_power]
        .prop_filter("non-zero", |m| !m.is_zero())
        .boxed()
}

/// A modulus and two operands: 0, 1, `m - 1`, `m + 1` or any value of up to
/// one limb more than `m` (reduced by `elem`).
fn operands(modulus: BoxedStrategy<Mag>) -> BoxedStrategy<(Mag, Mag, Mag)> {
    modulus
        .prop_flat_map(|m| {
            let near = move |m: &Mag| {
                let m = m.to_rug();
                prop_oneof![
                    1 => Just(Mag::default()),
                    1 => Just(Mag::from_u64(1)),
                    1 => Just(Mag::from_rug(&(m.clone() - 1))),
                    1 => Just(Mag::from_rug(&(m.clone() + 1))),
                    4 => strategy::natural_up_to(m.significant_digits::<u64>() + 1),
                ]
            };
            (Just(m.clone()), near(&m), near(&m))
        })
        .boxed()
}

fn ring(m: &Mag) -> ModRing {
    ModRing::new(m.to_natural())
}

/// `op` on the residues of `a` and `b`, as a value.
fn binary(
    m: &Mag,
    a: &Mag,
    b: &Mag,
    op: fn(&ModRing, &super::Residue, &super::Residue) -> super::Residue,
) -> Natural {
    let r = ring(m);
    r.to_natural(&op(&r, &r.elem(&a.to_natural()), &r.elem(&b.to_natural())))
}

differential! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    fn elem_round_trip((m, a, _b) in operands(small_modulus())) {
        ours: { let r = ring(&m); r.to_natural(&r.elem(&a.to_natural())) },
        gmp: a.to_rug().modulo(&m.to_rug()),
    }

    fn elem_integer((m, _a, _b) in operands(small_modulus()), a in strategy::integer()) {
        ours: { let r = ring(&m); r.to_natural(&r.elem_integer(&a.to_integer())) },
        gmp: a.to_rug().modulo(&m.to_rug()),
    }

    fn add((m, a, b) in operands(small_modulus())) {
        ours: binary(&m, &a, &b, ModRing::add),
        gmp: (a.to_rug() % m.to_rug() + b.to_rug() % m.to_rug()).modulo(&m.to_rug()),
    }

    fn sub((m, a, b) in operands(small_modulus())) {
        ours: binary(&m, &a, &b, ModRing::sub),
        gmp: (a.to_rug() - b.to_rug()).modulo(&m.to_rug()),
    }

    fn neg((m, a, _b) in operands(small_modulus())) {
        ours: { let r = ring(&m); r.to_natural(&r.neg(&r.elem(&a.to_natural()))) },
        gmp: (-a.to_rug()).modulo(&m.to_rug()),
    }

    fn mul((m, a, b) in operands(small_modulus())) {
        ours: binary(&m, &a, &b, ModRing::mul),
        gmp: (a.to_rug() * b.to_rug()).modulo(&m.to_rug()),
    }

    fn square((m, a, _b) in operands(small_modulus())) {
        ours: { let r = ring(&m); r.to_natural(&r.square(&r.elem(&a.to_natural()))) },
        gmp: a.to_rug().square().modulo(&m.to_rug()),
    }

    fn pow((m, a, _b) in operands(small_modulus()), e in strategy::natural_up_to(12)) {
        ours: { let r = ring(&m); r.to_natural(&r.pow(&r.elem(&a.to_natural()), &e.to_natural())) },
        gmp: a.to_rug().pow_mod(&e.to_rug(), &m.to_rug()).unwrap(),
    }

    fn inv((m, a, _b) in operands(small_modulus())) {
        ours: { let r = ring(&m); r.inv(&r.elem(&a.to_natural())).map(|x| r.to_natural(&x)) },
        gmp: a.to_rug().invert(&m.to_rug()).ok(),
    }

    fn is_zero((m, a, _b) in operands(small_modulus())) {
        ours: { let r = ring(&m); r.is_zero(&r.elem(&a.to_natural())) },
        gmp: a.to_rug().is_divisible(&m.to_rug()),
    }
}

differential! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    fn large_elem_round_trip((m, a, _b) in operands(large_modulus())) {
        ours: { let r = ring(&m); r.to_natural(&r.elem(&a.to_natural())) },
        gmp: a.to_rug().modulo(&m.to_rug()),
    }

    fn large_mul((m, a, b) in operands(large_modulus())) {
        ours: binary(&m, &a, &b, ModRing::mul),
        gmp: (a.to_rug() * b.to_rug()).modulo(&m.to_rug()),
    }

    fn large_square((m, a, _b) in operands(large_modulus())) {
        ours: { let r = ring(&m); r.to_natural(&r.square(&r.elem(&a.to_natural()))) },
        gmp: a.to_rug().square().modulo(&m.to_rug()),
    }

    fn large_pow((m, a, _b) in operands(large_modulus()), e in strategy::natural_up_to(1)) {
        ours: { let r = ring(&m); r.to_natural(&r.pow(&r.elem(&a.to_natural()), &e.to_natural())) },
        gmp: a.to_rug().pow_mod(&e.to_rug(), &m.to_rug()).unwrap(),
    }

    fn large_inv((m, a, _b) in operands(large_modulus())) {
        ours: { let r = ring(&m); r.inv(&r.elem(&a.to_natural())).map(|x| r.to_natural(&x)) },
        gmp: a.to_rug().invert(&m.to_rug()).ok(),
    }
}

/// A chain of operations stays consistent: residues built by different
/// routes compare equal, so the representation is unique.
#[test]
fn representation_is_unique() {
    for m in [
        "1000000007",
        "340282366920938463463374607431768211457",
        "18446744073709551616",
    ] {
        let n: Natural = m.parse().unwrap();
        let r = ModRing::new(n.clone());
        let x = r.elem(&Natural::from(12345u64));
        let y = r.mul(&r.add(&x, &r.one()), &r.sub(&x, &r.one()));
        let z = r.sub(&r.square(&x), &r.one());
        assert_eq!(y, z);
        assert_eq!(r.elem(&n), r.zero());
        assert_eq!(r.add(&r.elem(&(&n - 1u64)), &r.one()), r.zero());
        assert_eq!(r.pow(&x, &Natural::ZERO), r.one());
        assert_eq!(r.modulus(), &n);
    }
}

#[test]
fn modulus_one() {
    let r = ModRing::new(Natural::one());
    assert_eq!(r.one(), r.zero());
    let x = r.elem(&Natural::from(5u64));
    assert!(r.is_zero(&x));
    assert_eq!(r.inv(&x), Some(r.zero()));
    assert_eq!(r.to_natural(&r.pow(&x, &Natural::ZERO)), Natural::ZERO);
}

/// `(t R^-1 mod n, the result of redc_n)` for an odd `n` and `t < n R`; the
/// result is `r + carry R`, which must be `< 2n`.
fn run_redc_n(n: &[Limb], t: &[Limb]) -> (Integer, Integer) {
    let k = n.len();
    let (nv, tv) = (
        Mag::from_limbs(n.to_vec()).to_rug(),
        Mag::from_limbs(t.to_vec()).to_rug(),
    );
    let r_big = Integer::from(1) << (64 * u32::try_from(k).unwrap());
    let ninv = &r_big - nv.clone().invert(&r_big).unwrap();
    let mut ninv_full = Mag::from_rug(&ninv).limbs().to_vec();
    ninv_full.resize(k, 0);
    let mut r = vec![0xdead_beef; k];
    let mut u = vec![0xdead_beef; 4 * k];
    let mut work = vec![0xdead_beef; redc_n_scratch_len(k)];
    let carry = redc_n(&mut r, t, n, &ninv_full, &mut u, &mut work);
    let got =
        Mag::from_limbs(r).to_rug() + (Integer::from(carry) << (64 * u32::try_from(k).unwrap()));
    assert!(got < Integer::from(2 * &nv), "redc_n: result >= 2n");
    let want = tv * r_big.invert(&nv).unwrap() % &nv;
    (want, got % nv)
}

/// Modulus lengths for `redc_n`: small ones, and lengths around the powers of
/// two from which `q n` is taken modulo `B^m - 1`.
fn redc_lens() -> BoxedStrategy<usize> {
    let m = tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two().max(4);
    prop_oneof![
        1..=12usize,
        Just(m / 2),
        Just(m / 2 + 1),
        Just(m * 3 / 4),
        Just(m - 1),
        Just(m),
        Just(m + 1)
    ]
    .boxed()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// `t` random below `n R`, with a zero low half (`q = 0`'s neighbour: the
    /// low halves sum to 0, not `R`), or `n R - 1`.
    #[test]
    fn redc_n_matches(
        (n, t) in redc_lens().prop_flat_map(|k| (strategy::raw_limbs(k), strategy::raw_limbs(2 * k), 0..4u8))
            .prop_map(|(mut n, t, shape)| {
                let k = n.len();
                n[0] |= 1;
                n[k - 1] = n[k - 1].max(1);
                let nv = Mag::from_limbs(n.clone()).to_rug();
                let nr = nv << (64 * u32::try_from(k).unwrap());
                let tv = match shape {
                    0 => Integer::from(&nr - 1),
                    1 => (Mag::from_limbs(t[k..].to_vec()).to_rug() << (64 * u32::try_from(k).unwrap())) % &nr,
                    _ => Mag::from_limbs(t).to_rug() % &nr,
                };
                let mut t = Mag::from_rug(&tv).limbs().to_vec();
                t.resize(2 * k, 0);
                (n, t)
            }),
    ) {
        let (want, got) = run_redc_n(&n, &t);
        prop_assert_eq!(got, want);
    }
}

/// All-ones moduli (`B^k - 1`) and `t = n R - 1` at every wrap shape, and the
/// wrap rule itself.
#[test]
fn redc_n_extremes() {
    let m = tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two().max(4);
    for k in [1, 2, 3, m / 2 + 1, m * 3 / 4, m] {
        let n = vec![Limb::MAX; k];
        let nr = Mag::from_limbs(n.clone()).to_rug() << (64 * u32::try_from(k).unwrap());
        for tv in [Integer::from(&nr - 1), Integer::from(1), Integer::from(0)] {
            let mut t = Mag::from_rug(&tv).limbs().to_vec();
            t.resize(2 * k, 0);
            let (want, got) = run_redc_n(&n, &t);
            assert_eq!(got, want, "k = {k}");
        }
    }
    if MUL_WRAP_NTT_THRESHOLD <= crate::testkit::TESTED_SIZE_MAX {
        assert!(redc_wraps(m * 3 / 4) && redc_wraps(m));
        assert_eq!(redc_wraps(m / 2), m / 2 >= MUL_WRAP_NTT_THRESHOLD);
        assert!(redc_wraps(m + 1));
    }
}

#[test]
fn forms() {
    let ninv_full_len = |limbs: usize| {
        let ring = ModRing::new(Natural::power_of_two(64 * limbs as u64) - 1u64);
        let Form::Montgomery { ninv_full, .. } = &ring.form else {
            panic!("odd modulus not in Montgomery form");
        };
        ninv_full.len()
    };
    if REDC_N_THRESHOLD <= crate::testkit::TESTED_SIZE_MAX {
        assert_eq!(ninv_full_len(REDC_N_THRESHOLD), REDC_N_THRESHOLD);
    }
    if REDC_N_THRESHOLD > 1 {
        let short = crate::testkit::tested(REDC_N_THRESHOLD - 1);
        assert_eq!(ninv_full_len(short), 0);
    }
    assert!(matches!(
        ModRing::new(Natural::from(7u64)).form,
        Form::Montgomery { .. }
    ));
    let even = ModRing::new(Natural::power_of_two(200));
    assert!(matches!(even.form, Form::Plain));
}

#[test]
#[should_panic(expected = "zero modulus")]
fn zero_modulus_panics() {
    let _ = ModRing::new(Natural::ZERO);
}

/// Times one Montgomery product with each REDC variant, to place
/// `REDC_N_THRESHOLD`. Run with
/// `cargo test --release -p cas-num tune_redc -- --ignored --nocapture`.
#[test]
#[ignore = "benchmark"]
fn tune_redc() {
    use std::time::Instant;
    for k in [32, 64, 96, 128, 144, 160, 192, 256, 384] {
        let n = (Natural::power_of_two(64 * k as u64) - 1u64) / 3u64 * 2u64 + 1u64;
        assert_eq!(n.limb_len(), k);
        let time = |ring: &ModRing| {
            let mut x = ring.elem(&(&n >> 3));
            let start = Instant::now();
            let reps = 200_000 / u32::try_from(k).unwrap();
            for _ in 0..reps {
                x = ring.square(&x);
            }
            std::hint::black_box(&x);
            start.elapsed().as_secs_f64() * 1e9 / f64::from(reps)
        };
        let mut one = ModRing::new(n.clone());
        let mut full = ModRing::new(n.clone());
        let r = Natural::power_of_two(64 * k as u64);
        let inv = crate::Integer::from(&n)
            .mod_inverse(&crate::Integer::from(&r))
            .unwrap()
            .unsigned_abs();
        let inv = super::padded(&(&r - &inv), k).into_vec();
        if let Form::Montgomery { ninv_full, .. } = &mut one.form {
            ninv_full.clear();
        }
        if let Form::Montgomery { ninv_full, .. } = &mut full.form {
            *ninv_full = inv;
        }
        println!(
            "k = {k:4}: redc_1 {:9.0} ns, redc_n {:9.0} ns",
            time(&one),
            time(&full)
        );
    }
}

/// Each width takes over exactly where it saves multiplications:
/// `bits / (w + 1) + 2^(w - 1)` drops below the narrower width's count.
#[test]
#[allow(clippy::cast_precision_loss, reason = "counts below 2^20")]
fn window_widths() {
    let cost = |bits: u64, w: u64| bits as f64 / (w + 1) as f64 + (1u64 << (w - 1)) as f64;
    for bits in 1..20_000u64 {
        let w = window_bits(bits);
        assert!(cost(bits, w) <= cost(bits, w + 1), "{bits} bits: {w}");
        if w > 1 {
            assert!(cost(bits, w) < cost(bits, w - 1), "{bits} bits: {w}");
        }
    }
    assert_eq!(window_bits(6), 1);
    assert_eq!(window_bits(7), 2);
    assert_eq!(window_bits(4096), 7);
    assert_eq!(window_bits(4609), 8);
}
