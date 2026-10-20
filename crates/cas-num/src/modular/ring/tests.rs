//! `ModRing` checked against GMP, on both sides of `REDC_N_THRESHOLD` and
//! for odd (Montgomery) and even (plain) moduli.

use proptest::prelude::*;
use rug::Integer;

use super::{Form, ModRing, REDC_N_THRESHOLD, Redc, RingWork, redc_n, redc_wrap_len, window_bits};
use crate::kernel::mul::split::split_wrap_len;
use crate::kernel::mul::{
    MUL_WRAP_NTT_THRESHOLD, MUL_WRAP_SPLIT_THRESHOLD, wrap_len, wraps_cyclically,
};
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

/// `len` deterministic pseudo-random limbs (`SplitMix64` from `seed`).
fn splitmix_limbs(len: usize, seed: u64) -> Vec<Limb> {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        })
        .collect()
}

/// Modulus lengths where `ModRing` changes how REDC takes its products, with
/// their neighbours: the split of `B^m - 1` from its threshold (or REDC's,
/// if higher); `ninv` prepared from the length whose double has a cyclic
/// transform length; `n` prepared once the length's own transform length is
/// cyclic, and the jump of that length past a power of two. Thresholds beyond
/// test sizes are left out.
fn redc_choice_lens() -> Vec<usize> {
    use crate::testkit::TESTED_SIZE_MAX;
    let mut lens = Vec::new();
    let split = MUL_WRAP_SPLIT_THRESHOLD.max(REDC_N_THRESHOLD);
    if split <= TESTED_SIZE_MAX {
        lens.extend([split - 1, split, split + 1]);
    }
    if MUL_WRAP_NTT_THRESHOLD <= TESTED_SIZE_MAX {
        let cyclic = MUL_WRAP_NTT_THRESHOLD.next_power_of_two().max(4);
        lens.extend([
            cyclic / 2 - 1,
            cyclic / 2,
            cyclic / 2 + 1,
            cyclic * 3 / 4,
            cyclic * 3 / 4 + 1,
            cyclic,
            cyclic + 1,
        ]);
    }
    lens.retain(|&k| (2..=TESTED_SIZE_MAX).contains(&k));
    lens.sort_unstable();
    lens.dedup();
    lens
}

/// Products and squares at [`redc_choice_lens`], against GMP, for the moduli
/// `B^k - 1` (all ones), `2^(64k - 1) + 1`, `2^(64k - 1) - 1` and an odd
/// random one, of the operands `m - 1`, `B^k - 1` reduced, `2^(64k - 1)`
/// reduced and a random value; and a power of the random value.
#[test]
fn redc_choices_match() {
    for k in redc_choice_lens() {
        let bits = 64 * u32::try_from(k).unwrap();
        let b = Integer::from(1) << bits;
        let half = Integer::from(1) << (bits - 1);
        let random = Mag::from_limbs(splitmix_limbs(k, k as u64)).to_rug() | 1u32;
        for m in [
            Integer::from(&b - 1u32),
            Integer::from(&half + 1u32),
            Integer::from(&half - 1u32),
            random,
        ] {
            let ring = ModRing::new(Mag::from_rug(&m).to_natural());
            let random = Mag::from_limbs(splitmix_limbs(k, 7 + k as u64)).to_rug();
            let values = [
                Integer::from(&m - 1u32),
                Integer::from(&b - 1u32) % &m,
                Integer::from(&half % &m),
                random % &m,
            ];
            let elem = |v: &Integer| ring.elem(&Mag::from_rug(v).to_natural());
            let value = |x: &super::Residue| ring.to_natural(x).to_rug();
            for (i, x) in values.iter().enumerate() {
                let rx = elem(x);
                let y = &values[(i + 1) % values.len()];
                let want = Integer::from(x * y) % &m;
                assert_eq!(value(&ring.mul(&rx, &elem(y))), want, "mul, k = {k}");
                let want = Integer::from(x * x) % &m;
                assert_eq!(value(&ring.square(&rx)), want, "square, k = {k}");
            }
            let (x, e) = (&values[3], Integer::from(0xf7u32));
            let want = x.clone().pow_mod(&e, &m).unwrap();
            let got = ring.pow(&elem(x), &Mag::from_rug(&e).to_natural());
            assert_eq!(value(&got), want, "pow, k = {k}");
        }
    }
}

/// `RingWork`'s in-place operations, each writing over a residue that holds
/// another value, against `ModRing`'s allocating ones.
fn in_place_matches(m: &Mag, a: &Mag, b: &Mag) -> Result<(), TestCaseError> {
    let r = ring(m);
    let w = RingWork::new(&r);
    let (x, y) = (r.elem(&a.to_natural()), r.elem(&b.to_natural()));
    let mut out = r.sub(&y, &x);
    w.add(&mut out, &x, &y);
    prop_assert_eq!(&out, &r.add(&x, &y));
    w.sub(&mut out, &x, &y);
    prop_assert_eq!(&out, &r.sub(&x, &y));
    w.neg(&mut out, &x);
    prop_assert_eq!(&out, &r.neg(&x));
    w.mul(&mut out, &x, &y);
    prop_assert_eq!(&out, &r.mul(&x, &y));
    w.square(&mut out, &x);
    prop_assert_eq!(&out, &r.square(&x));
    // The product buffers are reused: a second product is still right.
    w.mul(&mut out, &y, &x);
    prop_assert_eq!(&out, &r.mul(&x, &y));
    out.clone_from(&y);
    prop_assert_eq!(&out, &y);
    Ok(())
}

/// `RingWork::lucas_double` from residues of `a` and `b` modulo `m` (made
/// odd and at least 3), against the doubling formulas evaluated by GMP:
/// `(B - C, B - q A)` for `A = a^2`, `B = b^2`, `C = (a - b)^2`, then
/// `(B - q A, B - q A - q (B - C))` if `up`.
fn lucas_matches(m: &Integer, a: &Integer, b: &Integer, q: i64, up: bool) {
    let n = Integer::from(m | 1u32).max(Integer::from(3));
    let (x, y) = (Integer::from(a % &n), Integer::from(b % &n));
    let r = ModRing::new(Mag::from_rug(&n).to_natural());
    let w = RingWork::new(&r);
    let elem = |v: &Integer| r.elem(&Mag::from_rug(v).to_natural());
    let (mut u, mut u1) = (elem(&x), elem(&y));
    w.lucas_double(&mut u, &mut u1, q, up);
    let md = |v: Integer| v.modulo(&n);
    let sq = |v: &Integer| Integer::from(v * v);
    let even = md(sq(&y) - sq(&Integer::from(&x - &y)));
    let odd = md(sq(&y) - Integer::from(q) * sq(&x));
    let want = if up {
        let next = md(&odd - Integer::from(q) * &even);
        (odd, next)
    } else {
        (even, odd)
    };
    let got = (r.to_natural(&u).to_rug(), r.to_natural(&u1).to_rug());
    assert_eq!(got, want, "n = {n}, q = {q}, up = {up}");
    // Canonical: the same residue as built from the value.
    assert_eq!(u, elem(&want.0));
    assert_eq!(u1, elem(&want.1));
}

/// The Lucas step's extremes: `|Q| = 255`, residues `n - 1`, and moduli
/// `B^k - 1` (top limb all ones) and `B^(k-1) + 1` (top limb 1, the
/// largest normalising shift).
#[test]
fn lucas_double_extremes() {
    for k in 1..=12u32 {
        let full = (Integer::from(1) << (64 * k)) - 1;
        let low = (Integer::from(1) << (64 * (k - 1))) + 1;
        for n in [full, low] {
            let top = Integer::from(&n - 1);
            for q in [-255, -1, 1, 2, 255] {
                for up in [false, true] {
                    lucas_matches(&n, &top, &top, q, up);
                    lucas_matches(&n, &top, &Integer::new(), q, up);
                    lucas_matches(&n, &Integer::new(), &top, q, up);
                }
            }
        }
    }
}

proptest! {
    #[test]
    fn lucas_double_matches(
        (m, a, b) in operands(small_modulus()),
        q in -255..=255i64,
        up in any::<bool>(),
    ) {
        lucas_matches(&m.to_rug(), &a.to_rug(), &b.to_rug(), q, up);
    }

    #[test]
    fn large_lucas_double_matches(
        (m, a, b) in operands(large_modulus()),
        q in -255..=255i64,
        up in any::<bool>(),
    ) {
        lucas_matches(&m.to_rug(), &a.to_rug(), &b.to_rug(), q, up);
    }

    #[test]
    fn ring_work_matches((m, a, b) in operands(small_modulus())) {
        in_place_matches(&m, &a, &b)?;
    }

    #[test]
    fn large_ring_work_matches((m, a, b) in operands(large_modulus())) {
        in_place_matches(&m, &a, &b)?;
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

/// `(t R^-1 mod n, the result of redc_n)` for an odd `n` and `t < n R`, with
/// `q n` taken each way that applies (the full product, and modulo `B^m - 1`
/// for the ring's choice of `m`, the shortest transform length and the
/// split's length) and `q` by `mullo` or, where the transform applies, from
/// `ninv` prepared at the shortest transform length from `2k`. Every way must
/// give the same result `r + carry R`, which must be `< 2n`.
fn run_redc_n(n: &[Limb], t: &[Limb]) -> (Integer, Integer) {
    let k = n.len();
    let mut wraps = vec![0, redc_wrap_len(k), wrap_len(k), split_wrap_len(k)];
    wraps.sort_unstable();
    wraps.dedup();
    let mut ninv_lens = vec![0];
    if wraps_cyclically(wrap_len(2 * k)) {
        ninv_lens.push(wrap_len(2 * k));
    }
    let mut results = Vec::new();
    for &wrap in &wraps {
        for &ninv_len in &ninv_lens {
            results.push((
                wrap,
                ninv_len,
                run_redc_n_with(n, t, &Redc::with(n, wrap, ninv_len)),
            ));
        }
    }
    results.push((0, 0, run_redc_n_with(n, t, &Redc::new(n))));
    for (wrap, ninv_len, r) in &results {
        assert_eq!(
            r, &results[0].2,
            "redc_n with wrap {wrap} and ninv_len {ninv_len} differs, k = {k}"
        );
    }
    results.swap_remove(0).2
}

/// [`run_redc_n`] for one choice of products.
fn run_redc_n_with(n: &[Limb], t: &[Limb], redc: &Redc) -> (Integer, Integer) {
    let k = n.len();
    let (nv, tv) = (
        Mag::from_limbs(n.to_vec()).to_rug(),
        Mag::from_limbs(t.to_vec()).to_rug(),
    );
    let r_big = Integer::from(1) << (64 * u32::try_from(k).unwrap());
    let ninv = &r_big - nv.clone().invert(&r_big).unwrap();
    assert_eq!(
        Mag::from_limbs(redc.ninv.clone()).to_rug(),
        ninv,
        "-n^-1 mod R"
    );
    let mut r = vec![0xdead_beef; k];
    let mut u = vec![0xdead_beef; 4 * k];
    let mut work = vec![0xdead_beef; redc.scratch_len(k)];
    let carry = redc_n(&mut r, t, n, redc, &mut u, &mut work);
    let got =
        Mag::from_limbs(r).to_rug() + (Integer::from(carry) << (64 * u32::try_from(k).unwrap()));
    assert!(got < Integer::from(2 * &nv), "redc_n: result >= 2n");
    let want = tv * r_big.invert(&nv).unwrap() % &nv;
    (want, got % nv)
}

/// Modulus lengths for `redc_n`: small ones, lengths around the powers of
/// two from which `q n` is taken modulo `B^m - 1` by a cyclic transform, and
/// lengths around the split's threshold, odd and even, including those
/// rounded up for three halvings.
fn redc_lens() -> BoxedStrategy<usize> {
    let m = tested(MUL_WRAP_NTT_THRESHOLD).next_power_of_two().max(4);
    let s = tested(MUL_WRAP_SPLIT_THRESHOLD).max(2);
    prop_oneof![
        1..=12usize,
        Just(m / 2),
        Just(m / 2 + 1),
        Just(m * 3 / 4),
        Just(m - 1),
        Just(m),
        Just(m + 1),
        s - 1..=s + 2,
        Just(2 * s + 1),
        Just(4 * s + 3),
        Just(8 * s + 5),
        Just(8 * s + 8),
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
    // Wrap lengths of 2^j and 3 2^j: m, 3m/4 (for k from 5m/8 + 1), m/2.
    for k in [
        1,
        2,
        3,
        m / 2 + 1,
        m * 5 / 8 + 1,
        m * 3 / 4,
        m * 3 / 4 + 1,
        m,
    ] {
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
        for k in [
            m / 2,
            m / 2 + 1,
            m * 5 / 8 + 1,
            m * 3 / 4,
            m * 3 / 4 + 1,
            m,
            m + 1,
        ] {
            assert_eq!(
                wraps_cyclically(redc_wrap_len(k)),
                wrap_len(k) >= MUL_WRAP_NTT_THRESHOLD,
                "k = {k}"
            );
        }
        assert!(redc_wrap_len(m) == m && redc_wrap_len(m + 1) == wrap_len(m + 1));
        assert_eq!(wrap_len(m * 5 / 8 + 1), m * 3 / 4);
    }
}

#[test]
fn forms() {
    let ninv_len = |limbs: usize| {
        let ring = ModRing::new(Natural::power_of_two(64 * limbs as u64) - 1u64);
        let Form::Montgomery { redc, .. } = &ring.form else {
            panic!("odd modulus not in Montgomery form");
        };
        redc.as_ref().map_or(0, |r| r.ninv.len())
    };
    if REDC_N_THRESHOLD <= crate::testkit::TESTED_SIZE_MAX {
        assert_eq!(ninv_len(REDC_N_THRESHOLD), REDC_N_THRESHOLD);
    }
    if REDC_N_THRESHOLD > 1 {
        let short = crate::testkit::tested(REDC_N_THRESHOLD - 1);
        assert_eq!(ninv_len(short), 0);
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
        if let Form::Montgomery { redc, .. } = &mut one.form {
            *redc = None;
        }
        if let Form::Montgomery { redc, .. } = &mut full.form {
            *redc = Some(Redc::new(n.limbs()));
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
