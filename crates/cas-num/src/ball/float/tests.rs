//! `Float` against MPFR: rounded operations must give MPFR's value and ternary
//! value in every mode; exact ones must give the exact rational result.
#![allow(
    clippy::many_single_char_names,
    clippy::needless_pass_by_value,
    clippy::float_cmp,
    reason = "test names: operands a, b, x, y and MPFR's f, g; results passed as built; \
              conversions to f64 must be exact"
)]

use std::cmp::Ordering;

use proptest::prelude::*;
use rug::Float as MpFloat;

use super::{Float, Round};
use crate::testkit::float::{
    Dyadic, div_round, dyadic, dyadic_with, far_apart_pair, huge_shift, rug_round,
};
use crate::testkit::{Signed, ToRug, check_same, strategy};
use crate::{Integer, Limb, Rational};

/// Precisions at limb boundaries, small ones, and random ones.
fn prec() -> BoxedStrategy<u64> {
    prop_oneof![
        2 => proptest::sample::select(&[1u64, 2, 3, 53, 63, 64, 65, 127, 128, 129, 191, 192, 193, 256][..]),
        1 => 1..600u64,
    ]
    .boxed()
}

fn round_mode() -> impl Strategy<Value = Round> {
    proptest::sample::select(&Round::ALL[..])
}

/// [`dyadic_with`] without zero (mapped to one), for divisors.
fn nonzero(max_limbs: usize, exps: std::ops::RangeInclusive<i64>) -> BoxedStrategy<Dyadic> {
    dyadic_with(max_limbs, exps)
        .prop_map(|d| {
            if d.is_zero() {
                Dyadic::new(&rug::Integer::from(1), &rug::Integer::new())
            } else {
                d
            }
        })
        .boxed()
}

/// The exact MPFR value of `d`.
fn mp(d: &Dyadic) -> MpFloat {
    d.to_rug().expect("strategy stays in MPFR's range")
}

/// Fails unless ours equals MPFR's `(value, ternary)`.
fn check_round(
    ours: (Float, Ordering),
    expected: (MpFloat, Ordering),
) -> Result<(), TestCaseError> {
    let want = Dyadic::from_rug(&expected.0).expect("finite");
    let got = ours.0.to_rug();
    prop_assert_eq!(&got, &want, "value");
    prop_assert_eq!(ours.1, expected.1, "ternary of {:?}", got);
    Ok(())
}

/// A `k`-limb mantissa (top bit set) of exponent `e`: all ones (shape 0),
/// the top and lowest bits (1), the top bit (2), a top limb of ones over
/// limbs holding only their top bit (3), or the random `limbs` (4).
fn shaped(k: usize, shape: u8, limbs: [Limb; 4], neg: bool, e: i64) -> Float {
    let top = 1 << (Limb::BITS - 1);
    let mut m: Vec<Limb> = match shape {
        0 => vec![Limb::MAX; k],
        1 => {
            let mut v = vec![0; k];
            v[0] |= 1;
            v
        }
        2 => vec![0; k],
        3 => {
            let mut v = vec![top; k];
            v[k - 1] = Limb::MAX;
            v
        }
        _ => limbs[..k].to_vec(),
    };
    m[k - 1] |= top;
    let bits = i64::try_from(64 * k).unwrap();
    let mag = Integer::from_parts(neg, crate::Natural::from_limbs(&m));
    Float::from_integer_exp(&mag, &Integer::from(e - bits))
}

/// `x` scaled by `2^k`, which commutes with every operation tested here.
fn scaled(x: &Float, k: &Integer) -> Float {
    x.mul_2exp(k)
}

proptest! {
    #[test]
    fn dyadic_round_trip(a in dyadic()) {
        let x = a.to_ours();
        prop_assert!(x.is_canonical());
        prop_assert_eq!(x.to_rug(), a);
    }

    /// `from_raw` on arbitrary limbs, with zero limbs at both ends.
    #[test]
    fn from_raw_normalises(
        limbs in proptest::collection::vec(strategy::limb(), 0..8),
        neg: bool,
        e in -500i64..500,
    ) {
        let x = Float::from_raw(neg, &limbs, Integer::from(e));
        let m = crate::testkit::Mag::from_limbs(limbs).to_rug();
        let m = if neg { -m } else { m };
        prop_assert_eq!(x.to_rug(), Dyadic::new(&m, &rug::Integer::from(e)));
    }

    #[test]
    fn accessors(a in dyadic()) {
        let x = a.to_ours();
        let f = mp(&a);
        prop_assert_eq!(x.is_zero(), f.is_zero());
        prop_assert_eq!(x.sign(), f.cmp0().unwrap());
        if !x.is_zero() {
            prop_assert_eq!(x.exponent().to_i64(), f.get_exp().map(i64::from));
        }
        prop_assert_eq!(x.bit_length(), u64::from(a.mant_bits()));
        prop_assert_eq!(x.neg().to_rug(), Dyadic::from_rug(&(-f.clone())).unwrap());
        prop_assert_eq!(x.abs().to_rug(), Dyadic::from_rug(&f.clone().abs()).unwrap());
    }

    #[test]
    fn exact_ops(a in dyadic_with(6, -500..=500), b in dyadic_with(6, -500..=500)) {
        let (x, y) = (a.to_ours(), b.to_ours());
        let (p, q) = (a.to_rational(), b.to_rational());
        prop_assert_eq!(x.add(&y).to_rational().to_rug(), rug::Rational::from(&p + &q));
        prop_assert_eq!(x.sub(&y).to_rational().to_rug(), rug::Rational::from(&p - &q));
        prop_assert_eq!(x.mul(&y).to_rational().to_rug(), rug::Rational::from(&p * &q));
        prop_assert!(x.add(&y).is_canonical() && x.sub(&y).is_canonical() && x.mul(&y).is_canonical());
    }

    #[test]
    fn comparisons(a in dyadic(), b in dyadic()) {
        let (x, y) = (a.to_ours(), b.to_ours());
        let (f, g) = (mp(&a), mp(&b));
        prop_assert_eq!(x.cmp(&y), f.partial_cmp(&g).unwrap());
        prop_assert_eq!(x.cmp_abs(&y), f.cmp_abs(&g).unwrap());
        prop_assert_eq!(x == y, f == g);
        prop_assert_eq!(x.cmp(&x), Ordering::Equal);
    }

    #[test]
    fn round_matches_mpfr(a in dyadic_with(12, -300..=300), p in prec(), r in round_mode()) {
        let f = mp(&a);
        check_round(a.to_ours().round(p, r), MpFloat::with_val_round(u32::try_from(p).unwrap(), &f, rug_round(r)))?;
    }

    #[test]
    fn add_sub_match_mpfr(a in dyadic(), b in dyadic(), p in prec(), r in round_mode()) {
        let (x, y) = (a.to_ours(), b.to_ours());
        let (f, g) = (mp(&a), mp(&b));
        let pp = u32::try_from(p).unwrap();
        check_round(x.add_round(&y, p, r), MpFloat::with_val_round(pp, &f + &g, rug_round(r)))?;
        check_round(x.sub_round(&y, p, r), MpFloat::with_val_round(pp, &f - &g, rug_round(r)))?;
    }

    /// The destination and assign forms (CAS-91) agree with the value
    /// forms, whatever the destination held before (zero, short, long
    /// enough to live on the heap, a boxed exponent).
    #[test]
    fn in_place_forms_agree(
        a in dyadic(),
        b in dyadic(),
        old in prop_oneof![dyadic(), dyadic_with(12, -5..=5)],
        huge_old: bool,
        p in prec(),
        r in round_mode(),
    ) {
        type Value = fn(&Float, &Float, u64, Round) -> (Float, Ordering);
        type From = fn(&mut Float, &Float, &Float, u64, Round) -> Ordering;
        type Assign = fn(&mut Float, &Float, u64, Round) -> Ordering;
        let (x, y) = (a.to_ours(), b.to_ours());
        let fresh = || {
            let o = old.to_ours();
            if huge_old { o.mul_2exp(&(Integer::from(i64::MAX) * Integer::from(4))) } else { o }
        };
        let forms: [(&str, Value, From, Assign); 3] = [
            ("add", Float::add_round, Float::add_round_from, Float::add_round_assign),
            ("sub", Float::sub_round, Float::sub_round_from, Float::sub_round_assign),
            ("mul", Float::mul_round, Float::mul_round_from, Float::mul_round_assign),
        ];
        for (name, value, from, assign) in forms {
            let want = value(&x, &y, p, r);
            let mut d = fresh();
            let ord = from(&mut d, &x, &y, p, r);
            prop_assert_eq!(&d, &want.0, "{}_round_from", name);
            prop_assert_eq!(ord, want.1, "{}_round_from ternary", name);
            prop_assert!(d.is_canonical());
            let mut z = x.clone();
            let ord = assign(&mut z, &y, p, r);
            prop_assert_eq!(&z, &want.0, "{}_round_assign", name);
            prop_assert_eq!(ord, want.1, "{}_round_assign ternary", name);
        }
        if !y.is_zero() {
            let want = x.div_round(&y, p, r);
            let mut d = fresh();
            prop_assert_eq!(d.div_round_from(&x, &y, p, r), want.1, "div_round_from ternary");
            prop_assert_eq!(&d, &want.0, "div_round_from");
        }
        let ax = x.abs();
        let want = ax.sqrt_round(p, r);
        let mut d = fresh();
        prop_assert_eq!(d.sqrt_round_from(&ax, p, r), want.1, "sqrt_round_from ternary");
        prop_assert_eq!(&d, &want.0, "sqrt_round_from");
    }

    #[test]
    fn mul_matches_mpfr(
        a in dyadic_with(8, -(1 << 27)..=(1 << 27)),
        b in dyadic_with(8, -(1 << 27)..=(1 << 27)),
        p in prec(),
        r in round_mode(),
    ) {
        let (f, g) = (mp(&a), mp(&b));
        let pp = u32::try_from(p).unwrap();
        check_round(a.to_ours().mul_round(&b.to_ours(), p, r), MpFloat::with_val_round(pp, &f * &g, rug_round(r)))?;
    }

    /// One operand partly or entirely below the other's rounding position.
    #[test]
    fn far_apart_match_mpfr((a, b) in far_apart_pair(1 << 20), p in prec(), r in round_mode()) {
        let (x, y) = (a.to_ours(), b.to_ours());
        let (f, g) = (mp(&a), mp(&b));
        let pp = u32::try_from(p).unwrap();
        check_round(x.add_round(&y, p, r), MpFloat::with_val_round(pp, &f + &g, rug_round(r)))?;
        check_round(x.sub_round(&y, p, r), MpFloat::with_val_round(pp, &f - &g, rug_round(r)))?;
        check_round(y.sub_round(&x, p, r), MpFloat::with_val_round(pp, &g - &f, rug_round(r)))?;
    }

    /// One-limb mantissas at every gap around the `u128` window's edges (63,
    /// 64, 126, 127 bits), where bits start to fall off.
    #[test]
    fn one_limb_window_edges(
        x in strategy::top_limb(),
        y in strategy::top_limb(),
        gap in 0i64..200,
        nx: bool,
        ny: bool,
        p in 1u64..=64,
        r in round_mode(),
    ) {
        let lead = |m: Limb| m | 1 << 63;
        let mk = |m: Limb, neg: bool, e: i64| {
            let m = Integer::from(lead(m));
            Float::from_integer_exp(&if neg { -m } else { m }, &Integer::from(e))
        };
        let (a, b) = (mk(x, nx, 0), mk(y, ny, -gap));
        let (f, g) = (mp(&a.to_rug()), mp(&b.to_rug()));
        let pp = u32::try_from(p).unwrap();
        check_round(a.add_round(&b, p, r), MpFloat::with_val_round(pp, &f + &g, rug_round(r)))?;
        check_round(a.sub_round(&b, p, r), MpFloat::with_val_round(pp, &f - &g, rug_round(r)))?;
        check_round(b.sub_round(&a, p, r), MpFloat::with_val_round(pp, &g - &f, rug_round(r)))?;
        check_round(a.mul_round(&b, p, r), MpFloat::with_val_round(pp, &f * &g, rug_round(r)))?;
    }

    /// Mantissas of one to four limbs in adversarial shapes at gaps around
    /// every limb boundary, against MPFR: the fixed windows of `small`
    /// (CAS-91, CAS-99), where bits start to fall off, a carry shifts the
    /// lowest bit out, a difference cancels, or a quotient is exact.
    #[test]
    fn small_window_edges(
        (kx, ky) in (1usize..=4, 1usize..=4),
        (sx, sy) in (0u8..5, 0u8..5),
        (rx, ry) in (any::<[Limb; 4]>(), any::<[Limb; 4]>()),
        limb_gap in 0i64..=5,
        bit_gap in prop_oneof![-2i64..=2, Just(63), Just(-63)],
        (nx, ny) in (any::<bool>(), any::<bool>()),
        p in prop_oneof![1u64..=256, Just(64), Just(127), Just(128), Just(129), Just(192), Just(256)],
        r in round_mode(),
    ) {
        let (a, b) = (shaped(kx, sx, rx, nx, 0), shaped(ky, sy, ry, ny, -(64 * limb_gap + bit_gap)));
        let (f, g) = (mp(&a.to_rug()), mp(&b.to_rug()));
        let pp = u32::try_from(p).unwrap();
        check_round(a.add_round(&b, p, r), MpFloat::with_val_round(pp, &f + &g, rug_round(r)))?;
        check_round(a.sub_round(&b, p, r), MpFloat::with_val_round(pp, &f - &g, rug_round(r)))?;
        check_round(b.sub_round(&a, p, r), MpFloat::with_val_round(pp, &g - &f, rug_round(r)))?;
        check_round(a.mul_round(&b, p, r), MpFloat::with_val_round(pp, &f * &g, rug_round(r)))?;
        check_round(a.div_round(&b, p, r), MpFloat::with_val_round(pp, &f / &g, rug_round(r)))?;
        let (aa, fa) = (a.abs(), f.clone().abs());
        check_round(aa.sqrt_round(p, r), MpFloat::with_val_round(pp, fa.sqrt_ref(), rug_round(r)))?;
        // Squares are exact roots; their neighbours sit next to a root.
        let sq = aa.mul(&aa);
        check_round(sq.sqrt_round(p, r), MpFloat::with_val_round(pp, mp(&sq.to_rug()).sqrt_ref(), rug_round(r)))?;
    }

    /// Values near powers of two, where a rounding carry or a cancellation
    /// moves the exponent: `2^k ± small`, `2^k - 2^j`.
    #[test]
    fn near_powers_of_two(
        k in -200i64..200,
        j in 1u32..300,
        neg: bool,
        p in prec(),
        r in round_mode(),
    ) {
        let x = Float::from_integer_exp(&Integer::ONE, &Integer::from(k));
        let small = Float::from_integer_exp(&if neg { -Integer::ONE } else { Integer::ONE }, &Integer::from(k - i64::from(j)));
        let (f, g) = (mp(&x.to_rug()), mp(&small.to_rug()));
        let pp = u32::try_from(p).unwrap();
        check_round(x.add_round(&small, p, r), MpFloat::with_val_round(pp, &f + &g, rug_round(r)))?;
        check_round(x.sub_round(&small, p, r), MpFloat::with_val_round(pp, &f - &g, rug_round(r)))?;
        let ones = x.sub(&Float::from_integer_exp(&Integer::ONE, &Integer::from(k - i64::from(j))));
        check_round(ones.round(p, r), MpFloat::with_val_round(pp, &mp(&ones.to_rug()), rug_round(r)))?;
        check_round(ones.mul_round(&ones, p, r), MpFloat::with_val_round(pp, mp(&ones.to_rug()).square_ref(), rug_round(r)))?;
    }

    /// Exponents near and beyond ±2^63: scaling both operands by 2^E scales
    /// the result and keeps the ternary value.
    #[test]
    fn huge_exponents(a in dyadic_with(4, -300..=300), b in dyadic_with(4, -300..=300), k in huge_shift(), p in prec(), r in round_mode()) {
        let k = Signed::from_rug(&k).to_integer();
        let (x, y) = (a.to_ours(), b.to_ours());
        let (xs, ys) = (scaled(&x, &k), scaled(&y, &k));
        let (s, o) = x.add_round(&y, p, r);
        prop_assert_eq!(xs.add_round(&ys, p, r), (scaled(&s, &k), o));
        let (d, o) = x.sub_round(&y, p, r);
        prop_assert_eq!(xs.sub_round(&ys, p, r), (scaled(&d, &k), o));
        let (m, o) = x.mul_round(&y, p, r);
        let k2 = -&k;
        prop_assert_eq!(xs.mul_round(&scaled(&y, &k2), p, r), (m, o));
        prop_assert_eq!(xs.add(&ys), scaled(&x.add(&y), &k));
        prop_assert_eq!(xs.cmp(&ys), x.cmp(&y));
        prop_assert_eq!(xs.round(p, r), {
            let (v, o) = x.round(p, r);
            (scaled(&v, &k), o)
        });
    }

    #[test]
    fn from_rational_matches_mpfr(q in strategy::rational(), p in prec(), r in round_mode()) {
        let q = Rational::new(q.0.to_integer(), q.1.to_integer());
        let rq = q.to_rug();
        let pp = u32::try_from(p).unwrap();
        check_round(Float::from_rational(&q, p, r), MpFloat::with_val_round(pp, &rq, rug_round(r)))?;
    }

    #[test]
    fn to_f64_matches_mpfr(a in prop_oneof![
        dyadic_with(3, -1200..=1200),
        dyadic_with(2, -1080..=-1015),
        dyadic_with(2, 1015..=1030),
    ]) {
        let x = a.to_ours();
        let want = mp(&a).to_f64();
        let got = x.to_f64();
        prop_assert!(got.to_bits() == want.to_bits() || (got == 0.0 && want == 0.0), "{} vs {}", got, want);
    }

    #[test]
    fn from_f64_exact(bits: u64) {
        let v = f64::from_bits(bits);
        match Float::from_f64(v) {
            None => prop_assert!(!v.is_finite()),
            Some(x) => {
                prop_assert_eq!(x.to_rug(), Dyadic::from_rug(&MpFloat::with_val(53, v)).unwrap());
                prop_assert_eq!(x.to_f64(), if v == 0.0 { 0.0 } else { v });
            }
        }
    }

    /// Equal values have equal hashes and compare equal, however built.
    #[test]
    fn structural_eq_is_value_eq(a in dyadic(), pad in 0usize..3, shift in 0u32..64) {
        use std::hash::{BuildHasher, BuildHasherDefault};
        let x = a.to_ours();
        let (m, e) = x.to_integer_exp();
        // The same value with extra zero bits below the mantissa.
        let k = 64 * pad as u64 + u64::from(shift);
        let (neg, mag) = m.into_parts();
        let y = Float::from_raw(neg, (&mag << k).limbs(), e - Integer::from(k));
        prop_assert_eq!(&x, &y);
        let h = BuildHasherDefault::<std::collections::hash_map::DefaultHasher>::default();
        prop_assert_eq!(h.hash_one(&x), h.hash_one(&y));
    }
}

/// `1 ± 2^-10^12` at 64 bits: no exact intermediate (it would take 125 GB, far
/// beyond the test's memory cap), and the right value and ternary in every
/// mode.
#[test]
fn far_apart_add_without_huge_intermediate() {
    let one = Float::from(1);
    let tiny = Float::from_integer_exp(&Integer::ONE, &Integer::from(-1_000_000_000_000i64));
    let ulp = Float::from_integer_exp(&Integer::ONE, &Integer::from(-63));
    for r in Round::ALL {
        let (s, o) = one.add_round(&tiny, 64, r);
        let up = matches!(r, Round::Up | Round::Ceil);
        assert_eq!(s, if up { one.add(&ulp) } else { one.clone() }, "{r:?}");
        assert_eq!(
            o,
            if up {
                Ordering::Greater
            } else {
                Ordering::Less
            },
            "{r:?}"
        );
        let (d, o) = one.sub_round(&tiny, 64, r);
        let down = matches!(r, Round::Down | Round::Floor);
        // 1 - 2^-64 is the representable value just below 1.
        let below = one.sub(&Float::from_integer_exp(&Integer::ONE, &Integer::from(-64)));
        assert_eq!(d, if down { below } else { one.clone() }, "{r:?}");
        assert_eq!(
            o,
            if down {
                Ordering::Less
            } else {
                Ordering::Greater
            },
            "{r:?}"
        );
    }
}

/// Exponents straddling the `i64` edge of `Integer`'s inline form.
#[test]
fn exponent_at_i64_edge() {
    let edge = Integer::from(i64::MAX);
    let x = Float::from_integer_exp(&Integer::from(3), &(&edge - &Integer::from(2)));
    let y = x.mul_2exp(&Integer::from(5));
    assert!(y.exponent().to_i64().is_none());
    assert_eq!(y.mul_2exp(&Integer::from(-5)), x);
    let (s, o) = y.add_round(&y, 2, Round::Nearest);
    assert_eq!(o, Ordering::Equal);
    assert_eq!(s.exponent(), &(y.exponent() + &Integer::ONE));
    let z = Float::from_integer_exp(&Integer::ONE, &(-&edge - &Integer::from(10)));
    let (p, _) = z.mul_round(&z, 10, Round::Nearest);
    assert!(*p.exponent() < &(-&edge * &Integer::from(2)) + &Integer::from(10));
    assert_eq!(z.to_f64(), 0.0);
    assert_eq!(y.to_f64(), f64::INFINITY);
}

/// The carry case of the fixed window whose only sticky information is the
/// bit the carry shifts out: `x = 0.1…1 1000…0`, `y = 0.1000…0 0…01` 64
/// bits lower, at 128 bits the sum is inexact by `y`'s lowest bit alone.
#[test]
fn small_window_carry_out() {
    let x = shaped(2, 3, [0; 4], false, 0);
    let y = shaped(2, 1, [0; 4], false, -64);
    let (f, g) = (mp(&x.to_rug()), mp(&y.to_rug()));
    for r in Round::ALL {
        // `&f + &g` stays exact until rounded once, at 128 bits.
        let want = MpFloat::with_val_round(128, &f + &g, rug_round(r));
        check_round(x.add_round(&y, 128, r), want).unwrap();
    }
    assert_ne!(x.add_round(&y, 128, Round::Down).1, Ordering::Equal);
}

#[test]
fn small_cases() {
    assert_eq!(
        Float::ZERO.round(1, Round::Up),
        (Float::ZERO, Ordering::Equal)
    );
    assert_eq!(Float::from(3).bit_length(), 2);
    assert_eq!(Float::from(3).mant().len(), 1);
    assert_eq!(
        Float::from(-6).to_rational(),
        Rational::from(Integer::from(-6))
    );
    let (x, o) = Float::from(3).round(1, Round::Nearest);
    assert_eq!((x, o), (Float::from(4), Ordering::Greater));
    let (x, o) = Float::from(5).round(2, Round::Nearest);
    assert_eq!((x, o), (Float::from(4), Ordering::Less));
    assert!(Float::from(-1) < Float::ZERO && Float::ZERO < Float::from(1));
    assert_eq!(format!("{:?}", Float::from(-3)), "-0x0.c000000000000000p2");
    check_same(
        Float::from(7),
        Dyadic::new(&rug::Integer::from(7), &rug::Integer::new()),
    )
    .unwrap();
}

/// `|true - approx|` checked exactly: `approx` is a lower bound of the
/// magnitude within `APPROX_ULPS` units of `2^(exp(approx) - prec)`.
fn check_approx_quotient(q: &Float, exact: &Rational, prec: u64) -> Result<(), TestCaseError> {
    prop_assert!(q.bit_length() <= prec, "more than prec bits");
    let qv = q.to_rational();
    prop_assert!(qv.abs() <= exact.abs(), "not a lower bound: {:?}", q);
    prop_assert!(q.is_zero() || q.sign() == exact.sign());
    if !q.is_zero() {
        let ulp = Float::from_integer_exp(&Integer::ONE, &(q.exponent() - &Integer::from(prec)))
            .to_rational();
        let bound = &ulp * &Rational::from(Integer::from(super::APPROX_ULPS));
        prop_assert!(
            &(exact.abs()) - &qv.abs() < bound,
            "beyond the bound: {:?}",
            q
        );
    }
    Ok(())
}

/// For a square root `s` of `a`: `s^2 <= a < (s + u)^2` with
/// `u = APPROX_ULPS × 2^(exp(s) - prec)`, and similarly `y` for `1 / sqrt(a)`
/// (`y^2 a <= 1 < (y + u)^2 a`).
fn check_approx_root(
    r: &Float,
    a: &Rational,
    prec: u64,
    inverse: bool,
) -> Result<(), TestCaseError> {
    prop_assert!(r.bit_length() <= prec && !r.is_negative());
    let rv = r.to_rational();
    let ulp = Float::from_integer_exp(&Integer::ONE, &(r.exponent() - &Integer::from(prec)))
        .to_rational();
    let hi = &rv + &(&ulp * &Rational::from(Integer::from(super::APPROX_ULPS)));
    let sq = |x: &Rational| x * x;
    if inverse {
        prop_assert!(&sq(&rv) * a <= Rational::ONE, "not a lower bound");
        prop_assert!(&sq(&hi) * a > Rational::ONE, "beyond the bound");
    } else {
        prop_assert!(sq(&rv) <= *a, "not a lower bound");
        prop_assert!(sq(&hi) > *a, "beyond the bound");
    }
    Ok(())
}

proptest! {
    #[test]
    fn div_matches_mpfr(
        a in dyadic_with(8, -(1 << 27)..=(1 << 27)),
        b in dyadic_with(8, -(1 << 27)..=(1 << 27)),
        p in prec(),
        r in round_mode(),
    ) {
        let (x, y) = (a.to_ours(), b.to_ours());
        match x.checked_div_round(&y, p, r) {
            None => prop_assert!(y.is_zero()),
            Some(ours) => {
                let pp = u32::try_from(p).unwrap();
                check_round(ours, div_round(&a, &b, pp, rug_round(r)))?;
            }
        }
    }

    /// Exact quotients and quotients exactly halfway between two `p`-bit
    /// values: `a = b h` for an odd `h` of `p + 1` bits (or of at most `p`).
    #[test]
    fn div_exact_and_halfway(
        b in nonzero(4, -300..=300),
        h in strategy::top_limb(),
        p in 1u64..=63,
        exact: bool,
        r in round_mode(),
    ) {
        let y = b.to_ours();
        let bits = if exact { p } else { p + 1 };
        let h = (h >> (64 - bits)) | 1 | (1 << (bits - 1));
        let hf = Float::from_integer_exp(&Integer::from(h), &Integer::from(-7));
        let x = y.mul(&hf);
        let pp = u32::try_from(p).unwrap();
        let ours = x.div_round(&y, p, r);
        if exact {
            prop_assert_eq!(ours.1, Ordering::Equal);
        }
        check_round(ours, div_round(&x.to_rug(), &b, pp, rug_round(r)))?;
    }

    #[test]
    fn sqrt_matches_mpfr(a in dyadic_with(8, -(1 << 27)..=(1 << 27)), p in prec(), r in round_mode()) {
        let x = a.to_ours();
        match x.checked_sqrt_round(p, r) {
            None => prop_assert!(x.is_negative()),
            Some(ours) => {
                let pp = u32::try_from(p).unwrap();
                check_round(ours, MpFloat::with_val_round(pp, mp(&a).sqrt_ref(), rug_round(r)))?;
            }
        }
    }

    /// Perfect squares, including roots of `p + 1` bits (halfway at `p`).
    #[test]
    fn sqrt_exact_and_halfway(
        s in strategy::top_limb(),
        p in 1u64..=63,
        exact: bool,
        e in -200i64..200,
        r in round_mode(),
    ) {
        let bits = if exact { p } else { p + 1 };
        let s = (s >> (64 - bits)) | 1 | (1 << (bits - 1));
        let root = Float::from_integer_exp(&Integer::from(s), &Integer::from(e));
        let x = root.mul(&root);
        let pp = u32::try_from(p).unwrap();
        let ours = x.sqrt_round(p, r);
        if exact {
            prop_assert_eq!(&ours, &(root.clone(), Ordering::Equal));
        }
        check_round(ours, MpFloat::with_val_round(pp, mp(&x.to_rug()).sqrt_ref(), rug_round(r)))?;
    }

    #[test]
    fn div_approx_within_bound(a in dyadic_with(8, -500..=500), b in dyadic_with(8, -500..=500), p in prec()) {
        let (x, y) = (a.to_ours(), b.to_ours());
        match x.div_approx(&y, p) {
            None => prop_assert!(y.is_zero()),
            Some(q) => check_approx_quotient(&q, &(&x.to_rational() / &y.to_rational()), p)?,
        }
    }

    #[test]
    fn sqrt_approx_within_bound(a in dyadic_with(8, -500..=500), p in prec()) {
        let x = a.to_ours();
        match (x.sqrt_approx(p), x.rsqrt_approx(p)) {
            (None, None) => prop_assert!(x.is_negative()),
            (Some(s), None) => {
                prop_assert!(x.is_zero() && s.is_zero());
            }
            (Some(s), Some(y)) => {
                let av = x.to_rational();
                check_approx_root(&s, &av, p, false)?;
                check_approx_root(&y, &av, p, true)?;
            }
            (None, Some(_)) => prop_assert!(false, "rsqrt of a negative"),
        }
    }

    /// Exponents near and beyond ±2^63: scaling commutes with division and
    /// (by even powers) with the square root.
    #[test]
    fn div_sqrt_huge_exponents(a in dyadic_with(4, -300..=300), b in nonzero(4, -300..=300), k in huge_shift(), p in prec(), r in round_mode()) {
        let k = Signed::from_rug(&k).to_integer();
        let (x, y) = (a.to_ours().abs(), b.to_ours());
        let (q, o) = x.div_round(&y, p, r);
        prop_assert_eq!(x.mul_2exp(&k).div_round(&y, p, r), (q.mul_2exp(&k), o));
        let k2 = &k * &Integer::from(2);
        let (s, o) = x.sqrt_round(p, r);
        prop_assert_eq!(x.mul_2exp(&k2).sqrt_round(p, r), (s.mul_2exp(&k), o));
        let k1 = &k2 + &Integer::ONE;
        let lhs = x.mul_2exp(&k1).sqrt_round(p, r);
        let rhs = x.mul_2exp(&Integer::ONE).sqrt_round(p, r);
        prop_assert_eq!(lhs, (rhs.0.mul_2exp(&k), rhs.1));
        if let Some(q) = x.div_approx(&y, p) {
            prop_assert_eq!(x.mul_2exp(&k).div_approx(&y, p), Some(q.mul_2exp(&k)));
        }
    }
}

/// The fixed-size square root just above `i64::MIN`, where the exponent
/// still takes it but `e - 128 N` does not fit: the general path answers
/// (CAS-113: the subtraction overflowed), and scaling by an even power of
/// two still commutes with the root.
#[test]
fn sqrt_next_to_the_smallest_i64_exponent() {
    for x in [Float::from(3), Float::from(6), Float::from(7)] {
        for d in [0i64, 1, 2, 3, 255, 256, 257, 511, 512, 513, 641, 700] {
            // An even shift taking x's exponent to about i64::MIN + d.
            let e = i128::try_from(x.exponent()).expect("small");
            let s = (i128::from(i64::MIN) + i128::from(d) - e) & !1;
            let (s, half) = (Integer::from(s), Integer::from(s / 2));
            for p in [10, 64, 100, 150, 200, 250, 300] {
                for r in [Round::Down, Round::Up, Round::Nearest] {
                    let (q, o) = x.sqrt_round(p, r);
                    assert_eq!(
                        x.mul_2exp(&s).sqrt_round(p, r),
                        (q.mul_2exp(&half), o),
                        "{x:?} 2^{s} at {p}"
                    );
                }
            }
        }
    }
}

#[test]
fn div_sqrt_edge_cases() {
    assert_eq!(
        Float::from(1).checked_div_round(&Float::ZERO, 10, Round::Nearest),
        None
    );
    assert_eq!(Float::from(1).div_approx(&Float::ZERO, 10), None);
    assert_eq!(Float::from(-4).checked_sqrt_round(10, Round::Nearest), None);
    assert_eq!(Float::from(-4).sqrt_approx(10), None);
    assert_eq!(Float::ZERO.rsqrt_approx(10), None);
    assert_eq!(
        Float::ZERO.sqrt_round(5, Round::Up),
        (Float::ZERO, Ordering::Equal)
    );
    assert_eq!(
        Float::from(9).sqrt_round(2, Round::Down),
        (Float::from(3), Ordering::Equal)
    );
    assert_eq!(
        Float::from(6).div_round(&Float::from(3), 1, Round::Down),
        (Float::from(2), Ordering::Equal)
    );
    assert!(
        std::panic::catch_unwind(|| Float::from(1).div_round(&Float::ZERO, 10, Round::Down))
            .is_err()
    );
}
