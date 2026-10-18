//! `Ball` containment (CAS-DOC-5 §11): the result of every operation must
//! contain the exact result at points sampled from the inputs (endpoints,
//! midpoint, interior), and the MPFR bracket of the exact result on exact
//! inputs; results at 64, 128 and 1000 bits must overlap; exact results must
//! stay exact; and the predicates must agree with exact rational reasoning.
#![allow(
    clippy::many_single_char_names,
    clippy::needless_pass_by_value,
    reason = "balls a, b, points x, y, values p, q; results passed as built"
)]

use std::cmp::Ordering;

use proptest::prelude::*;
use rug::Float as MpFloat;

use super::{Ball, RealBalls};
use crate::algebra::Ring;
use crate::ball::{Float, Mag};
use crate::testkit::float::{GUARD, bracket, check_encloses, dyadic_with, huge_shift};
use crate::testkit::{Signed, ToRug};
use crate::{Integer, Rational};

/// Precisions: small, limb edges, random up to 1000 bits.
fn prec() -> BoxedStrategy<u64> {
    prop_oneof![
        2 => proptest::sample::select(&[2u64, 10, 53, 64, 65, 128, 200][..]),
        1 => 1..1000u64,
    ]
    .boxed()
}

/// A ball: dyadic midpoint, radius 0 (often), relative to the midpoint, or
/// absolute.
fn ball() -> BoxedStrategy<Ball> {
    let mid = dyadic_with(4, -200..=200).prop_map(|d| d.to_ours());
    let rel = (0i64..120).prop_map(Integer::from);
    let abs = dyadic_with(1, -250..=150).prop_map(|d| Mag::from_float(&d.to_ours()));
    prop_oneof![
        2 => mid.clone().prop_map(Ball::exact),
        2 => (mid.clone(), rel).prop_map(|(m, k)| {
            let r = Mag::from_float(&m).mul_2exp(&-k);
            Ball::new(m, r)
        }),
        1 => (mid, abs).prop_map(|(m, r)| Ball::new(m, r)),
    ]
    .boxed()
}

fn rat(x: &Float) -> Rational {
    x.to_rational()
}

fn rad(b: &Ball) -> Rational {
    b.rad().to_float().expect("finite").to_rational()
}

/// Points of a finite ball: midpoint, both endpoints and two interior points.
fn points(b: &Ball) -> Vec<Rational> {
    let (m, r) = (rat(b.mid()), rad(b));
    let third = Rational::new(Integer::ONE, Integer::from(3));
    let other = Rational::new(Integer::from(-5), Integer::from(7));
    vec![
        m.clone(),
        &m + &r,
        &m - &r,
        &m + &(&r * &third),
        &m + &(&r * &other),
    ]
}

/// Fails unless `result` contains `exact`.
fn check_contains(what: &str, result: &Ball, exact: &Rational) -> Result<(), TestCaseError> {
    prop_assert!(
        result.contains_rational(exact),
        "{}: {:?} misses {}",
        what,
        result,
        exact
    );
    Ok(())
}

/// The MPFR oracle on exact inputs: `result` against the bracket of the true
/// value computed at `prec + GUARD` bits.
fn check_oracle(
    result: &Ball,
    prec: u64,
    f: impl Fn(u32, rug::float::Round) -> MpFloat,
) -> Result<(), TestCaseError> {
    if !result.is_finite() {
        return Ok(());
    }
    let b = bracket(u32::try_from(prec).unwrap(), f);
    check_encloses(&rat(result.mid()).to_rug(), &rad(result).to_rug(), &b)
}

fn mp(x: &Float) -> MpFloat {
    x.to_rug().to_rug().expect("in MPFR's range")
}

proptest! {
    #[test]
    fn arithmetic_contains_sampled_points(a in ball(), b in ball(), p in prec()) {
        let (sum, diff, prod, quot) = (a.add(&b, p), a.sub(&b, p), a.mul(&b, p), a.div(&b, p));
        for x in points(&a) {
            for y in points(&b) {
                check_contains("add", &sum, &(&x + &y))?;
                check_contains("sub", &diff, &(&x - &y))?;
                check_contains("mul", &prod, &(&x * &y))?;
                if !y.is_zero() {
                    check_contains("div", &quot, &(&x / &y))?;
                }
            }
        }
    }

    /// Square roots of sampled points: the result must meet the MPFR bracket
    /// of `sqrt(x)` at 2048 bits, far narrower than any radius here.
    #[test]
    fn sqrt_contains_sampled_points(a in ball(), p in prec()) {
        let s = a.sqrt(p);
        if !s.is_finite() {
            // Only for a ball reaching below 0.
            let below = &rat(a.mid()) - &rad(&a);
            prop_assert!(below.is_negative(), "indeterminate sqrt of {:?}", a);
            return Ok(());
        }
        for x in points(&a) {
            prop_assert!(!x.is_negative());
            let xr = MpFloat::with_val(4096, &x.to_rug());
            let b = bracket(2048 - GUARD, |q, r| MpFloat::with_val_round(q, xr.sqrt_ref(), r).0);
            check_encloses(&rat(s.mid()).to_rug(), &rad(&s).to_rug(), &b)?;
        }
    }

    #[test]
    fn oracle_on_exact_inputs(a in dyadic_with(4, -200..=200), b in dyadic_with(4, -200..=200), p in prec()) {
        let (x, y) = (Ball::exact(a.to_ours()), Ball::exact(b.to_ours()));
        let (f, g) = (mp(x.mid()), mp(y.mid()));
        check_oracle(&x.add(&y, p), p, |q, r| MpFloat::with_val_round(q, &f + &g, r).0)?;
        check_oracle(&x.sub(&y, p), p, |q, r| MpFloat::with_val_round(q, &f - &g, r).0)?;
        check_oracle(&x.mul(&y, p), p, |q, r| MpFloat::with_val_round(q, &f * &g, r).0)?;
        if !g.is_zero() {
            check_oracle(&x.div(&y, p), p, |q, r| MpFloat::with_val_round(q, &f / &g, r).0)?;
        }
        if !f.is_sign_negative() {
            check_oracle(&x.sqrt(p), p, |q, r| MpFloat::with_val_round(q, f.sqrt_ref(), r).0)?;
        }
    }

    /// The same expression at 64, 128 and 1000 bits gives pairwise
    /// overlapping balls: a wrong error bound shows up here without an oracle.
    #[test]
    fn cross_precision_overlap(a in ball(), b in ball(), n in -6i64..6) {
        let expr = |p: u64| {
            let s = a.mul(&b, p).add(&a, p);
            let q = s.div(&b.add(&Ball::from_integer(&Integer::from(3)), p), p);
            (q.sqrt(p), a.pow(n, p), s)
        };
        let rs: Vec<_> = [64u64, 128, 1000].into_iter().map(expr).collect();
        for i in 0..3 {
            for j in i + 1..3 {
                prop_assert!(rs[i].0.overlaps(&rs[j].0), "sqrt quotient: {:?} vs {:?}", rs[i].0, rs[j].0);
                prop_assert!(rs[i].1.overlaps(&rs[j].1), "pow: {:?} vs {:?}", rs[i].1, rs[j].1);
                prop_assert!(rs[i].2.overlaps(&rs[j].2), "mul add: {:?} vs {:?}", rs[i].2, rs[j].2);
            }
        }
    }

    /// Exact inputs whose results fit the precision give exact results.
    #[test]
    fn exact_results_stay_exact(x in -1_000_000i64..1_000_000, y in -1_000_000i64..1_000_000, k in 0i64..4) {
        let (a, b) = (Ball::from_integer(&Integer::from(x)), Ball::from_integer(&Integer::from(y)));
        for p in [64u64, 200] {
            prop_assert!(a.add(&b, p).is_exact() && a.sub(&b, p).is_exact() && a.mul(&b, p).is_exact());
            prop_assert!(a.pow(k, p).is_exact());
            let prod = a.mul(&b, p);
            if y != 0 {
                let q = prod.div(&b, p);
                prop_assert!(q.is_exact() && q == a, "{:?} / {:?} = {:?}", prod, b, q);
            }
            let sq = a.mul(&a, p).sqrt(p);
            prop_assert!(sq.is_exact() && sq == a.abs());
        }
    }

    /// The predicates against exact rational reasoning on the endpoints.
    #[test]
    fn predicates(a in ball(), b in ball(), q in crate::testkit::strategy::rational()) {
        let q = Rational::new(q.0.to_integer(), q.1.to_integer());
        let (m, r, n, s) = (rat(a.mid()), rad(&a), rat(b.mid()), rad(&b));
        let d = (&m - &n).abs();
        prop_assert_eq!(a.contains_rational(&q), (&q - &m).abs() <= r);
        prop_assert_eq!(a.overlaps(&b), d <= &r + &s);
        prop_assert_eq!(a.contains_ball(&b), &d + &s <= r);
        for x in points(&a) {
            prop_assert!(a.contains_rational(&x));
            if let Some(sg) = a.sign() {
                prop_assert_eq!(x.sign(), sg);
            }
            for y in points(&b) {
                if let Some(o) = a.cmp(&b) {
                    prop_assert_eq!(x.cmp(&y), o);
                }
            }
        }
        // sign decides whenever 0 is outside the ball.
        prop_assert_eq!(a.sign().is_some(), r < m.abs() || (r.is_zero() && m.is_zero()));
    }

    #[test]
    fn union_and_intersection(a in ball(), b in ball(), p in prec()) {
        let u = a.union(&b, p);
        prop_assert!(u.contains_ball(&a) && u.contains_ball(&b), "union {:?}", u);
        match a.intersection(&b, p) {
            None => prop_assert!(!a.overlaps(&b)),
            Some(i) => {
                for x in points(&a).into_iter().chain(points(&b)) {
                    if a.contains_rational(&x) && b.contains_rational(&x) {
                        prop_assert!(i.contains_rational(&x), "intersection {:?} misses {}", i, x);
                    }
                }
            }
        }
    }

    #[test]
    fn rounding_keeps_the_enclosure(a in ball(), bits in 1u64..200) {
        let r = a.round_to_accuracy(bits);
        prop_assert!(r.contains_ball(&a), "{:?} rounded to {:?}", a, r);
        prop_assert!(r.mid().bit_length() <= bits);
        if a.is_exact() && a.mid().bit_length() <= bits {
            prop_assert_eq!(&r, &a);
        }
    }

    /// Exponents near and beyond ±2^63: scaling both inputs scales sums,
    /// differences and square roots (by even powers), and cancels in
    /// quotients.
    #[test]
    fn huge_exponents(a in ball(), b in ball(), k in huge_shift(), p in prec()) {
        let k = Signed::from_rug(&k).to_integer();
        let (ak, bk) = (a.mul_2exp(&k), b.mul_2exp(&k));
        prop_assert_eq!(ak.add(&bk, p), a.add(&b, p).mul_2exp(&k));
        prop_assert_eq!(ak.sub(&bk, p), a.sub(&b, p).mul_2exp(&k));
        prop_assert_eq!(ak.div(&bk, p), a.div(&b, p));
        prop_assert_eq!(ak.mul(&b.mul_2exp(&-&k), p), a.mul(&b, p));
        let k2 = &k * &Integer::from(2);
        prop_assert_eq!(a.mul_2exp(&k2).sqrt(p), a.sqrt(p).mul_2exp(&k));
        prop_assert_eq!(ak.overlaps(&bk), a.overlaps(&b));
        prop_assert_eq!(ak.cmp(&bk), a.cmp(&b));
    }

    /// The ring laws hold exactly on exact integer balls at a precision that
    /// keeps every result exact.
    #[test]
    fn ring_laws_on_exact_balls(x in -10_000i64..10_000, y in -10_000i64..10_000, z in -10_000i64..10_000) {
        let rb = RealBalls::new(4096);
        let [a, b, c] = [x, y, z].map(|v| Ball::from_integer(&Integer::from(v)));
        prop_assert_eq!(rb.add(&rb.add(&a, &b), &c), rb.add(&a, &rb.add(&b, &c)));
        prop_assert_eq!(rb.mul(&rb.mul(&a, &b), &c), rb.mul(&a, &rb.mul(&b, &c)));
        prop_assert_eq!(rb.mul(&a, &rb.add(&b, &c)), rb.add(&rb.mul(&a, &b), &rb.mul(&a, &c)));
        prop_assert_eq!(rb.sub(&a, &b), rb.add(&a, &rb.neg(&b)));
        prop_assert_eq!(rb.pow(&a, 5), rb.mul(&rb.square(&rb.square(&a)), &a));
        prop_assert!(rb.is_zero(&rb.sub(&a, &a)));
        prop_assert_eq!(rb.elem_integer(&Integer::from(x)), a);
    }
}

#[test]
fn indeterminate_cases() {
    let p = 64;
    let zero_ish = Ball::new(Float::from(1), Mag::from_integer(&Integer::from(2)));
    assert_eq!(
        Ball::from_integer(&Integer::ONE).div(&zero_ish, p),
        Ball::INDETERMINATE
    );
    assert_eq!(
        Ball::from_integer(&Integer::ONE).div(&Ball::ZERO, p),
        Ball::INDETERMINATE
    );
    assert_eq!(zero_ish.sqrt(p), Ball::INDETERMINATE);
    assert_eq!(
        Ball::from_integer(&Integer::from(-4)).sqrt(p),
        Ball::INDETERMINATE
    );
    assert_eq!(Ball::ZERO.sqrt(p), Ball::ZERO);
    // A ball touching 0 from above: sqrt([1 ± 1]) is defined.
    let touching = Ball::new(Float::from(1), Mag::ONE);
    assert!(touching.sqrt(p).is_finite());
    assert!(touching.sqrt(p).contains(&Float::ZERO));
    let ind = Ball::INDETERMINATE;
    assert_eq!(ind.add(&touching, p), ind);
    assert_eq!(ind.sign(), None);
    assert_eq!(ind.cmp(&touching), None);
    assert!(ind.contains(&Float::from(12345)));
    assert!(
        ind.overlaps(&touching) && ind.contains_ball(&touching) && !touching.contains_ball(&ind)
    );
    assert_eq!(ind.rel_accuracy_bits(), i64::MIN);
    assert_eq!(
        Ball::from_integer(&Integer::from(5)).rel_accuracy_bits(),
        i64::MAX
    );
    assert_eq!(Ball::ZERO.sign(), Some(Ordering::Equal));
    assert_eq!(touching.sign(), None);
}

#[test]
fn from_rational_encloses() {
    let third = Rational::new(Integer::ONE, Integer::from(3));
    for p in [1u64, 2, 53, 64, 300] {
        let b = Ball::from_rational(&third, p);
        assert!(b.contains_rational(&third) && !b.is_exact(), "{p}: {b:?}");
        assert!(b.rel_accuracy_bits() >= i64::try_from(p).unwrap());
    }
    let half = Rational::new(Integer::ONE, Integer::from(2));
    assert!(Ball::from_rational(&half, 1).is_exact());
}

/// The tightness the formulas promise in easy cases: a few ulps of radius
/// for one operation on exact inputs, and sqrt of a narrow ball.
#[test]
fn radii_are_tight() {
    let p = 100;
    let x = Ball::from_rational(&Rational::new(Integer::from(2), Integer::from(3)), p);
    let y = Ball::from_rational(&Rational::new(Integer::from(5), Integer::from(7)), p);
    for r in [x.add(&y, p), x.mul(&y, p), x.div(&y, p), x.sqrt(p)] {
        assert!(r.rel_accuracy_bits() >= 96, "{r:?}");
    }
}
