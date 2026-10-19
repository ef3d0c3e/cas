//! Values with enormous exponents (CAS-93): representable and cheap, but
//! the operations that need them exactly must refuse instead of building
//! mantissas as long as an exponent gap, and the predicates must stay exact
//! without such sums. Every test here would exhaust memory if they didn't.
#![allow(
    clippy::many_single_char_names,
    reason = "e, a, s, c, x as in the comments"
)]

use std::cmp::Ordering;

use crate::ball::{Ball, EXACT_BITS_MAX, Float, Mag};
use crate::{Integer, Rational};

/// `2^e`.
fn p2(e: &Integer) -> Float {
    Float::from(1i64).mul_2exp(e)
}

fn big_e() -> Integer {
    Integer::from(10u64).pow(30)
}

#[test]
fn checked_forms_refuse() {
    let big = p2(&big_e());
    let one = Float::from(1i64);
    assert_eq!(big.checked_add(&one), None);
    assert_eq!(big.checked_sub(&one), None);
    assert_eq!(one.checked_sub(&big), None);
    assert_eq!(big.checked_to_rational(), None);
    assert_eq!(p2(&-big_e()).checked_to_rational(), None);
    // Cheap ones still work: products, scalings, rounded sums, comparisons.
    assert_eq!(
        big.mul(&big).exponent(),
        &(&(&big_e() * &Integer::from(2)) + &Integer::ONE)
    );
    assert!(big.add_round(&one, 64, crate::ball::Round::Up).0 > big);
    assert_eq!(big.cmp(&one), Ordering::Greater);
    // Just inside the limit.
    let edge = p2(&Integer::from(EXACT_BITS_MAX - 2));
    assert!(
        edge.checked_add(&one)
            .is_some_and(|s| s.bit_length() == EXACT_BITS_MAX - 1)
    );
    assert!(
        p2(&Integer::from(EXACT_BITS_MAX + 2))
            .checked_add(&one)
            .is_none()
    );
}

#[test]
#[should_panic(expected = "Float::add: exact sum beyond EXACT_BITS_MAX bits")]
fn exact_add_panics_past_the_limit() {
    let _ = p2(&big_e()).add(&Float::from(1i64));
}

#[test]
#[should_panic(expected = "Float::to_rational: value beyond EXACT_BITS_MAX bits")]
fn rational_panics_past_the_limit() {
    let _ = Rational::from(&p2(&big_e()));
}

#[test]
fn predicates_stay_exact_across_huge_gaps() {
    let e = big_e();
    let big = p2(&e);
    let r = Mag::power_of_two(&e);
    let r_val = r.to_float().unwrap();
    assert_eq!(r_val, big);
    // [2^e ± 2^e] holds exactly [0, 2^(e+1)].
    let b = Ball::new(big.clone(), r.clone());
    let eps = p2(&Integer::from(-5));
    assert!(b.contains(&Float::ZERO));
    assert!(b.contains(&eps));
    assert!(!b.contains(&eps.neg()));
    assert!(b.contains(&big.mul_2exp(&Integer::ONE)));
    assert!(
        !b.contains(
            &big.mul_2exp(&Integer::ONE)
                .add_round(&eps, 64, crate::ball::Round::Up)
                .0
        )
    );
    // Another ball touching or missing 0 from below.
    let s = Mag::power_of_two(&Integer::from(-6));
    let touch = Ball::new(Float::from_f64(-1.0 / 64.0).unwrap(), s.clone());
    let miss = Ball::new(Float::from_f64(-1.0 / 32.0).unwrap(), s.clone());
    assert!(b.overlaps(&touch) && touch.overlaps(&b));
    assert!(!b.overlaps(&miss) && !miss.overlaps(&b));
    assert_eq!(b.cmp(&miss), Some(Ordering::Greater));
    assert_eq!(b.cmp(&touch), None);
    // contains_ball: [2^-6 ± 2^-6] = [0, 2^-5] lies inside; shifted down it doesn't.
    let inside = Ball::new(Float::from_f64(1.0 / 64.0).unwrap(), s.clone());
    let across = Ball::new(Float::from_f64(1.0 / 128.0).unwrap(), s);
    assert!(b.contains_ball(&inside));
    assert!(!b.contains_ball(&across));
    // A tiny ball against a huge one, the other way round.
    // [2^-e ± 2^-(e-1)] = [-2^-e, 3 · 2^-e].
    let tiny = Ball::new(p2(&-e.clone()), Mag::power_of_two(&-(&e - &Integer::ONE)));
    assert!(tiny.contains(&Float::ZERO) && !tiny.contains(&big));
    assert!(!tiny.overlaps(&Ball::exact(Float::from(1i64))));
}

#[test]
fn functions_of_tiny_and_huge_arguments() {
    let e = big_e();
    let tiny = Ball::exact(p2(&-e.clone()));
    // atanh, sin_pi, cos_pi of 2^-(10^30): no exact 1 - x or x mod 2 shift.
    let a = tiny.atanh(64);
    assert!(
        a.contains(tiny.mid()) && a.rel_accuracy_bits() >= 63,
        "{a:?}"
    );
    let s = tiny.sin_pi(64);
    // π 2^-e has exponent 2 - e.
    assert!(
        s.mid().exponent() == &(Integer::from(2) - &e) && s.rel_accuracy_bits() >= 60,
        "{s:?}"
    );
    let c = tiny.cos_pi(64);
    assert!(c.contains(&Float::from(1i64)) && c.rel_accuracy_bits() >= 60);
    // The other functions at 2^-(10^30): relative accuracy, no panic.
    for (name, f) in [
        ("exp", Ball::exp as fn(&Ball, u64) -> Ball),
        ("expm1", Ball::expm1),
        ("log1p", Ball::log1p),
        ("sin", Ball::sin),
        ("cos", Ball::cos),
        ("atan", Ball::atan),
        ("sinh", Ball::sinh),
        ("tanh", Ball::tanh),
        ("asinh", Ball::asinh),
    ] {
        let v = f(&tiny, 64);
        assert!(
            v.is_finite() && v.rel_accuracy_bits() >= 60,
            "{name}: {v:?}"
        );
    }
    // exp(10^30) and its log, as before.
    let x = Ball::from_integer(&e);
    let ex = x.exp(64);
    assert!(ex.is_finite() && ex.log(160).contains_integer(&e));
    // Printing falls back to hex beyond i64 decimal exponents.
    assert!(
        ex.to_decimal(5, crate::ball::DecimalStyle::Digits)
            .starts_with("[0x1.")
    );
}
