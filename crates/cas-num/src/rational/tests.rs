//! `Rational` checked against GMP's `mpq`. Every result goes through `ToRug`,
//! which also asserts the canonical form (reduced, positive denominator).

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use proptest::prelude::*;
use rug::Rational as Gmpq;

use super::Rational;
use crate::testkit::strategy::{self, rational, rational_pair};
use crate::testkit::{Signed, ToRug, all_forms, differential};
use crate::{Integer, ParseErrorKind};

type Frac = (Signed, Signed);

fn ours((n, d): &Frac) -> Rational {
    Rational::new(n.to_integer(), d.to_integer())
}

fn gmp((n, d): &Frac) -> Gmpq {
    Gmpq::from((n.to_rug(), d.to_rug()))
}

fn hash_of(x: &Rational) -> u64 {
    let mut h = DefaultHasher::new();
    x.hash(&mut h);
    h.finish()
}

/// The integer nearest `x`, halves to even, computed from GMP's floor.
fn gmp_round_ties_even(x: &Gmpq) -> rug::Integer {
    let f = x.clone().floor();
    let diff = Gmpq::from(x - &f);
    let half = Gmpq::from((1, 2));
    let f = f.numer().clone();
    let up = match diff.cmp(&half) {
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Equal => f.is_odd(),
    };
    if up { f + 1 } else { f }
}

fn radix() -> impl Strategy<Value = u32> {
    prop_oneof![2 => Just(10u32), 3 => 2..=36u32]
}

differential! {
    /// Any numerator and denominator, unreduced and with any signs.
    fn new_reduces(n in strategy::integer(), d in strategy::integer()) {
        ours: Rational::checked_new(n.to_integer(), d.to_integer()),
        gmp: (!d.mag().is_zero()).then(|| Gmpq::from((n.to_rug(), d.to_rug()))),
    }

    fn add((a, b) in rational_pair()) {
        ours: all_forms!(ours(&a), +, +=, ours(&b)),
        gmp: gmp(&a) + gmp(&b),
    }

    fn sub((a, b) in rational_pair()) {
        ours: all_forms!(ours(&a), -, -=, ours(&b)),
        gmp: gmp(&a) - gmp(&b),
    }

    fn mul((a, b) in rational_pair()) {
        ours: all_forms!(ours(&a), *, *=, ours(&b)),
        gmp: gmp(&a) * gmp(&b),
    }

    fn div((a, b) in rational_pair()) {
        ours: {
            let (x, y) = (ours(&a), ours(&b));
            let q = x.checked_div(&y);
            if let Some(q) = &q {
                prop_assert_eq!(&all_forms!(x, /, /=, y), q);
            }
            q
        },
        gmp: {
            let d = gmp(&b);
            (d != 0).then(|| gmp(&a) / d)
        },
    }

    fn neg_abs_recip(a in rational()) {
        ours: {
            let x = ours(&a);
            prop_assert_eq!(-&x, -x.clone());
            (-&x, (x.abs(), x.checked_recip()))
        },
        gmp: {
            let g = gmp(&a);
            let recip = (g != 0).then(|| g.clone().recip());
            (-g.clone(), (g.abs(), recip))
        },
    }

    fn cmp((a, b) in rational_pair()) {
        ours: {
            let (x, y) = (ours(&a), ours(&b));
            prop_assert_eq!(x == y, x.cmp(&y).is_eq());
            prop_assert_eq!(y.cmp(&x), x.cmp(&y).reverse());
            x.cmp(&y)
        },
        gmp: gmp(&a).cmp(&gmp(&b)),
    }

    fn rounding(a in rational()) {
        ours: {
            let x = ours(&a);
            ((x.floor(), x.ceil()), (x.trunc(), (x.round(), x.round_ties_even())))
        },
        gmp: {
            let g = gmp(&a);
            let int = |q: Gmpq| q.numer().clone();
            (
                (int(g.clone().floor()), int(g.clone().ceil())),
                (int(g.clone().trunc()), (int(g.clone().round()), gmp_round_ties_even(&g))),
            )
        },
    }

    fn pow(a in rational(), e in -12i32..=12) {
        ours: {
            let x = ours(&a);
            if x.is_zero() && e < 0 { x } else { x.pow(i64::from(e)) }
        },
        gmp: {
            let g = gmp(&a);
            if g == 0 && e < 0 { g } else { rug::ops::Pow::pow(g, e) }
        },
    }

    fn to_string_radix(a in rational(), r in radix()) {
        ours: ours(&a).to_string_radix(r),
        gmp: gmp(&a).to_string_radix(i32::try_from(r).unwrap()),
    }

    /// GMP's output, unreduced by a factor `k`, in either case and with an
    /// optional `+`, parses back to the value.
    fn from_str_radix(a in rational(), r in radix(), k in 1i64..=12, upper in any::<bool>(), plus in any::<bool>()) {
        ours: {
            let g = gmp(&a);
            let r32 = i32::try_from(r).unwrap();
            let num = rug::Integer::from(g.numer() * k).to_string_radix(r32);
            let den = rug::Integer::from(g.denom() * k).to_string_radix(r32);
            let mut s = format!("{}{num}/{den}", if plus && g >= 0 { "+" } else { "" });
            if upper {
                s = s.to_uppercase();
            }
            Rational::from_str_radix(&s, r).unwrap()
        },
        gmp: gmp(&a),
    }

    fn display(a in rational()) {
        ours: {
            let x = ours(&a);
            prop_assert_eq!(format!("{x:?}"), x.to_string());
            prop_assert!(x.to_string().parse::<Rational>().unwrap() == x);
            x.to_string()
        },
        gmp: gmp(&a).to_string(),
    }
}

proptest! {
    /// Equal values (built from different unreduced fractions) hash equally.
    #[test]
    fn eq_and_hash((n, d) in rational(), k in strategy::integer()) {
        prop_assume!(!k.mag().is_zero());
        let x = ours(&(n.clone(), d.clone()));
        let k = k.to_integer();
        let y = Rational::new(n.to_integer() * &k, d.to_integer() * &k);
        prop_assert!(x == y);
        prop_assert_eq!(hash_of(&x), hash_of(&y));
    }

    /// Width, fill, alignment and sign flags as `rug` formats.
    #[test]
    fn format_flags(a in rational()) {
        let (x, g) = (ours(&a), gmp(&a));
        macro_rules! same {
            ($($fmt:literal),*) => {$(
                prop_assert_eq!(format!($fmt, x), format!($fmt, g), "format {}", $fmt);
            )*};
        }
        same!("{}", "{:+}", "{:>60}", "{:<60}", "{:^60}", "{:*^60}", "{:060}", "{:+060}");
    }

    #[test]
    fn integer_round_trip(n in strategy::integer()) {
        let n = n.to_integer();
        let x = Rational::from(n.clone());
        prop_assert!(x.is_integer());
        prop_assert!(Integer::try_from(&x).unwrap() == n);
    }
}

#[test]
fn edge_cases() {
    let q = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
    assert_eq!(q(0, -7), Rational::ZERO);
    assert_eq!(q(6, -4), q(-3, 2));
    assert_eq!(
        q(i64::MIN, -1).numer().to_rug(),
        rug::Integer::from(1u64 << 63)
    );
    assert_eq!(q(i64::MIN, i64::MIN), Rational::ONE);
    assert!(Rational::checked_new(Integer::ONE, Integer::ZERO).is_none());
    assert!(Rational::ONE.checked_div(&Rational::ZERO).is_none());
    assert!(Rational::ZERO.checked_recip().is_none());
    assert!(Integer::try_from(q(1, 2)).is_err());
    assert_eq!(q(-7, 2).round(), Integer::from(-4));
    assert_eq!(q(-5, 2).round_ties_even(), Integer::from(-2));
    assert_eq!(q(-3, 2).to_string(), "-3/2");
    assert_eq!(q(4, 2).to_string(), "2");
    assert_eq!(Rational::default(), Rational::ZERO);
    assert_eq!(q(-2, 3).pow(-3), q(-27, 8));
    assert_eq!(Rational::ZERO.pow(0), Rational::ONE);

    let parse = |s: &str| s.parse::<Rational>().map_err(|e| e.kind());
    assert_eq!(parse(""), Err(ParseErrorKind::Empty));
    assert_eq!(parse("/2"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1/"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1/-2"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1/+2"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1/2/3"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("1 /2"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("-/2"), Err(ParseErrorKind::InvalidDigit));
    assert_eq!(parse("3/0"), Err(ParseErrorKind::ZeroDenominator));
    assert_eq!(parse("-0/5"), Ok(Rational::ZERO));
    assert_eq!(parse("-6/4"), Ok(q(-3, 2)));
    assert_eq!(parse("+12"), Ok(q(12, 1)));
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn division_by_zero_panics() {
    let _ = Rational::ONE / Rational::ZERO;
}

#[test]
#[should_panic(expected = "attempt to divide by zero")]
fn zero_to_negative_power_panics() {
    let _ = Rational::ZERO.pow(-1);
}

#[test]
fn simplest_between() {
    let q = |n: i64, d: i64| Rational::new(n.into(), d.into());
    for (lo, hi, simplest) in [
        (q(3, 10), q(4, 10), q(1, 3)),
        (q(1, 3), q(1, 3), q(1, 3)),
        (q(-1, 2), q(1, 3), q(0, 1)),
        (q(5, 2), q(7, 2), q(3, 1)),
        (q(3, 1), q(7, 2), q(3, 1)),
        (q(-7, 2), q(-5, 2), q(-3, 1)),
        // Numbers that round to 0.333333333333 at 12 digits.
        (
            q(3_333_333_333_325, 10_000_000_000_000),
            q(3_333_333_333_335, 10_000_000_000_000),
            q(1, 3),
        ),
        // ... and to 0.7 at 12 digits: no fraction simpler than itself.
        (
            q(6_999_999_999_995, 10_000_000_000_000),
            q(7_000_000_000_005, 10_000_000_000_000),
            q(7, 10),
        ),
    ] {
        assert_eq!(
            Rational::simplest_between(&lo, &hi),
            Some(simplest),
            "[{lo}, {hi}]"
        );
    }
    assert_eq!(Rational::simplest_between(&q(1, 2), &q(1, 3)), None);
}

proptest! {
    /// Against a search of the denominators in order: the interval holds
    /// the result, no smaller denominator has a fraction in it, and none
    /// of the same denominator is smaller in magnitude.
    #[test]
    fn simplest_between_has_the_least_denominator(
        (a, b, c, d) in (-60i64..60, 1i64..40, -60i64..60, 1i64..40)
    ) {
        let (x, y) = (Rational::new(a.into(), b.into()), Rational::new(c.into(), d.into()));
        let (lo, hi) = if x <= y { (x, y) } else { (y, x) };
        let s = Rational::simplest_between(&lo, &hi).expect("not empty");
        prop_assert!(lo <= s && s <= hi);
        let found = (1i64..)
            .find_map(|den| {
                let den_q = Rational::from(Integer::from(den));
                let first = (&lo * &den_q).ceil();
                let last = (&hi * &den_q).floor();
                (first <= last).then(|| {
                    // The least magnitude: 0 if it is there, else the end nearest 0.
                    let zero = Integer::from(0);
                    let num = if first <= zero && zero <= last {
                        zero
                    } else if first.is_positive() {
                        first
                    } else {
                        last
                    };
                    Rational::new(num, den.into())
                })
            })
            .expect("a denominator");
        prop_assert_eq!(s, found);
    }
}
