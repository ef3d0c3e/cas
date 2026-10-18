//! `Complex` against GMP rationals part by part, the field identities, the
//! Gaussian division and GCD, and the printed form. The ring laws of
//! `GaussianRationals` and `GaussianIntegers` are in `algebra::tests`.
#![allow(
    clippy::many_single_char_names,
    reason = "the identities use the textbook names (a, b, c, z, w, q, r, g)"
)]

use proptest::prelude::*;

use super::{Complex, GaussianIntegers};
use crate::algebra::EuclideanDomain;
use crate::testkit::{ToRug, check_same, strategy};
use crate::{Integer, ParseErrorKind, Rational};

/// A rational with parts of at most `limbs` limbs.
fn rational_up_to(limbs: usize) -> BoxedStrategy<Rational> {
    (strategy::signed_up_to(limbs), strategy::signed_up_to(limbs))
        .prop_filter_map("zero denominator", |(n, d)| {
            Rational::checked_new(n.to_integer(), d.to_integer())
        })
        .boxed()
}

/// Complex rationals: general, real, purely imaginary, and units.
pub(crate) fn complex() -> BoxedStrategy<Complex> {
    let r = || rational_up_to(3);
    prop_oneof![
        6 => (r(), r()).prop_map(|(a, b)| Complex::new(a, b)),
        1 => r().prop_map(Complex::from),
        1 => r().prop_map(|b| Complex::new(Rational::ZERO, b)),
        1 => (0u8..4).prop_map(super::unit),
    ]
    .boxed()
}

/// Gaussian integers, including zero, units and pure parts.
pub(crate) fn gaussian() -> BoxedStrategy<Complex> {
    let i = || strategy::signed_up_to(3).prop_map(|x| Rational::from(x.to_integer()));
    prop_oneof![
        6 => (i(), i()).prop_map(|(a, b)| Complex::new(a, b)),
        1 => i().prop_map(Complex::from),
        1 => (0u8..4).prop_map(super::unit),
    ]
    .boxed()
}

fn rug_parts(z: &Complex) -> (rug::Rational, rug::Rational) {
    (z.re().to_rug(), z.im().to_rug())
}

proptest! {
    #[test]
    fn arithmetic_matches_gmp(z in complex(), w in complex()) {
        let ((a, b), (c, d)) = (rug_parts(&z), rug_parts(&w));
        let mul = |x: &rug::Rational, y: &rug::Rational| rug::Rational::from(x * y);
        check_same(&z + &w, (rug::Rational::from(&a + &c), rug::Rational::from(&b + &d)))?;
        check_same(&z - &w, (rug::Rational::from(&a - &c), rug::Rational::from(&b - &d)))?;
        check_same(&z * &w, (mul(&a, &c) - mul(&b, &d), mul(&a, &d) + mul(&b, &c)))?;
        check_same(-&z, (-a.clone(), -b.clone()))?;
        check_same(z.conj(), (a.clone(), -b.clone()))?;
        check_same(z.norm(), mul(&a, &a) + mul(&b, &b))?;
        let n = mul(&c, &c) + mul(&d, &d);
        let expected = (n != 0).then(|| {
            ((mul(&a, &c) + mul(&b, &d)) / &n, (mul(&b, &c) - mul(&a, &d)) / &n)
        });
        check_same(z.checked_div(&w), expected)?;
    }

    /// Every owned/borrowed form agrees with `&z op &w`.
    #[test]
    fn operator_forms(z in complex(), w in complex()) {
        let sum = crate::testkit::all_forms!(z.clone(), +, +=, w.clone());
        prop_assert_eq!(sum, &z + &w);
        let prod = crate::testkit::all_forms!(z.clone(), *, *=, w.clone());
        prop_assert_eq!(prod, &z * &w);
        if !w.is_zero() {
            let q = crate::testkit::all_forms!(z.clone(), /, /=, w.clone());
            prop_assert_eq!(q, &z / &w);
        }
    }

    #[test]
    fn field_identities(z in complex(), w in complex(), v in complex()) {
        prop_assert_eq!(&z * &(&w + &v), &(&z * &w) + &(&z * &v), "distributive");
        prop_assert_eq!((&z * &w).conj(), &z.conj() * &w.conj(), "conj multiplicative");
        prop_assert_eq!((&z * &w).norm(), &z.norm() * &w.norm(), "norm multiplicative");
        prop_assert_eq!(&z * &z.conj(), Complex::from(z.norm()), "z conj(z) = norm");
        match z.checked_recip() {
            None => prop_assert!(z.is_zero()),
            Some(r) => {
                prop_assert_eq!(&z * &r, Complex::ONE, "inverse");
                prop_assert_eq!(&(&w / &z) * &z, w.clone(), "division");
                prop_assert_eq!(z.pow(-3), r.pow(3));
            }
        }
        prop_assert_eq!(z.pow(5), &(&z * &z) * &(&(&z * &z) * &z));
        prop_assert_eq!(z.pow(0), Complex::ONE);
        prop_assert!((&z * &w).is_canonical() && (&z + &w).is_canonical());
    }

    /// `a = q b + r` with `norm(r) <= norm(b) / 2`, and `q, r` in ℤ[i].
    #[test]
    fn gaussian_division(a in gaussian(), b in gaussian()) {
        prop_assume!(!b.is_zero());
        let (q, r) = GaussianIntegers.div_rem(&a, &b);
        prop_assert!(q.is_gaussian_integer() && r.is_gaussian_integer());
        prop_assert_eq!(&(&q * &b) + &r, a.clone());
        prop_assert!(r.norm() * Rational::from(Integer::from(2)) <= b.norm());
    }

    /// The GCD of `g x` and `g y` is a multiple of `g`, divides both, and its
    /// norm divides both norms.
    #[test]
    fn gaussian_gcd(g in gaussian(), x in gaussian(), y in gaussian()) {
        let (a, b) = (&g * &x, &g * &y);
        let d = GaussianIntegers.gcd(&a, &b);
        let (n, _) = GaussianIntegers.normal_form(&d);
        prop_assert_eq!(&n, &d, "normal form");
        if d.is_zero() {
            prop_assert!(a.is_zero() && b.is_zero());
        } else {
            prop_assert!((&a / &d).is_gaussian_integer() && (&b / &d).is_gaussian_integer());
            prop_assert!(g.is_zero() || (&d / &g).is_gaussian_integer(), "g divides the gcd");
            let (na, nb, nd) = (a.norm(), b.norm(), d.norm());
            prop_assert!((&na / &nd).is_integer() && (&nb / &nd).is_integer());
            prop_assert!(n.re().is_positive() && !n.im().is_negative());
        }
    }

    #[test]
    fn print_parse_round_trip(z in complex()) {
        let s = z.to_string();
        prop_assert_eq!(s.parse::<Complex>(), Ok(z.clone()), "{}", s);
        prop_assert_eq!(format!("{z:?}"), s);
    }
}

#[test]
fn printed_forms() {
    let q = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
    let cases = [
        (Complex::ZERO, "0"),
        (Complex::ONE, "1"),
        (Complex::I, "I"),
        (-Complex::I, "-I"),
        (Complex::new(q(0, 1), q(-1, 2)), "-1/2*I"),
        (Complex::new(q(3, 1), q(1, 1)), "3 + I"),
        (Complex::new(q(1, 2), q(-3, 4)), "1/2 - 3/4*I"),
        (Complex::new(q(-2, 1), q(5, 1)), "-2 + 5*I"),
    ];
    for (z, s) in cases {
        assert_eq!(z.to_string(), s);
        assert_eq!(s.parse::<Complex>(), Ok(z));
    }
    // Other spellings of the same values parse and reduce.
    assert_eq!("1 + 0*I".parse::<Complex>(), Ok(Complex::ONE));
    assert_eq!(
        "2/4*I".parse::<Complex>(),
        Ok(Complex::new(q(0, 1), q(1, 2)))
    );
}

#[test]
fn parse_errors() {
    let kind = |s: &str| s.parse::<Complex>().unwrap_err().kind();
    assert_eq!(kind(""), ParseErrorKind::Empty);
    for s in [
        "1 + -2*I",
        "1+2*I",
        "*I",
        "1 + *I",
        "2I",
        "1 + 2*I + 3*I",
        "i",
        "1 +  I",
        "--I",
        "1 - +I",
    ] {
        assert_eq!(kind(s), ParseErrorKind::InvalidDigit, "{s:?}");
    }
    assert_eq!(kind("1/0 + I"), ParseErrorKind::ZeroDenominator);
}

#[test]
fn division_by_zero() {
    assert_eq!(Complex::ONE.checked_div(&Complex::ZERO), None);
    assert_eq!(Complex::ZERO.checked_recip(), None);
    assert!(std::panic::catch_unwind(|| Complex::ONE / Complex::ZERO).is_err());
    assert!(std::panic::catch_unwind(|| Complex::ZERO.pow(-1)).is_err());
    assert_eq!(Complex::ZERO.pow(0), Complex::ONE);
}

/// The normal forms of the associates of `2 + i` are all `2 + i`.
#[test]
fn associates_share_a_normal_form() {
    let z = Complex::new(Rational::from(Integer::from(2)), Rational::ONE);
    for k in 0..4 {
        let a = &z * &super::unit(k);
        let (n, u) = GaussianIntegers.normal_form(&a);
        assert_eq!(n, z);
        assert_eq!(&a * &u, z);
    }
    // 5 = (2 + i)(2 - i) and 3 + i = (1 + i)(2 - i), so the gcd is 2 - i up to
    // a unit; its first-quadrant associate is (2 - i) i = 1 + 2i.
    let five = Complex::from(5);
    let w = Complex::new(Rational::from(Integer::from(3)), Rational::ONE);
    let g = GaussianIntegers.gcd(&five, &w);
    assert_eq!(
        g,
        Complex::new(Rational::ONE, Rational::from(Integer::from(2)))
    );
}
