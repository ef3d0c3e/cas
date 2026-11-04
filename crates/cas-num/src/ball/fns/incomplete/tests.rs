//! Tests of `Ball::gamma_p`, `Ball::gamma_q` and `Ball::beta_i` (CAS-468):
//! 50-digit references (Legendre's continued fraction and Numerical
//! Recipes' `betacf` in Python's `decimal` at 90 digits, with `log Γ` by
//! Stirling's series; giac's doubles and the closed forms agree), balls,
//! the domain, and `Γ(a, x) = Q(a, x) Γ(a)` against MPFR's `gamma_inc`.

use proptest::prelude::*;
use rug::Float as MpFloat;
use rug::float::Round as MpRound;

use crate::ball::{Ball, Float, Mag};
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses};
use crate::{Integer, Rational};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn exact(r: &Rational) -> Ball {
    Ball::from_rational(r, 256)
}

/// Whether `b` holds the decimal `reference` (50 digits, so within
/// `10⁻⁴⁹` of it relatively) and has `want` bits.
fn holds(b: &Ball, reference: &str, want: i64) -> Result<(), String> {
    let r = MpFloat::with_val(512, MpFloat::parse(reference).unwrap());
    let mid = MpFloat::with_val(4096, b.mid().to_rational().to_rug());
    let rad = MpFloat::with_val(4096, b.rad().to_float().unwrap().to_rational().to_rug());
    let slack = MpFloat::with_val(512, r.abs_ref()) * MpFloat::with_val(512, 1e-49);
    let gap = MpFloat::with_val(4096, &mid - &r).abs();
    if gap > rad + slack {
        return Err(format!("{b:?} misses {reference}"));
    }
    if b.rel_accuracy_bits() < want {
        return Err(format!("{b:?}: {} bits", b.rel_accuracy_bits()));
    }
    Ok(())
}

#[test]
fn gamma_values() {
    let cases: [(&str, Rational, Rational, &str); 10] = [
        (
            "P",
            q(5, 2),
            q(3, 2),
            "3.00014164121372490900198448375331778846520914019585e-1",
        ),
        (
            "Q",
            q(5, 2),
            q(3, 2),
            "6.99985835878627509099801551624668221153479085980415e-1",
        ),
        (
            "P",
            q(1, 3),
            q(10, 1),
            "9.99996559781516976120849071869077646872754361199009e-1",
        ),
        (
            "Q",
            q(1, 3),
            q(100, 1),
            "6.40321467929021216959579737374031238609132697157907e-46",
        ),
        (
            "Q",
            q(7, 1),
            q(200, 1),
            "1.26797548256697582011420469256507110278887271267867e-76",
        ),
        (
            "P",
            q(100, 1),
            q(90, 1),
            "1.58220989186430168104969699670910531699823345743347e-1",
        ),
        (
            "P",
            q(3, 10),
            q(1, 1000),
            "1.40242458924867359534928716286394738503584478730885e-1",
        ),
        (
            "Q",
            q(3, 10),
            q(1, 1000),
            "8.59757541075132640465071283713605261496415521269115e-1",
        ),
        (
            "Q",
            q(5, 2),
            q(30, 1),
            "1.21545697771830389483485299983775770353223862630395e-11",
        ),
        (
            "P",
            q(1, 1),
            q(2, 1),
            "8.64664716763387308106000505027515596592368454090424e-1",
        ),
    ];
    let mut wrong = Vec::new();
    for prec in [53, 128] {
        let want = i64::try_from(prec).unwrap() - 3;
        for (f, a, x, reference) in &cases {
            let (a, x) = (exact(a), exact(x));
            let b = if *f == "P" {
                Ball::gamma_p(&a, &x, prec)
            } else {
                Ball::gamma_q(&a, &x, prec)
            };
            if let Err(e) = holds(&b, reference, want) {
                wrong.push(format!("{f}({a:?}, {x:?}) at {prec}: {e}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn beta_values() {
    let mut wrong = Vec::new();
    let betas: [(Rational, Rational, Rational, &str); 5] = [
        (
            q(5, 2),
            q(3, 2),
            q(3, 10),
            "8.89437231706655993544962677753584934052880057156755e-2",
        ),
        (
            q(2, 1),
            q(3, 1),
            q(9, 10),
            "9.96300000000000000000000000000000000000000000000000e-1",
        ),
        (
            q(1, 2),
            q(1, 3),
            q(99, 100),
            "8.46158403063435878082629336362436842319912538748412e-1",
        ),
        (
            q(50, 1),
            q(60, 1),
            q(1, 2),
            "8.30907293901669414336975456958316912284840704428194e-1",
        ),
        (
            q(1, 2),
            q(1, 2),
            q(1, 4),
            "3.33333333333333333333333333333333333333333333333333e-1",
        ),
    ];
    for prec in [53, 128] {
        let want = i64::try_from(prec).unwrap() - 3;
        for (a, b, x, reference) in &betas {
            let v = Ball::beta_i(&exact(a), &exact(b), &exact(x), prec);
            if let Err(e) = holds(&v, reference, want) {
                wrong.push(format!("I_{x}({a}, {b}) at {prec}: {e}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Balls hold the values at their points, and `P + Q` holds 1.
#[test]
fn balls_and_domain() {
    let a = Ball::new(
        Float::from(Integer::from(5)).mul_2exp(&-Integer::ONE),
        Mag::power_of_two(&Integer::from(-20)),
    );
    let x = exact(&q(3, 2));
    let p = Ball::gamma_p(&a, &x, 64);
    holds(
        &p,
        "3.00014164121372490900198448375331778846520914019585e-1",
        0,
    )
    .unwrap();
    assert!(p.rel_accuracy_bits() > 15, "{p:?}");
    let sum = p.add(&Ball::gamma_q(&a, &x, 64), 64);
    assert!(sum.contains(&Float::from(Integer::ONE)), "{sum:?}");
    let i = Ball::beta_i(
        &exact(&q(2, 1)),
        &exact(&q(3, 1)),
        &Ball::new(
            Float::from(Integer::from(9)).mul_2exp(&-Integer::from(4)),
            Mag::power_of_two(&Integer::from(-30)),
        ),
        64,
    );
    // I_x(2, 3) at 9/16: 1 − (1−x)^3 (1 + 3x), by the binomial sum.
    let y = q(7, 16);
    let expected = Rational::from(Integer::ONE)
        - &(&(&y * &y) * &y) * &(Rational::from(Integer::ONE) + &(q(3, 1) * &q(9, 16)));
    assert!(i.contains_rational(&expected), "{i:?} misses {expected}");
    // The domain.
    let zero = Ball::ZERO;
    assert!(!Ball::gamma_p(&zero, &x, 53).is_finite());
    assert!(!Ball::gamma_q(&exact(&q(1, 2)), &exact(&q(-1, 2)), 53).is_finite());
    assert!(!Ball::beta_i(&exact(&q(1, 2)), &exact(&q(1, 2)), &exact(&q(3, 2)), 53).is_finite());
    // The ends.
    assert_eq!(Ball::gamma_p(&exact(&q(1, 2)), &zero, 53), Ball::ZERO);
    assert!(Ball::gamma_q(&exact(&q(1, 2)), &zero, 53).contains(&Float::from(Integer::ONE)));
    assert!(
        Ball::beta_i(&exact(&q(1, 2)), &exact(&q(3, 1)), &exact(&q(1, 1)), 53)
            .contains(&Float::from(Integer::ONE))
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// `Q(a, x) Γ(a)` holds MPFR's `Γ(a, x)`.
    #[test]
    fn upper_against_mpfr(an in 1i64..400, ad in 1i64..8, xn in 1i64..2000, xd in 1i64..8, prec in 20u64..200) {
        let (a, x) = (q(an, ad), q(xn, xd));
        let ours = Ball::gamma_q(&exact(&a), &exact(&x), prec).mul(&exact(&a).gamma(prec + 16), prec);
        prop_assume!(ours.is_finite());
        let (am, xm) = (MpFloat::with_val(4096, a.to_rug()), MpFloat::with_val(4096, x.to_rug()));
        let br = bracket(u32::try_from(prec).unwrap(), |p, r: MpRound| {
            MpFloat::with_val_round(p, am.gamma_inc_ref(&xm), r).0
        });
        let mid = ours.mid().to_rational().to_rug();
        let rad = ours.rad().to_float().unwrap().to_rational().to_rug();
        check_encloses(&mid, &rad, &br)?;
    }
}
