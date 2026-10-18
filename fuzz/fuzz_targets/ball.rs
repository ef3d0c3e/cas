//! Ball arithmetic against exact `Rational` evaluation on dyadic inputs
//! (CAS-DOC-5 §11): every result ball must contain the exact result at the
//! input midpoints (which lie in the input balls), and division by a ball
//! containing 0 must give the indeterminate ball.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Rational;
use cas_num::ball::{Ball, Float, Mag, Round};
use cas_num_fuzz::{check, dyadic, gmp_rat};
use libfuzzer_sys::fuzz_target;

/// The exact float of a dyadic rational (at most 16 limbs, the square of an
/// 8-limb input, so 2048 bits hold it).
fn exact(q: &Rational) -> Float {
    let (f, o) = Float::from_rational(q, 2048, Round::Nearest);
    assert_eq!(o, std::cmp::Ordering::Equal, "dyadic inputs are exact");
    f
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok(prec), Ok((a, ag)), Ok((b, bg)), Ok((r, _)), Ok(radii)) = (
        u.int_in_range(0u8..=4),
        u.int_in_range(1u64..=512),
        dyadic(&mut u, 8),
        dyadic(&mut u, 8),
        dyadic(&mut u, 2),
        u.int_in_range(0u8..=3),
    ) else {
        return;
    };
    let rad = Mag::from_float(&exact(&r));
    let ball =
        |q: &Rational, with: bool| Ball::new(exact(q), if with { rad.clone() } else { Mag::ZERO });
    let (x, y) = (ball(&a, radii & 1 == 1), ball(&b, radii & 2 == 2));
    let (result, value) = match op {
        0 => (x.add(&y, prec), Some(&a + &b)),
        1 => (x.sub(&y, prec), Some(&a - &b)),
        2 => (x.mul(&y, prec), Some(&a * &b)),
        3 => (x.div(&y, prec), a.checked_div(&b)),
        _ => {
            // sqrt(a^2) = |a|, a case with an exact answer.
            let sq = &a * &a;
            (
                ball(&sq, radii & 1 == 1).sqrt(prec),
                Some(Rational::from(a.numer().abs()) / Rational::from(a.denom().clone())),
            )
        }
    };
    let gmp = match op {
        0 => Some(ag.clone() + &bg),
        1 => Some(ag.clone() - &bg),
        2 => Some(ag.clone() * &bg),
        3 => (bg != 0).then(|| ag.clone() / &bg),
        _ => Some(ag.clone().abs()),
    };
    check("exact reference", value.as_ref().map(gmp_rat), gmp);
    match value {
        Some(v) => assert!(
            result.contains_rational(&v),
            "op {op} at {prec}: {result:?} misses {v}"
        ),
        None => assert_eq!(
            result,
            Ball::INDETERMINATE,
            "division by a ball containing 0"
        ),
    }
    if op == 3 && y.contains(&Float::ZERO) {
        assert_eq!(result, Ball::INDETERMINATE);
    }
});
