#![allow(
    clippy::many_single_char_names,
    reason = "a, b, g, s, t mirror the notation"
)]
//! `Natural` GCD and LCM checked against GMP, and the Lehmer machinery checked
//! on its own.

use proptest::prelude::*;
use rug::Integer;

use super::{Step, reduce};
use crate::Natural;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

fn pair() -> BoxedStrategy<(Mag, Mag)> {
    prop_oneof![
        4 => strategy::gcd_pair(DEFAULT_MAX_LIMBS),
        1 => strategy::unbalanced_pair(DEFAULT_MAX_LIMBS),
    ]
    .boxed()
}

/// Matrix and division steps `reduce` takes on `(a, b)`.
fn steps(a: &Natural, b: &Natural) -> (usize, usize) {
    let (mut x, mut y) = super::ordered(a, b);
    let (mut matrices, mut divisions) = (0, 0);
    reduce(&mut x, &mut y, |step| match step {
        Step::Matrix(_) => matrices += 1,
        Step::Quotient(_) => divisions += 1,
    });
    (matrices, divisions)
}

/// Uniformly random limbs with a top limb of at least `2^32`.
fn dense(len: usize) -> impl Strategy<Value = Natural> {
    proptest::collection::vec(any::<u64>(), len).prop_map(move |mut limbs| {
        limbs[len - 1] |= 1 << 32;
        Natural::from_limbs(&limbs)
    })
}

fn fibonacci(n: u32) -> Natural {
    Mag::from_rug(&Integer::from(Integer::fibonacci(n))).to_natural()
}

differential! {
    fn gcd((a, b) in pair()) {
        ours: {
            let (x, y) = (a.to_natural(), b.to_natural());
            let g = x.gcd(&y);
            prop_assert_eq!(&g, &y.gcd(&x));
            g
        },
        gmp: a.to_rug().gcd(&b.to_rug()),
    }

    fn lcm((a, b) in pair()) {
        ours: a.to_natural().lcm(&b.to_natural()),
        gmp: a.to_rug().lcm(&b.to_rug()),
    }

    /// `g = s a + t b` for an integer `t`, i.e. `b` divides `s a - g`.
    fn gcd_cofactor((a, b) in pair()) {
        ours: {
            let (g, s) = a.to_natural().gcd_cofactor(&b.to_natural());
            let rest = s.to_rug() * a.to_rug() - g.to_rug();
            let b = b.to_rug();
            prop_assert!(
                if b == 0 { rest == 0 } else { rest.is_divisible(&b) },
                "s = {:?} is not a Bézout coefficient", s
            );
            g
        },
        gmp: a.to_rug().gcd(&b.to_rug()),
    }
}

proptest! {
    /// On dense operands nearly every step is a Lehmer matrix; the full
    /// division is the exception.
    #[test]
    fn lehmer_does_the_work(a in dense(12), b in dense(12)) {
        let (matrices, divisions) = steps(&a, &b);
        prop_assert!(matrices >= 8 * divisions, "{} matrices, {} divisions", matrices, divisions);
        prop_assert!(matrices >= 12, "{} matrices", matrices);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Operands of hundreds of limbs, where division steps take the
    /// divide-and-conquer path.
    #[test]
    fn large_operands((a, b, c) in (dense(200), dense(150), dense(120))) {
        let (x, y) = (&a * &c, &b * &c);
        let g = x.gcd(&y);
        prop_assert!(g.to_rug() == x.to_rug().gcd(&y.to_rug()));
        let (g2, s) = x.gcd_cofactor(&y);
        prop_assert!(g2 == g);
        prop_assert!((s.to_rug() * x.to_rug() - g.to_rug()).is_divisible(&y.to_rug()));
    }
}

#[test]
fn fibonacci_worst_case_uses_matrices() {
    let (a, b) = (fibonacci(5001), fibonacci(5000));
    assert!(a.gcd(&b).is_one());
    let (matrices, divisions) = steps(&a, &b);
    assert!(
        matrices >= 100 && divisions <= matrices / 8,
        "{matrices} matrices, {divisions} divisions"
    );
}

#[test]
fn small_cases() {
    let n = |x: u64| Natural::from(x);
    assert_eq!(n(0).gcd(&n(0)), n(0));
    assert_eq!(n(0).gcd(&n(12)), n(12));
    assert_eq!(n(12).gcd(&n(0)), n(12));
    assert_eq!(n(12).gcd(&n(18)), n(6));
    assert_eq!(n(0).lcm(&n(5)), n(0));
    assert_eq!(n(4).lcm(&n(6)), n(12));
    let big = Natural::power_of_two(200);
    assert_eq!(
        big.gcd(&Natural::power_of_two(130)),
        Natural::power_of_two(130)
    );
    assert_eq!(big.gcd(&n(0)), big);
    assert_eq!(n(0).gcd(&big), big);
}
