//! The Jacobi symbol checked against GMP.

use proptest::prelude::*;

use crate::Integer;
use crate::testkit::{Mag, ToRug, differential, strategy};

/// An odd positive modulus, single-limb or multi-limb.
fn odd() -> BoxedStrategy<Mag> {
    strategy::natural_up_to(6)
        .prop_map(|m| Mag::from_rug(&(m.to_rug() | 1)))
        .boxed()
}

differential! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    fn jacobi(a in strategy::integer(), n in odd()) {
        ours: a.to_integer().jacobi(&Integer::from(n.to_natural())),
        gmp: a.to_rug().jacobi(&n.to_rug()),
    }

    /// Operands sharing a factor give 0.
    fn jacobi_common_factor(a in strategy::signed_up_to(3), n in odd(), k in odd()) {
        ours: (a.to_integer() * Integer::from(k.to_natural())).jacobi(&(Integer::from(n.to_natural()) * Integer::from(k.to_natural()))),
        gmp: (a.to_rug() * k.to_rug()).jacobi(&(n.to_rug() * k.to_rug())),
    }
}

#[test]
fn small_exhaustive() {
    for n in (1..200i64).step_by(2) {
        for a in -200..200i64 {
            let ours = Integer::from(a).jacobi(&Integer::from(n));
            let gmp = rug::Integer::from(a).jacobi(&rug::Integer::from(n));
            assert_eq!(ours, gmp, "({a} / {n})");
        }
    }
}

#[test]
#[should_panic(expected = "odd positive modulus")]
fn even_modulus_panics() {
    let _ = Integer::ONE.jacobi(&Integer::from(4));
}

#[test]
#[should_panic(expected = "odd positive modulus")]
fn negative_modulus_panics() {
    let _ = Integer::ONE.jacobi(&Integer::from(-3));
}
