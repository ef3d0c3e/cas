//! `RealBalls`: real balls at a fixed precision as a ring object (§5.1).

use super::Ball;
use crate::algebra::Ring;
use crate::{Integer, Natural};

/// The real numbers as balls at a working precision of `prec` bits.
///
/// Implements [`Ring`] only. It is deliberately not a `EuclideanDomain` or a
/// `Field`: algorithms that branch on exact zero tests (GCD, elimination,
/// degree normalisation) would go silently wrong on a ball that merely
/// contains 0. `is_zero` means the exact zero ball.
///
/// ```compile_fail
/// use cas_num::algebra::Field;
/// use cas_num::ball::RealBalls;
///
/// fn needs_a_field<F: Field>(_: &F) {}
/// needs_a_field(&RealBalls::new(64));
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct RealBalls {
    prec: u64,
}

impl RealBalls {
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn new(prec: u64) -> RealBalls {
        assert!(prec >= 1, "RealBalls: precision 0");
        RealBalls { prec }
    }

    /// The working precision in bits.
    #[must_use]
    pub fn prec(&self) -> u64 {
        self.prec
    }

    /// `a / b`; indeterminate if `b` contains 0.
    #[must_use]
    pub fn div(&self, a: &Ball, b: &Ball) -> Ball {
        a.div(b, self.prec)
    }

    /// `1 / a`; indeterminate if `a` contains 0.
    #[must_use]
    pub fn inv(&self, a: &Ball) -> Ball {
        a.inv(self.prec)
    }

    /// `sqrt(a)`; indeterminate if `a` reaches below 0.
    #[must_use]
    pub fn sqrt(&self, a: &Ball) -> Ball {
        a.sqrt(self.prec)
    }
}

impl Ring for RealBalls {
    type Elem = Ball;

    fn zero(&self) -> Ball {
        Ball::ZERO
    }

    fn one(&self) -> Ball {
        Ball::from_integer(&Integer::ONE)
    }

    fn add(&self, a: &Ball, b: &Ball) -> Ball {
        a.add(b, self.prec)
    }

    fn sub(&self, a: &Ball, b: &Ball) -> Ball {
        a.sub(b, self.prec)
    }

    fn neg(&self, a: &Ball) -> Ball {
        a.neg()
    }

    fn mul(&self, a: &Ball, b: &Ball) -> Ball {
        a.mul(b, self.prec)
    }

    /// The exact ball of `n`.
    fn elem_integer(&self, n: &Integer) -> Ball {
        Ball::from_integer(n)
    }

    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }

    fn pow(&self, a: &Ball, e: u64) -> Ball {
        if let Ok(e) = i64::try_from(e) {
            return a.pow(e, self.prec);
        }
        // Beyond i64: the generic square-and-multiply.
        let mut r = self.one();
        for i in (0..u64::BITS - e.leading_zeros()).rev() {
            r = self.square(&r);
            if (e >> i) & 1 == 1 {
                r = self.mul(&r, a);
            }
        }
        r
    }
}
