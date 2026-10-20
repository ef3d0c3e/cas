//! Exact numbers and infinities, as the arithmetic natives see them.

use cas_num::{Complex, Integer, Rational};

use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The value of a number atom: real values as rationals, so real arithmetic
/// never pays for an imaginary part.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Num {
    Real(Rational),
    Complex(Complex),
}

impl Num {
    /// The value of `e`, if it is a number atom.
    pub(super) fn of(e: &Expr) -> Option<Num> {
        match e.view() {
            View::Complex(z) => Some(Num::Complex(z.clone())),
            _ => e.to_rational().map(Num::Real),
        }
    }

    pub(super) fn zero() -> Num {
        Num::Real(Rational::from(0))
    }

    pub(super) fn one() -> Num {
        Num::Real(Rational::from(1))
    }

    pub(super) fn is_zero(&self) -> bool {
        match self {
            Num::Real(q) => q.is_zero(),
            Num::Complex(z) => z.is_zero(),
        }
    }

    pub(super) fn is_one(&self) -> bool {
        match self {
            Num::Real(q) => q.is_one(),
            Num::Complex(z) => z.is_one(),
        }
    }

    pub(super) fn real(&self) -> Option<&Rational> {
        match self {
            Num::Real(q) => Some(q),
            Num::Complex(_) => None,
        }
    }

    fn complex(&self) -> Complex {
        match self {
            Num::Real(q) => Complex::from(q.clone()),
            Num::Complex(z) => z.clone(),
        }
    }

    /// The canonical form of a complex value: a real one is real.
    fn from_complex(z: Complex) -> Num {
        if z.is_real() {
            Num::Real(z.re().clone())
        } else {
            Num::Complex(z)
        }
    }

    pub(super) fn add(&self, other: &Num) -> Num {
        match (self, other) {
            (Num::Real(a), Num::Real(b)) => Num::Real(a + b),
            _ => Num::from_complex(&self.complex() + &other.complex()),
        }
    }

    pub(super) fn mul(&self, other: &Num) -> Num {
        match (self, other) {
            (Num::Real(a), Num::Real(b)) => Num::Real(a * b),
            _ => Num::from_complex(&self.complex() * &other.complex()),
        }
    }

    /// `self^n`; `None` for zero to a negative power.
    pub(super) fn pow(&self, n: i64) -> Option<Num> {
        if self.is_zero() && n < 0 {
            return None;
        }
        Some(match self {
            Num::Real(q) => Num::Real(q.pow(n)),
            Num::Complex(z) => Num::from_complex(z.pow(n)),
        })
    }

    /// An estimate of the bits of `self^n`, to refuse powers too large to
    /// hold.
    pub(super) fn power_bits(&self, n: i64) -> u64 {
        let bits = |q: &Rational| q.numer().bit_length().max(q.denom().bit_length());
        let base = match self {
            Num::Real(q) => bits(q),
            Num::Complex(z) => bits(z.re()).max(bits(z.im())) + 1,
        };
        base.saturating_mul(n.unsigned_abs())
    }

    /// `self / |self|` when it is exact: the sign of a real number, or a
    /// complex number whose modulus is rational. `None` for zero and for a
    /// complex number with an irrational modulus.
    pub(super) fn unit(&self) -> Option<Num> {
        match self {
            Num::Real(q) if q.is_zero() => None,
            Num::Real(q) => Some(Num::Real(Rational::from(if q.is_negative() {
                -1
            } else {
                1
            }))),
            Num::Complex(z) => {
                let modulus = rational_sqrt(&z.norm())?;
                Some(Num::from_complex(Complex::new(
                    z.re() / &modulus,
                    z.im() / &modulus,
                )))
            }
        }
    }

    pub(super) fn to_expr(&self) -> Expr {
        match self {
            Num::Real(q) => Expr::from(q),
            Num::Complex(z) => Expr::from(z.clone()),
        }
    }
}

/// The square root of a non-negative rational, if it is rational.
fn rational_sqrt(q: &Rational) -> Option<Rational> {
    let root = |n: &Integer| {
        let (r, rem) = n.sqrt_rem();
        rem.is_zero().then_some(r)
    };
    Some(Rational::new(root(q.numer())?, root(q.denom())?))
}

/// An infinite or indeterminate quantity (CAS-DOC-4 R14).
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Infinity {
    /// `DirectedInfinity[d]`: the direction `d` has modulus 1. `Infinity` is
    /// the direction 1.
    Directed(Num),
    /// `ComplexInfinity`, without a direction.
    Complex,
    Indeterminate,
}

impl Infinity {
    /// What `e` is, if it is an infinity or `Indeterminate`: `Infinity`,
    /// `ComplexInfinity`, `DirectedInfinity[d]` with a number `d`, `DirectedInfinity[]`
    /// or `Indeterminate`.
    pub(super) fn of(e: &Expr) -> Option<Infinity> {
        match e.view() {
            View::Symbol(Symbol::INFINITY) => Some(Infinity::Directed(Num::one())),
            View::Symbol(Symbol::COMPLEX_INFINITY) => Some(Infinity::Complex),
            View::Symbol(Symbol::INDETERMINATE) => Some(Infinity::Indeterminate),
            View::Apply { head, args } if head.as_symbol() == Some(Symbol::DIRECTED_INFINITY) => {
                match args {
                    [] => Some(Infinity::Complex),
                    [d] => Num::of(d).map(Infinity::Directed),
                    _ => None,
                }
            }
            _ => None,
        }
    }

    /// The sum of two infinite terms.
    pub(super) fn add(self, other: Infinity) -> Infinity {
        match (self, other) {
            (Infinity::Directed(a), Infinity::Directed(b)) if a == b => Infinity::Directed(a),
            _ => Infinity::Indeterminate,
        }
    }

    /// The product of two infinite factors.
    pub(super) fn mul(self, other: Infinity) -> Infinity {
        match (self, other) {
            (Infinity::Indeterminate, _) | (_, Infinity::Indeterminate) => Infinity::Indeterminate,
            (Infinity::Directed(a), Infinity::Directed(b)) => Infinity::Directed(a.mul(&b)),
            _ => Infinity::Complex,
        }
    }

    /// The canonical expression: `Infinity` for the direction 1,
    /// `DirectedInfinity[d]` for another, `ComplexInfinity`,
    /// `Indeterminate`.
    pub(super) fn to_expr(&self) -> Expr {
        match self {
            Infinity::Directed(d) if d.is_one() => Symbol::INFINITY.into(),
            Infinity::Directed(d) => Expr::apply(&Symbol::DIRECTED_INFINITY.into(), &[d.to_expr()]),
            Infinity::Complex => Symbol::COMPLEX_INFINITY.into(),
            Infinity::Indeterminate => Symbol::INDETERMINATE.into(),
        }
    }
}
