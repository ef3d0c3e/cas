//! Numbers and infinities, as the arithmetic natives see them: exact ones,
//! and balls (CAS-DOC-5 §9.1), which absorb exact ones.

use std::cmp::Ordering;

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Complex, Integer, Rational};

use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The value of a number atom: real values as rationals, so real arithmetic
/// never pays for an imaginary part, and balls, real or complex. The balls
/// are boxed, which keeps a `Num` as small as a `Complex`: exact arithmetic
/// moves no more than it did before them.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Num {
    Real(Rational),
    Complex(Complex),
    Ball(Box<Ball>),
    /// Its imaginary part is not exactly 0.
    ComplexBall(Box<ComplexBall>),
}

/// The guard bits of ball arithmetic: an operation runs at this many bits
/// more than its result keeps.
const GUARD: u64 = 16;

/// The least precision ball arithmetic keeps, a double's.
const MIN_BITS: u64 = 53;

/// The bits a ball result keeps (CAS-DOC-5 §9.1): the least accuracy of the
/// inexact balls among `operands`, since the least accurate bounds the
/// result, or the longest midpoint of exact ones; at least [`MIN_BITS`].
/// `None` without a ball.
fn ball_bits(operands: &[&Num]) -> Option<u64> {
    let mut inexact: Option<u64> = None;
    let mut exact: Option<u64> = None;
    let mut see = |b: &Ball| {
        if b.is_exact() {
            exact = Some(exact.unwrap_or(0).max(b.mid().bit_length()));
        } else {
            let accuracy = u64::try_from(b.rel_accuracy_bits()).unwrap_or(0);
            inexact = Some(inexact.map_or(accuracy, |p| p.min(accuracy)));
        }
    };
    for n in operands {
        match n {
            Num::Ball(b) => see(b),
            Num::ComplexBall(z) => {
                see(z.re());
                see(z.im());
            }
            Num::Real(_) | Num::Complex(_) => {}
        }
    }
    inexact.or(exact).map(|bits| bits.max(MIN_BITS))
}

/// A complex ball rounded to `bits`.
fn round_complex(z: &ComplexBall, bits: u64) -> ComplexBall {
    ComplexBall::new(
        z.re().round_to_accuracy(bits),
        z.im().round_to_accuracy(bits),
    )
}

impl Num {
    /// The value of `e`, if it is a number atom.
    pub(super) fn of(e: &Expr) -> Option<Num> {
        match e.view() {
            View::Complex(z) => Some(Num::Complex(z.clone())),
            View::Ball(b) => Some(Num::Ball(Box::new(b.clone()))),
            View::ComplexBall(z) => Some(Num::ComplexBall(Box::new(z.clone()))),
            _ => e.to_rational().map(Num::Real),
        }
    }

    pub(super) fn zero() -> Num {
        Num::Real(Rational::from(0))
    }

    pub(super) fn one() -> Num {
        Num::Real(Rational::from(1))
    }

    /// Exactly 0, as `0.0` is too.
    pub(super) fn is_zero(&self) -> bool {
        match self {
            Num::Real(q) => q.is_zero(),
            Num::Complex(z) => z.is_zero(),
            Num::Ball(b) => b.is_exact() && b.mid().is_zero(),
            Num::ComplexBall(_) => false,
        }
    }

    /// Exactly 1 and exact: a ball coefficient stays, so `1.0 x` keeps the
    /// mark of its inexactness, as in Mathematica.
    pub(super) fn is_one(&self) -> bool {
        match self {
            Num::Real(q) => q.is_one(),
            Num::Complex(z) => z.is_one(),
            Num::Ball(_) | Num::ComplexBall(_) => false,
        }
    }

    pub(super) fn real(&self) -> Option<&Rational> {
        match self {
            Num::Real(q) => Some(q),
            _ => None,
        }
    }

    /// Whether this is a ball, real or complex.
    pub(super) fn is_ball(&self) -> bool {
        matches!(self, Num::Ball(_) | Num::ComplexBall(_))
    }

    fn complex(&self) -> Complex {
        match self {
            Num::Real(q) => Complex::from(q.clone()),
            Num::Complex(z) => z.clone(),
            Num::Ball(_) | Num::ComplexBall(_) => unreachable!("an exact number"),
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

    /// The canonical form of a complex ball: one whose imaginary part is
    /// exactly 0 is a real ball.
    fn from_complex_ball(z: ComplexBall) -> Num {
        if z.is_real() {
            Num::Ball(Box::new(z.re().clone()))
        } else {
            Num::ComplexBall(Box::new(z))
        }
    }

    /// A real value as a ball at `prec` bits.
    fn to_ball(&self, prec: u64) -> Option<Ball> {
        match self {
            Num::Real(q) => Some(Ball::from_rational(q, prec)),
            Num::Ball(b) => Some(Ball::clone(b)),
            Num::Complex(_) | Num::ComplexBall(_) => None,
        }
    }

    /// The value as a complex ball at `prec` bits.
    fn to_complex_ball(&self, prec: u64) -> ComplexBall {
        match self {
            Num::Real(q) => ComplexBall::from_real(Ball::from_rational(q, prec)),
            Num::Complex(z) => ComplexBall::from_complex(z, prec),
            Num::Ball(b) => ComplexBall::from_real(Ball::clone(b)),
            Num::ComplexBall(z) => ComplexBall::clone(z),
        }
    }

    /// `real(a, b, bits)` on real balls or `complex(a, b, bits)` on complex
    /// ones, at the guarded precision of [`ball_bits`], rounded back.
    fn ball_op(
        &self,
        other: &Num,
        real: impl Fn(&Ball, &Ball, u64) -> Ball,
        complex: impl Fn(&ComplexBall, &ComplexBall, u64) -> ComplexBall,
    ) -> Num {
        let bits = ball_bits(&[self, other]).expect("a ball");
        let working = bits + GUARD;
        if let (Some(a), Some(b)) = (self.to_ball(working), other.to_ball(working)) {
            return Num::Ball(Box::new(real(&a, &b, working).round_to_accuracy(bits)));
        }
        let (a, b) = (
            self.to_complex_ball(working),
            other.to_complex_ball(working),
        );
        Num::from_complex_ball(round_complex(&complex(&a, &b, working), bits))
    }

    /// Exactly 0 or 1, and exact: identities of ball arithmetic too, which
    /// must leave a ball as it is. Ball arithmetic by an exact 1 would still
    /// grow the radius, so `Times[ball, x]` would not be a fixpoint.
    fn is_exact(&self, one: bool) -> bool {
        match self {
            Num::Real(q) => {
                if one {
                    q.is_one()
                } else {
                    q.is_zero()
                }
            }
            _ => false,
        }
    }

    pub(super) fn add(&self, other: &Num) -> Num {
        match (self, other) {
            (Num::Real(a), Num::Real(b)) => Num::Real(a + b),
            _ if self.is_ball() || other.is_ball() => {
                if self.is_exact(false) {
                    return other.clone();
                }
                if other.is_exact(false) {
                    return self.clone();
                }
                self.ball_op(other, Ball::add, ComplexBall::add)
            }
            _ => Num::from_complex(&self.complex() + &other.complex()),
        }
    }

    pub(super) fn mul(&self, other: &Num) -> Num {
        match (self, other) {
            (Num::Real(a), Num::Real(b)) => Num::Real(a * b),
            _ if self.is_ball() || other.is_ball() => {
                if self.is_exact(true) {
                    return other.clone();
                }
                if other.is_exact(true) {
                    return self.clone();
                }
                self.ball_op(other, Ball::mul, ComplexBall::mul)
            }
            _ => Num::from_complex(&self.complex() * &other.complex()),
        }
    }

    /// `self^n`; `None` for zero to a negative power. A ball containing 0
    /// to a negative power is the indeterminate ball, which becomes
    /// `Indeterminate`.
    pub(super) fn pow(&self, n: i64) -> Option<Num> {
        if self.is_zero() && n < 0 {
            return None;
        }
        Some(match self {
            Num::Real(q) => Num::Real(q.pow(n)),
            Num::Complex(z) => Num::from_complex(z.pow(n)),
            Num::Ball(b) => {
                let bits = ball_bits(&[self]).expect("a ball");
                Num::Ball(Box::new(b.pow(n, bits + GUARD).round_to_accuracy(bits)))
            }
            Num::ComplexBall(z) => {
                let bits = ball_bits(&[self]).expect("a ball");
                Num::from_complex_ball(round_complex(&z.pow_int(n, bits + GUARD), bits))
            }
        })
    }

    /// `self^exponent` when either is a ball: `x^y`, the real root of odd
    /// order in the Real domain for a negative `x` and an exact `p/q` of odd
    /// `q` (O7), and the principal complex power when the real one has no
    /// value. `None` for 0 to a negative or 0th power, which the caller
    /// decides.
    pub(super) fn ball_power(&self, exponent: &Num, real_domain: bool) -> Option<Num> {
        let bits = ball_bits(&[self, exponent])?;
        let working = bits + GUARD;
        if self.is_zero() && !matches!(exponent.sign(), Some(Ordering::Greater)) {
            return None;
        }
        if let (Some(x), Some(y)) = (self.to_ball(working), exponent.to_ball(working)) {
            if let (true, Some(q)) = (real_domain, exponent.real()) {
                let odd = q.denom().unsigned_abs().bit(0);
                if odd && !q.is_integer() && x.sign() == Some(Ordering::Less) {
                    let root = x.root(u64::try_from(q.denom()).ok()?, working);
                    let power = root.pow(q.numer().to_i64()?, working);
                    return Some(Num::Ball(Box::new(power.round_to_accuracy(bits))));
                }
            }
            let power = x.pow_ball(&y, working);
            if power.is_finite() || x.sign() != Some(Ordering::Less) {
                return Some(Num::Ball(Box::new(power.round_to_accuracy(bits))));
            }
        }
        let (x, y) = (
            self.to_complex_ball(working),
            exponent.to_complex_ball(working),
        );
        Some(Num::from_complex_ball(round_complex(
            &x.pow(&y, working),
            bits,
        )))
    }

    /// Whether this is a real value known to be negative.
    pub(super) fn is_negative(&self) -> bool {
        self.sign() == Some(Ordering::Less)
    }

    /// The bits ball arithmetic on this value keeps ([`ball_bits`]); the
    /// least for an exact value.
    pub(super) fn precision(&self) -> u64 {
        ball_bits(&[self]).unwrap_or(MIN_BITS)
    }

    /// The sign of a real value, if known.
    fn sign(&self) -> Option<Ordering> {
        match self {
            Num::Real(q) => Some(q.cmp(&Rational::from(0))),
            Num::Ball(b) => b.sign(),
            Num::Complex(_) | Num::ComplexBall(_) => None,
        }
    }

    /// An estimate of the bits of `self^n`, to refuse powers too large to
    /// hold.
    pub(super) fn power_bits(&self, n: i64) -> u64 {
        let bits = |q: &Rational| q.numer().bit_length().max(q.denom().bit_length());
        let base = match self {
            Num::Real(q) => bits(q),
            Num::Complex(z) => bits(z.re()).max(bits(z.im())) + 1,
            // A ball's power keeps its precision.
            Num::Ball(_) | Num::ComplexBall(_) => return 0,
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
            Num::Ball(b) => match b.sign()? {
                Ordering::Greater => Some(Num::one()),
                Ordering::Less => Some(Num::Real(Rational::from(-1))),
                Ordering::Equal => None,
            },
            Num::ComplexBall(_) => None,
        }
    }

    /// The expression of the value; `Indeterminate` for a ball that is not
    /// finite, which is no value.
    pub(super) fn to_expr(&self) -> Expr {
        match self {
            Num::Real(q) => Expr::from(q),
            Num::Complex(z) => Expr::from(z.clone()),
            Num::Ball(b) if b.is_finite() => Expr::from(&**b),
            Num::ComplexBall(z) if z.is_finite() => Expr::from(&**z),
            Num::Ball(_) | Num::ComplexBall(_) => Symbol::INDETERMINATE.into(),
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
