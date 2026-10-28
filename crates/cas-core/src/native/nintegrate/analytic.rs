//! Integrands on complex balls (CAS-DOC-26 §3, CAS-323): an integrand is
//! compiled once per precision into a tree of ball operations ([`Node`]),
//! its parts free of the variable into balls by `N`'s values. A tree is
//! evaluated at a complex ball in one of two [`Mode`]s:
//!
//! - **Analytic**: the analytic continuation of the real integrand,
//!   holomorphic on the ball wherever an answer comes.
//!   - Each pole and branch cut is kept away from its argument's ball:
//!     - 0 from a negative integer power's base;
//!     - `(−∞, 0]` from `Log`'s argument and from the base of another power;
//!     - `{iy : |y| ≥ 1}` from `ArcTan` and `ArcSinh`;
//!     - `(−∞, −1] ∪ [1, ∞)` from `ArcSin`, `ArcCos` and `ArcTanh`;
//!     - `(−∞, 1]` from `ArcCosh`;
//!     - the zeros of the denominators of `Tan` and `Tanh`.
//!   - A composition is holomorphic where each part is, on the enclosure of
//!     its argument.
//!   - `Abs`, `Sign` and `HeavisideTheta` give nothing, until
//!     [`Node::resolved`] continues them from a real segment by the sign of
//!     their argument there: `±u`, `±1`, 1 or 0.
//! - **Real**: the values of the integrand itself, enclosed, for the crude
//!   bound of §4.
//!   - The functions are the ball functions, which cover both sides of a
//!     cut.
//!   - `Abs` is the modulus, and `Sign` and `HeavisideTheta` are their
//!     values, or `[−1, 1]` and `[0, 1]` where the argument's sign is not
//!     known.
//!
//! A ball that is not finite, as at a pole, gives nothing in either mode.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the tree: f in x, b^n, the ball z and the value v"
)]

use std::cmp::Ordering;

use cas_num::ball::{Ball, ComplexBall, Float, Mag};
use cas_num::{Integer, Rational};

use crate::context::Context;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;

/// The largest integer exponent compiled.
const MAX_EXPONENT: i64 = 1 << 16;

/// A switch of the real integrand.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Switch {
    Abs,
    Sign,
    Step,
}

/// A function of one argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Function {
    Exp,
    Log,
    Sin,
    Cos,
    Tan,
    Sinh,
    Cosh,
    Tanh,
    ArcTan,
    ArcSin,
    ArcCos,
    ArcSinh,
    ArcCosh,
    ArcTanh,
}

impl Function {
    fn of(s: Symbol) -> Option<Function> {
        Some(match s {
            Symbol::EXP => Function::Exp,
            Symbol::LOG => Function::Log,
            Symbol::SIN => Function::Sin,
            Symbol::COS => Function::Cos,
            Symbol::TAN => Function::Tan,
            Symbol::SINH => Function::Sinh,
            Symbol::COSH => Function::Cosh,
            Symbol::TANH => Function::Tanh,
            Symbol::ARC_TAN => Function::ArcTan,
            Symbol::ARC_SIN => Function::ArcSin,
            Symbol::ARC_COS => Function::ArcCos,
            Symbol::ARC_SINH => Function::ArcSinh,
            Symbol::ARC_COSH => Function::ArcCosh,
            Symbol::ARC_TANH => Function::ArcTanh,
            _ => return None,
        })
    }
}

/// How a tree is evaluated (the module's two modes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Mode {
    Analytic,
    Real,
}

/// An integrand as ball operations.
#[derive(Clone, Debug)]
pub(super) enum Node {
    Variable,
    Constant(ComplexBall),
    Sum(Vec<Node>),
    Product(Vec<Node>),
    /// `b^n` for an integer `n`.
    IntegerPower(Box<Node>, i64),
    /// The principal `b^w` for a constant `w`.
    Power(Box<Node>, ComplexBall),
    /// `b^g = e^(g log b)`, principal.
    Exponential(Box<Node>, Box<Node>),
    Apply(Function, Box<Node>),
    Switch(Switch, Box<Node>),
}

/// The ball of `f`, free of the variable, by `N`'s values; `None` when it
/// is no number.
fn constant(ctx: &Context, f: &Expr, prec: u64) -> Option<ComplexBall> {
    let value = approximation(ctx, f, prec).ok()?;
    let ball = match value.view() {
        View::Ball(b) => ComplexBall::from_real(b.clone()),
        View::ComplexBall(z) => z.clone(),
        _ => return None,
    };
    ball.is_finite().then_some(ball)
}

fn real(q: &Rational, prec: u64) -> ComplexBall {
    ComplexBall::from_real(Ball::from_rational(q, prec))
}

fn q(n: i64) -> Rational {
    Rational::from(Integer::from(n))
}

/// Whether `z` may meet the ray `(−∞, c]` (`left`) or `[c, ∞)`.
fn meets_ray(z: &ComplexBall, c: i64, left: bool, prec: u64) -> bool {
    if !z.im().contains_rational(&Rational::ZERO) {
        return false;
    }
    let side = z
        .re()
        .sub(&Ball::from_integer(&Integer::from(c)), prec)
        .sign();
    if left {
        side != Some(Ordering::Greater)
    } else {
        side != Some(Ordering::Less)
    }
}

/// Whether `z` may meet `{iy : |y| ≥ 1}`.
fn meets_imaginary_rays(z: &ComplexBall, prec: u64) -> bool {
    z.re().contains_rational(&Rational::ZERO)
        && z.im()
            .abs()
            .sub(&Ball::from_integer(&Integer::from(1)), prec)
            .sign()
            != Some(Ordering::Less)
}

impl Node {
    /// `f` in the variable `x` as a tree at `prec` bits; `None` for a head
    /// without a ball function here, or a part free of `x` that is no
    /// number.
    pub(super) fn compile(ctx: &Context, f: &Expr, x: &Expr, prec: u64) -> Option<Node> {
        if f == x {
            return Some(Node::Variable);
        }
        if !f.contains(x) {
            return constant(ctx, f, prec).map(Node::Constant);
        }
        let head = f.head()?.as_symbol()?;
        let args = f.args();
        let compile = |g: &Expr| Node::compile(ctx, g, x, prec).map(Box::new);
        let all = || {
            args.iter()
                .map(|a| Node::compile(ctx, a, x, prec))
                .collect::<Option<Vec<_>>>()
        };
        Some(match (head, args) {
            (Symbol::PLUS, _) => Node::Sum(all()?),
            (Symbol::TIMES, _) => Node::Product(all()?),
            (Symbol::POWER, [b, g]) if b.as_symbol() == Some(Symbol::E) => {
                Node::Apply(Function::Exp, compile(g)?)
            }
            (Symbol::POWER, [b, n]) if !n.contains(x) => {
                let integer = n
                    .to_rational()
                    .filter(Rational::is_integer)
                    .and_then(|k| k.numer().to_i64());
                match integer {
                    Some(k) if k.abs() <= MAX_EXPONENT => Node::IntegerPower(compile(b)?, k),
                    Some(_) => return None,
                    None => Node::Power(compile(b)?, constant(ctx, n, prec)?),
                }
            }
            (Symbol::POWER, [b, g]) => Node::Exponential(compile(b)?, compile(g)?),
            (s, [g]) => {
                let switch = match s {
                    Symbol::ABS => Some(Switch::Abs),
                    Symbol::SIGN => Some(Switch::Sign),
                    Symbol::HEAVISIDE_THETA => Some(Switch::Step),
                    _ => None,
                };
                match (Function::of(s), switch) {
                    (Some(function), _) => Node::Apply(function, compile(g)?),
                    (None, Some(switch)) => Node::Switch(switch, compile(g)?),
                    (None, None) => return None,
                }
            }
            _ => return None,
        })
    }

    /// Whether the tree has a switch.
    pub(super) fn has_switch(&self) -> bool {
        match self {
            Node::Variable | Node::Constant(_) => false,
            Node::Sum(parts) | Node::Product(parts) => parts.iter().any(Node::has_switch),
            Node::IntegerPower(b, _) | Node::Power(b, _) | Node::Apply(_, b) => b.has_switch(),
            Node::Exponential(b, g) => b.has_switch() || g.has_switch(),
            Node::Switch(..) => true,
        }
    }

    /// The tree with each switch continued from the real segment
    /// `segment`, by the sign of its argument there; `None` where an
    /// argument is not real there, or its sign is not known.
    pub(super) fn resolved(&self, segment: &ComplexBall, prec: u64) -> Option<Node> {
        let each = |parts: &[Node]| {
            parts
                .iter()
                .map(|p| p.resolved(segment, prec))
                .collect::<Option<Vec<_>>>()
        };
        Some(match self {
            Node::Variable | Node::Constant(_) => self.clone(),
            Node::Sum(parts) => Node::Sum(each(parts)?),
            Node::Product(parts) => Node::Product(each(parts)?),
            Node::IntegerPower(b, n) => {
                Node::IntegerPower(Box::new(b.resolved(segment, prec)?), *n)
            }
            Node::Power(b, w) => Node::Power(Box::new(b.resolved(segment, prec)?), w.clone()),
            Node::Exponential(b, g) => Node::Exponential(
                Box::new(b.resolved(segment, prec)?),
                Box::new(g.resolved(segment, prec)?),
            ),
            Node::Apply(f, g) => Node::Apply(*f, Box::new(g.resolved(segment, prec)?)),
            Node::Switch(s, u) => {
                let u = u.resolved(segment, prec)?;
                let value = u.eval(segment, Mode::Analytic, prec)?;
                if !value.is_real() {
                    return None;
                }
                let positive = match value.re().sign()? {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal => return None,
                };
                let sign = if positive { 1 } else { -1 };
                match s {
                    Switch::Abs if positive => u,
                    Switch::Abs => Node::Product(vec![Node::Constant(real(&q(-1), prec)), u]),
                    Switch::Sign => Node::Constant(real(&q(sign), prec)),
                    Switch::Step => Node::Constant(real(&q(i64::from(positive)), prec)),
                }
            }
        })
    }

    /// The tree at `z` in `mode` (the module's rules); `None` where the
    /// mode's conditions fail or the value is not finite.
    pub(super) fn eval(&self, z: &ComplexBall, mode: Mode, prec: u64) -> Option<ComplexBall> {
        let analytic = mode == Mode::Analytic;
        let value = match self {
            Node::Variable => z.clone(),
            Node::Constant(c) => c.clone(),
            Node::Sum(parts) => {
                let mut sum = ComplexBall::ZERO;
                for p in parts {
                    sum = sum.add(&p.eval(z, mode, prec)?, prec);
                }
                sum
            }
            Node::Product(parts) => {
                let mut product = real(&q(1), prec);
                for p in parts {
                    product = product.mul(&p.eval(z, mode, prec)?, prec);
                }
                product
            }
            Node::IntegerPower(b, n) => {
                let v = b.eval(z, mode, prec)?;
                if *n < 0 && v.may_be_zero() {
                    return None;
                }
                v.pow_int(*n, prec)
            }
            Node::Power(b, w) => {
                let v = b.eval(z, mode, prec)?;
                if analytic && meets_ray(&v, 0, true, prec) {
                    return None;
                }
                v.pow(w, prec)
            }
            Node::Exponential(b, g) => {
                let v = b.eval(z, mode, prec)?;
                if analytic && meets_ray(&v, 0, true, prec) {
                    return None;
                }
                g.eval(z, mode, prec)?.mul(&v.log(prec), prec).exp(prec)
            }
            Node::Apply(f, g) => apply(*f, &g.eval(z, mode, prec)?, analytic, prec)?,
            Node::Switch(s, u) => {
                if analytic {
                    return None;
                }
                switch(*s, &u.eval(z, mode, prec)?, prec)?
            }
        };
        value.is_finite().then_some(value)
    }
}

/// `f(v)`, its cuts kept away from `v` when `analytic`.
fn apply(f: Function, v: &ComplexBall, analytic: bool, prec: u64) -> Option<ComplexBall> {
    let away = match f {
        Function::Log => !meets_ray(v, 0, true, prec),
        Function::ArcTan | Function::ArcSinh => !meets_imaginary_rays(v, prec),
        Function::ArcSin | Function::ArcCos | Function::ArcTanh => {
            !meets_ray(v, -1, true, prec) && !meets_ray(v, 1, false, prec)
        }
        Function::ArcCosh => !meets_ray(v, 1, true, prec),
        _ => true,
    };
    if analytic && !away {
        return None;
    }
    Some(match f {
        Function::Exp => v.exp(prec),
        Function::Log => v.log(prec),
        Function::Sin => v.sin(prec),
        Function::Cos => v.cos(prec),
        Function::Tan => {
            let c = v.cos(prec);
            if c.may_be_zero() {
                return None;
            }
            v.sin(prec).div(&c, prec)
        }
        Function::Sinh => v.sinh(prec),
        Function::Cosh => v.cosh(prec),
        Function::Tanh => {
            let c = v.cosh(prec);
            if c.may_be_zero() {
                return None;
            }
            v.sinh(prec).div(&c, prec)
        }
        Function::ArcTan => v.atan(prec),
        Function::ArcSin => v.asin(prec),
        Function::ArcCos => v.acos(prec),
        Function::ArcSinh => v.asinh(prec),
        Function::ArcCosh => v.acosh(prec),
        Function::ArcTanh => v.atanh(prec),
    })
}

/// A switch of the real integrand at `v` (the Real mode).
fn switch(s: Switch, v: &ComplexBall, prec: u64) -> Option<ComplexBall> {
    if s == Switch::Abs {
        return Some(ComplexBall::from_real(if v.is_real() {
            v.re().abs()
        } else {
            v.abs(prec)
        }));
    }
    if !v.is_real() {
        return None;
    }
    let (mid, rad) = match (s, v.re().sign()) {
        (Switch::Sign, Some(o)) => (q(o as i64), 0),
        (Switch::Sign, None) => (q(0), 1),
        (_, Some(Ordering::Greater)) => (q(1), 0),
        (_, Some(Ordering::Less)) => (q(0), 0),
        (_, Some(Ordering::Equal)) => (Rational::new(Integer::from(1), Integer::from(2)), 0),
        (_, None) => (Rational::new(Integer::from(1), Integer::from(2)), 1),
    };
    let ball = Ball::from_rational(&mid, prec);
    Some(ComplexBall::from_real(if rad == 0 {
        ball
    } else if s == Switch::Sign {
        Ball::new(Float::ZERO, Mag::from_integer(&Integer::from(1)))
    } else {
        ball.add_error(&Mag::from_integer(&Integer::from(1)).mul_2exp(&Integer::from(-1)))
    }))
}

#[cfg(test)]
mod tests;
