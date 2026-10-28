//! Right sides of differential equations (CAS-DOC-35 §3, CAS-391): `f(t,
//! y₁, …, yₙ)` compiled once into a tree of real ball operations, its parts
//! free of the variables into balls by `N`'s values, evaluated on balls,
//! for the a priori enclosure, and on jets, for the Taylor coefficients
//! (§4).
//!
//! The heads are those of `NIntegrate`'s integrands but its switches:
//! `Abs`, `Sign` and steps are refused, as the method needs a right side
//! analytic along the path. A complex constant is refused too. A ball or a
//! jet outside a function's domain, as the logarithm of a ball meeting 0,
//! gives nothing, and the step that met it is shortened, or the call
//! stays.

use cas_num::ball::{Ball, Jet};
use cas_num::{Integer, Rational};

use crate::context::Context;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;

/// The largest integer exponent compiled.
const MAX_EXPONENT: i64 = 1 << 16;

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

    fn ball(self, x: &Ball, prec: u64) -> Ball {
        match self {
            Function::Exp => x.exp(prec),
            Function::Log => x.log(prec),
            Function::Sin => x.sin(prec),
            Function::Cos => x.cos(prec),
            Function::Tan => x.tan(prec),
            Function::Sinh => x.sinh(prec),
            Function::Cosh => x.cosh(prec),
            Function::Tanh => x.tanh(prec),
            Function::ArcTan => x.atan(prec),
            Function::ArcSin => x.asin(prec),
            Function::ArcCos => x.acos(prec),
            Function::ArcSinh => x.asinh(prec),
            Function::ArcCosh => x.acosh(prec),
            Function::ArcTanh => x.atanh(prec),
        }
    }

    fn jet(self, x: &Jet, prec: u64) -> Jet {
        match self {
            Function::Exp => x.exp(prec),
            Function::Log => x.log(prec),
            Function::Sin => x.sin_cos(prec).0,
            Function::Cos => x.sin_cos(prec).1,
            Function::Tan => x.tan(prec),
            Function::Sinh => x.sinh_cosh(prec).0,
            Function::Cosh => x.sinh_cosh(prec).1,
            Function::Tanh => x.tanh(prec),
            Function::ArcTan => x.atan(prec),
            Function::ArcSin => x.asin(prec),
            Function::ArcCos => x.acos(prec),
            Function::ArcSinh => x.asinh(prec),
            Function::ArcCosh => x.acosh(prec),
            Function::ArcTanh => x.atanh(prec),
        }
    }
}

/// A right side's tree.
#[derive(Clone, Debug)]
pub(super) enum Node {
    /// `t` for 0, `y_k` for `k`.
    Variable(usize),
    Constant(Ball),
    Sum(Vec<Node>),
    Product(Vec<Node>),
    IntegerPower(Box<Node>, i64),
    /// `b^w` for a real ball `w`, `b` positive.
    Power(Box<Node>, Ball),
    /// `b^g = e^(g log b)`, `b` positive.
    Exponential(Box<Node>, Box<Node>),
    Apply(Function, Box<Node>),
}

/// The real ball of `f`, free of the variables, by `N`'s values.
fn constant(ctx: &Context, f: &Expr, prec: u64) -> Option<Ball> {
    let value = approximation(ctx, f, prec).ok()?;
    match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => value.to_rational().map(|q| Ball::from_rational(&q, prec)),
    }
}

/// Whether the ball is above 0, as a power's base must be.
fn positive(b: &Ball) -> bool {
    b.sign() == Some(std::cmp::Ordering::Greater)
}

/// A ball's value, `None` where it is not finite.
fn finite(b: Ball) -> Option<Ball> {
    b.is_finite().then_some(b)
}

/// A jet's value, `None` where a coefficient is not finite.
fn finite_jet(j: Jet) -> Option<Jet> {
    j.coeffs().iter().all(Ball::is_finite).then_some(j)
}

fn half() -> Rational {
    Rational::new(Integer::ONE, Integer::from(2))
}

impl Node {
    /// `f` in `variables`, `t` first, as a tree at `prec` bits; `None` for
    /// a head without a ball function here, or a part free of the
    /// variables that is no real number.
    pub(super) fn compile(ctx: &Context, f: &Expr, variables: &[Expr], prec: u64) -> Option<Node> {
        if let Some(k) = variables.iter().position(|v| v == f) {
            return Some(Node::Variable(k));
        }
        if !variables.iter().any(|v| f.contains(v)) {
            return constant(ctx, f, prec).map(Node::Constant);
        }
        let head = f.head()?.as_symbol()?;
        let args = f.args();
        let compile = |g: &Expr| Node::compile(ctx, g, variables, prec).map(Box::new);
        let all = || {
            args.iter()
                .map(|a| Node::compile(ctx, a, variables, prec))
                .collect::<Option<Vec<_>>>()
        };
        let free = |g: &Expr| !variables.iter().any(|v| g.contains(v));
        Some(match (head, args) {
            (Symbol::PLUS, _) => Node::Sum(all()?),
            (Symbol::TIMES, _) => Node::Product(all()?),
            (Symbol::POWER, [b, g]) if b.as_symbol() == Some(Symbol::E) => {
                Node::Apply(Function::Exp, compile(g)?)
            }
            (Symbol::POWER, [b, n]) if free(n) => {
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
            (s, [g]) => Node::Apply(Function::of(s)?, compile(g)?),
            _ => return None,
        })
    }

    /// The value at the balls `at`, `t` first; `None` outside a domain.
    pub(super) fn ball(&self, at: &[Ball], prec: u64) -> Option<Ball> {
        Some(match self {
            Node::Variable(k) => at.get(*k)?.clone(),
            Node::Constant(c) => c.clone(),
            Node::Sum(parts) => {
                let mut acc = Ball::ZERO;
                for p in parts {
                    acc = acc.add(&p.ball(at, prec)?, prec);
                }
                acc
            }
            Node::Product(parts) => {
                let mut acc = Ball::from_integer(&Integer::ONE);
                for p in parts {
                    acc = acc.mul(&p.ball(at, prec)?, prec);
                }
                acc
            }
            Node::IntegerPower(b, n) => finite(b.ball(at, prec)?.pow(*n, prec))?,
            Node::Power(b, w) => {
                let base = b.ball(at, prec)?;
                if !positive(&base) {
                    return None;
                }
                finite(base.log(prec).mul(w, prec).exp(prec))?
            }
            Node::Exponential(b, g) => {
                let base = b.ball(at, prec)?;
                if !positive(&base) {
                    return None;
                }
                let exponent = g.ball(at, prec)?;
                finite(base.log(prec).mul(&exponent, prec).exp(prec))?
            }
            Node::Apply(f, g) => finite(f.ball(&g.ball(at, prec)?, prec))?,
        })
    }

    /// The jet at the jets `at`, `t` first, all of one order; `None`
    /// outside a domain.
    pub(super) fn jet(&self, at: &[Jet], prec: u64) -> Option<Jet> {
        let n = at.first()?.len();
        Some(match self {
            Node::Variable(k) => at.get(*k)?.clone(),
            Node::Constant(c) => Jet::constant(c.clone(), n),
            Node::Sum(parts) => {
                let mut acc = Jet::constant(Ball::ZERO, n);
                for p in parts {
                    acc = acc.add(&p.jet(at, prec)?, prec);
                }
                acc
            }
            Node::Product(parts) => {
                let mut acc = Jet::constant(Ball::from_integer(&Integer::ONE), n);
                for p in parts {
                    acc = acc.mul(&p.jet(at, prec)?, prec);
                }
                acc
            }
            Node::IntegerPower(b, k) => finite_jet(b.jet(at, prec)?.powi(*k, prec))?,
            Node::Power(b, w) => {
                let base = b.jet(at, prec)?;
                if !positive(&base.coeffs()[0]) {
                    return None;
                }
                if w.is_exact() && w.mid().to_rational() == half() {
                    finite_jet(base.sqrt(prec))?
                } else {
                    finite_jet(base.pow(w, prec))?
                }
            }
            Node::Exponential(b, g) => {
                let base = b.jet(at, prec)?;
                if !positive(&base.coeffs()[0]) {
                    return None;
                }
                let exponent = g.jet(at, prec)?;
                finite_jet(exponent.mul(&base.log(prec), prec).exp(prec))?
            }
            Node::Apply(f, g) => finite_jet(f.jet(&g.jet(at, prec)?, prec))?,
        })
    }
}

#[cfg(test)]
mod tests;
