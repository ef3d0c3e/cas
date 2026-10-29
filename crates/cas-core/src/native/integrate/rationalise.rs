//! Substitutions that make the integrand a rational function of a new
//! variable (CAS-DOC-40), tried after the heuristics of CAS-199:
//!
//! - rational functions of exponentials `e^(a x + b)`, `a` rational, and of
//!   the hyperbolic functions, written with exponentials first (§2,
//!   CAS-433): `u = e^(g x)` for `g` the gcd of the `a`s, so that
//!   `e^(a x + b) = e^b u^(a/g)` with an integer `a/g`, and
//!   `dx = du/(g u)`. giac's `int(1/(exp(x)+exp(-x)),x)` is
//!   `atan(exp(x))`, here from `1/(u² + 1)`.
//! - rational functions `F(sin v, cos v)` of one linear argument
//!   `v = a x + b`, `tan`, `cot`, `sec` and `csc` written with them (§3,
//!   CAS-434): by Bioche's rules on `ω = F dv`, checked with `Together`,
//!   `u = cos v` when `ω` is unchanged by `v → −v`, `u = sin v` by
//!   `v → π − v`, `u = tan v` by `v → π + v`; otherwise Weierstrass's
//!   `t = tan(v/2)`. The powers of the other function are even, so a square
//!   root of `1 − u²` or `1 + u²` put for it cancels. `ArcTan[Tan[w]]` goes
//!   back as `w`: giac's `int(tan(x)^2,x)` is `tan(x)-x`.
//!
//! The rational function is integrated as an integrand of its own, a level
//! deeper; an antiderivative with an `Integrate` left fails, and the
//! substitution goes back into the antiderivative.

use cas_num::{Integer, Rational};

use super::integral_at;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Found = Result<Option<Expr>, EvalError>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// The new variable of a substitution `depth` deep.
fn fresh(depth: u8) -> Expr {
    Expr::from(Symbol::new(&format!("Integrate$u{depth}")))
}

/// An antiderivative of `f` by one of the module's substitutions.
pub(super) fn integral(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr, depth: u8) -> Found {
    if let Some(g) = exponentials(ev, e, f, x, depth)? {
        return Ok(Some(g));
    }
    trigonometric(ev, e, f, x, depth)
}

/// The rational function `r` of `u` integrated whole, `None` otherwise.
fn rational_integral(ev: &mut Eval<'_>, e: &Expr, r: &Expr, u: &Expr, depth: u8) -> Found {
    let Some(antiderivative) = integral_at(ev, e, r, u, depth + 1)? else {
        return Ok(None);
    };
    Ok((!antiderivative.contains(&Expr::from(Symbol::INTEGRATE))).then_some(antiderivative))
}

/// The gcd of rationals: the numerators', each over the lcm `l` of the
/// denominators, over `l`.
fn gcd(rs: &[Rational]) -> Rational {
    let l = rs.iter().fold(Integer::ONE, |l, r| l.lcm(r.denom()));
    let g = rs.iter().fold(Integer::ZERO, |g, r| {
        g.gcd(&(r.numer() * &l.div_exact(r.denom())))
    });
    Rational::new(g, l)
}

/// The exponent `w` of `Exp[w]` or `E^w`.
fn exponent(g: &Expr) -> Option<&Expr> {
    match (head(g), g.args()) {
        (Some(Symbol::EXP), [w]) => Some(w),
        (Some(Symbol::POWER), [b, w]) if b.as_symbol() == Some(Symbol::E) => Some(w),
        _ => None,
    }
}

/// The slope of `w` in `x`, a rational not 0, when `w` is linear so.
fn rational_slope(ev: &mut Eval<'_>, w: &Expr, x: &Expr) -> Result<Option<Rational>, EvalError> {
    let a = ev.eval(&Expr::call(Symbol::D, [w.clone(), x.clone()]))?;
    Ok(a.to_rational().filter(|a| !a.is_zero()))
}

/// Whether `x` is in `f` only in exponentials of linear arguments with
/// rational slopes, through sums, products and powers free of `x`; the
/// slopes collected.
fn exponential_parts(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    slopes: &mut Vec<Rational>,
) -> Result<bool, EvalError> {
    if !f.contains(x) {
        return Ok(true);
    }
    if let Some(w) = exponent(f) {
        let Some(a) = rational_slope(ev, w, x)? else {
            return Ok(false);
        };
        slopes.push(a);
        return Ok(true);
    }
    match (head(f), f.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), parts) => {
            for part in parts {
                if !exponential_parts(ev, part, x, slopes)? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
        (Some(Symbol::POWER), [base, n]) if !n.contains(x) => {
            exponential_parts(ev, base, x, slopes)
        }
        _ => Ok(false),
    }
}

/// `f` with each `e^(a x + b)` written `e^b u^(a/g)`.
fn in_powers(ev: &mut Eval<'_>, f: &Expr, x: &Expr, g: &Rational, u: &Expr) -> Found {
    if !f.contains(x) {
        return Ok(Some(f.clone()));
    }
    if let Some(w) = exponent(f) {
        let Some(a) = rational_slope(ev, w, x)? else {
            return Ok(None);
        };
        // w is linear: its constant is its value at 0.
        let b = ev.eval(&w.substitute(&[(x.clone(), Expr::from(0))]))?;
        let k = &a / g;
        return Ok(Some(Expr::from(Symbol::E).pow(b) * u.pow(Expr::from(k))));
    }
    let Some(h) = f.head() else {
        return Ok(None);
    };
    let mut args = Vec::with_capacity(f.args().len());
    for part in f.args() {
        let Some(part) = in_powers(ev, part, x, g, u)? else {
            return Ok(None);
        };
        args.push(part);
    }
    Ok(Some(Expr::apply(h, &args)))
}

/// Rational functions of exponentials and hyperbolic functions (§2): `u =
/// e^(g x)`.
fn exponentials(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr, depth: u8) -> Found {
    let hyperbolic = f
        .distinct()
        .any(|g| matches!(head(g), Some(Symbol::SINH | Symbol::COSH | Symbol::TANH)));
    let f = if hyperbolic {
        ev.eval(&Expr::call(Symbol::HYPERBOLIC_TO_EXP, [f.clone()]))?
    } else {
        f.clone()
    };
    let mut slopes = Vec::new();
    if !exponential_parts(ev, &f, x, &mut slopes)? || slopes.is_empty() {
        return Ok(None);
    }
    let g = gcd(&slopes);
    let u = fresh(depth);
    if f.contains(&u) {
        return Ok(None);
    }
    let Some(powers) = in_powers(ev, &f, x, &g, &u)? else {
        return Ok(None);
    };
    // dx = du/(g u).
    let rational = ev.eval(&(powers / (Expr::from(g.clone()) * u.clone())))?;
    if rational.contains(x) {
        return Ok(None);
    }
    let Some(antiderivative) = rational_integral(ev, e, &rational, &u, depth)? else {
        return Ok(None);
    };
    // Back as the integrand wrote its exponentials: giac's atan(exp(x)).
    let w = Expr::from(g) * x.clone();
    let back = if f.distinct().any(|part| head(part) == Some(Symbol::EXP)) {
        Expr::call(Symbol::EXP, [w])
    } else {
        Expr::from(Symbol::E).pow(w)
    };
    Ok(Some(antiderivative.substitute(&[(u, back)])))
}

/// Whether `x` is in `f` only in sines, cosines, tangents and their
/// reciprocals, through sums, products and powers free of `x`; their
/// arguments collected.
fn trigonometric_parts(f: &Expr, x: &Expr, arguments: &mut Vec<Expr>) -> bool {
    if !f.contains(x) {
        return true;
    }
    match (head(f), f.args()) {
        (
            Some(Symbol::SIN | Symbol::COS | Symbol::TAN | Symbol::COT | Symbol::SEC | Symbol::CSC),
            [v],
        ) => {
            arguments.push(v.clone());
            true
        }
        (Some(Symbol::PLUS | Symbol::TIMES), parts) => {
            parts.iter().all(|p| trigonometric_parts(p, x, arguments))
        }
        (Some(Symbol::POWER), [base, n]) if !n.contains(x) => {
            trigonometric_parts(base, x, arguments)
        }
        _ => false,
    }
}

/// `f` with its trigonometric functions of `v` written with `s = sin v`
/// and `c = cos v`.
fn in_sines(f: &Expr, v: &Expr, s: &Expr, c: &Expr) -> Expr {
    let one = || Expr::from(1);
    match (head(f), f.args()) {
        (Some(function), [w]) if w == v => match function {
            Symbol::SIN => return s.clone(),
            Symbol::COS => return c.clone(),
            Symbol::TAN => return s.clone() / c.clone(),
            Symbol::COT => return c.clone() / s.clone(),
            Symbol::SEC => return one() / c.clone(),
            Symbol::CSC => return one() / s.clone(),
            _ => {}
        },
        _ => {}
    }
    match f.head() {
        Some(h) if f.contains(v) => {
            let args: Vec<Expr> = f.args().iter().map(|a| in_sines(a, v, s, c)).collect();
            Expr::apply(h, &args)
        }
        _ => f.clone(),
    }
}

/// Whether `a` and `b` are the same rational function of their
/// generators.
fn same(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    let difference = ev.eval(&Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]))?;
    Ok(difference == Expr::from(0))
}

/// Whether `r` has a root, a power with an exponent that is not an
/// integer, of something with `u`.
fn has_root(r: &Expr, u: &Expr) -> bool {
    r.distinct().any(|part| {
        head(part) == Some(Symbol::POWER)
            && part.args()[0].contains(u)
            && !part.args()[1].to_rational().is_some_and(|q| q.is_integer())
    })
}

/// `ArcTan[Tan[w]]` as `w`: they differ by a constant on each interval where
/// `Tan[w]` is continuous, as antiderivatives may.
fn untangled(e: &Expr) -> Expr {
    if let (Some(Symbol::ARC_TAN), [t]) = (head(e), e.args())
        && let (Some(Symbol::TAN), [w]) = (head(t), t.args())
    {
        return w.clone();
    }
    match e.head() {
        Some(h) if !e.args().is_empty() => {
            let args: Vec<Expr> = e.args().iter().map(untangled).collect();
            Expr::apply(h, &args)
        }
        _ => e.clone(),
    }
}

/// The substitutions of §3.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rule {
    Cosine,
    Sine,
    Tangent,
    HalfTangent,
}

/// Bioche's rule for `ω = F(s, c) dv`, or Weierstrass's.
fn rule(ev: &mut Eval<'_>, big_f: &Expr, s: &Expr, c: &Expr) -> Result<Rule, EvalError> {
    let flip = |e: &Expr, sine: bool, cosine: bool| {
        let mut rules = Vec::new();
        if sine {
            rules.push((s.clone(), -s.clone()));
        }
        if cosine {
            rules.push((c.clone(), -c.clone()));
        }
        e.substitute(&rules)
    };
    // v → −v: F(−s, c) = −F; v → π − v: F(s, −c) = −F; v → π + v: F(−s, −c) = F.
    if same(ev, &flip(big_f, true, false), &-big_f.clone())? {
        return Ok(Rule::Cosine);
    }
    if same(ev, &flip(big_f, false, true), &-big_f.clone())? {
        return Ok(Rule::Sine);
    }
    if same(ev, &flip(big_f, true, true), big_f)? {
        return Ok(Rule::Tangent);
    }
    Ok(Rule::HalfTangent)
}

/// Rational functions of sines and cosines of one linear argument (§3).
fn trigonometric(ev: &mut Eval<'_>, e: &Expr, f: &Expr, x: &Expr, depth: u8) -> Found {
    let mut arguments = Vec::new();
    if !trigonometric_parts(f, x, &mut arguments) {
        return Ok(None);
    }
    let Some(v) = arguments.first().cloned() else {
        return Ok(None);
    };
    if arguments.iter().any(|w| *w != v) {
        return Ok(None);
    }
    let a = ev.eval(&Expr::call(Symbol::D, [v.clone(), x.clone()]))?;
    if a.contains(x) || a == Expr::from(0) {
        return Ok(None);
    }
    let s = Expr::from(Symbol::new(&format!("Integrate$s{depth}")));
    let c = Expr::from(Symbol::new(&format!("Integrate$c{depth}")));
    let u = fresh(depth);
    if [&s, &c, &u].iter().any(|y| f.contains(y)) {
        return Ok(None);
    }
    let big_f = ev.eval(&in_sines(f, &v, &s, &c))?;
    if big_f.contains(x) {
        return Ok(None);
    }
    let one = || Expr::from(1);
    let half = || Expr::from(Rational::new(Integer::ONE, Integer::from(2)));
    let root = |q: Expr| q.pow(half());
    let rule = rule(ev, &big_f, &s, &c)?;
    let (r, back) = match rule {
        // du = −s dv, s² = 1 − u².
        Rule::Cosine => (
            (-big_f / s.clone()).substitute(&[(s.clone(), root(one() - u.pow(2))), (c, u.clone())]),
            Expr::call(Symbol::COS, [v.clone()]),
        ),
        // du = c dv, c² = 1 − u².
        Rule::Sine => (
            (big_f / c.clone()).substitute(&[(c.clone(), root(one() - u.pow(2))), (s, u.clone())]),
            Expr::call(Symbol::SIN, [v.clone()]),
        ),
        // s = u c, c² = 1/(1 + u²), dv = du/(1 + u²).
        Rule::Tangent => {
            let c_u = one() / root(one() + u.pow(2));
            (
                big_f.substitute(&[(s, u.clone() * c_u.clone()), (c, c_u)]) / (one() + u.pow(2)),
                Expr::call(Symbol::TAN, [v.clone()]),
            )
        }
        // s = 2t/(1 + t²), c = (1 − t²)/(1 + t²), dv = 2 dt/(1 + t²).
        Rule::HalfTangent => {
            let denominator = one() + u.pow(2);
            (
                big_f.substitute(&[
                    (s, Expr::from(2) * u.clone() / denominator.clone()),
                    (c, (one() - u.pow(2)) / denominator.clone()),
                ]) * Expr::from(2)
                    / denominator,
                Expr::call(Symbol::TAN, [v.clone() * half()]),
            )
        }
    };
    let r = ev.eval(&Expr::call(Symbol::TOGETHER, [r]))?;
    if has_root(&r, &u) || r.contains(x) {
        return Ok(None);
    }
    let Some(antiderivative) = rational_integral(ev, e, &r, &u, depth)? else {
        return Ok(None);
    };
    // dv = a dx.
    let result = antiderivative.substitute(&[(u, back)]) / a;
    Ok(Some(untangled(&ev.eval(&result)?)))
}

#[cfg(test)]
mod tests;
