//! Substitutions that make the integrand a rational function of a new
//! variable (CAS-DOC-40), tried after the heuristics of CAS-199:
//!
//! - rational functions of exponentials `e^(a x + b)`, `a` rational, and of
//!   the hyperbolic functions, written with exponentials first (§2,
//!   CAS-433): `u = e^(g x)` for `g` the gcd of the `a`s, so that
//!   `e^(a x + b) = e^b u^(a/g)` with an integer `a/g`, and
//!   `dx = du/(g u)`. giac's `int(1/(exp(x)+exp(-x)),x)` is
//!   `atan(exp(x))`, here from `1/(u² + 1)`.
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
    exponentials(ev, e, f, x, depth)
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
    let Some(antiderivative) = integral_at(ev, e, &rational, &u, depth + 1)? else {
        return Ok(None);
    };
    if antiderivative.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(None);
    }
    // Back as the integrand wrote its exponentials: giac's atan(exp(x)).
    let w = Expr::from(g) * x.clone();
    let back = if f.distinct().any(|part| head(part) == Some(Symbol::EXP)) {
        Expr::call(Symbol::EXP, [w])
    } else {
        Expr::from(Symbol::E).pow(w)
    };
    Ok(Some(antiderivative.substitute(&[(u, back)])))
}

#[cfg(test)]
mod tests;
