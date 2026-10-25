//! `Limit` (CAS-DOC-10 §4, §6, CAS-205): the limit of `f` as `x` tends to
//! `a`, by series. For a side, `x = a ± t` or `±1/t` with `t → 0⁺`, and the
//! leading term `c t^q` of the series of `f` in `t` (CAS-203) decides: 0
//! for `q > 0`, `c` for `q = 0`, and for `q < 0` `Infinity` times the sign
//! of `c`, `ComplexInfinity` when it has none that the assumptions or its
//! ball tell (giac's unsigned `infinity` for `a x` at infinity). `log t` in
//! `c` tends to `-Infinity`: `c` a polynomial in it, its leading term in
//! `log t` beats any constant and loses to any power of `t` (`x ln x` at
//! `0⁺` is 0); `c` at `q = 0`, the limit of `c` as `log t → −∞`.
//!
//! Two sides: the same value, that value; infinities of opposite signs,
//! `ComplexInfinity` (giac's `infinity` for `1/x` at 0); anything else
//! stays, where giac answers `undef`. At infinity there is one side.
//!
//! A leading coefficient that is a number `N` cannot tell from 0 stays the
//! call, with a warning: a 0 that evaluation does not see would give a
//! wrong limit (§6).
//!
//! The forms: `Limit[f, x, a]`, `Limit[f, x, a, dir]` (giac's
//! `limit(f, x, a, dir)`), `Limit[f, x == a]` and `Limit[f, x == a, dir]`
//! (giac's `limit(f, x=a, dir)`), and Mathematica's `Limit[f, x -> a]`;
//! `dir` 1 from the right, -1 from the left, 0 for both sides.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of limits: f as x tends to a, c t^q, from a side d"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::approximation;
use super::series::{Expander, Point, Series, expansion, fresh, log_t};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The orders tried for a leading term, each until one shows.
const ORDERS: [i64; 4] = [2, 4, 8, 16];

/// How many logarithms deep a leading coefficient is followed (§4.2).
const MAX_LOGARITHMS: u32 = 4;

/// A limit from one side.
#[derive(Clone, Debug, PartialEq)]
enum Value {
    Finite(Expr),
    /// `Infinity` (`true`) or `-Infinity`.
    Infinite(bool),
    Unsigned,
}

impl Value {
    fn expr(self) -> Expr {
        match self {
            Value::Finite(e) => e,
            Value::Infinite(true) => Expr::from(Symbol::INFINITY),
            Value::Infinite(false) => -Expr::from(Symbol::INFINITY),
            Value::Unsigned => Expr::from(Symbol::COMPLEX_INFINITY),
        }
    }
}

/// Why a side has no limit here.
enum Failure {
    /// No series, or its terms undecided: the call stays.
    Unknown,
    /// A leading coefficient that `N` cannot tell from 0.
    Undecided(Expr),
}

/// `f`, `x`, the point and the direction of a call.
fn parts(e: &Expr) -> Option<(&Expr, Expr, Point, i64)> {
    let direction = |d: Option<&Expr>| -> Option<i64> {
        match d {
            None => Some(0),
            Some(d) => d.as_i64().filter(|d| (-1..=1).contains(d)),
        }
    };
    let equation = |s: &Expr| -> Option<(Expr, Expr)> {
        match (s.head()?.as_symbol()?, s.args()) {
            (Symbol::EQUATION | Symbol::EQUAL | Symbol::RULE, [x, a]) => {
                Some((x.clone(), a.clone()))
            }
            _ => None,
        }
    };
    let (f, spec) = e.args().split_first()?;
    let (x, a, dir) = match spec {
        [s] => {
            let (x, a) = equation(s)?;
            (x, a, direction(None)?)
        }
        [s, d] if equation(s).is_some() => {
            let (x, a) = equation(s)?;
            (x, a, direction(Some(d))?)
        }
        [x, a] => (x.clone(), a.clone(), direction(None)?),
        [x, a, d] => (x.clone(), a.clone(), direction(Some(d))?),
        _ => return None,
    };
    x.as_symbol()?;
    Some((f, x, Point::of(&a)?, dir))
}

pub(crate) fn limit(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x, point, dir)) = parts(e) else {
        return Ok(None);
    };
    let t = fresh();
    if f.contains(&t) || f.contains(&log_t()) {
        return Ok(None);
    }
    let infinite = matches!(point, Point::Infinity { .. });
    let side = |ev: &mut Eval<'_>, left: bool| -> Result<Result<Value, Failure>, EvalError> {
        let local = ev.eval(&f.substitute(&[(x.clone(), point.variable(&t, left))]))?;
        one_side(ev, &local, &t, 0)
    };
    let value = if infinite || dir == 1 {
        side(ev, false)?
    } else if dir == -1 {
        side(ev, true)?
    } else {
        match (side(ev, false)?, side(ev, true)?) {
            (Ok(right), Ok(left)) => both(ev, right, left)?,
            (Err(failure), _) | (_, Err(failure)) => Err(failure),
        }
    };
    match value {
        Ok(v) => {
            let result = v.expr();
            ev.step("Limit", e, &result);
            Ok(Some(result))
        }
        Err(Failure::Undecided(c)) => {
            ev.warn(Warning::Native {
                message: "the limit's leading coefficient cannot be told from 0".into(),
                expr: c,
            });
            Ok(None)
        }
        Err(Failure::Unknown) => Ok(None),
    }
}

/// Both sides: the same, or opposite infinities.
fn both(ev: &mut Eval<'_>, right: Value, left: Value) -> Result<Result<Value, Failure>, EvalError> {
    Ok(match (right, left) {
        (Value::Finite(r), Value::Finite(l)) => {
            let difference = ev.eval(&Expr::call(Symbol::TOGETHER, [r.clone() - l]))?;
            if difference == Expr::from(0) {
                Ok(Value::Finite(r))
            } else {
                Err(Failure::Unknown)
            }
        }
        (Value::Infinite(r), Value::Infinite(l)) if r == l => Ok(Value::Infinite(r)),
        (Value::Infinite(_) | Value::Unsigned, Value::Infinite(_) | Value::Unsigned) => {
            Ok(Value::Unsigned)
        }
        _ => Err(Failure::Unknown),
    })
}

/// The limit of `local` as `t → 0⁺`, `depth` logarithms in.
fn one_side(
    ev: &mut Eval<'_>,
    local: &Expr,
    t: &Expr,
    depth: u32,
) -> Result<Result<Value, Failure>, EvalError> {
    let Some(series) = leading(ev, local, t)? else {
        return Ok(Err(Failure::Unknown));
    };
    let c = series.coeffs[0].clone();
    let q = Rational::new(
        Integer::from(series.valuation),
        Integer::from(series.ramification),
    );
    let l = log_t();
    if c.contains(&l) {
        return logarithmic(ev, &c, &q, depth);
    }
    if !decided(ev, &c) {
        return Ok(Err(Failure::Undecided(c)));
    }
    Ok(Ok(match q.cmp(&Rational::ZERO) {
        Ordering::Greater => Value::Finite(Expr::from(0)),
        Ordering::Equal => Value::Finite(c),
        Ordering::Less => {
            let mut expander = Expander::new(ev, t.clone());
            match expander.sign(&c) {
                Some(Ordering::Greater) => Value::Infinite(true),
                Some(Ordering::Less) => Value::Infinite(false),
                _ => Value::Unsigned,
            }
        }
    }))
}

/// The series of `local` from its first nonzero term, to a few orders.
fn leading(ev: &mut Eval<'_>, local: &Expr, t: &Expr) -> Result<Option<Series>, EvalError> {
    for order in ORDERS {
        let Some(s) = expansion(ev, local, t, &Rational::from(Integer::from(order)))? else {
            return Ok(None);
        };
        if !s.coeffs.is_empty() {
            return Ok(Some(s));
        }
    }
    Ok(None)
}

/// Whether `c` is told from 0: a number exactly, a numeric expression by
/// its ball; any other, as one with parameters, is taken as it is.
fn decided(ev: &Eval<'_>, c: &Expr) -> bool {
    if c.is_number() {
        return true;
    }
    let Ok(value) = approximation(ev.context(), c, 128) else {
        return true;
    };
    match value.view() {
        View::Ball(b) => !b.contains_rational(&Rational::ZERO),
        View::ComplexBall(z) => !z.may_be_zero(),
        _ => true,
    }
}

/// `c t^q` with `log t` in `c`, `log t → −∞` (§4.2).
fn logarithmic(
    ev: &mut Eval<'_>,
    c: &Expr,
    q: &Rational,
    depth: u32,
) -> Result<Result<Value, Failure>, EvalError> {
    let l = log_t();
    if q.is_zero() {
        if depth >= MAX_LOGARITHMS {
            return Ok(Err(Failure::Unknown));
        }
        // The limit of c as log t → −∞: log t = −1/s, s → 0⁺.
        let s = fresh_logarithm();
        let local = ev.eval(&c.substitute(&[(l, -(s.pow(-1)))]))?;
        return one_side(ev, &local, &s, depth + 1);
    }
    // c a polynomial in log t: its leading term in log t, of sign
    // sign(c_m) (−1)^m.
    let degree = ev.eval(&Expr::call(Symbol::EXPONENT, [c.clone(), l.clone()]))?;
    let Some(m) = degree.as_i64().filter(|m| *m >= 0) else {
        return Ok(Err(Failure::Unknown));
    };
    let lead = ev.eval(&Expr::call(
        Symbol::COEFFICIENT,
        [c.clone(), l.clone(), Expr::from(m)],
    ))?;
    let rest = ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [c.clone() - lead.clone() * l.pow(m)],
    ))?;
    let lower = rest == Expr::from(0) || {
        let k = ev.eval(&Expr::call(Symbol::EXPONENT, [rest, l.clone()]))?;
        k.as_i64().is_some_and(|k| k < m)
    };
    if lead.contains(&l) || !lower {
        return Ok(Err(Failure::Unknown));
    }
    if q.is_positive() {
        return Ok(Ok(Value::Finite(Expr::from(0))));
    }
    let mut expander = Expander::new(ev, fresh());
    Ok(Ok(match expander.sign(&lead) {
        Some(s) => {
            let positive = (s == Ordering::Greater) == (m % 2 == 0);
            Value::Infinite(positive)
        }
        None => Value::Unsigned,
    }))
}

/// The variable of a limit in `log t`.
fn fresh_logarithm() -> Expr {
    Expr::from(Symbol::new("Limit$s"))
}

#[cfg(test)]
mod tests;
