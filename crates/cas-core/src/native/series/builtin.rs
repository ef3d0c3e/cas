//! `Series` (CAS-DOC-10 §1, §3, CAS-204): the expansion of `f` at `x = a`
//! to the order `n`, as giac's `series` and `taylor` write it: the terms in
//! `(x − a)^k` for `k ≤ n`, in `(1/x)^k` at infinity, and the remainder
//! `(x − a)^m OrderSize[x − a]` at the next exponent `m` of the series
//! (`n + 1` but for Puiseux exponents), `OrderSize` an inert head. A Laurent
//! polynomial in `x − a` whose terms all come by the order is itself,
//! without a remainder (giac's `series(1/x+x,x=0,3)`). `log` of a vanishing
//! argument keeps `Log[x − a]` in a coefficient.
//!
//! The forms: `Series[f, {x, a, n}]` (Mathematica's), `Series[f, x == a, n]`
//! from giac's `series(f, x=a, n)`, `Series[f, x, a, n]`; `n` 5 when left
//! out, and `a` 0 for `Series[f, x]`. A variable that is not a symbol, an
//! order that is no nonnegative integer, or `f` without such a series
//! (`exp(1/x)` at 0) leaves the call as it is.

use cas_num::{Integer, Rational};

use super::{Series, expansion, log_t};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's default order.
const DEFAULT_ORDER: i64 = 5;

/// The orders a call may ask, at most: the coefficients cost their square.
const MAX_ORDER: i64 = 64;

/// Where the expansion is: at a point, or at `+∞` or `−∞`.
pub(crate) enum Point {
    At(Expr),
    Infinity { negative: bool },
}

impl Point {
    /// `a`, read: an infinity is one of the two ends.
    pub(crate) fn of(a: &Expr) -> Option<Point> {
        let negative_infinity = Expr::call(
            Symbol::TIMES,
            [Expr::from(-1), Expr::from(Symbol::INFINITY)],
        );
        if a.as_symbol() == Some(Symbol::INFINITY) || is_directed(a, 1) {
            return Some(Point::Infinity { negative: false });
        }
        if *a == negative_infinity || is_directed(a, -1) {
            return Some(Point::Infinity { negative: true });
        }
        if a.contains(&Expr::from(Symbol::INFINITY))
            || a.contains(&Expr::from(Symbol::DIRECTED_INFINITY))
            || a.contains(&Expr::from(Symbol::COMPLEX_INFINITY))
        {
            return None;
        }
        Some(Point::At(a.clone()))
    }

    /// `x` in terms of `t → 0`: `a + t`, `a − t` from the left, `±1/t`.
    pub(crate) fn variable(&self, t: &Expr, left: bool) -> Expr {
        match self {
            Point::At(a) if left => a.clone() - t.clone(),
            Point::At(a) => a.clone() + t.clone(),
            Point::Infinity { negative } => {
                let inverse = t.pow(-1);
                if *negative { -inverse } else { inverse }
            }
        }
    }

    /// `t` in terms of `x`: `x − a`, `1/x` or `−1/x`.
    pub(crate) fn local(&self, x: &Expr) -> Expr {
        match self {
            Point::At(a) => x.clone() - a.clone(),
            Point::Infinity { negative } => {
                let inverse = x.pow(-1);
                if *negative { -inverse } else { inverse }
            }
        }
    }
}

/// Whether `a` is `DirectedInfinity[s]`.
fn is_directed(a: &Expr, s: i64) -> bool {
    a.head() == Some(&Expr::from(Symbol::DIRECTED_INFINITY))
        && a.args().len() == 1
        && a.args()[0].as_i64() == Some(s)
}

/// The variable of the expansions: `t`.
pub(crate) fn fresh() -> Expr {
    Expr::from(Symbol::new("Series$t"))
}

/// `f`, `x`, the point and the order of a call.
fn parts(e: &Expr) -> Option<(&Expr, Expr, Point, i64)> {
    let point = |x: &Expr, a: &Expr, n: Option<&Expr>| -> Option<(Expr, Point, i64)> {
        x.as_symbol()?;
        let n = match n {
            Some(n) => n.as_i64().filter(|n| (0..=MAX_ORDER).contains(n))?,
            None => DEFAULT_ORDER,
        };
        Some((x.clone(), Point::of(a)?, n))
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
    let (x, point, n) = match spec {
        [l] if l.head() == Some(&Expr::from(Symbol::LIST)) => match l.args() {
            [x, a, n] => point(x, a, Some(n))?,
            [x, a] => point(x, a, None)?,
            _ => return None,
        },
        [s] if equation(s).is_some() => {
            let (x, a) = equation(s)?;
            point(&x, &a, None)?
        }
        [s, n] if equation(s).is_some() => {
            let (x, a) = equation(s)?;
            point(&x, &a, Some(n))?
        }
        [x] => point(x, &Expr::from(0), None)?,
        [x, a] => point(x, a, None)?,
        [x, a, n] => point(x, a, Some(n))?,
        _ => return None,
    };
    Some((f, x, point, n))
}

/// The terms of `s`, a series in `t` written in `u` (`x − a` or `1/x`), up
/// to the exponent `n`, and the exponent of the remainder.
fn terms(s: &Series, u: &Expr, n: i64) -> (Vec<Expr>, Rational) {
    let d = s.ramification;
    let last = n * d;
    let mut out = Vec::new();
    for (k, c) in s.coeffs.iter().enumerate() {
        let units = s.valuation + i64::try_from(k).expect("small");
        if units > last {
            break;
        }
        if *c == Expr::from(0) {
            continue;
        }
        let exponent = Rational::new(Integer::from(units), Integer::from(d));
        out.push(c.clone() * u.pow(Expr::from(exponent)));
    }
    let next = Rational::new(
        Integer::from(last.max(s.valuation - 1) + 1),
        Integer::from(d),
    );
    (out, next)
}

/// `Σ c t^e` for the terms of `s` up to the exponent `n`, `log t` as
/// `Log[t]`: what `f` is when its series ends by `n`.
fn truncation(s: &Series, t: &Expr, n: i64) -> Expr {
    let (terms, _) = terms(s, t, n);
    let sum: Expr = terms.into_iter().sum();
    sum.substitute(&[(log_t(), Expr::call(Symbol::LOG, [t.clone()]))])
}

pub(crate) fn series(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((f, x, point, n)) = parts(e) else {
        return Ok(None);
    };
    let t = fresh();
    if f.contains(&t) {
        return Ok(None);
    }
    let local = ev.eval(&f.substitute(&[(x.clone(), point.variable(&t, false))]))?;
    let order = Rational::from(Integer::from(n + 1));
    let Some(s) = expansion(ev, &local, &t, &order)? else {
        return Ok(None);
    };
    let u = point.local(&x);
    let (mut out, next) = terms(&s, &u, n);
    // A series that ends by `n` (giac's `series(1/x+x,x=0,3)`) is itself.
    let rest = local - truncation(&s, &t, n);
    let ends = ev.eval(&Expr::call(Symbol::TOGETHER, [rest]))? == Expr::from(0);
    if !ends {
        out.push(u.pow(Expr::from(next)) * Expr::call(Symbol::ORDER_SIZE, [u.clone()]));
    }
    let log = Expr::call(Symbol::LOG, [u]);
    let result: Expr = out.into_iter().sum();
    let result = result.substitute(&[(log_t(), log)]);
    ev.step("Series", e, &result);
    Ok(Some(result))
}
