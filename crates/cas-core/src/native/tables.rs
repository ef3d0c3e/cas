//! Variation and sign tables (CAS-582): Xcas's `tabvar` and `tabsign`,
//! ports of giac's (misc.cc `step_func_` :8628, `_tabvar` :9152,
//! `_tabsign` :9322).
//!
//! `TabVar[f]`, `[f, x]`, `[f, x, a, b]`, `[f, x = a..b]` (`x` the symbol
//! `x` when not given, the whole line by default), with `diff = 1` (no
//! second derivative) or `diff = 2` (the default) last; `TabSign` with the
//! same ranges. The answer is giac's table, a list of rows: `x` and the
//! points with `" "` between them; for `tabvar` the rows `y' = f′` (`f′`
//! factored), `y = f` and `(y')'` (without `diff = 1`), for `tabsign` the
//! row `f`. A point's cells are values or one-sided limits; an interval's
//! are `"+"`, `"-"`, arrows `"↑"`, `"↓"`, convexity `"+ (∪)"`, `"- (∩)"`,
//! or `"X"` where the function is undefined. giac's steps:
//!
//! - a periodic `f` (sines, cosines, tangents of `a x + b`) is studied on
//!   one period: from `a`, up to `b`, or `[−T/2, T/2]` on the whole line;
//!   then with `f''`'s row, an even or odd `f` on `[−c, c]` from 0;
//! - the points: the singular values of `f` (`domain` in mode 1) in the
//!   closure of its domain, the zeros of `f′` and of `f''`'s numerator
//!   (`Solve`; a periodic one repeated across the range) in the domain,
//!   the ends of the range and of the domain that belong to the domain, an
//!   infinite end of the range where the domain reaches it, the zeros of
//!   the arguments of `abs` and `sign`; sorted, each once;
//! - between two points the signs at the midpoint (one beyond a finite
//!   point next to an infinite one, 0 between two infinite ones), `"X"`
//!   when the midpoint is outside the domain; a zero of `f′` there makes
//!   the table undefined;
//! - an inner singular point gives two columns, `f`'s limits from either
//!   side and `"||"` in the derivatives' rows; another point one column,
//!   `f′` as the pair of its one-sided limits when they differ, `f''` by
//!   its two-sided limit (`infinity` for infinite sides that differ);
//! - a column equal to the one before goes (giac's loop: of three equal
//!   columns one repeat stays).
//!
//! `undef` for infinitely many critical points: a periodic `f′` or `f''`
//! on an unbounded range, a root that is not a number, `f''` zero with
//! `abs` or `sign` in `f`. giac's `fsolve` pass (numeric roots `Solve`
//! misses, rounded to 3 digits) is not ported: such a `Solve` leaves the
//! call, as do parameters other than `x`, a domain `Solve` leaves,
//! undecided comparisons, more than [`MAX_POINTS`] points, and giac's
//! parametric `tabvar([f, g], t)`. The labels giac prints `y'=(…)` and
//! `(y')'` are here `"y'" = …` and the string `"(y')'"`.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: f, x, a, b, p, q, t, v for points and values"
)]

use std::cmp::Ordering;

use super::function_study::domain;
use super::optimization::{period, repeated};
use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most points a table has: each takes up to two columns, each column
/// a limit per row (C0019).
const MAX_POINTS: usize = 32;

/// giac's `do_inflex_tabsign`: 0, 1, 2.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mode {
    /// `diff = 1`: no second derivative's row.
    First,
    /// `tabvar`'s default: with the row of `f''`.
    Second,
    /// `tabsign`: the row of `f`.
    Sign,
}

/// Why a table is not made.
enum Halt {
    /// The call stays.
    Leave,
    /// giac's `undef` (`step_func_` returning 0).
    Undefined,
    Eval(EvalError),
}

impl From<EvalError> for Halt {
    fn from(e: EvalError) -> Self {
        Halt::Eval(e)
    }
}

type Step<T> = Result<T, Halt>;

/// A point of the line or an infinite end.
#[derive(Clone, PartialEq, Debug)]
enum Point {
    Low,
    High,
    At(Expr),
}

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn is_zero(e: &Expr) -> bool {
    e.to_rational().is_some_and(|q| q.is_zero())
}

/// `TabVar[…]`.
pub(super) fn tabvar(ev: &mut Eval<'_>, e: &Expr) -> Native {
    run(ev, e, false)
}

/// `TabSign[…]`.
pub(super) fn tabsign(ev: &mut Eval<'_>, e: &Expr) -> Native {
    run(ev, e, true)
}

fn run(ev: &mut Eval<'_>, e: &Expr, sign: bool) -> Native {
    let Some((f, x, low, high, mode)) = arguments(e, sign) else {
        return Ok(None);
    };
    let minus_infinity = ev.eval(&-Expr::from(Symbol::INFINITY))?;
    let mut study = Study {
        ev,
        x,
        minus_infinity,
        df: Expr::from(Symbol::TRUE),
    };
    match study.table(f, low, high, mode) {
        Ok(table) => Ok(Some(table)),
        Err(Halt::Leave) => Ok(None),
        Err(Halt::Undefined) => Ok(Some(Expr::from(Symbol::INDETERMINATE))),
        Err(Halt::Eval(e)) => Err(e),
    }
}

/// `diff = k`, `k` 1 or 2.
fn option(o: &Expr) -> Option<i64> {
    match (head(o), o.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [d, k]) if d.as_symbol() == Some(Symbol::D) => {
            k.as_i64().filter(|k| (1..=2).contains(k))
        }
        _ => None,
    }
}

type Arguments<'e> = (&'e Expr, Expr, Option<Expr>, Option<Expr>, Mode);

/// `f`, `x`, the ends given and the mode.
fn arguments(e: &Expr, sign: bool) -> Option<Arguments<'_>> {
    let (f, mut rest) = e.args().split_first()?;
    if head(f) == Some(Symbol::LIST) {
        return None;
    }
    let mut mode = if sign { Mode::Sign } else { Mode::Second };
    while let Some((last, before)) = rest.split_last()
        && let Some(k) = option(last)
    {
        if !sign {
            mode = if k == 1 { Mode::First } else { Mode::Second };
        }
        rest = before;
    }
    let variable = |v: &Expr| v.as_symbol().is_some_and(|s| !s.is_builtin());
    let (x, low, high) = match rest {
        [] => (Expr::symbol("x"), None, None),
        [v] if variable(v) => (v.clone(), None, None),
        [r] => {
            let (Some(Symbol::EQUATION | Symbol::EQUAL), [v, range]) = (head(r), r.args()) else {
                return None;
            };
            let (Some(Symbol::RANGE), [a, b]) = (head(range), range.args()) else {
                return None;
            };
            if !variable(v) {
                return None;
            }
            (v.clone(), Some(a.clone()), Some(b.clone()))
        }
        [v, a, b] if variable(v) => (v.clone(), Some(a.clone()), Some(b.clone())),
        _ => return None,
    };
    Some((f, x, low, high, mode))
}

/// The closure of a domain (giac's `strict2large`): `<` as `<=`, `>` as
/// `>=`, `≠` true.
fn closure(e: &Expr) -> Expr {
    match head(e) {
        Some(Symbol::GREATER) => Expr::call(Symbol::GREATER_EQUAL, e.args().to_vec()),
        Some(Symbol::LESS) => Expr::call(Symbol::LESS_EQUAL, e.args().to_vec()),
        Some(Symbol::UNEQUAL) => Expr::from(Symbol::TRUE),
        Some(s @ (Symbol::AND | Symbol::OR | Symbol::LIST)) => {
            Expr::call(s, e.args().iter().map(closure).collect::<Vec<_>>())
        }
        _ => e.clone(),
    }
}

/// The ends a domain names: the sides of its comparisons free of `x`.
fn ends(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    match head(e) {
        Some(
            Symbol::GREATER
            | Symbol::GREATER_EQUAL
            | Symbol::LESS
            | Symbol::LESS_EQUAL
            | Symbol::UNEQUAL
            | Symbol::EQUAL,
        ) => {
            for side in e.args() {
                if !side.contains(x) && !out.contains(side) {
                    out.push(side.clone());
                }
            }
        }
        Some(Symbol::AND | Symbol::OR | Symbol::LIST) => {
            for part in e.args() {
                ends(part, x, out);
            }
        }
        _ => {}
    }
}

/// Whether `f` has a function of `heads`.
fn has(f: &Expr, heads: &[Symbol]) -> bool {
    f.distinct()
        .any(|s| head(s).is_some_and(|h| heads.contains(&h)))
}

/// Functions with jumps: a value there is not a limit.
const JUMPS: [Symbol; 6] = [
    Symbol::SIGN,
    Symbol::FLOOR,
    Symbol::CEILING,
    Symbol::ROUND,
    Symbol::HEAVISIDE_THETA,
    Symbol::PIECEWISE,
];

struct Study<'a, 'e> {
    ev: &'a mut Eval<'e>,
    x: Expr,
    minus_infinity: Expr,
    /// The domain of `f`, as `domain` gives it.
    df: Expr,
}

impl Study<'_, '_> {
    fn eval(&mut self, e: &Expr) -> Step<Expr> {
        Ok(self.ev.eval(e)?)
    }

    fn call<const N: usize>(&mut self, f: Symbol, args: [Expr; N]) -> Step<Expr> {
        self.eval(&Expr::call(f, args))
    }

    fn expr(&self, p: &Point) -> Expr {
        match p {
            Point::Low => self.minus_infinity.clone(),
            Point::High => Expr::from(Symbol::INFINITY),
            Point::At(a) => a.clone(),
        }
    }

    /// An end of the range: `None` the infinite one, otherwise a real
    /// number or an infinity.
    fn end(&mut self, given: Option<Expr>, infinite: Point) -> Step<Point> {
        let Some(given) = given else {
            return Ok(infinite);
        };
        let v = self.eval(&given)?;
        if v == Expr::from(Symbol::INFINITY) {
            Ok(Point::High)
        } else if v == self.minus_infinity {
            Ok(Point::Low)
        } else if self.real(&v)? {
            Ok(Point::At(v))
        } else {
            Err(Halt::Leave)
        }
    }

    /// Whether `v` is a real number: constants only, a finite ball.
    fn real(&mut self, v: &Expr) -> Step<bool> {
        if !v.symbols().iter().all(|s| s.is_builtin())
            || [
                Symbol::INFINITY,
                Symbol::COMPLEX_INFINITY,
                Symbol::DIRECTED_INFINITY,
                Symbol::INDETERMINATE,
            ]
            .iter()
            .any(|s| v.contains(&Expr::from(*s)))
        {
            return Ok(false);
        }
        let n = self.call(Symbol::N, [v.clone()])?;
        Ok(matches!(n.view(), View::Ball(b) if b.is_finite()))
    }

    /// The sign of a real number; the call stays when it is not decided.
    fn sign(&mut self, e: &Expr) -> Step<Ordering> {
        let d = self.eval(e)?;
        if let Some(q) = d.to_rational() {
            return Ok(q.cmp(&cas_num::Rational::ZERO));
        }
        if let Some(s) = real_sign(self.ev.context(), &d) {
            return Ok(s);
        }
        let s = self.call(Symbol::SIMPLIFY, [d])?;
        if is_zero(&s) {
            return Ok(Ordering::Equal);
        }
        real_sign(self.ev.context(), &s).ok_or(Halt::Leave)
    }

    fn cmp(&mut self, a: &Point, b: &Point) -> Step<Ordering> {
        Ok(match (a, b) {
            (Point::Low, Point::Low) | (Point::High, Point::High) => Ordering::Equal,
            (Point::Low, _) | (_, Point::High) => Ordering::Less,
            (_, Point::Low) | (Point::High, _) => Ordering::Greater,
            (Point::At(a), Point::At(b)) => self.sign(&(a.clone() - b.clone()))?,
        })
    }

    fn within(&mut self, p: &Point, low: &Point, high: &Point) -> Step<bool> {
        Ok(self.cmp(p, low)? != Ordering::Less && self.cmp(p, high)? != Ordering::Greater)
    }

    /// Whether `p` satisfies the conditions `set` (`x` the whole line, a
    /// list any of them), giac's `in_domain`.
    fn holds(&mut self, set: &Expr, p: &Point) -> Step<bool> {
        if *set == self.x {
            return Ok(true);
        }
        if head(set) == Some(Symbol::LIST) {
            for part in set.args() {
                if self.holds(part, p)? {
                    return Ok(true);
                }
            }
            return Ok(false);
        }
        // An infinite point is no excluded value.
        let set = match p {
            Point::At(_) => set.clone(),
            _ => closure(set),
        };
        let at = set.substitute(&[(self.x.clone(), self.expr(p))]);
        Ok(self.eval(&at)?.as_symbol() == Some(Symbol::TRUE))
    }

    fn in_domain(&mut self, p: &Point) -> Step<bool> {
        let df = self.df.clone();
        self.holds(&df, p)
    }

    /// Whether `e` is 0, after `Together` and then `Simplify`.
    fn vanishes(&mut self, e: &Expr) -> Step<bool> {
        let v = self.eval(e)?;
        if is_zero(&v) {
            return Ok(true);
        }
        let t = self.call(Symbol::TOGETHER, [v])?;
        if is_zero(&t) {
            return Ok(true);
        }
        Ok(is_zero(&self.call(Symbol::SIMPLIFY, [t])?))
    }

    /// giac's `try_limit_undef`: the limit of `g` at `p` from the right
    /// (`1`), the left (`-1`) or both sides (`0`); a value where `g` has
    /// one and no jump function; `infinity` for an unsigned infinity,
    /// `undef` for none or a complex one.
    fn limit(&mut self, g: &Expr, p: &Point, dir: i64) -> Step<Expr> {
        let v = if let Point::At(a) = p {
            if !has(g, &JUMPS) {
                let v = self.eval(&g.substitute(&[(self.x.clone(), a.clone())]))?;
                if self.real(&v)? {
                    return self.normal(v);
                }
            }
            let mut args = vec![g.clone(), self.x.clone(), a.clone()];
            if dir != 0 {
                args.push(Expr::from(dir));
            }
            self.eval(&Expr::call(Symbol::LIMIT, args))?
        } else {
            let at = self.expr(p);
            self.call(Symbol::LIMIT, [g.clone(), self.x.clone(), at])?
        };
        if v == Expr::from(Symbol::INFINITY)
            || v == self.minus_infinity
            || v == Expr::from(Symbol::COMPLEX_INFINITY)
        {
            return Ok(v);
        }
        if self.real(&v)? {
            return self.normal(v);
        }
        if head(&v) != Some(Symbol::LIMIT)
            && (v.contains(&Expr::from(Symbol::COMPLEX_INFINITY))
                || v.contains(&Expr::from(Symbol::DIRECTED_INFINITY)))
        {
            return Ok(Expr::from(Symbol::COMPLEX_INFINITY));
        }
        Ok(Expr::from(Symbol::INDETERMINATE))
    }

    /// A value as giac's `recursive_normal` leaves it: `Together`.
    fn normal(&mut self, v: Expr) -> Step<Expr> {
        if v.to_rational().is_some() {
            return Ok(v);
        }
        self.call(Symbol::TOGETHER, [v])
    }

    /// Whether a value is positive (giac's `is_strictly_positive`); `None`
    /// for `undef`.
    fn positive(&mut self, v: &Expr) -> Step<Option<bool>> {
        if *v == Expr::from(Symbol::INDETERMINATE) {
            return Ok(None);
        }
        if *v == Expr::from(Symbol::INFINITY) {
            return Ok(Some(true));
        }
        if *v == self.minus_infinity || *v == Expr::from(Symbol::COMPLEX_INFINITY) {
            return Ok(Some(false));
        }
        Ok(Some(self.sign(v)? == Ordering::Greater))
    }

    /// Whether two limits are one value.
    fn same(&mut self, a: &Expr, b: &Expr) -> Step<bool> {
        if a == b {
            return Ok(true);
        }
        Ok(self.real(a)?
            && self.real(b)?
            && self.sign(&(a.clone() - b.clone()))? == Ordering::Equal)
    }

    /// The real zeros of `g` in `[low, high]`: `Solve`'s, a periodic `g`'s
    /// repeated across the range; `undef` for infinitely many.
    fn roots(&mut self, g: &Expr, low: &Point, high: &Point) -> Step<Vec<Point>> {
        if !g.contains(&self.x) {
            return if is_zero(g) {
                Err(Halt::Undefined)
            } else {
                Ok(Vec::new())
            };
        }
        let periodic = period(self.ev, g, &self.x)?.is_some();
        let bounded = match (low, high) {
            (Point::At(a), Point::At(b)) => Some((a.clone(), b.clone())),
            _ => None,
        };
        if periodic && bounded.is_none() {
            return Err(Halt::Undefined);
        }
        // giac's solve takes the numerator.
        let t = self.call(Symbol::TOGETHER, [g.clone()])?;
        let n = self.call(Symbol::NUMERATOR, [t])?;
        if !n.contains(&self.x) {
            return if is_zero(&n) {
                Err(Halt::Undefined)
            } else {
                Ok(Vec::new())
            };
        }
        let equation = Expr::call(Symbol::EQUATION, [n, Expr::from(0)]);
        let solved = self.call(Symbol::SOLVE, [equation, self.x.clone()])?;
        if head(&solved) != Some(Symbol::LIST) {
            return Err(Halt::Leave);
        }
        let mut found = Vec::new();
        for r in solved.args() {
            if !r.symbols().iter().all(|s| s.is_builtin()) {
                return Err(Halt::Undefined);
            }
            if self.real(r)? {
                found.push(r.clone());
            }
        }
        if periodic && let Some((a, b)) = &bounded {
            found = repeated(self.ev, g, &self.x, found, a, b)?.ok_or(Halt::Leave)?;
        }
        let mut out = Vec::new();
        for r in found {
            let p = Point::At(r);
            if self.within(&p, low, high)? {
                out.push(p);
            }
        }
        Ok(out)
    }

    /// The test point between `cur` and `next`; `None` when the midpoint
    /// of two finite points is outside the domain.
    fn test_point(&mut self, cur: &Point, next: &Point) -> Step<Option<Point>> {
        let t = match (cur, next) {
            (Point::At(a), Point::At(b)) => {
                let m = if a.to_rational().is_some() && b.to_rational().is_some() {
                    self.eval(&((a.clone() + b.clone()) / Expr::from(2)))?
                } else {
                    let (p, q) = (self.approximation(a)?, self.approximation(b)?);
                    self.eval(&((Expr::from(p) + Expr::from(q)) / Expr::from(2)))?
                };
                let m = Point::At(m);
                return Ok(self.in_domain(&m)?.then_some(m));
            }
            (Point::Low, Point::At(b)) => self.eval(&(b.clone() - Expr::from(1)))?,
            (Point::At(a), Point::High) => self.eval(&(a.clone() + Expr::from(1)))?,
            _ => Expr::from(0),
        };
        Ok(Some(Point::At(t)))
    }

    /// The midpoint of `N`'s ball for a real number, a rational.
    fn approximation(&mut self, e: &Expr) -> Step<cas_num::Rational> {
        let n = self.call(Symbol::N, [e.clone()])?;
        match n.view() {
            View::Ball(b) if b.is_finite() => Ok(b.mid().to_rational()),
            _ => Err(Halt::Leave),
        }
    }

    /// The points, sorted, each once.
    fn sorted(&mut self, points: Vec<Point>) -> Step<Vec<Point>> {
        let mut out: Vec<Point> = Vec::new();
        for p in points {
            let mut at = out.len();
            let mut seen = false;
            for (i, q) in out.iter().enumerate() {
                match self.cmp(&p, &q.clone())? {
                    Ordering::Equal => {
                        seen = true;
                        break;
                    }
                    Ordering::Less => {
                        at = i;
                        break;
                    }
                    Ordering::Greater => {}
                }
            }
            if !seen {
                out.insert(at, p);
            }
            if out.len() > MAX_POINTS {
                return Err(Halt::Leave);
            }
        }
        Ok(out)
    }

    /// giac's `step_func_`, the table part.
    #[allow(
        clippy::too_many_lines,
        reason = "giac's steps in order, as the module documents them"
    )]
    fn table(&mut self, f: &Expr, low: Option<Expr>, high: Option<Expr>, mode: Mode) -> Step<Expr> {
        let x = self.x.clone();
        if !f.contains(&x)
            || f.symbols()
                .iter()
                .any(|s| !s.is_builtin() && Expr::from(*s) != x)
        {
            return Err(Halt::Leave);
        }
        let mut xmin = self.end(low, Point::Low)?;
        let mut xmax = self.end(high, Point::High)?;
        if self.cmp(&xmin, &xmax)? != Ordering::Less {
            return Err(Halt::Leave);
        }
        // One period.
        if let Some(tau) = period(self.ev, f, &x)? {
            let t = self.eval(&(Expr::from(tau) * Expr::from(Symbol::PI)))?;
            let wide = match (&xmin, &xmax) {
                (Point::At(a), Point::At(b)) => {
                    self.sign(&(b.clone() - a.clone() - t.clone()))? == Ordering::Greater
                }
                _ => true,
            };
            if wide {
                match (&xmin, &xmax) {
                    (Point::At(a), _) => xmax = Point::At(self.eval(&(a.clone() + t))?),
                    (_, Point::At(b)) => xmin = Point::At(self.eval(&(b.clone() - t))?),
                    _ => {
                        let half = self.eval(&(t / Expr::from(2)))?;
                        xmin = Point::At(self.eval(&-half.clone())?);
                        xmax = Point::At(half);
                    }
                }
            }
        }
        // Even or odd, with f''.
        if mode == Mode::Second {
            let symmetric = match (&xmin, &xmax) {
                (Point::Low, Point::High) => true,
                (Point::At(a), Point::At(b)) => {
                    self.sign(&(a.clone() + b.clone()))? == Ordering::Equal
                }
                _ => false,
            };
            if symmetric {
                let g = self.eval(&f.substitute(&[(x.clone(), -x.clone())]))?;
                if self.vanishes(&(g.clone() - f.clone()))? || self.vanishes(&(g + f.clone()))? {
                    xmin = Point::At(Expr::from(0));
                }
            }
        }
        // The domain and the singular values.
        let Some(df) = domain(self.ev, f, &x, false)? else {
            return Err(Halt::Leave);
        };
        self.df = df.clone();
        let large = closure(&df);
        let Some(singular) = domain(self.ev, f, &x, true)? else {
            return Err(Halt::Leave);
        };
        let mut sing = Vec::new();
        for v in singular.args() {
            if self.real(v)? {
                let p = Point::At(v.clone());
                if self.holds(&large, &p)? && self.within(&p, &xmin, &xmax)? {
                    sing.push(p);
                }
            }
        }
        // The derivatives and the critical points.
        let f1 = if mode == Mode::Sign {
            f.clone()
        } else {
            let d = Expr::call(Symbol::D, [f.clone(), x.clone()]);
            self.call(Symbol::FACTOR, [d])?
        };
        let f2 = if mode == Mode::Second {
            Some(self.call(Symbol::D, [f1.clone(), x.clone()])?)
        } else {
            None
        };
        // f'' of |u| is 0 but not structurally: giac's solve gives x.
        let mut inflexions = Vec::new();
        if let Some(f2) = &f2 {
            let t = self.call(Symbol::TOGETHER, [f2.clone()])?;
            if is_zero(&t) {
                if has(f, &[Symbol::ABS, Symbol::SIGN]) {
                    return Err(Halt::Undefined);
                }
            } else {
                inflexions = self.roots(&t, &xmin, &xmax)?;
            }
        }
        let mut critical = self.roots(&f1, &xmin, &xmax)?;
        critical.extend(inflexions);
        let mut points = Vec::new();
        if self.in_domain(&xmin)? {
            points.push(xmin.clone());
        }
        for end in [Point::Low, Point::High] {
            if (xmin == end || xmax == end) && self.in_domain(&end)? {
                sing.push(end);
            }
        }
        points.extend(sing.iter().cloned());
        for p in critical {
            if self.in_domain(&p)? {
                points.push(p);
            }
        }
        if self.in_domain(&xmax)? {
            points.push(xmax.clone());
        }
        let mut named = Vec::new();
        ends(&df, &x, &mut named);
        for e in named {
            if self.real(&e)? {
                let p = Point::At(e);
                if self.within(&p, &xmin, &xmax)? && self.in_domain(&p)? {
                    points.push(p);
                }
            }
        }
        let kinks: Vec<Expr> = f
            .distinct()
            .filter(|s| matches!(head(s), Some(Symbol::ABS | Symbol::SIGN)) && s.contains(&x))
            .filter_map(|s| s.args().first().cloned())
            .collect();
        for u in kinks {
            points.extend(self.roots(&u, &xmin, &xmax)?);
        }
        let points = self.sorted(points)?;
        let Some(first) = points.first() else {
            return Err(Halt::Leave);
        };
        // The first column.
        let mut rx = vec![x.clone(), self.expr(first)];
        let mut rd = vec![
            if mode == Mode::Sign {
                f.clone()
            } else {
                Expr::call(Symbol::EQUATION, [Expr::string("y'"), f1.clone()])
            },
            self.limit(&f1, first, 1)?,
        ];
        let mut ry = vec![
            Expr::call(Symbol::EQUATION, [Expr::symbol("y"), f.clone()]),
            self.limit(f, first, 1)?,
        ];
        let mut rdd = Vec::new();
        if let Some(f2) = &f2 {
            rdd = vec![Expr::string("(y')'"), self.limit(f2, first, 1)?];
        }
        let s = Expr::string;
        // The intervals and the points after them.
        for i in 1..points.len() {
            let (cur, next) = (&points[i - 1], &points[i]);
            rx.push(s(" "));
            let (dfx, d2) = match self.test_point(cur, next)? {
                Some(t) => (
                    self.limit(&f1, &t, 0)?,
                    match &f2 {
                        Some(f2) => self.limit(f2, &t, 0)?,
                        None => Expr::from(Symbol::INDETERMINATE),
                    },
                ),
                None => (
                    Expr::from(Symbol::INDETERMINATE),
                    Expr::from(Symbol::INDETERMINATE),
                ),
            };
            if is_zero(&dfx) {
                return Err(Halt::Undefined);
            }
            let (arrow, plus) = match self.positive(&dfx)? {
                None => ("X", "X"),
                Some(true) => ("↑", "+"),
                Some(false) => ("↓", "-"),
            };
            ry.push(s(arrow));
            rd.push(s(plus));
            if f2.is_some() {
                rdd.push(s(match self.positive(&d2)? {
                    None => "X",
                    Some(true) => "+ (∪)",
                    Some(false) => "- (∩)",
                }));
            }
            let mut inner = false;
            if i < points.len() - 1 {
                for q in &sing {
                    if self.cmp(next, q)? == Ordering::Equal {
                        inner = true;
                        break;
                    }
                }
            }
            if inner {
                for dir in [-1, 1] {
                    rx.push(self.expr(next));
                    ry.push(self.limit(f, next, dir)?);
                    rd.push(s("||"));
                    if f2.is_some() {
                        rdd.push(s("||"));
                    }
                }
            } else {
                rx.push(self.expr(next));
                ry.push(self.limit(f, next, -1)?);
                let l = self.limit(&f1, next, -1)?;
                let mut cell = l.clone();
                if matches!(next, Point::At(_)) {
                    let r = self.limit(&f1, next, 1)?;
                    if !self.same(&l, &r)? {
                        cell = list(vec![l, r]);
                    }
                }
                rd.push(cell);
                if let Some(f2) = &f2 {
                    rdd.push(self.limit(f2, next, 0)?);
                }
            }
        }
        let mut rows = match mode {
            Mode::Sign => vec![rx, rd],
            _ => vec![rx, rd, ry],
        };
        if f2.is_some() {
            rows.push(rdd);
        }
        // A column equal to the one before goes; the next is not compared
        // with it (giac's loop).
        let mut columns: Vec<Vec<Expr>> = (0..rows[0].len())
            .map(|j| rows.iter().map(|r| r[j].clone()).collect())
            .collect();
        let mut j = 1;
        while j < columns.len() {
            if columns[j] == columns[j - 1] {
                columns.remove(j);
            }
            j += 1;
        }
        Ok(list(
            (0..rows.len())
                .map(|k| list(columns.iter().map(|c| c[k].clone()).collect()))
                .collect(),
        ))
    }
}

#[cfg(test)]
mod tests;
