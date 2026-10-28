//! Inequalities (CAS-DOC-19, CAS-277): `Solve[ineq, x]` for a real `x` by
//! sign charts, as giac's list of conditions.
//!
//! `l op r` is `f op 0` with `f = Together[l − r] = N/D`. The real roots of
//! `N` (zeros) and of `D` (poles), by [`roots`], sorted by the sign of their
//! differences as `Sign` decides it (assumptions, balls, approximation),
//! cut the line into pieces on which `f` keeps its sign: an exact test
//! point decides each (the midpoint, or the first point less 1 and the last
//! plus 1), a zero holds where `op` admits equality, a pole never. Adjacent
//! pieces join. So `x > a` with `a` a parameter is `x > a`: one point, and
//! `f` at `a ± 1` is `±1`.
//!
//! The answer lists the intervals in increasing order: `x < c`, `x <= c`,
//! `x > c`, `x >= c`, `a < x < b` (an `And` of two comparisons, printed as a
//! chain, giac's `((x>a) and (x<b))`), a point `c` alone, `x` for the whole
//! line; no interval, no element. What is not decided stays.
//!
//! Systems (CAS-278): a list or an `And` of conditions (a chain `a < x <= b`
//! is one) intersects their sets, an `Or` unites them, a `Not` takes the
//! complement; an equation among them is charted as `f = 0`. Each set is
//! disjoint intervals in increasing order, ends compared as the sign
//! chart's points are.
//!
//! Functions (CAS-279): absolute values, powers, exponentials and
//! logarithms of `x`, nothing periodic. Radicands of even roots (`u >= 0`,
//! `u > 0` under a negative power; odd roots of negatives are real, O7,
//! but in the Complex domain) and logarithms' arguments (`u > 0`) make a
//! domain,
//! itself solved so, whose ends are cuts of the chart, where `f` is
//! evaluated, and outside of which nothing holds. The poles of functions'
//! arguments are poles of the chart (`exp(1/x)` at 0). Zeros: of a
//! product, its factors'; with absolute values, piecewise between the
//! zeros and poles of their arguments, which are cuts; else `Solve`'s,
//! radicals and inversions of functions (CAS-213) among them.

use std::cmp::Ordering;

use super::super::rounding::real_sign;
use super::roots;
use crate::assume::Domain;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// Whether `e` is an inequality `Solve` hands here: a comparison, or a
/// list, `And`, `Or` or `Not` of conditions (equations among them) with
/// at least one comparison (CAS-278).
pub(super) fn is_inequality(e: &Expr) -> bool {
    condition(e) && has_comparison(e)
}

/// Whether `e` is a condition on `x`: a comparison, an equation, `True`,
/// `False`, or a connective of conditions.
fn condition(e: &Expr) -> bool {
    if relation(e).is_some() || equation(e).is_some() {
        return true;
    }
    match e.as_symbol() {
        Some(Symbol::TRUE | Symbol::FALSE) => return true,
        Some(_) => return false,
        None => {}
    }
    matches!(
        e.head().and_then(Expr::as_symbol),
        Some(Symbol::AND | Symbol::OR | Symbol::NOT | Symbol::LIST)
    ) && e.args().iter().all(condition)
}

/// Whether a comparison occurs in the condition `e`.
fn has_comparison(e: &Expr) -> bool {
    relation(e).is_some() || (equation(e).is_none() && e.args().iter().any(has_comparison))
}

/// The sides of an equation, `Equation[l, r]` or `Equal[l, r]`.
fn equation(e: &Expr) -> Option<(&Expr, &Expr)> {
    match (e.head()?.as_symbol()?, e.args()) {
        (Symbol::EQUATION | Symbol::EQUAL, [l, r]) => Some((l, r)),
        _ => None,
    }
}

/// The relation and the sides of a comparison other than an equation.
fn relation(e: &Expr) -> Option<(Symbol, &Expr, &Expr)> {
    match (e.head()?.as_symbol()?, e.args()) {
        (
            s @ (Symbol::LESS
            | Symbol::LESS_EQUAL
            | Symbol::GREATER
            | Symbol::GREATER_EQUAL
            | Symbol::UNEQUAL),
            [l, r],
        ) => Some((s, l, r)),
        _ => None,
    }
}

/// Whether a sign satisfies the relation `op` against 0, an equation's
/// `Equal` among them.
fn holds(op: Symbol, sign: Ordering) -> bool {
    match op {
        Symbol::LESS => sign == Ordering::Less,
        Symbol::LESS_EQUAL => sign != Ordering::Greater,
        Symbol::GREATER => sign == Ordering::Greater,
        Symbol::GREATER_EQUAL => sign != Ordering::Less,
        Symbol::UNEQUAL => sign != Ordering::Equal,
        _ => sign == Ordering::Equal,
    }
}

/// An end of an interval: a real, held or not.
#[derive(Clone, Debug)]
struct End {
    at: Expr,
    closed: bool,
}

/// A real interval, `None` ends infinite.
#[derive(Clone, Debug)]
struct Interval {
    lo: Option<End>,
    hi: Option<End>,
}

/// The solutions of `ineq` in `x` as giac's list of conditions; `None`
/// when undecided.
pub(super) fn solve(
    ev: &mut Eval<'_>,
    e: &Expr,
    ineq: &Expr,
    x: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let Some(set) = set(ev, e, ineq, x)? else {
        return Ok(None);
    };
    let conditions: Vec<Expr> = set.iter().map(|i| render(i, x)).collect();
    ev.eval(&Expr::call(Symbol::LIST, conditions)).map(Some)
}

/// The solution set of the condition `c`, intervals in increasing order:
/// a comparison or an equation by its sign chart, a list or `And` the
/// intersection of its members', `Or` the union, `Not` the complement
/// (CAS-278).
fn set(
    ev: &mut Eval<'_>,
    e: &Expr,
    c: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Interval>>, EvalError> {
    if let Some((op, l, r)) = relation(c) {
        return one(ev, e, op, l, r, x);
    }
    if let Some((l, r)) = equation(c) {
        return one(ev, e, Symbol::EQUAL, l, r, x);
    }
    match c.as_symbol() {
        Some(Symbol::TRUE) => return Ok(Some(vec![Interval { lo: None, hi: None }])),
        Some(Symbol::FALSE) => return Ok(Some(Vec::new())),
        _ => {}
    }
    let mut members = Vec::with_capacity(c.args().len());
    for m in c.args() {
        let Some(s) = set(ev, e, m, x)? else {
            return Ok(None);
        };
        members.push(s);
    }
    match c.head().and_then(Expr::as_symbol) {
        Some(Symbol::AND | Symbol::LIST) => {
            let mut all = vec![Interval { lo: None, hi: None }];
            for m in &members {
                let Some(next) = intersection(ev, &all, m)? else {
                    return Ok(None);
                };
                all = next;
            }
            Ok(Some(all))
        }
        Some(Symbol::OR) => union(ev, members.into_iter().flatten().collect()),
        Some(Symbol::NOT) => match members.as_slice() {
            [m] => Ok(Some(complement(m))),
            _ => Ok(None),
        },
        _ => Ok(None),
    }
}

/// Lower ends ordered: −∞ first; at one value, a closed end first.
fn lower(
    ev: &mut Eval<'_>,
    a: Option<&End>,
    b: Option<&End>,
) -> Result<Option<Ordering>, EvalError> {
    Ok(match (a, b) {
        (None, None) => Some(Ordering::Equal),
        (None, Some(_)) => Some(Ordering::Less),
        (Some(_), None) => Some(Ordering::Greater),
        (Some(a), Some(b)) => match compare(ev, &a.at, &b.at)? {
            Some(Ordering::Equal) => Some(b.closed.cmp(&a.closed)),
            other => other,
        },
    })
}

/// Upper ends ordered: +∞ last; at one value, an open end first.
fn upper(
    ev: &mut Eval<'_>,
    a: Option<&End>,
    b: Option<&End>,
) -> Result<Option<Ordering>, EvalError> {
    Ok(match (a, b) {
        (None, None) => Some(Ordering::Equal),
        (None, Some(_)) => Some(Ordering::Greater),
        (Some(_), None) => Some(Ordering::Less),
        (Some(a), Some(b)) => match compare(ev, &a.at, &b.at)? {
            Some(Ordering::Equal) => Some(a.closed.cmp(&b.closed)),
            other => other,
        },
    })
}

/// Whether the interval from `lo` to `hi` holds a point.
fn nonempty(
    ev: &mut Eval<'_>,
    lo: Option<&End>,
    hi: Option<&End>,
) -> Result<Option<bool>, EvalError> {
    let (Some(lo), Some(hi)) = (lo, hi) else {
        return Ok(Some(true));
    };
    Ok(compare(ev, &lo.at, &hi.at)?.map(|o| match o {
        Ordering::Less => true,
        Ordering::Equal => lo.closed && hi.closed,
        Ordering::Greater => false,
    }))
}

/// `a ∩ b` of two sets of disjoint intervals in increasing order.
fn intersection(
    ev: &mut Eval<'_>,
    a: &[Interval],
    b: &[Interval],
) -> Result<Option<Vec<Interval>>, EvalError> {
    let mut out = Vec::new();
    for p in a {
        for q in b {
            // The greater lower end and the lesser upper end.
            let Some(lo_order) = lower(ev, p.lo.as_ref(), q.lo.as_ref())? else {
                return Ok(None);
            };
            let lo = if lo_order == Ordering::Less {
                &q.lo
            } else {
                &p.lo
            };
            let Some(hi_order) = upper(ev, p.hi.as_ref(), q.hi.as_ref())? else {
                return Ok(None);
            };
            let hi = if hi_order == Ordering::Greater {
                &q.hi
            } else {
                &p.hi
            };
            match nonempty(ev, lo.as_ref(), hi.as_ref())? {
                Some(true) => out.push(Interval {
                    lo: lo.clone(),
                    hi: hi.clone(),
                }),
                Some(false) => {}
                None => return Ok(None),
            }
        }
    }
    union(ev, out)
}

/// The union of intervals, sorted by their lower ends and joined where
/// they overlap or meet at a point one of them holds.
fn union(ev: &mut Eval<'_>, intervals: Vec<Interval>) -> Result<Option<Vec<Interval>>, EvalError> {
    let mut sorted: Vec<Interval> = Vec::with_capacity(intervals.len());
    for i in intervals {
        // By insertion.
        let mut at = sorted.len();
        for (k, s) in sorted.iter().enumerate() {
            match lower(ev, i.lo.as_ref(), s.lo.as_ref())? {
                Some(Ordering::Less) => {
                    at = k;
                    break;
                }
                Some(_) => {}
                None => return Ok(None),
            }
        }
        sorted.insert(at, i);
    }
    let mut out: Vec<Interval> = Vec::with_capacity(sorted.len());
    for i in sorted {
        if let Some(last) = out.last_mut() {
            let meets = match (&last.hi, &i.lo) {
                (None, _) | (_, None) => true,
                (Some(h), Some(l)) => match compare(ev, &l.at, &h.at)? {
                    Some(Ordering::Less) => true,
                    Some(Ordering::Equal) => h.closed || l.closed,
                    Some(Ordering::Greater) => false,
                    None => return Ok(None),
                },
            };
            if meets {
                let Some(order) = upper(ev, last.hi.as_ref(), i.hi.as_ref())? else {
                    return Ok(None);
                };
                if order == Ordering::Less {
                    last.hi = i.hi;
                }
                continue;
            }
        }
        out.push(i);
    }
    Ok(Some(out))
}

/// The complement of a set of disjoint intervals in increasing order: the
/// gaps between them, each end held where the interval's is not.
fn complement(set: &[Interval]) -> Vec<Interval> {
    let flip = |end: &End| End {
        at: end.at.clone(),
        closed: !end.closed,
    };
    let mut out = Vec::new();
    // The gap's lower end, while one is open: `Some(None)` from −∞.
    let mut lo: Option<Option<End>> = Some(None);
    for i in set {
        if let (Some(start), Some(end)) = (lo.take(), &i.lo) {
            out.push(Interval {
                lo: start,
                hi: Some(flip(end)),
            });
        }
        lo = i.hi.as_ref().map(|h| Some(flip(h)));
    }
    if let Some(start) = lo {
        out.push(Interval {
            lo: start,
            hi: None,
        });
    }
    out
}

/// `a`'s order against `b`, by their difference's sign.
fn compare(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<Ordering>, EvalError> {
    if a == b {
        return Ok(Some(Ordering::Equal));
    }
    let gap = ev.eval(&Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]))?;
    Ok(real_sign(ev.context(), &gap))
}

/// `f`'s sign at `x = p`, `f` put together there.
fn sign_at(ev: &mut Eval<'_>, f: &Expr, x: &Expr, p: &Expr) -> Result<Option<Ordering>, EvalError> {
    let value = ev.eval(&Expr::call(
        Symbol::TOGETHER,
        [f.substitute(&[(x.clone(), p.clone())])],
    ))?;
    Ok(real_sign(ev.context(), &value))
}

/// What a point of a sign chart is, by precedence where points meet: a
/// cut, where `f` is evaluated (an end of the domain, a kink of an absolute
/// value); a zero of `f`; a pole, of `f` or of a function's argument.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Kind {
    Cut,
    Zero,
    Pole,
}

/// The real zeros of an expression in `x`, and the cuts where its sign may
/// change besides, each maybe repeated.
#[derive(Default)]
struct Crossings {
    zeros: Vec<Expr>,
    cuts: Vec<Expr>,
}

/// The solutions of `l op r`, `op` a comparison or `Equal`, by a sign
/// chart.
fn one(
    ev: &mut Eval<'_>,
    e: &Expr,
    op: Symbol,
    l: &Expr,
    r: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Interval>>, EvalError> {
    if !charted(l, x) || !charted(r, x) {
        return Ok(None);
    }
    let Some(domain) = domain(ev, e, l, r, x)? else {
        return Ok(None);
    };
    let Some(singular) = singular(ev, e, l, r, x)? else {
        return Ok(None);
    };
    let f = ev.eval(&Expr::call(Symbol::TOGETHER, [l.clone() - r.clone()]))?;
    let n = ev.eval(&Expr::call(Symbol::NUMERATOR, [f.clone()]))?;
    let d = ev.eval(&Expr::call(Symbol::DENOMINATOR, [f.clone()]))?;
    let Some(above) = crossings(ev, e, &n, x)? else {
        return Ok(None);
    };
    let Some(below) = crossings(ev, e, &d, x)? else {
        return Ok(None);
    };
    // A root that is `x` itself: every x, as for f ≡ 0, which `Together`
    // already took.
    if above
        .zeros
        .iter()
        .chain(&below.zeros)
        .any(|p| p.contains(x))
    {
        return Ok(None);
    }
    let ends: Vec<Expr> = domain
        .iter()
        .flat_map(|i| [&i.lo, &i.hi])
        .flatten()
        .map(|end| end.at.clone())
        .collect();
    let all = (above.zeros.into_iter().map(|p| (p, Kind::Zero)))
        .chain(
            below
                .zeros
                .into_iter()
                .chain(singular)
                .map(|p| (p, Kind::Pole)),
        )
        .chain(
            above
                .cuts
                .into_iter()
                .chain(below.cuts)
                .chain(ends)
                .map(|p| (p, Kind::Cut)),
        );
    let Some(points) = sorted(ev, all)? else {
        return Ok(None);
    };
    // Whether each piece holds, outside the domain never, else by `f`'s
    // sign at a point of it: piece i lies below point i.
    let mut pieces = Vec::with_capacity(points.len() + 1);
    for i in 0..=points.len() {
        let test = between(
            i.checked_sub(1).map(|j| &points[j].0),
            points.get(i).map(|(p, _)| p),
        );
        let holding = match member(ev, &domain, &test)? {
            Some(true) => sign_at(ev, &f, x, &test)?.map(|s| holds(op, s)),
            other => other,
        };
        let Some(holding) = holding else {
            return Ok(None);
        };
        pieces.push(holding);
    }
    // And each point: a pole never, a zero where `op` admits equality, a
    // cut by `f`'s sign there; outside the domain never.
    let mut held = Vec::with_capacity(points.len());
    for (p, kind) in &points {
        let holding = match kind {
            Kind::Pole => Some(false),
            _ => match member(ev, &domain, p)? {
                Some(true) if *kind == Kind::Zero => Some(holds(op, Ordering::Equal)),
                Some(true) => sign_at(ev, &f, x, p)?.map(|s| holds(op, s)),
                other => other,
            },
        };
        let Some(holding) = holding else {
            return Ok(None);
        };
        held.push(holding);
    }
    let points: Vec<Expr> = points.into_iter().map(|(p, _)| p).collect();
    Ok(Some(join(&points, &pieces, &held)))
}

/// Whether the sign chart applies to `e` (CAS-279): every part of it in
/// `x` a sum, a product, a power, an absolute value, an exponential or a
/// logarithm, continuous where defined. Periodic functions are not.
fn charted(e: &Expr, x: &Expr) -> bool {
    e.distinct().all(|s| {
        s == x
            || !s.contains(x)
            || matches!(
                (s.head().and_then(Expr::as_symbol), s.args()),
                (Some(Symbol::PLUS | Symbol::TIMES), _)
                    | (Some(Symbol::POWER), [_, _])
                    | (Some(Symbol::ABS | Symbol::EXP | Symbol::LOG), [_])
            )
    })
}

/// Where `l` and `r` are real (CAS-279): a power `u^q` of something in
/// `x`, `q` a rational but no integer, needs `u >= 0`, `u > 0` for
/// `q < 0`, unless `q`'s denominator is odd in the Real domain, where odd
/// roots of negatives are real (O7; giac's `x^(1/3)` is not); a logarithm
/// `log(u)` needs `u > 0`. The whole line when nothing does.
fn domain(
    ev: &mut Eval<'_>,
    e: &Expr,
    l: &Expr,
    r: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Interval>>, EvalError> {
    let real = ev.context().domain() == Domain::Real;
    let mut conditions: Vec<Expr> = Vec::new();
    for s in l.distinct().chain(r.distinct()) {
        let condition = match (s.head().and_then(Expr::as_symbol), s.args()) {
            (Some(Symbol::POWER), [u, q]) if u.contains(x) => match q.to_rational() {
                // An odd root of a negative is real in the Real domain (O7).
                Some(q) if q.is_integer() || real && q.denom().unsigned_abs().bit(0) => continue,
                Some(q) => {
                    let op = if q.is_negative() {
                        Symbol::GREATER
                    } else {
                        Symbol::GREATER_EQUAL
                    };
                    Expr::call(op, [u.clone(), Expr::from(0)])
                }
                None => continue,
            },
            (Some(Symbol::LOG), [u]) if u.contains(x) => {
                Expr::call(Symbol::GREATER, [u.clone(), Expr::from(0)])
            }
            _ => continue,
        };
        if !conditions.contains(&condition) {
            conditions.push(condition);
        }
    }
    if conditions.is_empty() {
        return Ok(Some(vec![Interval { lo: None, hi: None }]));
    }
    set(ev, e, &Expr::call(Symbol::AND, conditions), x)
}

/// The poles of `l`, of `r` and of their functions' arguments (CAS-279),
/// where the sign may change without a zero: `exp(1/x)` at 0.
fn singular(
    ev: &mut Eval<'_>,
    e: &Expr,
    l: &Expr,
    r: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut arguments: Vec<&Expr> = vec![l, r];
    for s in l.distinct().chain(r.distinct()) {
        match (s.head().and_then(Expr::as_symbol), s.args()) {
            (Some(Symbol::ABS | Symbol::EXP | Symbol::LOG), [u]) => arguments.push(u),
            (Some(Symbol::POWER), [b, q]) if !q.to_rational().is_some_and(|q| q.is_integer()) => {
                arguments.extend([b, q]);
            }
            _ => {}
        }
    }
    let mut poles = Vec::new();
    for u in arguments {
        if !u.contains(x) {
            continue;
        }
        let together = ev.eval(&Expr::call(Symbol::TOGETHER, [u.clone()]))?;
        let d = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
        let Some(c) = crossings(ev, e, &d, x)? else {
            return Ok(None);
        };
        poles.extend(c.zeros);
    }
    Ok(Some(poles))
}

/// The crossings of `f` in `x` (CAS-279): of a product, its factors'; of
/// `u^q` with `q > 0` and of `|u|`, `u`'s; of an exponential, none; with
/// absolute values of something in `x` in it, by [`piecewise`]; else the
/// roots `Solve` finds. `None` when they are not found.
fn crossings(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
) -> Result<Option<Crossings>, EvalError> {
    if !f.contains(x) {
        return Ok(Some(Crossings::default()));
    }
    match (f.head().and_then(Expr::as_symbol), f.args()) {
        (Some(Symbol::TIMES), factors) => {
            let mut all = Crossings::default();
            for g in factors {
                let Some(c) = crossings(ev, e, g, x)? else {
                    return Ok(None);
                };
                all.zeros.extend(c.zeros);
                all.cuts.extend(c.cuts);
            }
            return Ok(Some(all));
        }
        (Some(Symbol::POWER), [u, q]) if q.to_rational().is_some_and(|q| q.is_positive()) => {
            return crossings(ev, e, u, x);
        }
        (Some(Symbol::ABS), [u]) => return crossings(ev, e, u, x),
        (Some(Symbol::EXP), [_]) => return Ok(Some(Crossings::default())),
        _ => {}
    }
    if f.distinct().any(|s| absolute(s, x).is_some()) {
        return piecewise(ev, e, f, x);
    }
    Ok(roots(ev, e, f, x, true)?.map(|zeros| Crossings {
        zeros,
        cuts: Vec::new(),
    }))
}

/// The argument of `|u|`, `u` in `x`.
fn absolute<'a>(s: &'a Expr, x: &Expr) -> Option<&'a Expr> {
    match (s.head()?.as_symbol()?, s.args()) {
        (Symbol::ABS, [u]) if u.contains(x) => Some(u),
        _ => None,
    }
}

/// The crossings of `f`, which holds absolute values of expressions in
/// `x`: the zeros and poles of their arguments cut the line, and on each
/// piece, each `|u|` taken as `u` or `−u` by `u`'s sign at a point of it,
/// `f`'s zeros inside the piece are `f`'s. The cuts are crossings too, as
/// `f` may be 0 along a piece (`|x| − x` for `x > 0`).
fn piecewise(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    x: &Expr,
) -> Result<Option<Crossings>, EvalError> {
    let arguments: Vec<Expr> = f
        .distinct()
        .filter_map(|s| absolute(s, x))
        .cloned()
        .collect();
    let mut cuts = Vec::new();
    for u in arguments {
        let u = ev.eval(&Expr::call(Symbol::TOGETHER, [u]))?;
        for part in [Symbol::NUMERATOR, Symbol::DENOMINATOR] {
            let p = ev.eval(&Expr::call(part, [u.clone()]))?;
            let Some(c) = crossings(ev, e, &p, x)? else {
                return Ok(None);
            };
            cuts.extend(c.zeros.into_iter().chain(c.cuts));
        }
    }
    let Some(cuts) = sorted(ev, cuts.into_iter().map(|c| (c, Kind::Cut)))? else {
        return Ok(None);
    };
    let cuts: Vec<Expr> = cuts.into_iter().map(|(c, _)| c).collect();
    let mut zeros = Vec::new();
    for i in 0..=cuts.len() {
        let (lo, hi) = (i.checked_sub(1).map(|j| &cuts[j]), cuts.get(i));
        let test = between(lo, hi);
        // The outermost absolute values first, until none is left.
        let mut g = f.clone();
        loop {
            let found: Vec<Expr> = g
                .distinct()
                .filter(|s| absolute(s, x).is_some())
                .cloned()
                .collect();
            if found.is_empty() {
                break;
            }
            let mut rules = Vec::with_capacity(found.len());
            for a in found {
                let u = a.args()[0].clone();
                let Some(sign) = sign_at(ev, &u, x, &test)? else {
                    return Ok(None);
                };
                rules.push((a, if sign == Ordering::Less { -u } else { u }));
            }
            g = g.substitute(&rules);
        }
        let g = ev.eval(&g)?;
        if g.distinct().any(|s| absolute(s, x).is_some()) {
            return Ok(None);
        }
        let Some(c) = crossings(ev, e, &g, x)? else {
            return Ok(None);
        };
        for z in c.zeros {
            let after = match lo {
                Some(lo) => compare(ev, lo, &z)?,
                None => Some(Ordering::Less),
            };
            let before = match hi {
                Some(hi) => compare(ev, &z, hi)?,
                None => Some(Ordering::Less),
            };
            match (after, before) {
                (Some(Ordering::Less), Some(Ordering::Less)) => zeros.push(z),
                (Some(_), Some(_)) => {}
                _ => return Ok(None),
            }
        }
    }
    Ok(Some(Crossings { zeros, cuts }))
}

/// `points` in increasing order, equal ones merged into the greatest kind;
/// `None` when an order is not decided.
fn sorted(
    ev: &mut Eval<'_>,
    points: impl IntoIterator<Item = (Expr, Kind)>,
) -> Result<Option<Vec<(Expr, Kind)>>, EvalError> {
    let mut out: Vec<(Expr, Kind)> = Vec::new();
    for (p, kind) in points {
        // By insertion.
        let mut at = Some(out.len());
        for (i, (q, q_kind)) in out.iter_mut().enumerate() {
            match compare(ev, &p, q)? {
                Some(Ordering::Less) => {
                    at = Some(i);
                    break;
                }
                Some(Ordering::Equal) => {
                    *q_kind = (*q_kind).max(kind);
                    at = None;
                    break;
                }
                Some(Ordering::Greater) => {}
                None => return Ok(None),
            }
        }
        if let Some(at) = at {
            out.insert(at, (p, kind));
        }
    }
    Ok(Some(out))
}

/// A point of the open piece from `lo` to `hi`, `None` infinite: the
/// midpoint, or the end ∓ 1, or 0.
fn between(lo: Option<&Expr>, hi: Option<&Expr>) -> Expr {
    match (lo, hi) {
        (None, None) => Expr::from(0),
        (None, Some(hi)) => hi.clone() - Expr::from(1),
        (Some(lo), None) => lo.clone() + Expr::from(1),
        (Some(lo), Some(hi)) => (lo.clone() + hi.clone()) / Expr::from(2),
    }
}

/// Whether the point `p` lies in `set`.
fn member(ev: &mut Eval<'_>, set: &[Interval], p: &Expr) -> Result<Option<bool>, EvalError> {
    let end = Some(End {
        at: p.clone(),
        closed: true,
    });
    let point = [Interval {
        lo: end.clone(),
        hi: end,
    }];
    Ok(intersection(ev, set, &point)?.map(|s| !s.is_empty()))
}

/// The intervals of the pieces and points that hold, adjacent ones joined:
/// piece `i` lies between point `i − 1` and point `i`.
fn join(points: &[Expr], pieces: &[bool], held: &[bool]) -> Vec<Interval> {
    let mut out = Vec::new();
    let mut lo: Option<Option<End>> = None;
    // Along the line: piece 0, point 0, piece 1, …, point k−1, piece k.
    for i in 0..=points.len() {
        if pieces[i] {
            if lo.is_none() {
                lo = Some(i.checked_sub(1).map(|j| End {
                    at: points[j].clone(),
                    closed: false,
                }));
            }
        } else if let Some(start) = lo.take() {
            out.push(Interval {
                lo: start,
                hi: Some(End {
                    at: points[i - 1].clone(),
                    closed: held[i - 1],
                }),
            });
        }
        let Some(p) = points.get(i) else {
            break;
        };
        if held[i] {
            if lo.is_none() {
                lo = Some(Some(End {
                    at: p.clone(),
                    closed: true,
                }));
            }
        } else if let Some(start) = lo.take() {
            out.push(Interval {
                lo: start,
                hi: Some(End {
                    at: p.clone(),
                    closed: false,
                }),
            });
        }
    }
    if let Some(start) = lo {
        out.push(Interval {
            lo: start,
            hi: None,
        });
    }
    out
}

/// An interval as giac's condition on `x`.
fn render(i: &Interval, x: &Expr) -> Expr {
    let below = |end: &End| {
        if end.closed {
            Symbol::LESS_EQUAL
        } else {
            Symbol::LESS
        }
    };
    match (&i.lo, &i.hi) {
        (None, None) => x.clone(),
        (None, Some(hi)) => Expr::call(below(hi), [x.clone(), hi.at.clone()]),
        (Some(lo), None) => {
            let above = if lo.closed {
                Symbol::GREATER_EQUAL
            } else {
                Symbol::GREATER
            };
            Expr::call(above, [x.clone(), lo.at.clone()])
        }
        (Some(lo), Some(hi)) if lo.closed && hi.closed && lo.at == hi.at => lo.at.clone(),
        (Some(lo), Some(hi)) => Expr::call(
            Symbol::AND,
            [
                Expr::call(below(lo), [lo.at.clone(), x.clone()]),
                Expr::call(below(hi), [x.clone(), hi.at.clone()]),
            ],
        ),
    }
}

#[cfg(test)]
mod tests;
