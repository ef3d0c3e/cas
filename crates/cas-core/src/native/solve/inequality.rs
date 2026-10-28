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
//! complement; an equation among them is its real roots, points. Each set
//! is disjoint intervals in increasing order, ends compared as the sign
//! chart's points are.

use std::cmp::Ordering;

use super::super::rounding::real_sign;
use super::roots;
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

/// Whether a sign satisfies the relation `op` against 0.
fn holds(op: Symbol, sign: Ordering) -> bool {
    match op {
        Symbol::LESS => sign == Ordering::Less,
        Symbol::LESS_EQUAL => sign != Ordering::Greater,
        Symbol::GREATER => sign == Ordering::Greater,
        Symbol::GREATER_EQUAL => sign != Ordering::Less,
        _ => sign != Ordering::Equal,
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
/// a comparison by its sign chart, an equation by its real roots, a list
/// or `And` the intersection of its members', `Or` the union, `Not` the
/// complement (CAS-278).
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
        return points(ev, e, l, r, x);
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

/// The real roots of `l = r` as points.
fn points(
    ev: &mut Eval<'_>,
    e: &Expr,
    l: &Expr,
    r: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Interval>>, EvalError> {
    let f = ev.eval(&(l.clone() - r.clone()))?;
    let Some(found) = roots(ev, e, &f, x, true)? else {
        return Ok(None);
    };
    if found.iter().any(|p| p.contains(x)) {
        // f ≡ 0: every x.
        return Ok(Some(vec![Interval { lo: None, hi: None }]));
    }
    let points = found
        .into_iter()
        .map(|p| Interval {
            lo: Some(End {
                at: p.clone(),
                closed: true,
            }),
            hi: Some(End {
                at: p,
                closed: true,
            }),
        })
        .collect();
    union(ev, points)
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

/// The solutions of `l op r` by a sign chart.
fn one(
    ev: &mut Eval<'_>,
    e: &Expr,
    op: Symbol,
    l: &Expr,
    r: &Expr,
    x: &Expr,
) -> Result<Option<Vec<Interval>>, EvalError> {
    let f = ev.eval(&Expr::call(Symbol::TOGETHER, [l.clone() - r.clone()]))?;
    let whole = || vec![Interval { lo: None, hi: None }];
    if !f.contains(x) {
        return Ok(
            real_sign(ev.context(), &f).map(|s| if holds(op, s) { whole() } else { Vec::new() })
        );
    }
    let n = ev.eval(&Expr::call(Symbol::NUMERATOR, [f.clone()]))?;
    let d = ev.eval(&Expr::call(Symbol::DENOMINATOR, [f.clone()]))?;
    // The critical points: zeros, then poles.
    let Some(zeros) = roots(ev, e, &n, x, true)? else {
        return Ok(None);
    };
    let poles = if d.contains(x) {
        let Some(p) = roots(ev, e, &d, x, true)? else {
            return Ok(None);
        };
        p
    } else {
        Vec::new()
    };
    // A root that is `x` itself: every x, as for f ≡ 0, which `Together`
    // already took.
    if zeros.iter().chain(&poles).any(|p| p.contains(x)) {
        return Ok(None);
    }
    let mut points: Vec<(Expr, bool)> = Vec::new();
    for (p, zero) in zeros
        .into_iter()
        .map(|p| (p, true))
        .chain(poles.into_iter().map(|p| (p, false)))
    {
        // Sorted by insertion, equal points merged (a pole wins).
        let mut at = points.len();
        for (i, (q, q_zero)) in points.iter_mut().enumerate() {
            match compare(ev, &p, q)? {
                Some(Ordering::Less) => {
                    at = i;
                    break;
                }
                Some(Ordering::Equal) => {
                    *q_zero = *q_zero && zero;
                    at = usize::MAX;
                    break;
                }
                Some(Ordering::Greater) => {}
                None => return Ok(None),
            }
        }
        if at != usize::MAX {
            points.insert(at, (p, zero));
        }
    }
    // Whether each piece, and each point, holds: piece i lies below point i.
    let mut pieces = Vec::with_capacity(points.len() + 1);
    for i in 0..=points.len() {
        let test = match (
            i.checked_sub(1).map(|j| &points[j].0),
            points.get(i).map(|(p, _)| p),
        ) {
            (None, None) => Expr::from(0),
            (None, Some(hi)) => hi.clone() - Expr::from(1),
            (Some(lo), None) => lo.clone() + Expr::from(1),
            (Some(lo), Some(hi)) => (lo.clone() + hi.clone()) / Expr::from(2),
        };
        let Some(sign) = sign_at(ev, &f, x, &test)? else {
            return Ok(None);
        };
        pieces.push(holds(op, sign));
    }
    let held: Vec<bool> = points
        .iter()
        .map(|&(_, zero)| zero && holds(op, Ordering::Equal))
        .collect();
    Ok(Some(join(&points, &pieces, &held)))
}

/// The intervals of the pieces and points that hold, adjacent ones joined:
/// piece `i` lies between point `i − 1` and point `i`.
fn join(points: &[(Expr, bool)], pieces: &[bool], held: &[bool]) -> Vec<Interval> {
    let mut out = Vec::new();
    let mut lo: Option<Option<End>> = None;
    // Along the line: piece 0, point 0, piece 1, …, point k−1, piece k.
    for i in 0..=points.len() {
        if pieces[i] {
            if lo.is_none() {
                lo = Some(i.checked_sub(1).map(|j| End {
                    at: points[j].0.clone(),
                    closed: false,
                }));
            }
        } else if let Some(start) = lo.take() {
            out.push(Interval {
                lo: start,
                hi: Some(End {
                    at: points[i - 1].0.clone(),
                    closed: held[i - 1],
                }),
            });
        }
        let Some((p, _)) = points.get(i) else {
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
