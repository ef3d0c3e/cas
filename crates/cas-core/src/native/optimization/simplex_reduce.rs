//! giac's `simplex_reduce` (CAS-606), a port of misc.cc `simplex_reduce`
//! and `_simplex_reduce` (:5878-6017): maximize `c·x` over `A x ≤ b`,
//! `x ≥ 0`, on giac's own tableau.
//!
//! - `SimplexReduce[A, b, c]` makes the tableau `[A | I | b ; −c | 0 | 0]`;
//!   `SimplexReduce[M]` takes the tableau `M` as given. The answer is
//!   giac's `optimum, bfs, M'` as `Sequence[opt, List[…], List[List[…]…]]`:
//!   the optimum, the basic feasible solution over every column but the
//!   last, and the reduced tableau; an unbounded program is
//!   `Sequence[Infinity, List[], M']`, as `LPSolve`'s `Infinity`.
//! - giac's rules, which differ from [`super::simplex`]: the entering
//!   column is the most negative cost (the first of the least), the leaving
//!   row the least ratio `b_i / a_iJ` over `a_iJ > 0`, the lowest row on a
//!   tie; after a zero ratio the entering column is the first negative
//!   cost (Bland's), never before. Columns that are a unit vector in the
//!   constraint rows are made basic first by clearing their cost. There is
//!   no first phase, and the raw tableau is returned.
//! - Entries are exact or decimal as giac's gens: exact with exact stays
//!   exact, a decimal makes the result decimal. Decimals are read as
//!   doubles; an entry not a finite number leaves the call.
//! - Departures: giac pivots from the tableau even with a negative `b`,
//!   whose basis is then infeasible and its answer no optimum; ours warns
//!   and leaves the call. giac's dimension error (fewer columns than rows
//!   plus one) and any other shape leave the call.
//! - Limits (C0019): a tableau of more cells than the `terms` limit leaves
//!   the call before any work; pivots times cells past `match_steps`, or an
//!   exact entry past `number_bits` bits, warns and leaves it. The host's
//!   interrupt is polled at each pivot.

#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    reason = "giac's names (m, nr, nc, I, J) and its index loops, kept to keep its tie rules"
)]

use std::cmp::Ordering;
use std::sync::atomic::Ordering as Atomic;

use cas_num::Rational;
use cas_num::ball::{Ball, Float, Round};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// An entry of the tableau, as giac's gen: exact or a double.
#[derive(Clone, Debug)]
enum Entry {
    Exact(Rational),
    Float(f64),
}

/// The nearest `f64`.
fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

impl Entry {
    fn zero() -> Entry {
        Entry::Exact(Rational::from(0))
    }

    fn double(&self) -> f64 {
        match self {
            Entry::Exact(q) => to_f64(q),
            Entry::Float(f) => *f,
        }
    }

    /// The order of two entries, as numbers; `None` with a NaN.
    fn cmp(&self, other: &Entry) -> Option<Ordering> {
        match (self, other) {
            (Entry::Exact(a), Entry::Exact(b)) => Some(a.cmp(b)),
            _ => self.double().partial_cmp(&other.double()),
        }
    }

    fn is(&self, n: i64) -> bool {
        self.cmp(&Entry::Exact(Rational::from(n))) == Some(Ordering::Equal)
    }

    fn above(&self, other: &Entry) -> bool {
        self.cmp(other) == Some(Ordering::Greater)
    }

    fn sub(&self, other: &Entry) -> Entry {
        match (self, other) {
            (Entry::Exact(a), Entry::Exact(b)) => Entry::Exact(a - b),
            _ => Entry::Float(self.double() - other.double()),
        }
    }

    fn mul(&self, other: &Entry) -> Entry {
        match (self, other) {
            (Entry::Exact(a), Entry::Exact(b)) => Entry::Exact(a * b),
            _ => Entry::Float(self.double() * other.double()),
        }
    }

    /// `self / other`, `other` not zero.
    fn div(&self, other: &Entry) -> Entry {
        match (self, other) {
            (Entry::Exact(a), Entry::Exact(b)) => Entry::Exact(a / b),
            _ => Entry::Float(self.double() / other.double()),
        }
    }

    /// Whether an exact entry has more than `bits` bits a part.
    fn past(&self, bits: u64) -> bool {
        match self {
            Entry::Exact(q) => q.numer().bit_length() > bits || q.denom().bit_length() > bits,
            Entry::Float(_) => false,
        }
    }

    /// The expression; `None` for a double not finite.
    fn expr(&self) -> Option<Expr> {
        match self {
            Entry::Exact(q) => Some(Expr::from(q.clone())),
            Entry::Float(f) => Float::from_f64(*f).map(|x| Expr::from(Ball::exact(x))),
        }
    }
}

/// An entry of the input: a number, or a finite decimal.
fn entry(e: &Expr) -> Option<Entry> {
    match e.view() {
        View::Int(n) => Some(Entry::Exact(Rational::from(n))),
        View::Integer(n) => Some(Entry::Exact(Rational::from(n))),
        View::Rational(q) => Some(Entry::Exact(q.clone())),
        View::Ball(b) if b.is_finite() => {
            Some(Entry::Float(b.mid().to_f64())).filter(|x| x.double().is_finite())
        }
        _ => None,
    }
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// The entries of a list.
fn vector(e: &Expr) -> Option<Vec<Entry>> {
    if !is_list(e) {
        return None;
    }
    e.args().iter().map(entry).collect()
}

/// The rows of a matrix, a list of lists of entries, rectangular and not
/// empty (giac's `ckmatrix`).
fn matrix(e: &Expr) -> Option<Vec<Vec<Entry>>> {
    if !is_list(e) || e.args().is_empty() {
        return None;
    }
    let m: Vec<Vec<Entry>> = e.args().iter().map(vector).collect::<Option<_>>()?;
    let nc = m[0].len();
    (nc > 0 && m.iter().all(|r| r.len() == nc)).then_some(m)
}

/// giac's tableau `[A | I | b ; −c | 0 | 0]`.
fn tableau(a: &Expr, b: &Expr, c: &Expr) -> Option<Vec<Vec<Entry>>> {
    let a = matrix(a)?;
    let (b, c) = (vector(b)?, vector(c)?);
    let (n, k) = (a.len(), a[0].len());
    if b.len() != n || c.len() != k {
        return None;
    }
    let mut m = Vec::with_capacity(n + 1);
    for (i, row) in a.into_iter().enumerate() {
        let mut r = row;
        r.extend((0..n).map(|j| Entry::Exact(Rational::from(i64::from(i == j)))));
        r.push(b[i].clone());
        m.push(r);
    }
    let mut cost: Vec<Entry> = c.iter().map(|x| Entry::zero().sub(x)).collect();
    cost.extend((0..=n).map(|_| Entry::zero()));
    m.push(cost);
    Some(m)
}

/// Why the reduction stopped before its end.
enum Stop {
    Limits,
    Interrupted,
}

/// The answer: the optimum (`None` unbounded), the basic feasible solution
/// and the reduced tableau.
type Reduced = (Option<Entry>, Vec<Entry>, Vec<Vec<Entry>>);

/// `row -= a·pivot`.
fn eliminate(row: &mut [Entry], a: &Entry, pivot: &[Entry]) {
    for (x, y) in row.iter_mut().zip(pivot) {
        *x = x.sub(&a.mul(y));
    }
}

/// The basic feasible solution of an optimal tableau: a column of cost 0
/// with one entry in the constraint rows, a 1, takes that row's right side,
/// up to one per constraint row; the others are 0.
fn solution(m: &[Vec<Entry>]) -> Vec<Entry> {
    let (nr, nc) = (m.len(), m[0].len());
    let mut bfs = vec![Entry::zero(); nc - 1];
    let mut counter = nr - 1;
    for i in 0..nc - 1 {
        if counter == 0 {
            break;
        }
        let mut cur = None;
        if m[nr - 1][i].is(0) {
            for j in 0..nr - 1 {
                if m[j][i].is(0) {
                    continue;
                }
                if cur.is_some() {
                    cur = None;
                    break;
                }
                cur = Some(j);
            }
        }
        if let Some(j) = cur
            && m[j][i].is(1)
        {
            counter -= 1;
            bfs[i] = m[j][nc - 1].clone();
        }
    }
    bfs
}

/// giac's `simplex_reduce` (misc.cc:5878-5982) on a tableau of `nr ≥ 1`
/// rows and `nc ≥ nr + 1` columns.
fn reduce(ev: &Eval<'_>, mut m: Vec<Vec<Entry>>) -> Result<Reduced, Stop> {
    let (nr, nc) = (m.len(), m[0].len());
    let limits = ev.context().limits();
    let cells = u64::try_from(nr * nc).unwrap_or(u64::MAX);
    // Step 0: columns that are a unit vector in the constraint rows become
    // basic, their cost cleared.
    let mut last = vec![Entry::zero(); nr - 1];
    for i in 0..nc - 1 {
        let mut cur = None;
        for j in 0..nr - 1 {
            if m[j][i].is(0) {
                continue;
            }
            if !m[j][i].is(1) || cur.is_some() {
                cur = None;
                break;
            }
            cur = Some(j);
        }
        if let Some(j) = cur {
            last[j] = m[nr - 1][i].clone();
        }
    }
    for i in 0..nr - 1 {
        if !last[i].is(0) {
            let pivot = m[i].clone();
            eliminate(&mut m[nr - 1], &last[i], &pivot);
        }
    }
    let mut first = false;
    let mut steps = 0u64;
    loop {
        let mut col = None;
        let mut mc = Entry::zero();
        for j in 0..nc - 1 {
            if mc.above(&m[nr - 1][j]) {
                col = Some(j);
                mc = m[nr - 1][j].clone();
                if first {
                    break;
                }
            }
        }
        let Some(jj) = col else {
            let opt = m[nr - 1][nc - 1].clone();
            let bfs = solution(&m);
            return Ok((Some(opt), bfs, m));
        };
        let mut row: Option<(usize, Entry)> = None;
        for i in 0..nr - 1 {
            if m[i][jj].above(&Entry::zero()) {
                let r = m[i][nc - 1].div(&m[i][jj]);
                if row.as_ref().is_none_or(|(_, mc)| mc.above(&r)) {
                    row = Some((i, r));
                }
            }
        }
        let Some((ii, mc)) = row else {
            return Ok((None, Vec::new(), m));
        };
        if mc.is(0) {
            first = true;
        }
        if ev.context().interrupt_flag().load(Atomic::Relaxed) {
            return Err(Stop::Interrupted);
        }
        steps = steps.saturating_add(cells);
        if steps > limits.match_steps {
            return Err(Stop::Limits);
        }
        let p = m[ii][jj].clone();
        for x in &mut m[ii] {
            *x = x.div(&p);
        }
        let pivot = m[ii].clone();
        for i in 0..nr {
            if i != ii {
                let a = m[i][jj].clone();
                eliminate(&mut m[i], &a, &pivot);
            }
        }
        if m.iter().flatten().any(|x| x.past(limits.number_bits)) {
            return Err(Stop::Limits);
        }
    }
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &'static str) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

fn list(row: &[Entry]) -> Option<Expr> {
    let items: Vec<Expr> = row.iter().map(Entry::expr).collect::<Option<_>>()?;
    Some(Expr::call(Symbol::LIST, items))
}

/// `SimplexReduce[A, b, c]` and `SimplexReduce[M]`.
pub(in crate::native) fn simplex_reduce(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let m = match e.args() {
        [a, b, c] => tableau(a, b, c),
        [m] => matrix(m),
        _ => None,
    };
    let Some(m) = m else {
        return Ok(None);
    };
    let (nr, nc) = (m.len(), m[0].len());
    if nc < nr + 1 {
        return Ok(None);
    }
    if m[..nr - 1].iter().any(|r| Entry::zero().above(&r[nc - 1])) {
        warn(
            ev,
            e,
            "simplex_reduce: a negative right side, no feasible basis",
        );
        return Ok(None);
    }
    let limits = ev.context().limits();
    if nr.saturating_mul(nc) > limits.terms
        || m.iter().flatten().any(|x| x.past(limits.number_bits))
    {
        warn(ev, e, "simplex_reduce: past the limits");
        return Ok(None);
    }
    let (opt, bfs, m) = match reduce(ev, m) {
        Ok(r) => r,
        Err(Stop::Interrupted) => return Err(EvalError::Interrupted),
        Err(Stop::Limits) => {
            warn(ev, e, "simplex_reduce: past the limits");
            return Ok(None);
        }
    };
    let opt = match opt {
        Some(x) => x.expr(),
        None => Some(Expr::from(Symbol::INFINITY)),
    };
    let rows: Option<Vec<Expr>> = m.iter().map(|r| list(r)).collect();
    let (Some(opt), Some(bfs), Some(rows)) = (opt, list(&bfs), rows) else {
        return Ok(None);
    };
    let result = Expr::call(Symbol::SEQUENCE, [opt, bfs, Expr::call(Symbol::LIST, rows)]);
    ev.step("simplex_reduce", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
