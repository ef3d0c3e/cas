//! Transportation problems (CAS-598): `TPSolve`, Xcas's `tpsolve`.
//!
//! `TPSolve[supply, demand, cost]`:
//!
//! - `supply` and `demand` are lists of nonnegative integers, `cost` a
//!   matrix of their sizes with nonnegative rational entries, or with one
//!   symbol `M` (a forbidden route, priced `1/ε`) among them. Anything else
//!   (decimals, negatives, `2*M`, two symbols, a wrong shape, empty lists)
//!   leaves the call.
//! - The algorithm is giac's (optimization.cc `_tpsolve`, :3253): a dummy
//!   row or column balances the totals, `tprob::north_west_corner` gives a
//!   first plan, degenerate cells held by `ε`, and `tprob::modi` improves
//!   it along `tprob::stepping_stone_path` until no reduced cost is
//!   negative. The answer is giac's `cost, plan`: the dummy dropped, the
//!   cost `Σ P·X` over the real cells.
//! - Departures: exact rationals throughout, `ε = 1/(10¹² (m′n′ + 1))` on
//!   the sizes with the dummy; the potentials from a walk of the basis tree
//!   (`u₀ = 0`); each cell rounded to the nearest integer at the end, for
//!   giac's `epsilon2zero` and `exact`. A plan that uses an `M` cell, a
//!   basis that does not span, an empty path or more than [`MAX_PIVOTS`]
//!   steps warns that the problem has no feasible solution and leaves the
//!   call, as giac errs.
//! - Caps, before any work: at most [`MAX_CELLS`] cells with the dummy,
//!   and each supply, demand, numerator and denominator of at most
//!   [`MAX_BITS`] bits.

#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    reason = "giac's names (supply s, demand d, plan X, costs P, potentials U, V) and its \
              index loops, kept to keep its tie rules"
)]

use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// Most cells of the balanced problem.
const MAX_CELLS: usize = 1024;
/// Most MODI steps.
const MAX_PIVOTS: usize = 4 * MAX_CELLS;
/// Most bits of a supply, a demand, or a cost's numerator or denominator.
const MAX_BITS: u64 = 1024;

/// A problem with no feasible solution, as far as this method sees.
struct Infeasible;

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// A nonnegative integer of at most [`MAX_BITS`] bits.
fn quantity(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Int(n) if n >= 0 => Some(Rational::from(n)),
        View::Integer(n) if !n.is_negative() && n.bit_length() <= MAX_BITS => {
            Some(Rational::from(n))
        }
        _ => None,
    }
}

/// A nonnegative rational of at most [`MAX_BITS`] bits a part.
fn price(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Rational(q)
            if !q.is_negative()
                && q.numer().bit_length() <= MAX_BITS
                && q.denom().bit_length() <= MAX_BITS =>
        {
            Some(q.clone())
        }
        _ => quantity(e),
    }
}

/// The north-west corner plan (giac's `tprob::north_west_corner`).
fn north_west(s: &[Rational], d: &[Rational], eps: &Rational) -> Vec<Vec<Rational>> {
    let (m, n) = (s.len(), d.len());
    let tiny = Rational::new(Integer::from(1), Integer::from(10).pow(12));
    let kept = |x: &Rational| {
        if x.abs() >= tiny {
            x.clone()
        } else {
            Rational::from(0)
        }
    };
    let mut x = vec![vec![Rational::from(0); n]; m];
    let (mut i, mut j) = (0, 0);
    while i < m && j < n {
        let v = (0..i).fold(Rational::from(0), |t, k| t + kept(&x[k][j]));
        let u = (0..j).fold(Rational::from(0), |t, k| t + kept(&x[i][k]));
        let a = (&s[i] - &u).min(&d[j] - &v);
        x[i][j] = a.clone();
        let k = i + j;
        if &u + &a == s[i] {
            i += 1;
        }
        if &v + &a == d[j] {
            j += 1;
        }
        if i < m && j < n && i + j - k == 2 {
            x[i - 1][j] = eps.clone();
        }
    }
    x
}

/// The stepping-stone path from `path[0]` (giac's
/// `tprob::stepping_stone_path`), or `false` and `path` as it came.
fn stepping_stone(path: &mut Vec<(usize, usize)>, x: &[Vec<Rational>]) -> bool {
    let (row, column) = path[path.len() - 1];
    let (m, n) = (x.len(), x[0].len());
    if path.len() > 1 && path[0].1 == column {
        return true;
    }
    let horizontal = path.len() % 2 == 1;
    for i in 0..(if horizontal { n } else { m }) {
        let count = path
            .iter()
            .filter(|&&(a, b)| if horizontal { b == i } else { a == i })
            .count();
        let (r, c) = if horizontal { (row, i) } else { (i, column) };
        if count < 2 && !x[r][c].is_zero() {
            path.push((r, c));
            if stepping_stone(path, x) {
                return true;
            }
            path.pop();
        }
    }
    false
}

/// The potentials `u`, `v` with `u₀ = 0` and `u_i + v_j = P_ij` on the
/// basis, by a walk of its tree.
fn potentials(
    p: &[Vec<Rational>],
    x: &[Vec<Rational>],
) -> Result<(Vec<Rational>, Vec<Rational>), Infeasible> {
    let (m, n) = (x.len(), x[0].len());
    let mut u: Vec<Option<Rational>> = vec![None; m];
    let mut v: Vec<Option<Rational>> = vec![None; n];
    u[0] = Some(Rational::from(0));
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..m {
            for j in 0..n {
                if x[i][j].is_zero() {
                    continue;
                }
                match (&u[i], &v[j]) {
                    (Some(ui), None) => {
                        v[j] = Some(&p[i][j] - ui);
                        changed = true;
                    }
                    (None, Some(vj)) => {
                        u[i] = Some(&p[i][j] - vj);
                        changed = true;
                    }
                    _ => {}
                }
            }
        }
    }
    let u = u
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(Infeasible)?;
    let v = v
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(Infeasible)?;
    Ok((u, v))
}

/// The optimal plan from `x` (giac's `tprob::modi`).
fn modi(p: &[Vec<Rational>], x: &mut [Vec<Rational>], eps: &Rational) -> Result<(), Infeasible> {
    let (m, n) = (x.len(), x[0].len());
    let half = Rational::new(Integer::from(1), Integer::from(2));
    for _ in 0..MAX_PIVOTS {
        let (u, v) = potentials(p, x)?;
        let mut least = Rational::from(0);
        let mut best = None;
        for i in 0..m {
            for j in 0..n {
                if x[i][j].is_zero() {
                    let c = &p[i][j] - &u[i] - &v[j];
                    if least > c {
                        least = c;
                        best = Some((i, j));
                    }
                }
            }
        }
        let Some(best) = best else {
            return Ok(());
        };
        let mut path = vec![best];
        if !stepping_stone(&mut path, x) {
            return Err(Infeasible);
        }
        let mut d = x[path[1].0][path[1].1].clone();
        for &(a, b) in path.iter().skip(3).step_by(2) {
            if x[a][b] < d {
                d = x[a][b].clone();
            }
        }
        for (k, &(a, b)) in path.iter().enumerate() {
            let mut y = if k % 2 == 1 {
                &x[a][b] - &d
            } else {
                &x[a][b] + &d
            };
            let has_zero = path.iter().any(|&(r, c)| x[r][c].is_zero());
            if (!y.is_zero() && half > y) || (y.is_zero() && has_zero) {
                y = eps.clone();
            }
            x[a][b] = y;
        }
    }
    Err(Infeasible)
}

/// The cost matrix, rows of `n` prices, `None` for the one symbol `M`.
fn costs(cost: &Expr, n: usize) -> Option<Vec<Vec<Option<Rational>>>> {
    let mut forbidden: Option<&Expr> = None;
    let mut p = Vec::with_capacity(cost.args().len());
    for row in cost.args() {
        if !is_list(row) || row.args().len() != n {
            return None;
        }
        let mut cells = Vec::with_capacity(n);
        for c in row.args() {
            if let Some(q) = price(c) {
                cells.push(Some(q));
            } else if c.as_symbol().is_some_and(|x| !x.is_builtin())
                && forbidden.is_none_or(|f| f == c)
            {
                forbidden = Some(c);
                cells.push(None);
            } else {
                return None;
            }
        }
        p.push(cells);
    }
    Some(p)
}

/// `TPSolve[supply, demand, cost]`.
pub(in crate::native) fn tp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [supply, demand, cost] = e.args() else {
        return Ok(None);
    };
    if !is_list(supply) || !is_list(demand) || !is_list(cost) {
        return Ok(None);
    }
    let (Some(mut s), Some(mut d)) = (
        supply
            .args()
            .iter()
            .map(quantity)
            .collect::<Option<Vec<_>>>(),
        demand
            .args()
            .iter()
            .map(quantity)
            .collect::<Option<Vec<_>>>(),
    ) else {
        return Ok(None);
    };
    let (m, n) = (s.len(), d.len());
    if m == 0 || n == 0 || cost.args().len() != m {
        return Ok(None);
    }
    let Some(mut p) = costs(cost, n) else {
        return Ok(None);
    };
    let ts = s.iter().fold(Rational::from(0), |t, q| t + q);
    let td = d.iter().fold(Rational::from(0), |t, q| t + q);
    let cells = (m + usize::from(td > ts)) * (n + usize::from(ts > td));
    if cells > MAX_CELLS {
        return Ok(None);
    }
    let zero = || Some(Rational::from(0));
    if ts > td {
        d.push(&ts - &td);
        for row in &mut p {
            row.push(zero());
        }
    } else if td > ts {
        s.push(&td - &ts);
        p.push(vec![zero(); n]);
    }
    let eps = Rational::new(
        Integer::from(1),
        Integer::from(10).pow(12) * Integer::from(cells as u64 + 1),
    );
    let big = eps.recip();
    let priced: Vec<Vec<Rational>> = p
        .iter()
        .map(|row| {
            row.iter()
                .map(|c| c.clone().unwrap_or_else(|| big.clone()))
                .collect()
        })
        .collect();
    let mut x = north_west(&s, &d, &eps);
    let solved = modi(&priced, &mut x, &eps);
    let mut total = Rational::from(0);
    let mut plan = Vec::with_capacity(m);
    let mut feasible = solved.is_ok();
    for i in 0..m {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            let r = x[i][j].round_ties_even();
            if !r.is_zero() {
                match &p[i][j] {
                    Some(q) => total += q * &Rational::from(&r),
                    None => feasible = false,
                }
            }
            row.push(Expr::from(r));
        }
        plan.push(Expr::call(Symbol::LIST, row));
    }
    if !feasible {
        ev.warn(Warning::Native {
            message: "the problem has no feasible solution".to_owned().into(),
            expr: e.clone(),
        });
        return Ok(None);
    }
    let result = Expr::call(
        Symbol::SEQUENCE,
        [Expr::from(total), Expr::call(Symbol::LIST, plan)],
    );
    ev.step("tpsolve", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
