//! The simplex method on rationals (CAS-DOC-38 §4, CAS-417): minimize `c·x`
//! over `x ≥ 0` and rows `a·x ≤ b`, `= b` or `≥ b`, exactly.
//!
//! - Rows with `b < 0` are negated. A `≤` row takes a slack variable, a `≥`
//!   row a surplus and an artificial one, an `=` row an artificial one;
//!   the slacks and artificials make the first basis.
//! - The first phase minimizes the artificials' sum: above 0 the program
//!   has no feasible point. Artificials left in the basis at 0 are pivoted
//!   out, or their rows, redundant, dropped.
//! - The second phase minimizes `c·x`, artificials out. A column that would
//!   enter with no positive entry makes the program unbounded.
//! - Bland's rule, the entering and leaving variables of the least index,
//!   keeps degenerate pivots from cycling.
//! - [`Work`] counts the iterations of a whole call, every program it
//!   solves, as giac's `icount` (lpsolve.cc:392-401): each search for an
//!   entering column, the one that finds the optimum too, in both phases.
//!   Past the user's `lp_iterationlimit` the method stops
//!   ([`Outcome::Limit`]), past [`MAX_PIVOTS`] it gives up; it polls the
//!   host's interrupt at each iteration, and a tableau of more cells than
//!   the evaluator's `terms` limit gives up before any work.
//! - The optimal tableau is kept ([`Final`]) for the cuts of
//!   [`super::lpsolve::relaxation`].

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the simplex method: c·x, n variables, m rows, i and j"
)]

use std::sync::atomic::{AtomicBool, Ordering};

use cas_num::Rational;

/// The most simplex iterations of one call, every relaxation of a branch
/// and bound together (was the most pivots of one program, also `100_000`).
pub(super) const MAX_PIVOTS: usize = 100_000;

/// A row's relation.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Relation {
    AtMost,
    Equal,
    AtLeast,
}

/// What a program has.
#[derive(Debug, PartialEq, Eq)]
pub(super) enum Outcome {
    /// A best point and its value.
    Optimal(Vec<Rational>, Rational),
    Unbounded,
    Infeasible,
    /// Past [`MAX_PIVOTS`], or past the limits.
    GaveUp,
    /// Past the user's `lp_iterationlimit`.
    Limit,
    /// The host asked to stop.
    Interrupted,
}

/// The work of one call: giac's simplex iteration count (`icount`,
/// lpsolve.cc:395-401) over every program it solves, its cap, the tableau
/// cells allowed, and the host's interrupt.
pub(super) struct Work<'a> {
    /// The iterations so far.
    pub(super) iterations: usize,
    /// The most iterations: the user's limit, at most [`MAX_PIVOTS`], or
    /// [`MAX_PIVOTS`].
    limit: usize,
    /// Whether `limit` is the user's.
    user: bool,
    /// The most cells of a tableau.
    cells: usize,
    interrupt: &'a AtomicBool,
}

impl<'a> Work<'a> {
    /// The work of a call with the user's iteration limit `user` (0 for
    /// none), tableaux of at most `cells` cells, and the host's
    /// `interrupt`.
    pub(super) fn new(user: usize, cells: usize, interrupt: &'a AtomicBool) -> Work<'a> {
        Work {
            iterations: 0,
            limit: if user > 0 {
                user.min(MAX_PIVOTS)
            } else {
                MAX_PIVOTS
            },
            user: user > 0,
            cells,
            interrupt,
        }
    }

    /// Whether the host asked to stop.
    pub(super) fn interrupted(&self) -> bool {
        self.interrupt.load(Ordering::Relaxed)
    }

    /// One iteration more; the outcome that stops the call past the limit
    /// (giac: `++icount; if (limit>0 && icount>limit) break;`) or on an
    /// interrupt.
    fn tick(&mut self) -> Result<(), Outcome> {
        if self.interrupted() {
            return Err(Outcome::Interrupted);
        }
        self.iterations += 1;
        if self.iterations > self.limit {
            return Err(if self.user {
                Outcome::Limit
            } else {
                Outcome::GaveUp
            });
        }
        Ok(())
    }
}

/// The optimal tableau of a program: its rows `x_B + Σ a_k x_k = b`, each
/// with its basic column, over the program's columns, then one slack or
/// surplus per inequality, then the artificials, the right side last.
pub(super) struct Final {
    pub(super) rows: Vec<Vec<Rational>>,
    pub(super) basis: Vec<usize>,
    /// The first artificial column, after the program's and the slacks.
    pub(super) real: usize,
    /// The program's row of each slack or surplus column, in order.
    pub(super) slack_rows: Vec<usize>,
}

/// A row `a·x (relation) b`.
pub(super) type Row = (Vec<Rational>, Relation, Rational);

/// A tableau in canonical form: each row's basic variable has a unit
/// column; the last entry of a row is its right side.
struct Tableau {
    rows: Vec<Vec<Rational>>,
    basis: Vec<usize>,
    columns: usize,
}

impl Tableau {
    fn pivot(&mut self, r: usize, c: usize) {
        let p = self.rows[r][c].clone();
        for v in &mut self.rows[r] {
            *v = &*v / &p;
        }
        let pivot_row = self.rows[r].clone();
        for (i, row) in self.rows.iter_mut().enumerate() {
            if i == r || row[c].is_zero() {
                continue;
            }
            let factor = row[c].clone();
            for (v, w) in row.iter_mut().zip(&pivot_row) {
                *v = &*v - &(&factor * w);
            }
        }
        self.basis[r] = c;
    }

    /// The reduced cost of column `j` for `cost`.
    fn reduced(&self, cost: &[Rational], j: usize) -> Rational {
        let mut r = cost[j].clone();
        for (row, &b) in self.rows.iter().zip(&self.basis) {
            r = &r - &(&cost[b] * &row[j]);
        }
        r
    }

    /// Minimizes `cost` over the columns `allowed` to enter; `Ok(false)`
    /// for an unbounded program, `Err` with the outcome that stops the
    /// call.
    fn minimize(
        &mut self,
        cost: &[Rational],
        allowed: &[bool],
        work: &mut Work<'_>,
    ) -> Result<bool, Outcome> {
        let rhs = self.columns;
        loop {
            work.tick()?;
            let entering = (0..self.columns).find(|&j| {
                allowed[j] && !self.basis.contains(&j) && self.reduced(cost, j).is_negative()
            });
            let Some(j) = entering else {
                return Ok(true);
            };
            let mut leaving: Option<(usize, Rational)> = None;
            for (i, row) in self.rows.iter().enumerate() {
                if !row[j].is_positive() {
                    continue;
                }
                let ratio = &row[rhs] / &row[j];
                let better = match &leaving {
                    None => true,
                    Some((k, best)) => {
                        ratio < *best || ratio == *best && self.basis[i] < self.basis[*k]
                    }
                };
                if better {
                    leaving = Some((i, ratio));
                }
            }
            let Some((i, _)) = leaving else {
                return Ok(false);
            };
            self.pivot(i, j);
        }
    }
}

/// The first tableau of `n` variables and the rows, their right sides made
/// nonnegative: the program's row of each slack or surplus too.
fn tableau(n: usize, rows: &[Row]) -> (Tableau, Vec<usize>) {
    let rows: Vec<Row> = rows
        .iter()
        .map(|(a, rel, b)| {
            if b.is_negative() {
                let flipped = match rel {
                    Relation::AtMost => Relation::AtLeast,
                    Relation::Equal => Relation::Equal,
                    Relation::AtLeast => Relation::AtMost,
                };
                (a.iter().map(|v| -v.clone()).collect(), flipped, -b.clone())
            } else {
                (a.clone(), *rel, b.clone())
            }
        })
        .collect();
    // Columns: the variables, then one slack or surplus per inequality,
    // then one artificial per `≥` or `=` row.
    let slacks = rows
        .iter()
        .filter(|(_, r, _)| *r != Relation::Equal)
        .count();
    let artificials = rows
        .iter()
        .filter(|(_, r, _)| *r != Relation::AtMost)
        .count();
    let columns = n + slacks + artificials;
    let mut table = Tableau {
        rows: Vec::with_capacity(rows.len()),
        basis: Vec::with_capacity(rows.len()),
        columns,
    };
    let mut slack_rows = Vec::with_capacity(slacks);
    let (mut slack, mut artificial) = (n, n + slacks);
    for (k, (a, rel, b)) in rows.iter().enumerate() {
        let mut row = vec![Rational::ZERO; columns + 1];
        row[..n].clone_from_slice(a);
        row[columns] = b.clone();
        if *rel != Relation::Equal {
            row[slack] = if *rel == Relation::AtMost {
                Rational::ONE
            } else {
                -Rational::ONE
            };
            slack_rows.push(k);
            slack += 1;
        }
        if *rel == Relation::AtMost {
            table.basis.push(slack - 1);
        } else {
            row[artificial] = Rational::ONE;
            table.basis.push(artificial);
            artificial += 1;
        }
        table.rows.push(row);
    }
    (table, slack_rows)
}

/// The first phase: the artificials' sum minimized, then those left in the
/// basis at 0 pivoted out, or their rows, redundant, dropped. `Ok(false)`
/// for no feasible point, `Err` with the outcome that stops the call.
fn phase_one(table: &mut Tableau, real: usize, work: &mut Work<'_>) -> Result<bool, Outcome> {
    let columns = table.columns;
    let artificial = |j: usize| j >= real;
    if !table.basis.iter().any(|&b| artificial(b)) {
        return Ok(true);
    }
    let cost: Vec<Rational> = (0..columns)
        .map(|j| {
            if artificial(j) {
                Rational::ONE
            } else {
                Rational::ZERO
            }
        })
        .collect();
    table.minimize(&cost, &vec![true; columns], work)?;
    let sum = table
        .rows
        .iter()
        .zip(&table.basis)
        .filter(|&(_, &b)| artificial(b))
        .fold(Rational::ZERO, |s, (row, _)| &s + &row[columns]);
    if sum.is_positive() {
        return Ok(false);
    }
    let mut i = 0;
    while i < table.rows.len() {
        if artificial(table.basis[i]) {
            if let Some(j) = (0..real).find(|&j| !table.rows[i][j].is_zero()) {
                table.pivot(i, j);
            } else {
                table.rows.remove(i);
                table.basis.remove(i);
                continue;
            }
        }
        i += 1;
    }
    Ok(true)
}

/// Minimizes `c·x` over `x ≥ 0` and the rows, with the optimal tableau.
pub(super) fn minimize(
    c: &[Rational],
    rows: &[Row],
    work: &mut Work<'_>,
) -> (Outcome, Option<Final>) {
    let n = c.len();
    // At most one slack and one artificial per row.
    let slacks = rows
        .iter()
        .filter(|(_, r, _)| *r != Relation::Equal)
        .count();
    let cells = rows.len().saturating_mul(n + slacks + rows.len() + 1);
    if cells > work.cells {
        return (Outcome::GaveUp, None);
    }
    let (mut table, slack_rows) = tableau(n, rows);
    let real = n + slack_rows.len();
    match phase_one(&mut table, real, work) {
        Err(stop) => return (stop, None),
        Ok(false) => return (Outcome::Infeasible, None),
        Ok(true) => {}
    }
    let columns = table.columns;
    let mut cost = vec![Rational::ZERO; columns];
    cost[..n].clone_from_slice(c);
    let allowed: Vec<bool> = (0..columns).map(|j| j < real).collect();
    match table.minimize(&cost, &allowed, work) {
        Err(stop) => (stop, None),
        Ok(false) => (Outcome::Unbounded, None),
        Ok(true) => {
            let mut x = vec![Rational::ZERO; n];
            for (row, &b) in table.rows.iter().zip(&table.basis) {
                if b < n {
                    x[b] = row[columns].clone();
                }
            }
            let value = c
                .iter()
                .zip(&x)
                .fold(Rational::ZERO, |s, (a, v)| &s + &(a * v));
            let last = Final {
                rows: table.rows,
                basis: table.basis,
                real,
                slack_rows,
            };
            (Outcome::Optimal(x, value), Some(last))
        }
    }
}

#[cfg(test)]
mod tests;
