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
//!   keeps degenerate pivots from cycling; past [`MAX_PIVOTS`] the method
//!   gives up.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the simplex method: c·x, n variables, m rows, i and j"
)]

use cas_num::Rational;

/// The most pivots of one program.
const MAX_PIVOTS: usize = 100_000;

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
    /// Past [`MAX_PIVOTS`].
    GaveUp,
}

/// A row `a·x (relation) b`.
pub(super) type Row = (Vec<Rational>, Relation, Rational);

/// A tableau in canonical form: each row's basic variable has a unit
/// column; the last entry of a row is its right side.
struct Tableau {
    rows: Vec<Vec<Rational>>,
    basis: Vec<usize>,
    columns: usize,
    pivots: usize,
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
        self.pivots += 1;
    }

    /// The reduced cost of column `j` for `cost`.
    fn reduced(&self, cost: &[Rational], j: usize) -> Rational {
        let mut r = cost[j].clone();
        for (row, &b) in self.rows.iter().zip(&self.basis) {
            r = &r - &(&cost[b] * &row[j]);
        }
        r
    }

    /// Minimizes `cost` over the columns `allowed` to enter; `Some(false)`
    /// for an unbounded program, `None` past the pivots.
    fn minimize(&mut self, cost: &[Rational], allowed: &[bool]) -> Option<bool> {
        let rhs = self.columns;
        loop {
            if self.pivots > MAX_PIVOTS {
                return None;
            }
            let entering = (0..self.columns).find(|&j| {
                allowed[j] && !self.basis.contains(&j) && self.reduced(cost, j).is_negative()
            });
            let Some(j) = entering else {
                return Some(true);
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
                return Some(false);
            };
            self.pivot(i, j);
        }
    }
}

/// The first tableau of `n` variables and the rows, their right sides made
/// nonnegative: the number of slacks and surpluses too.
fn tableau(n: usize, rows: &[Row]) -> (Tableau, usize) {
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
        pivots: 0,
    };
    let (mut slack, mut artificial) = (n, n + slacks);
    for (a, rel, b) in &rows {
        let mut row = vec![Rational::ZERO; columns + 1];
        row[..n].clone_from_slice(a);
        row[columns] = b.clone();
        if *rel != Relation::Equal {
            row[slack] = if *rel == Relation::AtMost {
                Rational::ONE
            } else {
                -Rational::ONE
            };
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
    (table, slacks)
}

/// The first phase: the artificials' sum minimized, then those left in the
/// basis at 0 pivoted out, or their rows, redundant, dropped. `Some(false)`
/// for no feasible point, `None` past the pivots.
fn phase_one(table: &mut Tableau, real: usize) -> Option<bool> {
    let columns = table.columns;
    let artificial = |j: usize| j >= real;
    if !table.basis.iter().any(|&b| artificial(b)) {
        return Some(true);
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
    table.minimize(&cost, &vec![true; columns])?;
    let sum = table
        .rows
        .iter()
        .zip(&table.basis)
        .filter(|&(_, &b)| artificial(b))
        .fold(Rational::ZERO, |s, (row, _)| &s + &row[columns]);
    if sum.is_positive() {
        return Some(false);
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
    Some(true)
}

/// Minimizes `c·x` over `x ≥ 0` and the rows.
pub(super) fn minimize(c: &[Rational], rows: &[Row]) -> Outcome {
    let n = c.len();
    let (mut table, slacks) = tableau(n, rows);
    let real = n + slacks;
    match phase_one(&mut table, real) {
        None => return Outcome::GaveUp,
        Some(false) => return Outcome::Infeasible,
        Some(true) => {}
    }
    let columns = table.columns;
    let mut cost = vec![Rational::ZERO; columns];
    cost[..n].clone_from_slice(c);
    let allowed: Vec<bool> = (0..columns).map(|j| j < real).collect();
    match table.minimize(&cost, &allowed) {
        None => Outcome::GaveUp,
        Some(false) => Outcome::Unbounded,
        Some(true) => {
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
            Outcome::Optimal(x, value)
        }
    }
}

#[cfg(test)]
mod tests;
