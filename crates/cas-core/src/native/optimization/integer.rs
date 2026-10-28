//! Integer programs (CAS-DOC-38 §4, CAS-418): giac's `assume = integer` or
//! `lp_integer` (integers), `nonnegint` (nonnegative integers) and
//! `lp_binary` (0 or 1), by branch and bound on the continuous
//! relaxations.
//!
//! - A relaxation's best point with an integer variable at a fraction `v`
//!   branches into `x ≤ ⌊v⌋` and `x ≥ ⌈v⌉`, depth first, the first such
//!   variable first.
//! - A relaxation no better than the best integer point found prunes its
//!   branch; an infeasible one ends it.
//! - Past [`MAX_NODES`] relaxations, or a first relaxation that is
//!   unbounded (the integer program is then unbounded or infeasible, which
//!   this does not tell apart), the method gives up and the call stays.

use cas_num::Rational;

use super::lpsolve::{Program, relaxation};
use super::simplex::Outcome;
use crate::expr::Expr;

/// The most relaxations of one program.
const MAX_NODES: usize = 10_000;

/// A variable's ends, `None` for an infinite one.
type Bounds = Vec<(Option<Rational>, Option<Rational>)>;

/// The variables that must be integers under giac's option `k`, their
/// bounds narrowed for `nonnegint` and `lp_binary`; `None` for another
/// option.
pub(super) fn kind(k: &Expr, bounds: &mut Bounds) -> Option<Vec<bool>> {
    let is = |name: &str| *k == Expr::symbol(name);
    if is("integer") || is("lp_integer") {
        return Some(vec![true; bounds.len()]);
    }
    let (low, high) = if is("nonnegint") {
        (Rational::ZERO, None)
    } else if is("lp_binary") {
        (Rational::ZERO, Some(Rational::ONE))
    } else {
        return None;
    };
    for (l, h) in bounds.iter_mut() {
        if l.as_ref().is_none_or(|l| *l < low) {
            *l = Some(low.clone());
        }
        if let Some(high) = &high
            && h.as_ref().is_none_or(|h| h > high)
        {
            *h = Some(high.clone());
        }
    }
    Some(vec![true; bounds.len()])
}

/// Whether `value` is better than `best`.
fn better(value: &Rational, best: &Rational, maximum: bool) -> bool {
    if maximum { value > best } else { value < best }
}

/// The best integer point of the program.
pub(super) fn branch_and_bound(p: &Program) -> Outcome {
    let mut best: Option<(Vec<Rational>, Rational)> = None;
    let mut stack: Vec<Bounds> = vec![p.bounds.clone()];
    let mut nodes = 0;
    while let Some(bounds) = stack.pop() {
        nodes += 1;
        if nodes > MAX_NODES {
            return Outcome::GaveUp;
        }
        let empty = bounds
            .iter()
            .any(|(l, h)| matches!((l, h), (Some(l), Some(h)) if l > h));
        if empty {
            continue;
        }
        let (x, value) = match relaxation(p, &bounds) {
            Outcome::Optimal(x, value) => (x, value),
            Outcome::Infeasible => continue,
            Outcome::Unbounded | Outcome::GaveUp => return Outcome::GaveUp,
        };
        if let Some((_, b)) = &best
            && !better(&value, b, p.maximum)
        {
            continue;
        }
        let fraction = x
            .iter()
            .zip(&p.integer)
            .position(|(v, &integer)| integer && !v.is_integer());
        let Some(j) = fraction else {
            best = Some((x, value));
            continue;
        };
        let (floor, ceil) = (Rational::from(x[j].floor()), Rational::from(x[j].ceil()));
        let mut down = bounds.clone();
        if down[j].1.as_ref().is_none_or(|h| *h > floor) {
            down[j].1 = Some(floor);
        }
        let mut up = bounds;
        if up[j].0.as_ref().is_none_or(|l| *l < ceil) {
            up[j].0 = Some(ceil);
        }
        // The lower branch first.
        stack.push(up);
        stack.push(down);
    }
    match best {
        Some((x, value)) => Outcome::Optimal(x, value),
        None => Outcome::Infeasible,
    }
}
