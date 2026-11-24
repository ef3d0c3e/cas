//! Integer programs (CAS-DOC-38 §4, CAS-418, CAS-600): giac's
//! `assume = integer` or `lp_integer` (integers), `lp_nonnegint` or
//! `nonnegint` (nonnegative integers) and `lp_binary` (0 or 1), each
//! variable's own by `lp_integervariables` and `lp_binaryvariables` (the
//! others continuous), by branch and bound on the continuous relaxations.
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

use super::lpsolve::{Program, relaxation, tighten};
use super::simplex::Outcome;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most relaxations of one program.
const MAX_NODES: usize = 10_000;

/// A variable's ends, `None` for an infinite one.
type Bounds = Vec<(Option<Rational>, Option<Rational>)>;

/// Applies giac's assumption `assume = k` to every variable: `integer` and
/// `lp_integer` make it an integer, `lp_nonnegint` (`nonnegint`) a
/// nonnegative integer, `lp_binary` an integer in `0..1`, and
/// `lp_nonnegative` nonnegative, its bounds intersected with those
/// (lpsolve.cc:2672-2689); `false` for another value.
pub(super) fn kind(k: &Expr, bounds: &mut Bounds, integer: &mut [bool]) -> bool {
    let zero = Some(&Rational::ZERO);
    // `integer` is a name of [types], not a builtin.
    let (whole, low, high) = if *k == Expr::symbol("integer") {
        (true, None, None)
    } else {
        match k.as_symbol() {
            Some(Symbol::LP_INTEGER) => (true, None, None),
            Some(Symbol::LP_NON_NEG_INT) => (true, zero, None),
            Some(Symbol::LP_BINARY) => (true, zero, Some(&Rational::ONE)),
            Some(Symbol::LP_NON_NEGATIVE) => (false, zero, None),
            _ => return false,
        }
    };
    for (b, is) in bounds.iter_mut().zip(integer) {
        tighten(b, low, high);
        *is |= whole;
    }
    true
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
