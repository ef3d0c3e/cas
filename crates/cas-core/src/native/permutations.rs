//! Permutations (CAS-DOC-45 §2, CAS-463): giac's permutations, lists that
//! hold `0, …, n − 1` once each, and its cycles, lists of distinct
//! nonnegative integers, `[a, b, c]` taking `a` to `b`, `b` to `c` and `c`
//! to `a`. A cycle's permutation has the size `1 +` its largest entry, as
//! giac's `cycle2perm([0,2])` is `[2,1,0]`.
//!
//! - `PermutationQ`, giac's `is_permu`: `True` for a permutation.
//! - `PermutationCycles`, `permu2cycles`: the cycles of length 2 or more,
//!   each from its least entry, ordered by it; when `n − 1` is fixed, the
//!   cycle `[n − 1]` comes first to record the size, as giac's
//!   `permu2cycles([0,2,1,3])` is `[[3],[1,2]]`. `CyclesPermutation`,
//!   `cycles2permu`, gives the permutation back; `CyclePermutation`,
//!   `cycle2perm`, is one cycle's.
//! - `InversePermutation` and `InverseCycle`, `perminv` and `cycleinv`.
//! - `PermutationProduct[p, q]`, `p1op2`: `(p ∘ q)(i) = p(q(i))`, giac's
//!   `p1op2([1,2,0],[2,1,0])` being `[0,2,1]`; `PermutationCycleProduct`,
//!   `CyclePermutationProduct` and `CycleProduct`, `p1oc2`, `c1op2` and
//!   `c1oc2`, take a cycle's permutation in place of either. The shorter
//!   permutation fixes the points it lacks.
//! - `PermutationOrder`, `permuorder`: the lcm of the cycle lengths.
//!   `PermutationSignature`, `signature`: `(−1)^(n − c)` for `c` cycles,
//!   fixed points among them.
//! - `PermutationMatrix`, `permu2mat`: 1 at `(i, p(i))`.
//!
//! Arguments of another form leave the call; `PermutationQ` says `False`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of permutations: p and q, cycles c and d, the call e"
)]

use cas_num::Integer;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The entries of a list of nonnegative integers that fit an index.
fn naturals(e: &Expr) -> Option<Vec<usize>> {
    if e.head().and_then(Expr::as_symbol) != Some(Symbol::LIST) {
        return None;
    }
    e.args()
        .iter()
        .map(|x| x.as_i64().and_then(|k| usize::try_from(k).ok()))
        .collect()
}

/// Whether the entries are distinct.
fn distinct(v: &[usize]) -> bool {
    let mut seen = v.to_vec();
    seen.sort_unstable();
    seen.windows(2).all(|w| w[0] != w[1])
}

/// `e` as a permutation.
fn permutation(e: &Expr) -> Option<Vec<usize>> {
    let p = naturals(e)?;
    (p.iter().all(|&i| i < p.len()) && distinct(&p)).then_some(p)
}

/// `e` as a cycle: distinct entries, one at least.
fn cycle(e: &Expr) -> Option<Vec<usize>> {
    let c = naturals(e)?;
    (!c.is_empty() && distinct(&c)).then_some(c)
}

/// The permutation of disjoint cycles, of the size `1 +` their largest
/// entry; `None` for cycles that meet.
fn from_cycles(cycles: &[Vec<usize>]) -> Option<Vec<usize>> {
    let n = cycles.iter().flatten().max().map_or(0, |m| m + 1);
    let all: Vec<usize> = cycles.iter().flatten().copied().collect();
    if !distinct(&all) {
        return None;
    }
    let mut p: Vec<usize> = (0..n).collect();
    for c in cycles {
        for (i, &a) in c.iter().enumerate() {
            p[a] = c[(i + 1) % c.len()];
        }
    }
    Some(p)
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn indices(p: &[usize]) -> Expr {
    list(
        p.iter()
            .map(|&i| Expr::from(i64::try_from(i).expect("an index"))),
    )
}

/// The cycles of `p`, every one, each from its least entry, ordered by it.
fn all_cycles(p: &[usize]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; p.len()];
    let mut cycles = Vec::new();
    for start in 0..p.len() {
        if seen[start] {
            continue;
        }
        let mut c = Vec::new();
        let mut i = start;
        while !seen[i] {
            seen[i] = true;
            c.push(i);
            i = p[i];
        }
        cycles.push(c);
    }
    cycles
}

/// `(p ∘ q)(i) = p(q(i))`, the shorter fixing the points it lacks.
fn compose(p: &[usize], q: &[usize]) -> Vec<usize> {
    let n = p.len().max(q.len());
    let at = |r: &[usize], i: usize| r.get(i).copied().unwrap_or(i);
    (0..n).map(|i| at(p, at(q, i))).collect()
}

fn one_argument(e: &Expr) -> Option<&Expr> {
    match e.args() {
        [a] => Some(a),
        _ => None,
    }
}

fn answered(ev: &mut Eval<'_>, name: &'static str, e: &Expr, result: Expr) -> Native {
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// `PermutationQ[p]`.
pub(super) fn permutation_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e) else {
        return Ok(None);
    };
    Ok(Some(Expr::from(if permutation(p).is_some() {
        Symbol::TRUE
    } else {
        Symbol::FALSE
    })))
}

/// `PermutationCycles[p]`.
pub(super) fn permutation_cycles(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e).and_then(permutation) else {
        return Ok(None);
    };
    let n = p.len();
    let mut cycles: Vec<Vec<usize>> = all_cycles(&p).into_iter().filter(|c| c.len() > 1).collect();
    if n > 0 && p[n - 1] == n - 1 {
        cycles.insert(0, vec![n - 1]);
    }
    let result = list(cycles.iter().map(|c| indices(c)));
    answered(ev, "permu2cycles", e, result)
}

/// `CyclesPermutation[{c₁, c₂, …}]`.
pub(super) fn cycles_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(arg) = one_argument(e) else {
        return Ok(None);
    };
    if arg.head().and_then(Expr::as_symbol) != Some(Symbol::LIST) {
        return Ok(None);
    }
    let Some(cycles) = arg.args().iter().map(cycle).collect::<Option<Vec<_>>>() else {
        return Ok(None);
    };
    let Some(p) = from_cycles(&cycles) else {
        return Ok(None);
    };
    answered(ev, "cycles2permu", e, indices(&p))
}

/// `CyclePermutation[c]`.
pub(super) fn cycle_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(c) = one_argument(e).and_then(cycle) else {
        return Ok(None);
    };
    let p = from_cycles(&[c]).expect("one cycle");
    answered(ev, "cycle2perm", e, indices(&p))
}

/// `InversePermutation[p]`.
pub(super) fn inverse_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e).and_then(permutation) else {
        return Ok(None);
    };
    let mut inverse = vec![0; p.len()];
    for (i, &j) in p.iter().enumerate() {
        inverse[j] = i;
    }
    answered(ev, "perminv", e, indices(&inverse))
}

/// `InverseCycle[c]`: the cycle backwards.
pub(super) fn inverse_cycle(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(mut c) = one_argument(e).and_then(cycle) else {
        return Ok(None);
    };
    c.reverse();
    answered(ev, "cycleinv", e, indices(&c))
}

/// The permutation of an argument that is a permutation, or a cycle when
/// `is_cycle`.
fn operand(e: &Expr, is_cycle: bool) -> Option<Vec<usize>> {
    if is_cycle {
        from_cycles(&[cycle(e)?])
    } else {
        permutation(e)
    }
}

/// The products: `first` and `second` say which operands are cycles.
fn product(ev: &mut Eval<'_>, e: &Expr, first: bool, second: bool, name: &'static str) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let (Some(p), Some(q)) = (operand(a, first), operand(b, second)) else {
        return Ok(None);
    };
    answered(ev, name, e, indices(&compose(&p, &q)))
}

/// `PermutationProduct[p, q]`.
pub(super) fn permutation_product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    product(ev, e, false, false, "p1op2")
}

/// `PermutationCycleProduct[p, c]`.
pub(super) fn permutation_cycle_product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    product(ev, e, false, true, "p1oc2")
}

/// `CyclePermutationProduct[c, p]`.
pub(super) fn cycle_permutation_product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    product(ev, e, true, false, "c1op2")
}

/// `CycleProduct[c, d]`.
pub(super) fn cycle_product(ev: &mut Eval<'_>, e: &Expr) -> Native {
    product(ev, e, true, true, "c1oc2")
}

/// `PermutationOrder[p]`.
pub(super) fn permutation_order(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e).and_then(permutation) else {
        return Ok(None);
    };
    let order = all_cycles(&p).iter().fold(Integer::ONE, |l, c| {
        let k = Integer::from(i64::try_from(c.len()).expect("a length"));
        let g = l.gcd(&k);
        (&l * &k).div_exact(&g)
    });
    answered(ev, "permuorder", e, Expr::from(order))
}

/// `PermutationSignature[p]`.
pub(super) fn permutation_signature(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e).and_then(permutation) else {
        return Ok(None);
    };
    let even = (p.len() - all_cycles(&p).len()) % 2 == 0;
    answered(ev, "signature", e, Expr::from(if even { 1 } else { -1 }))
}

/// `PermutationMatrix[p]`.
pub(super) fn permutation_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = one_argument(e).and_then(permutation) else {
        return Ok(None);
    };
    if p.is_empty() {
        return Ok(None);
    }
    let rows = p
        .iter()
        .map(|&j| list((0..p.len()).map(|k| Expr::from(i64::from(k == j)))));
    answered(ev, "permu2mat", e, list(rows))
}

#[cfg(test)]
mod tests;
