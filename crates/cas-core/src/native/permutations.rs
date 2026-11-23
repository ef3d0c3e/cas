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
//! - `NextPermutation` and `PreviousPermutation`, `nextperm` and `prevperm`
//!   (CAS-561, permu.cc): the next and previous permutation in
//!   lexicographic order (`std::next_permutation`), `undef` after the last
//!   and before the first.
//! - `PermutationGroup[p, q]`, `groupermu`, also of the list `[p, q]`: the
//!   group `p` and `q` generate, giac's order (permu.cc `groupermu`): `p`
//!   and `q`, the shorter fixing the points it lacks, then breadth first
//!   `p ∘ s` and `q ∘ s` for each `s` found, new ones appended; `p = q`
//!   stays twice. More than lists.rs's `MAX_ELEMENTS` entries in all leave
//!   the call (C0019).
//!
//! Arguments of another form leave the call; `PermutationQ` says `False`.
//! For the three of CAS-561, giac's error, Bad Argument Type for what is
//! not a list (of two lists for `groupermu`) and Bad Argument Value for a
//! list that is not a permutation, is a warning.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of permutations: p and q, cycles c and d, the call e"
)]

use cas_num::Integer;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::lists::MAX_ELEMENTS;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

const BAD_VALUE: &str = "bad argument value";
const BAD_TYPE: &str = "bad argument type";
const SIZE_LIMIT: &str = "calculation size limit exceeded";

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

/// giac's error for `e`: the call stays, with `why` as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// `std::next_permutation`, or `prev_permutation` when `!next`: whether
/// there was one.
fn step(p: &mut [usize], next: bool) -> bool {
    let after = |a: usize, b: usize| if next { a < b } else { a > b };
    let Some(i) = (1..p.len()).rev().find(|&i| after(p[i - 1], p[i])) else {
        p.reverse();
        return false;
    };
    let j = (i..p.len())
        .rev()
        .find(|&j| after(p[i - 1], p[j]))
        .expect("p[i] follows");
    p.swap(i - 1, j);
    p[i..].reverse();
    true
}

/// `NextPermutation[p]` and `PreviousPermutation[p]`.
fn neighbour(ev: &mut Eval<'_>, e: &Expr, next: bool, name: &'static str) -> Native {
    let Some(arg) = one_argument(e) else {
        return refuse(ev, e, BAD_TYPE);
    };
    if arg.head().and_then(Expr::as_symbol) != Some(Symbol::LIST) {
        return refuse(ev, e, BAD_TYPE);
    }
    let Some(mut p) = permutation(arg) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let result = if step(&mut p, next) {
        indices(&p)
    } else {
        Expr::from(Symbol::INDETERMINATE)
    };
    answered(ev, name, e, result)
}

/// `NextPermutation[p]`.
pub(super) fn next_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    neighbour(ev, e, true, "nextperm")
}

/// `PreviousPermutation[p]`.
pub(super) fn previous_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    neighbour(ev, e, false, "prevperm")
}

/// `PermutationGroup[p, q]` or `PermutationGroup[{p, q}]`.
pub(super) fn permutation_group(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let is_list = |x: &Expr| x.head().and_then(Expr::as_symbol) == Some(Symbol::LIST);
    let (p, q) = match e.args() {
        [p, q] => (p, q),
        [pq] if is_list(pq) && pq.args().len() == 2 => (&pq.args()[0], &pq.args()[1]),
        _ => return refuse(ev, e, BAD_TYPE),
    };
    if !is_list(p) || !is_list(q) {
        return refuse(ev, e, BAD_TYPE);
    }
    let (Some(mut a), Some(mut b)) = (permutation(p), permutation(q)) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let n = a.len().max(b.len());
    a.extend(a.len()..n);
    b.extend(b.len()..n);
    let fits = |count: usize| {
        count
            .checked_mul(n)
            .and_then(|total| i64::try_from(total).ok())
            .is_some_and(|total| total <= MAX_ELEMENTS)
    };
    if !fits(2) {
        return refuse(ev, e, SIZE_LIMIT);
    }
    // giac's `est_dans` looks through every element found so far.
    let mut found: std::collections::HashSet<Vec<usize>> = [a.clone(), b.clone()].into();
    let mut s = vec![a.clone(), b.clone()];
    let mut first = 0;
    while first < s.len() {
        let end = s.len();
        for j in first..end {
            for g in [&a, &b] {
                let r = compose(g, &s[j]);
                if !found.contains(&r) {
                    if !fits(s.len() + 1) {
                        return refuse(ev, e, SIZE_LIMIT);
                    }
                    found.insert(r.clone());
                    s.push(r);
                }
            }
        }
        first = end;
    }
    answered(ev, "groupermu", e, list(s.iter().map(|r| indices(r))))
}

#[cfg(test)]
mod tests;
