//! Integer programs (CAS-DOC-38 §4, CAS-418, CAS-600, CAS-602, CAS-603):
//! giac's `assume = integer` or `lp_integer` (integers), `lp_nonnegint` or
//! `nonnegint` (nonnegative integers) and `lp_binary` (0 or 1), each
//! variable's own by `lp_integervariables` and `lp_binaryvariables` (the
//! others continuous), by giac's branch and bound on the continuous
//! relaxations (lpsolve.cc:1678-1890), minimizing (`−c` for a maximum,
//! :1626-1629), every bound, gap and point exact:
//!
//! - Each relaxation takes giac's Gomory mixed-integer cuts, up to
//!   `lp_maxcuts` (5) of them, re-solved from scratch
//!   ([`super::lpsolve::relaxation`]).
//! - The active nodes are the relaxations with an integer variable at a
//!   fraction. `lp_nodeselect` picks one (:1704-1749): `lp_depthfirst` the
//!   newest, `lp_breadthfirst` the oldest, `lp_bestlocalbound` (the
//!   default) the least bound and then the deepest, `lp_bestprojection` the
//!   least `bound + (incumbent − root bound)·infeasibility/root
//!   infeasibility`, the incumbent 0 while there is none (giac's), then the
//!   least bound, and `lp_hybrid` depth first until an incumbent, then the
//!   least bound. Where giac draws among equals with `_rand` (:1741-1745)
//!   ours takes the oldest node.
//! - `lp_varselect` picks the variable (:1750-1789): `lp_firstfractional`
//!   and `lp_lastfractional` by index, `lp_mostfractional` the fraction
//!   nearest 1/2, the lowest index among equals (giac's `_rand`, :831), and
//!   `lp_pseudocost` (the default) the best pseudocost score (:197-213,
//!   doubles as giac's), the most fractional while a fractional variable
//!   scores 0, that is until both its branches were taken.
//! - The node branches into `x ≥ ⌈v⌉` and `x ≤ ⌊v⌋`, the up branch first
//!   when diving, else down then up (:1799-1866). A child better than the
//!   incumbent becomes it when integer, else an active node, and then giac's
//!   rounding heuristic (:844-933) may find a better incumbent
//!   (`lp_heuristic`, on). A new incumbent prunes the active nodes no better
//!   than it.
//! - The search stops when the relative gap between the incumbent and the
//!   least active bound is at most `lp_gaptolerance` (0; :1790-1797), past
//!   `lp_nodelimit` child relaxations (:1696), past `lp_iterationlimit`
//!   simplex iterations (:1700, [`super::simplex::Work`]), or at a child
//!   deeper than `lp_depthlimit` (:1803-1809); a limit of 0 is none, and
//!   the limits are at most [`MAX_NODES`] and
//!   [`super::simplex::MAX_PIVOTS`]. Stopped so, the answer is the
//!   incumbent, or none, with a warning that names the option.
//! - Past [`MAX_NODES`] or the simplex's iterations with no limit of the
//!   user's, or a first relaxation that is unbounded (the integer program
//!   is then unbounded or infeasible, which this does not tell apart; a
//!   child's region is within it, so only the first can be), the method
//!   gives up and the call stays.

use cas_num::Rational;
use cas_num::ball::{Float, Round};

use super::lpsolve::{Program, relaxation, tighten};
use super::simplex::{Outcome, Relation, Work};
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most child relaxations of one program, and of a user's
/// `lp_nodelimit` (was the most relaxations, the first one too, also
/// `10_000`).
const MAX_NODES: usize = 10_000;

/// giac's `LP_SCORE_FACTOR` (lpsolve.h:35).
const SCORE_FACTOR: f64 = 0.1667;

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

/// giac's node selection, `lp_nodeselect`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum NodeSelect {
    DepthFirst,
    BreadthFirst,
    BestLocalBound,
    BestProjection,
    Hybrid,
}

/// giac's branching variable selection, `lp_varselect`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum VarSelect {
    FirstFractional,
    LastFractional,
    MostFractional,
    PseudoCost,
}

/// The settings of branch and bound, giac's defaults (`lp_settings`,
/// lpsolve.cc:238-256; the selections :1681, :1751); a limit of 0 is none.
#[derive(Clone, Debug)]
pub(super) struct Settings {
    pub(super) node_select: NodeSelect,
    pub(super) var_select: VarSelect,
    pub(super) max_cuts: usize,
    pub(super) heuristic: bool,
    pub(super) depth_limit: usize,
    pub(super) node_limit: usize,
    pub(super) iteration_limit: usize,
    pub(super) gap_tolerance: Rational,
}

impl Default for Settings {
    fn default() -> Settings {
        Settings {
            node_select: NodeSelect::BestLocalBound,
            var_select: VarSelect::PseudoCost,
            max_cuts: 5,
            heuristic: true,
            depth_limit: 0,
            node_limit: 0,
            iteration_limit: 0,
            gap_tolerance: Rational::ZERO,
        }
    }
}

/// The option that stopped a search before its end, giac's name.
pub(super) type Stop = Option<&'static str>;

/// The nearest `f64`.
fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// The fractional part `q − ⌊q⌋`, giac's `fracpart` (lpsolve.cc:908-910).
pub(super) fn fraction(q: &Rational) -> Rational {
    q - &Rational::from(q.floor())
}

/// An active node: its bounds and depth, its relaxation's point and value
/// (minimized), its integer variables at a fraction with that fraction,
/// and their infeasibility `Σ min(f, 1 − f)` (lpsolve.cc:818-837).
struct Node {
    bounds: Bounds,
    depth: usize,
    x: Vec<Rational>,
    value: Rational,
    fractional: Vec<(usize, Rational)>,
    infeasibility: Rational,
}

impl Node {
    fn new(p: &Program, bounds: Bounds, depth: usize, x: Vec<Rational>, value: &Rational) -> Node {
        let fractional: Vec<(usize, Rational)> = x
            .iter()
            .zip(&p.integer)
            .enumerate()
            .filter(|(_, (_, is))| **is)
            .map(|(i, (v, _))| (i, fraction(v)))
            .filter(|(_, f)| !f.is_zero())
            .collect();
        let infeasibility = fractional.iter().fold(Rational::ZERO, |s, (_, f)| {
            let g = &Rational::ONE - f;
            &s + if *f < g { f } else { &g }
        });
        Node {
            bounds,
            depth,
            x,
            value: minimized(p, value),
            fractional,
            infeasibility,
        }
    }

    /// The variable whose fraction is nearest 1/2, the lowest index among
    /// equals where giac draws (lpsolve.cc:826-837).
    fn most_fractional(&self) -> usize {
        let mut best: Option<(usize, Rational)> = None;
        for (i, f) in &self.fractional {
            let g = &Rational::ONE - f;
            let away = if *f < g { f.clone() } else { g };
            if best.as_ref().is_none_or(|(_, b)| away > *b) {
                best = Some((*i, away));
            }
        }
        best.expect("a fractional variable").0
    }

    /// The fraction of variable `j`, 0 for an integer.
    fn fraction_of(&self, j: usize) -> Rational {
        self.fractional
            .iter()
            .find(|(i, _)| *i == j)
            .map_or(Rational::ZERO, |(_, f)| f.clone())
    }
}

/// `value` minimized: negated for a maximum.
fn minimized(p: &Program, value: &Rational) -> Rational {
    if p.maximum {
        -value.clone()
    } else {
        value.clone()
    }
}

/// giac's pseudocosts of each variable, down (0) and up (1), with the
/// branches they average (lpsolve.cc:197-213).
struct Pseudocosts {
    cost: Vec<[f64; 2]>,
    branches: Vec<[u32; 2]>,
}

impl Pseudocosts {
    /// `update_pseudocost(delta, fr, dir)`.
    fn update(&mut self, j: usize, delta: f64, fr: f64, dir: usize) {
        if fr == 0.0 {
            return;
        }
        let sigma = self.cost[j][dir] * f64::from(self.branches[j][dir])
            + delta / if dir == 0 { fr } else { 1.0 - fr };
        self.branches[j][dir] += 1;
        self.cost[j][dir] = sigma / f64::from(self.branches[j][dir]);
    }

    /// `score(fr)`: 0 until both branches were taken.
    fn score(&self, j: usize, fr: f64) -> f64 {
        if self.branches[j][0] == 0 || self.branches[j][1] == 0 {
            return 0.0;
        }
        let (low, high) = (fr * self.cost[j][0], (1.0 - fr) * self.cost[j][1]);
        (1.0 - SCORE_FACTOR) * low.min(high) + SCORE_FACTOR * low.max(high)
    }

    /// giac's branching variable (lpsolve.cc:1750-1789): the best score
    /// under `lp_pseudocost`, else, or while a variable scores 0, by
    /// `select`.
    fn branching(&self, node: &Node, select: VarSelect) -> usize {
        let mut j = None;
        if select == VarSelect::PseudoCost {
            let mut best = 0.0;
            for (i, f) in &node.fractional {
                let score = self.score(*i, to_f64(f));
                // One variable scoring 0 falls back to the most fractional.
                if score == 0.0 {
                    j = None;
                    break;
                }
                if score > best {
                    (j, best) = (Some(*i), score);
                }
            }
        }
        j.unwrap_or_else(|| match select {
            VarSelect::FirstFractional => node.fractional[0].0,
            VarSelect::LastFractional => node.fractional[node.fractional.len() - 1].0,
            VarSelect::MostFractional | VarSelect::PseudoCost => node.most_fractional(),
        })
    }
}

/// giac's node selection (lpsolve.cc:1704-1749): the index of the active
/// node to branch on; `best` the incumbent's value, `root_value` and
/// `root_infeasibility` the first relaxation's.
fn select(
    select: NodeSelect,
    dive: bool,
    active: &[Node],
    best: Option<&Rational>,
    root_value: &Rational,
    root_infeasibility: &Rational,
) -> usize {
    if select == NodeSelect::BreadthFirst {
        0
    } else if dive {
        active.len() - 1
    } else if select == NodeSelect::BestProjection {
        // giac's quirk: the incumbent counts 0 while there is none.
        let best = best.map_or(Rational::ZERO, Clone::clone);
        let slope = &(&best - root_value) / root_infeasibility;
        // The least (projection, bound), the oldest among equals.
        let quality = |n: &Node| (&n.value + &(&slope * &n.infeasibility), n.value.clone());
        let mut k = 0;
        let mut least = quality(&active[0]);
        for (i, n) in active.iter().enumerate().skip(1) {
            let q = quality(n);
            if q < least {
                (k, least) = (i, q);
            }
        }
        k
    } else {
        // The least bound, then the deepest, the oldest among equals.
        let mut k = 0;
        for (i, n) in active.iter().enumerate().skip(1) {
            let m = &active[k];
            if n.value < m.value || n.value == m.value && n.depth > m.depth {
                k = i;
            }
        }
        k
    }
}

/// The relaxation of a node: `Ok(None)` when it has no point (giac's
/// empty range, lpsolve.cc:607-608, takes no iteration), `Err` with the
/// outcome that ends the search.
fn solve(
    p: &Program,
    bounds: Bounds,
    depth: usize,
    work: &mut Work<'_>,
) -> Result<Option<Node>, Outcome> {
    let empty = bounds
        .iter()
        .any(|(l, h)| matches!((l, h), (Some(l), Some(h)) if l > h));
    if empty {
        return Ok(None);
    }
    match relaxation(p, &bounds, work, p.settings.max_cuts) {
        Outcome::Optimal(x, value) => Ok(Some(Node::new(p, bounds, depth, x, &value))),
        Outcome::Infeasible => Ok(None),
        other => Err(other),
    }
}

/// The violated rows of the program at `x` and the sum of their squared
/// violations, giac's `violated_constraints` (lpsolve.cc:1105-1117, its
/// norm's square, exact); ours also counts the program's bounds, which
/// giac's leaves out, so that the heuristic never answers a point outside
/// them.
fn violation(p: &Program, x: &[Rational]) -> (usize, Rational) {
    let mut count = 0;
    let mut squares = Rational::ZERO;
    let mut add = |d: Rational| {
        count += 1;
        squares = &squares + &(&d * &d);
    };
    for (a, relation, b) in &p.rows {
        let d = &a
            .iter()
            .zip(x)
            .fold(Rational::ZERO, |s, (a, v)| &s + &(a * v))
            - b;
        let kept = match relation {
            Relation::AtMost => !d.is_positive(),
            Relation::AtLeast => !d.is_negative(),
            Relation::Equal => d.is_zero(),
        };
        if !kept {
            add(d);
        }
    }
    for (v, (low, high)) in x.iter().zip(&p.bounds) {
        if let Some(low) = low
            && v < low
        {
            add(low - v);
        } else if let Some(high) = high
            && v > high
        {
            add(v - high);
        }
    }
    (count, squares)
}

/// giac's rounding heuristic (lpsolve.cc:844-880): the node's fractional
/// variables rounded to the nearest integer (halves away from zero, giac's
/// `_round`), then while a row is violated the one flip, by one toward the
/// other neighbour, that lowers the violation most (count, then sum of
/// squares), never the last one flipped back; an integer point that keeps
/// every row and bound and its value (minimized), or `None`.
fn rounding(p: &Program, node: &Node) -> Option<(Vec<Rational>, Rational)> {
    let mut point = node.x.clone();
    // Each variable rounded, and whether it went up.
    let mut up: Vec<(usize, bool)> = Vec::with_capacity(node.fractional.len());
    for (i, _) in &node.fractional {
        point[*i] = Rational::from(node.x[*i].round());
        up.push((*i, point[*i] >= node.x[*i]));
    }
    let step = |up: bool| if up { -Rational::ONE } else { Rational::ONE };
    let mut current = violation(p, &point);
    let mut last = None;
    while current.0 > 0 {
        let mut best = current.clone();
        let mut flip = None;
        for (k, &(i, is_up)) in up.iter().enumerate() {
            if last == Some(i) {
                continue;
            }
            let mut trial = point.clone();
            trial[i] = &trial[i] + &step(is_up);
            let v = violation(p, &trial);
            if v < best {
                flip = Some(k);
                best = v;
            }
        }
        let Some(k) = flip else {
            break;
        };
        let (i, is_up) = up[k];
        point[i] = &point[i] + &step(is_up);
        up[k].1 = !is_up;
        last = Some(i);
        current = best;
    }
    if current.0 > 0 {
        return None;
    }
    let value = p
        .objective
        .iter()
        .zip(&point)
        .fold(p.constant.clone(), |s, (a, v)| &s + &(a * v));
    let value = minimized(p, &value);
    Some((point, value))
}

/// The best integer point of the program, giac's branch and bound
/// (lpsolve.cc:1671-1890), and the option that stopped it early, if one
/// did.
#[allow(
    clippy::too_many_lines,
    reason = "giac's loop (lpsolve.cc:1695-1890), one step after another"
)]
pub(super) fn branch_and_bound(p: &Program, work: &mut Work<'_>) -> (Outcome, Stop) {
    let s = &p.settings;
    // The incumbent: a point and its value, minimized.
    let mut incumbent: Option<(Vec<Rational>, Rational)> = None;
    let finish = |incumbent: Option<(Vec<Rational>, Rational)>, stop: Stop| match incumbent {
        Some((x, value)) => (Outcome::Optimal(x, minimized(p, &value)), stop),
        None => (Outcome::Infeasible, stop),
    };
    // A relaxation's outcome that ends the search.
    let end = |outcome: Outcome, incumbent: Option<(Vec<Rational>, Rational)>| match outcome {
        Outcome::Limit => finish(incumbent, Some("lp_iterationlimit")),
        Outcome::Interrupted => (Outcome::Interrupted, None),
        // A child's region is within the first's: only the first can be
        // unbounded, and the method then gives up.
        _ => (Outcome::GaveUp, None),
    };
    let root = match solve(p, p.bounds.clone(), 0, work) {
        Ok(Some(root)) => root,
        Ok(None) => return (Outcome::Infeasible, None),
        Err(outcome) => return end(outcome, None),
    };
    if root.fractional.is_empty() {
        return finish(Some((root.x, root.value)), None);
    }
    let node_limit = if s.node_limit > 0 {
        s.node_limit.min(MAX_NODES)
    } else {
        MAX_NODES
    };
    let mut pseudocosts = Pseudocosts {
        cost: vec![[0.0; 2]; p.xs.len()],
        branches: vec![[0; 2]; p.xs.len()],
    };
    let (root_value, root_infeasibility) = (root.value.clone(), root.infeasibility.clone());
    let mut active = vec![root];
    // giac's `subproblems_examined`: the child relaxations.
    let mut nodes = 0;
    while !active.is_empty() {
        if work.interrupted() {
            return (Outcome::Interrupted, None);
        }
        if nodes >= node_limit {
            if s.node_limit > 0 {
                return finish(incumbent, Some("lp_nodelimit"));
            }
            return (Outcome::GaveUp, None);
        }
        if s.iteration_limit > 0 && work.iterations > s.iteration_limit {
            return finish(incumbent, Some("lp_iterationlimit"));
        }
        let dive = s.node_select == NodeSelect::DepthFirst
            || s.node_select == NodeSelect::Hybrid && incumbent.is_none();
        let bound = active
            .iter()
            .map(|n| &n.value)
            .min()
            .expect("an active node")
            .clone();
        let best = incumbent.as_ref().map(|(_, v)| v);
        let k = select(
            s.node_select,
            dive,
            &active,
            best,
            &root_value,
            &root_infeasibility,
        );
        let node = &active[k];
        let j = pseudocosts.branching(node, s.var_select);
        if let Some((_, best)) = &incumbent {
            let gap = if best.is_zero() {
                -bound.clone()
            } else {
                &(best - &bound) / &best.abs()
            };
            if gap <= s.gap_tolerance {
                let stop = s.gap_tolerance.is_positive().then_some("lp_gaptolerance");
                return finish(incumbent, stop);
            }
        }
        let fr = to_f64(&node.fraction_of(j));
        let mut updated = false;
        let mut children = Vec::new();
        for i in 0..2 {
            let depth = node.depth + 1;
            if s.depth_limit > 0 && depth > s.depth_limit {
                return finish(incumbent, Some("lp_depthlimit"));
            }
            let mut bounds = node.bounds.clone();
            if i == usize::from(!dive) {
                let ceil = Rational::from(node.x[j].ceil());
                tighten(&mut bounds[j], Some(&ceil), None);
            } else {
                let floor = Rational::from(node.x[j].floor());
                tighten(&mut bounds[j], None, Some(&floor));
            }
            nodes += 1;
            let child = match solve(p, bounds, depth, work) {
                Ok(Some(child)) => child,
                Ok(None) => continue,
                Err(outcome) => return end(outcome, incumbent),
            };
            // giac's direction: down (0) and up (1) when diving, else
            // always 1, its quirk.
            let dir = if dive { 1 - i } else { 1 };
            let delta = (to_f64(&child.value) - to_f64(&node.value)).abs();
            pseudocosts.update(j, delta, fr, dir);
            if incumbent.as_ref().is_some_and(|(_, v)| *v <= child.value) {
                continue;
            }
            if child.fractional.is_empty() {
                incumbent = Some((child.x, child.value));
                updated = true;
                continue;
            }
            if s.heuristic
                && let Some((x, value)) = rounding(p, &child)
                && incumbent.as_ref().is_none_or(|(_, v)| value < *v)
            {
                incumbent = Some((x, value));
                updated = true;
            }
            children.push(child);
        }
        active.remove(k);
        active.extend(children);
        if let Some((_, best)) = &incumbent
            && updated
        {
            active.retain(|n| n.value < *best);
        }
    }
    finish(incumbent, None)
}
