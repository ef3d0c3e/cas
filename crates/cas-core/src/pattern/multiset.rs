//! The arguments of an `Orderless` head, matched as a multiset (CAS-DOC-4
//! §3.2–3.3, CAS-57).
//!
//! The subject's arguments are sorted, so equal ones are adjacent. The
//! argument nodes take subjects in the order of `Mode::order`: literals
//! first, then nodes with a fixed head, which consider only subjects with that
//! head, then blanks, then sequences. A node whose value is known, a literal
//! or a variable bound before, takes only a subject equal to it. Under `Flat`
//! a node may also take a group of several subjects, bound wrapped in the
//! head. The last node takes whatever is left, so one trailing sequence is
//! linear; several sequences enumerate the ways of splitting what is left,
//! which the step budget bounds.
//!
//! Equal subjects are interchangeable: among equal free subjects a choice
//! takes the first ones, so no two choices give the same bindings. The
//! choices of a node come in a fixed order, which defines the first match:
//! by size, smallest first, then taking as many as possible of the earliest
//! subjects.
//!
//! In `Times`, lists do not commute (R13): once every node has its subjects,
//! the lists must have been taken in the order of the nodes.

use super::{Matcher, Mode, Node, Rest, fixed_head, sequence};
use crate::context::EvalError;
use crate::eval::is_list_like;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The state of one multiset match.
struct Multiset<'a> {
    mode: &'a Mode,
    nodes: &'a [Node],
    subjects: &'a [Expr],
    /// The node that took each subject.
    taker: Vec<Option<usize>>,
    /// How many subjects no node has taken.
    free: usize,
}

impl Multiset<'_> {
    /// The free subjects that `keep` accepts, in runs of equal ones.
    fn runs(&self, keep: impl Fn(&Expr) -> bool) -> Vec<Vec<usize>> {
        let mut runs: Vec<Vec<usize>> = Vec::new();
        for (j, s) in self.subjects.iter().enumerate() {
            if self.taker[j].is_some() || !keep(s) {
                continue;
            }
            match runs.last_mut() {
                Some(run) if self.subjects[run[0]] == *s => run.push(j),
                _ => runs.push(vec![j]),
            }
        }
        runs
    }

    fn mark(&mut self, chosen: &[usize], node: usize) {
        for &j in chosen {
            self.taker[j] = Some(node);
        }
        self.free -= chosen.len();
    }

    fn unmark(&mut self, chosen: &[usize]) {
        for &j in chosen {
            self.taker[j] = None;
        }
        self.free += chosen.len();
    }

    /// `chosen` as what a node binds: the subject itself, or a group in the
    /// head.
    fn item(&self, chosen: &[usize]) -> Expr {
        if let [j] = chosen {
            return self.subjects[*j].clone();
        }
        let group: Vec<Expr> = chosen.iter().map(|&j| self.subjects[j].clone()).collect();
        Expr::apply(self.mode.head.as_ref().expect("a Flat head"), &group)
    }

    /// The first free subjects equal to `value`, or under `Flat` to the
    /// arguments of `value` when it is an application of the head: the only
    /// choice of a node whose value is known.
    fn find(&self, value: &Expr, group: bool) -> Option<Vec<usize>> {
        let parts = match value.view() {
            crate::expr::View::Apply { head, args }
                if group && Some(head) == self.mode.head.as_ref() && args.len() >= 2 =>
            {
                args
            }
            _ => std::slice::from_ref(value),
        };
        let mut chosen: Vec<usize> = Vec::with_capacity(parts.len());
        for part in parts {
            let j = (0..self.subjects.len()).find(|&j| {
                self.taker[j].is_none() && !chosen.contains(&j) && self.subjects[j] == *part
            })?;
            chosen.push(j);
        }
        chosen.sort_unstable();
        Some(chosen)
    }

    /// In `Times`, whether the lists were taken in the order of the nodes
    /// (R13).
    fn lists_in_order(&self) -> bool {
        let mut last = 0;
        for (s, taker) in self.subjects.iter().zip(&self.taker) {
            if is_list_like(s) {
                let t = taker.expect("every subject is taken");
                if t < last {
                    return false;
                }
                last = t;
            }
        }
        true
    }
}

/// The ways to choose `size` subjects from runs of equal ones, among equal
/// ones the first: as the counts taken from each run, most from the first
/// runs first. In terms of the subjects, the choices come in increasing
/// order of their indices.
struct Choices<'r> {
    runs: &'r [Vec<usize>],
    counts: Vec<usize>,
    started: bool,
    size: usize,
}

impl<'r> Choices<'r> {
    fn new(runs: &'r [Vec<usize>], size: usize) -> Choices<'r> {
        Choices {
            runs,
            counts: vec![0; runs.len()],
            started: false,
            size,
        }
    }

    /// Fills the runs from `from` on with `n` subjects, as many as possible
    /// in the first ones; `false` if they cannot hold `n`.
    fn fill(&mut self, from: usize, mut n: usize) -> bool {
        for j in from..self.runs.len() {
            self.counts[j] = self.runs[j].len().min(n);
            n -= self.counts[j];
        }
        n == 0
    }

    /// The next choice, as subject indices in increasing order.
    fn next(&mut self) -> Option<Vec<usize>> {
        if !self.started {
            self.started = true;
            if !self.fill(0, self.size) {
                return None;
            }
            return Some(self.indices());
        }
        // The last run that can give one subject to the runs after it.
        let (mut after, mut room) = (0, 0);
        for j in (0..self.runs.len()).rev() {
            if self.counts[j] > 0 && room > after {
                self.counts[j] -= 1;
                self.fill(j + 1, after + 1);
                return Some(self.indices());
            }
            after += self.counts[j];
            room += self.runs[j].len();
        }
        None
    }

    fn indices(&self) -> Vec<usize> {
        self.runs
            .iter()
            .zip(&self.counts)
            .flat_map(|(run, &c)| run[..c].iter().copied())
            .collect()
    }
}

impl Matcher<'_, '_, '_> {
    /// Matches the argument nodes of an `Orderless` head against the
    /// subject's arguments as a multiset, then the rest.
    pub(super) fn multiset(
        &mut self,
        mode: &Mode,
        nodes: &[Node],
        subjects: &[Expr],
        rest: Rest<'_>,
    ) -> Result<bool, EvalError> {
        let mut ms = Multiset {
            mode,
            nodes,
            subjects,
            taker: vec![None; subjects.len()],
            free: subjects.len(),
        };
        self.take(&mut ms, 0, rest)
    }

    /// Gives subjects to the node at position `k` of the order, then to
    /// those after it, then calls the rest.
    fn take(&mut self, ms: &mut Multiset<'_>, k: usize, rest: Rest<'_>) -> Result<bool, EvalError> {
        if ms.free < ms.mode.least[k] || ms.mode.most[k].is_some_and(|most| ms.free > most) {
            return Ok(false);
        }
        let Some(&i) = ms.mode.order.get(k) else {
            if ms.mode.times && !ms.lists_in_order() {
                return Ok(false);
            }
            return rest(self);
        };
        let node = &ms.nodes[i];
        let last = k + 1 == ms.mode.order.len();
        // What the nodes after this one leave it.
        let spare = ms.free - ms.mode.least[k + 1];
        if let Some((slot, min, head)) = sequence(node) {
            let runs = ms.runs(|s| head.is_none_or(|h| s.head() == Some(h)));
            let sizes = if last { ms.free..=ms.free } else { min..=spare };
            for size in sizes {
                let mut choices = Choices::new(&runs, size);
                while let Some(chosen) = choices.next() {
                    self.step()?;
                    if size < min {
                        break;
                    }
                    let run: Vec<Expr> = chosen.iter().map(|&j| ms.subjects[j].clone()).collect();
                    let mark = self.trail.len();
                    let bound = match slot {
                        Some(slot) => self.bind(slot, Expr::apply(&Symbol::SEQUENCE.into(), &run)),
                        None => true,
                    };
                    ms.mark(&chosen, i);
                    let found = bound && self.take(ms, k + 1, rest)?;
                    ms.unmark(&chosen);
                    if found {
                        return Ok(true);
                    }
                    self.undo(mark);
                }
            }
            return Ok(false);
        }
        let (inner, default) = match node {
            Node::Optional { inner, default } => (inner.as_ref(), default.as_ref()),
            _ => (node, None),
        };
        let group = ms.mode.groups[i];
        if let Some(value) = self.known(inner) {
            // Its value is known: only an equal subject, or group, will do.
            if let Some(chosen) = ms.find(&value, group) {
                if !last || chosen.len() == ms.free {
                    let item = ms.item(&chosen);
                    ms.mark(&chosen, i);
                    let found = self.node(inner, &item, &mut |m| m.take(ms, k + 1, rest))?;
                    ms.unmark(&chosen);
                    if found {
                        return Ok(true);
                    }
                }
            }
        } else {
            let most = if group { spare } else { spare.min(1) };
            let sizes = if last { ms.free..=ms.free } else { 1..=most };
            let fixed = fixed_head(inner);
            for size in sizes {
                if size == 0 || size > 1 && !group {
                    continue;
                }
                // A group has the head of the application.
                let runs = ms.runs(|s| match fixed {
                    Some(h) if size == 1 => s.head() == Some(h),
                    Some(h) => Some(h) == ms.mode.head.as_ref(),
                    None => true,
                });
                let mut choices = Choices::new(&runs, size);
                while let Some(chosen) = choices.next() {
                    let item = ms.item(&chosen);
                    ms.mark(&chosen, i);
                    let found = self.node(inner, &item, &mut |m| m.take(ms, k + 1, rest))?;
                    ms.unmark(&chosen);
                    if found {
                        return Ok(true);
                    }
                }
            }
        }
        match default {
            Some(default) => self.missing(inner, default, &mut |m| m.take(ms, k + 1, rest)),
            None => Ok(false),
        }
    }

    /// The value `node` must have, when it is known before matching: a
    /// literal, or a variable bound already.
    fn known(&self, node: &Node) -> Option<Expr> {
        match node {
            Node::Literal(e) => Some(e.clone()),
            Node::Bind { slot, .. } => self.bindings[*slot].clone(),
            Node::Condition { inner, .. } | Node::Test { inner, .. } => self.known(inner),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests;
