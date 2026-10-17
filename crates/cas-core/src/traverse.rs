//! Structural traversal: parts and positions, visitors, queries and
//! substitution.
//!
//! Nothing here evaluates or normalises: rebuilt applications are interned
//! with [`Expr::apply`] exactly as they are. Every traversal uses an explicit
//! stack, so any depth works in constant native stack.
//!
//! Expressions are DAGs: a subexpression shared `k` times occurs `k` times in
//! the tree. [`Expr::pre_order`], [`Expr::post_order`] and
//! [`Expr::positions`] visit the tree, every occurrence. The rest visit each
//! distinct subexpression once, memoised by node identity (`==` and `Hash`
//! are O(1)), so they take time linear in the number of distinct
//! subexpressions however much is shared.
//!
//! Positions follow the `Part` convention: index `0` is the head and `i ≥ 1`
//! the `i`-th argument.
//!
//! ```
//! use cas_core::build::{sin, symbols};
//! use cas_core::Expr;
//!
//! let [x, y] = symbols(["x", "y"]);
//! let e = sin(&x) + 2 * x.pow(2);
//! assert_eq!(e.to_string(), "2*x^2 + Sin[x]");
//! assert_eq!(e.part(&[1, 2, 1]), Some(&x));
//! assert_eq!(e.positions(|s| *s == x), [vec![1, 2, 1], vec![2, 1]]);
//! assert_eq!(e.substitute(&[(x.clone(), &y + 1)]).to_string(), "2*(1 + y)^2 + Sin[1 + y]");
//! assert!(e.free_of(&y));
//! ```

use std::ops::Index;

use foldhash::fast::FixedState;
use hashbrown::{HashMap, HashSet};

use crate::expr::Expr;
use crate::symbol::Symbol;

type Memo<'a> = HashMap<&'a Expr, Expr, FixedState>;

impl Expr {
    /// The part at `position`: `0` is the head and `i ≥ 1` the `i`-th
    /// argument, so `&[]` is `self` and `&[2, 1]` the first argument of the
    /// second argument. `None` if the position does not exist.
    #[must_use]
    pub fn part(&self, position: &[usize]) -> Option<&Expr> {
        position.iter().try_fold(self, |e, &i| e.child(i))
    }

    /// Head (`0`) or argument (`i ≥ 1`) of an application.
    fn child(&self, i: usize) -> Option<&Expr> {
        match i {
            0 => self.head(),
            _ => self.args().get(i - 1),
        }
    }

    /// `self` with the part at `position` replaced by `new`; `None` if the
    /// position does not exist. The rest is shared with `self`.
    #[must_use]
    pub fn replace_part(&self, position: &[usize], new: Expr) -> Option<Expr> {
        let mut path = Vec::with_capacity(position.len());
        let mut e = self;
        for &i in position {
            path.push(e);
            e = e.child(i)?;
        }
        // Rebuild from the replaced part up to the root.
        Some(
            path.iter()
                .zip(position)
                .rev()
                .fold(new, |part, (&parent, &i)| with_child(parent, i, part)),
        )
    }

    /// The positions of the parts for which `pred` holds, in pre-order (a
    /// part before its subparts, the head before the arguments). Visits
    /// every occurrence: on a heavily shared expression the result itself
    /// can be exponentially long.
    pub fn positions(&self, mut pred: impl FnMut(&Expr) -> bool) -> Vec<Vec<usize>> {
        let mut found = Vec::new();
        let mut position = Vec::new();
        // `(e, len, index)`: `e` is at a position of length `len` whose
        // last index is `index`; the rest is the current `position`.
        let mut stack = vec![(self, 0, 0)];
        while let Some((e, len, index)) = stack.pop() {
            if len > 0 {
                position.truncate(len - 1);
                position.push(index);
            }
            if pred(e) {
                found.push(position.clone());
            }
            if let Some(head) = e.head() {
                let len = position.len() + 1;
                stack.extend(
                    e.args()
                        .iter()
                        .enumerate()
                        .rev()
                        .map(|(i, a)| (a, len, i + 1)),
                );
                stack.push((head, len, 0));
            }
        }
        found
    }

    /// Every subexpression, `self` first, then the head and the arguments,
    /// each in pre-order. Every occurrence is visited.
    #[must_use]
    pub fn pre_order(&self) -> PreOrder<'_> {
        PreOrder { stack: vec![self] }
    }

    /// Every subexpression in post-order: the head, then the arguments, then
    /// the application. Every occurrence is visited.
    #[must_use]
    pub fn post_order(&self) -> PostOrder<'_> {
        PostOrder {
            stack: vec![(self, false)],
        }
    }

    /// Each distinct subexpression once, in post-order (subexpressions
    /// before the expressions that contain them).
    #[must_use]
    pub fn distinct(&self) -> Distinct<'_> {
        Distinct {
            stack: vec![(self, false)],
            seen: HashSet::with_hasher(FixedState::default()),
        }
    }

    /// The symbols in `self`, heads included, each once, in canonical order.
    ///
    /// There are no binding constructs yet, so every symbol is free; once
    /// there are (Rewriting), bound variables are left out.
    #[must_use]
    pub fn symbols(&self) -> Vec<Symbol> {
        let mut symbols: Vec<Symbol> = self.distinct().filter_map(Expr::as_symbol).collect();
        symbols.sort_by_cached_key(|&s| Expr::from(s));
        symbols
    }

    /// Whether `sub` occurs nowhere in `self` (as a head, an argument or
    /// `self` itself).
    #[must_use]
    pub fn free_of(&self, sub: &Expr) -> bool {
        !self.contains(sub)
    }

    /// Whether `sub` occurs in `self`. Skips any subexpression that cannot
    /// contain it by its metadata (symbols, depth, size), so a symbol that
    /// does not occur is usually rejected in O(1).
    #[must_use]
    pub fn contains(&self, sub: &Expr) -> bool {
        let mut seen = HashSet::with_hasher(FixedState::default());
        let mut stack = vec![self];
        while let Some(e) = stack.pop() {
            if e == sub {
                return true;
            }
            if !may_occur_in(sub, e) || !seen.insert(e) {
                continue;
            }
            if let Some(head) = e.head() {
                stack.extend(e.args().iter().rev());
                stack.push(head);
            }
        }
        false
    }

    /// `head[f(arg)…]`, the same expression if `f` changes no argument.
    /// Atoms are returned unchanged.
    #[must_use]
    pub fn map_args(&self, f: impl FnMut(&Expr) -> Expr) -> Expr {
        let Some(head) = self.head() else {
            return self.clone();
        };
        let args: Vec<Expr> = self.args().iter().map(f).collect();
        rebuild(self, head, &args)
    }

    /// Replaces subexpressions top-down: where `f` gives `Some(r)`, the
    /// subexpression becomes `r` and is not descended into; elsewhere the
    /// head and arguments are replaced in turn. Unchanged subexpressions are
    /// returned as the same handle.
    ///
    /// `f` is called at most once per distinct subexpression, and the
    /// result is shared wherever the subexpression was, so the time is
    /// linear in the number of distinct subexpressions.
    ///
    /// ```
    /// use cas_core::build::symbols;
    ///
    /// let [x, y] = symbols(["x", "y"]);
    /// let e = (&x + &y).pow(2);
    /// let swapped = e.replace(|s| if *s == x { Some(y.clone()) } else if *s == y { Some(x.clone()) } else { None });
    /// assert_eq!(swapped.to_string(), "(y + x)^2"); // not re-sorted
    /// assert_eq!(e.replace(|_| None), e);
    /// ```
    #[must_use]
    pub fn replace(&self, mut f: impl FnMut(&Expr) -> Option<Expr>) -> Expr {
        rewrite(self, |e| match f(e) {
            Some(r) => Step::Done(r),
            None => Step::Descend,
        })
    }

    /// Replaces each occurrence of a `from` by its `to`, simultaneously and
    /// top-down, without replacing inside a result: `ReplaceAll` with
    /// literal rules, structurally. The first matching rule wins.
    #[must_use]
    pub fn substitute(&self, rules: &[(Expr, Expr)]) -> Expr {
        let mut table: HashMap<&Expr, &Expr, FixedState> =
            HashMap::with_capacity_and_hasher(rules.len(), FixedState::default());
        for (from, to) in rules {
            table.entry(from).or_insert(to);
        }
        rewrite(self, |e| {
            if let Some(&to) = table.get(e) {
                Step::Done(to.clone())
            } else if rules.iter().any(|(from, _)| may_occur_in(from, e)) {
                Step::Descend
            } else {
                Step::Done(e.clone())
            }
        })
    }

    /// Applies `f` bottom-up: to every subexpression after its head and
    /// arguments have been mapped, `self` last. Unchanged subexpressions are
    /// returned as the same handle; `f` is called once per distinct
    /// subexpression.
    ///
    /// ```
    /// use cas_core::build::symbols;
    /// use cas_core::Symbol;
    ///
    /// let [x, y] = symbols(["x", "y"]);
    /// // b^2 → b*b, inner squares first.
    /// let e = (x.pow(2) + y.pow(2)).pow(2);
    /// let expanded = e.map_bottom_up(|s| match (s.head().and_then(|h| h.as_symbol()), s.args()) {
    ///     (Some(Symbol::POWER), [b, n]) if n.as_i64() == Some(2) => b * b,
    ///     _ => s.clone(),
    /// });
    /// assert_eq!(expanded.to_string(), "(x*x + y*y)*(x*x + y*y)");
    /// ```
    #[must_use]
    pub fn map_bottom_up(&self, f: impl FnMut(&Expr) -> Expr) -> Expr {
        bottom_up(self, f)
    }
}

/// See [`Expr::map_bottom_up`].
fn bottom_up(root: &Expr, mut f: impl FnMut(&Expr) -> Expr) -> Expr {
    let mut memo = Memo::with_hasher(FixedState::default());
    // Post-order: the parts of `e` are mapped before `e`.
    for e in root.distinct() {
        let rebuilt = match e.head() {
            Some(head) => {
                let head = memo[head].clone();
                let args: Vec<Expr> = e.args().iter().map(|a| memo[a].clone()).collect();
                rebuild(e, &head, &args)
            }
            None => e.clone(),
        };
        let mapped = f(&rebuilt);
        memo.insert(e, mapped);
    }
    memo.remove(root).expect("the root is mapped")
}

/// `parent` with its part `i` (`0` the head) replaced by `part`.
fn with_child(parent: &Expr, i: usize, part: Expr) -> Expr {
    let mut head = parent.head().expect("an application").clone();
    let mut args = parent.args().to_vec();
    if i == 0 {
        head = part;
    } else {
        args[i - 1] = part;
    }
    Expr::apply(&head, &args)
}

/// What [`rewrite`] does with a subexpression.
enum Step {
    /// Replace it by this.
    Done(Expr),
    /// Rewrite its head and arguments.
    Descend,
}

/// Top-down rewriting, memoised: `step` is called once per distinct
/// subexpression reached. The stack is LIFO, so a subexpression's whole
/// subtree is finished before any later occurrence of it is popped.
fn rewrite(root: &Expr, mut step: impl FnMut(&Expr) -> Step) -> Expr {
    let mut memo = Memo::with_hasher(FixedState::default());
    let mut stack = vec![(root, false)];
    while let Some((e, expanded)) = stack.pop() {
        if expanded {
            let head = memo[e.head().expect("an application")].clone();
            let args: Vec<Expr> = e.args().iter().map(|a| memo[a].clone()).collect();
            let result = rebuild(e, &head, &args);
            memo.insert(e, result);
            continue;
        }
        if memo.contains_key(e) {
            continue;
        }
        match step(e) {
            Step::Done(r) => {
                memo.insert(e, r);
            }
            Step::Descend => match e.head() {
                Some(head) => {
                    stack.push((e, true));
                    stack.extend(e.args().iter().rev().map(|a| (a, false)));
                    stack.push((head, false));
                }
                None => {
                    memo.insert(e, e.clone());
                }
            },
        }
    }
    memo.remove(root).expect("the root is rewritten")
}

/// `head[args…]`, or `e` itself if that is what it already is.
fn rebuild(e: &Expr, head: &Expr, args: &[Expr]) -> Expr {
    if Some(head) == e.head() && args == e.args() {
        e.clone()
    } else {
        Expr::apply(head, args)
    }
}

/// `false` if `sub` certainly does not occur in `e`: it would need all its
/// symbols, and no more depth or size than `e` has. Saturated counts stay
/// ordered, so the test stays sound.
fn may_occur_in(sub: &Expr, e: &Expr) -> bool {
    let (s, m) = (sub.meta(), e.meta());
    s.symbols & m.symbols == s.symbols && s.depth <= m.depth && s.size <= m.size
}

/// The parts of an application: `e[0]` is the head and `e[i]` the `i`-th
/// argument, as in [`Expr::part`].
///
/// # Panics
///
/// If `e` is an atom or has fewer than `i` arguments.
impl Index<usize> for Expr {
    type Output = Expr;

    fn index(&self, i: usize) -> &Expr {
        self.child(i)
            .unwrap_or_else(|| panic!("no part {i}: {} arguments", self.args().len()))
    }
}

/// See [`Expr::pre_order`].
#[derive(Debug)]
pub struct PreOrder<'a> {
    stack: Vec<&'a Expr>,
}

impl<'a> Iterator for PreOrder<'a> {
    type Item = &'a Expr;

    fn next(&mut self) -> Option<&'a Expr> {
        let e = self.stack.pop()?;
        if let Some(head) = e.head() {
            self.stack.extend(e.args().iter().rev());
            self.stack.push(head);
        }
        Some(e)
    }
}

/// See [`Expr::post_order`].
#[derive(Debug)]
pub struct PostOrder<'a> {
    /// `(e, expanded)`: `e` is yielded once its parts have been.
    stack: Vec<(&'a Expr, bool)>,
}

impl<'a> Iterator for PostOrder<'a> {
    type Item = &'a Expr;

    fn next(&mut self) -> Option<&'a Expr> {
        loop {
            let (e, expanded) = self.stack.pop()?;
            match e.head() {
                Some(head) if !expanded => {
                    self.stack.push((e, true));
                    self.stack.extend(e.args().iter().rev().map(|a| (a, false)));
                    self.stack.push((head, false));
                }
                _ => return Some(e),
            }
        }
    }
}

/// See [`Expr::distinct`].
#[derive(Debug)]
pub struct Distinct<'a> {
    stack: Vec<(&'a Expr, bool)>,
    seen: HashSet<&'a Expr, FixedState>,
}

impl<'a> Iterator for Distinct<'a> {
    type Item = &'a Expr;

    fn next(&mut self) -> Option<&'a Expr> {
        loop {
            let (e, expanded) = self.stack.pop()?;
            if expanded {
                return Some(e);
            }
            // LIFO: an occurrence popped later than another is popped only
            // after the first one's subtree, and itself, were yielded.
            if !self.seen.insert(e) {
                continue;
            }
            match e.head() {
                Some(head) => {
                    self.stack.push((e, true));
                    self.stack.extend(e.args().iter().rev().map(|a| (a, false)));
                    self.stack.push((head, false));
                }
                None => return Some(e),
            }
        }
    }
}

#[cfg(test)]
mod tests;
