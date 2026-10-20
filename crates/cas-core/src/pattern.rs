//! Patterns and rules (CAS-DOC-4 §3, CAS-55): compiled once, matched by
//! backtracking.
//!
//! A pattern is an ordinary expression built from `Blank`, `BlankSequence`,
//! `BlankNullSequence`, `Pattern`, `Alternatives`, `Optional`, `Condition`,
//! `PatternTest` and `HoldPattern`; anything else matches itself. [`Pattern`]
//! compiles one into a tree of match nodes, numbers its variables (slots), and
//! records the literal symbols every match contains, as a Bloom filter in the
//! format of the node metadata, so a subject without them is rejected at once.
//!
//! The matcher backtracks through continuations: each node, once it matches,
//! calls the rest of the match, and a failure undoes the bindings made since.
//! So every way of matching is tried, in order, until one completes: shorter
//! sequences first, alternatives from the first, and a repeated variable must
//! bind the same expression, which is an O(1) check. Matching is syntactic:
//! `Flat` and `Orderless` heads are matched like the others until CAS-57.
//!
//! [`substitute`] instantiates a right-hand side: it replaces the variables,
//! but never inside a `Function` that binds them, so `f[a_] :> x -> a x`
//! keeps its parameter `x` (CAS-DOC-4 §3.1).

use std::fmt;

use crate::context::{Context, EvalError};
use crate::eval::{Eval, run};
use crate::expr::{Expr, View, symbol_bit};
use crate::symbol::Symbol;

/// Why an expression is not a pattern.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatternError {
    /// The construct at fault.
    pub part: Expr,
    /// Why.
    pub kind: PatternErrorKind,
}

/// What a [`PatternError`] is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PatternErrorKind {
    /// `Pattern[x, p]` whose `x` is not a symbol.
    NotAName,
    /// A construct with another number of arguments: `Blank[a, b]`,
    /// `Pattern[x]`, `Condition[p]`.
    Arity,
}

impl fmt::Display for PatternError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind {
            PatternErrorKind::NotAName => {
                write!(f, "a pattern variable must be a symbol: {:?}", self.part)
            }
            PatternErrorKind::Arity => {
                write!(
                    f,
                    "a pattern construct with the wrong arguments: {:?}",
                    self.part
                )
            }
        }
    }
}

impl std::error::Error for PatternError {}

/// A compiled pattern (CAS-DOC-4 §3.1).
///
/// ```
/// use cas_core::{Context, Expr, Pattern};
///
/// let p = Pattern::new(&Expr::parse_full_form("f[Pattern[x, Blank[]], Pattern[x, Blank[]]]").unwrap()).unwrap();
/// let ctx = Context::new();
/// assert!(ctx.matches(&Expr::parse_full_form("f[a, a]").unwrap(), &p).unwrap());
/// assert!(!ctx.matches(&Expr::parse_full_form("f[a, b]").unwrap(), &p).unwrap());
/// ```
#[derive(Clone, Debug)]
pub struct Pattern {
    root: Node,
    /// The variables, by slot.
    variables: Vec<Symbol>,
    /// The bits of the literal symbols every match contains.
    required: u64,
    /// The heads of the applications it mentions: matching depends on their
    /// attributes (CAS-56 recompiles when they change).
    heads: Vec<Symbol>,
}

#[derive(Clone, Debug)]
enum Node {
    /// Only itself.
    Literal(Expr),
    /// `_` or `_h`: one expression, with head `h`.
    Blank(Option<Expr>),
    /// `__` (at least 1) or `___` (at least 0), with head `h`: a run of
    /// arguments. Elsewhere, one expression.
    Sequence {
        min: usize,
        head: Option<Expr>,
    },
    /// `Pattern[x, p]`: `p`, bound to the slot of `x`.
    Bind {
        slot: usize,
        inner: Box<Node>,
    },
    /// `h[p1, …]`, its head matched by a node too.
    Apply {
        head: Box<Node>,
        args: Vec<Node>,
    },
    Alternatives(Vec<Node>),
    /// `Optional[p, d]`, as an argument: when it is missing, `p`'s variable
    /// takes the default, if there is one.
    Optional {
        inner: Box<Node>,
        default: Option<Expr>,
    },
    /// `Condition[p, test]`: `test`, with the variables bound, must evaluate
    /// to `True`.
    Condition {
        inner: Box<Node>,
        test: Expr,
    },
    /// `PatternTest[p, f]`: `f[subject]` must evaluate to `True`.
    Test {
        inner: Box<Node>,
        test: Expr,
    },
}

impl Pattern {
    /// Compiles `e`.
    ///
    /// # Errors
    ///
    /// A [`PatternError`] for a malformed construct.
    pub fn new(e: &Expr) -> Result<Pattern, PatternError> {
        Pattern::with_defaults(e, &|_| None)
    }

    /// Compiles `e`, where `x_.` as an argument of a user head `f` takes
    /// `defaults(f)` when it is missing (`Context::set_default`, CAS-56).
    pub(crate) fn with_defaults(
        e: &Expr,
        defaults: &dyn Fn(Symbol) -> Option<Expr>,
    ) -> Result<Pattern, PatternError> {
        let mut compiler = Compiler {
            variables: Vec::new(),
            heads: Vec::new(),
            defaults,
        };
        let root = compiler.node(e, None)?;
        let mut required = 0;
        root.required(&mut required);
        Ok(Pattern {
            root,
            variables: compiler.variables,
            required,
            heads: compiler.heads,
        })
    }

    /// The variables, in the order of their slots.
    #[must_use]
    pub fn variables(&self) -> &[Symbol] {
        &self.variables
    }

    /// The heads whose attributes matching depends on.
    #[must_use]
    pub fn heads(&self) -> &[Symbol] {
        &self.heads
    }

    /// Matches `subject`: `accept` receives each complete set of bindings,
    /// in order, until it accepts one. Returns whether one was accepted.
    pub(crate) fn match_with(
        &self,
        ev: &mut Eval<'_>,
        subject: &Expr,
        accept: Accept<'_>,
    ) -> Result<bool, EvalError> {
        if (subject.meta().symbols & self.required) != self.required {
            return Ok(false);
        }
        let mut m = Matcher {
            ev,
            variables: &self.variables,
            bindings: vec![None; self.variables.len()],
            trail: Vec::new(),
        };
        m.node(&self.root, subject, &mut |m| {
            // A variable left unbound, in an alternative not taken, is the
            // empty sequence.
            let bound: Vec<Expr> = m
                .bindings
                .iter()
                .map(|b| {
                    b.clone()
                        .unwrap_or_else(|| Expr::apply(&Symbol::SEQUENCE.into(), &[]))
                })
                .collect();
            accept(m.ev, &bound)
        })
    }

    /// Whether `subject` matches.
    pub(crate) fn matches(&self, ev: &mut Eval<'_>, subject: &Expr) -> Result<bool, EvalError> {
        self.match_with(ev, subject, &mut |_, _| Ok(true))
    }

    /// What a subject must be for this pattern to match it, read from the
    /// compiled nodes so that it agrees with the matcher.
    pub(crate) fn filter(&self) -> Filter {
        let mut root = &self.root;
        while let Node::Bind { inner, .. }
        | Node::Condition { inner, .. }
        | Node::Test { inner, .. } = root
        {
            root = inner;
        }
        let Node::Apply { head, args } = root else {
            return Filter::default();
        };
        let (mut min, mut max) = (0, Some(0));
        for a in args {
            let (least, most) = match (a, sequence(a)) {
                (Node::Optional { .. }, _) => (0, Some(1)),
                (_, Some((_, least, _))) => (least, None),
                (_, None) => (1, Some(1)),
            };
            min += least;
            max = max.zip(most).map(|(m, n)| m + n);
        }
        Filter {
            head: match head.as_ref() {
                Node::Literal(h) => Some(h.clone()),
                _ => None,
            },
            min,
            max,
            first: args.first().map_or(First::Any, first_argument),
        }
    }
}

/// What a subject must be for a pattern to match it (CAS-DOC-4 §5.1): its
/// head, its number of arguments, and its first argument, or that
/// argument's head, when the pattern fixes them. Checked before matching,
/// it skips most rules of a head in O(1).
///
/// Matching is syntactic: CAS-57 widens the filter for `Flat` and
/// `Orderless` heads, whose arguments may group or move.
#[derive(Clone, Debug, Default)]
pub(crate) struct Filter {
    head: Option<Expr>,
    min: usize,
    /// `None`: any number from `min`.
    max: Option<usize>,
    first: First,
}

#[derive(Clone, Debug, Default)]
enum First {
    #[default]
    Any,
    /// The first argument has this head.
    Head(Expr),
    /// The first argument is this.
    Is(Expr),
}

impl Filter {
    /// Whether a match of `e` is possible.
    pub(crate) fn admits(&self, e: &Expr) -> bool {
        let args = e.args();
        self.head.as_ref().is_none_or(|h| e.head() == Some(h))
            && args.len() >= self.min
            && self.max.is_none_or(|max| args.len() <= max)
            && match &self.first {
                First::Any => true,
                First::Head(h) => args.first().and_then(Expr::head) == Some(h),
                First::Is(x) => args.first() == Some(x),
            }
    }
}

/// What the first argument of a subject must be for `node`, the first
/// argument node, to match it. Only nodes that always take exactly the
/// first argument fix anything.
fn first_argument(node: &Node) -> First {
    match node {
        Node::Literal(e) => First::Is(e.clone()),
        Node::Blank(Some(h))
        | Node::Sequence {
            min: 1..,
            head: Some(h),
        } => First::Head(h.clone()),
        Node::Bind { inner, .. } | Node::Condition { inner, .. } | Node::Test { inner, .. } => {
            first_argument(inner)
        }
        Node::Apply { head, .. } => match head.as_ref() {
            Node::Literal(h) => First::Head(h.clone()),
            _ => First::Any,
        },
        Node::Blank(None)
        | Node::Sequence { .. }
        | Node::Alternatives(_)
        | Node::Optional { .. } => First::Any,
    }
}

/// What receives a complete set of bindings, by slot, and says whether it
/// accepts them.
type Accept<'a> = &'a mut dyn FnMut(&mut Eval<'_>, &[Expr]) -> Result<bool, EvalError>;

struct Compiler<'d> {
    variables: Vec<Symbol>,
    heads: Vec<Symbol>,
    /// The defaults of user heads.
    defaults: &'d dyn Fn(Symbol) -> Option<Expr>,
}

impl Compiler<'_> {
    fn slot(&mut self, x: Symbol) -> usize {
        if let Some(i) = self.variables.iter().position(|&v| v == x) {
            return i;
        }
        self.variables.push(x);
        self.variables.len() - 1
    }

    /// The node for `e`, an argument at `position` of `head` when it is one:
    /// that decides the default of `Optional[p]`.
    fn node(&mut self, e: &Expr, position: Option<(Symbol, usize)>) -> Result<Node, PatternError> {
        if !e.has_pattern() {
            return Ok(Node::Literal(e.clone()));
        }
        let View::Apply { head, args } = e.view() else {
            return Ok(Node::Literal(e.clone()));
        };
        let arity = || PatternError {
            part: e.clone(),
            kind: PatternErrorKind::Arity,
        };
        let head_of = |args: &[Expr]| -> Result<Option<Expr>, PatternError> {
            match args {
                [] => Ok(None),
                [h] => Ok(Some(h.clone())),
                _ => Err(arity()),
            }
        };
        Ok(match head.as_symbol() {
            Some(Symbol::BLANK) => Node::Blank(head_of(args)?),
            Some(Symbol::BLANK_SEQUENCE) => Node::Sequence {
                min: 1,
                head: head_of(args)?,
            },
            Some(Symbol::BLANK_NULL_SEQUENCE) => Node::Sequence {
                min: 0,
                head: head_of(args)?,
            },
            Some(Symbol::PATTERN) => {
                let [x, p] = args else {
                    return Err(arity());
                };
                let Some(x) = x.as_symbol() else {
                    return Err(PatternError {
                        part: e.clone(),
                        kind: PatternErrorKind::NotAName,
                    });
                };
                let slot = self.slot(x);
                Node::Bind {
                    slot,
                    inner: Box::new(self.node(p, position)?),
                }
            }
            Some(Symbol::ALTERNATIVES) => Node::Alternatives(
                args.iter()
                    .map(|a| self.node(a, position))
                    .collect::<Result<_, _>>()?,
            ),
            Some(Symbol::OPTIONAL) => match args {
                [p] => Node::Optional {
                    inner: Box::new(self.node(p, position)?),
                    default: position
                        .and_then(|(h, i)| default(h, i).or_else(|| (self.defaults)(h))),
                },
                [p, d] => Node::Optional {
                    inner: Box::new(self.node(p, position)?),
                    default: Some(d.clone()),
                },
                _ => return Err(arity()),
            },
            Some(Symbol::CONDITION) => {
                let [p, test] = args else {
                    return Err(arity());
                };
                Node::Condition {
                    inner: Box::new(self.node(p, position)?),
                    test: test.clone(),
                }
            }
            Some(Symbol::PATTERN_TEST) => {
                let [p, test] = args else {
                    return Err(arity());
                };
                Node::Test {
                    inner: Box::new(self.node(p, position)?),
                    test: test.clone(),
                }
            }
            Some(Symbol::HOLD_PATTERN) => {
                let [p] = args else {
                    return Err(arity());
                };
                self.node(p, position)?
            }
            symbol => {
                if let Some(s) = symbol.filter(|s| !s.is_builtin() && !self.heads.contains(s)) {
                    self.heads.push(s);
                }
                let mut nodes = Vec::with_capacity(args.len());
                for (i, a) in args.iter().enumerate() {
                    nodes.push(self.node(a, symbol.map(|s| (s, i + 1)))?);
                }
                Node::Apply {
                    head: Box::new(self.node(head, None)?),
                    args: nodes,
                }
            }
        })
    }
}

/// The default of `Optional[p]` as argument `i` of a builtin `head`: 0 in
/// sums, 1 in products and as an exponent (CAS-DOC-4 §3).
fn default(head: Symbol, i: usize) -> Option<Expr> {
    match (head, i) {
        (Symbol::PLUS, _) => Some(Expr::from(0)),
        (Symbol::TIMES, _) | (Symbol::POWER, 2) => Some(Expr::from(1)),
        _ => None,
    }
}

impl Node {
    /// Adds the bits of the literal symbols every match of this node
    /// contains: not those under alternatives or optional parts.
    fn required(&self, bits: &mut u64) {
        match self {
            Node::Literal(e) => *bits |= e.meta().symbols,
            Node::Blank(Some(h))
            | Node::Sequence {
                head: Some(h),
                min: 1..,
            } => {
                if let Some(s) = h.as_symbol() {
                    *bits |= symbol_bit(s);
                }
            }
            Node::Bind { inner, .. } | Node::Condition { inner, .. } | Node::Test { inner, .. } => {
                inner.required(bits);
            }
            Node::Apply { head, args } => {
                head.required(bits);
                for a in args {
                    a.required(bits);
                }
            }
            Node::Blank(_)
            | Node::Sequence { .. }
            | Node::Alternatives(_)
            | Node::Optional { .. } => {}
        }
    }
}

/// What the rest of a match does once a node has matched.
type Rest<'r> = &'r mut dyn FnMut(&mut Matcher<'_, '_, '_>) -> Result<bool, EvalError>;

struct Matcher<'e, 'c, 'p> {
    ev: &'e mut Eval<'c>,
    variables: &'p [Symbol],
    bindings: Vec<Option<Expr>>,
    /// The slots bound so far, in order, to undo on backtracking.
    trail: Vec<usize>,
}

impl Matcher<'_, '_, '_> {
    fn undo(&mut self, mark: usize) {
        for slot in self.trail.drain(mark..) {
            self.bindings[slot] = None;
        }
    }

    /// Binds `slot` to `value`, or checks it is bound to it already.
    fn bind(&mut self, slot: usize, value: Expr) -> bool {
        if let Some(bound) = &self.bindings[slot] {
            return *bound == value;
        }
        self.bindings[slot] = Some(value);
        self.trail.push(slot);
        true
    }

    /// The bindings so far, as substitutions.
    fn substitutions(&self) -> Vec<(Symbol, Expr)> {
        self.variables
            .iter()
            .zip(&self.bindings)
            .filter_map(|(&x, b)| b.clone().map(|b| (x, b)))
            .collect()
    }

    fn holds(&mut self, test: &Expr) -> Result<bool, EvalError> {
        let value = self.ev.eval(test)?;
        Ok(value.as_symbol() == Some(Symbol::TRUE))
    }

    /// Matches `node` against `subject`, then the rest. Each node takes a
    /// level of the recursion limit while it is on the stack: the rest of
    /// the match runs inside it.
    fn node(&mut self, node: &Node, subject: &Expr, rest: Rest<'_>) -> Result<bool, EvalError> {
        self.ev.enter()?;
        let result = self.node_inner(node, subject, rest);
        self.ev.leave();
        result
    }

    fn node_inner(
        &mut self,
        node: &Node,
        subject: &Expr,
        rest: Rest<'_>,
    ) -> Result<bool, EvalError> {
        match node {
            Node::Literal(e) => {
                if e == subject {
                    rest(self)
                } else {
                    Ok(false)
                }
            }
            Node::Blank(head) | Node::Sequence { head, .. } => {
                if head.as_ref().is_none_or(|h| subject.head() == Some(h)) {
                    rest(self)
                } else {
                    Ok(false)
                }
            }
            Node::Bind { slot, inner } => {
                let mark = self.trail.len();
                if self.bind(*slot, subject.clone()) && self.node(inner, subject, rest)? {
                    return Ok(true);
                }
                self.undo(mark);
                Ok(false)
            }
            Node::Apply { head, args } => {
                let View::Apply {
                    head: subject_head,
                    args: subject_args,
                } = subject.view()
                else {
                    return Ok(false);
                };
                self.node(head, subject_head, &mut |m| {
                    m.args(args, subject_args, rest)
                })
            }
            Node::Alternatives(alternatives) => {
                for a in alternatives {
                    let mark = self.trail.len();
                    if self.node(a, subject, rest)? {
                        return Ok(true);
                    }
                    self.undo(mark);
                }
                Ok(false)
            }
            Node::Optional { inner, .. } => self.node(inner, subject, rest),
            Node::Condition { inner, test } => self.node(inner, subject, &mut |m| {
                let test = substitute(test, &m.substitutions());
                if m.holds(&test)? { rest(m) } else { Ok(false) }
            }),
            Node::Test { inner, test } => self.node(inner, subject, &mut |m| {
                let test = Expr::apply(test, std::slice::from_ref(subject));
                if m.holds(&test)? { rest(m) } else { Ok(false) }
            }),
        }
    }

    /// Matches the argument nodes against the subject's arguments, then the
    /// rest. A sequence takes the fewest arguments first.
    fn args(
        &mut self,
        nodes: &[Node],
        subjects: &[Expr],
        rest: Rest<'_>,
    ) -> Result<bool, EvalError> {
        let Some((first, more)) = nodes.split_first() else {
            return if subjects.is_empty() {
                rest(self)
            } else {
                Ok(false)
            };
        };
        if let Some((slot, min, head)) = sequence(first) {
            for len in min..=subjects.len() {
                let run = &subjects[..len];
                if head.is_some_and(|h| run.iter().any(|s| s.head() != Some(h))) {
                    continue;
                }
                let mark = self.trail.len();
                let bound = match slot {
                    Some(slot) => self.bind(slot, Expr::apply(&Symbol::SEQUENCE.into(), run)),
                    None => true,
                };
                if bound && self.args(more, &subjects[len..], rest)? {
                    return Ok(true);
                }
                self.undo(mark);
            }
            return Ok(false);
        }
        if let Node::Optional { inner, default } = first {
            if let Some((subject, tail)) = subjects.split_first() {
                let mark = self.trail.len();
                if self.node(inner, subject, &mut |m| m.args(more, tail, rest))? {
                    return Ok(true);
                }
                self.undo(mark);
            }
            // Missing: the variable takes the default.
            let Some(default) = default else {
                return Ok(false);
            };
            let mark = self.trail.len();
            let bound = match inner.as_ref() {
                Node::Bind { slot, .. } => self.bind(*slot, default.clone()),
                _ => true,
            };
            if bound && self.args(more, subjects, rest)? {
                return Ok(true);
            }
            self.undo(mark);
            return Ok(false);
        }
        let Some((subject, tail)) = subjects.split_first() else {
            return Ok(false);
        };
        self.node(first, subject, &mut |m| m.args(more, tail, rest))
    }
}

/// `(slot, least length, head)` when `node` matches a run of arguments:
/// `__`, `___`, or one of them bound to a variable.
fn sequence(node: &Node) -> Option<(Option<usize>, usize, Option<&Expr>)> {
    match node {
        Node::Sequence { min, head } => Some((None, *min, head.as_ref())),
        Node::Bind { slot, inner } => match inner.as_ref() {
            Node::Sequence { min, head } => Some((Some(*slot), *min, head.as_ref())),
            _ => None,
        },
        _ => None,
    }
}

/// `e` with each symbol of `map` replaced by its value, except inside a
/// `Function` that binds it as a parameter (CAS-DOC-4 §3.1). Parts without
/// any of the symbols are kept as they are, by their metadata.
pub(crate) fn substitute(e: &Expr, map: &[(Symbol, Expr)]) -> Expr {
    if map.is_empty() {
        return e.clone();
    }
    e.replace(|part| {
        if map.iter().all(|(s, _)| !part.may_contain(*s)) {
            return Some(part.clone());
        }
        if let Some(s) = part.as_symbol() {
            return map.iter().find(|(x, _)| *x == s).map(|(_, v)| v.clone());
        }
        let parameters = function_parameters(part)?;
        let free: Vec<(Symbol, Expr)> = map
            .iter()
            .filter(|(x, _)| !parameters.contains(x))
            .cloned()
            .collect();
        if free.len() == map.len() {
            return None;
        }
        let [params, body] = part.args() else {
            return None;
        };
        let body = substitute(body, &free);
        Some(Expr::apply(
            &Symbol::FUNCTION.into(),
            &[params.clone(), body],
        ))
    })
}

/// The parameters of `Function[x, body]` or `Function[List[x, …], body]`.
fn function_parameters(e: &Expr) -> Option<Vec<Symbol>> {
    if e.head().and_then(Expr::as_symbol) != Some(Symbol::FUNCTION) {
        return None;
    }
    let [params, _] = e.args() else {
        return None;
    };
    match params.view() {
        View::Symbol(x) => Some(vec![x]),
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::LIST) => {
            args.iter().map(Expr::as_symbol).collect()
        }
        _ => None,
    }
}

/// `function[args…]` for a `Function` head: its body with the parameters
/// replaced (`Function[List[x, y], body]`, or `Function[x, body]`), or with
/// `Slot[n]` replaced by the `n`-th argument (`Function[body]`). `None` when
/// there are fewer arguments than parameters or slots, or the function is
/// malformed. The body is not evaluated here.
pub(crate) fn apply_function(function: &Expr, args: &[Expr]) -> Option<Expr> {
    match function.args() {
        [body] => {
            let mut missing = false;
            let result = body.replace(|part| {
                // Slots inside a nested pure function are its own; a function
                // with named parameters has none.
                if part.head().and_then(Expr::as_symbol) == Some(Symbol::FUNCTION)
                    && part.args().len() == 1
                {
                    return Some(part.clone());
                }
                if part.head().and_then(Expr::as_symbol) != Some(Symbol::SLOT) {
                    return None;
                }
                let n = match part.args() {
                    [] => Some(1),
                    [n] => n.as_i64(),
                    _ => None,
                };
                match n.and_then(|n| usize::try_from(n).ok()).filter(|&n| n >= 1) {
                    Some(n) if n <= args.len() => Some(args[n - 1].clone()),
                    _ => {
                        missing = true;
                        Some(part.clone())
                    }
                }
            });
            (!missing).then_some(result)
        }
        [_, body] => {
            let parameters = function_parameters(function)?;
            if args.len() < parameters.len() {
                return None;
            }
            let map: Vec<(Symbol, Expr)> =
                parameters.into_iter().zip(args.iter().cloned()).collect();
            Some(substitute(body, &map))
        }
        _ => None,
    }
}

/// A rule, `lhs -> rhs` or `lhs :> rhs`, compiled once (CAS-DOC-4 §5.3): a
/// program keeps it and applies it without compiling it again.
///
/// A right-hand side `Condition[rhs, test]` applies only where `test`, with
/// the bindings, evaluates to `True`; otherwise the match goes on to the
/// next way of matching.
///
/// ```
/// use cas_core::{Context, Expr, Rule};
///
/// let read = |s: &str| Expr::parse_full_form(s).unwrap();
/// let rule = Rule::new(&read("f[Pattern[x, Blank[]]]"), &read("g[x, x]")).unwrap();
/// let ctx = Context::new();
/// let e = ctx.replace_all(&read("h[f[1], f[a]]"), &[rule]).unwrap();
/// assert_eq!(format!("{e:?}"), "h[g[1, 1], g[a, a]]");
/// ```
#[derive(Clone, Debug)]
pub struct Rule {
    lhs: Pattern,
    rhs: Expr,
}

impl Rule {
    /// The rule `lhs -> rhs`.
    ///
    /// # Errors
    ///
    /// A [`PatternError`] if `lhs` is not a pattern.
    pub fn new(lhs: &Expr, rhs: &Expr) -> Result<Rule, PatternError> {
        Ok(Rule {
            lhs: Pattern::new(lhs)?,
            rhs: rhs.clone(),
        })
    }

    /// The rule `lhs -> rhs` for a compiled `lhs`.
    pub(crate) fn from_pattern(lhs: Pattern, rhs: Expr) -> Rule {
        Rule { lhs, rhs }
    }

    /// The rule an expression `Rule[lhs, rhs]` or `RuleDelayed[lhs, rhs]`
    /// stands for; `None` for anything else.
    pub(crate) fn of(e: &Expr) -> Option<Result<Rule, PatternError>> {
        let View::Apply {
            head,
            args: [lhs, rhs],
        } = e.view()
        else {
            return None;
        };
        matches!(head.as_symbol(), Some(Symbol::RULE | Symbol::RULE_DELAYED))
            .then(|| Rule::new(lhs, rhs))
    }

    /// The left-hand side.
    #[must_use]
    pub fn lhs(&self) -> &Pattern {
        &self.lhs
    }

    /// The right-hand side, with the pattern variables.
    #[must_use]
    pub fn rhs(&self) -> &Expr {
        &self.rhs
    }

    /// The right-hand side, instantiated for the first match of `subject`
    /// whose condition holds; not evaluated.
    pub(crate) fn apply(
        &self,
        ev: &mut Eval<'_>,
        subject: &Expr,
    ) -> Result<Option<Expr>, EvalError> {
        let variables = self.lhs.variables();
        let mut result = None;
        let condition = match self.rhs.view() {
            View::Apply {
                head,
                args: [body, test],
            } if head.as_symbol() == Some(Symbol::CONDITION) => Some((body, test)),
            _ => None,
        };
        self.lhs.match_with(ev, subject, &mut |ev, bound| {
            let map: Vec<(Symbol, Expr)> = variables
                .iter()
                .copied()
                .zip(bound.iter().cloned())
                .collect();
            match condition {
                Some((body, test)) => {
                    let value = ev.eval(&substitute(test, &map))?;
                    if value.as_symbol() != Some(Symbol::TRUE) {
                        return Ok(false);
                    }
                    result = Some(substitute(body, &map));
                }
                None => result = Some(substitute(&self.rhs, &map)),
            }
            Ok(true)
        })?;
        Ok(result)
    }
}

/// The rules of a `Replace` argument: one rule, or a list of them. `None` if
/// anything is not a valid rule.
pub(crate) fn rules(e: &Expr) -> Option<Vec<Rule>> {
    if e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
        e.args().iter().map(|r| Rule::of(r)?.ok()).collect()
    } else {
        Some(vec![Rule::of(e)?.ok()?])
    }
}

/// The first rule that rewrites `e`, applied.
pub(crate) fn rewrite(
    ev: &mut Eval<'_>,
    rules: &[Rule],
    e: &Expr,
) -> Result<Option<Expr>, EvalError> {
    for rule in rules {
        if let Some(result) = rule.apply(ev, e)? {
            return Ok(Some(result));
        }
    }
    Ok(None)
}

/// `e` with every part that a rule rewrites replaced, from the top: a
/// rewritten part is not looked into again. Not evaluated.
pub(crate) fn replace_all(ev: &mut Eval<'_>, rules: &[Rule], e: &Expr) -> Result<Expr, EvalError> {
    let mut error = None;
    let result = e.replace(|part| {
        if error.is_some() {
            return Some(part.clone());
        }
        match rewrite(ev, rules, part) {
            Ok(r) => r,
            Err(err) => {
                error = Some(err);
                Some(part.clone())
            }
        }
    });
    error.map_or(Ok(result), Err)
}

impl Context {
    /// Whether `e` matches `pattern` (`MatchQ`). `e` is not evaluated first.
    ///
    /// # Errors
    ///
    /// [`EvalError`] from evaluating a condition or a test.
    pub fn matches(&self, e: &Expr, pattern: &Pattern) -> Result<bool, EvalError> {
        run(self.limits(), |recursion| {
            pattern.matches(&mut Eval::new(self, recursion), e)
        })
    }

    /// `e` with every part that a rule rewrites replaced, then evaluated
    /// (`ReplaceAll`).
    ///
    /// # Errors
    ///
    /// [`EvalError`] from conditions, tests or the evaluation.
    pub fn replace_all(&self, e: &Expr, rules: &[Rule]) -> Result<Expr, EvalError> {
        run(self.limits(), |recursion| {
            let mut ev = Eval::new(self, recursion);
            let replaced = replace_all(&mut ev, rules, e)?;
            ev.eval(&replaced)
        })
    }
}

#[cfg(test)]
mod tests;
