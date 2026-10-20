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
//! bind the same expression, which is an O(1) check.
//!
//! The attributes of a head, as they are when the pattern is compiled, decide
//! how its arguments match (§3.2–3.3, CAS-57): in order; under `Flat`, a run
//! of arguments grouped in the head; under `Orderless`, as a multiset
//! (`multiset`); and under `OneIdentity`, a lone expression stands for the
//! application when the other arguments may be missing, so `a_. x_` matches
//! `y`. One match attempt has a budget of steps,
//! [`Limits::match_steps`](crate::Limits::match_steps).
//!
//! [`substitute`] instantiates a right-hand side: it replaces the variables,
//! but never inside a `Function` that binds them, so `f[a_] :> x -> a x`
//! keeps its parameter `x` (CAS-DOC-4 §3.1).

mod multiset;

use std::fmt;

use crate::attributes::Attributes;
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
    /// `h[p1, …]`, its head matched by a node too, its arguments as the
    /// head's attributes say.
    Apply {
        head: Box<Node>,
        args: Vec<Node>,
        mode: Mode,
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
    /// Compiles `e` with the attributes of the builtins only: a user head
    /// matches its arguments in order. [`Pattern::with_context`] knows the
    /// attributes of user heads.
    ///
    /// # Errors
    ///
    /// A [`PatternError`] for a malformed construct.
    pub fn new(e: &Expr) -> Result<Pattern, PatternError> {
        Pattern::compile(
            e,
            &|s| s.builtin_attributes().unwrap_or(Attributes::NONE),
            &|_| None,
        )
    }

    /// Compiles `e` with the attributes and the defaults that heads have in
    /// `ctx` (CAS-DOC-4 §4): how a head's arguments match depends on its
    /// attributes as they are now.
    ///
    /// # Errors
    ///
    /// A [`PatternError`] for a malformed construct.
    pub fn with_context(e: &Expr, ctx: &Context) -> Result<Pattern, PatternError> {
        Pattern::compile(e, &|s| ctx.attributes(s), &|s| ctx.default_of(s))
    }

    /// Compiles `e`: heads have `attributes`, and `x_.` as an argument of a
    /// user head `f` takes `defaults(f)` when it is missing.
    fn compile(
        e: &Expr,
        attributes: &dyn Fn(Symbol) -> Attributes,
        defaults: &dyn Fn(Symbol) -> Option<Expr>,
    ) -> Result<Pattern, PatternError> {
        let mut compiler = Compiler {
            variables: Vec::new(),
            heads: Vec::new(),
            attributes,
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
    ///
    /// # Errors
    ///
    /// [`EvalError::MatchLimit`] after
    /// [`Limits::match_steps`](crate::Limits::match_steps) steps, and what
    /// evaluating a condition or a test gives.
    pub(crate) fn match_with(
        &self,
        ev: &mut Eval<'_>,
        subject: &Expr,
        accept: Accept<'_>,
    ) -> Result<bool, EvalError> {
        if (subject.meta().symbols & self.required) != self.required {
            return Ok(false);
        }
        let budget = ev.context().limits().match_steps;
        let mut m = Matcher {
            ev,
            variables: &self.variables,
            bindings: vec![None; self.variables.len()],
            trail: Vec::new(),
            steps: 0,
            budget,
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
        let Node::Apply { head, args, mode } = root else {
            return Filter::default();
        };
        // A lone expression may stand for the application.
        if mode.lone {
            return Filter::default();
        }
        let (min, max) = mode.arity();
        Filter {
            head: match head.as_ref() {
                Node::Literal(h) => Some(h.clone()),
                _ => None,
            },
            min,
            max,
            // Under `Orderless`, any argument may come first.
            first: if mode.orderless() {
                First::Any
            } else {
                args.first().map_or(First::Any, first_argument)
            },
        }
    }
}

/// What a subject must be for a pattern to match it (CAS-DOC-4 §5.1): its
/// head, its number of arguments, and its first argument, or that
/// argument's head, when the pattern fixes them. Checked before matching,
/// it skips most rules of a head in O(1). Under `Flat` any number of
/// arguments from the least may group, under `Orderless` any may come
/// first, and under `OneIdentity` the subject may not be an application of
/// the head at all.
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
        Node::Apply { head, mode, .. } => match head.as_ref() {
            Node::Literal(h) if !mode.lone => First::Head(h.clone()),
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
    attributes: &'d dyn Fn(Symbol) -> Attributes,
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
            symbol => self.application(head, args, symbol)?,
        })
    }

    /// The node for `head[args…]` that is no pattern construct, `symbol`
    /// its head when it is one.
    fn application(
        &mut self,
        head: &Expr,
        args: &[Expr],
        symbol: Option<Symbol>,
    ) -> Result<Node, PatternError> {
        if let Some(s) = symbol.filter(|s| !s.is_builtin() && !self.heads.contains(s)) {
            self.heads.push(s);
        }
        let mut nodes = Vec::with_capacity(args.len());
        for (i, a) in args.iter().enumerate() {
            nodes.push(self.node(a, symbol.map(|s| (s, i + 1)))?);
        }
        let attributes = symbol.map_or(Attributes::NONE, |s| (self.attributes)(s));
        let mode = Mode::new(symbol, attributes, &nodes);
        Ok(Node::Apply {
            head: Box::new(self.node(head, None)?),
            args: nodes,
            mode,
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

/// How the arguments of an application match: by the attributes its head
/// had when the pattern was compiled (CAS-DOC-4 §3.2–3.4, §4).
#[derive(Clone, Debug)]
struct Mode {
    /// The head, when it is a symbol: under `Flat`, a group of several
    /// subjects is bound wrapped in it.
    head: Option<Expr>,
    /// The head's attributes, of which `Flat` and `Orderless` count here.
    attributes: Attributes,
    /// `OneIdentity`, and every argument but one may be missing: a lone
    /// expression then stands for the application.
    lone: bool,
    /// `Times`, whose lists do not commute (R13): a multiset match takes
    /// them in order.
    times: bool,
    /// The order in which the arguments take subjects: as written, or for a
    /// multiset literals first, then fixed heads, blanks and sequences.
    order: Vec<usize>,
    /// Per argument: whether it may take a group of several subjects.
    groups: Vec<bool>,
    /// From each position of `order` on: the fewest subjects the arguments
    /// take, and the most (`None`: any number).
    least: Vec<usize>,
    most: Vec<Option<usize>>,
}

impl Mode {
    /// The mode of an application with a symbol head that has `attributes`,
    /// or with a compound head (`None`), whose arguments match in order.
    fn new(head: Option<Symbol>, attributes: Attributes, args: &[Node]) -> Mode {
        let attributes = if head.is_some() {
            attributes
        } else {
            Attributes::NONE
        };
        let flat = attributes.contains(Attributes::FLAT);
        let h = head.map(Expr::from);
        let groups: Vec<bool> = args
            .iter()
            .map(|a| flat && h.as_ref().is_some_and(|h| groupable(a, h)))
            .collect();
        let mut order: Vec<usize> = (0..args.len()).collect();
        if attributes.contains(Attributes::ORDERLESS) {
            order.sort_by_key(|&i| category(&args[i]));
        }
        let mut least = vec![0; args.len() + 1];
        let mut most = vec![Some(0); args.len() + 1];
        for k in (0..args.len()).rev() {
            let i = order[k];
            let (fewest, max) = takes(&args[i], groups[i]);
            least[k] = least[k + 1] + fewest;
            most[k] = most[k + 1].zip(max).map(|(a, b)| a + b);
        }
        Mode {
            lone: attributes.contains(Attributes::ONE_IDENTITY)
                && least[0] <= 1
                && most[0] != Some(0),
            times: head == Some(Symbol::TIMES),
            head: h,
            attributes,
            order,
            groups,
            least,
            most,
        }
    }

    /// The fewest and the most arguments an application of the head may
    /// have to match (`None`: any number).
    fn arity(&self) -> (usize, Option<usize>) {
        (self.least[0], self.most[0])
    }

    fn flat(&self) -> bool {
        self.attributes.contains(Attributes::FLAT)
    }

    fn orderless(&self) -> bool {
        self.attributes.contains(Attributes::ORDERLESS)
    }
}

/// The fewest and the most subjects an argument node takes: one, or under
/// `Flat` (`group`) any number from one; a sequence any number from its
/// least; an optional argument may take none.
fn takes(node: &Node, group: bool) -> (usize, Option<usize>) {
    if let Some((_, min, _)) = sequence(node) {
        return (min, None);
    }
    let most = if group { None } else { Some(1) };
    match node {
        Node::Optional {
            default: Some(_), ..
        } => (0, most),
        _ => (1, most),
    }
}

/// Where an argument comes in the order of a multiset match (CAS-DOC-4
/// §3.3): literals, then fixed heads, then blanks, then applications that a
/// lone expression may stand for (`x_^n_.`, which match anything but should
/// take their own head's subjects, left by the blanks), then sequences; an
/// optional argument where what it holds comes.
fn category(node: &Node) -> u8 {
    if sequence(node).is_some() {
        return 4;
    }
    match node {
        Node::Literal(_) => 0,
        Node::Optional { inner, .. } => category(inner),
        _ if stands_alone(node) => 3,
        _ if fixed_head(node).is_some() => 1,
        _ => 2,
    }
}

/// Whether `node` is an application that a lone expression may stand for.
fn stands_alone(node: &Node) -> bool {
    match node {
        Node::Apply { mode, .. } => mode.lone,
        Node::Bind { inner, .. } | Node::Condition { inner, .. } | Node::Test { inner, .. } => {
            stands_alone(inner)
        }
        _ => false,
    }
}

/// The head of every expression `node` matches, when it fixes one.
fn fixed_head(node: &Node) -> Option<&Expr> {
    match node {
        Node::Literal(e) => e.head(),
        Node::Blank(head) => head.as_ref(),
        Node::Bind { inner, .. } | Node::Condition { inner, .. } | Node::Test { inner, .. } => {
            fixed_head(inner)
        }
        Node::Apply { head, mode, .. } => match head.as_ref() {
            Node::Literal(h) if !mode.lone => Some(h),
            _ => None,
        },
        Node::Sequence { .. } | Node::Alternatives(_) | Node::Optional { .. } => None,
    }
}

/// Whether `node` may match an application of `h`: a group of several
/// subjects under `Flat`.
fn groupable(node: &Node, h: &Expr) -> bool {
    match node {
        Node::Literal(e) => e.head() == Some(h),
        Node::Blank(head) => head.as_ref().is_none_or(|g| g == h),
        Node::Sequence { .. } => false,
        Node::Bind { inner, .. }
        | Node::Condition { inner, .. }
        | Node::Test { inner, .. }
        | Node::Optional { inner, .. } => groupable(inner, h),
        Node::Apply { head, mode, .. } => match head.as_ref() {
            Node::Literal(g) => mode.lone || g == h,
            _ => true,
        },
        Node::Alternatives(alternatives) => alternatives.iter().any(|a| groupable(a, h)),
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
            Node::Apply { head, args, mode } => {
                // A lone expression standing for the application need not
                // have its head.
                if !mode.lone {
                    head.required(bits);
                }
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
    /// The steps taken in this attempt, and the most it may take.
    steps: u64,
    budget: u64,
}

impl Matcher<'_, '_, '_> {
    /// Counts a step: a node tried, or a way of splitting arguments.
    fn step(&mut self) -> Result<(), EvalError> {
        self.steps += 1;
        if self.steps > self.budget {
            Err(EvalError::MatchLimit)
        } else {
            Ok(())
        }
    }

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
        let value = self.ev.test(test)?;
        Ok(value.as_symbol() == Some(Symbol::TRUE))
    }

    /// Matches `node` against `subject`, then the rest. Each node takes a
    /// level of the recursion limit while it is on the stack: the rest of
    /// the match runs inside it.
    fn node(&mut self, node: &Node, subject: &Expr, rest: Rest<'_>) -> Result<bool, EvalError> {
        self.step()?;
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
            Node::Apply { head, args, mode } => {
                if let View::Apply {
                    head: subject_head,
                    args: subject_args,
                } = subject.view()
                {
                    if self.node(head, subject_head, &mut |m| {
                        m.arguments(mode, args, subject_args, rest)
                    })? {
                        return Ok(true);
                    }
                    // An application of the head matches only as one.
                    if mode.head.as_ref() == Some(subject_head) {
                        return Ok(false);
                    }
                }
                if mode.lone {
                    // `OneIdentity`: the subject is the only argument.
                    return self.arguments(mode, args, std::slice::from_ref(subject), rest);
                }
                Ok(false)
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

    /// Matches the argument nodes against the subject's arguments as `mode`
    /// says, then the rest.
    fn arguments(
        &mut self,
        mode: &Mode,
        nodes: &[Node],
        subjects: &[Expr],
        rest: Rest<'_>,
    ) -> Result<bool, EvalError> {
        if mode.orderless() {
            self.multiset(mode, nodes, subjects, rest)
        } else if mode.flat() {
            self.runs(mode, nodes, 0, subjects, rest)
        } else {
            self.args(nodes, subjects, rest)
        }
    }

    /// Matches the argument nodes against the subject's arguments, in order,
    /// then the rest. A sequence takes the fewest arguments first.
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
                self.step()?;
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
            return self.missing(inner, default, &mut |m| m.args(more, subjects, rest));
        }
        let Some((subject, tail)) = subjects.split_first() else {
            return Ok(false);
        };
        self.node(first, subject, &mut |m| m.args(more, tail, rest))
    }

    /// An optional argument left out: its variable, if any, takes `default`,
    /// then the rest.
    fn missing(&mut self, inner: &Node, default: &Expr, rest: Rest<'_>) -> Result<bool, EvalError> {
        self.step()?;
        let mark = self.trail.len();
        let bound = match inner {
            Node::Bind { slot, .. } => self.bind(*slot, default.clone()),
            _ => true,
        };
        if bound && rest(self)? {
            return Ok(true);
        }
        self.undo(mark);
        Ok(false)
    }

    /// The argument nodes of a `Flat` head that is not `Orderless`, from the
    /// `k`-th, against the subjects left: each node takes a run, grouped in
    /// the head when it has several subjects, shortest first; the last takes
    /// what is left (CAS-DOC-4 §3.2).
    fn runs(
        &mut self,
        mode: &Mode,
        nodes: &[Node],
        k: usize,
        subjects: &[Expr],
        rest: Rest<'_>,
    ) -> Result<bool, EvalError> {
        let Some(node) = nodes.get(k) else {
            return if subjects.is_empty() {
                rest(self)
            } else {
                Ok(false)
            };
        };
        // What the nodes after this one leave it.
        let Some(spare) = subjects.len().checked_sub(mode.least[k + 1]) else {
            return Ok(false);
        };
        let last = k + 1 == nodes.len();
        let lengths = |least: usize, most: usize| {
            if last {
                subjects.len()..=subjects.len()
            } else {
                least..=most
            }
        };
        if let Some((slot, min, head)) = sequence(node) {
            for len in lengths(min, spare) {
                self.step()?;
                let run = &subjects[..len];
                if len < min || head.is_some_and(|h| run.iter().any(|s| s.head() != Some(h))) {
                    continue;
                }
                let mark = self.trail.len();
                let bound = match slot {
                    Some(slot) => self.bind(slot, Expr::apply(&Symbol::SEQUENCE.into(), run)),
                    None => true,
                };
                if bound && self.runs(mode, nodes, k + 1, &subjects[len..], rest)? {
                    return Ok(true);
                }
                self.undo(mark);
            }
            return Ok(false);
        }
        let (inner, default) = match node {
            Node::Optional { inner, default } => (inner.as_ref(), default.as_ref()),
            _ => (node, None),
        };
        let most = if mode.groups[k] { spare } else { spare.min(1) };
        for len in lengths(1, most) {
            if len == 0 || len > 1 && !mode.groups[k] {
                continue;
            }
            let item = match &subjects[..len] {
                [one] => one.clone(),
                run => Expr::apply(mode.head.as_ref().expect("a Flat head"), run),
            };
            if self.node(inner, &item, &mut |m| {
                m.runs(mode, nodes, k + 1, &subjects[len..], rest)
            })? {
                return Ok(true);
            }
        }
        match default {
            Some(default) => self.missing(inner, default, &mut |m| {
                m.runs(mode, nodes, k + 1, subjects, rest)
            }),
            None => Ok(false),
        }
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
    /// The rule `lhs -> rhs`, with the attributes of the builtins only, as
    /// [`Pattern::new`].
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

    /// The rule `lhs -> rhs`, with the attributes and defaults that heads
    /// have in `ctx`, as [`Pattern::with_context`].
    ///
    /// # Errors
    ///
    /// A [`PatternError`] if `lhs` is not a pattern.
    pub fn with_context(lhs: &Expr, rhs: &Expr, ctx: &Context) -> Result<Rule, PatternError> {
        Ok(Rule {
            lhs: Pattern::with_context(lhs, ctx)?,
            rhs: rhs.clone(),
        })
    }

    /// The rule `lhs -> rhs` for a compiled `lhs`.
    pub(crate) fn from_pattern(lhs: Pattern, rhs: Expr) -> Rule {
        Rule { lhs, rhs }
    }

    /// The rule an expression `Rule[lhs, rhs]` or `RuleDelayed[lhs, rhs]`
    /// stands for, compiled in `ctx`; `None` for anything else.
    pub(crate) fn of(e: &Expr, ctx: &Context) -> Option<Result<Rule, PatternError>> {
        let View::Apply {
            head,
            args: [lhs, rhs],
        } = e.view()
        else {
            return None;
        };
        matches!(head.as_symbol(), Some(Symbol::RULE | Symbol::RULE_DELAYED))
            .then(|| Rule::with_context(lhs, rhs, ctx))
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
                    let value = ev.test(&substitute(test, &map))?;
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

/// The rules of a `Replace` argument, compiled in `ctx`: one rule, or a list
/// of them. `None` if anything is not a valid rule.
pub(crate) fn rules(e: &Expr, ctx: &Context) -> Option<Vec<Rule>> {
    if e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
        e.args().iter().map(|r| Rule::of(r, ctx)?.ok()).collect()
    } else {
        Some(vec![Rule::of(e, ctx)?.ok()?])
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
