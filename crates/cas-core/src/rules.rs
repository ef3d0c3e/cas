//! The rule store and dispatch (CAS-DOC-4 §5, CAS-56).
//!
//! A symbol's rules are its [`SymbolRules`]: an exact table for ground
//! left-hand sides, keyed by the left-hand side evaluated and normalised
//! (§5.1); patterned down-rules in definition order (R4); the up-rules it
//! carries for other heads (§5.2); sub-rules for `f[…][…]`; a native
//! function; and the default of `x_.` among its arguments. They live in a
//! persistent map of the context, so a definition copies O(log n) of it.
//!
//! Step 6 of the evaluator tries, in order (R10): the up-rules of each
//! argument's symbol, the head's native function, the system rules of a
//! builtin head (§5.4) or the rules of a user's, exact first, and for a
//! compound head the sub-rules of its innermost symbol. Each patterned rule
//! has a filter on the subject's head, number of arguments and first
//! argument, checked before matching.
//!
//! A left-hand side is stored with its arguments evaluated and normalised, as
//! the evaluator presents subjects. Changing the attributes or the default of
//! a symbol recompiles the rules that depend on them, which an index finds:
//! the rules of that symbol, and the rules whose left-hand sides mention it
//! as a head. Recompiling sorts and flattens the stored left-hand sides again
//! under the new attributes; it does not evaluate them again (§4).

use std::fmt;
use std::sync::{Arc, LazyLock};

use foldhash::fast::FixedState;
use hashbrown::hash_map::Entry;
use hashbrown::{HashMap, HashSet};

use crate::attributes::Attributes;
use crate::context::{Context, ContextError, EvalError};
use crate::eval::{Eval, Normal, holds, normalise, run};
use crate::expr::{Expr, View, symbol_bit};
use crate::pattern::{Filter, Pattern, PatternError, Rule, substitute};
use crate::symbol::Symbol;
use crate::symbol_map::SymbolMap;

/// The native function of a user's symbol (CAS-DOC-4 §6), tried after the
/// up-rules and before the symbol's rules. Like a builtin's
/// [`NativeFn`](crate::NativeFn), it returns the rewrite of the application
/// or `None`, and must be deterministic in the state and its input and change
/// nothing visible, since results are stamped and memoised.
pub type UserNative =
    Arc<dyn Fn(&mut Eval<'_>, &Expr) -> Result<Option<Expr>, EvalError> + Send + Sync>;

/// Why a definition is refused. The context is left as it was.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum DefineError {
    /// The symbol the definition attaches to is `PROTECTED`: a builtin, or a
    /// user symbol given that attribute.
    Protected(Symbol),
    /// A value for a symbol with up-rules, which it would hide (O6).
    HasUpRules(Symbol),
    /// An up-rule for a symbol with a value: the value replaces the symbol
    /// before its up-rules could apply (O6).
    HasValue(Symbol),
    /// The tag of an up-rule is not at level 1 of its left-hand side: an
    /// argument, the head of one, or the head a blank requires (`_f`).
    NotATag(Symbol),
    /// The left-hand side has no symbol to attach the definition to.
    NoHead(Expr),
    /// The left-hand side is not a valid pattern.
    Pattern(PatternError),
    /// Evaluating a side of the definition failed.
    Eval(EvalError),
}

impl fmt::Display for DefineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DefineError::Protected(s) => ContextError::Protected(*s).fmt(f),
            DefineError::HasUpRules(s) => ContextError::HasUpRules(*s).fmt(f),
            DefineError::HasValue(s) => {
                write!(f, "symbol {s:?} has a value, which would hide its up-rules")
            }
            DefineError::NotATag(s) => write!(
                f,
                "symbol {s:?} is not an argument of the left-hand side or the head of one"
            ),
            DefineError::NoHead(e) => write!(f, "{e:?} has no symbol to define"),
            DefineError::Pattern(e) => e.fmt(f),
            DefineError::Eval(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for DefineError {}

impl From<ContextError> for DefineError {
    fn from(e: ContextError) -> DefineError {
        match e {
            ContextError::Protected(s) => DefineError::Protected(s),
            ContextError::HasUpRules(s) => DefineError::HasUpRules(s),
        }
    }
}

impl From<EvalError> for DefineError {
    fn from(e: EvalError) -> DefineError {
        DefineError::Eval(e)
    }
}

impl From<PatternError> for DefineError {
    fn from(e: PatternError) -> DefineError {
        DefineError::Pattern(e)
    }
}

/// The rules of a context.
#[derive(Clone, Default)]
pub(crate) struct Rules {
    symbols: SymbolMap<SymbolRules>,
    /// For each symbol, the symbols whose rules depend on its attributes or
    /// default: those to recompile when they change. Usually the symbol
    /// itself and a few more, so a small table.
    dependents: SymbolMap<HashSet<Symbol, FixedState>>,
    /// How many symbols have up-rules, and a Bloom filter of them. With none,
    /// the common case, step 6 skips up-rules at once; otherwise it skips an
    /// application whose own filter has none of these bits, and looks up
    /// only the arguments the filter lets through.
    up_count: usize,
    up_bits: u64,
}

impl Rules {
    pub(crate) fn has_up_rules(&self, s: Symbol) -> bool {
        self.symbols.get(s).is_some_and(|r| !r.up.rules.is_empty())
    }

    /// Records that rules of `owner` depend on each symbol of `on`.
    fn depend(&mut self, on: &[Symbol], owner: Symbol) {
        for &s in on {
            if !self
                .dependents
                .get(s)
                .is_some_and(|owners| owners.contains(&owner))
            {
                self.dependents.entry(s, HashSet::default).insert(owner);
            }
        }
    }
}

impl fmt::Debug for Rules {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rules")
            .field("symbols", &self.symbols.len())
            .field("up", &self.up_count)
            .finish_non_exhaustive()
    }
}

/// The rules of one symbol.
#[derive(Clone, Default)]
pub(crate) struct SymbolRules {
    /// Ground left-hand sides, evaluated and normalised (`f[2] = 3`).
    exact: HashMap<Expr, Exact, FixedState>,
    down: RuleList,
    /// Rules for other heads that this symbol carries (§5.2).
    up: RuleList,
    /// Rules for `f[…][…]`, kept on the innermost head `f`.
    sub: RuleList,
    native: Option<UserNative>,
    /// The value of `x_.` as a missing argument of this head.
    default: Option<Expr>,
    /// The definitions made so far: the time of the next one.
    next: u64,
}

impl fmt::Debug for SymbolRules {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SymbolRules")
            .field("exact", &self.exact.len())
            .field("down", &self.down.rules.len())
            .field("up", &self.up.rules.len())
            .field("sub", &self.sub.rules.len())
            .field("native", &self.native.is_some())
            .field("default", &self.default)
            .finish_non_exhaustive()
    }
}

/// The right-hand side of an exact definition, with the times of the first
/// definition of its left-hand side, which orders the table, and of the
/// last, which wrote `rhs`.
#[derive(Clone, Debug)]
struct Exact {
    rhs: Expr,
    rank: u64,
    written: u64,
}

/// Patterned rules in the order they are tried.
#[derive(Clone, Debug, Default)]
struct RuleList {
    rules: Vec<Stored>,
    /// Where each key is in `rules`.
    index: HashMap<Expr, usize, FixedState>,
}

impl RuleList {
    /// Adds `stored` in place of the rule with the same key, or last (R4).
    fn put(&mut self, stored: Stored) {
        match self.index.entry(stored.key.clone()) {
            Entry::Occupied(o) => self.rules[*o.get()] = stored,
            Entry::Vacant(v) => {
                v.insert(self.rules.len());
                self.rules.push(stored);
            }
        }
    }

    /// The rules again, those that depend on `s` recompiled in `ctx`.
    /// Rules whose keys become equal collapse as if defined again in order:
    /// into the place of the first, with the last one written.
    fn recompiled(&self, ctx: &Context, s: Symbol) -> RuleList {
        let mut list = RuleList::default();
        for stored in &self.rules {
            let new = if stored.depends.contains(&s) {
                let lhs = renormalised(ctx, &stored.lhs);
                // The same constructs as before: compiling cannot fail.
                ctx.compile(lhs, stored.rule.rhs().clone(), stored.written)
                    .unwrap_or_else(|_| stored.clone())
            } else {
                stored.clone()
            };
            match list.index.entry(new.key.clone()) {
                Entry::Occupied(o) => {
                    let old = &mut list.rules[*o.get()];
                    if new.written > old.written {
                        *old = new;
                    }
                }
                Entry::Vacant(v) => {
                    v.insert(list.rules.len());
                    list.rules.push(new);
                }
            }
        }
        list
    }
}

/// A patterned rule as stored.
#[derive(Clone, Debug)]
struct Stored {
    /// The left-hand side, evaluated and normalised.
    lhs: Expr,
    /// The left-hand side and the right-hand side's condition, with the
    /// variables renamed by slot: rules with equal keys replace each other.
    key: Expr,
    rule: Rule,
    filter: Filter,
    /// The user heads in the left-hand side, on whose attributes and
    /// defaults the rule depends.
    depends: Vec<Symbol>,
    /// The time of the definition, among those of the symbol.
    written: u64,
}

/// The system rules of the builtins (§5.4), built on first use and shared by
/// every context. CAS-59 adds those of `D`.
static SYSTEM: LazyLock<Rules> = LazyLock::new(|| system(&[]));

/// The rules of `definitions`, delayed, on builtin heads. They are written in
/// normal form and not evaluated, so building them evaluates nothing.
///
/// # Panics
///
/// If a definition has no head or is not a valid pattern.
fn system(definitions: &[(Expr, Expr)]) -> Rules {
    let mut ctx = Context::new();
    for (lhs, rhs) in definitions {
        let (owner, kind) = attachment(lhs).expect("a system rule has a head");
        ctx.place(owner, kind, lhs.clone(), rhs.clone())
            .expect("a valid system rule");
    }
    ctx.rules().clone()
}

/// Step 6, first: the up-rules of the arguments' symbols, each symbol once.
/// Inlined, so that a context without up-rules pays one comparison, and an
/// application whose Bloom filter has none of the tags one more.
#[inline]
pub(crate) fn up(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    let rules = ev.context().rules();
    if rules.up_count == 0 || e.meta().symbols & rules.up_bits == 0 {
        return Ok(None);
    }
    up_rules(ev, e)
}

fn up_rules(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    let rules = ev.context().rules();
    let mut tried: Vec<Symbol> = Vec::new();
    for a in e.args() {
        let Some(tag) = a.as_symbol().or_else(|| head_symbol(a)) else {
            continue;
        };
        if symbol_bit(tag) & rules.up_bits == 0 || tried.contains(&tag) {
            continue;
        }
        tried.push(tag);
        if let Some(own) = rules.symbols.get(tag) {
            if let Some(rewritten) = first(ev, &own.up, e)? {
                return Ok(Some(rewritten));
            }
        }
    }
    Ok(None)
}

/// Step 6 for a symbol head, after a builtin's native function: a user's
/// native function, then the system rules of a builtin or the rules of a
/// user's symbol, exact first. Only a builtin head forces the system rules.
#[inline]
pub(crate) fn down(ev: &mut Eval<'_>, e: &Expr, head: Symbol) -> Result<Option<Expr>, EvalError> {
    let rules = if head.is_builtin() {
        &*SYSTEM
    } else {
        ev.context().rules()
    };
    match rules.symbols.get(head) {
        Some(own) => own_rules(own, ev, e),
        None => Ok(None),
    }
}

/// A symbol's native function, then its rules, exact ones first.
fn own_rules(own: &SymbolRules, ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    if let Some(native) = &own.native {
        if let Some(rewritten) = native(ev, e)? {
            return Ok(Some(rewritten));
        }
    }
    if let Some(exact) = own.exact.get(e) {
        return Ok(Some(exact.rhs.clone()));
    }
    first(ev, &own.down, e)
}

/// Step 6 for any other head: the sub-rules of its innermost symbol.
pub(crate) fn sub(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some(f) = innermost(e) else {
        return Ok(None);
    };
    if f.is_builtin() {
        sub_rules(&SYSTEM, ev, e, f)
    } else {
        sub_rules(ev.context().rules(), ev, e, f)
    }
}

/// The sub-rules of `f` in `rules`.
fn sub_rules(
    rules: &Rules,
    ev: &mut Eval<'_>,
    e: &Expr,
    f: Symbol,
) -> Result<Option<Expr>, EvalError> {
    match rules.symbols.get(f) {
        Some(own) => first(ev, &own.sub, e),
        None => Ok(None),
    }
}

/// The first rule of `list` that rewrites `e`, applied.
fn first(ev: &mut Eval<'_>, list: &RuleList, e: &Expr) -> Result<Option<Expr>, EvalError> {
    for stored in &list.rules {
        if stored.filter.admits(e) {
            if let Some(rewritten) = stored.rule.apply(ev, e)? {
                return Ok(Some(rewritten));
            }
        }
    }
    Ok(None)
}

fn head_symbol(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// The symbol at the bottom of the heads of `e`, as in `f` for `f[a][b]`.
fn innermost(e: &Expr) -> Option<Symbol> {
    let mut e = e;
    while let Some(head) = e.head() {
        e = head;
    }
    e.as_symbol()
}

/// What a definition makes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    /// The value of a symbol.
    Value,
    /// A down-rule, or an exact one, of the head.
    Down,
    /// A sub-rule of the innermost head.
    Sub,
}

/// `lhs` inside `HoldPattern` and `Condition`.
fn core(lhs: &Expr) -> &Expr {
    let mut e = lhs;
    loop {
        match (head_symbol(e), e.args()) {
            (Some(Symbol::HOLD_PATTERN), [p]) | (Some(Symbol::CONDITION), [p, _]) => e = p,
            _ => return e,
        }
    }
}

/// The symbol a definition of `lhs` attaches to, and what it makes: a symbol
/// gets a value, `f[…]` a rule of `f`, and `f[…][…]` a sub-rule of `f`,
/// through `HoldPattern` and `Condition`.
fn attachment(lhs: &Expr) -> Result<(Symbol, Kind), DefineError> {
    let core = core(lhs);
    match core.view() {
        View::Symbol(s) if core == lhs => Ok((s, Kind::Value)),
        View::Apply { head, .. } => match head.as_symbol() {
            Some(s) => Ok((s, Kind::Down)),
            None => innermost(head)
                .map(|s| (s, Kind::Sub))
                .ok_or_else(|| DefineError::NoHead(lhs.clone())),
        },
        _ => Err(DefineError::NoHead(lhs.clone())),
    }
}

/// The symbol an argument pattern of an up-rule ties it to (§5.2): an atomic
/// symbol, the head of an application, or the head a blank requires (`_f`,
/// `x_f`), through `Pattern`, `Condition`, `PatternTest`, `Optional` and
/// `HoldPattern`.
fn tag(p: &Expr) -> Option<Symbol> {
    match p.view() {
        View::Symbol(s) => Some(s),
        View::Apply { head, args } => match head.as_symbol()? {
            Symbol::BLANK | Symbol::BLANK_SEQUENCE | Symbol::BLANK_NULL_SEQUENCE => {
                args.first()?.as_symbol()
            }
            Symbol::PATTERN => tag(args.get(1)?),
            Symbol::CONDITION | Symbol::PATTERN_TEST | Symbol::OPTIONAL | Symbol::HOLD_PATTERN => {
                tag(args.first()?)
            }
            s => Some(s),
        },
        _ => None,
    }
}

/// `lhs` as subjects come to it: its arguments evaluated, except those its
/// head holds, and normalised under its head's attributes (§5.1). A compound
/// head is normalised the same way; `HoldPattern` keeps its part as written,
/// and `Condition` its test.
fn normalised(ev: &mut Eval<'_>, lhs: &Expr) -> Result<Expr, EvalError> {
    let View::Apply { head, args } = lhs.view() else {
        return Ok(lhs.clone());
    };
    match (head.as_symbol(), args) {
        (Some(Symbol::HOLD_PATTERN), _) => return Ok(lhs.clone()),
        (Some(Symbol::CONDITION), [p, test]) => {
            return Ok(Expr::apply(head, &[normalised(ev, p)?, test.clone()]));
        }
        _ => {}
    }
    let head = if head.is_atom() {
        head.clone()
    } else {
        normalised(ev, head)?
    };
    let attributes = head
        .as_symbol()
        .map_or(Attributes::NONE, |s| ev.context().attributes(s));
    let mut evaluated = Vec::with_capacity(args.len());
    for (i, a) in args.iter().enumerate() {
        evaluated.push(if holds(attributes, i) {
            a.clone()
        } else {
            ev.eval(a)?
        });
    }
    Ok(match normalise(&head, attributes, &evaluated)? {
        Normal::Changed(args) => Expr::apply(&head, &args),
        // A `Listable` head over lists: no subject ever has this form.
        Normal::Same | Normal::Threaded(_) => Expr::apply(&head, &evaluated),
    })
}

/// A stored left-hand side normalised again under the attributes of `ctx`,
/// as [`normalised`] would but without evaluating: sorted under `Orderless`
/// heads, flattened under `Flat` ones, `Sequence`s spliced, bottom up.
/// Held arguments stay as they are, and nothing is threaded.
fn renormalised(ctx: &Context, lhs: &Expr) -> Expr {
    match (head_symbol(lhs), lhs.args()) {
        (Some(Symbol::HOLD_PATTERN), _) => lhs.clone(),
        (Some(Symbol::CONDITION), [p, test]) => Expr::apply(
            &Symbol::CONDITION.into(),
            &[renormalised(ctx, p), test.clone()],
        ),
        _ => reevaluated(ctx, lhs),
    }
}

/// `e`, which was evaluated, normalised again bottom up (see
/// [`renormalised`]).
fn reevaluated(ctx: &Context, e: &Expr) -> Expr {
    let View::Apply { head, args } = e.view() else {
        return e.clone();
    };
    let new_head = reevaluated(ctx, head);
    let attributes = new_head
        .as_symbol()
        .map_or(Attributes::NONE, |s| ctx.attributes(s));
    let new_args: Vec<Expr> = args
        .iter()
        .enumerate()
        .map(|(i, a)| {
            if holds(attributes, i) {
                a.clone()
            } else {
                reevaluated(ctx, a)
            }
        })
        .collect();
    match normalise(
        &new_head,
        attributes.difference(Attributes::LISTABLE),
        &new_args,
    ) {
        Ok(Normal::Changed(args)) => Expr::apply(&new_head, &args),
        _ if new_head == *head && new_args == args => e.clone(),
        _ => Expr::apply(&new_head, &new_args),
    }
}

/// The user symbols that are heads in `e`.
fn user_heads(e: &Expr) -> Vec<Symbol> {
    let mut heads = Vec::new();
    for part in e.distinct() {
        if let Some(h) = head_symbol(part).filter(|h| !h.is_builtin()) {
            if !heads.contains(&h) {
                heads.push(h);
            }
        }
    }
    heads
}

/// The key of a rule: `lhs` with each variable renamed `Slot[i]` by its slot
/// where it is bound or tested, with the condition of the right-hand side.
/// Rules with equal keys have the same left-hand side and condition up to
/// the names of their variables (R4); a literal symbol that shares a
/// variable's name keeps it.
fn key(lhs: &Expr, rhs: &Expr, variables: &[Symbol]) -> Expr {
    let slots: Vec<(Symbol, Expr)> = variables
        .iter()
        .zip(1i64..)
        .map(|(&x, i)| (x, Expr::apply(&Symbol::SLOT.into(), &[i.into()])))
        .collect();
    let lhs = rename(lhs, &slots);
    match rhs.view() {
        View::Apply {
            head,
            args: [_, test],
        } if head.as_symbol() == Some(Symbol::CONDITION) => {
            Expr::apply(head, &[lhs, substitute(test, &slots)])
        }
        _ => lhs,
    }
}

fn rename(e: &Expr, slots: &[(Symbol, Expr)]) -> Expr {
    e.replace(|part| {
        if !part.has_pattern() {
            return Some(part.clone());
        }
        let View::Apply { head, args } = part.view() else {
            return Some(part.clone());
        };
        match (head.as_symbol(), args) {
            (Some(Symbol::PATTERN), [x, p]) => {
                let x = x
                    .as_symbol()
                    .and_then(|x| slots.iter().find(|(v, _)| *v == x))
                    .map_or_else(|| x.clone(), |(_, slot)| slot.clone());
                Some(Expr::apply(head, &[x, rename(p, slots)]))
            }
            (Some(Symbol::CONDITION), [p, test]) => Some(Expr::apply(
                head,
                &[rename(p, slots), substitute(test, slots)],
            )),
            _ => None,
        }
    })
}

impl SymbolRules {
    /// The time of a new definition.
    fn tick(&mut self) -> u64 {
        self.next += 1;
        self.next - 1
    }

    /// Adds `lhs = rhs` to the exact table, in place of a definition of the
    /// same left-hand side.
    fn insert_exact(&mut self, lhs: Expr, rhs: Expr) {
        let now = self.tick();
        match self.exact.entry(lhs) {
            Entry::Occupied(mut o) => {
                let x = o.get_mut();
                x.rhs = rhs;
                x.written = now;
            }
            Entry::Vacant(v) => {
                v.insert(Exact {
                    rhs,
                    rank: now,
                    written: now,
                });
            }
        }
    }

    /// These rules again, those that depend on `s` recompiled in `ctx`.
    fn recompiled(&self, ctx: &Context, s: Symbol) -> SymbolRules {
        let mut exact: Vec<(&Expr, &Exact)> = self.exact.iter().collect();
        exact.sort_unstable_by_key(|(_, x)| x.rank);
        let mut new = SymbolRules {
            down: self.down.recompiled(ctx, s),
            up: self.up.recompiled(ctx, s),
            sub: self.sub.recompiled(ctx, s),
            native: self.native.clone(),
            default: self.default.clone(),
            next: self.next,
            ..SymbolRules::default()
        };
        // Left-hand sides that become equal collapse as if defined again
        // in order: with the rank of the first, the value of the last.
        for (lhs, x) in exact {
            let lhs = if lhs.may_contain(s) {
                renormalised(ctx, lhs)
            } else {
                lhs.clone()
            };
            match new.exact.entry(lhs) {
                Entry::Occupied(mut o) => {
                    let old = o.get_mut();
                    if x.written > old.written {
                        old.rhs = x.rhs.clone();
                        old.written = x.written;
                    }
                }
                Entry::Vacant(v) => {
                    v.insert(x.clone());
                }
            }
        }
        new
    }
}

impl Context {
    /// `lhs = rhs` (`Set`): `rhs` is evaluated now (CAS-DOC-4 §5.3).
    ///
    /// - A symbol gets the value.
    /// - `f[…]` without patterns goes into the exact table of `f`, keyed by
    ///   its arguments evaluated and normalised: `f[1 + 1] = 3` is the rule
    ///   for `f[2]`.
    /// - `f[…]` with patterns is a down-rule of `f`, tried in definition
    ///   order, after the exact ones. A rule with the same left-hand side and
    ///   condition up to the names of its variables replaces the old one in
    ///   place.
    /// - `f[…][…]` is a sub-rule of `f`.
    ///
    /// `HoldPattern[lhs]` keeps `lhs` as written; `Condition[lhs, test]`
    /// applies only where `test` holds.
    ///
    /// ```
    /// use cas_core::{Context, Expr};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let mut ctx = Context::new();
    /// ctx.define(&read("fact[Plus[1, -1]]"), &read("1")).unwrap();
    /// ctx.define_delayed(
    ///     &read("fact[Pattern[n, Blank[]]]"),
    ///     &read("Times[n, fact[Plus[n, -1]]]"),
    /// )
    /// .unwrap();
    /// assert_eq!(ctx.eval(&read("fact[5]")).unwrap(), Expr::from(120));
    /// ```
    ///
    /// # Errors
    ///
    /// A [`DefineError`] if the symbol is protected, the left-hand side is
    /// not a pattern or has no symbol, the symbol has up-rules and `lhs` is
    /// the symbol, or evaluation fails.
    pub fn define(&mut self, lhs: &Expr, rhs: &Expr) -> Result<(), DefineError> {
        run(self.limits(), |recursion| {
            self.definition(recursion, lhs, rhs, true)
        })
    }

    /// `lhs := rhs` (`SetDelayed`): as [`Context::define`], but `rhs` is
    /// evaluated at each application.
    ///
    /// # Errors
    ///
    /// As [`Context::define`].
    pub fn define_delayed(&mut self, lhs: &Expr, rhs: &Expr) -> Result<(), DefineError> {
        run(self.limits(), |recursion| {
            self.definition(recursion, lhs, rhs, false)
        })
    }

    /// Delayed definitions in order, each as [`Context::define_delayed`] in
    /// the state those before it made, in one call: a library loads on one
    /// evaluation thread. All or nothing: on an error the context is as it
    /// was.
    ///
    /// # Errors
    ///
    /// The first [`DefineError`].
    pub fn define_many<'a>(
        &mut self,
        definitions: impl IntoIterator<Item = (&'a Expr, &'a Expr)>,
    ) -> Result<(), DefineError> {
        let definitions: Vec<(&Expr, &Expr)> = definitions.into_iter().collect();
        let before = self.clone();
        let result = run(self.limits(), |recursion| {
            definitions
                .iter()
                .try_for_each(|&(lhs, rhs)| self.definition(recursion, lhs, rhs, false))
        });
        if result.is_err() {
            *self = before;
        }
        result
    }

    /// `tag /: lhs = rhs` (`TagSet`, and `UpSet` for one tag): an up-rule
    /// carried by `tag`, which occurs at level 1 of `lhs` as an argument,
    /// the head of one, or the head a blank requires: `Power[phi, 2]`,
    /// `Length[vec[x___]]` or `g[Blank[f]]` (§5.2, O6). `rhs` is evaluated
    /// now. Up-rules are tried before anything else, so they extend the
    /// builtins, which are protected.
    ///
    /// ```
    /// use cas_core::{Context, Expr, Symbol};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let mut ctx = Context::new();
    /// let phi = Symbol::new("phi");
    /// ctx.define_up(phi, &read("Power[phi, 2]"), &read("Plus[phi, 1]")).unwrap();
    /// assert_eq!(ctx.eval(&read("Times[phi, phi]")).unwrap(), read("Plus[1, phi]"));
    /// ```
    ///
    /// # Errors
    ///
    /// [`DefineError::NotATag`], [`DefineError::HasValue`],
    /// [`DefineError::Protected`] for a protected tag, or the errors of
    /// [`Context::define`].
    pub fn define_up(&mut self, tag: Symbol, lhs: &Expr, rhs: &Expr) -> Result<(), DefineError> {
        run(self.limits(), |recursion| {
            self.up_definition(recursion, tag, lhs, rhs, true)
        })
    }

    /// `tag /: lhs := rhs` (`TagSetDelayed`, and `UpSetDelayed` for one
    /// tag): as [`Context::define_up`], but `rhs` is evaluated at each
    /// application.
    ///
    /// # Errors
    ///
    /// As [`Context::define_up`].
    pub fn define_up_delayed(
        &mut self,
        tag: Symbol,
        lhs: &Expr,
        rhs: &Expr,
    ) -> Result<(), DefineError> {
        run(self.limits(), |recursion| {
            self.up_definition(recursion, tag, lhs, rhs, false)
        })
    }

    /// Removes the value, the rules and the native function of `s`; its
    /// attributes and default stay.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected.
    pub fn clear(&mut self, s: Symbol) -> Result<(), ContextError> {
        self.unset(s)?;
        let Some(own) = self.rules().symbols.get(s) else {
            return Ok(());
        };
        let had_up = !own.up.rules.is_empty();
        let default = own.default.clone();
        let rules = self.rules_mut();
        match default {
            Some(default) => rules.symbols.insert(
                s,
                SymbolRules {
                    default: Some(default),
                    ..SymbolRules::default()
                },
            ),
            None => {
                rules.symbols.remove(s);
            }
        }
        if had_up {
            rules.up_count -= 1;
            if rules.up_count == 0 {
                rules.up_bits = 0;
            }
        }
        Ok(())
    }

    /// Gives `s` a native function, tried after the up-rules and before its
    /// rules (CAS-DOC-4 §6), in place of any it had.
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use cas_core::{Context, Expr, Symbol};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let mut ctx = Context::new();
    /// // `count[e]` is the number of arguments of `e`.
    /// ctx.define_native(Symbol::new("count"), Arc::new(|_, e| {
    ///     let [arg] = e.args() else { return Ok(None) };
    ///     Ok(Some(Expr::from(i64::try_from(arg.args().len()).unwrap())))
    /// }))
    /// .unwrap();
    /// assert_eq!(ctx.eval(&read("count[f[a, b, c]]")).unwrap(), Expr::from(3));
    /// ```
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected.
    pub fn define_native(&mut self, s: Symbol, native: UserNative) -> Result<(), ContextError> {
        self.check_unprotected(s)?;
        self.rules_mut()
            .symbols
            .entry(s, SymbolRules::default)
            .native = Some(native);
        Ok(())
    }

    /// Sets the value `x_.` takes as a missing argument of `s` (CAS-DOC-4
    /// §4), and recompiles the rules that mention `s`.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected.
    pub fn set_default(&mut self, s: Symbol, value: Expr) -> Result<(), ContextError> {
        self.check_unprotected(s)?;
        self.rules_mut()
            .symbols
            .entry(s, SymbolRules::default)
            .default = Some(value);
        self.recompile(s);
        Ok(())
    }

    /// A definition with `recursion` levels to evaluate in, its right-hand
    /// side evaluated `now` or kept.
    fn definition(
        &mut self,
        recursion: u32,
        lhs: &Expr,
        rhs: &Expr,
        now: bool,
    ) -> Result<(), DefineError> {
        let (owner, kind) = attachment(lhs)?;
        self.check_unprotected(owner)?;
        let (lhs, rhs) = {
            let mut ev = Eval::new(self, recursion);
            let rhs = if now { ev.eval(rhs)? } else { rhs.clone() };
            let lhs = match kind {
                Kind::Value => lhs.clone(),
                Kind::Down | Kind::Sub => normalised(&mut ev, lhs)?,
            };
            (lhs, rhs)
        };
        self.place(owner, kind, lhs, rhs)
    }

    /// Stores a definition of `owner`, its left-hand side normalised.
    fn place(
        &mut self,
        owner: Symbol,
        kind: Kind,
        lhs: Expr,
        rhs: Expr,
    ) -> Result<(), DefineError> {
        match kind {
            Kind::Value => self.set(owner, rhs)?,
            Kind::Down if !lhs.has_pattern() => {
                let depends = user_heads(&lhs);
                let rules = self.rules_mut();
                rules.depend(&depends, owner);
                rules
                    .symbols
                    .entry(owner, SymbolRules::default)
                    .insert_exact(lhs, rhs);
            }
            Kind::Down | Kind::Sub => {
                let now = self.rules().symbols.get(owner).map_or(0, |r| r.next);
                let stored = self.compile(lhs, rhs, now)?;
                let rules = self.rules_mut();
                rules.depend(&stored.depends, owner);
                let own = rules.symbols.entry(owner, SymbolRules::default);
                own.tick();
                if kind == Kind::Down {
                    own.down.put(stored);
                } else {
                    own.sub.put(stored);
                }
            }
        }
        Ok(())
    }

    fn up_definition(
        &mut self,
        recursion: u32,
        tag: Symbol,
        lhs: &Expr,
        rhs: &Expr,
        now: bool,
    ) -> Result<(), DefineError> {
        self.check_unprotected(tag)?;
        if self.value(tag).is_some() {
            return Err(DefineError::HasValue(tag));
        }
        let (lhs, rhs) = {
            let mut ev = Eval::new(self, recursion);
            let rhs = if now { ev.eval(rhs)? } else { rhs.clone() };
            (normalised(&mut ev, lhs)?, rhs)
        };
        if !core(&lhs).args().iter().any(|a| self::tag(a) == Some(tag)) {
            return Err(DefineError::NotATag(tag));
        }
        let now = self.rules().symbols.get(tag).map_or(0, |r| r.next);
        let stored = self.compile(lhs, rhs, now)?;
        let rules = self.rules_mut();
        rules.depend(&stored.depends, tag);
        let own = rules.symbols.entry(tag, SymbolRules::default);
        own.tick();
        let new_tag = own.up.rules.is_empty();
        own.up.put(stored);
        if new_tag {
            rules.up_count += 1;
            rules.up_bits |= symbol_bit(tag);
        }
        Ok(())
    }

    /// The rule `lhs :> rhs`, defined at time `written`, compiled with the
    /// defaults of this context.
    fn compile(&self, lhs: Expr, rhs: Expr, written: u64) -> Result<Stored, PatternError> {
        let defaults = |s: Symbol| self.rules().symbols.get(s).and_then(|r| r.default.clone());
        let pattern = Pattern::with_defaults(&lhs, &defaults)?;
        Ok(Stored {
            key: key(&lhs, &rhs, pattern.variables()),
            filter: pattern.filter(),
            // Literal parts too: their normal form depends on their heads.
            depends: user_heads(&lhs),
            rule: Rule::from_pattern(pattern, rhs),
            lhs,
            written,
        })
    }

    /// Recompiles the rules that depend on the attributes or the default of
    /// `s`, under those it has now.
    pub(crate) fn recompile(&mut self, s: Symbol) {
        let Some(owners) = self.rules().dependents.get(s) else {
            return;
        };
        let updated: Vec<(Symbol, SymbolRules)> = owners
            .iter()
            .filter_map(|&owner| {
                let own = self.rules().symbols.get(owner)?;
                Some((owner, own.recompiled(self, s)))
            })
            .collect();
        let rules = self.rules_mut();
        for (owner, own) in updated {
            for stored in own
                .down
                .rules
                .iter()
                .chain(&own.up.rules)
                .chain(&own.sub.rules)
            {
                rules.depend(&stored.depends, owner);
            }
            rules.symbols.insert(owner, own);
        }
    }
}

#[cfg(test)]
mod tests;
