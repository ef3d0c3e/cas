//! The evaluation context: everything semantic about expressions.

use std::fmt;
use std::num::NonZeroU64;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::attributes::Attributes;
use crate::expr::{Expr, Repr};
use crate::rules::Rules;
use crate::symbol::Symbol;
use crate::symbol_map::SymbolMap;

/// Values, rules, the attributes of user symbols, assumptions and limits.
///
/// Expressions do not depend on a context; evaluation always does. A
/// `Context` is a value: cloning it is O(1) and shares the state, mutation
/// copies whatever is shared (`&mut self`, copy on write), and a scope is a
/// new context ([`Context::assuming`]), so a child's changes never reach its
/// parent.
///
/// ```
/// use cas_core::{Attributes, Context, Expr, Symbol};
///
/// let (f, x) = (Symbol::new("f"), Symbol::new("x"));
/// let mut ctx = Context::new();
/// ctx.set_attributes(f, Attributes::ORDERLESS).unwrap();
/// assert!(ctx.set_attributes(Symbol::PLUS, Attributes::NONE).is_err());
///
/// let positive = Expr::apply(&Symbol::new("Greater").into(), &[x.into(), 0.into()]);
/// let child = ctx.assuming(positive.clone());
/// assert_eq!(child.assumptions(), &[positive]);
/// assert!(ctx.assumptions().is_empty());
/// assert_ne!(child.id(), ctx.id());
/// ```
#[derive(Clone)]
pub struct Context {
    state: Arc<State>,
    id: StateId,
}

/// The maps are persistent, so a change to a shared state copies O(log n)
/// of one map, never a whole map: a session that keeps each state it passes
/// through loads n definitions in O(n log n).
#[derive(Clone, Default)]
struct State {
    values: SymbolMap<Expr>,
    /// User symbols only; builtin attributes are fixed.
    attributes: SymbolMap<Attributes>,
    rules: Rules,
    assumptions: Arc<Vec<Expr>>,
    limits: Limits,
}

/// Identifies one context state. Every new state (a new context, each
/// mutation, each scope) takes the next id of one global counter, so ids are
/// unique across all contexts and threads, and two contexts have the same id
/// iff one is an unmodified clone of the other.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct StateId(NonZeroU64);

impl StateId {
    fn fresh() -> StateId {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        // 2^64 states are never reached; the check catches a wrap to zero.
        let id = NEXT.fetch_add(1, Ordering::Relaxed);
        StateId(NonZeroU64::new(id).expect("StateId overflow"))
    }

    fn get(self) -> u64 {
        self.0.get()
    }
}

impl Default for Context {
    fn default() -> Context {
        Context {
            state: Arc::default(),
            id: StateId::fresh(),
        }
    }
}

/// Evaluation limits (CAS-DOC-4 §2.1).
///
/// The defaults suit [`Context::eval`] on the caller's thread with a 2 MiB
/// stack, the least a Rust thread gets by default. [`Limits::deep`] runs
/// evaluation on a thread of its own with a large stack.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub struct Limits {
    /// Nested evaluations allowed: each level of an expression takes one, and
    /// so does each evaluation a native function or a rule starts. Reaching
    /// it is [`EvalError::RecursionLimit`], before the stack runs out as long
    /// as the stack has [`STACK_PER_LEVEL`](crate::STACK_PER_LEVEL) bytes per
    /// level.
    pub recursion: u32,
    /// Rewrites of one expression before giving up.
    pub iterations: u32,
    /// Steps of one match attempt, one pattern against one subject: each
    /// node tried and each way of splitting arguments counts one. Running
    /// out is [`EvalError::MatchLimit`], never "no match" (CAS-DOC-4 §3.3,
    /// R8).
    pub match_steps: u64,
    /// The stack of a thread of its own to evaluate on, in bytes; `0`
    /// evaluates on the caller's thread. On its own thread, evaluation stops
    /// at `stack / STACK_PER_LEVEL` levels, whatever `recursion` says.
    pub stack: usize,
}

impl Default for Limits {
    fn default() -> Limits {
        Limits {
            recursion: 512,
            iterations: 4096,
            match_steps: 1_000_000,
            stack: 0,
        }
    }
}

/// A rejected change to a context.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum ContextError {
    /// The symbol is `PROTECTED`: a builtin, or a user symbol given that
    /// attribute.
    Protected(Symbol),
    /// A value for a symbol with up-rules: it would replace the symbol before
    /// they could apply (CAS-DOC-4 §5.2, O6).
    HasUpRules(Symbol),
}

impl fmt::Display for ContextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ContextError::Protected(s) => write!(f, "symbol {s:?} is protected"),
            ContextError::HasUpRules(s) => {
                write!(f, "symbol {s:?} has up-rules, which a value would hide")
            }
        }
    }
}

impl std::error::Error for ContextError {}

/// Why evaluation stopped without a result.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum EvalError {
    /// More nested evaluations than [`Limits::recursion`].
    RecursionLimit,
    /// More rewrites of one expression than [`Limits::iterations`].
    IterationLimit,
    /// Lists of different shapes: lists of different lengths under a
    /// `Listable` head, or a list with something that is not a list
    /// (CAS-DOC-4 §2, O8).
    Shape,
    /// More steps in one match attempt than [`Limits::match_steps`].
    MatchLimit,
}

impl fmt::Display for EvalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            EvalError::RecursionLimit => "recursion limit exceeded",
            EvalError::IterationLimit => "iteration limit exceeded",
            EvalError::Shape => "lists of different shapes",
            EvalError::MatchLimit => "match step limit exceeded",
        })
    }
}

impl std::error::Error for EvalError {}

impl Context {
    /// A context with only the builtins.
    #[must_use]
    pub fn new() -> Context {
        Context::default()
    }

    /// The id of this context's state.
    #[must_use]
    pub fn id(&self) -> StateId {
        self.id
    }

    /// The state, for a change: a new id, and a private copy of whatever is
    /// shared.
    fn state_mut(&mut self) -> &mut State {
        self.id = StateId::fresh();
        Arc::make_mut(&mut self.state)
    }

    /// The fixed attributes of a builtin, else those set in this context.
    #[must_use]
    pub fn attributes(&self, s: Symbol) -> Attributes {
        s.builtin_attributes()
            .or_else(|| self.state.attributes.get(s).copied())
            .unwrap_or(Attributes::NONE)
    }

    /// Replaces the attributes of `s`, and recompiles the rules that depend
    /// on them (CAS-DOC-4 §4): those of `s`, and those whose patterns mention
    /// `s`.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected: always for builtins,
    /// whose attributes are fixed.
    pub fn set_attributes(
        &mut self,
        s: Symbol,
        attributes: Attributes,
    ) -> Result<(), ContextError> {
        self.check_unprotected(s)?;
        let before = self.attributes(s);
        let map = &mut self.state_mut().attributes;
        if attributes.is_empty() {
            map.remove(s);
        } else {
            map.insert(s, attributes);
        }
        // Protection changes nothing in matching.
        if before.difference(Attributes::PROTECTED) != attributes.difference(Attributes::PROTECTED)
        {
            self.recompile(s);
        }
        Ok(())
    }

    /// Removes `PROTECTED` from a user symbol; nothing else does.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is a builtin.
    pub fn unprotect(&mut self, s: Symbol) -> Result<(), ContextError> {
        if s.is_builtin() {
            return Err(ContextError::Protected(s));
        }
        let a = self.attributes(s);
        if a.contains(Attributes::PROTECTED) {
            let rest = a.difference(Attributes::PROTECTED);
            let map = &mut self.state_mut().attributes;
            if rest.is_empty() {
                map.remove(s);
            } else {
                map.insert(s, rest);
            }
        }
        Ok(())
    }

    /// The value of `s`, if set.
    #[must_use]
    pub fn value(&self, s: Symbol) -> Option<&Expr> {
        self.state.values.get(s)
    }

    /// Sets the value of `s`.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected, and
    /// [`ContextError::HasUpRules`] if it has up-rules.
    pub fn set(&mut self, s: Symbol, value: Expr) -> Result<(), ContextError> {
        self.check_unprotected(s)?;
        if self.rules().has_up_rules(s) {
            return Err(ContextError::HasUpRules(s));
        }
        self.state_mut().values.insert(s, value);
        Ok(())
    }

    /// Removes the value of `s`, if any.
    ///
    /// # Errors
    ///
    /// [`ContextError::Protected`] if `s` is protected.
    pub fn unset(&mut self, s: Symbol) -> Result<(), ContextError> {
        self.check_unprotected(s)?;
        if self.state.values.get(s).is_some() {
            self.state_mut().values.remove(s);
        }
        Ok(())
    }

    /// The rules of this state.
    pub(crate) fn rules(&self) -> &Rules {
        &self.state.rules
    }

    /// The rules, for a change: a new state.
    pub(crate) fn rules_mut(&mut self) -> &mut Rules {
        &mut self.state_mut().rules
    }

    pub(crate) fn check_unprotected(&self, s: Symbol) -> Result<(), ContextError> {
        if self.attributes(s).contains(Attributes::PROTECTED) {
            Err(ContextError::Protected(s))
        } else {
            Ok(())
        }
    }

    /// A child scope that also assumes `assumption`. `self` is unchanged;
    /// the child has its own id.
    #[must_use]
    pub fn assuming(&self, assumption: Expr) -> Context {
        let mut child = self.clone();
        Arc::make_mut(&mut child.state_mut().assumptions).push(assumption);
        child
    }

    /// The assumptions in scope, outermost first.
    #[must_use]
    pub fn assumptions(&self) -> &[Expr] {
        &self.state.assumptions
    }

    #[must_use]
    pub fn limits(&self) -> Limits {
        self.state.limits
    }

    pub fn set_limits(&mut self, limits: Limits) {
        self.state_mut().limits = limits;
    }

    /// Whether `e` was marked a fixpoint in this state
    /// ([`mark_fixpoint`](Self::mark_fixpoint)). Always `false` for atoms,
    /// which have no stamp; the evaluator handles them directly.
    ///
    /// A node holds one stamp, so marking it in another state forgets this
    /// one: a miss, never a wrong hit, since ids are unique.
    #[must_use]
    pub fn is_fixpoint(&self, e: &Expr) -> bool {
        match &e.0 {
            Repr::Node(n) => n.eval_stamp.load(Ordering::Relaxed) == self.id.get(),
            _ => false,
        }
    }

    /// Records that evaluating `e` in this state gives `e` itself. For
    /// evaluators: call it only for a true fixpoint (no limit hit, nothing
    /// aborted). No-op for atoms.
    pub fn mark_fixpoint(&self, e: &Expr) {
        // Relaxed: the node is immutable, so the stamp publishes nothing,
        // and a stale read only misses.
        if let Repr::Node(n) = &e.0 {
            n.eval_stamp.store(self.id.get(), Ordering::Relaxed);
        }
    }
}

impl fmt::Debug for Context {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Context")
            .field("id", &self.id)
            .field("values", &self.state.values.len())
            .field("attributes", &self.state.attributes.len())
            .field("rules", &self.state.rules)
            .field("assumptions", &self.state.assumptions)
            .field("limits", &self.state.limits)
            .finish()
    }
}

#[cfg(test)]
mod tests;
