//! Evaluation (CAS-DOC-4 §2): the loop behind [`Context::eval`].
//!
//! `eval(e)`:
//!
//! 1. A number or a string is itself. A symbol with a value evaluates to its
//!    value (`I` has the value `i`), and any other symbol is itself.
//! 2. A node stamped as a fixpoint of this state, or met before in this call
//!    (the memo), is done at once.
//! 3. The head is evaluated, and its attributes rule the rest.
//! 4. The arguments are evaluated left to right, except the held ones.
//! 5. Normalisation: `Sequence` arguments are spliced, nested applications
//!    of a `Flat` head flattened, a `Listable` head threaded over lists and
//!    the arguments of an `Orderless` head sorted. In `Times`, lists keep
//!    their order, after the other factors, since `*` does not commute on
//!    them (§6.2).
//! 6. Rewrite, by the first source that applies (R10): the up-rules of the
//!    arguments, the head's native function (the builtins' are in `native`,
//!    CAS-54, a user's in the context), the system rules or the down-rules
//!    of the head, exact ones first (`rules`, CAS-56), then for a compound
//!    head its sub-rules and the application of a `Function` (`pattern`,
//!    CAS-55).
//! 7. A change in step 5 or 6 goes back to step 1, counted against
//!    [`Limits::iterations`]. Otherwise the expression is a fixpoint: it is
//!    stamped, so evaluating it again in this state is O(1), and memoised
//!    for the rest of the call.
//!
//! The evaluator is recursive. Each nested evaluation counts against
//! [`Limits::recursion`], which stops it before the stack runs out as long
//! as each level has [`STACK_PER_LEVEL`] bytes of stack.
//!
//! [`Context::eval_with`] also reports warnings and, when asked, the steps
//! taken (§9, `trace`). The loop comes in two instances, with a tracer and
//! without one, which has no branch for it; the sources of rewrites test for
//! one once each. Recording returns before evaluation goes deeper, so a
//! level takes no more stack. A part whose evaluation warned is neither
//! stamped nor memoised: evaluating it again warns again.

use std::borrow::Cow;
use std::fmt;

use foldhash::fast::FixedState;
use hashbrown::HashMap;

use crate::attributes::Attributes;
use crate::build::flatten;
use crate::context::{Context, EvalError, Limits};
use crate::expr::{Expr, View};
use crate::rules;
use crate::symbol::Symbol;
use crate::trace::{
    EvalOptions, Origin, Report, RuleId, RuleKind, StepKind, TraceLevel, Tracer, Warning,
};

/// The stack one level of evaluation may take, in bytes. A test measures the
/// deepest path, in debug and release builds, and checks it stays below
/// this. So a thread with `recursion * STACK_PER_LEVEL` bytes of stack never
/// overflows, which is how [`Limits::deep`] is sized.
pub const STACK_PER_LEVEL: usize = 3072;

/// The stack a thread of its own keeps besides its levels: for the thread
/// itself, the start of evaluation and what a level deeper than the measured
/// paths takes.
const STACK_RESERVE: usize = 256 << 10;

/// A native function: an application rewritten in Rust (CAS-DOC-4 §6).
///
/// It receives the application after normalisation and returns its rewrite,
/// or `None` when it does not apply. It must be deterministic in the state
/// and its input and change nothing visible, since the result is stamped
/// and memoised. It may evaluate parts with [`Eval::eval`].
pub type NativeFn = fn(&mut Eval<'_>, &Expr) -> Result<Option<Expr>, EvalError>;

/// One call of [`Context::eval`]: what native functions receive to evaluate
/// parts, read the context, warn and report their steps.
///
/// It holds the call's memo, which maps each node evaluated so far to its
/// result, so a repeated subexpression is evaluated once per call, the
/// recursion depth, and what [`Context::eval_with`] records: the warnings
/// and the trace (§9). The local frames of procedures (the Language module)
/// will live here too.
pub struct Eval<'c> {
    context: &'c Context,
    recursion: u32,
    iterations: u32,
    depth: u32,
    memo: HashMap<Expr, Expr, FixedState>,
    /// Whether stamps and the memo are used: off only in the tests that
    /// check they change no result.
    caches: bool,
    /// What [`Context::eval_with`] reports, lent by it; [`Context::eval`]
    /// keeps nothing, and its call has nothing to drop for it.
    recording: Option<&'c mut Recording>,
    /// How many warnings were given, kept or not.
    warned: u32,
}

/// The warnings so far, each once, and the trace being recorded.
#[derive(Default)]
struct Recording {
    warnings: Vec<Warning>,
    tracer: Option<Tracer>,
}

impl Recording {
    /// The recording of the evaluation of `input` with a trace at `level`.
    fn new(level: TraceLevel, input: &Expr) -> Recording {
        Recording {
            warnings: Vec::new(),
            tracer: (level != TraceLevel::Off).then(|| Tracer::new(level, input)),
        }
    }

    /// The report, the trace closed.
    fn into_report(self) -> Report {
        Report {
            warnings: self.warnings,
            trace: self.tracer.map(Tracer::finish).unwrap_or_default(),
        }
    }
}

impl fmt::Debug for Eval<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Eval")
            .field("context", &self.context.id())
            .field("depth", &self.depth)
            .field("memo", &self.memo.len())
            .finish_non_exhaustive()
    }
}

/// What steps 3 to 6 make of an application.
enum Next {
    /// Something changed: evaluate this from step 1.
    Again(Expr),
    /// A fixpoint.
    Done(Expr),
}

impl<'c> Eval<'c> {
    /// A call in `context`, stopping at `recursion` nested evaluations.
    pub(crate) fn new(context: &'c Context, recursion: u32) -> Eval<'c> {
        Eval {
            context,
            recursion,
            iterations: context.limits().iterations,
            depth: 0,
            memo: HashMap::with_hasher(FixedState::default()),
            caches: true,
            recording: None,
            warned: 0,
        }
    }

    /// The context evaluation runs in.
    #[must_use]
    pub fn context(&self) -> &'c Context {
        self.context
    }

    /// Evaluates `e` as part of this call: one level deeper, with its memo.
    /// In a trace, the evaluation is a step of the native function or the
    /// condition that asked for it ([`StepKind::Evaluation`]).
    ///
    /// # Errors
    ///
    /// [`EvalError`] when a limit is reached or an operation has no result.
    pub fn eval(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        if self.tracing() {
            return self.eval_traced(e);
        }
        self.part::<false>(e)
    }

    /// [`Eval::eval`] in a trace, which records the evaluation as a step.
    #[inline(never)]
    fn eval_traced(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        self.trace_open(e, false);
        let result = self.part::<true>(e);
        self.trace_close_evaluation(&result);
        result
    }

    /// Evaluates a condition's test, as [`Eval::eval`] but in a trace only
    /// at [`TraceLevel::Verbose`].
    pub(crate) fn test(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        let Some(recording) = &mut self.recording else {
            return self.eval(e);
        };
        if recording.tracer.as_ref().map(Tracer::level) != Some(TraceLevel::Rules) {
            return self.eval(e);
        }
        let tracer = recording.tracer.take();
        let result = self.eval(e);
        if let Some(recording) = &mut self.recording {
            recording.tracer = tracer;
        }
        result
    }

    /// Evaluates the root of a call, or a part of the expression being
    /// evaluated: one level deeper, located where the caller is. The loop
    /// comes in two instances, with a tracer (`TRACED`) and without one,
    /// which has no branch for it; inlined, so it takes no frame of its own.
    #[allow(
        clippy::inline_always,
        reason = "no frame of its own on each level, in a debug build too"
    )]
    #[inline(always)]
    fn part<const TRACED: bool>(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        self.enter()?;
        let result = self.evaluate::<TRACED>(e);
        self.leave();
        result
    }

    /// Whether this evaluation records a trace, so that a native function
    /// may skip building what only a step would show.
    #[must_use]
    pub fn tracing(&self) -> bool {
        self.recording.as_ref().is_some_and(|r| r.tracer.is_some())
    }

    /// The trace being recorded, if any.
    fn tracer(&mut self) -> Option<&mut Tracer> {
        self.recording
            .as_deref_mut()
            .and_then(|r| r.tracer.as_mut())
    }

    /// Records a step a native function took itself (§9): `before` became
    /// `after`, as `note` says. It is a child of the native's own step;
    /// without a trace, nothing.
    pub fn step(&mut self, note: impl Into<Cow<'static, str>>, before: &Expr, after: &Expr) {
        if let Some(t) = self.tracer() {
            t.record(StepKind::Algorithm(note.into()), before, after, Vec::new());
        }
    }

    /// Gives `warning`, which the report of [`Context::eval_with`] lists
    /// once; [`Context::eval`] drops it. The part being evaluated is not
    /// stamped, so evaluating it again warns again.
    pub fn warn(&mut self, warning: Warning) {
        self.warned = self.warned.wrapping_add(1);
        if let Some(r) = &mut self.recording {
            if !r.warnings.contains(&warning) {
                r.warnings.push(warning);
            }
        }
    }

    /// What rewrote the subject of the attempt in progress, for the trace.
    pub(crate) fn trace_rewrote(&mut self, kind: impl FnOnce() -> StepKind) {
        if let Some(t) = self.tracer() {
            t.rewrote(kind());
        }
    }

    /// A source that did not rewrite the subject of the attempt in
    /// progress: a patterned rule tried (`counted`), or a native function.
    pub(crate) fn trace_failed(&mut self, counted: bool) {
        if let Some(t) = self.tracer() {
            t.failed(counted);
        }
    }

    #[inline(never)]
    fn trace_open(&mut self, subject: &Expr, attempt: bool) {
        if let Some(t) = self.tracer() {
            t.open(subject, attempt);
        }
    }

    #[inline(never)]
    fn trace_close_evaluation(&mut self, result: &Result<Expr, EvalError>) {
        if let Some(t) = self.tracer() {
            t.close_evaluation(result.as_ref().map_err(|_| ()));
        }
    }

    #[inline(never)]
    fn trace_close_attempt(&mut self, rewritten: &Result<Option<Expr>, EvalError>) {
        if let Some(t) = self.tracer() {
            t.close_attempt(rewritten.as_ref().map(Option::as_ref).map_err(|_| ()));
        }
    }

    #[inline(never)]
    fn trace_enter(&mut self, i: usize) {
        if let Some(t) = self.tracer() {
            t.enter(i);
        }
    }

    #[inline(never)]
    fn trace_leave(&mut self) {
        if let Some(t) = self.tracer() {
            t.leave();
        }
    }

    #[inline(never)]
    fn trace_value(&mut self, s: Symbol, value: &Expr) {
        let kind = if self.context.value(s).is_some() {
            StepKind::Rule(RuleId {
                symbol: s,
                kind: RuleKind::Value,
                index: 0,
                origin: Origin::User,
            })
        } else {
            StepKind::Native(s)
        };
        if let Some(t) = self.tracer() {
            t.record(kind, &s.into(), value, Vec::new());
        }
    }

    fn trace_next_id(&self) -> usize {
        self.recording
            .as_ref()
            .and_then(|r| r.tracer.as_ref())
            .map_or(0, Tracer::next_id)
    }

    #[inline(never)]
    fn trace_memo(&mut self, before: &Expr, after: &Expr) {
        if let Some(t) = self.tracer() {
            t.memo(before, after);
        }
    }

    #[inline(never)]
    fn trace_memoised(&mut self, passed: &[Expr], result: &Expr, first: usize) {
        if let Some(t) = self.tracer() {
            t.memoised(passed, result, first);
        }
    }

    #[inline(never)]
    fn trace_normalise(&mut self, head: &Expr, args: &[Expr], after: &Expr) {
        if let Some(t) = self.tracer() {
            if t.level() == TraceLevel::Verbose {
                t.record(
                    StepKind::Normalise,
                    &Expr::apply(head, args),
                    after,
                    Vec::new(),
                );
            }
        }
    }

    /// Takes one level of [`Limits::recursion`], for evaluation or for the
    /// recursion of a native, such as the matcher's; [`Eval::leave`] gives
    /// it back.
    pub(crate) fn enter(&mut self) -> Result<(), EvalError> {
        if self.depth >= self.recursion {
            return Err(EvalError::RecursionLimit);
        }
        self.depth += 1;
        Ok(())
    }

    pub(crate) fn leave(&mut self) {
        self.depth -= 1;
    }

    /// Steps 1 to 7.
    fn evaluate<const TRACED: bool>(&mut self, input: &Expr) -> Result<Expr, EvalError> {
        let first = if TRACED { self.trace_next_id() } else { 0 };
        let warned = self.warned;
        let mut e = input.clone();
        // The nodes passed on the way, memoised with the result.
        let mut passed = Vec::new();
        let mut rewrites = 0;
        let mut fixpoint = false;
        let result = loop {
            match e.view() {
                View::Symbol(s) => {
                    match self
                        .context
                        .value(s)
                        .cloned()
                        .or_else(|| crate::native::value(s))
                    {
                        Some(value) => {
                            if TRACED {
                                self.trace_value(s, &value);
                            }
                            e = value;
                            self.count(&mut rewrites)?;
                            continue;
                        }
                        None => break e,
                    }
                }
                View::Apply { .. } => {}
                _ => break e,
            }
            if self.caches {
                if self.context.is_fixpoint(&e) {
                    break e;
                }
                if let Some(result) = self.memo.get(&e) {
                    let result = result.clone();
                    if TRACED {
                        self.trace_memo(&e, &result);
                    }
                    break result;
                }
            }
            passed.push(e.clone());
            let (head, attributes, args) = self.parts::<TRACED>(&e)?;
            match self.normal::<TRACED>(&e, &head, attributes, &args)? {
                Next::Again(next) => {
                    e = next;
                    self.count(&mut rewrites)?;
                }
                Next::Done(result) => {
                    fixpoint = true;
                    break result;
                }
            }
        };
        if self.caches && self.warned == warned {
            if fixpoint {
                self.context.mark_fixpoint(&result);
            }
            if TRACED {
                self.trace_memoised(&passed, &result, first);
            }
            for p in passed {
                if p != result {
                    self.memo.insert(p, result.clone());
                }
            }
        }
        Ok(result)
    }

    fn count(&self, rewrites: &mut u32) -> Result<(), EvalError> {
        *rewrites += 1;
        if *rewrites > self.iterations {
            Err(EvalError::IterationLimit)
        } else {
            Ok(())
        }
    }

    /// Steps 3 and 4 on an application: the head evaluated, then the
    /// arguments it does not hold, located in a trace. Apart from steps 5
    /// and 6 ([`Eval::normal`]), so that in a debug build a level through
    /// the arguments has on its stack only what evaluating them takes, and a
    /// level through a rewrite only what rewriting takes.
    fn parts<const TRACED: bool>(
        &mut self,
        e: &Expr,
    ) -> Result<(Expr, Attributes, Vec<Expr>), EvalError> {
        let View::Apply { head, args } = e.view() else {
            unreachable!("only applications take steps");
        };
        if TRACED {
            self.trace_enter(0);
        }
        let new_head = self.part::<TRACED>(head);
        if TRACED {
            self.trace_leave();
        }
        let new_head = new_head?;
        let attributes = new_head
            .as_symbol()
            .map_or(Attributes::NONE, |s| self.context.attributes(s));
        let mut evaluated = Vec::with_capacity(args.len());
        for (i, a) in args.iter().enumerate() {
            evaluated.push(if holds(attributes, i) {
                a.clone()
            } else {
                if TRACED {
                    self.trace_enter(i + 1);
                }
                let value = self.part::<TRACED>(a);
                if TRACED {
                    self.trace_leave();
                }
                value?
            });
        }
        Ok((new_head, attributes, evaluated))
    }

    /// Steps 5 and 6 on `e`, whose head and arguments evaluated to `head`
    /// and `evaluated`: normalisation, then a rewrite.
    fn normal<const TRACED: bool>(
        &mut self,
        e: &Expr,
        head: &Expr,
        attributes: Attributes,
        evaluated: &[Expr],
    ) -> Result<Next, EvalError> {
        match normalise(head, attributes, evaluated)? {
            Normal::Same => {}
            Normal::Changed(args) => {
                let next = Expr::apply(head, &args);
                if TRACED {
                    self.trace_normalise(head, evaluated, &next);
                }
                return Ok(Next::Again(next));
            }
            Normal::Threaded(list) => {
                if TRACED {
                    self.trace_normalise(head, evaluated, &list);
                }
                return Ok(Next::Again(list));
            }
        }
        let e = if e.head() == Some(head) && evaluated == e.args() {
            e.clone()
        } else {
            Expr::apply(head, evaluated)
        };
        if TRACED {
            self.trace_open(&e, true);
        }
        let rewritten = self.rewrite(&e, head.as_symbol());
        if TRACED {
            self.trace_close_attempt(&rewritten);
        }
        Ok(match rewritten? {
            Some(rewritten) if rewritten != e => Next::Again(rewritten),
            _ => Next::Done(e),
        })
    }

    /// Step 6: the first source that rewrites `e`, an application with
    /// `head`. Inlined into the loop in release builds; apart in debug
    /// builds, whose frames would add up.
    #[allow(
        clippy::inline_always,
        reason = "inlined into the loop, which counting instructions measured (CAS-61)"
    )]
    #[cfg_attr(not(debug_assertions), inline(always))]
    fn rewrite(&mut self, e: &Expr, head: Option<Symbol>) -> Result<Option<Expr>, EvalError> {
        // 1. The up-rules of the arguments.
        if let Some(rewritten) = rules::up(self, e)? {
            return Ok(Some(rewritten));
        }
        if let Some(head) = head {
            // 2. A builtin's native function.
            if let Some(native) = native(head) {
                if let Some(rewritten) = native(self, e)? {
                    self.trace_rewrote(|| StepKind::Native(head));
                    return Ok(Some(rewritten));
                }
                self.trace_failed(false);
            }
            // 2 to 4. A user's native function, the system rules of a
            // builtin, the down-rules of a user's symbol.
            return rules::down(self, e, head);
        }
        // 5. Sub-rules, and a `Function` head applied.
        if let Some(rewritten) = rules::sub(self, e)? {
            return Ok(Some(rewritten));
        }
        if let View::Apply { head, args } = e.view() {
            if head.head().and_then(Expr::as_symbol) == Some(Symbol::FUNCTION) {
                let applied = crate::pattern::apply_function(head, args);
                if applied.is_some() {
                    self.trace_rewrote(|| StepKind::Native(Symbol::FUNCTION));
                }
                return Ok(applied);
            }
        }
        Ok(None)
    }
}

/// Whether argument `i` (from 0) of a head with `attributes` is held.
pub(crate) fn holds(attributes: Attributes, i: usize) -> bool {
    attributes.contains(if i == 0 {
        Attributes::HOLD_FIRST
    } else {
        Attributes::HOLD_REST
    })
}

/// The native function of a builtin (CAS-54), or one a test gave.
#[allow(
    clippy::inline_always,
    reason = "inlined into the loop, which counting instructions measured (CAS-61)"
)]
#[cfg_attr(not(debug_assertions), inline(always))]
fn native(s: Symbol) -> Option<NativeFn> {
    crate::native::get(s).or_else(|| test_native(s))
}

#[cfg(test)]
use tests::native as test_native;

#[cfg(not(test))]
fn test_native(_: Symbol) -> Option<NativeFn> {
    None
}

/// What step 5 makes of evaluated arguments.
pub(crate) enum Normal {
    Same,
    Changed(Vec<Expr>),
    /// A `Listable` head threaded over lists: the list of its applications.
    Threaded(Expr),
}

/// Step 5 on the evaluated arguments of `head`, whose attributes are
/// `attributes`. Inlined into the loop in release builds.
#[allow(
    clippy::inline_always,
    reason = "inlined into the loop, which counting instructions measured (CAS-61)"
)]
#[cfg_attr(not(debug_assertions), inline(always))]
pub(crate) fn normalise(
    head: &Expr,
    attributes: Attributes,
    args: &[Expr],
) -> Result<Normal, EvalError> {
    let sequence = Expr::from(Symbol::SEQUENCE);
    let mut args = Cow::Borrowed(args);
    if args.iter().any(|a| a.head() == Some(&sequence)) {
        args = Cow::Owned(flatten(&sequence, &args));
    }
    if attributes.contains(Attributes::FLAT) && args.iter().any(|a| a.head() == Some(head)) {
        args = Cow::Owned(flatten(head, &args));
    }
    if attributes.contains(Attributes::LISTABLE) && args.iter().any(is_list) {
        return thread(head, &args).map(Normal::Threaded);
    }
    if attributes.contains(Attributes::ORDERLESS) {
        if let Some(sorted) = sort_arguments(head, &args) {
            args = Cow::Owned(sorted);
        }
    }
    Ok(match args {
        Cow::Borrowed(_) => Normal::Same,
        Cow::Owned(args) => Normal::Changed(args),
    })
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// Whether `e` is a list, or a pattern that matches only lists: `_List`,
/// `x_List`, or a test `ListQ` or `MatrixQ`. In `Times`, those keep their
/// order, in patterns as in subjects, since `*` does not commute on them
/// (R13).
pub(crate) fn is_list_like(e: &Expr) -> bool {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::LIST), _) => true,
        (Some(Symbol::BLANK | Symbol::BLANK_SEQUENCE | Symbol::BLANK_NULL_SEQUENCE), [h]) => {
            h.as_symbol() == Some(Symbol::LIST)
        }
        (Some(Symbol::PATTERN), [_, p]) | (Some(Symbol::CONDITION), [p, _]) => is_list_like(p),
        (Some(Symbol::PATTERN_TEST), [p, test]) => {
            matches!(test.as_symbol(), Some(Symbol::LIST_Q | Symbol::MATRIX_Q)) || is_list_like(p)
        }
        _ => false,
    }
}

/// `head` threaded over `args`: `f[{a, b}, {c, d}]` is `{f[a, c], f[b, d]}`.
/// Every argument must be a list, all of one length: a scalar is not
/// broadcast (CAS-DOC-4 O8).
fn thread(head: &Expr, args: &[Expr]) -> Result<Expr, EvalError> {
    let len = args
        .iter()
        .find(|a| is_list(a))
        .map_or(0, |l| l.args().len());
    if !args.iter().all(|a| is_list(a) && a.args().len() == len) {
        return Err(EvalError::Shape);
    }
    let items: Vec<Expr> = (0..len)
        .map(|i| {
            let parts: Vec<Expr> = args.iter().map(|a| a.args()[i].clone()).collect();
            Expr::apply(head, &parts)
        })
        .collect();
    Ok(Expr::apply(&Symbol::LIST.into(), &items))
}

/// The arguments of an `Orderless` head in canonical order, if they are not
/// in it. In `Times`, lists, and patterns that match only lists, come after
/// the other factors, in their own order.
pub(crate) fn sort_arguments(head: &Expr, args: &[Expr]) -> Option<Vec<Expr>> {
    let sorted: Vec<Expr> =
        if head.as_symbol() == Some(Symbol::TIMES) && args.iter().any(is_list_like) {
            let (lists, mut factors): (Vec<Expr>, Vec<Expr>) =
                args.iter().cloned().partition(is_list_like);
            factors.sort();
            factors.extend(lists);
            factors
        } else {
            if args.windows(2).all(|w| w[0] <= w[1]) {
                return None;
            }
            let mut sorted = args.to_vec();
            sorted.sort();
            sorted
        };
    (sorted != args).then_some(sorted)
}

impl Context {
    /// Evaluates `e` (CAS-DOC-4 §2), on a thread of its own when the limits
    /// ask for one ([`Limits::stack`]).
    ///
    /// ```
    /// use cas_core::{Attributes, Context, Expr, Symbol};
    ///
    /// let mut ctx = Context::new();
    /// let (f, x) = (Symbol::new("f"), Symbol::new("x"));
    /// ctx.set_attributes(f, Attributes::ORDERLESS).unwrap();
    /// ctx.set(x, 2.into()).unwrap();
    /// let e = Expr::parse_full_form("f[y, x, Sequence[a, Hold[x]]]").unwrap();
    /// assert_eq!(format!("{:?}", ctx.eval(&e).unwrap()), "f[2, a, y, Hold[x]]");
    /// ```
    ///
    /// # Errors
    ///
    /// [`EvalError`] when a limit is reached or an operation has no result.
    ///
    /// # Panics
    ///
    /// If the evaluation thread cannot start.
    pub fn eval(&self, e: &Expr) -> Result<Expr, EvalError> {
        run(self.limits(), |recursion| {
            Eval::new(self, recursion).part::<false>(e)
        })
    }

    /// Evaluates `e` as [`Context::eval`] does, and reports what happened
    /// (CAS-DOC-4 §9): the warnings and, as `options` ask, the steps. The
    /// report comes also when evaluation fails, its trace ending with the
    /// steps in progress.
    ///
    /// ```
    /// use cas_core::{Context, EvalOptions, Expr, StepKind, Symbol, TraceLevel, Warning};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let mut ctx = Context::new();
    /// let lhs = read("f[PatternTest[Pattern[x, Blank[]], IntegerQ]]");
    /// ctx.define_delayed(&lhs, &read("Power[x, 2]")).unwrap();
    /// let options = EvalOptions { trace: TraceLevel::Rules };
    /// let (result, report) = ctx.eval_with(&read("Plus[f[3], f[a]]"), &options);
    /// assert_eq!(result.unwrap(), read("Plus[9, f[a]]"));
    /// assert_eq!(report.warnings, [Warning::NoDefinitionMatches(read("f[a]"))]);
    /// // `f[3]`, the first argument, became `3^2` by the rule of `f`, then 9.
    /// let [rule, power] = &report.trace.steps[..] else { panic!() };
    /// assert!(matches!(rule.kind, StepKind::Rule(id) if id.symbol == Symbol::new("f")));
    /// assert_eq!((&rule.position[..], &rule.after), (&[1][..], &read("Power[3, 2]")));
    /// assert_eq!((&power.kind, &power.after), (&StepKind::Native(Symbol::POWER), &read("9")));
    /// ```
    ///
    /// # Panics
    ///
    /// If the evaluation thread cannot start.
    pub fn eval_with(&self, e: &Expr, options: &EvalOptions) -> (Result<Expr, EvalError>, Report) {
        run(self.limits(), |recursion| {
            let mut recording = Recording::new(options.trace, e);
            let traced = recording.tracer.is_some();
            let mut ev = Eval::new(self, recursion);
            ev.recording = Some(&mut recording);
            let result = if traced {
                ev.part::<true>(e)
            } else {
                ev.part::<false>(e)
            };
            drop(ev);
            (result, recording.into_report())
        })
    }

    /// [`Context::eval`] without fixpoint stamps or memo: the same results,
    /// which tests check.
    #[cfg(test)]
    pub(crate) fn eval_uncached(&self, e: &Expr) -> Result<Expr, EvalError> {
        self.eval_uncached_with(e, TraceLevel::Off).0
    }

    /// The report of `dispatch` run on `e` as one attempt to rewrite it,
    /// traced at `level`: a test of what a dispatch site records.
    #[cfg(test)]
    pub(crate) fn traced_attempt(
        &self,
        e: &Expr,
        level: TraceLevel,
        dispatch: impl FnOnce(&mut Eval<'_>) -> Result<Option<Expr>, EvalError>,
    ) -> Report {
        let mut recording = Recording::new(level, e);
        let mut ev = Eval::new(self, self.limits().recursion);
        ev.recording = Some(&mut recording);
        ev.trace_open(e, true);
        let rewritten = dispatch(&mut ev);
        ev.trace_close_attempt(&rewritten);
        drop(ev);
        recording.into_report()
    }

    /// [`Context::eval_with`] without fixpoint stamps or memo, on the
    /// caller's thread.
    #[cfg(test)]
    pub(crate) fn eval_uncached_with(
        &self,
        e: &Expr,
        level: TraceLevel,
    ) -> (Result<Expr, EvalError>, Report) {
        let mut recording = Recording::new(level, e);
        let traced = recording.tracer.is_some();
        let mut ev = Eval::new(self, self.limits().recursion);
        ev.caches = false;
        ev.recording = Some(&mut recording);
        let result = if traced {
            ev.part::<true>(e)
        } else {
            ev.part::<false>(e)
        };
        drop(ev);
        (result, recording.into_report())
    }
}

/// Runs `f`, which evaluates, with the recursion limit of `limits`: on a
/// thread of its own when they ask for one, whose stack caps the limit.
///
/// # Panics
///
/// If the thread cannot start.
pub(crate) fn run<R: Send>(limits: Limits, f: impl FnOnce(u32) -> R + Send) -> R {
    if limits.stack == 0 {
        return f(limits.recursion);
    }
    let recursion = limits.recursion.min(levels(limits.stack));
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .name("cas-eval".into())
            .stack_size(limits.stack)
            .spawn_scoped(scope, move || f(recursion))
            .expect("the evaluation thread starts")
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
    })
}

/// The levels of evaluation a thread with `stack` bytes of stack holds.
fn levels(stack: usize) -> u32 {
    u32::try_from(stack.saturating_sub(STACK_RESERVE) / STACK_PER_LEVEL).unwrap_or(u32::MAX)
}

impl Limits {
    /// The limits for deep evaluation (CAS-DOC-4 §2.1), which the session
    /// uses: a thread of its own with 256 MiB of stack, and the recursion
    /// that stack holds, about 8.7·10⁴ levels.
    #[must_use]
    pub fn deep() -> Limits {
        let stack = 256 << 20;
        Limits {
            recursion: levels(stack),
            stack,
            ..Limits::default()
        }
    }
}

#[cfg(test)]
pub(crate) mod tests;
