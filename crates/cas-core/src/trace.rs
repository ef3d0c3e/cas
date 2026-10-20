//! The evaluation report (CAS-DOC-4 §9, CAS-61): the warnings an evaluation
//! gives, and the steps it takes, as a tree.
//!
//! [`Context::eval_with`](crate::Context::eval_with) evaluates with
//! [`EvalOptions`] and returns a [`Report`] beside the result, also when
//! evaluation fails. Plain [`Context::eval`](crate::Context::eval) runs the
//! instance of the evaluation loop without a tracer, which has no branch for
//! one.
//!
//! A [`Step`] says that `before` became `after` at a position, and what did
//! it, its [`StepKind`]: a rule, known by its [`RuleId`]; a native function;
//! a step a native reports itself ([`Eval::step`](crate::Eval::step)); the
//! memo; normalisation. Its children are the evaluations it asked for, a
//! native's or a condition's, each a [`StepKind::Evaluation`] whose own
//! children are located from what it evaluated. So a trace is a tree, kept
//! as a list of the steps in the order taken, which name their parent and
//! children by id: a trace as deep as evaluation goes is dropped, cloned and
//! printed without recursion. At [`TraceLevel::Verbose`] a trace is
//! complete: its top steps, replayed from the input, give the result.
//!
//! [`Trace::verify`] checks the steps, each of which claims `before =
//! after`, at exact sample points that satisfy the assumptions.

use std::borrow::Cow;
use std::fmt;
use std::panic::Location;

use foldhash::fast::FixedState;
use hashbrown::HashMap;

use crate::expr::Expr;
use crate::symbol::Symbol;

mod verify;

pub use verify::{Counterexample, Verification};

/// How much [`Context::eval_with`](crate::Context::eval_with) records.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TraceLevel {
    /// No steps: the report has the warnings only.
    #[default]
    Off,
    /// Rewrites, by rules and native functions, and the steps natives
    /// report.
    Rules,
    /// Also normalisations, the tests of conditions, and the rules tried
    /// that did not apply.
    Verbose,
}

/// What [`Context::eval_with`](crate::Context::eval_with) records.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct EvalOptions {
    /// How much of the steps.
    pub trace: TraceLevel,
}

/// What an evaluation reports beside its result.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    /// What went wrong without stopping evaluation, in the order it
    /// happened, each once.
    pub warnings: Vec<Warning>,
    /// The steps, none at [`TraceLevel::Off`].
    pub trace: Trace,
}

/// Something wrong that evaluation went past (CAS-DOC-4 O2).
///
/// An expression whose evaluation gave a warning is not stamped as a
/// fixpoint, so evaluating it again gives the warning again. Its message
/// writes out an expression of at most 1000 parts, and of a larger one its
/// head: an expression that shares its parts can have far more parts than
/// it takes memory.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Warning {
    /// A zero to a negative power, as `1/0`, which is `ComplexInfinity`.
    Infinite(Expr),
    /// An indeterminate form, `0^0`, `0·Infinity`, `Infinity - Infinity`
    /// or `Infinity^0`, which is `Indeterminate`.
    Indeterminate(Expr),
    /// A user's symbol with down-rules applied where none of them matches
    /// (CAS-DOC-4 §5.3): `f[1, 2]` for `f[x_] := …`.
    NoDefinitionMatches(Expr),
    /// A number applied as a function (CAS-DOC-7 §3.1), as `a(x + 1)` after
    /// `a := 2`, which stays as it is.
    NumberApplied(Expr),
    /// A native function's own: its message, and the expression.
    Native {
        /// What went wrong.
        message: Cow<'static, str>,
        /// Where.
        expr: Expr,
    },
}

/// Messages write out an expression of at most this many parts: a larger one
/// could take long, sharing parts that its text repeats.
const MAX_SHOWN: u32 = 1000;

/// An expression in a message: written out when small, else its head.
struct Shown<'e>(&'e Expr);

impl fmt::Display for Shown<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let e = self.0;
        if e.size() <= MAX_SHOWN {
            return write!(f, "{e}");
        }
        match e.head() {
            Some(head) if head.size() <= MAX_SHOWN => write!(f, "{head}[…]"),
            _ => f.write_str("…"),
        }
    }
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Warning::Infinite(e) => write!(f, "infinite expression {} encountered", Shown(e)),
            Warning::Indeterminate(e) => {
                write!(f, "indeterminate expression {} encountered", Shown(e))
            }
            Warning::NoDefinitionMatches(e) => match e.head() {
                Some(head) => write!(f, "no definition of {} matches {}", Shown(head), Shown(e)),
                None => write!(f, "no definition matches {}", Shown(e)),
            },
            Warning::NumberApplied(e) => match e.head() {
                Some(head) => write!(
                    f,
                    "the number {} applied as a function in {}",
                    Shown(head),
                    Shown(e)
                ),
                None => write!(f, "a number applied as a function in {}", Shown(e)),
            },
            Warning::Native { message, expr } => write!(f, "{message}: {}", Shown(expr)),
        }
    }
}

/// The steps of an evaluation, located from its input.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trace {
    /// Every step, children included, in the order taken: `steps[i].id` is
    /// `i`, and a step comes after its children.
    pub steps: Vec<Step>,
}

impl Trace {
    /// The steps of the evaluation itself, children of no other, in the
    /// order taken. At [`TraceLevel::Verbose`], replayed from the input,
    /// they give the result.
    pub fn top(&self) -> impl Iterator<Item = &Step> {
        self.steps.iter().filter(|s| s.parent.is_none())
    }

    /// The children of `step`, in the order taken.
    pub fn children<'t>(&'t self, step: &'t Step) -> impl Iterator<Item = &'t Step> {
        step.children.iter().map(|&c| &self.steps[c])
    }
}

/// One step of a [`Trace`]: `before` became `after` at `position`.
#[derive(Clone, Debug, PartialEq)]
pub struct Step {
    /// Its place in [`Trace::steps`], the order taken.
    pub id: usize,
    /// What did it.
    pub kind: StepKind,
    /// The path to the part rewritten, in the `Part` convention: 0 is the
    /// head, `i` the `i`-th argument, and the empty path the whole. Children
    /// are located from the expression their evaluation started at; a step
    /// a native reports, and an evaluation, have the empty path.
    pub position: Vec<usize>,
    /// The part before the step.
    pub before: Expr,
    /// The part after it.
    pub after: Expr,
    /// The step this is a child of, if any.
    pub parent: Option<usize>,
    /// The evaluations the step asked for, as [`StepKind::Evaluation`]
    /// steps, and the steps a native reported, in the order taken.
    pub children: Vec<usize>,
}

/// What took a [`Step`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum StepKind {
    /// A rule rewrote the part, or a symbol's value replaced it.
    Rule(RuleId),
    /// A native function did: a builtin's, a user's, or a `Function`
    /// applied (`Symbol::FUNCTION`), or the value of a builtin symbol.
    Native(Symbol),
    /// A step a native function reported itself, with its note.
    Algorithm(Cow<'static, str>),
    /// A part met before in this evaluation, which evaluated to `after`:
    /// the id of the first step that evaluation took, if the trace has it.
    Memo {
        /// The first step of the earlier evaluation.
        first: Option<usize>,
    },
    /// Normalisation: `Sequence` spliced, a `Flat` head flattened, a
    /// `Listable` head threaded, `Orderless` arguments sorted. Recorded at
    /// [`TraceLevel::Verbose`].
    Normalise,
    /// An evaluation that a native function or a condition asked for:
    /// `before` evaluated to `after`, and the children are its steps.
    /// Conditions are recorded at [`TraceLevel::Verbose`].
    Evaluation,
    /// Sources tried on `before` that did not rewrite it, recorded at
    /// [`TraceLevel::Verbose`], and after an error for the attempt in
    /// progress: how many patterned rules among them, and as children the
    /// evaluations they asked for. `after` is `before`.
    Attempts(u32),
}

/// Which rule a [`Step`] applied (CAS-DOC-4 §9).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RuleId {
    /// The symbol the rule is stored with: the head of a down-rule or an
    /// exact one, the tag of an up-rule, the innermost head of a sub-rule,
    /// the symbol of a value.
    pub symbol: Symbol,
    /// What kind of rule.
    pub kind: RuleKind,
    /// Its place in definition order: for a patterned rule, among those of
    /// its kind on `symbol` in the order they are tried, from 0; for an
    /// exact one, the time of its first definition among those of `symbol`;
    /// 0 for a value.
    pub index: u64,
    /// Where it comes from.
    pub origin: Origin,
}

/// What kind of rule a [`RuleId`] names (CAS-DOC-4 §5).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RuleKind {
    /// The value of a symbol.
    Value,
    /// A ground left-hand side, from the exact table.
    Exact,
    /// A patterned rule of the head.
    Down,
    /// A rule an argument's symbol carries.
    Up,
    /// A rule for `f[…][…]`.
    Sub,
}

/// Where a rule comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Origin {
    /// A user's definition.
    User,
    /// A system rule (CAS-DOC-4 §5.4): its stable name, such as
    /// `"D/product"`, and where in the sources it is written.
    System {
        /// The name.
        name: &'static str,
        /// Where it is written.
        location: &'static Location<'static>,
    },
}

/// The trace an evaluation records, for
/// [`Context::eval_with`](crate::Context::eval_with).
#[derive(Debug)]
pub(crate) struct Tracer {
    level: TraceLevel,
    /// Every step recorded, by id, and whether it was dropped: the steps of
    /// a source that did not rewrite, below the Verbose level.
    steps: Vec<Step>,
    dropped: Vec<bool>,
    /// The steps being collected: the root's, then those of each attempt to
    /// rewrite and each evaluation asked for, the innermost last.
    frames: Vec<Frame>,
    /// For each part memoised, the id of the first step its evaluation took.
    first: HashMap<Expr, usize, FixedState>,
    /// What rewrote the subject of the innermost attempt: set where it did,
    /// taken when the attempt closes.
    source: Option<StepKind>,
}

/// The steps collected for one expression, by id.
#[derive(Debug)]
struct Frame {
    /// The subject of an attempt, or what an evaluation evaluates.
    subject: Expr,
    attempt: bool,
    /// The path from `subject` to the part being evaluated.
    position: Vec<usize>,
    steps: Vec<usize>,
    /// In an attempt at [`TraceLevel::Verbose`]: the steps of the sources
    /// that did not rewrite, and how many of those were patterned rules.
    failed: Vec<usize>,
    unmatched: u32,
}

impl Frame {
    fn new(subject: &Expr, attempt: bool) -> Frame {
        Frame {
            subject: subject.clone(),
            attempt,
            position: Vec::new(),
            steps: Vec::new(),
            failed: Vec::new(),
            unmatched: 0,
        }
    }
}

impl Tracer {
    /// A trace at `level`, not `Off`, of the evaluation of `input`.
    pub(crate) fn new(level: TraceLevel, input: &Expr) -> Tracer {
        Tracer {
            level,
            steps: Vec::new(),
            dropped: Vec::new(),
            frames: vec![Frame::new(input, false)],
            first: HashMap::with_hasher(FixedState::default()),
            source: None,
        }
    }

    pub(crate) fn level(&self) -> TraceLevel {
        self.level
    }

    fn verbose(&self) -> bool {
        self.level == TraceLevel::Verbose
    }

    fn top(&mut self) -> &mut Frame {
        self.frames.last_mut().expect("the root frame stays")
    }

    /// Enters part `i` of the expression being evaluated.
    pub(crate) fn enter(&mut self, i: usize) {
        self.top().position.push(i);
    }

    /// Leaves the part entered last.
    pub(crate) fn leave(&mut self) {
        self.top().position.pop();
    }

    /// The id the next step will have.
    pub(crate) fn next_id(&self) -> usize {
        self.steps.len()
    }

    /// Records a step at the current position, the parent of `children`.
    pub(crate) fn record(
        &mut self,
        kind: StepKind,
        before: &Expr,
        after: &Expr,
        children: Vec<usize>,
    ) {
        let id = self.steps.len();
        for &c in &children {
            self.steps[c].parent = Some(id);
        }
        let top = self.frames.last_mut().expect("the root frame stays");
        self.steps.push(Step {
            id,
            kind,
            position: top.position.clone(),
            before: before.clone(),
            after: after.clone(),
            parent: None,
            children,
        });
        self.dropped.push(false);
        top.steps.push(id);
    }

    /// Drops `ids` and their descendants.
    fn drop_steps(&mut self, mut ids: Vec<usize>) {
        while let Some(id) = ids.pop() {
            self.dropped[id] = true;
            ids.extend_from_slice(&self.steps[id].children);
        }
    }

    /// A memo hit: `before`, met before, evaluated to `after`.
    pub(crate) fn memo(&mut self, before: &Expr, after: &Expr) {
        let first = self.first.get(before).copied();
        self.record(StepKind::Memo { first }, before, after, Vec::new());
    }

    /// Notes, for each part of `passed` memoised with `result`, the first
    /// step of its evaluation, which started when the next id was `first`.
    pub(crate) fn memoised(&mut self, passed: &[Expr], result: &Expr, first: usize) {
        if self.steps.len() > first {
            for p in passed.iter().filter(|p| *p != result) {
                self.first.insert(p.clone(), first);
            }
        }
    }

    /// Opens an attempt to rewrite `subject`, or an evaluation of it.
    pub(crate) fn open(&mut self, subject: &Expr, attempt: bool) {
        self.frames.push(Frame::new(subject, attempt));
    }

    /// What rewrote the subject of the innermost attempt.
    pub(crate) fn rewrote(&mut self, kind: StepKind) {
        self.source = Some(kind);
    }

    /// A source tried on the subject of the innermost attempt that did not
    /// rewrite it: a patterned rule (`counted`) or a native function. Its
    /// steps are kept for an [`StepKind::Attempts`] step at the Verbose
    /// level, and dropped otherwise.
    pub(crate) fn failed(&mut self, counted: bool) {
        let verbose = self.verbose();
        let top = self.top();
        if counted {
            top.unmatched += 1;
        }
        let steps = std::mem::take(&mut top.steps);
        if verbose {
            top.failed.extend(steps);
        } else {
            self.drop_steps(steps);
        }
    }

    /// Closes the innermost attempt, whose subject became `after`: `None`
    /// if nothing rewrote it, `Err` if evaluation failed in it.
    pub(crate) fn close_attempt(&mut self, after: Result<Option<&Expr>, ()>) {
        let frame = self.frames.pop().expect("an attempt is open");
        debug_assert!(frame.attempt, "an attempt closes");
        let source = self.source.take();
        let subject = &frame.subject;
        let attempts = StepKind::Attempts(frame.unmatched);
        match after {
            Ok(Some(after)) if after != subject => {
                if self.verbose() && (frame.unmatched > 0 || !frame.failed.is_empty()) {
                    self.record(attempts, subject, subject, frame.failed);
                }
                let kind = source.expect("what rewrote is known");
                self.record(kind, subject, after, frame.steps);
            }
            Ok(_) if self.verbose() => {
                let mut failed = frame.failed;
                failed.extend(frame.steps);
                if frame.unmatched > 0 || !failed.is_empty() {
                    self.record(attempts, subject, subject, failed);
                }
            }
            Ok(_) => {
                debug_assert!(frame.failed.is_empty(), "only Verbose keeps failures");
                self.drop_steps(frame.steps);
            }
            Err(()) => {
                let mut children = frame.failed;
                children.extend(frame.steps);
                self.record(attempts, subject, subject, children);
            }
        }
    }

    /// Closes the innermost evaluation, of its subject to `result`, or
    /// failed: then `after` is the subject.
    pub(crate) fn close_evaluation(&mut self, result: Result<&Expr, ()>) {
        let frame = self.frames.pop().expect("an evaluation is open");
        debug_assert!(!frame.attempt, "an evaluation closes");
        let after = result.unwrap_or(&frame.subject);
        if !frame.steps.is_empty() || *after != frame.subject {
            self.record(StepKind::Evaluation, &frame.subject, after, frame.steps);
        }
    }

    /// The trace. Frames still open close as after an error: an attempt as
    /// [`StepKind::Attempts`], an evaluation with `after` its subject. The
    /// steps dropped are left out and the others numbered again in order; a
    /// memo step whose earlier evaluation was dropped points nowhere.
    pub(crate) fn finish(mut self) -> Trace {
        while self.frames.len() > 1 {
            if self.frames.last().is_some_and(|f| f.attempt) {
                self.close_attempt(Err(()));
            } else {
                self.close_evaluation(Err(()));
            }
        }
        let mut new = vec![None; self.steps.len()];
        let mut kept = 0;
        for (n, &dropped) in new.iter_mut().zip(&self.dropped) {
            if !dropped {
                *n = Some(kept);
                kept += 1;
            }
        }
        let steps = self
            .steps
            .into_iter()
            .zip(self.dropped)
            .filter(|(_, dropped)| !dropped)
            .map(|(mut s, _)| {
                s.id = new[s.id].expect("a step kept");
                s.parent = s.parent.and_then(|p| new[p]);
                for c in &mut s.children {
                    *c = new[*c].expect("the children of a step kept are kept");
                }
                if let StepKind::Memo { first } = &mut s.kind {
                    *first = first.and_then(|f| new[f]);
                }
                s
            })
            .collect();
        Trace { steps }
    }
}

#[cfg(test)]
mod tests;
