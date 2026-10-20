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
//! 6. Rewrite, by the first source that applies: up-rules, the head's native
//!    function (the builtins' are in `native`, CAS-54), system rules,
//!    down-rules, sub-rules. The rules come with CAS-55 and CAS-56.
//! 7. A change in step 5 or 6 goes back to step 1, counted against
//!    [`Limits::iterations`]. Otherwise the expression is a fixpoint: it is
//!    stamped, so evaluating it again in this state is O(1), and memoised
//!    for the rest of the call.
//!
//! The evaluator is recursive. Each nested evaluation counts against
//! [`Limits::recursion`], which stops it before the stack runs out as long
//! as each level has [`STACK_PER_LEVEL`] bytes of stack.

use std::borrow::Cow;
use std::fmt;

use foldhash::fast::FixedState;
use hashbrown::HashMap;

use crate::attributes::Attributes;
use crate::build::flatten;
use crate::context::{Context, EvalError, Limits};
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

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
/// parts and read the context.
///
/// It holds the call's memo, which maps each node evaluated so far to its
/// result, so a repeated subexpression is evaluated once per call, and the
/// recursion depth. The evaluation report's warnings and tracer (CAS-61) and
/// the local frames of procedures (the Language module) will live here too.
pub struct Eval<'c> {
    context: &'c Context,
    recursion: u32,
    iterations: u32,
    depth: u32,
    memo: HashMap<Expr, Expr, FixedState>,
    /// Whether stamps and the memo are used: off only in the tests that
    /// check they change no result.
    caches: bool,
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
enum Step {
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
        }
    }

    /// The context evaluation runs in.
    #[must_use]
    pub fn context(&self) -> &'c Context {
        self.context
    }

    /// Evaluates `e` as part of this call: one level deeper, with its memo.
    ///
    /// # Errors
    ///
    /// [`EvalError`] when a limit is reached or an operation has no result.
    pub fn eval(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        if self.depth >= self.recursion {
            return Err(EvalError::RecursionLimit);
        }
        self.depth += 1;
        let result = self.evaluate(e);
        self.depth -= 1;
        result
    }

    /// Steps 1 to 7.
    fn evaluate(&mut self, input: &Expr) -> Result<Expr, EvalError> {
        let mut e = input.clone();
        // The nodes passed on the way, memoised with the result.
        let mut passed = Vec::new();
        let mut rewrites = 0;
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
                    break result.clone();
                }
            }
            passed.push(e.clone());
            match self.step(&e)? {
                Step::Again(next) => {
                    e = next;
                    self.count(&mut rewrites)?;
                }
                Step::Done(result) => {
                    if self.caches {
                        self.context.mark_fixpoint(&result);
                    }
                    break result;
                }
            }
        };
        if self.caches {
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

    /// Steps 3 to 6 on an application.
    fn step(&mut self, e: &Expr) -> Result<Step, EvalError> {
        let View::Apply { head, args } = e.view() else {
            unreachable!("only applications take steps");
        };
        let new_head = self.eval(head)?;
        let attributes = new_head
            .as_symbol()
            .map_or(Attributes::NONE, |s| self.context.attributes(s));
        let mut evaluated = Vec::with_capacity(args.len());
        for (i, a) in args.iter().enumerate() {
            let hold = if i == 0 {
                Attributes::HOLD_FIRST
            } else {
                Attributes::HOLD_REST
            };
            evaluated.push(if attributes.contains(hold) {
                a.clone()
            } else {
                self.eval(a)?
            });
        }
        match normalise(&new_head, attributes, &evaluated)? {
            Normal::Same => {}
            Normal::Changed(args) => return Ok(Step::Again(Expr::apply(&new_head, &args))),
            Normal::Threaded(list) => return Ok(Step::Again(list)),
        }
        let e = if new_head == *head && evaluated == args {
            e.clone()
        } else {
            Expr::apply(&new_head, &evaluated)
        };
        Ok(match self.rewrite(&e, new_head.as_symbol())? {
            Some(rewritten) if rewritten != e => Step::Again(rewritten),
            _ => Step::Done(e),
        })
    }

    /// Step 6: the first source that rewrites `e`, an application with
    /// `head`.
    fn rewrite(&mut self, e: &Expr, head: Option<Symbol>) -> Result<Option<Expr>, EvalError> {
        // 1. The up-rules of the arguments' heads: CAS-56.
        // 2. The head's native function.
        if let Some(native) = head.and_then(native) {
            if let Some(rewritten) = native(self, e)? {
                return Ok(Some(rewritten));
            }
        }
        // 3. System rules, 4. down-rules, 5. sub-rules and the application of
        // `Function` heads: CAS-55 and CAS-56.
        Ok(None)
    }
}

/// The native function of a builtin (CAS-54), or one a test gave.
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
enum Normal {
    Same,
    Changed(Vec<Expr>),
    /// A `Listable` head threaded over lists: the list of its applications.
    Threaded(Expr),
}

/// Step 5 on the evaluated arguments of `head`, whose attributes are
/// `attributes`.
fn normalise(head: &Expr, attributes: Attributes, args: &[Expr]) -> Result<Normal, EvalError> {
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
/// in it. In `Times`, lists come after the other factors, in their own
/// order.
pub(crate) fn sort_arguments(head: &Expr, args: &[Expr]) -> Option<Vec<Expr>> {
    let sorted: Vec<Expr> = if head.as_symbol() == Some(Symbol::TIMES) && args.iter().any(is_list) {
        let (lists, mut factors): (Vec<Expr>, Vec<Expr>) = args.iter().cloned().partition(is_list);
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
        let limits = self.limits();
        if limits.stack == 0 {
            return Eval::new(self, limits.recursion).eval(e);
        }
        let recursion = limits.recursion.min(levels(limits.stack));
        std::thread::scope(|scope| {
            std::thread::Builder::new()
                .name("cas-eval".into())
                .stack_size(limits.stack)
                .spawn_scoped(scope, || Eval::new(self, recursion).eval(e))
                .expect("the evaluation thread starts")
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
        })
    }

    /// [`Context::eval`] without fixpoint stamps or memo: the same results,
    /// which tests check.
    #[cfg(test)]
    pub(crate) fn eval_uncached(&self, e: &Expr) -> Result<Expr, EvalError> {
        let mut ev = Eval::new(self, self.limits().recursion);
        ev.caches = false;
        ev.eval(e)
    }
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
