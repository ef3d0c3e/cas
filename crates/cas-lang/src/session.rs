//! The session (CAS-DOC-7 §4): a program run one statement at a time
//! against a [`Context`].
//!
//! [`Session::run`] parses and lowers a text, then runs each statement in
//! order, in the state the ones before it left:
//!
//! - an expression is evaluated, and its result printed in input form;
//! - `x := e` evaluates `e` once and assigns it (`define`); `[a, b] := e`
//!   assigns the elements of a list in order, and `v[k] := e` a copy of `v`
//!   with an element replaced;
//! - `f(x) := …` is a delayed rule (`define_delayed`), and a definition of a
//!   builtin an up-rule on each user name at level 1 (`define_up_delayed`);
//! - `assume(…)` replaces the context with a scope that assumes the facts,
//!   and keeps it; `purge(…)` clears the names and forgets the assumptions
//!   on them, as in Xcas.
//!
//! Each statement gives an [`Output`]: the result unless `;` follows it or
//! it is a definition, and [`Message`]s: the parser's and lowering's
//! diagnostics, evaluation's warnings and errors. An evaluation error points
//! at the part of the statement that fails: the smallest part whose own
//! evaluation fails the same way.
//!
//! An evaluation can be stopped from another thread with the session's
//! [`Interrupt`]: the statement in progress ends with an error, and the
//! session goes on with the next one.

use std::ops::Range;

use cas_core::{
    Attributes, Context, DefineError, EvalError, EvalOptions, Expr, Interrupt, Limits, Symbol, View,
};

use crate::lower::{Fix, Severity, Statement, lower};
use crate::names::{Names, unescape_name};
use crate::parser::parse;
use crate::print::print;

/// Messages write out an expression of at most this many parts, as
/// cas-core's do: a larger one could take long, sharing parts that its text
/// repeats.
const MAX_WRITTEN: u32 = 1000;

/// A context, the names of the language, and the handle that interrupts
/// evaluation.
#[derive(Debug)]
pub struct Session {
    context: Context,
    names: Names,
    interrupt: Interrupt,
}

/// What one statement of a program gave.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Output {
    /// The text of the statement.
    pub range: Range<usize>,
    /// Its result, if it shows one: not when `;` follows, nor for
    /// definitions, `assume`, `purge` and `Null`.
    pub result: Option<Shown>,
    /// What there is to say about it, in order: the diagnostics of the
    /// parser and of lowering by position, then those of evaluation.
    pub messages: Vec<Message>,
}

/// A result as the session shows it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Shown {
    /// In input form, which reads back as the result.
    Input(String),
    /// In the evaluator's notation, where the language cannot write the
    /// result; a note says why.
    Evaluator(String),
}

impl Shown {
    /// The text, in either notation.
    #[must_use]
    pub fn text(&self) -> &str {
        match self {
            Shown::Input(text) | Shown::Evaluator(text) => text,
        }
    }
}

/// An error, a warning or a note about part of a program.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Message {
    /// How serious it is.
    pub level: Level,
    /// The text it is about.
    pub range: Range<usize>,
    /// What it says.
    pub text: String,
    /// A replacement that fixes it, from lowering.
    pub fix: Option<Fix>,
}

/// How serious a [`Message`] is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    /// The statement did not run, or stopped.
    Error,
    /// Evaluation went past something wrong: `1/0`.
    Warning,
    /// A remark: the product reading of `x(x + 1)`.
    Note,
}

impl Default for Session {
    fn default() -> Session {
        Session::new()
    }
}

impl Session {
    /// A session with the English names and the limits for deep evaluation
    /// ([`Limits::deep`]).
    #[must_use]
    pub fn new() -> Session {
        Session::with_names(Names::english().clone())
    }

    /// A session with the names of a locale ([`Names::with_locale`]).
    #[must_use]
    pub fn with_names(names: Names) -> Session {
        let mut context = Context::new();
        context.set_limits(Limits::deep());
        let interrupt = Interrupt::new();
        context.set_interrupt(interrupt.clone());
        Session {
            context,
            names,
            interrupt,
        }
    }

    /// The context, with what the statements so far defined and assumed.
    #[must_use]
    pub fn context(&self) -> &Context {
        &self.context
    }

    /// Sets the limits of evaluation, [`Limits::deep`] in a new session.
    pub fn set_limits(&mut self, limits: Limits) {
        self.context.set_limits(limits);
    }

    /// The names of the language in this session.
    #[must_use]
    pub fn names(&self) -> &Names {
        &self.names
    }

    /// A handle that stops the statement in progress from another thread:
    /// its evaluation ends with an error, and the next statement runs.
    #[must_use]
    pub fn interrupt(&self) -> Interrupt {
        self.interrupt.clone()
    }

    /// Runs the statements of `text` in order, and says what each gave.
    ///
    /// ```
    /// use cas_lang::Session;
    ///
    /// let mut session = Session::new();
    /// let outputs = session.run("f(x) := x^2 + 1\ny := f(3);\nf(y) - 1");
    /// let results: Vec<_> = outputs
    ///     .iter()
    ///     .map(|o| o.result.as_ref().map(|r| r.text()))
    ///     .collect();
    /// assert_eq!(results, [None, None, Some("100")]);
    /// ```
    pub fn run(&mut self, text: &str) -> Vec<Output> {
        let parse = parse(text);
        let lowered = lower(&parse, &self.names);
        let mut diagnostics: Vec<Message> = parse
            .errors()
            .iter()
            .map(|e| Message {
                level: Level::Error,
                range: e.range.clone(),
                text: e.to_string(),
                fix: None,
            })
            .chain(lowered.diagnostics.iter().map(|d| Message {
                level: match d.severity {
                    Severity::Error => Level::Error,
                    Severity::Note => Level::Note,
                },
                range: d.range.clone(),
                text: d.to_string(),
                fix: d.fix.clone(),
            }))
            .collect();
        diagnostics.sort_by_key(|m| m.range.start);
        let mut diagnostics = diagnostics.into_iter().peekable();
        let mut statements = lowered.statements.iter().peekable();
        let mut outputs = Vec::new();
        for node in parse.syntax().children() {
            let range =
                usize::from(node.text_range().start())..usize::from(node.text_range().end());
            let mut output = Output {
                range: range.clone(),
                result: None,
                messages: Vec::new(),
            };
            // A diagnostic before this statement, past the last one, is its own.
            while let Some(d) = diagnostics.next_if(|d| d.range.start < range.start) {
                outputs.push(Output {
                    range: d.range.clone(),
                    result: None,
                    messages: vec![d],
                });
            }
            while let Some(d) = diagnostics.next_if(|d| d.range.start < range.end) {
                output.messages.push(d);
            }
            if let Some(statement) = statements.next_if(|s| {
                range.start <= s.spans.range().start && s.spans.range().end <= range.end
            }) {
                self.execute(statement, &mut output);
            }
            outputs.push(output);
        }
        for d in diagnostics {
            outputs.push(Output {
                range: d.range.clone(),
                result: None,
                messages: vec![d],
            });
        }
        outputs
    }

    /// Runs one statement.
    fn execute(&mut self, statement: &Statement, out: &mut Output) {
        self.interrupt.reset();
        let e = &statement.expr;
        match (e.head().and_then(Expr::as_symbol), e.args()) {
            (Some(Symbol::SET), [target, _]) => self.assign(statement, target, out),
            (Some(Symbol::SET_DELAYED), [lhs, rhs]) => {
                let result = self.context.define_delayed(lhs, rhs);
                defined(statement.spans.range(), result, out);
            }
            (Some(Symbol::UP_SET_DELAYED), [lhs, rhs]) => {
                for tag in tags(lhs) {
                    let result = self.context.define_up_delayed(tag, lhs, rhs);
                    if !defined(statement.spans.range(), result, out) {
                        break;
                    }
                }
            }
            (Some(Symbol::ASSUME), facts) => {
                for fact in facts {
                    self.context = self.context.assuming(fact.clone());
                }
            }
            (Some(Symbol::CLEAR), names) => {
                for s in names.iter().filter_map(Expr::as_symbol) {
                    match self.context.clear(s) {
                        Ok(()) => self.context = self.context.forgetting(s),
                        Err(e) => out
                            .messages
                            .push(error(statement.spans.range(), e.to_string())),
                    }
                }
            }
            _ => {
                if let Some(value) = self.evaluate(statement, &[], out) {
                    if !statement.quiet {
                        self.show(&value, statement.spans.range(), out);
                    }
                }
            }
        }
    }

    /// `Set[target, value]`: the value evaluated once, then assigned to a
    /// name, or element by element to a list of names.
    fn assign(&mut self, statement: &Statement, target: &Expr, out: &mut Output) {
        let Some(value) = self.evaluate(statement, &[2], out) else {
            return;
        };
        let spans = &statement.spans;
        if let View::Symbol(s) = target.view() {
            if has_head(&statement.expr.args()[1], Symbol::REPLACE_AT)
                && has_head(&value, Symbol::REPLACE_AT)
            {
                let message = format!(
                    "{} is not a list, or the index is neither an integer nor a slice",
                    spelling(s)
                );
                out.messages.push(error(spans.find(&[1]), message));
                return;
            }
            if !self.assign_one(s, &value, spans.find(&[1]), out) {
                return;
            }
        } else {
            let names = target.args();
            let fits = has_head(&value, Symbol::LIST) && value.args().len() == names.len();
            if !fits {
                let message = format!("this is not a list of {} elements", names.len());
                out.messages.push(error(spans.find(&[2]), message));
                return;
            }
            for (i, (name, element)) in names.iter().zip(value.args()).enumerate() {
                let Some(s) = name.as_symbol() else {
                    continue;
                };
                if !self.assign_one(s, element, spans.find(&[1, i + 1]), out) {
                    return;
                }
            }
        }
        if !statement.quiet {
            self.show(&value, spans.range(), out);
        }
    }

    /// `s := value`, unless `value` refers to `s` while `s` has no value:
    /// evaluating it would not end (giac's "recursive definition").
    fn assign_one(
        &mut self,
        s: Symbol,
        value: &Expr,
        range: Range<usize>,
        out: &mut Output,
    ) -> bool {
        let symbol = Expr::from(s);
        if self.context.value(s).is_none() && !value.free_of(&symbol) {
            out.messages.push(error(
                range,
                format!(
                    "recursive definition: the value refers to {}, which has no value",
                    spelling(s)
                ),
            ));
            return false;
        }
        let result = self.context.define(&symbol, value);
        defined(range, result, out)
    }

    /// Evaluates the part of `statement` at `path`, with its warnings; on an
    /// error, says where it comes from.
    fn evaluate(&self, statement: &Statement, path: &[usize], out: &mut Output) -> Option<Expr> {
        let e = statement.expr.part(path)?;
        let (result, report) = self.context.eval_with(e, &EvalOptions::default());
        let range = statement.spans.find(path);
        for w in &report.warnings {
            out.messages.push(Message {
                level: Level::Warning,
                range: range.clone(),
                text: w.message(|e| self.written(e)),
                fix: None,
            });
        }
        match result {
            Ok(value) => Some(value),
            Err(e) => {
                let mut at = path.to_vec();
                if e != EvalError::Interrupted {
                    at.extend(self.failing(statement.expr.part(path)?, e));
                }
                out.messages
                    .push(error(statement.spans.find(&at), e.to_string()));
                None
            }
        }
    }

    /// The path in `e` to the smallest part whose own evaluation fails with
    /// `error`: down through the parts that are evaluated, while one of them
    /// fails alone. Those that evaluate are stamped, so trying them again is
    /// cheap.
    fn failing(&self, e: &Expr, error: EvalError) -> Vec<usize> {
        let mut path = Vec::new();
        let mut current = e.clone();
        'down: while let View::Apply { head, args } = current.view() {
            if self.context.eval(head) == Err(error) {
                path.push(0);
                current = head.clone();
                continue;
            }
            let attributes = head
                .as_symbol()
                .map_or(Attributes::NONE, |s| self.context.attributes(s));
            for (i, a) in args.iter().enumerate() {
                let held = if i == 0 {
                    Attributes::HOLD_FIRST
                } else {
                    Attributes::HOLD_REST
                };
                if !attributes.contains(held) && self.context.eval(a) == Err(error) {
                    path.push(i + 1);
                    current = a.clone();
                    continue 'down;
                }
            }
            break;
        }
        path
    }

    /// `e` for a message: in input form, or in the evaluator's notation
    /// where the language cannot write it; a large one by its head.
    fn written(&self, e: &Expr) -> String {
        if e.size() > MAX_WRITTEN {
            return match e.head() {
                Some(head) if head.size() <= MAX_WRITTEN => format!("{}(…)", self.written(head)),
                _ => "…".to_owned(),
            };
        }
        print(e, &self.names).unwrap_or_else(|_| e.to_string())
    }

    /// The result of a statement in input form; nothing for `Null`, and the
    /// evaluator's notation, with a note, for what the language cannot
    /// write.
    fn show(&self, value: &Expr, range: Range<usize>, out: &mut Output) {
        if value.as_symbol() == Some(Symbol::NULL) {
            return;
        }
        match print(value, &self.names) {
            Ok(text) => out.result = Some(Shown::Input(text)),
            Err(e) => {
                out.result = Some(Shown::Evaluator(value.to_string()));
                out.messages.push(Message {
                    level: Level::Note,
                    range,
                    text: format!("shown in the evaluator's notation: {e}"),
                    fix: None,
                });
            }
        }
    }
}

/// The user names an up-rule attaches to: those at level 1 of `lhs`, as
/// arguments or heads of arguments (CAS-DOC-4 §5.2).
fn tags(lhs: &Expr) -> Vec<Symbol> {
    let mut tags: Vec<Symbol> = Vec::new();
    for arg in lhs.args() {
        let symbol = match arg.view() {
            View::Symbol(s) => Some(s),
            View::Apply { head, .. } => head.as_symbol(),
            _ => None,
        };
        if let Some(s) = symbol.filter(|s| !s.is_builtin()) {
            if !tags.contains(&s) {
                tags.push(s);
            }
        }
    }
    tags
}

/// Reports a refused definition, about the text `range`; whether it was
/// made.
fn defined(range: Range<usize>, result: Result<(), DefineError>, out: &mut Output) -> bool {
    match result {
        Ok(()) => true,
        Err(e) => {
            out.messages.push(error(range, e.to_string()));
            false
        }
    }
}

fn has_head(e: &Expr, s: Symbol) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(s)
}

/// How the user writes the name of `s`.
fn spelling(s: Symbol) -> &'static str {
    unescape_name(s.name())
}

fn error(range: Range<usize>, text: String) -> Message {
    Message {
        level: Level::Error,
        range,
        text,
        fix: None,
    }
}

#[cfg(test)]
mod tests;
