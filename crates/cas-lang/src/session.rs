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
//!   on them, as in Xcas;
//! - `Digits := n` sets the working precision (CAS-151), and
//!   `all_trig_solutions := v`, `with_sqrt := v` or giac's calls
//!   `all_trig_solutions(v)`, `with_sqrt(v)` the context's flags, `v` 0 or
//!   1 (CAS-545);
//! - Xcas's settings of phase 1.6 (CAS-617, CAS-618) are set the same ways,
//!   `name := v` or `name(v)`, and `name()` reads one: `approx_mode`,
//!   `complex_mode`, `complex_variables` (the domain), `keep_algext`,
//!   `increasing_power` (0 or 1), `angle_radian` (0 degrees, 1 radians, 2
//!   grads), `epsilon` (a positive number), `eval_level` (a positive
//!   integer), `autosimplify` (0 none, 1 `regroup`, 2 `simplify`, or a
//!   command's name; it answers 1), and giac's `rur_certify`,
//!   `gbasis_max_pairs`, `gbasis_simult_primes`, answered as giac does
//!   (solve.cc). `keep_algext`, `eval_level`, `rur_certify` and the
//!   `gbasis_*` are stored only: this CAS evaluates to a fixpoint, and its
//!   Gröbner bases and rational univariate representations are exact.
//!   `cas_setup(l)`, a list of at least 7 elements or as many arguments,
//!   sets giac's elements 0 to 10 (`approx_mode`, `complex_variables`,
//!   `complex_mode`, `angle_radian`, -, `epsilon`, `Digits`, `[…,
//!   eval_level]`, `increasing_power`, `with_sqrt`, `all_trig_solutions`)
//!   and answers `l`; `cas_setup()` lists them;
//! - an expression's result goes through the `autosimplify` command when
//!   one is set, then `N` when `approx_mode` is 1.
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
use std::sync::Arc;

use cas_core::{
    Angle, Attributes, Context, DefineError, EvalError, EvalOptions, Expr, Interrupt, Limits,
    Symbol, View,
};
use cas_num::ball::Ball;

use crate::hook::Hook;
use crate::literal::decimal_balls;
use crate::lower::{Fix, Severity, Statement, lower};
use crate::names::{Names, Target, unescape_name};
use crate::parser::parse;
use crate::print::{Balls, PrintError, power_order, print_with};

/// Messages write out an expression of at most this many parts, as
/// cas-core's do: a larger one could take long, sharing parts that its text
/// repeats.
const MAX_WRITTEN: u32 = 1000;

/// The settings `cas_setup(l)` sets, at giac's positions in `l` (prog.cc):
/// `Null` for the one ignored.
const CAS_SETUP: [Symbol; 11] = [
    Symbol::APPROX_MODE,
    Symbol::COMPLEX_VARIABLES,
    Symbol::COMPLEX_MODE,
    Symbol::ANGLE_RADIAN,
    Symbol::NULL,
    Symbol::EPSILON,
    Symbol::DIGITS,
    Symbol::EVAL_LEVEL,
    Symbol::INCREASING_POWER,
    Symbol::WITH_SQRT,
    Symbol::ALL_TRIG_SOLUTIONS,
];

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

    /// A session with the names of a locale ([`Names::with_locale`]),
    /// which its context's language hook ([`Hook`]) uses too.
    #[must_use]
    pub fn with_names(names: Names) -> Session {
        let mut context = Context::new();
        context.set_limits(Limits::deep());
        context.set_language(Some(Arc::new(Hook::new(names.clone()))));
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

    /// The bits balls are made with, 53 in a new session (CAS-DOC-5 §6.4):
    /// decimal literals are read at it, and balls shown and compared at it.
    #[must_use]
    pub fn working_precision(&self) -> u64 {
        self.context.working_precision()
    }

    /// Sets the working precision, for the statements run after it.
    pub fn set_working_precision(&mut self, bits: u64) {
        self.context.set_working_precision(bits);
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

    /// How balls are shown: as the decimals they show at the working
    /// precision (D2).
    fn shown(&self) -> Balls {
        Balls::Shown {
            working_precision: self.context.working_precision(),
        }
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
                // Its decimals are made balls now, at the precision in force
                // (CAS-148).
                let statement = Statement {
                    expr: decimal_balls(&statement.expr, self.context.working_precision()),
                    quiet: statement.quiet,
                    spans: statement.spans.clone(),
                };
                self.execute(&statement, &mut output);
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
            // `all_trig_solutions(v)`, `with_sqrt(v)`: giac's form of
            // setting them (CAS-545), and of phase 1.6's (CAS-617, CAS-618).
            (Some(s), [_]) if is_setting(s) => {
                if let Some(value) = self.evaluate(statement, &[1], out) {
                    self.set(statement, s, &value, &[1], out);
                }
            }
            // `name()` reads a setting: the head alone would evaluate to
            // its value before a native could see the call.
            (Some(s), []) if is_setting(s) => {
                if !statement.quiet {
                    if let Ok(value) = self.context.eval(&Expr::from(s)) {
                        self.show(&value, statement.spans.range(), out);
                    }
                }
            }
            (Some(Symbol::CAS_SETUP), [_, ..]) => self.cas_setup(statement, out),
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
            _ => self.expression(statement, out),
        }
    }

    /// An expression statement: its value, through the `autosimplify`
    /// command when one is set, then `N` when `approx_mode` is 1
    /// (CAS-617).
    fn expression(&mut self, statement: &Statement, out: &mut Output) {
        let Some(mut value) = self.evaluate(statement, &[], out) else {
            return;
        };
        if statement.quiet {
            return;
        }
        let mut commands = Vec::new();
        commands.extend(self.context.autosimplify().cloned());
        if self.context.approx_mode() {
            commands.push(Expr::from(Symbol::N));
        }
        for command in commands {
            match self.context.eval(&Expr::apply(&command, &[value])) {
                Ok(next) => value = next,
                Err(e) => {
                    out.messages
                        .push(error(statement.spans.range(), e.to_string()));
                    return;
                }
            }
        }
        self.show(&value, statement.spans.range(), out);
    }

    /// `cas_setup(l)`, giac's (prog.cc): `l` a list of at least 7 elements,
    /// or as many arguments, sets the settings at giac's positions 0 to 10
    /// and is the answer; others are left to evaluation, where they stay.
    fn cas_setup(&mut self, statement: &Statement, out: &mut Output) {
        let args = statement.expr.args();
        let elements = if let [_] = args {
            let Some(value) = self.evaluate(statement, &[1], out) else {
                return;
            };
            if !has_head(&value, Symbol::LIST) || value.args().len() < 7 {
                return self.expression(statement, out);
            }
            value
        } else if args.len() >= 7 {
            let mut values = Vec::with_capacity(args.len());
            for i in 1..=args.len() {
                let Some(value) = self.evaluate(statement, &[i], out) else {
                    return;
                };
                values.push(value);
            }
            Expr::call(Symbol::LIST, values)
        } else {
            return self.expression(statement, out);
        };
        let mut context = self.context.clone();
        for (&s, element) in CAS_SETUP.iter().zip(elements.args()) {
            let value = match s {
                Symbol::NULL => continue,
                Symbol::EPSILON if has_head(element, Symbol::LIST) => match element.args() {
                    [first, ..] => first,
                    [] => continue,
                },
                Symbol::EVAL_LEVEL => match element.args() {
                    [_, _, _, level, ..] if has_head(element, Symbol::LIST) => level,
                    _ => continue,
                },
                _ => element,
            };
            match changed(&context, &self.names, s, value) {
                Ok((next, _)) => context = next,
                Err(message) => {
                    out.messages.push(error(statement.spans.range(), message));
                    return;
                }
            }
        }
        self.context = context;
        if !statement.quiet {
            self.show(&elements, statement.spans.range(), out);
        }
    }

    /// `name := v` or `name(v)` for a setting `s` (`Digits`, CAS-151; the
    /// flags of CAS-545; phase 1.6's, CAS-617, CAS-618): the context
    /// changed for the statements after it, and the answer shown; a value
    /// the setting does not take is an error at `path`, and nothing is set.
    fn set(
        &mut self,
        statement: &Statement,
        s: Symbol,
        value: &Expr,
        path: &[usize],
        out: &mut Output,
    ) {
        match changed(&self.context, &self.names, s, value) {
            Ok((context, answer)) => {
                self.context = context;
                if !statement.quiet {
                    self.show(&answer, statement.spans.range(), out);
                }
            }
            Err(message) => out
                .messages
                .push(error(statement.spans.find(path), message)),
        }
    }

    /// `Set[target, value]`: the value evaluated once, then assigned to a
    /// name, or element by element to a list of names.
    fn assign(&mut self, statement: &Statement, target: &Expr, out: &mut Output) {
        let Some(value) = self.evaluate(statement, &[2], out) else {
            return;
        };
        if let Some(s) = target
            .as_symbol()
            .filter(|&s| s == Symbol::DIGITS || is_setting(s))
        {
            self.set(statement, s, &value, &[2], out);
            return;
        }
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
                level: if w.is_note() {
                    Level::Note
                } else {
                    Level::Warning
                },
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
        self.text(e).unwrap_or_else(|_| e.to_string())
    }

    /// `e` in input form, its sums of powers of one symbol in the order of
    /// the context's `increasing_power` ([`power_order`], CAS-618).
    fn text(&self, e: &Expr) -> Result<String, PrintError> {
        let ordered = power_order(e, self.context.increasing_power());
        print_with(&ordered, &self.names, self.shown())
    }

    /// The result of a statement in input form; nothing for `Null`, and the
    /// evaluator's notation, with a note, for what the language cannot
    /// write.
    fn show(&self, value: &Expr, range: Range<usize>, out: &mut Output) {
        if value.as_symbol() == Some(Symbol::NULL) {
            return;
        }
        match self.text(value) {
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

/// Whether `s` is a setting giac sets by its call `s(v)` and reads by `s()`:
/// `all_trig_solutions`, `with_sqrt` (CAS-545) and phase 1.6's (CAS-617,
/// CAS-618).
fn is_setting(s: Symbol) -> bool {
    matches!(
        s,
        Symbol::ALL_TRIG_SOLUTIONS
            | Symbol::WITH_SQRT
            | Symbol::APPROX_MODE
            | Symbol::COMPLEX_MODE
            | Symbol::COMPLEX_VARIABLES
            | Symbol::ANGLE_RADIAN
            | Symbol::KEEP_ALGEXT
            | Symbol::INCREASING_POWER
            | Symbol::EPSILON
            | Symbol::EVAL_LEVEL
            | Symbol::AUTO_SIMPLIFY
            | Symbol::RUR_CERTIFY
            | Symbol::GBASIS_MAX_PAIRS
            | Symbol::GBASIS_SIMULT_PRIMES
    )
}

/// `context` with the setting `s` set from `value`, and the answer of
/// setting it: `value` itself, but `autosimplify`'s 1, and giac's answers
/// of `rur_certify` (the stored value, or `-n` for `n < 0`, unset) and of
/// `gbasis_max_pairs`, `gbasis_simult_primes` (the old value; solve.cc,
/// giac 2.0.0.21). An error message for a value `s` does not take.
fn changed(
    context: &Context,
    names: &Names,
    s: Symbol,
    value: &Expr,
) -> Result<(Context, Expr), String> {
    let name = names
        .spelling(&Target::Builtin(s))
        .unwrap_or_else(|| spelling(s));
    let mut context = context.clone();
    let flag = || match (value.as_i64(), value.as_symbol()) {
        (Some(0), _) | (_, Some(Symbol::FALSE)) => Ok(false),
        (Some(1), _) | (_, Some(Symbol::TRUE)) => Ok(true),
        _ => Err(format!("{name} takes 0 or 1, false or true")),
    };
    let integer = || {
        value
            .as_i64()
            .ok_or_else(|| format!("{name} takes an integer"))
    };
    let mut answer = value.clone();
    match s {
        Symbol::DIGITS => {
            let number_bits = context.limits().number_bits;
            let bits = value
                .as_i64()
                .and_then(|n| u64::try_from(n).ok())
                .filter(|&n| n >= 1)
                .map(cas_core::digits_bits)
                .filter(|&bits| bits <= number_bits)
                .ok_or_else(|| {
                    format!(
                        "Digits takes a positive integer, up to {}: the decimal digits balls are made and shown with",
                        cas_core::shown_digits(number_bits)
                    )
                })?;
            context.set_working_precision(bits);
        }
        Symbol::ALL_TRIG_SOLUTIONS => context.set_all_trig_solutions(flag()?),
        Symbol::WITH_SQRT => context.set_with_sqrt(flag()?),
        Symbol::APPROX_MODE => context.set_approx_mode(flag()?),
        Symbol::COMPLEX_MODE => context.set_complex_mode(flag()?),
        Symbol::COMPLEX_VARIABLES => context.set_complex_variables(flag()?),
        Symbol::KEEP_ALGEXT => context.set_keep_algext(flag()?),
        Symbol::INCREASING_POWER => context.set_increasing_power(flag()?),
        Symbol::ANGLE_RADIAN => context.set_angle(match value.as_i64() {
            Some(0) => Angle::Degree,
            Some(1) => Angle::Radian,
            Some(2) => Angle::Grad,
            _ => return Err(format!("{name} takes 0, 1 or 2: degrees, radians or grads")),
        }),
        Symbol::EPSILON => {
            let x = match value.view() {
                View::Ball(b) => Some(b.mid().to_f64()),
                _ => value
                    .to_rational()
                    .map(|q| Ball::from_rational(&q, 64).mid().to_f64()),
            };
            match x {
                Some(x) if x.is_finite() && x > 0.0 => context.set_epsilon(x),
                _ => return Err(format!("{name} takes a positive number")),
            }
        }
        Symbol::EVAL_LEVEL => match value.as_i64().and_then(|n| u64::try_from(n).ok()) {
            Some(n) if n >= 1 => context.set_eval_level(n),
            _ => return Err(format!("{name} takes a positive integer")),
        },
        Symbol::AUTO_SIMPLIFY => {
            context.set_autosimplify(match (value.as_i64(), value.as_symbol()) {
                (Some(0), _) => None,
                (Some(1), _) => Some(Expr::from(Symbol::REGROUP)),
                (Some(2), _) => Some(Expr::from(Symbol::SIMPLIFY)),
                (_, Some(command)) => Some(Expr::from(command)),
                _ => return Err(format!("{name} takes 0, 1, 2 or the name of a command")),
            });
            answer = Expr::from(1);
        }
        // solve.cc:8684: n ≥ 0 stores n - 1, n < 0 answers -n.
        Symbol::RUR_CERTIFY => {
            let n = integer()?;
            if n >= 0 {
                context.set_rur_certify(n - 1);
                answer = Expr::from(n - 1);
            } else {
                answer = Expr::from(n.saturating_neg());
            }
        }
        // solve.cc:9274.
        Symbol::GBASIS_MAX_PAIRS => {
            let n = integer()?;
            answer = Expr::from(context.gbasis_max_pairs());
            context.set_gbasis_max_pairs(if n < 256 { 1 << 30 } else { n });
        }
        // solve.cc:9289.
        Symbol::GBASIS_SIMULT_PRIMES => {
            let n = integer()?;
            answer = Expr::from(context.gbasis_simult_primes());
            context.set_gbasis_simult_primes(n.max(1));
        }
        _ => return Err(format!("{name} is not a setting")),
    }
    Ok((context, answer))
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
