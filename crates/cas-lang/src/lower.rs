//! Lowering: syntax trees to expressions (CAS-DOC-7 §3.1–§7, §9, §11).
//!
//! Each statement becomes one [`Expr`]: the expression to evaluate, or a
//! statement head that the session runs and evaluation leaves alone
//! (`Set`, `SetDelayed`, `UpSetDelayed`, `Assume`, `Clear`; CAS-DOC-4 R1).
//! Expressions are built structurally, in the order of the text: a chain of
//! one operator is flattened (`a + b - c` is one `Plus`), but nothing is
//! sorted or evaluated, so every part keeps its position, which [`Spans`]
//! maps back to the text.
//!
//! Lowering takes the decisions that are not syntax: whether `name(arg)` is
//! a call or a product (giac's rule, §3.1), what a name stands for (§5), the
//! patterns of definitions and rules (§6), and the errors of §11. A
//! statement with an error, the parser's or its own, is left out, and the
//! diagnostics say why.

use std::fmt;
use std::ops::Range;

use cas_core::{Attributes, Expr, Symbol, View, build};
use rowan::NodeOrToken;

use crate::literal::{NumberError, name, number_value, string_value, superscript_value};
use crate::names::{Names, SyntaxForm, Target, TypeName, escape_name};
use crate::parser::Parse;
use crate::syntax::{SyntaxKind as K, SyntaxNode, SyntaxToken};

/// The statements of a text as expressions, and what lowering found.
#[derive(Clone, Debug)]
pub struct Lowered {
    /// The statements without errors, in order.
    pub statements: Vec<Statement>,
    /// Errors and notes, ordered by where they start.
    pub diagnostics: Vec<Diagnostic>,
}

impl Lowered {
    /// Whether lowering found an error; notes do not count.
    #[must_use]
    pub fn has_errors(&self) -> bool {
        self.diagnostics
            .iter()
            .any(|d| d.severity == Severity::Error)
    }
}

/// One statement, as an expression.
#[derive(Clone, Debug)]
pub struct Statement {
    /// The expression to evaluate, or a statement head: `Set[x, e]`,
    /// `SetDelayed[f[Pattern[x, Blank[]]], body]`, `UpSetDelayed[lhs, rhs]`,
    /// `Assume[facts…]` or `Clear[names…]`.
    pub expr: Expr,
    /// Whether a `;` follows on its line: the session does not print the
    /// result.
    pub quiet: bool,
    /// Where the parts of `expr` come from.
    pub spans: Spans,
}

/// Where the parts of a lowered expression come from in the text.
///
/// Positions follow the `Part` convention: 0 is the head of an application,
/// 1 its first argument. A part that lowering made up, such as the `-1` of
/// `x - y`, has the text of the syntax it was made for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Spans {
    range: Range<usize>,
    head: Option<Box<Spans>>,
    args: Vec<Spans>,
}

impl Spans {
    fn leaf(range: Range<usize>) -> Spans {
        Spans {
            range,
            head: None,
            args: Vec::new(),
        }
    }

    /// The text the whole expression comes from.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// The text the part at `path` comes from, or that of its nearest
    /// enclosing part with a text of its own.
    ///
    /// ```
    /// use cas_lang::{Names, lower, parse};
    ///
    /// let text = "f(x) + 2 y";
    /// let lowered = lower(&parse(text), Names::english());
    /// let spans = &lowered.statements[0].spans;
    /// assert_eq!(&text[spans.find(&[2, 2])], "y");
    /// assert_eq!(&text[spans.find(&[1, 0])], "f");
    /// ```
    #[must_use]
    pub fn find(&self, path: &[usize]) -> Range<usize> {
        let mut spans = self;
        for &i in path {
            let part = if i == 0 {
                spans.head.as_deref()
            } else {
                spans.args.get(i - 1)
            };
            match part {
                Some(part) => spans = part,
                None => break,
            }
        }
        spans.range.clone()
    }
}

/// A problem or a remark about the text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Diagnostic {
    /// The bytes concerned.
    pub range: Range<usize>,
    /// An error leaves its statement out; a note does not.
    pub severity: Severity,
    /// What it is about.
    pub kind: LowerErrorKind,
    /// A replacement that fixes it, when there is an obvious one.
    pub fix: Option<Fix>,
}

/// How serious a [`Diagnostic`] is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Severity {
    /// The statement is left out.
    Error,
    /// A remark: the statement is kept.
    Note,
}

/// A replacement of part of the text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fix {
    /// The bytes to replace.
    pub range: Range<usize>,
    /// What to put instead.
    pub replacement: String,
}

/// What a [`Diagnostic`] is about.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum LowerErrorKind {
    /// `sin x`: a function takes its argument in parentheses (CAS-DOC-7 L3).
    FunctionWithoutParentheses(String),
    /// `sin^2(x)`: a power of a function, which reads as nothing sensible.
    PowerOfFunction(String),
    /// A name of the table assigned, defined or bound (CAS-DOC-7 L7).
    ReservedName(String),
    /// A builtin that the design names but that does not exist yet.
    NotAvailable(String),
    /// `call`, `assume` or `purge` somewhere it has no meaning.
    MisplacedForm(String),
    /// `x :: t` outside the left-hand side of a definition or a rule.
    TypeOutsidePattern,
    /// `xs...` outside the parameters of a definition and the calls of its
    /// body.
    MisplacedSpread,
    /// `x :: t` with a type the table does not name.
    UnknownType(String),
    /// A left-hand side that cannot be assigned or defined.
    InvalidTarget,
    /// `when` after an assignment: guards belong to definitions.
    GuardOnAssignment,
    /// A definition of a builtin without a user name to attach it to, such
    /// as `sin(x :: any) := 0`: builtins are protected.
    NoUserName(String),
    /// A number without an exact value.
    Number(NumberError),
    /// Statements of procedures, which come with stage 2 (CAS-122).
    NotYetSupported(&'static str),
    /// A note: `name(arg)` read as a product.
    ProductReading(String),
}

impl fmt::Display for LowerErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LowerErrorKind::FunctionWithoutParentheses(n) => {
                write!(f, "{n} takes its argument in parentheses: {n}(x)")
            }
            LowerErrorKind::PowerOfFunction(n) => {
                write!(
                    f,
                    "a power of {n} itself: write the power of the value, {n}(x)^2"
                )
            }
            LowerErrorKind::ReservedName(n) => {
                write!(f, "{n} is a name of the language: choose another one")
            }
            LowerErrorKind::NotAvailable(n) => write!(f, "{n} is not available yet"),
            LowerErrorKind::MisplacedForm(n) => match n.as_str() {
                "call" => f.write_str("call is used as call(f, arguments…)"),
                _ => write!(f, "{n}(…) is a statement of its own"),
            },
            LowerErrorKind::TypeOutsidePattern => f.write_str(
                "`::` gives the type of a parameter: it belongs to the left-hand side of a \
                 definition or a rule",
            ),
            LowerErrorKind::MisplacedSpread => f.write_str(
                "`...` belongs to the parameters of a definition, and to the calls of its body",
            ),
            LowerErrorKind::UnknownType(t) => write!(f, "unknown type {t}"),
            LowerErrorKind::InvalidTarget => f.write_str("this cannot be assigned or defined"),
            LowerErrorKind::GuardOnAssignment => {
                f.write_str("`when` belongs to definitions: f(x) := … when …")
            }
            LowerErrorKind::NoUserName(n) => write!(
                f,
                "{n} is a builtin: a definition of it needs a name of yours on its left-hand side"
            ),
            LowerErrorKind::Number(e) => e.fmt(f),
            LowerErrorKind::NotYetSupported(what) => write!(f, "{what} are not supported yet"),
            LowerErrorKind::ProductReading(reading) => write!(f, "read as {reading}"),
        }
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.kind.fmt(f)
    }
}

/// Lowers the statements of a parsed text, with the names of `names`.
///
/// ```
/// use cas_lang::{Names, lower, parse};
///
/// let lowered = lower(&parse("f(x) := x^2\ny := f(3);"), Names::english());
/// let forms: Vec<_> = lowered.statements.iter().map(|s| format!("{:?}", s.expr)).collect();
/// assert_eq!(forms, ["SetDelayed[f[Pattern[x, Blank[]]], Power[x, 2]]", "Set[y, f[3]]"]);
/// assert!(lowered.statements[1].quiet);
/// ```
#[must_use]
pub fn lower(parse: &Parse, names: &Names) -> Lowered {
    let mut lowerer = Lowerer {
        names,
        diagnostics: Vec::new(),
        failed: false,
        mode: Mode::Value,
        sequences: Vec::new(),
    };
    let mut statements = Vec::new();
    for node in parse.syntax().children() {
        let range = range(&node);
        let broken = node.kind() == K::ErrorNode
            || node.descendants().any(|d| d.kind() == K::ErrorNode)
            || parse
                .errors()
                .iter()
                .any(|e| range.start <= e.range.start && e.range.start <= range.end);
        if broken {
            continue;
        }
        lowerer.failed = false;
        let (expr, spans) = lowerer.statement(&node);
        if !lowerer.failed {
            statements.push(Statement {
                expr,
                quiet: followed_by_semicolon(&node),
                spans,
            });
        }
    }
    lowerer.diagnostics.sort_by_key(|d| d.range.start);
    Lowered {
        statements,
        diagnostics: lowerer.diagnostics,
    }
}

/// An expression and where its parts come from.
type Lowering = (Expr, Spans);

/// What names mean in the part being lowered.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mode {
    /// An expression to evaluate.
    Value,
    /// The arguments of a definition: plain names are parameters (§6).
    Parameters,
    /// The left-hand side of a rule or of a definition of a builtin: names
    /// are literal, and `::` marks the variables.
    Literal,
}

struct Lowerer<'a> {
    names: &'a Names,
    diagnostics: Vec<Diagnostic>,
    /// Whether the statement being lowered has an error.
    failed: bool,
    mode: Mode,
    /// The variable-argument parameters of the definition whose body is being
    /// lowered: there, `xs` is the list of the arguments and `xs...` splices
    /// them.
    sequences: Vec<String>,
}

impl Lowerer<'_> {
    fn error(&mut self, range: Range<usize>, kind: LowerErrorKind, fix: Option<Fix>) {
        self.failed = true;
        self.diagnostics.push(Diagnostic {
            range,
            severity: Severity::Error,
            kind,
            fix,
        });
    }

    /// Reports `kind` and returns a placeholder: the statement is left out.
    fn fail(&mut self, range: Range<usize>, kind: LowerErrorKind) -> Lowering {
        self.error(range.clone(), kind, None);
        placeholder(range)
    }

    /// Leaves the statement out without a diagnostic: the parser reported it.
    fn fail_silently(&mut self, range: Range<usize>) -> Lowering {
        self.failed = true;
        placeholder(range)
    }

    /// Runs `lower` with names meaning what `mode` says.
    fn in_mode<T>(&mut self, mode: Mode, lower: impl FnOnce(&mut Self) -> T) -> T {
        let saved = std::mem::replace(&mut self.mode, mode);
        let result = lower(self);
        self.mode = saved;
        result
    }

    fn statement(&mut self, node: &SyntaxNode) -> Lowering {
        match node.kind() {
            K::ExprStmt => match node.first_child() {
                Some(expr) => match self.statement_form(&expr) {
                    Some(lowered) => lowered,
                    None => self.expr(&expr),
                },
                None => self.fail_silently(range(node)),
            },
            K::Definition => self.definition(node),
            _ => self.fail(
                range(node),
                LowerErrorKind::NotYetSupported("local variables, loops and `return`"),
            ),
        }
    }

    /// `assume(…)` and `purge(…)`, the statements the session runs itself.
    fn statement_form(&mut self, node: &SyntaxNode) -> Option<Lowering> {
        let form = self.syntax_form(node)?;
        let args = arguments(node);
        let range = range(node);
        Some(match form {
            SyntaxForm::Assume => {
                let facts = args.iter().map(|a| self.expr(a)).collect();
                apply(Symbol::ASSUME, facts, range)
            }
            SyntaxForm::Purge => {
                let names = args
                    .iter()
                    .map(|a| {
                        if a.kind() == K::Name {
                            self.assignable(a)
                        } else {
                            self.fail(self::range(a), LowerErrorKind::InvalidTarget)
                        }
                    })
                    .collect();
                apply(Symbol::CLEAR, names, range)
            }
            SyntaxForm::Call => return None,
        })
    }

    /// The syntax form an application calls, as in `assume(x > 0)`.
    fn syntax_form(&self, node: &SyntaxNode) -> Option<SyntaxForm> {
        if node.kind() != K::Application {
            return None;
        }
        match self.target(&node.first_child()?)? {
            Target::Syntax(form) => Some(*form),
            _ => None,
        }
    }

    /// What a [`Name`](K::Name) node stands for, if it is in the table.
    fn target(&self, node: &SyntaxNode) -> Option<&Target> {
        if node.kind() != K::Name {
            return None;
        }
        self.names.target(&name_of(node))
    }

    // Definitions and assignments (CAS-DOC-7 §4).

    fn definition(&mut self, node: &SyntaxNode) -> Lowering {
        let mut children = node.children();
        let (Some(lhs), Some(rhs)) = (children.next(), children.next()) else {
            return self.fail_silently(range(node));
        };
        let guard = children.next().and_then(|w| w.first_child());
        let range = range(node);
        let assignment = matches!(lhs.kind(), K::Name | K::ListExpr | K::IndexExpr);
        if assignment {
            if let Some(guard) = &guard {
                self.error(self::range(guard), LowerErrorKind::GuardOnAssignment, None);
            }
        }
        match lhs.kind() {
            K::Name => {
                let target = self.assignable(&lhs);
                let value = self.expr(&rhs);
                apply(Symbol::SET, vec![target, value], range)
            }
            K::ListExpr => {
                let targets = lhs
                    .children()
                    .map(|item| {
                        if item.kind() == K::Name {
                            self.assignable(&item)
                        } else {
                            self.fail(self::range(&item), LowerErrorKind::InvalidTarget)
                        }
                    })
                    .collect();
                let list = apply(Symbol::LIST, targets, self::range(&lhs));
                let value = self.expr(&rhs);
                apply(Symbol::SET, vec![list, value], range)
            }
            K::IndexExpr => {
                let Some(base) = lhs.first_child().filter(|b| b.kind() == K::Name) else {
                    return self.fail(self::range(&lhs), LowerErrorKind::InvalidTarget);
                };
                let target = self.assignable(&base);
                let mut args = vec![target.clone()];
                args.extend(arguments(&lhs).iter().map(|i| self.expr(i)));
                args.push(self.expr(&rhs));
                let copy = apply(Symbol::REPLACE_AT, args, self::range(&rhs));
                apply(Symbol::SET, vec![target, copy], range)
            }
            _ => self.rule_definition(&lhs, &rhs, guard.as_ref(), range),
        }
    }

    /// A name to assign or bind: a user symbol, or an error for a name of
    /// the table.
    fn assignable(&mut self, node: &SyntaxNode) -> Lowering {
        let n = name_of(node);
        let range = range(node);
        if self.names.is_reserved(&n) {
            return self.fail(range, LowerErrorKind::ReservedName(n));
        }
        (Expr::from(user_symbol(&n)), Spans::leaf(range))
    }

    /// `f(params) := body`, a rule on `f`, or a definition of a builtin
    /// (`phi^2 := phi + 1`), an up-rule on the user names in it.
    fn rule_definition(
        &mut self,
        lhs: &SyntaxNode,
        rhs: &SyntaxNode,
        guard: Option<&SyntaxNode>,
        range: Range<usize>,
    ) -> Lowering {
        match self.defined_head(lhs) {
            Defined::User => {
                self.sequences.clear();
                let pattern = self.in_mode(Mode::Parameters, |l| l.defined_application(lhs));
                let body = self.body(rhs, guard);
                self.sequences.clear();
                apply(Symbol::SET_DELAYED, vec![pattern, body], range)
            }
            Defined::Builtin(spelling) => {
                self.sequences.clear();
                let pattern = self.in_mode(Mode::Literal, |l| l.expr(lhs));
                if !pattern.0.args().iter().any(has_user_tag) && !self.failed {
                    self.error(self::range(lhs), LowerErrorKind::NoUserName(spelling), None);
                }
                let body = self.body(rhs, guard);
                self.sequences.clear();
                apply(Symbol::UP_SET_DELAYED, vec![pattern, body], range)
            }
            Defined::Invalid => self.fail(self::range(lhs), LowerErrorKind::InvalidTarget),
            Defined::Failed => placeholder(range),
        }
    }

    /// What the left-hand side of a definition defines.
    fn defined_head(&mut self, lhs: &SyntaxNode) -> Defined {
        match lhs.kind() {
            K::Application => {
                let mut callee = lhs.first_child();
                while let Some(c) = callee.as_ref().filter(|c| c.kind() == K::Application) {
                    callee = c.first_child();
                }
                let Some(callee) = callee.filter(|c| c.kind() == K::Name) else {
                    return Defined::Invalid;
                };
                let n = name_of(&callee);
                match self.names.target(&n) {
                    None => Defined::User,
                    Some(Target::Builtin(s)) if !is_value(*s) => Defined::Builtin(n),
                    Some(Target::Planned(_)) => {
                        self.error(range(&callee), LowerErrorKind::NotAvailable(n), None);
                        Defined::Failed
                    }
                    Some(_) => Defined::Invalid,
                }
            }
            K::BinaryExpr
            | K::ImplicitProduct
            | K::PowerExpr
            | K::PrefixExpr
            | K::Factorial
            | K::Comparison
            | K::Equation => Defined::Builtin(lhs.text().to_string()),
            _ => Defined::Invalid,
        }
    }

    /// The left-hand side of a definition of a user function: the head and
    /// its arguments as patterns, `f(a)(b)` for sub-rules.
    fn defined_application(&mut self, node: &SyntaxNode) -> Lowering {
        let Some(callee) = node.first_child() else {
            return self.fail_silently(range(node));
        };
        let head = if callee.kind() == K::Application {
            self.defined_application(&callee)
        } else {
            let n = name_of(&callee);
            (Expr::from(user_symbol(&n)), Spans::leaf(range(&callee)))
        };
        let args = arguments(node).iter().map(|a| self.parameter(a)).collect();
        apply_to(head, args, range(node))
    }

    /// An argument of a definition (§6): a pattern, `n = 1` for a default,
    /// `xs...` for the remaining arguments.
    fn parameter(&mut self, node: &SyntaxNode) -> Lowering {
        if node.kind() == K::Equation {
            let mut sides = node.children();
            let (Some(variable), Some(default)) = (sides.next(), sides.next()) else {
                return self.fail_silently(range(node));
            };
            let pattern = self.expr(&variable);
            let default = self.in_mode(Mode::Value, |l| l.expr(&default));
            return apply(Symbol::OPTIONAL, vec![pattern, default], range(node));
        }
        self.expr(node)
    }

    /// The right-hand side of a definition, with its guard.
    fn body(&mut self, rhs: &SyntaxNode, guard: Option<&SyntaxNode>) -> Lowering {
        let value = self.expr(rhs);
        match guard {
            Some(guard) => {
                let test = self.expr(guard);
                let range = range(rhs).start..range(guard).end;
                apply(Symbol::CONDITION, vec![value, test], range)
            }
            None => value,
        }
    }

    // Expressions (CAS-DOC-7 §9).

    fn expr(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        match node.kind() {
            K::Literal => self.literal(node),
            K::Name => self.name(node),
            K::Paren => match node.first_child() {
                Some(inner) => self.expr(&inner),
                None => self.fail_silently(range),
            },
            K::BinaryExpr => self.binary(node),
            K::ImplicitProduct => self.product(node),
            K::PrefixExpr => self.prefix(node),
            K::PowerExpr => self.power(node, false),
            K::Equation => {
                let sides = node.children().map(|n| self.expr(&n)).collect();
                apply(Symbol::EQUATION, sides, range)
            }
            K::Comparison => self.comparison(node),
            K::RangeExpr => {
                let ends = node.children().map(|n| self.expr(&n)).collect();
                apply(Symbol::RANGE, ends, range)
            }
            K::Application => self.application(node).0,
            K::IndexExpr => {
                let base = match node.first_child() {
                    Some(b) => self.expr(&b),
                    None => self.fail_silently(range.clone()),
                };
                let mut parts = vec![base];
                parts.extend(arguments(node).iter().map(|i| self.expr(i)));
                apply(Symbol::AT, parts, range)
            }
            K::Factorial => {
                let operand = match node.first_child() {
                    Some(n) => self.expr(&n),
                    None => self.fail_silently(range.clone()),
                };
                apply(Symbol::FACTORIAL, vec![operand], range)
            }
            K::Spread => self.spread(node),
            K::Typed => self.typed(node),
            K::ListExpr => {
                let items = self.items(&node.children().collect::<Vec<_>>());
                apply(Symbol::LIST, items, range)
            }
            K::SetExpr => {
                let items = self.items(&node.children().collect::<Vec<_>>());
                apply(Symbol::FINITE_SET, items, range)
            }
            K::Lambda => self.lambda(node),
            K::IfExpr => self.if_expr(node),
            K::BlockExpr => self.fail(range, LowerErrorKind::NotYetSupported("blocks")),
            _ => self.fail_silently(range),
        }
    }

    fn literal(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let Some(token) = node.first_token() else {
            return self.fail_silently(range);
        };
        match token.kind() {
            K::Number => match number_value(token.text()) {
                Ok(q) => (Expr::from(q), Spans::leaf(range)),
                Err(e) => self.fail(range, LowerErrorKind::Number(e)),
            },
            K::String => match string_value(token.text()) {
                Some(s) => (Expr::string(&s), Spans::leaf(range)),
                None => self.fail_silently(range),
            },
            _ => self.fail_silently(range),
        }
    }

    /// A name: what the table says, or a user symbol; a pattern variable
    /// among the parameters of a definition; the list of the arguments for a
    /// variable-argument parameter in a body.
    fn name(&mut self, node: &SyntaxNode) -> Lowering {
        let n = name_of(node);
        let range = range(node);
        match self.mode {
            Mode::Parameters if !self.names.is_reserved(&n) => {
                let variable = Expr::from(user_symbol(&n));
                (pattern(&variable, Symbol::BLANK), Spans::leaf(range))
            }
            Mode::Value if self.sequences.contains(&n) => {
                let list = Expr::apply(&Symbol::LIST.into(), &[user_symbol(&n).into()]);
                (list, Spans::leaf(range))
            }
            _ => self.name_value(&n, range),
        }
    }

    /// What the name `n` stands for as a value.
    fn name_value(&mut self, n: &str, range: Range<usize>) -> Lowering {
        match self.names.target(n) {
            Some(Target::Builtin(s)) => (Expr::from(*s), Spans::leaf(range)),
            Some(Target::Imaginary) => (build::i(), Spans::leaf(range)),
            Some(Target::Planned(_)) => self.fail(range, LowerErrorKind::NotAvailable(n.into())),
            Some(Target::Syntax(_)) => self.fail(range, LowerErrorKind::MisplacedForm(n.into())),
            None => (Expr::from(user_symbol(n)), Spans::leaf(range)),
        }
    }

    /// `+`, `-`, `*`, `/` and the logical and element-wise operators. A chain
    /// of one associative operator is one node.
    fn binary(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        match operator(node) {
            Some(K::Plus | K::Minus) => {
                let mut terms = Vec::new();
                self.terms(node, &mut terms);
                apply(Symbol::PLUS, terms, range)
            }
            Some(K::Star | K::Slash) => self.product(node),
            Some(op @ (K::AndKw | K::OrKw | K::DotStar)) => {
                let head = match op {
                    K::AndKw => Symbol::AND,
                    K::OrKw => Symbol::OR,
                    _ => Symbol::ELEMENTWISE_TIMES,
                };
                let mut operands = Vec::new();
                self.chain(node, op, &mut operands);
                apply(head, operands, range)
            }
            Some(K::DotSlash) => {
                let operands = node.children().map(|n| self.expr(&n)).collect();
                apply(Symbol::ELEMENTWISE_DIVIDE, operands, range)
            }
            _ => self.fail_silently(range),
        }
    }

    /// The terms of a sum, `a - b` giving `a` and `-b`.
    fn terms(&mut self, node: &SyntaxNode, out: &mut Vec<Lowering>) {
        let op = operator(node);
        if node.kind() == K::BinaryExpr && matches!(op, Some(K::Plus | K::Minus)) {
            let mut sides = node.children();
            if let (Some(lhs), Some(rhs)) = (sides.next(), sides.next()) {
                self.terms(&lhs, out);
                let term = self.expr(&rhs);
                out.push(if op == Some(K::Minus) {
                    negate(term)
                } else {
                    term
                });
                return;
            }
        }
        out.push(self.expr(node));
    }

    /// The operands of a chain of `op`.
    fn chain(&mut self, node: &SyntaxNode, op: K, out: &mut Vec<Lowering>) {
        if node.kind() == K::BinaryExpr && operator(node) == Some(op) {
            for side in node.children() {
                self.chain(&side, op, out);
            }
        } else {
            out.push(self.expr(node));
        }
    }

    /// A product, `*`, `/` and implicit products flattened into one `Times`.
    fn product(&mut self, node: &SyntaxNode) -> Lowering {
        let mut factors = Vec::new();
        self.factors(node, &mut factors);
        apply(Symbol::TIMES, factors, range(node))
    }

    fn factors(&mut self, node: &SyntaxNode, out: &mut Vec<Lowering>) {
        match node.kind() {
            K::BinaryExpr if operator(node) == Some(K::Star) => {
                for side in node.children() {
                    self.factors(&side, out);
                }
            }
            K::BinaryExpr if operator(node) == Some(K::Slash) => {
                let mut sides = node.children();
                if let (Some(lhs), Some(rhs)) = (sides.next(), sides.next()) {
                    self.factors(&lhs, out);
                    let divisor = self.expr(&rhs);
                    let range = range(&rhs);
                    out.push(apply(
                        Symbol::POWER,
                        vec![divisor, (Expr::from(-1), Spans::leaf(range.clone()))],
                        range,
                    ));
                }
            }
            K::ImplicitProduct => self.implicit_factors(node, out),
            K::Application if self.mode == Mode::Value => {
                let ((expr, spans), product) = self.application(node);
                if product && !expr.is_atom() {
                    out.extend(expr.args().iter().cloned().zip(spans.args));
                } else {
                    out.push((expr, spans));
                }
            }
            _ => out.push(self.expr(node)),
        }
    }

    /// The factors of an implicit product, where a function name without
    /// parentheses is an error with a fix (§3.2): `sin x`, `sin^2(x)`.
    fn implicit_factors(&mut self, node: &SyntaxNode, out: &mut Vec<Lowering>) {
        let children: Vec<SyntaxNode> = node.children().collect();
        for (i, child) in children.iter().enumerate() {
            let next = children.get(i + 1);
            if let Some(f) = self.function_name(child) {
                let fix = next.map(|next| Fix {
                    range: range(child).start..range(next).end,
                    replacement: format!("{}({})", child.text(), unparenthesised(next)),
                });
                self.error(
                    range(child),
                    LowerErrorKind::FunctionWithoutParentheses(f),
                    fix,
                );
                out.push(placeholder(range(child)));
            } else if let Some(f) = self.function_power(child) {
                let fix = next.and_then(|next| {
                    let exponent = exponent_text(child)?;
                    Some(Fix {
                        range: range(child).start..range(next).end,
                        replacement: format!("{f}({})^{exponent}", unparenthesised(next)),
                    })
                });
                self.error(range(child), LowerErrorKind::PowerOfFunction(f), fix);
                out.push(placeholder(range(child)));
            } else {
                self.factors(child, out);
            }
        }
    }

    /// The function of the table that a power raises, as in `sin^2`.
    fn function_power(&self, node: &SyntaxNode) -> Option<String> {
        if node.kind() != K::PowerExpr {
            return None;
        }
        self.function_name(&node.first_child()?)
    }

    /// The spelling of a function of the table that `node` names, as in
    /// `sin`; `None` for other names, constants and other nodes.
    fn function_name(&self, node: &SyntaxNode) -> Option<String> {
        match self.target(node)? {
            Target::Builtin(s) if !is_value(*s) => Some(node.text().to_string()),
            Target::Planned(_) => Some(node.text().to_string()),
            _ => None,
        }
    }

    /// `-x`, `+x`, `√x`, `not p`, and the directed infinities.
    fn prefix(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let Some(operand) = node.first_child() else {
            return self.fail_silently(range);
        };
        let op = operator(node);
        // `+infinity` and `-infinity` are whole forms (§5).
        if matches!(op, Some(K::Plus | K::Minus))
            && matches!(
                self.target(&operand),
                Some(Target::Builtin(Symbol::COMPLEX_INFINITY))
            )
        {
            let expr = if op == Some(K::Plus) {
                Expr::from(Symbol::INFINITY)
            } else {
                Expr::apply(&Symbol::DIRECTED_INFINITY.into(), &[Expr::from(-1)])
            };
            return (expr, Spans::leaf(range));
        }
        let value = self.expr(&operand);
        match op {
            Some(K::Minus) => negate(value),
            Some(K::Plus) => value,
            Some(K::NotKw) => apply(Symbol::NOT, vec![value], range),
            Some(root @ (K::Sqrt | K::Cbrt | K::FourthRoot)) => {
                let degree = match root {
                    K::Sqrt => 2,
                    K::Cbrt => 3,
                    _ => 4,
                };
                let exponent = Expr::from(cas_core::Rational::new(1.into(), degree.into()));
                apply(
                    Symbol::POWER,
                    vec![value, (exponent, Spans::leaf(range.clone()))],
                    range,
                )
            }
            _ => self.fail_silently(range),
        }
    }

    /// `b^e`, `b.^e` and `b²`. With `checked`, the caller has dealt with a
    /// function as the base.
    fn power(&mut self, node: &SyntaxNode, checked: bool) -> Lowering {
        let range = range(node);
        let Some(base) = node.first_child() else {
            return self.fail_silently(range);
        };
        if !checked {
            if let Some(f) = self.function_name(&base) {
                return self.fail(range, LowerErrorKind::PowerOfFunction(f));
            }
        }
        let base = self.expr(&base);
        let op = significant_tokens(node).next();
        match op.as_ref().map(SyntaxToken::kind) {
            Some(K::Superscript) => {
                let token = op.expect("a superscript");
                let Some(n) = superscript_value(token.text()) else {
                    return self.fail_silently(range);
                };
                let exponent = (Expr::from(n), Spans::leaf(token_range(&token)));
                apply(Symbol::POWER, vec![base, exponent], range)
            }
            Some(kind @ (K::Caret | K::DotCaret)) => {
                let Some(exponent) = node.children().nth(1) else {
                    return self.fail_silently(range);
                };
                let exponent = self.expr(&exponent);
                let head = if kind == K::Caret {
                    Symbol::POWER
                } else {
                    Symbol::ELEMENTWISE_POWER
                };
                apply(head, vec![base, exponent], range)
            }
            _ => self.fail_silently(range),
        }
    }

    /// A chain of comparisons: one test, or `And` of the tests pairwise
    /// (`0 < x < 1` is `0 < x and x < 1`).
    fn comparison(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let operands: Vec<SyntaxNode> = node.children().collect();
        let ops: Vec<K> = significant_tokens(node).map(|t| t.kind()).collect();
        let lowered: Vec<Lowering> = operands.iter().map(|n| self.expr(n)).collect();
        let mut tests: Vec<Lowering> = Vec::new();
        for (i, op) in ops.iter().enumerate() {
            let (Some(lhs), Some(rhs)) = (lowered.get(i), lowered.get(i + 1)) else {
                return self.fail_silently(range);
            };
            let head = match op {
                K::EqEq => Symbol::EQUAL,
                K::Neq => Symbol::UNEQUAL,
                K::Lt => Symbol::LESS,
                K::Le => Symbol::LESS_EQUAL,
                K::Gt => Symbol::GREATER,
                K::Ge => Symbol::GREATER_EQUAL,
                _ => Symbol::ELEMENT,
            };
            let span = lhs.1.range.start..rhs.1.range.end;
            tests.push(apply(head, vec![lhs.clone(), rhs.clone()], span));
        }
        if tests.len() == 1 {
            tests.pop().expect("one test")
        } else {
            apply(Symbol::AND, tests, range)
        }
    }

    /// An application (§3.1): a call, or a product. Returns whether it reads
    /// as a product.
    fn application(&mut self, node: &SyntaxNode) -> (Lowering, bool) {
        let range = range(node);
        let Some(callee) = node.first_child() else {
            return (self.fail_silently(range), false);
        };
        let args = arguments(node);
        if self.mode != Mode::Value {
            // In a pattern, an application is structure: its head is literal.
            let head = self.head(&callee);
            let args = args.iter().map(|a| self.argument(a)).collect();
            return (apply_to(head, args, range), false);
        }
        match callee.kind() {
            K::Name => {
                let n = name_of(&callee);
                let target = self.names.target(&n).cloned();
                match target {
                    Some(Target::Syntax(SyntaxForm::Call)) => {
                        let Some((head, rest)) = args.split_first() else {
                            return (self.fail(range, LowerErrorKind::MisplacedForm(n)), false);
                        };
                        let head = self.expr(head);
                        let rest = rest.iter().map(|a| self.argument(a)).collect();
                        (apply_to(head, rest, range), false)
                    }
                    Some(Target::Syntax(_)) => {
                        (self.fail(range, LowerErrorKind::MisplacedForm(n)), false)
                    }
                    Some(ref t) if is_constant(t) && args.len() == 1 => {
                        (self.reading_product(&callee, &args[0], range, false), true)
                    }
                    _ if args.len() == 1 && occurs_in(&n, &args[0]) => {
                        (self.reading_product(&callee, &args[0], range, true), true)
                    }
                    Some(Target::Planned(_)) => {
                        (self.fail(range, LowerErrorKind::NotAvailable(n)), false)
                    }
                    Some(Target::Builtin(Symbol::REPLACE_ALL)) => {
                        let head = (
                            Expr::from(Symbol::REPLACE_ALL),
                            Spans::leaf(self::range(&callee)),
                        );
                        let mut lowered = Vec::new();
                        for (i, a) in args.iter().enumerate() {
                            lowered.push(if i == 0 { self.expr(a) } else { self.rules(a) });
                        }
                        (apply_to(head, lowered, range), false)
                    }
                    _ => {
                        let head = self.expr(&callee);
                        let args = args.iter().map(|a| self.argument(a)).collect();
                        (apply_to(head, args, range), false)
                    }
                }
            }
            K::Application => {
                let (inner, product) = self.application(&callee);
                if product && args.len() == 1 {
                    let factor = self.expr(&args[0]);
                    (apply(Symbol::TIMES, vec![inner, factor], range), true)
                } else {
                    let args = args.iter().map(|a| self.argument(a)).collect();
                    (apply_to(inner, args, range), false)
                }
            }
            _ => {
                let head = self.expr(&callee);
                let args = args.iter().map(|a| self.argument(a)).collect();
                (apply_to(head, args, range), false)
            }
        }
    }

    /// `name(arg)` read as `name*arg`, with a note when the reading comes
    /// from the name occurring in its argument.
    fn reading_product(
        &mut self,
        callee: &SyntaxNode,
        arg: &SyntaxNode,
        range: Range<usize>,
        note: bool,
    ) -> Lowering {
        if note {
            let reading = format!("{}*({})", callee.text(), arg.text());
            self.diagnostics.push(Diagnostic {
                range: range.clone(),
                severity: Severity::Note,
                kind: LowerErrorKind::ProductReading(reading),
                fix: None,
            });
        }
        let factor = self.expr(callee);
        let argument = self.expr(arg);
        apply(Symbol::TIMES, vec![factor, argument], range)
    }

    /// The head of an application in a pattern: a name is literal.
    fn head(&mut self, callee: &SyntaxNode) -> Lowering {
        if callee.kind() == K::Name {
            let n = name_of(callee);
            self.name_value(&n, range(callee))
        } else {
            self.expr(callee)
        }
    }

    /// An argument of a call: `xs...` splices a variable-argument
    /// parameter; in a pattern, `xs...` is one.
    fn argument(&mut self, node: &SyntaxNode) -> Lowering {
        if node.kind() == K::Spread {
            return self.spread(node);
        }
        if self.mode == Mode::Parameters {
            return self.parameter(node);
        }
        self.expr(node)
    }

    /// `xs...`: in a pattern, a variable for any number of arguments; in a
    /// body, the arguments it stands for, spliced.
    fn spread(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let Some(variable) = node.first_child().filter(|v| v.kind() == K::Name) else {
            return self.fail(range, LowerErrorKind::MisplacedSpread);
        };
        let n = name_of(&variable);
        match self.mode {
            Mode::Parameters | Mode::Literal => {
                if self.names.is_reserved(&n) {
                    return self.fail(range, LowerErrorKind::ReservedName(n));
                }
                self.sequences.push(n.clone());
                let variable = Expr::from(user_symbol(&n));
                (
                    pattern(&variable, Symbol::BLANK_NULL_SEQUENCE),
                    Spans::leaf(range),
                )
            }
            Mode::Value if self.sequences.contains(&n) => {
                (Expr::from(user_symbol(&n)), Spans::leaf(range))
            }
            Mode::Value => self.fail(range, LowerErrorKind::MisplacedSpread),
        }
    }

    /// `x :: t`, a pattern variable with a type.
    fn typed(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        if self.mode == Mode::Value {
            return self.fail(range, LowerErrorKind::TypeOutsidePattern);
        }
        let mut children = node.children();
        let (Some(variable), Some(ty)) = (children.next(), children.next()) else {
            return self.fail_silently(range);
        };
        if variable.kind() != K::Name {
            return self.fail(self::range(&variable), LowerErrorKind::InvalidTarget);
        }
        let n = name_of(&variable);
        if self.names.is_reserved(&n) {
            return self.fail(self::range(&variable), LowerErrorKind::ReservedName(n));
        }
        let type_name = name_of(&ty);
        let Some(ty) = self.names.type_name(&type_name) else {
            return self.fail(self::range(&ty), LowerErrorKind::UnknownType(type_name));
        };
        let blank = pattern(&Expr::from(user_symbol(&n)), Symbol::BLANK);
        let expr = match type_test(ty) {
            None => blank,
            Some(test) => Expr::apply(&Symbol::PATTERN_TEST.into(), &[blank, test.into()]),
        };
        (expr, Spans::leaf(range))
    }

    /// The rules of `replace` and `subst` (§6): an equation is a rule whose
    /// left-hand side is a pattern; a list holds several.
    fn rules(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        match node.kind() {
            K::Equation => {
                let mut sides = node.children();
                let (Some(lhs), Some(rhs)) = (sides.next(), sides.next()) else {
                    return self.fail_silently(range);
                };
                let lhs = self.in_mode(Mode::Literal, |l| l.expr(&lhs));
                let rhs = self.expr(&rhs);
                apply(Symbol::RULE_DELAYED, vec![lhs, rhs], range)
            }
            K::ListExpr => {
                let rules = node.children().map(|n| self.rules(&n)).collect();
                apply(Symbol::LIST, rules, range)
            }
            _ => self.expr(node),
        }
    }

    /// The items of a list or a set; `xs...` splices in a body.
    fn items(&mut self, nodes: &[SyntaxNode]) -> Vec<Lowering> {
        nodes.iter().map(|n| self.argument(n)).collect()
    }

    /// `x -> body`, `(x, y) -> body`: a `Function` of its parameters.
    fn lambda(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let mut children = node.children();
        let (Some(params), Some(body)) = (children.next(), children.next()) else {
            return self.fail_silently(range);
        };
        let names: Vec<SyntaxNode> = if params.kind() == K::Name {
            vec![params.clone()]
        } else {
            params.children().collect()
        };
        let parameters = names.iter().map(|n| self.assignable(n)).collect();
        let parameters = apply(Symbol::LIST, parameters, self::range(&params));
        let body = self.expr(&body);
        apply(Symbol::FUNCTION, vec![parameters, body], range)
    }

    /// `if c then a elif d then b else e end`: nested `If`s. In stage 1 a
    /// branch is one expression.
    fn if_expr(&mut self, node: &SyntaxNode) -> Lowering {
        let range = range(node);
        let nodes: Vec<SyntaxNode> = node.children().collect();
        let has_else = node.children_with_tokens().any(|e| e.kind() == K::ElseKw);
        let mut branches = Vec::new();
        let mut i = 0;
        while i + 1 < nodes.len() && nodes[i].kind() != K::Body {
            branches.push((nodes[i].clone(), nodes[i + 1].clone()));
            i += 2;
        }
        let mut result = if has_else {
            nodes.get(i).map(|b| self.branch(b))
        } else {
            None
        };
        for (condition, body) in branches.iter().rev() {
            let condition = self.expr(condition);
            let value = self.branch(body);
            let span = condition.1.range.start..range.end;
            let mut args = vec![condition, value];
            args.extend(result.take());
            result = Some(apply(Symbol::IF, args, span));
        }
        result.unwrap_or_else(|| self.fail_silently(range))
    }

    fn branch(&mut self, body: &SyntaxNode) -> Lowering {
        let statements: Vec<SyntaxNode> = body.children().collect();
        match statements.as_slice() {
            [] => (Expr::from(Symbol::NULL), Spans::leaf(range(body))),
            [statement] if statement.kind() == K::ExprStmt => match statement.first_child() {
                Some(expr) => self.expr(&expr),
                None => self.fail_silently(range(body)),
            },
            _ => self.fail(
                range(body),
                LowerErrorKind::NotYetSupported("branches of several statements"),
            ),
        }
    }
}

/// What a definition defines.
enum Defined {
    /// A function of the user's.
    User,
    /// A builtin, spelled as in the text: an up-rule.
    Builtin(String),
    /// Nothing that can be defined.
    Invalid,
    /// Something reported already.
    Failed,
}

/// The name a [`Name`](K::Name) node stands for, as [`name`] gives it.
fn name_of(node: &SyntaxNode) -> String {
    name(&node.text().to_string()).into_owned()
}

/// The arguments of an application or an index.
fn arguments(node: &SyntaxNode) -> Vec<SyntaxNode> {
    node.children()
        .find(|c| c.kind() == K::ArgList)
        .map(|args| args.children().collect())
        .unwrap_or_default()
}

fn range(node: &SyntaxNode) -> Range<usize> {
    let r = node.text_range();
    usize::from(r.start())..usize::from(r.end())
}

fn token_range(token: &SyntaxToken) -> Range<usize> {
    let r = token.text_range();
    usize::from(r.start())..usize::from(r.end())
}

/// The tokens of `node` itself that are not trivia: its operators.
fn significant_tokens(node: &SyntaxNode) -> impl Iterator<Item = SyntaxToken> {
    node.children_with_tokens()
        .filter_map(NodeOrToken::into_token)
        .filter(|t| !t.kind().is_trivia() && t.kind() != K::Newline)
}

/// The operator of a binary or prefix node.
fn operator(node: &SyntaxNode) -> Option<K> {
    significant_tokens(node).next().map(|t| t.kind())
}

/// `head[args…]`, built structurally.
fn apply(head: Symbol, args: Vec<Lowering>, range: Range<usize>) -> Lowering {
    let (exprs, spans): (Vec<Expr>, Vec<Spans>) = args.into_iter().unzip();
    let expr = Expr::apply(&head.into(), &exprs);
    (
        expr,
        Spans {
            range,
            head: None,
            args: spans,
        },
    )
}

/// `head[args…]` with a head from the text.
fn apply_to(head: Lowering, args: Vec<Lowering>, range: Range<usize>) -> Lowering {
    let (exprs, spans): (Vec<Expr>, Vec<Spans>) = args.into_iter().unzip();
    let expr = Expr::apply(&head.0, &exprs);
    (
        expr,
        Spans {
            range,
            head: Some(Box::new(head.1)),
            args: spans,
        },
    )
}

/// `-value`: the negated number for a number, `Times[-1, value]` otherwise.
fn negate((expr, spans): Lowering) -> Lowering {
    if expr.is_number() {
        return (-&expr, spans);
    }
    let range = spans.range.clone();
    apply(
        Symbol::TIMES,
        vec![(Expr::from(-1), Spans::leaf(range.clone())), (expr, spans)],
        range,
    )
}

/// `Pattern[variable, blank[]]`.
fn pattern(variable: &Expr, blank: Symbol) -> Expr {
    let blank = Expr::apply(&blank.into(), &[]);
    Expr::apply(&Symbol::PATTERN.into(), &[variable.clone(), blank])
}

/// The symbol of a user's name: `$`-escaped when capitalised (§5).
fn user_symbol(n: &str) -> Symbol {
    Symbol::new(&escape_name(n))
}

/// What stands for a statement or an expression that has an error.
fn placeholder(range: Range<usize>) -> Lowering {
    (Expr::from(Symbol::NULL), Spans::leaf(range))
}

/// Whether a builtin is a value rather than a function: a constant, a
/// logical value or an infinity.
pub(crate) fn is_value(s: Symbol) -> bool {
    s.builtin_attributes()
        .is_some_and(|a| a.contains(Attributes::CONSTANT))
        || matches!(
            s,
            Symbol::TRUE
                | Symbol::FALSE
                | Symbol::NULL
                | Symbol::INDETERMINATE
                | Symbol::INFINITY
                | Symbol::COMPLEX_INFINITY
        )
}

/// Whether a target is a constant, before `(…)` a factor (§3.1).
fn is_constant(target: &Target) -> bool {
    match target {
        Target::Imaginary => true,
        Target::Builtin(s) => is_value(*s),
        _ => false,
    }
}

/// Whether `n` names something in `arg` other than the head of a call
/// (giac's rule, §3.1): `x` in `x + 1`, not `f` in `f(y)`.
fn occurs_in(n: &str, arg: &SyntaxNode) -> bool {
    arg.descendants()
        .filter(|d| d.kind() == K::Name)
        .filter(|d| {
            !d.parent()
                .is_some_and(|p| p.kind() == K::Application && p.first_child().as_ref() == Some(d))
        })
        .any(|d| name_of(&d) == n)
}

/// Whether a left-hand side holds a user's name at level 1, where an
/// up-rule attaches (CAS-DOC-4 §5.2).
fn has_user_tag(arg: &Expr) -> bool {
    let symbol = match arg.view() {
        View::Symbol(s) => Some(s),
        View::Apply { head, .. } => head.as_symbol(),
        _ => None,
    };
    symbol.is_some_and(|s| !s.is_builtin())
}

/// The predicate a type tests with; `None` for `any`.
fn type_test(ty: TypeName) -> Option<Symbol> {
    Some(match ty {
        TypeName::Any => return None,
        TypeName::Integer => Symbol::INTEGER_Q,
        TypeName::Rational => Symbol::RATIONAL_Q,
        TypeName::Number => Symbol::NUMBER_Q,
        TypeName::Real => Symbol::REAL_Q,
        TypeName::Complex => Symbol::COMPLEX_Q,
        TypeName::List => Symbol::LIST_Q,
        TypeName::Matrix => Symbol::MATRIX_Q,
        TypeName::String => Symbol::STRING_Q,
        TypeName::Symbol => Symbol::SYMBOL_Q,
    })
}

/// The text of a factor, without its parentheses if it has some.
fn unparenthesised(node: &SyntaxNode) -> String {
    match node.kind() {
        K::Paren => node
            .first_child()
            .map_or_else(|| node.text().to_string(), |inner| inner.text().to_string()),
        _ => node.text().to_string(),
    }
}

/// The exponent of a power as plain text: `2` for `^2` and for `²`.
fn exponent_text(power: &SyntaxNode) -> Option<String> {
    let op = significant_tokens(power).next()?;
    if op.kind() == K::Superscript {
        superscript_value(op.text()).map(|n| n.to_string())
    } else {
        power.children().nth(1).map(|e| unparenthesised(&e))
    }
}

/// Whether a `;` follows `node` on its line.
fn followed_by_semicolon(node: &SyntaxNode) -> bool {
    let mut next = node.next_sibling_or_token();
    while let Some(element) = next {
        match element {
            NodeOrToken::Token(t) if t.kind().is_trivia() => next = t.next_sibling_or_token(),
            NodeOrToken::Token(t) => return t.kind() == K::Semicolon,
            NodeOrToken::Node(_) => return false,
        }
    }
    false
}

#[cfg(test)]
mod tests;
