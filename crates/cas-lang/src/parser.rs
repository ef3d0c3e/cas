//! The parser: text to a lossless syntax tree (CAS-DOC-7 §2–§4, §10).
//!
//! Hand-written: one method per rule of the grammar for statements and
//! blocks, and a Pratt loop over binding powers for the operators
//! (`grammar.rs`). It builds a `rowan` green tree as it goes, with
//! checkpoints to wrap an operand once its operator is seen. Every token of
//! the text is in the tree, trivia included: a node spans from its first to
//! its last significant token, and the trivia before a token belongs to the
//! node that is open when the token is reached.
//!
//! Parsing never fails. A token that has no place goes into an
//! [`ErrorNode`](SyntaxKind::ErrorNode), each problem is a [`SyntaxError`],
//! and the lexer's problems come along. Recovery stays local: a statement
//! that goes wrong is skipped to the end of its line, an argument to the next
//! comma or the closing bracket. Text that nests beyond [`MAX_DEPTH`] goes
//! whole into one error node, so that the parser's recursion stays within
//! any thread's stack.

use std::fmt;
use std::ops::Range;

use rowan::{Checkpoint, GreenNode, GreenNodeBuilder};

use crate::lexer::{LexErrorKind, Lexed, lex};
use crate::syntax::{SyntaxKind, SyntaxNode};

mod grammar;

/// How deeply expressions and blocks may nest: parentheses within
/// parentheses, operators within operators' operands, blocks within blocks.
///
/// Nested blocks take the most stack per level. On a 2 MiB stack in a debug
/// build they overflow between 1000 and 2000 levels, so this bound leaves at
/// least a fourfold margin on any thread.
pub const MAX_DEPTH: usize = 256;

/// The syntax tree of a text, with the problems found in it.
#[derive(Clone, Debug)]
pub struct Parse {
    green: GreenNode,
    errors: Vec<SyntaxError>,
}

impl Parse {
    /// The root of the tree, a [`Root`](SyntaxKind::Root) node. Its text is
    /// the whole text that was parsed.
    #[must_use]
    pub fn syntax(&self) -> SyntaxNode {
        SyntaxNode::new_root(self.green.clone())
    }

    /// The problems, from the lexer and the parser, ordered by where they
    /// start.
    #[must_use]
    pub fn errors(&self) -> &[SyntaxError] {
        &self.errors
    }
}

/// Parses `text` into a syntax tree.
///
/// ```
/// use cas_lang::{SyntaxKind, parse};
///
/// let tree = parse("f(x) := x^2 when x > 0\ny := f(3)").syntax();
/// let kinds: Vec<_> = tree.children().map(|n| n.kind()).collect();
/// assert_eq!(kinds, [SyntaxKind::Definition, SyntaxKind::Definition]);
/// assert_eq!(tree.text(), "f(x) := x^2 when x > 0\ny := f(3)");
/// ```
#[must_use]
pub fn parse(text: &str) -> Parse {
    let lexed = lex(text);
    let mut parser = Parser::new(&lexed);
    parser.root();
    parser.finish()
}

/// A problem in the text, found by the lexer or by the parser. `Display`
/// gives a message for the user.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SyntaxError {
    /// The bytes concerned; empty at the end of the text.
    pub range: Range<usize>,
    /// What is wrong.
    pub kind: SyntaxErrorKind,
}

/// What is wrong in a [`SyntaxError`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SyntaxErrorKind {
    /// Found by the lexer.
    Lex(LexErrorKind),
    /// Found by the parser.
    Parse(ParseErrorKind),
}

/// What the parser found wrong.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ParseErrorKind {
    /// An expression was expected here.
    ExpectedExpression,
    /// This token was expected here.
    Expected(SyntaxKind),
    /// A name was expected here: a parameter, a `for` variable, a local or a
    /// type.
    ExpectedName,
    /// Something follows a complete statement on its line.
    ExpectedSeparator,
    /// A token that has no place here.
    Unexpected(SyntaxKind),
    /// `a = b = c`: an equation has two sides, and does not chain with
    /// comparisons.
    ChainedEquation,
    /// `a..b..c`.
    ChainedRange,
    /// `(a, b)` outside the parameters of a lambda.
    Tuple,
    /// `x 2`: a number cannot continue an implicit product.
    NumberInProduct,
    /// Nesting beyond [`MAX_DEPTH`].
    TooDeep,
}

impl fmt::Display for ParseErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseErrorKind::ExpectedExpression => f.write_str("expected an expression"),
            ParseErrorKind::Expected(kind) => write!(f, "expected {}", kind.describe()),
            ParseErrorKind::ExpectedName => f.write_str("expected a name"),
            ParseErrorKind::ExpectedSeparator => {
                f.write_str("expected `;` or a line break after the statement")
            }
            ParseErrorKind::Unexpected(kind) => write!(f, "unexpected {}", kind.describe()),
            ParseErrorKind::ChainedEquation => f.write_str(
                "an equation has two sides: `a = b = c` and `a = b < c` need `and` between them",
            ),
            ParseErrorKind::ChainedRange => f.write_str("a range has two ends: `a..b..c`"),
            ParseErrorKind::Tuple => f.write_str("there are no tuples: write a list, [a, b]"),
            ParseErrorKind::NumberInProduct => {
                f.write_str("a number cannot follow a factor: write `*` before it")
            }
            ParseErrorKind::TooDeep => write!(f, "nested more than {MAX_DEPTH} deep"),
        }
    }
}

impl fmt::Display for SyntaxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            SyntaxErrorKind::Lex(kind) => kind.fmt(f),
            SyntaxErrorKind::Parse(kind) => kind.fmt(f),
        }
    }
}

impl std::error::Error for SyntaxError {}

/// The state of one parse.
struct Parser<'a> {
    lexed: &'a Lexed<'a>,
    /// The next token to add to the tree, trivia included.
    pos: usize,
    /// The next significant token: neither trivia nor a line break.
    next: usize,
    /// Whether a line break lies between `pos` and `next`.
    newline: bool,
    builder: GreenNodeBuilder<'static>,
    errors: Vec<SyntaxError>,
    /// Whether line breaks separate statements, innermost context last:
    /// they do in blocks, not inside brackets (CAS-DOC-7 §4).
    newlines: Vec<bool>,
    depth: usize,
    /// Set when the text nests too deeply; the rest of it is in an error
    /// node, and the parser sees the end of the text.
    stopped: bool,
}

impl<'a> Parser<'a> {
    fn new(lexed: &'a Lexed<'a>) -> Parser<'a> {
        let mut parser = Parser {
            lexed,
            pos: 0,
            next: 0,
            newline: false,
            builder: GreenNodeBuilder::new(),
            errors: Vec::new(),
            newlines: Vec::new(),
            depth: 0,
            stopped: false,
        };
        parser.seek();
        parser
    }

    /// Adds the rest of the text, closes the root and collects the errors.
    fn finish(mut self) -> Parse {
        self.add_tokens(self.lexed.len());
        self.builder.finish_node();
        let mut errors: Vec<SyntaxError> = self
            .lexed
            .errors()
            .iter()
            .map(|e| SyntaxError {
                range: e.range.clone(),
                kind: SyntaxErrorKind::Lex(e.kind),
            })
            .collect();
        errors.extend(self.errors);
        errors.sort_by_key(|e| e.range.start);
        Parse {
            green: self.builder.finish(),
            errors,
        }
    }

    /// Finds the next significant token after `pos`.
    fn seek(&mut self) {
        self.next = self.pos;
        self.newline = false;
        while self.next < self.lexed.len() {
            match self.lexed.kind(self.next) {
                SyntaxKind::Newline => self.newline = true,
                kind if kind.is_trivia() => {}
                _ => break,
            }
            self.next += 1;
        }
    }

    /// The kind of the next significant token, or [`SyntaxKind::Eof`].
    fn peek(&self) -> SyntaxKind {
        if self.stopped || self.next >= self.lexed.len() {
            SyntaxKind::Eof
        } else {
            self.lexed.kind(self.next)
        }
    }

    fn at(&self, kind: SyntaxKind) -> bool {
        self.peek() == kind
    }

    /// The kinds of the significant tokens from the next one on.
    fn lookahead(&self) -> impl Iterator<Item = SyntaxKind> + '_ {
        let end = if self.stopped { 0 } else { self.lexed.len() };
        (self.next..end)
            .map(|i| self.lexed.kind(i))
            .filter(|&kind| kind != SyntaxKind::Newline && !kind.is_trivia())
    }

    /// Whether a line break that separates statements comes before the next
    /// significant token.
    fn newline_before(&self) -> bool {
        self.newline && self.newlines.last().copied().unwrap_or(true)
    }

    /// The next token's kind if it continues what is being parsed: on the
    /// same line, or after a line break when it cannot start a statement
    /// (CAS-DOC-7 §4). Otherwise [`SyntaxKind::Eof`].
    fn peek_continuing(&self) -> SyntaxKind {
        let kind = self.peek();
        if self.newline_before() && starts_statement(kind) {
            SyntaxKind::Eof
        } else {
            kind
        }
    }

    /// Adds the tokens before `end` to the tree.
    fn add_tokens(&mut self, end: usize) {
        while self.pos < end {
            let kind = self.lexed.kind(self.pos);
            self.builder.token(kind.into(), self.lexed.text(self.pos));
            self.pos += 1;
        }
    }

    /// Adds the trivia and line breaks before the next significant token.
    fn flush_trivia(&mut self) {
        self.add_tokens(self.next);
    }

    /// Adds the next significant token, and the trivia before it.
    fn bump(&mut self) {
        if self.at(SyntaxKind::Eof) {
            return;
        }
        self.add_tokens(self.next + 1);
        self.seek();
    }

    fn start(&mut self, kind: SyntaxKind) {
        self.flush_trivia();
        self.builder.start_node(kind.into());
    }

    fn checkpoint(&mut self) -> Checkpoint {
        self.flush_trivia();
        self.builder.checkpoint()
    }

    /// Starts a node of `kind` at `checkpoint`, around what was added since.
    fn start_at(&mut self, checkpoint: Checkpoint, kind: SyntaxKind) {
        self.builder.start_node_at(checkpoint, kind.into());
    }

    fn finish_node(&mut self) {
        self.builder.finish_node();
    }

    /// A node of `kind` holding the next token.
    fn single(&mut self, kind: SyntaxKind) {
        self.start(kind);
        self.bump();
        self.finish_node();
    }

    /// Adds the next token if it is `kind`; reports it missing otherwise.
    fn expect(&mut self, kind: SyntaxKind) -> bool {
        let found = self.at(kind);
        if found {
            self.bump();
        } else {
            self.error(ParseErrorKind::Expected(kind));
        }
        found
    }

    /// Reports `kind` at the next significant token, or at the end.
    fn error(&mut self, kind: ParseErrorKind) {
        if self.stopped {
            return;
        }
        let range = if self.next < self.lexed.len() {
            self.lexed.range(self.next)
        } else {
            let end = self.lexed.source().len();
            end..end
        };
        self.errors.push(SyntaxError {
            range,
            kind: SyntaxErrorKind::Parse(kind),
        });
    }

    /// Runs `parse` in a context where line breaks separate statements or
    /// not.
    fn with_newlines<T>(&mut self, separate: bool, parse: impl FnOnce(&mut Self) -> T) -> T {
        self.newlines.push(separate);
        let result = parse(self);
        self.newlines.pop();
        result
    }

    /// Enters one level of nesting, or stops the parse when it is one too
    /// many: the rest of the text goes into an error node.
    fn enter(&mut self) -> bool {
        if self.depth < MAX_DEPTH {
            self.depth += 1;
            return true;
        }
        if !self.stopped {
            self.error(ParseErrorKind::TooDeep);
            self.flush_trivia();
            self.builder.start_node(SyntaxKind::ErrorNode.into());
            self.add_tokens(self.lexed.len());
            self.builder.finish_node();
            self.stopped = true;
            self.seek();
        }
        false
    }

    fn leave(&mut self) {
        self.depth -= 1;
    }
}

/// Whether a token of `kind` starts an expression.
fn starts_expression(kind: SyntaxKind) -> bool {
    matches!(
        kind,
        SyntaxKind::Number
            | SyntaxKind::String
            | SyntaxKind::Ident
            | SyntaxKind::LParen
            | SyntaxKind::LBracket
            | SyntaxKind::LBrace
            | SyntaxKind::Minus
            | SyntaxKind::Plus
            | SyntaxKind::Sqrt
            | SyntaxKind::Cbrt
            | SyntaxKind::FourthRoot
            | SyntaxKind::NotKw
            | SyntaxKind::IfKw
            | SyntaxKind::BeginKw
            | SyntaxKind::Error
    )
}

/// Whether a token of `kind` starts a statement; one that does not continues
/// the statement before it across a line break (CAS-DOC-7 §4).
fn starts_statement(kind: SyntaxKind) -> bool {
    starts_expression(kind)
        || matches!(
            kind,
            SyntaxKind::LocalKw
                | SyntaxKind::ReturnKw
                | SyntaxKind::BreakKw
                | SyntaxKind::ContinueKw
                | SyntaxKind::ForKw
                | SyntaxKind::WhileKw
        )
}

#[cfg(test)]
mod tests;
