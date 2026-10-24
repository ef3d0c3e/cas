//! The rules of the grammar (CAS-DOC-7 §2): one method per rule for
//! statements and blocks, and a Pratt loop over the binding powers of §3 for
//! the operators.

use rowan::Checkpoint;

use super::{ParseErrorKind as E, Parser, starts_expression, starts_statement};
use crate::syntax::SyntaxKind as K;

/// Binding powers, higher binding tighter (CAS-DOC-7 §3). An infix operator
/// applies when its left power is at least the context's minimum, and parses
/// its right operand with its right power as the minimum: a left power below
/// the right one makes it left-associative.
mod bp {
    pub(super) const OR: (u8, u8) = (1, 2);
    pub(super) const AND: (u8, u8) = (3, 4);
    /// The operand of `not`: everything that binds tighter than `and`.
    pub(super) const NOT: u8 = 5;
    pub(super) const COMPARE: (u8, u8) = (7, 8);
    pub(super) const RANGE: (u8, u8) = (9, 10);
    pub(super) const SUM: (u8, u8) = (11, 12);
    pub(super) const PRODUCT: (u8, u8) = (13, 14);
    pub(super) const IMPLICIT: (u8, u8) = (15, 16);
    /// The operand of `-`, `+` and the roots: powers and postfix forms, not
    /// implicit products (`-2x` is `(-2)·x`).
    pub(super) const PREFIX: u8 = 17;
    /// Right-associative: `2^3^2` is `2^(3^2)`.
    pub(super) const POWER: (u8, u8) = (20, 19);
}

/// The binding powers of the operators that make a
/// [`BinaryExpr`](K::BinaryExpr).
fn binary(kind: K) -> Option<(u8, u8)> {
    Some(match kind {
        K::OrKw => bp::OR,
        K::AndKw => bp::AND,
        K::Plus | K::Minus => bp::SUM,
        K::Star | K::Slash | K::DotStar | K::DotSlash => bp::PRODUCT,
        _ => return None,
    })
}

fn is_comparison(kind: K) -> bool {
    matches!(
        kind,
        K::EqEq | K::Neq | K::Lt | K::Le | K::Gt | K::Ge | K::InKw
    )
}

/// Whether a token of `kind` starts a factor of an implicit product: not a
/// sign, which makes a sum (`a -b` is `a - b`), and not a block.
fn starts_factor(kind: K) -> bool {
    matches!(
        kind,
        K::Ident
            | K::Number
            | K::String
            | K::LParen
            | K::LBracket
            | K::LBrace
            | K::Sqrt
            | K::Cbrt
            | K::FourthRoot
            | K::Error
    )
}

/// Tokens that close an enclosing construct, where recovery stops.
fn is_boundary(kind: K) -> bool {
    matches!(
        kind,
        K::RParen
            | K::RBracket
            | K::RBrace
            | K::Semicolon
            | K::EndKw
            | K::ThenKw
            | K::ElifKw
            | K::ElseKw
            | K::DoKw
            | K::Eof
    )
}

impl Parser<'_> {
    pub(super) fn root(&mut self) {
        self.builder.start_node(K::Root.into());
        self.with_newlines(true, |p| p.statements(&[]));
    }

    /// Statements up to one of `terminators` or the end, separated by `;` or
    /// line breaks. Inside brackets, a closing bracket ends them too: it
    /// belongs to the enclosing construct, which reports what is missing.
    fn statements(&mut self, terminators: &[K]) {
        loop {
            while self.at(K::Semicolon) {
                self.bump();
            }
            let kind = self.peek();
            if kind == K::Eof || terminators.contains(&kind) || self.at_enclosing_closer() {
                break;
            }
            if starts_statement(kind) {
                self.statement();
                self.statement_end(terminators);
            } else {
                self.error(E::Unexpected(kind));
                self.single(K::ErrorNode);
            }
        }
    }

    /// Checks that a statement ends here, or skips to the end of its line.
    fn statement_end(&mut self, terminators: &[K]) {
        let ended = |p: &Self| {
            let kind = p.peek();
            kind == K::Semicolon
                || kind == K::Eof
                || p.newline_before()
                || terminators.contains(&kind)
                || p.at_enclosing_closer()
        };
        if ended(self) {
            return;
        }
        self.error(E::ExpectedSeparator);
        self.start(K::ErrorNode);
        while !ended(self) {
            self.bump();
        }
        self.finish_node();
    }

    /// Whether the next token closes a bracket around the block being
    /// parsed.
    fn at_enclosing_closer(&self) -> bool {
        matches!(self.peek(), K::RParen | K::RBracket | K::RBrace) && self.newlines.contains(&false)
    }

    /// The statements of a block, up to one of `terminators`. An empty block
    /// is an empty node: the line breaks in it go to the enclosing node.
    fn body(&mut self, terminators: &[K]) {
        if !self.enter() {
            return;
        }
        let next = self.peek();
        if next == K::Eof || terminators.contains(&next) || self.at_enclosing_closer() {
            self.builder.start_node(K::Body.into());
        } else {
            self.start(K::Body);
        }
        self.with_newlines(true, |p| p.statements(terminators));
        self.finish_node();
        self.leave();
    }

    fn statement(&mut self) {
        match self.peek() {
            K::LocalKw => self.local_stmt(),
            K::ReturnKw => {
                self.start(K::ReturnStmt);
                self.bump();
                if starts_expression(self.peek_continuing()) {
                    self.expr();
                }
                self.finish_node();
            }
            K::BreakKw => self.single(K::BreakStmt),
            K::ContinueKw => self.single(K::ContinueStmt),
            K::ForKw => self.for_stmt(),
            K::WhileKw => self.while_stmt(),
            _ => {
                let checkpoint = self.checkpoint();
                self.expr();
                if self.at(K::ColonEq) {
                    self.start_at(checkpoint, K::Definition);
                    self.bump();
                    self.expr();
                    if self.at(K::WhenKw) {
                        self.start(K::WhenClause);
                        self.bump();
                        self.expr();
                        self.finish_node();
                    }
                } else {
                    self.start_at(checkpoint, K::ExprStmt);
                }
                self.finish_node();
            }
        }
    }

    fn local_stmt(&mut self) {
        self.start(K::LocalStmt);
        self.bump();
        loop {
            self.start(K::LocalItem);
            self.name();
            if self.at(K::ColonEq) {
                self.bump();
                self.expr();
            }
            self.finish_node();
            if !self.at(K::Comma) {
                break;
            }
            self.bump();
        }
        self.finish_node();
    }

    fn for_stmt(&mut self) {
        self.start(K::ForStmt);
        self.bump();
        self.name();
        self.expect(K::InKw);
        self.expr();
        self.expect(K::DoKw);
        self.body(&[K::EndKw]);
        self.expect(K::EndKw);
        self.finish_node();
    }

    fn while_stmt(&mut self) {
        self.start(K::WhileStmt);
        self.bump();
        self.expr();
        self.expect(K::DoKw);
        self.body(&[K::EndKw]);
        self.expect(K::EndKw);
        self.finish_node();
    }

    /// A name, as a [`Name`](K::Name) node, or an error.
    fn name(&mut self) {
        if self.at(K::Ident) {
            self.single(K::Name);
        } else {
            self.error(E::ExpectedName);
        }
    }

    /// An expression: a lambda, or operators over operands.
    pub(super) fn expr(&mut self) {
        if self.at_lambda() {
            self.lambda();
        } else {
            self.expr_bp(0);
        }
    }

    /// Whether a lambda starts here: a name or parenthesised names, then
    /// `->`. The scan stops at the first token that cannot be a parameter,
    /// so it costs at most the parameters.
    fn at_lambda(&self) -> bool {
        let mut tokens = self.lookahead();
        match tokens.next() {
            Some(K::Ident) => tokens.next() == Some(K::Arrow),
            Some(K::LParen) => {
                let mut token = tokens.next();
                if token == Some(K::Ident) {
                    token = tokens.next();
                    while token == Some(K::Comma) {
                        if tokens.next() != Some(K::Ident) {
                            return false;
                        }
                        token = tokens.next();
                    }
                }
                token == Some(K::RParen) && tokens.next() == Some(K::Arrow)
            }
            _ => false,
        }
    }

    fn lambda(&mut self) {
        if !self.enter() {
            return;
        }
        self.start(K::Lambda);
        if self.at(K::Ident) {
            self.single(K::Name);
        } else {
            self.param_list();
        }
        self.expect(K::Arrow);
        self.expr();
        self.finish_node();
        self.leave();
    }

    /// `(x, y)`, the parameters of a lambda; [`at_lambda`](Self::at_lambda)
    /// has checked their shape.
    fn param_list(&mut self) {
        self.start(K::ParamList);
        self.bump();
        self.with_newlines(false, |p| {
            while p.at(K::Ident) {
                p.single(K::Name);
                if !p.at(K::Comma) {
                    break;
                }
                p.bump();
            }
        });
        self.expect(K::RParen);
        self.finish_node();
    }

    /// An expression whose operators bind at least as tightly as `min`.
    ///
    /// The forms other than binary operators each wrap what came before
    /// them, so each after the first is a level of nesting: `x² ² ²`,
    /// `a = b = c` (CAS-168). The first is free, so that `x^2 = 1` takes
    /// no more levels than `x^2`.
    fn expr_bp(&mut self, min: u8) {
        if !self.enter() {
            return;
        }
        let checkpoint = self.checkpoint();
        self.operand();
        let mut wrapped = false;
        let mut levels = 0;
        loop {
            let kind = self.peek_continuing();
            if let Some((left, right)) = binary(kind) {
                if left < min {
                    break;
                }
                self.start_at(checkpoint, K::BinaryExpr);
                self.bump();
                self.expr_bp(right);
                self.finish_node();
                continue;
            }
            let left = match kind {
                K::Caret | K::DotCaret | K::Superscript => bp::POWER.0,
                K::Eq | K::EqEq | K::Neq | K::Lt | K::Le | K::Gt | K::Ge | K::InKw => bp::COMPARE.0,
                K::DotDot => bp::RANGE.0,
                kind if starts_factor(kind) => bp::IMPLICIT.0,
                _ => break,
            };
            if left < min {
                break;
            }
            if wrapped {
                if !self.enter() {
                    break;
                }
                levels += 1;
            }
            wrapped = true;
            match kind {
                K::Caret | K::DotCaret => {
                    self.start_at(checkpoint, K::PowerExpr);
                    self.bump();
                    self.expr_bp(bp::POWER.1);
                    self.finish_node();
                }
                K::Superscript => {
                    self.start_at(checkpoint, K::PowerExpr);
                    self.bump();
                    self.finish_node();
                }
                K::Eq => self.equation(checkpoint),
                K::DotDot => self.range(checkpoint),
                kind if is_comparison(kind) => self.comparison(checkpoint),
                _ => self.implicit(checkpoint),
            }
        }
        for _ in 0..levels {
            self.leave();
        }
        self.leave();
    }

    /// A prefix operator with its operand, or a postfix expression.
    fn operand(&mut self) {
        match self.peek() {
            K::Minus | K::Plus | K::Sqrt | K::Cbrt | K::FourthRoot => {
                self.start(K::PrefixExpr);
                self.bump();
                self.expr_bp(bp::PREFIX);
                self.finish_node();
            }
            K::NotKw => {
                self.start(K::PrefixExpr);
                self.bump();
                self.expr_bp(bp::NOT);
                self.finish_node();
            }
            _ => self.postfix_expr(),
        }
    }

    /// `a = b`; a second `=` or a comparison after it is an error.
    fn equation(&mut self, checkpoint: Checkpoint) {
        self.start_at(checkpoint, K::Equation);
        self.bump();
        self.expr_bp(bp::COMPARE.1);
        let next = self.peek_continuing();
        if next == K::Eq || is_comparison(next) {
            self.error(E::ChainedEquation);
        }
        self.finish_node();
    }

    /// `a < b <= c`, one node for the whole chain.
    fn comparison(&mut self, checkpoint: Checkpoint) {
        self.start_at(checkpoint, K::Comparison);
        while is_comparison(self.peek_continuing()) {
            self.bump();
            self.expr_bp(bp::COMPARE.1);
        }
        if self.peek_continuing() == K::Eq {
            self.error(E::ChainedEquation);
        }
        self.finish_node();
    }

    /// `a..b`; a second `..` is an error.
    fn range(&mut self, checkpoint: Checkpoint) {
        self.start_at(checkpoint, K::RangeExpr);
        self.bump();
        self.expr_bp(bp::RANGE.1);
        if self.peek_continuing() == K::DotDot {
            self.error(E::ChainedRange);
        }
        self.finish_node();
    }

    /// Factors side by side, one node for the whole product.
    fn implicit(&mut self, checkpoint: Checkpoint) {
        self.start_at(checkpoint, K::ImplicitProduct);
        while starts_factor(self.peek_continuing()) {
            if self.at(K::Number) {
                self.error(E::NumberInProduct);
            }
            self.expr_bp(bp::IMPLICIT.1);
        }
        self.finish_node();
    }

    /// A primary and its postfix forms: applications, indexes, `!`, `...` and
    /// `::`. Each form nests the ones before it, so each is a level of
    /// nesting (`x!!…!` is as deep as it is long).
    fn postfix_expr(&mut self) {
        let checkpoint = self.checkpoint();
        let mut callable = self.primary();
        let mut levels = 0;
        loop {
            let kind = self.peek_continuing();
            let postfix = matches!(kind, K::LBracket | K::Bang | K::Ellipsis | K::ColonColon)
                || (kind == K::LParen && callable);
            if !postfix || !self.enter() {
                break;
            }
            levels += 1;
            match kind {
                K::LParen if callable => {
                    self.start_at(checkpoint, K::Application);
                    self.args(K::RParen);
                    self.finish_node();
                }
                K::LBracket => {
                    self.start_at(checkpoint, K::IndexExpr);
                    self.args(K::RBracket);
                    self.finish_node();
                    callable = true;
                }
                kind @ (K::Bang | K::Ellipsis) => {
                    let node = if kind == K::Bang {
                        K::Factorial
                    } else {
                        K::Spread
                    };
                    self.start_at(checkpoint, node);
                    self.bump();
                    self.finish_node();
                    callable = false;
                }
                K::ColonColon => {
                    self.start_at(checkpoint, K::Typed);
                    self.bump();
                    self.name();
                    self.finish_node();
                    callable = false;
                }
                _ => unreachable!("not a postfix form"),
            }
        }
        for _ in 0..levels {
            self.leave();
        }
    }

    /// A primary expression. Returns whether `(…)` after it applies it:
    /// after a name or a parenthesised lambda. After a number or another
    /// parenthesised expression, `(…)` is a factor of an implicit product
    /// (`2(x + 1)`), so that `x^2(x + 1)` is `x^2·(x + 1)` (CAS-DOC-7 §3.1).
    fn primary(&mut self) -> bool {
        match self.peek() {
            K::Number | K::String => {
                self.single(K::Literal);
                false
            }
            K::Ident => {
                self.single(K::Name);
                true
            }
            K::LParen => self.paren(),
            K::LBracket => {
                self.delimited(K::ListExpr, K::RBracket);
                false
            }
            K::LBrace => {
                self.delimited(K::SetExpr, K::RBrace);
                false
            }
            K::IfKw => {
                self.if_expr();
                false
            }
            K::BeginKw => {
                self.block_expr();
                false
            }
            // The lexer has reported it.
            K::Error => {
                self.single(K::ErrorNode);
                false
            }
            _ => {
                self.error(E::ExpectedExpression);
                false
            }
        }
    }

    /// `(e)`. Returns whether `e` is a lambda.
    fn paren(&mut self) -> bool {
        self.start(K::Paren);
        self.bump();
        let (lambda, reported) = self.with_newlines(false, |p| {
            let lambda = p.at_lambda();
            if starts_expression(p.peek()) {
                p.expr();
            } else {
                p.error(E::ExpectedExpression);
            }
            if p.at(K::Comma) {
                p.error(E::Tuple);
                while p.at(K::Comma) {
                    p.bump();
                    if starts_expression(p.peek()) {
                        p.expr();
                    }
                }
            }
            let junk = !p.at(K::RParen) && !is_boundary(p.peek());
            if junk {
                p.error(E::Expected(K::RParen));
                p.skip_to(K::RParen);
            }
            (lambda, junk)
        });
        if self.at(K::RParen) {
            self.bump();
        } else if !reported {
            self.error(E::Expected(K::RParen));
        }
        self.finish_node();
        lambda
    }

    /// Puts the tokens before a comma, `closer` or a boundary into an error
    /// node.
    fn skip_to(&mut self, closer: K) {
        self.start(K::ErrorNode);
        while !self.at(K::Comma) && !self.at(closer) && !is_boundary(self.peek()) {
            self.bump();
        }
        self.finish_node();
    }

    /// A list or a set: `kind`, delimited by its opening token and `closer`.
    fn delimited(&mut self, kind: K, closer: K) {
        self.start(kind);
        self.bump();
        self.with_newlines(false, |p| p.items(closer));
        self.expect(closer);
        self.finish_node();
    }

    /// The arguments of an application or an index, up to `closer`.
    fn args(&mut self, closer: K) {
        self.start(K::ArgList);
        self.bump();
        self.with_newlines(false, |p| p.items(closer));
        self.expect(closer);
        self.finish_node();
    }

    /// Expressions separated by commas up to `closer`, a trailing comma
    /// allowed. Junk is skipped to the next comma or the closer, with one
    /// error.
    fn items(&mut self, closer: K) {
        while !self.at(closer) && !is_boundary(self.peek()) {
            if starts_expression(self.peek()) {
                self.expr();
                if !self.at(K::Comma) && !self.at(closer) && !is_boundary(self.peek()) {
                    self.error(E::Expected(closer));
                    self.skip_to(closer);
                }
            } else {
                self.error(E::ExpectedExpression);
                if !self.at(K::Comma) {
                    self.skip_to(closer);
                }
            }
            if self.at(K::Comma) {
                self.bump();
            }
        }
    }

    fn if_expr(&mut self) {
        const BRANCH_END: &[K] = &[K::ElifKw, K::ElseKw, K::EndKw];
        self.start(K::IfExpr);
        self.bump();
        self.expr();
        self.expect(K::ThenKw);
        self.body(BRANCH_END);
        while self.at(K::ElifKw) {
            self.bump();
            self.expr();
            self.expect(K::ThenKw);
            self.body(BRANCH_END);
        }
        if self.at(K::ElseKw) {
            self.bump();
            self.body(&[K::EndKw]);
        }
        self.expect(K::EndKw);
        self.finish_node();
    }

    fn block_expr(&mut self) {
        self.start(K::BlockExpr);
        self.bump();
        self.body(&[K::EndKw]);
        self.expect(K::EndKw);
        self.finish_node();
    }
}
