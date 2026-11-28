//! What is at an offset, for completion and signature help (CAS-DOC-60 §1):
//! the token there, what kind of spot it is (a name, a type after `::`, a
//! number, a string, a comment) and the call whose arguments it is in.
//!
//! An offset past the text is the text's end, one inside a character its
//! start. At a boundary between two tokens the left one is taken when typing
//! extends it (a name, a keyword, a number, a superscript, a line comment),
//! else the right one. A query walks the ancestors of at most three tokens
//! and never scans the tree.

use std::ops::Range;

use rowan::{NodeOrToken, TextRange, TextSize, TokenAtOffset};

use crate::lexer::LexErrorKind;
use crate::parser::{Parse, SyntaxErrorKind};
use crate::syntax::{SyntaxKind, SyntaxNode, SyntaxToken};

/// What is at an offset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cursor {
    /// The offset asked for, clamped to the text and rounded down to a
    /// character's start.
    pub offset: usize,
    /// What kind of spot it is.
    pub kind: CursorKind,
    /// The call whose arguments hold the offset, if any.
    pub call: Option<Call>,
}

/// The kind of spot at an offset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CursorKind {
    /// Inside a comment.
    Comment,
    /// Inside a string.
    String,
    /// On a number or a superscript.
    Number,
    /// On a name or a keyword (its range), or where one would start (an empty
    /// range at the offset).
    Name {
        /// The name's bytes.
        range: Range<usize>,
    },
    /// As [`Name`](CursorKind::Name), after `::`.
    Type {
        /// The type name's bytes.
        range: Range<usize>,
    },
    /// Anything else: inside an operator such as `:=`, or on a character that
    /// starts no token.
    Other,
}

/// A call around an offset: `f(…)`, never an index `v[…]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Call {
    /// What is called: the name, or the application for `f(a)(b)`.
    pub callee: Range<usize>,
    /// The arguments, from `(` to `)` or to where they stop.
    pub arguments: Range<usize>,
    /// The argument the offset is in, from 0: the commas before it.
    pub argument: usize,
}

/// The token at `offset`: the one it is strictly inside, else at a boundary
/// the left token when typing extends it (a name, a keyword, a number, a
/// superscript, a line comment), else the right one. `None` for the empty
/// text.
#[must_use]
pub fn token_at(parse: &Parse, offset: usize) -> Option<SyntaxToken> {
    pick(locate(&parse.syntax(), offset).1)
}

/// What is at `offset` (CAS-DOC-60 §1).
#[must_use]
pub fn cursor(parse: &Parse, offset: usize) -> Cursor {
    let (at, around) = locate(&parse.syntax(), offset);
    let offset = usize::from(at);
    let Some(token) = pick(around.clone()) else {
        return Cursor {
            offset,
            kind: CursorKind::Name {
                range: offset..offset,
            },
            call: None,
        };
    };
    let last = last_significant(&token, at);
    Cursor {
        offset,
        kind: kind_at(parse, &token, at, last.as_ref()),
        call: call_at(around, at, last.as_ref()),
    }
}

/// `offset` within the text, rounded down to a character's start, and the
/// tokens there: the one holding it strictly inside, or the ones that end
/// and start there (only one at the text's start or end). Each level of the
/// tree is binary-searched, as `token_at_offset` does not.
fn locate(root: &SyntaxNode, offset: usize) -> (TextSize, TokenAtOffset<SyntaxToken>) {
    let len = root.text_range().end();
    let at = TextSize::try_from(offset).map_or(len, |o| o.min(len));
    if at == len {
        let last = at
            .checked_sub(TextSize::from(1))
            .and_then(|byte| holding(root, byte));
        return (at, last.map_or(TokenAtOffset::None, TokenAtOffset::Single));
    }
    let Some(right) = holding(root, at) else {
        return (at, TokenAtOffset::None);
    };
    let start = right.text_range().start();
    let text = right.text();
    let back = (0..=usize::from(at - start))
        .rev()
        .find(|&i| text.is_char_boundary(i))
        .unwrap_or(0);
    let at = start + TextSize::of(&text[..back]);
    if at > start {
        return (at, TokenAtOffset::Single(right));
    }
    match right.prev_token() {
        Some(left) => (at, TokenAtOffset::Between(left, right)),
        None => (at, TokenAtOffset::Single(right)),
    }
}

/// The token holding the byte at `byte`, inside the text.
fn holding(root: &SyntaxNode, byte: TextSize) -> Option<SyntaxToken> {
    root.covering_element(TextRange::at(byte, TextSize::from(1)))
        .into_token()
}

/// The token at a character boundary by [`token_at`]'s rule, from the
/// tokens there. An unterminated string or comment runs to the text's end,
/// so it is never the left token of a boundary inside the text.
fn pick(around: TokenAtOffset<SyntaxToken>) -> Option<SyntaxToken> {
    match around {
        TokenAtOffset::None => None,
        TokenAtOffset::Single(token) => Some(token),
        TokenAtOffset::Between(left, right) => {
            Some(if extends(left.kind()) { left } else { right })
        }
    }
}

/// Whether typing right after a token of `kind` extends it.
fn extends(kind: SyntaxKind) -> bool {
    is_name(kind)
        || matches!(
            kind,
            SyntaxKind::Number | SyntaxKind::Superscript | SyntaxKind::LineComment
        )
}

fn is_name(kind: SyntaxKind) -> bool {
    kind == SyntaxKind::Ident || kind.is_keyword()
}

/// Neither trivia nor a line break.
fn significant(kind: SyntaxKind) -> bool {
    !kind.is_trivia() && kind != SyntaxKind::Newline
}

/// The last significant token that ends at or before `at`, from `token`,
/// which holds `at`, backwards.
fn last_significant(token: &SyntaxToken, at: TextSize) -> Option<SyntaxToken> {
    let mut next = Some(token.clone());
    while let Some(t) = next {
        if t.text_range().end() <= at && significant(t.kind()) {
            return Some(t);
        }
        next = t.prev_token();
    }
    None
}

/// The significant token before `token`.
fn previous_significant(token: &SyntaxToken) -> Option<SyntaxToken> {
    let mut next = token.prev_token();
    while let Some(t) = next {
        if significant(t.kind()) {
            return Some(t);
        }
        next = t.prev_token();
    }
    None
}

/// Whether the string or comment at `range` is unterminated: the lexer's
/// `error` starts where it does and ends at the text's end.
fn unterminated(parse: &Parse, range: TextRange, error: LexErrorKind) -> bool {
    let start = usize::from(range.start());
    let end = usize::from(range.end());
    let errors = parse.errors();
    let first = errors.partition_point(|e| e.range.start < start);
    errors[first..]
        .iter()
        .take_while(|e| e.range.start == start)
        .any(|e| e.kind == SyntaxErrorKind::Lex(error) && e.range.end == end)
}

fn span(range: TextRange) -> Range<usize> {
    usize::from(range.start())..usize::from(range.end())
}

/// The kind of spot at `at` on `token`; `last` is the last significant
/// token ending at or before `at`.
fn kind_at(
    parse: &Parse,
    token: &SyntaxToken,
    at: TextSize,
    last: Option<&SyntaxToken>,
) -> CursorKind {
    let range = token.text_range();
    let after_start = range.start() < at;
    let inside = after_start && at < range.end();
    match token.kind() {
        SyntaxKind::LineComment if after_start => CursorKind::Comment,
        SyntaxKind::BlockComment
            if inside
                || (after_start
                    && unterminated(parse, range, LexErrorKind::UnterminatedComment)) =>
        {
            CursorKind::Comment
        }
        SyntaxKind::String
            if inside
                || (after_start
                    && unterminated(parse, range, LexErrorKind::UnterminatedString)) =>
        {
            CursorKind::String
        }
        SyntaxKind::Number | SyntaxKind::Superscript => CursorKind::Number,
        kind if is_name(kind) => {
            let range = span(range);
            if previous_significant(token).is_some_and(|t| t.kind() == SyntaxKind::ColonColon) {
                CursorKind::Type { range }
            } else {
                CursorKind::Name { range }
            }
        }
        SyntaxKind::Error => CursorKind::Other,
        SyntaxKind::Whitespace
        | SyntaxKind::Newline
        | SyntaxKind::LineComment
        | SyntaxKind::BlockComment
        | SyntaxKind::String => empty(at, last),
        _ if inside => CursorKind::Other,
        _ => empty(at, last),
    }
}

/// Where a name would start: a type after `::`, else a name.
fn empty(at: TextSize, last: Option<&SyntaxToken>) -> CursorKind {
    let o = usize::from(at);
    if last.is_some_and(|t| t.kind() == SyntaxKind::ColonColon) {
        CursorKind::Type { range: o..o }
    } else {
        CursorKind::Name { range: o..o }
    }
}

/// The innermost call whose arguments hold `at`, among the ancestors of the
/// tokens on both sides of `at` and of `last`.
fn call_at(
    around: TokenAtOffset<SyntaxToken>,
    at: TextSize,
    last: Option<&SyntaxToken>,
) -> Option<Call> {
    let (left, right) = match around {
        TokenAtOffset::None => (None, None),
        TokenAtOffset::Single(token) => (Some(token), None),
        TokenAtOffset::Between(left, right) => (Some(left), Some(right)),
    };
    let mut best: Option<(TextSize, SyntaxNode)> = None;
    for token in [left, right, last.cloned()].into_iter().flatten() {
        // The first along a chain of ancestors opens last.
        let found = token
            .parent_ancestors()
            .find_map(|node| opener(&node, at, last).map(|open| (open, node)));
        if let Some((open, node)) = found {
            if best.as_ref().is_none_or(|(start, _)| open > *start) {
                best = Some((open, node));
            }
        }
    }
    let (_, arguments) = best?;
    let application = arguments.parent()?;
    let commas = |node: &SyntaxNode| {
        node.children_with_tokens()
            .take_while(|e| e.text_range().start() < at)
            .map(|e| match e {
                NodeOrToken::Token(t) => usize::from(t.kind() == SyntaxKind::Comma),
                NodeOrToken::Node(n) if n.kind() == SyntaxKind::ErrorNode => n
                    .children_with_tokens()
                    .take_while(|e| e.text_range().start() < at)
                    .filter(|e| e.kind() == SyntaxKind::Comma)
                    .count(),
                NodeOrToken::Node(_) => 0,
            })
            .sum::<usize>()
    };
    let argument = commas(&arguments);
    Some(Call {
        callee: span(application.first_child_or_token()?.text_range()),
        arguments: span(arguments.text_range()),
        argument,
    })
}

/// The start of `node`'s `(` when it is the arguments of an application
/// that hold `at`: after the `(`, and before the `)` or, when there is none,
/// at most at their end or with `last` inside them.
fn opener(node: &SyntaxNode, at: TextSize, last: Option<&SyntaxToken>) -> Option<TextSize> {
    if node.kind() != SyntaxKind::ArgList
        || node.parent().map(|p| p.kind()) != Some(SyntaxKind::Application)
    {
        return None;
    }
    let open = node
        .children_with_tokens()
        .find(|e| significant(e.kind()))?
        .into_token()
        .filter(|t| t.kind() == SyntaxKind::LParen)?;
    if open.text_range().end() > at {
        return None;
    }
    let range = node.text_range();
    // The `)` is added last, right after the arguments.
    let close = node
        .last_child_or_token()
        .and_then(NodeOrToken::into_token)
        .filter(|t| t.kind() == SyntaxKind::RParen);
    let holds = match close {
        Some(close) => close.text_range().start() >= at,
        None => at <= range.end() || last.is_some_and(|t| range.contains_range(t.text_range())),
    };
    holds.then(|| open.text_range().start())
}

#[cfg(test)]
mod tests;
