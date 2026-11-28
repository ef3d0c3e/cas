//! Fixes for the lexer's and the parser's errors (CAS-DOC-60 §5): one
//! replacement each, given only where the intended text is plain.
//!
//! Runs on every parse, so it stays linear in the text: the errors that fix
//! together (the bars of a line, a run of `&`, the primes of `f'(x)`) are
//! grouped in one pass over the tokens, and line starts are found once.

use std::collections::HashMap;
use std::ops::Range;

use rowan::{GreenNode, TextSize};

use super::{Fix, ParseErrorKind as E, SyntaxError, SyntaxErrorKind, starts_statement};
use crate::lexer::{LexErrorKind as L, Lexed};
use crate::syntax::{SyntaxKind as K, SyntaxNode};

/// Fills the fix of each error in `errors`, the errors of `lexed` parsed
/// into `green`.
pub(super) fn fill(lexed: &Lexed<'_>, green: &GreenNode, errors: &mut [SyntaxError]) {
    if errors.is_empty() {
        return;
    }
    let mut fixes = Fixes {
        lexed,
        text: lexed.source(),
        root: SyntaxNode::new_root(green.clone()),
        grouped: None,
        lines: None,
        line_break: None,
        indent: None,
    };
    for e in errors.iter_mut() {
        e.fix = match e.kind {
            SyntaxErrorKind::Lex(kind) => fixes.lex(e.range.clone(), kind),
            SyntaxErrorKind::Parse(kind) => fixes.parse(e.range.clone(), kind),
        };
    }
}

/// What the fixes of one text share, computed when first needed.
struct Fixes<'a> {
    lexed: &'a Lexed<'a>,
    text: &'a str,
    root: SyntaxNode,
    /// The fixes of bars, ampersands and apostrophes, by where the error
    /// starts.
    grouped: Option<HashMap<usize, Fix>>,
    /// Where each line starts.
    lines: Option<Vec<usize>>,
    /// The text's first line break, found with the line starts.
    line_break: Option<&'a str>,
    /// The last line start asked for, and its indentation.
    indent: Option<(usize, &'a str)>,
}

impl<'a> Fixes<'a> {
    fn lex(&mut self, range: Range<usize>, kind: L) -> Option<Fix> {
        let text = self.text;
        match kind {
            L::Bar | L::Ampersand | L::Apostrophe => self.grouped().get(&range.start).cloned(),
            L::Hash => Some(fix(range, "//")),
            L::LoneColon => {
                let rest = &text[range.end..];
                let after = rest.trim_start_matches(is_space);
                let spaces = range.end + rest.len() - after.len();
                if after.starts_with('=') {
                    Some(fix(range.start..spaces + 1, ":="))
                } else if after.is_empty() || after.starts_with(['\n', '\r']) {
                    Some(fix(range, ";"))
                } else if after.starts_with(';') {
                    Some(fix(range, ""))
                } else {
                    None
                }
            }
            L::LeadingPoint => Some(fix(range.start..range.start, "0")),
            L::MisplacedUnderscore => {
                let mut digits = String::new();
                for c in text[range.clone()].chars() {
                    if c != '_' || !digits.ends_with('_') {
                        digits.push(c);
                    }
                }
                let digits = digits.trim_end_matches('_').to_owned();
                Some(fix(range, &digits))
            }
            L::UnterminatedString => {
                let at = text[range.start..]
                    .find(['\n', '\r'])
                    .map_or(text.len(), |i| range.start + i);
                Some(fix(at..at, "\""))
            }
            L::UnterminatedComment => {
                // The nesting left open at the end, as the lexer counts it.
                let mut depth = 0usize;
                let mut chars = text[range.clone()].chars().peekable();
                while let Some(c) = chars.next() {
                    match (c, chars.peek()) {
                        ('/', Some('*')) => {
                            chars.next();
                            depth += 1;
                        }
                        ('*', Some('/')) => {
                            chars.next();
                            depth = depth.saturating_sub(1);
                        }
                        _ => {}
                    }
                }
                let end = text.len();
                Some(fix(end..end, &"*/".repeat(depth.max(1))))
            }
            L::InvalidEscape => {
                let escape = &text[range.clone()];
                (escape.len() > 1 && !escape[1..].starts_with('u'))
                    .then(|| fix(range.start..range.start, "\\"))
            }
            _ => None,
        }
    }

    fn parse(&mut self, range: Range<usize>, kind: E) -> Option<Fix> {
        let at = range.start..range.start;
        match kind {
            E::Expected(closer) => {
                let word = match closer {
                    K::RParen => ")",
                    K::RBracket => "]",
                    K::RBrace => "}",
                    K::Arrow => "->",
                    K::ThenKw => return Some(padded(self.text, at, "then")),
                    K::DoKw => return Some(padded(self.text, at, "do")),
                    K::InKw => return Some(padded(self.text, at, "in")),
                    K::EndKw => return Some(padded(self.text, at, "end")),
                    _ => return None,
                };
                Some(fix(at, word))
            }
            E::ExpectedSeparator => self.separator(range.start),
            E::Unexpected(K::RParen | K::RBracket | K::RBrace | K::EndKw) => Some(fix(range, "")),
            E::Tuple => {
                let comma = self.token_at(range.start)?;
                let paren = comma.parent()?;
                let closed = paren
                    .last_child_or_token()
                    .is_some_and(|last| last.kind() == K::RParen);
                if paren.kind() != K::Paren || !closed {
                    return None;
                }
                let span = paren.text_range();
                let inner = &self.text[usize::from(span.start()) + 1..usize::from(span.end()) - 1];
                Some(fix(
                    usize::from(span.start())..usize::from(span.end()),
                    &format!("[{inner}]"),
                ))
            }
            E::NumberInProduct => Some(fix(at, "*")),
            _ => None,
        }
    }

    /// A line break, with the indentation of the statement's line, in place
    /// of the whitespace before the token at `start`, when it can start a
    /// statement.
    fn separator(&mut self, start: usize) -> Option<Fix> {
        let i = token_index(self.lexed, start)?;
        if !starts_statement(self.lexed.kind(i)) {
            return None;
        }
        let from = if i > 0 && self.lexed.kind(i - 1) == K::Whitespace {
            self.lexed.range(i - 1).start
        } else {
            start
        };
        // The statement is the node before the error node that holds the
        // token.
        let statement = self
            .token_at(start)
            .and_then(|t| t.parent())
            .filter(|n| n.kind() == K::ErrorNode)
            .and_then(|n| n.prev_sibling())
            .map_or(start, |n| usize::from(n.text_range().start()));
        let indent = self.indent(statement);
        let line_break = self.line_break.unwrap_or("\n");
        Some(fix(from..start, &format!("{line_break}{indent}")))
    }

    /// The indentation of the line holding `offset`.
    fn indent(&mut self, offset: usize) -> &'a str {
        let lexed = self.lexed;
        let line_break = &mut self.line_break;
        let lines = self.lines.get_or_insert_with(|| {
            let breaks: Vec<usize> = (0..lexed.len())
                .filter(|&i| lexed.kind(i) == K::Newline)
                .collect();
            *line_break = breaks.first().map(|&i| lexed.text(i));
            let mut starts = vec![0];
            starts.extend(breaks.iter().map(|&i| lexed.range(i).end));
            starts
        });
        let line = lines[lines.partition_point(|&s| s <= offset) - 1];
        if let Some((start, indent)) = self.indent
            && start == line
        {
            return indent;
        }
        let rest = &self.text[line..];
        let indent = &rest[..rest.len() - rest.trim_start_matches(is_space).len()];
        self.indent = Some((line, indent));
        indent
    }

    /// The token of the tree that starts at or covers `offset`.
    fn token_at(&self, offset: usize) -> Option<crate::syntax::SyntaxToken> {
        let offset = TextSize::from(u32::try_from(offset).ok()?);
        self.root.token_at_offset(offset).right_biased()
    }

    fn grouped(&mut self) -> &HashMap<usize, Fix> {
        let (lexed, text) = (self.lexed, self.text);
        self.grouped.get_or_insert_with(|| group(lexed, text))
    }
}

/// The fixes of the errors that go together, in one pass over the tokens:
/// the bars of each line, runs of `&`, and `f'(x)`.
fn group(lexed: &Lexed<'_>, text: &str) -> HashMap<usize, Fix> {
    let is = |i: usize, s: &str| lexed.kind(i) == K::Error && lexed.text(i) == s;
    let mut fixes = HashMap::new();
    let mut bars = Vec::new();
    let mut i = 0;
    while i < lexed.len() {
        let kind = lexed.kind(i);
        if kind == K::Newline {
            bar_pairs(lexed, text, &bars, &mut fixes);
            bars.clear();
        } else if is(i, "|") {
            bars.push(i);
        } else if is(i, "&") {
            let mut end = i + 1;
            while end < lexed.len() && is(end, "&") {
                end += 1;
            }
            if end - i <= 2 {
                let span = lexed.range(i).start..lexed.range(end - 1).end;
                let and = padded(text, span, "and");
                for k in i..end {
                    fixes.insert(lexed.range(k).start, and.clone());
                }
            }
            i = end;
            continue;
        } else if kind == K::Ident {
            // `f''(x)`: the name, primes, then `(`, one name and `)`, with
            // nothing between them.
            let mut end = i + 1;
            while end < lexed.len() && is(end, "'") {
                end += 1;
            }
            let primes = end - i - 1;
            let shape = [K::LParen, K::Ident, K::RParen];
            if primes > 0
                && end + 3 <= lexed.len()
                && (0..3).all(|k| lexed.kind(end + k) == shape[k])
            {
                let (f, x) = (lexed.text(i), lexed.text(end + 1));
                let mut replacement = format!("diff({f}({x})");
                for _ in 0..primes {
                    replacement.push_str(", ");
                    replacement.push_str(x);
                }
                replacement.push(')');
                let diff = fix(lexed.range(i).start..lexed.range(end + 2).end, &replacement);
                for k in i + 1..end {
                    fixes.insert(lexed.range(k).start, diff.clone());
                }
            }
            i = end;
            continue;
        }
        i += 1;
    }
    bar_pairs(lexed, text, &bars, &mut fixes);
    fixes
}

/// The fixes of the bars of one line: paired left to right, each pair
/// `|e|` becomes `abs(e)` when the count is even and every pair holds text
/// with balanced brackets; otherwise none.
fn bar_pairs(lexed: &Lexed<'_>, text: &str, bars: &[usize], fixes: &mut HashMap<usize, Fix>) {
    if bars.is_empty() || bars.len() % 2 != 0 {
        return;
    }
    let mut pairs = Vec::new();
    for pair in bars.chunks(2) {
        let (open, close) = (lexed.range(pair[0]), lexed.range(pair[1]));
        let inner = &text[open.end..close.start];
        if inner.trim().is_empty() || !balanced(lexed, pair[0] + 1..pair[1]) {
            return;
        }
        pairs.push(fix(open.start..close.end, &format!("abs({inner})")));
    }
    for (pair, abs) in bars.chunks(2).zip(pairs) {
        fixes.insert(lexed.range(pair[0]).start, abs.clone());
        fixes.insert(lexed.range(pair[1]).start, abs);
    }
}

/// Whether the tokens `tokens` have their `()`, `[]` and `{}` balanced.
fn balanced(lexed: &Lexed<'_>, tokens: Range<usize>) -> bool {
    let mut open = Vec::new();
    for i in tokens {
        match lexed.kind(i) {
            K::LParen => open.push(K::RParen),
            K::LBracket => open.push(K::RBracket),
            K::LBrace => open.push(K::RBrace),
            kind @ (K::RParen | K::RBracket | K::RBrace) if open.pop() != Some(kind) => {
                return false;
            }
            _ => {}
        }
    }
    open.is_empty()
}

/// `word` in place of `at`, with a space before it after a character that is
/// not whitespace and one after it before such a character.
fn padded(text: &str, at: Range<usize>, word: &str) -> Fix {
    let mut replacement = String::new();
    if text[..at.start]
        .chars()
        .next_back()
        .is_some_and(|c| !c.is_whitespace())
    {
        replacement.push(' ');
    }
    replacement.push_str(word);
    if text[at.end..]
        .chars()
        .next()
        .is_some_and(|c| !c.is_whitespace())
    {
        replacement.push(' ');
    }
    Fix {
        range: at,
        replacement,
    }
}

/// The index of the token that covers `offset`.
fn token_index(lexed: &Lexed<'_>, offset: usize) -> Option<usize> {
    let (mut low, mut high) = (0, lexed.len());
    while low < high {
        let mid = low + (high - low) / 2;
        if lexed.range(mid).end <= offset {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    (low < lexed.len()).then_some(low)
}

/// Whitespace within a line.
fn is_space(c: char) -> bool {
    c.is_whitespace() && c != '\n' && c != '\r'
}

fn fix(range: Range<usize>, replacement: &str) -> Fix {
    Fix {
        range,
        replacement: replacement.to_owned(),
    }
}
