//! The cas-lang analysis queries (CAS-DOC-60 §7) on arbitrary text, parsed
//! once, at most 1024 offsets per input. T3: the line index's round trips,
//! monotonic positions, line count and clamping, in both encodings. T4: every
//! fix of a parse error or of a lowering diagnostic lies within the text on
//! character boundaries, and the text it gives parses without a panic. T5:
//! at every sampled offset, boundaries or not and past the end, the cursor's
//! offset is the clamped one, its ranges lie within the text on character
//! boundaries, a name's range holds the offset, a call's argument is at most
//! the commas of its arguments, and there is a token iff the text is not
//! empty, holding the offset. T6 adds its properties (symbols).
#![no_main]

use cas_lang::analysis::cursor::{CursorKind, cursor, token_at};
use cas_lang::analysis::lines::{Encoding, LineIndex, Position};
use cas_lang::{Names, lower, parse};
use libfuzzer_sys::fuzz_target;

/// The most offsets an input's queries run on.
const CAP: usize = 1024;

/// `0..=len` when it has at most `CAP` values, else `CAP` evenly spaced ones.
fn spaced(len: usize) -> Vec<usize> {
    if len < CAP {
        (0..=len).collect()
    } else {
        (0..CAP).map(|k| k * len / (CAP - 1)).collect()
    }
}

/// `offset` (at most the length) rounded down to a character boundary.
fn floor(text: &str, mut offset: usize) -> usize {
    while !text.is_char_boundary(offset) {
        offset -= 1;
    }
    offset
}

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let text: &str = &text;
    let parse = parse(text);
    let index = LineIndex::new(text);
    let bytes = text.as_bytes();
    let len = text.len();

    // Fixes: within the text, on boundaries, and the fixed text parses.
    let lowered = lower(&parse, Names::english());
    let fixes = parse.errors().iter().filter_map(|e| e.fix.as_ref());
    for fix in fixes.chain(lowered.diagnostics.iter().filter_map(|d| d.fix.as_ref())) {
        let r = &fix.range;
        let inside = r.start <= r.end && r.end <= len;
        assert!(inside && text.is_char_boundary(r.start) && text.is_char_boundary(r.end), "{fix:?}");
        let fixed = format!("{}{}{}", &text[..r.start], fix.replacement, &text[r.end..]);
        let _ = cas_lang::parse(&fixed);
    }
    assert_eq!(index.len(), len);

    // The sampled character boundaries, in order.
    let boundaries: Vec<usize> = {
        let all: Vec<usize> = (0..=len).filter(|&o| text.is_char_boundary(o)).collect();
        if all.len() <= CAP {
            all
        } else {
            let mut sampled: Vec<usize> = spaced(len).into_iter().map(|o| floor(text, o)).collect();
            sampled.dedup();
            sampled
        }
    };
    let breaks = bytes
        .iter()
        .enumerate()
        .filter(|&(i, &b)| b == b'\n' || (b == b'\r' && bytes.get(i + 1) != Some(&b'\n')))
        .count();
    assert_eq!(index.line_count() as usize, breaks + 1);
    let line_count = index.line_count();
    // Between the `\r` and `\n` of a `\r\n`.
    let in_crlf = |o: usize| o > 0 && bytes[o - 1] == b'\r' && bytes.get(o) == Some(&b'\n');

    for encoding in [Encoding::Utf8, Encoding::Utf16] {
        // Every byte offset, and past the end, uses the character's start.
        for o in spaced(len) {
            assert_eq!(index.position(o, encoding), index.position(floor(text, o), encoding));
        }
        for o in [len + 1, usize::MAX] {
            assert_eq!(index.position(o, encoding), index.position(len, encoding));
        }

        let mut last = Position { line: 0, character: 0 };
        for &o in &boundaries {
            let p = index.position(o, encoding);
            assert!(last <= p, "position not monotonic at {o}");
            last = p;
            let back = index.offset(p, encoding);
            assert_eq!(back, if in_crlf(o) { o - 1 } else { o });
            assert_eq!(index.position(back, encoding), p);

            let near = [
                Position { line: p.line, character: p.character.saturating_add(1) },
                Position { line: p.line, character: p.character.saturating_add(2) },
                Position { line: p.line, character: p.character.saturating_add(3) },
                Position { line: p.line.saturating_add(1), character: 0 },
                Position { line: line_count.saturating_add(1), character: 0 },
            ];
            for q in near {
                let b = index.offset(q, encoding);
                assert!(b <= len && text.is_char_boundary(b) && !in_crlf(b), "{q:?} -> {b}");
                assert!(index.position(b, encoding) <= q, "{q:?} -> {b}");
            }
        }
    }

    // The cursor and the token at sampled offsets, boundaries or not.
    let on_text = |r: &std::ops::Range<usize>| {
        let inside = r.start <= r.end && r.end <= len;
        inside && text.is_char_boundary(r.start) && text.is_char_boundary(r.end)
    };
    let raw = spaced(len).into_iter().chain([len + 1, usize::MAX]);
    for o in boundaries.iter().copied().chain(raw) {
        let c = cursor(&parse, o);
        let at = floor(text, o.min(len));
        assert_eq!(c.offset, at, "{o}");
        if let CursorKind::Name { range } | CursorKind::Type { range } = &c.kind {
            assert!(on_text(range) && range.start <= at && at <= range.end, "{o}: {c:?}");
        }
        if let Some(call) = &c.call {
            assert!(on_text(&call.callee) && on_text(&call.arguments), "{o}: {c:?}");
            let commas = bytes[call.arguments.clone()].iter().filter(|&&b| b == b',').count();
            assert!(call.argument <= commas, "{o}: {c:?}");
        }
        let token = token_at(&parse, o);
        assert_eq!(token.is_some(), len > 0, "{o}");
        if let Some(token) = token {
            let r = token.text_range();
            assert!(usize::from(r.start()) <= at && at <= usize::from(r.end()), "{o}: {token:?}");
        }
    }
});
