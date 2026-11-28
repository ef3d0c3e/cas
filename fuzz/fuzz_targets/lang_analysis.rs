//! The cas-lang analysis queries (CAS-DOC-60 §7) on arbitrary text, parsed
//! once, at most 1024 offsets per input. T3: the line index's round trips,
//! monotonic positions, line count and clamping, in both encodings. T4-T6 add
//! their properties (fixes, cursor, symbols).
#![no_main]

use cas_lang::analysis::lines::{Encoding, LineIndex, Position};
use cas_lang::parse;
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
    let _parse = parse(text);
    let index = LineIndex::new(text);
    let bytes = text.as_bytes();
    let len = text.len();
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
});
