//! The cas-lang analysis queries (CAS-DOC-60 §7) on arbitrary text, parsed
//! once, at most 1024 offsets per input. T3: the line index's round trips,
//! monotonic positions, line count and clamping, in both encodings. T4: every
//! fix of a parse error or of a lowering diagnostic lies within the text on
//! character boundaries, and the text it gives parses without a panic. T5:
//! at the sampled boundaries, the inner bytes of the first 64 multi-byte
//! characters, `len + 1` and `usize::MAX`, the cursor's offset is the
//! clamped one, its ranges lie within the text on character boundaries, a
//! name's range holds the offset, a call's argument is at most the commas of
//! its arguments, and there is a token iff the text is not empty, holding the
//! offset. T6: the symbols are built once; their occurrences are sorted,
//! disjoint `Name` nodes on boundaries, each in its symbol's definitions or
//! uses (and nothing else is listed), its text keyed as `name` gives it: a
//! user symbol's name, else the type or the reserved target (English names)
//! of the symbol's kind; every definition, use and scope lies within the text
//! on boundaries, the document is scope 0 and every other scope lies in its
//! parent, after it, disjoint from its siblings; `at` at the sampled offsets
//! is the occurrence whose range is the token's.
#![no_main]

use std::collections::HashMap;

use cas_lang::analysis::cursor::{CursorKind, cursor, token_at};
use cas_lang::analysis::lines::{Encoding, LineIndex, Position};
use cas_lang::analysis::symbols::{ScopeKind, SymbolKind, symbols};
use cas_lang::{Names, SyntaxKind, Target, lower, parse};
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

    // Symbols, built once: occurrences are `Name` nodes, sorted, disjoint,
    // each listed in its symbol and keyed by its text.
    let names = Names::english();
    let table = symbols(&parse, names);
    let name_nodes: Vec<(usize, usize)> = parse
        .syntax()
        .descendants()
        .filter(|n| n.kind() == SyntaxKind::Name)
        .map(|n| (usize::from(n.text_range().start()), usize::from(n.text_range().end())))
        .collect();
    let occurrences = table.occurrences();
    let by_range: HashMap<(usize, usize), usize> =
        occurrences.iter().enumerate().map(|(k, o)| ((o.range.start, o.range.end), k)).collect();
    for pair in occurrences.windows(2) {
        assert!(pair[0].range.end <= pair[1].range.start, "{pair:?}");
    }
    for o in occurrences {
        let r = (o.range.start, o.range.end);
        assert!(on_text(&o.range) && name_nodes.binary_search(&r).is_ok(), "{o:?}");
        let info = &table.symbols()[o.symbol.0];
        let listed = if o.definition {
            info.definitions.iter().any(|d| d.name == o.range)
        } else {
            info.uses.contains(&o.range)
        };
        assert!(listed, "{o:?} {info:?}");
        let key = cas_lang::name(&text[o.range.clone()]);
        let keyed = match (info.kind, names.target(&key)) {
            (SymbolKind::Type(ty), _) => names.type_name(&key) == Some(ty),
            (kind, None) => kind.is_user() && info.name == key,
            (SymbolKind::Builtin(s), Some(Target::Builtin(t))) => s == *t,
            (SymbolKind::Planned, Some(Target::Planned(_))) => true,
            (SymbolKind::Imaginary, Some(Target::Imaginary)) => true,
            (SymbolKind::Syntax(f), Some(Target::Syntax(g))) => f == *g,
            _ => false,
        };
        assert!(keyed, "{o:?} {key} {info:?}");
    }
    // Each symbol lists only its occurrences, with their flags.
    let mut count = 0;
    for (id, info) in table.symbols().iter().enumerate() {
        let definitions = info.definitions.iter().map(|d| (&d.name, true));
        for (r, definition) in definitions.chain(info.uses.iter().map(|u| (u, false))) {
            let o = by_range.get(&(r.start, r.end)).map(|&k| &occurrences[k]);
            let same = o.is_some_and(|o| o.symbol.0 == id && o.definition == definition);
            assert!(same, "{r:?} {info:?}");
            count += 1;
        }
        assert!(info.definitions.iter().all(|d| on_text(&d.node)), "{info:?}");
    }
    assert_eq!(count, occurrences.len());
    // Scopes nest: the document first, each other one in its parent, after
    // it, and disjoint from its siblings.
    let scopes = table.scopes();
    let document = &scopes[0];
    assert!(document.kind == ScopeKind::Document && document.parent.is_none());
    assert_eq!(document.range, 0..len);
    for (k, s) in scopes.iter().enumerate().skip(1) {
        let p = s.parent.map_or(usize::MAX, |p| p.0);
        assert!(p < k && on_text(&s.range), "{s:?}");
        let outer = &scopes[p].range;
        assert!(outer.start <= s.range.start && s.range.end <= outer.end, "{s:?}");
    }
    let mut siblings: Vec<usize> = (1..scopes.len()).collect();
    siblings.sort_by_key(|&k| (scopes[k].parent.map(|p| p.0), scopes[k].range.start));
    for pair in siblings.windows(2) {
        let (a, b) = (&scopes[pair[0]], &scopes[pair[1]]);
        assert!(a.parent != b.parent || a.range.end <= b.range.start, "{a:?} {b:?}");
    }

    // Raw non-boundary offsets only inside the first 64 multi-byte characters.
    let inner = text.char_indices().filter(|(_, ch)| ch.len_utf8() > 1).take(64);
    let raw = inner.flat_map(|(i, ch)| i + 1..i + ch.len_utf8()).chain([len + 1, usize::MAX]);
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
        let covering = token.as_ref().and_then(|t| {
            let r = t.text_range();
            by_range.get(&(usize::from(r.start()), usize::from(r.end())))
        });
        assert_eq!(table.at(o), covering.map(|&k| &occurrences[k]), "{o}");
        if let Some(token) = token {
            let r = token.text_range();
            assert!(usize::from(r.start()) <= at && at <= usize::from(r.end()), "{o}: {token:?}");
        }
    }
});
