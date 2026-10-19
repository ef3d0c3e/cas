//! The surface-language parser on arbitrary text: lossless, total, with its
//! errors in order and within the text.
#![no_main]

use cas_lang::{SyntaxKind, parse};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let parse = parse(&text);
    let root = parse.syntax();
    assert_eq!(root.kind(), SyntaxKind::Root);
    assert_eq!(root.text().to_string(), text);
    let mut last = 0;
    for e in parse.errors() {
        assert!(last <= e.range.start && e.range.start <= e.range.end && e.range.end <= text.len());
        assert!(text.is_char_boundary(e.range.start) && text.is_char_boundary(e.range.end));
        last = e.range.start;
    }
    // Children cover their parent exactly, in order.
    for node in root.descendants() {
        let mut end = node.text_range().start();
        for child in node.children_with_tokens() {
            assert_eq!(child.text_range().start(), end);
            end = child.text_range().end();
        }
        assert_eq!(end, node.text_range().end());
    }
});
