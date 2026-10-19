//! The surface-language lexer on arbitrary text: lossless, total, and every
//! literal token's value readable without panicking.
#![no_main]

use cas_lang::{SyntaxKind, lex, name, number_value, string_value, superscript_value};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let lexed = lex(&text);
    let mut end = 0;
    for i in 0..lexed.len() {
        let range = lexed.range(i);
        assert_eq!(range.start, end);
        assert!(range.start < range.end);
        end = range.end;
        let token = lexed.text(i);
        let alone = lex(token);
        assert!(alone.len() == 1 && alone.kind(0) == lexed.kind(i), "{token:?}");
        match lexed.kind(i) {
            SyntaxKind::Number => {
                let _ = number_value(token);
            }
            SyntaxKind::String => {
                let _ = string_value(token);
            }
            SyntaxKind::Superscript => {
                let _ = superscript_value(token);
            }
            SyntaxKind::Ident => assert!(!name(token).is_empty()),
            SyntaxKind::Error => {
                assert!(lexed.errors().iter().any(|e| e.range == range));
            }
            _ => {}
        }
    }
    assert_eq!(end, text.len());
    let mut last = 0;
    for e in lexed.errors() {
        assert!(last <= e.range.start && e.range.start < e.range.end && e.range.end <= end);
        assert!(text.is_char_boundary(e.range.start) && text.is_char_boundary(e.range.end));
        last = e.range.start;
    }
});
