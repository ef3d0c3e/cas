use proptest::prelude::*;

use super::{Encoding, LineIndex, Position};

const U8: Encoding = Encoding::Utf8;
const U16: Encoding = Encoding::Utf16;

fn pos(line: u32, character: u32) -> Position {
    Position { line, character }
}

/// Checks `position` at `offset` in both encodings, as `(line, utf8, utf16)`.
fn check_position(index: &LineIndex, offset: usize, line: u32, utf8: u32, utf16: u32) {
    assert_eq!(
        index.position(offset, U8),
        pos(line, utf8),
        "utf-8 at {offset}"
    );
    assert_eq!(
        index.position(offset, U16),
        pos(line, utf16),
        "utf-16 at {offset}"
    );
}

#[test]
fn ascii() {
    let index = LineIndex::new("ab\ncd");
    assert_eq!(index.len(), 5);
    assert_eq!(index.line_count(), 2);
    check_position(&index, 0, 0, 0, 0);
    check_position(&index, 2, 0, 2, 2);
    check_position(&index, 3, 1, 0, 0);
    check_position(&index, 5, 1, 2, 2);
    for e in [U8, U16] {
        assert_eq!(index.offset(pos(0, 1), e), 1);
        assert_eq!(index.offset(pos(1, 1), e), 4);
    }
}

#[test]
fn wide_characters() {
    // é is 2 bytes, € 3, 𝄞 4 (a surrogate pair in UTF-16).
    let index = LineIndex::new("aé€𝄞b");
    check_position(&index, 1, 0, 1, 1);
    check_position(&index, 3, 0, 3, 2);
    check_position(&index, 6, 0, 6, 3);
    check_position(&index, 10, 0, 10, 5);
    check_position(&index, 11, 0, 11, 6);
    // Inside a character: its start.
    check_position(&index, 2, 0, 1, 1);
    check_position(&index, 4, 0, 3, 2);
    check_position(&index, 5, 0, 3, 2);
    check_position(&index, 7, 0, 6, 3);
    check_position(&index, 9, 0, 6, 3);
    assert_eq!(index.offset(pos(0, 2), U16), 3);
    assert_eq!(index.offset(pos(0, 3), U16), 6);
    // Between the halves of the pair: rounded down.
    assert_eq!(index.offset(pos(0, 4), U16), 6);
    assert_eq!(index.offset(pos(0, 5), U16), 10);
    assert_eq!(index.offset(pos(0, 6), U16), 11);
    // Inside a UTF-8 sequence: rounded down.
    assert_eq!(index.offset(pos(0, 2), U8), 1);
    assert_eq!(index.offset(pos(0, 5), U8), 3);
    assert_eq!(index.offset(pos(0, 8), U8), 6);
    assert_eq!(index.offset(pos(0, 10), U8), 10);
}

#[test]
fn crlf() {
    let index = LineIndex::new("a\r\nb");
    assert_eq!(index.line_count(), 2);
    check_position(&index, 1, 0, 1, 1);
    // On the `\n` of a `\r\n`: the end of the line's content.
    check_position(&index, 2, 0, 1, 1);
    check_position(&index, 3, 1, 0, 0);
    for e in [U8, U16] {
        // Past the content, never inside the break.
        assert_eq!(index.offset(pos(0, 2), e), 1);
        assert_eq!(index.offset(pos(0, 9), e), 1);
        assert_eq!(index.offset(pos(1, 0), e), 3);
    }
}

#[test]
fn lone_cr() {
    let index = LineIndex::new("a\rb\r");
    assert_eq!(index.line_count(), 3);
    check_position(&index, 1, 0, 1, 1);
    check_position(&index, 2, 1, 0, 0);
    check_position(&index, 4, 2, 0, 0);
    assert_eq!(index.offset(pos(0, 5), U16), 1);
    assert_eq!(index.offset(pos(1, 1), U8), 3);
}

#[test]
fn empty() {
    let index = LineIndex::new("");
    assert_eq!(index.len(), 0);
    assert_eq!(index.line_count(), 1);
    for o in [0, 1, usize::MAX] {
        check_position(&index, o, 0, 0, 0);
    }
    for e in [U8, U16] {
        assert_eq!(index.offset(pos(0, 0), e), 0);
        assert_eq!(index.offset(pos(0, 7), e), 0);
        assert_eq!(index.offset(pos(u32::MAX, u32::MAX), e), 0);
    }
}

#[test]
fn empty_last_line() {
    let index = LineIndex::new("a\n");
    assert_eq!(index.line_count(), 2);
    check_position(&index, 2, 1, 0, 0);
    for e in [U8, U16] {
        assert_eq!(index.offset(pos(1, 0), e), 2);
        assert_eq!(index.offset(pos(1, 3), e), 2);
    }
}

#[test]
fn past_end() {
    let index = LineIndex::new("ab\nc€");
    check_position(&index, 6, 1, 4, 2);
    check_position(&index, 100, 1, 4, 2);
    check_position(&index, usize::MAX, 1, 4, 2);
    for e in [U8, U16] {
        // A line at or past the count: the end.
        assert_eq!(index.offset(pos(2, 0), e), 7);
        assert_eq!(index.offset(pos(u32::MAX, 0), e), 7);
        // A column past the content: the content's end.
        assert_eq!(index.offset(pos(0, 3), e), 2);
        assert_eq!(index.offset(pos(1, u32::MAX), e), 7);
    }
}

#[test]
fn send_sync() {
    fn check<T: Send + Sync>() {}
    check::<LineIndex>();
}

/// Text with every kind of line break and characters of 1 to 4 bytes.
fn text() -> impl Strategy<Value = String> {
    let piece = prop_oneof![
        Just("\n".to_owned()),
        Just("\r".to_owned()),
        Just("\r\n".to_owned()),
        Just("é".to_owned()),
        Just("€".to_owned()),
        Just("𝄞".to_owned()),
        "[a-z ]{1,3}",
        any::<char>().prop_map(String::from),
    ];
    prop::collection::vec(piece, 0..40).prop_map(|p| p.concat())
}

proptest! {
    #[test]
    fn round_trips(text in text()) {
        let index = LineIndex::new(&text);
        let bytes = text.as_bytes();
        let breaks = bytes
            .iter()
            .enumerate()
            .filter(|&(i, &b)| b == b'\n' || (b == b'\r' && bytes.get(i + 1) != Some(&b'\n')))
            .count();
        prop_assert_eq!(index.line_count() as usize, breaks + 1);
        for e in [U8, U16] {
            let mut last = pos(0, 0);
            for o in (0..=text.len()).filter(|&o| text.is_char_boundary(o)) {
                let p = index.position(o, e);
                prop_assert!(last <= p, "not monotonic at {}", o);
                last = p;
                let back = index.offset(p, e);
                let crlf_lf = o > 0 && bytes[o - 1] == b'\r' && bytes.get(o) == Some(&b'\n');
                prop_assert_eq!(back, if crlf_lf { o - 1 } else { o });
                prop_assert_eq!(index.position(back, e), p);
            }
        }
    }
}
