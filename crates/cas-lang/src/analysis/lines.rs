//! Lines and columns (CAS-DOC-60 §2). Line breaks are `\r\n`, `\n`, and a
//! `\r` not followed by `\n`, as the lexer has them; a line's content runs
//! from its start to its break. Columns count bytes ([`Encoding::Utf8`]) or
//! UTF-16 code units ([`Encoding::Utf16`]) from the line's start.
//!
//! The index keeps the line starts, the break widths and the non-ASCII
//! characters, so a conversion is a binary search plus a walk over one
//! line's non-ASCII characters, and the text is not borrowed. Offsets are
//! kept as `usize`; a line or column that does not fit a `u32` (texts past
//! 4 GiB, out of scope) saturates instead of panicking.

/// How columns are counted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    /// Bytes of UTF-8.
    Utf8,
    /// UTF-16 code units: 2 for a character above U+FFFF, else 1.
    Utf16,
}

/// A line and a column in it, both from 0; the column counts in an
/// [`Encoding`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Position {
    /// The line, from 0.
    pub line: u32,
    /// The column from the line's start, from 0.
    pub character: u32,
}

/// A non-ASCII character: its byte offset and its length in UTF-8 (2 to 4).
#[derive(Clone, Copy, Debug)]
struct Wide {
    start: usize,
    len: u8,
}

impl Wide {
    fn end(self) -> usize {
        self.start + usize::from(self.len)
    }

    /// Its length in UTF-16 code units.
    fn len_utf16(self) -> usize {
        if self.len == 4 { 2 } else { 1 }
    }
}

/// Byte offsets of a text to and from [`Position`]s, built once in O(n).
#[derive(Clone, Debug)]
pub struct LineIndex {
    /// The text's length in bytes.
    len: usize,
    /// Each line's start; the first is 0.
    starts: Vec<usize>,
    /// Each line's break width in bytes: 1 or 2, and 0 for the last line.
    breaks: Vec<u8>,
    /// The non-ASCII characters, in text order.
    wide: Vec<Wide>,
}

impl LineIndex {
    /// The index of `text`.
    #[must_use]
    pub fn new(text: &str) -> LineIndex {
        let bytes = text.as_bytes();
        let mut starts = vec![0];
        let mut breaks = Vec::new();
        let mut wide = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            let b = bytes[i];
            if b.is_ascii() {
                let width = match b {
                    b'\r' if bytes.get(i + 1) == Some(&b'\n') => 2,
                    b'\n' | b'\r' => 1,
                    _ => 0,
                };
                if width == 0 {
                    i += 1;
                } else {
                    i += usize::from(width);
                    breaks.push(width);
                    starts.push(i);
                }
            } else {
                // A leading byte: its count of leading ones is the length.
                let len = b.leading_ones().clamp(2, 4) as u8;
                wide.push(Wide { start: i, len });
                i += usize::from(len);
            }
        }
        breaks.push(0);
        LineIndex {
            len: text.len(),
            starts,
            breaks,
            wide,
        }
    }

    /// The text's length in bytes.
    #[allow(clippy::len_without_is_empty)]
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// The number of lines: 1 plus the number of line breaks.
    #[must_use]
    pub fn line_count(&self) -> u32 {
        saturate(self.starts.len())
    }

    /// The position of `offset`. Past the end, the end; inside a character,
    /// its start; on the `\n` of a `\r\n`, the end of the line's content.
    #[must_use]
    pub fn position(&self, offset: usize, encoding: Encoding) -> Position {
        let offset = self.floor(offset.min(self.len));
        let line = self.starts.partition_point(|&s| s <= offset) - 1;
        let start = self.starts[line];
        let offset = offset.min(self.content_end(line));
        let bytes = offset - start;
        let character = match encoding {
            Encoding::Utf8 => bytes,
            Encoding::Utf16 => {
                let first = self.wide.partition_point(|w| w.start < start);
                let saved: usize = self.wide[first..]
                    .iter()
                    .take_while(|w| w.start < offset)
                    .map(|w| usize::from(w.len) - w.len_utf16())
                    .sum();
                bytes - saved
            }
        };
        Position {
            line: saturate(line),
            character: saturate(character),
        }
    }

    /// The offset of `position`. A line at or past [`line_count`] gives the
    /// text's end; a column past the line's content its end; a column inside
    /// a character that character's start. Never inside a line break.
    ///
    /// [`line_count`]: LineIndex::line_count
    #[must_use]
    pub fn offset(&self, position: Position, encoding: Encoding) -> usize {
        let line = position.line as usize;
        if line >= self.starts.len() {
            return self.len;
        }
        let start = self.starts[line];
        let end = self.content_end(line);
        let character = position.character as usize;
        match encoding {
            Encoding::Utf8 => self.floor(start.saturating_add(character).min(end)),
            Encoding::Utf16 => {
                let first = self.wide.partition_point(|w| w.start < start);
                let mut at = start;
                let mut left = character;
                for &w in self.wide[first..].iter().take_while(|w| w.start < end) {
                    let gap = w.start - at;
                    if left <= gap {
                        return at + left;
                    }
                    left -= gap;
                    at = w.start;
                    if left < w.len_utf16() {
                        return at;
                    }
                    left -= w.len_utf16();
                    at = w.end();
                }
                at + left.min(end - at)
            }
        }
    }

    /// The end of `line`'s content: the start of its break.
    fn content_end(&self, line: usize) -> usize {
        let next = self.starts.get(line + 1).copied().unwrap_or(self.len);
        next - usize::from(self.breaks[line])
    }

    /// `offset` (at most the length) rounded down to a character boundary.
    fn floor(&self, offset: usize) -> usize {
        let i = self.wide.partition_point(|w| w.start <= offset);
        match i.checked_sub(1).map(|i| self.wide[i]) {
            Some(w) if offset < w.end() => w.start,
            _ => offset,
        }
    }
}

/// `n` as a `u32`, saturating.
fn saturate(n: usize) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

#[cfg(test)]
mod tests;
