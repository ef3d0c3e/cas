//! Queries over source text and its tree for editors (CAS-DOC-60). Every
//! query takes and gives byte offsets (`usize`); an offset past the text is
//! the text's end, one inside a character is rounded down to its start, and
//! no query panics on any offset. Only [`lines`] converts to lines and
//! columns, at the edge.
//!
//! - [`lines`]: [`LineIndex`](lines::LineIndex), byte offsets to and from
//!   [`Position`](lines::Position)s (line and column) in UTF-8 or UTF-16
//!   ([`Encoding`](lines::Encoding)), as LSP 3.17 counts them.
//! - [`cursor`]: what is at an offset ([`Cursor`](cursor::Cursor)): the token
//!   ([`token_at`](cursor::token_at)), the kind of spot
//!   ([`CursorKind`](cursor::CursorKind)) and the call around it
//!   ([`Call`](cursor::Call)).
//! - [`symbols`]: the names of a text ([`Symbols`](symbols::Symbols)), each
//!   [`Occurrence`](symbols::Occurrence) resolved to a
//!   [`SymbolInfo`](symbols::SymbolInfo) (a reserved target, a type, a
//!   global, a parameter or a local) with its definitions and uses, and the
//!   [`Scope`](symbols::Scope)s that bind them, without lowering.

pub mod cursor;
pub mod lines;
pub mod symbols;
