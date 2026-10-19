# cas-lang

The surface language of the CAS: one math-like grammar for students and programs. The design is
in Lific page `CAS-DOC-7`; the implementation plan is `CAS-PLAN-11`.

## Layout

| Module | Contents |
|---|---|
| `syntax` | `SyntaxKind`: the token kinds (the parser adds node kinds). Aliases share a kind with what they stand for: `×` is a `Star`, `≤` a `Le`, `∈` the keyword `in`, each with its own text. `is_trivia` (whitespace and comments; line breaks are not trivia, the parser decides which ones end a statement), `is_keyword`, `keyword`. |
| `lexer` | `lex`: text to `Lexed` tokens, lossless (the token texts concatenate back to the text, so the language server can map every byte) and total (a character that starts no token becomes an `Error` token). `LexError` carries a byte range and a `LexErrorKind` whose message names the replacement where there is one (`'` → `diff`, `\|` → `abs`, `&` → `and`, `#` → `//`). Numbers take `_` between digits, `0x`/`0o`/`0b` only before a digit of their base (`0x^2` is `0·x^2`), a fraction only before a digit (`2.*x` is `2 .* x`, `1..2` a range) and an exponent only before a digit (`2e` is `2·e`). `π` and `ℯ` never join a longer name (`2πr` is `2·π·r`), `·` is a product even between letters, block comments nest, and superscript digits with an optional sign are an exponent token (`x²`, `x⁻¹`). |
| `literal` | The values of literal tokens. `number_value` is exact (`0.1` is `1/10`, `1.5e-3` is `3/2000`), with exponents bounded by `MAX_DECIMAL_EXPONENT` (10⁶) so a short literal can't become gigabytes; `string_value` replaces the escapes `\"`, `\\`, `\n`, `\t` and `\u{…}`; `superscript_value` gives the exponent as an `Integer`; `name` puts a name in Unicode normal form C (borrowed when it already is) and spells `π`, `ℯ`, `∞` as `pi`, `e`, `inf`. Each accepts exactly the texts the lexer gives as one error-free token of its kind. |

## Testing

Lexer tests are in `src/lexer/tests.rs`: golden token tables for trivia and nested comments, names
and keywords (case-sensitive, `π` standing alone, combining accents), every number form and the
malformed ones, strings and every bad escape, every operator and Unicode alias, superscripts, and
the error tokens with their messages. Proptests over arbitrary text and over text built from the
characters that start and end tokens check that the tokens cover the input exactly, that each token
lexes alone as one token of the same kind, and that errors are ordered and on character boundaries.
Literal tests are in `src/literal/tests.rs`: value tables, the exponent bound, and proptests that
decimals with separators and exponents equal a digit-at-a-time reference, that `0x`/`0o`/`0b`
literals read back any `u64`, that any string escaped in every allowed way reads back, and that
superscripts read back any `i64`. The `lang_lex` fuzz target (`fuzz/`) runs the same checks on
arbitrary bytes.
