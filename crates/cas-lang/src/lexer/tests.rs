use proptest::prelude::*;

use super::lex;
use crate::LexErrorKind as E;
use crate::syntax::SyntaxKind as K;

/// The tokens of `text`, whitespace left out.
fn tokens(text: &str) -> Vec<(K, &str)> {
    lex(text)
        .iter()
        .filter(|(k, _)| *k != K::Whitespace)
        .collect()
}

/// Checks the tokens of `text`, whitespace left out, and that there are no
/// errors.
fn check(text: &str, expected: &[(K, &str)]) {
    check_errors(text, expected, &[]);
}

/// Checks the tokens of `text`, whitespace left out, and its errors as
/// `(covered text, kind)`.
fn check_errors(text: &str, expected: &[(K, &str)], errors: &[(&str, E)]) {
    assert_eq!(tokens(text), expected, "{text:?}");
    let lexed = lex(text);
    let found: Vec<_> = lexed
        .errors()
        .iter()
        .map(|e| (&text[e.range.clone()], e.kind))
        .collect();
    assert_eq!(found, errors, "{text:?}");
}

#[test]
fn trivia() {
    let all: Vec<_> = lex("a \t b\nc\r\nd\re").iter().collect();
    assert_eq!(
        all,
        [
            (K::Ident, "a"),
            (K::Whitespace, " \t "),
            (K::Ident, "b"),
            (K::Newline, "\n"),
            (K::Ident, "c"),
            (K::Newline, "\r\n"),
            (K::Ident, "d"),
            (K::Newline, "\r"),
            (K::Ident, "e"),
        ]
    );
    assert_eq!(
        lex("\n\n").iter().collect::<Vec<_>>(),
        [(K::Newline, "\n"), (K::Newline, "\n")]
    );
    check(
        "x // note\ny",
        &[
            (K::Ident, "x"),
            (K::LineComment, "// note"),
            (K::Newline, "\n"),
            (K::Ident, "y"),
        ],
    );
    check(
        "a /* b /* c */ d */ e",
        &[
            (K::Ident, "a"),
            (K::BlockComment, "/* b /* c */ d */"),
            (K::Ident, "e"),
        ],
    );
    check("/**/x", &[(K::BlockComment, "/**/"), (K::Ident, "x")]);
    check("\u{feff}x\u{a0}y", &[(K::Ident, "x"), (K::Ident, "y")]);
    check_errors(
        "/* a /* b */",
        &[(K::BlockComment, "/* a /* b */")],
        &[("/* a /* b */", E::UnterminatedComment)],
    );
    assert!(lex("").is_empty());
    assert!(K::Whitespace.is_trivia() && K::BlockComment.is_trivia());
    assert!(!K::Newline.is_trivia() && !K::Ident.is_trivia());
}

#[test]
fn names_and_keywords() {
    check(
        "x x_1 max_iter θ αβ _tmp Δx ℝ",
        &[
            (K::Ident, "x"),
            (K::Ident, "x_1"),
            (K::Ident, "max_iter"),
            (K::Ident, "θ"),
            (K::Ident, "αβ"),
            (K::Ident, "_tmp"),
            (K::Ident, "Δx"),
            (K::Ident, "ℝ"),
        ],
    );
    // `π` and `ℯ` never join a name; `∞` is not a letter.
    check("2πr", &[(K::Number, "2"), (K::Ident, "π"), (K::Ident, "r")]);
    check(
        "aπb ℯx ∞",
        &[
            (K::Ident, "a"),
            (K::Ident, "π"),
            (K::Ident, "b"),
            (K::Ident, "ℯ"),
            (K::Ident, "x"),
            (K::Ident, "∞"),
        ],
    );
    // The middle dot is a product, also between letters.
    check("a·b", &[(K::Ident, "a"), (K::Star, "·"), (K::Ident, "b")]);
    // A combining accent stays in the name.
    check("e\u{301}", &[(K::Ident, "e\u{301}")]);
    let keywords = "if then elif else end for in do while begin return break continue \
                    local when and or not";
    let kinds: Vec<_> = tokens(keywords).into_iter().map(|(k, _)| k).collect();
    assert_eq!(
        kinds,
        [
            K::IfKw,
            K::ThenKw,
            K::ElifKw,
            K::ElseKw,
            K::EndKw,
            K::ForKw,
            K::InKw,
            K::DoKw,
            K::WhileKw,
            K::BeginKw,
            K::ReturnKw,
            K::BreakKw,
            K::ContinueKw,
            K::LocalKw,
            K::WhenKw,
            K::AndKw,
            K::OrKw,
            K::NotKw,
        ]
    );
    assert!(kinds.iter().all(|k| k.is_keyword()));
    assert!(!K::Ident.is_keyword() && !K::Plus.is_keyword());
    // Keywords are case-sensitive and whole words.
    check(
        "If iff ends",
        &[(K::Ident, "If"), (K::Ident, "iff"), (K::Ident, "ends")],
    );
    check("x ∈ S", &[(K::Ident, "x"), (K::InKw, "∈"), (K::Ident, "S")]);
}

#[test]
fn numbers() {
    check(
        "42 1_000_000 0x1f 0xFF 0o17 0b101 007 0x1e5",
        &[
            (K::Number, "42"),
            (K::Number, "1_000_000"),
            (K::Number, "0x1f"),
            (K::Number, "0xFF"),
            (K::Number, "0o17"),
            (K::Number, "0b101"),
            (K::Number, "007"),
            (K::Number, "0x1e5"),
        ],
    );
    check(
        "0.1 1.5e-3 2E10 1e+5 6.02e23 1_0.0_1e1_0",
        &[
            (K::Number, "0.1"),
            (K::Number, "1.5e-3"),
            (K::Number, "2E10"),
            (K::Number, "1e+5"),
            (K::Number, "6.02e23"),
            (K::Number, "1_0.0_1e1_0"),
        ],
    );
    // A prefix needs a digit of its base after it, an exponent a digit after
    // its sign.
    check(
        "0x^2",
        &[
            (K::Number, "0"),
            (K::Ident, "x"),
            (K::Caret, "^"),
            (K::Number, "2"),
        ],
    );
    check(
        "0b 0o8",
        &[
            (K::Number, "0"),
            (K::Ident, "b"),
            (K::Number, "0"),
            (K::Ident, "o8"),
        ],
    );
    check("2e", &[(K::Number, "2"), (K::Ident, "e")]);
    check(
        "2e-x",
        &[
            (K::Number, "2"),
            (K::Ident, "e"),
            (K::Minus, "-"),
            (K::Ident, "x"),
        ],
    );
    check("1e5x", &[(K::Number, "1e5"), (K::Ident, "x")]);
    check("2x", &[(K::Number, "2"), (K::Ident, "x")]);
    // The point needs a digit after it.
    check(
        "2.*x",
        &[(K::Number, "2"), (K::DotStar, ".*"), (K::Ident, "x")],
    );
    check(
        "1..2",
        &[(K::Number, "1"), (K::DotDot, ".."), (K::Number, "2")],
    );
    check(
        "1./2 1.^2",
        &[
            (K::Number, "1"),
            (K::DotSlash, "./"),
            (K::Number, "2"),
            (K::Number, "1"),
            (K::DotCaret, ".^"),
            (K::Number, "2"),
        ],
    );
}

#[test]
fn malformed_numbers() {
    check_errors(
        ".5 .5e3",
        &[(K::Number, ".5"), (K::Number, ".5e3")],
        &[(".5", E::LeadingPoint), (".5e3", E::LeadingPoint)],
    );
    check_errors(
        "1.",
        &[(K::Number, "1"), (K::Error, ".")],
        &[(".", E::LoneDot)],
    );
    check_errors(
        "1_",
        &[(K::Number, "1_")],
        &[("1_", E::MisplacedUnderscore)],
    );
    check_errors(
        "1__0",
        &[(K::Number, "1__0")],
        &[("1__0", E::MisplacedUnderscore)],
    );
    check_errors(
        "2_x",
        &[(K::Number, "2_"), (K::Ident, "x")],
        &[("2_", E::MisplacedUnderscore)],
    );
    check_errors(
        "1.5_e3",
        &[(K::Number, "1.5_e3")],
        &[("5_", E::MisplacedUnderscore)],
    );
    check_errors(
        "0b102 0o78",
        &[(K::Number, "0b102"), (K::Number, "0o78")],
        &[
            ("2", E::InvalidDigit { radix: 2 }),
            ("8", E::InvalidDigit { radix: 8 }),
        ],
    );
}

#[test]
fn strings() {
    check(
        r#""abc" "" "a\"b" "\\" "\n\t" "\u{3c0}" "two
lines""#,
        &[
            (K::String, r#""abc""#),
            (K::String, r#""""#),
            (K::String, r#""a\"b""#),
            (K::String, r#""\\""#),
            (K::String, r#""\n\t""#),
            (K::String, r#""\u{3c0}""#),
            (K::String, "\"two\nlines\""),
        ],
    );
    for bad in [
        r"\q",
        r"\u{110000}",
        r"\u{}",
        r"\u{1234567}",
        r"\u{d800}",
        r"\u",
    ] {
        let text = format!("\"{bad}\"");
        check_errors(&text, &[(K::String, &text)], &[(bad, E::InvalidEscape)]);
    }
    // `\u41` is the bad escape `\u` and the text `41`.
    check_errors(
        r#""\u41""#,
        &[(K::String, r#""\u41""#)],
        &[(r"\u", E::InvalidEscape)],
    );
    check_errors(
        "\"abc",
        &[(K::String, "\"abc")],
        &[("\"abc", E::UnterminatedString)],
    );
    // Errors are ordered by where they start.
    check_errors(
        r#""a\q"#,
        &[(K::String, r#""a\q"#)],
        &[
            (r#""a\q"#, E::UnterminatedString),
            (r"\q", E::InvalidEscape),
        ],
    );
    check_errors(
        "\"\\",
        &[(K::String, "\"\\")],
        &[("\"\\", E::UnterminatedString), ("\\", E::InvalidEscape)],
    );
}

#[test]
fn operators_and_aliases() {
    check(
        "+ - * / ^ .* ./ .^ ! = == != < <= > >= .. ... :: := -> ( ) [ ] { } , ;",
        &[
            (K::Plus, "+"),
            (K::Minus, "-"),
            (K::Star, "*"),
            (K::Slash, "/"),
            (K::Caret, "^"),
            (K::DotStar, ".*"),
            (K::DotSlash, "./"),
            (K::DotCaret, ".^"),
            (K::Bang, "!"),
            (K::Eq, "="),
            (K::EqEq, "=="),
            (K::Neq, "!="),
            (K::Lt, "<"),
            (K::Le, "<="),
            (K::Gt, ">"),
            (K::Ge, ">="),
            (K::DotDot, ".."),
            (K::Ellipsis, "..."),
            (K::ColonColon, "::"),
            (K::ColonEq, ":="),
            (K::Arrow, "->"),
            (K::LParen, "("),
            (K::RParen, ")"),
            (K::LBracket, "["),
            (K::RBracket, "]"),
            (K::LBrace, "{"),
            (K::RBrace, "}"),
            (K::Comma, ","),
            (K::Semicolon, ";"),
        ],
    );
    check(
        "− × · ⋅ ÷ ≠ ≤ ≥ → √ ∛ ∜",
        &[
            (K::Minus, "−"),
            (K::Star, "×"),
            (K::Star, "·"),
            (K::Star, "⋅"),
            (K::Slash, "÷"),
            (K::Neq, "≠"),
            (K::Le, "≤"),
            (K::Ge, "≥"),
            (K::Arrow, "→"),
            (K::Sqrt, "√"),
            (K::Cbrt, "∛"),
            (K::FourthRoot, "∜"),
        ],
    );
    // Longest match: `n!=3` is `n != 3`, as the language server points out.
    check("n!=3", &[(K::Ident, "n"), (K::Neq, "!="), (K::Number, "3")]);
    check(
        "n! = 3",
        &[
            (K::Ident, "n"),
            (K::Bang, "!"),
            (K::Eq, "="),
            (K::Number, "3"),
        ],
    );
    check(
        "a<-b",
        &[
            (K::Ident, "a"),
            (K::Lt, "<"),
            (K::Minus, "-"),
            (K::Ident, "b"),
        ],
    );
    check("===", &[(K::EqEq, "=="), (K::Eq, "=")]);
    check(
        "√2x xs...",
        &[
            (K::Sqrt, "√"),
            (K::Number, "2"),
            (K::Ident, "x"),
            (K::Ident, "xs"),
            (K::Ellipsis, "..."),
        ],
    );
    check_errors(
        "....",
        &[(K::Ellipsis, "..."), (K::Error, ".")],
        &[(".", E::LoneDot)],
    );
}

#[test]
fn superscripts() {
    check(
        "x² x⁻¹ (a+b)¹⁰ x⁺³ x⁰⁴⁵⁶⁷⁸⁹",
        &[
            (K::Ident, "x"),
            (K::Superscript, "²"),
            (K::Ident, "x"),
            (K::Superscript, "⁻¹"),
            (K::LParen, "("),
            (K::Ident, "a"),
            (K::Plus, "+"),
            (K::Ident, "b"),
            (K::RParen, ")"),
            (K::Superscript, "¹⁰"),
            (K::Ident, "x"),
            (K::Superscript, "⁺³"),
            (K::Ident, "x"),
            (K::Superscript, "⁰⁴⁵⁶⁷⁸⁹"),
        ],
    );
    check(
        "x²y",
        &[(K::Ident, "x"), (K::Superscript, "²"), (K::Ident, "y")],
    );
    check_errors(
        "x⁻ x²⁻",
        &[
            (K::Ident, "x"),
            (K::Superscript, "⁻"),
            (K::Ident, "x"),
            (K::Superscript, "²"),
            (K::Superscript, "⁻"),
        ],
        &[
            ("⁻", E::MissingSuperscriptDigits),
            ("⁻", E::MissingSuperscriptDigits),
        ],
    );
}

#[test]
fn error_tokens() {
    check_errors(
        "f'(x)",
        &[
            (K::Ident, "f"),
            (K::Error, "'"),
            (K::LParen, "("),
            (K::Ident, "x"),
            (K::RParen, ")"),
        ],
        &[("'", E::Apostrophe)],
    );
    check_errors(
        "|x| a && b # c",
        &[
            (K::Error, "|"),
            (K::Ident, "x"),
            (K::Error, "|"),
            (K::Ident, "a"),
            (K::Error, "&"),
            (K::Error, "&"),
            (K::Ident, "b"),
            (K::Error, "#"),
            (K::Ident, "c"),
        ],
        &[
            ("|", E::Bar),
            ("|", E::Bar),
            ("&", E::Ampersand),
            ("&", E::Ampersand),
            ("#", E::Hash),
        ],
    );
    check_errors(
        "a : b . c",
        &[
            (K::Ident, "a"),
            (K::Error, ":"),
            (K::Ident, "b"),
            (K::Error, "."),
            (K::Ident, "c"),
        ],
        &[(":", E::LoneColon), (".", E::LoneDot)],
    );
    check_errors(
        "@\u{7}\u{301}",
        &[(K::Error, "@"), (K::Error, "\u{7}"), (K::Error, "\u{301}")],
        &[
            ("@", E::UnexpectedChar('@')),
            ("\u{7}", E::UnexpectedChar('\u{7}')),
            ("\u{301}", E::UnexpectedChar('\u{301}')),
        ],
    );
    // Every kind has a message, and the hints name the replacement.
    assert!(E::Apostrophe.to_string().contains("diff"));
    assert!(E::Bar.to_string().contains("abs"));
    assert!(E::Ampersand.to_string().contains("and"));
    assert!(E::Hash.to_string().contains("//"));
    assert_eq!(
        E::UnexpectedChar('@').to_string(),
        "unexpected character '@'"
    );
}

/// Pieces that start, end or join tokens, mixed with any character, so that
/// concatenations hit sequences such as `*/x`, `0x1`, `2e-3` and `\u{41}`.
fn fragment() -> impl Strategy<Value = String> {
    const FRAGMENTS: &[&str] = &[
        "0", "1", "9", "_", ".", "..", "...", "e", "E", "e-", "x", "0x", "0b", "0o", "a", "f", "π",
        "ℯ", "∞", "·", "\"", "\\", "\\u{", "41", "}", "{", "/", "*", "/*", "*/", "//", "\n", "\r",
        "\r\n", " ", "\t", "²", "⁻", "⁺", "-", "+", ":", "=", "!", "<", ">", "'", "|", "#", "é",
        "\u{301}", "\u{feff}", "√", "×", "≤", "→", "∈", "(", ")", "[", "]", ",", ";", "if", "end",
    ];
    prop_oneof![
        4 => prop::sample::select(FRAGMENTS).prop_map(String::from),
        1 => any::<char>().prop_map(String::from),
    ]
}

/// The tokens cover the text exactly; errors are ordered, in bounds and on
/// character boundaries; every error token has an error over it; and each
/// token lexes alone as one token of the same kind.
fn check_lossless(text: &str) {
    let lexed = lex(text);
    let mut end = 0;
    for i in 0..lexed.len() {
        let range = lexed.range(i);
        assert_eq!(range.start, end, "{text:?}");
        assert!(range.start < range.end, "{text:?}");
        end = range.end;
        let alone = lex(lexed.text(i));
        assert_eq!(alone.len(), 1, "{text:?}: {:?}", lexed.text(i));
        assert_eq!(
            alone.kind(0),
            lexed.kind(i),
            "{text:?}: {:?}",
            lexed.text(i)
        );
        if lexed.kind(i) == K::Error {
            assert!(lexed.errors().iter().any(|e| e.range == range), "{text:?}");
        }
    }
    assert_eq!(end, text.len());
    assert_eq!(lexed.iter().map(|(_, t)| t).collect::<String>(), text);
    let mut last = 0;
    for e in lexed.errors() {
        assert!(
            last <= e.range.start && e.range.start < e.range.end,
            "{text:?}"
        );
        assert!(e.range.end <= text.len(), "{text:?}");
        assert!(text.is_char_boundary(e.range.start) && text.is_char_boundary(e.range.end));
        last = e.range.start;
    }
}

proptest! {
    #[test]
    fn lossless(fragments in prop::collection::vec(fragment(), 0..40)) {
        check_lossless(&fragments.concat());
    }

    #[test]
    fn lossless_any_text(text in any::<String>()) {
        check_lossless(&text);
    }
}
