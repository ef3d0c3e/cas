use std::borrow::Cow;
use std::collections::HashSet;

use cas_core::Symbol;
use proptest::prelude::*;

use super::{
    Names, NamesError, NamesErrorKind as E, SyntaxForm, Target, TypeName, escape_name,
    unescape_name,
};
use crate::lexer::lex;
use crate::syntax::SyntaxKind;

fn planned(name: &str) -> Target {
    Target::Planned(name.into())
}

#[test]
fn english_names() {
    let names = Names::english();
    let cases = [
        ("sin", Target::Builtin(Symbol::SIN)),
        ("exp", Target::Builtin(Symbol::EXP)),
        ("ln", Target::Builtin(Symbol::LOG)),
        ("log", Target::Builtin(Symbol::LOG)),
        ("pi", Target::Builtin(Symbol::PI)),
        ("e", Target::Builtin(Symbol::E)),
        ("hold", Target::Builtin(Symbol::HOLD)),
        ("i", Target::Imaginary),
        ("call", Target::Syntax(SyntaxForm::Call)),
        ("assume", Target::Syntax(SyntaxForm::Assume)),
        ("purge", Target::Syntax(SyntaxForm::Purge)),
        ("diff", planned("D")),
        ("int", planned("Integrate")),
        ("integrate", planned("Integrate")),
        ("evalf", planned("N")),
        ("subst", Target::Builtin(Symbol::REPLACE_ALL)),
        ("replace", Target::Builtin(Symbol::REPLACE_ALL)),
        ("undef", Target::Builtin(Symbol::INDETERMINATE)),
        ("infinity", Target::Builtin(Symbol::COMPLEX_INFINITY)),
        ("inf", Target::Builtin(Symbol::INFINITY)),
        ("true", Target::Builtin(Symbol::TRUE)),
        ("sqrt", Target::Builtin(Symbol::SQRT)),
        ("abs", planned("Abs")),
    ];
    for (name, target) in cases {
        assert_eq!(names.target(name), Some(&target), "{name}");
        assert!(names.is_reserved(name));
    }
    assert_eq!(names.spelling(&Target::Builtin(Symbol::LOG)), Some("ln"));
    assert_eq!(names.spelling(&planned("Integrate")), Some("integrate"));
    assert_eq!(
        names.spelling(&Target::Builtin(Symbol::LENGTH)),
        Some("size")
    );
    assert_eq!(names.spelling(&Target::Builtin(Symbol::PLUS)), None);
    // User names, capitalised or not, and the aliases the lexer resolves.
    for name in ["x", "E", "D", "Sin", "N", "π", "∞", "Log"] {
        assert_eq!(names.target(name), None, "{name}");
        assert!(!names.is_reserved(name));
    }
    assert_eq!(names.type_name("integer"), Some(TypeName::Integer));
    assert_eq!(names.type_name("any"), Some(TypeName::Any));
    assert_eq!(names.type_spelling(TypeName::Matrix), Some("matrix"));
    // Types are a namespace of their own.
    assert_eq!(names.type_name("sin"), None);
    assert!(!names.is_reserved("integer"));
}

/// Every builtin has a surface name or is listed as internal, never both.
#[test]
fn every_builtin_is_named_or_internal() {
    let table = Names::english();
    let named: HashSet<Symbol> = table
        .names()
        .filter_map(|(_, t)| match t {
            Target::Builtin(s) => Some(*s),
            _ => None,
        })
        .collect();
    let internal: HashSet<Symbol> = table.internal().iter().copied().collect();
    for s in Symbol::builtins() {
        assert!(
            named.contains(&s) != internal.contains(&s),
            "{s:?} must have a name in names/en.txt or be listed in [internal], not both"
        );
    }
}

/// A planned builtin moves to `[builtins]` when cas-core defines it; until
/// then no builtin has its internal name.
#[test]
fn planned_builtins_are_not_defined() {
    let defined: HashSet<&str> = Symbol::builtins().map(Symbol::name).collect();
    for (name, target) in Names::english().names() {
        if let Target::Planned(internal) = target {
            assert!(!defined.contains(&**internal), "{name}: {internal}");
        }
    }
}

/// Every name is one identifier token in normal form C, and every target's
/// spelling is one of its names.
fn check_round_trip(names: &Names) {
    for (name, target) in names.names() {
        let lexed = lex(name);
        assert!(
            lexed.len() == 1 && lexed.kind(0) == SyntaxKind::Ident && lexed.errors().is_empty(),
            "{name}"
        );
        assert_eq!(crate::name(name), name);
        let spelling = names.spelling(target).unwrap();
        assert_eq!(names.target(spelling), Some(target), "{name}");
    }
    for (name, ty) in names.type_names() {
        let spelling = names.type_spelling(ty).unwrap();
        assert_eq!(names.type_name(spelling), Some(ty), "{name}");
    }
}

const TEST_LOCALE: &str = "
# A test locale: new spellings for a defined, a planned and a syntax target,
# a type, and an English name repeated to make it the spelling.
[builtins]
Sin: sinus
Log: logn ln
D: deriver derivee
Expand: developper

[syntax]
call: appel

[types]
integer: entier

[constants]
imaginary: i
";

#[test]
fn locales() {
    check_round_trip(Names::english());
    let names = Names::with_locale(TEST_LOCALE).unwrap();
    check_round_trip(&names);
    assert_eq!(names.target("sinus"), Some(&Target::Builtin(Symbol::SIN)));
    assert_eq!(names.target("sin"), Some(&Target::Builtin(Symbol::SIN)));
    assert_eq!(names.spelling(&Target::Builtin(Symbol::SIN)), Some("sinus"));
    assert_eq!(names.spelling(&Target::Builtin(Symbol::LOG)), Some("logn"));
    assert_eq!(names.target("derivee"), Some(&planned("D")));
    assert_eq!(names.spelling(&planned("D")), Some("deriver"));
    assert_eq!(
        names.target("appel"),
        Some(&Target::Syntax(SyntaxForm::Call))
    );
    assert_eq!(names.type_name("entier"), Some(TypeName::Integer));
    assert_eq!(names.type_name("integer"), Some(TypeName::Integer));
    assert_eq!(names.type_spelling(TypeName::Integer), Some("entier"));
    // Untouched targets keep their English spelling.
    assert_eq!(names.spelling(&Target::Builtin(Symbol::COS)), Some("cos"));
    // The English table itself is unchanged.
    assert_eq!(
        Names::english().spelling(&Target::Builtin(Symbol::SIN)),
        Some("sin")
    );
}

#[test]
fn malformed_locales() {
    let cases = [
        ("Sin: sinus", 1, E::OutsideSection),
        ("[words]\n", 1, E::UnknownSection("[words]".into())),
        ("[builtins", 1, E::UnknownSection("[builtins".into())),
        ("[planned]\n", 1, E::SectionNotInLocale("[planned]".into())),
        (
            "[internal]\n",
            1,
            E::SectionNotInLocale("[internal]".into()),
        ),
        ("[builtins]\nSin sinus", 2, E::MissingColon),
        ("[builtins]\nSin:", 2, E::NoNames),
        (
            "[builtins]\n\n# x\nSni: sinus",
            4,
            E::UnknownTarget("Sni".into()),
        ),
        // Internal builtins have no names.
        ("[builtins]\nPlus: plus", 2, E::UnknownTarget("Plus".into())),
        (
            "[types]\nreal: reel\nnatural: naturel",
            3,
            E::UnknownTarget("natural".into()),
        ),
        ("[constants]\npi: pie", 2, E::UnknownTarget("pi".into())),
        ("[syntax]\nlet: soit", 2, E::UnknownTarget("let".into())),
        ("[builtins]\nSin: 2x", 2, E::NotAName("2x".into())),
        ("[builtins]\nSin: if", 2, E::NotAName("if".into())),
        ("[builtins]\nSin: π", 2, E::NotAName("π".into())),
        (
            "[builtins]\nSin: e\u{301}",
            2,
            E::NotAName("e\u{301}".into()),
        ),
        ("[builtins]\nSin: a-b", 2, E::NotAName("a-b".into())),
        ("[builtins]\nCos: sin", 2, E::Conflict("sin".into())),
        ("[builtins]\nSin: s\nCos: s", 3, E::Repeated("s".into())),
        ("[builtins]\nSin: s t s", 2, E::Repeated("s".into())),
        ("[builtins]\nSin: s\nSin: t", 3, E::Repeated("Sin".into())),
        ("[types]\nreal: integer", 2, E::Conflict("integer".into())),
    ];
    for (text, line, kind) in cases {
        assert_eq!(
            Names::with_locale(text).unwrap_err(),
            NamesError { line, kind },
            "{text:?}"
        );
    }
    let error = Names::with_locale("[builtins]\nSni: s").unwrap_err();
    assert_eq!(error.to_string(), "line 2: unknown target Sni");
    // An empty overlay is the English table.
    let empty = Names::with_locale("").unwrap();
    assert_eq!(empty.names().count(), Names::english().names().count());
}

#[test]
fn malformed_base_tables() {
    let load = |text: &str| Names::empty().load(text, false).unwrap_err();
    assert_eq!(
        load("[planned]\nSin: sinus"),
        NamesError {
            line: 2,
            kind: E::AlreadyBuiltin("Sin".into())
        }
    );
    assert_eq!(
        load("[planned]\nlower: f"),
        NamesError {
            line: 2,
            kind: E::UnknownTarget("lower".into())
        }
    );
    assert_eq!(
        load("[internal]\nPlus Plus"),
        NamesError {
            line: 2,
            kind: E::Repeated("Plus".into())
        }
    );
    assert_eq!(
        load("[internal]\nNotABuiltin"),
        NamesError {
            line: 2,
            kind: E::UnknownTarget("NotABuiltin".into())
        }
    );
}

#[test]
fn escaping() {
    assert!(matches!(escape_name("x"), Cow::Borrowed("x")));
    assert!(matches!(escape_name("_X"), Cow::Borrowed("_X")));
    assert!(matches!(escape_name("Δx"), Cow::Borrowed("Δx")));
    assert_eq!(escape_name("E"), "$E");
    assert_eq!(escape_name("Delta"), "$Delta");
    assert_eq!(unescape_name("$E"), "E");
    assert_eq!(unescape_name("x"), "x");
    assert_eq!(unescape_name("Sin"), "Sin");
    // No escaped name is a builtin's.
    for s in Symbol::builtins() {
        assert_ne!(escape_name(s.name()), s.name(), "{s:?}");
    }
}

proptest! {
    #[test]
    fn escaping_round_trips(name in "[A-Za-z_θΔ][A-Za-z0-9_θΔ]{0,8}") {
        let escaped = escape_name(&name);
        prop_assert_eq!(unescape_name(&escaped), name.as_str());
        prop_assert_eq!(
            escaped.starts_with('$'),
            name.starts_with(|c: char| c.is_ascii_uppercase())
        );
    }
}
