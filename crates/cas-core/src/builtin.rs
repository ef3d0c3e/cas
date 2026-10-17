//! The builtin symbols: fixed ids `0..COUNT` and fixed attributes.
//!
//! The symbol table is created with these names first, in this order, so a
//! builtin's id is its index here and its constant is known at compile time.
//! Builtin attributes never change (CAS-DOC-3, D3). Ids are not persisted, so
//! builtins can be added anywhere in the list.

use crate::attributes::Attributes;
use crate::symbol::Symbol;

pub(crate) struct Builtin {
    pub(crate) name: &'static str,
    pub(crate) attributes: Attributes,
}

macro_rules! builtins {
    ($($konst:ident = $name:literal [$($attr:ident),*];)*) => {
        #[allow(non_camel_case_types, clippy::upper_case_acronyms)]
        #[repr(u32)]
        enum Id { $($konst,)* }

        /// The builtin symbols. Every builtin is `PROTECTED`.
        impl Symbol {
            $(
                #[doc = concat!("`", $name, "`: ", stringify!($($attr),*), ".")]
                pub const $konst: Symbol = Symbol::from_id(Id::$konst as u32);
            )*
        }

        pub(crate) const COUNT: usize = [$($name),*].len();

        pub(crate) static BUILTINS: [Builtin; COUNT] = [$(
            Builtin {
                name: $name,
                attributes: Attributes::NONE $(.union(Attributes::$attr))*,
            },
        )*];
    };
}

builtins! {
    // Arithmetic.
    PLUS = "Plus" [FLAT, ORDERLESS, ONE_IDENTITY, LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    TIMES = "Times" [FLAT, ORDERLESS, ONE_IDENTITY, LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    POWER = "Power" [ONE_IDENTITY, LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Elementary functions.
    EXP = "Exp" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    LOG = "Log" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SIN = "Sin" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    COS = "Cos" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    TAN = "Tan" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Constants.
    PI = "Pi" [CONSTANT, PROTECTED];
    E = "E" [CONSTANT, PROTECTED];
    // Structure.
    LIST = "List" [PROTECTED];
    SEQUENCE = "Sequence" [PROTECTED];
    HOLD = "Hold" [HOLD_ALL, PROTECTED];
    RULE = "Rule" [PROTECTED];
    RULE_DELAYED = "RuleDelayed" [HOLD_REST, PROTECTED];
    // Patterns.
    BLANK = "Blank" [PROTECTED];
    BLANK_SEQUENCE = "BlankSequence" [PROTECTED];
    BLANK_NULL_SEQUENCE = "BlankNullSequence" [PROTECTED];
    PATTERN = "Pattern" [HOLD_FIRST, PROTECTED];
    ALTERNATIVES = "Alternatives" [PROTECTED];
    OPTIONAL = "Optional" [PROTECTED];
    CONDITION = "Condition" [HOLD_ALL, PROTECTED];
    PATTERN_TEST = "PatternTest" [HOLD_REST, PROTECTED];
    HOLD_PATTERN = "HoldPattern" [HOLD_ALL, PROTECTED];
}

pub(crate) fn get(s: Symbol) -> Option<&'static Builtin> {
    BUILTINS.get(s.id() as usize)
}

/// Whether an application with head `s` is a pattern construct, so that an
/// expression containing one may match more than itself. Over-approximates:
/// the matcher adds heads here as it supports them.
pub(crate) fn is_pattern_head(s: Symbol) -> bool {
    matches!(
        s,
        Symbol::BLANK
            | Symbol::BLANK_SEQUENCE
            | Symbol::BLANK_NULL_SEQUENCE
            | Symbol::PATTERN
            | Symbol::ALTERNATIVES
            | Symbol::OPTIONAL
            | Symbol::CONDITION
            | Symbol::PATTERN_TEST
    )
}
