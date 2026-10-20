//! The builtin symbols: fixed ids `0..COUNT` and fixed attributes.
//!
//! The symbol table is created with these names first, in this order, so a
//! builtin's id is its index here and its constant is known at compile time.
//! Builtin attributes never change (CAS-DOC-3, D3). Ids are not persisted, so
//! builtins can be added anywhere in the list.
//!
//! A name is an ASCII capital, then ASCII letters and digits. The language
//! stores a user's name that starts with a capital with a `$` in front (`E`
//! becomes `$E`), so it never meets a builtin's, now or later (CAS-DOC-7 §5).

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
    TIMES = "Times" [FLAT, ORDERLESS, ONE_IDENTITY, NUMERIC_FUNCTION, PROTECTED];
    POWER = "Power" [ONE_IDENTITY, NUMERIC_FUNCTION, PROTECTED];
    // Elementary functions.
    EXP = "Exp" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    LOG = "Log" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SIN = "Sin" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    COS = "Cos" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    TAN = "Tan" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SQRT = "Sqrt" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    FACTORIAL = "Factorial" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Element-wise operators, next to `Times` and `Power` (CAS-DOC-4 §6.2).
    ELEMENTWISE_TIMES = "ElementwiseTimes" [FLAT, ONE_IDENTITY, LISTABLE, PROTECTED];
    ELEMENTWISE_DIVIDE = "ElementwiseDivide" [LISTABLE, PROTECTED];
    ELEMENTWISE_POWER = "ElementwisePower" [LISTABLE, PROTECTED];
    // Constants.
    PI = "Pi" [CONSTANT, PROTECTED];
    E = "E" [CONSTANT, PROTECTED];
    // The imaginary unit's symbol: it evaluates to the complex number
    // `0 + 1·I`, which the builder and the language produce directly
    // (CAS-51).
    I = "I" [CONSTANT, PROTECTED];
    // Numbers: the head of a complex number in full form, `Complex[re, im]`.
    COMPLEX = "Complex" [PROTECTED];
    // Logical values, infinities and nothing (CAS-DOC-4 §6.1, R14).
    TRUE = "True" [PROTECTED];
    FALSE = "False" [PROTECTED];
    NULL = "Null" [PROTECTED];
    INDETERMINATE = "Indeterminate" [PROTECTED];
    INFINITY = "Infinity" [PROTECTED];
    COMPLEX_INFINITY = "ComplexInfinity" [PROTECTED];
    DIRECTED_INFINITY = "DirectedInfinity" [PROTECTED];
    // Equations, comparisons and logic (CAS-DOC-7 §9). `Equation` is the
    // language's `=`, which never becomes a boolean; `Equal` is `==`.
    EQUATION = "Equation" [PROTECTED];
    EQUAL = "Equal" [PROTECTED];
    UNEQUAL = "Unequal" [PROTECTED];
    LESS = "Less" [PROTECTED];
    LESS_EQUAL = "LessEqual" [PROTECTED];
    GREATER = "Greater" [PROTECTED];
    GREATER_EQUAL = "GreaterEqual" [PROTECTED];
    ELEMENT = "Element" [PROTECTED];
    AND = "And" [FLAT, ONE_IDENTITY, HOLD_ALL, PROTECTED];
    OR = "Or" [FLAT, ONE_IDENTITY, HOLD_ALL, PROTECTED];
    NOT = "Not" [PROTECTED];
    IF = "If" [HOLD_REST, PROTECTED];
    // Structure.
    LIST = "List" [PROTECTED];
    SEQUENCE = "Sequence" [PROTECTED];
    HOLD = "Hold" [HOLD_ALL, PROTECTED];
    RULE = "Rule" [PROTECTED];
    RULE_DELAYED = "RuleDelayed" [HOLD_REST, PROTECTED];
    FUNCTION = "Function" [HOLD_ALL, PROTECTED];
    RANGE = "Range" [PROTECTED];
    FINITE_SET = "FiniteSet" [ORDERLESS, PROTECTED];
    // Indexing from 0, and a copy with an element replaced (CAS-DOC-7 §7).
    AT = "At" [PROTECTED];
    REPLACE_AT = "ReplaceAt" [PROTECTED];
    REPLACE_ALL = "ReplaceAll" [PROTECTED];
    // Statements: the session runs them, evaluation leaves them alone
    // (CAS-DOC-4 R1, CAS-DOC-7 §4).
    SET = "Set" [HOLD_FIRST, PROTECTED];
    SET_DELAYED = "SetDelayed" [HOLD_ALL, PROTECTED];
    UP_SET_DELAYED = "UpSetDelayed" [HOLD_ALL, PROTECTED];
    ASSUME = "Assume" [HOLD_ALL, PROTECTED];
    CLEAR = "Clear" [HOLD_ALL, PROTECTED];
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
    // Type tests of pattern variables, `n :: integer` (CAS-DOC-7 §6).
    INTEGER_Q = "IntegerQ" [PROTECTED];
    RATIONAL_Q = "RationalQ" [PROTECTED];
    NUMBER_Q = "NumberQ" [PROTECTED];
    REAL_Q = "RealQ" [PROTECTED];
    COMPLEX_Q = "ComplexQ" [PROTECTED];
    LIST_Q = "ListQ" [PROTECTED];
    MATRIX_Q = "MatrixQ" [PROTECTED];
    STRING_Q = "StringQ" [PROTECTED];
    SYMBOL_Q = "SymbolQ" [PROTECTED];
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
