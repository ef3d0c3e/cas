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
    ARC_SIN = "ArcSin" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_COS = "ArcCos" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_TAN = "ArcTan" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Hyperbolic functions (CAS-DOC-15 §2).
    SINH = "Sinh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    COSH = "Cosh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    TANH = "Tanh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_SINH = "ArcSinh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_COSH = "ArcCosh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_TANH = "ArcTanh" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Reciprocal trigonometric functions, rewritten on evaluation
    // (CAS-DOC-15 §3).
    SEC = "Sec" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    CSC = "Csc" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    COT = "Cot" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_SEC = "ArcSec" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_CSC = "ArcCsc" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARC_COT = "ArcCot" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Complex parts (CAS-DOC-15 §4).
    RE = "Re" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    IM = "Im" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    CONJUGATE = "Conjugate" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ARG = "Arg" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Rounding, sign and order (CAS-DOC-15 §5).
    FLOOR = "Floor" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    CEILING = "Ceiling" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ROUND = "Round" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SIGN = "Sign" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    MAX = "Max" [FLAT, ORDERLESS, ONE_IDENTITY, NUMERIC_FUNCTION, PROTECTED];
    MIN = "Min" [FLAT, ORDERLESS, ONE_IDENTITY, NUMERIC_FUNCTION, PROTECTED];
    // Special functions (CAS-DOC-15 §6–§8).
    GAMMA = "Gamma" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    BETA = "Beta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ZETA = "Zeta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ERF = "Erf" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ERFC = "Erfc" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // The step and the delta (CAS-DOC-18 §2).
    HEAVISIDE_THETA = "HeavisideTheta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    DIRAC_DELTA = "DiracDelta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SQRT = "Sqrt" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    FACTORIAL = "Factorial" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Calculus (CAS-DOC-4 §8): the derivative, as system rules, and the
    // derivative of a function, `Derivative[n][f]`; integrals (CAS-DOC-9).
    D = "D" [PROTECTED];
    DERIVATIVE = "Derivative" [PROTECTED];
    INTEGRATE = "Integrate" [PROTECTED];
    // Series and limits (CAS-DOC-10): `OrderSize[u]`, giac's `order_size`,
    // a bounded function of `u`, inert.
    SERIES = "Series" [PROTECTED];
    LIMIT = "Limit" [PROTECTED];
    ORDER_SIZE = "OrderSize" [PROTECTED];
    // Differential equations (CAS-DOC-16).
    DSOLVE = "DSolve" [PROTECTED];
    // Laplace transforms (CAS-DOC-18).
    LAPLACE_TRANSFORM = "LaplaceTransform" [PROTECTED];
    INVERSE_LAPLACE_TRANSFORM = "InverseLaplaceTransform" [PROTECTED];
    // Equations (CAS-DOC-11): the real solutions, the complex ones, and
    // numeric roots.
    SOLVE = "Solve" [PROTECTED];
    CSOLVE = "CSolve" [PROTECTED];
    FSOLVE = "FSolve" [PROTECTED];
    // Linear algebra (CAS-DOC-12), next to the products of lists.
    TRANSPOSE = "Transpose" [PROTECTED];
    IDENTITY_MATRIX = "IdentityMatrix" [PROTECTED];
    TR = "Tr" [PROTECTED];
    DOT = "Dot" [PROTECTED];
    CROSS = "Cross" [PROTECTED];
    DET = "Det" [PROTECTED];
    INVERSE = "Inverse" [PROTECTED];
    ROW_REDUCE = "RowReduce" [PROTECTED];
    MATRIX_RANK = "MatrixRank" [PROTECTED];
    NULL_SPACE = "NullSpace" [PROTECTED];
    COLUMN_SPACE = "ColumnSpace" [PROTECTED];
    LIN_SOLVE = "LinSolve" [PROTECTED];
    CHARACTERISTIC_POLYNOMIAL = "CharacteristicPolynomial" [PROTECTED];
    EIGENVALUES = "Eigenvalues" [PROTECTED];
    EIGENVECTORS = "Eigenvectors" [PROTECTED];
    // Sums and products (CAS-DOC-13).
    SUM = "Sum" [PROTECTED];
    PRODUCT = "Product" [PROTECTED];
    BINOMIAL = "Binomial" [LISTABLE, PROTECTED];
    // Number theory (CAS-DOC-17).
    PRIME_Q = "PrimeQ" [LISTABLE, PROTECTED];
    NEXT_PRIME = "NextPrime" [LISTABLE, PROTECTED];
    PREVIOUS_PRIME = "PreviousPrime" [LISTABLE, PROTECTED];
    INTEGER_FACTOR = "IntegerFactor" [LISTABLE, PROTECTED];
    INTEGER_FACTOR_LIST = "IntegerFactorList" [LISTABLE, PROTECTED];
    DIVISORS = "Divisors" [LISTABLE, PROTECTED];
    EULER_PHI = "EulerPhi" [LISTABLE, PROTECTED];
    POWER_MOD = "PowerMod" [LISTABLE, PROTECTED];
    EXTENDED_GCD = "ExtendedGCD" [PROTECTED];
    CHINESE_REMAINDER = "ChineseRemainder" [PROTECTED];
    JACOBI_SYMBOL = "JacobiSymbol" [LISTABLE, PROTECTED];
    LEGENDRE_SYMBOL = "LegendreSymbol" [LISTABLE, PROTECTED];
    BERNOULLI_B = "BernoulliB" [LISTABLE, PROTECTED];
    PRIME = "Prime" [LISTABLE, PROTECTED];
    PRIME_PI = "PrimePi" [LISTABLE, PROTECTED];
    // Statistics (CAS-DOC-22): Xcas's `mean`, `variance`, `stddev` and
    // `stddevp`. giac's `variance` and `stddev` are the population ones and
    // Mathematica's `Variance` and `StandardDeviation` the sample ones, hence
    // the names; `stddevp` is the sample deviation.
    MEAN = "Mean" [PROTECTED];
    POPULATION_VARIANCE = "PopulationVariance" [PROTECTED];
    POPULATION_STANDARD_DEVIATION = "PopulationStandardDeviation" [PROTECTED];
    STANDARD_DEVIATION = "StandardDeviation" [PROTECTED];
    // Its order statistics: `median`, `quantile`, `quartiles`, `quartile1`
    // and `quartile3`.
    MEDIAN = "Median" [PROTECTED];
    QUANTILE = "Quantile" [PROTECTED];
    QUARTILES = "Quartiles" [PROTECTED];
    FIRST_QUARTILE = "FirstQuartile" [PROTECTED];
    THIRD_QUARTILE = "ThirdQuartile" [PROTECTED];
    // Two variables: `covariance`, the population one as giac's, and
    // `correlation`.
    POPULATION_COVARIANCE = "PopulationCovariance" [PROTECTED];
    CORRELATION = "Correlation" [PROTECTED];
    // Simplification (CAS-DOC-14).
    SIMPLIFY = "Simplify" [PROTECTED];
    TRIG_EXPAND = "TrigExpand" [PROTECTED];
    TRIG_REDUCE = "TrigReduce" [PROTECTED];
    TRIG_TO_EXP = "TrigToExp" [PROTECTED];
    EXP_TO_TRIG = "ExpToTrig" [PROTECTED];
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
    // The heads atoms have for typed blanks, as in Mathematica (CAS-150):
    // `Blank[Integer]` matches `2`. With `Complex`, for complex numbers and
    // complex balls; `Real` is that of real balls.
    INTEGER = "Integer" [PROTECTED];
    RATIONAL = "Rational" [PROTECTED];
    REAL = "Real" [PROTECTED];
    STRING = "String" [PROTECTED];
    SYMBOL = "Symbol" [PROTECTED];
    // A ball `Ball[c, r]` (CAS-DOC-5 §6.5), also the head of a ball in full
    // form. It threads over lists itself, a scalar radius spread over them.
    BALL = "Ball" [PROTECTED];
    // A decimal literal of the language, `Decimal[digits, exponent]` for
    // `digits · 10^exponent` with its digits as written (CAS-148): the ball
    // of it at the working precision. The session makes the balls of a
    // statement's literals when it runs the statement.
    DECIMAL = "Decimal" [PROTECTED];
    // `Exact[e]`: each ball of `e` as the simplest fraction that matches the
    // digits it shows (CAS-DOC-5 §6.4, D4).
    EXACT = "Exact" [PROTECTED];
    // `N[e, d]`: `e` with balls for its numbers, to `d` digits or those of
    // the working precision, raised until they are reached (CAS-DOC-5 §9.2,
    // CAS-84).
    N = "N" [PROTECTED];
    // The decimal digits the working precision shows, Xcas's `Digits`: its
    // value is read from the context; the session sets it (CAS-151).
    DIGITS = "Digits" [PROTECTED];
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
    // Predicates that always decide (CAS-DOC-4 §6.1).
    SAME_Q = "SameQ" [PROTECTED];
    FREE_Q = "FreeQ" [PROTECTED];
    // Signs, which assumptions decide (CAS-DOC-4 §7).
    POSITIVE = "Positive" [LISTABLE, PROTECTED];
    NEGATIVE = "Negative" [LISTABLE, PROTECTED];
    NON_NEGATIVE = "NonNegative" [LISTABLE, PROTECTED];
    NON_POSITIVE = "NonPositive" [LISTABLE, PROTECTED];
    ABS = "Abs" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Integer division, Euclidean, and parity (CAS-122).
    REMAINDER = "Remainder" [LISTABLE, PROTECTED];
    QUOTIENT = "Quotient" [LISTABLE, PROTECTED];
    EVEN_Q = "EvenQ" [LISTABLE, PROTECTED];
    ODD_Q = "OddQ" [LISTABLE, PROTECTED];
    // Polynomials, phase 1 (CAS-DOC-8 §10, §15 detail 1; CAS-140): Xcas's
    // `expand`, `normal`, `gcd`, `lcm`, `quo`, `rem`, `quorem`, `degree`,
    // `coeff`, `lcoeff`, `numer`, `denom`, `content`, `primpart`,
    // `resultant` and `discriminant`.
    EXPAND = "Expand" [LISTABLE, PROTECTED];
    TOGETHER = "Together" [LISTABLE, PROTECTED];
    GCD = "GCD" [FLAT, ORDERLESS, ONE_IDENTITY, LISTABLE, PROTECTED];
    LCM = "LCM" [FLAT, ORDERLESS, ONE_IDENTITY, LISTABLE, PROTECTED];
    POLYNOMIAL_QUOTIENT = "PolynomialQuotient" [PROTECTED];
    POLYNOMIAL_REMAINDER = "PolynomialRemainder" [PROTECTED];
    POLYNOMIAL_QUOTIENT_REMAINDER = "PolynomialQuotientRemainder" [PROTECTED];
    EXPONENT = "Exponent" [LISTABLE, PROTECTED];
    COEFFICIENT = "Coefficient" [LISTABLE, PROTECTED];
    LEADING_COEFFICIENT = "LeadingCoefficient" [LISTABLE, PROTECTED];
    NUMERATOR = "Numerator" [LISTABLE, PROTECTED];
    DENOMINATOR = "Denominator" [LISTABLE, PROTECTED];
    CONTENT = "Content" [LISTABLE, PROTECTED];
    PRIMITIVE_PART = "PrimitivePart" [LISTABLE, PROTECTED];
    RESULTANT = "Resultant" [PROTECTED];
    DISCRIMINANT = "Discriminant" [LISTABLE, PROTECTED];
    // Factoring (CAS-DOC-8 §4, CAS-185, CAS-187): Xcas's `factor`,
    // `collect`, `factors` and `sqrfree`.
    FACTOR = "Factor" [LISTABLE, PROTECTED];
    FACTOR_RATIONAL = "FactorRational" [LISTABLE, PROTECTED];
    FACTOR_LIST = "FactorList" [LISTABLE, PROTECTED];
    FACTOR_SQUARE_FREE = "FactorSquareFree" [LISTABLE, PROTECTED];
    // Gröbner bases (CAS-DOC-21): Xcas's `gbasis`, `greduce` and
    // `eliminate`, and the monomial orders, Mathematica's `MonomialOrder`
    // values: `plex`, `tdeg` and `revlex`.
    GROEBNER_BASIS = "GroebnerBasis" [PROTECTED];
    GROEBNER_REDUCE = "GroebnerReduce" [PROTECTED];
    ELIMINATE = "Eliminate" [PROTECTED];
    LEXICOGRAPHIC = "Lexicographic" [PROTECTED];
    DEGREE_LEXICOGRAPHIC = "DegreeLexicographic" [PROTECTED];
    DEGREE_REVERSE_LEXICOGRAPHIC = "DegreeReverseLexicographic" [PROTECTED];
    // The number sets of `Element` (CAS-DOC-4 §7).
    INTEGERS = "Integers" [PROTECTED];
    RATIONALS = "Rationals" [PROTECTED];
    REALS = "Reals" [PROTECTED];
    COMPLEXES = "Complexes" [PROTECTED];
    // Structure.
    LIST = "List" [PROTECTED];
    SEQUENCE = "Sequence" [PROTECTED];
    HOLD = "Hold" [HOLD_ALL, PROTECTED];
    RULE = "Rule" [PROTECTED];
    RULE_DELAYED = "RuleDelayed" [HOLD_REST, PROTECTED];
    FUNCTION = "Function" [HOLD_ALL, PROTECTED];
    // A function whose body is evaluated once, as it is made (CAS-DOC-7 §4).
    UNAPPLY = "Unapply" [PROTECTED];
    RANGE = "Range" [PROTECTED];
    FINITE_SET = "FiniteSet" [ORDERLESS, PROTECTED];
    LENGTH = "Length" [PROTECTED];
    PART = "Part" [PROTECTED];
    MAP = "Map" [PROTECTED];
    APPLY = "Apply" [PROTECTED];
    // Indexing from 0, and a copy with an element replaced (CAS-DOC-7 §7).
    AT = "At" [PROTECTED];
    REPLACE_AT = "ReplaceAt" [PROTECTED];
    REPLACE_ALL = "ReplaceAll" [PROTECTED];
    REPLACE = "Replace" [PROTECTED];
    REPLACE_REPEATED = "ReplaceRepeated" [PROTECTED];
    MATCH_Q = "MatchQ" [PROTECTED];
    // `#n` of a pure function without named parameters (CAS-55).
    SLOT = "Slot" [PROTECTED];
    // Procedures (CAS-DOC-7 §8): a block with its locals, run statement by
    // statement in a frame of their values; a group of statements in the
    // enclosing block (a branch, the body of a loop); loops and jumps.
    BLOCK = "Block" [HOLD_ALL, PROTECTED];
    STATEMENTS = "Statements" [HOLD_ALL, PROTECTED];
    FOR = "For" [HOLD_ALL, PROTECTED];
    WHILE = "While" [HOLD_ALL, PROTECTED];
    RETURN = "Return" [HOLD_ALL, PROTECTED];
    BREAK = "Break" [PROTECTED];
    CONTINUE = "Continue" [PROTECTED];
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
            | Symbol::HOLD_PATTERN
    )
}
