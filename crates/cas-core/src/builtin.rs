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
    // The lower incomplete gamma function (CAS-DOC-46 §3); the upper is
    // `Gamma` with two arguments, the incomplete beta `Beta` with three.
    LOWER_GAMMA = "LowerGamma" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ZETA = "Zeta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ERF = "Erf" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    ERFC = "Erfc" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // The exponential, sine, cosine and logarithmic integrals (CAS-DOC-31
    // §3.1), Xcas's `Ei`, `Si`, `Ci` and `Li`; `LogIntegral` is rewritten
    // as `ExpIntegralEi` of the logarithm on evaluation.
    EXP_INTEGRAL_EI = "ExpIntegralEi" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SIN_INTEGRAL = "SinIntegral" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    COS_INTEGRAL = "CosIntegral" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    LOG_INTEGRAL = "LogIntegral" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // The digamma function and its derivatives (CAS-DOC-31 §3.2), Xcas's
    // `Psi(x)` and `Psi(x, n)`, in giac's argument order: not listable, as
    // the orthogonal polynomials, for `Psi([1,2], 1)`.
    POLY_GAMMA = "PolyGamma" [NUMERIC_FUNCTION, PROTECTED];
    // Lambert's W on its branches 0 and −1 (CAS-DOC-31 §3.3), Xcas's
    // `LambertW(x)` and `LambertW(x, k)`; not listable, as `PolyGamma`.
    LAMBERT_W = "LambertW" [NUMERIC_FUNCTION, PROTECTED];
    // The Bessel and Airy functions (CAS-DOC-34 §3), the order first, as
    // giac's.
    BESSEL_J = "BesselJ" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    BESSEL_Y = "BesselY" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    BESSEL_I = "BesselI" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    BESSEL_K = "BesselK" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    AIRY_AI = "AiryAi" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    AIRY_BI = "AiryBi" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    AIRY_AI_PRIME = "AiryAiPrime" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    AIRY_BI_PRIME = "AiryBiPrime" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // The step and the delta (CAS-DOC-18 §2).
    HEAVISIDE_THETA = "HeavisideTheta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    DIRAC_DELTA = "DiracDelta" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // The box, giac's `rect` (CAS-DOC-41 §1).
    UNIT_BOX = "UnitBox" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SQRT = "Sqrt" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    SURD = "Surd" [PROTECTED];
    FACTORIAL = "Factorial" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Calculus (CAS-DOC-4 §8): the derivative, as system rules, and the
    // derivative of a function, `Derivative[n][f]`; integrals (CAS-DOC-9).
    D = "D" [PROTECTED];
    DERIVATIVE = "Derivative" [PROTECTED];
    INTEGRATE = "Integrate" [PROTECTED];
    // Risch integration of elementary functions (CAS-585), giac's `risch`.
    RISCH = "Risch" [PROTECTED];
    // Series and limits (CAS-DOC-10): `OrderSize[u]`, giac's `order_size`,
    // a bounded function of `u`, inert.
    SERIES = "Series" [PROTECTED];
    LIMIT = "Limit" [PROTECTED];
    ORDER_SIZE = "OrderSize" [PROTECTED];
    // Differential equations (CAS-DOC-16).
    DSOLVE = "DSolve" [PROTECTED];
    // Liouvillian solutions by Kovacic's algorithm (CAS-584), Xcas's.
    KOVACICSOLS = "KovacicSols" [PROTECTED];
    // Numerical solutions of initial value problems (CAS-DOC-35).
    NDSOLVE = "NDSolve" [PROTECTED];
    BVP_SOLVE = "BVPSolve" [PROTECTED];
    // Signal processing (CAS-521): giac's window functions.
    BARTLETT_HANN_WINDOW = "BartlettHannWindow" [PROTECTED];
    BLACKMAN_HARRIS_WINDOW = "BlackmanHarrisWindow" [PROTECTED];
    BLACKMAN_WINDOW = "BlackmanWindow" [PROTECTED];
    BOHMAN_WINDOW = "BohmanWindow" [PROTECTED];
    COSINE_WINDOW = "CosineWindow" [PROTECTED];
    GAUSSIAN_WINDOW = "GaussianWindow" [PROTECTED];
    HAMMING_WINDOW = "HammingWindow" [PROTECTED];
    HANN_POISSON_WINDOW = "HannPoissonWindow" [PROTECTED];
    HANN_WINDOW = "HannWindow" [PROTECTED];
    PARZEN_WINDOW = "ParzenWindow" [PROTECTED];
    POISSON_WINDOW = "PoissonWindow" [PROTECTED];
    RIEMANN_WINDOW = "RiemannWindow" [PROTECTED];
    TRIANGLE_WINDOW = "TriangleWindow" [PROTECTED];
    TUKEY_WINDOW = "TukeyWindow" [PROTECTED];
    WELCH_WINDOW = "WelchWindow" [PROTECTED];
    // First-order RC filters (CAS-522).
    LOW_PASS = "LowPass" [PROTECTED];
    HIGH_PASS = "HighPass" [PROTECTED];
    // Correlations of signals (CAS-523).
    CROSS_CORRELATION = "CrossCorrelation" [PROTECTED];
    AUTO_CORRELATION = "AutoCorrelation" [PROTECTED];
    // Thresholds of lists (CAS-524).
    THRESHOLD = "Threshold" [PROTECTED];
    // Strings, Xcas's (CAS-611, CAS-613).
    CHAR = "Char" [PROTECTED];
    ASC = "Asc" [PROTECTED];
    ORD = "Ord" [PROTECTED];
    IN_STRING = "InString" [PROTECTED];
    CONTAINS = "Contains" [PROTECTED];
    MID = "Mid" [PROTECTED];
    LEFT = "Left" [PROTECTED];
    RIGHT = "Right" [PROTECTED];
    LEVENSHTEIN = "Levenshtein" [PROTECTED];
    HAMMING = "Hamming" [PROTECTED];
    // Strings, Xcas's (CAS-612).
    SPLIT = "Split" [PROTECTED];
    STRING_JOIN = "StringJoin" [PROTECTED];
    TRIM = "Trim" [PROTECTED];
    STRIP = "Strip" [PROTECTED];
    // Expressions as text and back, Xcas's (CAS-614).
    TO_STRING = "ToString" [PROTECTED];
    TO_EXPRESSION = "ToExpression" [PROTECTED];
    CAT = "Cat" [PROTECTED];
    // LaTeX, Xcas's (CAS-615).
    LA_TE_X = "LaTeX" [PROTECTED];
    // MathML, Xcas's (CAS-616).
    EXPORT_MATH_ML = "ExportMathML" [PROTECTED];
    // Laplace transforms (CAS-DOC-18).
    LAPLACE_TRANSFORM = "LaplaceTransform" [PROTECTED];
    INVERSE_LAPLACE_TRANSFORM = "InverseLaplaceTransform" [PROTECTED];
    // Modular arithmetic (CAS-DOC-43): giac's `a % n`.
    MODULAR = "Modular" [PROTECTED];
    // Integer matrices (CAS-DOC-42).
    HERMITE_DECOMPOSITION = "HermiteDecomposition" [PROTECTED];
    SMITH_DECOMPOSITION = "SmithDecomposition" [PROTECTED];
    LATTICE_REDUCE = "LatticeReduce" [PROTECTED];
    // Fourier transforms (CAS-DOC-41).
    FOURIER_TRANSFORM = "FourierTransform" [PROTECTED];
    INVERSE_FOURIER_TRANSFORM = "InverseFourierTransform" [PROTECTED];
    // Z-transforms and recurrences (CAS-DOC-24): Xcas's `ztrans` and
    // `Kronecker`.
    ZTRANSFORM = "ZTransform" [PROTECTED];
    INVERSE_ZTRANSFORM = "InverseZTransform" [PROTECTED];
    KRONECKER_DELTA = "KroneckerDelta" [LISTABLE, PROTECTED];
    RSOLVE = "RSolve" [PROTECTED];
    SEQ_SOLVE = "SeqSolve" [PROTECTED];
    // Residues and Fourier coefficients (CAS-DOC-25).
    RESIDUE = "Residue" [PROTECTED];
    FOURIER_COS_COEFFICIENT = "FourierCosCoefficient" [PROTECTED];
    FOURIER_SIN_COEFFICIENT = "FourierSinCoefficient" [PROTECTED];
    FOURIER_COEFFICIENT = "FourierCoefficient" [PROTECTED];
    // Certified numerical integration (CAS-DOC-26).
    NINTEGRATE = "NIntegrate" [PROTECTED];
    // Areas by numeric rules and their options (CAS-588, CAS-589).
    AREA = "Area" [PROTECTED];
    LEFT_RECTANGLE = "LeftRectangle" [PROTECTED];
    RIGHT_RECTANGLE = "RightRectangle" [PROTECTED];
    MIDDLE_POINT = "MiddlePoint" [PROTECTED];
    TRAPEZOID = "Trapezoid" [PROTECTED];
    SIMPSON = "Simpson" [PROTECTED];
    ROMBERG_TRAPEZOID = "RombergTrapezoid" [PROTECTED];
    ROMBERG_MIDPOINT = "RombergMidpoint" [PROTECTED];
    GAUSS15 = "Gauss15" [PROTECTED];
    // Interpolation and orthogonal polynomials (CAS-DOC-27).
    INTERPOLATING_POLYNOMIAL = "InterpolatingPolynomial" [PROTECTED];
    SPLINE = "Spline" [PROTECTED];
    THIELE_INTERPOLATION = "ThieleInterpolation" [PROTECTED];
    // Sampling and intervals (CAS-593).
    LIN_SPACE = "LinSpace" [PROTECTED];
    LINEAR_INTERPOLATE = "LinearInterpolate" [PROTECTED];
    INTERVAL_TO_CENTER = "IntervalToCenter" [PROTECTED];
    CENTER_TO_INTERVAL = "CenterToInterval" [PROTECTED];
    // Rational and trigonometric interpolation (CAS-594).
    RATIONAL_INTERPOLATION = "RationalInterpolation" [PROTECTED];
    TRIGONOMETRIC_INTERPOLATION = "TrigonometricInterpolation" [PROTECTED];
    // Numerical derivatives (CAS-595).
    NUMERICAL_DERIVATIVE = "NumericalDerivative" [PROTECTED];
    N_DERIV = "NDeriv" [PROTECTED];
    // B-splines and spline fits (CAS-596).
    B_SPLINE = "BSpline" [PROTECTED];
    SPLINE_FIT = "SplineFit" [PROTECTED];
    // Polynomial fits (CAS-597).
    POLYNOMIAL_FIT = "PolynomialFit" [PROTECTED];
    // Minimax approximation (CAS-599).
    MINIMAX_APPROXIMATION = "MinimaxApproximation" [PROTECTED];
    LEGENDRE_P = "LegendreP" [PROTECTED];
    HERMITE_H = "HermiteH" [PROTECTED];
    LAGUERRE_L = "LaguerreL" [PROTECTED];
    CHEBYSHEV_T = "ChebyshevT" [PROTECTED];
    CHEBYSHEV_U = "ChebyshevU" [PROTECTED];
    // List operations (CAS-DOC-28). `Table` holds its body, evaluated at
    // each value of the index (CAS-DOC-33 §2).
    TABLE = "Table" [HOLD_ALL, PROTECTED];
    MAKE_LIST = "MakeList" [PROTECTED];
    SORT = "Sort" [PROTECTED];
    REVERSE = "Reverse" [PROTECTED];
    ROTATE = "Rotate" [PROTECTED];
    JOIN = "Join" [PROTECTED];
    APPEND = "Append" [PROTECTED];
    PREPEND = "Prepend" [PROTECTED];
    FIRST = "First" [PROTECTED];
    REST = "Rest" [PROTECTED];
    FLATTEN = "Flatten" [PROTECTED];
    DIMENSIONS = "Dimensions" [PROTECTED];
    SELECT = "Select" [PROTECTED];
    REJECT = "Reject" [PROTECTED];
    COUNT = "Count" [PROTECTED];
    MEMBER_INDEX = "MemberIndex" [PROTECTED];
    ZIP = "Zip" [PROTECTED];
    ACCUMULATE = "Accumulate" [PROTECTED];
    // Xcas's list commands (CAS-555 to CAS-557).
    FIND_POSITIONS = "FindPositions" [PROTECTED];
    FIRST_POSITION = "FirstPosition" [PROTECTED];
    COUNT_EQUAL = "CountEqual" [PROTECTED];
    COUNT_LESS = "CountLess" [PROTECTED];
    COUNT_GREATER = "CountGreater" [PROTECTED];
    DIFFERENCES = "Differences" [PROTECTED];
    ACCUMULATE_HEAD_TAIL = "AccumulateHeadTail" [PROTECTED];
    LENGTHS = "Lengths" [PROTECTED];
    MAKE_VECTOR = "MakeVector" [PROTECTED];
    MAKE_SEQUENCE = "MakeSequence" [PROTECTED];
    NEW_LIST = "NewList" [PROTECTED];
    ENUMERATE = "Enumerate" [PROTECTED];
    COPY = "Copy" [PROTECTED];
    // Xcas's list commands (CAS-558, CAS-559).
    APPLY_EACH = "ApplyEach" [PROTECTED];
    FOLD_LEFT = "FoldLeft" [PROTECTED];
    FOLD_RIGHT = "FoldRight" [PROTECTED];
    MULTIPLY = "Multiply" [PROTECTED];
    INSERT_AT = "InsertAt" [PROTECTED];
    SUBSTITUTE_OPERANDS = "SubstituteOperands" [PROTECTED];
    // Xcas's heaps (CAS-560).
    HEAPIFY = "Heapify" [PROTECTED];
    HEAP_PUSH = "HeapPush" [PROTECTED];
    HEAP_POP = "HeapPop" [PROTECTED];
    // Xcas's expression commands (CAS-562 to CAS-567).
    POWER_EXPAND = "PowerExpand" [PROTECTED];
    COMMON_DENOMINATOR = "CommonDenominator" [PROTECTED];
    RAW_NUMERATOR = "RawNumerator" [LISTABLE, PROTECTED];
    RAW_DENOMINATOR = "RawDenominator" [LISTABLE, PROTECTED];
    CANCEL_PAIR = "CancelPair" [PROTECTED];
    REGROUP = "Regroup" [PROTECTED];
    SUM_OF_SQUARES = "SumOfSquares" [PROTECTED];
    LOG_B = "LogB" [LISTABLE, PROTECTED];
    LOG2 = "Log2" [LISTABLE, PROTECTED];
    LOG10 = "Log10" [LISTABLE, PROTECTED];
    EXP10 = "Exp10" [LISTABLE, PROTECTED];
    RADIANS = "Radians" [LISTABLE, PROTECTED];
    DEGREES = "Degrees" [LISTABLE, PROTECTED];
    NORMALIZE = "Normalize" [PROTECTED];
    ALGEBRAIC_SUBSTITUTE = "AlgebraicSubstitute" [PROTECTED];
    // Xcas's expression commands, second part (CAS-568 to CAS-571).
    IDENTIFIERS = "Identifiers" [PROTECTED];
    KERNELS = "Kernels" [PROTECTED];
    ALGEBRAIC_VARIABLES = "AlgebraicVariables" [PROTECTED];
    REORDER = "Reorder" [PROTECTED];
    EVALUATE = "Evaluate" [PROTECTED];
    EVALUATE_BOOLEAN = "EvaluateBoolean" [PROTECTED];
    EPSILON_TO_ZERO = "EpsilonToZero" [HOLD_ALL, PROTECTED];
    INFINITE_Q = "InfiniteQ" [PROTECTED];
    INDETERMINATE_Q = "IndeterminateQ" [PROTECTED];
    COMPARE = "Compare" [PROTECTED];
    EQUATION_TO_DIFFERENCE = "EquationToDifference" [PROTECTED];
    EQUATION_TO_LIST = "EquationToList" [PROTECTED];
    LIST_TO_EQUATIONS = "ListToEquations" [PROTECTED];
    EQUATIONS_TO_LIST = "EquationsToList" [PROTECTED];
    // Polynomial tools (CAS-DOC-29).
    APART = "Apart" [PROTECTED];
    COMPLEX_APART = "ComplexApart" [PROTECTED];
    PROPER_FRACTION = "ProperFraction" [PROTECTED];
    NUMERATOR_DENOMINATOR = "NumeratorDenominator" [PROTECTED];
    COEFFICIENT_LIST = "CoefficientList" [PROTECTED];
    FROM_COEFFICIENTS = "FromCoefficients" [PROTECTED];
    TRAILING_COEFFICIENT = "TrailingCoefficient" [PROTECTED];
    VALUATION = "Valuation" [PROTECTED];
    COUNT_ROOTS = "CountRoots" [PROTECTED];
    REAL_ROOTS = "RealRoots" [PROTECTED];
    // Roots (CAS-538, CAS-539).
    REAL_ROOT_ISOLATION = "RealRootIsolation" [PROTECTED];
    POSITIVE_ROOT_ISOLATION = "PositiveRootIsolation" [PROTECTED];
    POSITIVE_ROOT_UPPER_BOUND = "PositiveRootUpperBound" [PROTECTED];
    POSITIVE_ROOT_LOWER_BOUND = "PositiveRootLowerBound" [PROTECTED];
    STURM_SEQUENCE = "SturmSequence" [PROTECTED];
    RATIONAL_ROOTS = "RationalRoots" [PROTECTED];
    COMPLEX_RATIONAL_ROOTS = "ComplexRationalRoots" [PROTECTED];
    COMPLEX_ROOTS = "ComplexRoots" [PROTECTED];
    // Polynomial tools II (CAS-DOC-48).
    POLYNOMIAL_ROOTS = "PolynomialRoots" [PROTECTED];
    FROM_ROOTS = "FromRoots" [PROTECTED];
    HORNER_EVAL = "HornerEval" [PROTECTED];
    TAYLOR_SHIFT = "TaylorShift" [PROTECTED];
    ROOT_MULTIPLICITIES = "RootMultiplicities" [PROTECTED];
    POLYNOMIAL_EXTENDED_GCD = "PolynomialExtendedGCD" [PROTECTED];
    POLYNOMIAL_BEZOUT = "PolynomialBezout" [PROTECTED];
    POLYNOMIAL_CHINESE_REMAINDER = "PolynomialChineseRemainder" [PROTECTED];
    SERIES_DIVISION = "SeriesDivision" [PROTECTED];
    PADE_APPROXIMANT = "PadeApproximant" [PROTECTED];
    COMPLEX_FACTOR = "ComplexFactor" [PROTECTED];
    COMPLETE_SQUARE = "CompleteSquare" [PROTECTED];
    DIGIT_POLYNOMIAL = "DigitPolynomial" [PROTECTED];
    BERLEKAMP_MASSEY = "BerlekampMassey" [PROTECTED];
    // Polynomial tools III (CAS-540, CAS-543).
    ROOT_OF = "RootOf" [PROTECTED];
    FACTORED_FROM_ROOTS = "FactoredFromRoots" [PROTECTED];
    HOMOGENIZE = "Homogenize" [PROTECTED];
    TOTAL_DEGREE = "TotalDegree" [PROTECTED];
    FACTOR_XN = "FactorXn" [PROTECTED];
    BEZOUT_MATRIX = "BezoutMatrix" [PROTECTED];
    IRREDUCIBLE_POLYNOMIAL_Q = "IrreduciblePolynomialQ" [PROTECTED];
    // Matrix decompositions and functions (CAS-DOC-30).
    LU_DECOMPOSITION = "LUDecomposition" [PROTECTED];
    QR_DECOMPOSITION = "QRDecomposition" [PROTECTED];
    CHOLESKY_DECOMPOSITION = "CholeskyDecomposition" [PROTECTED];
    JORDAN_DECOMPOSITION = "JordanDecomposition" [PROTECTED];
    MINIMAL_POLYNOMIAL = "MinimalPolynomial" [PROTECTED];
    MATRIX_EXP = "MatrixExp" [PROTECTED];
    MATRIX_POWER = "MatrixPower" [PROTECTED];
    DIAGONAL_MATRIX = "DiagonalMatrix" [PROTECTED];
    HILBERT_MATRIX = "HilbertMatrix" [PROTECTED];
    VANDERMONDE_MATRIX = "VandermondeMatrix" [PROTECTED];
    NORM = "Norm" [PROTECTED];
    FROBENIUS_NORM = "FrobeniusNorm" [PROTECTED];
    // Matrix decompositions II (CAS-DOC-47).
    HESSENBERG_DECOMPOSITION = "HessenbergDecomposition" [PROTECTED];
    COMPANION_MATRIX = "CompanionMatrix" [PROTECTED];
    SYLVESTER_MATRIX = "SylvesterMatrix" [PROTECTED];
    HADAMARD_PRODUCT = "HadamardProduct" [PROTECTED];
    LDL_DECOMPOSITION = "LDLDecomposition" [PROTECTED];
    POSITIVE_DEFINITE_Q = "PositiveDefiniteQ" [PROTECTED];
    INERTIA = "Inertia" [PROTECTED];
    MATRIX_NORM = "MatrixNorm" [PROTECTED];
    L1_NORM = "L1Norm" [PROTECTED];
    L_INFINITY_NORM = "LInfinityNorm" [PROTECTED];
    MAX_NORM = "MaxNorm" [PROTECTED];
    CONDITION_NUMBER = "ConditionNumber" [PROTECTED];
    GRAM_SCHMIDT = "GramSchmidt" [PROTECTED];
    SINGULAR_VALUES = "SingularValues" [PROTECTED];
    CONJUGATE_GRADIENT = "ConjugateGradient" [PROTECTED];
    // Matrix decompositions III (CAS-DOC-51).
    SCHUR_DECOMPOSITION = "SchurDecomposition" [PROTECTED];
    SINGULAR_VALUE_DECOMPOSITION = "SingularValueDecomposition" [PROTECTED];
    JACOBI_LINEAR_SOLVE = "JacobiLinearSolve" [PROTECTED];
    GAUSS_SEIDEL_LINEAR_SOLVE = "GaussSeidelLinearSolve" [PROTECTED];
    // Equations (CAS-DOC-11): the real solutions, the complex ones, and
    // numeric roots.
    SOLVE = "Solve" [PROTECTED];
    CSOLVE = "CSolve" [PROTECTED];
    FSOLVE = "FSolve" [PROTECTED];
    CFSOLVE = "CFSolve" [PROTECTED];
    // giac's damped Newton iteration (CAS-591).
    NEWTON = "Newton" [PROTECTED];
    NSOLVE = "NSolve" [PROTECTED];
    // FSolve's methods, options (CAS-590).
    BISECTION_SOLVER = "BisectionSolver" [PROTECTED];
    FALSE_POSITION_SOLVER = "FalsePositionSolver" [PROTECTED];
    BRENT_SOLVER = "BrentSolver" [PROTECTED];
    NEWTON_SOLVER = "NewtonSolver" [PROTECTED];
    SECANT_SOLVER = "SecantSolver" [PROTECTED];
    STEFFENSEN_SOLVER = "SteffensenSolver" [PROTECTED];
    DNEWTON_SOLVER = "DNewtonSolver" [PROTECTED];
    HYBRID_SOLVER = "HybridSolver" [PROTECTED];
    HYBRIDS_SOLVER = "HybridsSolver" [PROTECTED];
    HYBRIDJ_SOLVER = "HybridjSolver" [PROTECTED];
    HYBRIDSJ_SOLVER = "HybridsjSolver" [PROTECTED];
    NEWTONJ_SOLVER = "NewtonjSolver" [PROTECTED];
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
    // giac's matrix constructors, slices and row operations (CAS-531,
    // CAS-532, CAS-533).
    MATRIX = "Matrix" [PROTECTED];
    MAKE_MATRIX = "MakeMatrix" [PROTECTED];
    LIST_TO_MATRIX = "ListToMatrix" [PROTECTED];
    MATRIX_TO_LIST = "MatrixToList" [PROTECTED];
    BLOCK_MATRIX = "BlockMatrix" [PROTECTED];
    JORDAN_BLOCK = "JordanBlock" [PROTECTED];
    BLOCK_DIAGONAL_MATRIX = "BlockDiagonalMatrix" [PROTECTED];
    ROW_DIMENSION = "RowDimension" [PROTECTED];
    COLUMN_DIMENSION = "ColumnDimension" [PROTECTED];
    MATRIX_ROW = "MatrixRow" [PROTECTED];
    MATRIX_COLUMN = "MatrixColumn" [PROTECTED];
    SUB_MATRIX = "SubMatrix" [PROTECTED];
    DELETE_ROWS = "DeleteRows" [PROTECTED];
    DELETE_COLUMNS = "DeleteColumns" [PROTECTED];
    MULTIPLY_ROW = "MultiplyRow" [PROTECTED];
    MULTIPLY_ROW_ADD = "MultiplyRowAdd" [PROTECTED];
    SCALE_ROW_ADD = "ScaleRowAdd" [PROTECTED];
    SCALE_ROW = "ScaleRow" [PROTECTED];
    ROW_ADD = "RowAdd" [PROTECTED];
    ROW_SWAP = "RowSwap" [PROTECTED];
    COLUMN_SWAP = "ColumnSwap" [PROTECTED];
    ROW_SPACE = "RowSpace" [PROTECTED];
    COL_SPACE = "ColSpace" [PROTECTED];
    INTERSECTION_BASIS = "IntersectionBasis" [PROTECTED];
    SYSTEM_MATRIX = "SystemMatrix" [PROTECTED];
    // giac's minors, echelon forms and pivots (CAS-530, CAS-528).
    DET_MINOR = "DetMinor" [PROTECTED];
    ADJOINT_MATRIX = "AdjointMatrix" [PROTECTED];
    ROW_ECHELON = "RowEchelon" [PROTECTED];
    PIVOT = "Pivot" [PROTECTED];
    JORDAN_MATRIX = "JordanMatrix" [PROTECTED];
    KEEP_PIVOT = "KeepPivot" [PROTECTED];
    // giac's sums of squares and rational Jordan forms (CAS-535).
    GAUSS_SQUARES = "GaussSquares" [PROTECTED];
    RATIONAL_JORDAN_FORM = "RationalJordanForm" [PROTECTED];
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
    ECM_FACTOR = "EcmFactor" [LISTABLE, PROTECTED];
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
    // Xcas's integer functions (CAS-DOC-55, CAS-548 to CAS-553).
    CHINESE_REMAINDER_LISTS = "ChineseRemainderLists" [PROTECTED];
    LIST_GCD = "ListGCD" [PROTECTED];
    FRACTION_MOD = "FractionMod" [PROTECTED];
    INTEGER_COMPOSITIONS = "IntegerCompositions" [PROTECTED];
    BERNOULLI_MOD = "BernoulliMod" [PROTECTED];
    IS_PRIME = "IsPrime" [PROTECTED];
    // Permutations (CAS-DOC-45 §2): giac's lists of 0, …, n − 1.
    PERMUTATION_Q = "PermutationQ" [PROTECTED];
    PERMUTATION_CYCLES = "PermutationCycles" [PROTECTED];
    CYCLES_PERMUTATION = "CyclesPermutation" [PROTECTED];
    CYCLE_PERMUTATION = "CyclePermutation" [PROTECTED];
    INVERSE_PERMUTATION = "InversePermutation" [PROTECTED];
    INVERSE_CYCLE = "InverseCycle" [PROTECTED];
    PERMUTATION_PRODUCT = "PermutationProduct" [PROTECTED];
    PERMUTATION_CYCLE_PRODUCT = "PermutationCycleProduct" [PROTECTED];
    CYCLE_PERMUTATION_PRODUCT = "CyclePermutationProduct" [PROTECTED];
    CYCLE_PRODUCT = "CycleProduct" [PROTECTED];
    PERMUTATION_ORDER = "PermutationOrder" [PROTECTED];
    PERMUTATION_SIGNATURE = "PermutationSignature" [PROTECTED];
    PERMUTATION_MATRIX = "PermutationMatrix" [PROTECTED];
    // giac's nextperm, prevperm and groupermu (CAS-561).
    NEXT_PERMUTATION = "NextPermutation" [PROTECTED];
    PREVIOUS_PERMUTATION = "PreviousPermutation" [PROTECTED];
    PERMUTATION_GROUP = "PermutationGroup" [PROTECTED];
    // Bits, digits and falling factorials (CAS-DOC-45 §3).
    BIT_AND = "BitAnd" [LISTABLE, PROTECTED];
    BIT_OR = "BitOr" [LISTABLE, PROTECTED];
    BIT_XOR = "BitXor" [LISTABLE, PROTECTED];
    HAMMING_DISTANCE = "HammingDistance" [LISTABLE, PROTECTED];
    CONVERT = "Convert" [PROTECTED];
    // convert's option `interval` (CAS-592).
    INTERVAL = "Interval" [PROTECTED];
    FACTORIAL_POWER = "FactorialPower" [LISTABLE, PROTECTED];
    // The discrete Fourier transform (CAS-DOC-45 §4).
    DISCRETE_FOURIER_TRANSFORM = "DiscreteFourierTransform" [PROTECTED];
    INVERSE_DISCRETE_FOURIER_TRANSFORM = "InverseDiscreteFourierTransform" [PROTECTED];
    CONVOLUTION = "Convolution" [PROTECTED];
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
    // Statistics II (CAS-DOC-50).
    FREQUENCIES = "Frequencies" [PROTECTED];
    CUMULATED_FREQUENCIES = "CumulatedFrequencies" [PROTECTED];
    CLASSES = "Classes" [PROTECTED];
    COVARIANCE_CORRELATION = "CovarianceCorrelation" [PROTECTED];
    ROOT_MEAN_SQUARE = "RootMeanSquare" [PROTECTED];
    MOVING_AVERAGE = "MovingAverage" [PROTECTED];
    // Hypothesis tests (CAS-DOC-52).
    Z_TEST = "ZTest" [PROTECTED];
    T_TEST = "TTest" [PROTECTED];
    CHI_SQUARE_TEST = "ChiSquareTest" [PROTECTED];
    KOLMOGOROV_CDF = "KolmogorovCDF" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    KOLMOGOROV_SMIRNOV_TEST = "KolmogorovSmirnovTest" [PROTECTED];
    WILCOXON_DISTRIBUTION = "WilcoxonDistribution" [PROTECTED];
    WILCOXON_STATISTIC = "WilcoxonStatistic" [PROTECTED];
    WILCOXON_TEST = "WilcoxonTest" [PROTECTED];
    // Regressions: `linear_regression`, `polynomial_regression`,
    // `exponential_regression`, `logarithmic_regression` and
    // `power_regression`.
    LINEAR_REGRESSION = "LinearRegression" [PROTECTED];
    POLYNOMIAL_REGRESSION = "PolynomialRegression" [PROTECTED];
    EXPONENTIAL_REGRESSION = "ExponentialRegression" [PROTECTED];
    LOGARITHMIC_REGRESSION = "LogarithmicRegression" [PROTECTED];
    POWER_REGRESSION = "PowerRegression" [PROTECTED];
    // Discrete distributions (CAS-DOC-23 §2): the binomial's mass is
    // `Binomial`'s three-argument form, Xcas's `binomial(n, p, k)`.
    BINOMIAL_CDF = "BinomialCDF" [PROTECTED];
    BINOMIAL_ICDF = "BinomialICDF" [PROTECTED];
    POISSON_PDF = "PoissonPDF" [PROTECTED];
    POISSON_CDF = "PoissonCDF" [PROTECTED];
    POISSON_ICDF = "PoissonICDF" [PROTECTED];
    GEOMETRIC_PDF = "GeometricPDF" [PROTECTED];
    GEOMETRIC_CDF = "GeometricCDF" [PROTECTED];
    GEOMETRIC_ICDF = "GeometricICDF" [PROTECTED];
    NEGATIVE_BINOMIAL_PDF = "NegativeBinomialPDF" [PROTECTED];
    NEGATIVE_BINOMIAL_CDF = "NegativeBinomialCDF" [PROTECTED];
    NEGATIVE_BINOMIAL_ICDF = "NegativeBinomialICDF" [PROTECTED];
    // The inverse of `Erf`, for the normal distribution's inverse (CAS-DOC-23
    // §3); giac has none.
    INVERSE_ERF = "InverseErf" [LISTABLE, NUMERIC_FUNCTION, PROTECTED];
    // Continuous distributions in closed form (CAS-DOC-23 §3).
    NORMAL_PDF = "NormalPDF" [PROTECTED];
    NORMAL_CDF = "NormalCDF" [PROTECTED];
    NORMAL_ICDF = "NormalICDF" [PROTECTED];
    EXPONENTIAL_PDF = "ExponentialPDF" [PROTECTED];
    EXPONENTIAL_CDF = "ExponentialCDF" [PROTECTED];
    EXPONENTIAL_ICDF = "ExponentialICDF" [PROTECTED];
    UNIFORM_PDF = "UniformPDF" [PROTECTED];
    UNIFORM_CDF = "UniformCDF" [PROTECTED];
    UNIFORM_ICDF = "UniformICDF" [PROTECTED];
    // Distributions with degrees of freedom (CAS-DOC-23 §4).
    STUDENT_T_PDF = "StudentTPDF" [PROTECTED];
    STUDENT_T_CDF = "StudentTCDF" [PROTECTED];
    STUDENT_T_ICDF = "StudentTICDF" [PROTECTED];
    CHI_SQUARE_PDF = "ChiSquarePDF" [PROTECTED];
    CHI_SQUARE_CDF = "ChiSquareCDF" [PROTECTED];
    CHI_SQUARE_ICDF = "ChiSquareICDF" [PROTECTED];
    FISHER_PDF = "FisherPDF" [PROTECTED];
    FISHER_CDF = "FisherCDF" [PROTECTED];
    FISHER_ICDF = "FisherICDF" [PROTECTED];
    // Random numbers (CAS-DOC-33): the generator's seed, and uniform
    // numbers.
    SEED_RANDOM = "SeedRandom" [PROTECTED];
    RANDOM = "Random" [PROTECTED];
    RANDOM_INTEGER = "RandomInteger" [PROTECTED];
    // Random variates of the distributions (CAS-DOC-33 §4).
    RANDOM_NORMAL = "RandomNormal" [PROTECTED];
    RANDOM_EXPONENTIAL = "RandomExponential" [PROTECTED];
    RANDOM_BINOMIAL = "RandomBinomial" [PROTECTED];
    RANDOM_POISSON = "RandomPoisson" [PROTECTED];
    RANDOM_GEOMETRIC = "RandomGeometric" [PROTECTED];
    RANDOM_CHI_SQUARE = "RandomChiSquare" [PROTECTED];
    RANDOM_STUDENT_T = "RandomStudentT" [PROTECTED];
    RANDOM_FISHER = "RandomFisher" [PROTECTED];
    // The gamma, beta, Weibull and Cauchy distributions (CAS-DOC-46 §4).
    GAMMA_PDF = "GammaPDF" [PROTECTED];
    GAMMA_CDF = "GammaCDF" [PROTECTED];
    GAMMA_ICDF = "GammaICDF" [PROTECTED];
    BETA_PDF = "BetaPDF" [PROTECTED];
    BETA_CDF = "BetaCDF" [PROTECTED];
    BETA_ICDF = "BetaICDF" [PROTECTED];
    WEIBULL_PDF = "WeibullPDF" [PROTECTED];
    WEIBULL_CDF = "WeibullCDF" [PROTECTED];
    WEIBULL_ICDF = "WeibullICDF" [PROTECTED];
    CAUCHY_PDF = "CauchyPDF" [PROTECTED];
    CAUCHY_CDF = "CauchyCDF" [PROTECTED];
    CAUCHY_ICDF = "CauchyICDF" [PROTECTED];
    RANDOM_GAMMA = "RandomGamma" [PROTECTED];
    RANDOM_BETA = "RandomBeta" [PROTECTED];
    RANDOM_WEIBULL = "RandomWeibull" [PROTECTED];
    // Random vectors, matrices and polynomials, permutations and samples
    // (CAS-DOC-33 §5, §6).
    RANDOM_VECTOR = "RandomVector" [PROTECTED];
    RANDOM_MATRIX = "RandomMatrix" [PROTECTED];
    RANDOM_POLYNOMIAL = "RandomPolynomial" [PROTECTED];
    RANDOM_PERMUTATION = "RandomPermutation" [PROTECTED];
    RANDOM_SAMPLE = "RandomSample" [PROTECTED];
    // Simplification (CAS-DOC-14).
    SIMPLIFY = "Simplify" [PROTECTED];
    TRIG_EXPAND = "TrigExpand" [PROTECTED];
    TRIG_REDUCE = "TrigReduce" [PROTECTED];
    TRIG_TO_EXP = "TrigToExp" [PROTECTED];
    EXP_TO_TRIG = "ExpToTrig" [PROTECTED];
    // Rewriting forms (CAS-DOC-36): exponentials, powers and logarithms
    // (§2).
    HYPERBOLIC_TO_EXP = "HyperbolicToExp" [PROTECTED];
    EXP_TO_POWER = "ExpToPower" [PROTECTED];
    POWER_TO_EXP = "PowerToExp" [PROTECTED];
    EXP_EXPAND = "ExpExpand" [PROTECTED];
    LOG_EXPAND = "LogExpand" [PROTECTED];
    LOG_COMBINE = "LogCombine" [PROTECTED];
    // Trigonometric substitutions (CAS-DOC-36 §3).
    HALF_TAN = "HalfTan" [PROTECTED];
    HALF_TAN_HYPERBOLIC_TO_EXP = "HalfTanHyperbolicToExp" [PROTECTED];
    TAN_TO_SIN_COS = "TanToSinCos" [PROTECTED];
    SIN_TO_COS_TAN = "SinToCosTan" [PROTECTED];
    COS_TO_SIN_TAN = "CosToSinTan" [PROTECTED];
    TAN_TO_SIN_COS_2 = "TanToSinCos2" [PROTECTED];
    TAN_TO_COS_SIN_2 = "TanToCosSin2" [PROTECTED];
    // Pythagorean forms (CAS-DOC-36 §4).
    TRIG_TO_COS = "TrigToCos" [PROTECTED];
    TRIG_TO_SIN = "TrigToSin" [PROTECTED];
    TRIG_TO_TAN = "TrigToTan" [PROTECTED];
    // giac's `tcollect` (CAS-DOC-36 §5).
    TRIG_COLLECT = "TrigCollect" [PROTECTED];
    // Inverse trigonometric functions (CAS-DOC-36 §6).
    ARC_COS_TO_ARC_SIN = "ArcCosToArcSin" [PROTECTED];
    ARC_SIN_TO_ARC_COS = "ArcSinToArcCos" [PROTECTED];
    ARC_COS_TO_ARC_TAN = "ArcCosToArcTan" [PROTECTED];
    ARC_SIN_TO_ARC_TAN = "ArcSinToArcTan" [PROTECTED];
    ARC_TAN_TO_ARC_COS = "ArcTanToArcCos" [PROTECTED];
    ARC_TAN_TO_ARC_SIN = "ArcTanToArcSin" [PROTECTED];
    ARC_TRIG_TO_LOG = "ArcTrigToLog" [PROTECTED];
    // Complex parts (CAS-DOC-36 §7).
    COMPLEX_EXPAND = "ComplexExpand" [PROTECTED];
    MULTIPLY_CONJUGATE = "MultiplyConjugate" [PROTECTED];
    MULTIPLY_COMPLEX_CONJUGATE = "MultiplyComplexConjugate" [PROTECTED];
    // Combinations and the Heaviside forms (CAS-DOC-36 §8).
    COMBINE = "Combine" [PROTECTED];
    TRANSCENDENTAL_SIMPLIFY = "TranscendentalSimplify" [PROTECTED];
    SIGN_TO_HEAVISIDE = "SignToHeaviside" [PROTECTED];
    HEAVISIDE_TO_SIGN = "HeavisideToSign" [PROTECTED];
    // Steps, deltas and periodic functions in a variable (CAS-570,
    // CAS-572, CAS-573).
    LINEAR_STEP = "LinearStep" [PROTECTED];
    SIMPLIFY_DIRAC = "SimplifyDirac" [PROTECTED];
    PERIODIC = "Periodic" [PROTECTED];
    // Integration by parts, values at bounds and derivatives of functions
    // (CAS-574, CAS-577, CAS-578).
    INTEGRATE_BY_PARTS_DV = "IntegrateByPartsDV" [PROTECTED];
    INTEGRATE_BY_PARTS_U = "IntegrateByPartsU" [PROTECTED];
    PREVAL = "Preval" [PROTECTED];
    FUNCTION_DIFF = "FunctionDiff" [PROTECTED];
    // Series reversion, Riemann sums and WZ certificates (CAS-576,
    // CAS-579, CAS-580).
    REVERT = "Revert" [PROTECTED];
    SUM_RIEMANN = "SumRiemann" [PROTECTED];
    WZ_CERTIFICATE = "WZCertificate" [PROTECTED];
    // The calculus of variations (CAS-575).
    EULER_LAGRANGE = "EulerLagrange" [PROTECTED];
    JACOBI_EQUATION = "JacobiEquation" [PROTECTED];
    CONJUGATE_EQUATION = "ConjugateEquation" [PROTECTED];
    // Function study: domains and convexity (CAS-581, CAS-583).
    FUNCTION_DOMAIN = "FunctionDomain" [PROTECTED];
    CONVEX = "Convex" [PROTECTED];
    // Variation and sign tables (CAS-582).
    TAB_VAR = "TabVar" [PROTECTED];
    TAB_SIGN = "TabSign" [PROTECTED];
    // Vector calculus (CAS-DOC-37): the differential operators (§2).
    GRAD = "Grad" [PROTECTED];
    DIVERGENCE = "Divergence" [PROTECTED];
    CURL = "Curl" [PROTECTED];
    LAPLACIAN = "Laplacian" [PROTECTED];
    HESSIAN = "Hessian" [PROTECTED];
    // Potentials (CAS-DOC-37 §3).
    POTENTIAL = "Potential" [PROTECTED];
    VECTOR_POTENTIAL = "VectorPotential" [PROTECTED];
    // Curves (CAS-DOC-37 §4).
    CURVATURE = "Curvature" [PROTECTED];
    FRENET = "Frenet" [PROTECTED];
    EVOLUTE = "Evolute" [PROTECTED];
    ARC_LENGTH = "ArcLength" [PROTECTED];
    // Implicit differentiation (CAS-DOC-37 §5).
    IMPLICIT_D = "ImplicitD" [PROTECTED];
    // Optimization (CAS-DOC-38).
    EXTREMA = "Extrema" [PROTECTED];
    MINIMIZE = "Minimize" [PROTECTED];
    MAXIMIZE = "Maximize" [PROTECTED];
    F_MIN = "FMin" [PROTECTED];
    F_MAX = "FMax" [PROTECTED];
    LP_SOLVE = "LPSolve" [PROTECTED];
    TP_SOLVE = "TPSolve" [PROTECTED];
    NLP_SOLVE = "NLPSolve" [PROTECTED];
    // giac's simplex_reduce (CAS-606).
    SIMPLEX_REDUCE = "SimplexReduce" [PROTECTED];
    // giac's find_minimum (CAS-608).
    FIND_MINIMUM = "FindMinimum" [PROTECTED];
    // giac's frank_wolfe (CAS-609).
    FRANK_WOLFE = "FrankWolfe" [PROTECTED];
    // giac's simulated_annealing (CAS-610).
    SIMULATED_ANNEALING = "SimulatedAnnealing" [PROTECTED];
    // LPSolve's options (CAS-600, CAS-601): no native.
    LP_ASSUME = "LPAssume" [PROTECTED];
    LP_MAXIMIZE = "LPMaximize" [PROTECTED];
    LP_NON_NEGATIVE = "LPNonNegative" [PROTECTED];
    LP_NON_NEG_INT = "LPNonNegInt" [PROTECTED];
    LP_INTEGER = "LPInteger" [PROTECTED];
    LP_BINARY = "LPBinary" [PROTECTED];
    LP_INTEGER_VARIABLES = "LPIntegerVariables" [PROTECTED];
    LP_BINARY_VARIABLES = "LPBinaryVariables" [PROTECTED];
    LP_METHOD = "LPMethod" [PROTECTED];
    LP_SIMPLEX = "LPSimplex" [PROTECTED];
    LP_INTERIOR_POINT = "LPInteriorPoint" [PROTECTED];
    LP_VERBOSE = "LPVerbose" [PROTECTED];
    LP_PRESOLVE = "LPPresolve" [PROTECTED];
    LP_NODE_SELECT = "LPNodeSelect" [PROTECTED];
    LP_DEPTH_FIRST = "LPDepthFirst" [PROTECTED];
    LP_BREADTH_FIRST = "LPBreadthFirst" [PROTECTED];
    LP_BEST_LOCAL_BOUND = "LPBestLocalBound" [PROTECTED];
    LP_BEST_PROJECTION = "LPBestProjection" [PROTECTED];
    LP_HYBRID = "LPHybrid" [PROTECTED];
    LP_VAR_SELECT = "LPVarSelect" [PROTECTED];
    LP_FIRST_FRACTIONAL = "LPFirstFractional" [PROTECTED];
    LP_LAST_FRACTIONAL = "LPLastFractional" [PROTECTED];
    LP_MOST_FRACTIONAL = "LPMostFractional" [PROTECTED];
    LP_PSEUDO_COST = "LPPseudoCost" [PROTECTED];
    LP_MAX_CUTS = "LPMaxCuts" [PROTECTED];
    LP_HEURISTIC = "LPHeuristic" [PROTECTED];
    LP_DEPTH_LIMIT = "LPDepthLimit" [PROTECTED];
    LP_NODE_LIMIT = "LPNodeLimit" [PROTECTED];
    LP_ITERATION_LIMIT = "LPIterationLimit" [PROTECTED];
    LP_TIME_LIMIT = "LPTimeLimit" [PROTECTED];
    LP_GAP_TOLERANCE = "LPGapTolerance" [PROTECTED];
    // NLPSolve's options (CAS-604): no native.
    NLP_INITIAL_POINT = "NLPInitialPoint" [PROTECTED];
    NLP_ITERATION_LIMIT = "NLPIterationLimit" [PROTECTED];
    NLP_NON_NEGATIVE = "NLPNonNegative" [PROTECTED];
    NLP_PRECISION = "NLPPrecision" [PROTECTED];
    NLP_MAXIMIZE = "NLPMaximize" [PROTECTED];
    NLP_PRESOLVE = "NLPPresolve" [PROTECTED];
    NLP_METHOD = "NLPMethod" [PROTECTED];
    NLP_INTEGER_VARIABLES = "NLPIntegerVariables" [PROTECTED];
    NLP_INTEGER = "NLPInteger" [PROTECTED];
    NLP_BINARY_VARIABLES = "NLPBinaryVariables" [PROTECTED];
    NLP_BINARY = "NLPBinary" [PROTECTED];
    NLP_TOLERANCE = "NLPTolerance" [PROTECTED];
    NLP_VERBOSE = "NLPVerbose" [PROTECTED];
    NLP_NON_NEG_INT = "NLPNonNegInt" [PROTECTED];
    // Number theory II (CAS-DOC-39).
    QUOTIENT_REMAINDER = "QuotientRemainder" [PROTECTED];
    SYMMETRIC_MOD = "SymmetricMod" [PROTECTED];
    PSEUDOPRIME_Q = "PseudoprimeQ" [PROTECTED];
    CYCLOTOMIC = "Cyclotomic" [PROTECTED];
    MULTIPLICATIVE_ORDER = "MultiplicativeOrder" [PROTECTED];
    PRIMITIVE_ROOT = "PrimitiveRoot" [PROTECTED];
    CONTINUED_FRACTION = "ContinuedFraction" [PROTECTED];
    FROM_CONTINUED_FRACTION = "FromContinuedFraction" [PROTECTED];
    INTEGER_SOLVE_LINEAR = "IntegerSolveLinear" [PROTECTED];
    PRIME_SUM_OF_SQUARES = "PrimeSumOfSquares" [PROTECTED];
    INTEGER_SOLVE = "IntegerSolve" [PROTECTED];
    // Element-wise operators, next to `Times` and `Power` (CAS-DOC-4 §6.2).
    ELEMENTWISE_TIMES = "ElementwiseTimes" [FLAT, ONE_IDENTITY, LISTABLE, PROTECTED];
    ELEMENTWISE_DIVIDE = "ElementwiseDivide" [LISTABLE, PROTECTED];
    ELEMENTWISE_POWER = "ElementwisePower" [LISTABLE, PROTECTED];
    // Constants.
    PI = "Pi" [CONSTANT, PROTECTED];
    E = "E" [CONSTANT, PROTECTED];
    // Euler's constant γ, Xcas's `euler_gamma` (CAS-DOC-31 §3.2).
    EULER_GAMMA = "EulerGamma" [CONSTANT, PROTECTED];
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
    // Xcas's `all_trig_solutions` and `with_sqrt`, flags of the context
    // read as 0 or 1; the session sets them (CAS-545).
    ALL_TRIG_SOLUTIONS = "AllTrigSolutions" [PROTECTED];
    WITH_SQRT = "WithSqrt" [PROTECTED];
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
    PIECEWISE = "Piecewise" [PROTECTED];
    PIECEWISE_MAP = "PiecewiseMap" [PROTECTED];
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
    // The gcd algorithms by name (CAS-526): Xcas's `heugcd`, `modgcd`,
    // `psrgcd`, `ezgcd` and `icontent`.
    HEURISTIC_GCD = "HeuristicGCD" [PROTECTED];
    MODULAR_GCD = "ModularGCD" [PROTECTED];
    SUBRESULTANT_GCD = "SubresultantGCD" [PROTECTED];
    EZ_GCD = "EZGCD" [PROTECTED];
    INTEGER_CONTENT = "IntegerContent" [LISTABLE, PROTECTED];
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
    // values: `plex`, `tdeg` and `revlex`. Xcas's `in_ideal`, `syzygy`,
    // `gbasis_vars` and `gbasis`'s option `rur` (CAS-541).
    GROEBNER_BASIS = "GroebnerBasis" [PROTECTED];
    GROEBNER_REDUCE = "GroebnerReduce" [PROTECTED];
    ELIMINATE = "Eliminate" [PROTECTED];
    IDEAL_MEMBER_Q = "IdealMemberQ" [PROTECTED];
    SYZYGIES = "Syzygies" [PROTECTED];
    GROEBNER_VARIABLES = "GroebnerVariables" [PROTECTED];
    RATIONAL_UNIVARIATE_REPRESENTATION = "RationalUnivariateRepresentation" [PROTECTED];
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
    // Set operations (CAS-DOC-32): giac's infix `union`, `intersect` and
    // `minus`, and `is_included`.
    UNION = "Union" [FLAT, ORDERLESS, PROTECTED];
    INTERSECTION = "Intersection" [FLAT, ORDERLESS, PROTECTED];
    COMPLEMENT = "Complement" [PROTECTED];
    INCLUDED_Q = "IncludedQ" [PROTECTED];
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
