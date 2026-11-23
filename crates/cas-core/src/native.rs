//! The native functions of the builtins (CAS-DOC-4 §6, CAS-54): arithmetic
//! and the infinities, comparisons and predicates, logic, structure, and
//! what the language lowers to (CAS-DOC-7 §9, CAS-123): indexing, sets,
//! factorials and `Unapply`.
//!
//! Each is a [`NativeFn`] in a table indexed by builtin, so dispatching to
//! one is an array read. Each receives the application after the
//! evaluator normalised it, and returns its rewrite or `None`.

#![allow(
    clippy::unnecessary_wraps,
    reason = "natives share the `NativeFn` signature, which can fail"
)]

mod arith;
mod ball;
mod bessel;
mod bits;
mod calculus_tools;
mod dft;
mod dsolve;
mod elementary;
mod erf;
mod expint;
mod expressions;
mod factorial;
mod fourier;
mod fourier_transform;
mod groebner;
mod incomplete;
mod index;
mod integer;
mod integrate;
mod interpolation;
mod lambert;
mod laplace;
mod lattice;
mod limit;
mod linear;
mod lists;
mod logic;
mod matrix;
mod modular;
mod ndsolve;
mod nintegrate;
mod number;
mod number_theory;
mod numeric;
mod optimization;
mod orthogonal;
mod parts;
mod permutations;
pub(crate) mod piecewise;
mod polygamma;
mod polynomial;
mod polytools;
mod procedure;
mod radical;
mod random;
mod replace;
mod residue;
mod rewrite;
mod rounding;
mod rsolve;
mod series;
mod sets;
mod signal;
mod simplify;
mod solve;
mod special;
mod statistics;
mod step;
mod structure;
mod sum;
mod trig;
mod vector;
mod zeta;
mod ztransform;

pub(crate) use numeric::approximation;

use crate::build;
use crate::builtin;
use crate::context::Context;
use crate::eval::NativeFn;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The native functions of the builtins, by id: dispatching to one is an
/// array read (CAS-DOC-4 §6).
static NATIVES: [Option<NativeFn>; builtin::COUNT] = {
    let entries: &[(Symbol, NativeFn)] = &[
        (Symbol::PLUS, arith::plus),
        (Symbol::TIMES, arith::times),
        (Symbol::POWER, arith::power),
        (Symbol::SQRT, arith::sqrt),
        (Symbol::SURD, arith::surd),
        (Symbol::EXP, elementary::exp),
        (Symbol::LOG, elementary::log),
        (Symbol::SIN, elementary::sin),
        (Symbol::COS, elementary::cos),
        (Symbol::TAN, elementary::tan),
        (Symbol::ARC_SIN, elementary::arc_sin),
        (Symbol::ARC_COS, elementary::arc_cos),
        (Symbol::ARC_TAN, elementary::arc_tan),
        (Symbol::SINH, elementary::sinh),
        (Symbol::COSH, elementary::cosh),
        (Symbol::TANH, elementary::tanh),
        (Symbol::ARC_SINH, elementary::arc_sinh),
        (Symbol::ARC_COSH, elementary::arc_cosh),
        (Symbol::ARC_TANH, elementary::arc_tanh),
        (Symbol::SEC, elementary::reciprocal),
        (Symbol::CSC, elementary::reciprocal),
        (Symbol::COT, elementary::reciprocal),
        (Symbol::ARC_SEC, elementary::reciprocal),
        (Symbol::ARC_CSC, elementary::reciprocal),
        (Symbol::ARC_COT, elementary::reciprocal),
        (Symbol::RE, parts::part),
        (Symbol::IM, parts::part),
        (Symbol::CONJUGATE, parts::part),
        (Symbol::ARG, parts::arg),
        (Symbol::FLOOR, rounding::round),
        (Symbol::CEILING, rounding::round),
        (Symbol::ROUND, rounding::round),
        (Symbol::SIGN, rounding::sign),
        (Symbol::MAX, rounding::extremum),
        (Symbol::MIN, rounding::extremum),
        (Symbol::GAMMA, factorial::gamma),
        (Symbol::BETA, factorial::beta),
        (Symbol::LOWER_GAMMA, incomplete::lower_gamma),
        (Symbol::ZETA, zeta::zeta),
        (Symbol::ERF, erf::erf),
        (Symbol::ERFC, erf::erf),
        (Symbol::EXP_INTEGRAL_EI, expint::integral),
        (Symbol::SIN_INTEGRAL, expint::integral),
        (Symbol::COS_INTEGRAL, expint::integral),
        (Symbol::LOG_INTEGRAL, expint::log_integral),
        (Symbol::POLY_GAMMA, polygamma::polygamma),
        (Symbol::LAMBERT_W, lambert::lambert_w),
        (Symbol::NDSOLVE, ndsolve::ndsolve),
        (Symbol::BVP_SOLVE, ndsolve::bvp::bvp_solve),
        (Symbol::BARTLETT_HANN_WINDOW, signal::window),
        (Symbol::BLACKMAN_HARRIS_WINDOW, signal::window),
        (Symbol::BLACKMAN_WINDOW, signal::window),
        (Symbol::BOHMAN_WINDOW, signal::window),
        (Symbol::COSINE_WINDOW, signal::window),
        (Symbol::GAUSSIAN_WINDOW, signal::window),
        (Symbol::HAMMING_WINDOW, signal::window),
        (Symbol::HANN_POISSON_WINDOW, signal::window),
        (Symbol::HANN_WINDOW, signal::window),
        (Symbol::PARZEN_WINDOW, signal::window),
        (Symbol::POISSON_WINDOW, signal::window),
        (Symbol::RIEMANN_WINDOW, signal::window),
        (Symbol::TRIANGLE_WINDOW, signal::window),
        (Symbol::TUKEY_WINDOW, signal::window),
        (Symbol::WELCH_WINDOW, signal::window),
        (Symbol::LOW_PASS, signal::filter),
        (Symbol::HIGH_PASS, signal::filter),
        (Symbol::CROSS_CORRELATION, signal::correlation),
        (Symbol::AUTO_CORRELATION, signal::correlation),
        (Symbol::THRESHOLD, signal::threshold),
        (Symbol::BESSEL_J, bessel::bessel),
        (Symbol::BESSEL_Y, bessel::bessel),
        (Symbol::BESSEL_I, bessel::bessel),
        (Symbol::BESSEL_K, bessel::bessel),
        (Symbol::AIRY_AI, bessel::airy),
        (Symbol::AIRY_BI, bessel::airy),
        (Symbol::AIRY_AI_PRIME, bessel::airy),
        (Symbol::AIRY_BI_PRIME, bessel::airy),
        (Symbol::UNION, sets::union),
        (Symbol::INTERSECTION, sets::intersection),
        (Symbol::COMPLEMENT, sets::complement),
        (Symbol::INCLUDED_Q, sets::included_q),
        (Symbol::HEAVISIDE_THETA, step::heaviside),
        (Symbol::DIRAC_DELTA, step::dirac),
        (Symbol::UNIT_BOX, step::unit_box),
        (Symbol::DIRECTED_INFINITY, arith::directed_infinity),
        (Symbol::EQUAL, logic::equal_native),
        (Symbol::UNEQUAL, logic::unequal),
        (Symbol::LESS, logic::less),
        (Symbol::LESS_EQUAL, logic::less_equal),
        (Symbol::GREATER, logic::greater),
        (Symbol::GREATER_EQUAL, logic::greater_equal),
        (Symbol::SAME_Q, logic::same_q),
        (Symbol::FREE_Q, logic::free_q),
        (Symbol::INTEGER_Q, logic::integer_q),
        (Symbol::RATIONAL_Q, logic::rational_q),
        (Symbol::REAL_Q, logic::real_q),
        (Symbol::COMPLEX_Q, logic::complex_q),
        (Symbol::NUMBER_Q, logic::number_q),
        (Symbol::LIST_Q, logic::list_q),
        (Symbol::MATRIX_Q, logic::matrix_q),
        (Symbol::STRING_Q, logic::string_q),
        (Symbol::SYMBOL_Q, logic::symbol_q),
        (Symbol::POSITIVE, logic::sign_predicate),
        (Symbol::NEGATIVE, logic::sign_predicate),
        (Symbol::NON_NEGATIVE, logic::sign_predicate),
        (Symbol::NON_POSITIVE, logic::sign_predicate),
        (Symbol::ELEMENT, logic::element),
        (Symbol::ABS, arith::abs),
        (Symbol::REMAINDER, integer::remainder),
        (Symbol::QUOTIENT, integer::quotient),
        (Symbol::EVEN_Q, integer::even_q),
        (Symbol::ODD_Q, integer::odd_q),
        (Symbol::ELEMENTWISE_TIMES, linear::elementwise_times),
        (Symbol::ELEMENTWISE_DIVIDE, linear::elementwise_divide),
        (Symbol::ELEMENTWISE_POWER, linear::elementwise_power),
        (Symbol::FACTORIAL, factorial::factorial),
        (Symbol::AND, logic::and),
        (Symbol::OR, logic::or),
        (Symbol::NOT, logic::not),
        (Symbol::IF, logic::if_native),
        (Symbol::PIECEWISE, piecewise::piecewise_native),
        (Symbol::PIECEWISE_MAP, piecewise::piecewise_map),
        (Symbol::LENGTH, structure::length),
        (Symbol::PART, structure::part),
        (Symbol::MAP, structure::map),
        (Symbol::APPLY, structure::apply),
        (Symbol::FINITE_SET, structure::finite_set),
        (Symbol::UNAPPLY, structure::unapply),
        (Symbol::AT, index::at),
        (Symbol::REPLACE_AT, index::replace_at),
        (Symbol::BLOCK, procedure::block),
        (Symbol::BALL, ball::ball),
        (Symbol::DECIMAL, ball::decimal),
        (Symbol::EXACT, ball::exact),
        (Symbol::EXPAND, polynomial::expand),
        (Symbol::TOGETHER, polynomial::together),
        (Symbol::GCD, polynomial::gcd),
        (Symbol::LCM, polynomial::lcm),
        (Symbol::HEURISTIC_GCD, polynomial::heuristic_gcd),
        (Symbol::MODULAR_GCD, polynomial::modular_gcd),
        (Symbol::SUBRESULTANT_GCD, polynomial::subresultant_gcd),
        (Symbol::EZ_GCD, polynomial::ez_gcd),
        (Symbol::INTEGER_CONTENT, polynomial::integer_content),
        (Symbol::POLYNOMIAL_QUOTIENT, polynomial::polynomial_quotient),
        (
            Symbol::POLYNOMIAL_REMAINDER,
            polynomial::polynomial_remainder,
        ),
        (
            Symbol::POLYNOMIAL_QUOTIENT_REMAINDER,
            polynomial::polynomial_quotient_remainder,
        ),
        (Symbol::EXPONENT, polynomial::exponent),
        (Symbol::COEFFICIENT, polynomial::coefficient),
        (Symbol::LEADING_COEFFICIENT, polynomial::leading_coefficient),
        (Symbol::NUMERATOR, polynomial::numerator),
        (Symbol::DENOMINATOR, polynomial::denominator),
        (Symbol::CONTENT, polynomial::content),
        (Symbol::PRIMITIVE_PART, polynomial::primitive_part),
        (Symbol::RESULTANT, polynomial::resultant),
        (Symbol::DISCRIMINANT, polynomial::discriminant),
        (Symbol::FACTOR, polynomial::factor),
        (Symbol::GROEBNER_BASIS, groebner::groebner_basis),
        (Symbol::GROEBNER_REDUCE, groebner::groebner_reduce),
        (Symbol::ELIMINATE, groebner::eliminate),
        (Symbol::IDEAL_MEMBER_Q, groebner::ideal_member_q),
        (Symbol::SYZYGIES, groebner::syzygies),
        (Symbol::GROEBNER_VARIABLES, groebner::groebner_variables),
        (Symbol::INTEGRATE, integrate::integrate),
        (Symbol::SERIES, series::series),
        (Symbol::LIMIT, limit::limit),
        (Symbol::DSOLVE, dsolve::dsolve),
        (Symbol::LAPLACE_TRANSFORM, laplace::laplace_transform),
        (
            Symbol::FOURIER_TRANSFORM,
            fourier_transform::fourier_transform,
        ),
        (
            Symbol::INVERSE_FOURIER_TRANSFORM,
            fourier_transform::inverse_fourier_transform,
        ),
        (
            Symbol::HERMITE_DECOMPOSITION,
            lattice::hermite_decomposition,
        ),
        (Symbol::SMITH_DECOMPOSITION, lattice::smith_decomposition),
        (Symbol::MODULAR, modular::modular_value),
        (Symbol::LATTICE_REDUCE, lattice::lattice_reduce),
        (
            Symbol::INVERSE_LAPLACE_TRANSFORM,
            laplace::inverse_laplace_transform,
        ),
        (Symbol::ZTRANSFORM, ztransform::z_transform),
        (Symbol::INVERSE_ZTRANSFORM, ztransform::inverse_z_transform),
        (Symbol::KRONECKER_DELTA, ztransform::kronecker_delta),
        (Symbol::RSOLVE, rsolve::rsolve),
        (Symbol::SEQ_SOLVE, rsolve::seqsolve),
        (Symbol::RESIDUE, residue::residue),
        (
            Symbol::FOURIER_COS_COEFFICIENT,
            fourier::fourier_cos_coefficient,
        ),
        (
            Symbol::FOURIER_SIN_COEFFICIENT,
            fourier::fourier_sin_coefficient,
        ),
        (Symbol::FOURIER_COEFFICIENT, fourier::fourier_coefficient),
        (Symbol::NINTEGRATE, nintegrate::nintegrate),
        (
            Symbol::INTERPOLATING_POLYNOMIAL,
            interpolation::interpolating_polynomial,
        ),
        (Symbol::SPLINE, interpolation::spline),
        (
            Symbol::THIELE_INTERPOLATION,
            interpolation::thiele_interpolation,
        ),
        (Symbol::LEGENDRE_P, orthogonal::legendre_p),
        (Symbol::HERMITE_H, orthogonal::hermite_h),
        (Symbol::LAGUERRE_L, orthogonal::laguerre_l),
        (Symbol::CHEBYSHEV_T, orthogonal::chebyshev_t),
        (Symbol::CHEBYSHEV_U, orthogonal::chebyshev_u),
        (Symbol::TABLE, lists::table),
        (Symbol::MAKE_LIST, lists::make_list),
        (Symbol::SORT, lists::sort),
        (Symbol::REVERSE, lists::reverse),
        (Symbol::ROTATE, lists::rotate),
        (Symbol::JOIN, lists::join),
        (Symbol::APPEND, lists::append),
        (Symbol::PREPEND, lists::prepend),
        (Symbol::FIRST, lists::first),
        (Symbol::REST, lists::rest),
        (Symbol::FLATTEN, lists::flatten),
        (Symbol::DIMENSIONS, lists::dimensions),
        (Symbol::SELECT, lists::select),
        (Symbol::REJECT, lists::reject),
        (Symbol::COUNT, lists::count),
        (Symbol::MEMBER_INDEX, lists::member_index),
        (Symbol::ZIP, lists::zip),
        (Symbol::ACCUMULATE, lists::accumulate),
        (Symbol::FIND_POSITIONS, lists::xcas::find_positions),
        (Symbol::FIRST_POSITION, lists::xcas::first_position),
        (Symbol::COUNT_EQUAL, lists::xcas::count_equal),
        (Symbol::COUNT_LESS, lists::xcas::count_less),
        (Symbol::COUNT_GREATER, lists::xcas::count_greater),
        (Symbol::DIFFERENCES, lists::xcas::differences),
        (
            Symbol::ACCUMULATE_HEAD_TAIL,
            lists::xcas::accumulate_head_tail,
        ),
        (Symbol::LENGTHS, lists::xcas::lengths),
        (Symbol::MAKE_VECTOR, lists::xcas::make_vector),
        (Symbol::MAKE_SEQUENCE, lists::xcas::make_sequence),
        (Symbol::NEW_LIST, lists::xcas::new_list),
        (Symbol::ENUMERATE, lists::xcas::enumerate),
        (Symbol::COPY, lists::xcas::copy),
        (Symbol::APPLY_EACH, lists::xcas::apply_each),
        (Symbol::FOLD_LEFT, lists::xcas::fold_left),
        (Symbol::FOLD_RIGHT, lists::xcas::fold_right),
        (Symbol::MULTIPLY, lists::xcas::multiply),
        (Symbol::INSERT_AT, lists::xcas::insert_at),
        (
            Symbol::SUBSTITUTE_OPERANDS,
            lists::xcas::substitute_operands,
        ),
        (Symbol::HEAPIFY, lists::xcas::heapify),
        (Symbol::HEAP_PUSH, lists::xcas::heap_push),
        (Symbol::HEAP_POP, lists::xcas::heap_pop),
        (Symbol::POWER_EXPAND, rewrite::power_expand),
        (Symbol::COMMON_DENOMINATOR, polynomial::common_denominator),
        (Symbol::RAW_NUMERATOR, polynomial::raw_numerator),
        (Symbol::RAW_DENOMINATOR, polynomial::raw_denominator),
        (Symbol::CANCEL_PAIR, polynomial::cancel_pair),
        (Symbol::REGROUP, simplify::regroup),
        (Symbol::SUM_OF_SQUARES, arith::sum_of_squares),
        (Symbol::LOG_B, elementary::log_b),
        (Symbol::LOG2, elementary::log_b),
        (Symbol::LOG10, elementary::log_b),
        (Symbol::EXP10, elementary::exp10),
        (Symbol::RADIANS, elementary::angle),
        (Symbol::DEGREES, elementary::angle),
        (Symbol::NORMALIZE, matrix::constructors::normalize),
        (Symbol::ALGEBRAIC_SUBSTITUTE, replace::algebraic_substitute),
        (Symbol::IDENTIFIERS, expressions::identifiers),
        (Symbol::KERNELS, expressions::kernels),
        (
            Symbol::ALGEBRAIC_VARIABLES,
            expressions::algebraic_variables,
        ),
        (Symbol::REORDER, expressions::reorder),
        (Symbol::EVALUATE, expressions::evaluate),
        (Symbol::EVALUATE_BOOLEAN, expressions::evaluate_boolean),
        (Symbol::EPSILON_TO_ZERO, expressions::epsilon_to_zero),
        (Symbol::INFINITE_Q, expressions::infinite_q),
        (Symbol::INDETERMINATE_Q, expressions::indeterminate_q),
        (Symbol::COMPARE, expressions::compare),
        (
            Symbol::EQUATION_TO_DIFFERENCE,
            expressions::equation_to_difference,
        ),
        (Symbol::EQUATION_TO_LIST, expressions::equation_to_list),
        (Symbol::LIST_TO_EQUATIONS, expressions::list_to_equations),
        (Symbol::EQUATIONS_TO_LIST, expressions::equations_to_list),
        (Symbol::APART, polytools::apart),
        (Symbol::COMPLEX_APART, polytools::complex_apart),
        (Symbol::PROPER_FRACTION, polytools::proper_fraction),
        (
            Symbol::NUMERATOR_DENOMINATOR,
            polytools::numerator_denominator,
        ),
        (Symbol::COEFFICIENT_LIST, polytools::coefficient_list),
        (Symbol::FROM_COEFFICIENTS, polytools::from_coefficients),
        (
            Symbol::TRAILING_COEFFICIENT,
            polytools::trailing_coefficient,
        ),
        (Symbol::VALUATION, polytools::valuation),
        (Symbol::COUNT_ROOTS, polytools::count_roots),
        (Symbol::REAL_ROOTS, polytools::real_roots),
        (
            Symbol::REAL_ROOT_ISOLATION,
            polytools::roots::real_root_isolation,
        ),
        (
            Symbol::POSITIVE_ROOT_ISOLATION,
            polytools::roots::positive_root_isolation,
        ),
        (
            Symbol::POSITIVE_ROOT_UPPER_BOUND,
            polytools::roots::positive_root_upper_bound,
        ),
        (
            Symbol::POSITIVE_ROOT_LOWER_BOUND,
            polytools::roots::positive_root_lower_bound,
        ),
        (Symbol::STURM_SEQUENCE, polytools::roots::sturm_sequence),
        (Symbol::RATIONAL_ROOTS, polytools::roots::rational_roots),
        (
            Symbol::COMPLEX_RATIONAL_ROOTS,
            polytools::roots::complex_rational_roots,
        ),
        (Symbol::COMPLEX_ROOTS, polytools::roots::complex_roots),
        (
            Symbol::POLYNOMIAL_ROOTS,
            polytools::two::polynomial_roots_native,
        ),
        (Symbol::FROM_ROOTS, polytools::two::from_roots),
        (Symbol::HORNER_EVAL, polytools::two::horner_eval),
        (Symbol::TAYLOR_SHIFT, polytools::two::taylor_shift),
        (
            Symbol::ROOT_MULTIPLICITIES,
            polytools::two::root_multiplicities,
        ),
        (
            Symbol::POLYNOMIAL_EXTENDED_GCD,
            polytools::two::polynomial_extended_gcd,
        ),
        (Symbol::POLYNOMIAL_BEZOUT, polytools::two::polynomial_bezout),
        (
            Symbol::POLYNOMIAL_CHINESE_REMAINDER,
            polytools::two::polynomial_chinese_remainder,
        ),
        (Symbol::SERIES_DIVISION, polytools::two::series_division),
        (Symbol::PADE_APPROXIMANT, polytools::two::pade_approximant),
        (Symbol::COMPLEX_FACTOR, polytools::two::complex_factor),
        (Symbol::COMPLETE_SQUARE, polytools::two::complete_square),
        (Symbol::DIGIT_POLYNOMIAL, polytools::two::digit_polynomial),
        (Symbol::BERLEKAMP_MASSEY, polytools::two::berlekamp_massey),
        (Symbol::ROOT_OF, polytools::three::root_of),
        (
            Symbol::FACTORED_FROM_ROOTS,
            polytools::three::factored_from_roots,
        ),
        (Symbol::HOMOGENIZE, polytools::three::homogenize),
        (Symbol::TOTAL_DEGREE, polytools::three::total_degree),
        (Symbol::FACTOR_XN, polytools::three::factor_xn),
        (
            Symbol::BEZOUT_MATRIX,
            matrix::decompose::forms::bezout_matrix,
        ),
        (
            Symbol::IRREDUCIBLE_POLYNOMIAL_Q,
            polytools::three::irreducible_polynomial_q,
        ),
        (
            Symbol::LU_DECOMPOSITION,
            matrix::decompose::lu_decomposition,
        ),
        (
            Symbol::QR_DECOMPOSITION,
            matrix::decompose::qr_decomposition,
        ),
        (
            Symbol::CHOLESKY_DECOMPOSITION,
            matrix::decompose::cholesky_decomposition,
        ),
        (
            Symbol::JORDAN_DECOMPOSITION,
            matrix::decompose::jordan_decomposition,
        ),
        (
            Symbol::MINIMAL_POLYNOMIAL,
            matrix::decompose::minimal_polynomial,
        ),
        (Symbol::MATRIX_EXP, matrix::decompose::matrix_exp),
        (Symbol::MATRIX_POWER, matrix::decompose::matrix_power),
        (
            Symbol::DIAGONAL_MATRIX,
            matrix::constructors::diagonal_matrix,
        ),
        (Symbol::HILBERT_MATRIX, matrix::constructors::hilbert_matrix),
        (
            Symbol::VANDERMONDE_MATRIX,
            matrix::constructors::vandermonde_matrix,
        ),
        (Symbol::NORM, matrix::constructors::norm),
        (Symbol::FROBENIUS_NORM, matrix::constructors::frobenius_norm),
        (
            Symbol::HESSENBERG_DECOMPOSITION,
            matrix::decompose::forms::hessenberg_decomposition,
        ),
        (
            Symbol::COMPANION_MATRIX,
            matrix::decompose::forms::companion_matrix,
        ),
        (
            Symbol::SYLVESTER_MATRIX,
            matrix::decompose::forms::sylvester_matrix,
        ),
        (
            Symbol::HADAMARD_PRODUCT,
            matrix::decompose::forms::hadamard_product,
        ),
        (
            Symbol::LDL_DECOMPOSITION,
            matrix::decompose::forms::ldl_decomposition,
        ),
        (
            Symbol::POSITIVE_DEFINITE_Q,
            matrix::decompose::forms::positive_definite_q,
        ),
        (Symbol::INERTIA, matrix::decompose::forms::inertia),
        (Symbol::MATRIX_NORM, matrix::decompose::norms::matrix_norm),
        (Symbol::L1_NORM, matrix::decompose::norms::l1_norm),
        (
            Symbol::L_INFINITY_NORM,
            matrix::decompose::norms::l_infinity_norm,
        ),
        (Symbol::MAX_NORM, matrix::decompose::norms::max_norm),
        (
            Symbol::CONDITION_NUMBER,
            matrix::decompose::norms::condition_number,
        ),
        (Symbol::GRAM_SCHMIDT, matrix::decompose::norms::gram_schmidt),
        (
            Symbol::SINGULAR_VALUES,
            matrix::decompose::norms::singular_values_native,
        ),
        (
            Symbol::CONJUGATE_GRADIENT,
            matrix::decompose::norms::conjugate_gradient,
        ),
        (
            Symbol::SCHUR_DECOMPOSITION,
            matrix::decompose::unitary::schur_decomposition,
        ),
        (
            Symbol::SINGULAR_VALUE_DECOMPOSITION,
            matrix::decompose::unitary::singular_value_decomposition,
        ),
        (
            Symbol::JACOBI_LINEAR_SOLVE,
            matrix::decompose::iterative::jacobi_linear_solve,
        ),
        (
            Symbol::GAUSS_SEIDEL_LINEAR_SOLVE,
            matrix::decompose::iterative::gauss_seidel_linear_solve,
        ),
        (Symbol::PRIME_Q, number_theory::prime_q),
        (Symbol::NEXT_PRIME, number_theory::next_prime),
        (Symbol::PREVIOUS_PRIME, number_theory::previous_prime),
        (Symbol::INTEGER_FACTOR, number_theory::integer_factor),
        (
            Symbol::INTEGER_FACTOR_LIST,
            number_theory::integer_factor_list,
        ),
        (Symbol::ECM_FACTOR, number_theory::ecm_factor),
        (Symbol::DIVISORS, number_theory::divisors),
        (Symbol::EULER_PHI, number_theory::euler_phi),
        (Symbol::POWER_MOD, number_theory::power_mod),
        (Symbol::EXTENDED_GCD, number_theory::extended_gcd),
        (Symbol::CHINESE_REMAINDER, number_theory::chinese_remainder),
        (Symbol::JACOBI_SYMBOL, number_theory::jacobi_symbol),
        (Symbol::LEGENDRE_SYMBOL, number_theory::jacobi_symbol),
        (Symbol::BERNOULLI_B, number_theory::bernoulli_b),
        (Symbol::PRIME, number_theory::prime),
        (Symbol::PRIME_PI, number_theory::prime_pi_native),
        (
            Symbol::CHINESE_REMAINDER_LISTS,
            number_theory::integers::chinese_remainder_lists,
        ),
        (Symbol::LIST_GCD, number_theory::integers::list_gcd),
        (Symbol::FRACTION_MOD, number_theory::integers::fraction_mod),
        (
            Symbol::INTEGER_COMPOSITIONS,
            number_theory::integers::integer_compositions,
        ),
        (
            Symbol::BERNOULLI_MOD,
            number_theory::integers::bernoulli_mod_native,
        ),
        (Symbol::IS_PRIME, number_theory::integers::is_prime),
        (Symbol::PERMUTATION_Q, permutations::permutation_q),
        (Symbol::PERMUTATION_CYCLES, permutations::permutation_cycles),
        (Symbol::CYCLES_PERMUTATION, permutations::cycles_permutation),
        (Symbol::CYCLE_PERMUTATION, permutations::cycle_permutation),
        (
            Symbol::INVERSE_PERMUTATION,
            permutations::inverse_permutation,
        ),
        (Symbol::INVERSE_CYCLE, permutations::inverse_cycle),
        (
            Symbol::PERMUTATION_PRODUCT,
            permutations::permutation_product,
        ),
        (
            Symbol::PERMUTATION_CYCLE_PRODUCT,
            permutations::permutation_cycle_product,
        ),
        (
            Symbol::CYCLE_PERMUTATION_PRODUCT,
            permutations::cycle_permutation_product,
        ),
        (Symbol::CYCLE_PRODUCT, permutations::cycle_product),
        (Symbol::PERMUTATION_ORDER, permutations::permutation_order),
        (
            Symbol::PERMUTATION_SIGNATURE,
            permutations::permutation_signature,
        ),
        (Symbol::PERMUTATION_MATRIX, permutations::permutation_matrix),
        (Symbol::NEXT_PERMUTATION, permutations::next_permutation),
        (
            Symbol::PREVIOUS_PERMUTATION,
            permutations::previous_permutation,
        ),
        (Symbol::PERMUTATION_GROUP, permutations::permutation_group),
        (Symbol::BIT_AND, bits::bit_and),
        (Symbol::BIT_OR, bits::bit_or),
        (Symbol::BIT_XOR, bits::bit_xor),
        (Symbol::HAMMING_DISTANCE, bits::hamming_distance),
        (Symbol::CONVERT, bits::convert),
        (Symbol::FACTORIAL_POWER, bits::factorial_power),
        (
            Symbol::DISCRETE_FOURIER_TRANSFORM,
            dft::discrete_fourier_transform,
        ),
        (
            Symbol::INVERSE_DISCRETE_FOURIER_TRANSFORM,
            dft::inverse_discrete_fourier_transform,
        ),
        (Symbol::CONVOLUTION, dft::convolution),
        (Symbol::MEAN, statistics::mean),
        (Symbol::POPULATION_VARIANCE, statistics::population_variance),
        (
            Symbol::POPULATION_STANDARD_DEVIATION,
            statistics::population_standard_deviation,
        ),
        (Symbol::STANDARD_DEVIATION, statistics::standard_deviation),
        (Symbol::MEDIAN, statistics::median),
        (Symbol::QUANTILE, statistics::quantile),
        (Symbol::QUARTILES, statistics::quartiles),
        (Symbol::FIRST_QUARTILE, statistics::first_quartile),
        (Symbol::THIRD_QUARTILE, statistics::third_quartile),
        (
            Symbol::POPULATION_COVARIANCE,
            statistics::population_covariance,
        ),
        (Symbol::CORRELATION, statistics::correlation),
        (Symbol::FREQUENCIES, statistics::descriptive::frequencies),
        (
            Symbol::CUMULATED_FREQUENCIES,
            statistics::descriptive::cumulated_frequencies,
        ),
        (Symbol::CLASSES, statistics::descriptive::classes),
        (
            Symbol::COVARIANCE_CORRELATION,
            statistics::descriptive::covariance_correlation,
        ),
        (
            Symbol::ROOT_MEAN_SQUARE,
            statistics::descriptive::root_mean_square,
        ),
        (
            Symbol::MOVING_AVERAGE,
            statistics::descriptive::moving_average,
        ),
        (Symbol::Z_TEST, statistics::hypothesis::location::z_test),
        (Symbol::T_TEST, statistics::hypothesis::location::t_test),
        (
            Symbol::CHI_SQUARE_TEST,
            statistics::hypothesis::chisquare::chi_square_test,
        ),
        (
            Symbol::KOLMOGOROV_CDF,
            statistics::hypothesis::kolmogorov::kolmogorov_cdf_native,
        ),
        (
            Symbol::KOLMOGOROV_SMIRNOV_TEST,
            statistics::hypothesis::kolmogorov::kolmogorov_smirnov_test,
        ),
        (
            Symbol::WILCOXON_DISTRIBUTION,
            statistics::hypothesis::wilcoxon::wilcoxon_distribution,
        ),
        (
            Symbol::WILCOXON_STATISTIC,
            statistics::hypothesis::wilcoxon::wilcoxon_statistic,
        ),
        (
            Symbol::WILCOXON_TEST,
            statistics::hypothesis::wilcoxon::wilcoxon_test,
        ),
        (
            Symbol::LINEAR_REGRESSION,
            statistics::regression::linear_regression,
        ),
        (
            Symbol::POLYNOMIAL_REGRESSION,
            statistics::regression::polynomial_regression,
        ),
        (
            Symbol::EXPONENTIAL_REGRESSION,
            statistics::regression::exponential_regression,
        ),
        (
            Symbol::LOGARITHMIC_REGRESSION,
            statistics::regression::logarithmic_regression,
        ),
        (
            Symbol::POWER_REGRESSION,
            statistics::regression::power_regression,
        ),
        (
            Symbol::BINOMIAL_CDF,
            statistics::distributions::discrete::binomial_cdf,
        ),
        (
            Symbol::BINOMIAL_ICDF,
            statistics::distributions::discrete::binomial_icdf,
        ),
        (
            Symbol::POISSON_PDF,
            statistics::distributions::discrete::poisson_pdf,
        ),
        (
            Symbol::POISSON_CDF,
            statistics::distributions::discrete::poisson_cdf,
        ),
        (
            Symbol::POISSON_ICDF,
            statistics::distributions::discrete::poisson_icdf,
        ),
        (
            Symbol::GEOMETRIC_PDF,
            statistics::distributions::discrete::geometric_pdf,
        ),
        (
            Symbol::GEOMETRIC_CDF,
            statistics::distributions::discrete::geometric_cdf,
        ),
        (
            Symbol::GEOMETRIC_ICDF,
            statistics::distributions::discrete::geometric_icdf,
        ),
        (
            Symbol::NEGATIVE_BINOMIAL_PDF,
            statistics::distributions::discrete::negative_binomial_pdf,
        ),
        (
            Symbol::NEGATIVE_BINOMIAL_CDF,
            statistics::distributions::discrete::negative_binomial_cdf,
        ),
        (
            Symbol::NEGATIVE_BINOMIAL_ICDF,
            statistics::distributions::discrete::negative_binomial_icdf,
        ),
        (
            Symbol::INVERSE_ERF,
            statistics::distributions::inversion::inverse_erf,
        ),
        (
            Symbol::NORMAL_PDF,
            statistics::distributions::continuous::normal_pdf,
        ),
        (
            Symbol::NORMAL_CDF,
            statistics::distributions::continuous::normal_cdf,
        ),
        (
            Symbol::NORMAL_ICDF,
            statistics::distributions::continuous::normal_icdf,
        ),
        (
            Symbol::EXPONENTIAL_PDF,
            statistics::distributions::continuous::exponential_pdf,
        ),
        (
            Symbol::EXPONENTIAL_CDF,
            statistics::distributions::continuous::exponential_cdf,
        ),
        (
            Symbol::EXPONENTIAL_ICDF,
            statistics::distributions::continuous::exponential_icdf,
        ),
        (
            Symbol::UNIFORM_PDF,
            statistics::distributions::continuous::uniform_pdf,
        ),
        (
            Symbol::UNIFORM_CDF,
            statistics::distributions::continuous::uniform_cdf,
        ),
        (
            Symbol::UNIFORM_ICDF,
            statistics::distributions::continuous::uniform_icdf,
        ),
        (
            Symbol::STUDENT_T_PDF,
            statistics::distributions::degrees::student_pdf,
        ),
        (
            Symbol::STUDENT_T_CDF,
            statistics::distributions::degrees::student_cdf,
        ),
        (
            Symbol::STUDENT_T_ICDF,
            statistics::distributions::degrees::student_icdf,
        ),
        (
            Symbol::CHI_SQUARE_PDF,
            statistics::distributions::degrees::chi_square_pdf,
        ),
        (
            Symbol::CHI_SQUARE_CDF,
            statistics::distributions::degrees::chi_square_cdf,
        ),
        (
            Symbol::CHI_SQUARE_ICDF,
            statistics::distributions::degrees::chi_square_icdf,
        ),
        (
            Symbol::FISHER_PDF,
            statistics::distributions::degrees::fisher_pdf,
        ),
        (
            Symbol::FISHER_CDF,
            statistics::distributions::degrees::fisher_cdf,
        ),
        (
            Symbol::FISHER_ICDF,
            statistics::distributions::degrees::fisher_icdf,
        ),
        (Symbol::SEED_RANDOM, random::seed_random),
        (Symbol::RANDOM, random::random),
        (Symbol::RANDOM_INTEGER, random::random_integer),
        (Symbol::RANDOM_NORMAL, random::variates::random_normal),
        (
            Symbol::RANDOM_EXPONENTIAL,
            random::variates::random_exponential,
        ),
        (Symbol::RANDOM_BINOMIAL, random::variates::random_binomial),
        (Symbol::RANDOM_POISSON, random::variates::random_poisson),
        (Symbol::RANDOM_GEOMETRIC, random::variates::random_geometric),
        (
            Symbol::RANDOM_CHI_SQUARE,
            random::variates::random_chi_square,
        ),
        (Symbol::RANDOM_STUDENT_T, random::variates::random_student_t),
        (Symbol::RANDOM_FISHER, random::variates::random_fisher),
        (
            Symbol::GAMMA_PDF,
            statistics::distributions::shapes::gamma_pdf,
        ),
        (
            Symbol::GAMMA_CDF,
            statistics::distributions::shapes::gamma_cdf,
        ),
        (
            Symbol::GAMMA_ICDF,
            statistics::distributions::shapes::gamma_icdf,
        ),
        (
            Symbol::BETA_PDF,
            statistics::distributions::shapes::beta_pdf,
        ),
        (
            Symbol::BETA_CDF,
            statistics::distributions::shapes::beta_cdf,
        ),
        (
            Symbol::BETA_ICDF,
            statistics::distributions::shapes::beta_icdf,
        ),
        (
            Symbol::WEIBULL_PDF,
            statistics::distributions::shapes::weibull_pdf,
        ),
        (
            Symbol::WEIBULL_CDF,
            statistics::distributions::shapes::weibull_cdf,
        ),
        (
            Symbol::WEIBULL_ICDF,
            statistics::distributions::shapes::weibull_icdf,
        ),
        (
            Symbol::CAUCHY_PDF,
            statistics::distributions::shapes::cauchy_pdf,
        ),
        (
            Symbol::CAUCHY_CDF,
            statistics::distributions::shapes::cauchy_cdf,
        ),
        (
            Symbol::CAUCHY_ICDF,
            statistics::distributions::shapes::cauchy_icdf,
        ),
        (Symbol::RANDOM_GAMMA, random::variates::random_gamma),
        (Symbol::RANDOM_BETA, random::variates::random_beta),
        (Symbol::RANDOM_WEIBULL, random::variates::random_weibull),
        (Symbol::RANDOM_VECTOR, random::collections::random_vector),
        (Symbol::RANDOM_MATRIX, random::collections::random_matrix),
        (
            Symbol::RANDOM_POLYNOMIAL,
            random::collections::random_polynomial,
        ),
        (
            Symbol::RANDOM_PERMUTATION,
            random::collections::random_permutation,
        ),
        (Symbol::RANDOM_SAMPLE, random::collections::random_sample),
        (Symbol::SOLVE, solve::solve),
        (Symbol::CSOLVE, solve::complex_solve),
        (Symbol::FSOLVE, solve::fsolve),
        (Symbol::CFSOLVE, solve::fsolve::complex::cfsolve),
        (Symbol::TRANSPOSE, matrix::transpose),
        (Symbol::IDENTITY_MATRIX, matrix::identity_matrix),
        (Symbol::TR, matrix::tr),
        (Symbol::DOT, matrix::dot),
        (Symbol::CROSS, matrix::cross),
        (Symbol::DET, matrix::det),
        (Symbol::INVERSE, matrix::inverse),
        (Symbol::ROW_REDUCE, matrix::row_reduce),
        (Symbol::MATRIX_RANK, matrix::matrix_rank),
        (Symbol::NULL_SPACE, matrix::null_space),
        (Symbol::COLUMN_SPACE, matrix::column_space),
        (Symbol::LIN_SOLVE, matrix::lin_solve),
        (
            Symbol::CHARACTERISTIC_POLYNOMIAL,
            matrix::characteristic_polynomial,
        ),
        (Symbol::EIGENVALUES, matrix::eigenvalues),
        (Symbol::EIGENVECTORS, matrix::eigenvectors),
        (Symbol::MATRIX, matrix::constructors::matrix),
        (Symbol::MAKE_MATRIX, matrix::constructors::make_matrix),
        (Symbol::LIST_TO_MATRIX, matrix::constructors::list_to_matrix),
        (Symbol::MATRIX_TO_LIST, matrix::constructors::matrix_to_list),
        (Symbol::BLOCK_MATRIX, matrix::constructors::block_matrix),
        (Symbol::JORDAN_BLOCK, matrix::constructors::jordan_block),
        (
            Symbol::BLOCK_DIAGONAL_MATRIX,
            matrix::constructors::block_diagonal_matrix,
        ),
        (Symbol::ROW_DIMENSION, matrix::slicing::row_dimension),
        (Symbol::COLUMN_DIMENSION, matrix::slicing::column_dimension),
        (Symbol::MATRIX_ROW, matrix::slicing::matrix_row),
        (Symbol::MATRIX_COLUMN, matrix::slicing::matrix_column),
        (Symbol::SUB_MATRIX, matrix::slicing::sub_matrix),
        (Symbol::DELETE_ROWS, matrix::slicing::delete_rows),
        (Symbol::DELETE_COLUMNS, matrix::slicing::delete_columns),
        (Symbol::MULTIPLY_ROW, matrix::slicing::multiply_row),
        (Symbol::MULTIPLY_ROW_ADD, matrix::slicing::multiply_row_add),
        (Symbol::SCALE_ROW_ADD, matrix::slicing::scale_row_add),
        (Symbol::SCALE_ROW, matrix::slicing::scale_row),
        (Symbol::ROW_ADD, matrix::slicing::row_add),
        (Symbol::ROW_SWAP, matrix::slicing::row_swap),
        (Symbol::COLUMN_SWAP, matrix::slicing::column_swap),
        (Symbol::ROW_SPACE, matrix::row_space),
        (Symbol::COL_SPACE, matrix::col_space),
        (Symbol::INTERSECTION_BASIS, matrix::intersection_basis),
        (Symbol::SYSTEM_MATRIX, matrix::system_matrix),
        (Symbol::DET_MINOR, matrix::echelon::det_minor),
        (Symbol::ADJOINT_MATRIX, matrix::echelon::adjoint_matrix),
        (Symbol::ROW_ECHELON, matrix::echelon::row_echelon),
        (Symbol::PIVOT, matrix::echelon::pivot),
        (Symbol::JORDAN_MATRIX, matrix::echelon::jordan_matrix),
        (
            Symbol::GAUSS_SQUARES,
            matrix::decompose::forms::gauss_squares,
        ),
        (
            Symbol::RATIONAL_JORDAN_FORM,
            matrix::decompose::rational_jordan_form,
        ),
        (Symbol::SIMPLIFY, simplify::simplify),
        (Symbol::TRIG_EXPAND, trig::trig_expand),
        (Symbol::TRIG_REDUCE, trig::trig_reduce),
        (Symbol::TRIG_TO_EXP, trig::trig_to_exp),
        (Symbol::EXP_TO_TRIG, trig::exp_to_trig),
        (Symbol::HYPERBOLIC_TO_EXP, rewrite::hyperbolic_to_exp),
        (Symbol::EXP_TO_POWER, rewrite::exp_to_power),
        (Symbol::POWER_TO_EXP, rewrite::power_to_exp),
        (Symbol::EXP_EXPAND, rewrite::exp_expand),
        (Symbol::LOG_EXPAND, rewrite::log_expand),
        (Symbol::LOG_COMBINE, rewrite::log_combine),
        (Symbol::HALF_TAN, rewrite::half_tan),
        (
            Symbol::HALF_TAN_HYPERBOLIC_TO_EXP,
            rewrite::half_tan_hyperbolic_to_exp,
        ),
        (Symbol::TAN_TO_SIN_COS, rewrite::tan_to_sin_cos),
        (Symbol::SIN_TO_COS_TAN, rewrite::sin_to_cos_tan),
        (Symbol::COS_TO_SIN_TAN, rewrite::cos_to_sin_tan),
        (Symbol::TAN_TO_SIN_COS_2, rewrite::tan_to_sin_cos_2),
        (Symbol::TAN_TO_COS_SIN_2, rewrite::tan_to_cos_sin_2),
        (Symbol::TRIG_TO_COS, rewrite::trig_to_cos),
        (Symbol::TRIG_TO_SIN, rewrite::trig_to_sin),
        (Symbol::TRIG_TO_TAN, rewrite::trig_to_tan),
        (Symbol::TRIG_COLLECT, rewrite::trig_collect),
        (Symbol::ARC_COS_TO_ARC_SIN, rewrite::arc_cos_to_arc_sin),
        (Symbol::ARC_SIN_TO_ARC_COS, rewrite::arc_sin_to_arc_cos),
        (Symbol::ARC_COS_TO_ARC_TAN, rewrite::arc_cos_to_arc_tan),
        (Symbol::ARC_SIN_TO_ARC_TAN, rewrite::arc_sin_to_arc_tan),
        (Symbol::ARC_TAN_TO_ARC_COS, rewrite::arc_tan_to_arc_cos),
        (Symbol::ARC_TAN_TO_ARC_SIN, rewrite::arc_tan_to_arc_sin),
        (Symbol::ARC_TRIG_TO_LOG, rewrite::arc_trig_to_log),
        (Symbol::COMPLEX_EXPAND, rewrite::complex_expand),
        (Symbol::MULTIPLY_CONJUGATE, rewrite::multiply_conjugate),
        (
            Symbol::MULTIPLY_COMPLEX_CONJUGATE,
            rewrite::multiply_complex_conjugate,
        ),
        (Symbol::COMBINE, rewrite::combine),
        (
            Symbol::TRANSCENDENTAL_SIMPLIFY,
            rewrite::transcendental_simplify,
        ),
        (Symbol::SIGN_TO_HEAVISIDE, rewrite::sign_to_heaviside),
        (Symbol::HEAVISIDE_TO_SIGN, rewrite::heaviside_to_sign),
        (Symbol::LINEAR_STEP, step::linear_step),
        (Symbol::SIMPLIFY_DIRAC, step::simplify_dirac),
        (Symbol::PERIODIC, piecewise::periodic),
        (Symbol::INTEGRATE_BY_PARTS_DV, calculus_tools::by_parts_dv),
        (Symbol::INTEGRATE_BY_PARTS_U, calculus_tools::by_parts_u),
        (Symbol::PREVAL, calculus_tools::preval),
        (Symbol::FUNCTION_DIFF, calculus_tools::function_diff),
        (Symbol::REVERT, series::revert),
        (Symbol::SUM_RIEMANN, sum::sum_riemann),
        (Symbol::WZ_CERTIFICATE, sum::wz_certificate),
        (Symbol::GRAD, vector::grad),
        (Symbol::DIVERGENCE, vector::divergence),
        (Symbol::CURL, vector::curl),
        (Symbol::LAPLACIAN, vector::laplacian),
        (Symbol::HESSIAN, vector::hessian),
        (Symbol::POTENTIAL, vector::potential::potential),
        (
            Symbol::VECTOR_POTENTIAL,
            vector::potential::vector_potential,
        ),
        (Symbol::CURVATURE, vector::curves::curvature),
        (Symbol::FRENET, vector::curves::frenet),
        (Symbol::EVOLUTE, vector::curves::evolute),
        (Symbol::ARC_LENGTH, vector::curves::arc_length),
        (Symbol::IMPLICIT_D, vector::implicit::implicit_d),
        (Symbol::EXTREMA, optimization::extrema),
        (Symbol::MINIMIZE, optimization::global::minimize),
        (Symbol::MAXIMIZE, optimization::global::maximize),
        (Symbol::F_MIN, optimization::global::f_min),
        (Symbol::F_MAX, optimization::global::f_max),
        (Symbol::LP_SOLVE, optimization::lpsolve::lp_solve),
        (Symbol::NLP_SOLVE, optimization::nlpsolve::nlp_solve),
        (
            Symbol::QUOTIENT_REMAINDER,
            number_theory::residues::quotient_remainder,
        ),
        (
            Symbol::SYMMETRIC_MOD,
            number_theory::residues::symmetric_mod,
        ),
        (
            Symbol::PSEUDOPRIME_Q,
            number_theory::residues::pseudoprime_q,
        ),
        (Symbol::CYCLOTOMIC, number_theory::residues::cyclotomic),
        (
            Symbol::MULTIPLICATIVE_ORDER,
            number_theory::residues::multiplicative_order,
        ),
        (
            Symbol::PRIMITIVE_ROOT,
            number_theory::residues::primitive_root,
        ),
        (
            Symbol::CONTINUED_FRACTION,
            number_theory::fractions::continued_fraction,
        ),
        (
            Symbol::FROM_CONTINUED_FRACTION,
            number_theory::fractions::from_continued_fraction,
        ),
        (
            Symbol::INTEGER_SOLVE_LINEAR,
            number_theory::diophantine::integer_solve_linear,
        ),
        (
            Symbol::PRIME_SUM_OF_SQUARES,
            number_theory::diophantine::prime_sum_of_squares,
        ),
        (
            Symbol::INTEGER_SOLVE,
            number_theory::diophantine::integer_solve,
        ),
        (Symbol::SUM, sum::sum),
        (Symbol::PRODUCT, sum::product),
        (Symbol::BINOMIAL, factorial::binomial),
        (Symbol::FACTOR_RATIONAL, polynomial::factor_rational),
        (Symbol::FACTOR_LIST, polynomial::factor_list),
        (Symbol::FACTOR_SQUARE_FREE, polynomial::factor_square_free),
        (Symbol::N, numeric::n),
        (Symbol::REPLACE_ALL, replace::replace_all),
        (Symbol::REPLACE, replace::replace),
        (Symbol::REPLACE_REPEATED, replace::replace_repeated),
        (Symbol::MATCH_Q, replace::match_q),
    ];
    let mut table: [Option<NativeFn>; builtin::COUNT] = [None; builtin::COUNT];
    let mut i = 0;
    while i < entries.len() {
        table[entries[i].0.id() as usize] = Some(entries[i].1);
        i += 1;
    }
    table
};

/// The native function of the builtin `s`, if it has one.
pub(crate) fn get(s: Symbol) -> Option<NativeFn> {
    NATIVES.get(s.id() as usize).copied().flatten()
}

/// The value of a builtin symbol: `I` is the complex number `i` (CAS-51).
pub(crate) fn value(context: &Context, s: Symbol) -> Option<Expr> {
    match s {
        Symbol::I => Some(build::i()),
        // CAS-151.
        Symbol::DIGITS => {
            let digits = crate::print::shown_digits(context.working_precision());
            Some(Expr::from(i64::try_from(digits).unwrap_or(i64::MAX)))
        }
        // CAS-545.
        Symbol::ALL_TRIG_SOLUTIONS => Some(Expr::from(i64::from(context.all_trig_solutions()))),
        Symbol::WITH_SQRT => Some(Expr::from(i64::from(context.with_sqrt()))),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
