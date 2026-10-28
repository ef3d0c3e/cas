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
mod dsolve;
mod elementary;
mod erf;
mod expint;
mod factorial;
mod fourier;
mod groebner;
mod index;
mod integer;
mod integrate;
mod interpolation;
mod laplace;
mod limit;
mod linear;
mod lists;
mod logic;
mod matrix;
mod nintegrate;
mod number;
mod number_theory;
mod numeric;
mod orthogonal;
mod parts;
mod polygamma;
mod polynomial;
mod polytools;
mod procedure;
mod radical;
mod replace;
mod residue;
mod rounding;
mod rsolve;
mod series;
mod simplify;
mod solve;
mod special;
mod statistics;
mod step;
mod structure;
mod sum;
mod trig;
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
        (Symbol::ZETA, zeta::zeta),
        (Symbol::ERF, erf::erf),
        (Symbol::ERFC, erf::erf),
        (Symbol::EXP_INTEGRAL_EI, expint::integral),
        (Symbol::SIN_INTEGRAL, expint::integral),
        (Symbol::COS_INTEGRAL, expint::integral),
        (Symbol::LOG_INTEGRAL, expint::log_integral),
        (Symbol::POLY_GAMMA, polygamma::polygamma),
        (Symbol::HEAVISIDE_THETA, step::heaviside),
        (Symbol::DIRAC_DELTA, step::dirac),
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
        (Symbol::INTEGRATE, integrate::integrate),
        (Symbol::SERIES, series::series),
        (Symbol::LIMIT, limit::limit),
        (Symbol::DSOLVE, dsolve::dsolve),
        (Symbol::LAPLACE_TRANSFORM, laplace::laplace_transform),
        (
            Symbol::INVERSE_LAPLACE_TRANSFORM,
            laplace::inverse_laplace_transform,
        ),
        (Symbol::ZTRANSFORM, ztransform::z_transform),
        (Symbol::INVERSE_ZTRANSFORM, ztransform::inverse_z_transform),
        (Symbol::KRONECKER_DELTA, ztransform::kronecker_delta),
        (Symbol::RSOLVE, rsolve::rsolve),
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
        (Symbol::APART, polytools::apart),
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
        (Symbol::PRIME_Q, number_theory::prime_q),
        (Symbol::NEXT_PRIME, number_theory::next_prime),
        (Symbol::PREVIOUS_PRIME, number_theory::previous_prime),
        (Symbol::INTEGER_FACTOR, number_theory::integer_factor),
        (
            Symbol::INTEGER_FACTOR_LIST,
            number_theory::integer_factor_list,
        ),
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
        (Symbol::SOLVE, solve::solve),
        (Symbol::CSOLVE, solve::complex_solve),
        (Symbol::FSOLVE, solve::fsolve),
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
        (Symbol::SIMPLIFY, simplify::simplify),
        (Symbol::TRIG_EXPAND, trig::trig_expand),
        (Symbol::TRIG_REDUCE, trig::trig_reduce),
        (Symbol::TRIG_TO_EXP, trig::trig_to_exp),
        (Symbol::EXP_TO_TRIG, trig::exp_to_trig),
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
        _ => None,
    }
}

#[cfg(test)]
mod tests;
