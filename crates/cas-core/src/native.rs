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
mod factorial;
mod groebner;
mod index;
mod integer;
mod integrate;
mod laplace;
mod limit;
mod linear;
mod logic;
mod matrix;
mod number;
mod number_theory;
mod numeric;
mod parts;
mod polynomial;
mod procedure;
mod radical;
mod replace;
mod rounding;
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
