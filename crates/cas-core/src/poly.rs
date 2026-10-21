//! Expressions as polynomials and rational functions (CAS-DOC-8 §2,
//! CAS-139): the generators of a list of evaluated expressions, their
//! numerators and denominators over ℤ, or ℤ\[i\] when a number is complex,
//! in cas-poly's sparse form, and the way back through the builder.
//!
//! - **Generators** are what the arithmetic of `Plus`, `Times` and `Power`
//!   to an integer leaves: symbols, constants (`Pi`, `E`), applications
//!   (`Sin[x]`, `f[x]`), and powers that are not polynomial operations.
//! - **Roots** `b^(p/q)` are `b^n (b^(1/q))^r` with `n = trunc(p/q)`, as giac
//!   writes them: `x^(3/2)` is `x · x^(1/2)` and `x^(2/3)` is `(x^(1/3))^2`,
//!   so roots of one order share their generator, `x` stays one of its own,
//!   and roots of different orders stay unrelated (P1). Roots of rationals,
//!   `2^(1/2)`, are generators too, in phase 1.
//! - **Exponentials of one base**, `b^(c u)` with rational `c`, share the
//!   generator `b^(g u)`, `g` the gcd of the `c`: `E^(2x)` with `E^(x/2)` is
//!   `(E^(x/2))^4` (P1). Evaluation adds the exponents of powers of one base
//!   (`E · E^x` is `E^(1 + x)`), so an exponent is taken as an exact
//!   constant and a sum of rational multiples of units, a basis over ℚ
//!   (ours): sums split, as `b^(s + t) = b^s b^t`; a multiple of a sum is
//!   distributed; and a complex multiple `(a + b i) u` is `a u + b (i u)`.
//! - `Exp[c u]` is grouped likewise but never split, as evaluation never
//!   merges `Exp`s and `expand` leaves `exp(a + b)` (P6); `Exp[x]` and `E^x`
//!   stay apart.
//! - A power whose integer exponent does not fit 64 bits stays a generator
//!   of its own.
//! - Infinities and `Indeterminate` are refused: evaluation's arithmetic of
//!   them is not a polynomial ring's.
//! - The generators are in decreasing canonical order (CAS-DOC-3 §4.1):
//!   variable `x₀` of the polynomials is the last in that order, the main
//!   variable of CAS-DOC-8 §15 detail 5, most significant in cas-poly's
//!   lexicographic order.
//!
//! The way back, built then evaluated, keeps the value. Without roots it
//! also converts to the same rational function; with them evaluation brings
//! `(x^(1/2))^2` back to `x`, a generator of its own.
//!
//! [`expand`] reads expressions another way, as giac's `expand` does (P6):
//! only positive powers multiply out, and every other power is a generator,
//! so that each denominator stays with its own terms.

#[cfg(test)]
mod tests;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_num::{Complex, Integer, Natural, Rational};
use cas_poly::{Cancel, Cancelled, Numerators, Sparse, SparseRing};
use foldhash::fast::FixedState;

use crate::{Expr, Symbol, View};

type HashMap<K, V> = hashbrown::HashMap<K, V, FixedState>;
type HashSet<K> = hashbrown::HashSet<K, FixedState>;

/// The numbers of the converted expressions, which decide the coefficient
/// domain: ℤ or ℚ, or ℤ\[i\] or ℚ(i) when one is complex. The numerators and
/// denominators are over ℤ or ℤ\[i\] in every case.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum CoefficientDomain {
    Integers,
    Rationals,
    GaussianIntegers,
    GaussianRationals,
}

impl CoefficientDomain {
    /// Whether a number is complex.
    #[must_use]
    pub fn is_complex(self) -> bool {
        matches!(
            self,
            CoefficientDomain::GaussianIntegers | CoefficientDomain::GaussianRationals
        )
    }
}

/// A rational function: a numerator and a denominator, the denominator's
/// leading coefficient in normal form (positive over ℤ).
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Fraction {
    /// Over ℤ.
    Integers {
        num: Sparse<Integer>,
        den: Sparse<Integer>,
    },
    /// Over ℤ\[i\].
    GaussianIntegers {
        num: Sparse<Complex>,
        den: Sparse<Complex>,
    },
}

/// Why expressions have no polynomial form.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PolyError {
    /// A ball coefficient: the exact algorithms refuse balls (CAS-DOC-8 §9).
    Ball(Expr),
    /// An infinity or `Indeterminate`, anywhere: evaluation's arithmetic of
    /// them (`∞ + x = ∞`) is not a polynomial ring's.
    Infinite(Expr),
    /// A result past the limit on terms, an exponent past `2^62`, or
    /// expressions nested past [`MAX_DEPTH`] polynomial operations.
    TooLarge,
    /// A power whose coefficients would pass the limit on bits (CAS-156).
    Bits,
    /// A denominator that is the zero polynomial, which evaluation did not
    /// see: it never expands.
    DivisionByZero,
    /// The interrupt was raised.
    Interrupted,
}

impl fmt::Display for PolyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolyError::Ball(e) => write!(f, "{e} is not exact"),
            PolyError::Infinite(e) => write!(f, "{e} is not finite"),
            PolyError::TooLarge => f.write_str("too large a polynomial"),
            PolyError::Bits => f.write_str("too large a coefficient"),
            PolyError::DivisionByZero => f.write_str("a division by zero"),
            PolyError::Interrupted => f.write_str("interrupted"),
        }
    }
}

impl std::error::Error for PolyError {}

impl From<Cancelled> for PolyError {
    fn from(Cancelled: Cancelled) -> PolyError {
        PolyError::Interrupted
    }
}

/// The largest exponent a conversion makes, below cas-poly's `2^63`, so that
/// one product of two more cannot reach it.
const MAX_EXPONENT: u64 = 1 << 62;

/// The deepest nesting of polynomial operations converted, within the stack
/// of a test thread.
pub const MAX_DEPTH: usize = 1024;

/// Expressions as rational functions of one list of generators.
#[derive(Clone, Debug)]
pub struct Polynomials {
    generators: Vec<Expr>,
    domain: CoefficientDomain,
    fractions: Vec<Fraction>,
}

/// What a conversion may make (CAS-DOC-8 §15 detail 3): the terms of each
/// product, estimated before it is taken, and the bits of the coefficients
/// of each power, as `Limits::terms` and `Limits::number_bits` bound them
/// in the builtins (CAS-156).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Bounds {
    /// The terms of a product, bounded before it is taken.
    pub terms: usize,
    /// The bits of a power's coefficients, estimated before it is taken.
    pub bits: u64,
}

impl Bounds {
    /// No bounds.
    pub const NONE: Bounds = Bounds {
        terms: usize::MAX,
        bits: u64::MAX,
    };

    /// At most `terms` terms a product, and no bound on bits.
    #[must_use]
    pub const fn terms(terms: usize) -> Bounds {
        Bounds {
            terms,
            bits: u64::MAX,
        }
    }
}

/// What a conversion takes as a generator, beyond the kernels of §2.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Reading {
    /// Nothing more: balls are refused.
    Exact,
    /// Balls.
    Balls,
    /// Balls, every power whose exponent is not an exact positive number,
    /// and every `Exp`: the expressions are polynomials, for [`expand`].
    Terms,
}

impl Polynomials {
    /// The rational functions of evaluated expressions, over their common
    /// generators, within `bounds`; `cancel` is polled in the long
    /// products.
    ///
    /// # Errors
    ///
    /// [`PolyError::Ball`] for a ball coefficient, [`PolyError::Infinite`]
    /// for an infinity or `Indeterminate`, [`PolyError::TooLarge`] past
    /// the bound on terms, an exponent of `2^62` or [`MAX_DEPTH`],
    /// [`PolyError::Bits`] past the bound on bits,
    /// [`PolyError::DivisionByZero`] for a zero denominator,
    /// [`PolyError::Interrupted`] when `cancel` says so.
    pub fn new(
        exprs: &[Expr],
        bounds: Bounds,
        cancel: &dyn Cancel,
    ) -> Result<Polynomials, PolyError> {
        Polynomials::convert(exprs, bounds, Reading::Exact, cancel)
    }

    /// As [`Polynomials::new`], each ball a generator of its own instead of
    /// an error: for arithmetic, which works with balls (CAS-DOC-8 §9). The
    /// way back multiplies them out in evaluation, so `(x + 1.5)^2` expands
    /// to `x^2 + 3.0 x + 2.25`.
    ///
    /// # Errors
    ///
    /// As [`Polynomials::new`], but for balls.
    pub fn with_balls(
        exprs: &[Expr],
        bounds: Bounds,
        cancel: &dyn Cancel,
    ) -> Result<Polynomials, PolyError> {
        Polynomials::convert(exprs, bounds, Reading::Balls, cancel)
    }

    fn convert(
        exprs: &[Expr],
        bounds: Bounds,
        reading: Reading,
        cancel: &dyn Cancel,
    ) -> Result<Polynomials, PolyError> {
        for e in exprs {
            if let Some(infinite) = infinite(e) {
                return Err(PolyError::Infinite(infinite));
            }
        }
        let kernels = Kernels::scan(exprs, reading)?;
        let generators = kernels.generators();
        let domain = kernels.domain();
        let fractions = if domain.is_complex() {
            convert::<GaussianIntegers>(exprs, &generators, reading, bounds, cancel)?
                .into_iter()
                .map(|(num, den)| Fraction::GaussianIntegers { num, den })
                .collect()
        } else {
            convert::<Integers>(exprs, &generators, reading, bounds, cancel)?
                .into_iter()
                .map(|(num, den)| Fraction::Integers { num, den })
                .collect()
        };
        Ok(Polynomials {
            generators: generators.list,
            domain,
            fractions,
        })
    }

    /// The generators: `generators()[k]` is the variable `x_k`, in
    /// decreasing canonical order.
    #[must_use]
    pub fn generators(&self) -> &[Expr] {
        &self.generators
    }

    /// The coefficient domain, from the numbers of the expressions.
    #[must_use]
    pub fn domain(&self) -> CoefficientDomain {
        self.domain
    }

    /// The variable of the generator `g`, if it is one.
    #[must_use]
    pub fn variable(&self, g: &Expr) -> Option<usize> {
        self.generators.iter().position(|h| h == g)
    }

    /// The rational functions, one an expression, in their order.
    #[must_use]
    pub fn fractions(&self) -> &[Fraction] {
        &self.fractions
    }

    /// The expression of the polynomial `p` in the generators, built, not
    /// evaluated: a sum of products of a coefficient and powers of
    /// generators.
    #[must_use]
    pub fn polynomial_expr<E: Coefficient>(&self, p: &Sparse<E>) -> Expr {
        self.terms_expr(p, None)
    }

    /// [`Polynomials::polynomial_expr`] of `p / d` for a constant `d`, each
    /// coefficient divided.
    fn terms_expr<E: Coefficient>(&self, p: &Sparse<E>, d: Option<&E>) -> Expr {
        sum(self.term_exprs(p, d))
    }

    /// The terms of `p / d`, built, for a constant `d`.
    fn term_exprs<E: Coefficient>(&self, p: &Sparse<E>, d: Option<&E>) -> Vec<Expr> {
        p.terms()
            .map(|(exps, c)| self.term_expr(&exps, c, d))
            .collect()
    }

    /// The term `c / d` times the generators to the powers `exps`, built.
    fn term_expr<E: Coefficient>(&self, exps: &[u64], c: &E, d: Option<&E>) -> Expr {
        let c = match d {
            Some(d) => c.over(d),
            None => c.to_expr(),
        };
        let mut factors = Vec::new();
        if c != Expr::from(1) {
            factors.push(c);
        }
        for (generator, &e) in self.generators.iter().zip(exps) {
            match e {
                0 => {}
                1 => factors.push(generator.clone()),
                _ => factors.push(generator.pow(Expr::from(Integer::from(e)))),
            }
        }
        product(factors)
    }

    /// The expression of a rational function, built, not evaluated: over a
    /// constant denominator, a polynomial with its coefficients divided;
    /// otherwise the numerator times the denominator to the power -1.
    #[must_use]
    pub fn fraction_expr(&self, fraction: &Fraction) -> Expr {
        match fraction {
            Fraction::Integers { num, den } => self.ratio_expr(num, den),
            Fraction::GaussianIntegers { num, den } => self.ratio_expr(num, den),
        }
    }

    /// [`Polynomials::fraction_expr`] of `num / den`.
    #[must_use]
    pub fn ratio_expr<E: Coefficient>(&self, num: &Sparse<E>, den: &Sparse<E>) -> Expr {
        if is_constant(den) {
            let d = &den.coeffs()[0];
            return self.terms_expr(num, (!d.is_one()).then_some(d));
        }
        Expr::call(
            Symbol::TIMES,
            [self.polynomial_expr(num), self.polynomial_expr(den).pow(-1)],
        )
    }

    /// `num / den` as `expand` writes it (CAS-DOC-8 P6): over a constant, a
    /// polynomial with its coefficients divided; otherwise each term of the
    /// numerator over the denominator, as giac writes
    /// `x^2/(x^2-2*x+1) + 2*x/(x^2-2*x+1) + 1/(x^2-2*x+1)`, so that over a
    /// monomial evaluation merges the powers (`x/y + 1/y`). Built, not
    /// evaluated.
    #[must_use]
    pub fn expanded_expr<E: Coefficient>(&self, num: &Sparse<E>, den: &Sparse<E>) -> Expr {
        if is_constant(den) {
            return self.ratio_expr(num, den);
        }
        let inverse = self.polynomial_expr(den).pow(-1);
        sum(self
            .term_exprs(num, None)
            .into_iter()
            .map(|t| Expr::call(Symbol::TIMES, [t, inverse.clone()]))
            .collect())
    }
}

/// The sum of `terms`: 0 for none, the term for one.
fn sum(mut terms: Vec<Expr>) -> Expr {
    match terms.len() {
        0 => Expr::from(0),
        1 => terms.pop().expect("one term"),
        _ => Expr::call(Symbol::PLUS, terms),
    }
}

/// The product of `factors`: 1 for none, the factor for one.
fn product(mut factors: Vec<Expr>) -> Expr {
    match factors.len() {
        0 => Expr::from(1),
        1 => factors.pop().expect("one factor"),
        _ => Expr::call(Symbol::TIMES, factors),
    }
}

/// `e` expanded as giac's `expand` writes it (CAS-DOC-8 P6), built, not
/// evaluated.
///
/// - Products and positive powers multiply out, of sums inside products
///   too: `(1 + 1/(x + 1)) (x + 2)` is `x + 2 + x/(x + 1) + 2/(x + 1)`.
/// - A power with a negative exponent is the denominator of its terms, its
///   base expanded, and the denominators of one term multiply into one, so
///   that `(x + 1)/(x (x + 2))` is `x/(x^2 + 2 x) + 1/(x^2 + 2 x)`.
/// - Terms are never brought over a common denominator: `1/x + 1/(x + 1)`
///   stays.
/// - Numbers stay coefficients: `(x + 1)^2/(2 y)` is `x^2/(2 y) + x/y +
///   1/(2 y)`. Balls are generators.
/// - Other powers and `Exp` are generators, `exp(a + b)` among them.
///
/// giac also writes `1/(x + 1)/(x - 1)` as `1/((x + 1) (x - 1))`, where
/// `1/((x + 1) (x - 1))` is `1/(x^2 - 1)`; evaluation makes them one
/// expression, which is the second.
///
/// # Errors
///
/// As [`Polynomials::with_balls`].
pub fn expand(e: &Expr, bounds: Bounds, cancel: &dyn Cancel) -> Result<Expr, PolyError> {
    let p = Polynomials::convert(std::slice::from_ref(e), bounds, Reading::Terms, cancel)?;
    match &p.fractions[0] {
        Fraction::Integers { num, den } => p.expansion(num, den, bounds, cancel),
        Fraction::GaussianIntegers { num, den } => p.expansion(num, den, bounds, cancel),
    }
}

/// The base `b` and the exponent `-r` of a generator `b^r` with a negative
/// exact `r`, which [`expand`] takes as a denominator.
fn inverse(g: &Expr) -> Option<(Expr, Rational)> {
    if let View::Apply { head, args: [b, r] } = g.view()
        && head.as_symbol() == Some(Symbol::POWER)
        && let Some(r) = r.to_rational()
        && r.is_negative()
    {
        return Some((b.clone(), -r));
    }
    None
}

impl Polynomials {
    /// [`expand`] of `num / den` for a constant `den`, read as
    /// [`Reading::Terms`] reads: the terms with the same powers of the
    /// inverses `b^(-r)` go over the expansion of their product of `b^r`.
    fn expansion<E: Coefficient>(
        &self,
        num: &Sparse<E>,
        den: &Sparse<E>,
        bounds: Bounds,
        cancel: &dyn Cancel,
    ) -> Result<Expr, PolyError> {
        debug_assert!(is_constant(den), "only numbers have denominators");
        let d = &den.coeffs()[0];
        let d = (!d.is_one()).then_some(d);
        let inverses: Vec<Option<(Expr, Rational)>> = self.generators.iter().map(inverse).collect();
        let mut groups: BTreeMap<Vec<u64>, Vec<Expr>> = BTreeMap::new();
        for (mut exps, c) in num.terms() {
            let mut powers = vec![0; exps.len()];
            for ((e, power), inverse) in exps.iter_mut().zip(&mut powers).zip(&inverses) {
                if inverse.is_some() {
                    *power = std::mem::take(e);
                }
            }
            groups
                .entry(powers)
                .or_default()
                .push(self.term_expr(&exps, c, d));
        }
        let mut terms = Vec::with_capacity(num.term_count());
        for (powers, group) in groups {
            if powers.iter().all(|&k| k == 0) {
                terms.extend(group);
                continue;
            }
            let factors = powers
                .iter()
                .zip(&inverses)
                .filter(|&(&k, _)| k > 0)
                .map(|(&k, inverse)| {
                    let (b, r) = inverse.as_ref().expect("a power of an inverse");
                    let r = r * &Rational::from(Integer::from(k));
                    if r.is_one() {
                        b.clone()
                    } else {
                        b.pow(Expr::from(r))
                    }
                })
                .collect();
            let over = expand(&product(factors), bounds, cancel)?.pow(-1);
            terms.extend(group.into_iter().map(|t| {
                if t == Expr::from(1) {
                    over.clone()
                } else {
                    Expr::call(Symbol::TIMES, [t, over.clone()])
                }
            }));
        }
        Ok(sum(terms))
    }
}

/// The coefficients of the polynomials of a conversion: integers and
/// Gaussian integers.
pub trait Coefficient: Clone {
    /// The number as an expression.
    fn to_expr(&self) -> Expr;

    /// `self / d` as an expression, for a non-zero `d`.
    fn over(&self, d: &Self) -> Expr;

    fn is_one(&self) -> bool;
}

impl Coefficient for Integer {
    fn to_expr(&self) -> Expr {
        Expr::from(self)
    }

    fn over(&self, d: &Integer) -> Expr {
        Expr::from(Rational::new(self.clone(), d.clone()))
    }

    fn is_one(&self) -> bool {
        Integer::is_one(self)
    }
}

impl Coefficient for Complex {
    fn to_expr(&self) -> Expr {
        Expr::from(self)
    }

    fn over(&self, d: &Complex) -> Expr {
        Expr::from(self / d)
    }

    fn is_one(&self) -> bool {
        Complex::is_one(self)
    }
}

/// The first infinity or `Indeterminate` that `e` contains.
fn infinite(e: &Expr) -> Option<Expr> {
    [
        Symbol::INFINITY,
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
        Symbol::INDETERMINATE,
    ]
    .into_iter()
    .map(Expr::from)
    .find(|s| e.contains(s))
}

/// Whether `p` is a constant.
fn is_constant<E>(p: &Sparse<E>) -> bool {
    p.degrees().iter().all(|&d| d == 0)
}

/// How an exponential is written: `b^e` or `Exp[e]`, which evaluation keeps
/// apart.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
enum Form {
    Power,
    Exp,
}

/// An exponential's group: its form, its base (`E` for `Exp`), and the
/// unit `u` of its exponent `c u`.
type Group = (Form, Expr, Expr);

/// What the arithmetic of the expressions leaves, found by one scan.
struct Kernels {
    /// Generators as they are.
    plain: BTreeSet<Expr>,
    /// Roots `b^(1/q)`: the base and the order.
    roots: BTreeSet<(Expr, u64)>,
    /// Exponentials: the rational multiples `c` of each group's unit.
    exponentials: BTreeMap<Group, Vec<Rational>>,
    /// Whether a number is not an integer, and whether one is complex.
    fractional: bool,
    complex: bool,
}

/// What a power is to the conversion.
enum Power {
    /// `b^n` for an integer `n` of 64 bits: arithmetic on `b`.
    Integer(i64),
    /// `b^(p/q)`: `b^n (b^(1/q))^r`.
    Root { n: i64, r: i64, q: u64 },
    /// `b^(c u)` for a rational `c` and a non-numeric `u`.
    Exponential(Rational, Expr),
    /// Anything else: a generator, this one.
    Plain(Expr),
}

/// `p / q` truncated and its remainder, for `q > 0`: `p = n q + r` with `r`
/// of the sign of `p`.
fn truncated(p: i64, q: i64) -> (i64, i64) {
    (p / q, p % q)
}

/// `(r, z₀)` with `z = r z₀`: `r` rational, `z₀` a Gaussian integer with
/// coprime parts whose first non-zero part is positive, 1 for a real `z`.
/// So `-2i x` and `-i x` are rational multiples of one unit, `i x`.
fn rational_part(z: &Complex) -> (Rational, Complex) {
    if z.is_real() {
        return (z.re().clone(), Complex::from(1));
    }
    let (re, im) = (z.re(), z.im());
    let mut r = Rational::new(re.numer().gcd(im.numer()), re.denom().lcm(im.denom()));
    let negative = if re.is_zero() {
        im.is_negative()
    } else {
        re.is_negative()
    };
    if negative {
        r = -r;
    }
    let z0 = z / &Complex::from(r.clone());
    (r, z0)
}

/// The product of the exact number `c` and the factors `rest`, the number
/// first unless it is 1.
fn times(c: &Complex, rest: &[Expr]) -> Expr {
    let mut factors = Vec::with_capacity(rest.len() + 1);
    if !c.is_one() {
        factors.push(Expr::from(c.clone()));
    }
    factors.extend(rest.iter().cloned());
    match factors.len() {
        1 => factors.pop().expect("one factor"),
        _ => Expr::call(Symbol::TIMES, factors),
    }
}

/// The leading exact number of a product, and its other factors.
fn coefficient(e: &Expr) -> (Complex, Vec<Expr>) {
    if let View::Apply { head, args } = e.view()
        && head.as_symbol() == Some(Symbol::TIMES)
        && let [first, rest @ ..] = args
        && let Some(c) = first.to_complex()
    {
        return (c, rest.to_vec());
    }
    (Complex::from(1), vec![e.clone()])
}

/// The rational multiple `c` and the unit `u` of an exponent `c u`: the
/// leading number of a product, its rational part.
fn split_exponent(e: &Expr) -> (Rational, Expr) {
    let (z, rest) = coefficient(e);
    let (r, z0) = rational_part(&z);
    (r, times(&z0, &rest))
}

/// `g u`, with `g` folded into `u`'s leading number, or into `u` itself
/// when it is one (`2 I` is the number `2 I`).
fn scaled(g: &Rational, unit: &Expr) -> Expr {
    let g = Complex::from(g.clone());
    if let Some(z) = unit.to_complex() {
        return Expr::from(&g * &z);
    }
    let (z, rest) = coefficient(unit);
    times(&(&g * &z), &rest)
}

/// An exponent as an exact constant and a sum of rational multiples of
/// units, a multiple of a sum distributed: `-(x - 1)` is `1 - x`. A complex
/// multiple `(a + b i) u` is `a u + b (i u)`, so the units are a basis over
/// ℚ: evaluation adds exponents, and `(1 + i) z + (1 - i) z` is `2 z`.
fn linear_terms(e: &Expr) -> (Complex, Vec<(Rational, Expr)>) {
    let mut constant = Complex::from(0);
    // The non-numeric factors of each term, and their complex multiple.
    let mut terms: BTreeMap<Expr, Complex> = BTreeMap::new();
    let mut work: Vec<(Complex, Expr)> = vec![(Complex::from(1), e.clone())];
    while let Some((c, e)) = work.pop() {
        if let Some(z) = e.to_complex() {
            constant = &constant + &(&c * &z);
            continue;
        }
        if e.head().and_then(Expr::as_symbol) == Some(Symbol::PLUS) {
            work.extend(e.args().iter().map(|t| (c.clone(), t.clone())));
            continue;
        }
        let (k, rest) = coefficient(&e);
        let c = &c * &k;
        if let [sum] = &rest[..]
            && sum.head().and_then(Expr::as_symbol) == Some(Symbol::PLUS)
        {
            work.push((c, sum.clone()));
            continue;
        }
        let entry = terms
            .entry(times(&Complex::from(1), &rest))
            .or_insert_with(|| Complex::from(0));
        *entry = &*entry + &c;
    }
    let mut units: BTreeMap<Expr, Rational> = BTreeMap::new();
    for (rest, c) in terms {
        if !c.re().is_zero() {
            units.insert(rest.clone(), c.re().clone());
        }
        if !c.im().is_zero() {
            let (_, factors) = coefficient(&rest);
            units.insert(times(&Complex::I, &factors), c.im().clone());
        }
    }
    (constant, units.into_iter().map(|(u, c)| (c, u)).collect())
}

/// What `Power[b, e]` is to the conversion: its factors, one a term of the
/// exponent's [`linear_terms`], as `b^(s + t) = b^s b^t`. Evaluation merges
/// `b^s b^t` into `b^(s + t)`, so this keeps a round trip's generators
/// (ours, CAS-139). `Exp` is not split: evaluation keeps `Exp[a] Exp[b]`
/// apart, and `expand` leaves `exp(a + b)` (CAS-DOC-8 P6).
fn powers(b: &Expr, e: &Expr) -> Vec<Power> {
    // A number to an integer power, or to a rational one past 1, is what
    // evaluation leaves only past `Limits::number_bits`: a generator, not a
    // number to compute (CAS-156), as `3^4782969` at a point `Trace::verify`
    // draws for `(3 3^(1/3))^(x^7)` took a minute.
    if b.is_number()
        && e.to_rational()
            .is_some_and(|r| r.is_integer() || r.abs() > Rational::from(1))
    {
        return vec![Power::Plain(b.pow(e.clone()))];
    }
    let (constant, terms) = match e.to_complex() {
        Some(z) => (z, Vec::new()),
        // A ball.
        None if e.is_number() => return vec![Power::Plain(b.pow(e.clone()))],
        None => linear_terms(e),
    };
    let mut factors = Vec::with_capacity(terms.len() + 2);
    // b^(r + s i) = b^r (b^i)^s: the imaginary unit is a unit like x.
    if !constant.re().is_zero() {
        factors.push(power(b, &Expr::from(constant.re().clone())));
    }
    if !constant.im().is_zero() {
        factors.push(Power::Exponential(
            constant.im().clone(),
            Expr::from(Complex::I),
        ));
    }
    for (c, unit) in terms {
        factors.push(if unit.is_number() {
            // A ball: a generator of its own.
            Power::Plain(b.pow(Expr::call(Symbol::TIMES, [Expr::from(c), unit])))
        } else {
            Power::Exponential(c, unit)
        });
    }
    factors
}

/// Whether `e` is an exact positive number, an exponent that [`expand`]
/// multiplies out.
fn positive(e: &Expr) -> bool {
    e.to_rational().is_some_and(|r| r.is_positive())
}

/// What `b^e` is to the conversion, for an `e` that is not a sum.
fn power(b: &Expr, e: &Expr) -> Power {
    if let Some(n) = e.as_i64() {
        return Power::Integer(n);
    }
    match e.view() {
        View::Rational(ratio) => match (ratio.numer().to_i64(), ratio.denom().to_i64()) {
            (Some(num), Some(den)) => {
                let (n, rem) = truncated(num, den);
                Power::Root {
                    n,
                    r: rem,
                    q: den.unsigned_abs(),
                }
            }
            _ => Power::Plain(b.pow(e.clone())),
        },
        _ if e.is_number() => Power::Plain(b.pow(e.clone())),
        _ => {
            let (c, u) = split_exponent(e);
            Power::Exponential(c, u)
        }
    }
}

impl Kernels {
    /// Scans the polynomial operations of `exprs`, every distinct
    /// subexpression once.
    fn scan(exprs: &[Expr], reading: Reading) -> Result<Kernels, PolyError> {
        let mut k = Kernels {
            plain: BTreeSet::new(),
            roots: BTreeSet::new(),
            exponentials: BTreeMap::new(),
            fractional: false,
            complex: false,
        };
        let mut seen: HashSet<Expr> = HashSet::default();
        let mut stack: Vec<Expr> = exprs.to_vec();
        while let Some(e) = stack.pop() {
            if !seen.insert(e.clone()) {
                continue;
            }
            match e.view() {
                View::Int(_) | View::Integer(_) => {}
                View::Rational(_) => k.fractional = true,
                View::Complex(z) => {
                    k.complex = true;
                    k.fractional |= !z.is_gaussian_integer();
                }
                View::Ball(_) | View::ComplexBall(_) => {
                    if reading == Reading::Exact {
                        return Err(PolyError::Ball(e.clone()));
                    }
                    k.plain.insert(e.clone());
                }
                View::Apply { head, args } => match (head.as_symbol(), args) {
                    (Some(Symbol::PLUS | Symbol::TIMES), _) => stack.extend(args.iter().cloned()),
                    (Some(Symbol::POWER), [_, exponent])
                        if reading == Reading::Terms && !positive(exponent) =>
                    {
                        k.plain.insert(e.clone());
                    }
                    (Some(Symbol::POWER), [b, exponent]) => {
                        for p in powers(b, exponent) {
                            match p {
                                Power::Integer(_) => stack.push(b.clone()),
                                Power::Root { n, q, .. } => {
                                    k.roots.insert((b.clone(), q));
                                    if n != 0 {
                                        stack.push(b.clone());
                                    }
                                }
                                Power::Exponential(c, u) => {
                                    k.exponentials
                                        .entry((Form::Power, b.clone(), u))
                                        .or_default()
                                        .push(c);
                                }
                                Power::Plain(g) => {
                                    k.plain.insert(g);
                                }
                            }
                        }
                    }
                    (Some(Symbol::EXP), [x]) if reading != Reading::Terms && !x.is_number() => {
                        let (c, u) = split_exponent(x);
                        k.exponentials
                            .entry((Form::Exp, Expr::from(Symbol::E), u))
                            .or_default()
                            .push(c);
                    }
                    _ => {
                        k.plain.insert(e.clone());
                    }
                },
                _ => {
                    k.plain.insert(e.clone());
                }
            }
        }
        Ok(k)
    }

    fn domain(&self) -> CoefficientDomain {
        match (self.complex, self.fractional) {
            (false, false) => CoefficientDomain::Integers,
            (false, true) => CoefficientDomain::Rationals,
            (true, false) => CoefficientDomain::GaussianIntegers,
            (true, true) => CoefficientDomain::GaussianRationals,
        }
    }

    /// The generators, in decreasing canonical order, and where each
    /// kernel goes.
    fn generators(&self) -> Generators {
        let mut all: BTreeSet<Expr> = self.plain.clone();
        let mut roots = HashMap::default();
        for (b, q) in &self.roots {
            let root = b.pow(Expr::from(Rational::new(Integer::ONE, Integer::from(*q))));
            roots.insert((b.clone(), *q), root.clone());
            all.insert(root);
        }
        let mut exponentials = HashMap::default();
        for (group, multiples) in &self.exponentials {
            let g = multiples
                .iter()
                .fold(Rational::from(0), |g, c| rational_gcd(&g, c));
            let (form, base, unit) = group;
            let exponent = scaled(&g, unit);
            let generator = match form {
                Form::Power => base.pow(exponent),
                Form::Exp => Expr::call(Symbol::EXP, [exponent]),
            };
            exponentials.insert(group.clone(), (generator.clone(), g));
            all.insert(generator);
        }
        let list: Vec<Expr> = all.into_iter().rev().collect();
        let index = list
            .iter()
            .enumerate()
            .map(|(k, e)| (e.clone(), k))
            .collect();
        Generators {
            list,
            index,
            roots,
            exponentials,
        }
    }
}

/// A bound on the terms of `a b`, whose degree in each variable is at most
/// `degrees` (ours): the least of the products of terms; the monomials
/// within those degrees, `∏ (d_i + 1)`; and the monomials of total degree
/// at most `D = tdeg a + tdeg b` in the `n` variables of the product,
/// `C(n + D, n)`. So a dense product is judged by its result, not by its
/// products of terms: `(1 + x + y + z + t)^20` squared makes 113 million
/// products of terms, and 135751 terms. Saturating.
fn product_terms<E>(a: &Sparse<E>, b: &Sparse<E>, degrees: &[u64]) -> usize {
    let pairs = a.term_count().saturating_mul(b.term_count());
    let present: Vec<u64> = degrees.iter().copied().filter(|&d| d > 0).collect();
    let boxed = present.iter().fold(1usize, |acc, &d| {
        acc.saturating_mul(usize::try_from(d).unwrap_or(usize::MAX).saturating_add(1))
    });
    let bound = pairs.min(boxed);
    let total = total_degree(a).saturating_add(total_degree(b));
    // C(D + k, k) for k = 1, …, n, each step exact, until it passes the
    // other bounds.
    let mut simplex: u128 = 1;
    for k in 1..=present.len() as u128 {
        simplex = simplex.saturating_mul(u128::from(total) + k) / k;
        if simplex > bound as u128 {
            return bound;
        }
    }
    usize::try_from(simplex).map_or(bound, |s| s.min(bound))
}

/// The largest total degree of the terms of `p`.
fn total_degree<E>(p: &Sparse<E>) -> u64 {
    p.terms()
        .map(|(exps, _)| exps.iter().fold(0u64, |acc, &e| acc.saturating_add(e)))
        .max()
        .unwrap_or(0)
}

/// The gcd of two rationals, non-negative: the gcd of the numerators over
/// the lcm of the denominators; `gcd(0, c) = |c|`.
fn rational_gcd(a: &Rational, b: &Rational) -> Rational {
    let num = a.numer().gcd(b.numer());
    let den = a.denom().lcm(b.denom());
    Rational::new(num, den)
}

/// The generators and the way to them.
struct Generators {
    /// In decreasing canonical order: `list[k]` is `x_k`.
    list: Vec<Expr>,
    index: HashMap<Expr, usize>,
    /// `(b, q)` to `b^(1/q)`.
    roots: HashMap<(Expr, u64), Expr>,
    /// A group to its generator `b^(g u)` and `g`.
    exponentials: HashMap<Group, (Expr, Rational)>,
}

impl Generators {
    fn var(&self, generator: &Expr) -> usize {
        self.index[generator]
    }
}

/// A numerator and a denominator over `N`.
type Frac<E> = (Sparse<E>, Sparse<E>);

/// The numerator rings of the conversion, ℤ and ℤ\[i\].
trait Numbers: Numerators {
    /// `(n, d)` with `z = n / d`, `d` a positive natural; `z` is real over ℤ.
    fn split_number(z: &Complex) -> (Self::Elem, Natural);

    /// Bits enough for `log2 |c|`, 0 for a unit: what a power of `c`
    /// grows by, each factor.
    fn bits(c: &Self::Elem) -> u64;
}

impl Numbers for Integers {
    fn split_number(z: &Complex) -> (Integer, Natural) {
        debug_assert!(z.is_real(), "complex numbers make the domain ℤ[i]");
        Integers::split(z.re())
    }

    fn bits(c: &Integer) -> u64 {
        match c.bit_length() {
            ..=1 => 0,
            n => n,
        }
    }
}

impl Numbers for GaussianIntegers {
    fn split_number(z: &Complex) -> (Complex, Natural) {
        GaussianIntegers::split(z)
    }

    fn bits(c: &Complex) -> u64 {
        // |a + b i| < 2 max(|a|, |b|); the units are ±1 and ±i.
        let bits = |r: &Rational| r.numer().bit_length();
        match bits(c.re()).max(bits(c.im())) {
            ..=1 if c.re().is_zero() || c.im().is_zero() => 0,
            n => n + 1,
        }
    }
}

/// The rational functions of `exprs` over `N`.
fn convert<N: Numbers>(
    exprs: &[Expr],
    generators: &Generators,
    reading: Reading,
    bounds: Bounds,
    cancel: &dyn Cancel,
) -> Result<Vec<Frac<N::Elem>>, PolyError> {
    let mut converter = Converter::<N> {
        ring: SparseRing::new(N::default(), generators.list.len()),
        generators,
        reading,
        max_terms: bounds.terms,
        max_bits: bounds.bits,
        cancel,
        memo: HashMap::default(),
        depth: 0,
    };
    exprs.iter().map(|e| converter.fraction(e)).collect()
}

/// The conversion of expressions to rational functions over `N`, memoised
/// by subexpression.
struct Converter<'a, N: Numbers> {
    ring: SparseRing<N>,
    generators: &'a Generators,
    reading: Reading,
    max_terms: usize,
    max_bits: u64,
    cancel: &'a dyn Cancel,
    memo: HashMap<Expr, Frac<N::Elem>>,
    depth: usize,
}

impl<N: Numbers> Converter<'_, N> {
    fn fraction(&mut self, e: &Expr) -> Result<Frac<N::Elem>, PolyError> {
        if let Some(f) = self.memo.get(e) {
            return Ok(f.clone());
        }
        if self.depth >= MAX_DEPTH {
            return Err(PolyError::TooLarge);
        }
        if self.cancel.cancelled() {
            return Err(PolyError::Interrupted);
        }
        self.depth += 1;
        let f = self.compute(e);
        self.depth -= 1;
        let f = f?;
        self.memo.insert(e.clone(), f.clone());
        Ok(f)
    }

    fn compute(&mut self, e: &Expr) -> Result<Frac<N::Elem>, PolyError> {
        if let Some(z) = e.to_complex() {
            return Ok(self.number(&z));
        }
        let View::Apply { head, args } = e.view() else {
            return Ok(self.generator(e));
        };
        match (head.as_symbol(), args) {
            (Some(Symbol::PLUS), _) => self.balanced(args, false),
            (Some(Symbol::TIMES), _) => self.balanced(args, true),
            (Some(Symbol::POWER), [_, exponent])
                if self.reading == Reading::Terms && !positive(exponent) =>
            {
                Ok(self.generator(e))
            }
            (Some(Symbol::POWER), [b, exponent]) => {
                let mut product = self.constant(1);
                for p in powers(b, exponent) {
                    let f = self.power(b, p)?;
                    product = self.mul(&product, &f)?;
                }
                Ok(product)
            }
            (Some(Symbol::EXP), [x]) if self.reading != Reading::Terms && !x.is_number() => {
                let (c, u) = split_exponent(x);
                self.exponential(&(Form::Exp, Expr::from(Symbol::E), u), &c)
            }
            _ => Ok(self.generator(e)),
        }
    }

    /// The sum of the fractions of `args`, or their product, in a balanced
    /// tree: a sum of `n` terms takes `O(n log n)` operations on terms,
    /// where adding each to the sum so far took `O(n²)` (CAS-157). The
    /// order changes no result: a sum's denominator is the lcm of the
    /// denominators, a product's their product, normalised either way.
    fn balanced(&mut self, args: &[Expr], product: bool) -> Result<Frac<N::Elem>, PolyError> {
        let mut parts = args
            .iter()
            .map(|a| self.fraction(a))
            .collect::<Result<Vec<_>, _>>()?;
        while parts.len() > 1 {
            let mut next = Vec::with_capacity(parts.len().div_ceil(2));
            let mut pairs = parts.into_iter();
            while let Some(a) = pairs.next() {
                next.push(match pairs.next() {
                    Some(b) if product => self.mul(&a, &b)?,
                    Some(b) => self.add(&a, &b)?,
                    None => a,
                });
            }
            parts = next;
        }
        Ok(parts
            .pop()
            .unwrap_or_else(|| self.constant(i64::from(product))))
    }

    /// One factor of a power of `b`.
    fn power(&mut self, b: &Expr, p: Power) -> Result<Frac<N::Elem>, PolyError> {
        match p {
            Power::Integer(n) => {
                let f = self.fraction(b)?;
                self.pow(&f, n)
            }
            Power::Root { n, r, q } => {
                let root = &self.generators.roots[&(b.clone(), q)];
                let g = self.generator(root);
                let g = self.pow(&g, r)?;
                if n == 0 {
                    return Ok(g);
                }
                let f = self.fraction(b)?;
                let f = self.pow(&f, n)?;
                self.mul(&f, &g)
            }
            Power::Exponential(c, u) => self.exponential(&(Form::Power, b.clone(), u), &c),
            Power::Plain(g) => Ok(self.generator(&g)),
        }
    }

    /// `b^(c u)` as a power of its group's generator `b^(g u)`: `c / g`.
    fn exponential(&mut self, group: &Group, c: &Rational) -> Result<Frac<N::Elem>, PolyError> {
        let (generator, g) = &self.generators.exponentials[group];
        let k = c / g;
        debug_assert!(k.is_integer(), "the gcd divides each multiple");
        let Some(k) = k.numer().to_i64() else {
            return Err(PolyError::TooLarge);
        };
        let f = self.generator(generator);
        self.pow(&f, k)
    }

    fn generator(&self, e: &Expr) -> Frac<N::Elem> {
        (self.ring.var(self.generators.var(e)), self.ring.one())
    }

    fn constant(&self, n: i64) -> Frac<N::Elem> {
        (
            self.ring
                .constant(self.ring.base().elem_integer(&Integer::from(n))),
            self.ring.one(),
        )
    }

    /// An exact number `n / d` over one natural denominator.
    fn number(&self, z: &Complex) -> Frac<N::Elem> {
        let (num, den) = N::split_number(z);
        let base = self.ring.base();
        (
            self.ring.constant(num),
            self.ring.constant(base.elem_integer(&Integer::from(den))),
        )
    }

    /// Whether a product of `a` and `b` would pass the limits: more terms
    /// than allowed, by [`product_terms`], or an exponent past
    /// [`MAX_EXPONENT`].
    fn too_large(&self, a: &Sparse<N::Elem>, b: &Sparse<N::Elem>) -> bool {
        let degrees: Vec<u64> = a
            .degrees()
            .iter()
            .zip(b.degrees())
            .map(|(x, y)| x.saturating_add(y))
            .collect();
        degrees.iter().any(|&d| d > MAX_EXPONENT) || product_terms(a, b, &degrees) > self.max_terms
    }

    fn product(
        &self,
        a: &Sparse<N::Elem>,
        b: &Sparse<N::Elem>,
    ) -> Result<Sparse<N::Elem>, PolyError> {
        if self.too_large(a, b) {
            return Err(PolyError::TooLarge);
        }
        Ok(self.ring.mul_with(a, b, self.cancel)?)
    }

    fn mul(&self, a: &Frac<N::Elem>, b: &Frac<N::Elem>) -> Result<Frac<N::Elem>, PolyError> {
        let num = self.product(&a.0, &b.0)?;
        let den = self.product(&a.1, &b.1)?;
        Ok(self.normal((num, den)))
    }

    fn add(&self, a: &Frac<N::Elem>, b: &Frac<N::Elem>) -> Result<Frac<N::Elem>, PolyError> {
        if a.1 == b.1 {
            return Ok((self.ring.add(&a.0, &b.0), a.1.clone()));
        }
        // Over the lcm of the denominators. The gcd of two that are not
        // constants goes dense in their variables: past `max_terms` in
        // degree, too large (CAS-156), as `1/(x^(10^8) - 1) + 1/(y + 1)`
        // ran out of memory.
        let (ea, eb) = (a.1.degrees(), b.1.degrees());
        let constant = |d: &[u64]| d.iter().all(|&e| e == 0);
        if !constant(&ea)
            && !constant(&eb)
            && ea
                .iter()
                .chain(&eb)
                .any(|&e| usize::try_from(e).map_or(true, |d| d >= self.max_terms))
        {
            return Err(PolyError::TooLarge);
        }
        let (_, da, db) = self.ring.gcd_cofactors_with(&a.1, &b.1, self.cancel)?;
        let num = self
            .ring
            .add(&self.product(&a.0, &db)?, &self.product(&b.0, &da)?);
        let den = self.product(&a.1, &db)?;
        Ok(self.normal((num, den)))
    }

    fn pow(&self, f: &Frac<N::Elem>, n: i64) -> Result<Frac<N::Elem>, PolyError> {
        let (num, den) = if n < 0 {
            if f.0.is_zero() {
                return Err(PolyError::DivisionByZero);
            }
            (f.1.clone(), f.0.clone())
        } else {
            f.clone()
        };
        let e = n.unsigned_abs();
        let power = |p: &Sparse<N::Elem>| -> Result<Sparse<N::Elem>, PolyError> {
            if p.degrees()
                .iter()
                .any(|&d| d.saturating_mul(e) > MAX_EXPONENT)
            {
                return Err(PolyError::TooLarge);
            }
            // The coefficients of p^e have about e times the bits of p's
            // and of its count of terms: past the bound, too large, before
            // any product (CAS-156), as `(3 3^(1/3))^4782969` at a point
            // of `Trace::verify` took a minute where evaluation left it.
            if e > 1 {
                let largest = p.coeffs().iter().map(N::bits).max().unwrap_or(0);
                // log2 of the count of terms, rounded up: what the
                // multinomial coefficients add, each factor.
                let spread =
                    u64::from(usize::BITS - p.term_count().saturating_sub(1).leading_zeros());
                if largest.saturating_add(spread).saturating_mul(e) > self.max_bits {
                    return Err(PolyError::Bits);
                }
            }
            let mut result = self.ring.one();
            for i in (0..u64::BITS - e.leading_zeros()).rev() {
                result = self.product(&result, &result)?;
                if (e >> i) & 1 == 1 {
                    result = self.product(&result, p)?;
                }
            }
            Ok(result)
        };
        Ok(self.normal((power(&num)?, power(&den)?)))
    }

    /// The fraction with the denominator's leading coefficient in normal
    /// form, the numerator scaled by the same unit.
    fn normal(&self, (num, den): Frac<N::Elem>) -> Frac<N::Elem> {
        let base = self.ring.base();
        let Some(lead) = den.coeffs().first() else {
            return (num, den);
        };
        let unit = base.normal_form(lead).1;
        if base.is_one(&unit) {
            return (num, den);
        }
        (self.ring.scale(&num, &unit), self.ring.scale(&den, &unit))
    }
}
