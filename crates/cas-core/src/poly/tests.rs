//! Tests of the conversion between expressions and polynomials: the
//! generators and domains of CAS-139's table, the rational functions, the
//! errors, and the round trip.

use std::sync::atomic::AtomicBool;
use std::time::{Duration, Instant};

use cas_num::algebra::{Integers, Ring};
use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Never, Sparse, SparseRing};
use proptest::prelude::*;

use super::{Bounds, CoefficientDomain, Fraction, PolyError, Polynomials, product_terms};
use crate::build::{exp, sin};
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The value of `text` in full form, evaluated.
fn eval(text: &str) -> Expr {
    Context::new()
        .eval(&read(text))
        .unwrap_or_else(|e| panic!("{text:?}: {e:?}"))
}

fn polys(texts: &[&str]) -> Polynomials {
    let exprs: Vec<Expr> = texts.iter().map(|t| eval(t)).collect();
    Polynomials::new(&exprs, Bounds::NONE, &Never).unwrap_or_else(|e| panic!("{texts:?}: {e}"))
}

/// The generators, in full form.
fn generators(texts: &[&str]) -> Vec<String> {
    polys(texts)
        .generators()
        .iter()
        .map(|g| format!("{g:?}"))
        .collect()
}

/// CAS-139's table: the generators of each input, and its domain.
#[test]
#[allow(clippy::too_many_lines, reason = "a table")]
fn generators_and_domains() {
    let cases: &[(&[&str], &[&str], CoefficientDomain)] = &[
        // sin(x)^2 - 1: sin(x).
        (
            &["Plus[Power[Sin[x], 2], -1]"],
            &["Sin[x]"],
            CoefficientDomain::Integers,
        ),
        // e^x with e^(2x): e^(2x) is (e^x)^2.
        (
            &["Plus[Power[E, x], Power[E, Times[2, x]]]"],
            &["Power[E, x]"],
            CoefficientDomain::Integers,
        ),
        // e^x with e^(x/2): e^x is (e^(x/2))^2.
        (
            &["Plus[Power[E, x], Power[E, Times[1/2, x]]]"],
            &["Power[E, Times[1/2, x]]"],
            CoefficientDomain::Integers,
        ),
        // Exp[x] and Exp[2x] alike; Exp[x] and E^x stay apart.
        (
            &["Plus[Exp[x], Exp[Times[2, x]]]"],
            &["Exp[x]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Plus[Exp[x], Power[E, x]]"],
            &["Exp[x]", "Power[E, x]"],
            CoefficientDomain::Integers,
        ),
        // x with sqrt(x): two generators, sqrt(x) first in canonical order.
        (
            &["Plus[x, Power[x, 1/2]]"],
            &["x", "Power[x, 1/2]"],
            CoefficientDomain::Integers,
        ),
        // sqrt(2) is a generator in phase 1, before x.
        (
            &["Times[Power[2, 1/2], x]"],
            &["x", "Power[2, 1/2]"],
            CoefficientDomain::Integers,
        ),
        (&["f[x]"], &["f[x]"], CoefficientDomain::Integers),
        (
            &["Power[x, n]"],
            &["Power[x, n]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Power[x, 1000000000000000000000000000000]"],
            &["Power[x, 1000000000000000000000000000000]"],
            CoefficientDomain::Integers,
        ),
        // Roots split as giac does: x^(3/2) = x x^(1/2), x^(2/3) =
        // (x^(1/3))^2; roots of different orders stay apart.
        (
            &["Plus[Power[x, 3/2], Power[x, 1/2], x]"],
            &["x", "Power[x, 1/2]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Plus[Power[x, 2/3], Power[x, 1/3]]"],
            &["Power[x, 1/3]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Plus[Power[x, 1/2], Power[x, 1/3]]"],
            &["Power[x, 1/2]", "Power[x, 1/3]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Power[x, -1/2]"],
            &["Power[x, 1/2]"],
            CoefficientDomain::Integers,
        ),
        // Symbolic exponents with rational multiples, as e^x.
        (
            &["Plus[Power[x, Times[3, n]], Power[x, Times[2, n]]]"],
            &["Power[x, n]"],
            CoefficientDomain::Integers,
        ),
        // Bases apart: 2^x and 4^x.
        (
            &["Plus[Power[2, x], Power[4, x]]"],
            &["Power[4, x]", "Power[2, x]"],
            CoefficientDomain::Integers,
        ),
        // Constants and several variables, the last in canonical order
        // first.
        (
            &["Plus[x, y, Times[Pi, z]]"],
            &["z", "y", "x", "Pi"],
            CoefficientDomain::Integers,
        ),
        // The domains.
        (&["Times[1/2, x]"], &["x"], CoefficientDomain::Rationals),
        (
            &["Times[Complex[0, 1], x]"],
            &["x"],
            CoefficientDomain::GaussianIntegers,
        ),
        (
            &["Plus[Times[1/2, x], Complex[0, 1]]"],
            &["x"],
            CoefficientDomain::GaussianRationals,
        ),
        (
            &["Times[Complex[1/2, 1], x]"],
            &["x"],
            CoefficientDomain::GaussianRationals,
        ),
        // A number inside a generator is not a coefficient.
        (&["Sin[1/2]"], &["Sin[1/2]"], CoefficientDomain::Integers),
        // Complex constants in exponents: E^I, and E^(1 + 2 I) = E (E^I)^2.
        (
            &["Plus[Power[E, Complex[0, 1]], x]"],
            &["x", "Power[E, Complex[0, 1]]"],
            CoefficientDomain::Integers,
        ),
        (
            &["Power[E, Complex[1, 2]]"],
            &["Power[E, Complex[0, 2]]", "E"],
            CoefficientDomain::Integers,
        ),
        (
            &["Plus[Power[E, Complex[1, 2]], Power[E, Complex[0, 1]]]"],
            &["Power[E, Complex[0, 1]]", "E"],
            CoefficientDomain::Integers,
        ),
    ];
    for (inputs, expected, domain) in cases {
        let p = polys(inputs);
        let found: Vec<String> = p.generators().iter().map(|g| format!("{g:?}")).collect();
        assert_eq!(found, *expected, "{inputs:?}");
        assert_eq!(p.domain(), *domain, "{inputs:?}");
    }
}

/// The integer polynomial of `terms` over the generators of `p`.
fn integer_poly(p: &Polynomials, terms: &[(&[u64], i64)]) -> Sparse<Integer> {
    let ring = SparseRing::new(Integers, p.generators().len());
    ring.from_terms(terms.iter().map(|(e, c)| (e.to_vec(), Integer::from(*c))))
}

fn integer_fraction(p: &Polynomials, k: usize) -> (&Sparse<Integer>, &Sparse<Integer>) {
    match &p.fractions()[k] {
        Fraction::Integers { num, den } => (num, den),
        Fraction::GaussianIntegers { .. } => panic!("a complex fraction"),
    }
}

#[test]
fn fractions() {
    // (x + 1)^2 = x^2 + 2x + 1, over 1.
    let p = polys(&["Power[Plus[x, 1], 2]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[2], 1), (&[1], 2), (&[0], 1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[0], 1)]));
    // x + 1/x = (x^2 + 1) / x.
    let p = polys(&["Plus[x, Power[x, -1]]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[2], 1), (&[0], 1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[1], 1)]));
    // 1/x + 1/x^2 = (x + 1) / x^2, over the lcm of the denominators.
    let p = polys(&["Plus[Power[x, -1], Power[x, -2]]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[1], 1), (&[0], 1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[2], 1)]));
    // x/2 - 1/3 = (3x - 2) / 6.
    let p = polys(&["Plus[Times[1/2, x], -1/3]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[1], 3), (&[0], -2)]));
    assert_eq!(den, &integer_poly(&p, &[(&[0], 6)]));
    // 1/(1 - x): the denominator's leading coefficient positive.
    let p = polys(&["Power[Plus[1, Times[-1, x]], -1]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[0], -1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[1], 1), (&[0], -1)]));
    // x^(3/2) = x · x^(1/2); x^(-3/2) = 1 / (x · x^(1/2)).
    let p = polys(&["Power[x, 3/2]", "Power[x, -3/2]"]);
    let (num, den) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[1, 1], 1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[0, 0], 1)]));
    let (num, den) = integer_fraction(&p, 1);
    assert_eq!(num, &integer_poly(&p, &[(&[0, 0], 1)]));
    assert_eq!(den, &integer_poly(&p, &[(&[1, 1], 1)]));
    // e^(2x) with e^(x/2) is the fourth power of e^(x/2).
    let p = polys(&["Plus[Power[E, Times[2, x]], Power[E, Times[1/2, x]]]"]);
    let (num, _) = integer_fraction(&p, 0);
    assert_eq!(num, &integer_poly(&p, &[(&[4], 1), (&[1], 1)]));
    // Several expressions share their generators: y is x₀.
    let p = polys(&[
        "Plus[Power[x, 2], Times[-1, Power[y, 2]]]",
        "Plus[x, Times[-1, y]]",
    ]);
    assert_eq!(generators(&["Plus[x, y]"]), ["y", "x"]);
    let (num, _) = integer_fraction(&p, 1);
    assert_eq!(num, &integer_poly(&p, &[(&[1, 0], -1), (&[0, 1], 1)]));
    // Complex: (1 + i) x over ℤ[i].
    let p = polys(&["Times[Complex[1, 1], x]"]);
    match &p.fractions()[0] {
        Fraction::GaussianIntegers { num, .. } => {
            assert_eq!(
                num.coeffs(),
                &[Complex::new(Rational::from(1), Rational::from(1))]
            );
        }
        Fraction::Integers { .. } => panic!("a complex number"),
    }
}

/// A long sum converts in near-linear work (CAS-157): 30276 monomials,
/// which adding each to the sum so far took minutes for here. The bound is
/// loose, for a debug build on a loaded machine.
#[test]
fn long_sums_convert_fast() {
    let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
    let terms: Vec<Expr> = (0..174u64)
        .flat_map(|i| {
            let (x, y) = (x.clone(), y.clone());
            (0..174u64).map(move |j| {
                Expr::call(
                    Symbol::TIMES,
                    [
                        Expr::from(Integer::from(i + j + 1)),
                        x.pow(Expr::from(Integer::from(i))),
                        y.pow(Expr::from(Integer::from(j))),
                    ],
                )
            })
        })
        .collect();
    let sum = Expr::call(Symbol::PLUS, terms);
    let start = Instant::now();
    let p = Polynomials::new(&[sum], Bounds::NONE, &Never).unwrap();
    let elapsed = start.elapsed();
    let (num, _) = integer_fraction(&p, 0);
    assert_eq!(num.term_count(), 174 * 174);
    assert!(elapsed < Duration::from_secs(20), "{elapsed:?}");
}

#[test]
fn errors() {
    let exprs = |texts: &[&str]| -> Vec<Expr> { texts.iter().map(|t| eval(t)).collect() };
    // A ball coefficient.
    let e = exprs(&["Times[x, N[3/2]]"]);
    assert!(matches!(
        Polynomials::new(&e, Bounds::NONE, &Never),
        Err(PolyError::Ball(_))
    ));
    // Infinities, even inside a generator.
    for text in [
        "Plus[x, Infinity]",
        "Power[E, ComplexInfinity]",
        "Times[x, Indeterminate]",
    ] {
        let e = exprs(&[text]);
        assert!(
            matches!(
                Polynomials::new(&e, Bounds::NONE, &Never),
                Err(PolyError::Infinite(_))
            ),
            "{text}"
        );
    }
    // The lcm of two denominators that are not constants takes a gcd that
    // goes dense in their variables (CAS-156): past the terms limit in
    // degree, too large, before any work.
    for text in [
        "Plus[Power[Plus[Power[x, 10000000], -1], -1], Power[Plus[y, 1], -1]]",
        "Plus[Power[Plus[Power[x, 10000000], -1], -1], Power[Plus[Power[x, 1000000], -1], -1]]",
    ] {
        let e = exprs(&[text]);
        assert_eq!(
            Polynomials::new(&e, Bounds::terms(1_000_000), &Never).err(),
            Some(PolyError::TooLarge),
            "{text}"
        );
    }
    // A power whose coefficients would pass the bound on bits (CAS-156),
    // before any product; a number to a power that evaluation left is a
    // generator, not a number to compute.
    let big = [read("Power[Times[3, Power[3, 1/3]], 100000]")];
    let bits = Bounds {
        terms: usize::MAX,
        bits: 4096,
    };
    assert_eq!(
        Polynomials::new(&big, bits, &Never).err(),
        Some(PolyError::Bits)
    );
    assert!(Polynomials::new(&big, Bounds::NONE, &Never).is_ok());
    let left = [read("Plus[x, Power[3, 1000000], Power[2, 3000001/3]]")];
    let p = Polynomials::new(&left, bits, &Never).unwrap();
    assert_eq!(
        p.generators(),
        [
            read("x"),
            read("Power[3, 1000000]"),
            read("Power[2, 3000001/3]")
        ]
    );
    // A constant denominator needs no gcd of polynomials.
    let e = exprs(&["Plus[Power[x, 10000000], Times[1/2, y]]"]);
    assert!(Polynomials::new(&e, Bounds::terms(1_000_000), &Never).is_ok());
    // Past the terms limit: (1 + x + y)^10 has 66 terms.
    let e = exprs(&["Power[Plus[1, x, y], 10]"]);
    assert_eq!(
        Polynomials::new(&e, Bounds::terms(30), &Never).err(),
        Some(PolyError::TooLarge)
    );
    assert!(Polynomials::new(&e, Bounds::terms(1000), &Never).is_ok());
    // Exponents past 2^62.
    let e = exprs(&["Power[Plus[1, x], 4611686018427387905]"]);
    assert_eq!(
        Polynomials::new(&e, Bounds::NONE, &Never).err(),
        Some(PolyError::TooLarge)
    );
    // A zero denominator that evaluation leaves: it never expands.
    let e = exprs(&[
        "Power[Plus[Power[Plus[1, x], 2], Times[-1, Plus[1, Times[2, x], Power[x, 2]]]], -1]",
    ]);
    assert_eq!(
        Polynomials::new(&e, Bounds::NONE, &Never).err(),
        Some(PolyError::DivisionByZero)
    );
    // The interrupt.
    let stop = AtomicBool::new(true);
    let e = exprs(&["Plus[x, 1]"]);
    assert_eq!(
        Polynomials::new(&e, Bounds::NONE, &stop).err(),
        Some(PolyError::Interrupted)
    );
}

/// The way back: built, then evaluated, to the canonical forms.
#[test]
fn back_to_expressions() {
    let ctx = Context::new();
    let cases = [
        ("Power[Plus[x, 1], 2]", "Plus[1, Times[2, x], Power[x, 2]]"),
        (
            "Plus[x, Power[x, -1]]",
            "Times[Power[x, -1], Plus[1, Power[x, 2]]]",
        ),
        ("Power[x, 3/2]", "Power[x, 3/2]"),
        ("Power[x, -3/2]", "Power[x, -3/2]"),
        (
            "Plus[Power[E, Times[2, x]], Power[E, Times[1/2, x]]]",
            "Plus[Power[E, Times[1/2, x]], Power[E, Times[2, x]]]",
        ),
        ("Times[1/2, Plus[1, x]]", "Plus[1/2, Times[1/2, x]]"),
        (
            "Times[Complex[0, 1], Plus[x, 1]]",
            "Plus[Complex[0, 1], Times[Complex[0, 1], x]]",
        ),
        (
            "Power[Plus[Sin[x], 1], 2]",
            "Plus[1, Times[2, Sin[x]], Power[Sin[x], 2]]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let p = polys(&[input]);
        let back = ctx.eval(&p.fraction_expr(&p.fractions()[0])).unwrap();
        if format!("{back:?}") != expected {
            wrong.push(format!("{input}: {back:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The cases the round trip found, where evaluation adds exponents:
/// each converts, with its way back, to one rational function.
#[test]
fn exponent_algebra() {
    let ctx = Context::new();
    let cases = [
        // E E^x is E^(1 + x): sums in exponents split.
        "Plus[Power[E, -1], Power[E, x]]",
        // A multiple of a sum: E^(-(x + Exp[x])) against E^(x + Exp[x]).
        "Plus[Power[E, Times[-1, Plus[x, Exp[x]]]], x]",
        "Times[Power[E, Times[-1, Plus[-1, x]]], Plus[1, Times[-1, Power[E, Plus[-1, x]]]]]",
        // Complex constants and multiples: E E^I is E^(1 + I), and
        // (1 + I) z + (1 - I) z is 2 z.
        "Times[E, Plus[Power[E, Complex[0, 1]], x]]",
        "Power[E, Times[2, Plus[Complex[0, -1], x]]]",
        "Power[E, Plus[z, Times[-1, Plus[x, Times[Complex[-1, 1], z]]]]]",
        "Times[Power[E, Times[Complex[1, 1], z]], Plus[x, Power[E, Times[Complex[1, -1], z]]]]",
        "Power[E, Times[-2, Plus[Times[Complex[0, -1], Power[x, -1]], x]]]",
    ];
    for text in cases {
        let e = eval(text);
        let p = Polynomials::new(std::slice::from_ref(&e), Bounds::NONE, &Never).unwrap();
        let back = ctx.eval(&p.fraction_expr(&p.fractions()[0])).unwrap();
        let both = Polynomials::new(&[e.clone(), back.clone()], Bounds::NONE, &Never).unwrap();
        assert!(
            same(&both.fractions()[0], &both.fractions()[1]),
            "{text}: {e} and {back}"
        );
    }
    // E^((1 + I) z) is E^z E^(I z), over the basis z, I z.
    assert_eq!(
        generators(&["Power[E, Times[Complex[1, 1], z]]"]),
        ["Power[E, Times[Complex[0, 1], z]]", "Power[E, z]"]
    );
}

/// A random expression in `x`, `y` and `z`: sums of products of numbers,
/// variables, integer powers, and sometimes roots, exponentials and sines.
#[derive(Clone, Debug)]
enum Gen {
    Var(usize),
    Int(i64),
    Rat(i64, i64),
    Gauss(i64, i64),
    Sum(Vec<Gen>),
    Product(Vec<Gen>),
    Pow(Box<Gen>, i64),
    Root(Box<Gen>, i64, i64),
    Exp(Box<Gen>, i64),
    Sin(Box<Gen>),
}

/// Which kernels the random expressions may have.
#[derive(Clone, Copy)]
struct Kinds {
    roots: bool,
    /// Roots of odd order too, which O7 takes as real for negatives.
    odd_roots: bool,
    transcendental: bool,
}

fn gen_strategy(kinds: Kinds) -> impl Strategy<Value = Gen> {
    let leaf = prop_oneof![
        6 => (0usize..3).prop_map(Gen::Var),
        2 => (-5i64..6).prop_map(Gen::Int),
        1 => (-5i64..6, 1i64..5).prop_map(|(n, d)| Gen::Rat(n, d)),
        1 => (-3i64..4, -3i64..4).prop_map(|(a, b)| Gen::Gauss(a, b)),
    ];
    leaf.prop_recursive(3, 40, 4, move |inner| {
        let mut factors: Vec<(u32, BoxedStrategy<Gen>)> = vec![
            (4, inner.clone().boxed()),
            (
                2,
                (inner.clone(), -2i64..4)
                    .prop_map(|(a, n)| Gen::Pow(Box::new(a), n))
                    .boxed(),
            ),
        ];
        if kinds.roots {
            let orders = if kinds.odd_roots {
                vec![2i64, 3]
            } else {
                vec![2, 4]
            };
            factors.push((
                1,
                (inner.clone(), -3i64..4, prop::sample::select(orders))
                    .prop_map(|(a, p, q)| Gen::Root(Box::new(a), p, q))
                    .boxed(),
            ));
        }
        if kinds.transcendental {
            factors.push((
                1,
                (inner.clone(), -2i64..3)
                    .prop_map(|(a, c)| Gen::Exp(Box::new(a), c))
                    .boxed(),
            ));
            factors.push((1, inner.prop_map(|a| Gen::Sin(Box::new(a))).boxed()));
        }
        let factor = proptest::strategy::Union::new_weighted(factors);
        let product = prop::collection::vec(factor, 1..4).prop_map(Gen::Product);
        prop::collection::vec(product, 1..4).prop_map(Gen::Sum)
    })
}

fn build(g: &Gen) -> Expr {
    let vars = [Expr::symbol("x"), Expr::symbol("y"), Expr::symbol("z")];
    match g {
        Gen::Var(v) => vars[*v].clone(),
        Gen::Int(n) => Expr::from(*n),
        Gen::Rat(n, d) => Expr::from(Rational::new(Integer::from(*n), Integer::from(*d))),
        Gen::Gauss(a, b) => Expr::from(Complex::new(Rational::from(*a), Rational::from(*b))),
        Gen::Sum(terms) => terms.iter().map(build).sum(),
        Gen::Product(factors) => factors.iter().map(build).product(),
        Gen::Pow(a, n) => build(a).pow(*n),
        Gen::Root(a, p, q) => build(a).pow(Rational::new(Integer::from(*p), Integer::from(*q))),
        Gen::Exp(a, c) => {
            let e = Expr::from(Symbol::E);
            if *c == 0 {
                exp(build(a))
            } else {
                e.pow(Expr::from(*c) * build(a))
            }
        }
        Gen::Sin(a) => sin(build(a)),
    }
}

/// `num_a den_b = num_b den_a`: the same rational function.
fn same(a: &Fraction, b: &Fraction) -> bool {
    fn cross<E: Clone + PartialEq + std::fmt::Debug, R: Ring<Elem = E>>(
        ring: &SparseRing<R>,
        (na, da): (&Sparse<E>, &Sparse<E>),
        (nb, db): (&Sparse<E>, &Sparse<E>),
    ) -> bool {
        ring.mul(na, db) == ring.mul(nb, da)
    }
    match (a, b) {
        (Fraction::Integers { num: na, den: da }, Fraction::Integers { num: nb, den: db }) => {
            cross(&SparseRing::new(Integers, na.nvars()), (na, da), (nb, db))
        }
        (
            Fraction::GaussianIntegers { num: na, den: da },
            Fraction::GaussianIntegers { num: nb, den: db },
        ) => cross(
            &SparseRing::new(cas_num::algebra::GaussianIntegers, na.nvars()),
            (na, da),
            (nb, db),
        ),
        _ => false,
    }
}

/// The rational functions of `e` and of its way back, evaluated, over one
/// set of generators; `None` where the conversion refuses `e`.
fn there_and_back(ctx: &Context, e: &Expr) -> Result<Option<(Expr, Polynomials)>, TestCaseError> {
    let p = match Polynomials::new(std::slice::from_ref(e), Bounds::terms(100_000), &Never) {
        Ok(p) => p,
        // 0^-1 and the like, which evaluation leaves unexpanded, and the
        // infinities of 1/0.
        Err(PolyError::DivisionByZero | PolyError::TooLarge | PolyError::Infinite(_)) => {
            return Ok(None);
        }
        Err(error) => return Err(TestCaseError::fail(format!("{e}: {error}"))),
    };
    let back = ctx.eval(&p.fraction_expr(&p.fractions()[0])).unwrap();
    let both = Polynomials::new(&[e.clone(), back.clone()], Bounds::terms(100_000), &Never)
        .map_err(|error| TestCaseError::fail(format!("{back}: {error}")))?;
    Ok(Some((back, both)))
}

/// A numeric value as a complex ball: from a ball, a complex ball or an
/// exact number; `None` for anything else (an infinity, an unevaluated
/// part).
fn numeric(e: &Expr) -> Option<ComplexBall> {
    match e.view() {
        crate::View::Ball(b) => Some(ComplexBall::new(b.clone(), Ball::ZERO)),
        crate::View::ComplexBall(z) => Some(z.clone()),
        _ => {
            let z = e.to_complex()?;
            Some(ComplexBall::new(
                Ball::from_rational(z.re(), 128),
                Ball::from_rational(z.im(), 128),
            ))
        }
    }
}

/// Whether a generator is a root `b^(1/q)`.
fn has_roots(p: &Polynomials) -> bool {
    p.generators().iter().any(|g| {
        g.head().and_then(Expr::as_symbol) == Some(Symbol::POWER)
            && g.args()[1].to_rational().is_some_and(|r| !r.is_integer())
    })
}

/// Points where roots of order 2 and 3 of the variables are exact: sixth
/// powers of rationals.
fn sixth_powers(point: &[(i64, i64)]) -> Vec<(Expr, Expr)> {
    ["x", "y", "z"]
        .iter()
        .zip(point)
        .map(|(v, (n, d))| {
            let r = Rational::new(Integer::from(*n), Integer::from(*d));
            let sixth = (0..6).fold(Rational::from(1), |a, _| &a * &r);
            (Expr::symbol(v), Expr::from(sixth))
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Without roots, converted and back, an expression evaluates to the
    /// same rational function: both convert, over one set of generators,
    /// to equal fractions; and so does a second round trip.
    #[test]
    fn round_trips(g in gen_strategy(Kinds { roots: false, odd_roots: false, transcendental: true })) {
        let ctx = Context::new();
        let Ok(e) = ctx.eval(&build(&g)) else { return Ok(()); };
        let Some((back, both)) = there_and_back(&ctx, &e)? else { return Ok(()); };
        // Exponentials can make roots too: E^(5/2) is E^2 (E^(1/2)).
        if has_roots(&both) {
            return Ok(());
        }
        prop_assert!(same(&both.fractions()[0], &both.fractions()[1]), "{} and {}", e, back);
        let Some((again, both)) = there_and_back(&ctx, &back)? else { return Ok(()); };
        prop_assert!(same(&both.fractions()[0], &both.fractions()[1]), "{} and {}", back, again);
    }

    /// With roots, evaluation brings `(x^(1/2))^2` back to `x`, a generator
    /// of its own, and `((x + 1)^(1/2))^2` to `x + 1`, unexpanded. So the
    /// way back keeps the value, not the fraction: equal at points where
    /// the roots are exact.
    #[test]
    fn roots_keep_the_value(
        g in gen_strategy(Kinds { roots: true, odd_roots: true, transcendental: false }),
        point in prop::collection::vec((1i64..4, 1i64..3), 3),
    ) {
        let ctx = Context::new();
        let Ok(e) = ctx.eval(&build(&g)) else { return Ok(()); };
        let Some((back, _)) = there_and_back(&ctx, &e)? else { return Ok(()); };
        let values = sixth_powers(&point);
        let at = |f: &Expr| ctx.eval(&f.substitute(&values)).unwrap();
        let (u, v) = (at(&e), at(&back));
        if u.is_exact_number() && v.is_exact_number() {
            prop_assert_eq!(u, v, "{} and {}", e, back);
        }
    }

    /// Every kind of kernel: the values at points agree, to 20 digits, the
    /// balls of `N` overlapping. Odd roots included: evaluation keeps what
    /// O7 reads as a real odd root apart from the principal branch
    /// (CAS-154).
    #[test]
    fn values_agree_numerically(
        g in gen_strategy(Kinds { roots: true, odd_roots: true, transcendental: true }),
        point in prop::collection::vec((-4i64..5, 1i64..4), 3),
    ) {
        let ctx = Context::new();
        let Ok(e) = ctx.eval(&build(&g)) else { return Ok(()); };
        let Some((back, _)) = there_and_back(&ctx, &e)? else { return Ok(()); };
        let values: Vec<(Expr, Expr)> = ["x", "y", "z"]
            .iter()
            .zip(&point)
            .map(|(v, (n, d))| (Expr::symbol(v), Expr::from(Rational::new(Integer::from(*n), Integer::from(*d)))))
            .collect();
        let n = |f: &Expr| {
            let at = f.substitute(&values);
            ctx.eval(&Expr::call(Symbol::N, [at, Expr::from(20)])).ok()
        };
        let (Some(u), Some(v)) = (n(&e), n(&back)) else { return Ok(()); };
        if let (Some(u), Some(v)) = (numeric(&u), numeric(&v)) {
            prop_assert!(u.overlaps(&v), "{} and {} at {:?}: {:?} and {:?}", e, back, point, u, v);
        }
    }

    /// For polynomial expressions, the values at points agree.
    #[test]
    fn values_agree(g in gen_strategy(Kinds { roots: false, odd_roots: false, transcendental: false }), point in prop::collection::vec((-4i64..5, 1i64..4), 3)) {
        let ctx = Context::new();
        let Ok(e) = ctx.eval(&build(&g)) else { return Ok(()); };
        let Ok(p) = Polynomials::new(std::slice::from_ref(&e), Bounds::terms(100_000), &Never) else {
            return Ok(());
        };
        let back = ctx.eval(&p.fraction_expr(&p.fractions()[0])).unwrap();
        let values: Vec<(Expr, Expr)> = ["x", "y", "z"]
            .iter()
            .zip(&point)
            .map(|(v, (n, d))| (Expr::symbol(v), Expr::from(Rational::new(Integer::from(*n), Integer::from(*d)))))
            .collect();
        let at = |f: &Expr| ctx.eval(&f.substitute(&values)).unwrap();
        let (u, v) = (at(&e), at(&back));
        // A pole of one is a pole of the other; compare finite values.
        if u.is_exact_number() && v.is_exact_number() {
            prop_assert_eq!(u, v);
        }
    }

    /// The bound of a product's terms is at least its terms, and is them
    /// for products of dense polynomials, every monomial up to a total
    /// degree.
    #[test]
    fn product_terms_bound_the_product(
        a in prop::collection::vec((prop::collection::vec(0u64..4, 3), -3i64..4), 0..12),
        b in prop::collection::vec((prop::collection::vec(0u64..4, 3), -3i64..4), 0..12),
        (n, da, db) in (1usize..4, 0u64..5, 0u64..5),
    ) {
        let ring = SparseRing::new(Integers, 3);
        let poly = |t: &[(Vec<u64>, i64)]| ring.from_terms(t.iter().map(|(e, c)| (e.clone(), Integer::from(*c))));
        let (a, b) = (poly(&a), poly(&b));
        let bound = |a: &Sparse<Integer>, b: &Sparse<Integer>| {
            let degrees: Vec<u64> = a.degrees().iter().zip(b.degrees()).map(|(x, y)| x + y).collect();
            product_terms(a, b, &degrees)
        };
        prop_assert!(bound(&a, &b) >= ring.mul(&a, &b).term_count());
        // Dense: all monomials of total degree at most d in n variables.
        let dense = |d: u64| {
            let mut terms = Vec::new();
            for e0 in 0..=d {
                for e1 in 0..=d - e0 {
                    for e2 in 0..=d - e0 - e1 {
                        let e = [e0, e1, e2];
                        if e[n..].iter().all(|&k| k == 0) {
                            terms.push((e.to_vec(), Integer::from(1)));
                        }
                    }
                }
            }
            ring.from_terms(terms)
        };
        let (p, q) = (dense(da), dense(db));
        prop_assert_eq!(bound(&p, &q), ring.mul(&p, &q).term_count());
    }
}
