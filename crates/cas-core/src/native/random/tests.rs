//! Tests of random numbers (CAS-370 to CAS-373): seeds, draws taken anew
//! wherever they are met, `Table`'s held body, verification, the ranges of
//! `Random` and `RandomInteger`, the variates' supports, vectors, matrices
//! and polynomials, permutations and samples.

use cas_num::random::Generator;
use cas_num::{Integer, Rational};

use crate::{Context, EvalOptions, Expr, Symbol, TraceLevel, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// A number as a rational: itself, or an exact ball's midpoint.
fn value(e: &Expr) -> Rational {
    match e.view() {
        View::Ball(b) => b.mid().to_rational(),
        _ => e
            .to_rational()
            .unwrap_or_else(|| panic!("{e:?}: no number")),
    }
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// `SeedRandom[n]` gives `n` and seeds with it; `SeedRandom[]` gives the
/// seed it took; a non-integer leaves the call.
#[test]
fn seeds() {
    let ctx = Context::new();
    for (input, seed) in [
        ("SeedRandom[5]", Integer::from(5)),
        ("SeedRandom[-2]", Integer::from(-2)),
        (
            "SeedRandom[Power[2, 70]]",
            Integer::from(1i64 << 35) * Integer::from(1i64 << 35),
        ),
    ] {
        assert_eq!(eval(&ctx, input), Expr::from(seed.clone()), "{input}");
        assert_eq!(ctx.random().generator, Generator::seeded(&seed), "{input}");
    }
    let clock = eval(&ctx, "SeedRandom[]");
    let seed = clock.to_integer().expect("an integer seed");
    assert_eq!(ctx.random().generator, Generator::seeded(&seed));
    for input in ["SeedRandom[x]", "SeedRandom[1/2]", "SeedRandom[1, 2]"] {
        assert_eq!(eval(&ctx, input), read(input), "{input}");
    }
}

/// Each `Random[]` met is a draw of its own: in a list, in `Table`'s body,
/// after a seed met twice; a new context starts from one seed.
#[test]
fn draws_are_new() {
    let ctx = Context::new();
    let two = eval(&ctx, "List[Random[], Random[]]");
    assert_ne!(two.args()[0], two.args()[1], "{two:?}");
    let three = eval(&ctx, "Table[Random[10000], k, 1, 3]");
    let values: Vec<Rational> = three.args().iter().map(value).collect();
    assert!(
        values[0] != values[1] || values[1] != values[2],
        "{three:?}"
    );
    let seeded = eval(
        &ctx,
        "List[SeedRandom[5], Random[], SeedRandom[5], Random[]]",
    );
    assert_eq!(seeded.args()[1], seeded.args()[3], "{seeded:?}");
    assert_ne!(seeded.args()[1], two.args()[0]);
    let first = |ctx: &Context| eval(ctx, "Random[]");
    assert_eq!(first(&Context::new()), first(&Context::new()));
}

/// `N[Random[]]` draws once: the draw after it is the second of the seed.
#[test]
fn numeric_draws_once() {
    let ctx = Context::new();
    let with_n = eval(&ctx, "List[SeedRandom[5], N[Random[]], Random[]]");
    let without = eval(&ctx, "List[SeedRandom[5], Random[], Random[]]");
    assert_eq!(with_n.args()[2], without.args()[2]);
}

/// A step that draws is unchecked, and verifying leaves the context's
/// sequence as it was.
#[test]
fn verification() {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[5]");
    let options = EvalOptions {
        trace: TraceLevel::Verbose,
    };
    let (result, report) = ctx.eval_with(&read("Plus[Random[], 1]"), &options);
    assert!(result.is_ok());
    let verification = report.trace.verify(&ctx);
    assert!(verification.wrong.is_none(), "{verification:?}");
    assert!(!verification.unchecked.is_empty(), "{verification:?}");
    let next = eval(&ctx, "Random[]");
    let again = eval(&ctx, "List[SeedRandom[5], Random[], Random[]]");
    assert_eq!(next, again.args()[2]);
}

/// `Table` holds its body (CAS-370): evaluated at each value, an index
/// with a value bound anew, a value that holds the index given it again.
#[test]
fn table_holds_its_body() {
    let mut bound = Context::new();
    bound.set(Symbol::new("k"), Expr::from(7)).unwrap();
    assert_eq!(
        eval(&bound, "Table[Power[k, 2], k, 1, 3]"),
        read("List[1, 4, 9]")
    );
    let ctx = Context::new();
    let mut with_f = ctx.clone();
    with_f.set(Symbol::new("f"), read("Power[k, 2]")).unwrap();
    assert_eq!(eval(&with_f, "Table[f, k, 1, 3]"), read("List[1, 4, 9]"));
    assert_eq!(
        eval(&ctx, "Table[Integrate[Power[x, k], x], k, -1, 1]"),
        eval(&ctx, "List[Log[x], x, Times[1/2, Power[x, 2]]]")
    );
    assert_eq!(
        eval(&ctx, "Table[k, Equation[k, Range[Plus[1, 1], 4]]]"),
        read("List[2, 3, 4]")
    );
}

/// The draws of `forms` from the seed 1, `times` of each.
fn draws(forms: &[&str], times: usize) -> Vec<Vec<Expr>> {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[1]");
    forms
        .iter()
        .map(|form| (0..times).map(|_| eval(&ctx, form)).collect())
        .collect()
}

/// Whether `x` is an integer in `[low, high]`.
fn integer_in(x: &Rational, low: i64, high: i64) -> bool {
    x.is_integer() && q(low, 1) <= *x && *x <= q(high, 1)
}

/// Reals: exact balls in `[0, 1)`, and between the ends given.
#[test]
fn real_ranges() {
    let found = draws(
        &["Random[]", "Random[2, 5]", "Random[5, 2]", "Random[0, Pi]"],
        200,
    );
    for unit in &found[0] {
        let View::Ball(ball) = unit.view() else {
            panic!("{unit:?}: no ball");
        };
        assert!(ball.is_exact(), "{unit:?}");
        let fraction = value(unit);
        assert!(
            Rational::ZERO <= fraction && fraction < Rational::ONE,
            "{unit:?}"
        );
    }
    for drawn in found[1].iter().chain(&found[2]) {
        let real = value(drawn);
        assert!(q(2, 1) <= real && real <= q(5, 1), "{drawn:?}");
    }
    for drawn in &found[3] {
        let real = value(drawn);
        assert!(Rational::ZERO <= real && real < q(315, 100), "{drawn:?}");
    }
}

/// Integers: `[0, 6)`, each met; `(−6, 0]`; `randint`'s ranges.
#[test]
fn integer_ranges() {
    let found = draws(
        &[
            "Random[6]",
            "Random[-6]",
            "RandomInteger[6]",
            "RandomInteger[-2, 2]",
            "RandomInteger[2, -2]",
            "RandomInteger[0]",
        ],
        200,
    );
    let ranges = [(0, 5), (-5, 0), (1, 6), (-2, 2), (-2, 2), (0, 1)];
    for (drawn, (low, high)) in found.iter().zip(ranges) {
        let values: Vec<Rational> = drawn.iter().map(value).collect();
        for v in &values {
            assert!(integer_in(v, low, high), "{v} in [{low}, {high}]");
        }
        for target in low..=high {
            assert!(values.contains(&q(target, 1)), "{target} never drawn");
        }
    }
    let ctx = Context::new();
    let big = value(&eval(&ctx, "RandomInteger[Power[10, 30]]"));
    let top = Rational::from(Integer::from(10).pow(30));
    assert!(
        big.is_integer() && Rational::ONE <= big && big <= top,
        "{big}"
    );
}

/// Distinct integers, elements of lists, and samples at distinct places;
/// all of `1 … 49` is a permutation.
#[test]
fn choices() {
    let found = draws(
        &[
            "Random[3, 1, 49]",
            "Random[List[a, b, c]]",
            "Random[2, List[a, b, c, d]]",
        ],
        200,
    );
    for lotto in &found[0] {
        let picked: Vec<Rational> = lotto.args().iter().map(value).collect();
        assert_eq!(picked.len(), 3, "{lotto:?}");
        for (j, p) in picked.iter().enumerate() {
            assert!(integer_in(p, 1, 49), "{lotto:?}");
            assert!(!picked[..j].contains(p), "{lotto:?}");
        }
    }
    let letters = read("List[a, b, c]");
    for element in &found[1] {
        assert!(letters.args().contains(element), "{element:?}");
    }
    for pair in &found[2] {
        assert_eq!(pair.args().len(), 2, "{pair:?}");
        assert_ne!(pair.args()[0], pair.args()[1], "{pair:?}");
    }
    let ctx = Context::new();
    let all = eval(&ctx, "Random[49, 1, 49]");
    let mut sorted: Vec<Rational> = all.args().iter().map(value).collect();
    sorted.sort();
    let expected: Vec<Rational> = (1..=49).map(|i| q(i, 1)).collect();
    assert_eq!(sorted, expected);
}

/// Variates in their laws' supports (CAS-372): reals anywhere, reals
/// from 0, integers in their supports.
#[test]
fn variates() {
    let found = draws(
        &[
            "RandomNormal[]",
            "RandomNormal[1, 2]",
            "RandomStudentT[3]",
            "RandomExponential[2]",
            "RandomChiSquare[3]",
            "RandomFisher[2, 3]",
            "RandomBinomial[10, 1/2]",
            "RandomPoisson[3]",
            "RandomGeometric[1/3]",
        ],
        40,
    );
    for drawn in found[..6].iter().flatten() {
        assert!(matches!(drawn.view(), View::Ball(_)), "{drawn:?}");
    }
    for drawn in found[3..6].iter().flatten() {
        assert!(value(drawn) >= Rational::ZERO, "{drawn:?}");
    }
    for drawn in &found[6] {
        assert!(integer_in(&value(drawn), 0, 10), "{drawn:?}");
    }
    for drawn in &found[7] {
        let count = value(drawn);
        assert!(count.is_integer() && count >= Rational::ZERO, "{drawn:?}");
    }
    for drawn in &found[8] {
        let trials = value(drawn);
        assert!(trials.is_integer() && trials >= Rational::ONE, "{drawn:?}");
    }
}

/// Vectors, matrices and polynomials (CAS-373): giac's default entries in
/// `[−99, 99]`, `[0, m)`, a range's reals, a law's variates; monic
/// polynomials of the degree asked.
#[test]
fn vectors_matrices_polynomials() {
    let ctx = Context::new();
    let entries = |text: &str| -> Vec<Rational> {
        let drawn = eval(&ctx, text);
        drawn
            .args()
            .iter()
            .flat_map(|row| {
                if row.head() == Some(&Expr::from(Symbol::LIST)) {
                    row.args().to_vec()
                } else {
                    vec![row.clone()]
                }
            })
            .map(|x| value(&x))
            .collect()
    };
    for (input, count, low, high) in [
        ("RandomVector[5]", 5, -99, 99),
        ("RandomVector[4, 6]", 4, 0, 5),
        ("RandomVector[2, Binomial, 10, 1/2]", 2, 0, 10),
        ("RandomMatrix[2, 3]", 6, -99, 99),
        ("RandomMatrix[2]", 4, -99, 99),
        ("RandomMatrix[2, 2, 6]", 4, 0, 5),
    ] {
        let found = entries(input);
        assert_eq!(found.len(), count, "{input}");
        for x in &found {
            assert!(integer_in(x, low, high), "{input}: {x}");
        }
    }
    for x in entries("RandomVector[3, Range[5, 7]]") {
        assert!(q(5, 1) <= x && x <= q(7, 1), "{x}");
    }
    assert_eq!(entries("RandomVector[3, NormalPDF, 0, 1]").len(), 3);
    assert_eq!(eval(&ctx, "RandomVector[0]"), read("List[]"));
    assert_eq!(
        eval(&ctx, "Dimensions[RandomMatrix[2, 3]]"),
        read("List[2, 3]")
    );
    for (input, variable) in [
        ("RandomPolynomial[3, x]", "x"),
        ("RandomPolynomial[x, 3]", "x"),
        ("RandomPolynomial[3]", "x"),
        ("RandomPolynomial[3, y, 5]", "y"),
    ] {
        let p = eval(&ctx, input);
        let degree = eval(&ctx, &format!("Exponent[{p:?}, {variable}]"));
        assert_eq!(degree, Expr::from(3), "{input}: {p:?}");
        let leading = eval(&ctx, &format!("Coefficient[{p:?}, {variable}, 3]"));
        assert_eq!(leading, Expr::from(1), "{input}: {p:?}");
    }
    assert_eq!(eval(&ctx, "RandomPolynomial[0]"), Expr::from(1));
}

/// Permutations of `0 … n − 1` and of lists, and samples at distinct
/// places (CAS-373).
#[test]
fn permutations_and_samples() {
    let ctx = Context::new();
    let sorted = |e: Expr| {
        let mut items = e.args().to_vec();
        items.sort();
        items
    };
    for _ in 0..20 {
        let order = sorted(eval(&ctx, "RandomPermutation[5]"));
        let expected: Vec<Expr> = (0..5).map(Expr::from).collect();
        assert_eq!(order, expected);
        let letters = sorted(eval(&ctx, "RandomPermutation[List[a, b, c]]"));
        assert_eq!(letters, read("List[a, b, c]").args());
        let pair = eval(&ctx, "RandomSample[List[a, b, c, d], 2]");
        assert_eq!(pair.args().len(), 2, "{pair:?}");
        assert_ne!(pair.args()[0], pair.args()[1], "{pair:?}");
    }
    assert_eq!(eval(&ctx, "RandomPermutation[0]"), read("List[]"));
    assert_eq!(eval(&ctx, "RandomSample[List[a, b], 0]"), read("List[]"));
}

/// What stays: no number, a range or a list too small for the count, an
/// integer form given a non-integer.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Random[x]",
        "Random[1/2]",
        "Random[0]",
        "Random[x, 1]",
        "Random[List[]]",
        "Random[4, List[a, b, c]]",
        "Random[5, 1, 3]",
        "Random[1/2, 1, 3]",
        "Random[2, 1/2, 3]",
        "Random[1, 2, 3, 4]",
        "RandomInteger[1/2, 3]",
        "RandomInteger[x]",
        "RandomBinomial[10, 2]",
        "RandomExponential[-1]",
        "RandomPoisson[0]",
        "RandomNormal[0, -1]",
        "RandomNormal[1]",
        "RandomPoisson[x]",
        "RandomVector[x]",
        "RandomVector[-1]",
        "RandomVector[3, 0]",
        "RandomVector[3, f]",
        "RandomMatrix[2, x]",
        "RandomPolynomial[-1]",
        "RandomPolynomial[2, 3]",
        "RandomPermutation[-1]",
        "RandomPermutation[x]",
        "RandomSample[List[a, b], 3]",
        "RandomSample[x, 1]",
    ] {
        let e = read(input);
        assert_eq!(eval(&ctx, input).head(), e.head(), "{input}");
    }
}
