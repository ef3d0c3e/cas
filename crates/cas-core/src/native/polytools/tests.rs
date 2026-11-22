//! Tests of the polynomial tools (CAS-DOC-29): giac's partial fractions,
//! compared after putting the difference together, coefficient lists,
//! trailing terms, root counts and real roots in balls.

use cas_num::{Integer, Rational};

use crate::expr::View;
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a − b` is 0 once expanded and put together.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// giac's `cpartfrac`, `propfrac`, `f2nd` and `ratnormal`: the same value.
#[test]
fn complex_and_proper_fractions() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in [
        // giac's i/2/(x+i)+1/(2*i)/(x-i).
        (
            "ComplexApart[Power[Plus[Power[x, 2], 1], -1]]",
            "Plus[Times[Complex[0, 1/2], Power[Plus[x, Complex[0, 1]], -1]], Times[Complex[0, -1/2], Power[Plus[x, Complex[0, -1]], -1]]]",
        ),
        // giac's x-1/2/(x+i)-1/2/(x-i).
        (
            "ComplexApart[Times[Power[x, 3], Power[Plus[Power[x, 2], 1], -1]]]",
            "Plus[x, Times[-1/2, Power[Plus[x, Complex[0, 1]], -1]], Times[-1/2, Power[Plus[x, Complex[0, -1]], -1]]]",
        ),
        // A negative discriminant: roots ±i√2.
        (
            "ComplexApart[Power[Plus[Power[x, 2], 2], -1]]",
            "Plus[Times[Complex[0, -1/4], Power[2, 1/2], Power[Plus[x, Times[Complex[0, -1], Power[2, 1/2]]], -1]], Times[Complex[0, 1/4], Power[2, 1/2], Power[Plus[x, Times[Complex[0, 1], Power[2, 1/2]]], -1]]]",
        ),
        // giac's x-x/(x^2+1) and x^2-y*x+y^2+(-y^3+y)/(x+y).
        (
            "ProperFraction[Times[Power[x, 3], Power[Plus[Power[x, 2], 1], -1]]]",
            "Plus[x, Times[-1, x, Power[Plus[Power[x, 2], 1], -1]]]",
        ),
        (
            "ProperFraction[Times[Plus[Power[x, 3], y], Power[Plus[x, y], -1]]]",
            "Plus[Power[x, 2], Times[-1, y, x], Power[y, 2], Times[Plus[y, Times[-1, Power[y, 3]]], Power[Plus[x, y], -1]]]",
        ),
        ("ProperFraction[28/12]", "7/3"),
        // giac's [x^2+1,x], [1,2].
        (
            "NumeratorDenominator[Plus[x, Power[x, -1]]]",
            "List[Plus[Power[x, 2], 1], x]",
        ),
        ("NumeratorDenominator[1/2]", "List[1, 2]"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = read(expected);
        let ok = if expected.head() == Some(&Expr::from(Symbol::LIST)) {
            found == ctx.eval(&expected).unwrap()
        } else {
            same(&ctx, &found, &expected)
        };
        if !ok {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // x⁴ + 1 over ℚ(i): x² ∓ i stay whole, two terms.
    let terms = ctx
        .eval(&read("ComplexApart[Power[Plus[Power[x, 4], 1], -1]]"))
        .unwrap();
    assert_eq!(terms.args().len(), 2, "{terms:?}");
}

/// giac's `partfrac`: the same value, and no fraction left whole.
#[test]
fn partial_fractions() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in [
        (
            "Apart[Times[x, Power[Plus[2, Times[-3, x], Power[x, 2]], -1]]]",
            "Plus[Times[2, Power[Plus[x, -2], -1]], Times[-1, Power[Plus[x, -1], -1]]]",
        ),
        (
            "Apart[Times[Plus[1, Power[x, 2]], Power[Plus[Power[x, 3], Times[-1, x]], -1]]]",
            "Plus[Times[-1, Power[x, -1]], Power[Plus[x, -1], -1], Power[Plus[x, 1], -1]]",
        ),
        (
            "Apart[Times[Plus[Power[x, 3], 1], Power[Times[x, Power[Plus[x, -1], 2]], -1]]]",
            "Plus[1, Power[x, -1], Times[2, Power[Plus[x, -1], -2]], Power[Plus[x, -1], -1]]",
        ),
        (
            "Apart[Power[Plus[Power[x, 3], 1], -1]]",
            "Plus[Times[1/3, Power[Plus[x, 1], -1]], Times[1/3, Plus[2, Times[-1, x]], Power[Plus[1, Times[-1, x], Power[x, 2]], -1]]]",
        ),
        (
            "Apart[Times[Plus[Times[2, x], 3], Power[Times[Power[Plus[x, 1], 2], Plus[Power[x, 2], 1]], -1]]]",
            "Plus[Times[1/2, Power[Plus[x, 1], -2]], Times[3/2, Power[Plus[x, 1], -1]], Times[1/2, Plus[2, Times[-3, x]], Power[Plus[Power[x, 2], 1], -1]]]",
        ),
        (
            "Apart[Times[Plus[Power[x, 4], 1], Power[Plus[Power[x, 2], -1], -1]]]",
            "Plus[Power[x, 2], 1, Power[Plus[x, -1], -1], Times[-1, Power[Plus[x, 1], -1]]]",
        ),
        // Square roots, as giac's 1/(2*sqrt(2))/(x-sqrt(2))-….
        (
            "Apart[Power[Plus[Power[x, 2], -2], -1]]",
            "Plus[Times[Power[Times[2, Power[2, 1/2]], -1], Power[Plus[x, Times[-1, Power[2, 1/2]]], -1]], Times[-1, Power[Times[2, Power[2, 1/2]], -1], Power[Plus[x, Power[2, 1/2]], -1]]]",
        ),
        // Symbols as coefficients.
        (
            "Apart[Power[Times[Plus[x, Times[-1, a]], Plus[x, Times[-1, b]]], -1], x]",
            "Plus[Times[Power[Plus[a, Times[-1, b]], -1], Power[Plus[x, Times[-1, a]], -1]], Times[-1, Power[Plus[a, Times[-1, b]], -1], Power[Plus[x, Times[-1, b]], -1]]]",
        ),
        ("Apart[Plus[1, Power[x, 2]]]", "Plus[1, Power[x, 2]]"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        if found.contains(&Expr::from(Symbol::APART)) || !same(&ctx, &found, &read(expected)) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // Each term's denominator is a factor's power: x/(x²−3x+2) has two.
    let terms = ctx
        .eval(&read(
            "Apart[Times[x, Power[Plus[2, Times[-3, x], Power[x, 2]], -1]]]",
        ))
        .unwrap();
    assert_eq!(terms.head(), Some(&Expr::from(Symbol::PLUS)), "{terms:?}");
    assert_eq!(terms.args().len(), 2, "{terms:?}");
}

/// `symb2poly`, `poly2symb`, `tcoeff`, `valuation` as giac answers them.
#[test]
fn coefficients() {
    let ctx = Context::new();
    for (input, expected) in [
        (
            "CoefficientList[Plus[Power[x, 3], Times[-2, x], 1]]",
            "List[1, 0, -2, 1]",
        ),
        (
            "CoefficientList[Plus[Times[Power[x, 2], y], y], x]",
            "List[y, 0, y]",
        ),
        ("CoefficientList[5]", "List[5]"),
        (
            "FromCoefficients[List[2, 0, 1], t]",
            "Plus[1, Times[2, Power[t, 2]]]",
        ),
        (
            "FromCoefficients[List[1, 0, -2, 1]]",
            "Plus[1, Times[-2, x], Power[x, 3]]",
        ),
        (
            "TrailingCoefficient[Plus[Times[3, Power[x, 4]], Times[2, Power[x, 2]]]]",
            "2",
        ),
        (
            "Valuation[Plus[Times[3, Power[x, 4]], Times[2, Power[x, 2]]]]",
            "2",
        ),
        ("Valuation[0]", "Times[-1, Infinity]"),
    ] {
        assert_eq!(
            ctx.eval(&read(input)).unwrap(),
            ctx.eval(&read(expected)).unwrap(),
            "{input}"
        );
    }
    for input in ["CoefficientList[Power[x, -1]]", "CoefficientList[Sin[x]]"] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap().head(), e.head(), "{input}");
    }
}

/// `sturmab` and `realroot`: distinct real roots in an interval, and the
/// roots in balls with their multiplicities.
#[test]
fn real_roots() {
    let ctx = Context::new();
    for (input, expected) in [
        ("CountRoots[Plus[Power[x, 3], Times[-2, x]], x, -2, 2]", 3),
        ("CountRoots[Plus[Power[x, 2], 1], x, -5, 5]", 0),
        (
            "CountRoots[Times[Power[Plus[x, -1], 2], Plus[x, 2]], x, -3, 3]",
            2,
        ),
        ("CountRoots[Plus[Power[x, 2], -2], x, 0, Infinity]", 1),
        (
            "CountRoots[Plus[Power[x, 2], -2], x, Times[-1, Infinity], Infinity]",
            2,
        ),
        (
            "CountRoots[Plus[Power[x, 2], -2], x, DirectedInfinity[-1], 0]",
            1,
        ),
    ] {
        assert_eq!(
            ctx.eval(&read(input)).unwrap(),
            Expr::from(expected),
            "{input}"
        );
    }
    let q = |n: i64| Rational::from(Integer::from(n));
    // x^4 − 5x^2 + 4 = (x² − 1)(x² − 4): −2, −1, 1, 2.
    let found = ctx
        .eval(&read(
            "RealRoots[Plus[Power[x, 4], Times[-5, Power[x, 2]], 4]]",
        ))
        .unwrap();
    let roots: Vec<(Expr, Expr)> = found
        .args()
        .iter()
        .map(|pair| (pair.args()[0].clone(), pair.args()[1].clone()))
        .collect();
    assert_eq!(roots.len(), 4, "{found:?}");
    for ((ball, m), r) in roots.iter().zip([-2, -1, 1, 2]) {
        let View::Ball(b) = ball.view() else {
            panic!("{found:?}");
        };
        assert!(
            b.contains_rational(&q(r)) && *m == Expr::from(1),
            "{found:?}"
        );
    }
    // (x − 1)²(x + 2): −2 once, 1 twice.
    let found = ctx
        .eval(&read("RealRoots[Times[Power[Plus[x, -1], 2], Plus[x, 2]]]"))
        .unwrap();
    assert_eq!(found.args().len(), 2, "{found:?}");
    assert_eq!(found.args()[1].args()[1], Expr::from(2), "{found:?}");
}

/// giac's `VAS`, `VAS_positive`, `posubLMQ`, `poslbdLMQ`, `sturmseq`,
/// `sturm` with bounds, `rationalroot` and `crationalroot` (CAS-538,
/// CAS-539), on giac's examples.
#[test]
fn isolation_and_sequences() {
    let ctx = Context::new();
    let cubic = "Plus[Power[x, 3], Times[-7, x], 7]";
    let mut wrong = Vec::new();
    for (input, expected) in [
        (
            format!("RealRootIsolation[{cubic}]"),
            "List[List[-4, 0], List[1, 3/2], List[3/2, 2]]",
        ),
        (
            format!("PositiveRootIsolation[{cubic}]"),
            "List[List[1, 3/2], List[3/2, 2]]",
        ),
        (
            "RealRootIsolation[Times[Power[Plus[x, -1], 2], Plus[x, -2]]]".to_owned(),
            "List[1, List[1, 4]]",
        ),
        (
            "RealRootIsolation[Plus[Power[x, 3], Times[-1, x]]]".to_owned(),
            "List[List[-2, 0], 0, List[0, 2]]",
        ),
        (
            "RealRootIsolation[Plus[Power[y, 2], -2], y]".to_owned(),
            "List[List[-3, 0], List[0, 3]]",
        ),
        ("RealRootIsolation[5]".to_owned(), "List[]"),
        (format!("PositiveRootUpperBound[{cubic}]"), "4"),
        (format!("PositiveRootLowerBound[{cubic}]"), "1/2"),
        ("PositiveRootUpperBound[Plus[x, -2]]".to_owned(), "5"),
        ("PositiveRootLowerBound[Plus[x, -100]]".to_owned(), "32"),
        (
            "PositiveRootLowerBound[Plus[Power[x, 3], x, 1]]".to_owned(),
            "Infinity",
        ),
        (
            format!("SturmSequence[{cubic}]"),
            "List[1, List[List[1, 0, -7, 7], List[3, 0, -7], List[42, -63], 49], 1]",
        ),
        (
            "SturmSequence[Plus[Times[2, Power[x, 2]], -4]]".to_owned(),
            "List[2, List[List[1, 0, -2], List[2, 0], 8], 1]",
        ),
        (
            "SturmSequence[Power[Plus[Power[x, 2], -2], 2]]".to_owned(),
            "List[1, 1]",
        ),
        (
            "SturmSequence[Times[Power[Plus[x, -1], 2], Plus[x, 2]]]".to_owned(),
            "List[1, List[List[1, 2], 1], 1]",
        ),
        (
            "SturmSequence[Times[Plus[x, -1], Power[Plus[x, 2], -1]]]".to_owned(),
            "List[1, List[List[1, -1], 1], 1, List[List[1, 2], 1]]",
        ),
        ("SturmSequence[-5]".to_owned(), "List[-5, 1]"),
        (format!("SturmSequence[{cubic}, x, -5, 5]"), "3"),
        (format!("SturmSequence[{cubic}, x, 1, 5]"), "2"),
        (
            "SturmSequence[Plus[Power[x, 2], -1], x, -1, 1]".to_owned(),
            "1",
        ),
        (
            "SturmSequence[Plus[Power[x, 2], -1], x, 1, 2]".to_owned(),
            "0",
        ),
        (
            "SturmSequence[Plus[Power[x, 2], -4], x, -1, 1]".to_owned(),
            "-1",
        ),
        (
            "SturmSequence[Plus[Power[x, 2], -1], x, 0, 0]".to_owned(),
            "-1",
        ),
        (
            "RationalRoots[Plus[Times[2, Power[x, 2]], Times[-3, x], 1]]".to_owned(),
            "List[1, 1/2]",
        ),
        (
            "RationalRoots[Times[Power[Plus[x, -1], 2], Plus[Times[2, x], 3]]]".to_owned(),
            "List[-3/2, 1]",
        ),
        ("RationalRoots[Plus[Power[x, 2], -2]]".to_owned(), "List[]"),
        ("RationalRoots[3]".to_owned(), "List[]"),
        (
            "ComplexRationalRoots[Plus[Power[x, 2], 1]]".to_owned(),
            "List[Complex[0, -1], Complex[0, 1]]",
        ),
        (
            "ComplexRationalRoots[Times[Plus[Power[x, 2], 1], Plus[Times[2, x], -1]]]".to_owned(),
            "List[1/2, Complex[0, -1], Complex[0, 1]]",
        ),
    ] {
        let got = ctx.eval(&read(&input));
        if got.as_ref() != Ok(&read(expected)) {
            wrong.push(format!("{input}: {got:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `complexroot`: every root with its multiplicity, a linear
/// factor's exact, the others in balls; the rectangle keeps those in it.
#[test]
fn complex_roots() {
    let ctx = Context::new();
    let i = cas_num::Complex::new(Rational::ZERO, Rational::ONE);
    // (x − 1)²(x² + 1): ±i in balls, then 1 twice.
    let found = ctx
        .eval(&read(
            "ComplexRoots[Times[Power[Plus[x, -1], 2], Plus[Power[x, 2], 1]], 1/10]",
        ))
        .unwrap();
    assert_eq!(found.args().len(), 3, "{found:?}");
    assert_eq!(found.args()[2], read("List[1, 2]"), "{found:?}");
    let mut signs = Vec::new();
    for pair in &found.args()[..2] {
        let View::ComplexBall(z) = pair.args()[0].view() else {
            panic!("{found:?}");
        };
        assert_eq!(pair.args()[1], Expr::from(1), "{found:?}");
        let plus = z.contains_complex(&i);
        assert!(plus || z.contains_complex(&-&i), "{found:?}");
        signs.push(plus);
    }
    assert!(signs[0] != signs[1], "{found:?}");
    // The rectangle from 0 to 2 + 2i.
    let found = ctx
        .eval(&read(
            "ComplexRoots[Times[Power[Plus[x, -1], 2], Plus[x, 2]], 1/10, 0, Complex[2, 2]]",
        ))
        .unwrap();
    assert_eq!(found, read("List[List[1, 2]]"));
    let found = ctx
        .eval(&read(
            "ComplexRoots[Plus[Power[x, 2], 1], 1/1000, 0, Complex[2, 2]]",
        ))
        .unwrap();
    assert_eq!(found.args().len(), 1, "{found:?}");
    let View::ComplexBall(z) = found.args()[0].args()[0].view() else {
        panic!("{found:?}");
    };
    assert!(z.contains_complex(&i), "{found:?}");
}
