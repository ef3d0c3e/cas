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
