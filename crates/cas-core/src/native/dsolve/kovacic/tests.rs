#![allow(
    clippy::many_single_char_names,
    reason = "the equation's coefficients a, b, c and the polynomial's names"
)]

use crate::context::Context;
use crate::expr::Expr;
use crate::symbol::Symbol;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e:?}"))
}

/// `KovacicSols[{a, b, c}, v]`: a list of `count` solutions, each making
/// `a Y″ + b Y′ + c Y` simplify to 0, or, `exact` false, vanish at
/// `v = 23/7` within 10⁻¹⁰ (radicals `Simplify` leaves). The answer.
fn solved(abc: [&str; 3], v: &str, count: usize, exact: bool) -> Expr {
    let ctx = Context::new();
    let [a, b, c] = abc;
    let input = format!("KovacicSols[List[{a}, {b}, {c}], {v}]");
    let answer = eval(&ctx, &input);
    assert_eq!(
        answer.head(),
        Some(&Expr::from(Symbol::LIST)),
        "{input}: {answer:?}"
    );
    assert_eq!(answer.args().len(), count, "{input}: {answer:?}");
    for y in answer.args() {
        let residual = format!(
            "Plus[Times[{a}, D[D[{y:?}, {v}], {v}]], Times[{b}, D[{y:?}, {v}]], Times[{c}, {y:?}]]"
        );
        if exact {
            let zero = eval(&ctx, &format!("Simplify[{residual}]"));
            assert_eq!(zero, Expr::from(0), "{input}: {y:?}");
        } else {
            let value = format!("N[ReplaceAll[Expand[{residual}], Rule[{v}, 23/7]]]");
            for part in ["Re", "Im"] {
                let at = format!("Less[Abs[{part}[{value}]], 1/10000000000]");
                assert_eq!(eval(&ctx, &at), read("True"), "{input}: {y:?}");
            }
        }
    }
    answer
}

/// Cases 1 and 2 on giac's examples (its aide, probe rows of CAS-584):
/// every solution satisfies its equation.
#[test]
fn liouvillian_solutions() {
    // Case 2, a pole of order 2 at 0: x^(1/4) e^(±2√x).
    solved(
        ["Times[16, Power[x, 2]]", "0", "Plus[3, Times[-16, x]]"],
        "x",
        2,
        true,
    );
    // Case 1, one pole: (x² − 2x)/√(x − 1).
    solved(
        ["Times[36, Power[Plus[-1, x], 2]]", "0", "-27"],
        "x",
        1,
        true,
    );
    // Case 2, two poles: fourth roots of radicals, checked numerically.
    solved(
        [
            "Times[16, Power[Plus[-1, t], 2], Power[t, 2]]",
            "0",
            "Times[3, Plus[1, Times[-1, t], Power[t, 2]]]",
        ],
        "t",
        2,
        false,
    );
    // A y′ term: (x² + 6x + 12) e^(−x)/x².
    solved(["Power[x, 2]", "Power[x, 2]", "-6"], "x", 1, true);
    // No poles, degree 9 and 10.
    solved(["1", "0", "Plus[19, Times[-1, Power[x, 2]]]"], "x", 1, true);
    solved(["1", "0", "Plus[21, Times[-1, Power[x, 2]]]"], "x", 1, true);
    // Poles ±i: (x² − 1)/√(x² + 1), the conjugate powers combined.
    let y = solved(["Power[Plus[1, Power[x, 2]], 2]", "0", "3"], "x", 1, true);
    assert_eq!(
        y,
        read("List[Times[Plus[-1, Power[x, 2]], Power[Plus[1, Power[x, 2]], -1/2]]]")
    );
    // Poles ±√2 (giac answers {}): case 2.
    solved(["Power[Plus[-2, Power[x, 2]], 2]", "0", "-1"], "x", 2, true);
    // giac skips the family of degree 0 (1/x): x², and x³.
    let y = solved(["Power[x, 2]", "0", "-2"], "x", 1, true);
    assert_eq!(y, read("List[Power[x, 2]]"));
    let y = solved(["Power[x, 2]", "0", "-6"], "x", 1, true);
    assert_eq!(y, read("List[Power[x, 3]]"));
}

/// The equation forms, `Equation` and `Equal` (C0016), with the
/// variables or `y[x]`.
#[test]
fn equation_forms() {
    let ctx = Context::new();
    let expected = read("List[Times[x, Exp[Times[-1/2, Power[x, 2]]]]]");
    for input in [
        "KovacicSols[Equation[Derivative[2][y][x], Times[Plus[-3, Power[x, 2]], y[x]]]]",
        "KovacicSols[Equal[Derivative[2][y][x], Times[Plus[-3, Power[x, 2]], y[x]]], x, y]",
        "KovacicSols[Plus[Derivative[2][y][x], Times[Plus[3, Times[-1, Power[x, 2]]], y[x]]], y[x]]",
        "KovacicSols[List[1, 0, Plus[3, Times[-1, Power[x, 2]]]], x]",
    ] {
        assert_eq!(eval(&ctx, input), expected, "{input}");
    }
    let input = "KovacicSols[Equation[Plus[Times[Power[t, 2], Derivative[2][y][t]], \
                 Times[-2, y[t]]], 0], y[t]]";
    assert_eq!(eval(&ctx, input), read("List[Power[t, 2]]"));
}

/// No solution is `{}`; giac's `return false` (desolve.cc:2157, :2183) is 0.
#[test]
fn none_and_false() {
    let ctx = Context::new();
    for (input, answer) in [
        ("KovacicSols[List[1, 0, Times[-1, x]]]", "List[]"),
        (
            "KovacicSols[List[1, 0, Times[-1, Plus[1, Power[x, 3]]]]]",
            "List[]",
        ),
        ("KovacicSols[List[1, 0, -1]]", "List[]"),
        (
            "KovacicSols[List[Power[x, 4], 0, Times[-1, Plus[1, Power[x, 2]]]]]",
            "List[]",
        ),
        ("KovacicSols[List[1, 0, 1]]", "0"),
        ("KovacicSols[List[1, 0, 0]]", "0"),
        (
            "KovacicSols[List[Power[t, 2], t, Plus[-1/4, Power[t, 2]]], t]",
            "0",
        ),
    ] {
        assert_eq!(eval(&ctx, input), read(answer), "{input}");
    }
}

/// Case 3 (n = 4) with a y′ term: a root ω of the polynomial is y′/y, so
/// `ω′ + ω² + p ω + q = 0` holds modulo the polynomial, `ω′ = −F_t/F_ω`.
/// giac's polynomial, with its coefficients reversed (desolve.cc:2417),
/// fails this.
#[test]
fn case_three_riccati() {
    let ctx = Context::new();
    let a = "Times[48, t, Plus[1, t], Plus[-4, Times[5, t]]]";
    let b = "Times[8, Plus[16, Times[25, t]], Plus[-2, t]]";
    let c = "Times[-1, Plus[68, Times[5, t]]]";
    let f = eval(&ctx, &format!("KovacicSols[List[{a}, {b}, {c}], t]"));
    let w = read("omega_");
    assert_eq!(
        ctx.eval(&Expr::call(Symbol::EXPONENT, [f.clone(), w]))
            .unwrap(),
        Expr::from(4),
        "{f:?}"
    );
    let p = format!("Times[{b}, Power[{a}, -1]]");
    let q = format!("Times[{c}, Power[{a}, -1]]");
    let g = format!(
        "Plus[Times[-1, D[{f:?}, t]], \
         Times[D[{f:?}, omega_], Plus[Power[omega_, 2], Times[{p}, omega_], {q}]]]"
    );
    let remainder = format!(
        "Together[PolynomialRemainder[Numerator[Together[{g}]], \
         Numerator[Together[{f:?}]], omega_]]"
    );
    assert_eq!(eval(&ctx, &remainder), Expr::from(0), "{f:?}");
}

/// A list holding `/` reads as written (giac quotes it and answers 0).
#[test]
fn quotients_in_lists() {
    let ctx = Context::new();
    let input = "KovacicSols[List[1, 0, Times[-2, Power[x, -2]]]]";
    assert_eq!(eval(&ctx, input), read("List[Power[x, 2]]"));
}

/// giac's errors, a pole past the port and a degree past the limits: the
/// call stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // Bad Argument Value: not second order, a's 0, four entries.
        "KovacicSols[Plus[Derivative[1][y][x], y[x]]]",
        "KovacicSols[List[0, 1, 1]]",
        "KovacicSols[List[Power[x, 2], 0, -2, 1]]",
        // Not linear, not rational in x.
        "KovacicSols[Equation[Derivative[2][y][x], Power[y[x], 2]]]",
        "KovacicSols[Equation[Derivative[2][y][x], Times[Sin[x], y[x]]]]",
        // A pole of x³ − 2.
        "KovacicSols[List[Plus[-2, Power[x, 3]], 0, -1]]",
        // y″ = (x² − 2001) y: a family of degree 1000.
        "KovacicSols[List[1, 0, Plus[2001, Times[-1, Power[x, 2]]]]]",
    ] {
        let answer = eval(&ctx, input);
        assert_eq!(
            answer.head(),
            Some(&Expr::from(Symbol::KOVACICSOLS)),
            "{input}: {answer:?}"
        );
    }
}
