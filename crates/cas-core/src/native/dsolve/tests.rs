use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Where residuals are evaluated: the variables and the constants.
const POINT: &str = "List[Rule[x, 2/9], Rule[t, 5/7], Rule[a, 3/4], \
                     Rule[c_0, 7/5], Rule[c_1, -2/3], Rule[c_2, 5/11], \
                     Rule[c_3, 3/13], Rule[c_4, -9/7], Rule[c_5, 4/3]]";

/// `DSolve[template = 0, y[v]]`, `Y` standing for `y[v]` in the template,
/// and the template at the answer: 0 at [`POINT`]. The answer, and its
/// constants `c_0 … c_(order−1)`.
fn general(template: &str, v: &str, order: usize) -> Expr {
    let ctx = Context::new();
    let unknown = format!("y[{v}]");
    let input = format!(
        "DSolve[Equation[{}, 0], {unknown}]",
        template.replace('Y', &unknown)
    );
    let solution = ctx.eval(&read(&input)).unwrap();
    assert!(
        solution.head() != Some(&Expr::from(crate::symbol::Symbol::DSOLVE)),
        "{input}: stays"
    );
    // A nonlinear equation's branches each satisfy it.
    let branches = if solution.head() == Some(&Expr::from(crate::symbol::Symbol::LIST)) {
        solution.args().to_vec()
    } else {
        vec![solution.clone()]
    };
    assert!(!branches.is_empty(), "{input}: no branch");
    for branch in &branches {
        // Each part apart: the modulus of a complex ball about 0 is
        // unbounded.
        let residual = template.replace('Y', &format!("{branch:?}"));
        // Expanded first: Times[-1, Plus[…]] is not distributed, which
        // leaves integrals that stay standing where they would cancel.
        let value = format!("N[ReplaceAll[Expand[{residual}], {POINT}]]");
        for part in ["Re", "Im"] {
            let at = format!("Less[Abs[{part}[{value}]], 1/10000000000]");
            assert_eq!(
                ctx.eval(&read(&at)).unwrap(),
                read("True"),
                "{input}: {branch:?}"
            );
        }
        for k in 0..order {
            let c = read(&format!("c_{k}"));
            assert!(branch.contains(&c), "{input}: {branch:?} lacks c_{k}");
        }
    }
    solution
}

/// First-order linear equations (CAS-259): their solutions, some exactly.
#[test]
fn first_order_linear() {
    let cases = [
        ("Plus[D[Y, x], Times[-1, Y]]", "x"),
        ("Plus[D[Y, x], Y, Times[-1, x]]", "x"),
        ("Plus[D[Y, x], Times[-2, Y, Power[x, -1]]]", "x"),
        (
            "Plus[D[Y, x], Times[-1, Y, Power[x, -1]], Times[-1, x]]",
            "x",
        ),
        (
            "Plus[Times[2, x, Y], Times[Plus[1, Power[x, 2]], D[Y, x]]]",
            "x",
        ),
        ("Plus[D[Y, t], Times[2, t, Y]]", "t"),
        ("Plus[D[Y, x], Times[2, x, Y], Times[-1, x]]", "x"),
        ("Plus[D[Y, x], Times[-1, Sin[x]]]", "x"),
        ("Plus[D[Y, x], Times[-1, a, Y]]", "x"),
        ("Plus[Times[t, D[Y, t]], Times[-1, Y]]", "t"),
        // An integral that stays: giac's Ei.
        ("Plus[D[Y, x], Y, Times[-1, Power[x, -1]]]", "x"),
    ];
    for (template, v) in cases {
        general(template, v, 1);
    }
    let ctx = Context::new();
    for (input, expected) in [
        (
            "DSolve[Equation[D[y[x], x], y[x]], y[x]]",
            "Times[c_0, Power[E, x]]",
        ),
        (
            "DSolve[Equation[D[y[x], x], Times[2, y[x], Power[x, -1]]], y[x]]",
            "Times[c_0, Power[x, 2]]",
        ),
        (
            "DSolve[Equation[Plus[Times[2, x, y[x]], Times[Plus[1, Power[x, 2]], D[y[x], x]]], 0], y[x]]",
            "Times[c_0, Power[Plus[1, Power[x, 2]], -1]]",
        ),
        // giac's order, the variable first, and the unknown alone.
        (
            "DSolve[Equation[D[y[t], t], y[t]], t, y]",
            "Times[c_0, Power[E, t]]",
        ),
        (
            "DSolve[Equation[D[y[t], t], y[t]], y]",
            "Times[c_0, Power[E, t]]",
        ),
        // A constant of the input is skipped.
        (
            "DSolve[Equation[D[y[x], x], Times[c_0, y[x]]], y[x]]",
            "Times[c_1, Power[E, Times[c_0, x]]]",
        ),
    ] {
        assert_eq!(
            format!("{:?}", ctx.eval(&read(input)).unwrap()),
            expected,
            "{input}"
        );
    }
}

/// First-order equations that are not linear (CAS-260): lists of branches.
#[test]
fn first_order_nonlinear() {
    let cases = [
        // Separable.
        "Plus[D[Y, x], Times[-1, Power[Y, 2]]]",
        "Plus[D[Y, x], Times[-1, x, Power[Y, 2]]]",
        "Plus[D[Y, x], Times[-1, x, Power[Y, -1]]]",
        "Plus[D[Y, x], -1, Times[-1, Power[Y, 2]]]",
        "Plus[D[Y, x], Times[-1, x, Power[E, Y]]]",
        "Plus[D[Y, x], Times[-1, Power[E, Plus[x, Y]]]]",
        // Bernoulli.
        "Plus[D[Y, x], Y, Times[-1, Power[Y, 2]]]",
        "Plus[D[Y, x], Times[Y, Power[x, -1]], Times[-1, Power[x, 2], Power[Y, 3]]]",
        "Plus[D[Y, x], Times[-1, Plus[Power[x, 2], Power[Y, 2]], Power[x, -1], Power[Y, -1]]]",
        // Homogeneous.
        "Plus[D[Y, x], Times[Plus[Times[3, x, Y], Power[Y, 2]], Power[Plus[Power[x, 2], Times[x, Y]], -1]]]",
        // Exact, and exact once multiplied by x.
        "Plus[Times[3, Power[x, 2]], Y, Times[Plus[x, Times[2, Y]], D[Y, x]]]",
        "Plus[Times[2, Power[Y, 2]], Times[3, x], Times[2, x, Y, D[Y, x]]]",
    ];
    for template in cases {
        general(template, "x", 1);
    }
    // giac's constant, H(x) − G(y) = c_0: y′ = y² is −1/(x − c_0).
    let ctx = Context::new();
    let found = ctx
        .eval(&read("DSolve[Equation[D[y[x], x], Power[y[x], 2]], y[x]]"))
        .unwrap();
    let expected = ctx
        .eval(&read("List[Times[-1, Power[Plus[Times[-1, c_0], x], -1]]]"))
        .unwrap();
    assert_eq!(found, expected);
}

/// Linear equations with constant coefficients (CAS-261): the basis from
/// the characteristic roots, real forms of conjugate pairs.
#[test]
fn constant_coefficients() {
    let cases = [
        ("Plus[D[Y, x, x], Y]", 2),
        ("Plus[D[Y, x, x], Times[2, D[Y, x]], Y]", 2),
        ("Plus[D[Y, x, x], D[Y, x], Y]", 2),
        ("Plus[D[Y, x, x, x], Times[-1, Y]]", 3),
        ("Plus[D[Y, x, x], Times[-1, Y]]", 2),
        ("Plus[D[Y, x, x, x, x], Times[-1, Y]]", 4),
        ("Plus[D[Y, x, x], Times[Power[a, 2], Y]]", 2),
        ("Plus[D[Y, x, x], Times[-1, a, Y]]", 2),
        ("Plus[D[Y, x, x], D[Y, x]]", 2),
        ("D[Y, x, x]", 2),
        (
            "Plus[D[Y, x, x, x], Times[-3, D[Y, x, x]], Times[3, D[Y, x]], Times[-1, Y]]",
            3,
        ),
    ];
    for (template, order) in cases {
        general(template, "x", order);
    }
    let ctx = Context::new();
    for (input, expected) in [
        (
            "DSolve[Equation[Plus[D[y[x], x, x], y[x]], 0], y[x]]",
            "Plus[Times[c_0, Cos[x]], Times[c_1, Sin[x]]]",
        ),
        (
            "DSolve[Equation[D[y[x], x, x], 0], y[x]]",
            "Plus[c_0, Times[c_1, x]]",
        ),
    ] {
        assert_eq!(
            format!("{:?}", ctx.eval(&read(input)).unwrap()),
            expected,
            "{input}"
        );
    }
}

/// Right sides (CAS-262): undetermined coefficients, with the resonance's
/// power of `x`, and variation of parameters for the rest.
#[test]
fn right_sides() {
    let cases = [
        "Plus[D[Y, x, x], Times[-3, D[Y, x]], Times[2, Y], Times[-1, Power[E, x]]]",
        "Plus[D[Y, x, x], Times[4, Y], Times[-1, Sin[x]]]",
        "Plus[D[Y, x, x], Times[-1, x]]",
        "Plus[D[Y, x, x], Times[-1, Y], Times[-1, Power[x, 2]]]",
        "Plus[D[Y, x, x], Y, Times[-1, Cos[x]]]",
        "Plus[D[Y, x, x], Y, Times[-1, x, Power[E, x]]]",
        "Plus[D[Y, x, x], Times[-2, D[Y, x]], Y, Times[-1, Power[E, x]]]",
        "Plus[D[Y, x, x], Y, Times[-1, Sin[x]], Times[-1, x]]",
        "Plus[D[Y, x, x], Y, Times[-1, Power[E, Plus[1, Times[2, x]]]]]",
        "Plus[D[Y, x, x, x], Times[-1, D[Y, x]], Times[-1, Sin[x], Cos[x]]]",
        // Variation of parameters.
        "Plus[D[Y, x, x], Y, Times[-1, Tan[x]]]",
        "Plus[D[Y, x, x], Times[-1, Y], Times[-1, Power[x, -1]]]",
    ];
    for template in cases {
        let order = if template.contains("D[Y, x, x, x]") {
            3
        } else {
            2
        };
        general(template, "x", order);
    }
}

/// Euler's equations (CAS-263), some of which giac does not solve.
#[test]
fn euler_equations() {
    let cases = [
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[x, D[Y, x]], Times[-1, Y]]",
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[3, x, D[Y, x]], Y]",
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[x, D[Y, x]], Y]",
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[-2, Y], Times[-1, x]]",
        // Up to a common factor: y″ + y′/x − y/x² = 0.
        "Plus[D[Y, x, x], Times[D[Y, x], Power[x, -1]], Times[-1, Y, Power[x, -2]]]",
        "Plus[Times[Power[x, 3], D[Y, x, x, x]], Times[-3, Power[x, 2], D[Y, x, x]], Times[6, x, D[Y, x]], Times[-6, Y]]",
    ];
    for template in cases {
        let order = if template.contains("D[Y, x, x, x]") {
            3
        } else {
            2
        };
        general(template, "x", order);
    }
}

/// Conditions (CAS-264): the constants fixed, free where they are fewer
/// than the order, `{}` when no solution meets them.
#[test]
fn conditions() {
    let ctx = Context::new();
    let second = "Equation[Plus[D[y[x], x, x], y[x]], 0]";
    for (input, expected) in [
        (
            "DSolve[List[Equation[D[y[x], x], y[x]], Equation[y[0], 1]], y[x]]".to_owned(),
            "Power[E, x]",
        ),
        (
            format!(
                "DSolve[List[{second}, Equation[y[0], 0], Equation[Derivative[1][y][0], 1]], y[x]]"
            ),
            "Sin[x]",
        ),
        (
            format!(
                "DSolve[List[{second}, Equation[y[0], 0], Equation[y[Times[1/2, Pi]], 1]], y[x]]"
            ),
            "Sin[x]",
        ),
        (
            format!("DSolve[List[{second}, Equation[y[0], 0], Equation[y[Pi], 1]], y[x]]"),
            "List[]",
        ),
        (
            format!("DSolve[List[{second}, Equation[y[0], 1]], y[x]]"),
            "Plus[Cos[x], Times[c_1, Sin[x]]]",
        ),
        (
            "DSolve[List[Equation[D[y[x], x], Power[y[x], 2]], Equation[y[0], 1]], y[x]]"
                .to_owned(),
            "List[Times[-1, Power[Plus[-1, x], -1]]]",
        ),
        // The unknown alone: its argument from the equation, not y[0].
        (
            "DSolve[List[Equation[D[y[t], t], y[t]], Equation[y[0], 2]], y]".to_owned(),
            "Times[2, Power[E, t]]",
        ),
    ] {
        assert_eq!(
            format!("{:?}", ctx.eval(&read(&input)).unwrap()),
            expected,
            "{input}"
        );
    }
}

/// Higher orders without `y` (CAS-265): `p = y⁽ᵏ⁾` solved, then integrated.
#[test]
fn without_y() {
    let cases = [
        ("Plus[D[Y, x, x], Times[-1, Power[D[Y, x], 2]]]", 2),
        ("Plus[D[Y, x, x], Times[-1, Power[D[Y, x], 2]], -1]", 2),
        ("Plus[D[Y, x, x], Times[-1, D[Y, x], Tan[x]]]", 2),
        ("Plus[D[Y, x, x, x], Times[-1, D[Y, x, x]]]", 3),
        ("Plus[Times[x, D[Y, x, x]], Times[-1, D[Y, x]]]", 2),
    ];
    for (template, order) in cases {
        general(template, "x", order);
    }
}

/// What no method solves, or what is not an equation in `y`, stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // y at another argument, or not differentiated.
        "DSolve[Equation[D[y[x], x], y[Plus[x, 1]]], y[x]]",
        "DSolve[Equation[y[x], x], y[x]]",
        // Not an unknown function.
        "DSolve[Equation[D[y[x], x], y[x]], Sin[x]]",
    ] {
        // The arguments evaluate: D makes Derivative.
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            found.head(),
            Some(&Expr::from(crate::symbol::Symbol::DSOLVE)),
            "{input}: {found:?}"
        );
    }
}

/// Bessel's and Airy's equations (CAS-385): their solutions satisfy them,
/// and two are giac-free forms: `c_0 Ai(x) + c_1 Bi(x)` and `c_0 J_1(x) +
/// c_1 Y_1(x)`.
#[test]
fn bessel_and_airy() {
    let cases = [
        // Airy's y″ = k x y.
        "Plus[D[Y, x, x], Times[-1, x, Y]]",
        "Plus[D[Y, x, x], Times[-8, x, Y]]",
        "Plus[D[Y, x, x], Times[x, Y]]",
        // Bessel's, and the modified one.
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[x, D[Y, x]], Times[Plus[Power[x, 2], -1], Y]]",
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[x, D[Y, x]], \
         Times[Plus[Times[4, Power[x, 2]], -9], Y]]",
        "Plus[Times[Power[x, 2], D[Y, x, x]], Times[x, D[Y, x]], Times[-1, Plus[Power[x, 2], 4], Y]]",
        "Plus[D[Y, x, x], Times[Power[x, -1], D[Y, x]], Y]",
    ];
    for template in cases {
        general(template, "x", 2);
    }
    let ctx = Context::new();
    for (input, expected) in [
        (
            "DSolve[Equation[D[y[x], x, x], Times[x, y[x]]], y[x]]",
            "Plus[Times[c_0, AiryAi[x]], Times[c_1, AiryBi[x]]]",
        ),
        (
            "DSolve[Equation[Plus[Times[Power[x, 2], D[y[x], x, x]], Times[x, D[y[x], x]], \
             Times[Plus[Power[x, 2], -1], y[x]]], 0], y[x]]",
            "Plus[Times[c_0, BesselJ[1, x]], Times[c_1, BesselY[1, x]]]",
        ),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
}
