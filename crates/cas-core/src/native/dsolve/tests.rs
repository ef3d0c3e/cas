use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Where residuals are evaluated: the variables and the constants.
const POINT: &str = "List[Rule[x, 2/9], Rule[t, 5/7], Rule[a, 3/4], \
                     Rule[c_0, 7/5], Rule[c_1, -2/3], Rule[c_2, 5/11]]";

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
    let residual = template.replace('Y', &format!("{solution:?}"));
    let at = format!("Less[Abs[N[ReplaceAll[{residual}, {POINT}]]], 1/10000000000]");
    assert_eq!(
        ctx.eval(&read(&at)).unwrap(),
        read("True"),
        "{input}: {solution:?}"
    );
    for k in 0..order {
        let c = read(&format!("c_{k}"));
        assert!(solution.contains(&c), "{input}: {solution:?} lacks c_{k}");
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

/// What no method solves, or what is not an equation in `y`, stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // Airy's equation.
        "DSolve[Equation[D[y[x], x, x], Times[x, y[x]]], y[x]]",
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
