//! Tests of systems of differential equations (CAS-283 …): residuals at a
//! point, the constants as initial values, and some forms.

use crate::context::Context;
use crate::expr::Expr;
use crate::symbol::Symbol;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Where residuals are evaluated: the variable, a parameter and the
/// constants.
const POINT: &str = "List[Rule[x, 2/9], Rule[a, 3/4], \
                     Rule[c_0, 7/5], Rule[c_1, -2/3], Rule[c_2, 5/11], \
                     Rule[c_3, 3/13], Rule[c_4, -9/7], Rule[c_5, 4/3]]";

/// The unknowns' letters in templates, and their names.
const UNKNOWNS: [(char, &str); 3] = [('Y', "y"), ('Z', "z"), ('W', "w")];

/// `template` with each letter of [`UNKNOWNS`] as `values`' entry.
fn instance(template: &str, values: &[String]) -> String {
    let mut out = template.to_owned();
    for ((letter, _), value) in UNKNOWNS.iter().zip(values) {
        out = out.replace(*letter, value);
    }
    out
}

/// `DSolve[{template = 0, …}, {y[x], z[x], …}]`, `Y`, `Z` and `W` standing
/// for the unknowns applied: each template is 0 at the answer at
/// [`POINT`], and the constants are the initial values at 0, `orders[j]`
/// of them for the `j`-th unknown, in order. The components.
fn general(templates: &[&str], orders: &[usize]) -> Vec<Expr> {
    let ctx = Context::new();
    let applied: Vec<String> = UNKNOWNS[..orders.len()]
        .iter()
        .map(|(_, y)| format!("{y}[x]"))
        .collect();
    let equations: Vec<String> = templates
        .iter()
        .map(|t| format!("Equation[{}, 0]", instance(t, &applied)))
        .collect();
    let input = format!(
        "DSolve[List[{}], List[{}]]",
        equations.join(", "),
        applied.join(", ")
    );
    let answer = ctx.eval(&read(&input)).unwrap();
    let [solution] = answer.args() else {
        panic!("{input}: {answer:?}");
    };
    assert!(
        answer.head() == Some(&Expr::from(Symbol::LIST))
            && solution.head() == Some(&Expr::from(Symbol::LIST))
            && solution.args().len() == orders.len(),
        "{input}: {answer:?}"
    );
    let components: Vec<String> = solution.args().iter().map(|v| format!("{v:?}")).collect();
    for template in templates {
        let value = format!(
            "N[ReplaceAll[Expand[{}], {POINT}]]",
            instance(template, &components)
        );
        // Each part apart: the modulus of a complex ball about 0 is
        // unbounded.
        for part in ["Re", "Im"] {
            let small = format!("Less[Abs[{part}[{value}]], 1/10000000000]");
            assert_eq!(
                ctx.eval(&read(&small)).unwrap(),
                read("True"),
                "{input}: {template} at {answer:?}"
            );
        }
    }
    let mut k = 0;
    for (component, &order) in components.iter().zip(orders) {
        for m in 0..order {
            let derivative = if m == 0 {
                component.clone()
            } else {
                format!("D[{component}, {}]", vec!["x"; m].join(", "))
            };
            let gap =
                format!("Expand[Plus[ReplaceAll[{derivative}, Rule[x, 0]], Times[-1, c_{k}]]]");
            assert_eq!(
                ctx.eval(&read(&gap)).unwrap(),
                Expr::from(0),
                "{input}: {answer:?}, c_{k}"
            );
            k += 1;
        }
    }
    solution.args().to_vec()
}

/// First-order systems with constant coefficients (CAS-283).
#[test]
fn first_order() {
    let cases: [&[&str]; 8] = [
        // Rotation, giac's [[0,1],[-1,0]].
        &["Plus[D[Y, x], Times[-1, Z]]", "Plus[D[Z, x], Y]"],
        // Real eigenvalues ±√2, and 2 and −1.
        &[
            "Plus[D[Y, x], Times[-1, Y], Times[-1, Z]]",
            "Plus[D[Z, x], Times[-1, Y], Z]",
        ],
        &[
            "Plus[D[Y, x], Times[-3, Y], Times[2, Z]]",
            "Plus[D[Z, x], Times[-2, Y], Times[2, Z]]",
        ],
        // A double eigenvalue with one eigenvector.
        &[
            "Plus[D[Y, x], Times[-1, Y], Times[-1, Z]]",
            "Plus[D[Z, x], Times[-1, Z]]",
        ],
        // Complex eigenvalues 1 ± i.
        &[
            "Plus[D[Y, x], Times[-1, Y], Z]",
            "Plus[D[Z, x], Times[-1, Y], Times[-1, Z]]",
        ],
        // Three unknowns.
        &[
            "Plus[D[Y, x], Times[-1, Z]]",
            "Plus[D[Z, x], Times[-1, W]]",
            "Plus[D[W, x], Times[-1, Y]]",
        ],
        // A parameter.
        &[
            "Plus[D[Y, x], Times[-1, a, Y]]",
            "Plus[D[Z, x], Times[-1, Y]]",
        ],
        // Derivatives mixed in the equations.
        &[
            "Plus[D[Y, x], D[Z, x], Times[-2, Y]]",
            "Plus[D[Y, x], Times[-1, D[Z, x]]]",
        ],
    ];
    for templates in cases {
        let orders = vec![1; templates.len()];
        general(templates, &orders);
    }
}

/// Forcing terms with a transform (CAS-283).
#[test]
fn forced() {
    let cases: [&[&str]; 4] = [
        &[
            "Plus[D[Y, x], Times[-1, Z], Times[-1, x]]",
            "Plus[D[Z, x], Times[-1, Y]]",
        ],
        &[
            "Plus[D[Y, x], Times[-1, Z]]",
            "Plus[D[Z, x], Y, Times[-1, Power[E, x]]]",
        ],
        &[
            "Plus[D[Y, x], Times[-1, Z], Times[-1, x]]",
            "Plus[D[Z, x], Y, Times[-1, Sin[Times[2, x]]]]",
        ],
        &[
            "Plus[D[Y, x], Times[-1, Z], Times[-1, HeavisideTheta[Plus[x, -1]]]]",
            "Plus[D[Z, x], Y]",
        ],
    ];
    for templates in cases {
        general(templates, &[1, 1]);
    }
}

/// Higher orders, and an unknown of order 0 (CAS-283).
#[test]
fn higher_order() {
    general(
        &["Plus[D[Y, x, x], Z]", "Plus[D[Z, x], Times[-1, Y]]"],
        &[2, 1],
    );
    general(
        &[
            "Plus[D[Y, x, x], Times[-1, Z]]",
            "Plus[D[Z, x, x], Times[-1, Y]]",
        ],
        &[2, 2],
    );
    // y′ = z with y + z = 0: z is −y.
    let components = general(&["Plus[D[Y, x], Times[-1, Z]]", "Plus[Y, Z]"], &[1, 0]);
    let ctx = Context::new();
    assert_eq!(
        components,
        [
            ctx.eval(&read("Times[c_0, Power[E, Times[-1, x]]]"))
                .unwrap(),
            ctx.eval(&read("Times[-1, c_0, Power[E, Times[-1, x]]]"))
                .unwrap(),
        ]
    );
}

/// The answers' forms, giac's for the matrix form, and the call's forms.
#[test]
fn forms() {
    let ctx = Context::new();
    let rotation = "List[Equation[D[y[x], x], z[x]], Equation[D[z[x], x], Times[-1, y[x]]]]";
    let expected = "List[List[Plus[Times[c_0, Cos[x]], Times[c_1, Sin[x]]], \
                    Plus[Times[c_1, Cos[x]], Times[-1, c_0, Sin[x]]]]]";
    for call in [
        format!("DSolve[{rotation}, List[y[x], z[x]]]"),
        format!("DSolve[{rotation}, x, List[y, z]]"),
        format!("DSolve[{rotation}, List[y, z]]"),
    ] {
        assert_eq!(
            ctx.eval(&read(&call)).unwrap(),
            ctx.eval(&read(expected)).unwrap(),
            "{call}"
        );
    }
    for (input, expected) in [
        (
            "DSolve[List[Equation[D[y[x], x], Times[2, y[x]]], \
             Equation[D[z[x], x], Times[3, z[x]]]], List[y[x], z[x]]]",
            "List[List[Times[c_0, Power[E, Times[2, x]]], Times[c_1, Power[E, Times[3, x]]]]]",
        ),
        // A constant of the input is skipped.
        (
            "DSolve[List[Equation[D[y[x], x], c_0], Equation[D[z[x], x], y[x]]], \
             List[y[x], z[x]]]",
            "List[List[Plus[c_1, Times[c_0, x]], \
             Plus[c_2, Times[c_1, x], Times[1/2, c_0, Power[x, 2]]]]]",
        ),
    ] {
        assert_eq!(
            ctx.eval(&read(input)).unwrap(),
            ctx.eval(&read(expected)).unwrap(),
            "{input}"
        );
    }
}

/// What is not solved here stays (CAS-283).
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // Nonlinear.
        "DSolve[List[Equation[D[y[x], x], Power[z[x], 2]], Equation[D[z[x], x], y[x]]], \
         List[y[x], z[x]]]",
        // A variable coefficient.
        "DSolve[List[Equation[D[y[x], x], Times[x, z[x]]], Equation[D[z[x], x], y[x]]], \
         List[y[x], z[x]]]",
        // The highest derivatives' matrix singular.
        "DSolve[List[Equation[Plus[D[y[x], x], D[z[x], x]], y[x]], \
         Equation[Plus[D[y[x], x], D[z[x], x]], z[x]]], List[y[x], z[x]]]",
        // A forcing term without a transform.
        "DSolve[List[Equation[D[y[x], x], Plus[z[x], Tan[x]]], \
         Equation[D[z[x], x], Times[-1, y[x]]]], List[y[x], z[x]]]",
        // Fewer equations than unknowns.
        "DSolve[List[Equation[D[y[x], x], z[x]]], List[y[x], z[x]]]",
    ] {
        let answer = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            answer.head(),
            Some(&Expr::from(Symbol::DSOLVE)),
            "{input}: {answer:?}"
        );
    }
}
