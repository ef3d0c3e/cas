//! Tests of inequalities (CAS-277 …): the answers' forms.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(inequality, answer)` in full form, `Solve[inequality, x]` and the
/// answer both evaluated: the same expression.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (ineq, expected) in cases {
        let found = ctx.eval(&read(&format!("Solve[{ineq}, x]"))).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{ineq}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// One rational inequality by its sign chart (CAS-277).
#[test]
fn sign_charts() {
    check(
        &Context::new(),
        &[
            (
                "Greater[Power[x, 2], 4]",
                "List[Less[x, -2], Greater[x, 2]]",
            ),
            ("Less[Power[x, 2], 4]", "List[And[Less[-2, x], Less[x, 2]]]"),
            (
                "GreaterEqual[Power[x, 2], 4]",
                "List[LessEqual[x, -2], GreaterEqual[x, 2]]",
            ),
            ("Greater[Plus[Times[2, x], 1], 0]", "List[Greater[x, -1/2]]"),
            (
                "LessEqual[Plus[Times[-3, x], 6], 0]",
                "List[GreaterEqual[x, 2]]",
            ),
            (
                "GreaterEqual[Times[Plus[x, -1], Power[Plus[x, 2], -1]], 0]",
                "List[Less[x, -2], GreaterEqual[x, 1]]",
            ),
            (
                "Less[Times[Plus[x, -1], Power[Plus[x, 2], -1]], 0]",
                "List[And[Less[-2, x], Less[x, 1]]]",
            ),
            (
                "Greater[Plus[Power[x, 3], Times[-1, x]], 0]",
                "List[And[Less[-1, x], Less[x, 0]], Greater[x, 1]]",
            ),
            ("Greater[Plus[Power[x, 2], 1], 0]", "List[x]"),
            ("Less[Plus[Power[x, 2], 1], 0]", "List[]"),
            (
                "LessEqual[Plus[Power[x, 2], Times[-2, x], -3], 0]",
                "List[And[LessEqual[-1, x], LessEqual[x, 3]]]",
            ),
            (
                "Greater[Power[x, -1], 1]",
                "List[And[Less[0, x], Less[x, 1]]]",
            ),
            (
                "Less[Power[x, 2], 2]",
                "List[And[Less[Times[-1, Power[2, 1/2]], x], Less[x, Power[2, 1/2]]]]",
            ),
            ("Greater[x, a]", "List[Greater[x, a]]"),
            ("Unequal[x, 1]", "List[Less[x, 1], Greater[x, 1]]"),
            (
                "Greater[Power[Plus[x, -1], 2], 0]",
                "List[Less[x, 1], Greater[x, 1]]",
            ),
            ("GreaterEqual[Power[Plus[x, -1], 2], 0]", "List[x]"),
            ("LessEqual[Power[Plus[x, -1], 2], 0]", "List[1]"),
            (
                "Less[Plus[Power[x, 4], Times[-5, Power[x, 2]], 4], 0]",
                "List[And[Less[-2, x], Less[x, -1]], And[Less[1, x], Less[x, 2]]]",
            ),
            // Free of x.
            ("Greater[2, 1]", "List[x]"),
        ],
    );
    // A parameter's sign from the assumptions: a x > 1 with a < 0.
    let negative = Context::new().assuming(read("Less[a, 0]"));
    check(
        &negative,
        &[("Greater[Times[a, x], 1]", "List[Less[x, Power[a, -1]]]")],
    );
}

/// Systems and connectives (CAS-278).
#[test]
fn systems() {
    check(
        &Context::new(),
        &[
            (
                "List[Greater[x, 1], Less[x, 3]]",
                "List[And[Less[1, x], Less[x, 3]]]",
            ),
            (
                "And[Greater[Plus[Power[x, 2], Times[-1, x]], 0], Less[x, 5]]",
                "List[Less[x, 0], And[Less[1, x], Less[x, 5]]]",
            ),
            (
                "Or[Less[x, -1], Greater[x, 2]]",
                "List[Less[x, -1], Greater[x, 2]]",
            ),
            (
                "And[Less[-1, x], LessEqual[x, 3]]",
                "List[And[Less[-1, x], LessEqual[x, 3]]]",
            ),
            (
                "List[Less[Power[x, 2], 9], Unequal[x, 0]]",
                "List[And[Less[-3, x], Less[x, 0]], And[Less[0, x], Less[x, 3]]]",
            ),
            ("Not[Greater[x, 1]]", "List[LessEqual[x, 1]]"),
            ("List[Greater[x, 3], Less[x, 1]]", "List[]"),
            (
                "Or[Greater[Power[x, 2], 1], Less[x, 0]]",
                "List[Less[x, 0], Greater[x, 1]]",
            ),
            ("List[GreaterEqual[x, 0], Equal[Power[x, 2], 4]]", "List[2]"),
            // Intervals meeting at a point one holds join.
            ("Or[LessEqual[x, 1], GreaterEqual[x, 1]]", "List[x]"),
            (
                "Or[Less[x, 1], Greater[x, 1]]",
                "List[Less[x, 1], Greater[x, 1]]",
            ),
            ("Or[Less[x, 1], Equal[x, 1]]", "List[LessEqual[x, 1]]"),
            ("Not[Unequal[x, 2]]", "List[2]"),
            (
                "Not[And[Less[0, x], LessEqual[x, 1]]]",
                "List[LessEqual[x, 0], Greater[x, 1]]",
            ),
            ("List[Greater[x, 0], True]", "List[Greater[x, 0]]"),
            ("Or[Greater[x, 0], False]", "List[Greater[x, 0]]"),
            ("And[Greater[x, 1], Less[x, 1]]", "List[]"),
            ("And[GreaterEqual[x, 1], LessEqual[x, 1]]", "List[1]"),
        ],
    );
}
