//! Tests of polynomial systems by lexicographic Gröbner bases (CAS-289).

use crate::context::Context;
use crate::expr::Expr;
use crate::symbol::Symbol;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, answer)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `count` tuples solving `equations` in `unknowns`, each meeting them to
/// within 10⁻²⁰ part by part.
fn solves(head: &str, equations: &str, unknowns: &str, count: usize) {
    let ctx = Context::new();
    let input = format!("{head}[List[{equations}], List[{unknowns}]]");
    let answer = ctx.eval(&read(&input)).unwrap();
    assert_eq!(
        answer.head(),
        Some(&Expr::from(Symbol::LIST)),
        "{input}: {answer:?}"
    );
    assert_eq!(answer.args().len(), count, "{input}: {answer:?}");
    let names: Vec<&str> = unknowns.split(", ").collect();
    for tuple in answer.args() {
        let rules: Vec<String> = names
            .iter()
            .zip(tuple.args())
            .map(|(x, v)| format!("Rule[{x}, {v:?}]"))
            .collect();
        for eq in split_equations(equations) {
            let residual = format!(
                "N[ReplaceAll[Plus[{}, Times[-1, {}]], List[{}]], 30]",
                eq.0,
                eq.1,
                rules.join(", ")
            );
            for part in ["Re", "Im"] {
                let small = format!("Less[Abs[{part}[{residual}]], Power[10, -20]]");
                assert_eq!(
                    ctx.eval(&read(&small)).unwrap(),
                    read("True"),
                    "{input}: {tuple:?}, {}",
                    eq.0
                );
            }
        }
    }
}

/// The sides of the top-level `Equation[l, r]`s of a comma-separated list.
fn split_equations(equations: &str) -> Vec<(String, String)> {
    let e = read(&format!("Hold[{equations}]"));
    e.args()
        .iter()
        .map(|eq| {
            let [l, r] = eq.args() else {
                panic!("{eq:?}");
            };
            (format!("{l:?}"), format!("{r:?}"))
        })
        .collect()
}

/// Systems with rational solutions, as giac's, ascending.
#[test]
fn exact() {
    check(&[
        (
            "Solve[List[Equation[Times[x, y], 1], Equation[Times[y, z], 1], \
             Equation[Times[x, z], 1]], List[x, y, z]]",
            "List[List[-1, -1, -1], List[1, 1, 1]]",
        ),
        (
            "Solve[List[Equation[Plus[x, y, z], 1], \
             Equation[Plus[Power[x, 2], Power[y, 2], Power[z, 2]], 1], \
             Equation[Plus[Power[x, 3], Power[y, 3], Power[z, 3]], 1]], List[x, y, z]]",
            "List[List[0, 0, 1], List[0, 1, 0], List[1, 0, 0]]",
        ),
        // More equations than unknowns.
        (
            "Solve[List[Plus[Power[x, 2], -1], Plus[Power[y, 2], -1], Plus[x, Times[-1, y]]], \
             List[x, y]]",
            "List[List[-1, -1], List[1, 1]]",
        ),
        (
            "Solve[List[Equation[Plus[Power[x, 2], Power[y, 2]], 2], Equation[Times[x, y], 1], \
             Equation[Plus[x, z], 0]], List[x, y, z]]",
            "List[List[-1, -1, 1], List[1, 1, -1]]",
        ),
        // None: inconsistent, or no real solution.
        (
            "Solve[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2], \
             Equation[Times[x, y], 3]], List[x, y]]",
            "List[]",
        ),
        (
            "Solve[List[Equation[Plus[x, y, z], 0], \
             Equation[Plus[Times[x, y], Times[y, z], Times[x, z]], 0], \
             Equation[Times[x, y, z], 1]], List[x, y, z]]",
            "List[]",
        ),
    ]);
}

/// Algebraic solutions, real and complex.
#[test]
fn algebraic() {
    solves(
        "Solve",
        "Equation[Plus[Power[x, 2], Power[y, 2], Power[z, 2]], 1], Equation[x, y], \
         Equation[y, z]",
        "x, y, z",
        2,
    );
    solves(
        "Solve",
        "Equation[Plus[Power[x, 2], y, z], 1], Equation[Plus[x, Power[y, 2], z], 1], \
         Equation[Plus[x, y, Power[z, 2]], 1]",
        "x, y, z",
        5,
    );
    // The roots of z³ = 1 and their permutations.
    solves(
        "CSolve",
        "Equation[Plus[x, y, z], 0], \
         Equation[Plus[Times[x, y], Times[y, z], Times[x, z]], 0], \
         Equation[Times[x, y, z], 1]",
        "x, y, z",
        6,
    );
}

/// Infinitely many solutions, or parameters: the call stays.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "Solve[List[Equation[Times[x, y], 0], Equation[Times[x, z], 0]], List[x, y, z]]",
        "Solve[List[Equation[Plus[x, y, z], a], Equation[Times[x, y, z], 1], \
         Equation[Plus[x, y], 1]], List[x, y, z]]",
    ] {
        let answer = ctx.eval(&read(input)).unwrap();
        assert_eq!(
            answer.head(),
            Some(&Expr::from(Symbol::SOLVE)),
            "{input}: {answer:?}"
        );
    }
}
