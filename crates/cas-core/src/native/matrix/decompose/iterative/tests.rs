//! Tests of the iterative solvers (CAS-499): balls holding the solution,
//! narrow, and the calls that stay.

use cas_num::Rational;
use cas_num::ball::Mag;

use crate::expr::View;
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The call's balls hold `solution` and are narrower than `2⁻³⁰`.
fn holds(input: &str, solution: &[(i64, i64)]) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    assert_eq!(found.args().len(), solution.len(), "{input}: {found:?}");
    let narrow = Mag::power_of_two(&cas_num::Integer::from(-30));
    for (x, &(p, q)) in found.args().iter().zip(solution) {
        let View::Ball(b) = x.view() else {
            panic!("{input}: {found:?}");
        };
        assert!(
            b.contains_rational(&Rational::new(p.into(), q.into())),
            "{input}: {found:?}"
        );
        assert!(*b.rad() < narrow, "{input}: {found:?}");
    }
}

#[test]
fn solutions() {
    let a = "List[List[4, 1], List[1, 3]]";
    holds(
        &format!("JacobiLinearSolve[{a}, List[1, 2], Power[10, -12]]"),
        &[(1, 11), (7, 11)],
    );
    holds(
        &format!("GaussSeidelLinearSolve[{a}, List[1, 2], Power[10, -12]]"),
        &[(1, 11), (7, 11)],
    );
    // Over-relaxation, giac's leading ω, and a bound on the steps.
    holds(
        &format!("GaussSeidelLinearSolve[6/5, {a}, List[1, 2], Power[10, -12], 500]"),
        &[(1, 11), (7, 11)],
    );
    // Three unknowns and decimals: x = (1, −1, 2).
    holds(
        "JacobiLinearSolve[N[List[List[10, -1, 2], List[-1, 11, -1], List[2, -1, 10]], 30], \
         List[16, -14, 23], Power[10, -14]]",
        &[(1, 1), (-1, 1), (2, 1)],
    );
}

#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        // Not diagonally dominant: no bound.
        "JacobiLinearSolve[List[List[1, 2], List[3, 4]], List[1, 2], Power[10, -12]]",
        // Two steps do not settle.
        "GaussSeidelLinearSolve[List[List[4, 1], List[1, 3]], List[1, 2], Power[10, -12], 2]",
        // ω outside (0, 2), ε not positive, symbols.
        "GaussSeidelLinearSolve[2, List[List[4, 1], List[1, 3]], List[1, 2], Power[10, -12]]",
        "JacobiLinearSolve[List[List[4, 1], List[1, 3]], List[1, 2], 0]",
        "JacobiLinearSolve[List[List[4, a], List[1, 3]], List[1, 2], Power[10, -12]]",
        "JacobiLinearSolve[List[List[4, 1], List[1, 3]], List[1, 2, 3], Power[10, -12]]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
