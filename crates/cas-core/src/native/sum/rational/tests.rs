//! Tests of rational sums through `PolyGamma` (CAS-359): giac's finite
//! sums, the convergent infinite ones in closed form, and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same once put together.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        let difference = ctx
            .eval(&Expr::call(
                Symbol::TOGETHER,
                [Expr::call(
                    Symbol::EXPAND,
                    [found.clone() - expected.clone()],
                )],
            ))
            .unwrap();
        if found != expected && difference != Expr::from(0) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `sum(1/k,k,1,n)`, `sum(1/(k+1),k,0,n)`, `sum(1/(2*k+1),k,0,n)`,
/// `sum(1/k^2,k,1,n)` and `sum(1/(k+1/3),k,0,n)`.
#[test]
fn finite_sums() {
    check(&[
        (
            "Sum[Power[k, -1], k, 1, n]",
            "Plus[EulerGamma, PolyGamma[Plus[1, n]]]",
        ),
        (
            "Sum[Power[Plus[1, k], -1], k, 0, n]",
            "Plus[EulerGamma, PolyGamma[Plus[2, n]]]",
        ),
        (
            "Sum[Power[Plus[1, Times[2, k]], -1], k, 0, n]",
            "Plus[Times[1/2, PolyGamma[Plus[3/2, n]]], Times[1/2, EulerGamma], Log[2]]",
        ),
        (
            "Sum[Power[k, -2], k, 1, n]",
            "Plus[Times[1/6, Power[Pi, 2]], Times[-1, PolyGamma[Plus[1, n], 1]]]",
        ),
        (
            "Sum[Power[Plus[1/3, k], -1], k, 0, n]",
            "Plus[PolyGamma[Plus[4/3, n]], Times[-1, PolyGamma[1/3]]]",
        ),
        // A polynomial part beside the terms.
        (
            "Sum[Plus[k, Power[k, -1]], k, 1, n]",
            "Plus[Times[1/2, n, Plus[1, n]], EulerGamma, PolyGamma[Plus[1, n]]]",
        ),
    ]);
}

/// giac's `sum(1/(k*(2*k+1)),k,1,inf)`, `sum(1/(k*(k+1/2)),k,1,inf)`,
/// `sum(1/(2*k-1)^2,k,1,inf)` and `sum(1/(3*k+1)^2,k,0,inf)`.
#[test]
fn infinite_sums() {
    check(&[
        (
            "Sum[Power[Times[k, Plus[1, Times[2, k]]], -1], k, 1, Infinity]",
            "Plus[2, Times[-2, Log[2]]]",
        ),
        (
            "Sum[Power[Times[k, Plus[1/2, k]], -1], k, 1, Infinity]",
            "Plus[4, Times[-4, Log[2]]]",
        ),
        (
            "Sum[Power[Plus[-1, Times[2, k]], -2], k, 1, Infinity]",
            "Times[1/8, Power[Pi, 2]]",
        ),
        (
            "Sum[Power[Plus[1, Times[3, k]], -2], k, 0, Infinity]",
            "Times[1/9, PolyGamma[1/3, 1]]",
        ),
    ]);
}

/// An irreducible quadratic, a pole inside the bounds, and a divergent
/// rational sum stay.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "Sum[Power[Plus[1, Power[k, 2]], -1], k, 0, n]",
        "Sum[Power[Plus[-3, k], -1], k, 1, n]",
        "Sum[Power[Plus[1, Times[2, k]], -1], k, 0, Infinity]",
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert!(
            found.contains(&Expr::from(Symbol::SUM)) || found == Expr::from(Symbol::INFINITY),
            "{input}: {found:?}"
        );
    }
}
