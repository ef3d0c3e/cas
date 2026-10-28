//! Tests of the orthogonal polynomials (CAS-DOC-27 §5): giac's polynomials,
//! values at numbers, exactly and in balls, lists, and what stays.

use crate::expr::View;
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: the same expression.
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

/// giac's `legendre(3,x)`, `hermite(4,x)`, `laguerre(2,x,0)`,
/// `laguerre(2,x,b)`, `tchebyshev1(4,x)`, `tchebyshev2(3,x)`, degree 0.
#[test]
fn polynomials() {
    check(&[
        (
            "LegendreP[3, x]",
            "Plus[Times[5/2, Power[x, 3]], Times[-3/2, x]]",
        ),
        ("LegendreP[0, x]", "1"),
        ("LegendreP[1, y]", "y"),
        (
            "HermiteH[4, x]",
            "Plus[Times[16, Power[x, 4]], Times[-48, Power[x, 2]], 12]",
        ),
        ("HermiteH[0, x]", "1"),
        (
            "LaguerreL[2, x]",
            "Plus[1, Times[-2, x], Times[1/2, Power[x, 2]]]",
        ),
        (
            "LaguerreL[2, x, b]",
            "Plus[Times[1/2, Power[b, 2]], Times[-1, b, x], Times[3/2, b], Times[1/2, Power[x, 2]], Times[-2, x], 1]",
        ),
        (
            "ChebyshevT[4, x]",
            "Plus[Times[8, Power[x, 4]], Times[-8, Power[x, 2]], 1]",
        ),
        (
            "ChebyshevU[3, x]",
            "Plus[Times[8, Power[x, 3]], Times[-4, x]]",
        ),
        ("ChebyshevU[0, x]", "1"),
    ]);
}

/// Values: giac's `legendre(3,1/2)` is `-7/16`; `T_n(cos θ) = cos(nθ)` at
/// `cos θ = 1/2`; a ball keeps its accuracy at a high degree.
#[test]
fn values() {
    check(&[
        ("LegendreP[3, 1/2]", "-7/16"),
        ("ChebyshevT[6, 1/2]", "1"),
        ("ChebyshevT[3, 1/2]", "-1"),
        ("HermiteH[3, 2]", "40"),
    ]);
    let ctx = Context::new();
    let found = ctx.eval(&read("LegendreP[200, N[1/3, 30]]")).unwrap();
    let View::Ball(b) = found.view() else {
        panic!("{found:?}");
    };
    assert!(b.rel_accuracy_bits() >= 60, "{found:?}");
}

/// What stays: a symbolic, negative or fractional degree.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "LegendreP[n, x]",
        "HermiteH[-1, x]",
        "ChebyshevT[1/2, x]",
        "LaguerreL[2]",
    ] {
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
    }
}
