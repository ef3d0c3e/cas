//! Tests of `Sinc`, `UnitTriangle` and `Boxcar` (CAS-648): giac's values
//! (crates/cas-core/proto/wavelet/giac/p5.out), the ties, the arity guard.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
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

#[test]
fn sinc() {
    check(&[
        ("Sinc[0]", "1"),
        ("Sinc[Pi]", "0"),
        ("Sinc[Times[2, Pi]]", "0"),
        ("Sinc[Times[1/2, Pi]]", "Times[2, Power[Pi, -1]]"),
        ("Sinc[x]", "Times[Sin[x], Power[x, -1]]"),
        ("Sinc[List[0, Pi]]", "List[1, 0]"),
    ]);
    let ctx = Context::new();
    let zero = Expr::call(Symbol::N, [Expr::from(0)]);
    let at_zero = Expr::call(Symbol::SINC, [ctx.eval(&zero).unwrap()]);
    assert_eq!(ctx.eval(&at_zero).unwrap(), Expr::from(1));
}

#[test]
fn triangle() {
    check(&[
        ("UnitTriangle[0]", "1"),
        ("UnitTriangle[1/2]", "1/2"),
        ("UnitTriangle[-1/2]", "1/2"),
        ("UnitTriangle[-1/3]", "2/3"),
        ("UnitTriangle[1]", "0"),
        ("UnitTriangle[-1]", "0"),
        ("UnitTriangle[2]", "0"),
        (
            "UnitTriangle[x]",
            "Plus[Times[Plus[1, Times[-1, x]], Plus[HeavisideTheta[x], \
             Times[-1, HeavisideTheta[Plus[x, -1]]]]], Times[Plus[1, x], \
             Plus[HeavisideTheta[Times[-1, x]], Times[-1, \
             HeavisideTheta[Plus[Times[-1, x], -1]]]]]]",
        ),
    ]);
}

#[test]
fn boxcar() {
    check(&[
        ("Boxcar[1, 2, 3/2]", "1"),
        ("Boxcar[1, 2, 1]", "1/2"),
        ("Boxcar[1, 2, 2]", "1/2"),
        ("Boxcar[1, 2, 0]", "0"),
        ("Boxcar[2, 1, 3/2]", "-1"),
        ("Boxcar[1, 1, 1]", "0"),
        (
            "Boxcar[a, b, x]",
            "Plus[HeavisideTheta[Plus[x, Times[-1, a]]], \
             Times[-1, HeavisideTheta[Plus[x, Times[-1, b]]]]]",
        ),
    ]);
    // Not three arguments: the call stays, with giac's message.
    let ctx = Context::new();
    for input in ["Boxcar[1, 2]", "Boxcar[1, 2, 3, 4]"] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
