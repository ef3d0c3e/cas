//! Tests of the Hilbert transform form of `HilbertMatrix` (CAS-648):
//! giac's answers (crates/cas-core/proto/wavelet/giac/p5.out) by values at
//! two points, `undef`, the inert `Hilbert`, the forms that keep their
//! meaning.

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

/// `(input, expected)`: the two within `10⁻²⁰` at `v = 7/10` and `13/10`.
fn check_values(v: &str, cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = read(expected);
        for k in [7, 13] {
            let at = |f: &Expr| {
                f.substitute(&[(Expr::from(Symbol::new(v)), Expr::from(k) / Expr::from(10))])
            };
            let gap = Expr::call(
                Symbol::N,
                [
                    Expr::call(Symbol::ABS, [at(&found) - at(&expected)]),
                    Expr::from(30),
                ],
            );
            let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-20)]);
            if ctx.eval(&small).unwrap().as_symbol() != Some(Symbol::TRUE) {
                wrong.push(format!(
                    "{input} at {k}/10: {found:?}, expected {expected:?}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn oscillations() {
    check_values(
        "t",
        &[
            ("HilbertMatrix[Cos[t], t]", "Sin[t]"),
            ("HilbertMatrix[Sin[t], t]", "Times[-1, Cos[t]]"),
            ("HilbertMatrix[Cos[Times[2, t]], t]", "Sin[Times[2, t]]"),
            (
                "HilbertMatrix[Sin[Times[3, t]], t]",
                "Times[-1, Cos[Times[3, t]]]",
            ),
            ("HilbertMatrix[Times[2, Cos[t]], t]", "Times[2, Sin[t]]"),
            (
                "HilbertMatrix[Plus[Cos[t], Sin[Times[2, t]]], t]",
                "Plus[Sin[t], Times[-1, Cos[Times[2, t]]]]",
            ),
        ],
    );
    check_values("x", &[("HilbertMatrix[Sin[x]]", "Times[-1, Cos[x]]")]);
}

/// The residue rule (CAS-654); the forms from a sympy prototype.
#[test]
fn rational_functions() {
    check(&[(
        "HilbertMatrix[Power[Plus[1, Power[t, 2]], -1], t]",
        "Times[t, Power[Plus[1, Power[t, 2]], -1]]",
    )]);
    check_values(
        "t",
        &[
            (
                "HilbertMatrix[Power[Plus[1, Power[t, 4]], -1], t]",
                "Times[Power[2, 1/2], t, Plus[1, Power[t, 2]], Power[Times[2, Plus[1, Power[t, 4]]], -1]]",
            ),
            (
                "HilbertMatrix[Power[Plus[1, Power[t, 2]], -2], t]",
                "Times[t, Plus[3, Power[t, 2]], Power[Times[2, Power[Plus[1, Power[t, 2]], 2]], -1]]",
            ),
        ],
    );
    // Real poles, a polynomial part, a symbolic coefficient: not the rule's.
    check(&[
        (
            "HilbertMatrix[Power[Plus[1, Times[-1, Power[t, 2]]], -1], t]",
            "Hilbert[Power[Plus[1, Times[-1, Power[t, 2]]], -1], t]",
        ),
        (
            "HilbertMatrix[Times[Power[t, 2], Power[Plus[1, Power[t, 2]], -1]], t]",
            "Hilbert[Times[Power[t, 2], Power[Plus[1, Power[t, 2]], -1]], t]",
        ),
        (
            "HilbertMatrix[Power[Plus[Power[a, 2], Power[t, 2]], -1], t]",
            "Hilbert[Power[Plus[Power[a, 2], Power[t, 2]], -1], t]",
        ),
    ]);
}

#[test]
fn constants_polynomials_and_inert() {
    check(&[
        ("HilbertMatrix[1, t]", "0"),
        ("HilbertMatrix[a, t]", "0"),
        ("HilbertMatrix[Power[t, 2], t]", "Indeterminate"),
        ("HilbertMatrix[t, t]", "Indeterminate"),
        ("HilbertMatrix[x]", "Indeterminate"),
        (
            "HilbertMatrix[Exp[Times[-1, Power[t, 2]]], t]",
            "Hilbert[Exp[Times[-1, Power[t, 2]]], t]",
        ),
        ("HilbertMatrix[f[t], t]", "Hilbert[f[t], t]"),
        // The matrix and the list form keep their meaning.
        ("HilbertMatrix[1]", "List[List[1]]"),
        ("HilbertMatrix[2]", "List[List[1, 1/2], List[1/2, 1/3]]"),
    ]);
    // No variable, or a list with one: the call stays.
    let ctx = Context::new();
    for input in ["HilbertMatrix[Sin[t], 2]", "HilbertMatrix[List[1, 2], t]"] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
