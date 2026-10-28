//! Tests of `ZTransform` and `KroneckerDelta` (CAS-308): the table, the
//! product rules, unknown sequences and their shifts, what stays.

use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, answer)` in full form, both evaluated, the answer put
/// together: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(&format!("Together[{expected}]"))).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn kronecker_delta() {
    let ctx = Context::new();
    for (input, value) in [
        ("KroneckerDelta[0]", "1"),
        ("KroneckerDelta[3]", "0"),
        ("KroneckerDelta[-1/2]", "0"),
        ("KroneckerDelta[n]", "KroneckerDelta[n]"),
        ("KroneckerDelta[List[0, 1]]", "List[1, 0]"),
    ] {
        assert_eq!(ctx.eval(&read(input)).unwrap(), read(value), "{input}");
    }
}

/// The table, as giac's answers put together (CAS-308).
#[test]
fn table() {
    check(&[
        ("ZTransform[1, n, z]", "Times[z, Power[Plus[z, -1], -1]]"),
        ("ZTransform[c, n, z]", "Times[c, z, Power[Plus[z, -1], -1]]"),
        ("ZTransform[n, n, z]", "Times[z, Power[Plus[z, -1], -2]]"),
        (
            "ZTransform[Power[n, 2], n, z]",
            "Times[z, Plus[z, 1], Power[Plus[z, -1], -3]]",
        ),
        (
            "ZTransform[Power[a, n], n, z]",
            "Times[z, Power[Plus[z, Times[-1, a]], -1]]",
        ),
        (
            "ZTransform[Times[n, Power[a, n]], n, z]",
            "Times[a, z, Power[Plus[z, Times[-1, a]], -2]]",
        ),
        (
            "ZTransform[Plus[Power[2, n], Power[3, n]], n, z]",
            "Plus[Times[z, Power[Plus[z, -2], -1]], Times[z, Power[Plus[z, -3], -1]]]",
        ),
        (
            "ZTransform[Power[2, Plus[n, 1]], n, z]",
            "Times[2, z, Power[Plus[z, -2], -1]]",
        ),
        (
            "ZTransform[Sin[n], n, z]",
            "Times[z, Sin[1], Power[Plus[Power[z, 2], Times[-2, z, Cos[1]], 1], -1]]",
        ),
        (
            "ZTransform[Cos[Times[w, n]], n, z]",
            "Times[z, Plus[z, Times[-1, Cos[w]]], \
             Power[Plus[Power[z, 2], Times[-2, z, Cos[w]], 1], -1]]",
        ),
        (
            "ZTransform[Power[Factorial[n], -1], n, z]",
            "Exp[Power[z, -1]]",
        ),
        (
            "ZTransform[Binomial[n, 2], n, z]",
            "Times[z, Power[Plus[z, -1], -3]]",
        ),
        ("ZTransform[KroneckerDelta[n], n, z]", "1"),
        (
            "ZTransform[KroneckerDelta[Plus[n, -2]], n, z]",
            "Power[z, -2]",
        ),
        ("ZTransform[KroneckerDelta[Plus[n, 2]], n, z]", "0"),
        // This CAS's step is 1/2 at 0.
        (
            "ZTransform[HeavisideTheta[Plus[n, -2]], n, z]",
            "Plus[Times[1/2, Power[z, -2]], Times[Power[z, -2], Power[Plus[z, -1], -1]]]",
        ),
        (
            "ZTransform[Times[Power[n, 3], Power[2, n]], n, z]",
            "Times[2, z, Plus[Power[z, 2], Times[8, z], 4], Power[Plus[z, -2], -4]]",
        ),
    ]);
}

/// Products, unknown sequences and their shifts (CAS-308).
#[test]
fn rules() {
    check(&[
        // Scaling: aⁿ/n! is e^(a/z).
        (
            "ZTransform[Times[Power[2, n], Power[Factorial[n], -1]], n, z]",
            "Exp[Times[2, Power[z, -1]]]",
        ),
        // Sifting by a Kronecker delta.
        (
            "ZTransform[Times[Power[n, 2], KroneckerDelta[Plus[n, -3]]], n, z]",
            "Times[9, Power[z, -3]]",
        ),
        // A product of sines reduced.
        (
            "ZTransform[Times[Sin[n], Cos[n]], n, z]",
            "ZTransform[Times[1/2, Sin[Times[2, n]]], n, z]",
        ),
        // u(n + 2) = z² (U − u(0) − u(1)/z).
        (
            "ZTransform[u[Plus[n, 2]], n, z]",
            "Times[Power[z, 2], Plus[ZTransform[u[n], n, z], Times[-1, u[0]], \
             Times[-1, u[1], Power[z, -1]]]]",
        ),
    ]);
    let ctx = Context::new();
    for input in [
        "ZTransform[u[n], n, z]",
        "ZTransform[Log[Plus[n, 1]], n, z]",
        "ZTransform[u[Plus[n, -1]], n, z]",
        "ZTransform[1, 2, z]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}

/// `(input, answer)`: the same expression, both evaluated, the answer as
/// it stands.
fn same(cases: &[(&str, &str)]) {
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

/// The inverse by residues, Kronecker deltas from the pole at 0, conjugate
/// pairs as sines and cosines (CAS-309).
#[test]
fn inverse() {
    same(&[
        (
            "InverseZTransform[Times[z, Power[Plus[z, -1], -1]], z, n]",
            "1",
        ),
        (
            "InverseZTransform[Times[z, Power[Plus[z, Times[-1, a]], -1]], z, n]",
            "Power[a, n]",
        ),
        (
            "InverseZTransform[Times[z, Power[Plus[z, -1], -2]], z, n]",
            "n",
        ),
        (
            "InverseZTransform[Power[Plus[z, -2], -1], z, n]",
            "Plus[Times[1/2, Power[2, n]], Times[-1/2, KroneckerDelta[n]]]",
        ),
        (
            "InverseZTransform[Times[Power[z, 2], Power[Plus[z, -1], -1], Power[Plus[z, -2], -1]], \
             z, n]",
            "Plus[Times[2, Power[2, n]], -1]",
        ),
        (
            "InverseZTransform[Times[z, Power[Plus[Power[z, 2], 1], -1]], z, n]",
            "Sin[Times[1/2, n, Pi]]",
        ),
        (
            "InverseZTransform[Times[Plus[Power[z, 2], z], Power[Plus[z, -1], -3]], z, n]",
            "Power[n, 2]",
        ),
        ("InverseZTransform[1, z, n]", "KroneckerDelta[n]"),
        (
            "InverseZTransform[Power[z, -1], z, n]",
            "KroneckerDelta[Plus[n, -1]]",
        ),
        (
            "InverseZTransform[Exp[Power[z, -1]], z, n]",
            "Power[Factorial[n], -1]",
        ),
        (
            "InverseZTransform[Exp[Times[2, Power[z, -1]]], z, n]",
            "Times[Power[2, n], Power[Factorial[n], -1]]",
        ),
        // giac refuses this one: the primitive sixth roots of unity.
        (
            "InverseZTransform[Times[z, Power[Plus[Power[z, 2], Times[-1, z], 1], -1]], z, n]",
            "Times[2, Power[3, -1/2], Sin[Times[1/3, n, Pi]]]",
        ),
    ]);
    let ctx = Context::new();
    for input in [
        // Not a one-sided transform: it grows like z.
        "InverseZTransform[Times[Power[z, 2], Power[Plus[z, -1], -1]], z, n]",
        "InverseZTransform[Log[z], z, n]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
