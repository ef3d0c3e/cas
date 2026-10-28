//! Tests of the Laplace transform (CAS-274) and its inverse (CAS-275).

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form: the difference `Together` makes 0.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let gap = format!("Together[Plus[{found:?}, Times[-1, {expected}]]]");
        if ctx.eval(&read(&gap)).unwrap() != Expr::from(0) {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn table_and_rules() {
    check(
        &Context::new(),
        &[
            ("LaplaceTransform[1, t, s]", "Power[s, -1]"),
            ("LaplaceTransform[t, t, s]", "Power[s, -2]"),
            (
                "LaplaceTransform[Power[t, 3], t, s]",
                "Times[6, Power[s, -4]]",
            ),
            (
                "LaplaceTransform[Power[E, Times[-2, t]], t, s]",
                "Power[Plus[s, 2], -1]",
            ),
            (
                "LaplaceTransform[Sin[Times[3, t]], t, s]",
                "Times[3, Power[Plus[Power[s, 2], 9], -1]]",
            ),
            (
                "LaplaceTransform[Cos[Times[3, t]], t, s]",
                "Times[s, Power[Plus[Power[s, 2], 9], -1]]",
            ),
            (
                "LaplaceTransform[Times[t, Power[E, t]], t, s]",
                "Power[Plus[s, -1], -2]",
            ),
            (
                "LaplaceTransform[Times[Power[E, t], Sin[t]], t, s]",
                "Power[Plus[Power[Plus[s, -1], 2], 1], -1]",
            ),
            (
                "LaplaceTransform[Sinh[Times[2, t]], t, s]",
                "Times[2, Power[Plus[Power[s, 2], -4], -1]]",
            ),
            (
                "LaplaceTransform[Times[t, Sin[t]], t, s]",
                "Times[2, s, Power[Plus[Power[s, 2], 1], -2]]",
            ),
            (
                "LaplaceTransform[Power[Sin[t], 2], t, s]",
                "Times[2, Power[Times[s, Plus[Power[s, 2], 4]], -1]]",
            ),
            (
                "LaplaceTransform[Power[Plus[t, 1], 2], t, s]",
                "Plus[Times[2, Power[s, -3]], Times[2, Power[s, -2]], Power[s, -1]]",
            ),
            (
                "LaplaceTransform[Sin[Plus[t, 1]], t, s]",
                "Times[Plus[Cos[1], Times[s, Sin[1]]], Power[Plus[Power[s, 2], 1], -1]]",
            ),
            // giac's is 0.
            (
                "LaplaceTransform[Power[t, 1/2], t, s]",
                "Times[1/2, Power[Pi, 1/2], Power[s, -3/2]]",
            ),
            (
                "LaplaceTransform[Power[t, n], t, s]",
                "Times[Gamma[Plus[n, 1]], Power[s, Times[-1, Plus[n, 1]]]]",
            ),
            (
                "LaplaceTransform[HeavisideTheta[Plus[t, -1]], t, s]",
                "Times[Power[E, Times[-1, s]], Power[s, -1]]",
            ),
            (
                "LaplaceTransform[DiracDelta[Plus[t, -2]], t, s]",
                "Power[E, Times[-2, s]]",
            ),
            (
                "LaplaceTransform[Times[t, HeavisideTheta[Plus[t, -1]]], t, s]",
                "Times[Power[E, Times[-1, s]], Plus[1, s], Power[s, -2]]",
            ),
            (
                "LaplaceTransform[Times[Sin[Plus[t, -1]], HeavisideTheta[Plus[t, -1]]], t, s]",
                "Times[Power[E, Times[-1, s]], Power[Plus[Power[s, 2], 1], -1]]",
            ),
            // laplace(f) in x, into x.
            ("LaplaceTransform[1]", "Power[x, -1]"),
            (
                "LaplaceTransform[Power[E, Times[-1, x]]]",
                "Power[Plus[x, 1], -1]",
            ),
        ],
    );
}

/// Unknown functions stay, and their derivatives transform (CAS-274).
#[test]
fn derivatives() {
    let ctx = Context::new();
    check(
        &ctx,
        &[
            (
                "LaplaceTransform[D[f[t], t], t, s]",
                "Plus[Times[s, LaplaceTransform[f[t], t, s]], Times[-1, f[0]]]",
            ),
            (
                "LaplaceTransform[D[f[t], t, t], t, s]",
                "Plus[Times[Power[s, 2], LaplaceTransform[f[t], t, s]], Times[-1, s, f[0]], \
                 Times[-1, Derivative[1][f][0]]]",
            ),
            (
                "LaplaceTransform[Plus[D[f[t], t], f[t]], t, s]",
                "Plus[Times[Plus[s, 1], LaplaceTransform[f[t], t, s]], Times[-1, f[0]]]",
            ),
        ],
    );
    let stays = read("LaplaceTransform[f[t], t, s]");
    assert_eq!(ctx.eval(&stays).unwrap(), stays);
    let stays = read("LaplaceTransform[Power[E, Power[t, 2]], t, s]");
    assert_eq!(ctx.eval(&stays).unwrap(), stays);
}
