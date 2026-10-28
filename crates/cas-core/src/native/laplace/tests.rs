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

/// `(input, expected)` in full form, equal by value at `t = 13/5`,
/// `a = 3/4`, each part apart.
fn check_values(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let gap = format!(
            "N[ReplaceAll[Plus[{found:?}, Times[-1, {expected}]], \
             List[Rule[t, 13/5], Rule[a, 3/4]]]]"
        );
        let close = ["Re", "Im"].iter().all(|part| {
            let test = format!("Less[Abs[{part}[{gap}]], 1/10000000000]");
            ctx.eval(&read(&test)).unwrap() == read("True")
        });
        if !close {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The inverse (CAS-275): giac's table, by value.
#[test]
fn inverses() {
    check_values(&[
        ("InverseLaplaceTransform[Power[s, -1], s, t]", "1"),
        ("InverseLaplaceTransform[Power[s, -2], s, t]", "t"),
        (
            "InverseLaplaceTransform[Power[Plus[s, 2], -1], s, t]",
            "Power[E, Times[-2, t]]",
        ),
        (
            "InverseLaplaceTransform[Times[3, Power[Plus[Power[s, 2], 9], -1]], s, t]",
            "Sin[Times[3, t]]",
        ),
        (
            "InverseLaplaceTransform[Times[s, Power[Plus[Power[s, 2], 9], -1]], s, t]",
            "Cos[Times[3, t]]",
        ),
        (
            "InverseLaplaceTransform[Power[Plus[Power[s, 2], Times[-2, s], 2], -1], s, t]",
            "Times[Power[E, t], Sin[t]]",
        ),
        (
            "InverseLaplaceTransform[Power[Times[Plus[s, 1], Plus[s, 2]], -1], s, t]",
            "Plus[Power[E, Times[-1, t]], Times[-1, Power[E, Times[-2, t]]]]",
        ),
        (
            "InverseLaplaceTransform[Times[Plus[s, 1], Power[Plus[Power[s, 2], Times[2, s], 5], -1]], s, t]",
            "Times[Power[E, Times[-1, t]], Cos[Times[2, t]]]",
        ),
        (
            "InverseLaplaceTransform[Power[Times[Power[s, 2], Plus[s, 1]], -1], s, t]",
            "Plus[t, -1, Power[E, Times[-1, t]]]",
        ),
        (
            "InverseLaplaceTransform[Power[Plus[s, -1], -3], s, t]",
            "Times[1/2, Power[t, 2], Power[E, t]]",
        ),
        (
            "InverseLaplaceTransform[Times[Power[E, Times[-2, s]], Power[s, -1]], s, t]",
            "HeavisideTheta[Plus[t, -2]]",
        ),
        (
            "InverseLaplaceTransform[Times[s, Power[Plus[Power[s, 2], 1], -2]], s, t]",
            "Times[1/2, t, Sin[t]]",
        ),
        (
            "InverseLaplaceTransform[Power[Plus[Power[s, 4], -1], -1], s, t]",
            "Plus[Times[1/4, Power[E, t]], Times[-1/4, Power[E, Times[-1, t]]], Times[-1/2, Sin[t]]]",
        ),
        (
            "InverseLaplaceTransform[Times[Power[E, Times[-1, s]], Power[Plus[Power[s, 2], 1], -1]], s, t]",
            "Times[Sin[Plus[t, -1]], HeavisideTheta[Plus[t, -1]]]",
        ),
        (
            "InverseLaplaceTransform[Power[Plus[s, Times[-1, a]], -1], s, t]",
            "Power[E, Times[a, t]]",
        ),
        (
            "InverseLaplaceTransform[Times[Plus[Times[2, s], 3], Power[Plus[Power[s, 2], Times[4, s], 13], -1]], s, t]",
            "Times[Power[E, Times[-2, t]], Plus[Times[2, Cos[Times[3, t]]], Times[-1/3, Sin[Times[3, t]]]]]",
        ),
        (
            "InverseLaplaceTransform[Power[Plus[Power[s, 3], 1], -1], s, t]",
            "Plus[Times[1/3, Power[E, Times[-1, t]]], \
             Times[-1/3, Cos[Times[1/2, Power[3, 1/2], t]], Power[E, Times[1/2, t]]], \
             Times[Power[3, -1/2], Sin[Times[1/2, Power[3, 1/2], t]], Power[E, Times[1/2, t]]]]",
        ),
        // giac's is wrong.
        (
            "InverseLaplaceTransform[Power[s, -1/2], s, t]",
            "Power[Times[Pi, t], -1/2]",
        ),
        // Round trips.
        (
            "InverseLaplaceTransform[LaplaceTransform[Power[t, 2], t, s], s, t]",
            "Power[t, 2]",
        ),
        (
            "InverseLaplaceTransform[LaplaceTransform[Times[Power[E, Times[-1, t]], Cos[t]], t, s], s, t]",
            "Times[Power[E, Times[-1, t]], Cos[t]]",
        ),
        (
            "InverseLaplaceTransform[LaplaceTransform[Times[t, Power[E, Times[2, t]]], t, s], s, t]",
            "Times[t, Power[E, Times[2, t]]]",
        ),
    ]);
    let ctx = Context::new();
    assert_eq!(
        ctx.eval(&read("InverseLaplaceTransform[1, s, t]")).unwrap(),
        read("DiracDelta[t]")
    );
    assert_eq!(
        ctx.eval(&read("InverseLaplaceTransform[Power[Plus[s, 1], -1]]"))
            .unwrap(),
        ctx.eval(&read("Power[E, Times[-1, x]]")).unwrap()
    );
}
