//! Tests of `TrigCollect` (CAS-399): giac's answers, values at a point,
//! and what stays.

use super::super::testing::{check, check_values};

#[test]
fn forms() {
    check(&[
        (
            "TrigCollect[Plus[Sin[x], Cos[x]]]",
            "Times[Power[2, 1/2], Cos[Plus[x, Times[-1/4, Pi]]]]",
        ),
        (
            "TrigCollect[Plus[Sin[x], Times[-1, Cos[x]]]]",
            "Times[Power[2, 1/2], Cos[Plus[x, Times[-3/4, Pi]]]]",
        ),
        (
            "TrigCollect[Plus[Times[-1, Sin[x]], Times[-1, Cos[x]]]]",
            "Times[Power[2, 1/2], Cos[Plus[x, Times[3/4, Pi]]]]",
        ),
        (
            "TrigCollect[Plus[Times[3, Sin[x]], Times[4, Cos[x]]]]",
            "Times[5, Cos[Plus[x, Times[-1, ArcTan[3/4]]]]]",
        ),
        (
            "TrigCollect[Times[Sin[x], Cos[x]]]",
            "Times[1/2, Sin[Times[2, x]]]",
        ),
        (
            "TrigCollect[Plus[Power[Cos[x], 2], Times[-1, Power[Sin[x], 2]]]]",
            "Cos[Times[2, x]]",
        ),
        (
            "TrigCollect[Plus[Sin[Times[2, x]], Cos[Times[2, x]], Sin[x]]]",
            "Plus[Sin[x], Times[Power[2, 1/2], Cos[Plus[Times[2, x], Times[-1/4, Pi]]]]]",
        ),
        (
            "TrigCollect[Plus[1, Sin[x], Cos[x]]]",
            "Plus[1, Times[Power[2, 1/2], Cos[Plus[x, Times[-1/4, Pi]]]]]",
        ),
        // What stays: a coefficient not known real, or with the angle's
        // symbol; angles apart; one function alone.
        (
            "TrigCollect[Plus[Times[a, Sin[x]], Cos[x]]]",
            "Plus[Times[a, Sin[x]], Cos[x]]",
        ),
        (
            "TrigCollect[Plus[Sin[Plus[1, x]], Cos[x]]]",
            "Plus[Sin[Plus[1, x]], Cos[x]]",
        ),
        ("TrigCollect[Sin[x]]", "Sin[x]"),
    ]);
    check_values(
        "TrigCollect",
        &[
            "Plus[Sin[x], Cos[x]]",
            "Plus[Times[3, Sin[y]], Times[-4, Cos[y]]]",
            "Plus[Sin[x], Times[Power[3, 1/2], Cos[x]]]",
            "Times[Sin[x], Power[Cos[x], 2]]",
            "Plus[Times[2, Sin[Times[2, x]]], Times[-5, Cos[Times[2, x]]], Cos[y]]",
            "Plus[Times[Pi, Sin[Plus[x, y]]], Times[E, Cos[Plus[x, y]]]]",
        ],
    );
}
