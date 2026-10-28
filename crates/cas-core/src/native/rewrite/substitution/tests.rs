//! Tests of the trigonometric substitutions (CAS-397): giac's forms, values
//! at a point, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn forms() {
    check(&[
        (
            "HalfTan[Sin[x]]",
            "Times[2, Tan[Times[1/2, x]], Power[Plus[1, Power[Tan[Times[1/2, x]], 2]], -1]]",
        ),
        (
            "HalfTan[Cos[x]]",
            "Times[Plus[1, Times[-1, Power[Tan[Times[1/2, x]], 2]]], \
             Power[Plus[1, Power[Tan[Times[1/2, x]], 2]], -1]]",
        ),
        (
            "HalfTan[Tan[x]]",
            "Times[2, Tan[Times[1/2, x]], Power[Plus[1, Times[-1, Power[Tan[Times[1/2, x]], 2]]], -1]]",
        ),
        // `tan(2x/2)` is `tan(x)` once evaluated.
        (
            "HalfTan[Sin[Times[2, x]]]",
            "Times[2, Tan[x], Power[Plus[1, Power[Tan[x], 2]], -1]]",
        ),
        ("TanToSinCos[Tan[x]]", "Times[Sin[x], Power[Cos[x], -1]]"),
        ("SinToCosTan[Sin[x]]", "Times[Cos[x], Tan[x]]"),
        ("CosToSinTan[Cos[x]]", "Times[Sin[x], Power[Tan[x], -1]]"),
        (
            "TanToSinCos2[Tan[x]]",
            "Times[Sin[Times[2, x]], Power[Plus[1, Cos[Times[2, x]]], -1]]",
        ),
        (
            "TanToCosSin2[Tan[x]]",
            "Times[Plus[1, Times[-1, Cos[Times[2, x]]]], Power[Sin[Times[2, x]], -1]]",
        ),
        // What stays: the other functions, and the hyperbolic ones but
        // in `HalfTanHyperbolicToExp`.
        ("TanToSinCos[Sin[x]]", "Sin[x]"),
        ("SinToCosTan[Plus[Cos[x], Tan[x]]]", "Plus[Cos[x], Tan[x]]"),
        ("HalfTan[Sinh[x]]", "Sinh[x]"),
        ("HalfTan[List[1, x]]", "List[1, x]"),
    ]);
    for f in [
        "HalfTan",
        "HalfTanHyperbolicToExp",
        "TanToSinCos",
        "SinToCosTan",
        "CosToSinTan",
        "TanToSinCos2",
        "TanToCosSin2",
    ] {
        check_values(
            f,
            &[
                "Sin[x]",
                "Cos[Times[2, y]]",
                "Tan[Plus[x, y]]",
                "Times[Sin[x], Cos[y], Tan[x]]",
                "Plus[Power[Sin[y], 2], Cot[x], Sec[y]]",
                "Plus[Sinh[x], Cos[Sinh[y]]]",
            ],
        );
    }
}
