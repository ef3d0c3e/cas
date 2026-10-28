//! Tests of the inverse trigonometric conversions (CAS-400): giac's forms,
//! values at points where each identity holds, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn forms() {
    check(&[
        (
            "ArcCosToArcSin[ArcCos[x]]",
            "Plus[Times[1/2, Pi], Times[-1, ArcSin[x]]]",
        ),
        // giac's example: the sum is π/2.
        (
            "ArcCosToArcSin[Plus[ArcCos[x], ArcSin[x]]]",
            "Times[1/2, Pi]",
        ),
        (
            "ArcSinToArcCos[ArcSin[x]]",
            "Plus[Times[1/2, Pi], Times[-1, ArcCos[x]]]",
        ),
        (
            "ArcSinToArcTan[ArcSin[x]]",
            "ArcTan[Times[x, Power[Plus[1, Times[-1, Power[x, 2]]], -1/2]]]",
        ),
        (
            "ArcCosToArcTan[ArcCos[x]]",
            "Plus[Times[1/2, Pi], \
             Times[-1, ArcTan[Times[x, Power[Plus[1, Times[-1, Power[x, 2]]], -1/2]]]]]",
        ),
        (
            "ArcTanToArcSin[ArcTan[x]]",
            "ArcSin[Times[x, Power[Plus[1, Power[x, 2]], -1/2]]]",
        ),
        (
            "ArcTanToArcCos[ArcTan[x]]",
            "Plus[Times[1/2, Pi], Times[-1, ArcCos[Times[x, Power[Plus[1, Power[x, 2]], -1/2]]]]]",
        ),
        (
            "ArcTrigToLog[ArcSin[x]]",
            "Times[Complex[0, -1], Log[Plus[Times[Complex[0, 1], x], \
             Power[Plus[1, Times[-1, Power[x, 2]]], 1/2]]]]",
        ),
        // What stays: the other functions.
        ("ArcSinToArcCos[ArcCos[x]]", "ArcCos[x]"),
        ("ArcTanToArcSin[ArcSinh[x]]", "ArcSinh[x]"),
    ]);
    // Within (−1, 1): x = 2/9 and y/3 = 7/9.
    for f in ["ArcCosToArcTan", "ArcSinToArcTan"] {
        check_values(
            f,
            &[
                "ArcSin[x]",
                "ArcCos[Times[-1, x]]",
                "Plus[ArcSin[Times[1/3, y]], Times[2, ArcCos[Times[1/3, y]]]]",
                "ArcCos[Times[-1/3, y]]",
            ],
        );
    }
    // Everywhere: past ±1 too.
    for f in ["ArcCosToArcSin", "ArcSinToArcCos"] {
        check_values(
            f,
            &[
                "ArcSin[x]",
                "ArcCos[Times[-1, x]]",
                "Plus[ArcSin[y], Times[2, ArcCos[y]]]",
                "ArcCos[Times[-1, y]]",
            ],
        );
    }
    for f in ["ArcTanToArcCos", "ArcTanToArcSin"] {
        check_values(
            f,
            &["ArcTan[x]", "ArcTan[Times[-1, y]]", "Times[y, ArcTan[y]]"],
        );
    }
    // Off the branch cuts: asin and acos on [−1, 1], atan on the real line.
    check_values(
        "ArcTrigToLog",
        &[
            "ArcSin[x]",
            "ArcSin[Times[-1/3, y]]",
            "ArcCos[x]",
            "ArcCos[Times[-1/3, y]]",
            "ArcTan[x]",
            "ArcTan[y]",
            "ArcTan[Times[-1, y]]",
            "Plus[ArcTan[x], ArcSin[x], ArcCos[x]]",
        ],
    );
}
