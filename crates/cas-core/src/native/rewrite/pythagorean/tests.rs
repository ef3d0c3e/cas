//! Tests of the Pythagorean forms (CAS-398): giac's answers, values at a
//! point, lists and relations, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn forms() {
    check(&[
        (
            "TrigToCos[Plus[Power[Sin[x], 4], Power[Sin[x], 2]]]",
            "Plus[2, Times[-3, Power[Cos[x], 2]], Power[Cos[x], 4]]",
        ),
        (
            "TrigToCos[Power[Sin[x], 3]]",
            "Plus[Sin[x], Times[-1, Power[Cos[x], 2], Sin[x]]]",
        ),
        ("TrigToCos[Tan[x]]", "Times[Sin[x], Power[Cos[x], -1]]"),
        (
            "TrigToCos[Plus[Power[Sin[Times[2, x]], 2], Power[Sin[x], 2]]]",
            "Plus[2, Times[-1, Power[Cos[Times[2, x]], 2]], Times[-1, Power[Cos[x], 2]]]",
        ),
        (
            "TrigToSin[Plus[Power[Cos[x], 4], Power[Cos[x], 2]]]",
            "Plus[2, Times[-3, Power[Sin[x], 2]], Power[Sin[x], 4]]",
        ),
        (
            "TrigToTan[Times[Sin[x], Cos[x]]]",
            "Times[Tan[x], Power[Plus[1, Power[Tan[x], 2]], -1]]",
        ),
        ("TrigToTan[Sin[x]]", "Times[Cos[x], Tan[x]]"),
        (
            "TrigToTan[Times[Power[Sin[x], 2], Power[Cos[x], -2]]]",
            "Power[Tan[x], 2]",
        ),
        (
            "TrigToTan[Plus[Power[Cos[x], 4], Power[Sin[x], 2]]]",
            "Times[Plus[1, Power[Tan[x], 2], Power[Tan[x], 4]], \
             Power[Plus[1, Times[2, Power[Tan[x], 2]], Power[Tan[x], 4]], -1]]",
        ),
        // Entry by entry.
        (
            "TrigToCos[List[Power[Sin[x], 2], 1]]",
            "List[Plus[1, Times[-1, Power[Cos[x], 2]]], 1]",
        ),
        (
            "TrigToCos[Equal[Power[Sin[x], 2], 0]]",
            "Equal[Plus[1, Times[-1, Power[Cos[x], 2]]], 0]",
        ),
        // Powers inside products and quotients too.
        (
            "TrigToCos[Power[Plus[Sin[x], Cos[x]], 3]]",
            "Plus[Sin[x], Times[3, Cos[x]], Times[2, Power[Cos[x], 2], Sin[x]], \
             Times[-2, Power[Cos[x], 3]]]",
        ),
        (
            "TrigToCos[Times[Power[Sin[x], 2], Power[Plus[1, Cos[x]], -1]]]",
            "Plus[1, Times[-1, Cos[x]]]",
        ),
        // What stays: an odd power of the other function.
        ("TrigToSin[Cos[x]]", "Cos[x]"),
        ("TrigToCos[Sin[x]]", "Sin[x]"),
    ]);
    for f in ["TrigToCos", "TrigToSin", "TrigToTan"] {
        check_values(
            f,
            &[
                "Plus[Power[Sin[x], 4], Power[Cos[x], 3]]",
                "Times[Power[Sin[y], -2], Tan[y]]",
                "Plus[Power[Tan[x], 2], Times[Sin[x], Cos[y]]]",
                "Power[Plus[Sin[x], Cos[x]], 3]",
                "Times[x, Power[Cos[Times[2, y]], -3]]",
            ],
        );
    }
}
