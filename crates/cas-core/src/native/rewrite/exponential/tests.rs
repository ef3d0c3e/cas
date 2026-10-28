//! Tests of the exponential, power and logarithm forms (CAS-396): giac's
//! forms where they are plain, values at a point, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn hyperbolic() {
    check(&[
        (
            "HyperbolicToExp[Sinh[x]]",
            "Times[1/2, Plus[Exp[x], Times[-1, Exp[Times[-1, x]]]]]",
        ),
        (
            "HyperbolicToExp[ArcSinh[x]]",
            "Log[Plus[x, Power[Plus[1, Power[x, 2]], 1/2]]]",
        ),
        // The circular functions are `TrigToExp`'s.
        ("HyperbolicToExp[Sin[x]]", "Sin[x]"),
    ]);
    // `ArcCosh` below 1 and below −1 too, where giac's form is not it.
    check_values(
        "HyperbolicToExp",
        &[
            "Sinh[x]",
            "Cosh[Times[2, y]]",
            "Tanh[Plus[x, y]]",
            "ArcSinh[x]",
            "ArcSinh[Times[-1, y]]",
            "ArcCosh[y]",
            "ArcCosh[x]",
            "ArcCosh[Times[-1, y]]",
            "ArcTanh[x]",
            "Times[Sinh[x], ArcSinh[y]]",
        ],
    );
}

#[test]
fn powers() {
    check(&[
        ("ExpToPower[Exp[Times[3, Log[x]]]]", "Power[x, 3]"),
        (
            "ExpToPower[Exp[Plus[Times[2, Log[x]], Times[3, Log[y]]]]]",
            "Times[Power[x, 2], Power[y, 3]]",
        ),
        ("ExpToPower[Exp[Times[x, Log[x]]]]", "Power[x, x]"),
        (
            "ExpToPower[Exp[Plus[1, Times[n, Log[x]]]]]",
            "Times[Exp[1], Power[x, n]]",
        ),
        ("ExpToPower[Power[E, Times[1/2, Log[x]]]]", "Power[x, 1/2]"),
        ("PowerToExp[Power[x, y]]", "Exp[Times[y, Log[x]]]"),
        ("PowerToExp[Power[2, x]]", "Exp[Times[x, Log[2]]]"),
        ("PowerToExp[Power[E, x]]", "Exp[x]"),
        // Rational exponents stay, as giac's do.
        ("PowerToExp[Power[x, 2]]", "Power[x, 2]"),
        ("PowerToExp[Power[x, 1/2]]", "Power[x, 1/2]"),
    ]);
    check_values(
        "ExpToPower",
        &[
            "Exp[Times[3, Log[x]]]",
            "Exp[Plus[1, Times[y, Log[x]], Times[-2, Log[y]]]]",
            "Power[E, Times[1/2, Log[x]]]",
        ],
    );
    check_values(
        "PowerToExp",
        &[
            "Power[x, y]",
            "Power[2, x]",
            "Power[Plus[1, x], Power[y, 1/2]]",
        ],
    );
}

#[test]
fn exponentials() {
    check(&[
        ("ExpExpand[Exp[Plus[x, y]]]", "Times[Exp[x], Exp[y]]"),
        ("ExpExpand[Exp[Times[3, x]]]", "Power[Exp[x], 3]"),
        ("ExpExpand[Exp[Times[-2, x]]]", "Power[Exp[x], -2]"),
        (
            "ExpExpand[Exp[Times[-1, I, x]]]",
            "Power[Exp[Times[I, x]], -1]",
        ),
        (
            "ExpExpand[Exp[Times[3/2, x]]]",
            "Power[Exp[Times[1/2, x]], 3]",
        ),
        (
            "ExpExpand[Exp[Plus[1, Times[3, x]]]]",
            "Times[Exp[1], Power[Exp[x], 3]]",
        ),
        ("ExpExpand[Power[E, Plus[x, y]]]", "Times[Exp[x], Exp[y]]"),
        // A primitive coefficient stays.
        ("ExpExpand[Exp[Times[1/2, x]]]", "Exp[Times[1/2, x]]"),
    ]);
    check_values(
        "ExpExpand",
        &[
            "Exp[Plus[x, Times[-2, y]]]",
            "Exp[Times[-6/5, x]]",
            "Exp[Times[Complex[2, -4], x]]",
            "Power[E, Plus[1, x, y]]",
        ],
    );
}

#[test]
fn logarithms_expanded() {
    check(&[
        ("LogExpand[Log[Times[x, y]]]", "Plus[Log[x], Log[y]]"),
        (
            "LogExpand[Log[Times[Power[x, 2], Power[y, -1]]]]",
            "Plus[Times[2, Log[Abs[x]]], Times[-1, Log[y]]]",
        ),
        ("LogExpand[Log[Power[x, 3]]]", "Times[3, Log[x]]"),
        ("LogExpand[Log[Power[x, 2/3]]]", "Times[2/3, Log[Abs[x]]]"),
        ("LogExpand[Log[Power[x, 1/2]]]", "Times[1/2, Log[x]]"),
        ("LogExpand[Log[Power[x, n]]]", "Times[n, Log[x]]"),
        ("LogExpand[Log[2/3]]", "Plus[Log[2], Times[-1, Log[3]]]"),
        (
            "LogExpand[Log[Times[2, x, Power[y, 3]]]]",
            "Plus[Log[2], Log[x], Times[3, Log[y]]]",
        ),
        ("LogExpand[Log[Times[Exp[x], y]]]", "Plus[x, Log[y]]"),
        // What stays: an integer, a negative factor, a power of a logarithm.
        ("LogExpand[Log[12]]", "Log[12]"),
        ("LogExpand[Log[Times[-2, x]]]", "Log[Times[-2, x]]"),
        ("LogExpand[Power[Log[x], 2]]", "Power[Log[x], 2]"),
    ]);
    check_values(
        "LogExpand",
        &[
            "Log[Times[x, Power[y, 3]]]",
            "Log[Times[3/5, Power[x, -2], Power[y, 1/2]]]",
            "Log[Power[Times[x, y], 1/3]]",
            "Times[y, Log[Times[Exp[x], y]]]",
        ],
    );
}

#[test]
fn logarithms_combined() {
    check(&[
        ("LogCombine[Plus[Log[x], Log[y]]]", "Log[Times[x, y]]"),
        (
            "LogCombine[Plus[Times[2, Log[x]], Times[-1, Log[y]]]]",
            "Log[Times[Power[x, 2], Power[y, -1]]]",
        ),
        (
            "LogCombine[Plus[Times[1/2, Log[x]], Times[1/3, Log[y]]]]",
            "Times[1/6, Log[Times[Power[x, 3], Power[y, 2]]]]",
        ),
        (
            "LogCombine[Plus[Times[-1, Log[x]], Times[-1, Log[y]]]]",
            "Times[-1, Log[Times[x, y]]]",
        ),
        ("LogCombine[Times[2, Log[x]]]", "Log[Power[x, 2]]"),
        (
            "LogCombine[Times[-2, Log[x]]]",
            "Times[-1, Log[Power[x, 2]]]",
        ),
        ("LogCombine[Times[3, Log[2]]]", "Log[8]"),
        ("LogCombine[Plus[Log[2], Log[3]]]", "Log[6]"),
        (
            "LogCombine[Plus[Times[a, Log[x]], Times[a, Log[y]]]]",
            "Times[a, Log[Times[x, y]]]",
        ),
        (
            "LogCombine[Times[Plus[Log[x], Log[y]], Power[z, -1]]]",
            "Times[Log[Times[x, y]], Power[z, -1]]",
        ),
        (
            "LogCombine[Equal[Plus[Log[x], Log[y]], 0]]",
            "Equal[Log[Times[x, y]], 0]",
        ),
        (
            "LogCombine[Plus[3, Log[x], Times[-1, Log[y]]]]",
            "Plus[3, Log[Times[x, Power[y, -1]]]]",
        ),
        // What stays: groups of one term, `-ln(x)`, a product of logarithms.
        (
            "LogCombine[Plus[Times[a, Log[x]], Log[y]]]",
            "Plus[Times[a, Log[x]], Log[y]]",
        ),
        ("LogCombine[Times[-1, Log[x]]]", "Times[-1, Log[x]]"),
        ("LogCombine[Times[Log[x], Log[y]]]", "Times[Log[x], Log[y]]"),
    ]);
    check_values(
        "LogCombine",
        &[
            "Plus[Log[x], Times[2, Log[y]]]",
            "Plus[Times[1/2, Log[x]], Times[-1/3, Log[y]], Log[5]]",
            "Plus[Times[y, Log[x]], Times[y, Log[3]], Log[y]]",
            "Times[-4, Log[x]]",
        ],
    );
}
