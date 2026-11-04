//! Tests of `Combine`, `TranscendentalSimplify` and the Heaviside forms
//! (CAS-402): giac's answers, values at a point, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn combinations() {
    check(&[
        ("Combine[Plus[Log[x], Log[y]], Log]", "Log[Times[x, y]]"),
        (
            "Combine[Plus[Times[2, Log[x]], Times[-1, Log[y]]], Log]",
            "Log[Times[Power[x, 2], Power[y, -1]]]",
        ),
        ("Combine[Times[Exp[x], Exp[y]], Exp]", "Exp[Plus[x, y]]"),
        ("Combine[Power[Exp[x], 2], Exp]", "Exp[Times[2, x]]"),
        (
            "Combine[Times[Exp[x], Power[Exp[y], 2]], Exp]",
            "Exp[Plus[x, Times[2, y]]]",
        ),
        (
            "Combine[Times[Sin[x], Cos[x]], trig]",
            "Times[1/2, Sin[Times[2, x]]]",
        ),
        (
            "Combine[Power[Sin[x], 2], Sin]",
            "Plus[1/2, Times[-1/2, Cos[Times[2, x]]]]",
        ),
        // Another class leaves the call.
        (
            "Combine[Plus[Log[x], Log[y]], foo]",
            "Combine[Plus[Log[x], Log[y]], foo]",
        ),
    ]);
}

/// Sums whose exponential forms have degrees near 150 in two variables,
/// found by the rewrite fuzz target: `Together`'s gcd over ℤ\[i\] took 25 s
/// to past 200 s by subresultants, where giac's `tsimplify` takes 50 ms
/// (CAS-494).
#[test]
fn transcendental_sums_of_high_degree() {
    check_values(
        "TranscendentalSimplify",
        &[
            "Plus[Times[-6, Power[Cos[Times[6, y]], 3]], Times[-3, Power[Sin[Times[2, y]], 3]], \
             Times[1/8, Power[Cos[Times[3, x]], 3], Tan[Times[3/4, y]]]]",
            "Plus[Times[-3/2, Cos[Times[5, y]]], Times[1/2, Power[Sin[Times[1/2, y]], 3]], \
             Times[1/6, Power[Cos[Times[3, x]], 3], Power[Tan[Times[3/4, y]], 2]]]",
            "Plus[Power[Cos[Times[5, y]], 2], Times[1/3, Power[Tan[Times[1/3, y]], 2]], \
             Times[1/8, Power[Cos[Times[3, x]], 3], Tan[Times[3/4, y]]]]",
            "Plus[Times[-3, Power[Cos[Times[5/2, y]], 2]], \
             Times[1/4, Power[Cos[Times[3, x]], 3], Power[Tan[Times[3/4, y]], 3]], \
             Times[1/3, Power[Tan[Times[4/3, y]], 3]]]",
        ],
    );
}

#[test]
fn transcendental() {
    check(&[
        (
            "TranscendentalSimplify[Plus[Exp[Times[2, x]], Times[2, Exp[x]], 1]]",
            "Plus[1, Times[2, Exp[x]], Power[Exp[x], 2]]",
        ),
        (
            "TranscendentalSimplify[Plus[Power[Sin[x], 2], Power[Cos[x], 2]]]",
            "1",
        ),
        (
            "TranscendentalSimplify[Plus[Power[Cosh[x], 2], Times[-1, Power[Sinh[x], 2]]]]",
            "1",
        ),
    ]);
    check_values(
        "TranscendentalSimplify",
        &[
            "Plus[Sin[x], Cos[Times[2, y]]]",
            "Tanh[x]",
            "Plus[ArcSinh[x], Exp[Times[-3, y]]]",
            "Times[Sin[x], Power[Cos[y], -1]]",
        ],
    );
}

#[test]
fn steps_and_signs() {
    check(&[
        (
            "SignToHeaviside[Sign[x]]",
            "Plus[-1, Times[2, HeavisideTheta[x]]]",
        ),
        (
            "HeavisideToSign[HeavisideTheta[x]]",
            "Times[1/2, Plus[1, Sign[x]]]",
        ),
    ]);
    check_values(
        "SignToHeaviside",
        &["Sign[Plus[x, -1]]", "Times[y, Sign[x]]", "Sign[Plus[x, y]]"],
    );
    check_values(
        "HeavisideToSign",
        &[
            "HeavisideTheta[Plus[x, -1]]",
            "HeavisideTheta[y]",
            "Times[x, HeavisideTheta[Plus[y, -7/3]]]",
        ],
    );
}
