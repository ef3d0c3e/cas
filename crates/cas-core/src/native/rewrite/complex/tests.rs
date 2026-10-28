//! Tests of the complex parts (CAS-401): giac's `evalc` answers, values at
//! a point, the conjugate multipliers, and what stays.

use super::super::testing::{check, check_values};

#[test]
fn expansions() {
    check(&[
        (
            "ComplexExpand[Power[Plus[x, Times[I, y]], -1]]",
            "Plus[Times[x, Power[Plus[Power[x, 2], Power[y, 2]], -1]], \
             Times[Complex[0, -1], y, Power[Plus[Power[x, 2], Power[y, 2]], -1]]]",
        ),
        (
            "ComplexExpand[Power[Plus[x, Times[I, y]], 2]]",
            "Plus[Power[x, 2], Times[-1, Power[y, 2]], Times[Complex[0, 2], x, y]]",
        ),
        (
            "ComplexExpand[Exp[Plus[x, Times[I, y]]]]",
            "Plus[Times[Exp[x], Cos[y]], Times[Complex[0, 1], Exp[x], Sin[y]]]",
        ),
        (
            "ComplexExpand[Sin[Plus[x, Times[I, y]]]]",
            "Plus[Times[Cosh[y], Sin[x]], Times[Complex[0, 1], Cos[x], Sinh[y]]]",
        ),
        (
            "ComplexExpand[Abs[Plus[x, Times[I, y]]]]",
            "Power[Plus[Power[x, 2], Power[y, 2]], 1/2]",
        ),
        (
            "ComplexExpand[Conjugate[Plus[x, Times[I, y]]]]",
            "Plus[x, Times[Complex[0, -1], y]]",
        ),
        // giac's example: re(x+i*y)+im(x*y) is x.
        (
            "ComplexExpand[Plus[Re[Plus[x, Times[I, y]]], Im[Times[x, y]]]]",
            "x",
        ),
        // Entry by entry; what stays: an unknown function.
        (
            "ComplexExpand[List[Plus[x, Times[I, y]], Conjugate[x]]]",
            "List[Plus[x, Times[I, y]], x]",
        ),
        ("ComplexExpand[f[x]]", "f[x]"),
    ]);
    check_values(
        "ComplexExpand",
        &[
            "Power[Plus[x, Times[I, y]], -2]",
            "Power[Plus[x, Times[I, y]], 3]",
            "Tan[Plus[x, Times[I, y]]]",
            "Tanh[Plus[y, Times[I, x]]]",
            "Cosh[Plus[x, Times[I, y]]]",
            "Cos[Plus[x, Times[-1, I, y]]]",
            "Sinh[Times[I, y]]",
            "Log[Plus[x, Times[I, y]]]",
            "Log[Times[-1, x]]",
            "Power[Times[-1, y], 1/3]",
            "Power[Plus[x, Times[I, y]], Plus[1/3, I]]",
            "Power[Plus[1, I], x]",
            "Times[Exp[Times[I, x]], Exp[Times[I, y]]]",
            "Power[x, Plus[y, I]]",
            "Abs[Plus[Exp[Times[I, x]], y]]",
        ],
    );
}

#[test]
fn conjugates() {
    check(&[
        (
            "MultiplyConjugate[Power[Plus[Power[3, 1/2], Times[-1, Power[2, 1/2]]], -1]]",
            "Plus[Power[2, 1/2], Power[3, 1/2]]",
        ),
        (
            "MultiplyComplexConjugate[Power[Plus[x, Times[I, y]], -1]]",
            "Times[Plus[x, Times[Complex[0, -1], y]], \
             Power[Plus[Power[x, 2], Power[y, 2]], -1]]",
        ),
        // What stays: no square root, and a real denominator.
        (
            "MultiplyConjugate[Power[Plus[1, x], -1]]",
            "Power[Plus[1, x], -1]",
        ),
        (
            "MultiplyComplexConjugate[Power[Plus[1, x], -1]]",
            "Power[Plus[1, x], -1]",
        ),
    ]);
    check_values(
        "MultiplyConjugate",
        &[
            "Power[Plus[Power[3, 1/2], Times[-1, Power[2, 1/2]]], -1]",
            "Times[x, Power[Plus[-2, Power[x, 1/2]], -1]]",
            "Power[Plus[1, Power[x, 1/2], Power[y, 1/2]], -1]",
            "Plus[Power[3, 1/2], Times[-1, Power[2, 1/2]]]",
            "Times[y, Power[Plus[Power[y, 1/2], Times[-1, x]], -1], Power[Plus[1, x], -1]]",
        ],
    );
    check_values(
        "MultiplyComplexConjugate",
        &[
            "Power[Plus[x, Times[I, y]], -1]",
            "Times[Plus[1, x], Power[Plus[2, Times[I, x]], -1]]",
            "Plus[x, Times[I, y]]",
            "Times[y, Power[Plus[Exp[Times[I, x]], 1], -1]]",
        ],
    );
}
