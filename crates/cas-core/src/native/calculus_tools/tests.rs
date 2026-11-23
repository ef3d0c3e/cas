//! Tests of `IntegrateByPartsDV`, `IntegrateByPartsU`, `Preval` and
//! `FunctionDiff` (CAS-574, CAS-577, CAS-578): giac's examples, the bounds,
//! the finishing step, and what stays.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
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

/// giac's `ibpdv` examples: one step, the bounds, and the finish.
#[test]
fn by_parts_dv() {
    check(&[
        (
            "IntegrateByPartsDV[Log[x], x]",
            "List[Times[x, Log[x]], -1]",
        ),
        (
            "IntegrateByPartsDV[Times[x, Log[x]], Times[1/2, Power[x, 2]]]",
            "List[Times[1/2, Power[x, 2], Log[x]], Times[-1/2, x]]",
        ),
        (
            "IntegrateByPartsDV[Log[x], x, x, 1, 3]",
            "List[Times[3, Log[3]], -1]",
        ),
        (
            "IntegrateByPartsDV[List[Times[x, Log[x]], -1], 0]",
            "Plus[Times[-1, x], Times[x, Log[x]]]",
        ),
        (
            "IntegrateByPartsDV[IntegrateByPartsDV[Log[x], x, x, 2, 3], 0, x, 2, 3]",
            "Plus[-1, Times[3, Log[3]], Times[-2, Log[2]]]",
        ),
        // Another variable.
        (
            "IntegrateByPartsDV[Log[y], y, y]",
            "List[Times[y, Log[y]], -1]",
        ),
        // Not a symbol, too few or too many arguments, a pair not of two.
        (
            "IntegrateByPartsDV[Log[x], x, 2]",
            "IntegrateByPartsDV[Log[x], x, 2]",
        ),
        ("IntegrateByPartsDV[Log[x]]", "IntegrateByPartsDV[Log[x]]"),
        (
            "IntegrateByPartsDV[List[1, 2, 3], 0]",
            "IntegrateByPartsDV[List[1, 2, 3], 0]",
        ),
    ]);
}

/// giac's `ibpu` examples.
#[test]
fn by_parts_u() {
    check(&[
        (
            "IntegrateByPartsU[Log[x], Log[x]]",
            "List[Times[x, Log[x]], -1]",
        ),
        (
            "IntegrateByPartsU[Log[x], Log[x], x, 1, 3]",
            "List[Times[3, Log[3]], -1]",
        ),
        (
            "IntegrateByPartsU[Times[x, Log[x]], Log[x]]",
            "List[Times[1/2, Power[x, 2], Log[x]], Times[-1/2, x]]",
        ),
        (
            "IntegrateByPartsU[List[Times[x, Log[x]], -1], 0]",
            "Plus[Times[-1, x], Times[x, Log[x]]]",
        ),
        (
            "IntegrateByPartsU[IntegrateByPartsU[Log[x], Log[x], x, 2, 3], 0, x, 2, 3]",
            "Plus[-1, Times[3, Log[3]], Times[-2, Log[2]]]",
        ),
        (
            "IntegrateByPartsU[Power[x, 2], 1]",
            "Times[1/3, Power[x, 3]]",
        ),
    ]);
}

/// `preval`: the default variable, another one, `x = a`, and a limit at a
/// removable singularity.
#[test]
fn preval() {
    check(&[
        ("Preval[Plus[Power[x, 2], -2], 2, 3]", "5"),
        ("Preval[Plus[Power[y, 2], -2], 2, 3, y]", "5"),
        ("Preval[Power[x, 2], Equation[x, 1], 2]", "3"),
        ("Preval[Integrate[y, y], 0, 1, y]", "1/2"),
        ("Preval[Times[Sin[x], Power[x, -1]], 0, Pi]", "-1"),
        (
            "Preval[Power[x, 2], 1, 2, 3]",
            "Preval[Power[x, 2], 1, 2, 3]",
        ),
        ("Preval[x, 1]", "Preval[x, 1]"),
    ]);
}

/// `function_diff`: builtins, sums and products of them, functions, lists,
/// and what stays.
#[test]
fn function_diff() {
    check(&[
        ("FunctionDiff[Sin]", "Function[List[x], Cos[x]]"),
        (
            "FunctionDiff[Log][List[x, y]]",
            "List[Power[x, -1], Power[y, -1]]",
        ),
        (
            "FunctionDiff[Plus[Sin, Cos]][t]",
            "Plus[Cos[t], Times[-1, Sin[t]]]",
        ),
        ("FunctionDiff[Times[2, Exp]][t]", "Times[2, Exp[t]]"),
        ("FunctionDiff[Function[List[y], Power[y, 3]]][2]", "12"),
        // A free x in the derivative: the parameter is x1.
        (
            "FunctionDiff[Function[List[y], Times[x, Power[y, 2]]]]",
            "Function[List[x1], Times[2, x1, x]]",
        ),
        // Unknown functions stay.
        ("FunctionDiff[f]", "FunctionDiff[f]"),
        ("FunctionDiff[Plus[Sin, id]]", "FunctionDiff[Plus[Sin, id]]"),
    ]);
}
