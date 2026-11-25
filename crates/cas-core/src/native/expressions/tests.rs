//! Tests of Xcas's expression commands, second part (CAS-568 to CAS-571):
//! giac's answers in full form, ours in canonical order, and what stays.

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

/// `lname`, `indets`, `lvar` and `algvar` (CAS-568): giac's answers, in
/// our order of the operands.
#[test]
fn identifiers_and_kernels() {
    check(&[
        ("Identifiers[Times[2, Exp[x], Sin[y]]]", "List[x, y]"),
        ("Identifiers[Plus[Times[x, y], z]]", "List[x, y, z]"),
        ("Identifiers[Plus[Pi, Times[Power[E, y], x]]]", "List[y, x]"),
        ("Identifiers[f[x, y]]", "List[f, x, y]"),
        ("Identifiers[3]", "List[]"),
        ("Kernels[Times[2, Exp[x], Sin[y]]]", "List[Exp[x], Sin[y]]"),
        (
            "Kernels[Plus[Log[x], Times[2, Exp[x], Sin[y]]]]",
            "List[Log[x], Exp[x], Sin[y]]",
        ),
        (
            "Kernels[Plus[Power[x, 1/2], Times[Power[x, 1/3], y]]]",
            "List[Power[x, 1/2], Power[x, 1/3], y]",
        ),
        ("Kernels[Power[Plus[1, x], -1]]", "List[x]"),
        ("Kernels[Plus[Power[x, 2], y]]", "List[x, y]"),
        ("Kernels[Power[x, n]]", "List[Power[x, n]]"),
        ("Kernels[List[x, Sin[y]]]", "List[x, Sin[y]]"),
        ("Kernels[Plus[1, Indeterminate]]", "List[]"),
        (
            "AlgebraicVariables[Plus[Power[x, 1/2], y]]",
            "List[List[y], List[x]]",
        ),
        (
            "AlgebraicVariables[Plus[Power[x, 1/2], Power[y, 1/2], z]]",
            "List[List[z], List[y], List[x]]",
        ),
        ("AlgebraicVariables[Plus[x, y]]", "List[List[x, y]]"),
        (
            "AlgebraicVariables[Power[Plus[x, Power[y, 1/2]], 1/2]]",
            "List[List[], List[x], List[y]]",
        ),
        (
            "AlgebraicVariables[Plus[Power[x, 1/7], y]]",
            "List[List[Power[x, 1/7], y]]",
        ),
    ]);
}

/// `reorder` (CAS-568) and `part` (CAS-570).
#[test]
fn reorder_and_part() {
    check(&[
        (
            "Reorder[Plus[Power[x, 2], Times[2, x], Power[y, 2]], List[y, x]]",
            "Plus[Power[x, 2], Times[2, x], Power[y, 2]]",
        ),
        (
            "Reorder[Power[Plus[1, x], 2], List[x]]",
            "Plus[1, Times[2, x], Power[x, 2]]",
        ),
        (
            "Reorder[Plus[y, Power[Plus[1, x], -1]], List[y, x]]",
            "Times[Plus[1, y, Times[x, y]], Power[Plus[1, x], -1]]",
        ),
        ("Reorder[3, List[x]]", "3"),
        ("Reorder[Plus[1, x], y]", "Reorder[Plus[1, x], y]"),
        ("Part[Times[Exp[x], Sin[x]], 1]", "Exp[x]"),
        ("Part[Plus[x, y, z], 3]", "z"),
    ]);
}

/// `eval`, `evalm` and `evalb` (CAS-569).
#[test]
fn evaluations() {
    check(&[
        ("Evaluate[Times[2, Sin[Pi]]]", "0"),
        ("Evaluate[Plus[x, x], 1]", "Times[2, x]"),
        (
            "Evaluate[Power[List[List[1, 2], List[3, 4]], 2]]",
            "List[List[7, 10], List[15, 22]]",
        ),
        ("EvaluateBoolean[Less[1, 2]]", "True"),
        ("EvaluateBoolean[Greater[2, 3]]", "False"),
        ("EvaluateBoolean[Equal[a, 2]]", "False"),
        ("EvaluateBoolean[Equation[1, 1]]", "True"),
        ("EvaluateBoolean[Equation[Plus[1, x], Plus[x, 1]]]", "True"),
        ("EvaluateBoolean[Less[x, Plus[1, x]]]", "True"),
        ("EvaluateBoolean[Greater[Pi, 3]]", "True"),
        ("EvaluateBoolean[Greater[Power[2, 1/2], 1]]", "True"),
        ("EvaluateBoolean[0]", "False"),
        ("EvaluateBoolean[2]", "True"),
        (
            "EvaluateBoolean[List[Less[1, 2], Less[2, 1]]]",
            "List[True, False]",
        ),
        ("EvaluateBoolean[Less[a, b]]", "EvaluateBoolean[Less[a, b]]"),
        (
            "EvaluateBoolean[Greater[Plus[Power[2, 1/2], Pi], a]]",
            "EvaluateBoolean[Greater[Plus[Power[2, 1/2], Pi], a]]",
        ),
    ]);
}

/// `epsilon2zero` (CAS-569): decimals held apart, giac's bound `10⁻¹²`
/// included, and a bound given.
#[test]
fn epsilon_to_zero() {
    check(&[
        ("EpsilonToZero[Plus[Decimal[1, -13], x, 5]]", "Plus[5, x]"),
        (
            "EpsilonToZero[List[Decimal[1, -12], Decimal[25, -1]]]",
            "List[0, Decimal[25, -1]]",
        ),
        ("EpsilonToZero[Decimal[11, -13]]", "Decimal[11, -13]"),
        (
            "EpsilonToZero[Plus[Decimal[1, -5], x], Decimal[1, -4]]",
            "x",
        ),
        (
            "EpsilonToZero[Plus[Decimal[1, -5], x], Decimal[1, -6]]",
            "Plus[Decimal[1, -5], x]",
        ),
        (
            "EpsilonToZero[Plus[Decimal[1, -13], Times[Decimal[2, -11], I]]]",
            "Times[Decimal[2, -11], I]",
        ),
        ("EpsilonToZero[x]", "x"),
    ]);
}

/// `isinf` and `isnan` (CAS-569).
#[test]
fn infinities() {
    check(&[
        ("InfiniteQ[DirectedInfinity[-1]]", "True"),
        ("InfiniteQ[Infinity]", "True"),
        ("InfiniteQ[ComplexInfinity]", "True"),
        ("InfiniteQ[1]", "False"),
        ("InfiniteQ[x]", "False"),
        ("InfiniteQ[Indeterminate]", "False"),
        ("IndeterminateQ[Indeterminate]", "True"),
        ("IndeterminateQ[0]", "False"),
        ("IndeterminateQ[Infinity]", "False"),
        ("IndeterminateQ[List[Indeterminate, 1]]", "True"),
        ("IndeterminateQ[List[1, Indeterminate]]", "False"),
    ]);
}

/// `compare` (CAS-569): giac's `islesscomplexthan`, types first.
#[test]
fn comparisons() {
    check(&[
        ("Compare[1, 2]", "1"),
        ("Compare[2, 1]", "0"),
        ("Compare[1, 1]", "0"),
        ("Compare[-1, 1]", "0"),
        ("Compare[1, -1]", "1"),
        ("Compare[N[1], 2]", "0"),
        ("Compare[2, N[1]]", "1"),
        ("Compare[Decimal[15, -1], Decimal[25, -1]]", "0"),
        ("Compare[\"ab\", \"cd\"]", "1"),
        ("Compare[\"ab\", \"ab\"]", "0"),
        ("Compare[1, \"ab\"]", "1"),
        ("Compare[x, y]", "1"),
        ("Compare[y, x]", "0"),
        ("Compare[1, x]", "1"),
        ("Compare[x, 1]", "0"),
        ("Compare[0, x]", "0"),
        ("Compare[x, 0]", "1"),
        ("Compare[1/2, 1]", "0"),
        ("Compare[1/2, 1/3]", "1"),
        ("Compare[Power[2, 40], 3]", "0"),
        ("Compare[3, Power[2, 40]]", "1"),
        ("Compare[Complex[1, 1], 2]", "0"),
        ("Compare[Complex[1, 2], Complex[2, 1]]", "0"),
        ("Compare[List[1, 2], List[1, 3]]", "1"),
        ("Compare[List[1, 2], List[3]]", "0"),
        ("Compare[Power[x, 2], x]", "0"),
        ("Compare[Sin[x], Sin[y]]", "1"),
        ("Compare[Plus[1, x], Plus[-1, x]]", "1"),
        ("Compare[Cos[x], Sin[x]]", "1"),
    ]);
}

/// `equal2diff`, `equal2list`, `list2exp` and `exp2list` (CAS-571), with
/// both `=` and `==`.
#[test]
fn equations_and_lists() {
    check(&[
        ("EquationToDifference[Equation[x, 2]]", "Plus[-2, x]"),
        (
            "EquationToDifference[Equal[x, Plus[1, y]]]",
            "Plus[x, Times[-1, Plus[1, y]]]",
        ),
        (
            "EquationToDifference[List[Equation[x, 1], Equation[y, 2]]]",
            "List[Plus[-1, x], Plus[-2, y]]",
        ),
        ("EquationToDifference[x]", "x"),
        ("List[EquationToList[Equation[x, 2]]]", "List[x, 2]"),
        (
            "List[EquationToList[Equal[Plus[1, x], y]]]",
            "List[Plus[1, x], y]",
        ),
        ("EquationToList[x]", "List[x, 0]"),
        (
            "EquationToList[List[Equation[x, 1], Equation[y, 2]]]",
            "List[x, 1, y, 2]",
        ),
        (
            "ListToEquations[List[List[3, 9], List[-1, 1]], List[x, y]]",
            "Or[And[Equation[x, 3], Equation[y, 9]], And[Equation[x, -1], Equation[y, 1]]]",
        ),
        (
            "ListToEquations[List[1, 2], x]",
            "Or[Equation[x, 1], Equation[x, 2]]",
        ),
        ("ListToEquations[List[1], x]", "Equation[x, 1]"),
        (
            "ListToEquations[List[List[1, 2]], List[x, y]]",
            "And[Equation[x, 1], Equation[y, 2]]",
        ),
        (
            "ListToEquations[List[1, 2], List[x, y]]",
            "And[Equation[x, 1], Equation[y, 2]]",
        ),
        ("ListToEquations[List[], x]", "False"),
        (
            "EquationsToList[Or[Equation[x, 2], Equation[x, 0]]]",
            "List[2, 0]",
        ),
        (
            "EquationsToList[Or[And[Equation[x, 3], Equation[y, 9]], And[Equation[x, -1], Equation[y, 1]]]]",
            "List[List[3, 9], List[-1, 1]]",
        ),
        ("EquationsToList[Equation[x, 2]]", "List[2]"),
        ("EquationsToList[Equal[x, 2]]", "List[2]"),
        (
            "EquationsToList[List[Equation[x, 1], Equation[x, 2]]]",
            "List[1, 2]",
        ),
        ("EquationsToList[0]", "List[]"),
        ("EquationsToList[x]", "List[x]"),
    ]);
}

/// `epsilon2zero`'s default bound: the context's `epsilon` (CAS-617).
#[test]
fn epsilon_to_zero_reads_the_context() {
    let mut ctx = Context::new();
    ctx.set_epsilon(1e-3);
    for (input, expected) in [
        ("EpsilonToZero[Plus[0x1p-14, x]]", "x"),
        ("EpsilonToZero[Plus[0x1p-7, x]]", "Plus[0x1p-7, x]"),
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        assert_eq!(found, ctx.eval(&read(expected)).unwrap(), "{input}");
    }
}
