//! Tests of the matrix basics (CAS-218): giac 1.9's `tran`, `idn`,
//! `trace`, `dot` and `cross`, and what stays.

use crate::trace::EvalOptions;
use crate::{Context, EvalError, Expr, Warning};

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

const M23: &str = "List[List[1, 2, 3], List[4, 5, 6]]";

#[test]
fn transposes() {
    check(&[
        (
            &format!("Transpose[{M23}]"),
            "List[List[1, 4], List[2, 5], List[3, 6]]",
        ),
        (
            "Transpose[List[1, 2, 3]]",
            "List[List[1], List[2], List[3]]",
        ),
        (
            "Transpose[List[List[a, b], List[c, d]]]",
            "List[List[a, c], List[b, d]]",
        ),
        ("Transpose[x]", "Transpose[x]"),
    ]);
}

#[test]
fn identities_and_traces() {
    check(&[
        (
            "IdentityMatrix[3]",
            "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]",
        ),
        ("IdentityMatrix[1]", "List[List[1]]"),
        ("IdentityMatrix[0]", "IdentityMatrix[0]"),
        ("IdentityMatrix[n]", "IdentityMatrix[n]"),
        ("Tr[List[List[1, 2], List[3, 4]]]", "5"),
        ("Tr[List[List[a, b], List[c, d]]]", "Plus[a, d]"),
        (&format!("Tr[{M23}]"), &format!("Tr[{M23}]")),
    ]);
}

#[test]
fn dot_products() {
    check(&[
        ("Dot[List[1, 2, 3], List[4, 5, 6]]", "32"),
        // The conjugate of the first, as giac's dot([i,0],[1,0]) is -i.
        ("Dot[List[I, 0], List[1, 0]]", "Times[-1, I]"),
        ("Dot[List[1, 0], List[I, 0]]", "I"),
        (
            "Dot[List[Plus[1, I], 2], List[3, Times[4, I]]]",
            "Plus[3, Times[5, I]]",
        ),
        // Real symbols are their own conjugates; of matrices, entrywise.
        (
            "Dot[List[a, b], List[c, d]]",
            "Plus[Times[a, c], Times[b, d]]",
        ),
        ("Dot[List[List[1, 2], List[3, 4]], IdentityMatrix[2]]", "5"),
    ]);
}

#[test]
fn cross_products() {
    check(&[
        ("Cross[List[1, 0, 0], List[0, 1, 0]]", "List[0, 0, 1]"),
        ("Cross[List[1, 2], List[3, 4]]", "List[0, 0, -2]"),
        (
            "Cross[List[a, b, c], List[x, y, z]]",
            "List[Plus[Times[b, z], Times[-1, c, y]], Plus[Times[c, x], Times[-1, a, z]], \
             Plus[Times[a, y], Times[-1, b, x]]]",
        ),
    ]);
}

/// Dimensions that do not fit are shape errors; a trace needs a square.
#[test]
fn shapes() {
    let ctx = Context::new();
    for input in [
        "Dot[List[1, 2], List[1, 2, 3]]",
        "Cross[List[1, 2, 3, 4], List[1, 2, 3, 4]]",
    ] {
        assert_eq!(ctx.eval(&read(input)), Err(EvalError::Shape), "{input}");
    }
    let (_, report) = ctx.eval_with(&read(&format!("Tr[{M23}]")), &EvalOptions::default());
    let messages: Vec<String> = report.warnings.iter().map(Warning::to_string).collect();
    assert!(
        messages.iter().any(|m| m.contains("a square matrix")),
        "{messages:?}"
    );
}

const M2: &str = "List[List[1, 2], List[3, 4]]";

/// Determinants (CAS-219): exact over the generators, numeric with balls.
#[test]
fn determinants() {
    check(&[
        (&format!("Det[{M2}]"), "-2"),
        (
            "Det[List[List[a, b], List[c, d]]]",
            "Plus[Times[a, d], Times[-1, b, c]]",
        ),
        (
            "Det[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 10]]]",
            "-3",
        ),
        ("Det[List[List[1, 2], List[2, 4]]]", "0"),
        ("Det[List[List[1/2, 1/3], List[1/4, 1/5]]]", "1/60"),
        ("Det[List[List[I, 1], List[1, I]]]", "-2"),
        ("Det[List[List[x, 1], List[1, x]]]", "Plus[-1, Power[x, 2]]"),
        ("Det[List[List[2]]]", "2"),
        // A row swap: [[0,1],[1,0]].
        ("Det[List[List[0, 1], List[1, 0]]]", "-1"),
    ]);
}

/// The value of `input`, each entry of a matrix through `Together`.
fn together(ctx: &Context, input: &str) -> Expr {
    let found = ctx.eval(&read(input)).unwrap();
    let rows = found
        .args()
        .iter()
        .map(|r| {
            let entries = r
                .args()
                .iter()
                .map(|x| {
                    ctx.eval(&Expr::call(crate::Symbol::TOGETHER, [x.clone()]))
                        .unwrap()
                })
                .collect::<Vec<_>>();
            Expr::call(crate::Symbol::LIST, entries)
        })
        .collect::<Vec<_>>();
    Expr::call(crate::Symbol::LIST, rows)
}

/// Inverses and negative powers (CAS-219).
#[test]
fn inverses() {
    check(&[
        (
            &format!("Inverse[{M2}]"),
            "List[List[-2, 1], List[3/2, -1/2]]",
        ),
        (
            &format!("Power[{M2}, -1]"),
            "List[List[-2, 1], List[3/2, -1/2]]",
        ),
        (
            &format!("Power[{M2}, -2]"),
            "List[List[11/2, -5/2], List[-15/4, 7/4]]",
        ),
        (
            "Inverse[List[List[1, 2, 3], List[0, 1, 4], List[5, 6, 0]]]",
            "List[List[-24, 18, 5], List[20, -15, -4], List[-5, 4, 1]]",
        ),
        ("Inverse[List[List[2]]]", "List[List[1/2]]"),
        (
            "Inverse[List[List[I, 0], List[0, 2]]]",
            "List[List[Times[-1, I], 0], List[0, 1/2]]",
        ),
    ]);
    // With parameters: the product with the matrix is the identity.
    let ctx = Context::new();
    for m in [
        "List[List[a, b], List[c, d]]",
        "List[List[1, a], List[a, 2]]",
        "List[List[x, 1, 0], List[0, x, 1], List[1, 0, x]]",
    ] {
        let product = together(&ctx, &format!("Times[{m}, Inverse[{m}]]"));
        let n = read(m).args().len();
        let identity = ctx.eval(&read(&format!("IdentityMatrix[{n}]"))).unwrap();
        assert_eq!(product, identity, "{m}");
    }
}

/// A singular matrix stays, with a warning.
#[test]
fn singular_matrices() {
    let ctx = Context::new();
    for input in [
        "Inverse[List[List[1, 2], List[2, 4]]]",
        "Power[List[List[1, 2], List[2, 4]], -1]",
        "Inverse[List[List[a, b], List[Times[2, a], Times[2, b]]]]",
    ] {
        let read = read(input);
        let (found, report) = ctx.eval_with(&read, &EvalOptions::default());
        assert_eq!(found.unwrap(), read, "{input}");
        let messages: Vec<String> = report.warnings.iter().map(Warning::to_string).collect();
        assert!(
            messages.iter().any(|m| m.contains("singular")),
            "{input}: {messages:?}"
        );
    }
}

/// Numeric entries: elimination on balls.
#[test]
fn numeric_matrices() {
    let ctx = Context::new();
    let value = |input: &str| -> f64 {
        let found = ctx.eval(&read(input)).unwrap();
        match found.view() {
            crate::View::Ball(b) => b.mid().to_f64(),
            _ => panic!("{input}: {found:?}"),
        }
    };
    let m = "List[List[Ball[3/2, 0], 2], List[3, Ball[5, 0]]]";
    assert!((value(&format!("Det[{m}]")) - 1.5).abs() < 1e-12);
    let inverse = ctx.eval(&read(&format!("Inverse[{m}]"))).unwrap();
    let entry = |i: usize, j: usize| match inverse.args()[i].args()[j].view() {
        crate::View::Ball(b) => b.mid().to_f64(),
        _ => panic!("{inverse:?}"),
    };
    // [[1.5, 2], [3, 5]]^-1 = [[10/3, -4/3], [-2, 1]].
    for (i, j, x) in [
        (0, 0, 10.0 / 3.0),
        (0, 1, -4.0 / 3.0),
        (1, 0, -2.0),
        (1, 1, 1.0),
    ] {
        assert!((entry(i, j) - x).abs() < 1e-12, "{inverse:?}");
    }
}

const M33: &str = "List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]]";

/// Row reductions (CAS-220): giac's rref, rank, ker and image.
#[test]
fn row_reductions() {
    check(&[
        (
            &format!("RowReduce[{M23}]"),
            "List[List[1, 0, -1], List[0, 1, 2]]",
        ),
        (
            "RowReduce[List[List[a, b], List[c, d]]]",
            "List[List[1, 0], List[0, 1]]",
        ),
        (
            "RowReduce[List[List[0, 1], List[0, 0]]]",
            "List[List[0, 1], List[0, 0]]",
        ),
        (
            &format!("RowReduce[{M33}]"),
            "List[List[1, 0, -1], List[0, 1, 2], List[0, 0, 0]]",
        ),
        ("MatrixRank[List[List[1, 2], List[2, 4]]]", "1"),
        (&format!("MatrixRank[{M33}]"), "2"),
        ("MatrixRank[List[List[0, 0], List[0, 0]]]", "0"),
        ("MatrixRank[List[List[a, b], List[c, d]]]", "2"),
        (
            "NullSpace[List[List[1, 2, 3], List[2, 4, 6]]]",
            "List[List[2, -1, 0], List[3, 0, -1]]",
        ),
        ("NullSpace[List[List[1, 0], List[0, 1]]]", "List[]"),
        (&format!("NullSpace[{M33}]"), "List[List[-1, 2, -1]]"),
        (
            "NullSpace[List[List[a, b], List[Times[2, a], Times[2, b]]]]",
            "List[List[Times[b, Power[a, -1]], -1]]",
        ),
        (
            "NullSpace[List[List[0, 0], List[0, 0]]]",
            "List[List[-1, 0], List[0, -1]]",
        ),
        (
            "ColumnSpace[List[List[1, 2, 3], List[2, 4, 6]]]",
            "List[List[1, 2]]",
        ),
        (
            "ColumnSpace[List[List[1, 0, 1], List[0, 1, 1]]]",
            "List[List[1, 0], List[0, 1]]",
        ),
        // Numeric entries.
        ("MatrixRank[List[List[Ball[1, 0], 2], List[2, 4]]]", "1"),
    ]);
}

/// Linear systems by `LinSolve` (CAS-220): one list.
#[test]
fn linear_solutions() {
    let two = "List[Equation[Plus[x, y], 1], Equation[Plus[x, Times[-1, y]], 3]]";
    check(&[
        (&format!("LinSolve[{two}, List[x, y]]"), "List[2, -1]"),
        (
            "LinSolve[List[Equation[Plus[x, y], 1], Equation[Plus[Times[2, x], Times[2, y]], 2]], List[x, y]]",
            "List[Plus[1, Times[-1, y]], y]",
        ),
        (
            "LinSolve[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2]], List[x, y]]",
            "List[]",
        ),
        (
            "LinSolve[List[List[1, 1], List[1, -1]], List[1, 3]]",
            "List[2, -1]",
        ),
        (
            "LinSolve[List[List[1, 2], List[3, 4]], List[5, 6]]",
            "List[-4, 9/2]",
        ),
        (
            "LinSolve[List[List[1, 1], List[1, 1]], List[1, 2]]",
            "List[]",
        ),
        // A family in matrix form stays, with a warning.
        (
            "LinSolve[List[List[1, 1], List[1, 1]], List[1, 1]]",
            "LinSolve[List[List[1, 1], List[1, 1]], List[1, 1]]",
        ),
    ]);
}
