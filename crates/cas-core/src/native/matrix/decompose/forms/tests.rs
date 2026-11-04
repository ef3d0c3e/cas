//! Tests of the Hessenberg form and the constructions (CAS-476): giac's
//! answers, the similarity, and the resultant.

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

#[test]
fn hessenberg() {
    let a = "List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 10]]";
    check(&[
        // giac's P and H.
        (
            &format!("HessenbergDecomposition[{a}]"),
            "List[List[List[1, 0, 0], List[0, 4/7, 1], List[0, 1, 0]], \
             List[List[1, 29/7, 2], List[7, 102/7, 8], List[0, 26/49, 3/7]]]",
        ),
        // P H P⁻¹ = A, with a symbol.
        (
            "Together[Times[Part[HessenbergDecomposition[List[List[a, 1, 2], List[1, 0, 1], \
             List[2, 1, 3]]], 1], Part[HessenbergDecomposition[List[List[a, 1, 2], List[1, 0, 1], \
             List[2, 1, 3]]], 2], Inverse[Part[HessenbergDecomposition[List[List[a, 1, 2], \
             List[1, 0, 1], List[2, 1, 3]]], 1]]]]",
            "List[List[a, 1, 2], List[1, 0, 1], List[2, 1, 3]]",
        ),
        // Already Hessenberg: P = I.
        (
            "HessenbergDecomposition[List[List[1, 2], List[3, 4]]]",
            "List[List[List[1, 0], List[0, 1]], List[List[1, 2], List[3, 4]]]",
        ),
    ]);
}

#[test]
fn constructions() {
    check(&[
        // giac's [[0,0,-5],[1,0,-2],[0,1,0]] and [[0,-3/2],[1,0]].
        (
            "CompanionMatrix[Plus[Power[x, 3], Times[2, x], 5], x]",
            "List[List[0, 0, -5], List[1, 0, -2], List[0, 1, 0]]",
        ),
        (
            "CompanionMatrix[Plus[Times[2, Power[x, 2]], 3], x]",
            "List[List[0, -3/2], List[1, 0]]",
        ),
        // giac's sylvester matrices; the determinant is the resultant.
        (
            "SylvesterMatrix[Plus[Power[x, 2], 1], Plus[x, 2], x]",
            "List[List[1, 0, 1], List[1, 2, 0], List[0, 1, 2]]",
        ),
        (
            "SylvesterMatrix[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]",
            "List[List[1, 0, 1, 1, 0], List[0, 1, 0, 1, 1], List[2, 0, 3, 0, 0], \
             List[0, 2, 0, 3, 0], List[0, 0, 2, 0, 3]]",
        ),
        (
            "Det[SylvesterMatrix[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]]",
            "Resultant[Plus[Power[x, 3], x, 1], Plus[Times[2, Power[x, 2]], 3], x]",
        ),
        // giac's [[5,12],[21,32]].
        (
            "HadamardProduct[List[List[1, 2], List[3, 4]], List[List[5, 6], List[7, 8]]]",
            "List[List[5, 12], List[21, 32]]",
        ),
        (
            "HadamardProduct[List[x, 2], List[y, 3]]",
            "List[Times[x, y], 6]",
        ),
    ]);
}

/// `L D Lᵀ`, positive definiteness and inertia (CAS-477).
#[test]
fn symmetric_matrices() {
    check(&[
        // giac's [0,1,2], [[1,0,0],[1/2,1,0],[1/4,1/4,1]], diag(4, 2, 37/8).
        (
            "LDLDecomposition[List[List[4, 2, 1], List[2, 3, 1], List[1, 1, 5]]]",
            "List[List[0, 1, 2], List[List[1, 0, 0], List[1/2, 1, 0], List[1/4, 1/4, 1]], \
             List[List[4, 0, 0], List[0, 2, 0], List[0, 0, 37/8]]]",
        ),
        // A zero pivot with an entry below: giac's 2 × 2 block, here the call.
        (
            "LDLDecomposition[List[List[0, 1], List[1, 0]]]",
            "LDLDecomposition[List[List[0, 1], List[1, 0]]]",
        ),
        // giac's true, and true for [[1,2],[2,1]], whose eigenvalues are 3
        // and −1: False here.
        ("PositiveDefiniteQ[List[List[2, 1], List[1, 2]]]", "True"),
        ("PositiveDefiniteQ[List[List[1, 2], List[2, 1]]]", "False"),
        ("PositiveDefiniteQ[List[List[1, 2], List[3, 4]]]", "False"),
        ("PositiveDefiniteQ[List[List[0, 0], List[0, 1]]]", "False"),
        (
            "PositiveDefiniteQ[List[List[a, 0], List[0, 1]]]",
            "PositiveDefiniteQ[List[List[a, 0], List[0, 1]]]",
        ),
        // giac's [1,1,0] twice; eigenvalues 2 and 0; the zero matrix.
        ("Inertia[List[List[1, 0], List[0, -1]]]", "List[1, 1, 0]"),
        ("Inertia[List[List[1, 2], List[2, 1]]]", "List[1, 1, 0]"),
        ("Inertia[List[List[1, 1], List[1, 1]]]", "List[1, 0, 1]"),
        ("Inertia[List[List[0, 0], List[0, 0]]]", "List[0, 0, 2]"),
        (
            "Inertia[List[List[2, -1, 0], List[-1, 2, -1], List[0, -1, 2]]]",
            "List[3, 0, 0]",
        ),
    ]);
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(crate::Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// `N` of both sides at 40 digits, entry by entry within `10⁻²⁰`: radicals
/// written in other forms.
fn check_values(cases: &[(&str, &str)]) {
    use crate::Symbol;
    let ctx = Context::new();
    let numeric = |e: &str| {
        ctx.eval(&Expr::call(Symbol::N, [read(e), Expr::from(40)]))
            .unwrap()
    };
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (found, want) = (numeric(input), numeric(expected));
        // Entry by entry, balls on both sides: no `N` of a difference near
        // 0, and no `Max` of gaps it cannot order.
        let (a, b) = (leaves(&found), leaves(&want));
        let near = a.len() == b.len()
            && a.iter().zip(&b).all(|(x, y)| {
                let gap = Expr::call(
                    Symbol::ABS,
                    [Expr::call(
                        Symbol::PLUS,
                        [x.clone(), Expr::from(-1) * y.clone()],
                    )],
                );
                let small = Expr::call(Symbol::LESS, [gap, read("Power[10, -20]")]);
                ctx.eval(&small).unwrap().as_symbol() == Some(Symbol::TRUE)
            });
        if !near {
            wrong.push(format!("{input}: {found:?}, expected {want:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Norms, condition numbers and Gram–Schmidt (CAS-478).
#[test]
fn norms() {
    let a = "List[List[1, -2], List[3, 4]]";
    check(&[
        // giac's 6, 7 and 6, 7; of a vector, 6 and 3.
        (&format!("MatrixNorm[{a}, 1]"), "6"),
        (&format!("MatrixNorm[{a}, Infinity]"), "7"),
        (&format!("MatrixNorm[{a}, frobenius]"), "Power[30, 1/2]"),
        (&format!("L1Norm[{a}]"), "6"),
        (&format!("LInfinityNorm[{a}]"), "7"),
        ("L1Norm[List[1, -2, 3]]", "6"),
        ("LInfinityNorm[List[1, -2, 3]]", "3"),
        // giac's 21.0 twice, exactly; a singular matrix.
        ("ConditionNumber[List[List[1, 2], List[3, 4]]]", "21"),
        (
            "ConditionNumber[List[List[1, 2], List[3, 4]], Infinity]",
            "21",
        ),
        ("ConditionNumber[List[List[1, 2], List[2, 4]]]", "Infinity"),
    ]);
    check_values(&[
        // giac's 5.11667273602 = √(15 + 5√5), and 14.9330343737.
        (
            &format!("MatrixNorm[{a}, 2]"),
            "Power[Plus[15, Times[5, Power[5, 1/2]]], 1/2]",
        ),
        (
            "ConditionNumber[List[List[1, 2], List[3, 4]], 2]",
            "Times[Power[Plus[15, Power[221, 1/2]], 1/2], Power[Plus[15, Times[-1, Power[221, 1/2]]], -1/2]]",
        ),
        // giac's [[1/√2,1/√2,0],[1/√6,−1/√6,2/√6]].
        (
            "GramSchmidt[List[List[1, 1, 0], List[1, 0, 1]]]",
            "List[List[Power[2, -1/2], Power[2, -1/2], 0], \
             List[Power[6, -1/2], Times[-1, Power[6, -1/2]], Times[2, Power[6, -1/2]]]]",
        ),
    ]);
    // A dependent row goes.
    check(&[("Length[GramSchmidt[List[List[1, 2], List[2, 4]]]]", "1")]);
}

/// Singular values and conjugate gradients (CAS-479).
#[test]
fn singular_values_and_conjugate_gradients() {
    check(&[
        // giac's [1/11,7/11], exactly; a 3 × 3 positive definite system.
        (
            "ConjugateGradient[List[List[4, 1], List[1, 3]], List[1, 2]]",
            "List[1/11, 7/11]",
        ),
        (
            "ConjugateGradient[List[List[2, -1, 0], List[-1, 2, -1], List[0, -1, 2]], List[1, 0, 1]]",
            "List[1, 1, 1]",
        ),
        // Not positive definite: pᵀ A p vanishes.
        (
            "ConjugateGradient[List[List[0, 1], List[1, 0]], List[1, 0]]",
            "ConjugateGradient[List[List[0, 1], List[1, 0]], List[1, 0]]",
        ),
        ("SingularValues[List[List[3, 0], List[0, 4]]]", "List[3, 4]"),
    ]);
    // giac's [0.365966190626,5.46498570422].
    check_values(&[(
        "SingularValues[List[List[1, 2], List[3, 4]]]",
        "List[Power[Plus[15, Times[-1, Power[221, 1/2]]], 1/2], Power[Plus[15, Power[221, 1/2]], 1/2]]",
    )]);
}
