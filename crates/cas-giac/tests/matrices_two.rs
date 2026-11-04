//! Matrix decompositions II against giac (CAS-DOC-47 §6, CAS-PLAN-41):
//! giac's exact answers read back are this CAS's; its decimals (`cond` of
//! a non-diagonal matrix, `svl`, the 2-norm) within a relative 10⁻¹¹;
//! `hessenberg` of a larger matrix by its identity `A = P H P⁻¹`, the form
//! not being unique. giac's `isposdef([[1,2],[2,1]])` (`true`, though an
//! eigenvalue is −1) and its `hadamard` of two vectors (the Hadamard bound
//! `|u| |v|`, not the entrywise product) are left out.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac's answers to `inputs`, read back, are this CAS's.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's sequence answers to `inputs`, bracketed to be read as lists,
/// are this CAS's.
fn check_sequences(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| {
            giac.answers(&format!("[{input}]"), input)
                .err()
                .map(|e| e.to_string())
        })
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// giac's decimals against this CAS's values, entry by entry, within a
/// relative 10⁻¹¹.
fn check_near(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let numeric = |e: Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e, Expr::from(30)]))
            .unwrap()
    };
    let mut wrong = Vec::new();
    for input in inputs {
        let (Ok(theirs), Ok(mine)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        let (t, o) = (
            leaves(&numeric(theirs.clone())),
            leaves(&numeric(mine.clone())),
        );
        let near = t.len() == o.len()
            && t.iter().zip(&o).all(|(a, b)| {
                let gap = Expr::call(
                    Symbol::ABS,
                    [Expr::call(
                        Symbol::PLUS,
                        [b.clone(), Expr::from(-1) * a.clone()],
                    )],
                );
                let scale = Expr::call(
                    Symbol::MAX,
                    [Expr::from(1), Expr::call(Symbol::ABS, [a.clone()])],
                );
                let bound = scale * Expr::from(10).pow(-11);
                ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
                    .ok()
                    .and_then(|e| e.as_symbol())
                    == Some(Symbol::TRUE)
            });
        if !near {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Hessenberg forms and constructions (CAS-476).
#[test]
fn forms_and_constructions() {
    check_sequences(&["hessenberg([[1,2,3],[4,5,6],[7,8,10]])"]);
    check(&[
        "companion(x^3+2*x+5,x)",
        "companion(2*x^2+3,x)",
        "sylvester(x^2+1,x+2,x)",
        "sylvester(x^3+x+1,2*x^2+3,x)",
        "hadamard([[1,2],[3,4]],[[5,6],[7,8]])",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    // A larger one: giac's and this CAS's P H P⁻¹ are A.
    let ctx = Context::new();
    let a = "[[2,1,0,3],[1,4,2,1],[5,0,1,2],[3,1,1,6]]";
    for d in [
        giac.value(&format!("[hessenberg({a})]")).unwrap(),
        ours(&format!("hessenberg({a})")).unwrap(),
    ] {
        let [p, h] = d.args() else {
            panic!("{d:?}");
        };
        let back = ctx
            .eval(&Expr::call(
                Symbol::TIMES,
                [
                    p.clone(),
                    h.clone(),
                    Expr::call(Symbol::INVERSE, [p.clone()]),
                ],
            ))
            .unwrap();
        assert_eq!(back, ours(a).unwrap(), "{d:?}");
    }
}

/// `L D Lᵀ`, positive definiteness, inertia (CAS-477).
#[test]
fn symmetric_matrices() {
    check_sequences(&["ldl([[4,2,1],[2,3,1],[1,1,5]])", "ldl([[2,1],[1,2]])"]);
    check(&[
        "isposdef([[2,1],[1,2]])",
        "isposdef([[2,-1,0],[-1,2,-1],[0,-1,2]])",
        "inertia([[1,0],[0,-1]])",
        "inertia([[2,-1,0],[-1,2,-1],[0,-1,2]])",
        "inertia([[0,1],[1,0]])",
    ]);
}

/// Norms, condition numbers, singular values, conjugate gradients
/// (CAS-478, CAS-479).
#[test]
fn norms_and_solutions() {
    check(&[
        "matrix_norm([[1,-2],[3,4]],1)",
        "matrix_norm([[1,-2],[3,4]],inf)",
        "l1norm([[1,-2],[3,4]])",
        "linfnorm([[1,-2],[3,4]])",
        "l1norm([1,-2,3])",
        "linfnorm([1,-2,3])",
        "cond([[1,2],[3,4]])",
        "cond([[1,2],[3,4]],inf)",
        "conjugate_gradient([[4,1],[1,3]],[1,2])",
        "conjugate_gradient([[2,-1,0],[-1,2,-1],[0,-1,2]],[1,0,1])",
    ]);
    check_near(&[
        "cond([[4,1],[1,3]])",
        "cond([[1,2],[3,4]],2)",
        "matrix_norm([[1,-2],[3,4]],2)",
        "svl([[1,2],[3,4]])",
        "svl([[3,0],[0,4]])",
        "gramschmidt([[1,1,0],[1,0,1]])",
    ]);
}
