//! giac's minors, echelon forms and pivots (CAS-530, CAS-528) against
//! giac 2.0: its answers read back are this CAS's, and the departures,
//! each with giac's answer as it prints it and this CAS's as decided.

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

#[test]
fn minors() {
    check(&[
        "det_minor([[1,2],[3,4]])",
        "det_minor([[1,2,3],[4,5,6],[7,8,10]])",
        "det_minor([[a,b],[c,d]])",
        "adjoint_matrix([[1,2],[3,4]])",
        "adjoint_matrix([[a,b],[c,d]])",
    ]);
}

#[test]
fn echelon_forms() {
    check(&[
        "ref([[1,2],[3,4]])",
        "ref([[2,1],[1,3]])",
        "ref([[3,1,-2],[3,2,2]])",
        "ref([[2,1,1,-1],[1,1,2,-1],[1,2,1,-4]])",
        "ref([[0,0],[0,3]])",
        "ref([[0,2],[0,3],[1,1]])",
        "ref([[a,b],[c,d]])",
        "rref([[1,2],[3,4]],keep_pivot)",
        "rref([[4,2],[2,4]],keep_pivot)",
        "rref([[0,2],[3,4]],keep_pivot)",
        "rref([[2,1],[1,3]],keep_pivot)",
        "rref([[2,4,6],[1,3,5]],keep_pivot)",
        "rref([[1,2,3],[4,5,6],[7,8,10]],keep_pivot)",
        "rref([[a,b],[c,d]],keep_pivot)",
        "rref([[2.0,1],[1,3]],keep_pivot)",
    ]);
}

#[test]
fn pivots() {
    check(&[
        "pivot([[1,2],[3,4]],1,1)",
        "pivot([[1,2,3],[4,5,6],[7,8,10]],1,2)",
        "pivot([[1,2],[3,4],[5,6]],0,1)",
        "egvl([[4,1],[-4,0]])",
        "eigVl([[4,1],[-4,0]])",
    ]);
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // The eigenvalues in this CAS's order, Eigenvalues', where giac's
        // are in its own.
        (
            "egvl([[1,2],[3,4]])",
            "[[(sqrt(33)+5)/2,0],[0,(-sqrt(33)+5)/2]]",
            "[[(-sqrt(33)+5)/2,0],[0,(sqrt(33)+5)/2]]",
        ),
        ("egvl([[1,1],[0,2]])", "[[2,0],[0,1]]", "[[1,0],[0,2]]"),
        // CAS-527: giac's inert capitalised forms stay unevaluated user
        // calls here; an alias would evaluate.
        (
            "Det([[1,2],[3,4]])",
            "det([[1,2],[3,4]])",
            "Det([[1,2],[3,4]])",
        ),
        (
            "Inverse([[1,2],[3,4]])",
            "inverse([[1,2],[3,4]])",
            "Inverse([[1,2],[3,4]])",
        ),
        (
            "Nullspace([[1,2],[2,4]])",
            "nullspace([[1,2],[2,4]])",
            "Nullspace([[1,2],[2,4]])",
        ),
        (
            "Rref([[1,2],[3,4]])",
            "rref([[1,2],[3,4]])",
            "Rref([[1,2],[3,4]])",
        ),
    ];
    let normal = |s: &str| s.replace("matrix[", "[").replace(' ', "");
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if normal(&printed) != normal(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
