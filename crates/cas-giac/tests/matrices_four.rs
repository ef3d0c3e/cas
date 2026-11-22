//! giac's minors, echelon forms and pivots (CAS-530, CAS-528), `gauss`
//! and `rat_jordan` (CAS-535) against giac 2.0: its answers read back are
//! this CAS's, and the departures, each with giac's answer as it prints it
//! and this CAS's as decided. `rat_jordan` returns a sequence, read in
//! brackets as a list.

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

#[test]
fn jordan_matrices() {
    check(&[
        "egvl([[1,2],[3,4]])",
        "egvl([[1,1],[0,2]])",
        "egvl([[1,0,1],[0,2,0],[0,0,3]])",
        "egvl([[1,2],[2,1]])",
        "egvl([[0,1],[-2,2]])",
        "egvl([[0,-1],[1,0]])",
        "egvl([[2,0,0],[0,1,0],[0,0,3]])",
        "egvl([[1,1,0],[0,1,0],[0,0,-2]])",
    ]);
}

#[test]
fn squares() {
    check(&[
        "gauss(x^2+2*x*y,[x,y])",
        "gauss(x^2+y^2+z^2+2*x*y+4*y*z,[x,y,z])",
        "gauss(x*y,[x,y])",
        "gauss(2*x^2-3*y^2,[x,y])",
        "gauss(2*x*y,[x,y])",
        "gauss(x*y+y*z,[x,y,z])",
        "gauss(x^2+x*y+y^2+z^2,[x,y,z])",
        "gauss(2*x^2+2*x*y+3*y^2,[x,y])",
        "gauss(x^2+2*a*x*y,[x,y])",
        "gauss(y*z+x*z,[x,y,z])",
    ]);
}

/// giac's `[rat_jordan(A)]`, read back, is this CAS's `rat_jordan(A)`.
#[test]
fn rational_jordan_forms() {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = [
        "rat_jordan([[1,1],[0,1]])",
        "rat_jordan([[1,2],[3,4]])",
        "rat_jordan([[0,-1],[1,0]])",
        "rat_jordan([[2,0],[0,3]])",
        "rat_jordan([[1,2,0],[3,4,0],[0,0,7]])",
    ]
    .iter()
    .filter_map(|input| {
        giac.answers(&format!("[{input}]"), input)
            .err()
            .map(|e| e.to_string())
    })
    .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac's P from its own kernels: the vectors of an eigenvalue
        // here are those of the reduced kernel, 1 at each free column.
        (
            "rat_jordan([[2,1],[0,3]])",
            "[[1,-1],[1,0]],[[3,0],[0,2]]",
            "[[[1,1],[1,0]],[[3,0],[0,2]]]",
        ),
        (
            "rat_jordan([[3,1],[0,2]])",
            "[[1,1],[0,-1]],[[3,0],[0,2]]",
            "[[[1,-1],[0,1]],[[3,0],[0,2]]]",
        ),
        (
            "rat_jordan([[-2,-2,1],[-2,1,-2],[1,-2,-2]])",
            "[[1,5,5],[-2,0,5],[1,-5,5]],[[3,0,0],[0,-3,0],[0,0,-3]]",
            "[[[1,2,-1],[-2,1,0],[1,0,1]],[[3,0,0],[0,-3,0],[0,0,-3]]]",
        ),
        // Already a rational Jordan form: P = I here, giac's P is not.
        (
            "rat_jordan([[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]])",
            "[[4,0,0,4],[0,4,-4,0],[0,0,4,0],[0,0,0,4]],[[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]]",
            "[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],[[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]]]",
        ),
        // Not a quadratic form: giac returns it, the call stays here.
        ("gauss(x^3+y^2,[x,y])", "x^3+y^2", "gauss(x^3+y^2,[x,y])"),
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
