//! Linear algebra against giac (CAS-DOC-12, CAS-PLAN-6): giac's answers,
//! its matrices `matrix[…]` read back as lists, are this CAS's.

use cas_giac::Giac;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

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

/// The basics of CAS-218.
#[test]
fn basics() {
    check(&[
        "tran([[1,2,3],[4,5,6]])",
        "transpose([[1,2],[3,4]])",
        "tran([1,2,3])",
        "idn(3)",
        "identity(2)",
        "trace([[1,2],[3,4]])",
        "trace([[1,2,3],[4,5,6],[7,8,9]])",
        "dot([1,2,3],[4,5,6])",
        "dot([i,0],[1,0])",
        "dot([1,0],[i,0])",
        "dot([1+i,2],[3,4*i])",
        "dot([[1,2],[3,4]],[[1,0],[0,1]])",
        "cross([1,0,0],[0,1,0])",
        "cross([1,2],[3,4])",
        "cross([1,2,3],[4,5,6])",
    ]);
}
