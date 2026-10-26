//! Linear algebra against giac (CAS-DOC-12, CAS-PLAN-6): giac's answers,
//! its matrices `matrix[…]` read back as lists, are this CAS's.

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

/// The parameters' values where entries with parameters are compared.
const PARAMETERS: [(&str, i64, i64); 5] = [
    ("a", 3, 7),
    ("b", 5, 11),
    ("c", -2, 13),
    ("d", 7, 5),
    ("x", 2, 9),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists part by part, entries the same, equal
/// by `==`, or equal by `==` with the parameters at their values.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let equal = |a: &Expr, b: &Expr| holds(ctx, &Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]));
    if a == b || equal(a, b) {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = PARAMETERS
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let at = |e: &Expr| ctx.eval(&e.substitute(&rules)).ok();
    match (at(a), at(b)) {
        (Some(a), Some(b)) => a == b || equal(&a, &b),
        _ => false,
    }
}

/// giac's answers and this CAS's, the same by [`same`].
fn check_values(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
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

/// Determinants and inverses (CAS-219).
#[test]
fn determinants_and_inverses() {
    check_values(&[
        "det([[1,2],[3,4]])",
        "det([[a,b],[c,d]])",
        "det([[1,2,3],[4,5,6],[7,8,10]])",
        "det([[x,1],[1,x]])",
        "det([[2]])",
        "inv([[1,2],[3,4]])",
        "inv([[a,b],[c,d]])",
        "inv([[1,2,3],[0,1,4],[5,6,0]])",
        "inv([[2]])",
        "[[1,2],[3,4]]^-1",
        "[[1,2],[3,4]]^-2",
    ]);
}
