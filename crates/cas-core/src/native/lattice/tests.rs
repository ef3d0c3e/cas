//! Tests of the integer matrices (CAS-445, CAS-446): giac's normal forms,
//! the identities `U A = H` and `U A V = S` with unimodular `U` and `V`,
//! singular and rectangular matrices, and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// Whether `det m` is ±1.
fn unimodular(ctx: &Context, m: &Expr) -> bool {
    let det = ctx.eval(&Expr::call(Symbol::DET, [m.clone()])).unwrap();
    det == Expr::from(1) || det == Expr::from(-1)
}

/// `m₁ m₂` (`Dot`).
fn dot(ctx: &Context, a: &Expr, b: &Expr) -> Expr {
    ctx.eval(&Expr::call(Symbol::DOT, [a.clone(), b.clone()]))
        .unwrap()
}

/// `HermiteDecomposition[A]` is `{U, H}` with `U A = H`, `U` unimodular,
/// `H` the expected.
fn check_hermite(a: &str, h: &str) {
    let ctx = Context::new();
    let found = eval(&ctx, &format!("HermiteDecomposition[{a}]"));
    let [u, found_h] = found.args() else {
        panic!("{a}: {found:?}");
    };
    assert_eq!(*found_h, eval(&ctx, h), "{a}: H");
    assert_eq!(dot(&ctx, u, &read(a)), *found_h, "{a}: U A");
    assert!(unimodular(&ctx, u), "{a}: U {u:?}");
}

/// `SmithDecomposition[A]` is `{U, S, V}` with `U A V = S`, `U` and `V`
/// unimodular, `S` the expected.
fn check_smith(a: &str, s: &str) {
    let ctx = Context::new();
    let found = eval(&ctx, &format!("SmithDecomposition[{a}]"));
    let [u, found_s, v] = found.args() else {
        panic!("{a}: {found:?}");
    };
    assert_eq!(*found_s, eval(&ctx, s), "{a}: S");
    let product = dot(&ctx, &dot(&ctx, u, &read(a)), v);
    assert_eq!(product, *found_s, "{a}: U A V");
    assert!(
        unimodular(&ctx, u) && unimodular(&ctx, v),
        "{a}: {u:?} {v:?}"
    );
}

#[test]
fn hermite_forms() {
    // giac's [[3,0,30],[0,12,0],[0,0,60]] and [[1,1],[0,2]].
    check_hermite(
        "List[List[9, -36, 30], List[-36, 192, -180], List[30, -180, 180]]",
        "List[List[3, 0, 30], List[0, 12, 0], List[0, 0, 60]]",
    );
    check_hermite(
        "List[List[2, 4], List[1, 3]]",
        "List[List[1, 1], List[0, 2]]",
    );
    // Singular, rectangular, zero rows last.
    check_hermite(
        "List[List[2, 4, 6], List[1, 2, 3]]",
        "List[List[1, 2, 3], List[0, 0, 0]]",
    );
    check_hermite(
        "List[List[0, 3], List[0, 5], List[4, 1]]",
        "List[List[4, 0], List[0, 1], List[0, 0]]",
    );
    check_hermite("List[List[-6]]", "List[List[6]]");
}

#[test]
fn smith_forms() {
    // giac's [[3,0,0],[0,12,0],[0,0,60]] and [[2,0],[0,4]].
    check_smith(
        "List[List[9, -36, 30], List[-36, 192, -180], List[30, -180, 180]]",
        "List[List[3, 0, 0], List[0, 12, 0], List[0, 0, 60]]",
    );
    check_smith(
        "List[List[2, 4], List[6, 8]]",
        "List[List[2, 0], List[0, 4]]",
    );
    // 2 ∤ 3: the corner becomes gcd(2, 3) = 1, and 6 the other.
    check_smith(
        "List[List[2, 0], List[0, 3]]",
        "List[List[1, 0], List[0, 6]]",
    );
    check_smith(
        "List[List[2, 4, 6], List[1, 2, 3]]",
        "List[List[1, 0, 0], List[0, 0, 0]]",
    );
    check_smith(
        "List[List[0, 0], List[0, 0]]",
        "List[List[0, 0], List[0, 0]]",
    );
}

#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "HermiteDecomposition[List[List[1, 1/2]]]",
        "HermiteDecomposition[List[List[1, 2], List[3]]]",
        "SmithDecomposition[x]",
        "SmithDecomposition[List[]]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
