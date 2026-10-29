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

/// `LatticeReduce[B]` is a basis of the same lattice (equal Hermite forms
/// of the rows that are not 0) that is LLL-reduced: `|μ_kj| ≤ 1/2` and
/// `B_k ≥ (3/4 − μ²_(k,k−1)) B_(k−1)`.
fn check_reduced(b: &str) -> Expr {
    let ctx = Context::new();
    let found = eval(&ctx, &format!("LatticeReduce[{b}]"));
    let reduced = super::integers(&found).unwrap_or_default();
    let input = super::integers(&read(b)).unwrap();
    let form = |m: &super::Matrix| {
        let mut h = m.clone();
        if !h.is_empty() {
            super::hermite(&mut h, &mut super::identity(m.len()));
        }
        h.retain(|r| r.iter().any(|x| !x.is_zero()));
        h
    };
    assert_eq!(
        form(&reduced),
        form(&input),
        "{b}: {found:?} is another lattice"
    );
    let (mu, norms) = super::gram_schmidt(&reduced);
    let half = cas_num::Rational::new(cas_num::Integer::ONE, cas_num::Integer::from(2));
    let delta = cas_num::Rational::new(cas_num::Integer::from(3), cas_num::Integer::from(4));
    for k in 1..reduced.len() {
        for (j, m) in mu[k].iter().enumerate().take(k) {
            assert!(m.abs() <= half, "{b}: {found:?}: μ_{k}{j}");
        }
        let bound = &(&delta - &(&mu[k][k - 1] * &mu[k][k - 1])) * &norms[k - 1];
        assert!(norms[k] >= bound, "{b}: {found:?}: Lovász at {k}");
    }
    found
}

#[test]
fn lattice_reduction() {
    let ctx = Context::new();
    assert_eq!(
        check_reduced("List[List[1, 0], List[1000, 1]]"),
        eval(&ctx, "List[List[1, 0], List[0, 1]]")
    );
    // giac's examples: lll([[1,1,1],[-1,0,2],[3,5,6]]) and the knapsack-like
    // basis, whose first reduced row is [1,-2,1,0] in giac.
    check_reduced("List[List[1, 1, 1], List[-1, 0, 2], List[3, 5, 6]]");
    let knapsack =
        check_reduced("List[List[1, 0, 0, 12345], List[0, 1, 0, 23456], List[0, 0, 1, 34567]]");
    // Its shortest vector has norm² 6: the first row's is within 2^2 of it.
    let first = super::integers(&knapsack).unwrap()[0].clone();
    let norm: i64 = first.iter().map(|x| x.to_i64().unwrap().pow(2)).sum();
    assert!(norm <= 24, "{knapsack:?}");
    // Dependent rows: the lattice of (1, 2); and none at all.
    assert_eq!(
        check_reduced("List[List[2, 4], List[1, 2]]"),
        eval(&ctx, "List[List[1, 2]]")
    );
    assert_eq!(
        eval(&ctx, "LatticeReduce[List[List[0, 0]]]"),
        eval(&ctx, "List[]")
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
