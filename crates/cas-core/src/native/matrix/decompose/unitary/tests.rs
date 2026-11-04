//! Tests of the unitary decompositions (CAS-DOC-51): `Q` unitary, `T` upper
//! triangular and `A = Q T Qᴴ`, exactly where the entries are rational and
//! at 40 digits where they have radicals.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The matrix product, its factors in their order.
fn dot(a: &Expr, b: &Expr) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])
}

/// `Mᴴ`.
fn adjoint(m: &Expr) -> Expr {
    Expr::call(
        Symbol::CONJUGATE,
        [Expr::call(Symbol::TRANSPOSE, [m.clone()])],
    )
}

/// Whether every entry of the matrix `m` is 0 within `10⁻²⁰` at 40 digits.
fn negligible(ctx: &Context, m: &Expr) -> bool {
    let value = eval(ctx, &Expr::call(Symbol::N, [m.clone(), Expr::from(40)]));
    value.args().iter().all(|row| {
        row.args().iter().all(|x| {
            let small = Expr::call(
                Symbol::LESS,
                [Expr::call(Symbol::ABS, [x.clone()]), read("Power[10, -20]")],
            );
            eval(ctx, &small).as_symbol() == Some(Symbol::TRUE)
        })
    })
}

/// The Schur form of `a`, checked: `Q` unitary, `T` upper triangular with
/// `diagonal` on its diagonal when given, `Q T Qᴴ = A`.
fn schur(ctx: &Context, a: &str, diagonal: Option<&[&str]>) -> (Expr, Expr) {
    let a = read(a);
    let found = eval(ctx, &Expr::call(Symbol::SCHUR_DECOMPOSITION, [a.clone()]));
    let [q, t] = found.args() else {
        panic!("{a:?}: {found:?}");
    };
    let n = t.args().len();
    for (i, row) in t.args().iter().enumerate() {
        for (k, x) in row.args().iter().enumerate() {
            assert!(k >= i || *x == Expr::from(0), "{a:?}: T = {t:?}");
        }
        if let Some(d) = diagonal {
            let want = eval(ctx, &read(d[i]));
            let gap = eval(ctx, &(row.args()[i].clone() - want));
            assert!(negligible(ctx, &read_matrix(gap)), "{a:?}: T = {t:?}");
        }
    }
    let ident = eval(
        ctx,
        &Expr::call(
            Symbol::IDENTITY_MATRIX,
            [Expr::from(i64::try_from(n).expect("small"))],
        ),
    );
    let unitary = eval(ctx, &(dot(&adjoint(q), q) - ident));
    assert!(negligible(ctx, &unitary), "{a:?}: Q = {q:?}");
    let back = eval(ctx, &(dot(&dot(q, t), &adjoint(q)) - a.clone()));
    assert!(negligible(ctx, &back), "{a:?}: Q = {q:?}, T = {t:?}");
    (q.clone(), t.clone())
}

/// The 1 × 1 matrix of `x`.
fn read_matrix(x: Expr) -> Expr {
    Expr::call(Symbol::LIST, [Expr::call(Symbol::LIST, [x])])
}

#[test]
fn schur_forms() {
    let ctx = Context::new();
    // Rational eigenvalues: exact, giac's T has 3 then 1 and the corner −2.
    let (q, t) = schur(&ctx, "List[List[1, 2], List[0, 3]]", None);
    let product = eval(&ctx, &dot(&dot(&q, &t), &adjoint(&q)));
    assert_eq!(
        eval(&ctx, &Expr::call(Symbol::TOGETHER, [product])),
        read("List[List[1, 2], List[0, 3]]")
    );
    // Radicals: (5 ± √33)/2.
    schur(&ctx, "List[List[1, 2], List[3, 4]]", None);
    // Symmetric: T diagonal, 3 and −1 in some order.
    let (_, t) = schur(&ctx, "List[List[1, 2], List[2, 1]]", None);
    let corner = t.args()[0].args()[1].clone();
    assert!(negligible(&ctx, &read_matrix(corner)), "T = {t:?}");
    // Complex eigenvalues ±i: the complex form.
    schur(&ctx, "List[List[0, -1], List[1, 0]]", None);
    // A Jordan block: not diagonalizable.
    schur(&ctx, "List[List[2, 1], List[0, 2]]", Some(&["2", "2"]));
    schur(
        &ctx,
        "List[List[2, 1, 0], List[-1, 4, 0], List[1, 1, 5]]",
        None,
    );
    // Three by three with distinct rational eigenvalues.
    schur(
        &ctx,
        "List[List[2, 0, 0], List[1, 3, 0], List[4, 5, 6]]",
        None,
    );
}

#[test]
fn schur_stays() {
    let ctx = Context::new();
    for input in [
        "SchurDecomposition[List[List[1, 2, 3], List[4, 5, 6]]]",
        "SchurDecomposition[x]",
        "SchurDecomposition[List[List[1, 2], List[3, 4]], 1]",
    ] {
        let e = read(input);
        assert_eq!(eval(&ctx, &e), e, "{input}");
    }
}
