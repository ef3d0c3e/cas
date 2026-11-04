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
        let jordan = eval(ctx, &Expr::call(Symbol::JORDAN_DECOMPOSITION, [a.clone()]));
        let qr = eval(
            ctx,
            &Expr::call(
                Symbol::QR_DECOMPOSITION,
                [jordan
                    .args()
                    .first()
                    .cloned()
                    .unwrap_or_else(|| Expr::from(0))],
            ),
        );
        panic!("{a:?}: {found:?}; Jordan {jordan:?}; its P's QR {qr:?}");
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

/// The SVD of `a`, checked: `U` and `V` unitary, `s` descending and equal
/// to `SingularValues` reversed, `U diag(s) Vᴴ = A`.
fn svd(ctx: &Context, a: &str) -> Expr {
    let a = read(a);
    let found = eval(
        ctx,
        &Expr::call(Symbol::SINGULAR_VALUE_DECOMPOSITION, [a.clone()]),
    );
    let [u, s, v] = found.args() else {
        let gram = eval(ctx, &dot(&adjoint(&a), &a));
        let vectors = eval(ctx, &Expr::call(Symbol::EIGENVECTORS, [gram.clone()]));
        let values = eval(ctx, &Expr::call(Symbol::EIGENVALUES, [gram.clone()]));
        let basis = eval(
            ctx,
            &Expr::call(
                Symbol::GRAM_SCHMIDT,
                [Expr::call(Symbol::TRANSPOSE, [vectors.clone()])],
            ),
        );
        panic!(
            "{a:?}: {found:?}; AᴴA {gram:?}, values {values:?}, vectors {vectors:?}, GramSchmidt {basis:?}"
        );
    };
    let n = s.args().len();
    let ident = eval(
        ctx,
        &Expr::call(
            Symbol::IDENTITY_MATRIX,
            [Expr::from(i64::try_from(n).expect("small"))],
        ),
    );
    for w in [u, v] {
        let unitary = eval(ctx, &(dot(&adjoint(w), w) - ident.clone()));
        assert!(negligible(ctx, &unitary), "{a:?}: {found:?}");
    }
    let ascending = eval(ctx, &Expr::call(Symbol::SINGULAR_VALUES, [a.clone()]));
    let reversed = eval(ctx, &Expr::call(Symbol::REVERSE, [ascending]));
    let gap = eval(ctx, &(s.clone() - reversed));
    assert!(
        negligible(ctx, &Expr::call(Symbol::LIST, [gap])),
        "{a:?}: s = {s:?}"
    );
    let diagonal = Expr::call(Symbol::DIAGONAL_MATRIX, [s.clone()]);
    let back = eval(ctx, &(dot(&dot(u, &diagonal), &adjoint(v)) - a.clone()));
    assert!(negligible(ctx, &back), "{a:?}: {found:?}");
    s.clone()
}

#[test]
fn singular_value_decompositions() {
    let ctx = Context::new();
    // giac's s = [3√5, √5].
    let s = svd(&ctx, "List[List[3, 0], List[4, 5]]");
    let gap = eval(
        &ctx,
        &(s.clone() - read("List[Times[3, Sqrt[5]], Sqrt[5]]")),
    );
    assert!(negligible(&ctx, &Expr::call(Symbol::LIST, [gap])), "{s:?}");
    // √(15 ± √221), giac's [5.46498570422,0.365966190626].
    svd(&ctx, "List[List[1, 2], List[3, 4]]");
    // Singular: a zero singular value, U completed.
    let s = svd(&ctx, "List[List[1, 2], List[2, 4]]");
    assert_eq!(s, eval(&ctx, &read("List[5, 0]")), "{s:?}");
    // A repeated singular value: an orthogonal matrix.
    svd(&ctx, "List[List[0, 1], List[1, 0]]");
    svd(&ctx, "List[List[2, 0, 0], List[0, 0, 3], List[0, 0, 0]]");
    // Complex entries.
    svd(&ctx, "List[List[1, I], List[0, 1]]");
}

/// Whether every entry of the matrix `m` is a ball holding 0, or 0.
fn holds_zero(ctx: &Context, m: &Expr) -> bool {
    use crate::expr::View;
    let m = eval(ctx, m);
    m.args().iter().all(|row| {
        row.args().iter().all(|x| match x.view() {
            View::Ball(b) => b.contains_rational(&cas_num::Rational::ZERO),
            View::ComplexBall(z) => z.may_be_zero(),
            _ => *x == Expr::from(0),
        })
    })
}

/// Balls (CAS-498): the identities hold as balls holding 0.
#[test]
fn ball_matrices() {
    let ctx = Context::new();
    for a in [
        "N[List[List[1, 2], List[3, 4]], 30]",
        "N[List[List[2, 1, 0], List[1, 3, 1], List[0, 1, 4]], 30]",
        "N[List[List[0, -1], List[1, 0]], 30]",
    ] {
        let a = eval(&ctx, &read(a));
        let found = eval(&ctx, &Expr::call(Symbol::SCHUR_DECOMPOSITION, [a.clone()]));
        let [q, t] = found.args() else {
            panic!("{a:?}: {found:?}");
        };
        let back = dot(&dot(q, t), &adjoint(q)) - a.clone();
        assert!(holds_zero(&ctx, &back), "{a:?}: {found:?}");
        let n = i64::try_from(q.args().len()).expect("small");
        let ident = Expr::call(Symbol::IDENTITY_MATRIX, [Expr::from(n)]);
        assert!(
            holds_zero(&ctx, &(dot(&adjoint(q), q) - ident.clone())),
            "{a:?}: {found:?}"
        );
        let found = eval(
            &ctx,
            &Expr::call(Symbol::SINGULAR_VALUE_DECOMPOSITION, [a.clone()]),
        );
        let [u, s, v] = found.args() else {
            panic!("{a:?}: {found:?}");
        };
        let diagonal = Expr::call(Symbol::DIAGONAL_MATRIX, [s.clone()]);
        let back = dot(&dot(u, &diagonal), &adjoint(v)) - a.clone();
        assert!(holds_zero(&ctx, &back), "{a:?}: {found:?}");
        // Eigenvectors: A P − P diag(λ) holds 0.
        let p = eval(&ctx, &Expr::call(Symbol::EIGENVECTORS, [a.clone()]));
        let values = Expr::call(Symbol::EIGENVALUES, [a.clone()]);
        let lambda = Expr::call(Symbol::DIAGONAL_MATRIX, [values]);
        let gap = dot(&a, &p) - dot(&p, &lambda);
        assert!(holds_zero(&ctx, &gap), "{a:?}: {p:?}");
    }
    // A double eigenvalue: not certified, the call stays.
    let a = eval(&ctx, &read("N[List[List[1, 0], List[0, 1]], 30]"));
    let found = eval(&ctx, &Expr::call(Symbol::EIGENVECTORS, [a]));
    assert_eq!(
        found.head(),
        Some(&Expr::from(Symbol::EIGENVECTORS)),
        "{found:?}"
    );
}
