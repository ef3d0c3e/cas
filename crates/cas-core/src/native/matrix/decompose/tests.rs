//! Tests of the decompositions and functions (CAS-DOC-30 §2, §3): their
//! defining identities, giac's answers where the form is canonical, and
//! what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The matrix product: `Times` applied as given, its factors in their
/// order (`Dot` of two matrices is their inner product).
fn dot(a: &Expr, b: &Expr) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])
}

/// Whether two matrices are equal entry by entry once put together.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = eval(ctx, &(a.clone() - b.clone()));
    difference.args().iter().all(|row| {
        row.args().iter().all(|x| {
            eval(
                ctx,
                &Expr::call(Symbol::TOGETHER, [Expr::call(Symbol::EXPAND, [x.clone()])]),
            ) == Expr::from(0)
        })
    })
}

const MATRICES: [&str; 5] = [
    "List[List[1, 2], List[3, 4]]",
    "List[List[0, 1], List[1, 1]]",
    "List[List[2, 1, 0], List[1, 3, 1], List[0, 1, 4]]",
    "List[List[0, 0, 1], List[1, 2, 3], List[2, 1, 0]]",
    "List[List[a, 1], List[1, b]]",
];

/// `A[p] = L U`, `L` unit lower triangular; giac's `lu([[1,2],[3,4]])`.
#[test]
fn lu() {
    let ctx = Context::new();
    for m in MATRICES {
        let a = read(m);
        let found = eval(&ctx, &Expr::call(Symbol::LU_DECOMPOSITION, [a.clone()]));
        let [p, l, u] = found.args() else {
            panic!("{m}: {found:?}");
        };
        let permuted = Expr::call(
            Symbol::LIST,
            p.args()
                .iter()
                .map(|i| {
                    Expr::call(
                        Symbol::PART,
                        [a.clone(), Expr::from(i.as_i64().unwrap() + 1)],
                    )
                })
                .collect::<Vec<_>>(),
        );
        assert!(same(&ctx, &dot(l, u), &permuted), "{m}: {found:?}");
    }
    assert_eq!(
        eval(&ctx, &read("LUDecomposition[List[List[1, 2], List[3, 4]]]")),
        eval(
            &ctx,
            &read("List[List[0, 1], List[List[1, 0], List[3, 1]], List[List[1, 2], List[0, -2]]]")
        )
    );
    assert_eq!(
        eval(&ctx, &read("LUDecomposition[List[List[0, 1], List[1, 1]]]")).args()[0],
        eval(&ctx, &read("List[1, 0]"))
    );
}

/// `Q R = A`, `Qᴴ Q = I`; giac's `qr([[3,0],[4,5]])`.
#[test]
fn qr() {
    let ctx = Context::new();
    for m in [
        "List[List[3, 0], List[4, 5]]",
        "List[List[1, 1], List[1, 0], List[0, 1]]",
        "List[List[1, 2], List[3, 4]]",
    ] {
        let a = read(m);
        let found = eval(&ctx, &Expr::call(Symbol::QR_DECOMPOSITION, [a.clone()]));
        let [q, r] = found.args() else {
            panic!("{m}: {found:?}");
        };
        assert!(same(&ctx, &dot(q, r), &a), "{m}: {found:?}");
        let n = r.args().len();
        let qtq = dot(&Expr::call(Symbol::TRANSPOSE, [q.clone()]), q);
        let ident = Expr::call(
            Symbol::IDENTITY_MATRIX,
            [Expr::from(i64::try_from(n).unwrap())],
        );
        assert!(same(&ctx, &qtq, &ident), "{m}: {found:?}");
    }
    assert_eq!(
        eval(&ctx, &read("QRDecomposition[List[List[3, 0], List[4, 5]]]")),
        eval(
            &ctx,
            &read("List[List[List[3/5, -4/5], List[4/5, 3/5]], List[List[5, 4], List[0, 3]]]")
        )
    );
}

/// `L Lᵀ = A`; giac's `cholesky([[4,2],[2,3]])` and `[[2,1],[1,2]]`.
#[test]
fn cholesky() {
    let ctx = Context::new();
    for (m, expected) in [
        (
            "List[List[4, 2], List[2, 3]]",
            "List[List[2, 0], List[1, Power[2, 1/2]]]",
        ),
        (
            "List[List[2, 1], List[1, 2]]",
            "List[List[Power[2, 1/2], 0], List[Times[1/2, Power[2, 1/2]], Times[1/2, Power[6, 1/2]]]]",
        ),
    ] {
        let found = eval(&ctx, &Expr::call(Symbol::CHOLESKY_DECOMPOSITION, [read(m)]));
        assert!(same(&ctx, &found, &read(expected)), "{m}: {found:?}");
        let back = dot(&found, &Expr::call(Symbol::TRANSPOSE, [found.clone()]));
        assert!(same(&ctx, &back, &read(m)), "{m}: {found:?}");
    }
}

/// `A = P J P⁻¹` with `J` in Jordan form; giac's `jordan([[1,1],[0,1]])`
/// and the nilpotent `3 × 3` block; the minimal polynomials.
#[test]
fn jordan_and_minimal_polynomial() {
    let ctx = Context::new();
    for m in [
        "List[List[1, 1], List[0, 1]]",
        "List[List[2, 1], List[1, 2]]",
        "List[List[0, 1, 0], List[0, 0, 1], List[0, 0, 0]]",
        "List[List[1, 2], List[3, 4]]",
        "List[List[2, 1, 0], List[0, 2, 0], List[0, 0, 3]]",
        "List[List[3, 1, 1], List[0, 3, 0], List[0, 0, 3]]",
    ] {
        let a = read(m);
        let found = eval(&ctx, &Expr::call(Symbol::JORDAN_DECOMPOSITION, [a.clone()]));
        let [p, j] = found.args() else {
            panic!("{m}: {found:?}");
        };
        let back = dot(&dot(p, j), &Expr::call(Symbol::INVERSE, [p.clone()]));
        assert!(same(&ctx, &back, &a), "{m}: {found:?}");
    }
    assert_eq!(
        eval(
            &ctx,
            &read("JordanDecomposition[List[List[1, 1], List[0, 1]]]")
        )
        .args()[1],
        eval(&ctx, &read("List[List[1, 1], List[0, 1]]"))
    );
    for (m, expected) in [
        (
            "MinimalPolynomial[List[List[1, 1], List[0, 1]], x]",
            "Plus[1, Times[-2, x], Power[x, 2]]",
        ),
        (
            "MinimalPolynomial[List[List[2, 0], List[0, 2]]]",
            "List[1, -2]",
        ),
        (
            "MinimalPolynomial[List[List[3, 1, 1], List[0, 3, 0], List[0, 0, 3]], x]",
            "Plus[9, Times[-6, x], Power[x, 2]]",
        ),
    ] {
        assert_eq!(eval(&ctx, &read(m)), eval(&ctx, &read(expected)), "{m}");
    }
}

/// giac's `exp` of matrices and `matpow`, symbolic `n` included.
#[test]
fn functions() {
    let ctx = Context::new();
    for (m, expected) in [
        (
            "MatrixExp[List[List[0, 1], List[0, 0]]]",
            "List[List[1, 1], List[0, 1]]",
        ),
        (
            "MatrixExp[List[List[1, 0], List[0, 2]]]",
            "List[List[Exp[1], 0], List[0, Exp[2]]]",
        ),
        (
            "MatrixExp[List[List[0, 1, 0], List[0, 0, 1], List[0, 0, 0]]]",
            "List[List[1, 1, 1/2], List[0, 1, 1], List[0, 0, 1]]",
        ),
        (
            "MatrixPower[List[List[1, 1], List[0, 1]], n]",
            "List[List[1, n], List[0, 1]]",
        ),
        (
            "MatrixPower[List[List[2, 1], List[0, 2]], n]",
            "List[List[Power[2, n], Times[n, Power[2, Plus[-1, n]]]], List[0, Power[2, n]]]",
        ),
        (
            "MatrixPower[List[List[1, 2], List[3, 4]], 2]",
            "List[List[7, 10], List[15, 22]]",
        ),
    ] {
        let found = eval(&ctx, &read(m));
        assert!(same(&ctx, &found, &read(expected)), "{m}: {found:?}");
    }
    // A^n of [[1,2],[3,4]] at n = 0 … 4 is the products.
    let symbolic = eval(&ctx, &read("MatrixPower[List[List[1, 2], List[3, 4]], n]"));
    for k in 0..=4i64 {
        let at = symbolic.substitute(&[(Expr::symbol("n"), Expr::from(k))]);
        let power = Expr::call(
            Symbol::POWER,
            [read("List[List[1, 2], List[3, 4]]"), Expr::from(k)],
        );
        assert!(same(&ctx, &at, &power), "n = {k}: {symbolic:?}");
    }
}

/// What stays: no matrix, not square, a dependent column for QR, a
/// matrix that is not symmetric or not positive for Cholesky.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "LUDecomposition[x]",
        "LUDecomposition[List[List[1, 2, 3], List[4, 5, 6]]]",
        "QRDecomposition[List[List[1, 2], List[2, 4]]]",
        "CholeskyDecomposition[List[List[1, 2], List[3, 4]]]",
        "CholeskyDecomposition[List[List[1, 2], List[2, 1]]]",
        "JordanDecomposition[List[List[1, 2, 3]]]",
    ] {
        let e = read(input);
        assert_eq!(eval(&ctx, &e).head(), e.head(), "{input}");
    }
}
