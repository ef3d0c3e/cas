//! Tests of `LeastSquares`, `LQDecomposition` and
//! `AscendingSingularValueDecomposition` (CAS-623): giac's exact answers,
//! the defining products, and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The matrix product: `Times` applied as given, its factors in order.
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

/// giac's answers (probe2, probe6), exact.
#[test]
fn least_squares() {
    let ctx = Context::new();
    let cases = [
        (
            "LeastSquares[List[List[1, 2], List[3, 4]], List[5, 11]]",
            "List[List[1], List[2]]",
        ),
        (
            "LeastSquares[List[List[1, 2], List[3, 4], List[5, 6]], List[1, 2, 4]]",
            "List[List[2/3], List[1/12]]",
        ),
        (
            "LeastSquares[List[List[1, 2], List[3, 4], List[5, 6]], List[List[1], List[2], List[4]]]",
            "List[List[2/3], List[1/12]]",
        ),
        (
            "LeastSquares[List[List[1, 2], List[3, 4]], List[List[5, -1], List[11, -1]]]",
            "List[List[1, 1], List[2, -1]]",
        ),
    ];
    for (input, expected) in cases {
        assert_eq!(
            eval(&ctx, &read(input)),
            eval(&ctx, &read(expected)),
            "{input}"
        );
    }
}

/// `A = L Q`, `L` lower triangular, `P` the identity of `A`'s columns.
#[test]
fn lq() {
    let ctx = Context::new();
    for a in [
        "List[List[1, 2], List[3, 4]]",
        "List[List[6, 12, 18], List[5, 14, 31], List[3, 8, 18]]",
        "List[List[4, 0, 0], List[8, -4, 3]]",
    ] {
        let a = read(a);
        let n = a.args()[0].args().len();
        let result = eval(&ctx, &Expr::call(Symbol::LQ_DECOMPOSITION, [a.clone()]));
        let [l, q, p] = result.args() else {
            panic!("{result:?}");
        };
        assert!(same(&ctx, &eval(&ctx, &dot(l, q)), &a), "{a:?}");
        for (i, row) in l.args().iter().enumerate() {
            for x in &row.args()[i + 1..] {
                assert_eq!(*x, Expr::from(0), "{a:?}");
            }
        }
        let identity = eval(
            &ctx,
            &Expr::call(
                Symbol::IDENTITY_MATRIX,
                [Expr::from(i64::try_from(n).unwrap())],
            ),
        );
        assert_eq!(*p, identity, "{a:?}");
    }
}

/// `A = U diag(s) Vᴴ` with `s` ascending (giac's `SVD`).
#[test]
fn ascending_svd() {
    let ctx = Context::new();
    for a in [
        "List[List[1, 4], List[4, 1]]",
        "List[List[3, 0], List[4, 5]]",
    ] {
        let a = read(a);
        let result = eval(
            &ctx,
            &Expr::call(Symbol::ASCENDING_SINGULAR_VALUE_DECOMPOSITION, [a.clone()]),
        );
        let [u, s, v] = result.args() else {
            panic!("{result:?}");
        };
        let ascending = eval(&ctx, &Expr::call(Symbol::SINGULAR_VALUES, [a.clone()]));
        let row = |x: &Expr| Expr::call(Symbol::LIST, [x.clone()]);
        assert!(same(&ctx, &row(s), &row(&ascending)), "{a:?}");
        let sigma = Expr::call(Symbol::DIAGONAL_MATRIX, [s.clone()]);
        let vh = Expr::call(
            Symbol::CONJUGATE,
            [Expr::call(Symbol::TRANSPOSE, [v.clone()])],
        );
        let product = eval(&ctx, &dot(&dot(u, &sigma), &vh));
        assert!(same(&ctx, &product, &a), "{a:?}");
    }
}

/// What stays: shapes that do not fit, an under-determined or singular
/// system, `LQ` of a tall or rank-deficient matrix.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "LeastSquares[List[List[1, 2], List[3, 4]], List[1, 2, 3]]",
        "LeastSquares[List[List[1, 2, 3], List[4, 5, 6]], List[1, 2]]",
        "LeastSquares[List[List[1, 2], List[2, 4]], List[1, 2]]",
        "LeastSquares[List[List[1, 2], List[3, 4]]]",
        "LQDecomposition[List[List[1, 2], List[3, 4], List[5, 6]]]",
        "LQDecomposition[List[List[1, 2], List[2, 4]]]",
        "LQDecomposition[x]",
    ] {
        let e = read(input);
        assert_eq!(eval(&ctx, &e), e, "{input}");
    }
}
