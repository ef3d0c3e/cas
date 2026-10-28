//! Matrix decompositions and functions (CAS-DOC-30 §5, CAS-348), one
//! check per input on a random integer matrix `A` of size 2 to 4:
//! - `LUDecomposition`: `L U = A[p]`;
//! - `QRDecomposition`, where it answers: `Q R = A`;
//! - `CholeskyDecomposition` of `B Bᵀ + I`: `L Lᵀ` is it;
//! - for a matrix with chosen integer eigenvalues, upper triangular and
//!   then conjugated by a unimodular matrix: `JordanDecomposition` gives
//!   `P J P⁻¹ = A`, the minimal polynomial annihilates `A`, and
//!   `MatrixPower[A, n]` at `n = 0 … 4` is the products.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

/// The matrix product: `Times`, as `Dot` of two matrices is their inner
/// product.
fn dot(a: &Expr, b: &Expr) -> Expr {
    Expr::call(Symbol::TIMES, [a.clone(), b.clone()])
}

fn matrix(rows: &[Vec<i64>]) -> Expr {
    Expr::call(
        Symbol::LIST,
        rows.iter()
            .map(|r| {
                Expr::call(
                    Symbol::LIST,
                    r.iter().map(|&x| Expr::from(x)).collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>(),
    )
}

/// Whether two matrices agree entry by entry once put together.
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

fn stayed(e: &Expr, head: &str) -> bool {
    e.head() == Some(&Expr::from(Symbol::new(head)))
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let n = u.int_in_range(2..=4usize)?;
    let mut entries = vec![vec![0i64; n]; n];
    for row in &mut entries {
        for x in row.iter_mut() {
            *x = u.int_in_range(-4..=4)?;
        }
    }
    let a = matrix(&entries);
    // LU.
    let found = eval(&ctx, &call("LUDecomposition", vec![a.clone()]));
    if let [p, l, upper] = found.args() {
        let permuted = Expr::call(
            Symbol::LIST,
            p.args()
                .iter()
                .map(|i| {
                    Expr::call(
                        Symbol::PART,
                        [a.clone(), Expr::from(i.as_i64().expect("an index") + 1)],
                    )
                })
                .collect::<Vec<_>>(),
        );
        assert!(
            same(&ctx, &dot(l, upper), &permuted),
            "LU of {a:?}: {found:?}"
        );
    }
    // QR.
    let found = eval(&ctx, &call("QRDecomposition", vec![a.clone()]));
    if !stayed(&found, "QRDecomposition") {
        let [q, r] = found.args() else {
            panic!("QR of {a:?}: {found:?}");
        };
        assert!(same(&ctx, &dot(q, r), &a), "QR of {a:?}: {found:?}");
    }
    // Cholesky of B Bᵀ + I.
    let spd = eval(
        &ctx,
        &(dot(&a, &Expr::call(Symbol::TRANSPOSE, [a.clone()]))
            + Expr::call(
                Symbol::IDENTITY_MATRIX,
                [Expr::from(i64::try_from(n).expect("small"))],
            )),
    );
    let l = eval(&ctx, &call("CholeskyDecomposition", vec![spd.clone()]));
    assert!(!stayed(&l, "CholeskyDecomposition"), "Cholesky of {spd:?}");
    let back = dot(&l, &Expr::call(Symbol::TRANSPOSE, [l.clone()]));
    assert!(same(&ctx, &back, &spd), "Cholesky of {spd:?}: {l:?}");
    // A matrix with integer eigenvalues: T upper triangular, S unimodular.
    let mut t = vec![vec![0i64; n]; n];
    for i in 0..n {
        t[i][i] = u.int_in_range(-2..=2)?;
        for j in i + 1..n {
            t[i][j] = u.int_in_range(-1..=1)?;
        }
    }
    let mut s = vec![vec![0i64; n]; n];
    for i in 0..n {
        s[i][i] = 1;
        for j in i + 1..n {
            s[i][j] = u.int_in_range(-1..=1)?;
        }
    }
    let (t, s) = (matrix(&t), matrix(&s));
    let b = eval(&ctx, &dot(&dot(&s, &t), &Expr::call(Symbol::INVERSE, [s])));
    let found = eval(&ctx, &call("JordanDecomposition", vec![b.clone()]));
    if !stayed(&found, "JordanDecomposition") {
        let [p, j] = found.args() else {
            panic!("Jordan of {b:?}: {found:?}");
        };
        let back = dot(&dot(p, j), &Expr::call(Symbol::INVERSE, [p.clone()]));
        assert!(same(&ctx, &back, &b), "Jordan of {b:?}: {found:?}");
        // The minimal polynomial annihilates B, by Horner on matrices.
        let coefficients = eval(&ctx, &call("MinimalPolynomial", vec![b.clone()]));
        let ident = Expr::call(
            Symbol::IDENTITY_MATRIX,
            [Expr::from(i64::try_from(n).expect("small"))],
        );
        let mut value = eval(&ctx, &(ident.clone() * Expr::from(0)));
        for c in coefficients.args() {
            value = eval(&ctx, &(dot(&value, &b) + c.clone() * ident.clone()));
        }
        let zero = eval(&ctx, &(ident * Expr::from(0)));
        assert!(
            same(&ctx, &value, &zero),
            "the minimal polynomial {coefficients:?} of {b:?}"
        );
        // MatrixPower at n = 0 … 4.
        let power = eval(
            &ctx,
            &call("MatrixPower", vec![b.clone(), Expr::symbol("n")]),
        );
        if !stayed(&power, "MatrixPower") {
            for k in 0..=4i64 {
                let at = power.substitute(&[(Expr::symbol("n"), Expr::from(k))]);
                let direct = Expr::call(Symbol::POWER, [b.clone(), Expr::from(k)]);
                // Blocks at 0 are generic in n: from the size on.
                if k >= i64::try_from(n).expect("small") || !entries_have_zero_eigenvalue(&ctx, &b)
                {
                    assert!(
                        same(&ctx, &at, &direct),
                        "MatrixPower[{b:?}, {k}]: {power:?}"
                    );
                }
            }
        }
    }
    Ok(())
}

/// Whether 0 is an eigenvalue of `b`.
fn entries_have_zero_eigenvalue(ctx: &Context, b: &Expr) -> bool {
    eval(ctx, &Expr::call(Symbol::DET, [b.clone()])) == Expr::from(0)
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
