//! Matrix decompositions III (CAS-DOC-51 §7, CAS-501), on a random integer
//! matrix `A`:
//! - `SchurDecomposition` (`A` up to 2 × 2 exact, 3 × 3 as balls): `T`
//!   upper triangular, `Qᴴ Q = I` and `Q T Qᴴ = A`;
//! - `SingularValueDecomposition` (the same sizes): `s` descending, `U`
//!   and `V` unitary, `U diag(s) Vᴴ = A`;
//! - `Eigenvalues` of `A` as balls, each `det(λ I − A)` holding 0;
//! - `Eigenvectors` of `A` as balls: each column `x` with `A x` parallel,
//!   `A x − (A x)_k x` holding 0 for its component `x_k = 1`;
//! - `JacobiLinearSolve` and `GaussSeidelLinearSolve` on `A` made strictly
//!   diagonally dominant: balls holding the exact solution;
//! - `SmithDecomposition` of a 2 × 2 matrix of polynomials in `x`: `U A V
//!   = D` diagonal, `d₁ | d₂`, `U` and `V` of constant determinants.
//!
//! Exactly, the identities are checked by `N` at 30 digits within 10⁻¹⁵;
//! with balls, by balls holding 0.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::Rational;
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn dot(a: &Expr, b: &Expr) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])
}

fn adjoint(m: &Expr) -> Expr {
    Expr::call(
        Symbol::CONJUGATE,
        [Expr::call(Symbol::TRANSPOSE, [m.clone()])],
    )
}

fn identity(n: usize) -> Expr {
    Expr::call(
        Symbol::IDENTITY_MATRIX,
        [Expr::from(i64::try_from(n).expect("small"))],
    )
}

fn entries(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(Symbol::LIST)) {
        e.args().iter().flat_map(entries).collect()
    } else {
        vec![e.clone()]
    }
}

/// Whether each entry of `m` is 0: a ball holding it, or within 10⁻¹⁵ at
/// 30 digits.
fn vanishes(ctx: &Context, m: &Expr, balls: bool) -> bool {
    let value = if balls {
        eval(ctx, m)
    } else {
        eval(ctx, &Expr::call(Symbol::N, [m.clone(), Expr::from(30)]))
    };
    entries(&value).iter().all(|x| {
        use cas_core::View;
        match x.view() {
            View::Ball(b) if balls => b.contains_rational(&Rational::ZERO),
            View::ComplexBall(z) if balls => z.may_be_zero(),
            _ if *x == Expr::from(0) => true,
            _ => {
                // `|Re| + |Im|`: `Abs` of a complex ball holding 0 stays.
                let size = Expr::call(
                    Symbol::PLUS,
                    [
                        Expr::call(Symbol::ABS, [Expr::call(Symbol::RE, [x.clone()])]),
                        Expr::call(Symbol::ABS, [Expr::call(Symbol::IM, [x.clone()])]),
                    ],
                );
                let small = Expr::call(Symbol::LESS, [size, Expr::from(10).pow(-15)]);
                !balls && eval(ctx, &small).as_symbol() == Some(Symbol::TRUE)
            }
        }
    })
}

/// For each column `x` of `p`, `A x − (A x)_k x`, `k` the component that
/// is exactly 1: 0 for an eigenvector, whatever its eigenvalue.
fn eigen_gaps(ctx: &Context, a: &Expr, p: &Expr) -> Expr {
    let n = p.args().len();
    let one = |x: &Expr| match x.view() {
        cas_core::View::Ball(b) => b.is_exact() && b.contains_rational(&Rational::ONE),
        _ => *x == Expr::from(1),
    };
    let mut gaps = Vec::with_capacity(n);
    for j in 0..n {
        let x: Vec<Expr> = p.args().iter().map(|row| row.args()[j].clone()).collect();
        let k = x
            .iter()
            .position(one)
            .unwrap_or_else(|| panic!("{p:?}: no 1 in column {j}"));
        let x = Expr::call(Symbol::LIST, x);
        let image = eval(ctx, &dot(a, &x));
        let scale = image.args()[k].clone();
        gaps.push(image - scale * x);
    }
    Expr::call(Symbol::LIST, gaps)
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn unitary_checks(ctx: &Context, a: &Expr, n: usize, balls: bool) {
    let found = eval(ctx, &Expr::call(Symbol::SCHUR_DECOMPOSITION, [a.clone()]));
    if is_list(&found) {
        let [q, t] = found.args() else {
            panic!("{a:?}: {found:?}");
        };
        for (i, row) in t.args().iter().enumerate() {
            for x in &row.args()[..i] {
                assert_eq!(*x, Expr::from(0), "{a:?}: T = {t:?}");
            }
        }
        assert!(
            vanishes(ctx, &(dot(&adjoint(q), q) - identity(n)), balls),
            "{a:?}: Q = {q:?}"
        );
        assert!(
            vanishes(ctx, &(dot(&dot(q, t), &adjoint(q)) - a.clone()), balls),
            "{a:?}: {found:?}"
        );
    }
    let found = eval(
        ctx,
        &Expr::call(Symbol::SINGULAR_VALUE_DECOMPOSITION, [a.clone()]),
    );
    if is_list(&found) {
        let [u, s, v] = found.args() else {
            panic!("{a:?}: {found:?}");
        };
        for w in [u, v] {
            assert!(
                vanishes(ctx, &(dot(&adjoint(w), w) - identity(n)), balls),
                "{a:?}: {found:?}"
            );
        }
        for pair in s.args().windows(2) {
            let order = Expr::call(Symbol::GREATER_EQUAL, [pair[0].clone(), pair[1].clone()]);
            let decided = eval(ctx, &order);
            assert_ne!(decided.as_symbol(), Some(Symbol::FALSE), "{a:?}: s = {s:?}");
        }
        let diagonal = Expr::call(Symbol::DIAGONAL_MATRIX, [s.clone()]);
        assert!(
            vanishes(
                ctx,
                &(dot(&dot(u, &diagonal), &adjoint(v)) - a.clone()),
                balls
            ),
            "{a:?}: {found:?}"
        );
    }
}

fn integer_matrix(u: &mut Unstructured<'_>, n: usize) -> Result<Vec<Vec<i64>>> {
    (0..n)
        .map(|_| (0..n).map(|_| u.int_in_range(-5..=5)).collect())
        .collect()
}

fn to_expr(m: &[Vec<i64>]) -> Expr {
    list(m.iter().map(|r| list(r.iter().map(|&x| Expr::from(x)))))
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let n = u.int_in_range(1..=3)?;
    let m = integer_matrix(u, n)?;
    let a = to_expr(&m);
    if n <= 2 {
        unitary_checks(&ctx, &a, n, false);
    }
    let balls = eval(&ctx, &Expr::call(Symbol::N, [a.clone(), Expr::from(30)]));
    unitary_checks(&ctx, &balls, n, true);
    // `Eigenvalues` (CAS-502): each `det(λ I − A)` holds 0.
    let values = eval(&ctx, &Expr::call(Symbol::EIGENVALUES, [balls.clone()]));
    if is_list(&values) {
        // `Det` of a ball matrix that may be singular stays: `det(λ I − A)`
        // from the characteristic polynomial.
        let c = eval(
            &ctx,
            &Expr::call(Symbol::CHARACTERISTIC_POLYNOMIAL, [balls.clone()]),
        );
        let dets = values.args().iter().map(|lambda| {
            let terms: Vec<Expr> = c
                .args()
                .iter()
                .enumerate()
                .map(|(i, ci)| {
                    ci.clone() * lambda.clone().pow(i64::try_from(n - i).expect("small"))
                })
                .collect();
            Expr::call(Symbol::PLUS, terms)
        });
        assert!(vanishes(&ctx, &list(dets), true), "{balls:?}: {values:?}");
    }
    let p = eval(&ctx, &Expr::call(Symbol::EIGENVECTORS, [balls.clone()]));
    if is_list(&p) {
        assert!(
            vanishes(&ctx, &eigen_gaps(&ctx, &balls, &p), true),
            "{balls:?}: {p:?}"
        );
    }
    // Diagonally dominant: the solvers' balls hold the solution.
    let mut d = m.clone();
    for (i, row) in d.iter_mut().enumerate() {
        let off: i64 = row
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, x)| x.abs())
            .sum();
        row[i] = off + 1 + i64::from(u.int_in_range(0..=3)?);
    }
    let b: Vec<i64> = (0..n)
        .map(|_| u.int_in_range(-9..=9))
        .collect::<Result<_>>()?;
    let (d, b) = (to_expr(&d), list(b.iter().map(|&x| Expr::from(x))));
    let exact = eval(&ctx, &dot(&Expr::call(Symbol::INVERSE, [d.clone()]), &b));
    let eps = Expr::from(10).pow(-12);
    for solver in [
        Symbol::JACOBI_LINEAR_SOLVE,
        Symbol::GAUSS_SEIDEL_LINEAR_SOLVE,
    ] {
        let found = eval(
            &ctx,
            &Expr::call(solver, [d.clone(), b.clone(), eps.clone()]),
        );
        assert!(is_list(&found), "{d:?} {b:?}: {found:?}");
        assert!(
            vanishes(&ctx, &list([found.clone() - exact.clone()]), true),
            "{d:?} {b:?}: {found:?}, {exact:?}"
        );
    }
    smith_checks(&ctx, u)
}

/// A random polynomial in `x` of degree at most 2.
fn polynomial(u: &mut Unstructured<'_>) -> Result<Expr> {
    let x = Expr::symbol("x");
    let terms: Vec<Expr> = (0..3)
        .map(|k| Ok(Expr::from(u.int_in_range(-3..=3)?) * x.clone().pow(k)))
        .collect::<Result<_>>()?;
    Ok(Expr::call(Symbol::PLUS, terms))
}

fn smith_checks(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let x = Expr::symbol("x");
    let mut rows = Vec::new();
    for _ in 0..2 {
        rows.push(list([polynomial(u)?, polynomial(u)?]));
    }
    let a = eval(ctx, &list(rows));
    // Integer or constant matrices take ℤ's form or stay: polynomials only.
    if !a.contains(&x) {
        return Ok(());
    }
    let found = eval(ctx, &Expr::call(Symbol::SMITH_DECOMPOSITION, [a.clone()]));
    let [left, d, right] = found.args() else {
        panic!("{a:?}: {found:?}");
    };
    let normal = |e: Expr| {
        eval(
            ctx,
            &Expr::call(Symbol::TOGETHER, [Expr::call(Symbol::EXPAND, [e])]),
        )
    };
    let product = dot(&dot(left, &a), right);
    let gap = normal(eval(ctx, &(product - d.clone())));
    assert!(
        entries(&gap)
            .iter()
            .all(|g| normal(g.clone()) == Expr::from(0)),
        "{a:?}: {found:?}"
    );
    assert_eq!(
        normal(d.args()[0].args()[1].clone()),
        Expr::from(0),
        "{found:?}"
    );
    assert_eq!(
        normal(d.args()[1].args()[0].clone()),
        Expr::from(0),
        "{found:?}"
    );
    let (d1, d2) = (d.args()[0].args()[0].clone(), d.args()[1].args()[1].clone());
    if d1 != Expr::from(0) {
        let r = eval(
            ctx,
            &Expr::call(Symbol::POLYNOMIAL_REMAINDER, [d2, d1, x.clone()]),
        );
        assert_eq!(normal(r), Expr::from(0), "{a:?}: {found:?}");
    }
    for m in [left, right] {
        let det = normal(Expr::call(Symbol::DET, [m.clone()]));
        assert!(
            det.to_rational().is_some_and(|q| !q.is_zero()),
            "{a:?}: det {m:?} = {det:?}"
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let _ = run(&mut Unstructured::new(data));
});
