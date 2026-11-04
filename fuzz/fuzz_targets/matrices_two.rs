//! Matrix decompositions II (CAS-DOC-47 §6, CAS-480), on a random integer
//! matrix `A` up to 5 × 5, against the target's own rational arithmetic:
//! - `HessenbergDecomposition`: `P H = A P`, `P` invertible, `H` zero below
//!   its subdiagonal;
//! - for `S = A + Aᵀ`: `L D Lᵀ = S` with `L` unit lower triangular and `D`
//!   diagonal when `LDLDecomposition` answers, `PositiveDefiniteQ` against
//!   the leading principal minors, `Inertia` summing to `n`, `{n, 0, 0}`
//!   exactly when positive definite, and the signs of `D` when it exists;
//! - `MatrixNorm` 1, ∞ and the Frobenius norm's square against the sums;
//!   `ConditionNumber` against `‖A‖₁ ‖A⁻¹‖₁` or `Infinity`;
//! - `ConjugateGradient` solving `(AᵀA + I) x = b`;
//! - `CompanionMatrix`'s characteristic polynomial the monic `p`, and
//!   `Det[SylvesterMatrix[p, q, x]]` the resultant.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

type Matrix = Vec<Vec<Rational>>;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn to_expr(m: &Matrix) -> Expr {
    list(m.iter().map(|r| list(r.iter().cloned().map(Expr::from))))
}

fn from_expr(e: &Expr) -> Option<Matrix> {
    if e.head() != Some(&Expr::from(Symbol::LIST)) {
        return None;
    }
    e.args()
        .iter()
        .map(|r| r.args().iter().map(Expr::to_rational).collect())
        .collect()
}

fn zero() -> Rational {
    Rational::from(Integer::ZERO)
}

fn one() -> Rational {
    Rational::from(Integer::ONE)
}

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    a.iter()
        .map(|row| {
            (0..b[0].len())
                .map(|j| {
                    row.iter()
                        .zip(b)
                        .fold(zero(), |s, (x, r)| &s + &(x * &r[j]))
                })
                .collect()
        })
        .collect()
}

fn transpose(a: &Matrix) -> Matrix {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j].clone()).collect())
        .collect()
}

/// The determinant by Gaussian elimination over ℚ.
fn det(m: &Matrix) -> Rational {
    let mut a = m.clone();
    let n = a.len();
    let mut d = one();
    for k in 0..n {
        let Some(p) = (k..n).find(|&i| !a[i][k].is_zero()) else {
            return zero();
        };
        if p != k {
            a.swap(p, k);
            d = -d;
        }
        d = &d * &a[k][k];
        let pivot = a[k].clone();
        for row in a.iter_mut().skip(k + 1) {
            let f = &row[k] / &pivot[k];
            for (x, y) in row.iter_mut().zip(&pivot) {
                *x = &*x - &(&f * y);
            }
        }
    }
    d
}

/// The inverse over ℚ by Gauss–Jordan, `None` when singular.
fn inverse(m: &Matrix) -> Option<Matrix> {
    let n = m.len();
    let mut a: Matrix = m
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut row = r.clone();
            row.extend((0..n).map(|j| if i == j { one() } else { zero() }));
            row
        })
        .collect();
    for k in 0..n {
        let p = (k..n).find(|&i| !a[i][k].is_zero())?;
        a.swap(p, k);
        let c = a[k][k].clone();
        for x in &mut a[k] {
            *x = &*x / &c;
        }
        let pivot = a[k].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i != k && !row[k].is_zero() {
                let f = row[k].clone();
                for (x, y) in row.iter_mut().zip(&pivot) {
                    *x = &*x - &(&f * y);
                }
            }
        }
    }
    Some(a.into_iter().map(|r| r[n..].to_vec()).collect())
}

fn one_norm(a: &Matrix) -> Rational {
    (0..a[0].len())
        .map(|j| a.iter().fold(zero(), |s, r| &s + &r[j].abs()))
        .max()
        .expect("a column")
}

fn infinity_norm(a: &Matrix) -> Rational {
    a.iter()
        .map(|r| r.iter().fold(zero(), |s, x| &s + &x.abs()))
        .max()
        .expect("a row")
}

fn call(s: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(s, args)
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let n = u.int_in_range(1..=5usize)?;
    let a: Matrix = (0..n)
        .map(|_| {
            (0..n)
                .map(|_| {
                    u.int_in_range(-5..=5i64)
                        .map(|v| Rational::from(Integer::from(v)))
                })
                .collect()
        })
        .collect::<Result<_>>()?;
    let ae = to_expr(&a);
    let what = format!("{ae:?}");
    match u.int_in_range(0..=4u8)? {
        0 => {
            let found = eval(&ctx, &call(Symbol::HESSENBERG_DECOMPOSITION, vec![ae]));
            let [p, h] = found.args() else {
                panic!("HessenbergDecomposition[{what}] is {found:?}");
            };
            let (p, h) = (
                from_expr(p).expect("rationals"),
                from_expr(h).expect("rationals"),
            );
            assert_eq!(product(&p, &h), product(&a, &p), "{what}: P H against A P");
            assert!(!det(&p).is_zero(), "{what}: P singular");
            for (i, row) in h.iter().enumerate() {
                for x in row.iter().take(i.saturating_sub(1)) {
                    assert!(x.is_zero(), "{what}: H {h:?} not Hessenberg");
                }
            }
        }
        1 => {
            let s: Matrix = product(&one_matrix(n), &a)
                .iter()
                .zip(&transpose(&a))
                .map(|(r, t)| r.iter().zip(t).map(|(x, y)| x + y).collect())
                .collect();
            let se = to_expr(&s);
            let what = format!("{se:?}");
            let minors = (1..=n).all(|k| {
                let sub: Matrix = s[..k].iter().map(|r| r[..k].to_vec()).collect();
                det(&sub).is_positive()
            });
            let definite = eval(&ctx, &call(Symbol::POSITIVE_DEFINITE_Q, vec![se.clone()]));
            assert_eq!(
                definite,
                Expr::from(if minors { Symbol::TRUE } else { Symbol::FALSE }),
                "PositiveDefiniteQ[{what}]"
            );
            let inertia = eval(&ctx, &call(Symbol::INERTIA, vec![se.clone()]));
            let counts: Vec<i64> = inertia
                .args()
                .iter()
                .map(|c| c.as_i64().expect("a count"))
                .collect();
            assert_eq!(
                counts.iter().sum::<i64>(),
                i64::try_from(n).expect("small"),
                "Inertia[{what}]"
            );
            assert_eq!(
                minors,
                counts[0] == i64::try_from(n).expect("small"),
                "Inertia[{what}]: {counts:?}"
            );
            let found = eval(&ctx, &call(Symbol::LDL_DECOMPOSITION, vec![se]));
            if let [_, l, d] = found.args() {
                let (l, d) = (
                    from_expr(l).expect("rationals"),
                    from_expr(d).expect("rationals"),
                );
                for i in 0..n {
                    assert!(
                        l[i][i] == one() && l[i][i + 1..].iter().all(Rational::is_zero),
                        "{what}: L {l:?}"
                    );
                    assert!(
                        (0..n).all(|j| j == i || d[i][j].is_zero()),
                        "{what}: D {d:?}"
                    );
                }
                assert_eq!(
                    product(&product(&l, &d), &transpose(&l)),
                    s,
                    "{what}: L D Lᵀ"
                );
                let positive = (0..n).filter(|&i| d[i][i].is_positive()).count();
                let negative = (0..n).filter(|&i| d[i][i].is_negative()).count();
                assert_eq!(
                    (
                        i64::try_from(positive).expect("small"),
                        i64::try_from(negative).expect("small")
                    ),
                    (counts[0], counts[1]),
                    "{what}: the signs of D against Inertia {counts:?}"
                );
            }
        }
        2 => {
            let norm = |k: Expr| eval(&ctx, &call(Symbol::MATRIX_NORM, vec![ae.clone(), k]));
            assert_eq!(
                norm(Expr::from(1)).to_rational(),
                Some(one_norm(&a)),
                "1-norm of {what}"
            );
            assert_eq!(
                norm(Expr::from(Symbol::INFINITY)).to_rational(),
                Some(infinity_norm(&a)),
                "∞-norm of {what}"
            );
            let frobenius = norm(Expr::symbol("frobenius"));
            let square = eval(&ctx, &Expr::call(Symbol::POWER, [frobenius, Expr::from(2)]));
            let sum = a.iter().flatten().fold(zero(), |s, x| &s + &(x * x));
            assert_eq!(square.to_rational(), Some(sum), "Frobenius norm of {what}");
            let cond = eval(&ctx, &call(Symbol::CONDITION_NUMBER, vec![ae]));
            match inverse(&a) {
                Some(b) => assert_eq!(
                    cond.to_rational(),
                    Some(&one_norm(&a) * &one_norm(&b)),
                    "ConditionNumber[{what}]"
                ),
                None => assert_eq!(
                    cond,
                    Expr::from(Symbol::INFINITY),
                    "ConditionNumber[{what}]"
                ),
            }
        }
        3 => {
            // (AᵀA + I) x = b, positive definite.
            let m: Matrix = product(&transpose(&a), &a)
                .into_iter()
                .enumerate()
                .map(|(i, mut r)| {
                    r[i] = &r[i] + &one();
                    r
                })
                .collect();
            let b: Vec<Rational> = (0..n)
                .map(|_| {
                    u.int_in_range(-5..=5i64)
                        .map(|v| Rational::from(Integer::from(v)))
                })
                .collect::<Result<_>>()?;
            let be = list(b.iter().cloned().map(Expr::from));
            let found = eval(
                &ctx,
                &call(Symbol::CONJUGATE_GRADIENT, vec![to_expr(&m), be]),
            );
            let x: Vec<Rational> = found
                .args()
                .iter()
                .map(|v| {
                    v.to_rational()
                        .unwrap_or_else(|| panic!("ConjugateGradient: {found:?}"))
                })
                .collect();
            let column: Matrix = x.into_iter().map(|v| vec![v]).collect();
            let image: Vec<Rational> = product(&m, &column)
                .into_iter()
                .map(|r| r[0].clone())
                .collect();
            assert_eq!(image, b, "ConjugateGradient[{:?}]", to_expr(&m));
        }
        _ => {
            // Companion and Sylvester matrices of rows of A as coefficients.
            let x = Expr::symbol("x");
            let poly = |row: &[Rational]| {
                Expr::call(
                    Symbol::PLUS,
                    row.iter()
                        .enumerate()
                        .map(|(k, c)| {
                            Expr::from(c.clone()) * x.pow(i64::try_from(k).expect("small"))
                        })
                        .chain([x.pow(i64::try_from(row.len()).expect("small"))])
                        .collect::<Vec<_>>(),
                )
            };
            let p = poly(&a[0]);
            let companion = call(Symbol::COMPANION_MATRIX, vec![p.clone(), x.clone()]);
            let char_poly = call(
                Symbol::CHARACTERISTIC_POLYNOMIAL,
                vec![companion, x.clone()],
            );
            let difference = eval(&ctx, &Expr::call(Symbol::EXPAND, [char_poly - p.clone()]));
            assert_eq!(difference, Expr::from(0), "CompanionMatrix of {p:?}");
            let q = poly(&a[n - 1][..n.min(3)]);
            let sylvester = call(
                Symbol::SYLVESTER_MATRIX,
                vec![p.clone(), q.clone(), x.clone()],
            );
            let d = eval(&ctx, &call(Symbol::DET, vec![sylvester]));
            let r = eval(
                &ctx,
                &call(Symbol::RESULTANT, vec![p.clone(), q.clone(), x]),
            );
            assert_eq!(d, r, "Det[SylvesterMatrix] of {p:?} and {q:?}");
        }
    }
    Ok(())
}

fn one_matrix(n: usize) -> Matrix {
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| if i == j { one() } else { zero() })
                .collect()
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
