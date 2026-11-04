//! Modular arithmetic II (CAS-DOC-44 §5, CAS-461), one check per input,
//! against the target's own arithmetic modulo `n`:
//! - a random matrix up to 5 × 5 modulo a small prime or composite:
//!   `Det` the integer determinant reduced; `Inverse` times the matrix the
//!   identity when the determinant is a unit, the call staying otherwise;
//!   `CharacteristicPolynomial` at two points `det(t I − m)`; for a prime,
//!   `RowReduce` the reduced echelon form, `MatrixRank` its rank, each
//!   `NullSpace` vector in the kernel and as many as the nullity, `LinSolve`
//!   solving, finding none, or staying for a family, and `Eigenvalues`
//!   roots of the characteristic polynomial whose sum is the trace and
//!   product the determinant;
//! - `Solve` of a random polynomial modulo a prime: the roots by brute
//!   force;
//! - `ChineseRemainder` of two modular values: the congruences, or none
//!   when they disagree;
//! - `PowerMod` of a polynomial modulo a prime and a monic polynomial:
//!   repeated multiplication.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::Integer;
use libfuzzer_sys::fuzz_target;

type Matrix = Vec<Vec<i64>>;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn modular(a: Expr, n: i64) -> Expr {
    Expr::call(Symbol::MODULAR, [a, Expr::from(n)])
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn matrix_expr(m: &Matrix) -> Expr {
    list(m.iter().map(|r| list(r.iter().map(|&x| Expr::from(x)))))
}

fn is(e: &Expr, s: Symbol) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(s)
}

/// An entry of an answer modulo `n`: a modular value or an integer.
fn residue(e: &Expr, n: i64) -> i64 {
    let a = match e.args() {
        [a, k] if is(e, Symbol::MODULAR) && k.as_i64() == Some(n) => a,
        _ => e,
    };
    a.as_i64()
        .unwrap_or_else(|| panic!("{e:?} is no residue modulo {n}"))
        .rem_euclid(n)
}

fn matrix_of(e: &Expr, n: i64) -> Matrix {
    assert!(is(e, Symbol::LIST), "{e:?} is no matrix");
    e.args()
        .iter()
        .map(|r| r.args().iter().map(|x| residue(x, n)).collect())
        .collect()
}

fn det(m: &[Vec<i64>], n: i64) -> i64 {
    if m.is_empty() {
        return 1 % n;
    }
    let mut total = 0;
    for (j, x) in m[0].iter().enumerate() {
        if *x == 0 {
            continue;
        }
        let minor: Matrix = m[1..]
            .iter()
            .map(|r| {
                r.iter()
                    .enumerate()
                    .filter(|&(c, _)| c != j)
                    .map(|(_, v)| *v)
                    .collect()
            })
            .collect();
        let term = x * det(&minor, n) % n;
        total = if j % 2 == 0 {
            total + term
        } else {
            total - term
        };
    }
    total.rem_euclid(n)
}

fn product(a: &[Vec<i64>], b: &[Vec<i64>], n: i64) -> Matrix {
    a.iter()
        .map(|row| {
            (0..b[0].len())
                .map(|j| {
                    row.iter()
                        .zip(b)
                        .map(|(x, r)| x * r[j])
                        .sum::<i64>()
                        .rem_euclid(n)
                })
                .collect()
        })
        .collect()
}

fn inverse(a: i64, p: i64) -> i64 {
    (1..p).find(|x| a * x % p == 1).expect("a unit")
}

/// The reduced row echelon form modulo the prime `p` and its pivots, among
/// the first `columns` columns.
fn rref(m: &[Vec<i64>], columns: usize, p: i64) -> (Matrix, Vec<usize>) {
    let mut a = m.to_vec();
    let mut pivots = Vec::new();
    for j in 0..columns {
        let r = pivots.len();
        let Some(k) = (r..a.len()).find(|&k| a[k][j] != 0) else {
            continue;
        };
        a.swap(r, k);
        let c = inverse(a[r][j], p);
        for x in &mut a[r] {
            *x = *x * c % p;
        }
        let pivot = a[r].clone();
        for (k, row) in a.iter_mut().enumerate() {
            if k != r && row[j] != 0 {
                let c = row[j];
                for (x, y) in row.iter_mut().zip(&pivot) {
                    *x = (*x - c * y).rem_euclid(p);
                }
            }
        }
        pivots.push(j);
    }
    (a, pivots)
}

fn random_matrix(u: &mut Unstructured<'_>, rows: usize, columns: usize, n: i64) -> Result<Matrix> {
    (0..rows)
        .map(|_| (0..columns).map(|_| u.int_in_range(0..=n - 1)).collect())
        .collect()
}

fn matrices(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let prime = u.arbitrary::<bool>()?;
    let n = if prime {
        *u.choose(&[2i64, 3, 5, 7, 11, 13, 101])?
    } else {
        *u.choose(&[4i64, 6, 8, 9, 12, 15])?
    };
    let rows = u.int_in_range(1..=5usize)?;
    let columns = if u.arbitrary()? {
        rows
    } else {
        u.int_in_range(1..=5usize)?
    };
    let a = random_matrix(u, rows, columns, n)?;
    let m = modular(matrix_expr(&a), n);
    let call = |s: Symbol, args: Vec<Expr>| Expr::call(s, args);
    if rows == columns {
        let d = det(&a, n);
        let found = eval(ctx, &call(Symbol::DET, vec![m.clone()]));
        assert_eq!(residue(&found, n), d, "Det of {a:?} modulo {n}: {found:?}");
        let inverse_call = call(Symbol::INVERSE, vec![m.clone()]);
        let found = eval(ctx, &inverse_call);
        if Integer::from(d).gcd(&Integer::from(n)) == Integer::ONE {
            let b = matrix_of(&found, n);
            let identity: Matrix = (0..rows)
                .map(|i| (0..rows).map(|j| i64::from(i == j)).collect())
                .collect();
            assert_eq!(
                product(&b, &a, n),
                identity,
                "{inverse_call:?} is {found:?}"
            );
        } else {
            assert!(is(&found, Symbol::INVERSE), "{inverse_call:?} is {found:?}");
        }
        let x = Expr::symbol("x");
        let poly = eval(
            ctx,
            &call(
                Symbol::CHARACTERISTIC_POLYNOMIAL,
                vec![m.clone(), x.clone()],
            ),
        );
        for _ in 0..2 {
            let t = u.int_in_range(0..=n - 1)?;
            let shifted: Matrix = (0..rows)
                .map(|i| (0..rows).map(|j| i64::from(i == j) * t - a[i][j]).collect())
                .collect();
            let at = eval(ctx, &poly.substitute(&[(x.clone(), Expr::from(t))]));
            assert_eq!(
                residue(&at, n),
                det(&shifted, n),
                "CharacteristicPolynomial of {a:?} modulo {n} is {poly:?}, at {t}"
            );
        }
    }
    if !prime {
        return Ok(());
    }
    let p = n;
    // The reduced form is unique.
    let (reduced, pivots) = rref(&a, columns, p);
    let found = eval(ctx, &call(Symbol::ROW_REDUCE, vec![m.clone()]));
    assert_eq!(
        matrix_of(&found, p),
        reduced,
        "RowReduce of {a:?} modulo {p}"
    );
    let rank = eval(ctx, &call(Symbol::MATRIX_RANK, vec![m.clone()]));
    assert_eq!(
        rank.as_i64(),
        i64::try_from(pivots.len()).ok(),
        "MatrixRank of {a:?} modulo {p}"
    );
    let kernel = eval(ctx, &call(Symbol::NULL_SPACE, vec![m.clone()]));
    assert_eq!(
        kernel.args().len(),
        columns - pivots.len(),
        "NullSpace of {a:?} modulo {p}: {kernel:?}"
    );
    for v in kernel.args() {
        let v: Matrix = v.args().iter().map(|x| vec![residue(x, p)]).collect();
        assert!(
            v.iter().any(|x| x[0] != 0),
            "NullSpace of {a:?} modulo {p}: {kernel:?}"
        );
        let image = product(&a, &v, p);
        assert!(
            image.iter().all(|r| r[0] == 0),
            "NullSpace of {a:?} modulo {p}: {kernel:?}"
        );
    }
    // A system a x = b.
    let b: Vec<i64> = (0..rows)
        .map(|_| u.int_in_range(0..=p - 1))
        .collect::<Result<_>>()?;
    let augmented: Matrix = a
        .iter()
        .zip(&b)
        .map(|(r, y)| r.iter().copied().chain([*y]).collect())
        .collect();
    let (_, all) = rref(&augmented, columns + 1, p);
    let solve_call = call(
        Symbol::LIN_SOLVE,
        vec![m.clone(), list(b.iter().map(|&y| Expr::from(y)))],
    );
    let found = eval(ctx, &solve_call);
    if all.contains(&columns) {
        assert_eq!(found, list([]), "{solve_call:?} is {found:?}");
    } else if pivots.len() < columns {
        assert!(is(&found, Symbol::LIN_SOLVE), "{solve_call:?} is {found:?}");
    } else {
        let x: Matrix = found.args().iter().map(|v| vec![residue(v, p)]).collect();
        let column: Matrix = b.iter().map(|&y| vec![y]).collect();
        assert_eq!(product(&a, &x, p), column, "{solve_call:?} is {found:?}");
    }
    if rows == columns {
        let eigen_call = call(Symbol::EIGENVALUES, vec![m]);
        let found = eval(ctx, &eigen_call);
        if !is(&found, Symbol::EIGENVALUES) {
            let values: Vec<i64> = found.args().iter().map(|v| residue(v, p)).collect();
            assert_eq!(values.len(), rows, "{eigen_call:?} is {found:?}");
            let trace: i64 = (0..rows).map(|i| a[i][i]).sum::<i64>().rem_euclid(p);
            let sum = values.iter().sum::<i64>().rem_euclid(p);
            let product = values.iter().fold(1, |acc, v| acc * v % p);
            assert_eq!(
                (sum, product),
                (trace, det(&a, p)),
                "{eigen_call:?} is {found:?}"
            );
        }
    }
    Ok(())
}

/// A polynomial over `𝔽ₚ`: its coefficients in `0..p` from the constant
/// up, without trailing zeros.
type Poly = Vec<i64>;

fn trim(mut f: Poly) -> Poly {
    while f.last() == Some(&0) {
        f.pop();
    }
    f
}

fn mul(f: &[i64], g: &[i64], p: i64) -> Poly {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let mut h = vec![0; f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        for (j, b) in g.iter().enumerate() {
            h[i + j] = (h[i + j] + a * b) % p;
        }
    }
    trim(h)
}

/// `f` modulo the monic `q`.
fn rem(f: &[i64], q: &[i64], p: i64) -> Poly {
    let mut r = trim(f.to_vec());
    while r.len() >= q.len() {
        let shift = r.len() - q.len();
        let c = *r.last().expect("not 0");
        for (x, b) in r[shift..].iter_mut().zip(q) {
            *x = (*x - c * b).rem_euclid(p);
        }
        r = trim(r);
    }
    r
}

fn poly_expr(f: &[i64]) -> Expr {
    let x = Expr::symbol("x");
    Expr::call(
        Symbol::PLUS,
        f.iter()
            .enumerate()
            .map(|(k, c)| Expr::from(*c) * x.pow(i64::try_from(k).expect("a degree")))
            .collect::<Vec<_>>(),
    )
}

fn value(f: &[i64], t: i64, p: i64) -> i64 {
    f.iter().rev().fold(0, |acc, c| (acc * t + c) % p)
}

fn polynomials(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let p = *u.choose(&[2i64, 3, 5, 7, 11, 13, 101])?;
    let d = u.int_in_range(1..=6usize)?;
    let mut f: Poly = (0..d)
        .map(|_| u.int_in_range(0..=p - 1))
        .collect::<Result<_>>()?;
    f.push(u.int_in_range(1..=p - 1)?);
    // Solve: the roots by brute force, ascending by symmetric residue.
    let symmetric = |r: i64| if 2 * r > p { r - p } else { r };
    let mut roots: Vec<i64> = (0..p)
        .filter(|&t| value(&f, t, p) == 0)
        .map(symmetric)
        .collect();
    roots.sort_unstable();
    let x = Expr::symbol("x");
    let solve_call = Expr::call(Symbol::SOLVE, [modular(poly_expr(&f), p), x.clone()]);
    let found = eval(ctx, &solve_call);
    let expected = list(roots.iter().map(|&r| modular(Expr::from(r), p)));
    assert_eq!(found, eval(ctx, &expected), "{solve_call:?}");
    // PowerMod against repeated multiplication.
    let e = u.int_in_range(1..=4usize)?;
    let mut q: Poly = (0..e)
        .map(|_| u.int_in_range(0..=p - 1))
        .collect::<Result<_>>()?;
    q.push(1);
    let k = u.int_in_range(0..=200i64)?;
    let mut expected = rem(&[1], &q, p);
    for _ in 0..k {
        expected = rem(&mul(&expected, &f, p), &q, p);
    }
    let power_call = Expr::call(
        Symbol::POWER_MOD,
        [
            poly_expr(&f),
            Expr::from(k),
            Expr::from(p),
            poly_expr(&q),
            x.clone(),
        ],
    );
    let found = eval(ctx, &power_call);
    let coefficients = eval(
        ctx,
        &Expr::call(Symbol::COEFFICIENT_LIST, [found.clone(), x]),
    );
    // From the leading coefficient down, as giac's `symb2poly`.
    let got: Poly = trim(
        coefficients
            .args()
            .iter()
            .rev()
            .map(|c| residue(c, p))
            .collect(),
    );
    assert_eq!(got, expected, "{power_call:?} is {found:?}");
    Ok(())
}

fn remainders(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let (n1, n2) = (u.int_in_range(2..=60i64)?, u.int_in_range(2..=60i64)?);
    let (a1, a2) = (u.int_in_range(0..=n1 - 1)?, u.int_in_range(0..=n2 - 1)?);
    let call = Expr::call(
        Symbol::CHINESE_REMAINDER,
        [modular(Expr::from(a1), n1), modular(Expr::from(a2), n2)],
    );
    let found = eval(ctx, &call);
    let lcm = n1
        / Integer::from(n1)
            .gcd(&Integer::from(n2))
            .to_i64()
            .expect("small")
        * n2;
    let solution = (0..lcm).find(|t| t % n1 == a1 && t % n2 == a2);
    match solution {
        Some(t) => assert_eq!(found, eval(ctx, &modular(Expr::from(t), lcm)), "{call:?}"),
        None => assert!(
            is(&found, Symbol::CHINESE_REMAINDER),
            "{call:?} is {found:?}"
        ),
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0..=2u8)? {
        0 => matrices(&ctx, u),
        1 => polynomials(&ctx, u),
        _ => remainders(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
