//! Integer matrices (CAS-DOC-42 §5, CAS-448), on a random matrix of up to
//! 4 × 4 small integers:
//! - `HermiteDecomposition`: `H` in Hermite form (giac's, the entries
//!   above a pivot `p` in `(−p/2, p/2]`), `U A = H`, `det U = ±1`,
//!   and the same `H` for `A`'s rows reversed, as it is unique;
//! - `SmithDecomposition`: `S` diagonal, nonnegative, `d₁ | d₂ | …`, `U A
//!   V = S`, `U` and `V` unimodular, and `d₁ ⋯ d_k` the gcd of `A`'s `k × k`
//!   minors;
//! - `LatticeReduce`: the same lattice (equal Hermite forms) and the LLL
//!   conditions, in exact rationals.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

type Matrix = Vec<Vec<Integer>>;

fn to_expr(m: &Matrix) -> Expr {
    Expr::call(
        Symbol::LIST,
        m.iter()
            .map(|r| Expr::call(Symbol::LIST, r.iter().cloned().map(Expr::from))),
    )
}

fn from_expr(e: &Expr) -> Option<Matrix> {
    e.args()
        .iter()
        .map(|r| r.args().iter().map(Expr::to_integer).collect())
        .collect()
}

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    a.iter()
        .map(|row| {
            (0..b[0].len())
                .map(|j| {
                    row.iter()
                        .zip(b)
                        .fold(Integer::ZERO, |s, (x, r)| &s + &(x * &r[j]))
                })
                .collect()
        })
        .collect()
}

/// The determinant by expansion along the first row.
fn det(m: &[Vec<Integer>]) -> Integer {
    if m.len() == 1 {
        return m[0][0].clone();
    }
    let mut total = Integer::ZERO;
    for (j, x) in m[0].iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        let minor: Vec<Vec<Integer>> = m[1..]
            .iter()
            .map(|r| {
                r.iter()
                    .enumerate()
                    .filter(|&(c, _)| c != j)
                    .map(|(_, v)| v.clone())
                    .collect()
            })
            .collect();
        let term = x * &det(&minor);
        total = if j % 2 == 0 {
            &total + &term
        } else {
            &total - &term
        };
    }
    total
}

fn unimodular(m: &Matrix) -> bool {
    let d = det(m);
    d == Integer::ONE || d == -Integer::ONE
}

/// Whether `h` is in Hermite form.
fn hermite_form(h: &Matrix) -> bool {
    let mut last: Option<usize> = None;
    let mut zero_seen = false;
    for (r, row) in h.iter().enumerate() {
        let Some(c) = row.iter().position(|x| !x.is_zero()) else {
            zero_seen = true;
            continue;
        };
        if zero_seen || last.is_some_and(|l| c <= l) || !h[r][c].is_positive() {
            return false;
        }
        // Above the pivot p, giac's symmetric remainders in (−p/2, p/2].
        for above in &h[..r] {
            let twice = &above[c] + &above[c];
            if twice > h[r][c] || twice <= -h[r][c].clone() {
                return false;
            }
        }
        for below in &h[r + 1..] {
            if !below[c].is_zero() {
                return false;
            }
        }
        last = Some(c);
    }
    true
}

/// The sets of `k` indices out of `0..n`.
fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![Vec::new()];
    }
    (0..n)
        .flat_map(|first| {
            subsets(n, k - 1)
                .into_iter()
                .filter(move |rest| rest.first().is_none_or(|&r| r > first))
                .map(move |mut rest| {
                    rest.insert(0, first);
                    rest
                })
        })
        .collect()
}

/// The gcd of `a`'s `k × k` minors.
fn minors_gcd(a: &Matrix, k: usize) -> Integer {
    let mut g = Integer::ZERO;
    for rows in subsets(a.len(), k) {
        for cols in subsets(a[0].len(), k) {
            let m: Vec<Vec<Integer>> = rows
                .iter()
                .map(|&r| cols.iter().map(|&c| a[r][c].clone()).collect())
                .collect();
            g = g.gcd(&det(&m));
        }
    }
    g
}

/// The Hermite form's rows that are not 0, by the CAS.
fn lattice(ctx: &Context, m: &Matrix) -> Matrix {
    if m.is_empty() {
        return Vec::new();
    }
    let e = Expr::call(Symbol::HERMITE_DECOMPOSITION, [to_expr(m)]);
    let found = ctx.eval(&e).unwrap_or_else(|err| panic!("{e:?}: {err}"));
    let mut h = from_expr(&found.args()[1]).expect("integers");
    h.retain(|r| r.iter().any(|x| !x.is_zero()));
    h
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let (m, n) = (u.int_in_range(1..=4usize)?, u.int_in_range(1..=4usize)?);
    let mut a: Matrix = Vec::with_capacity(m);
    for _ in 0..m {
        let row: Vec<Integer> = (0..n)
            .map(|_| u.int_in_range(-9..=9i64).map(Integer::from))
            .collect::<Result<_>>()?;
        a.push(row);
    }
    let input = to_expr(&a);
    // Hermite.
    let call = Expr::call(Symbol::HERMITE_DECOMPOSITION, [input.clone()]);
    let found = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    let (uu, h) = (
        from_expr(&found.args()[0]).unwrap(),
        from_expr(&found.args()[1]).unwrap(),
    );
    assert!(hermite_form(&h), "{call:?}: H {h:?}");
    assert_eq!(product(&uu, &a), h, "{call:?}: U A");
    assert!(unimodular(&uu), "{call:?}: U {uu:?}");
    let reversed: Matrix = a.iter().rev().cloned().collect();
    let again = ctx
        .eval(&Expr::call(
            Symbol::HERMITE_DECOMPOSITION,
            [to_expr(&reversed)],
        ))
        .unwrap();
    assert_eq!(
        from_expr(&again.args()[1]).unwrap(),
        h,
        "{call:?}: H of the rows reversed"
    );
    // Smith.
    let call = Expr::call(Symbol::SMITH_DECOMPOSITION, [input.clone()]);
    let found = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    let parts: Vec<Matrix> = found.args().iter().map(|p| from_expr(p).unwrap()).collect();
    let (uu, s, v) = (&parts[0], &parts[1], &parts[2]);
    assert_eq!(product(&product(uu, &a), v), *s, "{call:?}: U A V");
    assert!(unimodular(uu) && unimodular(v), "{call:?}: U, V");
    let diagonal: Vec<Integer> = (0..m.min(n)).map(|i| s[i][i].clone()).collect();
    for (i, row) in s.iter().enumerate() {
        for (j, x) in row.iter().enumerate() {
            assert!(i == j || x.is_zero(), "{call:?}: S {s:?}");
        }
    }
    let mut running = Integer::ONE;
    for (k, d) in diagonal.iter().enumerate() {
        assert!(!d.is_negative(), "{call:?}: S {s:?}");
        if k > 0 && !diagonal[k - 1].is_zero() {
            assert!(
                d.div_rem_euclid(&diagonal[k - 1]).1.is_zero(),
                "{call:?}: S {s:?} breaks the chain"
            );
        } else if k > 0 {
            assert!(d.is_zero(), "{call:?}: S {s:?}");
        }
        running = &running * d;
        assert_eq!(
            running,
            minors_gcd(&a, k + 1),
            "{call:?}: d₁⋯d_{} {s:?}",
            k + 1
        );
    }
    // LLL.
    let call = Expr::call(Symbol::LATTICE_REDUCE, [input]);
    let found = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    let reduced = from_expr(&found).unwrap();
    assert_eq!(
        lattice(&ctx, &reduced),
        lattice(&ctx, &a),
        "{call:?}: another lattice"
    );
    let rows: Vec<Vec<Rational>> = reduced
        .iter()
        .map(|r| r.iter().cloned().map(Rational::from).collect())
        .collect();
    let inner = |x: &[Rational], y: &[Rational]| {
        x.iter()
            .zip(y)
            .fold(Rational::ZERO, |s, (p, q)| &s + &(p * q))
    };
    let mut star: Vec<Vec<Rational>> = Vec::new();
    let mut norms: Vec<Rational> = Vec::new();
    let half = Rational::new(Integer::ONE, Integer::from(2));
    let delta = Rational::new(Integer::from(3), Integer::from(4));
    for (i, row) in rows.iter().enumerate() {
        let mut w = row.clone();
        let mut last_mu = Rational::ZERO;
        for j in 0..i {
            let mu = &inner(row, &star[j]) / &norms[j];
            assert!(mu.abs() <= half, "{call:?}: {reduced:?} not size-reduced");
            for (x, y) in w.iter_mut().zip(&star[j]) {
                *x = &*x - &(&mu * y);
            }
            last_mu = mu;
        }
        let b = inner(&w, &w);
        if i > 0 {
            let bound = &(&delta - &(&last_mu * &last_mu)) * &norms[i - 1];
            assert!(
                b >= bound,
                "{call:?}: {reduced:?} breaks Lovász's condition at {i}"
            );
        }
        assert!(!b.is_zero(), "{call:?}: {reduced:?} has dependent rows");
        norms.push(b);
        star.push(w);
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
