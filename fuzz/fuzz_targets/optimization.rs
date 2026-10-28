//! Optimization (CAS-DOC-38 §5, CAS-419), one check per input:
//! - a random linear program in 2 or 3 variables, each in `[−4, 4]` so
//!   that it is bounded, with up to 4 constraints of small integer
//!   coefficients, minimized or maximized: `LPSolve`'s optimum is the best
//!   vertex, found by solving every `n` of the constraints' and ranges'
//!   boundaries exactly, and its point is feasible with that value; `{}`
//!   when no vertex is feasible;
//! - the same with integer variables, against every integer point of the
//!   box;
//! - `Minimize` and `Maximize` of a random polynomial of degree at most 4
//!   on `[a, b]`: no value on a grid lies below the minimum or above the
//!   maximum.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// The box's half width.
const BOX: i64 = 4;

/// `A x = b` for a square `A`, exactly; `None` when singular.
fn solve(mut a: Vec<Vec<Rational>>, mut b: Vec<Rational>) -> Option<Vec<Rational>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, pivot);
        b.swap(col, pivot);
        for r in 0..n {
            if r == col || a[r][col].is_zero() {
                continue;
            }
            let factor = &a[r][col] / &a[col][col];
            for c in col..n {
                let delta = &factor * &a[col][c];
                a[r][c] = &a[r][c] - &delta;
            }
            let delta = &factor * &b[col];
            b[r] = &b[r] - &delta;
        }
    }
    Some((0..n).map(|i| &b[i] / &a[i][i]).collect())
}

/// Whether `x` meets every row `a·x ≤ b`.
fn feasible(rows: &[(Vec<Rational>, Rational)], x: &[Rational]) -> bool {
    rows.iter().all(|(a, b)| {
        let lhs = a
            .iter()
            .zip(x)
            .fold(Rational::ZERO, |s, (c, v)| &s + &(c * v));
        lhs <= *b
    })
}

fn dot(c: &[Rational], x: &[Rational]) -> Rational {
    c.iter()
        .zip(x)
        .fold(Rational::ZERO, |s, (a, v)| &s + &(a * v))
}

/// The best of `value` and `best`.
fn better(value: &Rational, best: Option<&Rational>, maximum: bool) -> bool {
    match best {
        None => true,
        Some(b) => {
            if maximum {
                value > b
            } else {
                value < b
            }
        }
    }
}

fn linear_program(ctx: &Context, u: &mut Unstructured<'_>, integer: bool) -> Result<()> {
    let n = u.int_in_range(2..=3usize)?;
    let names = ["x", "y", "z"];
    let xs: Vec<Expr> = names[..n].iter().map(|s| sym(s)).collect();
    let small = |u: &mut Unstructured<'_>| -> Result<i64> { u.int_in_range(-3..=3) };
    let objective: Vec<Rational> = (0..n)
        .map(|_| small(u).map(|c| q(c, 1)))
        .collect::<Result<_>>()?;
    let mut rows: Vec<(Vec<Rational>, Rational)> = Vec::new();
    for _ in 0..u.int_in_range(1..=4)? {
        let a: Vec<Rational> = (0..n)
            .map(|_| small(u).map(|c| q(c, 1)))
            .collect::<Result<_>>()?;
        let b = q(u.int_in_range(-6..=10)?, u.int_in_range(1..=2)?);
        rows.push((a, b));
    }
    let maximum: bool = u.arbitrary()?;
    // The boundaries: the constraints and the box.
    let mut boundaries = rows.clone();
    for j in 0..n {
        let mut e = vec![Rational::ZERO; n];
        e[j] = Rational::ONE;
        boundaries.push((e.clone(), q(BOX, 1)));
        e[j] = -Rational::ONE;
        boundaries.push((e, q(BOX, 1)));
    }
    let mut best: Option<(Vec<Rational>, Rational)> = None;
    if integer {
        let mut point = vec![-BOX; n];
        loop {
            let x: Vec<Rational> = point.iter().map(|&v| q(v, 1)).collect();
            if feasible(&rows, &x) {
                let value = dot(&objective, &x);
                if better(&value, best.as_ref().map(|(_, b)| b), maximum) {
                    best = Some((x, value));
                }
            }
            // The next point of the box.
            let Some(k) = (0..n).find(|&k| point[k] < BOX) else {
                break;
            };
            point[k] += 1;
            for v in &mut point[..k] {
                *v = -BOX;
            }
        }
    } else {
        let m = boundaries.len();
        let mut chosen = vec![0usize; n];
        // Every n of the m boundaries, as increasing indices.
        fn next(chosen: &mut [usize], m: usize) -> bool {
            let n = chosen.len();
            for i in (0..n).rev() {
                if chosen[i] < m - n + i {
                    chosen[i] += 1;
                    for k in i + 1..n {
                        chosen[k] = chosen[k - 1] + 1;
                    }
                    return true;
                }
            }
            false
        }
        for (i, c) in chosen.iter_mut().enumerate() {
            *c = i;
        }
        loop {
            let a: Vec<Vec<Rational>> = chosen.iter().map(|&i| boundaries[i].0.clone()).collect();
            let b: Vec<Rational> = chosen.iter().map(|&i| boundaries[i].1.clone()).collect();
            if let Some(x) = solve(a, b)
                && feasible(&boundaries, &x)
            {
                let value = dot(&objective, &x);
                if better(&value, best.as_ref().map(|(_, b)| b), maximum) {
                    best = Some((x, value));
                }
            }
            if !next(&mut chosen, m) {
                break;
            }
        }
    }
    // The call.
    let linear = |a: &[Rational]| -> Expr {
        Expr::call(
            Symbol::PLUS,
            a.iter()
                .zip(&xs)
                .map(|(c, x)| Expr::from(c.clone()) * x.clone()),
        )
    };
    let mut args = vec![
        linear(&objective),
        list(
            rows.iter()
                .map(|(a, b)| Expr::call(Symbol::LESS_EQUAL, [linear(a), Expr::from(b.clone())]))
                .collect(),
        ),
    ];
    for x in &xs {
        args.push(Expr::call(
            Symbol::EQUATION,
            [
                x.clone(),
                Expr::call(Symbol::RANGE, [Expr::from(-BOX), Expr::from(BOX)]),
            ],
        ));
    }
    if integer {
        args.push(Expr::call(
            Symbol::EQUATION,
            [Expr::from(Symbol::ASSUME), sym("integer")],
        ));
    }
    if maximum {
        args.push(Expr::from(Symbol::MAXIMIZE));
    }
    let call = Expr::call(Symbol::LP_SOLVE, args);
    let found = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    if found.head() == Some(&Expr::from(Symbol::LP_SOLVE)) {
        // Past its steps.
        return Ok(());
    }
    match best {
        None => assert_eq!(found, list(Vec::new()), "{call:?}: no point is feasible"),
        Some((_, value)) => {
            let [optimum, point] = found.args() else {
                panic!("{call:?}: {found:?}, expected {value:?}");
            };
            assert_eq!(
                optimum.to_rational(),
                Some(value.clone()),
                "{call:?}: {found:?}"
            );
            // By name: a variable with no coefficient in the objective comes
            // later in the answer.
            let x: Vec<Rational> = xs
                .iter()
                .map(|v| {
                    let eq = point
                        .args()
                        .iter()
                        .find(|eq| eq.args()[0] == *v)
                        .unwrap_or_else(|| panic!("{call:?}: {found:?} has no {v:?}"));
                    eq.args()[1].to_rational().expect("rational")
                })
                .collect();
            assert!(
                feasible(&boundaries, &x),
                "{call:?}: {found:?} is not feasible"
            );
            assert_eq!(dot(&objective, &x), value, "{call:?}: {found:?}");
            if integer {
                assert!(x.iter().all(Rational::is_integer), "{call:?}: {found:?}");
            }
        }
    }
    Ok(())
}

fn polynomial_optimum(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let x = sym("x");
    let degree = u.int_in_range(1..=4i64)?;
    let terms: Vec<Expr> = (0..=degree)
        .map(|k| {
            u.int_in_range(-4..=4i64)
                .map(|c| Expr::from(c) * x.clone().pow(k))
        })
        .collect::<Result<_>>()?;
    let f = Expr::call(Symbol::PLUS, terms);
    let a = q(u.int_in_range(-8..=4)?, 2);
    let b = &a + &q(u.int_in_range(1..=8)?, 2);
    let range = Expr::call(
        Symbol::EQUATION,
        [
            x.clone(),
            Expr::call(
                Symbol::RANGE,
                [Expr::from(a.clone()), Expr::from(b.clone())],
            ),
        ],
    );
    for (head, relation) in [
        (Symbol::MINIMIZE, Symbol::LESS),
        (Symbol::MAXIMIZE, Symbol::GREATER),
    ] {
        let call = Expr::call(head, [f.clone(), range.clone()]);
        let optimum = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
        if optimum.head() == Some(&Expr::from(head)) {
            continue;
        }
        for k in 0..=20 {
            let p = &a + &(&(&b - &a) * &q(k, 20));
            let value = f.substitute(&[(x.clone(), Expr::from(p.clone()))]);
            // f(p) < min, or f(p) > max, proved, is a mistake.
            let wrong = Expr::call(relation, [value - optimum.clone(), Expr::from(0)]);
            let told = ctx
                .eval(&wrong)
                .unwrap_or_else(|e| panic!("{wrong:?}: {e}"));
            assert!(
                told.as_symbol() != Some(Symbol::TRUE),
                "{call:?} is {optimum:?}, but f({p:?}) is beyond it"
            );
        }
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0u8..=2)? {
        0 => linear_program(&ctx, u, false),
        1 => linear_program(&ctx, u, true),
        _ => polynomial_optimum(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
