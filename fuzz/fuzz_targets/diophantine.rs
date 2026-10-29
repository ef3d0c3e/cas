//! Number theory II (CAS-DOC-39 §6, CAS-426), one check per input:
//! - `QuotientRemainder` and `SymmetricMod`: their defining relations;
//! - `MultiplicativeOrder` and `PrimitiveRoot` against brute force below
//!   2000; `PseudoprimeQ` against `PrimeQ`;
//! - `Cyclotomic[n]`: of degree `φ(n)`, and `Φ_n(2)` divides `2^n − 1`;
//! - `IntegerSolve` of random linear equations and systems: the answer is
//!   `x₀ + Σ tₖ cₖ`, each instance satisfies them, the `cₖ` are independent
//!   and as many as the solutions' dimension, and every solution in a box
//!   is an instance;
//! - `x y = n` and `x² + y² = n` against enumeration; Pell's `x² − D y² =
//!   1`: the fundamental solution, minimal where enumeration reaches it,
//!   and the next one;
//! - `FromContinuedFraction[ContinuedFraction[x]] = x` for rationals and
//!   quadratic irrationals, also with the period started one quotient late;
//!   `ContinuedFraction[x, n]` by balls against the exact quotients.
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

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(head: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(head, args)
}

fn rule(from: Expr, to: Expr) -> Expr {
    Expr::call(Symbol::RULE, [from, to])
}

fn stays(found: &Expr, head: Symbol) -> bool {
    found.head() == Some(&Expr::from(head))
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

fn order(a: i64, n: i64) -> i64 {
    if gcd(a, n) != 1 {
        return 0;
    }
    let mut x = a.rem_euclid(n);
    let mut k = 1;
    while x != 1 % n {
        x = x * a.rem_euclid(n) % n;
        k += 1;
    }
    k
}

fn phi(n: i64) -> i64 {
    (1..=n).filter(|&k| gcd(k, n) == 1).count() as i64
}

fn residues(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let a = u.int_in_range(-1_000_000..=1_000_000i64)?;
    let b = u.int_in_range(-1000..=1000i64)?;
    if b != 0 {
        let found = eval(
            ctx,
            &call(
                Symbol::QUOTIENT_REMAINDER,
                vec![Expr::from(a), Expr::from(b)],
            ),
        );
        let [qq, r] = found.args() else {
            panic!("QuotientRemainder[{a}, {b}]: {found:?}");
        };
        let (qq, r) = (qq.as_i64().unwrap(), r.as_i64().unwrap());
        assert!(
            a == b * qq + r && 0 <= r && r < b.abs(),
            "QuotientRemainder[{a}, {b}]: {found:?}"
        );
        let s = eval(
            ctx,
            &call(Symbol::SYMMETRIC_MOD, vec![Expr::from(a), Expr::from(b)]),
        );
        let s = s.as_i64().unwrap();
        let n = b.abs();
        assert!(
            (a - s) % n == 0 && -n < 2 * s && 2 * s <= n,
            "SymmetricMod[{a}, {b}] is {s}"
        );
    }
    let n = u.int_in_range(1..=2000i64)?;
    let k = u.int_in_range(-50..=50i64)?;
    let found = eval(
        ctx,
        &call(
            Symbol::MULTIPLICATIVE_ORDER,
            vec![Expr::from(k), Expr::from(n)],
        ),
    );
    assert_eq!(
        found,
        Expr::from(order(k, n)),
        "MultiplicativeOrder[{k}, {n}]"
    );
    let root = eval(ctx, &call(Symbol::PRIMITIVE_ROOT, vec![Expr::from(n)]));
    let expected = if n == 1 {
        Expr::from(0)
    } else {
        let f = phi(n);
        match (1..n).find(|&g| order(g, n) == f) {
            Some(g) => Expr::from(g),
            None => Expr::from(Symbol::INDETERMINATE),
        }
    };
    assert_eq!(root, expected, "PrimitiveRoot[{n}]");
    let pseudo = eval(ctx, &call(Symbol::PSEUDOPRIME_Q, vec![Expr::from(n)]));
    let prime = eval(ctx, &call(Symbol::PRIME_Q, vec![Expr::from(n)]));
    assert_eq!(
        pseudo == Expr::from(2),
        prime == Expr::from(Symbol::TRUE),
        "PseudoprimeQ[{n}] is {pseudo:?}, PrimeQ {prime:?}"
    );
    let m = u.int_in_range(1..=300i64)?;
    let coefficients = eval(ctx, &call(Symbol::CYCLOTOMIC, vec![Expr::from(m)]));
    assert_eq!(
        coefficients.args().len() as i64 - 1,
        phi(m),
        "Cyclotomic[{m}]"
    );
    let x = sym("x");
    let at_two = eval(
        ctx,
        &call(
            Symbol::REPLACE_ALL,
            vec![
                call(Symbol::CYCLOTOMIC, vec![Expr::from(m), x.clone()]),
                rule(x, Expr::from(2)),
            ],
        ),
    );
    let at_two = at_two.to_integer().unwrap();
    let power = (Integer::from(2).pow(u64::try_from(m).unwrap())) - Integer::ONE;
    assert!(
        power.div_rem_euclid(&at_two).1.is_zero(),
        "Φ_{m}(2) = {at_two:?} does not divide 2^{m} − 1"
    );
    Ok(())
}

/// Rows reduced by Gauss and Jordan over the rationals, their pivots
/// `(row, column)` among the first `columns` columns.
fn eliminate(m: &mut [Vec<Rational>], columns: usize) -> Vec<(usize, usize)> {
    let mut pivots = Vec::new();
    for col in 0..columns {
        let row = pivots.len();
        let Some(r) = (row..m.len()).find(|&r| !m[r][col].is_zero()) else {
            continue;
        };
        m.swap(row, r);
        for r in 0..m.len() {
            if r != row && !m[r][col].is_zero() {
                let f = &m[r][col] / &m[row][col];
                for c in 0..m[r].len() {
                    let v = &m[r][c] - &(&f * &m[row][c]);
                    m[r][c] = v;
                }
            }
        }
        pivots.push((row, col));
    }
    pivots
}

/// The rank of integer rows.
fn rank(rows: &[Vec<i64>]) -> usize {
    let width = rows.first().map_or(0, Vec::len);
    let mut m: Vec<Vec<Rational>> = rows
        .iter()
        .map(|r| r.iter().map(|&v| q(v, 1)).collect())
        .collect();
    eliminate(&mut m, width).len()
}

/// Whether `p − x₀` is an integer combination of the independent
/// `columns`.
fn reachable(x0: &[i64], columns: &[Vec<i64>], p: &[i64]) -> bool {
    let k = columns.len();
    let mut m: Vec<Vec<Rational>> = (0..x0.len())
        .map(|i| {
            let mut row: Vec<Rational> = columns.iter().map(|c| q(c[i], 1)).collect();
            row.push(q(p[i] - x0[i], 1));
            row
        })
        .collect();
    let pivots = eliminate(&mut m, k);
    m[pivots.len()..].iter().all(|row| row[k].is_zero())
        && pivots
            .iter()
            .all(|&(r, c)| (&m[r][k] / &m[r][c]).is_integer())
}

fn linear(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(2..=3usize)?;
    let m = u.int_in_range(1..=2usize)?;
    let names = ["x", "y", "z"];
    let xs: Vec<Expr> = names[..n].iter().map(|s| sym(s)).collect();
    let mut rows: Vec<(Vec<i64>, i64)> = Vec::new();
    for _ in 0..m {
        let a: Vec<i64> = (0..n)
            .map(|_| u.int_in_range(-6..=6i64))
            .collect::<Result<_>>()?;
        rows.push((a, u.int_in_range(-20..=20i64)?));
    }
    let equations: Vec<Expr> = rows
        .iter()
        .map(|(a, b)| {
            let lhs = Expr::call(
                Symbol::PLUS,
                a.iter().zip(&xs).map(|(c, x)| Expr::from(*c) * x.clone()),
            );
            Expr::call(Symbol::EQUATION, [lhs, Expr::from(*b)])
        })
        .collect();
    let input = call(Symbol::INTEGER_SOLVE, vec![list(equations), list(xs)]);
    let found = eval(ctx, &input);
    if stays(&found, Symbol::INTEGER_SOLVE) {
        return Ok(());
    }
    let holds = |values: &[i64]| {
        rows.iter()
            .all(|(a, b)| a.iter().zip(values).map(|(c, v)| c * v).sum::<i64>() == *b)
    };
    // The box's solutions.
    let mut in_box = Vec::new();
    let mut point = vec![-6i64; n];
    loop {
        if holds(&point) {
            in_box.push(point.clone());
        }
        let Some(k) = (0..n).find(|&k| point[k] < 6) else {
            break;
        };
        point[k] += 1;
        for v in &mut point[..k] {
            *v = -6;
        }
    }
    if found.args().is_empty() {
        assert!(
            in_box.is_empty(),
            "{input:?} is {{}}, but {in_box:?} solve it"
        );
        return Ok(());
    }
    // The answer at parameters t: x₀ at 0, x₀ + cₖ at the k-th unit.
    let parameters: Vec<Expr> = (0..n).map(|k| sym(&format!("n_{k}"))).collect();
    let count = parameters
        .iter()
        .take_while(|p| found.args().iter().any(|eq| eq.args()[1].contains(p)))
        .count();
    let at = |t: &[i64]| -> Vec<i64> {
        let rules = list(
            parameters
                .iter()
                .zip(t)
                .map(|(p, v)| rule(p.clone(), Expr::from(*v)))
                .collect(),
        );
        found
            .args()
            .iter()
            .map(|eq| {
                let v = eval(
                    ctx,
                    &call(
                        Symbol::REPLACE_ALL,
                        vec![eq.args()[1].clone(), rules.clone()],
                    ),
                );
                v.as_i64()
                    .unwrap_or_else(|| panic!("{input:?}: {found:?} gives {v:?}"))
            })
            .collect()
    };
    let x0 = at(&vec![0; n]);
    let columns: Vec<Vec<i64>> = (0..count)
        .map(|k| {
            let mut t = vec![0; n];
            t[k] = 1;
            at(&t).iter().zip(&x0).map(|(a, b)| a - b).collect()
        })
        .collect();
    for t in [[0, 0, 0], [1, -1, 2], [-2, 3, -1], [5, 0, 4]] {
        let values = at(&t[..n]);
        assert!(
            holds(&values),
            "{input:?}: {found:?} at {t:?} gives {values:?}"
        );
        let affine: Vec<i64> = (0..n)
            .map(|i| x0[i] + (0..count).map(|k| t[k] * columns[k][i]).sum::<i64>())
            .collect();
        assert_eq!(values, affine, "{input:?}: {found:?} is not affine");
    }
    let a: Vec<Vec<i64>> = rows.iter().map(|(a, _)| a.clone()).collect();
    assert_eq!(
        rank(&columns),
        count,
        "{input:?}: {found:?} has dependent parameters"
    );
    assert_eq!(
        count,
        n - rank(&a),
        "{input:?}: {found:?} has too few parameters"
    );
    for p in &in_box {
        assert!(
            reachable(&x0, &columns, p),
            "{input:?} is {found:?}, which misses {p:?}"
        );
    }
    Ok(())
}

/// `{x, y}` of an answer `{x = a, y = b}`.
fn pair(p: &Expr) -> (i64, i64) {
    (
        p.args()[0].args()[1].as_i64().unwrap(),
        p.args()[1].args()[1].as_i64().unwrap(),
    )
}

fn pell(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let (x, y) = (sym("x"), sym("y"));
    let d = u.int_in_range(2..=120i64)?;
    let lhs = x.clone().pow(2) - Expr::from(d) * y.clone().pow(2);
    let input = call(
        Symbol::INTEGER_SOLVE,
        vec![Expr::call(Symbol::EQUATION, [lhs, Expr::from(1)])],
    );
    let found = eval(ctx, &input);
    let square = (1..=11).any(|r| r * r == d);
    if square || stays(&found, Symbol::INTEGER_SOLVE) {
        assert!(square, "{input:?} stays");
        return Ok(());
    }
    let at = |k: i64| -> (Integer, Integer) {
        let values: Vec<Integer> = found
            .args()
            .iter()
            .map(|eq| {
                let e = call(
                    Symbol::EXPAND,
                    vec![call(
                        Symbol::REPLACE_ALL,
                        vec![eq.args()[1].clone(), rule(sym("n_0"), Expr::from(k))],
                    )],
                );
                let v = eval(ctx, &e);
                v.to_integer()
                    .unwrap_or_else(|| panic!("{input:?}: {found:?} gives {v:?} at {k}"))
            })
            .collect();
        let [a, b] = values.try_into().unwrap();
        (a, b)
    };
    let dd = Integer::from(d);
    let norm = |(a, b): &(Integer, Integer)| &(a * a) - &(&(&dd * b) * b);
    let (x1, y1) = at(1);
    assert!(
        x1.is_positive() && y1.is_positive() && norm(&(x1.clone(), y1.clone())) == Integer::ONE,
        "{input:?}: {found:?} gives ({x1:?}, {y1:?})"
    );
    // The fundamental solution is the least.
    if let Some(y1) = y1.to_i64().filter(|&v| v <= 20_000) {
        for b in 1..y1 {
            let square = Integer::from(1 + d * b * b);
            let (_, rest) = square.sqrt_rem();
            assert!(!rest.is_zero(), "{input:?}: y = {b} solves it before {y1}");
        }
    }
    // The next: (x₁ + y₁ √D)² = x₁² + D y₁² + 2 x₁ y₁ √D.
    let second = at(2);
    let expected = (
        &(&x1 * &x1) + &(&(&dd * &y1) * &y1),
        &(&Integer::from(2) * &x1) * &y1,
    );
    assert_eq!(second, expected, "{input:?}: {found:?} at 2");
    Ok(())
}

fn quadratic(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let (x, y) = (sym("x"), sym("y"));
    let kind = u.int_in_range(0..=2u8)?;
    if kind == 2 {
        return pell(ctx, u);
    }
    let n = u.int_in_range(1..=500i64)?;
    let product = kind == 0;
    let lhs = if product {
        x.clone() * y.clone()
    } else {
        x.clone().pow(2) + y.clone().pow(2)
    };
    let input = call(
        Symbol::INTEGER_SOLVE,
        vec![Expr::call(Symbol::EQUATION, [lhs, Expr::from(n)])],
    );
    let found = eval(ctx, &input);
    assert!(!stays(&found, Symbol::INTEGER_SOLVE), "{input:?} stays");
    let mut got: Vec<(i64, i64)> = found.args().iter().map(pair).collect();
    got.sort_unstable();
    let mut expected = Vec::new();
    for a in -n..=n {
        for b in -n..=n {
            if (product && a * b == n) || (!product && a * a + b * b == n) {
                expected.push((a, b));
            }
        }
    }
    assert_eq!(got, expected, "{input:?}: {found:?}");
    Ok(())
}

/// The first `n` quotients of an expansion with its period last.
fn unrolled(fraction: &Expr, n: usize) -> Vec<Expr> {
    let items = fraction.args();
    let (period, prefix) = items.split_last().unwrap();
    prefix
        .iter()
        .chain(period.args().iter().cycle())
        .take(n)
        .cloned()
        .collect()
}

fn fractions(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let rational: bool = u.arbitrary()?;
    let x = if rational {
        Expr::from(q(u.int_in_range(-1000..=1000)?, u.int_in_range(1..=500)?))
    } else {
        let d = u.int_in_range(2..=60i64)?;
        let a = q(u.int_in_range(-9..=9)?, u.int_in_range(1..=4)?);
        let b = q(
            u.int_in_range(1..=5)? * if u.arbitrary()? { -1 } else { 1 },
            u.int_in_range(1..=4)?,
        );
        Expr::from(a) + Expr::from(b) * Expr::from(d).pow(Expr::from(q(1, 2)))
    };
    let x = eval(ctx, &x);
    let fraction = eval(ctx, &call(Symbol::CONTINUED_FRACTION, vec![x.clone()]));
    let irrational = x.to_rational().is_none();
    assert!(
        !stays(&fraction, Symbol::CONTINUED_FRACTION),
        "ContinuedFraction[{x:?}] stays"
    );
    let back = |fraction: &Expr| {
        let back = eval(
            ctx,
            &call(Symbol::FROM_CONTINUED_FRACTION, vec![fraction.clone()]),
        );
        let difference = eval(ctx, &call(Symbol::EXPAND, vec![back.clone() - x.clone()]));
        assert_eq!(difference, Expr::from(0), "{x:?}: {fraction:?} is {back:?}");
    };
    back(&fraction);
    if !irrational {
        return Ok(());
    }
    // The period started one quotient late: {…, p₁, {p₂, …, p₁}}.
    let mut late = fraction.args().to_vec();
    let period = late.pop().unwrap();
    let mut rotated = period.args().to_vec();
    rotated.rotate_left(1);
    late.push(period.args()[0].clone());
    late.push(list(rotated));
    back(&list(late));
    // The numeric quotients, by balls.
    let n = u.int_in_range(1..=30usize)?;
    let numeric = eval(
        ctx,
        &call(
            Symbol::CONTINUED_FRACTION,
            vec![x.clone(), Expr::from(i64::try_from(n).unwrap())],
        ),
    );
    assert_eq!(
        numeric,
        list(unrolled(&fraction, n)),
        "ContinuedFraction[{x:?}, {n}] against {fraction:?}"
    );
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0u8..=3)? {
        0 => residues(&ctx, u),
        1 => linear(&ctx, u),
        2 => quadratic(&ctx, u),
        _ => fractions(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
