//! The hypothesis tests (CAS-DOC-52, CAS-509) on random small integer data:
//! - `WilcoxonDistribution[n]` and `[m, n]`: their lengths, a sum of 1,
//!   symmetric, and the rank-sum law the same for `[n, m]`;
//! - `WilcoxonStatistic`: `W⁺(l, μ) + W⁺(−l, −μ) = k(k+1)/2` over the `k`
//!   nonzero differences, and the two rank sums adding up to `N(N+1)/2`;
//! - the tests answering 0 or 1, never rejecting at a smaller level what
//!   they keep at a larger one;
//! - `ChiSquareTest` at least 0, and 0 for counts proportional to the
//!   probabilities;
//! - `KolmogorovSmirnovTest`: `D` in `[0, 1]`, the same for the samples
//!   swapped, 0 for a sample against itself;
//! - `KolmogorovCDF` in `[0, 1]` and not decreasing.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

fn list(values: &[i64]) -> Expr {
    Expr::call(
        Symbol::LIST,
        values.iter().map(|&x| Expr::from(x)).collect::<Vec<_>>(),
    )
}

fn sample(u: &mut Unstructured<'_>, least: usize) -> Result<Vec<i64>> {
    let n = u.int_in_range(least..=7)?;
    (0..n).map(|_| u.int_in_range(-6..=6)).collect()
}

fn rational(e: &Expr) -> Rational {
    e.to_rational()
        .unwrap_or_else(|| panic!("not rational: {e:?}"))
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// The 0 or 1 of a test.
fn answer(ctx: &Context, e: &Expr) -> i64 {
    let found = eval(ctx, e);
    match found.as_i64() {
        Some(a @ (0 | 1)) => a,
        _ => panic!("{e:?}: {found:?}"),
    }
}

/// A test at the levels 1/20 and 1/10: 0 or 1, the smaller level keeping
/// what the larger keeps.
fn levels(ctx: &Context, head: &str, args: &[Expr]) {
    let at = |alpha: Rational| {
        let mut a = args.to_vec();
        a.push(Expr::string("!="));
        a.push(Expr::from(alpha));
        answer(ctx, &call(head, a))
    };
    let (strict, loose) = (at(q(1, 20)), at(q(1, 10)));
    assert!(
        strict >= loose,
        "{head}{args:?}: {strict} at 5 %, {loose} at 10 %"
    );
}

fn laws(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(0..=12)?;
    let w = eval(ctx, &call("WilcoxonDistribution", vec![Expr::from(n)]));
    let p: Vec<Rational> = w.args().iter().map(rational).collect();
    let top = usize::try_from(n * (n + 1) / 2).expect("small");
    assert_eq!(p.len(), top + 1, "{w:?}");
    assert_eq!(p.iter().fold(Rational::ZERO, |s, x| &s + x), Rational::ONE);
    assert!(p.iter().eq(p.iter().rev()), "{w:?}");
    let m = u.int_in_range(0..=8)?;
    let k = u.int_in_range(0..=8)?;
    let a = eval(
        ctx,
        &call("WilcoxonDistribution", vec![Expr::from(m), Expr::from(k)]),
    );
    let b = eval(
        ctx,
        &call("WilcoxonDistribution", vec![Expr::from(k), Expr::from(m)]),
    );
    assert_eq!(a, b);
    let p: Vec<Rational> = a.args().iter().map(rational).collect();
    assert_eq!(p.len(), usize::try_from(m * k + 1).expect("small"), "{a:?}");
    assert_eq!(p.iter().fold(Rational::ZERO, |s, x| &s + x), Rational::ONE);
    assert!(p.iter().eq(p.iter().rev()), "{a:?}");
    Ok(())
}

fn statistics(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let l = sample(u, 1)?;
    let mu = u.int_in_range(-6..=6)?;
    let k = i64::try_from(l.iter().filter(|&&x| x != mu).count()).expect("small");
    let w = |l: &[i64], mu: i64| {
        rational(&eval(
            ctx,
            &call("WilcoxonStatistic", vec![list(l), Expr::from(mu)]),
        ))
    };
    let negated: Vec<i64> = l.iter().map(|x| -x).collect();
    assert_eq!(
        &w(&l, mu) + &w(&negated, -mu),
        q(k * (k + 1), 2),
        "{l:?} about {mu}"
    );
    levels(ctx, "WilcoxonTest", &[list(&l), Expr::from(mu)]);
    let other = sample(u, 1)?;
    let r = |a: &[i64], b: &[i64]| {
        rational(&eval(
            ctx,
            &call("WilcoxonStatistic", vec![list(a), list(b)]),
        ))
    };
    let total = i64::try_from(l.len() + other.len()).expect("small");
    assert_eq!(
        &r(&l, &other) + &r(&other, &l),
        q(total * (total + 1), 2),
        "{l:?}, {other:?}"
    );
    levels(ctx, "WilcoxonTest", &[list(&l), list(&other)]);
    // Kolmogorov–Smirnov.
    let d = |a: &[i64], b: &[i64]| {
        let found = eval(ctx, &call("KolmogorovSmirnovTest", vec![list(a), list(b)]));
        rational(&found.args()[1])
    };
    let gap = d(&l, &other);
    assert!(
        gap >= Rational::ZERO && gap <= Rational::ONE,
        "{l:?}, {other:?}: {gap:?}"
    );
    assert_eq!(gap, d(&other, &l));
    assert_eq!(d(&l, &l), Rational::ZERO);
    // Z and T tests.
    if l.len() >= 3 {
        let sigma = Expr::from(u.int_in_range(1..=5)?);
        levels(ctx, "ZTest", &[list(&l), Expr::from(mu), sigma.clone()]);
        if l.iter().any(|&x| x != l[0]) {
            levels(ctx, "TTest", &[list(&l), Expr::from(mu), sigma]);
        }
    }
    Ok(())
}

fn chi_square(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let counts: Vec<i64> = (0..u.int_in_range(1..=6)?)
        .map(|_| u.int_in_range(0..=20))
        .collect::<Result<_>>()?;
    if counts.iter().all(|&c| c == 0) {
        return Ok(());
    }
    let chi = rational(&eval(ctx, &call("ChiSquareTest", vec![list(&counts)])));
    assert!(chi >= Rational::ZERO, "{counts:?}: {chi:?}");
    // Proportional counts: 0.
    let c = u.int_in_range(1..=5)?;
    let total: i64 = counts.iter().sum();
    let probabilities = Expr::call(
        Symbol::LIST,
        counts
            .iter()
            .map(|&x| Expr::from(q(x, total)))
            .collect::<Vec<_>>(),
    );
    let scaled: Vec<i64> = counts.iter().map(|x| x * c).collect();
    if counts.iter().all(|&x| x > 0) {
        let chi = eval(
            ctx,
            &call("ChiSquareTest", vec![list(&scaled), probabilities]),
        );
        assert_eq!(chi, Expr::from(0), "{scaled:?}");
    }
    Ok(())
}

fn kolmogorov(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let a = q(u.int_in_range(1..=400)?, 100);
    let b = &a + &q(u.int_in_range(1..=100)?, 100);
    let value = |x: &Rational| {
        let e = Expr::call(
            Symbol::N,
            [
                call("KolmogorovCDF", vec![Expr::from(x.clone())]),
                Expr::from(30),
            ],
        );
        match eval(ctx, &e).view() {
            View::Ball(k) => k.clone(),
            _ => panic!("K({x:?})"),
        }
    };
    let (ka, kb) = (value(&a), value(&b));
    for k in [&ka, &kb] {
        assert!(
            k.mid().to_f64() >= -1e-20 && k.mid().to_f64() <= 1.0 + 1e-20,
            "{k:?}"
        );
    }
    assert_ne!(
        ka.cmp(&kb),
        Some(std::cmp::Ordering::Greater),
        "K({a:?}) = {ka:?}, K({b:?}) = {kb:?}"
    );
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0..=3)? {
        0 => laws(&ctx, u),
        1 => statistics(&ctx, u),
        2 => chi_square(&ctx, u),
        _ => kolmogorov(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let _ = run(&mut Unstructured::new(data));
});
