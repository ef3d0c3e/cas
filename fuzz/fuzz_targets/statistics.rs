//! The statistics builtins (CAS-DOC-22 §6, CAS-298) on random rational data
//! with frequencies, against the definitions computed here with exact
//! rationals:
//! - the mean, the population variance and covariance, the squares of the
//!   deviations (`stddev`, `stddevp`) and of the correlation;
//! - the median and quantiles, by the cumulative frequency of the sorted
//!   values;
//! - the linear regression, and the normal equations `Σ xʲ (y − p(x)) = 0`
//!   that the polynomial regression meets, or the call left when fewer
//!   distinct abscissas than coefficients make them singular;
//! - the exponential and power regressions of data on `y = b aˣ` and
//!   `y = b xᵃ`, which give back `a` and `b` exactly.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

/// At most this many values, small rationals.
const MAX_VALUES: usize = 8;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn value(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(-50..=50)?, u.int_in_range(1..=6)?))
}

fn list(values: &[Rational]) -> Expr {
    Expr::call(
        Symbol::LIST,
        values.iter().cloned().map(Expr::from).collect::<Vec<_>>(),
    )
}

/// The value of `head[args]`.
fn eval(ctx: &Context, head: Symbol, args: &[Expr]) -> Expr {
    let e = Expr::call(head, args.to_vec());
    ctx.eval(&e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The rational value of `head[args]`.
fn rational(ctx: &Context, head: Symbol, args: &[Expr]) -> Rational {
    let found = eval(ctx, head, args);
    found
        .to_rational()
        .unwrap_or_else(|| panic!("{head:?}{args:?} is {found:?}, not a rational"))
}

/// The rational value of the square of `head[args]`.
fn square(ctx: &Context, head: Symbol, args: &[Expr]) -> Rational {
    let found = eval(ctx, head, args);
    let squared = ctx.eval(&found.pow(2)).expect("a square");
    squared
        .to_rational()
        .unwrap_or_else(|| panic!("{head:?}{args:?} is {found:?}, squared {squared:?}"))
}

fn sum(values: impl Iterator<Item = Rational>) -> Rational {
    values.fold(Rational::ZERO, |s, v| &s + &v)
}

/// The values of positive frequency, ascending, with their cumulative
/// frequencies.
fn sorted(xs: &[Rational], fs: &[Rational]) -> Vec<(Rational, Rational)> {
    let mut pairs: Vec<(Rational, Rational)> = xs
        .iter()
        .cloned()
        .zip(fs.iter().cloned())
        .filter(|(_, f)| f.is_positive())
        .collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));
    let mut total = Rational::ZERO;
    for (_, f) in &mut pairs {
        total = &total + &*f;
        *f = total.clone();
    }
    pairs
}

/// The first value whose cumulative frequency reaches `target`, and
/// whether exactly.
fn reach(sorted: &[(Rational, Rational)], target: &Rational) -> (usize, bool) {
    let k = sorted
        .iter()
        .position(|(_, c)| c >= target)
        .expect("the total reaches every target");
    (k, sorted[k].1 == *target)
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(1..=MAX_VALUES)?;
    let xs = (0..n).map(|_| value(u)).collect::<Result<Vec<_>>>()?;
    let ys = (0..n).map(|_| value(u)).collect::<Result<Vec<_>>>()?;
    let fs = (0..n)
        .map(|_| Ok(q(u.int_in_range(0..=5)?, u.int_in_range(1..=3)?)))
        .collect::<Result<Vec<_>>>()?;
    let total = sum(fs.iter().cloned());
    let ctx = Context::new();
    let (x, y, f) = (list(&xs), list(&ys), list(&fs));
    let weighted = |g: &dyn Fn(usize) -> Rational| sum((0..n).map(|i| &fs[i] * &g(i)));

    if total.is_positive() {
        // Moments.
        let mean = &weighted(&|i| xs[i].clone()) / &total;
        let y_mean = &weighted(&|i| ys[i].clone()) / &total;
        let found = rational(&ctx, Symbol::MEAN, &[x.clone(), f.clone()]);
        assert_eq!(found, mean, "Mean[{xs:?}, {fs:?}]");
        let squares = weighted(&|i| {
            let d = &xs[i] - &mean;
            &d * &d
        });
        let variance = &squares / &total;
        let found = rational(&ctx, Symbol::POPULATION_VARIANCE, &[x.clone(), f.clone()]);
        assert_eq!(found, variance, "PopulationVariance[{xs:?}, {fs:?}]");
        let found = square(
            &ctx,
            Symbol::POPULATION_STANDARD_DEVIATION,
            &[x.clone(), f.clone()],
        );
        assert_eq!(
            found, variance,
            "PopulationStandardDeviation[{xs:?}, {fs:?}]"
        );
        if !total.is_one() {
            let sample = &squares / &(&total - &Rational::ONE);
            let found = square(&ctx, Symbol::STANDARD_DEVIATION, &[x.clone(), f.clone()]);
            assert_eq!(found, sample, "StandardDeviation[{xs:?}, {fs:?}]");
        }
        // Two variables.
        let sxy = weighted(&|i| &(&xs[i] - &mean) * &(&ys[i] - &y_mean));
        let syy = weighted(&|i| {
            let d = &ys[i] - &y_mean;
            &d * &d
        });
        let found = rational(
            &ctx,
            Symbol::POPULATION_COVARIANCE,
            &[x.clone(), y.clone(), f.clone()],
        );
        assert_eq!(
            found,
            &sxy / &total,
            "PopulationCovariance[{xs:?}, {ys:?}, {fs:?}]"
        );
        let r = eval(
            &ctx,
            Symbol::CORRELATION,
            &[x.clone(), y.clone(), f.clone()],
        );
        if squares.is_zero() || syy.is_zero() {
            assert_eq!(
                r,
                Expr::from(Symbol::INDETERMINATE),
                "Correlation, constant"
            );
        } else {
            let r2 = ctx.eval(&r.pow(2)).expect("a square").to_rational();
            assert_eq!(
                r2,
                Some(&(&sxy * &sxy) / &(&squares * &syy)),
                "Correlation² {r:?}"
            );
        }
        // Order statistics.
        let order = sorted(&xs, &fs);
        let (k, exactly) = reach(&order, &(&total / &q(2, 1)));
        let median = if exactly && k + 1 < order.len() {
            &(&order[k].0 + &order[k + 1].0) / &q(2, 1)
        } else {
            order[k].0.clone()
        };
        let found = rational(&ctx, Symbol::MEDIAN, &[x.clone(), f.clone()]);
        assert_eq!(found, median, "Median[{xs:?}, {fs:?}]");
        let p = q(u.int_in_range(0..=12)?, 12);
        let quantile = order[reach(&order, &(&p * &total)).0].0.clone();
        let found = rational(
            &ctx,
            Symbol::QUANTILE,
            &[x.clone(), f.clone(), Expr::from(p.clone())],
        );
        assert_eq!(found, quantile, "Quantile[{xs:?}, {fs:?}, {p:?}]");
    }

    // The linear regression, unweighted.
    let count = q(i64::try_from(n).expect("small"), 1);
    let x_mean = &sum(xs.iter().cloned()) / &count;
    let y_mean = &sum(ys.iter().cloned()) / &count;
    let sxx = sum(xs.iter().map(|v| {
        let d = v - &x_mean;
        &d * &d
    }));
    let sxy = sum(xs
        .iter()
        .zip(&ys)
        .map(|(a, b)| &(a - &x_mean) * &(b - &y_mean)));
    let line = eval(&ctx, Symbol::LINEAR_REGRESSION, &[x.clone(), y.clone()]);
    if sxx.is_zero() {
        let undefined = Expr::from(Symbol::INDETERMINATE);
        assert_eq!(
            line.args(),
            [undefined.clone(), undefined],
            "LinearRegression, constant"
        );
    } else {
        let a = &sxy / &sxx;
        let b = &y_mean - &(&a * &x_mean);
        let found: Vec<Option<Rational>> = line.args().iter().map(Expr::to_rational).collect();
        assert_eq!(
            found,
            [Some(a), Some(b)],
            "LinearRegression[{xs:?}, {ys:?}]"
        );
    }

    // The polynomial regression of a degree up to 3.
    let degree = u.int_in_range(0..=3usize.min(n - 1))?;
    let mut distinct = xs.clone();
    distinct.sort();
    distinct.dedup();
    let fit = eval(
        &ctx,
        Symbol::POLYNOMIAL_REGRESSION,
        &[
            x.clone(),
            y.clone(),
            Expr::from(i64::try_from(degree).expect("small")),
        ],
    );
    if distinct.len() <= degree {
        assert_eq!(
            fit.head(),
            Some(&Expr::from(Symbol::POLYNOMIAL_REGRESSION)),
            "singular"
        );
    } else {
        let coefficients: Vec<Rational> = fit
            .args()
            .iter()
            .map(|c| c.to_rational().unwrap_or_else(|| panic!("{fit:?}")))
            .collect();
        assert_eq!(coefficients.len(), degree + 1, "{fit:?}");
        let p = |v: &Rational| {
            coefficients
                .iter()
                .fold(Rational::ZERO, |acc, c| &(&acc * v) + c)
        };
        for j in 0..=degree {
            let residual = sum(xs.iter().zip(&ys).map(|(v, w)| {
                let power = (0..j).fold(Rational::ONE, |acc, _| &acc * v);
                &power * &(w - &p(v))
            }));
            assert!(residual.is_zero(), "normal equation {j} of {fit:?}");
        }
    }

    // Exact models: y = b aˣ at integer abscissas, y = b xᵃ at positive
    // ones.
    let a = q(u.int_in_range(1..=12)?, u.int_in_range(1..=12)?);
    let b = q(u.int_in_range(1..=12)?, u.int_in_range(1..=12)?);
    let m = u.int_in_range(2..=4usize)?;
    let abscissas: Vec<Rational> = (0..m)
        .map(|i| q(i64::try_from(i).expect("small"), 1))
        .collect();
    let pow = |base: &Rational, k: usize| (0..k).fold(Rational::ONE, |acc, _| &acc * base);
    let ordinates: Vec<Rational> = (0..m).map(|i| &b * &pow(&a, i)).collect();
    let found = eval(
        &ctx,
        Symbol::EXPONENTIAL_REGRESSION,
        &[list(&abscissas), list(&ordinates)],
    );
    let exact: Vec<Option<Rational>> = found.args().iter().map(Expr::to_rational).collect();
    assert_eq!(
        exact,
        [Some(a.clone()), Some(b.clone())],
        "ExponentialRegression: {found:?}"
    );
    let exponent = u.int_in_range(0..=3usize)?;
    let positive: Vec<Rational> = (1..=m)
        .map(|i| q(i64::try_from(i).expect("small"), 1))
        .collect();
    let ordinates: Vec<Rational> = positive.iter().map(|v| &b * &pow(v, exponent)).collect();
    let found = eval(
        &ctx,
        Symbol::POWER_REGRESSION,
        &[list(&positive), list(&ordinates)],
    );
    let exact: Vec<Option<Rational>> = found.args().iter().map(Expr::to_rational).collect();
    let k = q(i64::try_from(exponent).expect("small"), 1);
    assert_eq!(exact, [Some(k), Some(b)], "PowerRegression: {found:?}");
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
