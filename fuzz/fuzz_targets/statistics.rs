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
//!   `y = b xᵃ`, which give back `a` and `b` exactly;
//! - frequencies against the counts, cumulated ones ending at 1, classes
//!   counting their values, the root mean square, moving averages
//!   (CAS-DOC-50 §4, CAS-493);
//! - `Cluster` and `KMeans` (CAS-630, CAS-646) on small points or strings with
//!   random options (`output=part`, `count`, `index`, `index=` names): the
//!   parts are the data grouped by the labels, their number the count, and
//!   the silhouette rows' `k` decrease with values in `[−1, 1]` or `DBL_MIN`.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
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
    second(&ctx, u, &xs)
}

/// Statistics II (CAS-493) on the values `xs`.
fn second(ctx: &Context, u: &mut Unstructured<'_>, xs: &[Rational]) -> Result<()> {
    let n = xs.len();
    let total = q(i64::try_from(n).expect("small"), 1);
    let x = list(xs);
    let count = |v: &Rational| xs.iter().filter(|w| *w == v).count();
    let frequencies = eval(ctx, Symbol::FREQUENCIES, &[x.clone()]);
    let mut seen = 0usize;
    for pair in frequencies.args() {
        let [v, f] = pair.args() else {
            panic!("Frequencies[{xs:?}] is {frequencies:?}");
        };
        let (v, f) = (
            v.to_rational().expect("a value"),
            f.to_rational().expect("a share"),
        );
        let k = count(&v);
        assert_eq!(
            f,
            &q(i64::try_from(k).expect("small"), 1) / &total,
            "Frequencies[{xs:?}]"
        );
        seen += k;
    }
    assert_eq!(seen, n, "Frequencies[{xs:?}] is {frequencies:?}");
    let cumulated = eval(ctx, Symbol::CUMULATED_FREQUENCIES, &[x.clone()]);
    let last = cumulated
        .args()
        .last()
        .and_then(|p| p.args().get(1))
        .and_then(Expr::to_rational);
    assert_eq!(
        last,
        Some(q(1, 1)),
        "CumulatedFrequencies[{xs:?}] is {cumulated:?}"
    );
    // Classes from below the least value.
    let least = xs.iter().min().expect("a value").clone();
    let start = &least - &q(u.int_in_range(0..=4)?, 2);
    let width = q(u.int_in_range(1..=6)?, 2);
    let classes = eval(
        ctx,
        Symbol::CLASSES,
        &[x.clone(), Expr::from(start), Expr::from(width)],
    );
    let mut counted = 0i64;
    for row in classes.args() {
        let [range, c] = row.args() else {
            panic!("Classes[{xs:?}] is {classes:?}");
        };
        let [lo, hi] = range.args() else {
            panic!("Classes[{xs:?}] is {classes:?}");
        };
        let (lo, hi) = (
            lo.to_rational().expect("an end"),
            hi.to_rational().expect("an end"),
        );
        let inside = xs.iter().filter(|v| **v >= lo && **v < hi).count();
        let c = c.as_i64().expect("a count");
        assert_eq!(
            c,
            i64::try_from(inside).expect("small"),
            "Classes[{xs:?}]: {row:?}"
        );
        counted += c;
    }
    assert_eq!(
        counted,
        i64::try_from(n).expect("small"),
        "Classes[{xs:?}] is {classes:?}"
    );
    // The root mean square, squared.
    let squares = xs.iter().fold(q(0, 1), |acc, v| &acc + &(v * v));
    assert_eq!(
        square(ctx, Symbol::ROOT_MEAN_SQUARE, &[x.clone()]),
        &squares / &total,
        "RootMeanSquare[{xs:?}]"
    );
    // Moving averages.
    let k = u.int_in_range(1..=n)?;
    let means = eval(
        ctx,
        Symbol::MOVING_AVERAGE,
        &[x, Expr::from(i64::try_from(k).expect("small"))],
    );
    let expected: Vec<Rational> = xs
        .windows(k)
        .map(|w| {
            &w.iter().fold(q(0, 1), |acc, v| &acc + v) / &q(i64::try_from(k).expect("small"), 1)
        })
        .collect();
    let found: Vec<Option<Rational>> = means.args().iter().map(Expr::to_rational).collect();
    assert_eq!(
        found,
        expected.into_iter().map(Some).collect::<Vec<_>>(),
        "MovingAverage[{xs:?}, {k}]"
    );
    Ok(())
}

/// `Cluster` or `KMeans` on up to 7 points (2-D rationals or strings of
/// `a`, `b`, `é`) with random options: the parts, the count and the labels
/// agree; the silhouette rows are well formed.
fn clusters(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let m = u.int_in_range(1..=7)?;
    let strings = u.arbitrary::<bool>()?;
    let mut items = Vec::with_capacity(m);
    for _ in 0..m {
        items.push(if strings {
            let n = u.int_in_range(0..=4)?;
            let mut w = String::new();
            for _ in 0..n {
                w.push(*u.choose(&['a', 'b', 'é'])?);
            }
            format!("{w:?}")
        } else {
            format!(
                "List[{}, {}]",
                u.int_in_range(-3..=3)?,
                u.int_in_range(-3..=3)?
            )
        });
    }
    let head = if u.arbitrary::<bool>()? {
        "Cluster"
    } else {
        "KMeans"
    };
    let mut options = Vec::new();
    for _ in 0..u.int_in_range(0..=3)? {
        options.push(match u.int_in_range(0..=5)? {
            0 => format!("Equation[Count, {}]", u.int_in_range(1..=4)?),
            1 => format!("Equation[CountGreater, {}]", u.int_in_range(0..=8)?),
            2 => format!(
                "Equation[Type, {:?}]",
                u.choose(&["single", "complete", "average", "weighted", "ward"])?
            ),
            3 => format!(
                "Equation[FirstPosition, {}]",
                u.choose(&[
                    "\"silhouette\"",
                    "\"dunn\"",
                    "\"mcclain-rao\"",
                    "\"all\"",
                    "True",
                    "False",
                    "List[\"silhouette\", \"dunn\"]",
                    "\"foo\"",
                    "List[]",
                    // kmeans' cluster_crit names (T5.c3).
                    "\"calinski-harabasz\"",
                    "\"banfeld-raftery\"",
                    "\"davies-bouldin\"",
                    "\"ray-turi\"",
                    "\"scott-symons\"",
                    "\"pbm\"",
                    "\"ratkowsky-lance\"",
                    "\"det\"",
                    "\"log-det\"",
                    "\"log-ss\"",
                    "\"ksq-detW\"",
                    "\"ball-hall\"",
                    "\"trace-W\"",
                    "\"trace-WiB\"",
                    "List[\"pbm\", \"trace-WiB\"]",
                ])?
            ),
            4 => "Equation[Count, Range[2, 4]]".to_owned(),
            _ => "Equation[CountLess, 2]".to_owned(),
        });
    }
    let call = |output: &str| {
        let mut args = vec![format!("List[{}]", items.join(", "))];
        args.extend(options.iter().cloned());
        args.push(format!("Equation[output, {output}]"));
        let text = format!("{head}[{}]", args.join(", "));
        let e = Expr::parse_full_form(&text).unwrap_or_else(|err| panic!("{text}: {err}"));
        let found = ctx.eval(&e).unwrap_or_else(|err| panic!("{text}: {err}"));
        (text, (found.head() != e.head()).then_some(found))
    };
    let (text, labels) = call("List");
    let (_, parts) = call("Part");
    let (_, count) = call("Count");
    let (index_text, index) = call("FirstPosition");
    // cluster's silhouette rows: k decreasing, values in [-1, 1] (kmeans'
    // indices go up in k, with any value, "not defined" or ComplexInfinity).
    if let (Some(rows), "Cluster") = (index, head) {
        let mut last = usize::MAX;
        for row in rows.args() {
            let [k, v] = row.args() else {
                panic!("{index_text}: {rows:?}");
            };
            let k = k.as_i64().and_then(|k| usize::try_from(k).ok());
            let v = match v.view() {
                View::Ball(b) => b.mid().to_f64(),
                _ => panic!("{index_text}: {rows:?}"),
            };
            assert!(
                k.is_some_and(|k| k < last) && (-1.0..=1.0).contains(&v),
                "{index_text}: {rows:?}"
            );
            last = k.unwrap_or(0);
        }
    }
    // kmeans draws its centres anew per call: only cluster's three answers
    // agree with one another.
    let (Some(labels), Some(parts), Some(count), "Cluster") = (labels, parts, count, head) else {
        return Ok(());
    };
    let data = ctx
        .eval(&Expr::parse_full_form(&format!("List[{}]", items.join(", "))).expect("data"))
        .expect("data");
    let k = count.as_i64().and_then(|k| usize::try_from(k).ok());
    assert_eq!(k, Some(parts.args().len()), "{text}: {count:?} {parts:?}");
    for (c, group) in parts.args().iter().enumerate() {
        let expected: Vec<&Expr> = labels
            .args()
            .iter()
            .zip(data.args())
            .filter(|(l, _)| l.as_i64().and_then(|l| usize::try_from(l).ok()) == Some(c))
            .map(|(_, x)| x)
            .collect();
        assert!(
            !expected.is_empty() && group.args().iter().eq(expected),
            "{text}: {labels:?} {parts:?}"
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
    let _ = clusters(&mut Unstructured::new(data));
});
