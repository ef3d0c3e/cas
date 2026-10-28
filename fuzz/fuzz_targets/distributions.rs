//! The probability distributions (CAS-DOC-23 §7, CAS-305) at random
//! rational parameters, one family per input:
//! - binomial, geometric and negative binomial against the definitions
//!   computed here with exact rationals: the masses of a bounded support
//!   sum to 1, the cumulative function sums the masses, an interval's
//!   probability is a difference, and the inverse is the least `k` with
//!   `F(k) ≥ q`;
//! - Poisson's, whose masses hold `e^(−λ)`, by the CAS itself: the
//!   cumulative function is the expanded sum of the masses, and the
//!   inverse `k` has `F(k − 1) < q ≤ F(k)`;
//! - uniform and exponential: the cumulative function takes the inverse
//!   back, exactly;
//! - normal, Student, χ² and Fisher, and `erfinv`: `F(F⁻¹(q))` within
//!   `10⁻¹²` of `q`, the cumulative function increasing and in `[0, 1]`.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn int(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("small"))
}

/// A probability `k/d`, `0 ≤ k ≤ d`.
fn fraction(u: &mut Unstructured<'_>, positive: bool) -> Result<Rational> {
    let d = u.int_in_range(1..=12i64)?;
    let k = u.int_in_range(i64::from(positive)..=d)?;
    Ok(q(k, d))
}

/// The value of `head[args]`.
fn eval(ctx: &Context, head: &str, args: &[Expr]) -> Expr {
    let e = Expr::call(Symbol::new(head), args.to_vec());
    ctx.eval(&e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The rational value of `head[args]`.
fn rational(ctx: &Context, head: &str, args: &[Expr]) -> Rational {
    let found = eval(ctx, head, args);
    found
        .to_rational()
        .unwrap_or_else(|| panic!("{head}{args:?} is {found:?}, not a rational"))
}

/// Whether `claim` evaluates to `True`.
fn holds(ctx: &Context, claim: &Expr) -> bool {
    ctx.eval(claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether `a` is within `10⁻¹²` of `b`, numerically.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [a.clone() - b.clone(), Expr::from(20)],
        )],
    );
    holds(
        ctx,
        &Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-12)]),
    )
}

fn power(b: &Rational, k: usize) -> Rational {
    (0..k).fold(Rational::ONE, |acc, _| &acc * b)
}

fn binomial(n: usize, k: usize) -> Rational {
    (0..k).fold(Rational::ONE, |acc, i| {
        &(&acc * &q(i64::try_from(n - i).expect("small"), 1))
            / &q(i64::try_from(i + 1).expect("small"), 1)
    })
}

/// The masses of a discrete family on `low ..= low + masses.len() − 1`,
/// checked against the CAS: masses, cumulative sums, an interval, the
/// inverse at `q`.
fn discrete(
    ctx: &Context,
    names: [&str; 3],
    params: &[Expr],
    low: usize,
    masses: &[Rational],
    u: &mut Unstructured<'_>,
) -> Result<()> {
    let [pdf, cdf, icdf] = names;
    let args = |extra: &[Expr]| [params, extra].concat();
    let mut total = Rational::ZERO;
    let mut cumulative = Vec::with_capacity(masses.len());
    for (i, m) in masses.iter().enumerate() {
        let k = int(low + i);
        let found = rational(ctx, pdf, &args(&[k.clone()]));
        assert_eq!(found, *m, "{pdf}{params:?} at {}", low + i);
        total = &total + m;
        cumulative.push(total.clone());
        let found = rational(ctx, cdf, &args(&[k]));
        assert_eq!(found, total, "{cdf}{params:?} at {}", low + i);
    }
    // An interval [a, b] of the support.
    let a = u.int_in_range(0..=masses.len() - 1)?;
    let b = u.int_in_range(a..=masses.len() - 1)?;
    let before = if a == 0 {
        Rational::ZERO
    } else {
        cumulative[a - 1].clone()
    };
    let found = rational(ctx, cdf, &args(&[int(low + a), int(low + b)]));
    assert_eq!(
        found,
        &cumulative[b] - &before,
        "{cdf}{params:?} on [{a}, {b}]"
    );
    // The inverse at a q the listed masses reach.
    let reached = cumulative.last().expect("masses").clone();
    let p = fraction(u, false)?;
    let target = &p * &reached;
    let k = cumulative
        .iter()
        .position(|c| *c >= target)
        .expect("the last reaches it");
    let found = rational(ctx, icdf, &args(&[Expr::from(target.clone())]));
    let expected = if target.is_zero() {
        q(i64::try_from(low).expect("small"), 1)
    } else {
        q(i64::try_from(low + k).expect("small"), 1)
    };
    assert_eq!(found, expected, "{icdf}{params:?} at {target:?}");
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0..=8u8)? {
        0 => {
            // Binomial: the masses sum to 1.
            let n = u.int_in_range(0..=10usize)?;
            let p = fraction(u, false)?;
            let one_minus = &Rational::ONE - &p;
            let masses: Vec<Rational> = (0..=n)
                .map(|k| &(&binomial(n, k) * &power(&p, k)) * &power(&one_minus, n - k))
                .collect();
            let sum = masses.iter().fold(Rational::ZERO, |s, m| &s + m);
            assert_eq!(sum, Rational::ONE, "binomial({n}, {p:?}) sums to 1");
            let params = [int(n), Expr::from(p)];
            discrete(
                &ctx,
                ["Binomial", "BinomialCDF", "BinomialICDF"],
                &params,
                0,
                &masses,
                u,
            )?;
        }
        1 => {
            // Geometric, its first masses.
            let p = fraction(u, true)?;
            let one_minus = &Rational::ONE - &p;
            let count = u.int_in_range(1..=8usize)?;
            let masses: Vec<Rational> = (0..count).map(|i| &p * &power(&one_minus, i)).collect();
            let params = [Expr::from(p)];
            discrete(
                &ctx,
                ["GeometricPDF", "GeometricCDF", "GeometricICDF"],
                &params,
                1,
                &masses,
                u,
            )?;
        }
        2 => {
            // Negative binomial: failures before the n-th success.
            let n = u.int_in_range(1..=4usize)?;
            let p = fraction(u, true)?;
            let one_minus = &Rational::ONE - &p;
            let count = u.int_in_range(1..=8usize)?;
            let masses: Vec<Rational> = (0..count)
                .map(|k| &(&binomial(n + k - 1, k) * &power(&p, n)) * &power(&one_minus, k))
                .collect();
            let params = [int(n), Expr::from(p)];
            discrete(
                &ctx,
                [
                    "NegativeBinomialPDF",
                    "NegativeBinomialCDF",
                    "NegativeBinomialICDF",
                ],
                &params,
                0,
                &masses,
                u,
            )?;
        }
        3 => {
            // Poisson, by the CAS: the cumulative function is the sum of
            // the masses, and the inverse brackets q.
            let lambda = Expr::from(q(u.int_in_range(1..=20)?, u.int_in_range(1..=4)?));
            let k = u.int_in_range(0..=8usize)?;
            let masses: Vec<Expr> = (0..=k)
                .map(|j| eval(&ctx, "PoissonPDF", &[lambda.clone(), int(j)]))
                .collect();
            let f = eval(&ctx, "PoissonCDF", &[lambda.clone(), int(k)]);
            let difference = Expr::call(Symbol::EXPAND, [f - Expr::call(Symbol::PLUS, masses)]);
            assert_eq!(
                ctx.eval(&difference).ok(),
                Some(Expr::from(0)),
                "PoissonCDF[{lambda:?}, {k}]"
            );
            let target = Expr::from(fraction(u, true)?);
            let found = eval(&ctx, "PoissonICDF", &[lambda.clone(), target.clone()]);
            if let Some(k) = found
                .to_integer()
                .and_then(|k| usize::try_from(k.to_i64()?).ok())
            {
                let at = eval(&ctx, "PoissonCDF", &[lambda.clone(), int(k)]);
                assert!(
                    holds(
                        &ctx,
                        &Expr::call(Symbol::GREATER_EQUAL, [at, target.clone()])
                    ),
                    "F({k}) ≥ q"
                );
                if k > 0 {
                    let below = eval(&ctx, "PoissonCDF", &[lambda.clone(), int(k - 1)]);
                    assert!(
                        holds(&ctx, &Expr::call(Symbol::LESS, [below, target.clone()])),
                        "F({k} − 1) < q"
                    );
                }
            }
        }
        4 => {
            // Uniform and exponential: F(F⁻¹(q)) = q, exactly.
            let p = Expr::from(fraction(u, false)?);
            let a = Expr::from(q(u.int_in_range(-20..=20)?, u.int_in_range(1..=4)?));
            let width = Expr::from(q(u.int_in_range(1..=20)?, u.int_in_range(1..=4)?));
            let b = ctx.eval(&(a.clone() + width)).expect("a sum");
            let x = eval(&ctx, "UniformICDF", &[a.clone(), b.clone(), p.clone()]);
            let back = eval(&ctx, "UniformCDF", &[a, b, x]);
            assert_eq!(back, p, "UniformCDF of UniformICDF");
            let lambda = Expr::from(q(u.int_in_range(1..=20)?, u.int_in_range(1..=4)?));
            let p = Expr::from(fraction(u, false)?);
            if p != Expr::from(1) {
                let x = eval(&ctx, "ExponentialICDF", &[lambda.clone(), p.clone()]);
                let back = eval(&ctx, "ExponentialCDF", &[lambda, x]);
                assert_eq!(back, p, "ExponentialCDF of ExponentialICDF");
            }
        }
        5 => {
            // Normal: F(F⁻¹(q)) ≈ q.
            let p = Expr::from(fraction(u, true)?);
            if p != Expr::from(1) {
                let mu = Expr::from(q(u.int_in_range(-10..=10)?, u.int_in_range(1..=4)?));
                let sigma = Expr::from(q(u.int_in_range(1..=10)?, u.int_in_range(1..=4)?));
                let x = eval(&ctx, "NormalICDF", &[mu.clone(), sigma.clone(), p.clone()]);
                let back = eval(&ctx, "NormalCDF", &[mu, sigma, x]);
                assert!(
                    near(&ctx, &back, &p),
                    "NormalCDF of NormalICDF at {p:?}: {back:?}"
                );
            }
        }
        6 => {
            // erf(erfinv(y)) ≈ y.
            let y = Expr::from(q(u.int_in_range(-11..=11)?, 12));
            let x = Expr::call(Symbol::new("InverseErf"), [y.clone()]);
            let back = Expr::call(
                Symbol::new("Erf"),
                [Expr::call(Symbol::N, [x, Expr::from(20)])],
            );
            assert!(near(&ctx, &back, &y), "Erf of InverseErf at {y:?}");
        }
        _ => {
            // Student, χ² and Fisher: F(F⁻¹(q)) ≈ q, F increasing in [0, 1].
            let n = u.int_in_range(1..=6usize)?;
            let m = u.int_in_range(1..=6usize)?;
            let (names, params): ([&str; 2], Vec<Expr>) = match u.int_in_range(0..=2u8)? {
                0 => (["StudentTCDF", "StudentTICDF"], vec![int(n)]),
                1 => (["ChiSquareCDF", "ChiSquareICDF"], vec![int(n)]),
                _ => (["FisherCDF", "FisherICDF"], vec![int(m), int(n)]),
            };
            let [cdf, icdf] = names;
            let args = |extra: Expr| [params.as_slice(), &[extra]].concat();
            let p = Expr::from(fraction(u, true)?);
            if p != Expr::from(1) {
                let x = eval(&ctx, icdf, &args(p.clone()));
                let back = eval(&ctx, cdf, &args(x));
                assert!(
                    near(&ctx, &back, &p),
                    "{cdf}{params:?} of {icdf} at {p:?}: {back:?}"
                );
            }
            let x1 = Expr::from(q(u.int_in_range(0..=20)?, 4));
            let x2 = Expr::from(q(u.int_in_range(0..=20)?, 4));
            let (f1, f2) = (
                eval(&ctx, cdf, &args(x1.clone())),
                eval(&ctx, cdf, &args(x2.clone())),
            );
            let ordered =
                |a: &Expr, b: &Expr| Expr::call(Symbol::LESS_EQUAL, [a.clone(), b.clone()]);
            assert!(
                holds(&ctx, &ordered(&Expr::from(0), &f1))
                    && holds(&ctx, &ordered(&f1, &Expr::from(1))),
                "{cdf}{params:?} at {x1:?}: {f1:?}"
            );
            if holds(&ctx, &Expr::call(Symbol::LESS, [x1.clone(), x2.clone()])) {
                assert!(
                    holds(&ctx, &ordered(&f1, &f2)),
                    "{cdf}{params:?} increasing: {f1:?} {f2:?}"
                );
            }
        }
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
