//! Tests of random variables (CAS-624): `randvar`'s objects (the rows of
//! the prototype's fits and discrete forms, giac's answers made exact),
//! each warning, what stays without one, the supports of seeded draws of
//! `sample`, `randvector` and `randmultinomial`, a seed repeating a sample,
//! and the size refusals.

#![allow(clippy::cast_precision_loss, reason = "small test values")]

use cas_num::{Integer, Rational};

use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression,
/// with no warning.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        let found = found.unwrap();
        let expected = eval(&ctx, expected);
        if found != expected || !report.warnings.is_empty() {
            wrong.push(format!(
                "{input}: {found:?}, expected {expected:?}, {:?}",
                report.warnings
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head(), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

/// A number in `f64`: a ball's midpoint or a small integer.
fn float(e: &Expr) -> f64 {
    match e.view() {
        View::Ball(b) => b.mid().to_f64(),
        _ => e
            .as_i64()
            .unwrap_or_else(|| panic!("{e:?}: no small number")) as f64,
    }
}

/// The fits, as the prototype computes them (proto.out), and giac's
/// other objects.
#[test]
fn fits() {
    check(&[
        (
            "RandomVariable[GammaPDF, Equation[Mean, 12], Equation[PopulationVariance, 4]]",
            "GammaPDF[36, 3]",
        ),
        (
            "RandomVariable[Binomial, Equation[Mean, 18], Equation[PopulationStandardDeviation, 4]]",
            "RandomVariable[Binomial, 162, 1/9]",
        ),
        (
            "RandomVariable[UniformPDF, Equation[Mean, 10], Equation[PopulationStandardDeviation, 2]]",
            "UniformPDF[Plus[10, Times[-2, Sqrt[3]]], Plus[10, Times[2, Sqrt[3]]]]",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[PopulationStandardDeviation, 3/2]]",
            "NormalPDF[5, 3/2]",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[PopulationVariance, 2], \
             Equal[PopulationStandardDeviation, 1]]",
            "NormalPDF[5, 1]",
        ),
        (
            "RandomVariable[PoissonPDF, Equation[Mean, 3]]",
            "PoissonPDF[3]",
        ),
        (
            "RandomVariable[ExponentialPDF, Equation[Mean, 2]]",
            "ExponentialPDF[1/2]",
        ),
        (
            "RandomVariable[Exp, Equation[Mean, 2]]",
            "ExponentialPDF[1/2]",
        ),
        (
            "RandomVariable[GeometricPDF, Equation[Mean, 4]]",
            "GeometricPDF[1/4]",
        ),
        (
            "RandomVariable[GeometricPDF, Equation[PopulationVariance, 12]]",
            "GeometricPDF[1/4]",
        ),
        (
            "RandomVariable[BetaPDF, Equation[Mean, 2/5], Equation[PopulationVariance, 1/25]]",
            "BetaPDF[2, 3]",
        ),
        (
            "RandomVariable[Beta, Equation[Mean, 2/5], Equation[PopulationVariance, 1/25]]",
            "BetaPDF[2, 3]",
        ),
        (
            "RandomVariable[Gamma, Equation[Mean, 12], Equation[PopulationVariance, 4]]",
            "GammaPDF[36, 3]",
        ),
        (
            "RandomVariable[Binomial, Equation[Mean, 10], Equation[PopulationStandardDeviation, 2]]",
            "RandomVariable[Binomial, 17, 10/17]",
        ),
        (
            "RandomVariable[NegativeBinomialPDF, Equation[Mean, 6], \
             Equation[PopulationStandardDeviation, 3]]",
            "NegativeBinomialPDF[12, 2/3]",
        ),
        (
            "RandomVariable[UniformPDF, Equation[range, Range[1, 3]]]",
            "UniformPDF[1, 3]",
        ),
        (
            "RandomVariable[UniformPDF, Equation[range, List[1, 3]]]",
            "UniformPDF[1, 3]",
        ),
        (
            "RandomVariable[UniformPDF, Range[E, Pi]]",
            "UniformPDF[E, Pi]",
        ),
    ]);
}

/// A law named alone or with its parameters, and the objects themselves.
#[test]
fn positional() {
    check(&[
        ("RandomVariable[Together]", "NormalPDF[0, 1]"),
        ("RandomVariable[NormalPDF]", "NormalPDF[0, 1]"),
        ("RandomVariable[RandomNormal]", "NormalPDF[0, 1]"),
        ("RandomVariable[UniformPDF]", "UniformPDF[0, 1]"),
        ("RandomVariable[Exp]", "ExponentialPDF[1]"),
        ("RandomVariable[ExponentialPDF]", "ExponentialPDF[1]"),
        ("RandomVariable[Together, 5, 2]", "NormalPDF[5, 2]"),
        ("RandomVariable[Exp, 2]", "ExponentialPDF[2]"),
        ("RandomVariable[PoissonPDF, 3]", "PoissonPDF[3]"),
        ("RandomVariable[GeometricPDF, 1/4]", "GeometricPDF[1/4]"),
        ("RandomVariable[StudentTPDF, 3]", "StudentTPDF[3]"),
        ("RandomVariable[ChiSquarePDF, 3]", "ChiSquarePDF[3]"),
        ("RandomVariable[CauchyPDF, 0, 1]", "CauchyPDF[0, 1]"),
        ("RandomVariable[WeibullPDF, 2, 3]", "WeibullPDF[2, 3]"),
        ("RandomVariable[FisherPDF, 2, 3]", "FisherPDF[2, 3]"),
        (
            "RandomVariable[Together, mu, sigma]",
            "NormalPDF[mu, sigma]",
        ),
        ("RandomVariable[NormalPDF[5, 2]]", "NormalPDF[5, 2]"),
        ("RandomVariable[CauchyPDF[0, 1]]", "CauchyPDF[0, 1]"),
        (
            "RandomVariable[RandomVariable[List[1, 2]]]",
            "RandomVariable[List[1, 2]]",
        ),
        (
            "RandomVariable[Binomial, 10, 1/2]",
            "RandomVariable[Binomial, 10, 1/2]",
        ),
    ]);
}

/// The discrete forms: canonical objects, fixed points of `RandomVariable`.
#[test]
fn discrete_forms() {
    let support: Vec<i64> = (-10..=10).collect();
    let weights: Vec<String> = support
        .iter()
        .map(|&k| {
            let w = Rational::new(Integer::from(100 - k * k), Integer::from(100));
            if w.is_integer() {
                w.numer().to_string()
            } else {
                format!("{}/{}", w.numer(), w.denom())
            }
        })
        .collect();
    let points: Vec<String> = support.iter().map(i64::to_string).collect();
    let density = format!(
        "RandomVariable[List[{}], List[{}]]",
        weights.join(", "),
        points.join(", ")
    );
    check(&[
        (
            "RandomVariable[List[3, 1, 2, 5], List[alpha, beta, gamma, delta]]",
            "RandomVariable[List[3, 1, 2, 5], List[alpha, beta, gamma, delta]]",
        ),
        (
            "RandomVariable[List[List[\"apple\", 1/3], List[\"orange\", 1/4], List[\"pear\", 1/5], \
             List[\"plum\", 13/60]]]",
            "RandomVariable[List[1/3, 1/4, 1/5, 13/60], List[\"apple\", \"orange\", \"pear\", \"plum\"]]",
        ),
        (
            "RandomVariable[List[1, 2, 3]]",
            "RandomVariable[List[1, 2, 3]]",
        ),
        (
            "RandomVariable[List[0x1p-1, 0]]",
            "RandomVariable[List[0x1p-1, 0]]",
        ),
        (
            "RandomVariable[multinomial, List[1/2, 1/3, 1/6], List[a, b, c]]",
            "RandomVariable[List[1/2, 1/3, 1/6], List[a, b, c]]",
        ),
        (
            "RandomVariable[multinomial, List[1/2, 1/2]]",
            "RandomVariable[List[1/2, 1/2]]",
        ),
        (
            "RandomVariable[Function[List[k], Plus[1, Times[-1, Power[Times[k, 1/10], 2]]]], \
             Equation[range, Range[-10, 10]]]",
            &density,
        ),
        (
            "RandomVariable[Function[List[k], k], Equation[range, Range[1, 3]], List[a, b, c]]",
            "RandomVariable[List[1, 2, 3], List[a, b, c]]",
        ),
    ]);
}

/// Weibull's fit by Newton's method on giac's equation: giac's
/// `weibulld(3.08574940721, 13.9803128143)` for `mean=12.5, variance=1`
/// (the prototype's digits), from a decimal mean and an exact one.
#[test]
fn weibull() {
    let ctx = Context::new();
    for input in [
        "RandomVariable[WeibullPDF, Equation[Mean, 0x1.9p+3], Equation[PopulationVariance, 1]]",
        "RandomVariable[WeibullPDF, Equation[Mean, 25/2], Equation[PopulationVariance, 1]]",
    ] {
        let found = eval(&ctx, input);
        assert_eq!(
            found.head(),
            read("WeibullPDF[1, 1]").head(),
            "{input}: {found:?}"
        );
        let [k, lambda] = found.args() else {
            panic!("{input}: {found:?}");
        };
        assert!((float(k) - 3.085_749_407_212_56).abs() < 1e-9, "{k:?}");
        assert!(
            (float(lambda) - 13.980_312_814_270_7).abs() < 1e-9,
            "{lambda:?}"
        );
    }
}

/// giac's errors, each with giac's message.
#[test]
#[allow(clippy::too_many_lines, reason = "a table of inputs")]
fn warnings() {
    let ctx = Context::new();
    for (input, message) in [
        ("RandomVariable[FisherPDF, 2]", "Bad Argument Value"),
        ("RandomVariable[Binomial, 10]", "Bad Argument Value"),
        ("RandomVariable[List[1, -2, 3]]", "Bad Argument Value"),
        ("RandomVariable[List[0, 0]]", "Bad Argument Value"),
        ("RandomVariable[List[x, 1]]", "Bad Argument Value"),
        (
            "RandomVariable[Function[List[k], k], Equation[range, Range[0x1.8p+0, 3]]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[Function[List[k], k], Equation[range, Range[3, 1]]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[multinomial, Equation[Mean, 1]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[PoissonPDF, Equation[range, List[1, 2, 3]]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[PoissonPDF, Equation[Mean, 3], Equation[PopulationVariance, 2]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[GammaPDF, Equation[Mean, 3]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[UniformPDF, Equation[Mean, 3]]",
            "Bad Argument Value",
        ),
        (
            "RandomVariable[UniformPDF, Range[3, 1]]",
            "Invalid distribution parameter(s)",
        ),
        (
            "RandomVariable[BetaPDF, Equation[Mean, 1/2], Equation[PopulationVariance, 1]]",
            "Invalid distribution parameter(s)",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[PopulationVariance, -1]]",
            "Variance must be nonnegative",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[range, x]]",
            "Invalid range specification",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[foo, 2]]",
            "Unknown distribution parameter",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], \"x\"]",
            "Invalid distribution specification",
        ),
        (
            "RandomVariable[FisherPDF, Equation[Mean, 2], Equation[PopulationStandardDeviation, 1]]",
            "Specifying moments is not supported for this distribution",
        ),
        (
            "RandomVariable[CauchyPDF, Equation[Mean, 2], Equation[PopulationStandardDeviation, 1]]",
            "Specifying moments is not supported for this distribution",
        ),
        (
            "RandomVariable[List[1, 2], List[a, b, c]]",
            "Items and weights sizes do not match",
        ),
        (
            "RandomVariable[List[List[1, 2, 3], List[4, 5, 6]]]",
            "Expected a two-column matrix",
        ),
        (
            "RandomVariable[Together, Equation[Mean, 5], Equation[PopulationStandardDeviation, -2]]",
            "Standard deviation must be positive",
        ),
        (
            "RandomVariable[Function[List[k], k], Equation[range, Range[0, 10000000]]]",
            "too many weights",
        ),
        ("RandomMultinomial[List[1/2, 1/3]]", "Bad Argument Value"),
        (
            "RandomMultinomial[List[1/2, 1/2], List[a]]",
            "Bad Argument Value",
        ),
        ("RandomMultinomial[List[3/2, -1/2]]", "Bad Argument Value"),
        ("RandomMultinomial[List[1/2, 1/2], x]", "Bad Argument Value"),
        ("RandomSample[PoissonPDF[3], 10000000]", "too many draws"),
        (
            "RandomSample[Plus[PoissonPDF[3], Power[PoissonPDF[4], 2]], 600000]",
            "too many draws",
        ),
        (
            "RandomSample[PoissonPDF[3], -1]",
            "a count is a natural number",
        ),
    ] {
        assert_eq!(warning(&ctx, input), message, "{input}");
    }
}

/// Calls that stay without a warning.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "RandomVariable[]",
        "RandomVariable[PoissonPDF]",
        "RandomVariable[multinomial]",
        "RandomVariable[x]",
        "RandomVariable[x, 1, 2]",
        "RandomVariable[Together, Equation[Mean, m]]",
        "RandomVariable[Together, Equation[PopulationStandardDeviation, s]]",
        "RandomVariable[List[1, 2, 3]]",
        "RandomVariable[List[1/2, 1/2], List[a, b]]",
        "RandomVariable[Binomial, n, p]",
        "RandomMultinomial[x]",
        "RandomMultinomial[]",
        "RandomSample[x, 3]",
        "RandomSample[x]",
        "RandomSample[PoissonPDF[3], n]",
        "RandomSample[NormalPDF[mu, 1], 3]",
        "RandomVector[3, multinomial, List[1/2, 1/3]]",
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap().head(), e.head(), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// Seeded draws land in their supports.
#[test]
fn supports() {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[42]");
    let values = |input: &str| -> Vec<Expr> {
        let found = eval(&ctx, input);
        assert_eq!(found.head(), read("List[]").head(), "{input}: {found:?}");
        found.args().to_vec()
    };
    let named = values(
        "RandomSample[RandomVariable[List[3, 1, 2, 5], List[alpha, beta, gamma, delta]], 200]",
    );
    assert_eq!(named.len(), 200);
    for name in ["alpha", "beta", "gamma", "delta"] {
        assert!(named.contains(&read(name)), "{name}");
    }
    assert!(named.iter().all(|v| {
        ["alpha", "beta", "gamma", "delta"]
            .iter()
            .any(|n| *v == read(n))
    }));
    let indices = values("RandomSample[RandomVariable[List[1, 2, 3]], 100]");
    assert!(indices.iter().all(|i| matches!(i.as_i64(), Some(0..=2))));
    let zero = values("RandomSample[RandomVariable[List[0, 1, 0]], 50]");
    assert!(zero.iter().all(|i| i.as_i64() == Some(1)), "{zero:?}");
    let decimal = values("RandomSample[RandomVariable[List[0x1p-1, 0, 0x1p-1]], 50]");
    assert!(
        decimal.iter().all(|i| matches!(i.as_i64(), Some(0 | 2))),
        "{decimal:?}"
    );
    let poisson = values("RandomSample[PoissonPDF[3], 50]");
    assert!(poisson.iter().all(|i| i.as_i64().is_some_and(|i| i >= 0)));
    let binomial = values("RandomSample[RandomVariable[Binomial, 10, 1/2], 50]");
    assert!(binomial.iter().all(|i| matches!(i.as_i64(), Some(0..=10))));
    let uniform = values("RandomSample[UniformPDF[1, 3], 20]");
    assert!(uniform.iter().all(|x| (1.0..3.0).contains(&float(x))));
    let cauchy = values("RandomSample[CauchyPDF[0, 1], 5]");
    assert!(cauchy.iter().all(|x| float(x).is_finite()));
    let one = eval(&ctx, "RandomSample[NormalPDF[0, 1]]");
    assert!(float(&one).is_finite());
    let vector = values("RandomVector[5, RandomVariable[List[1, 2], List[a, b]]]");
    assert!(vector.iter().all(|v| *v == read("a") || *v == read("b")));
    let cauchy = values("RandomVector[3, CauchyPDF, 0, 1]");
    assert_eq!(cauchy.len(), 3);
    let elements = values("RandomVector[3, multinomial, List[1/2, 1/3, 1/6], List[a, b, c]]");
    assert!(
        elements
            .iter()
            .all(|v| ["a", "b", "c"].iter().any(|n| *v == read(n)))
    );
    // Without values, the counts of each index, n in all.
    let counts = values("RandomVector[20, multinomial, List[1/2, 1/3, 1/6]]");
    assert_eq!(counts.len(), 3);
    assert_eq!(counts.iter().map(|c| c.as_i64().unwrap()).sum::<i64>(), 20);
    for _ in 0..20 {
        let i = eval(&ctx, "RandomMultinomial[List[1/2, 1/3, 1/6]]");
        assert!(matches!(i.as_i64(), Some(0..=2)), "{i:?}");
        let v = eval(
            &ctx,
            "RandomMultinomial[List[1/2, 1/3, 1/6], List[a, b, c]]",
        );
        assert!(["a", "b", "c"].iter().any(|n| v == read(n)), "{v:?}");
        let d = eval(&ctx, "RandomMultinomial[List[0x1p-1, 0x1p-1]]");
        assert!(matches!(d.as_i64(), Some(0 | 1)), "{d:?}");
        let zero = eval(&ctx, "RandomMultinomial[List[0, 1]]");
        assert_eq!(zero.as_i64(), Some(1));
    }
}

/// Each occurrence of an object in an expression is drawn afresh: `X +
/// X²` with `X` in `{0, 2}` gives 2 or 4 only so.
#[test]
fn expressions() {
    let ctx = Context::new();
    eval(&ctx, "SeedRandom[3]");
    let found = eval(
        &ctx,
        "RandomSample[Plus[RandomVariable[List[1, 1], List[0, 2]], \
         Power[RandomVariable[List[1, 1], List[0, 2]], 2]], 64]",
    );
    let values: Vec<i64> = found.args().iter().map(|v| v.as_i64().unwrap()).collect();
    assert_eq!(values.len(), 64);
    assert!(
        values.iter().all(|v| [0, 2, 4, 6].contains(v)),
        "{values:?}"
    );
    assert!(values.iter().any(|v| [2, 4].contains(v)), "{values:?}");
}

/// A seed repeats a sample.
#[test]
fn seeds_repeat() {
    let ctx = Context::new();
    let twice = eval(
        &ctx,
        "List[SeedRandom[7], RandomSample[NormalPDF[0, 1], 5], \
         SeedRandom[7], RandomSample[NormalPDF[0, 1], 5]]",
    );
    assert_eq!(twice.args()[1], twice.args()[3]);
    let discrete = eval(
        &ctx,
        "List[SeedRandom[7], RandomSample[RandomVariable[List[1, 2, 3]], 10], \
         SeedRandom[7], RandomSample[RandomVariable[List[1, 2, 3]], 10]]",
    );
    assert_eq!(discrete.args()[1], discrete.args()[3]);
}
