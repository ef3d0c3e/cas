//! Tests of distribution fitting (CAS-629): `fitdistr`'s closed forms,
//! exact on exact data, its iterative fits against the roots of the score
//! equations (mpmath at 40 digits, from .run/1.7/T5/proto.py's equations
//! with the mean `Σx/N`), each warning, and the budget.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// `[1.2, 2.3, 0.7, 3.1, 1.9, 2.2]`.
const DECIMALS: &str = "N[List[6/5, 23/10, 7/10, 31/10, 19/10, 11/5]]";
/// `[0.2, 0.3, 0.7, 0.1, 0.9, 0.4]`.
const UNIT: &str = "N[List[1/5, 3/10, 7/10, 1/10, 9/10, 2/5]]";

/// Exact data give exact objects.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let eight = "List[1, 2, 3, 4, 5, 2, 3, 3]";
    let five = "List[1, 2, 3, 4, 5]";
    for (input, expected) in [
        (
            format!("FitDistribution[{five}, Together]"),
            "NormalPDF[3, Sqrt[5/2]]",
        ),
        (
            format!("FitDistribution[{five}, NORMALD]"),
            "NormalPDF[3, Sqrt[5/2]]",
        ),
        (
            format!("FitDistribution[{eight}, NormalPDF]"),
            "NormalPDF[23/8, Sqrt[87/56]]",
        ),
        (
            format!("FitDistribution[{eight}, PoissonPDF]"),
            "PoissonPDF[23/8]",
        ),
        (format!("FitDistribution[{five}, POISSON]"), "PoissonPDF[3]"),
        (
            format!("FitDistribution[{eight}, ExponentialPDF]"),
            "ExponentialPDF[8/23]",
        ),
        (
            format!("FitDistribution[{five}, Exp]"),
            "ExponentialPDF[1/3]",
        ),
        (
            format!("FitDistribution[{eight}, GeometricPDF]"),
            "GeometricPDF[8/23]",
        ),
        (
            "FitDistribution[List[2, 2, 2], Together]".to_owned(),
            "NormalPDF[2, 0]",
        ),
    ] {
        assert_eq!(eval(&ctx, &input), eval(&ctx, expected), "{input}");
    }
}

/// `found`'s head is `head` and its parameters balls within a relative
/// `10⁻¹²` of `expected`, their radii below it.
fn close(input: &str, found: &Expr, head: Symbol, expected: &[f64]) {
    assert_eq!(found.head(), Some(&Expr::from(head)), "{input}: {found:?}");
    assert_eq!(found.args().len(), expected.len(), "{input}: {found:?}");
    for (x, v) in found.args().iter().zip(expected) {
        let View::Ball(b) = x.view() else {
            panic!("{input}: {found:?}");
        };
        let radius = b.rad().to_float().map_or(f64::INFINITY, |r| r.to_f64());
        assert!(
            (b.mid().to_f64() - v).abs() <= 1e-12 * v.abs() && radius <= 1e-12 * v.abs(),
            "{input}: {found:?}, expected {expected:?}"
        );
    }
}

/// Decimal data, and the iterative fits on exact data, give certified balls.
#[test]
fn iterative_fits() {
    let ctx = Context::new();
    let five = "List[1, 2, 3, 4, 5]";
    for (input, head, expected) in [
        (
            format!("FitDistribution[{DECIMALS}, Together]"),
            Symbol::NORMAL_PDF,
            [1.9, 0.850_881_895_447_306],
        ),
        (
            format!("FitDistribution[{DECIMALS}, GammaPDF]"),
            Symbol::GAMMA_PDF,
            [4.916_025_249_443_367, 2.587_381_710_233_351],
        ),
        (
            format!("FitDistribution[{DECIMALS}, WeibullPDF]"),
            Symbol::WEIBULL_PDF,
            [2.711_528_234_976_74, 2.140_833_177_232_44],
        ),
        (
            format!("FitDistribution[{DECIMALS}, CauchyPDF]"),
            Symbol::CAUCHY_PDF,
            [2.076_582_676_865_552, 0.430_024_064_960_904_9],
        ),
        (
            format!("FitDistribution[{UNIT}, BetaPDF]"),
            Symbol::BETA_PDF,
            [1.158_398_650_961_8, 1.416_640_983_180_957_3],
        ),
        (
            format!("FitDistribution[{five}, Gamma]"),
            Symbol::GAMMA_PDF,
            [3.701_643_810_008_816_7, 1.233_881_270_002_939],
        ),
        (
            format!("FitDistribution[{five}, WeibullPDF]"),
            Symbol::WEIBULL_PDF,
            [2.293_806_670_712_317, 3.394_290_718_109_487],
        ),
        (
            format!("FitDistribution[{five}, CauchyPDF]"),
            Symbol::CAUCHY_PDF,
            [3.0, 1.062_017_918_685_918_7],
        ),
        (
            "FitDistribution[N[List[1/5, 1/2, 3/5, 9/10]], Beta]".to_owned(),
            Symbol::BETA_PDF,
            [1.794_813_870_078_337_1, 1.438_184_441_835_130_3],
        ),
        (
            "FitDistribution[N[List[1/2, 3/2, 2, 7/2, 4, 7, 19/2]], CauchyPDF]".to_owned(),
            Symbol::CAUCHY_PDF,
            [2.883_711_007_482_26, 1.687_891_017_701_232_2],
        ),
    ] {
        close(&input, &eval(&ctx, &input), head, &expected);
    }
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap(), e, "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

/// giac's errors, in its order.
#[test]
fn stays() {
    let ctx = Context::new();
    for (input, message) in [
        ("FitDistribution[List[1, 2, 3]]", "Bad Argument Type"),
        (
            "FitDistribution[List[1, 2, 3], Together, 1]",
            "Expected a sequence of two arguments",
        ),
        ("FitDistribution[3, Together]", "Expected a list"),
        (
            "FitDistribution[List[1], Together]",
            "At least two samples are required",
        ),
        (
            "FitDistribution[List[1, a], foo]",
            "Expected a real constant",
        ),
        (
            "FitDistribution[List[1, I], Together]",
            "Expected a real constant",
        ),
        (
            "FitDistribution[List[1, 2, 3], foo]",
            "Unrecognized probability distribution",
        ),
        (
            "FitDistribution[List[1, 2, 3], Binomial]",
            "Unrecognized probability distribution",
        ),
        (
            "FitDistribution[List[1, 2, 3], UniformPDF]",
            "Unrecognized probability distribution",
        ),
        (
            "FitDistribution[List[1, 2, 3], NegativeBinomialPDF]",
            "Unrecognized probability distribution",
        ),
        (
            "FitDistribution[List[1, 2, 3], RandomNormal]",
            "Unrecognized probability distribution",
        ),
        (
            "FitDistribution[List[1, 5/2], PoissonPDF]",
            "Expected a nonnegative integer",
        ),
        (
            "FitDistribution[List[1, -2], PoissonPDF]",
            "Expected a nonnegative integer",
        ),
        (
            "FitDistribution[List[1, -2], ExponentialPDF]",
            "Expected a positive number",
        ),
        (
            "FitDistribution[List[0, 2], GeometricPDF]",
            "Expected a positive integer",
        ),
        (
            "FitDistribution[List[1, -2], GammaPDF]",
            "Expected a positive number",
        ),
        (
            "FitDistribution[List[0, 2], GammaPDF]",
            "Expected a positive number",
        ),
        (
            "FitDistribution[List[1/2, 3/2], BetaPDF]",
            "Expected a number in [0,1]",
        ),
        (
            "FitDistribution[List[0, 1/2], BetaPDF]",
            "Unable to compute parameters for Beta distribution",
        ),
        (
            "FitDistribution[List[1/2, 1], BetaPDF]",
            "Unable to compute parameters for Beta distribution",
        ),
        (
            "FitDistribution[List[1/2, 1/2], BetaPDF]",
            "Unable to compute parameters for Beta distribution",
        ),
        (
            "FitDistribution[List[0, 0], WeibullPDF]",
            "Expected a nonzero variance",
        ),
        (
            "FitDistribution[List[2, 2, 2], WeibullPDF]",
            "Expected a nonzero variance",
        ),
        (
            "FitDistribution[List[1, -2], WeibullPDF]",
            "Expected a positive number",
        ),
        (
            "FitDistribution[List[0, 2], WeibullPDF]",
            "Expected a positive number",
        ),
    ] {
        assert_eq!(warning(&ctx, input), message, "{input}");
    }
}

/// No root found, or the budget spent: the call stays without a warning.
#[test]
fn stays_silently() {
    let quiet = |ctx: &Context, input: &str| {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap(), e, "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    };
    let ctx = Context::new();
    // Equal data: gamma's start m²/v divides by 0; Cauchy's scale starts at 0.
    quiet(&ctx, "FitDistribution[List[2, 2, 2], GammaPDF]");
    quiet(&ctx, "FitDistribution[List[2, 2, 2], CauchyPDF]");
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    // Six data: three scores at most, short of a Newton solve and its certificate.
    quiet(&small, &format!("FitDistribution[{DECIMALS}, CauchyPDF]"));
    quiet(&small, &format!("FitDistribution[{DECIMALS}, WeibullPDF]"));
}
