//! Tests of the probability distributions (CAS-DOC-23): masses,
//! cumulative probabilities and inverses, exact; what stays.

use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, answer)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input evaluates to itself.
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}

/// The binomial distribution: `binomial(n, p, k)`, giac's exact answers
/// (CAS-301).
#[test]
fn binomial() {
    check(&[
        ("Binomial[10, 3/10, 2]", "466948881/2000000000"),
        ("Binomial[10, 3/10, 11]", "0"),
        ("Binomial[10, 3/10, -1]", "0"),
        ("Binomial[10, 3/10, 5/2]", "0"),
        // The coefficient, as before.
        ("Binomial[10, 2]", "45"),
        // Symbols: expressions, as giac's.
        (
            "Binomial[10, p, 2]",
            "Times[45, Power[p, 2], Power[Plus[1, Times[-1, p]], 8]]",
        ),
        (
            "Binomial[10, 1/3, k]",
            "Times[Binomial[10, k], Power[1/3, k], Power[2/3, Plus[10, Times[-1, k]]]]",
        ),
        // p = 1: all of the mass at n, 0⁰ taken as 1.
        ("Binomial[3, 1, 3]", "1"),
        ("Binomial[3, 1, 2]", "0"),
        ("BinomialCDF[10, 3/10, 2]", "478478483/1250000000"),
        ("BinomialCDF[10, 3/10, 2, 4]", "1400846643/2000000000"),
        ("BinomialCDF[10, 3/10, 4, 2]", "0"),
        ("BinomialCDF[10, 1/2, 10]", "1"),
        ("BinomialCDF[10, 1/2, -1]", "0"),
        ("BinomialCDF[10, 1/2, 5/2]", "BinomialCDF[10, 1/2, 2]"),
        ("BinomialCDF[10, 1/2, Infinity]", "1"),
        ("BinomialICDF[10, 3/10, 1/2]", "3"),
        ("BinomialICDF[10, 1/2, 0]", "0"),
        ("BinomialICDF[10, 1/2, 1]", "10"),
        ("BinomialICDF[10, 1/2, 999/1000]", "9"),
        ("BinomialICDF[10, 0, 1]", "0"),
    ]);
}

/// Poisson's distribution: exact, where giac's answers are decimals
/// (CAS-301).
#[test]
fn poisson() {
    check(&[
        ("PoissonPDF[2, 3]", "Times[4/3, Exp[-2]]"),
        ("PoissonPDF[2, -1]", "0"),
        (
            "PoissonPDF[2, k]",
            "Times[Exp[-2], Power[2, k], Power[Factorial[k], -1]]",
        ),
        ("PoissonCDF[2, 3]", "Times[19/3, Exp[-2]]"),
        ("PoissonCDF[1, 1, 2]", "Times[3/2, Exp[-1]]"),
        ("PoissonCDF[1, -1]", "0"),
        ("PoissonICDF[2, 1/2]", "2"),
        ("PoissonICDF[1/2, 9/10]", "1"),
        ("PoissonICDF[2, 0]", "0"),
        ("PoissonICDF[2, 1]", "Infinity"),
    ]);
}

/// The geometric distribution, from 1, and the negative binomial, the
/// failures before the n-th success (CAS-301).
#[test]
fn geometric_and_negative_binomial() {
    check(&[
        ("GeometricPDF[1/3, 2]", "2/9"),
        ("GeometricPDF[1/3, 0]", "0"),
        (
            "GeometricPDF[1/3, k]",
            "Times[1/3, Power[2/3, Plus[-1, k]]]",
        ),
        ("GeometricCDF[1/3, 2]", "5/9"),
        ("GeometricCDF[1/3, 2, 4]", "38/81"),
        ("GeometricCDF[1/3, 0]", "0"),
        ("GeometricICDF[1/3, 1/2]", "2"),
        // (1/2)^2 = 1 − 3/4 exactly.
        ("GeometricICDF[1/2, 3/4]", "2"),
        ("GeometricICDF[1/3, 0]", "1"),
        ("GeometricICDF[1/3, 1]", "Infinity"),
        ("GeometricICDF[1, 1]", "1"),
        ("NegativeBinomialPDF[3, 1/2, 2]", "3/16"),
        ("NegativeBinomialCDF[3, 1/2, 2]", "1/2"),
        ("NegativeBinomialCDF[3, 1/2, 1, 2]", "3/8"),
        ("NegativeBinomialICDF[3, 1/2, 1/2]", "2"),
        ("NegativeBinomialICDF[3, 1, 1]", "0"),
    ]);
}

/// Parameters that are numbers out of range, symbols where a number is
/// needed, sums past the bound (CAS-301).
#[test]
fn discrete_stays() {
    stays(&[
        "Binomial[10, 3/2, 2]",
        "Binomial[-1, 1/2, 2]",
        "PoissonPDF[-1, 2]",
        "PoissonPDF[0, 2]",
        "GeometricPDF[0, 2]",
        "NegativeBinomialPDF[0, 1/2, 2]",
        "BinomialICDF[10, p, 1/2]",
        "BinomialICDF[10, 1/2, 3/2]",
        "BinomialICDF[10, 1/2, q]",
        "PoissonCDF[2, x]",
        "BinomialCDF[10000, 1/2, 5000]",
        "PoissonICDF[10000, 1/2]",
    ]);
}

/// Whether `claim` evaluates to `True`.
fn holds(ctx: &Context, claim: &str) -> bool {
    ctx.eval(&read(claim)).unwrap() == read("True")
}

/// `InverseErf`: its exact values, oddness, and balls that `Erf` takes
/// back to their argument (CAS-302).
#[test]
fn inverse_erf() {
    check(&[
        ("InverseErf[0]", "0"),
        ("InverseErf[1]", "Infinity"),
        ("InverseErf[-1]", "Times[-1, Infinity]"),
        ("InverseErf[-1/2]", "Times[-1, InverseErf[1/2]]"),
        ("InverseErf[List[0, 1]]", "List[0, Infinity]"),
    ]);
    stays(&["InverseErf[1/2]", "InverseErf[2]", "InverseErf[x]"]);
    let ctx = Context::new();
    for claim in [
        // To 30 digits, and erf takes them back.
        "Less[Abs[Plus[Erf[N[InverseErf[1/2], 30]], -1/2]], Power[10, -28]]",
        "Less[Abs[Plus[Erf[N[InverseErf[-99/100], 30]], 99/100]], Power[10, -28]]",
        "Less[Abs[Plus[Erf[N[InverseErf[1/1000000], 30]], -1/1000000]], Power[10, -34]]",
        // A decimal, at its precision.
        "Less[Abs[Plus[Erf[InverseErf[Decimal[9, -1]]], Times[-1, Decimal[9, -1]]]], \
         Power[10, -14]]",
        // erfinv(1/2) = 0.476936276204469873381…
        "Less[Abs[Plus[N[InverseErf[1/2], 30], Times[-1, Times[47693627620446987338, Power[10, -20]]]]], \
         Power[10, -19]]",
    ] {
        assert!(holds(&ctx, claim), "{claim}");
    }
}

/// The normal distribution, standard or not, its intervals by arity, and
/// its inverse through `erfinv` (CAS-303).
#[test]
fn normal() {
    check(&[
        ("NormalPDF[0]", "Power[Times[2, Pi], -1/2]"),
        (
            "NormalPDF[1, 2, x]",
            "Times[1/2, Power[Times[2, Pi], -1/2], \
             Exp[Times[-1/2, Power[Times[1/2, Plus[x, -1]], 2]]]]",
        ),
        ("NormalCDF[0]", "1/2"),
        ("NormalCDF[1]", "Times[1/2, Plus[1, Erf[Power[2, -1/2]]]]"),
        ("NormalCDF[Infinity]", "1"),
        ("NormalCDF[Times[-1, Infinity]]", "0"),
        ("NormalCDF[-1, 1]", "Erf[Power[2, -1/2]]"),
        ("NormalCDF[1, 2, 1]", "1/2"),
        ("NormalCDF[1, 2, -1, 3]", "Erf[Power[2, -1/2]]"),
        ("NormalICDF[1/2]", "0"),
        ("NormalICDF[1, 2, 1/2]", "1"),
        ("NormalICDF[0]", "Times[-1, Infinity]"),
        ("NormalICDF[1]", "Infinity"),
        (
            "NormalICDF[39/40]",
            "Times[Power[2, 1/2], InverseErf[19/20]]",
        ),
    ]);
    let ctx = Context::new();
    // The 97.5% quantile, 1.959963984540054…
    let claim = "Less[Abs[Plus[N[NormalICDF[39/40], 30], Times[-1959963984540054, Power[10, -15]]]], \
                 Power[10, -14]]";
    assert!(holds(&ctx, claim), "{claim}");
    stays(&[
        "NormalPDF[0, -1, 1]",
        "NormalPDF[0, 0, 1]",
        "NormalICDF[3/2]",
        "NormalCDF[1, 2, 3, 4, 5]",
    ]);
}

/// The exponential and uniform distributions, in closed form (CAS-303).
#[test]
fn exponential_and_uniform() {
    check(&[
        ("ExponentialPDF[2, 1]", "Times[2, Exp[-2]]"),
        ("ExponentialPDF[2, -1]", "0"),
        ("ExponentialPDF[l, x]", "Times[l, Exp[Times[-1, l, x]]]"),
        ("ExponentialCDF[2, 1]", "Plus[1, Times[-1, Exp[-2]]]"),
        (
            "ExponentialCDF[2, 1, 3]",
            "Plus[Exp[-2], Times[-1, Exp[-6]]]",
        ),
        ("ExponentialCDF[2, -1]", "0"),
        ("ExponentialICDF[2, 1/2]", "Times[-1/2, Log[1/2]]"),
        ("ExponentialICDF[2, 0]", "0"),
        ("ExponentialICDF[2, 1]", "Infinity"),
        ("UniformPDF[0, 2, 1]", "1/2"),
        ("UniformPDF[0, 2, 3]", "0"),
        ("UniformPDF[a, b, x]", "Power[Plus[b, Times[-1, a]], -1]"),
        ("UniformCDF[0, 2, 1]", "1/2"),
        ("UniformCDF[0, 2, 3]", "1"),
        ("UniformCDF[0, 2, -1]", "0"),
        ("UniformCDF[0, 2, 1/2, 3/2]", "1/2"),
        ("UniformCDF[0, 2, 3/2, 1/2]", "0"),
        ("UniformICDF[0, 2, 1/4]", "1/2"),
    ]);
    stays(&[
        "ExponentialPDF[-2, 1]",
        "UniformPDF[2, 0, 1]",
        "UniformICDF[0, 2, 2]",
    ]);
}

/// Student, χ² and Fisher: exact densities and cumulative probabilities,
/// certified inverses (CAS-304).
#[test]
fn degrees_of_freedom() {
    check(&[
        ("StudentTCDF[1, 1]", "3/4"),
        ("StudentTCDF[1, -1]", "1/4"),
        ("StudentTCDF[3, 0]", "1/2"),
        (
            "StudentTCDF[3, 1]",
            "Plus[2/3, Times[1/4, Power[3, 1/2], Power[Pi, -1]]]",
        ),
        ("StudentTCDF[3, Infinity]", "1"),
        ("ChiSquareCDF[4, 2]", "Plus[1, Times[-2, Exp[-1]]]"),
        ("ChiSquareCDF[2, 2]", "Plus[1, Times[-1, Exp[-1]]]"),
        ("ChiSquareCDF[3, -1]", "0"),
        ("FisherCDF[2, 4, 1]", "5/9"),
        ("FisherCDF[2, 4, 0]", "0"),
        ("FisherPDF[2, 4, 1]", "8/27"),
        ("ChiSquarePDF[2, 1]", "Times[1/2, Exp[-1/2]]"),
        ("ChiSquarePDF[4, 2]", "Times[1/2, Exp[-1]]"),
        ("ChiSquarePDF[3, -1]", "0"),
        (
            "StudentTPDF[1, x]",
            "Times[Power[Pi, -1], Power[Plus[1, Power[x, 2]], -1]]",
        ),
        ("StudentTICDF[3, 1/2]", "0"),
        ("StudentTICDF[3, 0]", "Times[-1, Infinity]"),
        ("ChiSquareICDF[2, 0]", "0"),
        ("FisherICDF[2, 3, 1]", "Infinity"),
    ]);
    // Outside `N` a ball beside an exact irrational, `Log[2]` or `Pi` in a
    // closed form, stays a sum: the claims compare under `N`.
    let ctx = Context::new();
    for claim in [
        // t(0.95, 3) = 2.353363434801826…, and its mirror.
        "Less[Abs[Plus[StudentTICDF[3, 95/100], Times[-2353363434801826, Power[10, -15]]]], Power[10, -12]]",
        "Less[Abs[Plus[StudentTICDF[3, 5/100], Times[2353363434801826, Power[10, -15]]]], Power[10, -12]]",
        // χ²(2) is E(1/2): its median is 2 ln 2.
        "Less[Abs[N[Plus[ChiSquareICDF[2, 1/2], Times[-2, Log[2]]], 20]], Power[10, -13]]",
        // F(2, 3)'s median is (3/2)(2^(2/3) − 1).
        "Less[Abs[N[Plus[FisherICDF[2, 3, 1/2], Times[-3/2, Plus[Power[2, 2/3], -1]]], 20]], \
         Power[10, -13]]",
        // The cumulative function takes an inverse back.
        "Less[Abs[N[Plus[StudentTCDF[3, StudentTICDF[3, 9/10]], -9/10], 20]], Power[10, -12]]",
        "Less[Abs[N[Plus[ChiSquareCDF[5, ChiSquareICDF[5, 1/3]], -1/3], 20]], Power[10, -12]]",
        // An exact root at a midpoint of the bisection: t(2) at 8/9 is 7/4,
        // which the bisection of [1, 3] reaches third and could not sign;
        // it gave the ball of [3/2, 2] (the fuzz target, CAS-305).
        "Less[Abs[N[Plus[StudentTICDF[2, 8/9], -7/4], 20]], Power[10, -14]]",
        "Less[Abs[N[Plus[StudentTICDF[2, 1/9], 7/4], 20]], Power[10, -14]]",
    ] {
        assert!(holds(&ctx, claim), "{claim}");
    }
    stays(&[
        "StudentTPDF[0, 1]",
        "StudentTPDF[-1, 1]",
        "StudentTPDF[5/2, 1]",
        "ChiSquareCDF[2, x]",
        "FisherICDF[2, 3, 3/2]",
    ]);
}
