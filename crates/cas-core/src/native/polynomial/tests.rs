//! Tests of the polynomial builtins: giac's answers where this CAS follows
//! them (CAS-DOC-8, CAS-140), the warnings, balls, the limit on terms, the
//! interrupt, the trace and its verification, and proptests that `Expand`
//! is cas-poly's product and `Together` keeps the value.

use std::time::{Duration, Instant};

use cas_num::algebra::{Integers, Ring};
use cas_num::{Integer, Rational};
use cas_poly::{Never, Sparse, SparseRing};
use proptest::prelude::*;

use crate::poly::{Fraction, Polynomials};
use crate::trace::{EvalOptions, Step, StepKind, TraceLevel};
use crate::{Attributes, Context, EvalError, Expr, Interrupt, Symbol, Warning};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// The value of `text` and the messages of its warnings.
fn eval(ctx: &Context, text: &str) -> (Expr, Vec<String>) {
    let (result, report) = ctx.eval_with(&read(text), &EvalOptions::default());
    let warnings = report.warnings.iter().map(Warning::to_string).collect();
    (result.unwrap_or_else(|e| panic!("{text}: {e}")), warnings)
}

/// `(input, result)` in full form, each without a warning.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let (found, warnings) = eval(&ctx, input);
        if format!("{found:?}") != *expected || !warnings.is_empty() {
            wrong.push(format!(
                "{input}: {found:?} {warnings:?}, expected {expected}"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `(input, warning)`: the call stays, with the warning.
fn check_refused(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, message) in cases {
        let (found, warnings) = eval(&ctx, input);
        let stays = found.head() == read(input).head();
        if !stays || !warnings.iter().any(|w| w.contains(message)) {
            wrong.push(format!(
                "{input}: {found:?} {warnings:?}, expected {message:?}"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac 1.9's `expand` and `normal` (CAS-DOC-8 §7, §8), in full form.
#[test]
fn expand_and_normal() {
    check(&[
        // expand((x+1)^2*(y-1)) = x^2*y-x^2+2*x*y-2*x+y-1.
        (
            "Expand[Times[Power[Plus[x, 1], 2], Plus[y, -1]]]",
            "Plus[-1, Times[-2, x], Times[-1, Power[x, 2]], y, Times[2, x, y], Times[Power[x, 2], y]]",
        ),
        // Inside function arguments: sin(x^2+2*x+1).
        (
            "Expand[Sin[Power[Plus[x, 1], 2]]]",
            "Sin[Plus[1, Times[2, x], Power[x, 2]]]",
        ),
        // Each term of the numerator over the expanded denominator.
        (
            "Expand[Times[Power[Plus[x, 1], 2], Power[Plus[x, -1], -2]]]",
            "Plus[Power[Plus[1, Times[-2, x], Power[x, 2]], -1], Times[2, x, Power[Plus[1, Times[-2, x], Power[x, 2]], -1]], Times[Power[x, 2], Power[Plus[1, Times[-2, x], Power[x, 2]], -1]]]",
        ),
        // Over a monomial: 1/2*x/y + 1/2/y.
        (
            "Expand[Times[Plus[x, 1], Power[Times[2, y], -1]]]",
            "Plus[Times[1/2, Power[y, -1]], Times[1/2, x, Power[y, -1]]]",
        ),
        ("Expand[Times[2, Plus[x, 1]]]", "Plus[2, Times[2, x]]"),
        // Algebra only: exp(a+b) stays.
        ("Expand[Exp[Plus[a, b]]]", "Exp[Plus[a, b]]"),
        (
            "Expand[Times[Power[Plus[a, b], 2], Sin[x]]]",
            "Plus[Times[Power[a, 2], Sin[x]], Times[2, a, b, Sin[x]], Times[Power[b, 2], Sin[x]]]",
        ),
        (
            "Expand[Power[Plus[x, 1], -2]]",
            "Power[Plus[1, Times[2, x], Power[x, 2]], -1]",
        ),
        ("Expand[5/4]", "5/4"),
        // An infinity is its own expansion and normal form.
        ("Expand[ComplexInfinity]", "ComplexInfinity"),
        ("Together[Times[-1, Infinity]]", "DirectedInfinity[-1]"),
        ("Together[Indeterminate]", "Indeterminate"),
        ("Expand[x]", "x"),
        (
            "Expand[List[Power[Plus[x, 1], 2], 3]]",
            "List[Plus[1, Times[2, x], Power[x, 2]], 3]",
        ),
        // normal(1/x+1/(x+1)) = (2*x+1)/(x^2+x).
        (
            "Together[Plus[Power[x, -1], Power[Plus[x, 1], -1]]]",
            "Times[Plus[1, Times[2, x]], Power[Plus[x, Power[x, 2]], -1]]",
        ),
        // (x^2-1)/(x^2+2*x+1) = (x-1)/(x+1).
        (
            "Together[Times[Plus[Power[x, 2], -1], Power[Plus[Power[x, 2], Times[2, x], 1], -1]]]",
            "Times[Plus[-1, x], Power[Plus[1, x], -1]]",
        ),
        // 1/(1-x) = -1/(x-1): the denominator's leading coefficient positive.
        (
            "Together[Power[Plus[1, Times[-1, x]], -1]]",
            "Times[-1, Power[Plus[-1, x], -1]]",
        ),
        // (x^2-a^2)/(x-a) = a+x.
        (
            "Together[Times[Plus[Power[x, 2], Times[-1, Power[a, 2]]], Power[Plus[x, Times[-1, a]], -1]]]",
            "Plus[a, x]",
        ),
        // A constant denominator stays in the coefficients: 1/2*x+1/3*y.
        (
            "Together[Plus[Times[1/2, x], Times[1/3, y]]]",
            "Plus[Times[1/2, x], Times[1/3, y]]",
        ),
        // 2/(4*x+4) = 1/(2*x+2); (x+1)/(2*x+2) = 1/2.
        (
            "Together[Times[2, Power[Plus[Times[4, x], 4], -1]]]",
            "Power[Plus[2, Times[2, x]], -1]",
        ),
        (
            "Together[Times[Plus[x, 1], Power[Plus[Times[2, x], 2], -1]]]",
            "1/2",
        ),
        (
            "Together[Plus[Power[Sin[x], 2], Times[2, Sin[x]], 1]]",
            "Plus[1, Times[2, Sin[x]], Power[Sin[x], 2]]",
        ),
        // Over ℤ[i] the denominator's leading coefficient is in the first
        // quadrant: (1+i)/(2x+4) = i/((1+i)x + 2 + 2i), where giac leaves
        // 1/((1-i)x + 2 - 2i).
        (
            "Together[Times[Complex[1/2, 1/2], Power[Plus[2, x], -1]]]",
            "Times[Complex[0, 1], Power[Plus[Complex[2, 2], Times[Complex[1, 1], x]], -1]]",
        ),
    ]);
}

/// `Expand` and `Together` stay out of the arguments a head holds, by the
/// context's attributes: `Hold`'s, and a user's `HoldAll`.
#[test]
fn held_arguments_stay() {
    let mut ctx = Context::new();
    ctx.set_attributes(Symbol::new("f"), Attributes::HOLD_ALL)
        .unwrap();
    for (input, expected) in [
        ("Expand[f[Power[Plus[x, 1], 2]]]", "f[Power[Plus[x, 1], 2]]"),
        (
            "Expand[Hold[Power[Plus[x, 1], 2]]]",
            "Hold[Power[Plus[x, 1], 2]]",
        ),
        (
            "Expand[g[Power[Plus[x, 1], 2]]]",
            "g[Plus[1, Times[2, x], Power[x, 2]]]",
        ),
        (
            "Together[f[Plus[Power[x, -1], 1]]]",
            "f[Plus[Power[x, -1], 1]]",
        ),
        (
            "Together[g[Plus[Power[x, -1], 1]]]",
            "g[Times[Power[x, -1], Plus[1, x]]]",
        ),
    ] {
        let (found, warnings) = eval(&ctx, input);
        assert_eq!(format!("{found:?}"), expected, "{input}");
        assert_eq!(warnings, Vec::<String>::new(), "{input}");
    }
}

/// giac's `expand` of fractions (P6), term by term: each power with a
/// negative exponent is the denominator of its terms, those of one term
/// multiplied into one and expanded; terms never go over a common
/// denominator.
#[test]
fn expand_term_by_term() {
    check(&[
        // expand(1/x+1/(x+1)) = 1/x+1/(x+1).
        (
            "Expand[Plus[Power[x, -1], Power[Plus[x, 1], -1]]]",
            "Plus[Power[x, -1], Power[Plus[1, x], -1]]",
        ),
        // (1+1/(x+1))*(x+2) = x+x/(x+1)+2+2/(x+1).
        (
            "Expand[Times[Plus[1, Power[Plus[x, 1], -1]], Plus[x, 2]]]",
            "Plus[2, x, Times[2, Power[Plus[1, x], -1]], Times[x, Power[Plus[1, x], -1]]]",
        ),
        // x/(1-y) = x/(-y+1): the denominator as written.
        (
            "Expand[Times[x, Power[Plus[1, Times[-1, y]], -1]]]",
            "Times[x, Power[Plus[1, Times[-1, y]], -1]]",
        ),
        // (x+1)/(x*(x+2)) = x/(x^2+2*x)+1/(x^2+2*x).
        (
            "Expand[Times[Plus[x, 1], Power[Times[x, Plus[x, 2]], -1]]]",
            "Plus[Power[Plus[Times[2, x], Power[x, 2]], -1], Times[x, Power[Plus[Times[2, x], Power[x, 2]], -1]]]",
        ),
        // (1+1/x)^2 = (1/x)^2+1+2/x; x^(-2)*(x+1)^2 = 1+2/x+1/x^2.
        (
            "Expand[Power[Plus[1, Power[x, -1]], 2]]",
            "Plus[1, Power[x, -2], Times[2, Power[x, -1]]]",
        ),
        (
            "Expand[Times[Power[x, -2], Power[Plus[x, 1], 2]]]",
            "Plus[1, Power[x, -2], Times[2, Power[x, -1]]]",
        ),
        // 1/(1+1/x) stays.
        (
            "Expand[Power[Plus[1, Power[x, -1]], -1]]",
            "Power[Plus[1, Power[x, -1]], -1]",
        ),
        // (x+1)^2/(2*y) = 1/2*x^2/y+x/y+1/2/y; 3/(4*x+4) stays.
        (
            "Expand[Times[Power[Plus[x, 1], 2], Power[Times[2, y], -1]]]",
            "Plus[Times[1/2, Power[y, -1]], Times[x, Power[y, -1]], Times[1/2, Power[x, 2], Power[y, -1]]]",
        ),
        (
            "Expand[Times[3, Power[Plus[Times[4, x], 4], -1]]]",
            "Times[3, Power[Plus[4, Times[4, x]], -1]]",
        ),
        // exp(-x)*(x+1) = x*exp(-x)+exp(-x); (x+1)/sqrt(x) = sqrt(x)+1/sqrt(x).
        (
            "Expand[Times[Power[E, Times[-1, x]], Plus[x, 1]]]",
            "Plus[Power[E, Times[-1, x]], Times[Power[E, Times[-1, x]], x]]",
        ),
        (
            "Expand[Times[Plus[x, 1], Power[x, -1/2]]]",
            "Plus[Power[x, -1/2], Power[x, 1/2]]",
        ),
        // (x+1)^(3/2) = x*sqrt(x+1)+sqrt(x+1).
        (
            "Expand[Power[Plus[x, 1], 3/2]]",
            "Plus[Power[Plus[1, x], 1/2], Times[x, Power[Plus[1, x], 1/2]]]",
        ),
        // 1/((x+1)^2*(y+1)) = 1/(x^2*y+x^2+2*x*y+2*x+y+1).
        (
            "Expand[Power[Times[Power[Plus[x, 1], 2], Plus[y, 1]], -1]]",
            "Power[Plus[1, Times[2, x], Power[x, 2], y, Times[2, x, y], Times[Power[x, 2], y]], -1]",
        ),
        // 1/((x+1)*(x-1)) = 1/(x^2-1), which is 1/(x+1)/(x-1) here.
        (
            "Expand[Power[Times[Plus[x, 1], Plus[x, -1]], -1]]",
            "Power[Plus[-1, Power[x, 2]], -1]",
        ),
        // (x+y)^2/(x-y) = x^2/(x-y)+y^2/(x-y)+2*x*y/(x-y).
        (
            "Expand[Times[Power[Plus[x, y], 2], Power[Plus[x, Times[-1, y]], -1]]]",
            "Plus[Times[Power[x, 2], Power[Plus[x, Times[-1, y]], -1]], Times[2, x, y, Power[Plus[x, Times[-1, y]], -1]], Times[Power[y, 2], Power[Plus[x, Times[-1, y]], -1]]]",
        ),
        // Ours: one denominator a term, expanded, where giac keeps
        // (1/(x+1))^2 from the square.
        (
            "Expand[Power[Plus[1, Power[Plus[x, 1], -1]], 2]]",
            "Plus[1, Times[2, Power[Plus[1, x], -1]], Power[Plus[1, Times[2, x], Power[x, 2]], -1]]",
        ),
        // Balls are coefficients: (x + 1.5)^2.
        (
            "Expand[Power[Plus[x, 3/2], 2]]",
            "Plus[9/4, Times[3, x], Power[x, 2]]",
        ),
    ]);
    let ctx = Context::new();
    let (e, warnings) = eval(&ctx, "Expand[Power[Plus[x, N[3/2]], 2]]");
    assert_eq!(warnings, Vec::<String>::new());
    let [constant, linear, square] = e.args() else {
        panic!("{e:?}");
    };
    assert!(constant.is_number() && !constant.is_exact_number(), "{e:?}");
    assert!(
        linear.args()[0].is_number() && !linear.args()[0].is_exact_number(),
        "{e:?}"
    );
    assert_eq!(*square, read("Power[x, 2]"));
}

/// giac's `gcd` and `lcm`, with §15 detail 5's normal forms.
#[test]
fn gcd_and_lcm() {
    check(&[
        // gcd(6x^2-6, 4x-4) = 2x-2, either sign.
        (
            "GCD[Plus[Times[6, Power[x, 2]], -6], Plus[Times[4, x], -4]]",
            "Plus[-2, Times[2, x]]",
        ),
        (
            "GCD[Plus[Times[-6, Power[x, 2]], 6], Plus[Times[-4, x], 4]]",
            "Plus[-2, Times[2, x]]",
        ),
        // Over ℚ, the gcd of the numerators: 10x-1, 3x-3, 2x-2.
        (
            "GCD[Plus[Power[x, 2], -1/100], Plus[x, -1/10]]",
            "Plus[-1, Times[10, x]]",
        ),
        (
            "GCD[Plus[Times[3/2, x], -3/2], Plus[Times[3, Power[x, 2]], -3]]",
            "Plus[-3, Times[3, x]]",
        ),
        (
            "GCD[Plus[Times[2, x], -2], Plus[Times[4/3, Power[x, 2]], -4/3]]",
            "Plus[-2, Times[2, x]]",
        ),
        // Several variables, three arguments, a number.
        (
            "GCD[Plus[Times[x, y], Times[-1, y]], Plus[Times[Power[x, 2], y], Times[-1, y]]]",
            "Plus[Times[-1, y], Times[x, y]]",
        ),
        (
            "GCD[Plus[Power[x, 2], -1], Plus[Power[x, 2], Times[2, x], 1], Plus[x, 1]]",
            "Plus[1, x]",
        ),
        ("GCD[4, Times[6, x]]", "2"),
        // Numbers as numbers: gcd(1/2, 1/3) = 1/6, gcd(2/3, 4/9) = 2/9.
        ("GCD[1/2, 1/3]", "1/6"),
        ("GCD[2/3, 4/9]", "2/9"),
        ("GCD[6, 4, 10]", "2"),
        ("GCD[Complex[0, 2], 4]", "2"),
        // Zeros: gcd(0, b) normalised, where giac returns b.
        ("GCD[0, Plus[Times[-2, x], 2]]", "Plus[-2, Times[2, x]]"),
        ("GCD[0, 0]", "0"),
        // lcm(6x^2-6, 4x-4) = 12x^2-12, expanded; lcm(x^2-1, x^2+2x+1).
        (
            "LCM[Plus[Times[6, Power[x, 2]], -6], Plus[Times[4, x], -4]]",
            "Plus[-12, Times[12, Power[x, 2]]]",
        ),
        (
            "LCM[Plus[Power[x, 2], -1], Plus[Power[x, 2], Times[2, x], 1]]",
            "Plus[-1, Times[-1, x], Power[x, 2], Power[x, 3]]",
        ),
        // lcm(x/2, x^2/3) = x^2/6; of numbers, lcm(1/2, 1/3) = 1.
        (
            "LCM[Times[1/2, x], Times[1/3, Power[x, 2]]]",
            "Times[1/6, Power[x, 2]]",
        ),
        ("LCM[1/2, 1/3]", "1"),
        ("LCM[4, 6]", "12"),
        // Normalised, where giac keeps signs: lcm(-4, 6) = -12 and
        // lcm(1-x, 1+x) = (-x+1)*(1+x).
        ("LCM[-4, 6]", "12"),
        ("GCD[-4, -6]", "2"),
        ("GCD[-4]", "4"),
        ("LCM[-6]", "6"),
        (
            "LCM[Plus[1, Times[-1, x]], Plus[1, x]]",
            "Plus[-1, Power[x, 2]]",
        ),
        ("LCM[Plus[1, Times[-1, x]], 0]", "0"),
    ]);
}

/// giac's `quo`, `rem`, `quorem`, over the fractions of the other
/// generators.
#[test]
fn division() {
    check(&[
        (
            "PolynomialQuotient[Plus[Power[x, 3], 1], Plus[x, 1]]",
            "Plus[1, Times[-1, x], Power[x, 2]]",
        ),
        ("PolynomialRemainder[Plus[Power[x, 3], 2], Plus[x, 1]]", "1"),
        (
            "PolynomialQuotientRemainder[Plus[Power[x, 3], 2], Plus[x, 1]]",
            "List[Plus[1, Times[-1, x], Power[x, 2]], 1]",
        ),
        // quo(x^2*y+1, x+y, x) = y*x-y^2; rem = y^3+1.
        (
            "PolynomialQuotient[Plus[Times[Power[x, 2], y], 1], Plus[x, y], x]",
            "Plus[Times[x, y], Times[-1, Power[y, 2]]]",
        ),
        (
            "PolynomialRemainder[Plus[Times[Power[x, 2], y], 1], Plus[x, y], x]",
            "Plus[1, Power[y, 3]]",
        ),
        // quo(2x^2+1, 3x+1) = 2/3*x-2/9, over ℚ.
        (
            "PolynomialQuotient[Plus[Times[2, Power[x, 2]], 1], Plus[Times[3, x], 1]]",
            "Plus[-2/9, Times[2/3, x]]",
        ),
        // quo(x^2, y*x+1, x) = 1/y*x-1/y^2; rem = 1/y^2.
        (
            "PolynomialQuotient[Power[x, 2], Plus[Times[y, x], 1], x]",
            "Plus[Times[-1, Power[y, -2]], Times[x, Power[y, -1]]]",
        ),
        (
            "PolynomialRemainder[Power[x, 2], Plus[Times[y, x], 1], x]",
            "Power[y, -2]",
        ),
        ("PolynomialQuotient[x, Power[x, 2]]", "0"),
        ("PolynomialRemainder[5, Plus[x, 1]]", "5"),
        (
            "PolynomialQuotientRemainder[Plus[Power[x, 2], 1], Plus[Power[x, 2], 1]]",
            "List[1, 0]",
        ),
    ]);
    // By 0, as Quotient of numbers.
    let ctx = Context::new();
    let (q, warnings) = eval(&ctx, "PolynomialQuotient[Plus[Power[x, 2], 1], 0]");
    assert_eq!(format!("{q:?}"), "ComplexInfinity");
    assert_eq!(warnings.len(), 1);
    let (r, _) = eval(&ctx, "PolynomialRemainder[Plus[Power[x, 2], 1], 0]");
    assert_eq!(format!("{r:?}"), "Indeterminate");
}

/// giac's `degree`, `coeff`, `lcoeff`, with §15 details 6 and 7.
#[test]
fn degrees_and_coefficients() {
    check(&[
        ("Exponent[Plus[Power[x, 3], Times[x, Power[y, 5]]], x]", "3"),
        ("Exponent[Plus[Power[x, 3], Times[x, Power[y, 5]]], y]", "5"),
        ("Exponent[Plus[Power[x, 3], x]]", "3"),
        ("Exponent[5, x]", "0"),
        // degree(0) is -infinity, where giac gives 0.
        ("Exponent[0, x]", "DirectedInfinity[-1]"),
        // In sin(x), where giac gives 1; constants don't count.
        ("Exponent[Plus[Power[Sin[x], 2], Sin[x]]]", "2"),
        ("Exponent[Times[Pi, Power[x, 2]]]", "2"),
        // Ball coefficients.
        ("Exponent[Times[N[3/2], Power[x, 2]], x]", "2"),
        ("Coefficient[Plus[Power[x, 3], Times[2, x], 5], x, 1]", "2"),
        ("Coefficient[Plus[Power[x, 3], Times[2, x], 5], x, 0]", "5"),
        ("Coefficient[Plus[Power[x, 3], Times[2, x], 5], 1]", "2"),
        (
            "Coefficient[Plus[Times[Power[x, 2], y], Times[3, x, Power[y, 2]]], x, 1]",
            "Times[3, Power[y, 2]]",
        ),
        ("Coefficient[Plus[Times[y, Power[x, 2]], 3], x, 2]", "y"),
        // The list of coefficients, as giac's poly1[1,3,0].
        (
            "Coefficient[Plus[Power[x, 2], Times[3, x]], x]",
            "List[1, 3, 0]",
        ),
        ("Coefficient[Plus[Power[x, 2], Times[3, x]], x, 7]", "0"),
        ("LeadingCoefficient[Plus[Times[3, Power[x, 2]], x]]", "3"),
        ("LeadingCoefficient[Plus[Times[-3, Power[x, 2]], x]]", "-3"),
        ("LeadingCoefficient[Plus[Times[y, Power[x, 2]], x], x]", "y"),
        (
            "LeadingCoefficient[Plus[Times[y, Power[x, 2]], x], y]",
            "Power[x, 2]",
        ),
        ("LeadingCoefficient[0]", "0"),
        ("LeadingCoefficient[5]", "5"),
    ]);
    let ctx = Context::new();
    let (c, _) = eval(&ctx, "Coefficient[Times[N[3/2], x], x, 1]");
    assert!(c.is_number() && !c.is_exact_number(), "{c:?}");
}

/// giac's `numer`, `denom`, `content`, `primpart`.
#[test]
fn numerators_and_contents() {
    check(&[
        ("Numerator[Times[x, Power[Plus[x, 1], -1]]]", "x"),
        ("Denominator[Times[x, Power[Plus[x, 1], -1]]]", "Plus[1, x]"),
        ("Numerator[Plus[1/2, x]]", "Plus[1, Times[2, x]]"),
        ("Denominator[Plus[1/2, x]]", "2"),
        ("Numerator[3/4]", "3"),
        ("Denominator[3/4]", "4"),
        ("Numerator[Plus[Times[x, Power[y, -1]], 1]]", "Plus[x, y]"),
        ("Denominator[Plus[Times[x, Power[y, -1]], 1]]", "y"),
        // Cancelled, the denominator's sign moved up.
        (
            "Numerator[Times[Plus[Power[x, 2], -1], Power[Plus[x, -1], -1]]]",
            "Plus[1, x]",
        ),
        (
            "Denominator[Times[Plus[Power[x, 2], -1], Power[Plus[x, -1], -1]]]",
            "1",
        ),
        (
            "Numerator[Times[x, Power[Plus[1, Times[-1, x]], -1]]]",
            "Times[-1, x]",
        ),
        (
            "Denominator[Times[x, Power[Plus[1, Times[-1, x]], -1]]]",
            "Plus[-1, x]",
        ),
        ("Content[Plus[Times[6, Power[x, 2]], Times[4, x]]]", "2"),
        ("Content[Plus[Times[-6, Power[x, 2]], Times[-4, x]]]", "2"),
        (
            "Content[Plus[Times[6, Power[x, 2], y], Times[4, x, y]], x]",
            "Times[2, y]",
        ),
        (
            "Content[Plus[Times[1/2, Power[x, 2]], Times[1/3, x]]]",
            "1/6",
        ),
        ("Content[0]", "0"),
        (
            "PrimitivePart[Plus[Times[6, Power[x, 2]], Times[4, x]]]",
            "Plus[Times[2, x], Times[3, Power[x, 2]]]",
        ),
        // The sign stays in the primitive part, as giac.
        (
            "PrimitivePart[Plus[Times[-6, Power[x, 2]], Times[-4, x]]]",
            "Plus[Times[-2, x], Times[-3, Power[x, 2]]]",
        ),
        (
            "PrimitivePart[Plus[Times[1/2, Power[x, 2]], Times[1/3, x]]]",
            "Plus[Times[2, x], Times[3, Power[x, 2]]]",
        ),
        (
            "PrimitivePart[Plus[Times[6, Power[x, 2], y], Times[4, x, y]], x]",
            "Plus[Times[2, x], Times[3, Power[x, 2]]]",
        ),
    ]);
}

/// giac's resultants and discriminants (CAS-137), but a constant's.
#[test]
fn resultants_and_discriminants() {
    check(&[
        (
            "Resultant[Plus[Power[x, 2], -2], Plus[Power[x, 3], -3], x]",
            "1",
        ),
        (
            "Resultant[Plus[Times[x, y], -1], Plus[Power[x, 2], Power[y, 2], -4], x]",
            "Plus[1, Times[-4, Power[y, 2]], Power[y, 4]]",
        ),
        (
            "Resultant[Plus[Power[x, 3], x], Plus[Power[x, 2], 3]]",
            "12",
        ),
        (
            "Resultant[Plus[Times[2, Power[x, 2]], 3], Plus[Times[4, Power[x, 3]], -1]]",
            "440",
        ),
        // A constant to the other's degree, where giac gives 1.
        ("Resultant[2, Plus[Power[x, 2], 1], x]", "4"),
        ("Resultant[0, Plus[x, 1], x]", "0"),
        // Two constants: the empty determinant, 1, but 0 for a zero.
        ("Resultant[2, 3, x]", "1"),
        ("Resultant[0, 3, x]", "0"),
        // Over ℚ: res(x/2 - 1, x - 3) = 1/2 b(2) = -1/2.
        ("Resultant[Plus[Times[1/2, x], -1], Plus[x, -3], x]", "-1/2"),
        (
            "Discriminant[Plus[Times[a, Power[x, 2]], Times[b, x], c], x]",
            "Plus[Power[b, 2], Times[-4, a, c]]",
        ),
        (
            "Discriminant[Plus[Power[x, 3], Times[p, x], q], x]",
            "Plus[Times[-4, Power[p, 3]], Times[-27, Power[q, 2]]]",
        ),
        ("Discriminant[Plus[Power[x, 4], 1]]", "256"),
        (
            "Discriminant[Plus[Times[2, Power[x, 2]], Times[3, x], 1]]",
            "1",
        ),
        // Degree 1 gives 1, where giac gives 1/3 for 3x+1.
        ("Discriminant[Plus[Times[3, x], 1]]", "1"),
        ("Discriminant[Plus[Times[1/2, Power[x, 2]], 1]]", "-2"),
    ]);
}

/// The calls that stay, with their warnings.
#[test]
fn refusals() {
    check_refused(&[
        // Two candidates: the variable must be given (P5).
        (
            "PolynomialQuotient[Plus[Times[Power[x, 2], y], 1], Plus[x, y]]",
            "the variable must be given: one of x, y",
        ),
        ("Exponent[Plus[x, y]]", "the variable must be given"),
        (
            "Exponent[Plus[x, Power[x, 1/2]]]",
            "the variable must be given",
        ),
        // Not a polynomial in the variable (§15 detail 6).
        ("Exponent[Sin[x], x]", "not a polynomial in x"),
        (
            "Exponent[Times[Power[x, 2], Power[Plus[x, 1], -1]], x]",
            "not a polynomial in x",
        ),
        (
            "Coefficient[Plus[x, Power[x, 1/2]], x, 1]",
            "not a polynomial in x",
        ),
        ("Exponent[Plus[x, 1], 2]", "cannot be a number"),
        ("Coefficient[Plus[x, 1], x, 1/2]", "an integer"),
        // Balls in the exact algorithms (P8).
        ("GCD[Times[N[2], x], x]", "must be exact"),
        (
            "Together[Times[x, Power[Plus[x, N[1]], -1]]]",
            "must be exact",
        ),
        (
            "PolynomialQuotient[Plus[Power[x, 2], N[1]], Plus[x, 1]]",
            "must be exact",
        ),
        // Polynomials only, for gcd and lcm.
        ("GCD[Power[x, -1], x]", "must be polynomials"),
        // A discriminant needs a degree.
        ("Discriminant[5, x]", "a degree of at least 1"),
        // Infinities, even inside a generator (CAS-139).
        ("Expand[f[Infinity]]", "not polynomials"),
        ("GCD[Plus[x, 1], f[Indeterminate]]", "not polynomials"),
    ]);
}

/// Past `Limits::terms`, which bounds the terms of each product before it
/// is taken, the call stays, with a warning that gives the limit (§15
/// detail 3). The bound is the result's for a dense product, not its
/// products of terms.
#[test]
fn the_limit_on_terms() {
    let input = "Expand[Power[Plus[1, x, y], 10]]";
    let with_limit = |terms| {
        let mut ctx = Context::new();
        let mut limits = ctx.limits();
        limits.terms = terms;
        ctx.set_limits(limits);
        eval(&ctx, input)
    };
    let (found, warnings) = with_limit(65);
    assert_eq!(found, read(input), "{found:?}");
    assert!(
        warnings.iter().any(|w| w.contains("could pass 65 terms")),
        "{warnings:?}"
    );
    // (1 + x + y)^5 squared: 441 products of terms, 66 terms, the bound of
    // the total degree.
    let (found, warnings) = with_limit(66);
    assert_eq!(found.args().len(), 66);
    assert_eq!(warnings, Vec::<String>::new());
    // Sparse: (x^9 + y^9 + 1)(x + y + 1) makes 9 terms, its products.
    let input = "Expand[Times[Plus[Power[x, 9], Power[y, 9], 1], Plus[x, y, 1]]]";
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.terms = 8;
    ctx.set_limits(limits);
    let (found, warnings) = eval(&ctx, input);
    assert_eq!(found.head(), Some(&Expr::from(Symbol::EXPAND)), "{found:?}");
    assert!(
        warnings.iter().any(|w| w.contains("could pass 8 terms")),
        "{warnings:?}"
    );
    limits.terms = 9;
    ctx.set_limits(limits);
    assert_eq!(eval(&ctx, input).0.args().len(), 9);
    // Within the default, a product of 1001-term powers: 10626 terms where
    // its products of terms pass 10^6.
    let (found, warnings) = eval(
        &Context::new(),
        "Expand[Times[Power[Plus[1, x, y, z, t], 10], Plus[Power[Plus[1, x, y, z, t], 10], 1]]]",
    );
    assert_eq!(found.args().len(), 10626);
    assert_eq!(warnings, Vec::<String>::new());
}

/// The interrupt stops a polynomial builtin, raised before it starts.
#[test]
fn the_interrupt_stops_them() {
    let mut ctx = Context::new();
    let interrupt = Interrupt::new();
    ctx.set_interrupt(interrupt.clone());
    interrupt.interrupt();
    for input in [
        "Expand[Power[Plus[1, x, y], 10]]",
        "GCD[Plus[Power[x, 2], -1], Plus[x, -1]]",
    ] {
        assert_eq!(
            ctx.eval(&read(input)),
            Err(EvalError::Interrupted),
            "{input}"
        );
    }
}

/// The interrupt stops a polynomial builtin inside its products: the last
/// squaring of `(1 + x + y + z + w + v)^40` takes `53130^2` products of
/// terms, minutes of work, and the interrupt comes 200 ms in.
#[test]
fn the_interrupt_stops_a_long_product() {
    let mut ctx = Context::new();
    let mut limits = ctx.limits();
    limits.terms = usize::MAX;
    ctx.set_limits(limits);
    let interrupt = Interrupt::new();
    ctx.set_interrupt(interrupt.clone());
    let raise = std::thread::spawn({
        let interrupt = interrupt.clone();
        move || {
            std::thread::sleep(Duration::from_millis(200));
            interrupt.interrupt();
        }
    });
    let start = Instant::now();
    let result = ctx.eval(&read("Expand[Power[Plus[1, x, y, z, w, v], 40]]"));
    raise.join().unwrap();
    assert_eq!(result, Err(EvalError::Interrupted));
    assert!(
        start.elapsed() < Duration::from_secs(30),
        "{:?}",
        start.elapsed()
    );
}

/// Every builtin's trace verifies: `Expand` and `Together` steps agree at
/// the sample points, and the formal ones are not checked.
#[test]
fn traces_verify() {
    let ctx = Context::new();
    let inputs = [
        "Expand[Times[Power[Plus[x, 1], 2], Plus[y, -1]]]",
        "Together[Plus[Power[x, -1], Power[Plus[x, 1], -1]]]",
        "GCD[Plus[Power[x, 2], -1], Plus[x, -1]]",
        "LCM[Plus[Power[x, 2], -1], Plus[x, -1]]",
        "PolynomialQuotient[Plus[Power[x, 3], 1], Plus[x, 1]]",
        "PolynomialRemainder[Plus[Power[x, 3], 2], Plus[x, 1]]",
        "PolynomialQuotientRemainder[Plus[Power[x, 3], 2], Plus[x, 1]]",
        "Exponent[Plus[Power[x, 3], x]]",
        "Coefficient[Plus[Power[x, 3], Times[2, x], 5], x, 1]",
        "LeadingCoefficient[Plus[Times[3, Power[x, 2]], x]]",
        "Numerator[Times[x, Power[Plus[x, 1], -1]]]",
        "Denominator[Times[x, Power[Plus[x, 1], -1]]]",
        "Content[Plus[Times[6, Power[x, 2], y], Times[4, x, y]], x]",
        "PrimitivePart[Plus[Times[6, Power[x, 2]], Times[4, x]]]",
        "Resultant[Plus[Times[x, y], -1], Plus[Power[x, 2], Power[y, 2], -4], x]",
        "Discriminant[Plus[Times[a, Power[x, 2]], Times[b, x], c], x]",
    ];
    for input in inputs {
        let options = EvalOptions {
            trace: TraceLevel::Rules,
        };
        let (result, report) = ctx.eval_with(&read(input), &options);
        let result = result.unwrap();
        let verification = report.trace.verify(&ctx);
        assert!(
            verification.wrong.is_none(),
            "{input}: {:?}",
            verification.wrong
        );
        // One step of the call's algorithm (CAS-DOC-8 §11.5).
        let algorithm: Vec<&Step> = report
            .trace
            .steps
            .iter()
            .filter(|s| matches!(s.kind, StepKind::Algorithm(_)))
            .collect();
        assert_eq!(algorithm.len(), 1, "{input}: {algorithm:?}");
        assert_eq!(algorithm[0].before.head(), read(input).head(), "{input}");
        assert_eq!(ctx.eval(&algorithm[0].after), Ok(result), "{input}");
    }
    // Expand's and Together's steps are checked.
    for input in &inputs[..2] {
        let options = EvalOptions {
            trace: TraceLevel::Rules,
        };
        let (_, report) = ctx.eval_with(&read(input), &options);
        assert!(report.trace.verify(&ctx).passed > 0, "{input}");
    }
}

/// Terms of `n` variables as exponent vectors and small coefficients.
type TermList = Vec<(Vec<u64>, i64)>;

fn terms(n: usize, len: usize) -> impl Strategy<Value = TermList> {
    prop::collection::vec((prop::collection::vec(0u64..4, n), -9i64..9), 0..=len)
}

/// The expression of a polynomial in `x`, `y`, `z`, built, not evaluated.
/// As [`terms`], at least one, each coefficient nonzero.
fn nonzero_terms(n: usize, len: usize) -> impl Strategy<Value = TermList> {
    prop::collection::vec(
        (
            prop::collection::vec(0u64..3, n),
            prop_oneof![-9i64..=-1, 1i64..=9],
        ),
        1..=len,
    )
}

/// Whether `e` is a sum of terms without sums.
fn expanded(e: &Expr) -> bool {
    let plus = Expr::from(Symbol::PLUS);
    let terms = if e.head() == Some(&plus) {
        e.args()
    } else {
        std::slice::from_ref(e)
    };
    terms.iter().all(|t| !t.contains(&plus))
}

/// Whether `e` is as `Together` writes it: a product of monomial factors,
/// at most one expanded sum, and at most one expanded sum to the power -1.
fn together_form(e: &Expr) -> bool {
    let plus = Expr::from(Symbol::PLUS);
    if e.head() == Some(&plus) {
        return expanded(e);
    }
    let factors = if e.head() == Some(&Expr::from(Symbol::TIMES)) {
        e.args()
    } else {
        std::slice::from_ref(e)
    };
    let (mut sums, mut inverses) = (0, 0);
    let fine = factors.iter().all(|f| {
        if f.head() == Some(&plus) {
            sums += 1;
            return expanded(f);
        }
        if let [b, k] = f.args()
            && f.head() == Some(&Expr::from(Symbol::POWER))
            && b.head() == Some(&plus)
        {
            inverses += 1;
            return *k == Expr::from(-1) && expanded(b);
        }
        !f.contains(&plus)
    });
    fine && sums <= 1 && inverses <= 1
}

fn poly_expr(t: &[(Vec<u64>, i64)]) -> Expr {
    let vars = [Expr::symbol("x"), Expr::symbol("y"), Expr::symbol("z")];
    t.iter()
        .map(|(e, c)| {
            e.iter()
                .zip(&vars)
                .map(|(&k, v)| v.pow(Expr::from(Integer::from(k))))
                .fold(Expr::from(*c), |acc, f| acc * f)
        })
        .sum()
}

fn sparse(t: &[(Vec<u64>, i64)]) -> Sparse<Integer> {
    SparseRing::new(Integers, 3).from_terms(t.iter().map(|(e, c)| (e.clone(), Integer::from(*c))))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `Expand` of a product is cas-poly's product of the factors.
    #[test]
    fn expand_is_the_product(a in terms(3, 6), b in terms(3, 6)) {
        let ctx = Context::new();
        let product = Expr::call(Symbol::EXPAND, [poly_expr(&a) * poly_expr(&b)]);
        let expanded = ctx.eval(&product).unwrap();
        let p = Polynomials::new(std::slice::from_ref(&expanded), usize::MAX, &Never).unwrap();
        // The generators of the expansion are among x, y, z: reorder its
        // terms into x, y, z.
        let (num, den) = match &p.fractions()[0] {
            Fraction::Integers { num, den } => (num.clone(), den.clone()),
            Fraction::GaussianIntegers { .. } => panic!("integer coefficients"),
        };
        prop_assert!(Integers.is_one(&den.coeffs()[0]) && den.term_count() == 1);
        // Expanded: a sum of products of numbers and powers of x, y, z.
        let terms_of = |e: &Expr| -> Vec<Expr> {
            if e.head() == Some(&Expr::from(Symbol::PLUS)) {
                e.args().to_vec()
            } else if *e == Expr::from(0) {
                Vec::new()
            } else {
                vec![e.clone()]
            }
        };
        let monomial = |t: &Expr| !t.contains(&Expr::from(Symbol::PLUS));
        prop_assert!(terms_of(&expanded).iter().all(monomial), "{}", expanded);
        prop_assert_eq!(terms_of(&expanded).len(), num.term_count());
        let terms_xyz: Vec<(Vec<u64>, Integer)> = num
            .terms()
            .map(|(e, c)| {
                let at = |v| p.variable(&Expr::symbol(v)).map_or(0, |k| e[k]);
                (vec![at("x"), at("y"), at("z")], c.clone())
            })
            .collect();
        let ring = SparseRing::new(Integers, 3);
        prop_assert_eq!(ring.from_terms(terms_xyz), ring.mul(&sparse(&a), &sparse(&b)));
    }

    /// `Together` of a sum of fractions keeps its value at rational points.
    #[test]
    fn together_keeps_the_value(
        fractions in prop::collection::vec((terms(3, 3), nonzero_terms(3, 3)), 1..4),
        point in prop::collection::vec((-4i64..5, 1i64..4), 3),
    ) {
        let ctx = Context::new();
        let sum: Expr = fractions
            .iter()
            .map(|(n, d)| poly_expr(n) * poly_expr(d).pow(-1))
            .sum();
        let Ok(sum) = ctx.eval(&sum) else { return Ok(()); };
        let Ok(together) = ctx.eval(&Expr::call(Symbol::TOGETHER, [sum.clone()])) else {
            return Ok(());
        };
        // A denominator whose terms cancel.
        let infinite = [Symbol::COMPLEX_INFINITY, Symbol::INDETERMINATE];
        if infinite.iter().any(|&s| sum.contains(&Expr::from(s))) {
            if infinite.iter().any(|&s| sum == Expr::from(s)) {
                prop_assert_eq!(together, sum);
            }
            return Ok(());
        }
        prop_assert!(together_form(&together), "{} is {}", sum, together);
        // Cancelled: the numerator and denominator are coprime.
        let p = Polynomials::new(std::slice::from_ref(&together), usize::MAX, &Never).unwrap();
        let Fraction::Integers { num, den } = &p.fractions()[0] else {
            panic!("integer coefficients");
        };
        let ring = SparseRing::new(Integers, p.generators().len());
        let g = ring.gcd(num, den);
        prop_assert!(g.degrees().iter().all(|&d| d == 0), "{} has {:?}", together, g);
        let values: Vec<(Expr, Expr)> = ["x", "y", "z"]
            .iter()
            .zip(&point)
            .map(|(v, (n, d))| (Expr::symbol(v), Expr::from(Rational::new(Integer::from(*n), Integer::from(*d)))))
            .collect();
        let at = |e: &Expr| ctx.eval(&e.substitute(&values)).unwrap();
        let (u, v) = (at(&sum), at(&together));
        if u.is_exact_number() && v.is_exact_number() {
            prop_assert_eq!(u, v, "{} and {}", sum, together);
        }
    }
}
