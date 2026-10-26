//! The polynomial builtins against giac (CAS-140, CAS-185; CAS-DOC-8 §4,
//! §7–§10, §15): giac's answers to each, read back; the departures, each as
//! decided; and a few hundred random inputs a run to `expand`, `normal` and
//! `gcd`.

use std::collections::BTreeMap;

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac 1.9's answers, as it prints them, to every builtin where this CAS
/// follows it: read back and evaluated, each is this CAS's value of the
/// input.
#[test]
#[allow(clippy::too_many_lines, reason = "a table")]
fn agreement() {
    let Some(giac) = giac() else {
        return;
    };
    let recorded = [
        // expand (P6), term by term.
        ("expand((x+1)^2*(y-1))", "x^2*y-x^2+2*x*y-2*x+y-1"),
        ("expand(sin((x+1)^2))", "sin(x^2+2*x+1)"),
        (
            "expand((x+1)^2/(x-1)^2)",
            "x^2/(x^2-2*x+1)+2*x/(x^2-2*x+1)+1/(x^2-2*x+1)",
        ),
        ("expand((x+1)/(2*y))", "1/2*x/y+1/2/y"),
        ("expand(2*(x+1))", "2*x+2"),
        ("expand(exp(a+b))", "exp(a+b)"),
        (
            "expand((a+b)^2*sin(x))",
            "a^2*sin(x)+b^2*sin(x)+2*a*b*sin(x)",
        ),
        ("expand(1/(x+1)^2)", "1/(x^2+2*x+1)"),
        ("expand(1/x+1/(x+1))", "1/x+1/(x+1)"),
        ("expand((1+1/(x+1))*(x+2))", "x+x/(x+1)+2+2/(x+1)"),
        ("expand(x/(1-y))", "x/(-y+1)"),
        ("expand((x+1)/(x*(x+2)))", "x/(x^2+2*x)+1/(x^2+2*x)"),
        ("expand((1+1/x)^2)", "(1/x)^2+1+2/x"),
        ("expand(x^(-2)*(x+1)^2)", "1+2/x+1/x^2"),
        ("expand(1/(1+1/x))", "1/(1+1/x)"),
        ("expand((x+1)^2/(2*y))", "1/2*x^2/y+x/y+1/2/y"),
        ("expand(3/(4*x+4))", "3/(4*x+4)"),
        ("expand(exp(-x)*(x+1))", "x*exp(-x)+exp(-x)"),
        ("expand((x+1)/sqrt(x))", "sqrt(x)+1/sqrt(x)"),
        ("expand((x+1)^(3/2))", "x*sqrt(x+1)+sqrt(x+1)"),
        ("expand(1/((x+1)^2*(y+1)))", "1/(x^2*y+x^2+2*x*y+2*x+y+1)"),
        ("expand(1/((x+1)*(x-1)))", "1/(x^2-1)"),
        ("expand((x+y)^2/(x-y))", "x^2/(x-y)+y^2/(x-y)+2*x*y/(x-y)"),
        // normal (P7).
        ("normal(1/x+1/(x+1))", "(2*x+1)/(x^2+x)"),
        ("normal((x^2-1)/(x^2+2*x+1))", "(x-1)/(x+1)"),
        ("normal(1/(1-x))", "-1/(x-1)"),
        ("normal((x^2-a^2)/(x-a))", "a+x"),
        ("normal(x/2+y/3)", "1/2*x+1/3*y"),
        ("normal(2/(4*x+4))", "1/(2*x+2)"),
        ("normal((x+1)/(2*x+2))", "1/2"),
        ("normal(sin(x)^2+2*sin(x)+1)", "sin(x)^2+2*sin(x)+1"),
        ("normal((x^3-1)/(x^2-1))", "(x^2+x+1)/(x+1)"),
        ("normal(1/(2*x)+1/(3*x))", "5/(6*x)"),
        ("normal(3*x/(6*x+6))", "x/(2*x+2)"),
        ("normal(1/(sqrt(2)-x))", "-1/(x-sqrt(2))"),
        // gcd and lcm (§15 detail 5); of numbers as numbers.
        ("gcd(6*x^2-6, 4*x-4)", "2*x-2"),
        ("gcd(-6*x^2+6, -4*x+4)", "2*x-2"),
        ("gcd(x^2-1/100, x-1/10)", "10*x-1"),
        ("gcd(3/2*x-3/2, 3*x^2-3)", "3*x-3"),
        ("gcd(x*y-y, x^2*y-y)", "x*y-y"),
        ("gcd(x^2-1, x^2+2*x+1, x+1)", "x+1"),
        ("gcd(4, 6*x)", "2"),
        ("gcd(1/2, 1/3)", "1/6"),
        ("gcd(2/3, 4/9)", "2/9"),
        ("gcd(6, 4, 10)", "2"),
        ("gcd(6*x*y-6*y, 4*x-4)", "2*x-2"),
        ("gcd(1/2*x+1/2, x^2-1)", "x+1"),
        ("lcm(x/2, x^2/3)", "1/2*x^2/3"),
        ("lcm(1/2, 1/3)", "1"),
        ("lcm(4, 6)", "12"),
        // quo, rem, quorem, over the fractions of the other variables.
        ("quo(x^3+1, x+1)", "x^2-x+1"),
        ("rem(x^3+2, x+1)", "1"),
        ("quorem(x^3+2, x+1)", "[x^2-x+1,1]"),
        ("quo(x^2*y+1, x+y, x)", "y*x-y^2"),
        ("rem(x^2*y+1, x+y, x)", "y^3+1"),
        ("quo(2*x^2+1, 3*x+1)", "2/3*x-2/9"),
        ("quo(x^2, y*x+1, x)", "1/y*x-1/y^2"),
        ("rem(x^2, y*x+1, x)", "1/y^2"),
        // degree, coeff, lcoeff.
        ("degree(x^3+x*y^5, x)", "3"),
        ("degree(x^3+x*y^5, y)", "5"),
        ("degree(x^3+x)", "3"),
        ("degree(pi*x^2)", "2"),
        ("coeff(x^3+2*x+5, x, 1)", "2"),
        ("coeff(x^3+2*x+5, x, 0)", "5"),
        ("coeff(x^2*y+3*x*y^2, x, 1)", "3*y^2"),
        ("coeff(x^2+3*x, x)", "poly1[1,3,0]"),
        ("coeff(x^2+3*x, x, 7)", "0"),
        ("lcoeff(3*x^2+x)", "3"),
        ("lcoeff(-3*x^2+x)", "-3"),
        ("lcoeff(y*x^2+x, x)", "y"),
        ("lcoeff(y*x^2+x, y)", "x^2"),
        // numer, denom.
        ("numer(x/(x+1))", "x"),
        ("denom(x/(x+1))", "x+1"),
        ("numer(1/2+x)", "2*x+1"),
        ("denom(1/2+x)", "2"),
        ("numer(3/4)", "3"),
        ("denom(3/4)", "4"),
        ("numer(x/y+1)", "x+y"),
        ("denom(x/y+1)", "y"),
        ("numer((x^2-1)/(x-1))", "x+1"),
        ("denom((x^2-1)/(x-1))", "1"),
        ("numer(x/(1-x))", "-x"),
        ("denom(x/(1-x))", "x-1"),
        // content, primpart.
        ("content(6*x^2+4*x)", "2"),
        ("content(-6*x^2-4*x)", "2"),
        ("content(6*x^2*y+4*x*y, x)", "2*y"),
        ("content(x^2/2+x/3)", "1/6"),
        ("primpart(6*x^2+4*x)", "3*x^2+2*x"),
        ("primpart(-6*x^2-4*x)", "-3*x^2-2*x"),
        ("primpart(x^2/2+x/3)", "3*x^2+2*x"),
        ("primpart(6*x^2*y+4*x*y, x)", "6*y/(2*y)*x^2+4*y/(2*y)*x"),
        // resultant, discriminant.
        ("resultant(x^2-2, x^3-3, x)", "1"),
        ("resultant(x*y-1, x^2+y^2-4, x)", "y^4-4*y^2+1"),
        ("resultant(x^3+x, x^2+3, x)", "12"),
        ("resultant(2*x^2+3, 4*x^3-1, x)", "440"),
        ("resultant(x/2-1, x-3, x)", "-1/2"),
        ("resultant(2, 3, x)", "1"),
        ("resultant(0, 3, x)", "0"),
        ("discriminant(a*x^2+b*x+c, x)", "-4*a*c+b^2"),
        ("discriminant(x^3+p*x+q, x)", "-4*p^3-27*q^2"),
        ("discriminant(x^4+1, x)", "256"),
        ("discriminant(2*x^2+3*x+1, x)", "1"),
        ("discriminant(x^2/2+1, x)", "-2"),
        // factor over ℚ (CAS-185), each factor with a positive leading
        // coefficient in the first variable, as giac's; quadratics in one
        // generator over square roots (CAS-187); collect over ℚ alone.
        ("factor(2*x^2-2)", "2*(x-1)*(x+1)"),
        ("factor(-x^2+1)", "-(x-1)*(x+1)"),
        ("factor(6*x+12)", "6*(x+2)"),
        ("factor(x^4-y^4)", "(x-y)*(x+y)*(x^2+y^2)"),
        ("factor(x^2*y+x*y^2)", "y*x*(x+y)"),
        ("factor(x*y-x)", "(y-1)*x"),
        ("factor(sin(x)^2-1)", "(sin(x)-1)*(sin(x)+1)"),
        ("factor((x^2-1)/(x^2+2*x+1))", "(x-1)/(x+1)"),
        ("factor((6*x^2-6)/(4*x+4))", "3*(x-1)/2"),
        ("factor(x^2/2-1/8)", "(2*x-1)*(2*x+1)/8"),
        ("factor((x+1)^3*(x-2)^2)", "(x-2)^2*(x+1)^3"),
        ("factor(x^6-1)", "(x-1)*(x+1)*(x^2-x+1)*(x^2+x+1)"),
        ("factor(x^4+x^2+1)", "(x^2-x+1)*(x^2+x+1)"),
        ("factor(x^4+4)", "(x^2-2*x+2)*(x^2+2*x+2)"),
        ("factor(1/(x^2-1))", "1/((x-1)*(x+1))"),
        ("factor(2/(x^2-1))", "2/((x-1)*(x+1))"),
        ("factor(x^2+1)", "x^2+1"),
        ("factor(x^2-2)", "(x-sqrt(2))*(x+sqrt(2))"),
        ("factor(x^2-3*x+1)", "(x+(-sqrt(5)-3)/2)*(x+(sqrt(5)-3)/2)"),
        ("factor(x^2-2*x-1)", "(x-sqrt(2)-1)*(x+sqrt(2)-1)"),
        ("factor(3*x^2-6)", "3*(x-sqrt(2))*(x+sqrt(2))"),
        ("factor(x^2-8)", "(x-2*sqrt(2))*(x+2*sqrt(2))"),
        ("factor(x^2-12)", "(x-2*sqrt(3))*(x+2*sqrt(3))"),
        ("factor((x^2-2)*(x+1))", "(x-sqrt(2))*(x+1)*(x+sqrt(2))"),
        ("factor((x^2-2)^2)", "(x-sqrt(2))^2*(x+sqrt(2))^2"),
        ("factor(2*x^2-4)", "2*(x-sqrt(2))*(x+sqrt(2))"),
        ("factor(x^2/2-1)", "(x-sqrt(2))*(x+sqrt(2))/2"),
        ("factor(-x^2+2)", "-(x-sqrt(2))*(x+sqrt(2))"),
        ("factor(x^2-x-1)", "(x+(-sqrt(5)-1)/2)*(x+(sqrt(5)-1)/2)"),
        (
            "factor((x^2-2)/(x^2-3))",
            "(x-sqrt(2))*(x+sqrt(2))/((x-sqrt(3))*(x+sqrt(3)))",
        ),
        ("factor(sin(x)^2-2)", "(sin(x)-sqrt(2))*(sin(x)+sqrt(2))"),
        ("factor(x^2*y-2*y)", "(x-sqrt(2))*(x+sqrt(2))*y"),
        (
            "factor(x^4-5*x^2+6)",
            "(x-sqrt(3))*(x-sqrt(2))*(x+sqrt(2))*(x+sqrt(3))",
        ),
        ("factor(y^2-2)", "(y-sqrt(2))*(y+sqrt(2))"),
        ("factor(x^4-4*x^2+1)", "x^4-4*x^2+1"),
        ("factor(x^3-2)", "x^3-2"),
        ("factor(x^4-10*x^2+1)", "x^4-10*x^2+1"),
        ("factor(exp(2*x)-1)", "exp(2*x)-1"),
        ("factor(x^(1/2)*x-x)", "x*(sqrt(x)-1)"),
        ("factor(x^2-x)", "x*(x-1)"),
        ("factor((x^2-1)/(x-1))", "x+1"),
        ("factor(-x)", "-x"),
        ("factor(x^2*y^2-1)", "(x*y-1)*(x*y+1)"),
        ("factor(x^2-y^2+x+y)", "(x+y)*(x-y+1)"),
        ("factor(2*x^2+4*x*y+2*y^2)", "2*(x+y)^2"),
        ("factor(3/4)", "3/4"),
        ("collect(x^2-2)", "x^2-2"),
        ("collect(x^2-4)", "(x-2)*(x+2)"),
        ("collect(2*x^2-2)", "2*(x-1)*(x+1)"),
        ("collect(x^6-1)", "(x-1)*(x+1)*(x^2-x+1)*(x^2+x+1)"),
        ("collect(x^4-y^4)", "(x-y)*(x+y)*(x^2+y^2)"),
        ("collect(sin(x)^2-1)", "(sin(x)-1)*(sin(x)+1)"),
        ("collect((x^2-1)/(x^2+2*x+1))", "(x-1)/(x+1)"),
        ("collect(-2*x^2+4)", "-2*(x^2-2)"),
        // Over ℚ(i) when a coefficient is complex (CAS-186).
        ("factor(x^2+(1+2*i)*x+2*i)", "(x+1)*(x+2*i)"),
        ("factor(i*x^2+i)", "i*(x+i)*(x-i)"),
        ("factor((2+i)*x+(1-2*i))", "(2+i)*(x-i)"),
        ("factor(x^2+i*x+2)", "(x+2*i)*(x-i)"),
        ("factor(2*i*x+2)", "2*i*(x-i)"),
        ("factor((x+i)^2*(x-i))", "(x+i)^2*(x-i)"),
        ("factor(x^2+2*x+1+i*(x+1))", "(x+1)*(x+1+i)"),
        ("factor(i*x^2-i)", "i*(x-1)*(x+1)"),
        ("factor(x^4+4*i)", "x^4+4*i"),
        (
            "factor(x^2+i*x+1)",
            "(x+((-i)*sqrt(5)+i)/2)*(x+(i*sqrt(5)+i)/2)",
        ),
        (
            "factor(x^3-i)",
            "(x+i)*(x+(-sqrt(3)-i)/2)*(x+(sqrt(3)-i)/2)",
        ),
        ("collect(x^2+(1+2*i)*x+2*i)", "(x+1)*(x+2*i)"),
        ("collect(i*x^2+i)", "i*(x+i)*(x-i)"),
        // factors.
        ("factors(x^2-1)", "[x-1,1,x+1,1]"),
        ("factors(12*x^4-12)", "[12,1,x-1,1,x+1,1,x^2+1,1]"),
        ("factors((x^2-1)/(x+2))", "[x-1,1,x+1,1,x+2,-1]"),
        ("factors((x^2-1)/(2*x+4))", "[x-1,1,x+1,1,2,-1,x+2,-1]"),
        ("factors(x^4-y^4)", "[x-y,1,x+y,1,x^2+y^2,1]"),
        ("factors(x^2/2-1/8)", "[2*x-1,1,2*x+1,1,8,-1]"),
        ("factors((x-1)^2/(x+1)^3)", "[x-1,2,x+1,-3]"),
        ("factors(1/(x^2-1))", "[x-1,-1,x+1,-1]"),
        ("factors(6)", "[6,1]"),
        ("factors(0)", "[0,1]"),
        ("factors(1)", "[1,1]"),
        ("factors(-1)", "[-1,1]"),
        ("factors(2/3)", "[2/3,1]"),
        ("factors(x)", "[x,1]"),
        // sqrfree, over ℤ and ℤ[i], a content's parts apart.
        ("sqrfree(x^3-x^2-x+1)", "(x+1)*(x-1)^2"),
        ("sqrfree(12*x^4-24*x^2+12)", "12*(x^2-1)^2"),
        ("sqrfree((x^2-1)/(x+1)^2)", "(x-1)/(x+1)"),
        ("sqrfree(x^2*y^2+2*x*y^2+y^2)", "y^2*(x+1)^2"),
        ("sqrfree(x^3/2-x^2-x/2+1)", "(x^3-2*x^2-x+2)/2"),
        ("sqrfree((x+1)^2/((x-1)^3*4))", "(x+1)^2/(4*(x-1)^3)"),
        ("sqrfree(i*x^2+2*i*x+i)", "i*(x+1)^2"),
        ("sqrfree(x^4-2*x^2+1)", "(x^2-1)^2"),
        ("sqrfree(6)", "6"),
    ];
    // giac's versions differ in prefixes of lists: `poly1[…]` for `[…]`.
    let normal = |s: &str| s.replace("poly1[", "[").replace("list[", "[");
    let mut wrong = Vec::new();
    for (input, answer) in recorded {
        let printed = giac.eval(input).unwrap();
        if normal(&printed) != normal(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if let Err(e) = giac.agrees(input) {
            wrong.push(e.to_string());
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // §15 detail 5: normalised, where giac returns an argument as it
        // is when the other is 0 or the same.
        ("gcd(0, -2*x+2)", "-2*x+2", "2*x-2"),
        ("gcd(-6*x, -6*x)", "-6*x", "6*x"),
        // A constant c to the other's degree; so a degree-1 discriminant is 1.
        ("resultant(2, x^2+1, x)", "1", "4"),
        ("discriminant(3*x+1, x)", "1/3", "1"),
        // §15 detail 7.
        ("degree(0, x)", "0", "-infinity"),
        // P5: the one variable, where giac assumes x.
        ("degree(y^3+y)", "0", "3"),
        ("quo(y^2+1, y+1)", "(y^2+1)/(y+1)", "y-1"),
        ("coeff(y^2+3*y, 1)", "0", "3"),
        // §15 detail 6: the degree in sin(x).
        ("degree(sin(x)^2+sin(x), sin(x))", "1", "2"),
        // §15 detail 5: lcm is expanded.
        ("lcm(x^2-1, x^2+2*x+1)", "(x-1)*(x^2+2*x+1)", "x^3+x^2-x-1"),
        ("lcm(6*x^2-6, 4*x-4)", "(3*x+3)*(4*x-4)", "12*x^2-12"),
        // Normalised, where giac keeps the signs.
        ("lcm(-4, 6)", "-12", "12"),
        ("lcm(1-x, 1+x)", "(-x+1)*(1+x)", "x^2-1"),
        // One expanded denominator a term, where giac keeps the factors
        // its tree has; evaluation makes 1/(x+1)/(x-1) 1/((x+1)*(x-1)).
        (
            "expand((1+1/(x+1))^2)",
            "(1/(x+1))^2+1+2/(x+1)",
            "1+2/(x+1)+1/(x^2+2*x+1)",
        ),
        ("expand(1/(x+1)/(x-1))", "1/((x+1)*(x-1))", "1/(x^2-1)"),
        // §15 detail 5: the main variable is the last in canonical order,
        // where giac's is the first in its order (pi before x).
        ("normal(1/(y-x))", "-1/(x-y)", "1/(y-x)"),
        ("numer((x-y)/(y-x+1))", "-x+y", "x-y"),
        ("denom((x-y)/(y-x+1))", "x-y-1", "-x+y+1"),
        ("normal(1/(pi-x))", "1/(pi-x)", "-1/(x-pi)"),
        // Over ℤ[i] the leading coefficient is in the first quadrant, where
        // giac leaves one that is neither real nor imaginary.
        (
            "normal((1+i)/(2*x+4))",
            "1/((1-i)*x+2-2*i)",
            "i/((1+i)*x+2+2*i)",
        ),
        // Quadratics split over square roots, but a non-monic one keeps
        // its leading coefficient apart, and neither quartics nor
        // quadratics in two variables split (CAS-187).
        (
            "factor(2*x^2-1)",
            "(x+sqrt(2)/2)*(2*x-sqrt(2))",
            "2*(x-sqrt(2)/2)*(x+sqrt(2)/2)",
        ),
        (
            "factor(4*x^2-8*x+1)",
            "(x+(sqrt(3)-2)/2)*(4*x-2*sqrt(3)-4)",
            "4*(x+(-sqrt(3)-2)/2)*(x+(sqrt(3)-2)/2)",
        ),
        (
            "factor(x^4+1)",
            "(x^2+(-sqrt(2))*x+1)*(x^2+sqrt(2)*x+1)",
            "x^4+1",
        ),
        ("factor(x^4-2)", "(x^2-sqrt(2))*(x^2+sqrt(2))", "x^4-2"),
        (
            "factor(x^2-2*y^2)",
            "(x+(-sqrt(2))*y)*(x+sqrt(2)*y)",
            "x^2-2*y^2",
        ),
        // The constant first, where giac lists -1 last.
        (
            "factors(-x^2+1)",
            "[x-1,1,x+1,1,-1,1]",
            "[-1,1,x-1,1,x+1,1]",
        ),
        // Over ℚ(i): factors in the order of their parts (CAS-186).
        ("factors(i*x^2+i)", "[i,1,x+i,1,x-i,1]", "[i,1,x-i,1,x+i,1]"),
        (
            "factors(x^2+(1+2*i)*x+2*i)",
            "[x+1,1,x+2*i,1]",
            "[x+2*i,1,x+1,1]",
        ),
        // Signs as factor's, where giac keeps -x+1.
        ("sqrfree(-2*x^2+4*x-2)", "-2*(-x+1)^2", "-2*(x-1)^2"),
    ];
    let normal = |s: &str| s.replace("poly1[", "[").replace("list[", "[");
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if normal(&printed) != normal(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The variables of the random inputs; none is a giac constant.
const VARIABLES: [&str; 3] = ["x", "y", "z"];

/// The terms of a polynomial in the first `n` of [`VARIABLES`]: exponent
/// vectors and coefficients.
type Terms = BTreeMap<Vec<u32>, i64>;

/// Up to `len` terms with exponents below 3 and coefficients in `-5..=5`,
/// some of them 0.
fn terms(n: usize, len: usize) -> impl Strategy<Value = Terms> {
    prop::collection::btree_map(prop::collection::vec(0u32..3, n), -5i64..=5, 0..=len)
}

/// As [`terms`], at least one and none 0, so a nonzero polynomial.
fn nonzero(n: usize, len: usize) -> impl Strategy<Value = Terms> {
    prop::collection::btree_map(
        prop::collection::vec(0u32..3, n),
        prop_oneof![-5i64..=-1, 1i64..=5],
        1..=len,
    )
}

/// The polynomial as Xcas text, in parentheses.
fn text(terms: &Terms) -> String {
    let parts: Vec<String> = terms
        .iter()
        .filter(|&(_, &c)| c != 0)
        .map(|(exps, c)| {
            let mut factors = vec![if *c < 0 {
                format!("({c})")
            } else {
                c.to_string()
            }];
            for (v, &k) in VARIABLES.iter().zip(exps) {
                match k {
                    0 => {}
                    1 => factors.push((*v).to_owned()),
                    _ => factors.push(format!("{v}^{k}")),
                }
            }
            factors.join("*")
        })
        .collect();
    if parts.is_empty() {
        "(0)".to_owned()
    } else {
        format!("({})", parts.join("+"))
    }
}

/// The terms in `x` alone: the other exponents dropped, like terms added,
/// those that cancel left out.
fn in_x(terms: Terms) -> Terms {
    let mut only = Terms::new();
    for (exps, c) in terms {
        let mut e = vec![0; exps.len()];
        e[0] = exps[0];
        *only.entry(e).or_default() += c;
    }
    only.retain(|_, c| *c != 0);
    only
}

/// `f` of giac's answer to `input`, read back, in a new context.
fn of_giac(giac: &Giac, input: &str, f: Symbol) -> Expr {
    let value = giac.value(input).unwrap_or_else(|e| panic!("{e}"));
    Context::new()
        .eval(&Expr::call(f, [value]))
        .unwrap_or_else(|e| panic!("{input}: {e}"))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `expand` of products, powers, a product over a polynomial, and a sum
    /// of two fractions, in three variables: giac's answer, term by term.
    #[test]
    fn random_expand(
        kind in 0..4,
        a in terms(3, 4),
        b in terms(3, 4),
        c in nonzero(3, 3),
        d in nonzero(3, 3),
        power in 2u32..5,
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let (a, b, c, d) = (text(&a), text(&b), text(&c), text(&d));
        let input = match kind {
            0 => format!("expand({a}*{b})"),
            1 => format!("expand({a}^{power})"),
            2 => format!("expand({a}*{b}/{c})"),
            _ => format!("expand({a}/{c}+{b}/{d})"),
        };
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }

    /// `normal` of a sum of two fractions: giac's answer in one variable;
    /// in two, the same once this CAS's `normal` gives it the sign of
    /// §15 detail 5.
    #[test]
    fn random_normal(
        univariate in any::<bool>(),
        a in terms(2, 3),
        b in terms(2, 3),
        c in nonzero(2, 3),
        d in nonzero(2, 3),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let strip = |t: Terms| if univariate { in_x(t) } else { t };
        let (c, d) = (strip(c), strip(d));
        // Coefficients that cancel to a zero denominator.
        if c.values().all(|&k| k == 0) || d.values().all(|&k| k == 0) {
            return Ok(());
        }
        let input = format!(
            "normal({}/{}+{}/{})",
            text(&strip(a)),
            text(&c),
            text(&strip(b)),
            text(&d),
        );
        if univariate {
            giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
        } else {
            let ours = ours(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
            prop_assert_eq!(of_giac(giac, &input, Symbol::TOGETHER), ours, "{}", input);
        }
    }

    /// `gcd` of two products with a common factor: giac's answer in one
    /// variable; in two, up to its sign, which giac takes from the inputs,
    /// as it does in one for two same arguments.
    #[test]
    fn random_gcd(
        univariate in any::<bool>(),
        a in nonzero(2, 3),
        b in nonzero(2, 3),
        g in nonzero(2, 3),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let strip = |t: Terms| if univariate { in_x(t) } else { t };
        let (a, b, g) = (strip(a), strip(b), strip(g));
        // giac returns a nonzero input as it is when the other is 0.
        if [&a, &b, &g].iter().any(|t| t.values().all(|&k| k == 0)) {
            return Ok(());
        }
        let same = a == b;
        let (a, b, g) = (text(&a), text(&b), text(&g));
        let input = format!("gcd({a}*{g},{b}*{g})");
        if univariate && !same {
            giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
        } else {
            let ours = ours(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
            let theirs = of_giac(giac, &input, Symbol::EXPAND);
            let negated = Context::new()
                .eval(&Expr::call(Symbol::EXPAND, [-theirs.clone()]))
                .unwrap();
            prop_assert!(ours == theirs || ours == negated, "{}: giac {}, ours {}", input, theirs, ours);
        }
    }
}
