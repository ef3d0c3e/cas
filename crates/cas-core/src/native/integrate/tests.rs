//! Tests of `Integrate` (CAS-194): giac's antiderivatives, compared as
//! values once evaluated; what stays; the trace and its verification; and a
//! proptest that the derivative of a polynomial's antiderivative is the
//! polynomial.

use proptest::prelude::*;

use crate::trace::{EvalOptions, Step, StepKind, Trace, TraceLevel};
use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a` and `b`, evaluated, are the same expression, or rational
/// functions of their generators whose difference is 0.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let difference = Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]);
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// `(input, expected)`, both evaluated: the same value.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx
            .eval(&read(input))
            .unwrap_or_else(|e| panic!("{input}: {e}"));
        let expected = ctx.eval(&read(expected)).unwrap();
        if !same(&ctx, &found, &expected) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answers (CAS-DOC-9 §2), as values.
#[test]
fn antiderivatives() {
    check(&[
        // integrate(x^2,x) = x^3/3.
        ("Integrate[Power[x, 2], x]", "Times[1/3, Power[x, 3]]"),
        // 3*x^3/3+2*x^2/2+x.
        (
            "Integrate[Plus[Times[3, Power[x, 2]], Times[2, x], 1], x]",
            "Plus[x, Power[x, 2], Power[x, 3]]",
        ),
        // ln(abs(x)) and -1/x.
        ("Integrate[Power[x, -1], x]", "Log[Abs[x]]"),
        ("Integrate[Power[x, -2], x]", "Times[-1, Power[x, -1]]"),
        // x^(n+1)/(n+1), even for a symbolic n.
        (
            "Integrate[Power[x, n], x]",
            "Times[Power[x, Plus[n, 1]], Power[Plus[n, 1], -1]]",
        ),
        // 3*x^(1/3)*x/4.
        ("Integrate[Power[x, 1/3], x]", "Times[3/4, Power[x, 4/3]]"),
        // a*x^2/2, and constants.
        ("Integrate[Times[a, x], x]", "Times[1/2, a, Power[x, 2]]"),
        ("Integrate[c, x]", "Times[c, x]"),
        ("Integrate[3, x]", "Times[3, x]"),
        ("Integrate[x, x]", "Times[1/2, Power[x, 2]]"),
        (
            "Integrate[Times[x, Power[Plus[y, 1], -1]], x]",
            "Times[1/2, Power[x, 2], Power[Plus[1, y], -1]]",
        ),
        // x^2/2+integrate(f(x),x), and a product that stays whole.
        (
            "Integrate[Plus[x, f[x]], x]",
            "Plus[Times[1/2, Power[x, 2]], Integrate[f[x], x]]",
        ),
        (
            "Integrate[Plus[x, Times[2, f[x]]], x]",
            "Plus[Times[1/2, Power[x, 2]], Integrate[Times[2, f[x]], x]]",
        ),
        // [1/2*x^2,1/3*x^3].
        (
            "Integrate[List[x, Power[x, 2]], x]",
            "List[Times[1/2, Power[x, 2]], Times[1/3, Power[x, 3]]]",
        ),
        // Polynomials that are not expanded; giac's (x+1)^3/3 is 1/3 more.
        (
            "Integrate[Power[Plus[x, 1], 2], x]",
            "Plus[x, Power[x, 2], Times[1/3, Power[x, 3]]]",
        ),
        (
            "Integrate[Times[x, Plus[x, y]], x]",
            "Plus[Times[1/3, Power[x, 3]], Times[1/2, Power[x, 2], y]]",
        ),
        (
            "Integrate[Times[Plus[x, 1], Plus[x, -1], Power[Plus[y, 2], -1]], x]",
            "Times[Plus[Times[1/3, Power[x, 3]], Times[-1, x]], Power[Plus[y, 2], -1]]",
        ),
        // Gaussian coefficients.
        (
            "Integrate[Times[Complex[0, 1], Power[x, 3]], x]",
            "Times[Complex[0, 1/4], Power[x, 4]]",
        ),
        // Definite integrals of polynomials: integrate(x^2,x,0,1) = 1/3.
        ("Integrate[Power[x, 2], x, 0, 1]", "1/3"),
        ("Integrate[Power[x, 2], List[x, 0, 1]]", "1/3"),
        (
            "Integrate[Plus[x, 1], x, a, b]",
            "Plus[Times[1/2, Power[b, 2]], b, Times[-1/2, Power[a, 2]], Times[-1, a]]",
        ),
        ("Integrate[Power[Plus[x, 1], 2], x, -1, 2]", "9"),
    ]);
}

/// What the methods of CAS-194 leave as it is.
#[test]
fn integrals_that_stay() {
    let ctx = Context::new();
    for input in [
        // The elementary functions come with CAS-195, rational functions with
        // CAS-196.
        "Integrate[Sin[x], x]",
        "Integrate[Exp[x], x]",
        "Integrate[Power[Plus[x, 1], -1], x]",
        "Integrate[f[x], x]",
        // giac's integrate(2*f(x),x) stays whole.
        "Integrate[Times[2, f[x]], x]",
        "Integrate[Power[x, x], x]",
        // A definite integral whose antiderivative is not a polynomial
        // waits for CAS-198.
        "Integrate[Power[x, -1], x, 1, 2]",
        // Not a symbol, or the wrong arguments.
        "Integrate[x, 2]",
        "Integrate[x, x, 1]",
        "Integrate[x]",
        "Integrate[x, List[x, 1]]",
    ] {
        // The arguments evaluate, the call stays.
        let e = read(input);
        let found = ctx.eval(&e).unwrap();
        assert_eq!(found.head(), e.head(), "{input}: {found:?}");
        assert_eq!(found.args().len(), e.args().len(), "{input}: {found:?}");
    }
}

/// One `"Integrate"` step a call, which `Trace::verify` checks through the
/// derivative of the antiderivative, as it checks the threading of lists.
#[test]
fn traces_verify() {
    let ctx = Context::new();
    for input in [
        "Integrate[Plus[Times[3, Power[x, 2]], Times[2, x], 1], x]",
        "Integrate[Power[x, -1], x]",
        "Integrate[Power[x, 1/2], x]",
        "Integrate[Times[a, Power[Plus[x, y], 2]], x]",
        "Integrate[Plus[x, f[x]], x]",
        "Integrate[List[x, Power[x, 2]], x]",
    ] {
        let options = EvalOptions {
            trace: TraceLevel::Rules,
        };
        let (result, report) = ctx.eval_with(&read(input), &options);
        result.unwrap();
        let verification = report.trace.verify(&ctx);
        assert!(
            verification.wrong.is_none(),
            "{input}: {:?}",
            verification.wrong
        );
        assert!(verification.passed > 0, "{input}");
        let steps: Vec<&Step> = report
            .trace
            .steps
            .iter()
            .filter(|s| matches!(s.kind, StepKind::Algorithm(_)))
            .collect();
        assert!(!steps.is_empty(), "{input}");
        assert_eq!(
            steps[0].before.head(),
            Some(&Expr::from(Symbol::INTEGRATE)),
            "{input}"
        );
    }
}

/// A wrong antiderivative fails the verification: the step's claim is its
/// derivative.
#[test]
fn a_wrong_antiderivative_is_caught() {
    let ctx = Context::new();
    let step = |after: &str| Step {
        id: 0,
        kind: StepKind::Algorithm("Integrate".into()),
        position: Vec::new(),
        before: read("Integrate[Power[x, 2], x]"),
        after: read(after),
        parent: None,
        children: Vec::new(),
    };
    let wrong = Trace {
        steps: vec![step("Times[1/2, Power[x, 3]]")],
    };
    assert!(wrong.verify(&ctx).wrong.is_some());
    // A constant more is right.
    let right = Trace {
        steps: vec![step("Plus[7, Times[1/3, Power[x, 3]]]")],
    };
    let verification = right.verify(&ctx);
    assert!(verification.wrong.is_none() && verification.passed == 1);
}

/// Terms in `x` and `y`, small exponents and coefficients.
fn terms() -> impl Strategy<Value = Vec<((u32, u32), i64)>> {
    prop::collection::vec(((0u32..5, 0u32..3), -9i64..9), 1..6)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// `D[Integrate[p, x], x]` is `p` for polynomials `p` in `x` and `y`.
    #[test]
    fn derivatives_of_antiderivatives(terms in terms()) {
        let ctx = Context::new();
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let p: Expr = terms
            .iter()
            .map(|&((i, j), c)| Expr::from(c) * x.pow(i64::from(i)) * y.pow(i64::from(j)))
            .sum();
        let p = ctx.eval(&p).unwrap();
        let integral = ctx.eval(&Expr::call(Symbol::INTEGRATE, [p.clone(), x.clone()])).unwrap();
        prop_assert!(!integral.contains(&Expr::from(Symbol::INTEGRATE)), "{integral:?}");
        let derivative = ctx.eval(&Expr::call(Symbol::D, [integral, x])).unwrap();
        prop_assert!(same(&ctx, &derivative, &p), "{derivative:?} against {p:?}");
    }
}
