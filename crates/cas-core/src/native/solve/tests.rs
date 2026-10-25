//! Tests of `Solve` and `CSolve` (CAS-210): giac 1.9's solutions of
//! polynomial and rational equations in one unknown, with parameters, in
//! both domains; the forms of an equation, poles, identities, and what
//! stays.

use crate::assume::Domain;
use crate::trace::EvalOptions;
use crate::{Context, Expr, Warning};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated in `ctx`: the same
/// expression.
fn check_in(ctx: &Context, cases: &[(&str, &str)]) {
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

fn check(cases: &[(&str, &str)]) {
    check_in(&Context::new(), cases);
}

/// The warnings of evaluating `input`, as messages.
fn warnings(input: &str) -> Vec<String> {
    let (_, report) = Context::new().eval_with(&read(input), &EvalOptions::default());
    report.warnings.iter().map(Warning::to_string).collect()
}

/// `(1 ± i√3)/2` or `(−1 ± i√3)/2`, as `c + d i √3` in full form.
fn sixth_root(c: &str, d: &str) -> String {
    format!("Plus[{c}, Times[{d}, I, Power[3, 1/2]]]")
}

/// giac's real solutions: ascending, each once.
#[test]
fn polynomial_equations() {
    check(&[
        // x^2-3x+2, x^2+1, x^4-5x^2+6, x^3-6x^2+11x-6.
        ("Solve[Plus[2, Times[-3, x], Power[x, 2]], x]", "List[1, 2]"),
        ("Solve[Equation[Plus[1, Power[x, 2]], 0], x]", "List[]"),
        (
            "Solve[Plus[6, Times[-5, Power[x, 2]], Power[x, 4]], x]",
            "List[Times[-1, Power[3, 1/2]], Times[-1, Power[2, 1/2]], Power[2, 1/2], Power[3, 1/2]]",
        ),
        (
            "Solve[Plus[-6, Times[11, x], Times[-6, Power[x, 2]], Power[x, 3]], x]",
            "List[1, 2, 3]",
        ),
        // x^2=2, 2x+3=x-1, repeated roots once, x^6-1, x^2+x+1.
        (
            "Solve[Equation[Power[x, 2], 2], x]",
            "List[Times[-1, Power[2, 1/2]], Power[2, 1/2]]",
        ),
        (
            "Solve[Equation[Plus[3, Times[2, x]], Plus[x, -1]], x]",
            "List[-4]",
        ),
        ("Solve[Plus[1, Times[-2, x], Power[x, 2]], x]", "List[1]"),
        ("Solve[Plus[-1, Power[x, 6]], x]", "List[-1, 1]"),
        ("Solve[Plus[1, x, Power[x, 2]], x]", "List[]"),
        // x^2-x-1: (1 - √5)/2 first.
        (
            "Solve[Plus[-1, Times[-1, x], Power[x, 2]], x]",
            "List[Plus[1/2, Times[-1/2, Power[5, 1/2]]], Plus[1/2, Times[1/2, Power[5, 1/2]]]]",
        ),
        // Equal, as Equation.
        ("Solve[Equal[Power[x, 2], 4], x]", "List[-2, 2]"),
    ]);
}

/// Rational equations: the denominator's zeros are no solutions.
#[test]
fn rational_equations() {
    check(&[
        // 1/x=2, (x-1)/(x^2-1), (x+1)/(x-2)=3, x/(x^2+1), 1/(x^2+1).
        ("Solve[Equation[Power[x, -1], 2], x]", "List[1/2]"),
        (
            "Solve[Times[Plus[-1, x], Power[Plus[-1, Power[x, 2]], -1]], x]",
            "List[]",
        ),
        (
            "Solve[Equation[Times[Plus[1, x], Power[Plus[-2, x], -1]], 3], x]",
            "List[7/2]",
        ),
        (
            "Solve[Times[x, Power[Plus[1, Power[x, 2]], -1]], x]",
            "List[0]",
        ),
        ("Solve[Power[Plus[1, Power[x, 2]], -1], x]", "List[]"),
    ]);
}

/// giac's complex solutions: the real ones first, ascending, then by real
/// part, the positive imaginary part first.
#[test]
fn complex_solutions() {
    check(&[
        (
            "CSolve[Equation[Plus[1, Power[x, 2]], 0], x]",
            "List[I, Times[-1, I]]",
        ),
        (
            "CSolve[Plus[1, x, Power[x, 2]], x]",
            &format!(
                "List[{}, {}]",
                sixth_root("-1/2", "1/2"),
                sixth_root("-1/2", "-1/2")
            ),
        ),
        (
            "CSolve[Plus[-1, Power[x, 6]], x]",
            &format!(
                "List[-1, 1, {}, {}, {}, {}]",
                sixth_root("-1/2", "1/2"),
                sixth_root("-1/2", "-1/2"),
                sixth_root("1/2", "1/2"),
                sixth_root("1/2", "-1/2")
            ),
        ),
        // x^4+5x^2+4: 2i before i.
        (
            "CSolve[Plus[4, Times[5, Power[x, 2]], Power[x, 4]], x]",
            "List[Times[2, I], Times[-2, I], I, Times[-1, I]]",
        ),
        // (x-2)(x^2+1)(x+3).
        (
            "CSolve[Times[Plus[-2, x], Plus[1, Power[x, 2]], Plus[3, x]], x]",
            "List[-3, 2, I, Times[-1, I]]",
        ),
        // A complex coefficient: every root, in Solve too, as giac.
        ("Solve[Plus[x, Times[-1, I]], x]", "List[I]"),
        (
            "Solve[Plus[-1, Times[2, I, x], Power[x, 2]], x]",
            "List[Times[-1, I]]",
        ),
        // x^2+i: √(−4i) = √2 − √2 i, exactly.
        (
            "CSolve[Plus[I, Power[x, 2]], x]",
            "List[Plus[Times[-1/2, Power[2, 1/2]], Times[1/2, I, Power[2, 1/2]]], \
             Plus[Times[1/2, Power[2, 1/2]], Times[-1/2, I, Power[2, 1/2]]]]",
        ),
    ]);
}

/// In the Complex domain, `Solve` is `CSolve`.
#[test]
fn the_complex_domain() {
    let mut ctx = Context::new();
    ctx.set_domain(Domain::Complex);
    check_in(
        &ctx,
        &[("Solve[Plus[1, Power[x, 2]], x]", "List[I, Times[-1, I]]")],
    );
}

/// Parameters: giac's formulas, `+√Δ` first, the sign of `Δ` from the
/// assumptions where they tell it.
#[test]
fn parameters() {
    check(&[
        // a*x+b, a*x+a, x^2-2ax+a^2.
        (
            "Solve[Plus[b, Times[a, x]], x]",
            "List[Times[-1, b, Power[a, -1]]]",
        ),
        ("Solve[Plus[a, Times[a, x]], x]", "List[-1]"),
        (
            "Solve[Plus[Power[a, 2], Times[-2, a, x], Power[x, 2]], x]",
            "List[a]",
        ),
        // x^2+y^2-1: the square 4 out of Δ = 4 - 4y^2.
        (
            "Solve[Plus[-1, Power[x, 2], Power[y, 2]], x]",
            "List[Power[Plus[1, Times[-1, Power[y, 2]]], 1/2], \
             Times[-1, Power[Plus[1, Times[-1, Power[y, 2]]], 1/2]]]",
        ),
        // x^2-a.
        (
            "Solve[Plus[Times[-1, a], Power[x, 2]], x]",
            "List[Power[a, 1/2], Times[-1, Power[a, 1/2]]]",
        ),
        // x^2+a^2: Δ = -4a^2 is not positive, so no real solution.
        ("Solve[Plus[Power[a, 2], Power[x, 2]], x]", "List[]"),
        // x^2-pi and x^2+pi: Δ's sign from its ball.
        (
            "Solve[Plus[Times[-1, Pi], Power[x, 2]], x]",
            "List[Times[-1, Power[Pi, 1/2]], Power[Pi, 1/2]]",
        ),
        ("Solve[Plus[Pi, Power[x, 2]], x]", "List[]"),
        // x^2-2√2x+2 = (x-√2)^2: one root.
        (
            "Solve[Plus[2, Times[-2, Power[2, 1/2], x], Power[x, 2]], x]",
            "List[Power[2, 1/2]]",
        ),
        // x^2*y-y: the factor y is no root.
        (
            "Solve[Plus[Times[-1, y], Times[Power[x, 2], y]], x]",
            "List[-1, 1]",
        ),
    ]);
    let ctx = Context::new().assuming(read("Greater[a, 0]"));
    check_in(&ctx, &[("Solve[Plus[a, Power[x, 2]], x]", "List[]")]);
}

/// An identity is solved by every `x`, a nonzero constant by none, as
/// giac's `solve(x=x,x)` is `[x]` and `solve(1,x)` is `[]`.
#[test]
fn identities_and_contradictions() {
    check(&[
        ("Solve[0, x]", "List[x]"),
        ("Solve[Equation[x, x], x]", "List[x]"),
        ("Solve[Equal[x, x], x]", "List[x]"),
        ("Solve[1, x]", "List[]"),
        ("Solve[a, x]", "List[]"),
        (
            "Solve[Equation[Plus[1, Power[x, 2]], Power[x, 2]], x]",
            "List[]",
        ),
        (
            "Solve[Equation[Plus[Times[-1, x], Power[x, 2]], Power[x, 2]], x]",
            "List[0]",
        ),
        ("CSolve[0, x]", "List[x]"),
    ]);
}

/// What these methods do not solve stays.
#[test]
fn what_stays() {
    check(&[
        (
            "Solve[Plus[-1, Times[-1, x], Power[x, 5]], x]",
            "Solve[Plus[-1, Times[-1, x], Power[x, 5]], x]",
        ),
        (
            "Solve[Plus[-1/2, Sin[x]], x]",
            "Solve[Plus[-1/2, Sin[x]], x]",
        ),
        ("Solve[Power[x, 2], 2]", "Solve[Power[x, 2], 2]"),
        ("Solve[Power[x, 2]]", "Solve[Power[x, 2]]"),
        (
            "Solve[List[Equation[Plus[x, y], 1]], List[x, y]]",
            "Solve[List[Equation[Plus[x, y], 1]], List[x, y]]",
        ),
    ]);
}

/// giac's warnings: an expression is solved as `e = 0`, and a factor of a
/// degree past 2 is not solved.
#[test]
fn warnings_say_why() {
    let expression = warnings("Solve[Plus[-1, Power[x, 2]], x]");
    assert!(
        expression
            .iter()
            .any(|m| m.contains("argument is not an equation")),
        "{expression:?}"
    );
    assert!(warnings("Solve[Equation[Power[x, 2], 1], x]").is_empty());
    let quintic = warnings("Solve[Equation[Power[x, 5], Plus[1, x]], x]");
    assert!(
        quintic.iter().any(|m| m.contains("a factor of degree 5")),
        "{quintic:?}"
    );
}
