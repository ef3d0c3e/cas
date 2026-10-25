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
        // x^3+a*x+1: parameters past degree 2, but for a binomial.
        (
            "Solve[Plus[1, Times[a, x], Power[x, 3]], x]",
            "Solve[Plus[1, Times[a, x], Power[x, 3]], x]",
        ),
        // Past the degree of numeric roots.
        (
            "Solve[Plus[-1, Times[-1, x], Power[x, 65]], x]",
            "Solve[Plus[-1, Times[-1, x], Power[x, 65]], x]",
        ),
        // x e^x = 1, giac's LambertW; sin x = cos x.
        (
            "Solve[Equation[Times[x, Exp[x]], 1], x]",
            "Solve[Equation[Times[x, Exp[x]], 1], x]",
        ),
        (
            "Solve[Equation[Sin[x], Cos[x]], x]",
            "Solve[Equation[Sin[x], Cos[x]], x]",
        ),
        ("Solve[Power[x, 2], 2]", "Solve[Power[x, 2], 2]"),
        ("Solve[Power[x, 2]]", "Solve[Power[x, 2]]"),
        // A nonlinear system (CAS-215), repeated unknowns.
        (
            "Solve[List[Equation[Times[x, y], 1], Equation[Plus[x, y], 2]], List[x, y]]",
            "Solve[List[Equation[Times[x, y], 1], Equation[Plus[x, y], 2]], List[x, y]]",
        ),
        (
            "Solve[List[Equation[Plus[x, y], 1]], List[x, x]]",
            "Solve[List[Equation[Plus[x, y], 1]], List[x, x]]",
        ),
    ]);
}

/// giac's warnings: an expression is solved as `e = 0`; and why a factor
/// of a degree past 2 is not solved.
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
    let cubic = warnings("Solve[Equation[Power[x, 3], Plus[1, Times[a, x]]], x]");
    assert!(
        cubic
            .iter()
            .any(|m| m.contains("a factor of degree 3 with parameters")),
        "{cubic:?}"
    );
    let high = warnings("Solve[Equation[Power[x, 65], Plus[1, x]], x]");
    assert!(
        high.iter()
            .any(|m| m.contains("not sought numerically past degree 64")),
        "{high:?}"
    );
}

/// The values of the solutions of `input`, approximated: real and
/// imaginary parts.
fn values(input: &str) -> Vec<(f64, f64)> {
    let ctx = Context::new();
    let found = ctx
        .eval(&Expr::call(crate::Symbol::N, [read(input)]))
        .unwrap();
    assert_eq!(
        found.head(),
        Some(&Expr::from(crate::Symbol::LIST)),
        "{input}: {found:?}"
    );
    found
        .args()
        .iter()
        .map(|v| match v.view() {
            crate::View::Ball(b) => (b.mid().to_f64(), 0.0),
            crate::View::ComplexBall(z) => (z.re().mid().to_f64(), z.im().mid().to_f64()),
            _ => panic!("{input}: {v:?} is no number"),
        })
        .collect()
}

/// `(input, values)`: the solutions, approximated, are the values, in
/// order, within `10⁻¹²`.
fn check_values(cases: &[(&str, &[(f64, f64)])]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = values(input);
        let close = found.len() == expected.len()
            && found
                .iter()
                .zip(expected.iter())
                .all(|(f, e)| (f.0 - e.0).abs() < 1e-12 && (f.1 - e.1).abs() < 1e-12);
        if !close {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Binomials and quadratics in a power, in radicals (CAS-211).
#[test]
fn binomials_and_powers() {
    check(&[
        // x^3-2, x^4-2, x^3-a, x^4-a: giac's radicals.
        ("Solve[Plus[-2, Power[x, 3]], x]", "List[Power[2, 1/3]]"),
        (
            "CSolve[Plus[-2, Power[x, 3]], x]",
            "List[Power[2, 1/3], \
             Times[Power[2, 1/3], Plus[-1/2, Times[1/2, I, Power[3, 1/2]]]], \
             Times[Power[2, 1/3], Plus[-1/2, Times[-1/2, I, Power[3, 1/2]]]]]",
        ),
        (
            "Solve[Plus[-2, Power[x, 4]], x]",
            "List[Times[-1, Power[2, 1/4]], Power[2, 1/4]]",
        ),
        ("Solve[Plus[2, Power[x, 4]], x]", "List[]"),
        (
            "Solve[Plus[2, Power[x, 3]], x]",
            "List[Times[-1, Power[2, 1/3]]]",
        ),
        (
            "Solve[Plus[Times[-1, a], Power[x, 3]], x]",
            "List[Power[a, 1/3]]",
        ),
        (
            "Solve[Plus[Times[-1, a], Power[x, 4]], x]",
            "List[Power[a, 1/4], Times[-1, Power[a, 1/4]]]",
        ),
        // x^4+1: (±1 ± i)/√2.
        (
            "CSolve[Plus[1, Power[x, 4]], x]",
            "List[Times[Plus[-1, I], Power[2, -1/2]], Times[Plus[-1, Times[-1, I]], Power[2, -1/2]], \
             Times[Plus[1, I], Power[2, -1/2]], Times[Plus[1, Times[-1, I]], Power[2, -1/2]]]",
        ),
        // x^4-10x^2+1: ±√2 ± √3, denested, as giac.
        (
            "Solve[Plus[1, Times[-10, Power[x, 2]], Power[x, 4]], x]",
            "List[Times[-1, Plus[Power[2, 1/2], Power[3, 1/2]]], \
             Times[-1, Plus[Times[-1, Power[2, 1/2]], Power[3, 1/2]]], \
             Plus[Times[-1, Power[2, 1/2]], Power[3, 1/2]], Plus[Power[2, 1/2], Power[3, 1/2]]]",
        ),
    ]);
    // x^6-3x^3+1: the cube roots of (3 ± √5)/2.
    check_values(&[(
        "Solve[Plus[1, Times[-3, Power[x, 3]], Power[x, 6]], x]",
        &[
            (0.725_562_630_246_326_3, 0.0),
            (1.378_240_772_489_210_2, 0.0),
        ],
    )]);
}

/// Cubics and quartics in radicals (CAS-211): Cardano's formula, the
/// trigonometric form, Ferrari's method.
#[test]
fn cubics_and_quartics() {
    let (re, im) = (0.341_163_901_914_009_66, 1.161_541_399_997_252);
    check_values(&[
        (
            "Solve[Plus[1, x, Power[x, 3]], x]",
            &[(-0.682_327_803_828_019_3, 0.0)],
        ),
        (
            "CSolve[Plus[1, x, Power[x, 3]], x]",
            &[(-0.682_327_803_828_019_3, 0.0), (re, im), (re, -im)],
        ),
        (
            "Solve[Plus[1, Times[-3, x], Power[x, 3]], x]",
            &[
                (-1.879_385_241_571_816_8, 0.0),
                (0.347_296_355_333_860_7, 0.0),
                (1.532_088_886_237_956, 0.0),
            ],
        ),
        // 2x^3-3x^2-1, shifted.
        (
            "Solve[Plus[-1, Times[-3, Power[x, 2]], Times[2, Power[x, 3]]], x]",
            &[(1.677_650_698_804_06, 0.0)],
        ),
        // x^4+4x-1: the resolvent m^3+m-2 has the root 1.
        (
            "Solve[Plus[-1, Times[4, x], Power[x, 4]], x]",
            &[
                (-1.663_251_938_771_469_4, 0.0),
                (0.249_038_376_398_374_33, 0.0),
            ],
        ),
        (
            "CSolve[Plus[-1, Times[4, x], Power[x, 4]], x]",
            &[
                (-1.663_251_938_771_469_4, 0.0),
                (0.249_038_376_398_374_33, 0.0),
                (0.707_106_781_186_547_5, 1.383_551_069_665_697_3),
                (0.707_106_781_186_547_5, -1.383_551_069_665_697_3),
            ],
        ),
    ]);
    // Cardano's real root is in radicals, not a ball.
    let found = Context::new()
        .eval(&read("Solve[Plus[1, x, Power[x, 3]], x]"))
        .unwrap();
    assert!(!found.args()[0].is_number(), "{found:?}");
}

/// Numeric roots of higher degrees (CAS-211), balls at the working
/// precision.
#[test]
fn numeric_roots() {
    let (a, b) = (-0.764_884_433_600_584_7, 0.352_471_546_031_726_25);
    let (c, d) = (0.181_232_444_469_875_38, 1.083_954_101_317_710_7);
    check_values(&[
        (
            "Solve[Plus[-1, Times[-1, x], Power[x, 5]], x]",
            &[(1.167_303_978_261_418_7, 0.0)],
        ),
        (
            "CSolve[Plus[-1, Times[-1, x], Power[x, 5]], x]",
            &[
                (1.167_303_978_261_418_7, 0.0),
                (a, b),
                (a, -b),
                (c, d),
                (c, -d),
            ],
        ),
        (
            "Solve[Plus[1, Times[-5, x], Power[x, 5]], x]",
            &[
                (-1.541_651_684_104_525, 0.0),
                (0.200_064_102_629_975_4, 0.0),
                (1.440_500_397_341_56, 0.0),
            ],
        ),
    ]);
    let found = Context::new()
        .eval(&read("Solve[Plus[-1, Times[-1, x], Power[x, 5]], x]"))
        .unwrap();
    assert!(
        matches!(found.args()[0].view(), crate::View::Ball(_)),
        "{found:?}"
    );
}

/// Linear systems (CAS-212), giac's tuples: one solution, a family in its
/// free unknowns, none; parameters; lists of one unknown.
#[test]
fn linear_systems() {
    let two = "List[Equation[Plus[x, y], 1], Equation[Plus[x, Times[-1, y]], 3]]";
    check(&[
        (&format!("Solve[{two}, List[x, y]]"), "List[List[2, -1]]"),
        // A free unknown no equation holds.
        (
            &format!("Solve[{two}, List[x, y, z]]"),
            "List[List[2, -1, z]]",
        ),
        (
            "Solve[List[Equation[Plus[x, y, z], 6], Equation[Plus[Times[2, x], Times[-1, y], z], 3], \
             Equation[Plus[x, Times[2, y], Times[-1, z]], 2]], List[x, y, z]]",
            "List[List[1, 2, 3]]",
        ),
        // A family: x + y = 1 twice over.
        (
            "Solve[List[Equation[Plus[x, y], 1], Equation[Plus[Times[2, x], Times[2, y]], 2]], List[x, y]]",
            "List[List[Plus[1, Times[-1, y]], y]]",
        ),
        (
            "Solve[Equation[Plus[x, y], 1], List[x, y]]",
            "List[List[Plus[1, Times[-1, y]], y]]",
        ),
        (
            "Solve[List[Equation[Plus[x, y, z], 1], Equation[Plus[x, Times[-1, y]], 0]], List[x, y, z]]",
            "List[List[Plus[1/2, Times[-1/2, z]], Plus[1/2, Times[-1/2, z]], z]]",
        ),
        // None.
        (
            "Solve[List[Equation[Plus[x, y], 1], Equation[Plus[x, y], 2]], List[x, y]]",
            "List[]",
        ),
        // Every point: 0 = 0.
        (
            "Solve[List[Equation[0, 0]], List[x, y]]",
            "List[List[x, y]]",
        ),
        // Parameters: ax + y = 1, x - y = b.
        (
            "Solve[List[Equation[Plus[Times[a, x], y], 1], Equation[Plus[x, Times[-1, y]], b]], List[x, y]]",
            "List[List[Times[Plus[1, b], Power[Plus[1, a], -1]], \
             Times[Plus[1, Times[-1, a, b]], Power[Plus[1, a], -1]]]]",
        ),
        // A complex coefficient: x + iy = 1, x = y.
        (
            "CSolve[List[Equation[Plus[x, Times[I, y]], 1], Equation[x, y]], List[x, y]]",
            "List[List[Plus[1/2, Times[-1/2, I]], Plus[1/2, Times[-1/2, I]]]]",
        ),
        // One unknown: a list of it, or of the equations.
        (
            "Solve[Plus[-1, Power[x, 2]], List[x]]",
            "List[List[-1], List[1]]",
        ),
        ("Solve[List[Equation[Plus[x, 1], 0]], x]", "List[-1]"),
        (
            "Solve[List[Equation[x, 1], Equation[Times[2, x], 2]], x]",
            "List[1]",
        ),
        ("Solve[List[Equation[x, 1], Equation[x, 2]], x]", "List[]"),
    ]);
}

/// Equations in functions of the unknown (CAS-213): substitution and
/// inversion to giac's principal solutions, square roots squared away and
/// the candidates checked.
#[test]
fn transcendental_equations() {
    check(&[
        ("Solve[Equation[Exp[x], 2], x]", "List[Log[2]]"),
        ("Solve[Equation[Log[x], 3], x]", "List[Exp[3]]"),
        ("Solve[Equation[Log[x], -1], x]", "List[Exp[-1]]"),
        ("Solve[Equation[Sin[x], 0], x]", "List[0, Pi]"),
        (
            "Solve[Equation[Sin[x], 1/2], x]",
            "List[Times[1/6, Pi], Times[5/6, Pi]]",
        ),
        (
            "Solve[Equation[Cos[x], 1/2], x]",
            "List[Times[-1/3, Pi], Times[1/3, Pi]]",
        ),
        ("Solve[Equation[Tan[x], 1], x]", "List[Times[1/4, Pi]]"),
        (
            "Solve[Equation[Cos[x], 0], x]",
            "List[Times[-1/2, Pi], Times[1/2, Pi]]",
        ),
        (
            "Solve[Equation[Sin[Times[2, x]], 1], x]",
            "List[Times[1/4, Pi]]",
        ),
        ("Solve[Equation[ArcTan[x], 1], x]", "List[Tan[1]]"),
        ("Solve[Equation[ArcSin[x], 1/2], x]", "List[Sin[1/2]]"),
        // No real solution; a complex one.
        ("Solve[Equation[Exp[x], -1], x]", "List[]"),
        ("CSolve[Equation[Exp[x], -1], x]", "List[Times[I, Pi]]"),
        ("Solve[Equation[Sin[x], 2], x]", "List[]"),
        ("Solve[Equation[Power[x, 1/2], -1], x]", "List[]"),
        // Roots: sqrt(x)=3, 1/sqrt(x)=2, x^(1/3)=2, sqrt(x^2+1)=2.
        ("Solve[Equation[Power[x, 1/2], 3], x]", "List[9]"),
        ("Solve[Equation[Power[x, -1/2], 2], x]", "List[1/4]"),
        ("Solve[Equation[Power[x, 1/3], 2], x]", "List[8]"),
        (
            "Solve[Equation[Power[Plus[1, Power[x, 2]], 1/2], 2], x]",
            "List[Times[-1, Power[3, 1/2]], Power[3, 1/2]]",
        ),
        // Squared away, the candidates checked: sqrt(x+1)=x-1 is 3, not 0.
        (
            "Solve[Equation[Power[Plus[x, 1], 1/2], Plus[x, -1]], x]",
            "List[3]",
        ),
        (
            "Solve[Equation[Plus[Power[Plus[x, -1], 1/2], 1], x], x]",
            "List[1, 2]",
        ),
        (
            "Solve[Equation[Plus[Power[x, 1/2], Power[Plus[x, 5], 1/2]], 5], x]",
            "List[4]",
        ),
        // Polynomial in one function: exp(2x)-3exp(x)+2, exp(x)+exp(-x)=3,
        // sin(x)^2=1/4.
        (
            "Solve[Plus[2, Times[-3, Exp[x]], Exp[Times[2, x]]], x]",
            "List[0, Log[2]]",
        ),
        (
            "Solve[Equation[Plus[Exp[x], Exp[Times[-1, x]]], 3], x]",
            "List[Log[Plus[3/2, Times[-1/2, Power[5, 1/2]]]], \
             Log[Plus[3/2, Times[1/2, Power[5, 1/2]]]]]",
        ),
        (
            "Solve[Equation[Power[Sin[x], 2], 1/4], x]",
            "List[Times[-1/6, Pi], Times[1/6, Pi], Times[5/6, Pi], Times[7/6, Pi]]",
        ),
        // 2^x = 8, ln(x^2-3) = 0, exp(x^2) = e.
        ("Solve[Equation[Power[2, x], 8], x]", "List[3]"),
        (
            "Solve[Equation[Log[Plus[-3, Power[x, 2]]], 0], x]",
            "List[-2, 2]",
        ),
        ("Solve[Equation[Exp[Power[x, 2]], E], x]", "List[-1, 1]"),
    ]);
}

/// The value of the number `input` evaluates to: real and imaginary
/// parts.
fn bare_value(input: &str) -> (f64, f64) {
    let found = Context::new().eval(&read(input)).unwrap();
    match found.view() {
        crate::View::Ball(b) => (b.mid().to_f64(), 0.0),
        crate::View::ComplexBall(z) => (z.re().mid().to_f64(), z.im().mid().to_f64()),
        _ => panic!("{input}: {found:?} is no number"),
    }
}

/// Numeric roots (CAS-214): a polynomial's real roots, the number alone
/// when there is one; sign changes bisected; Newton's iteration.
#[test]
fn numeric_solutions() {
    let close = |(a, b): (f64, f64), x: f64| (a - x).abs() < 1e-12 && b == 0.0;
    assert!(close(
        bare_value("FSolve[Plus[-1, Times[-1, x], Power[x, 5]], x]"),
        1.167_303_978_261_418_7
    ));
    assert!(close(
        bare_value("FSolve[Equation[Cos[x], x], x, 0]"),
        0.739_085_133_215_160_7
    ));
    assert!(close(
        bare_value("FSolve[Equation[Log[x], 1], x, 2]"),
        std::f64::consts::E
    ));
    assert!(close(
        bare_value("FSolve[Plus[-2, Power[x, 2]], x, 1]"),
        std::f64::consts::SQRT_2
    ));
    let r2 = std::f64::consts::SQRT_2;
    check_values(&[
        ("FSolve[Plus[-2, Power[x, 2]], x]", &[(-r2, 0.0), (r2, 0.0)]),
        (
            "FSolve[Plus[Times[-1, x], Power[x, 3]], x]",
            &[(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0)],
        ),
        ("FSolve[Plus[1, Power[x, 2]], x]", &[]),
        (
            "FSolve[Equation[Cos[x], x], x]",
            &[(0.739_085_133_215_160_7, 0.0)],
        ),
        (
            "FSolve[Equation[Cos[x], x], Equation[x, Range[0, 1]]]",
            &[(0.739_085_133_215_160_7, 0.0)],
        ),
        (
            "FSolve[Equation[Sin[x], 0], Equation[x, Range[2, 4]]]",
            &[(std::f64::consts::PI, 0.0)],
        ),
        (
            "FSolve[Equation[Exp[x], 3], x]",
            &[(1.098_612_288_668_109_8, 0.0)],
        ),
        (
            "FSolve[Equation[Times[x, Exp[x]], 1], x]",
            &[(0.567_143_290_409_783_8, 0.0)],
        ),
        (
            "FSolve[Equation[Tan[x], x], Equation[x, Range[4, 23/5]]]",
            &[(4.493_409_457_909_064, 0.0)],
        ),
        // No root in the bracket; a pole of tan is no root.
        (
            "FSolve[Equation[Power[x, 2], 2], Equation[x, Range[0, 1]]]",
            &[],
        ),
        ("FSolve[Tan[x], Equation[x, Range[1, 2]]]", &[]),
    ]);
    // Newton's iteration from 1 on x^2 + 1 finds nothing; parameters.
    check(&[
        (
            "FSolve[Plus[1, Power[x, 2]], x, 1]",
            "FSolve[Plus[1, Power[x, 2]], x, 1]",
        ),
        (
            "FSolve[Plus[Times[-1, a], Power[x, 2]], x]",
            "FSolve[Plus[Times[-1, a], Power[x, 2]], x]",
        ),
    ]);
}
