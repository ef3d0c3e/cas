//! Tests of giac's `period`, `linabs`, `lin`, `shift_phase`,
//! `non_recursive_normal` and `envelope` (CAS-640, CAS-641): giac's answers
//! as its `caseval` prints them, values where the forms differ.

use crate::{Context, Domain, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
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

/// `(input, expected)`: equal values at `x` = each of `points`, within
/// 10⁻⁹.
fn same_values(cases: &[(&str, &str)], points: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        for p in points {
            let gap = read(&format!(
                "Less[N[Abs[ReplaceAll[Plus[placeholderf, Times[-1, {expected}]], Rule[x, {p}]]]], 1/1000000000]"
            ));
            let gap = gap.substitute(&[(Expr::symbol("placeholderf"), found.clone())]);
            if ctx.eval(&gap).unwrap().as_symbol() != Some(Symbol::TRUE) {
                wrong.push(format!(
                    "{input} at x = {p}: {found:?}, expected {expected}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Inputs that stay as they are.
fn stay(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let found = ctx.eval(&read(input)).unwrap();
        if found != read(input) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn periods() {
    check(&[
        // giac's 2*pi, pi, pi*1/3, 2*pi*1/5*5, 2*pi*3, 12.
        ("Period[Sin[t], t]", "Times[2, Pi]"),
        ("Period[Power[Cos[x], 2]]", "Pi"),
        ("Period[Tan[Times[3, x]]]", "Times[1/3, Pi]"),
        (
            "Period[Times[Sin[Times[2, x]], Cos[Times[3, x]]]]",
            "Times[2, Pi]",
        ),
        (
            "Period[Plus[Sin[Times[1/3, x]], Tan[Times[1/2, x]]]]",
            "Times[6, Pi]",
        ),
        (
            "Period[Times[Sin[Times[1/2, Pi, x]], \
             Power[Plus[2, Cos[Plus[Times[1/3, Pi, x], Times[1/4, Pi]]]], -1]]]",
            "12",
        ),
        ("Period[Exp[Times[Complex[0, 1], x]]]", "Times[2, Pi]"),
        // giac's 2*pi*1/abs(a).
        (
            "Period[Sin[Times[a, x]]]",
            "Times[2, Pi, Power[Abs[a], -1]]",
        ),
        // giac's +infinity: not periodic, or not found so.
        ("Period[x]", "Infinity"),
        (
            "Period[Plus[Sin[x], Sin[Times[Power[2, 1/2], x]]], x]",
            "Infinity",
        ),
        ("Period[Abs[Sin[x]]]", "Infinity"),
        ("Period[Exp[x]]", "Infinity"),
        // Free of x: any period, giac's 0.
        ("Period[y]", "0"),
    ]);
    stay(&["Period[]", "Period[x, 2]", "Period[x, y, z]"]);
}

/// The eval fuzz target's input (Fuzz 45773886813): `Together` left this
/// as it is, which `periodic` took for a function of one argument and met
/// again without end, overflowing the stack (CAS-640).
#[test]
fn period_of_together_left_as_it_is() {
    let mut ctx = Context::new();
    ctx.set_domain(Domain::Complex);
    let ctx = ctx
        .assuming(read("Greater[x, 0]"))
        .assuming(read("Unequal[y, 0]"));
    let e = read("Period[s[-3, False, -3][k[ExpExpand[b, Infinity], x, Infinity]]]");
    assert_eq!(ctx.eval(&e), Ok(Expr::from(Symbol::INFINITY)));
}

#[test]
fn linearisations() {
    check(&[
        // giac's exp(x*n), exp(x+y), exp(6*x)+2*exp(4*x)+exp(2*x).
        ("ExpLinearize[Power[Exp[x], n]]", "Exp[Times[n, x]]"),
        ("ExpLinearize[Times[Exp[x], Exp[y]]]", "Exp[Plus[x, y]]"),
        (
            "ExpLinearize[Power[Plus[Power[Exp[x], 3], Exp[x]], 2]]",
            "Plus[Exp[Times[6, x]], Times[2, Exp[Times[4, x]]], Exp[Times[2, x]]]",
        ),
        // giac's sin(x)*exp(2*x), x*exp(2*x)+1, exp(x/2).
        (
            "ExpLinearize[Times[Sin[x], Power[Exp[x], 2]]]",
            "Times[Sin[x], Exp[Times[2, x]]]",
        ),
        (
            "ExpLinearize[Plus[Times[x, Power[Exp[x], 2]], 1]]",
            "Plus[Times[x, Exp[Times[2, x]]], 1]",
        ),
        ("ExpLinearize[Power[Exp[x], 1/2]]", "Exp[Times[1/2, x]]"),
        (
            "ExpLinearize[Cosh[x]]",
            "Plus[Times[1/2, Exp[x]], Times[1/2, Exp[Times[-1, x]]]]",
        ),
    ]);
    stay(&["ExpLinearize[]", "ExpLinearize[x, y]"]);
}

#[test]
fn abs_linearisations() {
    // giac's abs(x), 1/x*abs(x), 1/2*x+1/2*abs(x), abs(x-1),
    // 1/x*abs(x)+abs(x-1), sin(abs(x)).
    same_values(
        &[
            ("AbsLinearize[Times[x, Sign[x]]]", "Abs[x]"),
            ("AbsLinearize[Sign[x]]", "Times[Abs[x], Power[x, -1]]"),
            (
                "AbsLinearize[Times[x, HeavisideTheta[x]]]",
                "Plus[Times[1/2, x], Times[1/2, Abs[x]]]",
            ),
            (
                "AbsLinearize[Piecewise[Less[x, 1], Plus[1, Times[-1, x]], Plus[x, -1]]]",
                "Abs[Plus[x, -1]]",
            ),
            (
                "AbsLinearize[Plus[Times[Sign[Plus[x, -1]], Plus[x, -1]], Sign[x]], x]",
                "Plus[Times[Abs[x], Power[x, -1]], Abs[Plus[x, -1]]]",
            ),
            (
                "AbsLinearize[Plus[Times[Sin[x], HeavisideTheta[x]], \
                 Times[-1, Sin[x], HeavisideTheta[Times[-1, x]]]]]",
                "Sin[Abs[x]]",
            ),
        ],
        &["-7/3", "-2/3", "1/3", "2/3", "5/3", "8/3"],
    );
    check(&[("AbsLinearize[Times[x, Sign[x]]]", "Abs[x]")]);
    stay(&["AbsLinearize[]", "AbsLinearize[x, 2]"]);
}

#[test]
fn shifts_and_normals() {
    check(&[
        // giac's -cos((pi+2*x)/2), sin((pi+2*t)/2), -1/tan((pi+2*u)/2).
        (
            "ShiftPhase[Sin[x]]",
            "Times[-1, Cos[Times[1/2, Plus[Pi, Times[2, x]]]]]",
        ),
        (
            "ShiftPhase[Cos[t]]",
            "Sin[Times[1/2, Plus[Pi, Times[2, t]]]]",
        ),
        (
            "ShiftPhase[Tan[u]]",
            "Times[-1, Power[Tan[Times[1/2, Plus[Pi, Times[2, u]]]], -1]]",
        ),
        ("ShiftPhase[Plus[x, Exp[x]]]", "Plus[x, Exp[x]]"),
        // giac's 1 and x+exp((x^2-1)/(x-1))+1: the argument stays.
        (
            "NonRecursiveTogether[Plus[Times[x, Power[Plus[x, 1], -1]], Power[Plus[x, 1], -1]]]",
            "1",
        ),
        (
            "NonRecursiveTogether[Plus[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]], \
             Exp[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]]]",
            "Plus[1, x, Exp[Times[Plus[-1, Power[x, 2]], Power[Plus[-1, x], -1]]]]",
        ),
        (
            "NonRecursiveTogether[List[Plus[Power[x, -1], 1]]]",
            "List[Times[Plus[x, 1], Power[x, -1]]]",
        ),
    ]);
    stay(&["ShiftPhase[]", "NonRecursiveTogether[]"]);
}

#[test]
fn envelopes() {
    // The curve giac plots: 2 cos(t)^3, 2 sin(t)^3 at t = 1/3.
    let ctx = Context::new();
    for input in [
        "Envelope[Plus[y, Times[x, Tan[t]], Times[-2, Sin[t]]], t]",
        "Envelope[Plus[v, Times[u, Tan[t]], Times[-2, Sin[t]]], List[u, v, t]]",
        "Envelope[Equal[Plus[y, Times[x, Tan[t]]], Times[2, Sin[t]]], List[t]]",
    ] {
        let found = ctx.eval(&read(input)).unwrap();
        let [px, py] = found.args() else {
            panic!("{input}: {found:?}");
        };
        for (p, exact) in [
            (px, "Times[2, Power[Cos[1/3], 3]]"),
            (py, "Times[2, Power[Sin[1/3], 3]]"),
        ] {
            let gap = read(&format!(
                "Less[N[Abs[Plus[ReplaceAll[placeholderp, Rule[t, 1/3]], Times[-1, {exact}]]]], 1/1000000000]"
            ))
            .substitute(&[(Expr::symbol("placeholderp"), p.clone())]);
            assert_eq!(
                ctx.eval(&gap).unwrap().as_symbol(),
                Some(Symbol::TRUE),
                "{input}: {found:?}"
            );
        }
    }
    stay(&["Envelope[x]", "Envelope[x, 2]", "Envelope[x, List[a, b]]"]);
}
