//! Tests of `NDSolve` (CAS-391 to CAS-393): closed forms, systems, giac's
//! value, the curve, and what stays.

use cas_num::ball::{Ball, Float};

use crate::{Context, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// The ball of `e`, by `N`.
fn value(ctx: &Context, e: &Expr) -> Ball {
    let v = ctx
        .eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
        .unwrap();
    let View::Ball(b) = v.view() else {
        panic!("{e:?}: {v:?}");
    };
    b.clone()
}

/// `found` meets `expected` and has 30 bits.
fn holds(what: &str, found: &Expr, expected: &Ball) {
    let View::Ball(b) = found.view() else {
        panic!("{what}: {found:?}");
    };
    assert!(b.overlaps(expected), "{what}: {b:?} against {expected:?}");
    assert!(b.rel_accuracy_bits() >= 30, "{what}: {b:?}");
}

/// Equations with closed forms, forward and backward, in both of giac's
/// forms.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    for (input, expected) in [
        ("NDSolve[y, List[t, y], List[0, 1], 1]", "E"),
        (
            "NDSolve[Times[-2, t, y], List[t, y], List[0, 1], 3]",
            "Power[E, -9]",
        ),
        (
            "NDSolve[Times[y, Plus[1, Times[-1, y]]], List[t, y], List[0, 1/2], 5]",
            "Power[Plus[1, Power[E, -5]], -1]",
        ),
        ("NDSolve[y, Equation[t, Range[0, 1]], y, 1]", "E"),
        ("NDSolve[y, List[t, y], List[1, E], 0]", "1"),
        (
            "NDSolve[Plus[t, Times[-1, y]], List[t, y], List[0, 1], 2]",
            "Plus[1, Times[2, Power[E, -2]]]",
        ),
    ] {
        let found = eval(&ctx, input);
        assert_eq!(
            found.head(),
            Some(&Expr::from(Symbol::LIST)),
            "{input}: {found:?}"
        );
        holds(input, &found.args()[0], &value(&ctx, &read(expected)));
    }
}

/// Systems: the harmonic oscillator, and Bessel's equation from `t = 1`
/// against `J_0(5)`.
#[test]
fn systems() {
    let ctx = Context::new();
    let found = eval(
        &ctx,
        "NDSolve[List[z, Times[-1, y]], List[t, List[y, z]], List[0, List[0, 1]], 2]",
    );
    let pair = &found.args()[0];
    holds("sin 2", &pair.args()[0], &value(&ctx, &read("Sin[2]")));
    holds("cos 2", &pair.args()[1], &value(&ctx, &read("Cos[2]")));
    let found = eval(
        &ctx,
        "NDSolve[List[z, Plus[Times[-1, z, Power[t, -1]], Times[-1, y]]], List[t, List[y, z]], \
         List[1, List[BesselJ[0, 1], Times[-1, BesselJ[1, 1]]]], 5]",
    );
    let pair = &found.args()[0];
    let five = Ball::exact(Float::from(5i64));
    holds("J_0(5)", &pair.args()[0], &five.bessel_j(0, 64));
    holds("J_0′(5)", &pair.args()[1], &five.bessel_j(1, 64).neg());
}

/// giac's `odesolve(sin(t*y),[t,y],[0,1],2)`, `1.82241255674`, within
/// its integrator's error, which reached `2·10⁻⁹` on closed forms.
#[test]
fn giac_value() {
    let ctx = Context::new();
    let found = eval(&ctx, "NDSolve[Sin[Times[t, y]], List[t, y], List[0, 1], 2]");
    let View::Ball(b) = found.args()[0].view() else {
        panic!("{found:?}");
    };
    let giac = 1.822_412_556_74;
    assert!((b.mid().to_f64() - giac).abs() < 1e-8, "{b:?}");
}

/// The curve: `t0` first, `t1` last, each point `e^t`.
#[test]
fn curve() {
    let ctx = Context::new();
    let found = eval(&ctx, "NDSolve[y, List[t, y], List[0, 1], 1, curve]");
    let rows = found.args();
    assert!(rows.len() >= 2, "{found:?}");
    let View::Ball(t0) = rows[0].args()[0].view() else {
        panic!("{found:?}");
    };
    assert!(t0.overlaps(&Ball::ZERO), "{found:?}");
    for row in rows {
        let [t, ys] = row.args() else {
            panic!("{row:?}");
        };
        let expected = value(&ctx, &Expr::from(Symbol::E).pow(t.clone()));
        holds("e^t", &ys.args()[0], &expected);
    }
    let last = &rows[rows.len() - 1].args()[0];
    holds("t1", last, &value(&ctx, &Expr::from(1)));
}

/// What stays: a blow-up before `t1`, with a warning; a switch; a symbolic
/// start.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "NDSolve[Power[y, 2], List[t, y], List[0, 1], 2]",
        "NDSolve[Abs[y], List[t, y], List[0, 1], 1]",
        "NDSolve[y, List[t, y], List[0, a], 1]",
    ] {
        let found = eval(&ctx, input);
        assert_eq!(
            found.head(),
            Some(&Expr::from(Symbol::NDSOLVE)),
            "{input}: {found:?}"
        );
    }
}

/// Two-point boundary value problems by shooting (CAS-520): each value
/// holds the closed form's, the abscissae exact, the ends the given
/// values; resonance, giac's spline output and parameters stay.
#[test]
fn boundary_value_problems() {
    let ctx = Context::new();
    let cases = [
        // y″ = −y, y(0) = 0, y(1) = 1: sin x / sin 1.
        (
            "BVPSolve[Times[-1, y], List[Equation[x, Range[0, 1]], y], List[0, 1], 4]",
            "Times[Sin[x], Power[Sin[1], -1]]",
        ),
        // y″ = y, y(0) = 1, y(1) = 2: (sinh(1 − x) + 2 sinh x) / sinh 1.
        (
            "BVPSolve[y, List[Equation[x, Range[0, 1]], y], List[1, 2], 4]",
            "Times[Power[Sinh[1], -1], Plus[Sinh[Plus[1, Times[-1, x]]], Times[2, Sinh[x]]]]",
        ),
        // giac's example, y′ as diff(y(x), x) reads: x² + 16/x.
        (
            "BVPSolve[Times[1/8, Plus[32, Times[2, Power[x, 3]], \
             Times[-1, y, Derivative[1][y][x]]]], List[Equation[x, Range[1, 3]], y], \
             List[17, 43/3], 4]",
            "Plus[Power[x, 2], Times[16, Power[x, -1]]]",
        ),
    ];
    for (input, solution) in cases {
        let found = eval(&ctx, input);
        let rows = found.args();
        assert_eq!(rows.len(), 5, "{input}: {found:?}");
        for (k, row) in rows.iter().enumerate() {
            let [x, y] = row.args() else {
                panic!("{input}: {row:?}");
            };
            assert!(x.to_rational().is_some(), "{input}: {x:?}");
            let exact = read(solution).substitute(&[(Expr::symbol("x"), x.clone())]);
            if k == 0 || k == 4 {
                assert!(y.to_rational().is_some(), "{input}: end {y:?}");
                let gap = ctx
                    .eval(&Expr::call(Symbol::N, [exact - y.clone(), Expr::from(30)]))
                    .unwrap();
                let zero = match gap.view() {
                    View::Ball(b) => b.contains_rational(&cas_num::Rational::ZERO),
                    _ => gap.to_rational().is_some_and(|q| q.is_zero()),
                };
                assert!(zero, "{input}: end {y:?}, gap {gap:?}");
            } else {
                holds(input, y, &value(&ctx, &exact));
            }
        }
    }
    // y″ = −y on [0, π] with y(0) = 0, y(π) = 1 has no solution; giac's
    // spline output; a parameter.
    for input in [
        "BVPSolve[Times[-1, y], List[Equation[x, Range[0, Pi]], y], List[0, 1], 2]",
        "BVPSolve[Times[-1, y], List[Equation[x, Range[0, 1]], y], List[0, 1], 2, \
         Equation[output, spline]]",
        "BVPSolve[Times[-1, a, y], List[Equation[x, Range[0, 1]], y], List[0, 1], 2]",
    ] {
        assert_eq!(
            eval(&ctx, input).head().and_then(Expr::as_symbol),
            Some(Symbol::BVP_SOLVE),
            "{input}"
        );
    }
}
