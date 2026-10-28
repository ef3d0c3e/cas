//! Tests of the compiled right sides (CAS-391): values on balls against
//! `N`, jets against the ball functions' jets, refusals and domains.

use cas_num::Rational;
use cas_num::ball::{Ball, Jet};

use super::Node;
use crate::{Context, Expr, Symbol, View};

const PREC: u64 = 128;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn variables() -> Vec<Expr> {
    vec![Expr::symbol("t"), Expr::symbol("y"), Expr::symbol("z")]
}

fn ball(n: i64, d: i64) -> Ball {
    Ball::from_rational(
        &Rational::new(cas_num::Integer::from(n), cas_num::Integer::from(d)),
        PREC,
    )
}

/// The tree's value at `t = 1/2, y = 1/3, z = 5/4` holds `N` of the
/// expression there.
#[test]
fn values() {
    let ctx = Context::new();
    let at = [ball(1, 2), ball(1, 3), ball(5, 4)];
    for f in [
        "Plus[Times[y, t], Sin[y]]",
        "Times[Power[E, Times[-1, t]], Power[z, 1/2], ArcTan[y]]",
        "Plus[Power[y, 3], Times[-2, Power[z, -2]], Log[z], Pi]",
        "Power[z, t]",
        "Times[Cosh[y], Tanh[t], ArcSin[y], ArcCosh[z]]",
    ] {
        let e = read(f);
        let node = Node::compile(&ctx, &e, &variables(), PREC).unwrap_or_else(|| panic!("{f}"));
        let found = node.ball(&at, PREC).unwrap();
        let rules = read("List[Rule[t, 1/2], Rule[y, 1/3], Rule[z, 5/4]]");
        let expected = ctx
            .eval(&Expr::call(
                Symbol::N,
                [Expr::call(Symbol::REPLACE_ALL, [e, rules]), Expr::from(40)],
            ))
            .unwrap();
        let View::Ball(expected) = expected.view() else {
            panic!("{f}: {expected:?}");
        };
        assert!(
            found.overlaps(expected),
            "{f}: {found:?} against {expected:?}"
        );
    }
}

/// The jet of `e^y sin(t)` at `y = 1/3 + s`, `t = 1/2 + 2s` is the product
/// of the jets.
#[test]
fn jets() {
    let ctx = Context::new();
    let node = Node::compile(
        &ctx,
        &read("Times[Power[E, y], Sin[t]]"),
        &variables(),
        PREC,
    )
    .unwrap();
    let t = Jet::new(vec![
        ball(1, 2),
        ball(2, 1),
        Ball::ZERO,
        Ball::ZERO,
        Ball::ZERO,
    ]);
    let y = Jet::variable(ball(1, 3), 5);
    let z = Jet::constant(ball(5, 4), 5);
    let found = node.jet(&[t.clone(), y.clone(), z], PREC).unwrap();
    let expected = y.exp(PREC).mul(&t.sin_cos(PREC).0, PREC);
    for (a, b) in found.coeffs().iter().zip(expected.coeffs()) {
        assert!(a.overlaps(b), "{a:?} against {b:?}");
    }
}

/// Switches, unknown heads and complex constants are refused; a domain
/// left gives nothing.
#[test]
fn refusals_and_domains() {
    let ctx = Context::new();
    for f in ["Abs[y]", "f[y]", "Times[I, y]", "HeavisideTheta[t]"] {
        assert!(
            Node::compile(&ctx, &read(f), &variables(), PREC).is_none(),
            "{f}"
        );
    }
    let node = Node::compile(&ctx, &read("Log[y]"), &variables(), PREC).unwrap();
    assert!(
        node.ball(&[ball(0, 1), ball(-1, 1), ball(1, 1)], PREC)
            .is_none()
    );
    let node = Node::compile(&ctx, &read("Power[y, 1/2]"), &variables(), PREC).unwrap();
    let jets = [
        Jet::constant(ball(0, 1), 3),
        Jet::constant(ball(-2, 1), 3),
        Jet::constant(ball(1, 1), 3),
    ];
    assert!(node.jet(&jets, PREC).is_none());
}
