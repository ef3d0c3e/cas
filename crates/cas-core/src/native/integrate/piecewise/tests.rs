//! Tests of the piecewise splitting's steps (CAS-316): orders, zeros and
//! the integral under assumptions, apart, so that a failure names its
//! step.

use std::cmp::Ordering;

use super::{integral, order, zeros};
use crate::eval::{Eval, run};
use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `|x − c|` over `[0, 1]` for `0 < c < 1`: `0 < c < 1` told, the zero
/// `c`, and the integral `c²/2 + (1 − c)²/2`.
#[test]
fn a_symbolic_zero_between_the_bounds() {
    let ctx = Context::new().assuming(read("And[Greater[c, 0], Less[c, 1]]"));
    let (zero, one, c, x) = (
        Expr::from(0),
        Expr::from(1),
        Expr::symbol("c"),
        Expr::symbol("x"),
    );
    let e = read("Integrate[Abs[Plus[x, Times[-1, c]]], x, 0, 1]");
    let f = ctx.eval(&read("Abs[Plus[x, Times[-1, c]]]")).unwrap();
    let u = f.args()[0].clone();
    run(ctx.limits(), |recursion| {
        let mut ev = Eval::new(&ctx, recursion);
        assert_eq!(order(&mut ev, &zero, &c).unwrap(), Some(Ordering::Less));
        assert_eq!(order(&mut ev, &c, &one).unwrap(), Some(Ordering::Less));
        assert_eq!(order(&mut ev, &one, &c).unwrap(), Some(Ordering::Greater));
        let found = zeros(&mut ev, &e, &u, &x, [&zero, &one]).unwrap();
        assert_eq!(found, Some(vec![c.clone()]));
        let value = integral(&mut ev, &e, &f, &x, [&zero, &one]).unwrap();
        let value = value.unwrap_or_else(|| panic!("{e:?} stays"));
        let expected =
            read("Plus[Times[1/2, Power[c, 2]], Times[1/2, Power[Plus[1, Times[-1, c]], 2]]]");
        let difference = Expr::call(crate::Symbol::EXPAND, [value.clone() - expected]);
        assert_eq!(ctx.eval(&difference).unwrap(), Expr::from(0), "{value:?}");
    });
}

/// CAS-650: an `Abs` integral whose upper bound holds the level below, 30
/// deep, stays at once, not approximated again at each level (2^n). It
/// reaches `sign` (piecewise.rs) through `order(0, inner)`, the sign of the
/// held inner call; each level integrates in a variable of its own, so the
/// inner call is never rewritten and stays held, as its bound does.
#[test]
fn nested_bound() {
    let mut text = String::from("c");
    for i in 0..30 {
        text = format!("Integrate[Abs[x{i}], x{i}, 0, {text}]");
    }
    let ctx = Context::new();
    let e = read(&text);
    let start = std::time::Instant::now();
    let found = ctx.eval(&e).unwrap();
    let elapsed = start.elapsed();
    assert_eq!(found.head(), e.head());
    assert!(elapsed.as_secs_f64() < 2.0, "{elapsed:?}");
}
