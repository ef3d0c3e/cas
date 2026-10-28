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
