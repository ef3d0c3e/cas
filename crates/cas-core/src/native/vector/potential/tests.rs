//! Tests of the potentials (CAS-409): giac's examples, the definitions,
//! and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// `Simplify` of `e` is `expected`.
fn simplifies_to(ctx: &Context, e: Expr, expected: &str) -> bool {
    ctx.eval(&Expr::call(Symbol::SIMPLIFY, [e])).unwrap() == eval(ctx, expected)
}

#[test]
fn potentials() {
    let ctx = Context::new();
    for (field, vars, expected) in [
        (
            "List[Plus[Times[2, x, y], 3], Plus[Power[x, 2], Times[-4, z]], Times[-4, y]]",
            "List[x, y, z]",
            "Plus[Times[Power[x, 2], y], Times[3, x], Times[-4, y, z]]",
        ),
        (
            "List[Times[2, x, y], Plus[Power[x, 2], Times[3, Power[y, 2]]]]",
            "List[x, y]",
            "Plus[Times[Power[x, 2], y], Power[y, 3]]",
        ),
        (
            "List[y, Plus[x, z], y]",
            "List[x, y, z]",
            "Plus[Times[x, y], Times[y, z]]",
        ),
        (
            "List[Times[Cos[x], Exp[y]], Times[Sin[x], Exp[y]]]",
            "List[x, y]",
            "Times[Sin[x], Exp[y]]",
        ),
    ] {
        let u = eval(&ctx, &format!("Potential[{field}, {vars}]"));
        let difference = u.clone() - eval(&ctx, expected);
        assert!(
            simplifies_to(&ctx, difference, "0"),
            "Potential[{field}]: {u:?}"
        );
        // Its gradient is the field.
        let gradient = Expr::call(Symbol::GRAD, [u.clone(), read(vars)]);
        let difference = ctx
            .eval(&Expr::call(Symbol::SIMPLIFY, [gradient - read(field)]))
            .unwrap();
        assert!(
            difference.args().iter().all(|a| *a == Expr::from(0)),
            "grad Potential[{field}] = {u:?}: {difference:?}"
        );
    }
}

#[test]
fn vector_potentials() {
    let ctx = Context::new();
    for field in [
        "List[x, y, Times[-2, z]]",
        "List[Plus[Times[2, x, y], 3], Plus[Power[x, 2], Times[-4, z]], Times[-2, y, z]]",
        "List[Times[y, Exp[z]], Times[x, Sin[z]], Plus[x, y]]",
    ] {
        let a = eval(&ctx, &format!("VectorPotential[{field}, List[x, y, z]]"));
        assert_eq!(a.head(), Some(&Expr::from(Symbol::LIST)), "{field}: {a:?}");
        let curl = Expr::call(Symbol::CURL, [a.clone(), read("List[x, y, z]")]);
        let difference = ctx
            .eval(&Expr::call(Symbol::SIMPLIFY, [curl - read(field)]))
            .unwrap();
        assert!(
            difference.args().iter().all(|d| *d == Expr::from(0)),
            "curl VectorPotential[{field}] = {a:?}: {difference:?}"
        );
    }
    // giac's: any two differ by a gradient.
    assert_eq!(
        eval(
            &ctx,
            "VectorPotential[List[x, y, Times[-2, z]], List[x, y, z]]"
        ),
        eval(&ctx, "List[Times[y, z], Times[-1, x, z], 0]")
    );
}

/// What stays: a field that is no gradient, or no curl.
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Potential[List[y, Times[2, x]], List[x, y]]",
        "VectorPotential[List[x, y, z], List[x, y, z]]",
        "Potential[List[x, y], List[x]]",
        "VectorPotential[List[x, y], List[x, y]]",
    ] {
        let e = read(input);
        assert_eq!(eval(&ctx, input).head(), e.head(), "{input}");
    }
}
