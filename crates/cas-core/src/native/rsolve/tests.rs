//! Tests of `RSolve` (CAS-310): recurrences with constant coefficients,
//! alone and in systems, first order with coefficients in `n`, conditions,
//! the constants as the initial values, what stays. Answers compare by
//! their values at `n = 0 … 6`, closed forms coming in equal shapes.

use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// Whether `a` and `b` agree to `10⁻²⁰` at `n = 0 … 6`.
fn agree(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    (0..=6i64).all(|k| {
        let at = |e: &Expr| e.substitute(&[(read("n"), Expr::from(k))]);
        let claim = Expr::call(
            crate::symbol::Symbol::LESS,
            [
                Expr::call(
                    crate::symbol::Symbol::ABS,
                    [Expr::call(
                        crate::symbol::Symbol::N,
                        [at(a) - at(b), Expr::from(30)],
                    )],
                ),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&claim).unwrap() == read("True")
    })
}

/// Each input's answer `{u(n)}` against the expected sequence.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let ok = match found.args() {
            [u] if found.head() == Some(&read("List")) => agree(&ctx, u, &read(expected)),
            _ => false,
        };
        if !ok {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Constant coefficients, giac's answers (CAS-310).
#[test]
fn constant_coefficients() {
    check(&[
        (
            "RSolve[Equation[u[Plus[n, 1]], Plus[Times[2, u[n]], 1]], u[n], Equation[u[0], 1]]",
            "Plus[Power[2, Plus[n, 1]], -1]",
        ),
        (
            "RSolve[Equation[u[Plus[n, 1]], Plus[Times[2, u[n]], n]], u[n], Equation[u[0], 1]]",
            "Plus[Times[-1, n], Times[2, Power[2, n]], -1]",
        ),
        (
            "RSolve[Equation[Plus[u[Plus[n, 2]], Times[-3, u[Plus[n, 1]]], Times[2, u[n]]], 0], \
             u[n], Equation[u[0], 1], Equation[u[1], 3]]",
            "Plus[Times[2, Power[2, n]], -1]",
        ),
        (
            "RSolve[Equation[Plus[u[Plus[n, 2]], u[n]], 0], u[n], Equation[u[0], 1], \
             Equation[u[1], 0]]",
            "Cos[Times[1/2, n, Pi]]",
        ),
        (
            "RSolve[Equation[u[Plus[n, 1]], Plus[Times[2, u[n]], Power[2, n]]], u[n], \
             Equation[u[0], 0]]",
            "Times[1/2, n, Power[2, n]]",
        ),
        // A shift below n, moved up.
        (
            "RSolve[Equation[u[n], Plus[Times[2, u[Plus[n, -1]]], 1]], u[n], Equation[u[0], 0]]",
            "Plus[Power[2, n], -1]",
        ),
        // Fibonacci.
        (
            "RSolve[Equation[u[Plus[n, 2]], Plus[u[Plus[n, 1]], u[n]]], u[n], \
             Equation[u[0], 0], Equation[u[1], 1]]",
            "Times[Power[5, -1/2], Plus[Power[Times[1/2, Plus[1, Power[5, 1/2]]], n], \
             Times[-1, Power[Times[1/2, Plus[1, Times[-1, Power[5, 1/2]]]], n]]]]",
        ),
        // First order with coefficients in n: n!.
        (
            "RSolve[Equation[u[Plus[n, 1]], Times[Plus[n, 1], u[n]]], u[n], Equation[u[0], 1]]",
            "Factorial[n]",
        ),
    ]);
}

/// The constants are the initial values; a system; conditions none meets
/// (CAS-310).
#[test]
fn constants_and_systems() {
    let ctx = Context::new();
    // u(n + 1) = 2u(n): c_0 2ⁿ.
    let found = ctx
        .eval(&read(
            "RSolve[Equation[u[Plus[n, 1]], Times[2, u[n]]], u[n]]",
        ))
        .unwrap();
    let expected = ctx.eval(&read("List[Times[c_0, Power[2, n]]]")).unwrap();
    assert_eq!(found, expected);
    // u(n + 2) = 4u(n): its values at 0 and 1 are c_0 and c_1.
    let found = ctx
        .eval(&read(
            "RSolve[Equation[u[Plus[n, 2]], Times[4, u[n]]], u[n]]",
        ))
        .unwrap();
    let [u] = found.args() else {
        panic!("{found:?}");
    };
    for (k, c) in [(0, "c_0"), (1, "c_1")] {
        let at = ctx
            .eval(&Expr::call(
                crate::symbol::Symbol::EXPAND,
                [u.substitute(&[(read("n"), Expr::from(k))])],
            ))
            .unwrap();
        assert_eq!(at, read(c), "u({k}) in {found:?}");
    }
    // [u(n+1) = u(n) + v(n), v(n+1) = v(n)], u(0) = 0, v(0) = 1: {{n, 1}}.
    let found = ctx
        .eval(&read(
            "RSolve[List[Equation[u[Plus[n, 1]], Plus[u[n], v[n]]], Equation[v[Plus[n, 1]], v[n]]], \
             List[u[n], v[n]], Equation[u[0], 0], Equation[v[0], 1]]",
        ))
        .unwrap();
    assert_eq!(found, read("List[List[n, 1]]"));
    // Conditions none meets.
    let found = ctx
        .eval(&read(
            "RSolve[Equation[u[Plus[n, 1]], u[n]], u[n], Equation[u[0], 1], Equation[u[1], 2]]",
        ))
        .unwrap();
    assert_eq!(found, read("List[]"));
    for input in [
        // Not linear; second order with coefficients in n.
        "RSolve[Equation[u[Plus[n, 1]], Power[u[n], 2]], u[n]]",
        "RSolve[Equation[u[Plus[n, 2]], Times[n, u[n]]], u[n]]",
        "RSolve[Equation[u[Plus[n, 1]], u[n]], u[2]]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
