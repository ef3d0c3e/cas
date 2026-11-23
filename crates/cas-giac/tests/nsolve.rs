//! `nSolve` of a system against giac (CAS-591), a departure: giac's
//! Newton iteration restarts from random points, so its root depends on
//! its random state; ours is `FSolve`'s first real solution. giac's answer
//! with its generator seeded in the same call (C0014), and unseeded, and
//! ours each solve the equations.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The last value of giac's printed sequence, `srand`'s seed before it.
fn last(answer: &str) -> &str {
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => start = i + 1,
            _ => {}
        }
    }
    &answer[start..]
}

/// Whether `[x, y]` solves `x² + y² = 4`, `x y = 1` within `10⁻¹⁰`.
fn solves(ctx: &Context, tuple: &Expr) -> bool {
    let [x, y] = tuple.args() else {
        return false;
    };
    if tuple.head() != Some(&Expr::from(Symbol::LIST)) {
        return false;
    }
    let small = |e: Expr| {
        let size = Expr::call(Symbol::ABS, [e]);
        let bound = Expr::from(1) / Expr::from(10_000_000_000_i64);
        ctx.eval(&Expr::call(Symbol::LESS, [size, bound]))
            .ok()
            .and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
    };
    let sum = x.clone() * x.clone() + y.clone() * y.clone() - Expr::from(4);
    let product = x.clone() * y.clone() - Expr::from(1);
    small(sum) && small(product)
}

/// giac's roots of `[x^2+y^2=4,x*y=1]` seeded (42 gives
/// `[1.93185165258,0.517638090205]`, 7 `[0.517638090205,1.93185165258]`)
/// and unseeded (`[-0.517638090205,-1.93185165258]` from a fresh
/// instance), and ours, `[-1.93185165258,-0.517638090205]`: each a root.
#[test]
fn nsolve_systems() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let system = "nSolve([x^2+y^2=4,x*y=1],[x,y])";
    let mut wrong = Vec::new();
    for asked in [
        "srand(42);nSolve([x^2+y^2=4,x*y=1],[x,y])",
        "srand(7);nSolve([x^2+y^2=4,x*y=1],[x,y])",
        system,
    ] {
        let theirs = giac
            .eval(asked)
            .map_err(|e| e.to_string())
            .and_then(|t| ours(&last(&t).replace("list[", "[")).map_err(|e| e.to_string()));
        if !theirs.as_ref().is_ok_and(|t| solves(&ctx, t)) {
            wrong.push(format!("{asked}: giac's {theirs:?}"));
        }
    }
    match ours(system) {
        Ok(found) if solves(&ctx, &found) => {}
        found => wrong.push(format!("{system}: ours {found:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
