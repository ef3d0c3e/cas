//! `minimax` with `rand=` against giac (CAS-599), a departure: giac's
//! skewed starts are random, so its answer depends on its random state, and
//! its 2% stop leaves `x*exp(-x)` `[5,3]` on `0..10` far above the best.
//! giac's answer with its generator seeded in the same call (C0014), as
//! caseval prints it (C0007); ours seeded in the same input, its error
//! `4.35626e-5` (a dense sampling of the prototype's agrees) and
//! `|f − p/q| ≤ E` at 51 points. In its own file: giac's state is shared.

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

/// Whether `|a| < bound` numerically.
fn small(ctx: &Context, a: &Expr, bound: &Expr) -> bool {
    let size = Expr::call(Symbol::ABS, [Expr::call(Symbol::N, [a.clone()])]);
    ctx.eval(&Expr::call(Symbol::LESS, [size, bound.clone()]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// giac seeded with 42 (`0.000287`; 7 gives `0.000191`, unseeded in a
/// fresh instance `0.000256`), ours `4.35626e-5`.
#[test]
fn minimax_rand() {
    const GIAC: &str = "42,[0.000287329610301,(-1.95426626152e-06+1.00045709244*x\
        -0.384606606179*x^2+0.0574912095625*x^3-0.00393208610762*x^4\
        +0.000103144061133*x^5)/(1+0.622027669514*x+0.156597930714*x^2\
        +0.0399586267272*x^3)]";
    let input = "minimax(x*exp(-x),x=0..10,[5,3],rand=30)";
    let mut wrong = Vec::new();
    if let Some(giac) = giac() {
        let seeded = format!("srand(42);{input}");
        let text = giac.eval(&seeded).unwrap_or_default();
        if text != GIAC {
            wrong.push(format!("{seeded}: giac prints {text}, recorded {GIAC}"));
        }
    }
    let ctx = Context::new();
    let found = ours(&format!("[srand(42),{input}]")).ok();
    let answer = found.as_ref().and_then(|e| e.args().get(1));
    match answer.map(Expr::args) {
        Some([error, p]) => {
            let best = Expr::from(435_626) / Expr::from(10_000_000_000_i64);
            if !small(
                &ctx,
                &(error.clone() - best.clone()),
                &(best / Expr::from(1000)),
            ) {
                wrong.push(format!("{input}: error ours {error:?}"));
            }
            let x = Expr::symbol("x");
            let f = ours("x*exp(-x)").unwrap();
            let bound = error.clone() * (Expr::from(1) + Expr::from(10).pow(-6));
            for k in 0..=50 {
                let t = Expr::from(k) / Expr::from(5);
                let gap = (f.clone() - p.clone()).substitute(&[(x.clone(), t)]);
                if !small(&ctx, &gap, &bound) {
                    wrong.push(format!("{input}: |f - p/q| above ours {error:?} at {k}/5"));
                }
            }
        }
        _ => wrong.push(format!("{input}: ours {found:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
