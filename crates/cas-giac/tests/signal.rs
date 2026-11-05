//! Signal processing against giac (CAS-521 …): giac's answers on lists of
//! samples, read back, are this CAS's by value, entry by entry, within
//! `10⁻¹⁰` (giac's doubles; ours exact or balls).

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists part by part, the same, or numbers
/// within `10⁻¹⁰`.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [a.clone() - b.clone()])) else {
        return false;
    };
    let parts = match difference.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![difference],
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts.into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// giac's answers and this CAS's, the same by [`same`].
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The fifteen windows (CAS-521), on five and six samples, and with
/// parameters.
#[test]
fn windows() {
    let mut inputs = Vec::new();
    for w in [
        "bartlett_hann",
        "blackman_harris",
        "blackman",
        "bohman",
        "cosine",
        "gaussian",
        "hamming",
        "hann_poisson",
        "hann",
        "parzen",
        "poisson",
        "riemann",
        "triangle",
        "tukey",
        "welch",
    ] {
        inputs.push(format!("{w}_window([1,1,1,1,1])"));
        inputs.push(format!("{w}_window([1,2,3,4,5,6])"));
    }
    for input in [
        "tukey_window([1,1,1,1,1,1],0.4)",
        "gaussian_window([1,1,1,1,1],0.3)",
        "cosine_window([1,1,1,1,1],2)",
        "poisson_window([1,1,1,1,1],2)",
        "hann_poisson_window([1,1,1,1,1],2)",
        "triangle_window([1,1,1,1,1],1)",
        "triangle_window([1,1,1,1,1],-1)",
        "blackman_window([1,1,1,1,1],0.2)",
        "tukey_window([1,1,1,1,1],1)",
        "tukey_window([1,1,1,1,1],0)",
    ] {
        inputs.push(input.to_owned());
    }
    let inputs: Vec<&str> = inputs.iter().map(String::as_str).collect();
    check(&inputs);
}

/// First-order RC filters (CAS-522).
#[test]
fn rc_filters() {
    check(&[
        "lowpass([1,0,0,0,0],1000)",
        "lowpass([1,1,1,1],1000,10000)",
        "lowpass([1,2,3],1000)",
        "lowpass([5],1000)",
        "lowpass([],1000)",
        "highpass([1,1,1,1],1000,10000)",
        "highpass([1,0,0,0],1000,10000)",
        "highpass([3,1,4,1,5],2000)",
    ]);
}
