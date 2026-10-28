//! Z-transforms and recurrences against giac (CAS-DOC-24 §5, CAS-PLAN-18):
//! transforms as functions, at sample points of `z` and of the
//! parameters, since giac leaves `cos(1)^2+sin(1)^2` in its denominators;
//! sequences at `n = 0 … 6`, since giac writes `Heaviside(n-1)` and
//! complex exponentials where this CAS writes Kronecker deltas and powers.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// Whether two numbers agree to 10⁻²⁰, real and imaginary parts apart.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// `e` with each symbol of `values` put in.
fn at(e: &Expr, values: &[(&str, &str)]) -> Expr {
    let rules: Vec<(Expr, Expr)> = values
        .iter()
        .map(|(name, value)| (Expr::symbol(name), ours(value).expect("a number")))
        .collect();
    e.substitute(&rules)
}

/// Whether two transforms agree at the sample points.
fn same_function(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    [("7/3", "1/3", "2/7"), ("5/2", "3/4", "5/3")]
        .iter()
        .all(|&(z, a_value, w)| {
            let values = [("z", z), ("a", a_value), ("w", w), ("c", "3/5")];
            close(ctx, &at(a, &values), &at(b, &values))
        })
}

/// giac's transforms and this CAS's to `inputs`, as functions.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same_function(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The table (CAS-308).
#[test]
fn transforms() {
    check(&[
        "ztrans(1,n,z)",
        "ztrans(c,n,z)",
        "ztrans(n,n,z)",
        "ztrans(n^2,n,z)",
        "ztrans(a^n,n,z)",
        "ztrans(n*a^n,n,z)",
        "ztrans(2^n+3^n,n,z)",
        "ztrans(sin(n),n,z)",
        "ztrans(cos(w*n),n,z)",
        "ztrans(1/n!,n,z)",
        "ztrans(binomial(n,2),n,z)",
        "ztrans(n^3*2^n,n,z)",
        "ztrans(2^n/n!,n,z)",
    ]);
}

/// Whether this CAS's sequence for each of `inputs`, giac expressions in
/// `n`, is giac's at `n = 0 … 6`. giac puts in `n` itself, its
/// `Heaviside(0)` being 1 where this CAS's is 1/2.
fn check_sequences(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let Ok(sequence) = ours(input) else {
            wrong.push(format!("{input}: unreadable here"));
            continue;
        };
        for k in 0..=6i64 {
            let ours = ctx.eval(&sequence.substitute(&[(Expr::symbol("n"), Expr::from(k))]));
            let theirs = giac.value(&format!("subst({input},n={k})"));
            match (ours, theirs) {
                (Ok(o), Ok(t)) if close(&ctx, &o, &t) => {}
                (o, t) => wrong.push(format!("{input} at n = {k}: ours {o:?}, giac's {t:?}")),
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The inverse transform, by values (CAS-309).
#[test]
fn inverses() {
    check_sequences(&[
        "invztrans(z/(z-1),z,n)",
        "invztrans(z/(z-1)^2,z,n)",
        "invztrans(1/(z-2),z,n)",
        "invztrans(z^2/((z-1)*(z-2)),z,n)",
        "invztrans(z/(z^2+1),z,n)",
        "invztrans(1,z,n)",
        "invztrans(1/z,z,n)",
        "invztrans((z^2+z)/(z-1)^3,z,n)",
        "invztrans(z/(z-1/2),z,n)",
        "invztrans((3*z^2-z)/((z-1)*(z+1/2)^2),z,n)",
    ]);
}
