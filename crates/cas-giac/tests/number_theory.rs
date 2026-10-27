//! Number theory against giac (CAS-DOC-17, CAS-PLAN-11): giac's answers,
//! read back, are this CAS's, lists part by part.

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

/// Whether two answers are one: the same, or equal by `==`.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    a == b
        || ctx
            .eval(&Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]))
            .ok()
            .and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
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

/// Primes (CAS-267). giac threads `isprime` over a list with its
/// pseudoprime codes, `[2,2,0]`; here with `true` and `false`.
#[test]
fn primes() {
    check(&[
        "isprime(97)",
        "isprime(91)",
        "isprime(1)",
        "isprime(-7)",
        "isprime(2^61-1)",
        "isprime(2^64+13)",
        "nextprime(100)",
        "nextprime(-5)",
        "nextprime(2^64)",
        "prevprime(100)",
        "prevprime(3)",
    ]);
}

/// Factorisation and divisors (CAS-268).
#[test]
fn factorisation() {
    check(&[
        "ifactors(360)",
        "ifactors(-12)",
        "ifactors(1)",
        "ifactors(2^64+1)",
        "idivis(12)",
        "idivis(-12)",
        "euler(36)",
        "euler(1)",
        "euler(0)",
    ]);
}

/// `ifactor`: giac's product reads back multiplied out, so ours, held, is
/// compared by its value; its structure is `ifactors`'.
#[test]
fn held_factorisations() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "ifactor(360)",
        "ifactor(-12)",
        "ifactor(2^64+1)",
        "ifactor(97)",
    ] {
        let held = ours(input).ok();
        let value = held.as_ref().map(|h| match h.args() {
            [inner] if h.head() == Some(&Expr::from(Symbol::HOLD)) => ctx.eval(inner).ok(),
            _ => Some(h.clone()),
        });
        match (giac.value(input), value.flatten()) {
            (Ok(theirs), Some(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, _) => wrong.push(format!("{input}: ours {held:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
