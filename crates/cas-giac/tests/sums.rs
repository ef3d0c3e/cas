//! Sums and products against giac (CAS-DOC-13, CAS-PLAN-7): giac's `sum`
//! and `product`, read back, are this CAS's by value, the parameters and
//! the bounds at a point where they are symbols.

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

/// The symbols' values where answers are compared: `n` an integer, for
/// giac's factorials in it (`binomial(n,2)` is `n!/(2 (n-2)!)`), which
/// have no value at a fraction here.
const POINT: [(&str, i64, i64); 6] = [
    ("a", 3, 7),
    ("b", 5, 11),
    ("c", -2, 13),
    ("m", 4, 3),
    ("n", 17, 1),
    ("x", 2, 9),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: the same, equal by `==`, or equal by `==`
/// once approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let equal = |a: &Expr, b: &Expr| holds(ctx, &Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]));
    if a == b || equal(a, b) {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = POINT
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let at = |e: &Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
            .ok()
    };
    match (at(a), at(b)) {
        (Some(a), Some(b)) => a == b || equal(&a, &b),
        _ => false,
    }
}

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

/// Finite bounds, Karr's convention and polynomial antidifferences
/// (CAS-224).
#[test]
fn polynomial_sums() {
    check(&[
        "sum(k,k,1,n)",
        "sum(k^2,k,1,n)",
        "sum(k^3,k,1,n)",
        "sum(k^4,k,1,n)",
        "sum(k^2,k,0,n)",
        "sum(k^2,k,m,n)",
        "sum(a*k+b,k,1,n)",
        "sum(1,k,1,n)",
        "sum(k^2,k=1..n)",
        "sum(k,k,1,10)",
        "sum(k,k,5,1)",
        "sum(k,k,5,4)",
        "product(k,k,1,5)",
        "product(k,k,5,1)",
        "sum(k!,k,0,n)",
    ]);
}

/// Binomials and hypergeometric terms by Gosper's algorithm (CAS-225).
#[test]
fn hypergeometric_sums() {
    check(&[
        "binomial(5,2)",
        "comb(5,2)",
        "binomial(n,2)",
        "binomial(-3,2)",
        "sum(2^k,k,0,n)",
        "sum(x^k,k,0,n)",
        "sum(k*2^k,k,0,n)",
        "sum(k^2*2^k,k,0,n)",
        "sum((k+1)*2^k,k,0,n-1)",
        "sum(1/(k*(k+1)),k,1,n)",
        "sum(k/(k+1)!,k,0,n)",
        "sum(1/((2*k-1)*(2*k+1)),k,1,n)",
        "sum(1/(k*(k+1)*(k+2)),k,1,n)",
        "sum(1/(k^2-1),k,2,n)",
        "sum(k*k!,k,0,n)",
    ]);
}

/// Infinite sums (CAS-226): limits of antidifferences, known series and
/// divergence.
#[test]
fn infinite_sums() {
    check(&[
        "sum(1/2^k,k,0,inf)",
        "sum(1/(k*(k+1)),k,1,inf)",
        "sum(1/(k^2-1),k,2,inf)",
        "sum(k/2^k,k,1,inf)",
        "sum(1/k!,k,0,inf)",
        "sum(3^k/k!,k,0,inf)",
        "sum((-1)^k/k,k,1,inf)",
        "sum(1/k^2,k,1,inf)",
        "sum(1/k^4,k,1,inf)",
        "sum(1,k,1,inf)",
        "sum(1/k,k,1,inf)",
        "sum(1/k^2,k,2,inf)",
        "sum(1/k^6,k,1,inf)",
        "sum(2/k^2,k,1,inf)",
        "sum((-1)^k/k^2,k,1,inf)",
        "sum(1/k!,k,2,inf)",
        "sum((-1)^k/k!,k,0,inf)",
        "sum(x^k/k!,k,0,inf)",
        "sum(1/(2^k*k),k,1,inf)",
        "sum(1/(3^k*k),k,2,inf)",
        "sum(k^2/2^k,k,0,inf)",
        "sum(1/(k*(k+2)),k,1,inf)",
        "sum(1/(4*k^2-1),k,1,inf)",
        "sum((1/2)^k,k,1,inf)",
        "sum((-1/2)^k,k,0,inf)",
        "sum(k,k,1,inf)",
        "sum(-1,k,1,inf)",
    ]);
    // Odd zeta values stay, where giac has Zeta(3).
    let stays = ours("sum(1/k^3,k,1,inf)").unwrap();
    assert_eq!(stays.head(), Some(&Expr::from(Symbol::SUM)), "{stays:?}");
}
