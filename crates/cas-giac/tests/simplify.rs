//! Simplification against giac (CAS-DOC-14, CAS-PLAN-8): giac's `simplify`
//! and its rewritings, read back, are this CAS's by value, the symbols at a
//! point. giac's forms are not this CAS's (`sin(x)*cos(x)` is
//! `tan(x)/(tan(x)^2+1)` there), only their values are.

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

/// The symbols' values where answers are compared: `k` and `n` integers,
/// for factorials in them, which have no value at a fraction here.
const POINT: [(&str, i64, i64); 7] = [
    ("a", 3, 7),
    ("b", 5, 11),
    ("c", -2, 13),
    ("k", 5, 1),
    ("n", 17, 1),
    ("x", 2, 9),
    ("y", 7, 3),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists part by part, the same, equal by `==`,
/// or equal by `==` once approximated at [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
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

/// Rational normalization, inside functions and piece by piece (CAS-229).
#[test]
fn rational_simplification() {
    check(&[
        "simplify((x^2-1)/(x-1))",
        "simplify(x/(x^2+x))",
        "simplify((x^3-1)/(x^2+x+1))",
        "simplify((x+1)^2-x^2-2*x)",
        "simplify(1/(1/x+1/y))",
        "simplify((a^2-b^2)/(a+b))",
        "simplify(sqrt(x)*sqrt(x))",
        "simplify(sqrt(x^2))",
        "simplify(2^(1/3)*4^(1/3))",
        "simplify([(x^2-1)/(x-1),x/x])",
        "simplify(sin((x^2-1)/(x-1)))",
    ]);
}

/// Exponentials and logarithms (CAS-230). `ln(x^2)-2*ln(x)` under
/// `assume(x>0)` is the unit tests': an assumption would stay in giac's
/// shared session.
#[test]
fn exponentials_and_logarithms() {
    check(&[
        "simplify(exp(2*ln(x)))",
        "simplify(exp(x+ln(2)))",
        "simplify(exp(x)/exp(x-1))",
        "simplify(exp(x)*exp(y))",
        "simplify(ln(2)+ln(3))",
        "simplify(ln(4)-2*ln(2))",
        "simplify(ln(8)/ln(2))",
    ]);
}
