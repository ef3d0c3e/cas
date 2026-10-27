//! Simplification against giac (CAS-DOC-14, CAS-PLAN-8): giac's `simplify`
//! and its rewritings, read back, are this CAS's by value, the symbols at a
//! point. giac's forms are not this CAS's (`sin(x)*cos(x)` is
//! `tan(x)/(tan(x)^2+1)` there), only their values are.

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

/// The parts of a number, a complex ball's real and imaginary ones apart:
/// the modulus of a complex ball about 0 is unbounded.
fn parts(e: &Expr) -> Vec<Expr> {
    match e.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![e.clone()],
    }
}

/// Whether two values are one: lists part by part, the same, equal by `==`,
/// or equal by `==` or with each part of their difference within 10⁻¹⁰
/// once approximated at [`POINT`]: the values of `trig2exp` are complex
/// balls, which `==` does not tell apart.
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
    let close = |a: &Expr, b: &Expr| {
        let Ok(difference) = ctx.eval(&(a.clone() - b.clone())) else {
            return false;
        };
        let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
        parts(&difference).into_iter().all(|p| {
            let gap = Expr::call(Symbol::ABS, [p]);
            holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
        })
    };
    match (at(a), at(b)) {
        (Some(a), Some(b)) => a == b || equal(&a, &b) || close(&a, &b),
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

/// Compositions with the inverse functions, in evaluation (CAS-231).
#[test]
fn compositions() {
    check(&[
        "sin(asin(x))",
        "cos(asin(x))",
        "tan(asin(x))",
        "sin(acos(x))",
        "cos(acos(x))",
        "tan(acos(x))",
        "sin(atan(x))",
        "cos(atan(x))",
        "tan(atan(x))",
        "ln(exp(x))",
        "ln(exp(2))",
    ]);
}

/// Radicals (CAS-232): denominators rationalized, square roots denested.
#[test]
fn radicals() {
    check(&[
        "simplify(1/(1+sqrt(2)))",
        "simplify(1/(sqrt(3)-sqrt(2)))",
        "simplify(1/(1+sqrt(2)+sqrt(3)))",
        "simplify(sqrt(3+2*sqrt(2)))",
        "simplify(sqrt(5+2*sqrt(6)))",
        "simplify((sqrt(x)+1)*(sqrt(x)-1))",
        "simplify(sqrt(x)/x)",
    ]);
}

/// The trigonometric rewritings (CAS-233): giac's `texpand`, `tlin`,
/// `trig2exp` and `exp2trig`.
#[test]
fn trigonometric_rewritings() {
    check(&[
        "texpand(sin(2*x))",
        "texpand(cos(3*x))",
        "texpand(sin(x+y))",
        "texpand(tan(x+y))",
        "texpand(cos(2*x+y))",
        "tlin(sin(x)*cos(x))",
        "tlin(sin(x)^2)",
        "tlin(cos(x)^3)",
        "tlin(sin(x)*sin(y))",
        "tlin(sin(x)^2*cos(x))",
        "trig2exp(sin(x))",
        "trig2exp(tan(x))",
        "exp2trig(exp(i*x))",
        "exp2trig(trig2exp(cos(x)))",
    ]);
}

/// Trigonometric simplification (CAS-234).
#[test]
fn trigonometric_simplification() {
    check(&[
        "simplify(sin(x)^2+cos(x)^2)",
        "simplify(sin(2*x)/sin(x))",
        "simplify(tan(x)*cos(x))",
        "simplify(sin(x)/cos(x))",
        "simplify((1+tan(x)^2)*cos(x)^2)",
        "simplify(sin(x+y)-sin(x)*cos(y)-cos(x)*sin(y))",
        "simplify(1/sin(x)^2-1/tan(x)^2)",
        "simplify(cos(2*x)+2*sin(x)^2)",
        "simplify(sin(x)^2/(1-cos(x)))",
        "simplify(sin(x)^3+sin(x)*cos(x)^2)",
        "simplify(cos(x)/(1-sin(x)^2))",
        "simplify(sin(x)*tan(x)+cos(x))",
        "simplify(tan(2*x)*(1-tan(x)^2))",
        "simplify(sin(3*x)/sin(x))",
        "simplify(sin(x)*cos(y)+cos(x)*sin(y))",
    ]);
}

/// Factorials and binomials (CAS-235).
#[test]
fn factorials() {
    check(&[
        "simplify(n!/(n-1)!)",
        "simplify((n+1)!/n!)",
        "simplify((n+2)!/n!)",
        "simplify((2*n)!/(2*n-2)!)",
        "simplify(n!/(n+1)!)",
        "simplify(binomial(n,k)/binomial(n,k-1))",
        "simplify(binomial(n+1,k)/binomial(n,k))",
        "simplify(n!-n*(n-1)!)",
        "simplify((n+2)/(n+2)!)",
    ]);
    // giac keeps (n+2)/(n+2)!: here it is 1/(n+1)!.
    let found = ours("simplify((n+2)/(n+2)!)").unwrap();
    let expected = ours("1/(n+1)!").unwrap();
    assert_eq!(found, expected);
}
