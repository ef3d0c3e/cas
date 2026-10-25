//! Solutions against giac (CAS-DOC-11, CAS-PLAN-5): giac's `solve` and
//! `csolve`, read back, are this CAS's lists of solutions up to their
//! order, each of ours one of giac's by `==`, which takes overlapping balls
//! as equal, or, with parameters, at a point where they have values.

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

/// The parameters' values where solutions with parameters are compared.
const PARAMETERS: [(&str, i64, i64); 5] = [
    ("a", 3, 7),
    ("b", 5, 11),
    ("c", -2, 13),
    ("y", 1, 5),
    ("z", 2, 9),
];

/// Whether `a == b` evaluates to `True` in `ctx`.
fn equal(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    holds(ctx, &Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]))
}

/// Whether `proposition` evaluates to `True` in `ctx`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether numbers agree to the digits giac prints, its decimals being
/// balls at 53 bits: `|a − b| < 10⁻¹⁰ (1 + |b|)`.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(Symbol::ABS, [a.clone() - b.clone()]);
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(10_000_000_000_i64);
    holds(ctx, &Expr::call(Symbol::LESS, [difference, bound]))
}

/// Whether two solutions are one: the same, equal by `==`, or, with the
/// parameters at their values and approximated, equal by `==` or close;
/// tuples part by part.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b || equal(ctx, a, b) {
        return true;
    }
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let rules: Vec<(Expr, Expr)> = PARAMETERS
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let at = |e: &Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
            .ok()
    };
    match (at(a), at(b)) {
        (Some(a), Some(b)) => equal(ctx, &a, &b) || close(ctx, &a, &b),
        _ => false,
    }
}

/// giac's answer to `input` and this CAS's: lists of the same solutions,
/// in any order.
fn agrees(giac: &Giac, input: &str) -> Result<(), String> {
    agrees_with(giac, input, input)
}

/// giac's answer to `asked`, as `evalf(…)` where it writes `rootof`, and
/// this CAS's to `input`: lists of the same solutions, in any order.
fn agrees_with(giac: &Giac, asked: &str, input: &str) -> Result<(), String> {
    let theirs = giac.value(asked).map_err(|e| e.to_string())?;
    let ours = ours(input).map_err(|e| e.to_string())?;
    let differs = || format!("{input}: ours {ours:?}, giac's {theirs:?}");
    let list = Expr::from(Symbol::LIST);
    // fsolve's number alone.
    if theirs.head() != Some(&list) && ours.head() != Some(&list) {
        return if same(&Context::new(), &ours, &theirs) {
            Ok(())
        } else {
            Err(differs())
        };
    }
    if theirs.head() != Some(&list)
        || ours.head() != Some(&list)
        || theirs.args().len() != ours.args().len()
    {
        return Err(differs());
    }
    let ctx = Context::new();
    let mut left: Vec<&Expr> = theirs.args().iter().collect();
    for solution in ours.args() {
        let Some(k) = left.iter().position(|t| same(&ctx, solution, t)) else {
            return Err(differs());
        };
        left.swap_remove(k);
    }
    Ok(())
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| agrees(giac, input).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The equations of CAS-210: polynomial and rational, in one unknown.
#[test]
fn polynomial_and_rational_equations() {
    check(&[
        "solve(x^2-3*x+2,x)",
        "solve(x^2+1=0,x)",
        "csolve(x^2+1=0,x)",
        "solve(x^4-5*x^2+6,x)",
        "solve(a*x+b,x)",
        "solve(x^2=2,x)",
        "solve(1/x=2,x)",
        "solve((x-1)/(x^2-1)=0,x)",
        "solve(x^2+x+1,x)",
        "csolve(x^2+x+1,x)",
        "solve(x^3-6*x^2+11*x-6,x)",
        "solve(2*x+3=x-1,x)",
        "solve(x^2-2*a*x+a^2,x)",
        "solve(x^6-1,x)",
        "csolve(x^6-1,x)",
        "solve(x^2-x-1,x)",
        "solve(2*x^2-3,x)",
        "solve((x+1)/(x-2)=3,x)",
        "solve(x/(x^2+1),x)",
        "solve(1/(x^2+1),x)",
        "solve(x^3-x,x)",
    ]);
}

/// Parameters, constants and complex numbers among the coefficients.
#[test]
fn coefficients() {
    check(&[
        "solve(x^2-a,x)",
        "solve(x^2+b*x+c,x)",
        "solve(a*x^2+b*x+c,x)",
        "solve(x^2+y^2-1,x)",
        "solve(x*(x-a),x)",
        "solve((x-a)*(x-b),x)",
        "solve((x^2-2)*(x-a),x)",
        "solve(x^2-2*x-a,x)",
        "solve(x^2+a^2,x)",
        "csolve(x^2+a^2,x)",
        "solve(x^2*y-y,x)",
        "solve(sqrt(2)*x-1,x)",
        "solve(x^2-2*sqrt(2)*x+2,x)",
        "solve(x^2-pi,x)",
        "solve(x^2+pi,x)",
        "solve(exp(1)*x-1,x)",
        "solve(x-i,x)",
        "solve(x^2+2*i*x-1,x)",
        "csolve(x^2+i,x)",
        "csolve(x^2-2*x+5,x)",
        "csolve((x^2+1)*(x^2+4),x)",
        "csolve((x-2)*(x^2+1)*(x+3),x)",
        "solve(x^2-2.0,x)",
    ]);
}

/// Cubics, quartics and higher degrees (CAS-211): radicals where giac has
/// them, its `rootof` by value, numeric roots against its numbers.
#[test]
fn higher_degrees() {
    check(&[
        "solve(x^3-2,x)",
        "csolve(x^3-2,x)",
        "solve(x^4-2,x)",
        "csolve(x^4-2,x)",
        "csolve(x^4+1,x)",
        "solve(x^4-10*x^2+1,x)",
        "solve(x^6-3*x^3+1,x)",
        "csolve(x^6+x^3+1,x)",
        "solve(x^3-a,x)",
        "solve(x^4-a,x)",
        "solve(x^4+4*x-1,x)",
        "csolve(x^4+4*x-1,x)",
        "solve(x^5-x-1,x)",
        "csolve(x^5-x-1,x)",
        "solve(x^5-5*x+1,x)",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = [
        "solve(x^3+x+1,x)",
        "csolve(x^3+x+1,x)",
        "solve(x^3-3*x+1,x)",
        "csolve(x^3-3*x+1,x)",
        "solve(2*x^3-3*x^2-1,x)",
    ]
    .iter()
    .filter_map(|input| agrees_with(giac, &format!("evalf({input})"), input).err())
    .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Linear systems (CAS-212): giac's tuples, in any order.
#[test]
fn linear_systems() {
    check(&[
        "solve([x+y=1,x-y=3],[x,y])",
        "solve([x+y+z=6,2*x-y+z=3,x+2*y-z=2],[x,y,z])",
        "solve([x+y=1,2*x+2*y=2],[x,y])",
        "solve([x+y=1,x+y=2],[x,y])",
        "solve([a*x+y=1,x-y=b],[x,y])",
        "solve([x+y=a,x-y=b],[x,y])",
        "solve([x+y+z=1,x-y=0],[x,y,z])",
        "solve([x+y=1,x-y=3],[x,y,z])",
        "solve([x+y=1,x-y=3,2*x=4],[x,y])",
        "solve([x+y=1,x-y=3,2*x=5],[x,y])",
        "solve([2*x+4*y=6,3*x+6*y=9],[x,y])",
        "solve([1.5*x+y=2,x-y=1],[x,y])",
        "csolve([x+i*y=1,x-y=0],[x,y])",
        "solve(x+y=1,[x,y])",
        "solve([0=0],[x,y])",
        "solve(x^2-1,[x])",
        "solve([x+1=0],x)",
    ]);
}

/// Equations in functions of the unknown (CAS-213): giac's principal
/// solutions.
#[test]
fn transcendental_equations() {
    check(&[
        "solve(exp(x)=2,x)",
        "solve(ln(x)=3,x)",
        "solve(ln(x)=-1,x)",
        "solve(sin(x)=0,x)",
        "solve(sin(x)=1/2,x)",
        "solve(cos(x)=1/2,x)",
        "solve(tan(x)=1,x)",
        "solve(cos(x)=0,x)",
        "solve(sin(2*x)=1,x)",
        "solve(atan(x)=1,x)",
        "solve(asin(x)=1/2,x)",
        "solve(exp(x)=-1,x)",
        "csolve(exp(x)=-1,x)",
        "solve(sin(x)=2,x)",
        "solve(sqrt(x)=-1,x)",
        "solve(sqrt(x)=3,x)",
        "solve(1/sqrt(x)=2,x)",
        "solve(x^(1/3)=2,x)",
        "solve(sqrt(x^2+1)=2,x)",
        "solve(sqrt(x+1)=x-1,x)",
        "solve(sqrt(x-1)+1=x,x)",
        "solve(sqrt(x)+sqrt(x+5)=5,x)",
        "solve(exp(2*x)-3*exp(x)+2=0,x)",
        "solve(exp(x)+exp(-x)=3,x)",
        "solve(sin(x)^2=1/4,x)",
        "solve(2^x=8,x)",
        "solve(ln(x^2-3)=0,x)",
        "solve(exp(x^2)=e,x)",
    ]);
}

/// Two polynomial equations in two unknowns (CAS-215).
#[test]
fn polynomial_systems() {
    check(&[
        "solve([x^2+y^2=1,x-y=0],[x,y])",
        "solve([x*y=1,x+y=3],[x,y])",
        "solve([x*y=1,x+y=2],[x,y])",
        "solve([x^2+y^2=5,x*y=2],[x,y])",
        "solve([x^2-y=0,x+y=2],[x,y])",
        "solve([y=x^2,y=1],[x,y])",
        "solve([x*y=0,x+y=1],[x,y])",
        "solve([x^2=1,y^2=4],[x,y])",
        "solve([x^2+y^2=1,(x-1)^2+y^2=1],[x,y])",
        "solve([x^3+y=2,x-y=0],[x,y])",
        "solve([x^2+y^2=-1,x=y],[x,y])",
        "csolve([x^2+y^2=-1,x=y],[x,y])",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    let input = "solve([x^2+y=3,x+y^2=5],[x,y])";
    if let Err(wrong) = agrees_with(giac, &format!("evalf({input})"), input) {
        panic!("{wrong}");
    }
}

/// Numeric roots (CAS-214): giac's `fsolve`, by value.
#[test]
fn numeric_solutions() {
    check(&[
        "fsolve(x^5-x-1,x)",
        "fsolve(x^2-2,x)",
        "fsolve(x^3-x,x)",
        "fsolve(x^2+1,x)",
        "fsolve(cos(x)=x,x)",
        "fsolve(cos(x)=x,x,0)",
        "fsolve(cos(x)=x,x=0..1)",
        "fsolve(sin(x)=0,x=2..4)",
        "fsolve(sin(x)=0,x,3)",
        "fsolve(exp(x)=3,x)",
        "fsolve(x*exp(x)=1,x)",
        "fsolve(tan(x)=x,x=4..4.6)",
        "fsolve(x^2=2,x=0..1)",
        "fsolve(ln(x)=1,x,2)",
        "fsolve(x^2-2,x,1)",
    ]);
}

/// Identities and contradictions.
#[test]
fn identities() {
    check(&[
        "solve(x=x,x)",
        "solve(x^2-x=x^2,x)",
        "solve(x^2+1=x^2,x)",
        "solve(1,x)",
    ]);
}
