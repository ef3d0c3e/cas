//! Inequalities against giac (CAS-DOC-19, CAS-PLAN-13): the answers, lists
//! of conditions on `x`, are compared as sets, at every end of either
//! answer, the midpoints between and a point beyond each extreme.

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

/// The ends an answer names: the sides of its comparisons free of `x`, and
/// its points.
fn ends(answer: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    for part in answer.args() {
        if !part.contains(x) {
            out.push(part.clone());
        } else if part != x {
            ends(part, x, out);
        }
    }
}

/// Whether `v` meets the answer: some condition holds at it, `x` itself
/// every point, an expression free of `x` the point it is.
fn meets(ctx: &Context, answer: &Expr, x: &Expr, v: &Expr) -> Option<bool> {
    let mut any = false;
    for c in answer.args() {
        let test = if c == x {
            Expr::from(Symbol::TRUE)
        } else if !c.contains(x) {
            Expr::call(Symbol::EQUAL, [c.clone(), v.clone()])
        } else {
            c.substitute(&[(x.clone(), v.clone())])
        };
        match ctx.eval(&test).ok()?.as_symbol()? {
            Symbol::TRUE => any = true,
            Symbol::FALSE => {}
            _ => return None,
        }
    }
    Some(any)
}

/// Whether two answers are one set, sampled at the ends of both, the
/// midpoints and beyond.
fn same_set(ctx: &Context, a: &Expr, b: &Expr, x: &Expr) -> bool {
    let mut points = Vec::new();
    ends(a, x, &mut points);
    ends(b, x, &mut points);
    // Sorted by their approximations.
    let mut valued: Vec<(f64, Expr)> = points
        .into_iter()
        .filter_map(|p| Some((approximate(ctx, &p)?, p)))
        .collect();
    valued.sort_by(|p, q| p.0.total_cmp(&q.0));
    let mut samples: Vec<Expr> = valued.iter().map(|(_, p)| p.clone()).collect();
    for w in valued.windows(2) {
        samples.push((w[0].1.clone() + w[1].1.clone()) / Expr::from(2));
    }
    match (valued.first(), valued.last()) {
        (Some(lo), Some(hi)) => {
            samples.push(lo.1.clone() - Expr::from(1));
            samples.push(hi.1.clone() + Expr::from(1));
        }
        _ => samples.push(Expr::from(0)),
    }
    samples.iter().all(|v| {
        matches!(
            (meets(ctx, a, x, v), meets(ctx, b, x, v)),
            (Some(first), Some(second)) if first == second
        )
    })
}

/// A number's approximation as an `f64`, through `Floor` of it times 2^32.
fn approximate(ctx: &Context, n: &Expr) -> Option<f64> {
    let scaled = Expr::call(Symbol::FLOOR, [n.clone() * Expr::from(1_i64 << 32)]);
    // Below 2^31 in magnitude: an order between the ends, not their value.
    let k = ctx.eval(&scaled).ok()?.as_i64()?;
    #[allow(clippy::cast_precision_loss, reason = "an order, not a value")]
    Some(k as f64 / (1_i64 << 32) as f64)
}

/// giac's answers and this CAS's to `inputs`, in `x`, the same sets.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same_set(&ctx, &ours, &theirs, &x) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// One rational inequality (CAS-277).
#[test]
fn sign_charts() {
    check(&[
        "solve(x^2>4,x)",
        "solve(x^2<4,x)",
        "solve(x^2>=4,x)",
        "solve(2*x+1>0,x)",
        "solve(-3*x+6<=0,x)",
        "solve((x-1)/(x+2)>=0,x)",
        "solve((x-1)/(x+2)<0,x)",
        "solve(x^3-x>0,x)",
        "solve(x^2+1>0,x)",
        "solve(x^2+1<0,x)",
        "solve(x^2-2*x-3<=0,x)",
        "solve(1/x>1,x)",
        "solve(x^2<2,x)",
        "solve(x!=1,x)",
        "solve((x-1)^2>0,x)",
        "solve((x-1)^2>=0,x)",
        "solve((x-1)^2<=0,x)",
        "solve(x^4-5*x^2+4<0,x)",
    ]);
}

/// Systems and connectives (CAS-278).
#[test]
fn systems() {
    check(&[
        "solve([x>1,x<3],x)",
        "solve(x^2-x>0 and x<5,x)",
        "solve(x<-1 or x>2,x)",
        "solve(-1<x<=3,x)",
        "solve([x^2<9,x!=0],x)",
        "solve(not(x>1),x)",
        "solve([x>3,x<1],x)",
        "solve(x^2>1 or x<0,x)",
        "solve([x>=0,x^2=4],x)",
    ]);
}

/// Absolute values, radicals, exponentials and logarithms (CAS-279).
#[test]
fn functions() {
    check(&[
        "solve(abs(x-1)<2,x)",
        "solve(abs(x-1)>2,x)",
        "solve(exp(x)>1,x)",
        "solve(ln(x)<1,x)",
        "solve(sqrt(x)<2,x)",
        "solve(sqrt(x-1)>=1,x)",
        "solve(abs(x)<=x^2,x)",
        "solve(x^(1/3)<2,x)",
        "solve(1/sqrt(x)>1,x)",
        "solve(abs(x-1)+abs(x+1)<4,x)",
        "solve(ln(x)>0,x)",
        "solve(exp(-x)<2,x)",
        "solve(sqrt(x)>x,x)",
        "solve(sqrt(4-x^2)>=1,x)",
        "solve(ln(x^2-1)<0,x)",
        "solve(abs(2*x-3)<=1,x)",
        "solve(x*exp(x)>0,x)",
        "solve(sqrt(x+2)>x,x)",
        "solve(exp(x)-2>0,x)",
        "solve(2^x>8,x)",
        "solve(abs(x^2-4)<3,x)",
        "solve(sqrt(x)+sqrt(x-1)<3,x)",
    ]);
}
