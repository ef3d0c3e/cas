//! giac's `period`, `linabs`, `lin`, `shift_phase`, `non_recursive_normal`
//! and `envelope` of the gap ledger (CAS-640, CAS-641) against giac 2.0.
//! giac's answers read back are this CAS's: the same, or equal values at
//! sample points where the forms differ (`linabs` picks the smaller of
//! equal rewritings by size, ours and giac's sizes differ). The departure
//! carries giac's answer as it prints it and this CAS's as decided.

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

/// giac's answers to `inputs`, read back, are this CAS's.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whether `a − b` is within 10⁻⁹ of 0 once approximated.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::N,
        [Expr::call(Symbol::ABS, [a.clone() - b.clone()])],
    );
    let bound = Expr::from(1) / Expr::from(1_000_000_000_i64);
    ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// `e` with `x` replaced by `p`, evaluated.
fn at(ctx: &Context, e: &Expr, x: &str, p: &Expr) -> Expr {
    let e = e.substitute(&[(Expr::symbol(x), p.clone())]);
    ctx.eval(&e).unwrap_or(e)
}

/// giac's answers to `inputs`, read back, and this CAS's take the same
/// values at `x` = −7/3, −2/3, 1/3, 2/3, 5/3, 8/3.
fn same_values(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let points = [(-7, 3), (-2, 3), (1, 3), (2, 3), (5, 3), (8, 3)]
        .map(|(n, d)| Expr::from(n) / Expr::from(d));
    let mut wrong = Vec::new();
    for input in inputs {
        let (theirs, mine) = match (giac.value(input), ours(input)) {
            (Ok(t), Ok(m)) => (t, m),
            (t, m) => {
                wrong.push(format!("{input}: giac {t:?}, ours {m:?}"));
                continue;
            }
        };
        for p in &points {
            if !near(&ctx, &at(&ctx, &theirs, "x", p), &at(&ctx, &mine, "x", p)) {
                wrong.push(format!(
                    "{input} at x = {p:?}: giac {theirs:?}, ours {mine:?}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn periods() {
    check(&[
        "period(sin(x*pi/2)/(2+cos(x*pi/3+pi/4)))",
        "period(sin(x)+sin(sqrt(2)*x),x)",
        "period(sin(t),t)",
        "period(tan(3*x))",
        "period(cos(x)^2)",
        "period(sin(x)^2+cos(x))",
        "period(sin(2*x)*cos(3*x))",
        "period(x)",
        "period(exp(i*x))",
        "period(sin(x/3)+tan(x/2))",
        "period(sin(a*x))",
        "period(sin(t*a),t)",
        "period(sin(x)*cos(x))",
        "period(abs(sin(x)))",
        "period(floor(x)-x)",
        "period(sin(x)+x)",
        "period(cos(x)/(1+sin(x)^2))",
        "period(exp(sin(x)))",
        "period(sin(pi*x)+cos(2*pi*x/3))",
    ]);
}

#[test]
fn exponential_linearisations() {
    check(&[
        "lin(exp(x)^n)",
        "lin((exp(x)^3+exp(x))^2)",
        "lin(exp(x)*exp(y))",
        "lin(exp(x)^2*exp(y))",
        "lin(sin(x)*exp(x)^2)",
        "lin(exp(2*x)/exp(x))",
        "lin(x*exp(x)^2+1)",
        "lin(exp(x)^(1/2))",
        "lin(2^x)",
        "lin(sinh(x)^2)",
        "lin(exp(x)*exp(-x)+cosh(x))",
        "lin([exp(x)*exp(y),exp(x)^2])",
    ]);
}

#[test]
fn abs_linearisations() {
    same_values(&[
        "linabs(x*sign(x))",
        "linabs(sin(x)*Heaviside(x)-sin(x)*Heaviside(-x))",
        "linabs(abs(abs(x-2)+x*abs(x+1)+1))",
        "linabs(sign(x))",
        "linabs(Heaviside(x))",
        "linabs(x*Heaviside(x))",
        "linabs(piecewise(x<0,-x,x))",
        "linabs(piecewise(x<1,1-x,x-1))",
        "linabs(sign(x-1)*(x-1)+sign(x),x)",
        "linabs(abs(x-1)+abs(x+1))",
        "linabs(x*Heaviside(x-2))",
        "linabs(sign(x)*x^2)",
    ]);
}

#[test]
fn shifts_and_normals() {
    check(&[
        "shift_phase(sin(x))",
        "shift_phase(x+sin(x))",
        "shift_phase(cos(t))",
        "shift_phase(tan(u))",
        "shift_phase(exp(x)+cos(2*x))",
        "shift_phase(sin(x)+cos(x))",
        "shift_phase(1/cos(x))",
        "shift_phase(sin(2*x+1))",
        "shift_phase(sin(x)=cos(x))",
        "shift_phase(sin(sin(x)))",
        "non_recursive_normal(sin(x+x)+sin(2*x)+x+x)",
        "non_recursive_normal(sin(2*x)+sin(2*x)+x+x)",
        "non_recursive_normal((x^2-1)/(x-1)+exp((x^2-1)/(x-1)))",
        "non_recursive_normal(x/(x+1)+1/(x+1))",
        "non_recursive_normal(sin(x)/x+sin(x)/x)",
        "non_recursive_normal([x/(x+1)+1/(x+1),(x^2-1)/(x-1)])",
        "non_recursive_normal(sin((x^2-1)/(x-1))+1/x+1)",
    ]);
    // giac's x+abs(x+1): its normal takes the square root of the square.
    same_values(&["non_recursive_normal(sqrt(x^2+2*x+1)+x)"]);
}

/// Departure: giac's `envelope` answers the plot of the curve, printed
/// `Graphic_object`; ours is the curve `[x(t), y(t)]`, here `[k cos(t)³,
/// k sin(t)³]`, compared at `t = 1/3`.
#[test]
fn envelopes() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let t = Expr::from(1) / Expr::from(3);
    let mut wrong = Vec::new();
    for (input, k) in [
        ("envelope(y+x*tan(t)-2*sin(t),t)", 2),
        ("envelope(v+u*tan(t)-3*sin(t),[u,v,t])", 3),
        ("envelope(y+x*tan(t)=2*sin(t),[t])", 2),
    ] {
        let printed = giac.eval(input);
        if printed.as_deref().ok() != Some("Graphic_object") {
            wrong.push(format!("{input}: giac prints {printed:?}"));
        }
        let curve = |f: Symbol| Expr::from(k) * Expr::call(f, [t.clone()]).pow(Expr::from(3));
        match ours(input) {
            Ok(mine) => match mine.args() {
                [px, py]
                    if near(&ctx, &at(&ctx, px, "t", &t), &curve(Symbol::COS))
                        && near(&ctx, &at(&ctx, py, "t", &t), &curve(Symbol::SIN)) => {}
                _ => wrong.push(format!("{input}: ours {mine:?}")),
            },
            Err(e) => wrong.push(format!("{input}: {e}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
