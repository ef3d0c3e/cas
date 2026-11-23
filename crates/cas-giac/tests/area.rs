//! `area` against giac (CAS-588, CAS-589): giac's `area` with a number of
//! subdivisions and a rule answers a double, which this CAS's ball must
//! agree with to `10⁻¹⁰`, relatively; without them it answers the exact
//! integral.

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

/// Whether `a` is within `10⁻¹⁰ max(1, |b|)` of `b`.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(Symbol::N, [a.clone() - b.clone()])],
    );
    let scale = Expr::call(
        Symbol::MAX,
        [
            Expr::from(1),
            Expr::call(Symbol::ABS, [Expr::call(Symbol::N, [b.clone()])]),
        ],
    );
    let claim = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-10) * scale]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether `a − b` is 0 once expanded and put together.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// `e` with `a` replaced by `3`, evaluated.
fn at_three(ctx: &Context, e: &Expr) -> Option<Expr> {
    let rule = Expr::call(Symbol::RULE, [Expr::symbol("a"), Expr::from(3)]);
    ctx.eval(&Expr::call(Symbol::REPLACE_ALL, [e.clone(), rule]))
        .ok()
}

/// Every rule by name, giac's French names, Romberg's two forms by
/// depth (and depth 0, one panel), Simpson's odd counts, reversed bounds.
#[test]
fn rules() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "area(x^2,x=0..1,5,trapezoid)",
        "area(x^2,x=0..1,5,middle_point)",
        "area(x^2,x=0..1,5,simpson)",
        "area(x^2,x=0..1,5,rombergt)",
        "area(x^2,x=0..1,5,rombergm)",
        "area(x^2,x=0..1,5,gauss15)",
        "area(exp(x),x=0..1,4,trapezoid)",
        "area(exp(x),x=0..1,4,trapeze)",
        "area(exp(x),x=0..1,4,middle_point)",
        "area(exp(x),x=0..1,4,point_milieu)",
        "area(exp(x),x=0..1,4,left_rectangle)",
        "area(exp(x),x=0..1,4,rectangle_gauche)",
        "area(exp(x),x=0..1,4,right_rectangle)",
        "area(exp(x),x=0..1,4,rectangle_droit)",
        "area(exp(x),x=0..1,4,simpson)",
        "area(exp(x),x=0..1,4,gauss15)",
        "area(exp(x),x=0..1,1,rombergt)",
        "area(exp(x),x=0..1,2,rombergt)",
        "area(exp(x),x=0..1,3,rombergt)",
        "area(exp(x),x=0..1,1,rombergm)",
        "area(exp(x),x=0..1,2,rombergm)",
        "area(exp(x),x=0..1,3,rombergm)",
        "area(exp(x),x=0..1,3,simpson)",
        "area(sin(x),x=0..pi,10,simpson)",
        "area(exp(-x^2),x=0..2,6,trapezoid)",
        "area(sqrt(x),x=0..1,5,rombergm)",
        "area(x^2,x=1..0,5,trapezoid)",
        "area(1/x,x=1..100,4,gauss15)",
        "area(x^2,x=0..1,0,rombergt)",
        "area(x^2,x=0..1,0,rombergm)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Without a count the exact integral; no panels is `undef`, ours
/// `Indeterminate`.
#[test]
fn exact() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in ["area(exp(x),x=0..1)", "area(x^2,x=0..1)"] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    let input = "area(x^2,x=0..1,0,trapezoid)";
    if let Err(e) = giac.agrees(input) {
        wrong.push(format!("{e:?}"));
    }
    if ours(input).ok() != Some(Expr::from(Symbol::INDETERMINATE)) {
        wrong.push(format!("{input}: ours {:?}", ours(input)));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as libgiac
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    let printed = [
        // A count without a rule: giac falls through to its polygon
        // branch and answers an equation; ours stays unevaluated.
        (
            "area(exp(x),x=0..1,4)",
            "0=(1/2*exp(x)*im((-exp(x)) .. (-exp(x)+1))-2*im((-exp(x)) .. (-exp(x)+1)))",
        ),
        // gauss15 below 4 panels: giac answers a pair, its estimate
        // second; ours a ball around ln(100).
        (
            "area(1/x,x=1..100,1,gauss15)",
            "[17.9719306931,4.59372597601]",
        ),
        (
            "area(1/x,x=1..100,2,gauss15)",
            "[9.87315051625,4.60514490035]",
        ),
        (
            "area(1/x,x=1..100,3,gauss15)",
            "[6.59202964673,4.60517018599]",
        ),
        // A symbolic integrand: giac leaves the sum unsimplified, ours
        // the exact rule value.
        (
            "area(a*x,x=0..1,5,trapezoid)",
            "1/5*a*0.2+1/5*a*0.4+1/5*a*0.6+1/5*a*0.8+1/5*a/2",
        ),
        // A negative count: giac steps by (1-0)/-1 over [0,-1]; ours
        // stays unevaluated.
        ("area(x^2,x=0..1,-1,trapezoid)", "-0.5"),
        // Romberg on a symbolic integrand: giac answers garbage, ours a/2.
        ("area(a*x,x=0..1,2,rombergt)", "6.10593455919e+306"),
    ];
    for (input, answer) in printed {
        match giac.eval(input) {
            Ok(text) if text == answer => {}
            other => wrong.push(format!("{input}: giac gives {other:?}, recorded {answer}")),
        }
    }
    let area = ours("area").ok();
    let unevaluated = |input: &str| {
        let value = ours(input);
        if value.as_ref().ok().and_then(Expr::head) == area.as_ref() {
            None
        } else {
            Some(format!("{input}: ours {value:?} is not unevaluated"))
        }
    };
    wrong.extend(unevaluated("area(exp(x),x=0..1,4)"));
    wrong.extend(unevaluated("area(x^2,x=0..1,-1,trapezoid)"));
    let ln100 = Expr::call(Symbol::LOG, [Expr::from(100)]);
    for n in 1..=3 {
        let input = format!("area(1/x,x=1..100,{n},gauss15)");
        match ours(&input) {
            Ok(value) if near(&ctx, &value, &ln100) => {}
            other => wrong.push(format!("{input}: ours {other:?}, not ln(100)")),
        }
    }
    let half = Expr::symbol("a") / Expr::from(2);
    for input in [
        "area(a*x,x=0..1,5,trapezoid)",
        "area(a*x,x=0..1,2,rombergt)",
    ] {
        match ours(input) {
            Ok(value) if same(&ctx, &value, &half) => {}
            other => wrong.push(format!("{input}: ours {other:?}, not a/2")),
        }
    }
    let input = "area(a*x,x=0..1,5,trapezoid)";
    let three = |e: Result<Expr, _>| e.ok().and_then(|e| at_three(&ctx, &e));
    let target = Expr::from(3) / Expr::from(2);
    for (side, value) in [
        ("giac", three(giac.value(input))),
        ("ours", three(ours(input))),
    ] {
        if !value.is_some_and(|v| near(&ctx, &v, &target)) {
            wrong.push(format!("{input} at a = 3: {side}'s is not 1.5"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
