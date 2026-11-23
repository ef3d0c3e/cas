//! Function study against giac (CAS-581, CAS-583): `domain` in mode 0
//! compared as a set of points on the line, at every end of either answer,
//! the midpoints between and a point beyond each extreme, and by its shape
//! (a list or one condition); `domain` in mode 1 as sets of values;
//! `convex` condition by condition, each side's difference under
//! `Together`; `tabvar` and `tabsign` row by row past the labels
//! (CAS-582); and the departures, each with giac's answer.

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

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// An answer as a list of conditions.
fn conditions(e: &Expr) -> Expr {
    if is_list(e) {
        e.clone()
    } else {
        Expr::call(Symbol::LIST, [e.clone()])
    }
}

/// The ends an answer names: the sides of its comparisons free of `x`.
fn ends(answer: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    for part in answer.args() {
        if !part.contains(x) {
            out.push(part.clone());
        } else if part != x {
            ends(part, x, out);
        }
    }
}

/// Whether `v` meets the answer: some condition holds at it.
fn meets(ctx: &Context, answer: &Expr, x: &Expr, v: &Expr) -> Option<bool> {
    let mut any = false;
    for c in answer.args() {
        let test = if c == x {
            Expr::from(Symbol::TRUE)
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

/// A number's approximation as an `f64`, through `Floor` of it times 2^32.
fn approximate(ctx: &Context, n: &Expr) -> Option<f64> {
    let scaled = Expr::call(Symbol::FLOOR, [n.clone() * Expr::from(1_i64 << 32)]);
    let k = ctx.eval(&scaled).ok()?.as_i64()?;
    #[allow(clippy::cast_precision_loss, reason = "an order, not a value")]
    Some(k as f64 / (1_i64 << 32) as f64)
}

/// Whether two answers, lists of conditions, are one set, sampled at the
/// ends of both, the midpoints and beyond; equal answers settle it (ends
/// with parameters have no value to sample).
fn same_set(ctx: &Context, a: &Expr, b: &Expr, x: &Expr) -> bool {
    if a == b {
        return true;
    }
    let mut found = Vec::new();
    ends(a, x, &mut found);
    ends(b, x, &mut found);
    let mut valued: Vec<(f64, Expr)> = Vec::new();
    for p in found {
        if let Some(v) = approximate(ctx, &p)
            && !valued.iter().any(|(w, _)| (w - v).abs() < 1e-9)
        {
            valued.push((v, p));
        }
    }
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

/// giac's `domain` answers and this CAS's: the same shape and set, in `x`
/// (or `t` where the input names it).
fn check_domains(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let x = Expr::symbol(if input.contains(",t)") { "t" } else { "x" });
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours))
                if is_list(&theirs) == is_list(&ours)
                    && same_set(&ctx, &conditions(&ours), &conditions(&theirs), &x) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Domains, giac's examples and more (CAS-581).
#[test]
fn domains() {
    check_domains(&[
        "domain(log(x+1))",
        "domain(log(x)+asin(x))",
        "domain(1/x)",
        "domain(sqrt(x^2-1))",
        "domain(ln(x)/(x-2))",
        "domain(tan(x))",
        "domain(tan(2*x))",
        "domain(tan(x)+x)",
        "domain(acosh(x))",
        "domain(x^(1/3))",
        "domain(x^(2/3))",
        "domain(x^(-2))",
        "domain(sin(x))",
        "domain(x^2+1)",
        "domain(1/(x-1)+1/(x+2))",
        "domain(1/((x-1)*(x+2)))",
        "domain(1/(x^2-1))",
        "domain(1/(x^3-x))",
        "domain(1/sqrt(x))",
        "domain(x^(-1/2))",
        "domain(sqrt(x)/(x-1))",
        "domain(sqrt(x)/(x+1))",
        "domain(asin(x)/x)",
        "domain(atanh(x))",
        "domain(acos(2*x))",
        "domain(log(x^2))",
        "domain(log(t-1),t)",
        "domain(sqrt(4-x^2)/(x-1))",
        "domain(x^x)",
        "domain(exp(1/x))",
        "domain(atan(1/x))",
        "domain(2^x)",
        "domain(log10(x))",
        "domain(abs(x)/x)",
        "domain(Ei(x))",
        "domain(sqrt(x-1)+sqrt(3-x))",
        "domain(sqrt(x^2+1))",
        "domain(log(-x^2-1))",
        "domain(x/(x^2+1))",
        "domain(ln(x)/(x^2-4))",
        "domain(log(x)/(x-1)+sqrt(x-3))",
        "domain(sqrt(a-x))",
        "domain(1/(x-a))",
    ]);
}

/// Singular values, `domain` in mode 1, as sets (CAS-581).
#[test]
fn singular_values() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "domain(sqrt(x),x,1)",
        "domain(log(x),x,1)",
        "domain(1/x+1/(x-1),x,1)",
        "domain(asin(x),x,1)",
        "domain(acosh(x),x,1)",
        "domain(atanh(x),x,1)",
        "domain(sqrt(x^2-1),x,1)",
        "domain(1/x,1)",
        "domain(sin(x),x,1)",
        "domain(tan(x),x,1)",
        "domain(ln(x)/(x-2),1)",
        "domain(sqrt(x)/(x+1),1)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours))
                if is_list(&ours)
                    && theirs.args().len() == ours.args().len()
                    && ours.args().iter().all(|v| theirs.args().contains(v)) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// An answer to `convex` with each comparison `l op r` as `op` of
/// `Together[l − r]`.
fn normal(ctx: &Context, e: &Expr) -> Option<Expr> {
    if is_list(e) {
        let parts: Option<Vec<Expr>> = e.args().iter().map(|c| normal(ctx, c)).collect();
        return Some(Expr::call(Symbol::LIST, parts?));
    }
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (
            Some(
                op @ (Symbol::GREATER_EQUAL | Symbol::LESS_EQUAL | Symbol::GREATER | Symbol::LESS),
            ),
            [l, r],
        ) => {
            let difference = ctx
                .eval(&Expr::call(Symbol::TOGETHER, [l.clone() - r.clone()]))
                .ok()?;
            Some(Expr::call(op, [difference, Expr::from(0)]))
        }
        _ => Some(e.clone()),
    }
}

/// Convexity, giac's examples and more (CAS-583).
#[test]
fn convexity() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "convex(3*exp(x)+5*x^4-ln(x),x)",
        "convex(x^2+y^2+3*z^2-x*y+2*x*z+y*z,[x,y,z])",
        "convex(x1^2+2*x2^2+a*x3^2-2*x1*x2+2*x1*x3-6*x2*x3,[x1,x2,x3],simplify)",
        "convex(x1^2+2*x2^2+a*x3^2-2*x1*x2+2*x1*x3-6*x2*x3,[x1,x2,x3])",
        "convex(x1^3+2*x1^2+2*x1*x2+x2^2/2-8*x1-2*x2-8,[x1,x2],simplify)",
        "convex(x1^3+2*x1^2+2*x1*x2+x2^2/2-8*x1-2*x2-8,[x1,x2])",
        "convex(x^2,x)",
        "convex(-x^2,x)",
        "convex(x^3,x)",
        "convex(x^4,x)",
        "convex(sin(x),x)",
        "convex(x^4-x^2,x)",
        "convex(a*x^2,x)",
        "convex(log(x),x)",
        "convex(1/x,x)",
        "convex(sqrt(x),x)",
        "convex(exp(x),x)",
        "convex(x*log(x),x)",
        "convex(x^2,[x])",
        "convex(x*y,[x,y])",
        "convex(x^2+y^2,[x,y])",
        "convex(-x^2-y^2,[x,y])",
        "convex(x^2-y^2,[x,y])",
        "convex(x^2+x*y+y^2,[x,y])",
        "convex(x^2+2*x*y+y^2,[x,y])",
        "convex(x^2*y^2,[x,y])",
        "convex(x^4+y^2,[x,y])",
        "convex(x^4+y^4,[x,y])",
        "convex(exp(x)+y^2,[x,y])",
        "convex(1/x+y^2,[x,y])",
        "convex(sqrt(x)+y^2,[x,y])",
        "convex(x^3+y^2,[x,y],simplify)",
        "convex(a*x^2+y^2+x*y,[x,y])",
        "convex(a*x^2+b*y^2,[x,y])",
        "convex(x+y,[x,y])",
        "convex(x^2+y^2+z^2,[x,y,z])",
        "convex(x*y*z,[x,y,z])",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours))
                if normal(&ctx, &ours).is_some()
                    && normal(&ctx, &ours) == normal(&ctx, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac: giac's answer as it prints it, and
/// this CAS's as decided.
fn check_departures() -> Vec<String> {
    let Some(giac) = giac() else {
        return Vec::new();
    };
    let cases = [
        // giac rewrites asinh(x) as ln(x+sqrt(x^2+1)) and its solve of
        // x+sqrt(x^2+1)>0 finds nothing; asinh has no condition here.
        ("domain(asinh(x))", "[]", "x"),
        // giac assumes 0 <= x <= 2*pi for a periodic f; here only for
        // tan(c*x) alone, and sin(x) = 0 has its principal solutions.
        (
            "domain(1/sin(x))",
            "((x<>0) and (x<>pi) and (x<>(2*pi)))",
            "x != 0 and x != pi",
        ),
        // The same assumption; sin(x) >= 0 is not solved here.
        (
            "domain(sqrt(sin(x)))",
            "((x>=0) and (x<=pi))",
            "domain(sqrt(sin(x)))",
        ),
        // Variables u(t) with u'(t), giac's makevars, stay here.
        (
            "convex(u(t)^2+u'(t)^2,u(t))",
            "true",
            "convex(u(t)^2+u'(t)^2,u(t))",
        ),
        // Labels with primes (CAS-582): our parser has none; the row of
        // f' is labelled "y'" = f', the row of f'' "(y')'".
        (
            "tabvar(x^2-2*x,x)[1][0]",
            "y'=(2*(x-1))",
            "\"y'\" = 2*(x-1)",
        ),
        ("tabvar(x^2-2*x,x)[3][0]", "(y')'", "\"(y')'\""),
        // giac's parametric table (step_param_) stays here.
        (
            "tabvar([sin(2t),cos(3t)])",
            r#"[[t,0," ",pi/4," ",pi/3," ",pi/2],[x=(sin(2*t)),0,"↑",1,"↓",sqrt(3)/2,"↓",0],[x'=(2*cos(2*t)),2,"+",0,"-",-1,"-",-2],[y=(cos(3*t)),1,"↓",-sqrt(2)/2,"↓",-1,"↑",0],[y'=(-3*sin(3*t)),0,"-",-3*sqrt(2)/2,"-",0,"+",3],[x'*(y')'-y'*(x')',-18,"concav",-6*sqrt(2),"concav",-9,"concav",0]]"#,
            "tabvar([sin(2t),cos(3t)])",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!(
                "{input}: ours {:?} is not {decided} ({:?})",
                ours(input),
                ours(decided)
            ));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    wrong
}

#[test]
fn departures() {
    let wrong = check_departures();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `e` with each `Exp[u]` as `Power[E, u]`.
fn powers(e: &Expr) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    let args: Vec<Expr> = e.args().iter().map(powers).collect();
    match (e.head(), &args[..]) {
        (Some(h), [u]) if *h == Expr::from(Symbol::EXP) => {
            Expr::call(Symbol::POWER, [Expr::from(Symbol::E), u.clone()])
        }
        (Some(h), _) => Expr::apply(&powers(h), &args),
        (None, _) => e.clone(),
    }
}

/// Whether two cells of a table are one: the same string, lists cell by
/// cell, values whose difference `Simplify` makes 0.
fn same_cell(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    if is_list(a) && is_list(b) {
        return a.args().len() == b.args().len()
            && a.args()
                .iter()
                .zip(b.args())
                .all(|(p, q)| same_cell(ctx, p, q));
    }
    if a.as_str().is_some() || b.as_str().is_some() {
        return false;
    }
    // exp(-1) and giac's 1/exp(1) as powers of e.
    let difference = Expr::call(Symbol::SIMPLIFY, [powers(a) - powers(b)]);
    ctx.eval(&difference)
        .is_ok_and(|d| d.to_rational().is_some_and(|q| q.is_zero()))
}

/// giac's `tabvar` and `tabsign` tables and this CAS's: the number of
/// rows, each row past its label (`tail(t[k])`, giac's `y'=(…)` does not
/// read back), and the labels that do (`x`, `y = f`, `tabsign`'s `f`);
/// `undef` as one cell.
fn check_tables(inputs: &[&str]) -> Vec<String> {
    let Some(giac) = giac() else {
        return Vec::new();
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let Ok(table) = ours(input) else {
            wrong.push(format!("{input}: ours unreadable"));
            continue;
        };
        let mut parts = Vec::new();
        if is_list(&table) {
            parts.push(format!("size({input})"));
            for k in 0..table.args().len() {
                parts.push(format!("tail(({input})[{k}])"));
            }
            parts.push(format!("({input})[0][0]"));
            parts.push(format!(
                "({input})[{}][0]",
                if input.starts_with("tabsign") { 1 } else { 2 }
            ));
        } else {
            parts.push((*input).to_owned());
        }
        for part in parts {
            match (giac.value(&part), ours(&part)) {
                (Ok(theirs), Ok(ours)) if same_cell(&ctx, &theirs, &ours) => {}
                (theirs, ours) => wrong.push(format!("{part}: ours {ours:?}, giac's {theirs:?}")),
            }
        }
    }
    wrong
}

/// Variation and sign tables (CAS-582): giac's examples and more.
#[test]
fn tables() {
    let wrong = check_tables(&TABLES);
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

const TABLES: [&str; 31] = [
    "tabsign(x^2-x,x)",
    "tabsign((x-1)/(x+2),x)",
    "tabsign(x^2,x,-3,5)",
    "tabsign(x^3-x,x)",
    "tabsign(x^2-2,x)",
    "tabsign(ln(x),x)",
    "tabsign(abs(x)-1,x)",
    "tabvar(x^2-2*x,x)",
    "tabvar(x^2,x,-3,5)",
    "tabvar(sin(x),x)",
    "tabvar(1/x,x)",
    "tabvar(x^3-3*x,x)",
    "tabvar(x^2+x+1,x)",
    "tabvar(x/(x^2+1),x)",
    "tabvar(sqrt(x),x)",
    "tabvar(sqrt(4-x^2),x)",
    "tabvar(x*ln(x),x)",
    "tabvar(exp(-x^2),x)",
    "tabvar(abs(x-1),x)",
    "tabvar(sin(x),x,0,4*pi)",
    "tabvar(cos(x),x=0..pi)",
    "tabvar((x^2-1)/(x-2),x)",
    "tabvar(x^2/(x^2-1),x)",
    "tabvar(x*exp(x),x)",
    "tabvar(ln(x)/x,x)",
    "tabvar(x^4-2*x^2,x)",
    "tabvar(x+sin(x),x)",
    "tabvar(x^3,x,diff=1)",
    "tabvar(x^3-3*x,x=-2..3)",
    "tabvar(exp(x),x)",
    "tabvar(1/ln(x^2-1),x,diff=1)",
];
