//! giac's list, type, format, date and polynomial leftovers of the gap
//! ledger (CAS-644 to CAS-647) against giac 2.0: `back`, `range`, `shift`,
//! `suppress`, `at`, `rdiv`, `getType`, `subs`, `format`, `dayofweek`,
//! `convert(e, list)`, `roots`, `truncate`, `divis`, `choice`, and `fitdistr`'s
//! refusal of an exponential fit with a 0 datum. giac's answers read back
//! are this CAS's: the same, equal by `==`, or numbers within giac's 12
//! printed digits (C0100). `roots` lists the roots of giac's solver in its
//! order, ours those of `CSolve`: compared as multisets of `[root,
//! multiplicity]`, roots at a parameter's sample value (a departure in
//! order). The departures carry giac's answer as it prints it (`None` for
//! an error, where ours stays with giac's message as a warning) and this
//! CAS's as decided.

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

/// giac's answers to `inputs`, read back, are this CAS's: the same, equal
/// by `==`, or numbers within giac's 12 printed digits.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| {
            let error = giac.agrees(input).err()?;
            let close = match (giac.value(input), ours(input)) {
                (Ok(theirs), Ok(mine)) => near(&ctx, &theirs, &mine),
                _ => false,
            };
            (!close).then(|| error.to_string())
        })
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whether `a` and `b`, lists entry by entry, are numbers within 10⁻⁹ of
/// each other once approximated, the parameters `a` and `y` at 7/3.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| near(ctx, a, b));
    }
    let sample = Expr::from(7) / Expr::from(3);
    let gap = (a.clone() - b.clone()).substitute(&[
        (Expr::symbol("a"), sample.clone()),
        (Expr::symbol("y"), sample),
    ]);
    let gap = Expr::call(Symbol::N, [Expr::call(Symbol::ABS, [gap])]);
    let bound = Expr::from(1) / Expr::from(1_000_000_000_i64);
    ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

/// giac's `roots` and ours: the same `[root, multiplicity]` pairs in any
/// order, roots within 10⁻⁹.
fn same_roots(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let (theirs, mine) = match (giac.value(input), ours(input)) {
            (Ok(t), Ok(m)) => (t, m),
            (t, m) => {
                wrong.push(format!("{input}: giac {t:?}, ours {m:?}"));
                continue;
            }
        };
        let mut left: Vec<&Expr> = mine.args().iter().collect();
        let matched = theirs.args().len() == left.len()
            && theirs.args().iter().all(|pair| {
                let found = left.iter().position(|p| match (pair.args(), p.args()) {
                    ([r, k], [s, l]) => k == l && near(&ctx, r, s),
                    _ => false,
                });
                found.map(|i| left.remove(i)).is_some()
            });
        if !matched {
            wrong.push(format!("{input}: giac {theirs:?}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures: giac's printed answer (`None` for an error), and this CAS's
/// value, the same as `decided`'s; giac's does not agree.
fn departed(cases: &[(&str, Option<&str>, &str)]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for &(input, answer, decided) in cases {
        let printed = giac.eval(input);
        let recorded = match (answer, &printed) {
            (None, Ok(p)) => p.contains("rror"),
            (None, Err(_)) => true,
            (Some(a), Ok(p)) => p.replace(' ', "") == a,
            (Some(_), Err(_)) => false,
        };
        if !recorded {
            wrong.push(format!(
                "{input}: giac prints {printed:?}, recorded {answer:?}"
            ));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn lists() {
    check(&[
        "back(1,2,3)",
        "back([1,2,3])",
        "back(\"bonjour\")",
        "back([])",
        "back([[1,2],[3,4]])",
        "back(x)",
        "rdiv(3,5)",
        "rdiv(3.2,5.4)",
        "rdiv(x,y)",
        "rdiv(6,4)",
        "range(10)",
        "range(1,11)",
        "range(1,11,2)",
        "range(1,11,3)",
        "range(11,1,-1)",
        "range(0)",
        "range(1.5,4)",
        "range(1,2,0.25)",
        "range(1,3,1/2)",
        "range(3,1,-1/2)",
        "range(0,1,1/3)",
        "range(0,1,1/7)",
        "range(0,1,0.1)",
        "range(1,1)",
        "shift([0,1,2,3],2)",
        "shift([0,1,2,3])",
        "shift([0,1,2,3,4])",
        "shift([0,1,2,3],-2)",
        "shift([0,1,2,3],0)",
        "shift([0,1,2,3],5)",
        "shift([1,2],-5)",
        "shift([1,2,3],1.5)",
        "shift(5,1)",
        "shift(5)",
        "shift(-5,1)",
        "shift(\"abcd\",1)",
        "shift(\"abcd\")",
        "at([10,11,12],1)",
        "at([[1,2],[3,4]],[1,0])",
        "at(\"abc\",1)",
        "at([10,11,12],-1)",
        "at([[1,2],[3,4]],1)",
        "at([[1,2],[3,4]],[1])",
        "suppress([0,1,2,3],2)",
        "suppress(\"0123\",2)",
        "suppress([0,1,2,3],4)",
        "suppress([0,1,2,3],-1)",
        "suppress([0,1,2,3],[0,2])",
        "suppress(\"0123\",5)",
        "suppress([0,1,2,3],1..2)",
        "suppress([0,1,2,3],1,2)",
        "suppress(\"0123\",\"13\")",
        "subs(1/(4+x^2),x=2)",
        "subs(x^2,x=3)",
        "subs(x+y,[x=1,y=2])",
        "convert([1,2],list)",
        "convert(x^2+1,list)",
        "convert(matrix(2,2,[1,2,3,4]),list)",
        "convert(x,list)",
        "convert(\"ab\",list)",
    ]);
}

#[test]
fn types_formats_and_dates() {
    check(&[
        "getType(\"aze\")",
        "getType(2+i)",
        "getType(x->cos(2*x))",
        "getType(1.414)",
        "getType(x)",
        "getType([1,2])",
        "getType([[1,2],[3,4]])",
        "getType(x+1)",
        "getType(sin(x))",
        "getType(3)",
        "getType(1/2)",
        "getType(-3)",
        "getType(sqrt(2))",
        "getType(pi)",
        "getType(x=1)",
        "getType([])",
        "getType(sin)",
        "getType(e)",
        "getType(infinity)",
        "getType(undef)",
        "getType(i)",
        "getType(2.5)",
        "getType(true)",
        "getType(false)",
        "getType(-x)",
        "getType([[1],[2,3]])",
        "getType(1..2)",
        "getType(euler_gamma)",
        "format(9.3456,\"f3\")",
        "format(9.3456,\"s3\")",
        "format(9.3456,\"e3\")",
        "format(1234.5678,\"f2\")",
        "format(1234.5678,\"s2\")",
        "format(1234.5678,\"e2\")",
        "format(0.000123,\"e3\")",
        "format(3,\"f2\")",
        "format(pi,\"f5\")",
        "format(-2.5,\"f0\")",
        "format(1/3,\"f4\")",
        "format(123456789,\"s3\")",
        "format(2.5,\"f0\")",
        "format(0.5,\"f0\")",
        "format(1.0E-5,\"s3\")",
        "format(123456.0,\"s6\")",
        "format(100.0,\"s3\")",
        "format(2,\"e2\")",
        "format(-1234.5678,\"s2\")",
        "format(0.0,\"e2\")",
        "format(5,\"x3\")",
        "format(5,\"f\")",
        "format(7,\"s\")",
        "dayofweek(21,4,2014)",
        "dayofweek(15,10,1582)",
        "dayofweek(26,11,2026)",
        "dayofweek(1,1,2000)",
        "dayofweek(29,2,2001)",
        "dayofweek(4,10,1582)",
        "dayofweek(1,1,1)",
        "dayofweek(1,1,-1)",
    ]);
}

#[test]
fn polynomials() {
    check(&[
        "truncate((x^2+x)^2,3)",
        "truncate(1+x+x^2,1)",
        "truncate((1+x)^5,2)",
        "truncate((x+y)^3,2)",
        "truncate(sin(x),3)",
        "truncate(x*sin(x)+x^3,2)",
        "truncate(x^3/(1+y),2)",
        "truncate(x^3+y^4,[x,y],3)",
        "truncate(x^4+1)",
        "truncate(x^7+x^6+1)",
        "truncate(x^2+1,-1)",
        "truncate((x+1)^3,y,1)",
        "truncate((x+y)^2,[y],1)",
        "divis(x^2-1)",
        "divis(x^4-1)",
        "divis(t^4-1)",
        "divis([x^4-1,x^2-1])",
        "divis([x^4-1,t^2-1])",
        "divis(x^2)",
        "divis(2*x^2-2)",
        "divis(6)",
        "divis(x^2+1)",
        "divis((x-1)^2*(x+2))",
        "divis(0)",
        "divis(1/(x+1))",
        "divis(x^2-y^2)",
        "divis(x)",
        "divis(x^3)",
        "divis([6,x^2])",
    ]);
    same_roots(&[
        "roots(t^3-1,t)",
        "roots(x^5-2*x^4+x^3)",
        "roots(x^2+1)",
        "roots(x^2-2)",
        "roots(x^3-2)",
        "roots(x^2*y-1,x)",
        "roots(3)",
        "roots(x^4-1)",
        "roots(1/(x-1))",
        "roots((x-1)^2/(x+2))",
        "roots(x^2-3*x+2)",
        "roots(y^2-1,y)",
        "roots(x^3-x,x)",
        "roots(x^2-a)",
    ]);
}

/// Where giac errs, the call stays with a warning (`subs(x=2, e)`, a
/// `subst` with its arguments switched, stays without one); giac's
/// `format` of a non-real prints it, ours stays; `divis` lists the factor
/// −1 first, as our `factors` does.
#[test]
fn departures() {
    departed(&[
        ("range(1,2,0)", None, "range(1,2,0)"),
        ("range(x)", None, "range(x)"),
        ("range(5,1)", None, "range(5,1)"),
        ("range(-3)", None, "range(-3)"),
        ("range(2.5)", None, "range(2.5)"),
        ("subs(x=2,1/(4+x^2))", None, "subs(x=2,1/(4+x^2))"),
        ("at([10,11,12],3)", None, "at([10,11,12],3)"),
        ("at([[1,2],[3,4]],[1,5])", None, "at([[1,2],[3,4]],[1,5])"),
        ("suppress([0,1,2],x)", None, "suppress([0,1,2],x)"),
        ("truncate(1/(1+x),2)", None, "truncate(1/(1+x),2)"),
        ("truncate(x^2+1,1.5)", None, "truncate(x^2+1,1.5)"),
        ("dayofweek(1,13,2000)", None, "dayofweek(1,13,2000)"),
        ("dayofweek(32,1,2000)", None, "dayofweek(32,1,2000)"),
        ("dayofweek(1.5,1,2000)", None, "dayofweek(1.5,1,2000)"),
        ("dayofweek(1,1)", None, "dayofweek(1,1)"),
        (
            "fitdistr([0,1,2,3,4],exponential)",
            None,
            "fitdistr([0,1,2,3,4],exponential)",
        ),
        ("format(x,\"f3\")", Some("\"x\""), "format(x,\"f3\")"),
        (
            "format(2+i,\"f2\")",
            Some("\"2.00+i\""),
            "format(2+i,\"f2\")",
        ),
        (
            "divis(-x^2+1)",
            Some("[1,x-1,x+1,(x-1)*(x+1),-1,-x+1,-x-1,-(x-1)*(x+1)]"),
            "[1,-1,x-1,-(x-1),x+1,-(x+1),(x-1)*(x+1),-(x-1)*(x+1)]",
        ),
    ]);
}

/// Whether `e` is an item of the list `expected` (C0014: the seeded
/// generators differ, so draws are checked by property).
fn one_of(e: &Expr, expected: &str) -> bool {
    ours(expected).is_ok_and(|l| l.args().contains(e))
}

/// `choice(l)` is `rand(l)` (CAS-646): an element of `l` from giac and
/// ours, three after `srand(3)`; ours the same draw as `rand(l)` after the
/// same seed; anything else giac's "Bad Argument Value", ours stays.
#[test]
fn choice() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, items) in [
        ("choice([a,b,c])", "[a,b,c]"),
        (
            "choice([\"r\",\"r\",\"r\",\"b\",\"n\"])",
            "[\"r\",\"b\",\"n\"]",
        ),
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if one_of(&theirs, items) && one_of(&mine, items) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    let three = "[choice([a,b,c]),choice([a,b,c]),choice([a,b,c])]";
    let printed = giac.eval(&format!("srand(3);{three}"));
    let theirs = printed
        .as_deref()
        .ok()
        .and_then(|p| p.strip_prefix("3,"))
        .map(ours);
    let mine = ours(&format!("[srand(3),{three}]"));
    let mine = mine.as_ref().ok().and_then(|l| l.args().get(1));
    let drawn = |l: &Expr| l.args().len() == 3 && l.args().iter().all(|e| one_of(e, "[a,b,c]"));
    match (&theirs, mine) {
        (Some(Ok(t)), Some(m)) if drawn(t) && drawn(m) => {}
        _ => wrong.push(format!(
            "srand(3);{three}: ours {mine:?}, giac prints {printed:?}"
        )),
    }
    for (k, l) in [
        (5, "[1,2,3,4,5,6]"),
        (3, "[a,b,c]"),
        (42, "[\"r\",\"b\",\"n\"]"),
    ] {
        let chosen = ours(&format!("[srand({k}),choice({l})]"));
        let drawn = ours(&format!("[srand({k}),rand({l})]"));
        if chosen.is_err() || chosen != drawn {
            wrong.push(format!(
                "srand({k});choice({l}): ours {chosen:?}, rand {drawn:?}"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    departed(&[
        ("choice([])", None, "choice([])"),
        ("choice(5)", None, "choice(5)"),
        ("choice(1,2,3)", None, "choice(1,2,3)"),
    ]);
}
