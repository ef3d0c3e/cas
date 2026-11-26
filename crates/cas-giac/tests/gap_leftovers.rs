//! giac's expression, calculus, integer and box leftovers of the gap ledger
//! (CAS-640 to CAS-643) against giac 2.0: `root`, `atan2`, `has`, `fmod`,
//! `modf`, `convert(e, expln)`, `additionally`, `avgRC`, `addtable`,
//! `multinomial`, `order`, `maple_ifactors` and `box_constraints`. giac's
//! answers read back are this CAS's; the departures carry giac's answer as
//! it prints it and this CAS's as decided. `addtable` and `additionally`
//! act on the session, so their rows run statements in one session on both
//! sides, giac restarted around them.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use cas_lang::Session;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac's answers to `inputs`, read back, are this CAS's: the same, equal
/// by `==`, or numbers within giac's 12 printed digits (C0006).
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
/// each other once approximated.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| near(ctx, a, b));
    }
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
fn expressions() {
    check(&[
        "root(3,2)",
        "root(1/3,2)",
        "root(3,1.2)",
        "root(3.2,1.2)",
        "root(2,-4)",
        "root(3,8)",
        "root(4,16)",
        "root(2,9)",
        "root(2,x)",
        "root(n,x)",
        "root(1/2,x)",
        "root(2,-1.0)",
        "root(0,2)",
        "atan2(1,1)",
        "atan2(0,-1)",
        "atan2(-1,0)",
        "atan2(0,0)",
        "atan2(0,1)",
        "atan2(-1,-1)",
        "atan2(sqrt(3),1)",
        "atan2(1,-sqrt(3))",
        "atan2(-1/2,-sqrt(3)/2)",
        "atan2(1.0,2.0)",
        "atan2(y,x)",
        "has(x+y,x)",
        "has(x+y,n)",
        "has(sin(x^2),x)",
        "has(x,x)",
        "has(3,x)",
        "has([x,1],x)",
        "has(f(x),f)",
        "has(x+y,x+y)",
        "has(sin(x),sin)",
        "has(2*x,2)",
        "fmod(10,3)",
        "fmod(10.0,pi)",
        "fmod(-10.0,3.0)",
        "fmod(7.5,2)",
        "fmod(-7.5,2)",
        "fmod(7.5,-2)",
        "fmod(x,2)",
        "fmod(1/3,1/7)",
        "fmod(7,0)",
        "fmod(7.5,0)",
        "[modf(2.3)]",
        "[modf(-2.3)]",
        "[modf(5)]",
        "[modf(7/2)]",
        "[modf(-7/2)]",
        "evalf(subst(convert(sin(x),expln),x=1/3))",
        "evalf(subst(convert(cos(x)+exp(x),expln),x=1/3))",
    ]);
}

#[test]
fn calculus_integers_and_boxes() {
    check(&[
        "avgRC(f(x),x,h)",
        "avgRC(x^2,x,0)",
        "subst(avgRC(x^2,x),x=2)",
        "subst(avgRC(x^2,x,1/10),x=2)",
        "multinomial(10,[0.5,0.5],[3,7])",
        "multinomial(10,[0.2,0.3,0.5],[1,3,6])",
        "multinomial(10,[1/2,1/2],[3,7])",
        "multinomial(4,[1/2,1/3,1/6],[1,2,1])",
        "multinomial(2,[p,1-p],[1,1])",
        "multinomial(2,[1/3,1/3],[1,1])",
        "multinomial(3,[1/2,1/2],[1,2])",
        "multinomial(2.0,[1/2,1/2],[1,1])",
        "order(3 % 7)",
        "order(2 % 7)",
        "order(1 % 7)",
        "order(6 % 7)",
        "order(3 % 1000003)",
        "maple_ifactors(120)",
        "maple_ifactors(-120)",
        "maple_ifactors(-12)",
        "maple_ifactors(1000000007)",
        "maple_ifactors(0)",
        "maple_ifactors(1)",
        "maple_ifactors(-1)",
        "[box_constraints([x,y],[[-2,2],[-5,5]])]",
        "[box_constraints([x],[[0,1]])]",
        "[box_constraints([x,y,z],[[a,b],[1/2,3],[-1,1]])]",
    ]);
}

/// Where giac errs, the call stays with a warning; giac's `undef` for
/// `order` of an integer, its decimals for an exact `avgRC` step, and its
/// form of `atan2` on the axis `x = 0`, differ.
#[test]
fn departures() {
    departed(&[
        ("order(2 % 8)", None, "order(2 % 8)"),
        ("order(3)", Some("undef"), "order(3)"),
        (
            "multinomial(3,[1/2,1/2],[1,1])",
            None,
            "multinomial(3,[1/2,1/2],[1,1])",
        ),
        (
            "multinomial(2,[1/2],[1,1])",
            None,
            "multinomial(2,[1/2],[1,1])",
        ),
        ("maple_ifactors(2.5)", None, "maple_ifactors(2.5)"),
        ("maple_ifactors(x)", None, "maple_ifactors(x)"),
        ("modf(x)", None, "modf(x)"),
        ("modf(2+i)", None, "modf(2+i)"),
        (
            "box_constraints([x,y],[[0,1]])",
            None,
            "box_constraints([x,y],[[0,1]])",
        ),
        (
            "box_constraints([x],[[0]])",
            None,
            "box_constraints([x],[[0]])",
        ),
        (
            "avgRC(x^2,x,1/10)",
            Some("10.0*((x+0.1)^2-x^2)"),
            "((x+1/10)^2-x^2)/(1/10)",
        ),
        (
            "atan2(x,0)",
            Some("pi/2+(1-sign(x))*pi/2-2*pi*floor((pi/2+(1-sign(x))*pi/2)/(2*pi)+1/2)"),
            "sign(x)*pi/2",
        ),
    ]);
}

/// The values giac prints for a sequence of statements, in order.
fn values(answer: &str) -> Vec<&str> {
    let mut depth = 0i32;
    let mut start = 0;
    let mut values = Vec::new();
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => {
                values.push(&answer[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    values.push(&answer[start..]);
    values
}

/// This CAS's values of the statements of `input` that give one, giac's `;`
/// read as separate statements, run in one session.
fn session_values(input: &str) -> Vec<Expr> {
    let mut session = Session::new();
    session
        .run(&input.replace(';', "\n"))
        .iter()
        .filter_map(|output| {
            let Some(result) = &output.result else {
                assert!(output.messages.is_empty(), "{input}: {:?}", output.messages);
                return None;
            };
            Some(ours(result.text()).unwrap_or_else(|e| panic!("{input}: {e:?}")))
        })
        .collect()
}

/// Whether two values are one: the same, or equal by `==`.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    a == b
        || ctx
            .eval(&Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]))
            .ok()
            .and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
}

/// giac's `addtable(fourier, …)` and `additionally`, in a session: giac's
/// values, one by one; `true` where giac's `assume` answers its variable
/// and ours nothing, or giac's `addtable` answers 0 for a pair it holds
/// already (its table outlives `restart`): the last values only are
/// compared, which depend on the added fact or pair. The assumptions are on `w`, which no other row of
/// this binary uses: giac's session is shared by the threads.
#[test]
fn sessions() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, last_only) in [
        ("addtable(fourier,y(x),Y(s),x,s);fourier(y(x),x,s)", false),
        ("addtable(fourier,y(x),Y(s),x,s);fourier(y(t),t,u)", true),
        ("addtable(fourier,y(x),Y(s),x,s);fourier(2*y(x),x,s)", true),
        ("assume(w>-5);additionally(w>1);sign(w-1)", true),
        ("assume(w>-5);assume(w>1,additionally);sign(w-1)", true),
        ("assume(w>-5);additionally(w>1);abs(w-1)", true),
        ("assume(w>0);additionally(w<1);abs(w-1)", true),
    ] {
        giac.restart();
        let answer = giac
            .eval(input)
            .unwrap_or_else(|e| panic!("{input}: {e:?}"));
        let mut texts = values(&answer);
        let mut mine = session_values(input);
        if last_only {
            texts = texts.split_off(texts.len().saturating_sub(1));
            mine = mine.split_off(mine.len().saturating_sub(1));
        }
        let theirs: Vec<Expr> = texts
            .into_iter()
            .map(|v| ours(v).unwrap_or_else(|e| panic!("{input}: {e:?}")))
            .collect();
        if theirs.len() != mine.len() || theirs.iter().zip(&mine).any(|(a, b)| !same(&ctx, a, b)) {
            wrong.push(format!("{input}: giac's {answer}, ours {mine:?}"));
        }
    }
    giac.restart();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `laplace` does not consult the pair `addtable` stores, and
/// answers its integral (`1,integrate(g(t)*exp(-s*t),t,0,+infinity)`);
/// ours answers the pair, `1, G(s)`.
#[test]
fn laplace_table_departs() {
    let Some(giac) = giac() else {
        return;
    };
    let input = "addtable(laplace,g(t),G(s),t,s);laplace(g(t),t,s)";
    giac.restart();
    let answer = giac.eval(input).unwrap_or_default();
    giac.restart();
    let mut wrong = Vec::new();
    if answer.replace(' ', "") != "1,integrate(g(t)*exp(-s*t),t,0,+infinity)" {
        wrong.push(format!("{input}: giac prints {answer}"));
    }
    let mine = session_values(input);
    let decided = [ours("1").unwrap(), ours("G(s)").unwrap()];
    if mine != decided {
        wrong.push(format!("{input}: ours {mine:?}"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
