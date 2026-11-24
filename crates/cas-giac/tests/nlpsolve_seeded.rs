//! `nlpsolve` against giac (CAS-604, CAS-605), the rows whose giac answer draws
//! random numbers (C0014): giac seeded with `srand(42)` in the same call,
//! in a test binary of their own. Ours is checked against each program's
//! closed form first, giac's within its own error (C0006).

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

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists and equations part by part, the
/// same, or numbers within `1/scale`.
fn same(ctx: &Context, a: &Expr, b: &Expr, scale: i64) -> bool {
    if a == b {
        return true;
    }
    if a.head().is_some()
        && a.head() == b.head()
        && a.head()
            .and_then(Expr::as_symbol)
            .is_some_and(|h| matches!(h, Symbol::LIST | Symbol::EQUATION))
    {
        return a.args().len() == b.args().len()
            && a.args()
                .iter()
                .zip(b.args())
                .all(|(a, b)| same(ctx, a, b, scale));
    }
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [a.clone() - b.clone()])) else {
        return false;
    };
    let parts = match difference.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![difference],
    };
    let tolerance = Expr::from(1) / Expr::from(scale);
    parts.into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// giac's answer to `input` with its random generator seeded (42) first,
/// in the same call: its text, the last value of the printed sequence, or
/// the whole of an error.
fn seeded_text(giac: &Giac, input: &str) -> Option<String> {
    let answer = giac.eval(&format!("srand(42);{input}")).ok()?;
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => start = i + 1,
            _ => {}
        }
    }
    Some(answer[start..].replace("list[", "["))
}

/// `(input, closed form, giac's scale)`: giac's answer within `1/scale` of
/// the closed form, its observed error (T4 probes under seeds 42 and 1).
type Row<'a> = (&'a str, &'a str, i64);

const BOWL: &str = "[0,[x=1,y=2]]";

/// giac's options of `nlpsolve` on programs it solves (CAS-604): ours the
/// closed form to 10⁻¹⁰, certified; giac's within its error.
#[test]
fn options() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let rows: [Row<'_>; 14] = [
        ("nlpsolve((x-1)^2+(y-2)^2,[x,y]=0..3)", BOWL, 100_000),
        (
            "nlpsolve((x-1)^2+(y-2)^2,[x,y]=0..3,nlp_precision=1e-12)",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_tolerance=1e-3)",
            BOWL,
            1_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_verbose=2)",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_presolve=false)",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=\"automatic\")",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=[\"nm\",\"size\"=1.0])",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=\"de\")",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=\"cobyla\")",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve(-(x-1)^2-(y-2)^2,x=0..3,y=0..3,nlp_maximize=true)",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve(-(x-1)^2-(y-2)^2,x=0..3,y=0..3,nlp_maximize)",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,nlp_initialpoint=[x=0..1,y=0..1])",
            BOWL,
            100_000,
        ),
        (
            "nlpsolve(x^2+y^2,[x+y=1],nlp_method=\"cobyla\")",
            "[1/2,[x=1/2,y=1/2]]",
            100_000,
        ),
        (
            "nlpsolve(2-1/120*x1*x2*x3*x4*x5,[x1<=1,x2<=2,x3<=3,x4<=4,x5<=5],assume=nlp_nonnegative)",
            "[1,[x1=1,x2=2,x3=3,x4=4,x5=5]]",
            10_000,
        ),
    ];
    let mut wrong = Vec::new();
    for (input, closed, scale) in rows {
        let closed = ours(closed).expect("the closed form");
        let theirs = seeded_text(giac, input).and_then(|t| ours(&t).ok());
        match (theirs, ours(input)) {
            (Some(t), Ok(o))
                if same(&ctx, &o, &closed, 10_000_000_000) && same(&ctx, &t, &closed, scale) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where ours departs from giac (CAS-604), giac seeded and its answer
/// copied (C0007): giac's interior point finds no solution under an
/// inequality, with `assume=nlp_nonnegative` and with `nlp_method="ipt"`,
/// where ours answers the closed form; `nlp_iterationlimit=2` stops giac
/// short ("Failed to optimize at given precision" and its point), ours
/// leaves `NLPSolve[…]` with a warning.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let failed = "GIAC_ERROR: Failed to find a solution. Error: Bad Argument Value";
    let mut wrong = Vec::new();
    for (input, closed) in [
        ("nlpsolve((x-1)^2+(y-2)^2,[x+y<=2])", "[1/2,[x=1/2,y=3/2]]"),
        ("nlpsolve(x^2+y^2,[x+y>=1])", "[1/2,[x=1/2,y=1/2]]"),
        (
            "nlpsolve((x+1)^2+(y-2)^2,assume=nlp_nonnegative)",
            "[1,[x=0,y=2]]",
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=\"ipt\")",
            BOWL,
        ),
    ] {
        let closed = ours(closed).expect("the closed form");
        match (seeded_text(giac, input), ours(input)) {
            (Some(t), Ok(o)) if t == failed && same(&ctx, &o, &closed, 10_000_000_000) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    let short = "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_iterationlimit=2)";
    match (seeded_text(giac, short), ours(short)) {
        (Some(t), Ok(o))
            if t.starts_with("[\"Failed to optimize at given precision\",")
                && o.head().and_then(Expr::as_symbol) == Some(Symbol::NLP_SOLVE) => {}
        (t, o) => wrong.push(format!("{short}: ours {o:?}, giac's {t:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Integer and binary variables (CAS-605), where giac's `nlpsolve` fails
/// ("Failed to find a solution", seeded, as under seeds 7 and 1 and
/// unseeded, T4 probes), its answer copied (C0007): ours the known optimum
/// (`.run/1.5/check/known.py`, the T4.c2 prototype), exact at integer
/// points, the mixed program's continuous coordinate certified.
#[test]
fn integer_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let failed = "GIAC_ERROR: Failed to find a solution. Error: Bad Argument Value";
    let aide = "nlpsolve(x1^2/2+5x1+6x2,[6x1+x2^2<=23],[x1,x2]=1..5,";
    let cubic = "nlpsolve(4x1*x2*x3-x1*x2-x1*x3-x2*x3,";
    let lowest = "[-1,[x1=0,x2=1,x3=1]]";
    let rows = [
        (
            format!("{aide}assume=integer,maximize)"),
            "[63/2,[x1=3,x2=2]]",
        ),
        (
            format!("{aide}assume=nlp_integer,maximize)"),
            "[63/2,[x1=3,x2=2]]",
        ),
        (
            format!("{aide}nlp_integervariables=[x1,x2],maximize)"),
            "[63/2,[x1=3,x2=2]]",
        ),
        (
            format!("{aide}nlp_integervariables=[x1],maximize)"),
            "[39/2+6*sqrt(5),[x1=3,x2=sqrt(5)]]",
        ),
        (format!("{cubic}assume=nlp_binary)"), lowest),
        (
            format!("{cubic}assume=nlp_binary,maximize)"),
            "[1,[x1=1,x2=1,x3=1]]",
        ),
        (format!("{cubic}nlp_binaryvariables=[x1,x2,x3])"), lowest),
        (
            "nlpsolve((x-3)^2+(y+1)^2,assume=nlp_nonnegint)".to_owned(),
            "[1,[x=3,y=0]]",
        ),
        (
            "nlpsolve((x-3)^2+(y+1)^2,x=0..10,y=0..10,assume=integer)".to_owned(),
            "[1,[x=3,y=0]]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, known) in rows {
        let known = ours(known).expect("the known optimum");
        match (seeded_text(giac, &input), ours(&input)) {
            (Some(t), Ok(o)) if t == failed && same(&ctx, &o, &known, 10_000_000_000) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
