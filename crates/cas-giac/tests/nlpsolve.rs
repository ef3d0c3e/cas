//! `nlpsolve` against giac (CAS-604, CAS-605): its options matched by their names,
//! linear programs answered exactly as giac's `lpsolve` answers them
//! (optimization.cc:6080-6085), and the inputs where ours departs. The
//! rows whose giac answer draws random numbers are in `nlpsolve_seeded.rs`.

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

/// Linear programs (CAS-604): giac hands a linear objective under linear
/// constraints to `lpsolve`, ours to `LPSolve`; both exact, the optimum
/// unique (T4 prototype), the point in sorted order.
#[test]
fn linear_programs() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "nlpsolve(x+y,[x+2y<=4,3x+y<=6],x=0..10,y=0..10,maximize)",
        "nlpsolve(y+x,[2y+x<=4,y+3x<=6],assume=nlp_nonnegative,maximize)",
        "nlpsolve(b+a,[a+2b<=4,3a+b<=6],assume=nlp_nonnegative,maximize)",
        "nlpsolve(2x+3y,[x+y>=1],assume=nlp_nonnegative)",
        "nlpsolve(x-y,[x+y<=2.5,x>=0,y>=0],nlp_maximize=true)",
        "nlpsolve(y,[x+y>=1,x<=1/2],assume=nlp_nonnegative)",
        "nlpsolve(2x+y,[x+2y<=4.5,3x+y<=6.5],assume=nlp_nonnegint,maximize)",
        "nlpsolve(2x+y,[x+2y<=9/2,3x+y<=13/2],[x,y]=0..10,assume=integer,maximize)",
        "nlpsolve(2x+y,[x+2y<=9/2,3x+y<=13/2],assume=nlp_binary,maximize)",
        "nlpsolve(x+y,[x+2y<=4.5,3x+y<=6.5],assume=nlp_nonnegative,nlp_integervariables=[x],maximize)",
        "nlpsolve(2x+y,[x+y<=1.5],x=0..5,y=0..5,nlp_binaryvariables=[x,y],maximize)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(o)) if theirs == o => {}
            (theirs, o) => wrong.push(format!("{input}: ours {o:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A nonlinear program whose giac answer draws nothing from an initial
/// point: the same under seeds 42 and 1 and unseeded (T4 probes), good to
/// about 10⁻¹¹; ours certified, checked against the closed form `0` at
/// `(1, 2)` first.
#[test]
fn from_an_initial_point() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let input = "nlpsolve((x-1)^2+(y-2)^2,nlp_initialpoint=[x=5,y=5])";
    let closed = ours("[0,[x=1,y=2]]").expect("the closed form");
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(o))
            if same(&ctx, &o, &closed, 10_000_000_000) && same(&ctx, &theirs, &closed, 100_000) => {
        }
        (theirs, o) => panic!("{input}: ours {o:?}, giac's {theirs:?}"),
    }
}

/// An integer program giac answers (CAS-605): one variable, its points
/// enumerated (optimization.cc:5972-6102), the same under seeds 42 and 1
/// and unseeded (T4 probes); ours exact, `1/9` at `0`, giac's within its
/// error (C0006). The rows where giac fails are `nlpsolve_seeded.rs`'s.
#[test]
fn integer_programs() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let input = "nlpsolve((x-1/3)^2,x=-5..5,assume=integer)";
    let exact = ours("[1/9,[x=0]]").expect("the optimum");
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(o)) if o == exact && same(&ctx, &theirs, &exact, 100_000_000) => {}
        (theirs, o) => panic!("{input}: ours {o:?}, giac's {theirs:?}"),
    }
}

/// Where ours departs from giac (CAS-604), giac's answer copied (C0007):
/// an infeasible linear program is `LPSolve`'s `[]` where giac errs; an
/// unbounded one `LPSolve`'s `[+infinity, []]` where giac answers
/// `+infinity`; ours takes `nlp_nonnegative=true` and `nlp_verbose` and
/// `nlp_presolve` alone, which giac rejects; an unknown method string or a
/// method that is a symbol leaves ours `NLPSolve[…]` (with a warning) where
/// giac errs.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let invalid = "GIAC_ERROR: Invalid argument 3. Error: Bad Argument Value";
    let answered = [
        (
            "nlpsolve(x,[x>=1,x<=0])",
            "GIAC_ERROR: The problem is infeasible. Error: Bad Argument Value",
            "[]",
        ),
        (
            "nlpsolve(x+y,[x+y>=1],assume=nlp_nonnegative,maximize)",
            "+infinity",
            "[inf,[]]",
        ),
        (
            "nlpsolve((x+1)^2+(y-2)^2,x=0..3,y=0..3,nlp_nonnegative=true)",
            invalid,
            "[1,[x=0,y=2]]",
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_verbose)",
            invalid,
            "[0,[x=1,y=2]]",
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_presolve)",
            invalid,
            "[0,[x=1,y=2]]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, theirs, expected) in answered {
        let expected = ours(expected).expect("the expected answer");
        match (giac.eval(input), ours(input)) {
            (Ok(t), Ok(o)) if t == theirs && same(&ctx, &o, &expected, 10_000_000_000) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    for (input, theirs) in [
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=cobyla)",
            invalid,
        ),
        (
            "nlpsolve((x-1)^2+(y-2)^2,x=0..3,y=0..3,nlp_method=\"foo\")",
            "GIAC_ERROR: Method not supported. Error: Bad Argument Value",
        ),
    ] {
        match (giac.eval(input), ours(input)) {
            (Ok(t), Ok(o))
                if t == theirs && o.head().and_then(Expr::as_symbol) == Some(Symbol::NLP_SOLVE) => {
            }
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
