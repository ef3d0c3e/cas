//! Optimization against giac (CAS-DOC-38, CAS-PLAN-32): giac's answers to
//! its examples of `extrema`, `minimize`, `maximize`, `lpsolve` and `tpsolve`,
//! read back, are this CAS's by value, lists entry by entry.

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

/// giac's answers and this CAS's, the same by [`same`].
fn check(inputs: &[&str]) {
    check_within(inputs, 10_000_000_000);
}

/// giac's answers and this CAS's, the same by [`same`] with numbers within
/// `1/scale`.
fn check_within(inputs: &[&str], scale: i64) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs, scale) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Local extrema (CAS-415).
#[test]
fn extrema() {
    check(&[
        "extrema(x^3-3*x,x)",
        "extrema(x^2+y^2-2*x,[x,y])",
        "extrema(x^3+y^3-3*x*y,[x,y])",
        "extrema(x^4,x)",
        "extrema(x*y*z,x+y+z=1,[x,y,z])",
        "extrema(x^2+y^2,x*y=1,[x=0..inf,y=0..inf])",
        "extrema(sin(x),x=0..10)",
    ]);
}

/// Global optima (CAS-416).
#[test]
fn optima() {
    check(&[
        "minimize(x^2-2*x+3,x)",
        "maximize(sin(x),x=0..4)",
        "minimize(sin(x),x=0..10,point)",
        "maximize(sin(x),x=0..10,point)",
        "minimize(x^3-3*x,x=-3..3)",
        "minimize(x^2+y^2,x+y=1,[x,y])",
        "minimize(2x^2-y^2+6y,x^2+y^2<=16,[x,y])",
        "maximize(2x^2-y^2+6y,x^2+y^2<=16,[x,y])",
        "minimize(x^2-3x+y^2+3y+3,[x=2..4,y=-4..-2],point)",
    ]);
}

/// Linear programs (CAS-417, CAS-418).
#[test]
fn linear_programs() {
    check(&[
        "lpsolve(2x+y-z+4,[x<=1,y>=2,x+3y-z=2,2x-y+z<=8,-x+y<=5])",
        "lpsolve(-7x+2y,[4x-12y<=20,-x+3y<=3],x=-5..5,y=0..inf,maximize)",
        "lpsolve(-x-y,[y<=3x+1/2,y<=-5x+2],assume=lp_nonnegative)",
        "lpsolve(-6x+4y+z,[5x-10y<=20,2z-3y=6,-x+3y<=3],x=1..20,y=0..inf)",
        "lpsolve(-5x-7y,[7x+y<=35,-x+3y<=6],assume=integer)",
        "lpsolve(8x1+11x2+6x3+4x4,[5x1+7x2+4x3+3x4<=14],assume=lp_binary,maximize)",
        "lpsolve(x+y,[x>=1/3,y>=1/2],assume=integer)",
    ]);
}

/// giac's answer to `input` with its random generator seeded (42: every
/// input below then succeeds; seed 0 fails the first) first, in the
/// same call: `nlpsolve` starts from random points, so its answer (and
/// whether it finds one: `nlpsolve(x^2+y^2,[x+y=1],maximize=false)` fails
/// in about 1 order in 20 of this binary's tests) depends on the draws
/// the tests running before it made in the shared instance (T5).
fn seeded(giac: &Giac, input: &str) -> Option<Expr> {
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
    ours(&answer[start..].replace("list[", "[")).ok()
}

/// Nonlinear programs (CAS-519): giac's examples that its `nlpsolve`
/// solves, its answers good to about 10⁻⁷ (COBYLA's tolerance), ours
/// certified; giac seeded ([`seeded`]).
#[test]
fn nonlinear_programs() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "nlpsolve(sin(x1+x2)+(x1-x2)^2-1.5x1+2.5x2+1,x1=-1.5..4,x2=-3..3)",
        "nlpsolve(sin(x)*cos(y),[x=0..3,y=0..3])",
        "nlpsolve(x^2+y^2,[x+y=1],maximize=false)",
        "nlpsolve(x*y,[x+y=2],maximize)",
        "nlpsolve((x-2)^2,[x=0..1])",
    ] {
        match (seeded(giac, input), ours(input)) {
            (Some(theirs), Ok(ours)) if same(&ctx, &ours, &theirs, 100_000) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Transportation problems (CAS-598): giac's probe rows, each bracketed
/// as the answer is a sequence (C0023).
#[test]
fn transportation() {
    check(&[
        "[tpsolve([10,20],[15,15],[[1,2],[3,1]])]",
        "[tpsolve([10,20],[10,10],[[1,2],[3,1]])]",
        "[tpsolve([12,17,11],[10,10,10,10],[[500,750,300,450],[650,800,400,600],[400,700,500,550]])]",
        "[tpsolve([95,70,165,165],[195,150,30,45,75],[[15,M,45,M,0],[12,40,M,M,0],[0,15,25,25,0],[M,0,M,12,0]])]",
        "[tpsolve([1,1,1,1],[1,1,1,1],[[10,12,9,11],[5,10,7,8],[12,14,13,11],[8,15,11,9]])]",
        "[tpsolve([2,4,6,4],[6,1,2,1,3],[[2,2,M,8,6],[3,2,2,3,0],[9,5,8,6,8],[2,8,5,2,2]])]",
        "[tpsolve([6,5],[5,6],[[M,0],[5,0]])]",
        "[tpsolve([3,1],[2,1,1],[[2,0,2],[3,2,M]])]",
        "[tpsolve([3],[1,1,2,6],[[5,0,3,6]])]",
        "[tpsolve([3,3,1],[2,3,2],[[3,3,3],[2,0,0],[1,1,0]])]",
        "[tpsolve([5],[5,6,2,5,5],[[2,0,0,2,0]])]",
        "[tpsolve([3,1,4,4],[4,1,1,5],[[1,1,2,0],[2,0,2,2],[2,1,1,0],[0,1,2,0]])]",
        "[tpsolve([3],[2,1],[[0,1]])]",
        "[tpsolve([1,1,1,1,1],[1,1,1,1,1],[[0,2,0,3,2],[0,0,2,3,3],[0,3,0,3,2],[3,2,2,1,1],[0,2,2,0,2]])]",
        "[tpsolve([1,1,1,1,1],[1,1,1,1,1],[[1,0,0,0,0],[0,0,0,1,0],[0,1,1,1,0],[0,1,1,1,1],[1,0,0,1,1]])]",
        "[tpsolve([1,1,1,1,1,1],[1,1,1,1,1,1],[[0,1,1,0,M,M],[1,1,0,0,1,0],[1,M,1,1,0,0],[0,M,0,0,M,1],[0,1,1,0,1,M],[0,1,0,M,0,1]])]",
        "[tpsolve([1,1,1],[1,1,1],[[14,15,12],[10,15,3],[4,1,15]])]",
        "[tpsolve([14,8,15,10,9],[18,4,16],[[0,0,0],[1,0,1],[1,M,M],[1,0,1],[0,0,1]])]",
        "[tpsolve([8,13,17,13,17,8],[4,3,18,2,13,3],[[7,1,16,8,8,8],[7,M,12,11,1,M],[10,5,20,1,M,17],[19,M,0,14,17,14],[13,4,19,9,17,M],[12,3,M,6,11,5]])]",
    ]);
}

/// Transportation inputs ours leaves (CAS-598): giac finds no feasible
/// solution, or takes a decimal cost, a negative supply or `2*M`; ours
/// stays `TPSolve[…]`.
#[test]
fn transportation_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "tpsolve([6],[6],[[M]])",
        "tpsolve([6,4,1],[5,2,1,4],[[M,1,3,3],[M,1,2,2],[2,0,1,0]])",
        "tpsolve([10,20],[15,15],[[1.5,2],[3,1]])",
        "tpsolve([10,-20],[15,15],[[1,2],[3,1]])",
        "tpsolve([10,20],[15,15],[[2*M,2],[3,1]])",
    ] {
        let theirs = giac.value(input);
        match ours(input) {
            Ok(e) if e.head().and_then(Expr::as_symbol) == Some(Symbol::TP_SOLVE) => {}
            ours => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
