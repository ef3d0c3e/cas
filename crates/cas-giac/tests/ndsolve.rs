//! `NDSolve` against giac's `odesolve` (CAS-DOC-35 §5, CAS-PLAN-29):
//! giac integrates in doubles by a Runge–Kutta method, whose answers to
//! equations with closed forms err by up to about `2·10⁻⁹`; ours enclose
//! the solution, and the two agree within `10⁻⁸`, relative. giac 2.0.0
//! refuses the systems of its manual, `odesolve([z,-y],[t,[y,z]],…)`.

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

/// Whether two answers, lists of one number, agree within `10⁻⁸ (1 +
/// |b|)`.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    let (Some(a), Some(b)) = (
        (a.head() == Some(&list))
            .then(|| a.args().first())
            .flatten(),
        (b.head() == Some(&list))
            .then(|| b.args().first())
            .flatten(),
    ) else {
        return false;
    };
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(100_000_000_i64);
    let gap = Expr::call(Symbol::ABS, [a.clone() - b.clone()]);
    ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
        .ok()
        .and_then(|e| e.as_symbol())
        == Some(Symbol::TRUE)
}

#[test]
fn values() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "odesolve(sin(t*y),[t,y],[0,1],2)",
        "odesolve(y,t=0..1,y,1)",
        "odesolve(-2*t*y,t=0..3,y,1)",
        "odesolve(y*(1-y),t=0..5,y,1/2)",
        "odesolve(t^2-y,[t,y],[0,1],1)",
        "odesolve(y,[t,y],[1,1],0)",
        "odesolve(cos(t)*y+1,[t,y],[0,0],2)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if close(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whether two lists of points `[x, y]` agree within `10⁻⁴ (1 + |b|)`.
/// giac's `bvpsolve` integrates with a fixed step `h` and errs by about
/// `h⁴`: on `[1, 3]` by `4·10⁻²` at `n = 4` (`h = 1/2`) and by `4·10⁻⁵`
/// at its manual's `n = 20` (`h = 1/10`), so it is compared at `h = 1/20`
/// only (CAS-520).
fn close_points(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() != Some(&list) || b.head() != Some(&list) || a.args().len() != b.args().len() {
        return false;
    }
    a.args().iter().zip(b.args()).all(|(p, q)| {
        p.args().len() == 2
            && q.args().len() == 2
            && p.args().iter().zip(q.args()).all(|(u, v)| {
                let scale = Expr::from(1) + Expr::call(Symbol::ABS, [v.clone()]);
                let bound = scale / Expr::from(10_000_i64);
                let gap = Expr::call(Symbol::ABS, [u.clone() - v.clone()]);
                ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
                    .ok()
                    .and_then(|e| e.as_symbol())
                    == Some(Symbol::TRUE)
            })
    })
}

/// The closed form `solution` at the `n + 1` points from `a` to `b`, as
/// `bvpsolve` lists them: `evalf([[x_0, y(x_0)], …])`.
fn grid(a: i64, b: i64, n: i64, solution: fn(&str) -> String) -> String {
    let points: Vec<String> = (0..=n)
        .map(|k| {
            let point = format!("({a}+{k}*({b}-{a})/{n})");
            format!("[{point},{}]", solution(&point))
        })
        .collect();
    format!("evalf([{}])", points.join(","))
}

/// Boundary value problems (CAS-520): ours against the closed forms on
/// small grids, and giac's `bvpsolve`, by shooting in doubles, against
/// ours where its step is small enough for its error to be inside `10⁻⁴`.
#[test]
fn boundary_values() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    let sine: fn(&str) -> String = |x| format!("sin({x})/sin(1)");
    let hyperbolic: fn(&str) -> String = |x| format!("(sinh(1-{x})+2*sinh({x}))/sinh(1)");
    let rational: fn(&str) -> String = |x| format!("{x}^2+16/{x}");
    for (input, expected) in [
        ("bvpsolve(-y,[x=0..1,y],[0,1],4)", grid(0, 1, 4, sine)),
        ("bvpsolve(y,[x=0..1,y],[1,2],4)", grid(0, 1, 4, hyperbolic)),
        (
            "bvpsolve((32+2x^3-y*diff(y(x),x))/8,[x=1..3,y],[17,43/3],4)",
            grid(1, 3, 4, rational),
        ),
        (
            "bvpsolve(-y,[x=0..1,y],[0,1],2,output=list)",
            grid(0, 1, 2, sine),
        ),
    ] {
        match (ours(input), ours(&expected)) {
            (Ok(ours), Ok(closed)) if close_points(&ctx, &ours, &closed) => {}
            (ours, closed) => {
                wrong.push(format!(
                    "{input}: ours {ours:?}, the closed form's {closed:?}"
                ));
            }
        }
    }
    if let Some(giac) = giac() {
        for input in [
            "bvpsolve(-y,[x=0..1,y],[0,1],20)",
            "bvpsolve(y,[x=0..1,y],[1,2],20)",
            "bvpsolve((32+2x^3-y*diff(y(x),x))/8,[x=1..3,y],[17,43/3],40)",
            "bvpsolve(-y,[x=0..1,y],[0,1],20,output=list)",
        ] {
            match (giac.value(input), ours(input)) {
                (Ok(theirs), Ok(ours)) if close_points(&ctx, &ours, &theirs) => {}
                (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it (its doubles, compared within `close_points`), and the closed
/// form this CAS's answer is.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    let cases = [
        // CAS-520: y = x^2+16/x; giac's step 1/2 errs by 4·10⁻² at
        // x = 3/2, ours encloses the solution.
        (
            "bvpsolve((32+2x^3-y*diff(y(x),x))/8,[x=1..3,y],[17,43/3],4)",
            "[[1.0,17.0],[1.5,12.9604508108],[2.0,12.0235376441],\
             [2.5,12.6584840485],[3.0,14.3333333333]]",
            "[[1,17],[3/2,155/12],[2,12],[5/2,253/20],[3,43/3]]",
        ),
    ];
    for (input, answer, decided) in cases {
        let (Ok(printed), Ok(recorded), Ok(ours), Ok(decided)) = (
            giac.value(input),
            cas_giac::ours(answer),
            cas_giac::ours(input),
            cas_giac::ours(decided),
        ) else {
            wrong.push(format!("{input}: an answer does not read"));
            continue;
        };
        if !close_points(&ctx, &printed, &recorded) {
            wrong.push(format!(
                "{input}: giac gives {printed:?}, recorded {answer}"
            ));
        }
        if !close_points(&ctx, &ours, &decided) {
            wrong.push(format!("{input}: ours {ours:?} is not {decided}"));
        }
        if close_points(&ctx, &ours, &printed) {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
