//! `CFSolve`, `zeros`, `cZeros` and `Surd` against giac (CAS-536,
//! CAS-537). giac's `cfsolve` prints 12 digits; ours are certified balls,
//! compared with giac's within `10⁻⁹ (1 + |b|)` as sets (giac's order of
//! `cfsolve(x^4+1=0,x)` is irregular) and with closed forms first
//! (cookbook C0006). The orders that differ are departures.

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

/// Whether two numbers, or lists of them in order, agree within `10⁻⁹ (1 +
/// |b|)` in the real and the imaginary part of `N[a − b]`. Not by `Abs` of
/// the gap: a complex ball's modulus is `sqrt(re² + im²)`, indeterminate
/// when that ball holds 0, as it does when the two agree; a real ball's
/// `Abs` is `[|m| ± r]`, which `Less` decides against the bound.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if is_list(a) || is_list(b) {
        return is_list(a)
            && is_list(b)
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(u, v)| near(ctx, u, v));
    }
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = Expr::call(Symbol::N, [scale / Expr::from(1_000_000_000_i64)]);
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone()]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let size = Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]);
        ctx.eval(&Expr::call(Symbol::LESS, [size, bound.clone()]))
            .ok()
            .and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
    })
}

/// Whether two lists hold, in order, values whose differences `Expand` to
/// 0: closed forms written differently, `(-1 + i√3)/2` and `-1/2 + i√3/2`.
fn expands_same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    is_list(a)
        && is_list(b)
        && a.args().len() == b.args().len()
        && a.args().iter().zip(b.args()).all(|(u, v)| {
            ctx.eval(&Expr::call(Symbol::EXPAND, [u.clone() - v.clone()]))
                .ok()
                == Some(Expr::from(0))
        })
}

/// Whether two lists hold the same numbers or tuples in any order, each
/// matched once ([`near`]); two numbers, whether they are near.
fn same_set(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if !is_list(a) || !is_list(b) {
        return near(ctx, a, b);
    }
    if a.args().len() != b.args().len() {
        return false;
    }
    let mut used = vec![false; b.args().len()];
    a.args().iter().all(|u| {
        let found = b
            .args()
            .iter()
            .enumerate()
            .find(|&(j, v)| !used[j] && near(ctx, u, v))
            .map(|(j, _)| j);
        found.inspect(|&j| used[j] = true).is_some()
    })
}

/// giac's `cfsolve` against ours, as sets: polynomials, `exp`, `sin`, a
/// system, a start point.
#[test]
fn values() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "cfsolve(x^2+1=0,x)",
        "cfsolve(x^2+1,x)",
        "cfsolve(exp(x)=-1,x)",
        "cfsolve(x^3-1=0,x)",
        "cfsolve(x^2-2=0,x)",
        "cfsolve(x^4+1=0,x)",
        "cfsolve([x^2+y^2=-1,x-y=0],[x,y])",
        "cfsolve(sin(x)=2,x)",
        "cfsolve(x^2+1=0,x=1+i)",
        "cfsolve(x^3-1=0,x=-1+i)",
        "cfsolve(x*exp(x)=1,x)",
        "cfsolve(cos(x)=x,x)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same_set(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Ours against closed forms, in our order.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in [
        ("cfsolve(exp(x)=-1,x)", "[i*pi]"),
        (
            "cfsolve(sin(x)=2,x)",
            "[pi/2+i*ln(2+sqrt(3)),pi/2-i*ln(2+sqrt(3))]",
        ),
        (
            "cfsolve(x^4+1=0,x)",
            "[-sqrt(2)/2-i*sqrt(2)/2,-sqrt(2)/2+i*sqrt(2)/2,\
             sqrt(2)/2-i*sqrt(2)/2,sqrt(2)/2+i*sqrt(2)/2]",
        ),
        (
            "cfsolve([x^2+y^2=-1,x-y=0],[x=1,y=1])",
            "[i*sqrt(2)/2,i*sqrt(2)/2]",
        ),
        ("cfsolve(exp(x)=2,x)", "[ln(2)]"),
        ("cfsolve(x^3-1=0,x=1)", "1"),
    ] {
        match (ours(input), ours(expected)) {
            (Ok(found), Ok(closed)) if near(&ctx, &found, &closed) => {}
            (found, closed) => wrong.push(format!("{input}: ours {found:?}, {closed:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `zeros` and `cZeros`, `Solve`'s and `CSolve`'s names, and `surd`,
/// where giac's answer is ours.
#[test]
fn agreements() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "zeros(x^2+1,x)",
        "surd(-8,3)",
        "surd(8,3)",
        "surd(2,3)",
        "surd(-2,3)",
        "surd(-27,3)",
        "surd(x,3)",
    ] {
        if let Err(e) = giac.agrees(input) {
            wrong.push(format!("{e:?}"));
        }
    }
    // giac's `(-1 + i√3)/2` is our `-1/2 + i√3/2`, in the same order.
    let input = "cZeros(x^3-1,x)";
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(ours)) if expands_same(&ctx, &ours, &theirs) => {}
        (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
    }
    // giac's principal fourth root, `(√2/2 + i√2/2) 16^(1/4)`, by value.
    let input = "surd(-16,4)";
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => {}
        (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it prints
/// it and ours as decided, the same set in another order. `zeros` keeps
/// `Solve`'s order, giac's reverses it; `cZeros` keeps `CSolve`'s, as
/// giac's `csolve`; `cfsolve` orders by real then imaginary part. And
/// giac's `[undef]` where `CSolve` finds no root and ours is `[]`.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    let cases = [
        ("zeros(x^2-1,x)", "[1,-1]", "[-1,1]"),
        ("zeros(x^3-6*x^2+11*x-6,x)", "[3,2,1]", "[1,2,3]"),
        ("cZeros(x^2+1,x)", "[-i,i]", "[i,-i]"),
        (
            "zeros([x^2+y^2-1,x-y],[x,y])",
            "[[sqrt(2)/2,sqrt(2)/2],[-sqrt(2)/2,-sqrt(2)/2]]",
            "[[-sqrt(2)/2,-sqrt(2)/2],[sqrt(2)/2,sqrt(2)/2]]",
        ),
        (
            "cZeros([x^2+y^2-1,x-y],[x,y])",
            "[[sqrt(2)/2,sqrt(2)/2],[-sqrt(2)/2,-sqrt(2)/2]]",
            "[[-sqrt(2)/2,-sqrt(2)/2],[sqrt(2)/2,sqrt(2)/2]]",
        ),
        (
            "cfsolve(x^4+1=0,x)",
            "[0.707106781187+0.707106781187*i,-0.707106781187+0.707106781187*i,\
             -0.707106781187-0.707106781187*i,0.707106781187-0.707106781187*i]",
            "[-sqrt(2)/2-i*sqrt(2)/2,-sqrt(2)/2+i*sqrt(2)/2,\
             sqrt(2)/2-i*sqrt(2)/2,sqrt(2)/2+i*sqrt(2)/2]",
        ),
        (
            "cfsolve(sin(x)=2,x)",
            "[1.57079632679-1.31695789692*i,1.57079632679+1.31695789692*i]",
            "[pi/2+i*ln(2+sqrt(3)),pi/2-i*ln(2+sqrt(3))]",
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
        if !near(&ctx, &printed, &recorded) {
            wrong.push(format!(
                "{input}: giac gives {printed:?}, recorded {answer}"
            ));
        }
        if !near(&ctx, &ours, &decided) {
            wrong.push(format!("{input}: ours {ours:?} is not {decided:?}"));
        }
        if !same_set(&ctx, &ours, &printed) {
            wrong.push(format!(
                "{input}: ours {ours:?} and giac's {printed:?} differ"
            ));
        }
        if near(&ctx, &ours, &printed) {
            wrong.push(format!("{input} agrees now"));
        }
    }
    let input = "cfsolve(exp(x)=0,x)";
    match (giac.eval(input), ours(input)) {
        (Ok(theirs), Ok(found)) if theirs == "[undef]" && found == ours("[]").unwrap() => {}
        (theirs, found) => wrong.push(format!("{input}: ours {found:?}, giac's {theirs:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
