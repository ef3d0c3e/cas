//! `minimax` against giac (CAS-599): giac's converged rows, ours' error
//! within 2% of giac's and our `p` within twice giac's error of giac's at
//! 51 points; closed forms (C0006); departures, where giac stops early or
//! loses the alternation, each with giac's answer as caseval prints it
//! (C0007) and ours checked against the closed form or the known error,
//! `|f − p| ≤ E` at 51 points. Polynomials of balls are compared by value,
//! never by `CoefficientList`. The seeded `rand=` row is in
//! `tests/minimax_rand.rs`.

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

/// `{E, p}` of an answer.
fn parts(answer: &Expr) -> Option<(Expr, Expr)> {
    match answer.args() {
        [error, p] if answer.head() == Some(&Expr::from(Symbol::LIST)) => {
            Some((error.clone(), p.clone()))
        }
        _ => None,
    }
}

/// Whether `|a − b| < bound` numerically.
fn within(ctx: &Context, a: &Expr, b: &Expr, bound: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(Symbol::N, [a.clone() - b.clone()])],
    );
    let bound = Expr::call(Symbol::N, [bound.clone()]);
    let claim = Expr::call(Symbol::LESS, [gap, bound]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether `|a(t) − b(t)| < bound` at the 51 points `t` evenly spaced on
/// `[lo, hi]`.
fn within_at(ctx: &Context, a: &Expr, b: &Expr, (lo, hi): (i64, i64), bound: &Expr) -> bool {
    let x = Expr::symbol("x");
    (0..=50).all(|k| {
        let t = Expr::from(50 * lo + k * (hi - lo)) / Expr::from(50);
        let at = |e: &Expr| e.substitute(&[(x.clone(), t.clone())]);
        within(ctx, &at(a), &at(b), bound)
    })
}

/// giac's converged rows: ours' error within 2% of giac's, and our `p`
/// within twice giac's error of giac's at 51 points.
#[test]
fn oracle() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, range) in [
        ("minimax(exp(x),x=0..1,3)", (0, 1)),
        ("minimax(sin(x),x=0..1,[2,1])", (0, 1)),
        ("minimax(abs(x)*sqrt(abs(x)),x=-2..2,15)", (-2, 2)),
        ("minimax(abs(x),x=-1..1,3)", (-1, 1)),
        ("minimax(x^3,x=-1..1,2)", (-1, 1)),
    ] {
        let theirs = giac.value(input).ok().as_ref().and_then(parts);
        let found = ours(input).ok().as_ref().and_then(parts);
        let (Some((eg, pg)), Some((eo, po))) = (&theirs, &found) else {
            wrong.push(format!("{input}: giac {theirs:?}, ours {found:?}"));
            continue;
        };
        if !within(&ctx, eo, eg, &(eg.clone() / Expr::from(50))) {
            wrong.push(format!("{input}: error ours {eo:?}, giac {eg:?}"));
        }
        if !within_at(&ctx, po, pg, range, &(Expr::from(2) * eg.clone())) {
            wrong.push(format!("{input}: p ours {po:?}, giac {pg:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Closed forms: `E` within `10⁻⁶ E`, `p` within `10⁻⁸` at 51 points.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let level = "(1-(exp(1)-1)*(1-ln(exp(1)-1)))/2";
    let mut wrong = Vec::new();
    for (input, range, error, p) in [
        (
            "minimax(x^3,x=-1..1,2)",
            (-1, 1),
            "1/4".to_owned(),
            "3/4*x".to_owned(),
        ),
        (
            "minimax(abs(x),x=-1..1,2)",
            (-1, 1),
            "1/8".to_owned(),
            "x^2+1/8".to_owned(),
        ),
        (
            "minimax(abs(x),x=-1..1,3)",
            (-1, 1),
            "1/8".to_owned(),
            "x^2+1/8".to_owned(),
        ),
        (
            "minimax(exp(x),x=0..1,1)",
            (0, 1),
            level.to_owned(),
            format!("(exp(1)-1)*x+1-{level}"),
        ),
    ] {
        let found = ours(input).ok().as_ref().and_then(parts);
        let (Some((eo, po)), Ok(e), Ok(q)) = (&found, ours(&error), ours(&p)) else {
            wrong.push(format!("{input}: ours {found:?}"));
            continue;
        };
        let tiny = Expr::from(10).pow(-6) * e.clone();
        if !within(&ctx, eo, &e, &tiny) {
            wrong.push(format!("{input}: error ours {eo:?}, closed form {error}"));
        }
        if !within_at(&ctx, po, &q, range, &Expr::from(10).pow(-8)) {
            wrong.push(format!("{input}: p ours {po:?}, closed form {p}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures, giac's answers as caseval prints them (C0007):
/// - `abs(x)` degree 2: giac's exchange loses the alternation (from
///   Chebyshev's nodes it interpolates, `E = 0`), `0.2134`; the best is
///   `x^2+1/8`, error `1/8` ([`closed_forms`]).
/// - `exp(x)` degree 1: giac stops at its rule `max |e| ≤ 1.02 min |e|`,
///   `0.106657`; the best has error `0.1059334` ([`closed_forms`]).
/// - `min(1/cosh(3*sin(x)),sin(9*x/10))` degree 30: giac's exchange stalls
///   at `0.0947`; ours `0.0319105` (a 200001-point sampling gives
///   `0.0319148`), `|f − p| ≤ E (1 + 10⁻³)` at 51 points.
#[test]
fn departures() {
    const ABS: &str = "[0.213388376846,1.35980080546e-12-4.56102925069e-08*x+1.17157282181*x^2]";
    const EXP: &str = "[0.106657010751,0.894790178235+1.71828182846*x]";
    const MIN: &str = "[0.0946608433887,0.0630612016616+1.16353544765*x-0.222952158736*x^2\
        -1.59077630128*x^3-0.590950335264*x^4+1.45241324598*x^5+0.648034887642*x^6\
        -0.74612837184*x^7-0.270845561389*x^8+0.226980074708*x^9+0.0561969032923*x^10\
        -0.0419279781947*x^11-0.00573746843166*x^12+0.00462130154892*x^13\
        +0.000200144655075*x^14-0.000280157253701*x^15+8.35945141035e-06*x^16\
        +7.31046479372e-06*x^17-5.9359457628e-07*x^18-1.11558977681e-08*x^19\
        -5.23102869259e-10*x^20-4.07753931655e-11*x^21-4.10945028729e-12*x^22\
        -4.75926392933e-13*x^23-6.312103181e-14*x^24-8.94300160879e-15*x^25\
        -1.37612002065e-15*x^26-2.20675392682e-16*x^27-3.73553514777e-17*x^28\
        -6.51408404269e-18*x^29-1.17925973676e-18*x^30]";
    let f = "min(1/cosh(3*sin(x)),sin(9*x/10))";
    let min = format!("minimax({f},x=-3..4,30)");
    let mut wrong = Vec::new();
    if let Some(giac) = giac() {
        for (input, answer) in [
            ("minimax(abs(x),x=-1..1,2)", ABS),
            ("minimax(exp(x),x=0..1,1)", EXP),
            (min.as_str(), MIN),
        ] {
            let text = giac.eval(input).unwrap_or_default();
            if text != answer {
                wrong.push(format!("{input}: giac prints {text}, recorded {answer}"));
            }
        }
    }
    let ctx = Context::new();
    let found = ours(&min).ok().as_ref().and_then(parts);
    match (&found, ours(f)) {
        (Some((eo, po)), Ok(fx)) => {
            let best = Expr::from(319_105) / Expr::from(10_000_000);
            if !within(&ctx, eo, &best, &(best.clone() / Expr::from(1000))) {
                wrong.push(format!("{min}: error ours {eo:?}"));
            }
            let zero = Expr::from(0);
            // 10⁻³: at the ends |f − p| = E, and 53-bit N of degree 30 at 4 is wider than 10⁻⁶ E.
            let bound = eo.clone() * (Expr::from(1) + Expr::from(10).pow(-3));
            if !within_at(&ctx, &(fx - po.clone()), &zero, (-3, 4), &bound) {
                wrong.push(format!("{min}: |f - p| above ours {eo:?}"));
            }
        }
        _ => wrong.push(format!("{min}: ours {found:?}")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
