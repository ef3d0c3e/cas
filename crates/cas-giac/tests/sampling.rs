//! Sampling and intervals against giac (CAS-593): `linspace`,
//! `linear_interpolate`, `interval2center` and `center2interval`. giac
//! answers doubles where this CAS answers exactly, so lists and ranges are
//! compared element by element, numbers within `10⁻¹⁰` relatively.

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

/// Whether `a` is within `10⁻¹⁰ max(1, |b|)` of `b`, lists and ranges
/// element by element.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let compound = |e: &Expr| {
        matches!(
            e.head().and_then(Expr::as_symbol),
            Some(Symbol::LIST | Symbol::RANGE)
        )
    };
    if compound(a) || compound(b) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| near(ctx, x, y));
    }
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

/// giac's answers, exact or doubles; `center2interval` bracketed, since
/// it answers a sequence (C0023).
#[test]
fn answers() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "linspace(0,1,11)",
        "linspace(1,2,3)",
        "linspace(2,1,3)",
        "linspace(0.0,1.0,3)",
        "linear_interpolate([[1,2,6,9],[3,4,6,7]],1,9,1)",
        "linear_interpolate([[1,2,6,9],[3,4,6,7]],1,9,3)",
        "linear_interpolate([[0,2,3],[1,5,4]],0,3,1/2)",
        "linear_interpolate([[0,1],[0,1]],0,1,0.25)",
        "interval2center(2..5)",
        "interval2center(3)",
        "interval2center(1.5..2)",
        "interval2center([1..3,5..7])",
        "[center2interval([2,5,9],1)]",
        "[center2interval([2,5,8])]",
        "[center2interval([1,2,4])]",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where giac errs and this CAS leaves the call: giac's error as caseval
/// prints it, whitespace collapsed (C0007), and the head that stays.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        (
            "linspace(1,2,1)",
            "GIAC_ERROR: linspace(1,2,1) Error: Invalid dimension",
            Symbol::LIN_SPACE,
        ),
        (
            "linear_interpolate([[0,2,3],[1,5,4]],-1,4,1)",
            "GIAC_ERROR: Values out of range Error: Bad Argument Value",
            Symbol::LINEAR_INTERPOLATE,
        ),
        (
            "center2interval([3])",
            "GIAC_ERROR: center2interval([3]) Error: Bad Argument Value",
            Symbol::CENTER_TO_INTERVAL,
        ),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, head) in cases {
        let printed = giac.eval(input).unwrap();
        if collapsed(&printed) != collapsed(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let stays = ours(input)
            .ok()
            .and_then(|e| e.head().and_then(Expr::as_symbol));
        if stays != Some(head) {
            wrong.push(format!("{input}: ours {:?} does not stay", ours(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
