//! Xcas's wavelet transforms against giac (CAS-634): `dwt` and `idwt`
//! (giac's answers in `.run/1.7/T9/giac/p1.out` and `p2.out`, the prototype
//! `.run/1.7/T9/proto_dwt.py`, logs/T9.c1.md). giac's decimals read back are
//! ours within its printed digits (C0006, C0076): its rounding noise
//! (`2.2·10⁻¹⁶`) against ours of about `10⁻³⁰`; Haar on exact data is exact
//! here (`a + b√2`). giac's errors stay with its message as the warning
//! (C0056).

use cas_core::{Context, EvalOptions, Expr, Symbol, View, Warning};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, decimal_balls, lower, parse};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    match e.view() {
        View::Apply { .. } if e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) => {
            e.args().iter().flat_map(leaves).collect()
        }
        _ => vec![e.clone()],
    }
}

/// Whether `a` and `b` have the same shape and entries within giac's
/// printed digits (C0006).
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let numeric = |e: &Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
            .unwrap()
    };
    let (t, o) = (leaves(&numeric(a)), leaves(&numeric(b)));
    t.len() == o.len()
        && t.iter().zip(&o).all(|(a, b)| {
            let gap = Expr::call(
                Symbol::ABS,
                [Expr::call(
                    Symbol::PLUS,
                    [b.clone(), Expr::from(-1) * a.clone()],
                )],
            );
            let scale = Expr::call(
                Symbol::MAX,
                [Expr::from(1), Expr::call(Symbol::ABS, [a.clone()])],
            );
            let bound = scale * Expr::from(10).pow(-9);
            ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
                .ok()
                .and_then(|e| e.as_symbol())
                == Some(Symbol::TRUE)
        })
}

/// Each input: giac's answer, read back, is this CAS's value within giac's
/// error ([`near`]).
fn check_near(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if near(&ctx, &theirs, &mine) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answers (p1.out, p2.out): every family, members, centred,
/// inverse, padding, 2-D standard and `image`, and round trips.
#[test]
fn wavelets() {
    check_near(&[
        "dwt([1,2,3,4])",
        "dwt([1,2,3,4],4)",
        "dwt([1,2,3,4],\"daubechies\",4)",
        "dwt([1,2,3,4],\"haar\")",
        "dwt([1,2,3,4],\"bspline\")",
        "dwt([1,2,3,4,5,6,7,8],6)",
        "dwt([0.0,0.5,0.25,1.0,0.75,-0.25,-0.5,0.0])",
        "dwt([[15,25,5,70],[21,33,12,7],[64,18,50,10],[5,17,35,28]])",
        "dwt([[1,2],[3,4]])",
        "dwt([[1,2],[3,4]],image)",
        "dwt([1,2,3])",
        "dwt([1/2,1/3,1/4,1/5])",
        "dwt([1,2,3,4],center)",
        "dwt([1,2,3,4],\"haar\",center)",
        "dwt([1,2,3,4,5,6,7,8],center)",
        "dwt([1,2,3,4,5,6,7,8],\"bspline\",center)",
        "dwt([[1,2,3],[4,5,6]])",
        "dwt([[1,2,3,4],[4,5,6,7],[1,0,0,1],[2,2,1,1]],image)",
        "dwt([1,2,3,4,5,6,7,8],\"bspline\",309)",
        "idwt([1,2,3,4,5,6,7,8],\"bspline\",202,center)",
        "dwt([])",
        "dwt([5])",
        "dwt([1,2,3,4,5,6,7,8],\"daubechies\",20)",
        "dwt([1,2,3,4,5,6,7,8],\"bspline\",105)",
        "dwt([1,2,3,4,5,6,7,8],\"bspline\",208,center)",
        "dwt([1,2,3,4,5,6,7,8],\"bspline\",301)",
        "idwt([1,2,3,4,5,6,7,8],\"daubechies\",8,center)",
        "idwt([[1,2,3,4],[4,5,6,7],[1,0,0,1],[2,2,1,1]],image)",
        "idwt([[1,2,3,4],[4,5,6,7],[1,0,0,1],[2,2,1,1]],\"haar\",center)",
        "dwt([[1,2,3,4],[4,5,6,7],[1,0,0,1],[2,2,1,1]],\"bspline\",image,center)",
        "dwt([1.5,2,3,4],\"haar\")",
        "dwt([1,2,3,4],center,image)",
        "dwt([[1,2],[3,4],[5,6]],\"haar\")",
        "dwt([1,2,3,4,5],\"haar\")",
        "idwt(dwt([1,2,3,4],\"haar\"),\"haar\")",
        "idwt(dwt([1,2,3,4,5,6,7,8]))",
        "idwt(dwt([[15,25,5,70],[21,33,12,7],[64,18,50,10],[5,17,35,28]]))",
    ]);
}

/// Haar on exact data is exact here, giac's decimals above: ours against
/// the exact forms.
#[test]
fn exact_haar() {
    let rows = [
        ("dwt([1,2,3,4],\"haar\")", "[5,-2,-sqrt(2)/2,-sqrt(2)/2]"),
        (
            "dwt([1,2,3,4],\"haar\",center)",
            "[5,0,3*sqrt(2)/2,-sqrt(2)/2]",
        ),
        ("idwt(dwt([1,2,3,4],\"haar\"),\"haar\")", "[1,2,3,4]"),
        (
            "idwt(dwt([[1,2],[3,4],[5,6]],\"haar\"),\"haar\")",
            "[[1,2,0,0],[3,4,0,0],[5,6,0,0],[0,0,0,0]]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, form) in rows {
        match (ours(input), ours(form)) {
            (Ok(found), Ok(expected)) if found == expected => {}
            (found, expected) => wrong.push(format!("{input}: {found:?}, expected {expected:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `text` read as [`ours`] reads it, and its value with the warnings of
/// its evaluation.
fn evaluated(text: &str) -> Option<(Expr, Expr, Vec<Warning>)> {
    let parsed = parse(text);
    if !parsed.errors().is_empty() {
        return None;
    }
    let lowered = lower(&parsed, Names::english());
    let [statement] = &lowered.statements[..] else {
        return None;
    };
    let e = decimal_balls(&statement.expr, cas_core::DEFAULT_WORKING_PRECISION);
    let mut ctx = Context::new();
    ctx.set_language(Some(Hook::english()));
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    Some((e, found.ok()?, report.warnings))
}

/// Whitespace collapsed (C0047).
fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Departures that agree: giac's errors (p1.out, p2.out), as caseval prints
/// them, where ours stays with giac's message as the warning.
#[test]
fn errors_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let rows = [
        ("dwt([1,2,3,4],2)", value("Wavelet not implemented")),
        ("dwt([1,2,3,4],haar)", value("Invalid argument 2")),
        ("dwt([1,2,3,4],centre=true)", value("Invalid argument 2")),
        ("dwt([x,1,2,3])", value("Data must be numeric")),
        ("dwt([1,2,3,4],\"foo\")", value("Unsupported wavelet type")),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, error) in &rows {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
        let message = error
            .strip_prefix("GIAC_ERROR: ")
            .and_then(|m| m.split_once(". Error:"))
            .map_or("", |(m, _)| m);
        let stays = evaluated(input).is_some_and(|(e, found, warnings)| {
            found.head() == e.head()
                && warnings.iter().any(
                    |w| matches!(w, Warning::Native { message: m, .. } if collapsed(m) == message),
                )
        });
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {message:?}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
