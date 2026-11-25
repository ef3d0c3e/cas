//! Xcas's wavelet transforms against giac (CAS-634): `dwt` and `idwt`
//! (giac's answers in `crates/cas-core/proto/wavelet/giac/p1.cas` and `p2.cas`, the prototype
//! `crates/cas-core/proto/wavelet/proto_dwt.py`, logs/T9.c1.md). giac's decimals read back are
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

/// Each entry's real and imaginary parts, so [`near`] compares complex
/// values by parts (C0010, C0013).
fn by_parts(e: &Expr) -> Expr {
    let parts = leaves(e).into_iter().flat_map(|x| {
        [
            Expr::call(Symbol::RE, [x.clone()]),
            Expr::call(Symbol::IM, [x]),
        ]
    });
    Expr::call(Symbol::LIST, parts.collect::<Vec<_>>())
}

/// Each input: giac's answer, read back, is this CAS's value within giac's
/// error, entry by entry and part by part.
fn check_parts(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if near(&ctx, &by_parts(&theirs), &by_parts(&mine)) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's analytic signals (CAS-635, `crates/cas-core/proto/wavelet/giac/p3.cas`, the
/// prototypes `proto_emd.py` and `proto_analytic.py`, logs/T9.c3.md): `emd`
/// bracketed, its answer a sequence (C0023); the decimal 41-sample row
/// within giac's noise (`5.6·10⁻¹⁷` where its decimals differ from strtod's).
#[test]
fn analytic_signals() {
    check_parts(&[
        "[emd([0,1,0,-1,0,1,0,-1,0,1,0,-1,0])]",
        "[emd([0,0.289739961672,0.395600889145,0.805274305596,0.680499109408,\
         1.1440709714,0.772198812573,1.23498263843,0.637623918224,1.06498597541,\
         0.302062656876,0.681634550085,-0.156413647973,0.180943481355,\
         -0.625739692064,-0.316290261382,-0.987278893705,-0.693920605476,\
         -1.14516683232,-0.868443362401,-1.05036746099,-0.808727783336,\
         -0.713505671778,-0.542714521634,-0.20326848937,-0.14951370717,\
         0.369313524884,0.261077675498,0.876943694741,0.576535528022,\
         1.20619897585,0.71011413673,1.28443457025,0.622842159983,\
         1.09693515003,0.333563368693,0.690164036867,-0.0854759492027,\
         0.160342996143,-0.526428575903,-0.369837755626])]",
        "[emd([1,2,3,4])]",
        "[emd([1,2,1,3,1,2,1,3,1])]",
        "[emd([1,2,1,3,1,2,1,3,1],threshold=0.5)]",
        "[emd([1,2,1,3,1,2,1,3,1],limit=1)]",
        "[emd([1,2,1,3,1,2,1,3,1],extrema=6)]",
        "[emd([1,2,1,3,1,2,1,3,1],extrema=7)]",
        "[emd([1/2,2,1,3,1.5,2,1,3,1])]",
        "[emd([5])]",
        "hilbert([1,0,0,0])",
        "hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3])",
        "hilbert([1,2,3])",
        "hilbert([1,i,0])",
        "hilbert(3)",
        "instfreq(hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]))",
        "instphase(hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]))",
        "instfreq([1,i,-1,-i])",
        "instphase([1,i,-1,-i])",
        "instfreq([1,exp(i/2),exp(i)])",
        "instphase([1,exp(i),exp(2*i),exp(3*i)])",
        "instfreq([1,2])",
        "instfreq([1])",
        "instphase([])",
        "instphase([1.0,-1.0,1.0])",
    ]);
}

/// Whether `Simplify[a − b]` is 0 (C0032, `Simplify` for giac's
/// `cos²+sin²`), else whether `a` and `b` agree within [`near`] at `x = t =
/// 7/10` and `13/10`.
fn same_function(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::SIMPLIFY,
        [Expr::call(
            Symbol::PLUS,
            [a.clone(), Expr::from(-1) * b.clone()],
        )],
    );
    if ctx.eval(&gap).is_ok_and(|g| g == Expr::from(0)) {
        return true;
    }
    [7, 13].into_iter().all(|k| {
        let v = Expr::from(k) / Expr::from(10);
        let rule = |s: &str| Expr::call(Symbol::RULE, [Expr::from(Symbol::new(s)), v.clone()]);
        let at = |f: &Expr| {
            Expr::call(
                Symbol::REPLACE_ALL,
                [f.clone(), Expr::call(Symbol::LIST, [rule("x"), rule("t")])],
            )
        };
        near(ctx, &at(a), &at(b))
    })
}

/// giac's expression forms of `instfreq` and `instphase` (p3.out): its
/// `cos²+sin²` and `floor(1/2)` are ours closed.
#[test]
fn analytic_expressions() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "instfreq(exp(i*x^2),x)",
        "instphase(exp(i*x^2),x)",
        "instfreq(exp(2*i*x))",
        "instphase(exp(2*i*x))",
        "instfreq(cos(x)+i*sin(x),x)",
        "instphase(cos(t)+i*sin(t),t)",
        "instfreq(exp(i*t^3),t)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if same_function(&ctx, &mine, &theirs) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures (p3.out): `emd` of symbolic data stays with giac's message as
/// the warning; `emd`'s options giac refuses stay; giac's transform of an
/// expression, `hilbert(sin(t),t)` = `-cos(t)`, is not ported: it stays.
#[test]
fn analytic_departures() {
    let giac = giac();
    let mut wrong = Vec::new();
    let refused = [
        "emd([x,1,2])",
        "emd([1,2,1,3,1,2,1,3,1],extrema=1)",
        "emd([1,2,1,3,1,2,1,3,1],threshold=-1)",
        "emd([1,2,1,3,1,2,1,3,1],residue=0)",
    ];
    for input in refused {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if !printed.starts_with("GIAC_ERROR") {
                wrong.push(format!("{input}: giac prints {printed}, recorded an error"));
            }
        }
    }
    for input in [refused[0], "hilbert(sin(t),t)"] {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            let recorded = if input == refused[0] {
                "GIAC_ERROR: Data must be numeric Error: Bad Argument Value"
            } else {
                "-cos(t)"
            };
            if collapsed(&printed) != recorded {
                wrong.push(format!(
                    "{input}: giac prints {printed}, recorded {recorded}"
                ));
            }
        }
    }
    for input in refused.iter().chain(&["hilbert(sin(t),t)"]) {
        let stays = evaluated(input).is_some_and(|(e, found, warnings)| {
            found.head() == e.head()
                && (*input != refused[0]
                    || warnings.iter().any(|w| {
                        matches!(w, Warning::Native { message: m, .. } if m == "Data must be numeric")
                    }))
        });
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
