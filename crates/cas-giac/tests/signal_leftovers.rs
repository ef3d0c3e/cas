//! Xcas's pulses `sinc`, `tri`, `boxcar` and the Hilbert transform of an
//! expression, `hilbert(f,t)` and `hilbert(f)`, against giac (CAS-648;
//! giac's answers in `crates/cas-core/proto/wavelet/giac/p5.cas` and
//! `p5.out`). Exact answers, `undef` (C0081) and giac's inert
//! `Hilbert(f,t)` by `Giac::agrees`; expressions by `Simplify` of the
//! difference (C0032), else values within giac's digits at two points
//! (C0006). Departures: rows where our Fourier transforms leave a term giac
//! closes are ours inert, giac's lists of `sinc` and `tri` its oddities.
//! `stft`, `istft` and `hht(…, output=matrix)` (giac's answers in `p6` and
//! `p7`, checked by `proto_stft.py` and `proto_hht.py`): numbers within
//! giac's digits and FFT noise (C0084); giac's errors stay with its message;
//! departures: any function as `stft`'s window, `hht`'s plots and MIDI
//! ranges, and `bins=n` (giac's 161 rows).

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

/// Whether `Simplify[a − b]` is 0 (C0032), else whether `a` and `b` agree
/// within [`near`] at `x = t = s = a = 7/10` and `13/10`.
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
                [
                    f.clone(),
                    Expr::call(Symbol::LIST, [rule("x"), rule("t"), rule("s"), rule("a")]),
                ],
            )
        };
        near(ctx, &at(a), &at(b))
    })
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

/// Exact answers, `undef` and the inert `Hilbert(f,t)`: giac's answer read
/// back is ours (p5.out).
#[test]
fn exact_rows() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [
        "sinc(0)",
        "sinc(0.0)",
        "sinc(pi)",
        "sinc(2*pi)",
        "sinc(pi/2)",
        "tri(1/2)",
        "tri(-1/2)",
        "tri(-1/3)",
        "tri(0)",
        "tri(1)",
        "tri(-1)",
        "tri(2)",
        "boxcar(1,2,3/2)",
        "boxcar(1,2,1)",
        "boxcar(1,2,2)",
        "boxcar(1,2,0)",
        "boxcar(2,1,3/2)",
        "boxcar(0,1,1/2)",
        "boxcar(0,1,0)",
        "boxcar(1,1,1)",
        "boxcar(a,b,x)",
        "boxcar(1,2,x)",
        "boxcar(0,1,x)",
        "hilbert(1,t)",
        "hilbert(t^2,t)",
        "hilbert(t,t)",
        "hilbert(x)",
        "hilbert(exp(-t^2),t)",
        "Hilbert(exp(-t^2),t)",
        "hilbert(f(t),t)",
    ] {
        if let Err(e) = giac.agrees(input) {
            wrong.push(format!("{e:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Expressions and decimals (p5.out): `Simplify` of the difference, else
/// values at two points.
#[test]
fn expression_rows() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "sinc(x)",
        "sinc(t)",
        "sinc(pi*x)",
        "sinc(0.5)",
        "sinc(1.5)",
        "diff(sinc(x),x)",
        "tri(x)",
        "tri(t)",
        "tri(x-1)",
        "tri(0.25)",
        "boxcar(1,2,1.5)",
        "boxcar(-1/2,1/2,x)-rect(x)",
        "fourier(sinc(x),x,s)",
        "hilbert(cos(t),t)",
        "hilbert(sin(t),t)",
        "hilbert(sin(x))",
        "hilbert(cos(2*t),t)",
        "hilbert(cos(3*t),t)",
        "hilbert(sin(3*t),t)",
        "hilbert(2*cos(t),t)",
        "hilbert(cos(t)+sin(2*t),t)",
        "hilbert(sin(a*t),t)",
        "hilbert(cos(a*t),t)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if same_function(&ctx, &mine, &theirs) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's "Bad Argument Type" for `boxcar` without three arguments: ours
/// stays with the message as the warning.
#[test]
fn errors_stay() {
    let giac = giac();
    let mut wrong = Vec::new();
    for input in ["boxcar(1,2)", "boxcar(1,2,3,4)"] {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if !(printed.starts_with("GIAC_ERROR") && printed.contains("Bad Argument Type")) {
                wrong.push(format!("{input}: giac prints {printed}, recorded an error"));
            }
        }
        let stays = evaluated(input).is_some_and(|(e, found, warnings)| {
            found.head() == e.head()
                && warnings.iter().any(
                    |w| matches!(w, Warning::Native { message: m, .. } if m == "Bad Argument Type"),
                )
        });
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with the warning",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures (p5.out): giac closes the transform where our
/// `FourierTransform` or `InverseFourierTransform` leaves a term (rational
/// functions, whose transform holds `e^(−|ω|)`; `sinc`;
/// the box); ours is the inert `Hilbert(f,t)`.
#[test]
fn hilbert_departures() {
    let giac = giac();
    let mut wrong = Vec::new();
    let rows = [
        ("hilbert(1/(1+t^2),t)", "t/(t^2+1)"),
        ("hilbert(1/(4+t^2),t)", "t/(2*t^2+8)"),
        ("hilbert(t/(1+t^2),t)", "-1/(t^2+1)"),
        ("hilbert(x/(1+x^2),x)", "-1/(x^2+1)"),
        ("hilbert(1/(1+t^4),t)", "(t^3*sqrt(2)+t*sqrt(2))/(2*t^4+2)"),
        ("hilbert(sinc(x))", "(-cos(x)+1)/x"),
        ("hilbert(sinc(t),t)", "(-cos(t)+1)/t"),
        ("hilbert(sin(t)/t,t)", "(-cos(t)+1)/t"),
        ("hilbert(boxcar(0,1,x),x)", "(-ln(abs(-x+1))+ln(abs(x)))/pi"),
        (
            "hilbert(rect(x),x)",
            "(-ln(abs(-2*x+1)/2)+ln(abs(-2*x-1)/2))/pi",
        ),
    ];
    for (input, recorded) in rows {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != recorded {
                wrong.push(format!(
                    "{input}: giac prints {printed}, recorded {recorded}"
                ));
            }
        }
        let inert = evaluated(input).is_some_and(|(_, found, _)| {
            found.head().and_then(Expr::as_symbol) == Some(Symbol::HILBERT)
        });
        if !inert {
            wrong.push(format!(
                "{input}: ours {:?}, wanted Hilbert(f,t)",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures (p5.out): giac's `sinc` of a list is `[undef]` and its `tri`
/// of a list a `poly1`; ours maps them over the list.
#[test]
fn list_departures() {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, recorded, mine) in [
        ("sinc([0,pi,x])", "[undef]", "[1,0,sin(x)/x]"),
        ("tri([0,1/2])", "poly1[1/2]", "[1,1/2]"),
    ] {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != recorded {
                wrong.push(format!(
                    "{input}: giac prints {printed}, recorded {recorded}"
                ));
            }
        }
        if ours(input) != ours(mine) {
            wrong.push(format!("{input}: ours {:?}, wanted {mine}", ours(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `stft` and `istft` (CAS-648; giac's answers in `p6.cas` and `p6.out`,
/// checked by `proto_stft.py`): entries within giac's digits and its FFT
/// noise (C0084, C0006), complex ones by the modulus of the difference
/// (C0010).
#[test]
fn stft_rows() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let s = "[0.3,-0.1,0.7,0.2,-0.5,0.4,0.9,-0.6,0.1,0.8,-0.2,0.5,-0.7,0.6,0.2,-0.4]";
    let mut wrong = Vec::new();
    for input in [
        "stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3],2)".to_owned(),
        "stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3])".to_owned(),
        "stft([1,2,3,4,5,6,7,8],4)".to_owned(),
        "stft([1,2,3,4,5,6,7,8],2)".to_owned(),
        "stft([1,2,3,4,5,6,7,8,9,10],4)".to_owned(),
        "stft([1,2,3,4],4)".to_owned(),
        "stft([1,2,3,4,5,6,7,8],8)".to_owned(),
        "stft([1,2,3,4,5,6,7,8],4,hann_window)".to_owned(),
        format!("stft({s},4,blackman_window)"),
        format!("stft({s},8,bartlett_hann_window)"),
        format!("stft({s},2)"),
        "stft([1,2+i,3,4],2)".to_owned(),
        "stft([1/2,1/3,1/4,1/5],2)".to_owned(),
        "stft([1,2,3,4,5,6,7],2)".to_owned(),
        "stft([1,2,3,4,5,6,7,8],2,hann_window)".to_owned(),
        "istft(stft([1,2,3,4,5,6,7,8],2))".to_owned(),
        "istft(stft([1,2,3,4,5,6,7,8],4))".to_owned(),
        "istft(stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3],2))".to_owned(),
        "istft([[1,2],[3,4]])".to_owned(),
        "istft([[1+i,2],[3,4-2*i]])".to_owned(),
        "istft([[0.5,-1,2,0.25]])".to_owned(),
        "istft([[1,2,3]])".to_owned(),
        format!("istft(stft({s},4,blackman_window))"),
    ] {
        match (giac.value(&input), ours(&input)) {
            (Ok(theirs), Ok(mine)) if near(&ctx, &theirs, &mine) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `hht(…, output=matrix)` (CAS-648; giac's answers in `p7.cas` and
/// `p7.out`, checked by `proto_hht.py`): the matrices within giac's digits.
/// The top of the range lands one bin low as giac's does (`[[1,-1,1,-1]]`:
/// row 98 of 100); a range that is not positive or not numeric gives zeros.
#[test]
fn hht_rows() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let p = "[0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]";
    let q = "[0.1,0.4,0.2,-0.3,-0.1,0.2,0.6,0.1,-0.2]";
    let mut wrong = Vec::new();
    for input in [
        "hht([[1,0,-1,0,1,0,-1,0],[1,1,0,0,1,1,0,0]],output=matrix)".to_owned(),
        "hht([[1,0,-1,0]],output=matrix,samplerate=8,frequencies=0..4)".to_owned(),
        format!("hht([{p}],output=matrix)"),
        format!("hht([{p},{q}],output=matrix,frequencies=0.5..2.5)"),
        format!("hht([{p}],output=matrix,samplerate=100,ln)"),
        format!("hht([{p}],output=matrix,samplerate=100,ln=true,frequencies=10..40)"),
        "hht([0.1,0.9,-0.4,0.6,-0.8,0.3,0.7,-0.9,0.2,0.5,-0.6,0.8,-0.1,0.4,-0.7,0.2],output=matrix)"
            .to_owned(),
        "hht([[1,-1,1,-1]],output=matrix)".to_owned(),
        "hht([[1,-1,1,-1]],output=matrix,frequencies=0..4)".to_owned(),
        "hht([[1,0,-1,0]],output=matrix,samplerate=8,frequencies=0..1.99999)".to_owned(),
        "hht([[1,0,-1,0]],output=matrix,frequencies=a..b)".to_owned(),
        "hht([[1,0,-1,0]],output=matrix,frequencies=2..2)".to_owned(),
    ] {
        match (giac.value(&input), ours(&input)) {
            (Ok(theirs), Ok(mine)) if near(&ctx, &theirs, &mine) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's errors of `stft`, `istft` and `hht` (p6.out, p7.out): ours stays
/// with giac's message as the warning.
#[test]
fn transform_errors_stay() {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, message) in [
        (
            "stft([1,2,3,4,5,6,7,8],3)",
            "Window size must be a power of 2",
        ),
        (
            "stft([1,2,3,4,5,6,7,8],1)",
            "Window size should be at least 2",
        ),
        ("stft([1,2,3,4,5,6,7,8],2,4)", "Expected a window function"),
        ("stft([1,2,3,4,5,6,7,8],0)", "Expected a positive integer"),
        ("stft([1,2,3,4,5,6,7,8],2.0)", "Expected a positive integer"),
        (
            "stft([1,2,3,4,5,6,7,8],2,hann_window,1)",
            "Too many input arguments",
        ),
        ("stft([x,1,2,3],2)", "Bad Argument Value"),
        ("stft([])", "Invalid dimension"),
        ("stft(3)", "Bad Argument Type"),
        ("istft([[1,2],[3]])", "Bad Argument Type"),
        ("istft([[1,x],[3,4]])", "Bad Argument Value"),
        ("istft([1,2])", "Bad Argument Type"),
        ("hht(3)", "Bad Argument Type"),
        ("hht([1,2,3],output=matrix)", "No IMF found"),
        ("hht([[1,0,-1,0]],output=foo)", "Unknown output type"),
        (
            "hht([[1,0,-1,0]],ln,output=matrix)",
            "Log-frequency requires sample rate",
        ),
        (
            "hht([[1,0,-1,0]],output=matrix,samplerate=0)",
            "Sample rate must be a positive integer",
        ),
        (
            "hht([[1,0,-1,0]],output=matrix,ln=3,samplerate=8)",
            "Expected a boolean value",
        ),
        (
            "hht([[1,0,-1,0]],output=matrix,foo=3)",
            "Unknown optional argument",
        ),
    ] {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if !(printed.starts_with("GIAC_ERROR") && printed.contains(message)) {
                wrong.push(format!(
                    "{input}: giac prints {printed}, recorded {message}"
                ));
            }
        }
        let stays = evaluated(input).is_some_and(|(e, found, warnings)| {
            found.head() == e.head()
                && warnings
                    .iter()
                    .any(|w| matches!(w, Warning::Native { message: m, .. } if m == message))
        });
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {message}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures (p6.out, p7.out): giac's `stft` applies any function as the
/// window (ours takes the window functions); giac draws `hht`'s plot and
/// image (the default) and reads MIDI note ranges, which stay in ours; and
/// giac's `bins=n` gives 161 rows (it stores the option's tag), ours `n`.
#[test]
fn transform_departures() {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, recorded) in [
        ("stft([1,2,3,4,5,6,7,8],2,sin)", "[[-18.0376956654,"),
        ("hht([[1,0,-1,0]])", "Graphic_object"),
        ("hht([[1,2],[3,4]],output=plot)", "labels="),
        (
            "hht([[1,0,-1,0]],output=matrix,frequencies=[\"C4\",\"C5\"],samplerate=8)",
            "[[0,0,0,0],[0,0,0,0],",
        ),
    ] {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if !collapsed(&printed).starts_with(recorded) {
                wrong.push(format!(
                    "{input}: giac prints {printed}, recorded {recorded}…"
                ));
            }
        }
        let stays = evaluated(input)
            .is_some_and(|(e, found, warnings)| found.head() == e.head() && !warnings.is_empty());
        if !stays {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with a warning",
                evaluated(input)
            ));
        }
    }
    let bins = "hht([[1,0,-1,0,1,0,-1,0],[1,1,0,0,1,1,0,0]],output=matrix,bins=4)";
    if let Some(giac) = giac {
        let printed = giac.eval(&format!("size({bins})")).unwrap_or_default();
        if printed != "161" {
            wrong.push(format!("size({bins}): giac prints {printed}, recorded 161"));
        }
    }
    match ours(bins) {
        Ok(m) if m.args().len() == 4 && m.args().iter().all(|r| r.args().len() == 8) => {}
        found => wrong.push(format!("{bins}: ours {found:?}, wanted 4 rows of 8")),
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
