//! Xcas's sound clips against giac (CAS-631): `createwav`, `soundsec`,
//! `samplerate`, `bit_depth`, `channels`, `duration`, `channel_data` and
//! `about` (giac's answers in .run/1.7/giac/probe5.out, probe9.out,
//! probe10.out, probe11.out and .run/1.7/T7/dur.out; the prototype
//! .run/1.7/T7/proto.py, logs/T7.c1.md ## Rows). Exact answers against
//! giac's decimals within their error (C0006, C0032); several channels as a
//! list (C0023). giac's errors stay with its message as the warning (C0022,
//! C0056). Departures, giac's answer as caseval prints it (C0047) and ours:
//! a clip prints as the call `soundclip(rate, depth, [codes]…)`, an
//! ordinary call to `+`, `*`, `abs`, `[ ]`, `size`, `evalf` and `concat`;
//! giac keeps 47 fraction bits in its doubles (C0088). Session rows: the
//! session prints a clip as that call, which reads back, and `about(x)`
//! holds its facts (C0026). The clip edits (CAS-632): `set_channel_data`,
//! `stereo2mono`, `mixdown` and `splice` (.run/1.7/T8/probe1.out and the
//! prototype .run/1.7/T8/proto.py, logs/T8.c1.md).

use cas_core::{Context, Expr, Symbol, Warning};
use cas_core::{EvalOptions, View};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, Session, decimal_balls, lower, parse};

/// `createwav([0.5,-0.5,0.25,0])` as ours prints it.
const A: &str = "soundclip(44100, 16, [29204, -29204, 14602, 0])";
/// `createwav([0.5,-0.5])` and `createwav([1/2,-1/2])` as ours prints it.
const B: &str = "soundclip(44100,16,[29204,-29204])";

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
/// error ([`near`]), ours exact (R67).
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

/// Clips, their properties and channels, `soundsec` (probe5, probe9,
/// probe10, dur.out): giac's answers read back. Several channels are a
/// sequence in both: bracketed (C0023).
#[test]
#[allow(clippy::too_many_lines, reason = "one row per probe")]
fn values() {
    let s = "createwav([[0.5,-0.5,0.25,0],[0.1,0.2,0.3,0.4]])";
    let a = "createwav([0.5,-0.5,0.25,0])";
    let b = "createwav([0.1,0.2,0.3,0.4])";
    let z = "createwav([0,0,0,0])";
    let eight = "createwav([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])";
    let six = "createwav([0.5,-0.5,0.25,0,1,0.75]";
    let rows = [
        format!("channels({a})"),
        format!("bit_depth({a})"),
        format!("samplerate({a})"),
        format!("duration({a})"),
        format!("channel_data({a})"),
        format!("[channel_data({s})]"),
        format!("channel_data({s},matrix)"),
        format!("channel_data({s},right)"),
        format!("channel_data({s},left)"),
        format!("channel_data({six},samplerate=2),range=1..2)"),
        format!("channel_data({six},samplerate=2),range=1.0..2.5)"),
        format!("channel_data({six}),range=[1,3])"),
        format!("channels({s})"),
        "channel_data(createwav([1/2,-1/2,1/3]))".to_owned(),
        "channel_data(createwav([0.5,-0.5,0.25,0],bit_depth=8))".to_owned(),
        "channel_data(createwav([0.5,-0.5,0.25,0],bit_depth=24))".to_owned(),
        "channel_data(createwav([0.123456789,-0.5]))".to_owned(),
        "channel_data(createwav([2,-1,0.5]))".to_owned(),
        "channel_data(createwav([2,-1,0.5],normalize=-3))".to_owned(),
        "channel_data(createwav([0.2,-0.1,0.05],normalize=0))".to_owned(),
        "samplerate(createwav([0.1],samplerate=48000))".to_owned(),
        "duration(createwav([0.1,0.2,0.3],samplerate=2))".to_owned(),
        "[channel_data(createwav([0.1,0.2],channels=2))]".to_owned(),
        "size(soundsec(1))".to_owned(),
        "soundsec(1)[0..3]".to_owned(),
        "size(soundsec(1,22100))".to_owned(),
        "soundsec(1,22100)[0..3]".to_owned(),
        "soundsec(1/1000,8000)".to_owned(),
        "soundsec(0.5,4)".to_owned(),
        "soundsec(2,3)".to_owned(),
        "channel_data(createwav([0.5,-0.5],normalize=-1))".to_owned(),
        "channel_data(createwav([0.5,-0.5],normalize=1))".to_owned(),
        "channel_data(createwav([0.5,-0.5],normalize=-100))".to_owned(),
        "channel_data(createwav([0.2,-0.1],normalize=0))".to_owned(),
        "channel_data(createwav([1.5,-0.1],normalize=0))".to_owned(),
        "channel_data(createwav([0,0,0]))".to_owned(),
        "[channel_data(createwav([0.5,0.25],channels=3))]".to_owned(),
        "channel_data(createwav([0.5,0.25],duration=0.0001))".to_owned(),
        "channel_data(createwav(sin(2*pi*440*soundsec(1/1000,8000))))".to_owned(),
        "channel_data(createwav([0.5,0.25]),0)".to_owned(),
        "channel_data(createwav([[0.5,0.25],[0.1,0.2]]),1)".to_owned(),
        "channel_data(createwav([[0.5,0.25],[0.1,0.2]]),0)".to_owned(),
        format!("channel_data({six}),range=[4,10])"),
        format!("channel_data({six},samplerate=2),range=0.5..1.0)"),
        "size(channel_data(createwav(duration=0.5,channels=1,samplerate=4)))".to_owned(),
        "size(channel_data(createwav(duration=1,channels=1,samplerate=10)))".to_owned(),
        "size(channel_data(createwav(duration=0.25,channels=1,samplerate=10)))".to_owned(),
        "size(channel_data(createwav(length=3,duration=1,channels=1,samplerate=10)))".to_owned(),
        "size(channel_data(createwav([1,2],duration=1,samplerate=10)))".to_owned(),
        // The clip edits (CAS-632; probe5, probe9, .run/1.7/T8/probe1.out).
        format!("channel_data(set_channel_data({a},[0.9,0.8]))"),
        format!("channel_data(set_channel_data({a},[0.9,0.8],1))"),
        format!("channel_data(set_channel_data({a},[0.9,0.8,0.7,0.6,0.5],2))"),
        format!("channel_data(set_channel_data({a},[0.9,0.8],3))"),
        format!("channel_data(set_channel_data({a},[0.9,0.8],0,left))"),
        format!("channel_data(set_channel_data({a},[[0.9,0.8]],1))"),
        format!("channel_data(set_channel_data({s},[0.9],2,right),matrix)"),
        format!("channel_data(set_channel_data({s},[[0.9,0.8],[0.1,0.2]],1),matrix)"),
        format!("channel_data(set_channel_data({s},[0.9,0.8]),matrix)"),
        format!("channel_data(set_channel_data({a},[0.9,x],3))"),
        format!("channel_data(set_channel_data({a},[sqrt(2)/2,1/3,-2,pi]))"),
        format!("channel_data(set_channel_data({a},[]))"),
        format!("channel_data(set_channel_data({a},[0.9],0,0,7))"),
        format!("channel_data(set_channel_data({a},[1,1]))"),
        "channel_data(set_channel_data(createwav([0.5,-0.5],bit_depth=8),[0.9,-1]))".to_owned(),
        "channel_data(set_channel_data(createwav([0.5,-0.5],bit_depth=24),[0.9,-1]))".to_owned(),
        format!("channel_data(stereo2mono({s}))"),
        format!("channels(stereo2mono({s}))"),
        "channel_data(stereo2mono(createwav([0.5,0.25])))".to_owned(),
        "channel_data(stereo2mono(createwav([0.5,-0.25],bit_depth=8)))".to_owned(),
        "channel_data(stereo2mono(createwav([[0.5,0.25],[0.1,0.2],[0.3,-0.4]])))".to_owned(),
        "channel_data(stereo2mono(createwav([[0.7,0.25],[0.3,-0.2]],bit_depth=24)))".to_owned(),
        format!("channel_data(mixdown({a},1,{b},1))"),
        format!("channel_data(mixdown({a},[1,0],{b},[1,0]))"),
        format!("channel_data(mixdown({a},[1,0,0],{b},[1,0,0]))"),
        format!("channel_data(mixdown({a},0,{b},0))"),
        format!("channel_data(mixdown({a},2,{b},0))"),
        format!("channel_data(mixdown({a},[0,-6],{b},[0,-6]))"),
        format!("channel_data(mixdown({a},[0,0,-1],{b},[0,0,1]),matrix)"),
        format!("channels(mixdown({a},[0,0,-1],{b},[0,0,1]))"),
        format!("size(channel_data(mixdown({a},0.5,{b},0.5)))"),
        format!("size(channel_data(mixdown({a},0.0001,{b},0)))"),
        format!("size(channel_data(mixdown({a},0.0,{b},0)))"),
        format!("size(channel_data(mixdown({a},1/2,{b},0)))"),
        format!("size(channel_data(mixdown({a},sqrt(2),{b},0)))"),
        format!("size(channel_data(mixdown({a},0,{b},1e-5)))"),
        format!("channel_data(mixdown({a},0,{s},0),matrix)"),
        format!("channel_data(mixdown({a},[0,0,0.5],{s},0),matrix)"),
        format!("channel_data(mixdown({a},[0,0,-0.5],{b},[1,-3,0.5]),matrix)"),
        format!("channel_data(mixdown({s},[0,0,-0.5],{b},[1,-3,0.5]),matrix)"),
        format!("channel_data(mixdown({a},[0,0,0,2,2],{b},1))"),
        format!("channel_data(mixdown({a},[0,0,0,1,0],{b},[0,0,0,0,3]))"),
        format!("channel_data(mixdown(createwav([0.5,-0.5],bit_depth=8),0,{b},1))"),
        format!("bit_depth(mixdown(createwav([0.5,-0.5],bit_depth=8),0,{b},1))"),
        format!("channel_data(mixdown({a},0,{a},0,{a},0))"),
        "channel_data(mixdown(createwav([0.9,-0.9,0.9]),0,createwav([0.9,-0.9,0.1]),0))".to_owned(),
        format!("channel_data(mixdown({a},[],{b},0))"),
        format!("channel_data(mixdown({a},0,{b},0,{a},1))"),
        "channel_data(mixdown(createwav([0.5,0.25],bit_depth=24),0,createwav([0.1,0.3],bit_depth=24),0))".to_owned(),
        format!("channel_data(splice({a},{b}))"),
        format!("channel_data(splice({a},{b},2,1))"),
        format!("channel_data(splice({a},{b},2,0.5))"),
        format!("channel_data(splice({a},{b},2,2))"),
        format!("channel_data(splice({a},{b},2,1/2))"),
        format!("channel_data(splice({a},{b},2,0))"),
        format!("channel_data(splice({a},{b},2,-1))"),
        format!("channel_data(splice({a},{b},1))"),
        format!("channel_data(splice({a},{b},3))"),
        format!("channel_data(splice({a},{b},4))"),
        format!("channel_data(splice({a},{b},5))"),
        format!("channel_data(splice({a},{b},1.5))"),
        format!("channel_data(splice({a},{b},0.00004))"),
        format!("channel_data(splice({z},{b}))"),
        format!("channel_data(splice({z},{b},2,1))"),
        format!("channel_data(splice({a},{z}))"),
        format!("channel_data(splice({s},createwav([[0.1,0.2,0.3,0.4],[-0.5,-0.5,0.5,0.5]])),matrix)"),
        format!("channel_data(splice({s},createwav([[0.1,0.2,0.3,0.4],[0.5,0.5,0.5,0.5]])),matrix)"),
        format!("channel_data(splice({a},createwav([0.9])))"),
        format!("channel_data(splice({a},createwav([0.9]),3,1))"),
        format!("channel_data(splice({eight},createwav([0.9,-0.9]),6,1))"),
        format!("channel_data(splice({eight},createwav([0.9,-0.9]),7,0))"),
        "channel_data(splice(createwav([0.9,0.9,0.9,0.9]),createwav([0.9,0.9,0.9,0.9]),2,0))".to_owned(),
        "channel_data(splice(createwav([0.5,-0.5,0.25,0],bit_depth=24),createwav([0.1,0.2,0.3,0.4],bit_depth=24),2,1))".to_owned(),
        "channel_data(splice(createwav([0.5,-0.5,0.25,0],bit_depth=8),createwav([0.1,0.2,0.3,0.4],bit_depth=8)))".to_owned(),
    ];
    check_near(&rows.iter().map(String::as_str).collect::<Vec<_>>());
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

/// Whitespace removed: giac writes `[0.5,0.25]` in its messages, ours
/// `[0.5, 0.25]`.
fn bare(s: &str) -> String {
    s.split_whitespace().collect()
}

/// Whether ours of `input` stays, its head the input's, with a warning of
/// `expected` (whitespace aside).
fn stays_with(input: &str, expected: &str) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && warnings.iter().any(
                |w| matches!(w, Warning::Native { message, .. } if bare(message) == bare(expected)),
            )
    })
}

/// giac's errors that ours shares (probe5, probe9, probe10): the call
/// stays with giac's message as the warning (C0022, C0056).
#[test]
#[allow(clippy::too_many_lines, reason = "one row per error")]
fn errors_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let kind = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Type");
    let dim = |m: &str| format!("GIAC_ERROR: {m}. Error: Invalid dimension");
    let a = "createwav([0.5,-0.5,0.25,0])";
    let b = "createwav([0.1,0.2,0.3,0.4])";
    let slow = "createwav([0.1,0.2,0.3,0.4],samplerate=22050)";
    let same = "Clips should have the same bit depths, sample rates and channel counts";
    let s = "createwav([[0.5,-0.5,0.25,0],[0.1,0.2,0.3,0.4]])";
    let six = "createwav([0.5,-0.5,0.25,0,1,0.75])";
    let failed = "Failed to set channel data, [0.5,0.25] is not a valid sample value";
    let rate = "Invalid sample rate specification";
    let quality = "Invalid quality specification";
    let rows = [
        (
            format!("channel_data({s},2)"),
            value("Invalid channel number"),
        ),
        (
            "channel_data(createwav([0.5,-0.5,0.25,0],bit_depth=32))".to_owned(),
            value("Bit depth should be either 8, 16 or 24"),
        ),
        (
            "channel_data(createwav(duration=0.0001))".to_owned(),
            value("Number of channels not set"),
        ),
        (
            "channel_data(createwav(duration=0.0001,samplerate=20000))".to_owned(),
            value("Number of channels not set"),
        ),
        (
            "channel_data(createwav([x,1]))".to_owned(),
            value("Sample data should be numeric. Error: Bad Argument Value"),
        ),
        ("createwav([[0.5,0.25],[0.1]])".to_owned(), value(failed)),
        (
            "channel_data(createwav([[0.5,0.25],[0.1]]),matrix)".to_owned(),
            value(failed),
        ),
        (
            "createwav([0.5,0.25],samplerate=0)".to_owned(),
            value("Sample rate should be a positive integer"),
        ),
        (
            "createwav([0.5,0.25],samplerate=1.5)".to_owned(),
            value("Sample rate should be a positive integer"),
        ),
        (
            "createwav([],samplerate=8000)".to_owned(),
            value("Audio length not set"),
        ),
        (
            "samplerate([1,2])".to_owned(),
            "GIAC_ERROR: Invalid sound data. Error: Bad Argument Type".to_owned(),
        ),
        (
            "channel_data(createwav([0.5,0.25]),1)".to_owned(),
            value("Invalid channel number"),
        ),
        (
            format!("channel_data({s},right,range=1..2)"),
            value("Invalid slice bounds"),
        ),
        (
            format!("channel_data({six},range=4..10)"),
            value("Invalid slice bounds"),
        ),
        (
            "createwav(createwav([0.5,-0.5]))".to_owned(),
            value("Invalid argument 1"),
        ),
        (
            "channel_data(createwav(createwav([0.5,-0.5]),bit_depth=8))".to_owned(),
            value("Invalid argument 1"),
        ),
        // The clip edits (CAS-632; probe5, probe9, .run/1.7/T8/probe1.out).
        (
            format!("set_channel_data({a},[0.9,0.8],4)"),
            value("Invalid offset"),
        ),
        (
            format!("set_channel_data({a},[0.9,0.8],-1)"),
            value("Invalid offset"),
        ),
        (
            format!("set_channel_data({a},[0.9,0.8],1.0)"),
            kind("Expected an integer"),
        ),
        (
            format!("set_channel_data({a},[0.9,0.8],0,1)"),
            dim("Invalid channel number"),
        ),
        (
            format!("set_channel_data({a},[0.9,0.8],0,right)"),
            value("Invalid channel number"),
        ),
        (
            format!("set_channel_data({a},[0.9,0.8],0,x)"),
            kind("Invalid channel specification"),
        ),
        (
            format!("set_channel_data({a},0.9)"),
            kind("Invalid data specification"),
        ),
        (
            format!("set_channel_data({s},[[0.9,0.8],[0.1,0.2]],0,0)"),
            kind("Invalid data specification"),
        ),
        (
            format!("set_channel_data({a},[0.9,x])"),
            value("Failed to set channel data, x is not a valid sample value"),
        ),
        (
            format!("set_channel_data({a},[[0.9],[0.8,0.7]])"),
            value("Failed to set channel data, [0.9] is not a valid sample value"),
        ),
        (
            "set_channel_data([1,2],[0.9])".to_owned(),
            value("Expected audio clip or image"),
        ),
        (
            "stereo2mono([0.5,0.25])".to_owned(),
            kind("Invalid sound data"),
        ),
        (format!("stereo2mono({a},{b})"), kind("Invalid sound data")),
        ("stereo2mono(x)".to_owned(), kind("Invalid sound data")),
        (
            format!("mixdown({a},{b})"),
            value("Invalid offset specification"),
        ),
        (
            format!("mixdown({a},-6,{b},0)"),
            value("Invalid offset specification"),
        ),
        (
            format!("mixdown({a},x,{b},0)"),
            value("Invalid offset specification"),
        ),
        (
            format!("mixdown({a},[0.5,-1],{b},[0.5,1])"),
            value("Invalid gain specification"),
        ),
        (
            format!("mixdown({a},[0,0],{b},[0,2])"),
            value("Invalid gain specification"),
        ),
        (
            format!("mixdown({a},[0,x],{b},0)"),
            value("Invalid gain specification"),
        ),
        (
            format!("mixdown({a},[0,0,2],{b},1)"),
            value("Invalid pan specification"),
        ),
        (
            format!("mixdown({a},[0,0,0,-1],{b},1)"),
            value("Invalid fade-in specification"),
        ),
        (
            format!("mixdown({a},[0,0,0,0,x],{b},1)"),
            value("Invalid fade-out specification"),
        ),
        (
            format!("mixdown({a},[0,0,0,3,2],{b},1)"),
            dim("Fade-in/out is too long"),
        ),
        (
            format!("mixdown({a},0,[1,2],0)"),
            value("Invalid argument 3"),
        ),
        (
            format!("mixdown({a},0,{slow},0)"),
            value("All clips must have the same sample rate"),
        ),
        (
            format!("mixdown(createwav([0.1,0.2],channels=3),0,{b},0)"),
            value("Only mono and stereo clips are allowed"),
        ),
        (
            format!("splice({a},{b},0)"),
            value("Invalid splice position"),
        ),
        (
            format!("splice({a},{b},-3)"),
            value("Invalid splice position"),
        ),
        (format!("splice({a},{slow})"), value(same)),
        (format!("splice({s},{b},2)"), value(same)),
        (
            format!("splice({a},{b},x)"),
            kind("Crossfade length should be a real number"),
        ),
        (
            format!("splice({a},{b},2,x)"),
            value("Crossfade exponent should be a real number"),
        ),
        (
            format!("splice({a},[1,2])"),
            kind("First two arguments should be audio clips"),
        ),
        // resample (CAS-633, .run/1.7/T8/probe_rs.out).
        (format!("resample({a},0)"), value(rate)),
        (format!("resample({a},-3)"), value(rate)),
        (format!("resample({a},22050.0)"), value(rate)),
        (format!("resample({a},1/2)"), value(rate)),
        (format!("resample({a},22050,5)"), value(quality)),
        (format!("resample({a},22050,-1)"), value(quality)),
        (format!("resample({a},22050,2.0)"), value(quality)),
        ("resample(x)".to_owned(), kind("Invalid sound data")),
        ("resample([1,2,3])".to_owned(), kind("Invalid sound data")),
        ("resample(x,22050)".to_owned(), kind("Invalid sound data")),
        ("resample(x,0)".to_owned(), value(rate)),
        ("resample(x,22050,7)".to_owned(), value(quality)),
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
        if !stays_with(input, message) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {message:?}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whitespace collapsed (C0047).
fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Departures, giac's answer as caseval prints it (probe5, probe9, probe10,
/// dur.out) and ours: giac prints a clip as a sentence, ours as the call
/// `soundclip(…)`, which reads back (`string` and `about` too); giac's `type`
/// of a clip is 16, ours `DOM_SYMBOLIC`; giac's clip arithmetic (`+` splices,
/// `-` negates, `A - A` a clip of 4 samples, `*` and `abs` refused) where
/// ours sees an ordinary call: `+`, `-`, `*`, `abs`, `evalf` and `concat` stay, `A - A`
/// is 0; giac's `A[k]` is channel `k` and `size(A)` the length, ours `A[k]`
/// stays and `size(A)` is the call's 2 + channels; giac refuses `size = n`,
/// ours takes it as the length; giac's 47-bit doubles make 8 samples of
/// `duration = 0.001` (or `1/1000`) at 8000 Hz and 16 of 0.002, ours exactly 9
/// and 17 (C0088); giac's `resample` needs libsamplerate, which its Debian
/// build lacks, ours resamples (CAS-633).
#[test]
#[allow(clippy::too_many_lines, reason = "one row per departure")]
fn departures() {
    let s = "soundclip(44100,16,[29204,-29204,14602,0],[5841,11681,17522,23363])";
    let one = "soundclip(44100,16,[29204])";
    let mono = |n: usize| format!("a sound clip with {n} samples at 44100 Hz (16 bit, mono)");
    let stereo = "a sound clip with 2 samples at 44100 Hz (16 bit, stereo)";
    let refused = |m: &str| format!("GIAC_ERROR: {m} not redefined Error: Bad Argument Value");
    let argument = "GIAC_ERROR: Invalid argument 1. Error: Bad Argument Value";
    let spliced =
        "channel_data(splice(createwav([0.5,-0.5,0.25,0]),createwav([0.1,0.2,0.3,0.4],samplerate=";
    let no_lsr =
        "GIAC_ERROR: libsamplerate is required for resampling audio. Error: Bad Argument Value";
    let cases = [
        (
            "createwav([0.5,-0.5,0.25,0])".to_owned(),
            mono(4),
            A.to_owned(),
        ),
        (
            "string(createwav([0.5,-0.5,0.25,0]))".to_owned(),
            format!("\"{}\"", mono(4)),
            format!("\"{A}\""),
        ),
        (
            "about(createwav([0.5,-0.5,0.25,0]))".to_owned(),
            mono(4),
            A.to_owned(),
        ),
        ("createwav([1/2,-1/2])".to_owned(), mono(2), B.to_owned()),
        (
            "createwav([0.1,0.2],channels=2)".to_owned(),
            stereo.to_owned(),
            "soundclip(44100,16,[14602,29204],[14602,29204])".to_owned(),
        ),
        (
            "createwav([0.5,0.25],channels=3)".to_owned(),
            stereo.to_owned(),
            "soundclip(44100,16,[29204,14602],[29204,14602],[29204,14602])".to_owned(),
        ),
        (
            "about(createwav([0.5,0.25],channels=2))".to_owned(),
            stereo.to_owned(),
            "soundclip(44100,16,[29204,14602],[29204,14602])".to_owned(),
        ),
        (
            "type(createwav([0.5,-0.5,0.25,0]))".to_owned(),
            "16".to_owned(),
            "DOM_SYMBOLIC".to_owned(),
        ),
        (
            "createwav([0.5,-0.5])+createwav([0.25])".to_owned(),
            mono(3),
            format!("{B}+{one}"),
        ),
        (
            "-createwav([0.5,-0.5])".to_owned(),
            mono(2),
            format!("-{B}"),
        ),
        (
            "channel_data(2*createwav([0.25,-0.25]))".to_owned(),
            refused("Binary *"),
            format!("channel_data(2*{B})"),
        ),
        (
            "createwav([0.5,-0.5])-createwav([0.5,-0.5])".to_owned(),
            mono(4),
            "0".to_owned(),
        ),
        (
            "abs(createwav([0.5,-0.5]))".to_owned(),
            refused("Abs"),
            format!("abs({B})"),
        ),
        (
            "createwav([0.5,-0.5,0.25,0])[0]".to_owned(),
            "[0.891235351562,-0.891235351562,0.445617675781,0.0]".to_owned(),
            format!("{A}[0]"),
        ),
        (
            "createwav([[0.5,-0.5,0.25,0],[0.1,0.2,0.3,0.4]])[1]".to_owned(),
            "[0.178253173828,0.356475830078,0.534729003906,0.712982177734]".to_owned(),
            format!("{s}[1]"),
        ),
        (
            "size(createwav([0.5,-0.5,0.25,0]))".to_owned(),
            "4".to_owned(),
            "3".to_owned(),
        ),
        (
            "evalf(createwav([0.5,-0.5]))".to_owned(),
            mono(2),
            format!("evalf({B})"),
        ),
        (
            "concat(createwav([0.5,-0.5]),createwav([0.25]))".to_owned(),
            "0".to_owned(),
            format!("concat({B},{one})"),
        ),
        (
            "channel_data(createwav(size=3,channels=1))".to_owned(),
            argument.to_owned(),
            "[0,0,0]".to_owned(),
        ),
        (
            "channel_data(createwav(size=3))".to_owned(),
            argument.to_owned(),
            "channel_data(createwav(size=3))".to_owned(),
        ),
        (
            "createwav(duration=0.001,channels=1,samplerate=8000)".to_owned(),
            "a sound clip with 8 samples at 8000 Hz (16 bit, mono)".to_owned(),
            "soundclip(8000,16,[0,0,0,0,0,0,0,0,0])".to_owned(),
        ),
        (
            "channel_data(createwav(duration=0.001,channels=1,samplerate=8000))".to_owned(),
            "[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]".to_owned(),
            "[0,0,0,0,0,0,0,0,0]".to_owned(),
        ),
        (
            "size(channel_data(createwav(duration=0.002,channels=1,samplerate=8000)))".to_owned(),
            "16".to_owned(),
            "17".to_owned(),
        ),
        (
            "size(channel_data(createwav(duration=1/1000,channels=1,samplerate=8000)))".to_owned(),
            "8".to_owned(),
            "9".to_owned(),
        ),
        // splice: giac checks the clips' parameters only for the similarity
        // (a negative or absent exponent); ours always (CAS-632).
        (
            format!("{spliced}22050),2,1))"),
            "[0.891235351562,-0.891235351562,0.334197998047,0.445617675781,0.668426513672,0.891235351562]"
                .to_owned(),
            format!("{spliced}22050),2,1))"),
        ),
        (
            format!("{spliced}22050),2,0))"),
            "[0.891235351562,-0.891235351562,0.668395996094,0.445617675781,0.668426513672,0.891235351562]"
                .to_owned(),
            format!("{spliced}22050),2,0))"),
        ),
        // resample: giac's Debian build has no libsamplerate (CAS-633,
        // .run/1.7/T8/probe_rs.out); ours resamples (the prototype
        // .run/1.7/T8/resample.py), a clip of no samples stays.
        (
            "resample(createwav([0.5,-0.5,0.25,0]))".to_owned(),
            no_lsr.to_owned(),
            "soundclip(44100,16,[22068,-21986,7576,6577])".to_owned(),
        ),
        (
            "resample(createwav([0.5,-0.5,0.25,0]),22050)".to_owned(),
            no_lsr.to_owned(),
            "soundclip(22050,16,[4683,-1172])".to_owned(),
        ),
        (
            "resample(createwav([0.5,-0.5,0.25,0]),88200,0)".to_owned(),
            no_lsr.to_owned(),
            "soundclip(88200,16,[29203,29203,29203,0,-29203,-7301,14602,7301])".to_owned(),
        ),
        (
            "resample(createwav([0.5,-0.5,0.25,0]),88200,1)".to_owned(),
            no_lsr.to_owned(),
            "soundclip(88200,16,[29203,29203,29203,29203,-29203,-29203,14602,14602])".to_owned(),
        ),
        (
            "resample(createwav([0.5,-0.5,0.25,0]),48000,4)".to_owned(),
            no_lsr.to_owned(),
            "soundclip(48000,16,[27744,-26007,4832,11250])".to_owned(),
        ),
        (
            "resample(createwav([0.5,-0.5,0.25,0]),1)".to_owned(),
            no_lsr.to_owned(),
            "resample(soundclip(44100,16,[29204,-29204,14602,0]),1)".to_owned(),
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, mine) in &cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours(input);
        if found.as_ref().ok() != ours(mine).ok().as_ref() {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Clipping (CAS-632): ours warns once per call when any sample clips;
/// giac's `mix` and `splice` keep only the last sample's status
/// (signalprocessing.cc :1185, :1213), so `mixdown(A,0,B,0)`, which clips
/// its first and third samples (its value in [`values`]), warns in ours
/// and not in giac: a departure.
#[test]
fn clipping() {
    let a = "createwav([0.5,-0.5,0.25,0])";
    let b = "createwav([0.1,0.2,0.3,0.4])";
    let nine = "createwav([0.9,0.9,0.9,0.9])";
    let rows = [
        (format!("mixdown({a},0,{b},0)"), 1),
        (format!("mixdown({a},1,{b},1)"), 1),
        (format!("mixdown({a},[0,-6],{b},[0,-6])"), 0),
        (format!("splice({nine},{nine},2,0)"), 1),
        (format!("splice({a},{b})"), 0),
    ];
    let mut wrong = Vec::new();
    for (input, n) in &rows {
        let found = evaluated(input).map(|(_, _, warnings)| {
            warnings
                .iter()
                .filter(|w| matches!(w, Warning::Native { message, .. } if &**message == "clipping audio, the result will be distorted"))
                .count()
        });
        if found != Some(*n) {
            wrong.push(format!("{input}: {found:?} clipping warnings, wanted {n}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The results a session prints for `input`, statement by statement.
fn session_texts(input: &str) -> Vec<String> {
    Session::new()
        .run(input)
        .iter()
        .filter_map(|output| output.result.as_ref().map(|r| r.text().to_owned()))
        .collect()
}

/// The session (CAS-631): a clip prints as the call, which reads back to
/// the clip, and is its `string`; `about(x)` after `assume(x>0)` holds its
/// facts (C0026), where giac answers its assumption record (probe11, the
/// `assume` and the `purge` in one input: giac's session is shared).
#[test]
fn session() {
    let clip = "createwav([0.5,-0.5,0.25,0])";
    let mut wrong = Vec::new();
    let shown = session_texts(&format!("A := {clip}\nA"));
    if shown.last().map(String::as_str) != Some(A) {
        wrong.push(format!("A: ours {shown:?}"));
    }
    if ours(A).ok() != ours(clip).ok() {
        wrong.push(format!("{A} reads back as {:?}", ours(A)));
    }
    let string = session_texts(&format!("string({clip})"));
    let text = string.last().map(|t| ours(t).ok());
    if text != Some(ours(&format!("\"{A}\"")).ok()) {
        wrong.push(format!("string(A): ours {string:?}"));
    }
    let about = session_texts("assume(x > 0)\nabout(x)");
    if about.last().map(String::as_str) != Some("hold([x > 0])") {
        wrong.push(format!("about(x): ours {about:?}"));
    }
    if let Some(giac) = giac() {
        let input = "[assume(x>0),about(x),purge(x)][1]";
        let printed = giac.eval(input).unwrap_or_default();
        if collapsed(&printed) != "assume[[],[line[0,+infinity]],[0]]" {
            wrong.push(format!("{input}: giac prints {printed}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The codes of channel 0 of the clip `text` evaluates to.
fn codes(text: &str) -> Vec<i64> {
    let (_, found, _) = evaluated(text).unwrap_or_else(|| panic!("{text}: unreadable"));
    found
        .args()
        .get(2)
        .map(|c| c.args().iter().filter_map(Expr::as_i64).collect())
        .unwrap_or_default()
}

/// `soundclip(rate, depth, [codes])`.
fn mono_clip(rate: u32, depth: u32, codes: &[i64]) -> String {
    let codes: Vec<String> = codes.iter().map(i64::to_string).collect();
    format!("soundclip({rate},{depth},[{}])", codes.join(","))
}

/// `floor(range·x/2 + 1/2)` at 24 bits, `x = sin(2π·f·k/rate)`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "|range·x/2 + 1/2| < 2^23: the floor fits; k below 2^53"
)]
fn tone(f: f64, rate: f64, k: usize) -> i64 {
    let range = f64::from((1_u32 << 23) - 1);
    (range * 0.5 * (std::f64::consts::TAU * f * k as f64 / rate).sin() + 0.5).floor() as i64
}

/// Resampling's properties, ours only (CAS-633; the expected shapes and
/// tolerances from the prototype .run/1.7/T8/resample.py, resample.log,
/// with margin): `floor(len·r/sr)` samples; a constant clip stays constant
/// where the kernel lies inside it; at 48000 → 44100 Hz a 440 Hz tone of
/// amplitude 1/2 at 24 bits keeps its frequency and amplitude there
/// (measured 14, 1.8 and 1.4 codes off at the qualities 2, 3 and 4) and a
/// 23500 Hz tone, above the new Nyquist frequency, is attenuated to the
/// quality's stopband (measured at most 41, 3 and 1 codes).
#[test]
#[allow(
    clippy::cast_precision_loss,
    reason = "sample counts below 2^53: exact doubles"
)]
fn resample_properties() {
    let mut wrong = Vec::new();
    let zeros = |n: usize| vec![0_i64; n];
    for (len, sr, rate, n) in [
        (4, 44100, 88200, 8),
        (4, 44100, 22050, 2),
        (4, 44100, 30000, 2),
        (10, 48000, 44100, 9),
        (1000, 48000, 44100, 918),
        (3, 8000, 44100, 16),
        (1, 44100, 48000, 1),
    ] {
        for q in 0..=4 {
            let input = format!("resample({},{rate},{q})", mono_clip(sr, 16, &zeros(len)));
            let found = codes(&input).len();
            if found != n {
                wrong.push(format!(
                    "{len} samples {sr} -> {rate} q{q}: {found}, wanted {n}"
                ));
            }
        }
    }
    // (quality, zero crossings / cutoff, 440 Hz error, 23500 Hz peak) in codes.
    let qualities = [
        (2, 37.0 / 0.90, 20.0, 60),
        (3, 86.0 / 0.95, 4.0, 6),
        (4, 312.0 / 0.98, 3.0, 3),
    ];
    let rho = 44100.0 / 48000.0;
    let interior = |len: usize, reach: f64, k: usize| {
        let t = k as f64 / rho;
        t - reach / rho >= 0.0 && t + reach / rho <= (len - 1) as f64
    };
    let flat = mono_clip(48000, 16, &[10000; 900]);
    let low: Vec<i64> = (0..2048).map(|k| tone(440.0, 48000.0, k)).collect();
    let high: Vec<i64> = (0..2048).map(|k| tone(23500.0, 48000.0, k)).collect();
    let (low, high) = (mono_clip(48000, 24, &low), mono_clip(48000, 24, &high));
    for (q, reach, error, peak) in qualities {
        let out = codes(&format!("resample({flat},44100,{q})"));
        if out.len() != 826 || (0..826).any(|k| interior(900, reach, k) && out[k] != 10000) {
            wrong.push(format!("constant q{q}: {out:?}"));
        }
        let out = codes(&format!("resample({low},44100,{q})"));
        let off = (0..out.len())
            .filter(|&k| interior(2048, reach, k))
            .map(|k| (out[k] - tone(440.0, 44100.0, k)).abs())
            .max();
        if out.len() != 1881 || off.is_none_or(|d| d as f64 > error) {
            wrong.push(format!(
                "440 Hz q{q}: {} samples, {off:?} codes off",
                out.len()
            ));
        }
        let out = codes(&format!("resample({high},44100,{q})"));
        let top = (0..out.len())
            .filter(|&k| interior(2048, reach, k))
            .map(|k| out[k].abs())
            .max();
        if out.len() != 1881 || top.is_none_or(|t| t > peak) {
            wrong.push(format!(
                "23500 Hz q{q}: {} samples, peak {top:?}",
                out.len()
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
