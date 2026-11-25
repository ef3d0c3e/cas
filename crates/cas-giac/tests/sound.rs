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
//! holds its facts (C0026).

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
fn values() {
    let s = "createwav([[0.5,-0.5,0.25,0],[0.1,0.2,0.3,0.4]])";
    let a = "createwav([0.5,-0.5,0.25,0])";
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
fn errors_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let s = "createwav([[0.5,-0.5,0.25,0],[0.1,0.2,0.3,0.4]])";
    let six = "createwav([0.5,-0.5,0.25,0,1,0.75])";
    let failed = "Failed to set channel data, [0.5,0.25] is not a valid sample value";
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
/// and 17 (C0088).
#[test]
#[allow(clippy::too_many_lines, reason = "one row per departure")]
fn departures() {
    let s = "soundclip(44100,16,[29204,-29204,14602,0],[5841,11681,17522,23363])";
    let one = "soundclip(44100,16,[29204])";
    let mono = |n: usize| format!("a sound clip with {n} samples at 44100 Hz (16 bit, mono)");
    let stereo = "a sound clip with 2 samples at 44100 Hz (16 bit, stereo)";
    let refused = |m: &str| format!("GIAC_ERROR: {m} not redefined Error: Bad Argument Value");
    let argument = "GIAC_ERROR: Invalid argument 1. Error: Bad Argument Value";
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
