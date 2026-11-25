//! Tests of the sound clips (CAS-631): the rows of the prototype checked
//! against giac (crates/cas-core/proto/sound/clip/proto.py, logs/T7.c1.md ## Rows) that these
//! natives answer, their departures, giac's errors with their warnings,
//! `==`, `About`, and the limits.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr};

/// The double nearest `x`, as a hexadecimal float full form reads.
fn hex(x: f64) -> String {
    if x == 0.0 {
        return "0".to_owned();
    }
    let bits = x.to_bits();
    let exp = i64::try_from((bits >> 52) & 0x7ff).unwrap() - 1023;
    let frac = bits & ((1 << 52) - 1);
    format!("0x1.{frac:013x}p{exp}")
}

/// `text` with its decimals (`0.5`) as the language reads them: the
/// nearest doubles, which full form writes in hexadecimal.
fn decimals(text: &str) -> String {
    let b = text.as_bytes();
    let mut out = String::new();
    let mut i = 0;
    while i < b.len() {
        if b[i].is_ascii_digit() && (i == 0 || !b[i - 1].is_ascii_alphanumeric()) {
            let start = i;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
                i += 1;
            }
            let token = &text[start..i];
            if token.contains('.') {
                out.push_str(&hex(token.parse().unwrap()));
            } else {
                out.push_str(token);
            }
        } else {
            out.push(char::from(b[i]));
            i += 1;
        }
    }
    out
}

fn read(text: &str) -> Expr {
    let text = decimals(text);
    Expr::parse_full_form(&text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(ctx: &Context, cases: &[(String, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = eval(ctx, input);
        let expected = eval(ctx, expected);
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `input` with its arguments evaluated: the call as it stays.
fn stayed(ctx: &Context, input: &str) -> Expr {
    read(input).map_args(|a| ctx.eval(a).unwrap_or_else(|e| panic!("{input}: {e}")))
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap(), stayed(ctx, input), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

/// `CreateWav[args]`.
fn w(args: &str) -> String {
    format!("CreateWav[{args}]")
}

/// `ChannelData[CreateWav[args], options]`.
fn cd(args: &str, options: &str) -> String {
    if options.is_empty() {
        format!("ChannelData[{}]", w(args))
    } else {
        format!("ChannelData[{}, {options}]", w(args))
    }
}

const A: &str = "List[0.5, -0.5, 0.25, 0]";
const M: &str = "List[List[0.5, -0.5, 0.25, 0], List[0.1, 0.2, 0.3, 0.4]]";
const S6: &str = "List[0.5, -0.5, 0.25, 0, 1, 0.75]";
const A0: &str = "List[7301/8192, -7301/8192, 7301/16384, 0]";
const M1: &str = "List[5841/32768, 11681/32768, 8761/16384, 23363/32768]";

/// probe5's rows on one clip, its channels and slices.
#[test]
fn clips_and_channels() {
    let ctx = Context::new();
    let s6 = format!("{S6}, Equation[SampleRate, 2]");
    let both = format!("Sequence[{A0}, {M1}]");
    let rows = format!("List[{A0}, {M1}]");
    let cases: &[(String, &str)] = &[
        (format!("Channels[{}]", w(A)), "1"),
        (format!("BitDepth[{}]", w(A)), "16"),
        (format!("SampleRate[{}]", w(A)), "44100"),
        (format!("Duration[{}]", w(A)), "1/11025"),
        (cd(A, ""), A0),
        (cd(M, ""), &both),
        (cd(M, "Matrix"), &rows),
        (cd(M, "Right"), M1),
        (cd(M, "Left"), A0),
        (
            cd(&s6, "Equation[range, Range[1, 2]]"),
            "List[7301/32768, 0]",
        ),
        (
            cd(&s6, "Equation[range, Range[1.0, 2.5]]"),
            "List[7301/32768, 0, 7301/8192]",
        ),
        (
            cd(S6, "Equation[range, List[1, 3]]"),
            "List[-7301/16384, 7301/32768, 0]",
        ),
        (format!("Channels[{}]", w(M)), "2"),
        (
            w("List[1/2, -1/2]"),
            "SoundClip[44100, 16, List[29204, -29204]]",
        ),
        (
            cd("List[1/2, -1/2, 1/3]", ""),
            "List[7301/8192, -7301/8192, 19469/32768]",
        ),
        (
            cd(&format!("{A}, Equation[BitDepth, 8]"), ""),
            "List[113/128, -113/128, 57/128, 0]",
        ),
        (
            cd(&format!("{A}, Equation[BitDepth, 24]"), ""),
            "List[3738177/4194304, -3738177/4194304, 3738177/8388608, 0]",
        ),
        (
            cd("List[0.123456789, -0.5]", ""),
            "List[7211/32768, -7301/8192]",
        ),
        (
            cd("List[2, -1, 0.5]", ""),
            "List[7301/8192, -7301/16384, 7301/32768]",
        ),
        (
            cd("List[2, -1, 0.5], Equation[Normalize, -3]", ""),
            "List[23197/32768, -11599/32768, 5799/32768]",
        ),
        (
            cd("List[0.2, -0.1, 0.05], Equation[Normalize, 0]", ""),
            "List[6553/32768, -3277/32768, 819/16384]",
        ),
        (
            format!(
                "SampleRate[{}]",
                w("List[0.1], Equation[SampleRate, 48000]")
            ),
            "48000",
        ),
        (
            format!(
                "Duration[{}]",
                w("List[0.1, 0.2, 0.3], Equation[SampleRate, 2]")
            ),
            "3/2",
        ),
        (
            format!("Channels[{}]", w("List[0.1, 0.2], Equation[Channels, 2]")),
            "2",
        ),
        (
            cd("List[0.1, 0.2], Equation[Channels, 2]", ""),
            "Sequence[List[7301/16384, 7301/8192], List[7301/16384, 7301/8192]]",
        ),
    ];
    check(&ctx, cases);
}

/// probe9's and probe10's rows: normalisation, shapes, channels, slices.
#[test]
fn normalisation_shapes_slices() {
    let ctx = Context::new();
    let two = "List[List[0.5, 0.25], List[0.1, 0.2]]";
    let s6 = format!("{S6}, Equation[SampleRate, 2]");
    let sines = "List[0, Sin[Times[11/100, Pi]], Sin[Times[11/50, Pi]], \
                 Sin[Times[33/100, Pi]], Sin[Times[11/25, Pi]], Sin[Times[11/20, Pi]], \
                 Sin[Times[33/50, Pi]], Sin[Times[77/100, Pi]]]";
    let half = "List[7301/8192, 7301/16384]";
    let three = format!("Sequence[{half}, {half}, {half}]");
    let cases: &[(String, &str)] = &[
        (
            cd("List[0.5, -0.5], Equation[Normalize, -1]", ""),
            "List[7301/8192, -7301/8192]",
        ),
        (
            cd("List[0.5, -0.5], Equation[Normalize, 1]", ""),
            "List[1/2, -16383/32768]",
        ),
        (
            cd("List[0.5, -0.5], Equation[Normalize, -100]", ""),
            "List[0, 0]",
        ),
        (
            cd("List[0.2, -0.1], Equation[Normalize, 0]", ""),
            "List[6553/32768, -3277/32768]",
        ),
        (
            cd("List[1.5, -0.1], Equation[Normalize, 0]", ""),
            "List[32767/32768, -3277/32768]",
        ),
        (cd("List[0, 0, 0]", ""), "List[0, 0, 0]"),
        (cd("List[0.5, 0.25], Equation[Channels, 3]", ""), &three),
        (
            format!("Channels[{}]", w("List[0.5, 0.25], Equation[Channels, 3]")),
            "3",
        ),
        (
            format!(
                "Channels[{}]",
                w("List[List[0.5, 0.25], List[0.1, 0.2], List[0.3, 0.4]]")
            ),
            "3",
        ),
        (
            cd("Equation[Length, 3], Equation[Channels, 1]", ""),
            "List[0, 0, 0]",
        ),
        (
            cd("List[0.5, 0.25], Equation[Duration, 0.0001]", ""),
            "List[7301/8192, 7301/16384, 0, 0, 0]",
        ),
        (
            cd(sines, ""),
            "List[0, 313/1024, 18847/32768, 12725/16384, 7261/8192, 7301/8192, \
             12955/16384, 19553/32768]",
        ),
        (cd("List[0.5, 0.25]", "0"), half),
        (cd(two, "1"), "List[5841/32768, 11681/32768]"),
        (cd(two, "0"), half),
        (
            cd(S6, "Equation[range, List[4, 10]]"),
            "List[7301/8192, 21903/32768]",
        ),
        (
            cd(&s6, "Equation[range, Range[0.5, 1.0]]"),
            "List[-7301/16384]",
        ),
        // A mono clip with `matrix` answers the list (giac).
        (cd(A, "Matrix"), A0),
    ];
    check(&ctx, cases);
}

/// `SoundSeconds`, probe5:62-68.
#[test]
fn sound_seconds() {
    let ctx = Context::new();
    let cases = [
        ("Length[SoundSeconds[1]]".to_owned(), "44100"),
        ("Length[SoundSeconds[1, 22100]]".to_owned(), "22100"),
        (
            "SoundSeconds[1/1000, 8000]".to_owned(),
            "List[0, 1/8000, 1/4000, 3/8000, 1/2000, 1/1600, 3/4000, 7/8000]",
        ),
        ("SoundSeconds[0.5, 4]".to_owned(), "List[0, 1/4]"),
        (
            "SoundSeconds[2, 3]".to_owned(),
            "List[0, 1/3, 2/3, 1, 4/3, 5/3]",
        ),
    ];
    check(&ctx, &cases);
    for (input, first) in [
        ("SoundSeconds[1]", "List[0, 1/44100, 1/22050, 1/14700]"),
        (
            "SoundSeconds[1, 22100]",
            "List[0, 1/22100, 1/11050, 3/22100]",
        ),
    ] {
        let found = eval(&ctx, input);
        assert_eq!(&found.args()[..4], eval(&ctx, first).args(), "{input}");
    }
}

/// The duration rule exactly: giac's 47-bit doubles make 8 and 16 samples
/// of `duration = 0.001` and `0.002` at 8000 Hz, and of `1/1000` (C0088);
/// ours 9 and 17, departures. The other rows are giac's.
#[test]
fn durations() {
    let ctx = Context::new();
    let size = |args: &str| format!("Length[{}]", cd(args, ""));
    let cases = [
        (
            size("Equation[Duration, 0.002], Equation[Channels, 1], Equation[SampleRate, 8000]"),
            "17",
        ),
        (
            size("Equation[Duration, 1/1000], Equation[Channels, 1], Equation[SampleRate, 8000]"),
            "9",
        ),
        (
            cd(
                "Equation[Duration, 0.001], Equation[Channels, 1], Equation[SampleRate, 8000]",
                "",
            ),
            "List[0, 0, 0, 0, 0, 0, 0, 0, 0]",
        ),
        (
            size("Equation[Duration, 0.5], Equation[Channels, 1], Equation[SampleRate, 4]"),
            "3",
        ),
        (
            size("Equation[Duration, 1], Equation[Channels, 1], Equation[SampleRate, 10]"),
            "11",
        ),
        (
            size("Equation[Duration, 0.25], Equation[Channels, 1], Equation[SampleRate, 10]"),
            "3",
        ),
        (
            size(
                "Equation[Length, 3], Equation[Duration, 1], Equation[Channels, 1], \
                 Equation[SampleRate, 10]",
            ),
            "3",
        ),
        (
            size("List[1, 2], Equation[Duration, 1], Equation[SampleRate, 10]"),
            "11",
        ),
        // giac refuses `size = 3`: ours is the length, a departure.
        (
            w("Equation[Channels, 2], Equation[Length, 3]"),
            "SoundClip[44100, 16, List[0, 0, 0], List[0, 0, 0]]",
        ),
    ];
    check(&ctx, &cases);
}

/// giac's errors stay with its message.
#[test]
fn stays() {
    let ctx = Context::new();
    let ragged = eval(&ctx, "List[0.5, 0.25]");
    let ragged = format!("Failed to set channel data, {ragged} is not a valid sample value");
    let clip = w("List[0.5, -0.5]");
    let cases: &[(String, &str)] = &[
        (
            format!("ChannelData[{}, 2]", w(M)),
            "Invalid channel number",
        ),
        (
            w(&format!("{A}, Equation[BitDepth, 32]")),
            "Bit depth should be either 8, 16 or 24",
        ),
        (
            w("Equation[Duration, 0.0001]"),
            "Number of channels not set",
        ),
        (
            w("Equation[Duration, 0.0001], Equation[SampleRate, 20000]"),
            "Number of channels not set",
        ),
        (w("Equation[Length, 3]"), "Number of channels not set"),
        (w("List[x, 1]"), "Sample data should be numeric"),
        (w("List[List[0.5, 0.25], List[0.1]]"), &ragged),
        (
            w("List[0.5, 0.25], Equation[SampleRate, 0]"),
            "Sample rate should be a positive integer",
        ),
        (
            w("List[0.5, 0.25], Equation[SampleRate, 1.5]"),
            "Sample rate should be a positive integer",
        ),
        (
            w("List[], Equation[SampleRate, 8000]"),
            "Audio length not set",
        ),
        ("SampleRate[List[1, 2]]".to_owned(), "Invalid sound data"),
        ("Duration[List[1, 2]]".to_owned(), "Invalid sound data"),
        (
            "Channels[List[1, 2]]".to_owned(),
            "Expected audio clip or image",
        ),
        (
            "BitDepth[SoundClip[44100, 16, List[40000]]]".to_owned(),
            "Expected audio clip or image",
        ),
        (
            "ChannelData[List[1, 2]]".to_owned(),
            "Expected audio clip or image",
        ),
        (cd("List[0.5, 0.25]", "1"), "Invalid channel number"),
        (
            cd(M, "Right, Equation[range, Range[1, 2]]"),
            "Invalid slice bounds",
        ),
        (
            cd(S6, "Equation[range, Range[4, 10]]"),
            "Invalid slice bounds",
        ),
        (
            cd(S6, "Equation[range, List[1, 2, 3]]"),
            "Invalid range specification",
        ),
        (cd(S6, "\"x\""), "Invalid argument 2"),
        (w(&clip), "Invalid argument 1"),
        (
            w(&format!("{clip}, Equation[BitDepth, 8]")),
            "Invalid argument 1",
        ),
        (
            w("List[1], Equation[Channels, 0]"),
            "Number of channels should be a positive integer",
        ),
        (
            w("List[1], Equation[Length, 0]"),
            "Length should be a positive integer",
        ),
        (
            w("List[1], Equation[Normalize, x]"),
            "Normalization level should be a real number",
        ),
        (
            w("List[1], Equation[Duration, x]"),
            "Duration should be a real number",
        ),
        (
            w("Equation[Channels, 2], List[List[1], List[2], List[3]]"),
            "Number of channels does not match the number of data rows",
        ),
        (w("0"), "Expected a length in samples"),
        ("SoundSeconds[0]".to_owned(), "Invalid sound parameters"),
        (
            "SoundSeconds[1, 1/2]".to_owned(),
            "Invalid sound parameters",
        ),
        ("SoundSeconds[x]".to_owned(), "Invalid sound parameters"),
    ];
    for (input, message) in cases {
        assert_eq!(&warning(&ctx, input), message, "{input}");
    }
}

/// Option words alone, or another arity, stay silently.
#[test]
fn stays_silently() {
    let ctx = Context::new();
    for input in ["Channels[]", "SampleRate[1, 2]", "About[]", "ChannelData[]"] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap(), stayed(&ctx, input), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// giac's `==`: a clip equals itself only.
#[test]
fn equality() {
    let ctx = Context::new();
    let a = w(A);
    let b = w("List[0.5]");
    let cases = [
        (format!("Equal[{a}, {a}]"), "True"),
        (format!("Equal[{a}, {b}]"), "False"),
        (format!("Equal[{a}, 1]"), "False"),
        (format!("Equal[{a}, \"a\"]"), "False"),
    ];
    check(&ctx, &cases);
}

/// `About`: a clip with giac's description, a list mapped, a symbol's
/// assumptions, held.
#[test]
fn about() {
    let ctx = Context::new();
    for (args, message) in [
        (
            A,
            "Bit depth: 16\nChannels: 1 (mono)\nSample rate: 44100\nLength: 4 samples\n\
             Duration: 9.07029e-05 seconds\nLevel: -1.00015 dB",
        ),
        (
            "List[0.5, 0.25], Equation[Channels, 2]",
            "Bit depth: 16\nChannels: 2 (stereo)\nSample rate: 44100\nLength: 2 samples\n\
             Duration: 4.53515e-05 seconds\nLevel: -1.00015 dB",
        ),
        (
            "List[0, 0]",
            "Bit depth: 16\nChannels: 1 (mono)\nSample rate: 44100\nLength: 2 samples\n\
             Duration: 4.53515e-05 seconds\nLevel: -inf dB",
        ),
    ] {
        let clip = eval(&ctx, &w(args));
        let input = format!("About[{}]", w(args));
        let (found, report) = ctx.eval_with(&read(&input), &EvalOptions::default());
        assert_eq!(found.unwrap(), clip, "{input}");
        match &report.warnings[..] {
            [Warning::Native { message: m, .. }] => assert_eq!(m, message, "{input}"),
            other => panic!("{input}: {other:?}"),
        }
    }
    let assumed = ctx
        .assuming(read("Greater[x, 0]"))
        .assuming(read("Less[y, 1]"))
        .assuming(read("Less[x, 5]"));
    let cases = [
        (
            "About[x]".to_owned(),
            "Hold[List[Greater[x, 0], Less[x, 5]]]",
        ),
        ("About[z]".to_owned(), "z"),
        ("About[3]".to_owned(), "3"),
        (
            "About[List[y, 2]]".to_owned(),
            "List[Hold[List[Less[y, 1]]], 2]",
        ),
    ];
    for (input, expected) in &cases {
        let found = eval(&assumed, input);
        let expected = assumed.eval(&read(expected)).unwrap();
        assert_eq!(found, expected, "{input}");
    }
    assert_eq!(eval(&ctx, "About[x]"), read("x"));
}

/// Sizes past `Limits::terms` stay before anything is built.
#[test]
fn limits() {
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    for input in [
        w("Equation[Length, 21], Equation[Channels, 1]"),
        w("Equation[Length, 11], Equation[Channels, 2]"),
        w("Equation[Duration, 1000000], Equation[Channels, 1]"),
        w("Equation[Length, 10000000000000], Equation[Channels, 10000000000000]"),
    ] {
        assert_eq!(warning(&small, &input), "Too many samples", "{input}");
    }
    assert_eq!(warning(&small, "SoundSeconds[1, 21]"), "Too many records");
    check(&small, &[("Length[SoundSeconds[1, 20]]".to_owned(), "20")]);
}

const B: &str = "List[0.1, 0.2, 0.3, 0.4]";
const CLIPPING: &str = "clipping audio, the result will be distorted";

/// `SoundClip[44100, 16, codes…]`.
fn clip16(channels: &[&str]) -> String {
    format!("SoundClip[44100, 16, {}]", channels.join(", "))
}

/// The clip edits (CAS-632): the prototype's rows checked against giac
/// (crates/cas-core/proto/sound/edit/proto.py, logs/T8.c1.md), as codes.
#[test]
#[allow(clippy::too_many_lines, reason = "one row per probe")]
fn edits() {
    let ctx = Context::new();
    let (a, b, m) = (w(A), w(B), w(M));
    let a0 = "List[29204, -29204, 14602, 0]";
    let spliced = clip16(&["List[29204, -29204, 10180, 14602, 21903, 29204]"]);
    let first = clip16(&["List[29490, 26214, 14602, 0]"]);
    let eight = w("List[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]");
    let cases: &[(String, String)] = &[
        (
            format!("SetChannelData[{a}, List[0.9, 0.8]]"),
            first.clone(),
        ),
        (
            format!("SetChannelData[{a}, List[List[0.9, 0.8]]]"),
            first.clone(),
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 1]"),
            clip16(&["List[29204, 29490, 26214, 0]"]),
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8, 0.7, 0.6, 0.5], 2]"),
            clip16(&["List[29204, -29204, 29490, 26214]"]),
        ),
        (
            format!("SetChannelData[{m}, List[0.9], 2, Right]"),
            clip16(&[a0, "List[5841, 11681, 29490, 23363]"]),
        ),
        (
            format!("SetChannelData[{m}, List[List[0.9, 0.8], List[0.1, 0.2]], 1]"),
            clip16(&[
                "List[29204, 29490, 26214, 0]",
                "List[5841, 3277, 6553, 23363]",
            ]),
        ),
        (
            format!("SetChannelData[{m}, List[0.9, 0.8]]"),
            clip16(&[
                "List[29490, 26214, 14602, 0]",
                "List[29490, 26214, 17522, 23363]",
            ]),
        ),
        (
            format!("SetChannelData[{a}, List[0.9, x], 3]"),
            clip16(&["List[29204, -29204, 14602, 29490]"]),
        ),
        (
            format!("SetChannelData[{a}, List[Power[2, -1/2], 1/3, -2, Pi]]"),
            clip16(&["List[23170, 10922, -32767, 32767]"]),
        ),
        (
            format!("SetChannelData[{a}, List[0.9], 0, 0, 7]"),
            clip16(&["List[29490, -29204, 14602, 0]"]),
        ),
        (
            format!(
                "SetChannelData[{}, List[0.9, -1]]",
                w("List[0.5, -0.5], Equation[BitDepth, 8]")
            ),
            "SoundClip[44100, 8, List[114, -127]]".to_owned(),
        ),
        (format!("SetChannelData[{a}, List[]]"), clip16(&[a0])),
        (
            format!("Stereo2Mono[{m}]"),
            clip16(&["List[17522, -8761, 16062, 11681]"]),
        ),
        (
            format!("Stereo2Mono[{}]", w("List[0.5, 0.25]")),
            clip16(&["List[29203, 14602]"]),
        ),
        (
            format!(
                "Stereo2Mono[{}]",
                w("List[List[0.5, 0.25], List[0.1, 0.2], List[0.3, -0.4]]")
            ),
            clip16(&["List[17522, 973]"]),
        ),
        (
            format!("Mixdown[{a}, List[0, -6], {b}, List[0, -6]]"),
            clip16(&["List[18295, -7317, 18295, 14636]"]),
        ),
        (
            format!("Mixdown[{m}, List[0, 0, -0.5], {b}, List[1, -3, 0.5]]"),
            clip16(&[
                "List[26980, -25001, 17445, 5934, 7912]",
                "List[2235, 9245, 16255, 23265, 19100]",
            ]),
        ),
        (
            format!("Mixdown[{a}, List[0, 0, 0, 2, 2], {b}, 1]"),
            clip16(&["List[0, -7301, 21902, 21902, 29203]"]),
        ),
        (
            format!(
                "Mixdown[{}, 0, {b}, 1]",
                w("List[0.5, -0.5], Equation[BitDepth, 8]")
            ),
            clip16(&["List[28927, -21625, 14602, 21902, 29203]"]),
        ),
        (
            format!("Mixdown[{a}, List[0, -6, 0, 1, 1], {b}, List[1/44100, -3]]"),
            clip16(&["List[0, -9467, 17655, 15506, 20674]"]),
        ),
        (
            format!("Length[ChannelData[Mixdown[{a}, 0.5, {b}, 0.5]]]"),
            "22054".to_owned(),
        ),
        (
            format!("Length[ChannelData[Mixdown[{a}, 0.0001, {b}, 0]]]"),
            "8".to_owned(),
        ),
        (
            format!("Length[ChannelData[Mixdown[{a}, Sqrt[2], {b}, 0]]]"),
            "62371".to_owned(),
        ),
        (format!("Splice[{a}, {b}]"), spliced.clone()),
        (format!("Splice[{a}, {b}, 5]"), spliced.clone()),
        (format!("Splice[{a}, {b}, 2, -1]"), spliced.clone()),
        (format!("Splice[{a}, {b}, 1.5]"), spliced.clone()),
        (
            format!("Splice[{a}, {b}, 2, 1]"),
            clip16(&["List[29204, -29204, 10951, 14602, 21903, 29204]"]),
        ),
        (
            format!("Splice[{a}, {b}, 1]"),
            clip16(&["List[29204, -29204, 14602, 7301, 14602, 21903, 29204]"]),
        ),
        (
            format!("Splice[{a}, {b}, 3]"),
            clip16(&["List[29204, -8128, 14602, 21902, 29204]"]),
        ),
        (
            format!("Splice[{a}, {b}, 4]"),
            clip16(&["List[18704, 1125, 20869, 29203]"]),
        ),
        (
            format!("Splice[{}, {b}]", w("List[0, 0, 0, 0]")),
            clip16(&["List[0, 0, 0, 14602, 21903, 29204]"]),
        ),
        (
            format!(
                "Splice[{m}, {}]",
                w(&format!("List[{B}, List[-0.5, -0.5, 0.5, 0.5]]"))
            ),
            clip16(&[
                "List[29204, -29204, 0, 11681, 17522, 23363]",
                "List[5841, 11681, 0, -29203, 29204, 29204]",
            ]),
        ),
        (
            format!("Splice[{a}, {}, 3, 1]", w("List[0.9]")),
            clip16(&["List[29204, -29204]"]),
        ),
        (
            format!("Splice[{eight}, {}, 6, 1]", w("List[0.9, -0.9]")),
            clip16(&["List[3650, 7301, 10951, 14602]"]),
        ),
    ];
    let cases: Vec<(String, &str)> = cases.iter().map(|(i, x)| (i.clone(), x.as_str())).collect();
    check(&ctx, &cases);
}

/// A mix or a splice that clips answers with one warning per call (giac
/// reports the last sample's status only); one that does not, with none.
#[test]
fn edits_clipping() {
    let ctx = Context::new();
    let (a, b, m) = (w(A), w(B), w(M));
    let nine = w("List[0.9, 0.9, 0.9, 0.9]");
    let cases = [
        (
            format!("Mixdown[{a}, 1, {b}, 1]"),
            clip16(&["List[0, 32767, -14601, 32767, 29203]"]),
            true,
        ),
        (
            format!("Mixdown[{a}, 2, {b}, 0]"),
            clip16(&["List[7301, 14602, 32767, 1, 14602, 0]"]),
            true,
        ),
        (
            format!("Mixdown[{a}, List[0, 0, -1], {b}, List[0, 0, 1]]"),
            clip16(&["List[32767, -14601, 32767, 29203]"]),
            true,
        ),
        (
            format!("Mixdown[{a}, 0, {m}, 0]"),
            clip16(&[
                "List[32767, -32767, 20650, 0]",
                "List[24779, -12390, 22714, 16520]",
            ]),
            true,
        ),
        (
            format!("Splice[{nine}, {nine}, 2, 0]"),
            clip16(&["List[29204, 29204, 32767, 32767, 29204, 29204]"]),
            true,
        ),
        (
            format!("Mixdown[{a}, List[0, -6], {b}, List[0, -6]]"),
            clip16(&["List[18295, -7317, 18295, 14636]"]),
            false,
        ),
        (
            format!("Splice[{a}, {b}]"),
            clip16(&["List[29204, -29204, 10180, 14602, 21903, 29204]"]),
            false,
        ),
    ];
    for (input, expected, clips) in &cases {
        let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
        assert_eq!(found.unwrap(), eval(&ctx, expected), "{input}");
        let messages: Vec<String> = report
            .warnings
            .iter()
            .map(|w| match w {
                Warning::Native { message, .. } => message.to_string(),
                other => format!("{other:?}"),
            })
            .collect();
        let wanted: &[&str] = if *clips { &[CLIPPING] } else { &[] };
        assert_eq!(messages, wanted, "{input}");
    }
}

/// giac's errors of the clip edits: the call stays with its message.
#[test]
#[allow(clippy::too_many_lines, reason = "one row per error")]
fn edits_stay() {
    let ctx = Context::new();
    let (a, b, m) = (w(A), w(B), w(M));
    let slow = w(&format!("{B}, Equation[SampleRate, 22050]"));
    let same = "Clips should have the same bit depths, sample rates and channel counts";
    let cases: &[(String, &str)] = &[
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 4]"),
            "Invalid offset",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], -1]"),
            "Invalid offset",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 1.0]"),
            "Expected an integer",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 0, 1]"),
            "Invalid channel number",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 0, Right]"),
            "Invalid channel number",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, 0.8], 0, x]"),
            "Invalid channel specification",
        ),
        (
            format!("SetChannelData[{a}, 0.9]"),
            "Invalid data specification",
        ),
        (
            format!("SetChannelData[{m}, List[List[0.9, 0.8], List[0.1, 0.2]], 0, 0]"),
            "Invalid data specification",
        ),
        (
            format!("SetChannelData[{m}, List[List[0.9, 0.8, 0.7]]]"),
            "Invalid data specification",
        ),
        (
            format!("SetChannelData[{a}, List[0.9, x]]"),
            "Failed to set channel data, x is not a valid sample value",
        ),
        (
            "SetChannelData[List[1, 2], List[0.9]]".to_owned(),
            "Expected audio clip or image",
        ),
        (
            "Stereo2Mono[List[0.5, 0.25]]".to_owned(),
            "Invalid sound data",
        ),
        (format!("Stereo2Mono[{a}, {b}]"), "Invalid sound data"),
        (format!("Mixdown[{a}, {b}]"), "Invalid offset specification"),
        (
            format!("Mixdown[{a}, -6, {b}, 0]"),
            "Invalid offset specification",
        ),
        (
            format!("Mixdown[{a}, List[0.5, -1], {b}, List[0.5, 1]]"),
            "Invalid gain specification",
        ),
        (
            format!("Mixdown[{a}, List[0, x], {b}, 0]"),
            "Invalid gain specification",
        ),
        (
            format!("Mixdown[{a}, List[0, 0, 2], {b}, 1]"),
            "Invalid pan specification",
        ),
        (
            format!("Mixdown[{a}, List[0, 0, 0, -1], {b}, 1]"),
            "Invalid fade-in specification",
        ),
        (
            format!("Mixdown[{a}, List[0, 0, 0, 0, x], {b}, 1]"),
            "Invalid fade-out specification",
        ),
        (
            format!("Mixdown[{a}, List[0, 0, 0, 3, 2], {b}, 1]"),
            "Fade-in/out is too long",
        ),
        (
            format!("Mixdown[{a}, 0, List[1, 2], 0]"),
            "Invalid argument 3",
        ),
        (
            format!("Mixdown[{a}, 0, {slow}, 0]"),
            "All clips must have the same sample rate",
        ),
        (
            format!(
                "Mixdown[{}, 0, {b}, 0]",
                w("List[0.1, 0.2], Equation[Channels, 3]")
            ),
            "Only mono and stereo clips are allowed",
        ),
        (format!("Splice[{a}, {b}, 0]"), "Invalid splice position"),
        (format!("Splice[{a}, {b}, -3]"), "Invalid splice position"),
        (format!("Splice[{a}, {slow}]"), same),
        (format!("Splice[{a}, {slow}, 2, 1]"), same),
        (format!("Splice[{m}, {b}, 2]"), same),
        (
            format!("Splice[{a}, {b}, x]"),
            "Crossfade length should be a real number",
        ),
        (
            format!("Splice[{a}, {b}, 2, x]"),
            "Crossfade exponent should be a real number",
        ),
        (
            format!("Splice[{a}, List[1, 2]]"),
            "First two arguments should be audio clips",
        ),
    ];
    for (input, message) in cases {
        assert_eq!(&warning(&ctx, input), message, "{input}");
    }
}

/// giac's type and dimension errors of the clip edits stay silently, and
/// sizes past `Limits::terms` stay before anything is built.
#[test]
fn edits_silent_and_limits() {
    let ctx = Context::new();
    let (a, b) = (w(A), w(B));
    for input in [
        "Mixdown[]".to_owned(),
        format!("Mixdown[{a}]"),
        format!("Mixdown[{a}, 0, {b}]"),
        format!("Mixdown[{a}, List[0, 0, 0, 0, 0, 0], {b}, 1]"),
        format!("Splice[{a}]"),
        format!("Splice[{a}, {b}, 2, 0, 1]"),
        format!("Splice[{a}, {b}, 3/2]"),
        format!("SetChannelData[{a}]"),
    ] {
        let (found, report) = ctx.eval_with(&read(&input), &EvalOptions::default());
        assert_eq!(found.unwrap(), stayed(&ctx, &input), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    let fifteen = w("Equation[Length, 15], Equation[Channels, 1]");
    for input in [
        format!("Mixdown[{a}, 1000000.0, {b}, 0]"),
        format!("Mixdown[{a}, 1000000000000000000000000000000, {b}, 0]"),
        format!("Splice[{fifteen}, {fifteen}]"),
    ] {
        assert_eq!(warning(&small, &input), "Too many samples", "{input}");
    }
}

/// `SoundClip[rate, depth, codes…]`.
fn sc(rate: u32, depth: u32, channels: &[&str]) -> String {
    format!("SoundClip[{rate}, {depth}, {}]", channels.join(", "))
}

/// Resampling (CAS-633): the qualities 0 and 1 exactly, the others on
/// `A`, from the prototype (crates/cas-core/proto/sound/edit/resample.py, resample.log).
#[test]
#[allow(clippy::too_many_lines, reason = "one row per case")]
fn resample() {
    let ctx = Context::new();
    let a = sc(44100, 16, &["List[29204, -29204, 14602, 0]"]);
    let b = "List[5841, 11681, 17522, 23363]";
    let ab = sc(44100, 16, &["List[29204, -29204, 14602, 0]", b]);
    let one = sc(44100, 16, &["List[29204]"]);
    let eight = sc(8000, 8, &["List[64, -64, 32, 0]"]);
    let deep = sc(44100, 24, &["List[4194304, -4194304, 2097152, 0]"]);
    let cases: &[(String, String)] = &[
        (
            format!("Resample[{a}, 88200, 0]"),
            sc(
                88200,
                16,
                &["List[29203, 29203, 29203, 0, -29203, -7301, 14602, 7301]"],
            ),
        ),
        (
            format!("Resample[{a}, 88200, 1]"),
            sc(
                88200,
                16,
                &["List[29203, 29203, 29203, 29203, -29203, -29203, 14602, 14602]"],
            ),
        ),
        (
            format!("Resample[{a}, 22050, 0]"),
            sc(22050, 16, &["List[29203, -29203]"]),
        ),
        (
            format!("Resample[{a}, 22050, 1]"),
            sc(22050, 16, &["List[29203, -29203]"]),
        ),
        (
            format!("Resample[{a}, 30000, 0]"),
            sc(30000, 16, &["List[29203, 1752]"]),
        ),
        (
            format!("Resample[{a}, 30000, 1]"),
            sc(30000, 16, &["List[29203, 29203]"]),
        ),
        (
            format!("Resample[{a}, 132300, 0]"),
            sc(
                132_300,
                16,
                &[
                    "List[29203, 29203, 29203, 29203, 9734, -9734, -29203, -14602, 0, 14602, 9734, 4867]",
                ],
            ),
        ),
        (
            format!("Resample[{one}, 132300, 0]"),
            sc(132_300, 16, &["List[29203, 29203, 29203]"]),
        ),
        (
            format!("Resample[{one}, 132300, 1]"),
            sc(132_300, 16, &["List[29203, 29203, 29203]"]),
        ),
        (
            format!("Resample[{eight}, 20000, 0]"),
            sc(
                20000,
                8,
                &["List[64, 64, 64, 38, -13, -63, -25, 13, 25, 13]"],
            ),
        ),
        (
            format!("Resample[{eight}, 20000, 1]"),
            sc(
                20000,
                8,
                &["List[64, 64, 64, 64, 64, -63, -63, -63, 32, 32]"],
            ),
        ),
        (
            format!("Resample[{deep}, 88200, 0]"),
            sc(
                88200,
                24,
                &["List[4194304, 4194304, 4194304, 0, -4194303, -1048576, 2097152, 1048576]"],
            ),
        ),
        (
            format!("Resample[{ab}, 66150, 0]"),
            sc(
                66150,
                16,
                &[
                    "List[29203, 29203, 9734, -29203, 0, 9734]",
                    "List[5841, 5841, 7787, 11681, 15575, 19468]",
                ],
            ),
        ),
        (
            format!("Resample[{ab}, 66150, 1]"),
            sc(
                66150,
                16,
                &[
                    "List[29203, 29203, 29203, -29203, -29203, 14602]",
                    "List[5841, 5841, 5841, 11681, 11681, 17521]",
                ],
            ),
        ),
        (
            format!("Resample[{a}]"),
            sc(44100, 16, &["List[22068, -21986, 7576, 6577]"]),
        ),
        (
            format!("Resample[{a}, 44100, 2]"),
            sc(44100, 16, &["List[22068, -21986, 7576, 6577]"]),
        ),
        (
            format!("Resample[{a}, 22050]"),
            sc(22050, 16, &["List[4683, -1172]"]),
        ),
        (
            format!("Resample[{a}, 88200, 3]"),
            sc(
                88200,
                16,
                &["List[25574, -3008, -25563, -15281, 10986, 18701, 3556, -8673]"],
            ),
        ),
        (
            format!("Resample[{a}, 48000, 4]"),
            sc(48000, 16, &["List[27744, -26007, 4832, 11250]"]),
        ),
    ];
    let cases: Vec<(String, &str)> = cases.iter().map(|(i, o)| (i.clone(), o.as_str())).collect();
    check(&ctx, &cases);
}

/// Resampling's errors: giac's messages, ours for an empty result, odd
/// shapes silently, and the limits.
#[test]
fn resample_stays() {
    let ctx = Context::new();
    let a = sc(44100, 16, &["List[29204, -29204, 14602, 0]"]);
    let (rate, quality) = (
        "Invalid sample rate specification",
        "Invalid quality specification",
    );
    let cases: &[(String, &str)] = &[
        (format!("Resample[{a}, 0]"), rate),
        (format!("Resample[{a}, -3]"), rate),
        (format!("Resample[{a}, 22050.0]"), rate),
        (format!("Resample[{a}, 1/2]"), rate),
        (format!("Resample[{a}, 22050, 5]"), quality),
        (format!("Resample[{a}, 22050, -1]"), quality),
        (format!("Resample[{a}, 22050, 2.0]"), quality),
        ("Resample[x]".to_owned(), "Invalid sound data"),
        ("Resample[List[1, 2, 3]]".to_owned(), "Invalid sound data"),
        ("Resample[x, 22050]".to_owned(), "Invalid sound data"),
        ("Resample[x, 0]".to_owned(), rate),
        ("Resample[x, 22050, 7]".to_owned(), quality),
        (
            format!("Resample[{a}, 1]"),
            "The resampled clip would have no samples",
        ),
    ];
    for (input, message) in cases {
        assert_eq!(warning(&ctx, input), *message, "{input}");
    }
    for input in [
        "Resample[]".to_owned(),
        format!("Resample[{a}, 22050, 2, 1]"),
    ] {
        let (found, report) = ctx.eval_with(&read(&input), &EvalOptions::default());
        assert_eq!(found.unwrap(), stayed(&ctx, &input), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    let silent = |n: usize| format!("List[{}]", vec!["0"; n].join(", "));
    let two = sc(2, 16, &[&silent(4)]);
    let stereo = sc(2, 16, &[&silent(4), &silent(4)]);
    let long = sc(3000, 16, &[&silent(3000)]);
    for input in [
        format!("Resample[{a}, 21, 0]"),
        format!("Resample[{two}, 20, 0]"),
        format!("Resample[{stereo}, 6, 0]"),
        format!("Resample[{long}, 10, 4]"),
    ] {
        assert_eq!(warning(&small, &input), "Too many samples", "{input}");
    }
    check(
        &small,
        &[
            (
                format!("Resample[{long}, 10, 0]"),
                &sc(10, 16, &[&silent(10)]),
            ),
            (
                format!("Resample[{stereo}, 5, 0]"),
                &sc(5, 16, &[&silent(10), &silent(10)]),
            ),
        ],
    );
}
