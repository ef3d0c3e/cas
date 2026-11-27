//! Sound clips (CAS-631), Xcas's audio clips without the files: giac's
//! `createwav`, `soundsec`, `samplerate`, `bit_depth`, `channels`,
//! `duration`, `channel_data` and `about` (giac signalprocessing.cc
//! `_soundsec` :2069-2088, `_createwav` :2094-2190, `_channels` :2441,
//! `_bit_depth` :2455, `_samplerate` :2469, `_duration` :2480,
//! `_channel_data` :2491-2549, the clip's encoding :1292-1308; maple.cc
//! `_about` :323-366).
//!
//! The judgment: a clip is the inert call `SoundClip[rate, depth, {codes of
//! channel 0}, {codes of channel 1}, …]` (no native), giac's integer codes
//! `k`, `|k| ≤ 2^(b−1) − 1` at the bit depth `b` (8, 16 or 24), the sample
//! `k/2^(b−1)` exactly. It prints as a call everywhere (the session,
//! `string()`, the language hook), so what is printed reads back. `Plus`,
//! `Times`, `Part`, `Length` and `N` see a call like any other, departures
//! from giac's clip arithmetic: `A + B` stays, `A − A` is 0, `−A` stays,
//! `size(A)` is `2 +` the channels, `A[0]` stays (`At` selects in lists),
//! `N(A)` turns the codes to decimals (no longer a clip). `==` on clips is giac's: a well-formed
//! clip compares like a number (`logic.rs` `equal`), so two different clips,
//! or a clip and a number or a string, are not equal. The natives are pure:
//! nothing to declare volatile.
//!
//! `CreateWav` takes data as a list (one channel), a matrix (its rows the
//! channels) or `ChannelData = list`, and options `name = value` whose left
//! side is the builtin's symbol: `Channels`, `BitDepth` (16), `SampleRate`
//! (44100), `Length` (`size` or `length`; giac refuses `size = n`, a
//! departure), `Normalize` (by default the peak at −1 dBFS, one gain over
//! all channels; `normalize = 0` or more keeps the samples) and `Duration`
//! (`d` gives `1 + floor(d·sr)` samples, exactly: an exact `d` exactly, a
//! decimal by its midpoint's exact value, where giac's doubles, which keep
//! 47 fraction bits (C0088), make 8 samples of `duration = 0.001` at 8000
//! Hz and ours 9, a departure). The samples are encoded as giac does, in
//! `f64`: values to doubles (rationals correctly rounded, decimals by their
//! midpoint, other numeric constants approximated first), times the gain,
//! clamped to `[−1, 1]`, `floor(range·v + 1/2)` with `range = 2^(b−1) − 1`;
//! a silent clip normalised (gain ∞, `0·∞` NaN) gives code 0.
//!
//! `SoundSeconds[n, N]` (`soundsec`, `N` 44100 by default) is `[k/N : k <
//! n·N]`, exact, the count exact (decimals by their midpoint). `ChannelData`
//! answers exact samples: a channel from 0, `left`, `right`, `matrix`,
//! `range = [start, n]` in samples or `range = a .. b` in seconds (start
//! `floor(1/2 + sr·a)`, length `floor(1/2 + sr·(b − a))`); all channels are
//! a `Sequence` of lists, or a matrix with `matrix`. `About` answers a clip
//! with giac's description as a note, a list mapped, a symbol the
//! assumptions that mention it, held (`Hold[{facts}]`, printed
//! `hold([x > 0])`, so re-evaluation cannot decide them; the symbol itself
//! when none) and anything else itself.
//!
//! giac's errors stay with its message; sizes are checked against
//! `Limits::terms` before anything is built (C0019).
//!
//! The clip edits (CAS-632), `SetChannelData` (`set_channel_data`,
//! :2554-2602), `Stereo2Mono` (`stereo2mono`, :2339), `Mixdown`
//! (`mixdown`, :1791-1843) and `Splice` (`splice`, :1738-1773), live in
//! [`edits`]; `Resample` (`resample`, CAS-633, :2859-2881) in [`resample`].
//! A real argument not structurally numeric is refused before
//! `approximation`, which would evaluate it again (CAS-650).

use std::cmp::Ordering;

use cas_num::ball::{Float, Round};
use cas_num::{Integer, Rational};

use super::approximation;
use super::statistics::random_variables::{equation, range};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::language::Form;
use crate::symbol::Symbol;
use crate::trace::Warning;

mod edits;
mod resample;

pub(super) use edits::{mixdown, set_channel_data, splice, stereo2mono};
pub(super) use resample::resample;

type Native = Result<Option<Expr>, EvalError>;

const BAD_SOUND: &str = "Invalid sound data";
const NOT_CLIP: &str = "Expected audio clip or image";
const BAD_RANGE: &str = "Invalid range specification";

/// A well-formed clip.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(in crate::native) struct Clip {
    /// Samples per second, at least 1.
    pub(in crate::native) rate: i64,
    /// 8, 16 or 24.
    pub(in crate::native) depth: u32,
    /// The codes of each channel, at least one channel, all of one length
    /// of at least 1.
    pub(in crate::native) channels: Vec<Vec<i64>>,
}

impl Clip {
    fn len(&self) -> usize {
        self.channels[0].len()
    }

    /// The sample of code `k`, exactly.
    fn sample(&self, k: i64) -> Expr {
        Expr::from(Rational::new(
            Integer::from(k),
            Integer::from(1_i64 << (self.depth - 1)),
        ))
    }
}

/// `SoundClip[rate, depth, {codes}, …]` when well formed.
pub(in crate::native) fn clip(e: &Expr) -> Option<Clip> {
    if e.head().and_then(Expr::as_symbol) != Some(Symbol::SOUND_CLIP) {
        return None;
    }
    let [rate, depth, channels @ ..] = e.args() else {
        return None;
    };
    let rate = rate.as_i64().filter(|&r| r >= 1)?;
    let depth = match depth.as_i64()? {
        8 => 8,
        16 => 16,
        24 => 24,
        _ => return None,
    };
    let top = (1_i64 << (depth - 1)) - 1;
    let first = channels.first()?;
    let len = first.args().len();
    if len == 0 {
        return None;
    }
    let mut codes = Vec::with_capacity(channels.len());
    for c in channels {
        if !is_list(c) || c.args().len() != len {
            return None;
        }
        let row: Option<Vec<i64>> = c
            .args()
            .iter()
            .map(|k| k.as_i64().filter(|k| k.abs() <= top))
            .collect();
        codes.push(row?);
    }
    Some(Clip {
        rate,
        depth,
        channels: codes,
    })
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: impl Into<String>) -> Native {
    ev.warn(Warning::Native {
        message: message.into().into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// Whether `e` is the symbol `s` or a symbol named `name` (C0008).
fn is(e: &Expr, s: Symbol, name: &str) -> bool {
    e.as_symbol().is_some_and(|t| t == s || t.name() == name)
}

fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// A real number: its double as giac takes it, and its exact value when it
/// has one in reach (a decimal's midpoint).
struct Real {
    f: f64,
    exact: Option<Rational>,
}

/// `x` as a real number: exact numbers, decimals by their midpoint, other
/// numeric constants approximated first (C0028: no `to_rational` panics).
fn real(ev: &Eval<'_>, x: &Expr) -> Option<Real> {
    let from_ball = |e: &Expr| match e.view() {
        View::Ball(b) if b.is_finite() => Some(Real {
            f: b.mid().to_f64(),
            exact: b.mid().checked_to_rational(),
        }),
        _ => None,
    };
    match x.view() {
        View::Int(_) | View::Integer(_) | View::Rational(_) => {
            let q = x.to_rational()?;
            Some(Real {
                f: to_f64(&q),
                exact: Some(q),
            })
        }
        View::Ball(_) => from_ball(x),
        View::Complex(_) | View::ComplexBall(_) | View::String(_) => None,
        View::Symbol(_) | View::Apply { .. } => {
            if !x.is_numeric() {
                return None;
            }
            let a = approximation(ev.context(), x, 128).ok()?;
            from_ball(&a)
        }
    }
}

/// `e` as the language prints it, else in full form.
fn shown(ev: &Eval<'_>, e: &Expr) -> String {
    let ctx = ev.context();
    ctx.language()
        .and_then(|l| l.print(e, Form::Input, ctx, 1000).ok())
        .unwrap_or_else(|| e.to_string())
}

/// giac's encoding of the sample `s` at `depth` bits, in `f64`.
#[allow(
    clippy::cast_possible_truncation,
    reason = "|range·v + 1/2| < 2^24: the floor fits"
)]
fn encode(depth: u32, s: f64) -> i64 {
    if s.is_nan() {
        return 0;
    }
    let range = f64::from((1_u32 << (depth - 1)) - 1);
    let v = s.clamp(-1.0, 1.0);
    (range * v + 0.5).floor() as i64
}

fn positive_int(x: &Expr) -> Option<i64> {
    x.as_i64().filter(|&n| n >= 1)
}

/// `CreateWav[data, options…]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's options and checks, in its order"
)]
pub(super) fn create_wav(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    let mut data: Option<&Expr> = None;
    let mut length: i64 = 0;
    let mut rest = args;
    if let [n] = args
        && let Some(n) = n.to_integer()
    {
        if n.sign() != Ordering::Greater {
            return warn(ev, e, "Expected a length in samples");
        }
        let Some(n) = n.to_i64() else {
            return warn(ev, e, "Too many samples");
        };
        length = n;
        rest = &[];
    }
    let (mut depth, mut rate, mut channels) = (16_u32, 44100_i64, 0_i64);
    let mut norm = 1.0_f64;
    let mut duration: Option<Real> = None;
    for (i, x) in rest.iter().enumerate() {
        let i = i + 1;
        if let Some((lhs, v)) = equation(x) {
            let Some(name) = lhs.as_symbol() else {
                return warn(ev, e, format!("Invalid argument {i}"));
            };
            match name {
                Symbol::CHANNELS => match positive_int(v) {
                    Some(n) => channels = n,
                    None => return warn(ev, e, "Number of channels should be a positive integer"),
                },
                Symbol::BIT_DEPTH => match v.as_i64() {
                    Some(8) => depth = 8,
                    Some(16) => depth = 16,
                    Some(24) => depth = 24,
                    _ => return warn(ev, e, "Bit depth should be either 8, 16 or 24"),
                },
                Symbol::SAMPLE_RATE => match positive_int(v) {
                    Some(n) => rate = n,
                    None => return warn(ev, e, "Sample rate should be a positive integer"),
                },
                Symbol::LENGTH => match positive_int(v) {
                    Some(n) => length = n,
                    None => return warn(ev, e, "Length should be a positive integer"),
                },
                Symbol::NORMALIZE => match real(ev, v) {
                    Some(r) => norm = -r.f,
                    None => return warn(ev, e, "Normalization level should be a real number"),
                },
                Symbol::DURATION => match real(ev, v) {
                    Some(r) if r.exact.is_some() => duration = Some(r),
                    _ => return warn(ev, e, "Duration should be a real number"),
                },
                Symbol::CHANNEL_DATA => {
                    if !is_list(v) {
                        return warn(
                            ev,
                            e,
                            "Channel data should be either a matrix, vector, or string",
                        );
                    }
                    data = Some(v);
                }
                _ => return warn(ev, e, format!("Invalid argument {i}")),
            }
        } else if is_list(x) && data.is_none() {
            data = Some(x);
        } else {
            return warn(ev, e, format!("Invalid argument {i}"));
        }
    }
    let terms = ev.context().limits().terms;
    if let Some(d) = duration.as_ref().and_then(|d| d.exact.as_ref())
        && d.sign() != Ordering::Less
        && length == 0
    {
        let n = (&Rational::from(Integer::from(rate)) * d).floor();
        match n
            .to_i64()
            .filter(|&n| usize::try_from(n).is_ok_and(|n| n < terms))
        {
            Some(n) => length = n + 1,
            None => return warn(ev, e, "Too many samples"),
        }
    }
    let rows = data.map(Expr::args).unwrap_or_default();
    let matrix = data.is_some_and(|_| {
        !rows.is_empty()
            && rows.iter().all(is_list)
            && rows.iter().all(|r| r.args().len() == rows[0].args().len())
            && !rows[0].args().is_empty()
    });
    if data.is_some() {
        if length == 0 {
            let n = if matrix {
                rows[0].args().len()
            } else {
                rows.len()
            };
            length = i64::try_from(n).unwrap_or(i64::MAX);
        }
        if channels == 0 {
            channels = if matrix {
                i64::try_from(rows.len()).unwrap_or(i64::MAX)
            } else {
                1
            };
        } else if matrix && usize::try_from(channels).ok() != Some(rows.len()) {
            return warn(
                ev,
                e,
                "Number of channels does not match the number of data rows",
            );
        }
    }
    if channels == 0 {
        return warn(ev, e, "Number of channels not set");
    }
    if length == 0 {
        return warn(ev, e, "Audio length not set");
    }
    let (Ok(len), Ok(nc)) = (usize::try_from(length), usize::try_from(channels)) else {
        return warn(ev, e, "Too many samples");
    };
    if len.checked_mul(nc).is_none_or(|n| n > terms) {
        return warn(ev, e, "Too many samples");
    }
    let mut codes = vec![vec![0_i64; len]; nc];
    if data.is_some() {
        let mut ratio = 1.0_f64;
        if norm > 0.0 {
            let mut peak = 0.0_f64;
            let mut stack: Vec<&Expr> = rows.iter().collect();
            while let Some(y) = stack.pop() {
                if is_list(y) {
                    stack.extend(y.args());
                } else if let Some(r) = real(ev, y) {
                    peak = peak.max(r.f.abs());
                } else {
                    return warn(ev, e, "Sample data should be numeric");
                }
            }
            ratio = 10.0_f64.powf(-norm / 20.0) / peak;
        }
        for (c, out) in codes.iter_mut().enumerate() {
            let row = if matrix { rows[c].args() } else { rows };
            for (p, y) in row.iter().take(len).enumerate() {
                let Some(r) = (!is_list(y)).then(|| real(ev, y)).flatten() else {
                    let text = shown(ev, y);
                    return warn(
                        ev,
                        e,
                        format!("Failed to set channel data, {text} is not a valid sample value"),
                    );
                };
                out[p] = encode(depth, ratio * r.f);
            }
        }
    }
    let mut parts = Vec::with_capacity(nc + 2);
    parts.push(Expr::from(Integer::from(rate)));
    parts.push(Expr::from(Integer::from(depth)));
    parts.extend(codes.into_iter().map(|row| {
        list(
            row.into_iter()
                .map(|k| Expr::from(Integer::from(k)))
                .collect(),
        )
    }));
    Ok(Some(Expr::call(Symbol::SOUND_CLIP, parts)))
}

/// `SoundSeconds[n]`, `SoundSeconds[n, N]`.
pub(super) fn sound_seconds(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (n, rate) = match e.args() {
        [n] => (n, Some(44100)),
        [n, r] => (n, positive_int(r)),
        _ => return warn(ev, e, "Invalid sound parameters"),
    };
    let n = real(ev, n)
        .and_then(|r| r.exact)
        .filter(|q| q.sign() == Ordering::Greater);
    let (Some(n), Some(rate)) = (n, rate) else {
        return warn(ev, e, "Invalid sound parameters");
    };
    let count = (&n * &Rational::from(Integer::from(rate))).ceil();
    let terms = ev.context().limits().terms;
    let Some(count) = count
        .to_i64()
        .filter(|&c| usize::try_from(c).is_ok_and(|c| c <= terms))
    else {
        return warn(ev, e, "Too many records");
    };
    let rate = Integer::from(rate);
    Ok(Some(list(
        (0..count)
            .map(|k| Expr::from(Rational::new(Integer::from(k), rate.clone())))
            .collect(),
    )))
}

/// `SampleRate`, `BitDepth`, `Channels` and `Duration` of a clip.
pub(super) fn property(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(name) = e.head().and_then(Expr::as_symbol) else {
        return Ok(None);
    };
    let [x] = e.args() else {
        return Ok(None);
    };
    let Some(c) = clip(x) else {
        let message = match name {
            Symbol::SAMPLE_RATE | Symbol::DURATION => BAD_SOUND,
            _ => NOT_CLIP,
        };
        return warn(ev, e, message);
    };
    let int = |n: usize| Expr::from(Integer::from(n));
    Ok(Some(match name {
        Symbol::SAMPLE_RATE => Expr::from(Integer::from(c.rate)),
        Symbol::BIT_DEPTH => Expr::from(Integer::from(c.depth)),
        Symbol::CHANNELS => int(c.channels.len()),
        Symbol::DURATION => {
            Expr::from(Rational::new(Integer::from(c.len()), Integer::from(c.rate)))
        }
        _ => return Ok(None),
    }))
}

/// An integer, saturated to `i64`.
fn saturated(n: &Integer) -> i64 {
    n.to_i64()
        .unwrap_or(if n.is_negative() { i64::MIN } else { i64::MAX })
}

/// `ChannelData[clip, options…]`.
pub(super) fn channel_data(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, options @ ..] = e.args() else {
        return Ok(None);
    };
    let Some(c) = clip(x) else {
        return warn(ev, e, NOT_CLIP);
    };
    let len = i64::try_from(c.len()).unwrap_or(i64::MAX);
    let rate = Rational::from(Integer::from(c.rate));
    let half = Rational::new(Integer::from(1), Integer::from(2));
    let (mut channel, mut start, mut count, mut matrix) = (-1_i64, 0_i64, len, false);
    for (i, o) in options.iter().enumerate() {
        if let Some(n) = o.to_integer() {
            channel = saturated(&n);
        } else if is(o, Symbol::LEFT, "left") {
            channel = 0;
        } else if is(o, Symbol::RIGHT, "right") {
            channel = 1;
        } else if is(o, Symbol::MATRIX, "matrix") {
            matrix = true;
        } else if let Some((lhs, v)) = equation(o) {
            // `range` is giac's `range` builtin too (CAS-644).
            if !is(lhs, Symbol::RANGE, "range") && !is(lhs, Symbol::ARITHMETIC_RANGE, "range") {
                continue;
            }
            if is_list(v) {
                let [first, size] = v.args() else {
                    return warn(ev, e, BAD_RANGE);
                };
                let (Some(first), Some(size)) = (first.to_integer(), size.to_integer()) else {
                    return warn(ev, e, BAD_RANGE);
                };
                (start, count) = (saturated(&first), saturated(&size));
            } else if let Some((from, to)) = range(v) {
                let from = real(ev, from).and_then(|r| r.exact);
                let to = real(ev, to).and_then(|r| r.exact);
                let (Some(from), Some(to)) = (from, to) else {
                    return warn(ev, e, BAD_RANGE);
                };
                start = saturated(&(&half + &(&rate * &from)).floor());
                count = saturated(&(&half + &(&rate * &(&to - &from))).floor());
            } else {
                return warn(ev, e, BAD_RANGE);
            }
        } else {
            return warn(ev, e, format!("Invalid argument {}", i + 2));
        }
    }
    let count = count.min(len.saturating_sub(start));
    if start < 0 || start >= len || count <= 0 {
        return warn(ev, e, "Invalid slice bounds");
    }
    let (Ok(start), Ok(count)) = (usize::try_from(start), usize::try_from(count)) else {
        return warn(ev, e, "Invalid slice bounds");
    };
    let get = |codes: &[i64]| {
        list(
            codes[start..start + count]
                .iter()
                .map(|&k| c.sample(k))
                .collect(),
        )
    };
    if channel < 0 {
        let mut rows: Vec<Expr> = c.channels.iter().map(|ch| get(ch)).collect();
        if rows.len() == 1 {
            return Ok(rows.pop());
        }
        let head = if matrix {
            Symbol::LIST
        } else {
            Symbol::SEQUENCE
        };
        return Ok(Some(Expr::call(head, rows)));
    }
    match usize::try_from(channel)
        .ok()
        .and_then(|k| c.channels.get(k))
    {
        Some(ch) => Ok(Some(get(ch))),
        None => warn(ev, e, "Invalid channel number"),
    }
}

/// C++ ostream's default output of a double: `%g`, 6 significant digits.
fn g6(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_owned();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    let strip = |s: String| {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_owned()
        } else {
            s
        }
    };
    let sci = format!("{x:.5e}");
    let (mantissa, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    if !(-4..6).contains(&exp) {
        let sign = if exp < 0 { '-' } else { '+' };
        return format!("{}e{sign}{:02}", strip(mantissa.to_owned()), exp.abs());
    }
    let decimals = usize::try_from(5 - exp).unwrap_or(0);
    strip(format!("{x:.decimals$}"))
}

/// giac's `about` lines for a clip.
#[allow(
    clippy::cast_precision_loss,
    reason = "codes below 2^24 and lengths below Limits::terms: exact doubles"
)]
fn description(c: &Clip) -> String {
    let nc = c.channels.len();
    let kind = match nc {
        1 => " (mono)",
        2 => " (stereo)",
        _ => "",
    };
    let len = c.len();
    let scale = f64::from(1_u32 << (c.depth - 1));
    let peak = c
        .channels
        .iter()
        .flatten()
        .map(|&k| k as f64 / scale)
        .fold(0.0_f64, f64::max);
    let level = 20.0 * peak.log10();
    format!(
        "Bit depth: {}\nChannels: {nc}{kind}\nSample rate: {}\nLength: {len} samples\n\
         Duration: {} seconds\nLevel: {} dB",
        c.depth,
        c.rate,
        g6(len as f64 / c.rate as f64),
        g6(level)
    )
}

fn about_one(ev: &mut Eval<'_>, x: &Expr) -> Expr {
    if let Some(c) = clip(x) {
        let message = description(&c);
        ev.warn(Warning::Native {
            message: message.into(),
            expr: x.clone(),
        });
        return x.clone();
    }
    if is_list(x) {
        return list(x.args().iter().map(|y| about_one(ev, y)).collect());
    }
    if x.as_symbol().is_some() {
        let facts: Vec<Expr> = ev
            .context()
            .assumptions()
            .iter()
            .filter(|a| !a.free_of(x))
            .cloned()
            .collect();
        if !facts.is_empty() {
            // Held: re-evaluated under the same assumptions, each fact
            // would turn into True (C0026).
            return Expr::call(Symbol::HOLD, [list(facts)]);
        }
    }
    x.clone()
}

/// `About[x]`.
pub(super) fn about(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(about_one(ev, x)))
}

#[cfg(test)]
mod tests;
