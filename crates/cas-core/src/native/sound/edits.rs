//! Clip edits (CAS-632): giac's `set_channel_data`, `stereo2mono`,
//! `mixdown` and `splice` (giac signalprocessing.cc `_set_channel_data`
//! :2554-2602 with `audio_clip::set_channel_data`/`set_block` :1515-1544;
//! `_stereo2mono` :2339, `to_mono` :1223-1232; `_mixdown` :1791-1843,
//! `is_audio_position` :1777, `mix` :1164-1190; `_splice` :1738-1773,
//! `splice` :1194-1221, `similarity` :1130-1149; `get_sample`/`set_sample`
//! :1257, :1090). Each answers a new clip (giac's `set_channel_data` changes
//! its argument); the arithmetic is giac's, in `f64` and in its order:
//! samples decoded `k/2^(b−1)`, encoded as `CreateWav` does, so a decoded
//! top code encodes one lower (`0.5` normalised, `29204`, mixed with
//! silence gives `29203`).
//!
//! `SetChannelData[clip, data, offset, channel]`: `data` a list (written to
//! every channel when none is named) or a matrix, one row per channel; the
//! samples from `offset` (an integer, `0 ≤ offset < len`), those past the
//! end dropped and not read; `channel` an integer, `left` or `right`.
//! `Stereo2Mono[clip]`: the mean of the channels. `Mixdown[clip₁, parm₁,
//! …]`: mono and stereo clips of one rate into a clip of the largest depth
//! and channel count, of length `max(offset + len)`, each clip mixed in turn
//! into the encoded result; `parm` an offset or `[offset, gain, pan,
//! fade_in, fade_out]`, offsets and fades an integer in samples or a
//! decimal in seconds (`floor(1/2 + sr·x)`), the gain in dB at most 0, the
//! pan in `[−1, 1]` at the angle `θ = π(pan + 1)/4` (`cos θ` left, `sin θ`
//! right; none for a mono result), fades `cos²`. `Splice[a, b, n, p]`: `b`
//! after `a`, the last `n` samples of `a` (an integer, or a decimal in
//! seconds) crossfaded with `b`'s first by `((i + 1)/n)^p` (`p = 0` a plain
//! sum); `n` past `a`'s length, or absent, is giac's default
//! `floor(1/2 + sr·min(min(dur_a, dur_b)/2, 1/2))`, `p` negative or absent
//! `1/2 + |similarity|/2` (the cosine similarity of the overlap, the
//! geometric mean over the channels; a silent overlap's is NaN, which
//! encodes as code 0 except the last sample, `1^NaN = 1`). A `b` shorter
//! than the overlap leaves a result of `len_a − n + len_b` samples, giac's
//! copies past its end dropped.
//!
//! Departures: `Splice` refuses clips of other depths, rates or channel
//! counts always, giac only when it computes the similarity (with `p ≥ 0`
//! it mixes them, reading past the second clip's data); a result that
//! clips warns once per call, where giac reports only the last sample's
//! status. giac's errors stay with its message; its type and dimension
//! errors (arity, an exact non-integer overlap, a parm of more than five
//! entries) stay silently. Sizes are checked against `Limits::terms`
//! before anything is built (C0019).

use std::f64::consts::{FRAC_PI_2, PI};

use cas_num::Integer;

use super::{
    BAD_SOUND, Clip, NOT_CLIP, Native, clip, encode, is, is_list, list, real, saturated, shown,
    warn,
};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

const DATA: &str = "Invalid data specification";
const CHANNEL: &str = "Invalid channel number";
const SAME: &str = "Clips should have the same bit depths, sample rates and channel counts";
const CLIPPING: &str = "clipping audio, the result will be distorted";
const TOO_MANY: &str = "Too many samples";

/// giac's decoding of the code `k`.
#[allow(
    clippy::cast_precision_loss,
    reason = "codes below 2^24: exact doubles"
)]
fn decode(depth: u32, k: i64) -> f64 {
    k as f64 / f64::from(1_u32 << (depth - 1))
}

/// The count `x` as a double.
#[allow(
    clippy::cast_precision_loss,
    reason = "counts below Limits::terms, rates as giac's int to double"
)]
fn double(x: impl Into<i128>) -> f64 {
    x.into() as f64
}

/// `floor(1/2 + rate·x)` saturated to `i64` (giac's `int` overflows).
#[allow(
    clippy::cast_possible_truncation,
    reason = "a saturating float-to-int conversion, as intended"
)]
fn samples(rate: i64, x: f64) -> i64 {
    (0.5 + double(rate) * x).floor() as i64
}

fn count(n: usize) -> i128 {
    i128::try_from(n).unwrap_or(i128::MAX)
}

/// The clip as the call `SoundClip[rate, depth, {codes}…]`.
fn answer(c: Clip) -> Expr {
    let mut parts = Vec::with_capacity(c.channels.len() + 2);
    parts.push(Expr::from(Integer::from(c.rate)));
    parts.push(Expr::from(Integer::from(c.depth)));
    parts.extend(c.channels.into_iter().map(|row| {
        list(
            row.into_iter()
                .map(|k| Expr::from(Integer::from(k)))
                .collect(),
        )
    }));
    Expr::call(Symbol::SOUND_CLIP, parts)
}

fn clipping(ev: &mut Eval<'_>, e: &Expr) {
    ev.warn(Warning::Native {
        message: CLIPPING.to_owned().into(),
        expr: e.clone(),
    });
}

/// giac's `ckmatrix`: rows, all lists of one nonzero length.
fn is_matrix(rows: &[Expr]) -> bool {
    !rows.is_empty()
        && rows.iter().all(is_list)
        && rows.iter().all(|r| r.args().len() == rows[0].args().len())
        && !rows[0].args().is_empty()
}

/// `row` written into channel `chn` of `c` from `offset` (< the length),
/// the samples past the end neither written nor read; giac's message for
/// a sample that is not real.
fn write(
    ev: &Eval<'_>,
    c: &mut Clip,
    chn: usize,
    row: &[Expr],
    offset: usize,
) -> Result<(), String> {
    let n = c.len() - offset;
    let depth = c.depth;
    let Some(out) = c.channels.get_mut(chn) else {
        return Err(CHANNEL.to_owned());
    };
    for (p, y) in row.iter().take(n).enumerate() {
        let Some(r) = (!is_list(y)).then(|| real(ev, y)).flatten() else {
            let text = shown(ev, y);
            return Err(format!(
                "Failed to set channel data, {text} is not a valid sample value"
            ));
        };
        out[offset + p] = encode(depth, r.f);
    }
    Ok(())
}

/// `SetChannelData[clip, data, offset, channel]`.
pub(in crate::native) fn set_channel_data(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if !(2..=5).contains(&args.len()) {
        return Ok(None);
    }
    let Some(mut c) = clip(&args[0]) else {
        return warn(ev, e, NOT_CLIP);
    };
    let nc = c.channels.len();
    if !is_list(&args[1]) {
        return warn(ev, e, DATA);
    }
    let rows = args[1].args();
    let matrix = is_matrix(rows);
    if matrix && rows.len() != nc {
        return warn(ev, e, DATA);
    }
    let mut offset = 0;
    if let Some(o) = args.get(2) {
        let Some(o) = o.to_integer() else {
            return warn(ev, e, "Expected an integer");
        };
        match usize::try_from(saturated(&o)).ok().filter(|&o| o < c.len()) {
            Some(o) => offset = o,
            None => return warn(ev, e, "Invalid offset"),
        }
    }
    let mut channel = None;
    if let Some(ch) = args.get(3) {
        let k = if is(ch, Symbol::LEFT, "left") {
            0
        } else if is(ch, Symbol::RIGHT, "right") {
            1
        } else if let Some(n) = ch.to_integer() {
            match usize::try_from(saturated(&n)).ok().filter(|&k| k < nc) {
                Some(k) => k,
                None => return warn(ev, e, CHANNEL),
            }
        } else {
            return warn(ev, e, "Invalid channel specification");
        };
        if matrix {
            return warn(ev, e, DATA);
        }
        channel = Some(k);
    }
    let written = match channel {
        Some(k) => write(ev, &mut c, k, rows, offset),
        None if matrix => rows
            .iter()
            .enumerate()
            .try_for_each(|(k, row)| write(ev, &mut c, k, row.args(), offset)),
        None => (0..nc).try_for_each(|k| write(ev, &mut c, k, rows, offset)),
    };
    match written {
        Ok(()) => Ok(Some(answer(c))),
        Err(message) => warn(ev, e, message),
    }
}

/// `Stereo2Mono[clip]`.
pub(in crate::native) fn stereo2mono(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(c) = (match e.args() {
        [x] => clip(x),
        _ => None,
    }) else {
        return warn(ev, e, BAD_SOUND);
    };
    let nc = double(count(c.channels.len()));
    let mono = (0..c.len())
        .map(|i| {
            let tot = c
                .channels
                .iter()
                .fold(0.0, |tot, ch| tot + decode(c.depth, ch[i]));
            encode(c.depth, tot / nc)
        })
        .collect();
    Ok(Some(answer(Clip {
        rate: c.rate,
        depth: c.depth,
        channels: vec![mono],
    })))
}

/// One clip of a mixdown and its parameters, in samples.
struct Part {
    clip: Clip,
    offset: usize,
    gain: f64,
    pan: f64,
    fade_in: usize,
    fade_out: usize,
}

/// giac's `is_audio_position`: an integer `≥ 0` in samples, a real `≥ 0`
/// in seconds (`floor(1/2 + rate·x)`), saturated.
fn position(ev: &Eval<'_>, x: &Expr, rate: i64) -> Option<usize> {
    if let Some(n) = x.to_integer() {
        return (!n.is_negative()).then(|| usize::try_from(saturated(&n)).unwrap_or(usize::MAX));
    }
    let r = real(ev, x)?;
    (r.f >= 0.0).then(|| usize::try_from(samples(rate, r.f)).unwrap_or(usize::MAX))
}

/// giac's `mix` of `p` into `out` (long enough); whether a sample clipped.
fn mix(out: &mut Clip, p: &Part) -> bool {
    let nc = out.channels.len();
    let onc = p.clip.channels.len();
    let ratio = 10.0_f64.powf(p.gain / 20.0);
    let theta = PI * (p.pan + 1.0) / 4.0;
    let (pl, pr) = (if nc == 1 { 1.0 } else { theta.cos() }, theta.sin());
    let len = p.clip.len();
    let mut clipped = false;
    for i in 0..len {
        let mut r = ratio;
        // std::pow(x, 2) as gcc folds it: x·x.
        if i < p.fade_in {
            let k = (FRAC_PI_2 * (1.0 - double(count(i)) / double(count(p.fade_in)))).cos();
            r *= k * k;
        }
        if i >= len - p.fade_out {
            let k =
                (FRAC_PI_2 * (1.0 - double(count(len - i - 1)) / double(count(p.fade_out)))).cos();
            r *= k * k;
        }
        for c in 0..nc {
            let at = i + p.offset;
            let s1 = decode(out.depth, out.channels[c][at]);
            let gain = if c == 0 { pl } else { pr };
            let s2 = gain * r * decode(p.clip.depth, p.clip.channels[c.min(onc - 1)][i]);
            let s = s1 + s2;
            clipped |= !(-1.0..=1.0).contains(&s);
            out.channels[c][at] = encode(out.depth, s);
        }
    }
    clipped
}

/// `Mixdown[clip₁, parm₁, clip₂, parm₂, …]`.
pub(in crate::native) fn mixdown(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if args.is_empty() || args.len() % 2 != 0 {
        return Ok(None);
    }
    let (mut rate, mut depth, mut nc, mut len) = (0_i64, 0_u32, 0_usize, 0_usize);
    let mut parts = Vec::with_capacity(args.len() / 2);
    for (i, pair) in args.chunks_exact(2).enumerate() {
        let Some(c) = clip(&pair[0]) else {
            return warn(ev, e, format!("Invalid argument {}", 2 * i + 1));
        };
        if c.channels.len() > 2 {
            return warn(ev, e, "Only mono and stereo clips are allowed");
        }
        if rate == 0 {
            rate = c.rate;
        } else if rate != c.rate {
            return warn(ev, e, "All clips must have the same sample rate");
        }
        depth = depth.max(c.depth);
        nc = nc.max(c.channels.len());
        let parm = if is_list(&pair[1]) {
            if pair[1].args().len() > 5 {
                return Ok(None);
            }
            pair[1].args()
        } else {
            std::slice::from_ref(&pair[1])
        };
        let (mut offset, mut gain, mut pan, mut fade_in, mut fade_out) = (0, 0.0, 0.0, 0, 0);
        if let Some(x) = parm.first() {
            let Some(o) = position(ev, x, rate) else {
                return warn(ev, e, "Invalid offset specification");
            };
            offset = o;
        }
        if let Some(x) = parm.get(1) {
            match real(ev, x) {
                Some(r) if r.f <= 0.0 => gain = r.f,
                _ => return warn(ev, e, "Invalid gain specification"),
            }
        }
        if let Some(x) = parm.get(2) {
            match real(ev, x) {
                Some(r) if (-1.0..=1.0).contains(&r.f) => pan = r.f,
                _ => return warn(ev, e, "Invalid pan specification"),
            }
        }
        if let Some(x) = parm.get(3) {
            let Some(f) = position(ev, x, rate) else {
                return warn(ev, e, "Invalid fade-in specification");
            };
            fade_in = f;
        }
        if let Some(x) = parm.get(4) {
            let Some(f) = position(ev, x, rate) else {
                return warn(ev, e, "Invalid fade-out specification");
            };
            fade_out = f;
        }
        if fade_in.saturating_add(fade_out) > c.len() {
            return warn(ev, e, "Fade-in/out is too long");
        }
        len = len.max(offset.saturating_add(c.len()));
        parts.push(Part {
            clip: c,
            offset,
            gain,
            pan,
            fade_in,
            fade_out,
        });
    }
    let terms = ev.context().limits().terms;
    if len.checked_mul(nc).is_none_or(|n| n > terms) {
        return warn(ev, e, TOO_MANY);
    }
    let mut out = Clip {
        rate,
        depth,
        channels: vec![vec![0; len]; nc],
    };
    let mut clipped = false;
    for p in &parts {
        clipped |= mix(&mut out, p);
    }
    if clipped {
        clipping(ev, e);
    }
    Ok(Some(answer(out)))
}

/// giac's cosine similarity of `a` from `offset` and `b` from 0, the
/// geometric mean over the channels (clips of the same parameters).
fn similarity(a: &Clip, b: &Clip, offset: usize) -> f64 {
    let overlap = (a.len() - offset).min(b.len());
    let mut product = 1.0;
    for (ca, cb) in a.channels.iter().zip(&b.channels) {
        let (mut pp, mut na, mut nb) = (0.0, 0.0, 0.0);
        for (&ka, &kb) in ca[offset..offset + overlap].iter().zip(&cb[..overlap]) {
            let (sa, sb) = (decode(a.depth, ka), decode(b.depth, kb));
            pp += sa * sb;
            na += sa * sa;
            nb += sb * sb;
        }
        product *= pp / (na * nb).sqrt();
    }
    product.powf(1.0 / double(count(a.channels.len())))
}

/// `Splice[a, b, n, p]`.
pub(in crate::native) fn splice(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if !(2..=4).contains(&args.len()) {
        return Ok(None);
    }
    let (Some(a), Some(b)) = (clip(&args[0]), clip(&args[1])) else {
        return warn(ev, e, "First two arguments should be audio clips");
    };
    let len = a.len();
    let len_i = i64::try_from(len).unwrap_or(i64::MAX);
    let mut pos = -1_i64;
    if let Some(x) = args.get(2) {
        if let Some(n) = x.to_integer() {
            pos = len_i.saturating_sub(saturated(&n));
        } else if let Some(r) = real(ev, x) {
            if !matches!(x.view(), View::Ball(_)) {
                return Ok(None);
            }
            pos = len_i.saturating_sub(samples(a.rate, r.f));
        } else {
            return warn(ev, e, "Crossfade length should be a real number");
        }
    }
    let mut ex = -1.0;
    if let Some(x) = args.get(3) {
        let Some(r) = real(ev, x) else {
            return warn(ev, e, "Crossfade exponent should be a real number");
        };
        ex = r.f;
    }
    if pos < 0 {
        let da = double(count(len)) / double(a.rate);
        let db = double(count(b.len())) / double(b.rate);
        pos = len_i - samples(a.rate, (da.min(db) / 2.0).min(0.5));
    }
    let Some(pos) = usize::try_from(pos).ok().filter(|&p| p < len) else {
        return warn(ev, e, "Invalid splice position");
    };
    if a.depth != b.depth || a.rate != b.rate || a.channels.len() != b.channels.len() {
        return warn(ev, e, SAME);
    }
    let rl = pos + b.len();
    let terms = ev.context().limits().terms;
    if rl.checked_mul(a.channels.len()).is_none_or(|n| n > terms) {
        return warn(ev, e, TOO_MANY);
    }
    let ex = if ex < 0.0 {
        0.5 + 0.5 * similarity(&a, &b, pos).abs()
    } else {
        ex
    };
    let cflen = (len - pos).min(b.len());
    let pos = len - cflen;
    let mut clipped = false;
    let mut channels = Vec::with_capacity(a.channels.len());
    for (ca, cb) in a.channels.iter().zip(&b.channels) {
        let mut out = vec![0_i64; rl];
        let head = pos.min(rl);
        out[..head].copy_from_slice(&ca[..head]);
        if len < rl {
            let tail = (b.len() - cflen).min(rl - len);
            out[len..len + tail].copy_from_slice(&cb[cflen..cflen + tail]);
        }
        for i in 0..cflen {
            let from = decode(a.depth, ca[i + pos]);
            let to = decode(b.depth, cb[i]);
            let mixed = if ex == 0.0 {
                from + to
            } else {
                let cf = (double(count(i + 1)) / double(count(cflen))).powf(ex);
                from * (1.0 - cf) + to * cf
            };
            clipped |= mixed.abs() > 1.0;
            if let Some(slot) = out.get_mut(i + pos) {
                *slot = encode(a.depth, mixed);
            }
        }
        channels.push(out);
    }
    if clipped {
        clipping(ev, e);
    }
    Ok(Some(answer(Clip {
        rate: a.rate,
        depth: a.depth,
        channels,
    })))
}
