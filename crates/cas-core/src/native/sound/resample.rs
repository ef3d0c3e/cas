//! Resampling (CAS-633): giac's `resample(clip[, rate[, quality]])`
//! (giac signalprocessing.cc `_resample` :2859-2881, `audio_clip::resample`
//! :1981-2028), the rate 44100 and the quality 2 by default. The answer is
//! a clip at the new rate `r`, of the same depth and channels, of
//! `floor(len·r/sr)` samples, giac's count; a count of 0 stays with a
//! warning (giac's would ask libsamplerate for no frames). giac hands the
//! work to libsamplerate, which its Debian build lacks ("libsamplerate is
//! required for resampling audio"); ours computes it, with libsamplerate's
//! methods for the qualities 0 and 1 and a design of ours for 2 to 4.
//!
//! The output sample `k` sits at the input position `t = k·sr/r`, exact.
//! Quality 0 is libsamplerate's linear converter (`src_linear.c`), one frame
//! late: `y_k = X(t − 1)`, `X(s) = x_0` for `s ≤ 0`, else `x_i + f·(x_{i+1} −
//! x_i)`, `i = floor(s)`, `f = s − i`. Quality 1 is its zero-order hold
//! (`src_zoh.c`): `y_k = x_{max(0, floor(t) − 1)}`. Both are exact on the
//! decoded samples `k/2^(b−1)` and encoded exactly, `floor(range·v + 1/2)`
//! with `range = 2^(b−1) − 1`, so a decoded top code encodes one lower as
//! in giac's other edits. Positions past the last frame hold it (a clip of
//! one sample gives a constant), where libsamplerate gives no output there.
//!
//! Qualities 2, 3 and 4 are band-limited interpolation (J. O. Smith's
//! resample method), centred, without libsamplerate's lag and not its
//! tables: `y_k = Σ_j x_j·h(t − j)` over the input (zero outside it), `h(u) =
//! s·g(s·u)`, `s = c·min(1, r/sr)`, `g(x) = sinc(x)·I0(β√(1 −
//! (x/H)²))/I0(β)` for `|x| < H` (a Kaiser window, `β = 0.1102(A − 8.7)`),
//! `g` tabulated on `[0, H]` at `L` points per unit, once per quality, and
//! interpolated linearly; computed in `f64` and encoded as `CreateWav`
//! does. The designs, cutoff `c`, attenuation `A`, zero crossings `H`,
//! points per unit `L`: quality 2 `0.90`, 100 dB, 37, 512; quality 3 `0.95`,
//! 120 dB, 86, 2048; quality 4 `0.98`, 144 dB, 312, 8192. Measured on the
//! transform of the interpolated table (.run/1.7/T8/resample.py, its
//! images included), with `δ = 10^(−A/20)`: the passband deviation
//! `max |G − 1|` up to 80, 90 and 96 % of the lower Nyquist frequency is
//! `8.5e−6`, `9.2e−7` and `5.4e−8` (0.85, 0.92 and 0.86 `δ`); the stopband
//! from that Nyquist frequency on is attenuated 102.0, 121.7 and 146.6 dB;
//! the table's error is `1.57e−6`, `9.81e−8` and `6.13e−9`, 1.002, 1.0005
//! and 1.0000 times `π²/(24L²)` (the window's curvature adds to the sinc's
//! `π²/3`). The tables take 0.15, 1.4 and 20 MB. The designs first set,
//! `H = 29, 74, 232` and `L = 256, 1024, 4096`, missed: linear interpolation's droop
//! `sinc²(f/L)` alone is about `δ` at the passband edge with those `L`
//! (passband 1.9, 1.9 and 1.7 `δ`), and the transitions were too wide for
//! `A` at the Nyquist frequency (119.0 dB for the quality 3, 1.05 `δ` for
//! the quality 4).
//!
//! giac's errors stay with its message ("Invalid sample rate
//! specification" for a rate that is not an integer `≥ 1`, "Invalid quality
//! specification" for a quality that is not an integer from 0 to 4,
//! "Invalid sound data" for a non-clip), checked in its order; more than
//! three arguments stay silently. Sizes are checked before anything is
//! built (C0019): the rate and the result's samples against
//! `Limits::terms`, and the convolution's work, taps times outputs times
//! channels, against 1024 times it; the evaluation's interrupt is polled.

#![allow(
    clippy::many_single_char_names,
    reason = "x, h, t, s, k mirror the notation above"
)]

use std::f64::consts::PI;
use std::sync::OnceLock;
use std::sync::atomic::Ordering;

use super::edits::answer;
use super::{BAD_SOUND, Clip, Native, clip, encode, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;

const RATE: &str = "Invalid sample rate specification";
const QUALITY: &str = "Invalid quality specification";
const TOO_MANY: &str = "Too many samples";
const EMPTY: &str = "The resampled clip would have no samples";

/// The convolution's work allowed per term of `Limits::terms`.
const WORK_PER_TERM: f64 = 1024.0;

/// A windowed-sinc design.
struct Design {
    /// The cutoff, a fraction of the lower Nyquist frequency.
    cutoff: f64,
    /// The stopband attenuation in dB, which sets the window's `β`.
    attenuation: f64,
    /// The zero crossings `H` of each side of the kernel.
    zeros: u32,
    /// The table's points per unit `L`.
    per_unit: u32,
}

/// Qualities 2, 3 and 4.
const DESIGNS: [Design; 3] = [
    Design {
        cutoff: 0.90,
        attenuation: 100.0,
        zeros: 37,     // was 29
        per_unit: 512, // was 256
    },
    Design {
        cutoff: 0.95,
        attenuation: 120.0,
        zeros: 86,      // was 74
        per_unit: 2048, // was 1024
    },
    Design {
        cutoff: 0.98,
        attenuation: 144.0,
        zeros: 312,     // was 232
        per_unit: 8192, // was 4096
    },
];

static TABLES: [OnceLock<Vec<f64>>; 3] = [OnceLock::new(), OnceLock::new(), OnceLock::new()];

/// `I0(x)` by its series.
fn bessel_i0(x: f64) -> f64 {
    let y = x * x / 4.0;
    let (mut sum, mut term, mut k) = (1.0_f64, 1.0_f64, 1.0_f64);
    loop {
        term *= y / (k * k);
        sum += term;
        if term < 1e-17 * sum {
            return sum;
        }
        k += 1.0;
    }
}

/// `g` at `i/L`, `i = 0 … H·L`, then a zero (`g(H) = 0`) and a guard zero.
fn table(quality: usize) -> &'static [f64] {
    TABLES[quality].get_or_init(|| {
        let d = &DESIGNS[quality];
        let beta = 0.1102 * (d.attenuation - 8.7);
        let (h, l) = (f64::from(d.zeros), f64::from(d.per_unit));
        let scale = bessel_i0(beta);
        let mut t: Vec<f64> = (0..d.zeros * d.per_unit)
            .map(|i| {
                if i == 0 {
                    return 1.0;
                }
                let x = f64::from(i) / l;
                let r = x / h;
                (PI * x).sin() / (PI * x) * bessel_i0(beta * (1.0 - r * r).sqrt()) / scale
            })
            .collect();
        t.extend([0.0, 0.0]);
        t
    })
}

/// The linearly interpolated `g(x)`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "0 ≤ p < H·L, checked: the floor indexes the table"
)]
fn kernel(t: &[f64], per_unit: f64, end: f64, x: f64) -> f64 {
    let p = x.abs() * per_unit;
    let i = p.floor();
    if i >= end || i.is_nan() {
        return 0.0;
    }
    let k = i as usize;
    t[k] + (p - i) * (t[k + 1] - t[k])
}

/// libsamplerate's linear (quality 0) or zero-order hold (quality 1)
/// converter, exactly, `n` samples.
fn held(c: &Clip, rate: i64, linear: bool, n: usize) -> Vec<Vec<i64>> {
    let (sr, nsr) = (i128::from(c.rate), i128::from(rate));
    let last = c.len() - 1;
    let one = 1_i128 << (c.depth - 1);
    let range = one - 1;
    c.channels
        .iter()
        .map(|ch| {
            let x = |j: i128| i128::from(ch[usize::try_from(j).map_or(last, |j| j.min(last))]);
            (0..n)
                .map(|k| {
                    let k = i128::try_from(k).unwrap_or(i128::MAX);
                    // The sample is `v = num/(nsr·one)`.
                    let num = if linear {
                        let s = k * sr - nsr;
                        if s <= 0 {
                            x(0) * nsr
                        } else {
                            let (i, f) = (s / nsr, s % nsr);
                            x(i) * nsr + f * (x(i + 1) - x(i))
                        }
                    } else {
                        x((k * sr / nsr - 1).max(0)) * nsr
                    };
                    let code = (2 * range * num + nsr * one).div_euclid(2 * nsr * one);
                    i64::try_from(code).unwrap_or_default()
                })
                .collect()
        })
        .collect()
}

/// `x` as a double.
#[allow(
    clippy::cast_precision_loss,
    reason = "positions and counts as doubles, as intended"
)]
fn double(x: i128) -> f64 {
    x as f64
}

/// `x`, a nonnegative double, as an index.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a saturating float-to-int conversion, as intended"
)]
fn index(x: f64) -> usize {
    x as usize
}

/// The band-limited interpolation of design `d`, `n` samples.
fn filtered(
    ev: &Eval<'_>,
    c: &Clip,
    rate: i64,
    quality: usize,
    n: usize,
) -> Result<Vec<Vec<i64>>, EvalError> {
    let d = &DESIGNS[quality];
    let t = table(quality);
    let (sr, nsr) = (i128::from(c.rate), i128::from(rate));
    let s = d.cutoff * (double(nsr) / double(sr)).min(1.0);
    let half = f64::from(d.zeros) / s;
    let (per_unit, end) = (f64::from(d.per_unit), f64::from(d.zeros * d.per_unit));
    let one = f64::from(1_u32 << (c.depth - 1));
    let x: Vec<Vec<f64>> = c
        .channels
        .iter()
        .map(|ch| ch.iter().map(|&k| double(k.into()) / one).collect())
        .collect();
    let last = c.len() - 1;
    let mut out = vec![Vec::with_capacity(n); c.channels.len()];
    let mut w = Vec::new();
    for k in 0..n {
        if k % 1024 == 0 && ev.context().interrupt_flag().load(Ordering::Relaxed) {
            return Err(EvalError::Interrupted);
        }
        let num = i128::try_from(k).unwrap_or(i128::MAX) * sr;
        let pos = double(num) / double(nsr);
        let lo = index((pos - half).ceil().max(0.0));
        let hi = index((pos + half).floor()).min(last);
        w.clear();
        w.extend((lo..=hi).map(|j| {
            let j = i128::try_from(j).unwrap_or(i128::MAX);
            let u = double(num - j * nsr) / double(nsr);
            s * kernel(t, per_unit, end, s * u)
        }));
        for (ch, o) in x.iter().zip(&mut out) {
            let y = ch[lo..=hi].iter().zip(&w).fold(0.0, |y, (a, b)| y + a * b);
            o.push(encode(c.depth, y));
        }
    }
    Ok(out)
}

/// `Resample[clip, rate, quality]`.
pub(in crate::native) fn resample(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, rate, quality) = match e.args() {
        [x] => (x, None, None),
        [x, r] => (x, Some(r), None),
        [x, r, q] => (x, Some(r), Some(q)),
        _ => return Ok(None),
    };
    let rate = match rate.map(|r| r.as_i64().filter(|&r| r >= 1)) {
        None => 44100,
        Some(Some(r)) => r,
        Some(None) => return warn(ev, e, RATE),
    };
    let quality = match quality.map(|q| q.as_i64().filter(|q| (0..=4).contains(q))) {
        None => 2,
        Some(Some(q)) => q,
        Some(None) => return warn(ev, e, QUALITY),
    };
    let Some(c) = clip(x) else {
        return warn(ev, e, BAD_SOUND);
    };
    let terms = double(i128::try_from(ev.context().limits().terms).unwrap_or(i128::MAX));
    let nc = c.channels.len();
    let len = i128::try_from(c.len()).unwrap_or(i128::MAX);
    if double(rate.into()) > terms {
        return warn(ev, e, TOO_MANY);
    }
    let n = len * i128::from(rate) / i128::from(c.rate);
    if n < 1 {
        return warn(ev, e, EMPTY);
    }
    let samples = double(n) * double(i128::try_from(nc).unwrap_or(i128::MAX));
    let Ok(n) = usize::try_from(n) else {
        return warn(ev, e, TOO_MANY);
    };
    if samples > terms {
        return warn(ev, e, TOO_MANY);
    }
    let channels = if quality <= 1 {
        held(&c, rate, quality == 0, n)
    } else {
        let q = usize::try_from(quality - 2).unwrap_or_default();
        let d = &DESIGNS[q];
        let s = d.cutoff * (double(rate.into()) / double(c.rate.into())).min(1.0);
        let taps = (2.0 * (f64::from(d.zeros) / s).ceil() + 1.0).min(double(len));
        if taps * samples > WORK_PER_TERM * terms {
            return warn(ev, e, TOO_MANY);
        }
        filtered(ev, &c, rate, q, n)?
    };
    Ok(Some(answer(Clip {
        rate,
        depth: c.depth,
        channels,
    })))
}
