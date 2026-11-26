//! The Hilbert-Huang transform (CAS-648), giac's `hht` with
//! `output=matrix` (signalprocessing.cc:5391-5450 `inst_freq`, `inst_amp`,
//! `hht`; :5674-5791 `_hht`).
//!
//! `HilbertHuangTransform[data, options…]`: `data` a matrix whose rows are
//! the intrinsic mode functions, or a list of numbers decomposed first by
//! [`super::empirical_mode_decomposition`] with giac's parameters (no
//! limit, 2 extrema, threshold 0.1; giac's "No IMF found" when there is
//! none). Per IMF the analytic signal `z` by [`super::hilbert`], the
//! instantaneous frequencies `f_0 = 0`, `f_k = max(0, a arg(conj(z_(k−1))
//! z_k))` (`a = 1`, or `sr/(2π)` with a sample rate) and amplitudes `|z_k|`,
//! in doubles (an angle's imaginary part `+0` when its ball holds 0, as
//! giac's exact zeros give `arg = π` for a negative real). `fmin` and `fmax` are the extreme frequencies unless given;
//! with `ln`, `fmin = max(fmin, MIDI_NOTE_0)` and `r = ln(fmax/fmin)`, else
//! `r = fmax − fmin`; frequency `f` at sample `j` adds its amplitude to bin
//! `⌊(fbins − 1) d/r⌋` of column `j`, `d = ln(f/fmin)` or `f − fmin`, when
//! the bin is in `0 .. fbins − 1`; the answer is the `fbins × n` matrix of
//! the sums (exact 0 elsewhere), all zero when `r` is not positive. The bin
//! is taken in doubles in giac's order (C0085), and the top of the range,
//! `d = r`, lands one bin low as giac's does (observed in
//! `proto/wavelet/giac/p7.out`, not in the source).
//!
//! Options: `output=matrix` (required), `frequencies=a..b`,
//! `samplerate=n`, `ln` or `ln=bool` (needs a sample rate), `bins=n`
//! (`fbins`, 100 by default). giac's errors leave the call with its message
//! as the warning. Departures: giac's `bins=n` gives 161 bins whatever `n`
//! (it stores the option's tag, :5727), ours `n`; its plot and image
//! (`output=plot`, `output=image`, the default) and MIDI note ranges
//! (`frequencies=["C4","C5"]`) stay with a warning, as do its drawing
//! options and a matrix of more than `Limits::terms` entries (C0019).

use cas_num::Integer;
use cas_num::ball::{Ball, ComplexBall, Float};

use super::{dft, empirical_mode_decomposition, hilbert, is_list, list, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::approximation;
use crate::native::statistics::random_variables::{equation, range};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `MIDI_NOTE_0`, the frequency of MIDI note 0 in Hz.
const MIDI_NOTE_0: f64 = 8.175_798_915_64;

/// giac's frequency bins for a matrix.
const DEFAULT_BINS: usize = 100;

fn hold(ev: &mut Eval<'_>, e: &Expr, message: &'static str) -> Native {
    warn(ev, e, message);
    Ok(None)
}

fn named(e: &Expr, name: &str) -> bool {
    e.as_symbol().is_some_and(|s| s.name() == name)
}

/// A real number as a double: `None` for anything else.
fn double(ev: &Eval<'_>, x: &Expr) -> Option<f64> {
    if !x.is_numeric() {
        return None;
    }
    let n = approximation(ev.context(), x, 64).ok()?;
    let z = dft::ball(&n, 64)?;
    (z.im().is_exact() && *z.im() == Ball::ZERO)
        .then(|| z.re().mid().to_f64())
        .filter(|v| v.is_finite())
}

/// A positive integer option's value.
fn positive(v: &Expr) -> Option<usize> {
    v.as_i64()
        .filter(|n| *n >= 1 && *n <= i64::from(i32::MAX))
        .and_then(|n| usize::try_from(n).ok())
}

/// The options giac's matrix takes.
struct Options {
    /// `None` for an end that is not a real number: giac's range then is
    /// not positive.
    frequencies: Option<(Option<f64>, Option<f64>)>,
    rate: u32,
    log: bool,
    bins: usize,
}

/// The options, or giac's message for one it refuses, or ours for one that
/// stays.
fn options(ev: &Eval<'_>, given: &[Expr]) -> Result<Options, &'static str> {
    let mut o = Options {
        frequencies: None,
        rate: 0,
        log: false,
        bins: DEFAULT_BINS,
    };
    let mut matrix = false;
    for option in given {
        if option.as_symbol() == Some(Symbol::LOG) {
            o.log = true;
            continue;
        }
        let Some((k, v)) = equation(option) else {
            return Err("Unknown optional argument");
        };
        match k.as_symbol() {
            Some(Symbol::FREQUENCIES) => {
                if is_list(v) && v.args().len() == 2 {
                    return Err("MIDI note ranges stay: not built");
                }
                let Some((a, b)) = range(v) else {
                    return Err("Invalid frequency range specification");
                };
                o.frequencies = Some((double(ev, a), double(ev, b)));
            }
            Some(Symbol::SAMPLE_RATE) => {
                o.rate = positive(v)
                    .and_then(|n| u32::try_from(n).ok())
                    .ok_or("Sample rate must be a positive integer")?;
            }
            Some(Symbol::BINS) => {
                o.bins =
                    positive(v).ok_or("Number of frequency bins must be a positive integer")?;
            }
            Some(Symbol::LOG) => {
                o.log = match v.as_symbol() {
                    Some(Symbol::TRUE) => true,
                    Some(Symbol::FALSE) => false,
                    _ => return Err("Expected a boolean value"),
                };
            }
            _ if named(k, "output") => {
                if v.as_symbol() == Some(Symbol::MATRIX) || named(v, "matrix") {
                    matrix = true;
                } else if named(v, "plot")
                    || named(v, "image")
                    || v.as_symbol() == Some(Symbol::COLUMN_SPACE)
                {
                    matrix = false;
                } else {
                    return Err("Unknown output type");
                }
            }
            _ => return Err("Unknown optional argument"),
        }
    }
    if o.log && o.rate == 0 {
        return Err("Log-frequency requires sample rate");
    }
    if !matrix {
        return Err("the Hilbert spectrum's plot and image stay: output=matrix gives its matrix");
    }
    Ok(o)
}

/// The rows of a matrix: a nonempty list of lists of one length.
fn rows(data: &Expr) -> Option<&[Expr]> {
    let rows = data.args();
    let n = rows.first()?.args().len();
    (is_list(data) && n > 0 && rows.iter().all(|r| is_list(r) && r.args().len() == n))
        .then_some(rows)
}

/// `arg(conj(u) v)`, its imaginary part `+0` when its ball holds 0: giac's
/// exact zero there (its `arg` of a negative real is π).
fn step_angle(u: &ComplexBall, v: &ComplexBall, prec: u64) -> f64 {
    let p = u.conj().mul(v, prec);
    let im = if p.im().contains_integer(&Integer::ZERO) {
        0.0
    } else {
        p.im().mid().to_f64()
    };
    im.atan2(p.re().mid().to_f64())
}

/// `HilbertHuangTransform[data, options…]`, giac's `hht`.
#[allow(clippy::too_many_lines, reason = "giac's _hht and hht, in order")]
pub(in crate::native) fn hilbert_huang_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [data, given @ ..] = e.args() else {
        return Ok(None);
    };
    if !is_list(data) {
        return hold(ev, e, "Bad Argument Type");
    }
    let o = match options(ev, given) {
        Ok(o) => o,
        Err(message) => return hold(ev, e, message),
    };
    let imfs: Vec<Expr> = if let Some(rows) = rows(data) {
        rows.to_vec()
    } else {
        let emd = Expr::call(Symbol::EMPIRICAL_MODE_DECOMPOSITION, [data.clone()]);
        let Some(found) = empirical_mode_decomposition(ev, &emd)? else {
            return Ok(None);
        };
        let imfs = found
            .args()
            .first()
            .map(|l| l.args().to_vec())
            .unwrap_or_default();
        if imfs.is_empty() {
            return hold(ev, e, "No IMF found");
        }
        imfs
    };
    let n = imfs[0].args().len();
    if n == 0 {
        return Ok(None);
    }
    let fbins = o.bins;
    let limit = ev.context().limits().terms;
    if fbins.saturating_mul(n) > limit {
        warn(
            ev,
            e,
            format!("a matrix of more than {limit} entries (Limits::terms)"),
        );
        return Ok(None);
    }
    let a = if o.rate == 0 {
        1.0
    } else {
        std::f64::consts::FRAC_1_PI * 0.5 * f64::from(o.rate)
    };
    let prec = ev.context().working_precision();
    let mut f: Vec<Vec<f64>> = Vec::with_capacity(imfs.len());
    let mut amp: Vec<Vec<f64>> = Vec::with_capacity(imfs.len());
    for imf in &imfs {
        let call = Expr::call(Symbol::HILBERT_MATRIX, [imf.clone()]);
        let Some(z) = hilbert(ev, &call, imf)? else {
            return Ok(None);
        };
        let Some(z) = z
            .args()
            .iter()
            .map(|v| dft::ball(v, prec))
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(None);
        };
        let mut fi = Vec::with_capacity(n);
        fi.push(0.0);
        for w in z.windows(2) {
            fi.push(f64::max(0.0, a * step_angle(&w[0], &w[1], prec)));
        }
        f.push(fi);
        amp.push(
            z.iter()
                .map(|v| v.re().mid().to_f64().hypot(v.im().mid().to_f64()))
                .collect(),
        );
    }
    let lowest = || f.iter().flatten().copied().fold(f64::INFINITY, f64::min);
    let highest = || {
        f.iter()
            .flatten()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let (fmin, fmax) = match o.frequencies {
        None => (Some(lowest()), Some(highest())),
        Some(ends) => ends,
    };
    let mut sums: Vec<Option<f64>> = vec![None; fbins * n];
    if let (Some(fmin), Some(fmax)) = (fmin, fmax) {
        let fmin = if o.log { MIDI_NOTE_0.max(fmin) } else { fmin };
        let r = if o.log {
            (fmax / fmin).ln()
        } else {
            fmax - fmin
        };
        if r > 0.0 {
            #[allow(clippy::cast_precision_loss, reason = "fbins ≤ Limits::terms")]
            let top = (fbins - 1) as f64;
            for j in 0..n {
                for (fi, ai) in f.iter().zip(&amp) {
                    let d = if o.log {
                        (fi[j] / fmin).ln()
                    } else {
                        fi[j] - fmin
                    };
                    let mut h = (top * d / r).floor();
                    #[allow(clippy::float_cmp, reason = "giac's top of the range, exactly")]
                    if d == r && fbins >= 2 {
                        h = top - 1.0;
                    }
                    if h.is_finite() && h >= 0.0 && h < top + 1.0 {
                        #[allow(
                            clippy::cast_possible_truncation,
                            clippy::cast_sign_loss,
                            reason = "a bin in 0 .. fbins − 1"
                        )]
                        let cell = &mut sums[h as usize * n + j];
                        *cell = Some(cell.unwrap_or(0.0) + ai[j]);
                    }
                }
            }
        }
    }
    let number =
        |v: f64| Float::from_f64(v).map_or_else(|| Expr::from(0), |x| Expr::from(Ball::exact(x)));
    let result = list(
        sums.chunks(n)
            .map(|row| {
                list(
                    row.iter()
                        .map(|v| v.map_or_else(|| Expr::from(0), number))
                        .collect(),
                )
            })
            .collect(),
    );
    ev.step("hht", e, &result);
    Ok(Some(result))
}
