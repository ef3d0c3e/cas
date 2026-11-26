//! Short-time Fourier transforms (CAS-648), giac's `stft` and `istft`
//! (signalprocessing.cc:754-860, `stft`, `istft`, `_stft`, `_istft`).
//!
//! `ShortTimeFourierTransform[x]`, `[x, m]`, `[x, m, window]`: `x` a list of
//! numbers (the imaginary parts dropped with giac's warning "ignoring
//! imaginary part of input"), `m` a power of 2 (128 by default), `window` a
//! window function (`HammingWindow` by default). `x` is padded with zeros
//! to a multiple of `m`; frame `i < n/m − 1` packs the `2m` samples from
//! `i m` as `v_j = x_(i m + 2j) + i x_(i m + 2j + 1)`, windowed by the
//! window's weights on `m` samples, transformed (`u = fft(v)`), and
//! untangled into `r_j = a_j u_j + b_j conj(u_((m − j) mod m))` with `w_j =
//! i e^(−iπ j/m)`, `a_j = (1 − w_j)/2`, `b_j = (1 + w_j)/2`: one row per
//! frame, `[]` when there is none.
//!
//! `InverseShortTimeFourierTransform[y]`, `y` a numeric matrix of `N` rows
//! of `m`: with `w_j = i e^(iπ j/m)`, `a_j = (1 + w_j)/2`, `b_j = (1 −
//! w_j)/2`, row `i` gives `c_j = a_j y_ij + b_j conj(y_(i,(m − j) mod m))`,
//! `u = ifft(c)`, whose real and imaginary parts are added to `x_(i m +
//! 2j)` and `x_(i m + 2j + 1)`: `(N + 1) m` samples, overlapped and added
//! without normalising (giac's).
//!
//! Through the FFT core of [`super::dft`] in complex balls at the working
//! precision and a guard; the answers are numbers. giac's errors leave the
//! call with giac's message as the warning. Departures: the window is one
//! of the fifteen window functions, where giac takes any function (its
//! `stft(x, 2, sin)` applies `sin` to the frame); a window of more than
//! 2048 samples, and a padded list longer than `Limits::terms`, stay with a
//! warning (C0019).

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Integer, Rational};

use super::{dft, is_list, list, transform, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::approximation;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's default window size.
const DEFAULT_SIZE: usize = 128;

/// The longest window whose weights are evaluated (C0019).
const MAX_WINDOW: usize = dft::MAX_DIRECT;

/// The window functions `stft` takes.
const WINDOWS: [Symbol; 15] = [
    Symbol::BARTLETT_HANN_WINDOW,
    Symbol::BLACKMAN_HARRIS_WINDOW,
    Symbol::BLACKMAN_WINDOW,
    Symbol::BOHMAN_WINDOW,
    Symbol::COSINE_WINDOW,
    Symbol::GAUSSIAN_WINDOW,
    Symbol::HAMMING_WINDOW,
    Symbol::HANN_POISSON_WINDOW,
    Symbol::HANN_WINDOW,
    Symbol::PARZEN_WINDOW,
    Symbol::POISSON_WINDOW,
    Symbol::RIEMANN_WINDOW,
    Symbol::TRIANGLE_WINDOW,
    Symbol::TUKEY_WINDOW,
    Symbol::WELCH_WINDOW,
];

fn hold(ev: &mut Eval<'_>, e: &Expr, message: &'static str) -> Native {
    warn(ev, e, message);
    Ok(None)
}

/// The guard bits of a transform of length `m`, as `hilbert`'s.
fn precision(ev: &Eval<'_>, m: usize) -> u64 {
    ev.context().working_precision() + 2 * u64::from(usize::BITS - m.leading_zeros()) + 16
}

/// `(1 − w_j)/2` and `(1 + w_j)/2` for `w_j = i e^(∓iπ j/m)`, the sign `−`
/// for the transform: `w_j = sin(π j/m) + i cos(π j/m)`, or `−sin(π j/m) +
/// i cos(π j/m)`.
fn untangling(m: usize, inverse: bool, prec: u64) -> Vec<(ComplexBall, ComplexBall)> {
    let one = Ball::from_rational(&Rational::from(Integer::ONE), prec);
    let half = Integer::from(-1);
    let m_i = Integer::from(i64::try_from(m).expect("a length"));
    (0..m)
        .map(|j| {
            let t = Ball::from_rational(
                &Rational::new(
                    Integer::from(i64::try_from(j).expect("an index")),
                    m_i.clone(),
                ),
                prec,
            );
            let (s, c) = (t.sin_pi(prec), t.cos_pi(prec));
            let s = if inverse { s.neg() } else { s };
            // (1 − w)/2 and (1 + w)/2.
            let minus =
                ComplexBall::new(one.sub(&s, prec).mul_2exp(&half), c.neg().mul_2exp(&half));
            let plus = ComplexBall::new(one.add(&s, prec).mul_2exp(&half), c.mul_2exp(&half));
            (minus, plus)
        })
        .collect()
}

/// The weights of `window` on `m` samples, as balls.
fn weights(
    ev: &mut Eval<'_>,
    window: Symbol,
    m: usize,
    prec: u64,
) -> Result<Option<Vec<Ball>>, EvalError> {
    let ones = list(vec![Expr::from(1); m]);
    let found = ev.eval(&Expr::call(window, [ones]))?;
    if !is_list(&found) || found.args().len() != m {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(m);
    for w in found.args() {
        let Ok(n) = approximation(ev.context(), w, prec) else {
            return Ok(None);
        };
        match dft::ball(&n, prec) {
            Some(z) => out.push(z.re().clone()),
            None => return Ok(None),
        }
    }
    Ok(Some(out))
}

/// `ShortTimeFourierTransform[x, m, window]`, giac's `stft`.
pub(in crate::native) fn short_time_fourier_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if args.is_empty() {
        return Ok(None);
    }
    if args.len() > 3 {
        return hold(ev, e, "Too many input arguments");
    }
    let x = &args[0];
    if !is_list(x) {
        return hold(ev, e, "Bad Argument Type");
    }
    let window = match args.get(2) {
        None => Symbol::HAMMING_WINDOW,
        Some(w) => match w.as_symbol() {
            Some(s) if WINDOWS.contains(&s) => s,
            _ => return hold(ev, e, "Expected a window function"),
        },
    };
    let m = match args.get(1) {
        None => DEFAULT_SIZE,
        Some(s) => {
            let Some(m) = s
                .as_i64()
                .filter(|m| *m >= 1)
                .and_then(|m| usize::try_from(m).ok())
            else {
                return hold(ev, e, "Expected a positive integer");
            };
            if !m.is_power_of_two() {
                return hold(ev, e, "Window size must be a power of 2");
            }
            m
        }
    };
    let samples = x.args();
    if samples.is_empty() {
        return hold(ev, e, "Invalid dimension");
    }
    let prec = precision(ev, m);
    let Some(z) = samples
        .iter()
        .map(|v| dft::ball(v, prec))
        .collect::<Option<Vec<_>>>()
    else {
        return hold(ev, e, "Bad Argument Value");
    };
    if z.iter()
        .any(|v| !(v.im().is_exact() && *v.im() == Ball::ZERO))
    {
        warn(ev, e, "ignoring imaginary part of input");
    }
    let padded = samples.len().div_ceil(m).saturating_mul(m);
    let limit = ev.context().limits().terms;
    if padded > limit {
        warn(ev, e, format!("more than {limit} samples (Limits::terms)"));
        return Ok(None);
    }
    let frames = padded / m - 1;
    if frames == 0 {
        let result = list(Vec::new());
        ev.step("stft", e, &result);
        return Ok(Some(result));
    }
    if m < 2 {
        return hold(ev, e, "Window size should be at least 2");
    }
    if m > MAX_WINDOW {
        warn(ev, e, format!("a window of more than {MAX_WINDOW} samples"));
        return Ok(None);
    }
    let Some(w) = weights(ev, window, m, prec)? else {
        return Ok(None);
    };
    let mut real: Vec<Ball> = z.iter().map(|v| v.re().clone()).collect();
    real.resize(padded, Ball::ZERO);
    let ab = untangling(m, false, prec);
    let mut rows = Vec::with_capacity(frames);
    for i in 0..frames {
        let v: Vec<ComplexBall> = (0..m)
            .map(|j| {
                let k = i * m + 2 * j;
                ComplexBall::new(real[k].mul(&w[j], prec), real[k + 1].mul(&w[j], prec))
            })
            .collect();
        let u = transform(&v, false, prec);
        let r = (0..m)
            .map(|j| {
                let (a, b) = &ab[j];
                let found = a
                    .mul(&u[j], prec)
                    .add(&b.mul(&u[(m - j) % m].conj(), prec), prec);
                dft::number(found)
            })
            .collect();
        rows.push(list(r));
    }
    let result = list(rows);
    ev.step("stft", e, &result);
    Ok(Some(result))
}

/// `InverseShortTimeFourierTransform[y]`, giac's `istft`.
pub(in crate::native) fn inverse_short_time_fourier_transform(
    ev: &mut Eval<'_>,
    e: &Expr,
) -> Native {
    let [y] = e.args() else {
        return Ok(None);
    };
    let rows = y.args();
    let m = rows.first().map_or(0, |r| r.args().len());
    if !is_list(y) || m == 0 || rows.iter().any(|r| !is_list(r) || r.args().len() != m) {
        return hold(ev, e, "Bad Argument Type");
    }
    let n = (rows.len() + 1).saturating_mul(m);
    let limit = ev.context().limits().terms;
    if n > limit {
        warn(ev, e, format!("more than {limit} samples (Limits::terms)"));
        return Ok(None);
    }
    if !(m <= dft::MAX_DIRECT || (m.is_power_of_two() && m <= dft::MAX_RADIX_TWO)) {
        warn(
            ev,
            e,
            format!("rows of {m}, not a power of 2, above {}", dft::MAX_DIRECT),
        );
        return Ok(None);
    }
    let prec = precision(ev, m);
    let mut matrix = Vec::with_capacity(rows.len());
    for r in rows {
        let Some(z) = r
            .args()
            .iter()
            .map(|v| dft::ball(v, prec))
            .collect::<Option<Vec<_>>>()
        else {
            return hold(ev, e, "Bad Argument Value");
        };
        matrix.push(z);
    }
    let ab = untangling(m, true, prec);
    let mut x = vec![Ball::ZERO; n];
    for (i, r) in matrix.iter().enumerate() {
        // (1 + w)/2 and (1 − w)/2: `untangling`'s pair, swapped.
        let c: Vec<ComplexBall> = (0..m)
            .map(|j| {
                let (b, a) = &ab[j];
                a.mul(&r[j], prec)
                    .add(&b.mul(&r[(m - j) % m].conj(), prec), prec)
            })
            .collect();
        let u = transform(&c, true, prec);
        for (j, v) in u.iter().enumerate() {
            let k = i * m + 2 * j;
            x[k] = x[k].add(v.re(), prec);
            x[k + 1] = x[k + 1].add(v.im(), prec);
        }
    }
    let result = list(x.into_iter().map(Expr::from).collect());
    ev.step("istft", e, &result);
    Ok(Some(result))
}
