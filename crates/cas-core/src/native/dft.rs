//! The discrete Fourier transform and convolutions (CAS-DOC-45 §4,
//! CAS-465): Xcas's `fft`, `ifft` and `convolution`.
//!
//! - `DiscreteFourierTransform[l]`: `X_k = Σ_j x_j ω^(jk)`, `ω = e^(−2πi/n)`,
//!   giac's sign, for a list of numbers; `InverseDiscreteFourierTransform`
//!   with `ω⁻¹` and `1/n`. In complex balls at the working precision and a
//!   guard of `2 log₂ n + 16` bits, so each entry holds the exact
//!   transform: radix 2 for a power of 2 (Cooley and Tukey, iteratively),
//!   the direct sum up to [`MAX_DIRECT`] entries otherwise, `ω^(jk)` as
//!   `ω^(jk mod n)` from a table of `cos πt` and `sin πt` at `t = 2m/n`. An
//!   entry with a zero imaginary part that is exact is a real ball, as giac's
//!   `fft([1,2,3,4])` is `[10.0,-2.0+2.0*i,-2.0,-2.0-2.0*i]`.
//! - `DiscreteFourierTransform[l, w, p]`, giac's `fft(l, w, p)`: integer
//!   entries, `w` of order `n` modulo `p` (checked), the sums modulo `p` in
//!   `0, …, p − 1`, giac's `fft([1,2,3,4],2,5)` being `[0,2,3,4]`; the
//!   inverse with `w⁻¹` and `n⁻¹`. Another `w` leaves the call, with a
//!   warning.
//! - `Convolution[a, b]`: `c_k = Σ a_i b_(k−i)`, built and evaluated, so
//!   exact for exact entries (`[1,1,1,-3]` for giac's
//!   `convolution([1,2,3],[1,-1])`, which giac computes in decimals) and
//!   symbols alike, up to `Limits::terms` products.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of transforms: x of length n to X with roots w, butterflies t and u"
)]

use cas_num::ball::{Ball, ComplexBall};
use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The longest list transformed by the direct sum: `n²` products.
const MAX_DIRECT: usize = 1 << 11;

/// The longest list transformed modulo a prime, by the direct sum.
const MAX_MODULAR: usize = 1 << 12;

/// The longest power of 2 transformed.
const MAX_RADIX_TWO: usize = 1 << 20;

fn warn(ev: &mut Eval<'_>, e: &Expr, message: String) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// A number as a complex ball at `prec` bits.
fn ball(x: &Expr, prec: u64) -> Option<ComplexBall> {
    match x.view() {
        View::Ball(b) => Some(ComplexBall::from_real(b.clone())),
        View::ComplexBall(z) => Some(z.clone()),
        _ => x.to_complex().map(|z| ComplexBall::from_complex(&z, prec)),
    }
}

/// A complex ball as an expression, real when its imaginary part is an
/// exact 0.
fn number(z: ComplexBall) -> Expr {
    if z.im().is_exact() && *z.im() == Ball::ZERO {
        Expr::from(z.re().clone())
    } else {
        Expr::from(z)
    }
}

/// `ω^m = e^(∓2πi m/n)`, the sign `−` for the transform, for `m < n`.
fn roots(n: usize, inverse: bool, prec: u64) -> Vec<ComplexBall> {
    let n_i = i64::try_from(n).expect("a length");
    (0..n)
        .map(|m| {
            let t = Ball::from_rational(
                &Rational::new(
                    Integer::from(2 * i64::try_from(m).expect("an index")),
                    Integer::from(n_i),
                ),
                prec,
            );
            let (c, s) = (t.cos_pi(prec), t.sin_pi(prec));
            ComplexBall::new(c, if inverse { s } else { s.neg() })
        })
        .collect()
}

/// The transform of `x` with the roots `w`, by the direct sum.
fn direct(x: &[ComplexBall], w: &[ComplexBall], prec: u64) -> Vec<ComplexBall> {
    let n = x.len();
    (0..n)
        .map(|k| {
            x.iter()
                .enumerate()
                .fold(ComplexBall::from_real(Ball::ZERO), |sum, (j, xj)| {
                    sum.add(&xj.mul(&w[(j * k) % n], prec), prec)
                })
        })
        .collect()
}

/// The transform of `x`, of a length a power of 2, with the roots `w`:
/// bit-reversed order, then the butterflies.
fn radix_two(x: &[ComplexBall], w: &[ComplexBall], prec: u64) -> Vec<ComplexBall> {
    let n = x.len();
    let bits = n.trailing_zeros();
    let mut a: Vec<ComplexBall> = (0..n)
        .map(|i| {
            let r = if bits == 0 {
                0
            } else {
                i.reverse_bits() >> (usize::BITS - bits)
            };
            x[r].clone()
        })
        .collect();
    let mut size = 2;
    while size <= n {
        let step = n / size;
        for start in (0..n).step_by(size) {
            for j in 0..size / 2 {
                let t = a[start + j + size / 2].mul(&w[j * step], prec);
                let u = a[start + j].clone();
                a[start + j] = u.add(&t, prec);
                a[start + j + size / 2] = u.sub(&t, prec);
            }
        }
        size *= 2;
    }
    a
}

/// The numeric transform of the list `l`, or its inverse.
fn numeric(ev: &mut Eval<'_>, e: &Expr, l: &Expr, inverse: bool) -> Native {
    let n = l.args().len();
    let guard = 2 * u64::from(usize::BITS - n.leading_zeros()) + 16;
    let prec = ev.context().working_precision() + guard;
    let Some(x) = l
        .args()
        .iter()
        .map(|v| ball(v, prec))
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    if n == 0 {
        return Ok(None);
    }
    let w = roots(n, inverse, prec);
    let mut found = if n.is_power_of_two() && n <= MAX_RADIX_TWO {
        radix_two(&x, &w, prec)
    } else if n <= MAX_DIRECT {
        direct(&x, &w, prec)
    } else {
        warn(
            ev,
            e,
            format!("a length of {n}, not a power of 2, above {MAX_DIRECT}"),
        );
        return Ok(None);
    };
    if inverse {
        let scale = ComplexBall::from_real(Ball::from_rational(
            &Rational::new(
                Integer::ONE,
                Integer::from(i64::try_from(n).expect("a length")),
            ),
            prec,
        ));
        found = found.iter().map(|z| z.mul(&scale, prec)).collect();
    }
    let result = Expr::call(
        Symbol::LIST,
        found.into_iter().map(number).collect::<Vec<_>>(),
    );
    ev.step(if inverse { "ifft" } else { "fft" }, e, &result);
    Ok(Some(result))
}

/// Whether `w` has order exactly `n` modulo `p`.
fn of_order(w: &Integer, n: usize, p: &Integer) -> bool {
    let n_i = Integer::from(i64::try_from(n).expect("a length"));
    if w.pow_mod(&n_i, p) != Some(Integer::ONE.div_rem_euclid(p).1) {
        return false;
    }
    // Each prime q of n: w^(n/q) is not 1.
    let mut m = n;
    let mut q = 2;
    while m > 1 {
        if m % q == 0 {
            let e = Integer::from(i64::try_from(n / q).expect("a length"));
            if w.pow_mod(&e, p) == Some(Integer::ONE) {
                return false;
            }
            while m % q == 0 {
                m /= q;
            }
        }
        q += 1;
    }
    true
}

/// The transform of `l` modulo `p` with the root `w`, or its inverse.
fn modular(ev: &mut Eval<'_>, e: &Expr, l: &Expr, w: &Expr, p: &Expr, inverse: bool) -> Native {
    let (Some(x), Some(w), Some(p)) = (
        l.args()
            .iter()
            .map(Expr::to_integer)
            .collect::<Option<Vec<_>>>(),
        w.to_integer(),
        p.to_integer().filter(|p| *p >= Integer::from(2)),
    ) else {
        return Ok(None);
    };
    let n = x.len();
    if n == 0 || n > MAX_MODULAR {
        return Ok(None);
    }
    if !of_order(&w, n, &p) {
        warn(ev, e, format!("{w} is not of order {n} modulo {p}"));
        return Ok(None);
    }
    let n_i = Integer::from(i64::try_from(n).expect("a length"));
    let (w, scale) = if inverse {
        let (Some(w), Some(scale)) = (w.mod_inverse(&p), n_i.mod_inverse(&p)) else {
            warn(ev, e, format!("{n} is not invertible modulo {p}"));
            return Ok(None);
        };
        (w, scale)
    } else {
        (w, Integer::ONE)
    };
    // The powers w^m, m < n.
    let mut powers = Vec::with_capacity(n);
    let mut current = Integer::ONE;
    for _ in 0..n {
        powers.push(current.clone());
        current = (&current * &w).div_rem_euclid(&p).1;
    }
    let found = (0..n).map(|k| {
        let sum = x.iter().enumerate().fold(Integer::ZERO, |s, (j, xj)| {
            (&s + &(xj * &powers[(j * k) % n])).div_rem_euclid(&p).1
        });
        Expr::from((&sum * &scale).div_rem_euclid(&p).1)
    });
    let result = Expr::call(Symbol::LIST, found.collect::<Vec<_>>());
    ev.step(if inverse { "ifft" } else { "fft" }, e, &result);
    Ok(Some(result))
}

fn transform(ev: &mut Eval<'_>, e: &Expr, inverse: bool) -> Native {
    match e.args() {
        [l] if l.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) => {
            numeric(ev, e, l, inverse)
        }
        [l, w, p] if l.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) => {
            modular(ev, e, l, w, p, inverse)
        }
        _ => Ok(None),
    }
}

/// `DiscreteFourierTransform[l]` and `DiscreteFourierTransform[l, w, p]`.
pub(super) fn discrete_fourier_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    transform(ev, e, false)
}

/// `InverseDiscreteFourierTransform[l]` and its modular form.
pub(super) fn inverse_discrete_fourier_transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    transform(ev, e, true)
}

/// `Convolution[a, b]`.
pub(super) fn convolution(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let list = Some(Symbol::LIST);
    if a.head().and_then(Expr::as_symbol) != list || b.head().and_then(Expr::as_symbol) != list {
        return Ok(None);
    }
    let (a, b) = (a.args(), b.args());
    if a.is_empty() || b.is_empty() {
        return Ok(None);
    }
    if a.len().saturating_mul(b.len()) > ev.context().limits().terms {
        return Ok(None);
    }
    let entries = (0..a.len() + b.len() - 1).map(|k| {
        let terms = (0..a.len())
            .filter(|&i| k >= i && k - i < b.len())
            .map(|i| a[i].clone() * b[k - i].clone());
        Expr::call(Symbol::PLUS, terms.collect::<Vec<_>>())
    });
    let result = Expr::call(Symbol::LIST, entries.collect::<Vec<_>>());
    ev.step("convolution", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
