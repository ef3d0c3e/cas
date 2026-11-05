//! Kolmogorov's distribution (CAS-DOC-52 §4, CAS-505): `KolmogorovCDF`
//! (`kolmogorovd`), `K(x) = 1 − 2 Σ_{k≥1} (−1)^(k−1) e^(−2k²x²)`, the law of
//! `√n D` for large `n`.
//!
//! - **Balls**, `x > 0`: from `x ≥ 1` by the series above, alternating with
//!   decreasing terms, its tail below the first term left out; below 1 by
//!   Jacobi's form `K(x) = (√(2π)/x) Σ_{k≥1} e^(−(2k−1)² π²/(8x²))`, whose
//!   ratio of terms is below `e^(−π²)`, the tail below twice the first term
//!   left out.
//! - `K(x) = 0` for `x ≤ 0` (giac: `undef` at 0); an exact positive `x`
//!   stays, for `N`; anything else stays.

use std::cmp::Ordering;

use cas_num::Integer;
use cas_num::ball::{Ball, Mag};

use super::Native;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::number::{Num, Precision};
use crate::native::rounding::real_sign;

/// The most terms summed: far more than any precision asks.
const MAX_TERMS: u64 = 100_000;

/// An upper bound of `|x|` over the ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// `K(x)` for a ball `x > 0` at `prec` bits; `None` when a series does not
/// settle within [`MAX_TERMS`].
pub(in crate::native) fn kolmogorov_cdf(x: &Ball, prec: u64) -> Option<Ball> {
    let w = prec + 32;
    let small = Mag::power_of_two(&-Integer::from(w + 8));
    let one = Ball::from_integer(&Integer::from(1));
    let large = x.mid().to_f64() >= 1.0;
    let x2 = x.mul(x, w);
    if large {
        // 1 − 2 Σ (−1)^(k−1) e^(−2k²x²).
        let mut sum = Ball::from_integer(&Integer::from(0));
        for k in 1..=MAX_TERMS {
            let k2 = Ball::from_integer(&Integer::from(2 * k * k));
            let t = k2.mul(&x2, w).neg().exp(w);
            if upper(&t) < small {
                let tail = upper(&t);
                let sum = sum.add_error(&tail);
                return Some(one.sub(&sum.mul_2exp(&Integer::from(1)), w));
            }
            sum = if k % 2 == 1 {
                sum.add(&t, w)
            } else {
                sum.sub(&t, w)
            };
        }
        return None;
    }
    // (√(2π)/x) Σ e^(−(2k−1)² π²/(8x²)).
    let pi = Ball::const_pi(w);
    let scale = pi.mul(&pi, w).div(&x2.mul_2exp(&Integer::from(3)), w);
    let mut sum = Ball::from_integer(&Integer::from(0));
    let mut first: Option<Mag> = None;
    for k in 1..=MAX_TERMS {
        let odd = Ball::from_integer(&Integer::from((2 * k - 1) * (2 * k - 1)));
        let t = odd.mul(&scale, w).neg().exp(w);
        let bound = upper(&t);
        let lead = first.get_or_insert_with(|| bound.clone()).clone();
        if k > 1 && bound < lead.mul(&small) {
            let tail = bound.add(&bound);
            let sum = sum.add_error(&tail);
            let factor = pi.mul_2exp(&Integer::from(1)).sqrt(w).div(x, w);
            return Some(factor.mul(&sum, w));
        }
        sum = sum.add(&t, w);
    }
    None
}

/// `KolmogorovCDF[x]`.
pub(in crate::native) fn kolmogorov_cdf_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if let Some(Num::Ball(b)) = Num::of(x) {
        let bits = Num::Ball(b.clone()).precision(Precision::of(ev.context()));
        return Ok(match b.sign() {
            Some(Ordering::Greater) => kolmogorov_cdf(&b, bits)
                .filter(Ball::is_finite)
                .map(|k| Num::Ball(Box::new(k.round_to_accuracy(bits))).to_expr()),
            Some(_) => Some(Expr::from(0)),
            None => None,
        });
    }
    // Exact: 0 from 0 down; above, the call stays for `N`.
    if x.is_number() && real_sign(ev.context(), x).is_some_and(|s| s != Ordering::Greater) {
        return Ok(Some(Expr::from(0)));
    }
    Ok(None)
}
