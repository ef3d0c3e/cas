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
//!
//! `KolmogorovSmirnovTest` (`kolmogorovt`, CAS-506): of a sample against a
//! continuous law written as its density without the point (giac's
//! `normald(0,1)`, which stays: `NormalPDF[0, 1]`), `D = maxᵢ max(i/n −
//! F(x₍ᵢ₎), F(x₍ᵢ₎) − (i−1)/n)` and `K = D √n`; of two samples, `D =
//! sup |F₁ − F₂|` and `K = D √(n₁ n₂/(n₁ + n₂))`. `D` is exact as far as
//! comparisons decide, the answer giac's list `{"D=", D, "K=", K,
//! "1-kolmogorovd(K)=", 1 − K(K)}`.

#![allow(
    clippy::many_single_char_names,
    reason = "the statistics' names: D, K, n, m, F and the samples x and y"
)]

use std::cmp::Ordering;

use cas_num::ball::{Ball, Mag};
use cas_num::{Integer, Rational};

use super::{Native, compare, has_ball, settled, sorted};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::number::{Num, Precision};
use crate::native::rounding::real_sign;
use crate::native::statistics::{is_list, refuse};
use crate::symbol::Symbol;

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

/// The cumulative function of the continuous law written as its density
/// without the point, by its head.
fn cdf_of(law: &Expr) -> Option<Symbol> {
    let head = law.head()?.as_symbol()?;
    Some(match head {
        Symbol::NORMAL_PDF => Symbol::NORMAL_CDF,
        Symbol::UNIFORM_PDF => Symbol::UNIFORM_CDF,
        Symbol::EXPONENTIAL_PDF => Symbol::EXPONENTIAL_CDF,
        Symbol::STUDENT_T_PDF => Symbol::STUDENT_T_CDF,
        Symbol::CHI_SQUARE_PDF => Symbol::CHI_SQUARE_CDF,
        Symbol::FISHER_PDF => Symbol::FISHER_CDF,
        Symbol::GAMMA_PDF => Symbol::GAMMA_CDF,
        Symbol::BETA_PDF => Symbol::BETA_CDF,
        Symbol::WEIBULL_PDF => Symbol::WEIBULL_CDF,
        Symbol::CAUCHY_PDF => Symbol::CAUCHY_CDF,
        _ => return None,
    })
}

/// The largest of `candidates`; `None` when comparisons do not decide.
fn largest(ev: &mut Eval<'_>, candidates: Vec<Expr>) -> Result<Option<Expr>, EvalError> {
    let mut best: Option<Expr> = None;
    for c in candidates {
        let c = ev.eval(&c)?;
        best = Some(match best {
            None => c,
            Some(b) => match compare(ev, &c, &b)? {
                Some(Ordering::Greater) => c,
                Some(_) => b,
                None => return Ok(None),
            },
        });
    }
    Ok(best)
}

/// `D` of two sorted samples: the largest gap of their empirical laws,
/// exact.
fn two_sample(ev: &mut Eval<'_>, a: &[Expr], b: &[Expr]) -> Result<Option<Expr>, EvalError> {
    let (n, m) = (a.len(), b.len());
    let share = |k: usize, of: usize| {
        Rational::new(
            Integer::from(u64::try_from(k).unwrap_or(u64::MAX)),
            Integer::from(u64::try_from(of).unwrap_or(u64::MAX)),
        )
    };
    let (mut i, mut j) = (0, 0);
    let mut d = Rational::ZERO;
    while i < n && j < m {
        let v = match compare(ev, &a[i], &b[j])? {
            Some(Ordering::Greater) => b[j].clone(),
            Some(_) => a[i].clone(),
            None => return Ok(None),
        };
        // Past every value equal to `v` in both.
        while i < n {
            match compare(ev, &a[i], &v)? {
                Some(Ordering::Equal) => i += 1,
                Some(_) => break,
                None => return Ok(None),
            }
        }
        while j < m {
            match compare(ev, &b[j], &v)? {
                Some(Ordering::Equal) => j += 1,
                Some(_) => break,
                None => return Ok(None),
            }
        }
        let gap = (&share(i, n) - &share(j, m)).abs();
        if gap > d {
            d = gap;
        }
    }
    // Once one sample is spent, the gap only closes.
    Ok(Some(Expr::from(d)))
}

/// `KolmogorovSmirnovTest[l, law]` and `KolmogorovSmirnovTest[l₁, l₂]`,
/// giac's `kolmogorovt`: `{"D=", D, "K=", K, "1-kolmogorovd(K)=", 1 −
/// K(K)}`.
pub(in crate::native) fn kolmogorov_smirnov_test(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, other] = e.args() else {
        return Ok(None);
    };
    let Some(x) = sorted(ev, l)? else {
        return Ok(None);
    };
    let n = i64::try_from(x.len()).unwrap_or(i64::MAX);
    let (d, size) = if let Some(cdf) = cdf_of(other) {
        let mut candidates = Vec::with_capacity(2 * x.len());
        for (i, xi) in x.iter().enumerate() {
            let mut args = other.args().to_vec();
            args.push(xi.clone());
            let f = ev.eval(&Expr::call(cdf, args))?;
            let i = i64::try_from(i).unwrap_or(i64::MAX);
            candidates.push(Expr::from(Rational::new((i + 1).into(), n.into())) - f.clone());
            candidates.push(f - Expr::from(Rational::new(i.into(), n.into())));
        }
        let Some(d) = largest(ev, candidates)? else {
            return refuse(ev, e, "the gaps to the law are not ordered");
        };
        (d, Expr::from(n))
    } else if is_list(other) {
        let Some(y) = sorted(ev, other)? else {
            return Ok(None);
        };
        let Some(d) = two_sample(ev, &x, &y)? else {
            return refuse(ev, e, "the samples' values are not ordered");
        };
        let m = i64::try_from(y.len()).unwrap_or(i64::MAX);
        (d, Expr::from(Rational::new((n * m).into(), (n + m).into())))
    } else {
        return Ok(None);
    };
    // Decimal data: balls, not exact expressions of them.
    let inexact = has_ball(e);
    let d = settled(ev, d, inexact)?;
    let k = ev.eval(&(d.clone() * Expr::call(Symbol::SQRT, [size])))?;
    let k = settled(ev, k, inexact)?;
    let p = ev.eval(&(Expr::from(1) - Expr::call(Symbol::KOLMOGOROV_CDF, [k.clone()])))?;
    let p = settled(ev, p, inexact)?;
    let result = Expr::call(
        Symbol::LIST,
        [
            Expr::string("D="),
            d,
            Expr::string("K="),
            k,
            Expr::string("1-kolmogorovd(K)="),
            p,
        ],
    );
    ev.step("kolmogorovt", e, &result);
    Ok(Some(result))
}
