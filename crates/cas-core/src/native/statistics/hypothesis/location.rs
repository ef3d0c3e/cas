//! The Z and T tests on a mean (CAS-DOC-52 §2, CAS-503): `ZTest`
//! (`normalt`) and `TTest` (`studentt`).
//!
//! - `ZTest[d, μ₀, σ, alt, α]`, `TTest[d, μ₀, s, alt, α]`: the data `d` is
//!   `{x̄, n}`, a mean and a sample size as giac reads a list of two, or a
//!   list of at least three numbers. `TTest` of a data list may leave `s`
//!   out, the sample standard deviation then, as giac's "Estimated sigma
//!   from sample(s)".
//! - `z = (x̄ − μ₀) √n / σ`, its p-value by `NormalCDF`, or `StudentTCDF`
//!   with `n − 1` degrees of freedom; where giac's `{x̄, n}` form takes `n`,
//!   and its two-sided `studentt` the `1 − α` quantile.

use super::{Native, decide, is_string, options, positive, real, symmetric_p, written};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::statistics::{is_list, refuse};
use crate::symbol::Symbol;

/// Which law the statistic follows.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Law {
    Normal,
    Student,
}

/// `ZTest` or `TTest` of `e`.
fn location(ev: &mut Eval<'_>, e: &Expr, law: Law) -> Native {
    let [d, mu, rest @ ..] = e.args() else {
        return Ok(None);
    };
    if !is_list(d) || !real(ev, mu) {
        return Ok(None);
    }
    let values = d.args();
    if !values.iter().all(|x| real(ev, x)) {
        return Ok(None);
    }
    // The mean and size, and the data when given.
    let (mean, n, data) = match values {
        [mean, n] => match n.as_i64() {
            Some(k) if k >= 1 => (mean.clone(), k, None),
            _ => return refuse(ev, e, "the sample size is a positive integer"),
        },
        _ if values.len() >= 3 => {
            let k = i64::try_from(values.len()).unwrap_or(i64::MAX);
            let mean = ev.eval(&Expr::call(Symbol::MEAN, [d.clone()]))?;
            (mean, k, Some(d))
        }
        _ => return refuse(ev, e, "the data is {mean, size} or at least three values"),
    };
    // σ, given or estimated.
    let (sigma, rest) = match (rest, data) {
        ([s, rest @ ..], _) if !is_string(s) => (s.clone(), rest),
        (rest, Some(d)) if law == Law::Student => (
            ev.eval(&Expr::call(Symbol::STANDARD_DEVIATION, [d.clone()]))?,
            rest,
        ),
        _ => return refuse(ev, e, "the standard deviation is required"),
    };
    if !positive(ev, &sigma) {
        return refuse(ev, e, "the standard deviation is a positive number");
    }
    if law == Law::Student && n < 2 {
        return refuse(ev, e, "a T test takes at least two values");
    }
    let Some((alternative, alpha)) = options(ev, e, rest)? else {
        return Ok(None);
    };
    let size = Expr::from(n);
    let statistic =
        ev.eval(&((mean - mu.clone()) * Expr::call(Symbol::SQRT, [size.clone()]) / sigma))?;
    let df = Expr::from(n - 1);
    let cdf = |s: Expr| match law {
        Law::Normal => Expr::call(Symbol::NORMAL_CDF, [s]),
        Law::Student => Expr::call(Symbol::STUDENT_T_CDF, [df.clone(), s]),
    };
    let p = symmetric_p(alternative, &statistic, &cdf);
    let (what, name) = match law {
        Law::Normal => ("Z test", "z"),
        Law::Student => ("T test", "t"),
    };
    let what = format!("{what}, H1: μ {} μ₀", written(alternative));
    let mut shown = vec![(name, statistic)];
    if law == Law::Student {
        shown.push(("df", df.clone()));
    }
    decide(ev, e, &what, &shown, &p, &alpha)
}

/// `ZTest[d, μ₀, σ, alt, α]`, giac's `normalt`.
pub(in crate::native) fn z_test(ev: &mut Eval<'_>, e: &Expr) -> Native {
    location(ev, e, Law::Normal)
}

/// `TTest[d, μ₀, s, alt, α]`, giac's `studentt`.
pub(in crate::native) fn t_test(ev: &mut Eval<'_>, e: &Expr) -> Native {
    location(ev, e, Law::Student)
}
