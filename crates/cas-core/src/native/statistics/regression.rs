//! Regressions (CAS-DOC-22 §5, CAS-297): least squares, exact, of the
//! points `[x, y]` or of a matrix of two columns `[m]`.
//!
//! - `LinearRegression` (`linear_regression`): `{a, b}` for `y = a x + b`,
//!   `a = S_xy / S_xx` and `b = ȳ − a x̄`. Constant abscissas give
//!   `Indeterminate`, giac's `undef`.
//! - `PolynomialRegression[x, y, n]`, `PolynomialRegression[m, n]`
//!   (`polynomial_regression`): `{aₙ, …, a₀}`, from the normal equations
//!   `(VᵀV) c = Vᵀ y` of the Vandermonde matrix `V`. When they are
//!   singular, with fewer than `n + 1` distinct abscissas, the call stays;
//!   giac answers a parametric solution in `C_0`.
//! - `ExponentialRegression` (`exponential_regression`), `{a, b}` for
//!   `y = b aˣ`; `LogarithmicRegression` (`logarithmic_regression`), for
//!   `y = a ln x + b`; `PowerRegression` (`power_regression`), for
//!   `y = b xᵃ`. They are the linear regressions of `ln y` on `x`, of `y`
//!   on `ln x`, and of `ln y` on `ln x`. The logarithms of positive
//!   rationals are expanded over their factors (`ln 12 = 2 ln 2 + ln 3`),
//!   the slope is put together, and `e^(c ln p)` is `p^c`, so that exact
//!   data give exact answers: `exponential_regression([0,1,2],[3,6,12])` is
//!   `{2, 3}`, where giac answers an expression of logarithms. A value that
//!   is not positive under a logarithm leaves the call.
//! - `Logistic[x, x0, L, k]` (`logistic`, CAS-627), `L/(1 + e^(−k(x − x0)))`
//!   with giac's defaults `x0 = 0`, `L = k = 1`, evaluated as it stands:
//!   `logistic(0)` is `1/2`. No argument or more than four stay, with
//!   giac's Invalid dimension.
//! - `LogisticRegression[w, t0, P0]` (`logistic_regression`, CAS-627),
//!   giac's closed form (misc.cc `logistic_regression`): the productions
//!   `wᵢ` at times `t0, t0 + 1, …`, cumulated from `P0`, `Cᵢ`; the line
//!   `wᵢ/Cᵢ = a Cᵢ + b`, exact (balls for decimal data, as giac's), gives
//!   the total `U = −b/a`; the line of
//!   `ln(U/Cᵢ − 1)`, certified balls, on the times gives `b₂ x + c₂`. The
//!   answer is giac's `{U/(1 + e^(b₂x + c₂)), U b/2/(1 + cosh(b₂x + c₂)), U,
//!   U b/4, −c₂/b₂, r}` in the symbol `x`, `r` the correlation of `Cᵢ` and
//!   `wᵢ/Cᵢ`, exact as giac's. A slope `a` of 0 (constant quotients) gives
//!   giac's `{undef, undef, ∞, ∞, undef, undef}`, `ComplexInfinity` for ∞.
//!   Stay: a `t0` that is no real number (giac's Bad Argument Value), a
//!   negative `U` (giac's message), data that are not real numbers, a
//!   cumulated value that is not of known sign, and an `U/Cᵢ − 1` that is
//!   not positive, which giac answers with complex logarithms. Departure:
//!   the two-argument form, where giac guesses `P0` from 20 or more values,
//!   stays; so does giac's interval of times `t0..t1`.

use std::cmp::Ordering;

use cas_num::Integer;

use super::super::linear::dot;
use super::super::matrix::inverse_of;
use super::super::rounding::real_sign;
use super::pairs;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// Pollard rho's iterations on a factor under a logarithm. Exact data have
/// small factors; a composite it leaves unsplit is a factor of its own,
/// and the expansion is exact all the same.
const RHO_ITERATIONS: u64 = 1 << 12;

/// The points of `args`, `[x, y]` or `[m]`, without frequencies.
fn points(args: &[Expr]) -> Option<(Vec<Expr>, Vec<Expr>)> {
    let (x, y) = pairs(args)?;
    x.frequencies.is_none().then_some((x.values, y.values))
}

/// `Σ (uᵢ − ū)(vᵢ − v̄)`, expanded.
fn product_sum(u: &[Expr], u_mean: &Expr, v: &[Expr], v_mean: &Expr) -> Expr {
    let terms: Vec<Expr> = u
        .iter()
        .zip(v)
        .map(|(a, b)| (a - u_mean) * (b - v_mean))
        .collect();
    Expr::call(Symbol::EXPAND, [Expr::call(Symbol::PLUS, terms)])
}

/// The least-squares line `y = a x + b` through the points: `(a, b)`, `a`
/// put together and `b` expanded. A ball `a` stays a ball: `Together`
/// takes exact data only and would stay around it (CAS-627).
fn line(ev: &mut Eval<'_>, xs: &[Expr], ys: &[Expr]) -> Result<(Expr, Expr), EvalError> {
    let n = Expr::from(i64::try_from(xs.len()).unwrap_or(i64::MAX));
    let mean = |v: &[Expr]| Expr::call(Symbol::PLUS, v.to_vec()) / &n;
    let x_mean = ev.eval(&mean(xs))?;
    let y_mean = ev.eval(&mean(ys))?;
    let sxy = ev.eval(&product_sum(xs, &x_mean, ys, &y_mean))?;
    let sxx = ev.eval(&product_sum(xs, &x_mean, xs, &x_mean))?;
    let slope = ev.eval(&(sxy / sxx))?;
    let slope = if matches!(slope.view(), View::Ball(_)) {
        slope
    } else {
        ev.eval(&Expr::call(Symbol::TOGETHER, [slope]))?
    };
    let intercept = Expr::call(Symbol::EXPAND, [y_mean - &slope * x_mean]);
    Ok((slope, ev.eval(&intercept)?))
}

/// `ln v` for a positive `v`: over the factors of its numerator and
/// denominator when it is a rational; `None` when `v` is not positive.
fn log(ev: &mut Eval<'_>, v: &Expr) -> Result<Option<Expr>, EvalError> {
    if real_sign(ev.context(), v) != Some(Ordering::Greater) {
        return Ok(None);
    }
    let Some(q) = v.to_rational() else {
        return Ok(Some(ev.eval(&Expr::call(Symbol::LOG, [v.clone()]))?));
    };
    let mut terms = Vec::new();
    for (n, sign) in [(q.numer(), 1), (q.denom(), -1)] {
        let m = n.unsigned_abs();
        if m.is_one() {
            continue;
        }
        for f in m.factor_with_effort(RHO_ITERATIONS).factors() {
            let k = i64::try_from(f.exponent()).unwrap_or(i64::MAX) * sign;
            let base = Expr::from(Integer::from(f.base().clone()));
            terms.push(Expr::from(k) * Expr::call(Symbol::LOG, [base]));
        }
    }
    Ok(Some(ev.eval(&Expr::call(Symbol::PLUS, terms))?))
}

/// The logarithms of `values`; `None` when one is not positive.
fn logs(ev: &mut Eval<'_>, values: &[Expr]) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut out = Vec::with_capacity(values.len());
    for v in values {
        let Some(l) = log(ev, v)? else {
            return Ok(None);
        };
        out.push(l);
    }
    Ok(Some(out))
}

/// `c ln p`, `ln p` among them: `(c, p)`.
fn log_term(t: &Expr) -> Option<(Expr, Expr)> {
    let is_log = |l: &Expr| l.head() == Some(&Expr::from(Symbol::LOG)) && l.args().len() == 1;
    match t.args() {
        [p] if is_log(t) => Some((Expr::from(1), p.clone())),
        [c, l] if t.head() == Some(&Expr::from(Symbol::TIMES)) && c.is_number() && is_log(l) => {
            Some((c.clone(), l.args()[0].clone()))
        }
        _ => None,
    }
}

/// `e^x`, its exponent expanded and each term `c ln p` of it made `p^c`:
/// `e^(ln 2 − ln 3)` is `2/3`, which evaluation leaves as it is.
fn exp(ev: &mut Eval<'_>, x: &Expr) -> Result<Expr, EvalError> {
    let x = ev.eval(&Expr::call(Symbol::EXPAND, [x.clone()]))?;
    let terms = if x.head() == Some(&Expr::from(Symbol::PLUS)) {
        x.args().to_vec()
    } else {
        vec![x]
    };
    let mut factors = Vec::with_capacity(terms.len() + 1);
    let mut rest = Vec::new();
    for t in terms {
        match log_term(&t) {
            Some((c, p)) => factors.push(p.pow(c)),
            None => rest.push(t),
        }
    }
    factors.push(Expr::call(Symbol::EXP, [Expr::call(Symbol::PLUS, rest)]));
    ev.eval(&Expr::call(Symbol::TIMES, factors))
}

/// `{a, b}`, evaluated.
fn pair(ev: &mut Eval<'_>, a: Expr, b: Expr) -> Native {
    Ok(Some(ev.eval(&Expr::call(Symbol::LIST, [a, b]))?))
}

/// `LinearRegression[x, y]`, `LinearRegression[m]`.
pub(in super::super) fn linear_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys)) = points(e.args()) else {
        return Ok(None);
    };
    let (a, b) = line(ev, &xs, &ys)?;
    pair(ev, a, b)
}

/// `ExponentialRegression[x, y]`, `ExponentialRegression[m]`.
pub(in super::super) fn exponential_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys)) = points(e.args()) else {
        return Ok(None);
    };
    let Some(ys) = logs(ev, &ys)? else {
        return Ok(None);
    };
    let (a, b) = line(ev, &xs, &ys)?;
    let (a, b) = (exp(ev, &a)?, exp(ev, &b)?);
    pair(ev, a, b)
}

/// `LogarithmicRegression[x, y]`, `LogarithmicRegression[m]`.
pub(in super::super) fn logarithmic_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys)) = points(e.args()) else {
        return Ok(None);
    };
    let Some(xs) = logs(ev, &xs)? else {
        return Ok(None);
    };
    let (a, b) = line(ev, &xs, &ys)?;
    pair(ev, a, b)
}

/// `PowerRegression[x, y]`, `PowerRegression[m]`.
pub(in super::super) fn power_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((xs, ys)) = points(e.args()) else {
        return Ok(None);
    };
    let (Some(xs), Some(ys)) = (logs(ev, &xs)?, logs(ev, &ys)?) else {
        return Ok(None);
    };
    let (a, b) = line(ev, &xs, &ys)?;
    let b = exp(ev, &b)?;
    pair(ev, a, b)
}

/// `x^p`, with `x⁰ = 1` for every `x`, 0 included.
fn power(x: &Expr, p: usize) -> Expr {
    match p {
        0 => Expr::from(1),
        1 => x.clone(),
        _ => x.pow(i64::try_from(p).unwrap_or(i64::MAX)),
    }
}

/// `PolynomialRegression[x, y, n]`, `PolynomialRegression[m, n]`.
pub(in super::super) fn polynomial_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((degree, data)) = e.args().split_last() else {
        return Ok(None);
    };
    let Some((xs, ys)) = points(data) else {
        return Ok(None);
    };
    // Fewer points than coefficients: singular.
    let Some(n) = degree
        .to_integer()
        .and_then(|d| usize::try_from(d.to_i64()?).ok())
        .filter(|&n| n < xs.len())
    else {
        return Ok(None);
    };
    let mut power_sums = Vec::with_capacity(2 * n + 1);
    for p in 0..=2 * n {
        let terms: Vec<Expr> = xs.iter().map(|x| power(x, p)).collect();
        power_sums.push(ev.eval(&Expr::call(Symbol::PLUS, terms))?);
    }
    let normal = Expr::call(
        Symbol::LIST,
        (0..=n)
            .map(|j| Expr::call(Symbol::LIST, power_sums[j..=j + n].to_vec()))
            .collect::<Vec<_>>(),
    );
    let mut moments = Vec::with_capacity(n + 1);
    for j in 0..=n {
        let terms: Vec<Expr> = xs.iter().zip(&ys).map(|(x, y)| y * power(x, j)).collect();
        moments.push(ev.eval(&Expr::call(Symbol::PLUS, terms))?);
    }
    let Some(inverse) = inverse_of(ev, e, &normal)? else {
        return Ok(None);
    };
    // The matrix product: `Dot` is Xcas's `dot`, of two vectors.
    let product = dot(&inverse, &Expr::call(Symbol::LIST, moments))?;
    let solution = ev.eval(&product)?;
    let coefficients: Vec<Expr> = solution
        .args()
        .iter()
        .rev()
        .map(|c| Expr::call(Symbol::EXPAND, [c.clone()]))
        .collect();
    Ok(Some(ev.eval(&Expr::call(Symbol::LIST, coefficients))?))
}

/// The call stays, with giac's message as the warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// Whether `v` is a real number atom, exact or a ball.
fn real(v: &Expr) -> bool {
    matches!(
        v.view(),
        View::Int(_) | View::Integer(_) | View::Rational(_) | View::Ball(_)
    )
}

/// `Logistic[x, x0, L, k]`, the last three optional.
pub(in super::super) fn logistic(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    let [x, rest @ ..] = args else {
        return refuse(ev, e, "Invalid dimension");
    };
    if rest.len() > 3 {
        return refuse(ev, e, "Invalid dimension");
    }
    let x0 = rest.first().cloned().unwrap_or_else(|| Expr::from(0));
    let l = rest.get(1).cloned().unwrap_or_else(|| Expr::from(1));
    let k = rest.get(2).cloned().unwrap_or_else(|| Expr::from(1));
    let decay = Expr::call(Symbol::EXP, [-(k * (x - x0))]);
    Ok(Some(ev.eval(&(l / (1 + decay)))?))
}

/// `LogisticRegression[w, t0, P0]`.
pub(in super::super) fn logistic_regression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [data, t0, p0] = e.args() else {
        return Ok(None);
    };
    if data.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(None);
    }
    let values = data.args();
    if values.len() > ev.context().limits().terms {
        return refuse(ev, e, "too many values");
    }
    if !real(t0) {
        return refuse(ev, e, "Bad Argument Value");
    }
    if values.len() < 2 || !real(p0) || !values.iter().all(real) {
        return Ok(None);
    }
    // The cumulated productions and the quotients.
    let mut cumulated = Vec::with_capacity(values.len());
    let mut quotients = Vec::with_capacity(values.len());
    let mut total = p0.clone();
    for value in values {
        total = ev.eval(&(total + value))?;
        if !matches!(
            real_sign(ev.context(), &total),
            Some(Ordering::Less | Ordering::Greater)
        ) {
            return Ok(None);
        }
        quotients.push(ev.eval(&(value / &total))?);
        cumulated.push(total.clone());
    }
    let list = |v: &[Expr]| Expr::call(Symbol::LIST, v.to_vec());
    let correlation = ev.eval(&Expr::call(
        Symbol::CORRELATION,
        [list(&cumulated), list(&quotients)],
    ))?;
    let (slope, intercept) = line(ev, &cumulated, &quotients)?;
    let undef = Expr::from(Symbol::INDETERMINATE);
    let ctx = ev.context();
    if !matches!(
        real_sign(ctx, &slope),
        Some(Ordering::Less | Ordering::Greater)
    ) {
        // giac's −b/0: `infinity` (unsigned), or `undef` for 0/0.
        let total = match real_sign(ctx, &intercept) {
            Some(Ordering::Less | Ordering::Greater) => Expr::from(Symbol::COMPLEX_INFINITY),
            _ => undef.clone(),
        };
        let answer = [
            undef.clone(),
            undef.clone(),
            total.clone(),
            total,
            undef,
            correlation,
        ];
        return Ok(Some(ev.eval(&list(&answer))?));
    }
    let urr = ev.eval(&(-(&intercept / &slope)))?;
    if real_sign(ev.context(), &urr) == Some(Ordering::Less) {
        return refuse(
            ev,
            e,
            "Data not suited for logistic regression (instant production should be a \
             decreasing function of cumulated production)",
        );
    }
    // The line of ln(U/Cᵢ − 1) on the times.
    let mut times = Vec::with_capacity(values.len());
    let mut logs = Vec::with_capacity(values.len());
    for (i, cumul) in cumulated.iter().enumerate() {
        let ratio = ev.eval(&(&urr / cumul - 1))?;
        if real_sign(ev.context(), &ratio) != Some(Ordering::Greater) {
            return Ok(None);
        }
        let ln = Expr::call(Symbol::LOG, [ratio]);
        logs.push(ev.eval(&Expr::call(Symbol::N, [ln]))?);
        let i = i64::try_from(i).unwrap_or(i64::MAX);
        times.push(ev.eval(&(t0 + i))?);
    }
    let (b2, bt0) = line(ev, &times, &logs)?;
    let exponent = &b2 * Expr::symbol("x") + &bt0;
    let answer = [
        &urr / (1 + Expr::call(Symbol::EXP, [exponent.clone()])),
        &urr * &intercept / 2 / (1 + Expr::call(Symbol::COSH, [exponent])),
        urr.clone(),
        &urr * &intercept / 4,
        -(bt0 / b2),
        correlation,
    ];
    Ok(Some(ev.eval(&list(&answer))?))
}
