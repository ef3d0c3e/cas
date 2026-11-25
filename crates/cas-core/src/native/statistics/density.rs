//! Kernel density estimation (CAS-629): Xcas's `kernel_density` (`kde`)
//! (giac optimization.cc `select_bandwidth_dpi` :7168,
//! `select_bandwidth_dpi_bins` :7197, `kernel_density` :7218,
//! `_kernel_density` :7298; misc.cc quartiles :3134-3171).
//!
//! - **The data.** `KernelDensity[L, options…]`: a list of at least two
//!   real numbers (rationals, decimals, numeric constants), the Gaussian
//!   kernel estimate of their density. giac reads every number as a double;
//!   we keep rationals exact and take a decimal's midpoint (C0064).
//! - **The forms.** `exact` (or `output=exact`): `Σ exp(−(x−xᵢ)²/(2h²)) /
//!   (h n √(2π))` in the variable (`x`, or a bare symbol given as an
//!   option), exact on exact data with an exact `h`. `output=list`, the
//!   default: the estimate on `bins` points (100) of the grid `a … b`, giac's
//!   nearest-bin counts (`trunc((xᵢ − a)/d + 1/2)`, exact on exact data)
//!   convolved with the kernel's `2L + 1` weights, `L = min(bins − 1, ⌊1 +
//!   4h/d⌋)`. `piecewise` (or `output=piecewise`): giac's linear
//!   interpolation of that list, `Piecewise[x < g₀, 0, x < g₁, s₀(x − g₀) +
//!   y₀, …, x <= g_last, …, 0]`, its breakpoints exact on an exact grid. The
//!   last form given wins. giac convolves by FFT and so carries ~10⁻⁹
//!   relative noise (C0084); we sum directly, in certified balls at the
//!   working precision.
//! - **The bandwidth.** `bandwidth=h`, a positive real number; `bandwidth=
//!   select` or none: the direct plug-in rule, over the pairs of data for
//!   the exact form and for at most 1000 data (`select_bandwidth_dpi`), over
//!   the bin counts past them (`select_bandwidth_dpi_bins`, which sums
//!   ordered pairs; ported as it is); `bandwidth=gauss|normal|normald`:
//!   the rule of thumb `1.06 min(σ, IQR/1.34) n^(−1/5)` with giac's
//!   quartiles (the lower value, `v[⌈m/4⌉]`, `v[⌈3m/4⌉]`), which overrides
//!   any `h`. A selected bandwidth is a ball.
//! - **The grid.** `[min, max]` of the data when the bandwidth is selected
//!   by the plug-in rule, `[min − 3h, max + 3h]` otherwise; `range=a..b`,
//!   `range=[a, b]`, a bare `a..b` or `t=a..b` (which also names the
//!   variable) set it. `bins=n` sets its points.
//! - **Errors** leave the call with giac's message as the warning: a first
//!   argument no list, an option giac rejects (its text), fewer than two
//!   data, a datum no real number, `bins=1` (giac's FFT's "Invalid
//!   dimension"). giac's `interp=k` (`k > 1`), `spline=k` and `eval=v` are
//!   not ported: the call stays without a warning (`interp=1` is the
//!   default, accepted).
//! - **Limits** (`Limits::terms`), before any work: the exact form at most
//!   `terms` data, the bins times the data, the distinct pairs of the
//!   plug-in rule, the bins squared for its binned form, the bins times the
//!   kernel's `2L + 1` weights; past them the call stays with "too many
//!   points".

#![allow(
    clippy::many_single_char_names,
    reason = "a, b, d, h, n, x: giac's names for the grid, the bandwidth and the data"
)]

use std::cmp::Ordering;

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::numeric::approximation;
use crate::native::statistics::random_variables::{equation, range, sign};
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// giac's texts.
const LIST: &str = "Expected a list or a sequence of arguments";
const REAL: &str = "Expected a real constant";
const SMALL: &str = "Data too small";
const BANDWIDTH: &str = "Invalid bandwidth specification";
const BINS: &str = "Number of bins must be positive";
const RANGE: &str = "Invalid range specification";
const METHOD: &str = "Invalid method specification";
const INTERPOLATION: &str = "Invalid interpolation method specification";
const UNRECOGNIZED: &str = "Unrecognized option";
/// giac's FFT of one point (`bins=1`).
const DIMENSION: &str = "Invalid dimension";
/// Past `Limits::terms`.
const TOO_MANY: &str = "too many points";
/// giac's cut between the two plug-in rules.
const PAIRS_UP_TO: usize = 1000;

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// Whether `e` is the symbol `s` or a symbol named `name` (C0008).
fn is(e: &Expr, s: Symbol, name: &str) -> bool {
    e.as_symbol().is_some_and(|t| t == s || t.name() == name)
}

/// A real number: exact when `exact` is, and its ball at the working
/// precision.
#[derive(Clone, Debug)]
struct Real {
    exact: Option<Rational>,
    ball: Ball,
}

impl Real {
    fn rational(q: Rational, wp: u64) -> Real {
        Real {
            ball: Ball::from_rational(&q, wp),
            exact: Some(q),
        }
    }

    fn ball(ball: Ball) -> Real {
        Real { exact: None, ball }
    }

    /// `self + k other`.
    fn plus(&self, k: &Rational, other: &Real, wp: u64) -> Real {
        match (&self.exact, &other.exact) {
            (Some(a), Some(b)) => Real::rational(a + &(k * b), wp),
            _ => Real::ball(
                self.ball
                    .add(&Ball::from_rational(k, wp).mul(&other.ball, wp), wp),
            ),
        }
    }

    /// `k self`.
    fn scaled(&self, k: &Rational, wp: u64) -> Real {
        match &self.exact {
            Some(a) => Real::rational(a * k, wp),
            None => Real::ball(Ball::from_rational(k, wp).mul(&self.ball, wp)),
        }
    }

    /// The midpoint in `f64` (C0028).
    fn value(&self) -> f64 {
        self.ball.mid().to_f64()
    }

    fn expr(&self) -> Expr {
        match &self.exact {
            Some(q) => Expr::from(q),
            None => Expr::from(&self.ball),
        }
    }
}

fn ball(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

fn ratio(p: i64, q: i64) -> Rational {
    Rational::new(Integer::from(p), Integer::from(q))
}

/// A real number as giac's `is_real_number` reads it: a rational exactly,
/// a numeric expression by its approximation's midpoint; `None` otherwise.
fn real(ev: &Eval<'_>, x: &Expr) -> Option<Real> {
    let wp = ev.context().working_precision();
    if let Some(q) = x.to_rational() {
        return Some(Real::rational(q, wp));
    }
    if !x.is_numeric() {
        return None;
    }
    match approximation(ev.context(), x, wp).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(Real::ball(Ball::exact(b.mid().clone()))),
        _ => None,
    }
}

/// giac's `parse_interval`: two real ends, the right one larger.
fn interval(ev: &Eval<'_>, l: &Expr, r: &Expr) -> Option<(Real, Real)> {
    let (l, r) = (real(ev, l)?, real(ev, r)?);
    let larger = match (&l.exact, &r.exact) {
        (Some(p), Some(q)) => q > p,
        _ => r.value() > l.value(),
    };
    larger.then_some((l, r))
}

/// The output form.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Form {
    Exact,
    List,
    Piecewise,
}

/// The options, as giac's loop reads them (:7305-7363).
struct Options {
    variable: Expr,
    /// The given bandwidth, as given and as a number.
    bandwidth: Option<(Expr, Real)>,
    thumb: bool,
    bins: usize,
    ends: Option<(Real, Real)>,
    form: Form,
}

/// Why the options stop the call: giac's message, or `None` for an option
/// not ported (no warning).
type Stop = Option<&'static str>;

/// A positive integer option: its value (saturated), `None` for another.
fn positive(x: &Expr) -> Option<usize> {
    let n = x.to_integer()?;
    if n <= Integer::from(0) {
        return None;
    }
    Some(
        n.to_i64()
            .and_then(|v| usize::try_from(v).ok())
            .unwrap_or(usize::MAX),
    )
}

#[allow(
    clippy::too_many_lines,
    reason = "giac's options in its order, one arm each"
)]
fn options(ev: &Eval<'_>, given: &[Expr]) -> Result<Options, Stop> {
    let mut o = Options {
        variable: Expr::symbol("x"),
        bandwidth: None,
        thumb: false,
        bins: 100,
        ends: None,
        form: Form::List,
    };
    for x in given {
        if let Some((lhs, v)) = equation(x) {
            let name = lhs.as_symbol();
            if name == Some(Symbol::BANDWIDTH) {
                if is(v, Symbol::SELECT, "select") {
                    o.thumb = false;
                } else if is(v, Symbol::GAUSS_SQUARES, "gauss")
                    || is(v, Symbol::TOGETHER, "normal")
                    || is(v, Symbol::NORMAL_PDF, "normald")
                {
                    o.thumb = true;
                } else {
                    let h = real(ev, v).filter(|_| sign(ev, v) == Some(Ordering::Greater));
                    let Some(h) = h else {
                        return Err(Some(BANDWIDTH));
                    };
                    o.bandwidth = Some((v.clone(), h));
                }
            } else if name == Some(Symbol::BINS) {
                o.bins = positive(v).ok_or(Some(BINS))?;
            } else if name.is_some_and(|s| s.name() == "range") {
                let ends = if is_list(v) {
                    match v.args() {
                        [l, r] => interval(ev, l, r),
                        _ => None,
                    }
                } else {
                    range(v).and_then(|(l, r)| interval(ev, l, r))
                };
                o.ends = Some(ends.ok_or(Some(RANGE))?);
            } else if name.is_some_and(|s| matches!(s.name(), "output" | "Output")) {
                o.form = if is(v, Symbol::EXACT, "exact") {
                    Form::Exact
                } else if is(v, Symbol::PIECEWISE, "piecewise") {
                    Form::Piecewise
                } else if is(v, Symbol::LIST, "list") {
                    Form::List
                } else {
                    return Err(Some(METHOD));
                };
            } else if is(lhs, Symbol::INTERPOLATING_POLYNOMIAL, "interp")
                || is(lhs, Symbol::SPLINE, "spline")
            {
                let k = positive(v).ok_or(Some(INTERPOLATION))?;
                // giac's spline interpolation: not ported (interp=1 is ours).
                if k != 1 || is(lhs, Symbol::SPLINE, "spline") {
                    return Err(None);
                }
            } else if is(lhs, Symbol::EVALUATE, "eval") {
                return Err(None);
            } else if let Some(s) = name.filter(|s| !s.is_builtin()) {
                let ends = range(v).and_then(|(l, r)| interval(ev, l, r));
                o.ends = Some(ends.ok_or(Some(RANGE))?);
                o.variable = Expr::from(s);
            } else {
                return Err(Some(UNRECOGNIZED));
            }
        } else if is(x, Symbol::EXACT, "exact") {
            o.form = Form::Exact;
        } else if is(x, Symbol::PIECEWISE, "piecewise") {
            o.form = Form::Piecewise;
        } else if let Some(s) = x.as_symbol().filter(|s| !s.is_builtin()) {
            o.variable = Expr::from(s);
        } else if let Some((l, r)) = range(x) {
            o.ends = Some(interval(ev, l, r).ok_or(Some(RANGE))?);
        } else {
            return Err(Some(UNRECOGNIZED));
        }
    }
    Ok(o)
}

/// `KernelDensity[L, options…]`: Xcas's `kernel_density`.
pub(in crate::native) fn kernel_density(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((data, rest)) = e.args().split_first() else {
        return Ok(None);
    };
    if !is_list(data) {
        return refuse(ev, e, LIST);
    }
    let o = match options(ev, rest) {
        Ok(o) => o,
        Err(Some(why)) => return refuse(ev, e, why),
        Err(None) => return Ok(None),
    };
    let given = data.args();
    let n = given.len();
    if n < 2 {
        return refuse(ev, e, SMALL);
    }
    let Some(xs) = given
        .iter()
        .map(|x| real(ev, x))
        .collect::<Option<Vec<_>>>()
    else {
        return refuse(ev, e, REAL);
    };
    let terms = ev.context().limits().terms;
    let wp = ev.context().working_precision();
    let sorted = sorted(&xs);
    let sd = deviation(&xs, wp);
    let mut bandwidth = o.bandwidth.clone();
    if o.thumb {
        let h = thumb(&sorted, &sd, wp);
        // giac's rule of thumb at 0 leaves the plug-in rule to choose.
        bandwidth = (h.sign() == Some(Ordering::Greater)).then(|| (Expr::from(&h), Real::ball(h)));
    }
    if o.form == Form::Exact {
        if n > terms {
            return refuse(ev, e, TOO_MANY);
        }
        let h = if let Some((h, _)) = bandwidth {
            h
        } else {
            let groups = groups(&sorted);
            if !pairs_within(groups.len(), terms) {
                return refuse(ev, e, TOO_MANY);
            }
            let Some(h) = plug_in(&groups, n, &sd, wp) else {
                return Ok(None);
            };
            Expr::from(h)
        };
        return exact_form(ev, &o.variable, given, &h).map(Some);
    }
    if o.bins == 1 {
        return refuse(ev, e, DIMENSION);
    }
    if o.bins.checked_mul(n).is_none_or(|all| all > terms) {
        return refuse(ev, e, TOO_MANY);
    }
    let bandwidth = bandwidth.map(|(_, h)| h);
    let (a, b) = o.ends.clone().unwrap_or_else(|| {
        let (lo, hi) = (&sorted[0], &sorted[n - 1]);
        match &bandwidth {
            Some(h) => (lo.plus(&ratio(-3, 1), h, wp), hi.plus(&ratio(3, 1), h, wp)),
            None => (lo.clone(), hi.clone()),
        }
    });
    binned(ev, e, &o, &xs, &sorted, &sd, bandwidth, (&a, &b))
}

/// The list and piecewise forms on the grid `a … b` (:7233-7285).
#[allow(clippy::too_many_arguments, reason = "the call's parts so far")]
fn binned(
    ev: &mut Eval<'_>,
    e: &Expr,
    o: &Options,
    xs: &[Real],
    sorted: &[Real],
    sd: &Ball,
    bandwidth: Option<Real>,
    (a, b): (&Real, &Real),
) -> Native {
    let (n, bins) = (xs.len(), o.bins);
    let terms = ev.context().limits().terms;
    let wp = ev.context().working_precision();
    let span = b.plus(&ratio(-1, 1), a, wp);
    let d = span.scaled(&ratio(1, i64::try_from(bins - 1).unwrap_or(i64::MAX)), wp);
    if !(d.ball.is_finite() && d.value() > 0.0) {
        return Ok(None);
    }
    let counts = counts(xs, a, &d, bins, wp);
    let h = if let Some(h) = bandwidth {
        h
    } else {
        let h = if n <= PAIRS_UP_TO {
            let groups = groups(sorted);
            if !pairs_within(groups.len(), terms) {
                return refuse(ev, e, TOO_MANY);
            }
            plug_in(&groups, n, sd, wp)
        } else {
            if bins.checked_mul(bins).is_none_or(|all| all > terms) {
                return refuse(ev, e, TOO_MANY);
            }
            plug_in_bins(&counts, n, &d.ball, sd, wp)
        };
        let Some(h) = h else {
            return Ok(None);
        };
        Real::ball(h)
    };
    let reach = reach(&h, &d, bins);
    if bins
        .checked_mul(2 * reach + 1)
        .is_none_or(|all| all > terms)
    {
        return refuse(ev, e, TOO_MANY);
    }
    let values = convolved(&counts, n, &h.ball, &d.ball, reach, wp);
    if o.form == Form::List {
        return Ok(Some(Expr::call(
            Symbol::LIST,
            values.iter().map(Expr::from),
        )));
    }
    piecewise(ev, &o.variable, a, &d, &values).map(Some)
}

/// The data ascending: exactly when every datum is exact, by midpoint
/// otherwise (a total order either way).
fn sorted(xs: &[Real]) -> Vec<Real> {
    let mut v = xs.to_vec();
    if v.iter().all(|x| x.exact.is_some()) {
        v.sort_by(|p, q| p.exact.cmp(&q.exact));
    } else {
        v.sort_by(|p, q| p.value().total_cmp(&q.value()));
    }
    v
}

/// The sample standard deviation `√((Σx² − (Σx)²/n)/(n − 1))` (:7384),
/// its variance exact on exact data.
fn deviation(xs: &[Real], wp: u64) -> Ball {
    let n = xs.len() as u64;
    let exact: Option<Vec<&Rational>> = xs.iter().map(|x| x.exact.as_ref()).collect();
    if let Some(q) = exact {
        let s = q.iter().fold(Rational::from(0), |s, q| &s + *q);
        let s2 = q.iter().fold(Rational::from(0), |s, q| &s + &(*q * *q));
        let v = &(&s2 - &(&(&s * &s) / &Rational::from(n))) / &Rational::from(n - 1);
        return Ball::from_rational(&v, wp).sqrt(wp);
    }
    let mut s = Ball::ZERO;
    let mut s2 = Ball::ZERO;
    for x in xs {
        s.add_assign(&x.ball, wp);
        s2.add_assign(&x.ball.mul(&x.ball, wp), wp);
    }
    s2.sub(&s.mul(&s, wp).div(&ball(n), wp), wp)
        .div(&ball(n - 1), wp)
        .sqrt(wp)
}

/// giac's rule of thumb (:7385-7389): `1.06 min(σ, IQR/1.34) n^(−1/5)`,
/// the quartiles giac's lower values (misc.cc:3134-3171).
fn thumb(sorted: &[Real], sd: &Ball, wp: u64) -> Ball {
    let m = sorted.len();
    let q1 = &sorted[m.div_ceil(4) - 1];
    let q3 = &sorted[(3 * m).div_ceil(4) - 1];
    let iqr = q3
        .plus(&ratio(-1, 1), q1, wp)
        .ball
        .div(&Ball::from_rational(&ratio(134, 100), wp), wp);
    let least = if sd.mid().to_f64().total_cmp(&iqr.mid().to_f64()) == Ordering::Greater {
        iqr
    } else {
        sd.clone()
    };
    let fifth = ball(m as u64).pow_ball(&Ball::from_rational(&ratio(-1, 5), wp), wp);
    Ball::from_rational(&ratio(106, 100), wp)
        .mul(&least, wp)
        .mul(&fifth, wp)
}

/// The data's distinct values with their multiplicities: equal exact
/// values, or equal exact balls, adjacent in `sorted`.
fn groups(sorted: &[Real]) -> Vec<(Ball, u64)> {
    let mut out: Vec<(Real, u64)> = Vec::new();
    for x in sorted {
        if let Some((last, m)) = out.last_mut() {
            let same = match (&last.exact, &x.exact) {
                (Some(p), Some(q)) => p == q,
                (None, None) => {
                    last.ball.is_exact() && x.ball.is_exact() && last.ball.mid() == x.ball.mid()
                }
                _ => false,
            };
            if same {
                *m += 1;
                continue;
            }
        }
        out.push((x.clone(), 1));
    }
    out.into_iter().map(|(x, m)| (x.ball, m)).collect()
}

/// Whether the distinct pairs of `g` values are within `terms`.
fn pairs_within(g: usize, terms: usize) -> bool {
    g.checked_mul(g.saturating_sub(1))
        .is_some_and(|all| all / 2 <= terms)
}

/// The plug-in rule's two kernels of `t²`: the sixth and fourth
/// derivatives' sums (:7174, :7185).
fn sixth(t2: &Ball, wp: u64) -> Ball {
    // (2t²(t²(t² − 15) + 45) − 30) e^(−t²/2)
    let inner = t2.sub(&ball(15), wp).mul(t2, wp).add(&ball(45), wp);
    let poly = t2
        .mul(&inner, wp)
        .mul_2exp(&Integer::from(1))
        .sub(&ball(30), wp);
    poly.mul(&gauss(t2, wp), wp)
}

fn fourth(t2: &Ball, wp: u64) -> Ball {
    // (2t²(t² − 6) + 6) e^(−t²/2)
    let poly = t2
        .mul(&t2.sub(&ball(6), wp), wp)
        .mul_2exp(&Integer::from(1))
        .add(&ball(6), wp);
    poly.mul(&gauss(t2, wp), wp)
}

/// `e^(−t²/2)`.
fn gauss(t2: &Ball, wp: u64) -> Ball {
    t2.mul_2exp(&Integer::from(-1)).neg().exp(wp)
}

/// The final steps of both plug-in rules from the first stage's sum `s1`
/// and the second's at `g4` (:7179, :7189).
fn finish(g6: &Ball, s1: &Ball, n: u64, second: impl Fn(&Ball) -> Ball, wp: u64) -> Option<Ball> {
    let n = ball(n);
    let seventh = Ball::from_rational(&ratio(1, 7), wp);
    let g4 = g6.mul(
        &n.mul(&ball(6), wp).div(s1, wp).neg().pow_ball(&seventh, wp),
        wp,
    );
    let s2 = second(&g4);
    let root2 = ball(2).sqrt(wp);
    let fifth = Ball::from_rational(&ratio(1, 5), wp);
    let h = n
        .div(&root2.mul(&s2, wp), wp)
        .pow_ball(&fifth, wp)
        .mul(&g4, wp);
    (h.is_finite() && h.sign() == Some(Ordering::Greater)).then_some(h)
}

/// giac's `1.23044723 σ`.
fn first_scale(sd: &Ball, wp: u64) -> Ball {
    Ball::from_rational(&ratio(123_044_723, 100_000_000), wp).mul(sd, wp)
}

/// `select_bandwidth_dpi` (:7168-7190) over the distinct values: the sums
/// over pairs `i < j` and the diagonal constants `−15n`, `+3n`.
fn plug_in(groups: &[(Ball, u64)], n: usize, sd: &Ball, wp: u64) -> Option<Ball> {
    let n = n as u64;
    let pairs = |g: &Ball, f: &dyn Fn(&Ball, u64) -> Ball| {
        let mut s = Ball::ZERO;
        for (i, (u, mu)) in groups.iter().enumerate() {
            for (v, mv) in &groups[i + 1..] {
                let t = u.sub(v, wp).div(g, wp);
                s.add_assign(&f(&t.mul(&t, wp), wp).mul(&ball(mu * mv), wp), wp);
            }
            // The pairs within a value, t = 0.
            let within = mu * (mu - 1) / 2;
            if within > 0 {
                s.add_assign(&f(&Ball::ZERO, wp).mul(&ball(within), wp), wp);
            }
        }
        s
    };
    let g6 = first_scale(sd, wp);
    let s1 = pairs(&g6, &sixth).sub(&ball(15 * n), wp);
    finish(
        &g6,
        &s1,
        n,
        |g4| pairs(g4, &fourth).add(&ball(3 * n), wp),
        wp,
    )
}

/// `select_bandwidth_dpi_bins` (:7197-7215): the ordered sums over the bin
/// counts, `Σᵢ Σⱼ cᵢ cⱼ k(d(i − j)/g)`, through their autocorrelation.
fn plug_in_bins(counts: &[u64], n: usize, d: &Ball, sd: &Ball, wp: u64) -> Option<Ball> {
    let m = counts.len();
    let auto: Vec<u128> = (0..m)
        .map(|k| {
            counts[..m - k]
                .iter()
                .zip(&counts[k..])
                .map(|(&p, &q)| u128::from(p) * u128::from(q))
                .sum()
        })
        .collect();
    let sum = |g: &Ball, f: &dyn Fn(&Ball, u64) -> Ball| {
        let mut s = Ball::ZERO;
        for (k, &a) in auto.iter().enumerate() {
            if a == 0 {
                continue;
            }
            let t = d.mul(&ball(k as u64), wp).div(g, wp);
            let times = Ball::from_integer(&Integer::from(if k == 0 { a } else { 2 * a }));
            s.add_assign(&f(&t.mul(&t, wp), wp).mul(&times, wp), wp);
        }
        s
    };
    let g6 = first_scale(sd, wp);
    let s1 = sum(&g6, &sixth);
    finish(&g6, &s1, n as u64, |g4| sum(g4, &fourth), wp)
}

/// giac's nearest-bin counts (:7236-7240): `trunc((x − a)/d + 1/2)` when in
/// `0 … bins − 1`, exactly on exact data and grid.
fn counts(xs: &[Real], a: &Real, d: &Real, bins: usize, wp: u64) -> Vec<u64> {
    let mut c = vec![0u64; bins];
    let half = ratio(1, 2);
    for x in xs {
        let index = if let (Some(x), Some(a), Some(d)) = (&x.exact, &a.exact, &d.exact) {
            let t = (&(&(x - a) / d) + &half).trunc();
            t.to_i64().and_then(|t| usize::try_from(t).ok())
        } else {
            let u = x.ball.sub(&a.ball, wp).div(&d.ball, wp).mid().to_f64() + 0.5;
            index(u.trunc(), bins)
        };
        if let Some(i) = index.filter(|&i| i < bins) {
            c[i] += 1;
        }
    }
    c
}

/// A whole number `t` as an index below `below`; `None` outside (or NaN).
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    reason = "t is whole and checked in range"
)]
fn index(t: f64, below: usize) -> Option<usize> {
    (t >= 0.0 && t < below as f64).then_some(t as usize)
}

/// The kernel's reach `L = min(bins − 1, ⌊1 + 4h/d⌋)` (:7247), exactly
/// when both are exact.
fn reach(h: &Real, d: &Real, bins: usize) -> usize {
    let most = bins - 1;
    let wide = match (&h.exact, &d.exact) {
        (Some(h), Some(d)) => (&Rational::from(1) + &(&(h * &Rational::from(4)) / d))
            .floor()
            .to_i64()
            .and_then(|l| usize::try_from(l).ok()),
        _ => index((1.0 + 4.0 * h.value() / d.value()).floor(), most),
    };
    wide.map_or(most, |l| l.min(most))
}

/// The middle `bins` entries of the counts convolved with the kernel's
/// weights `e^(−(dk/h)²/2)/(n h √(2π))`, `k = −L … L` (:7248-7252), summed
/// directly (C0084).
fn convolved(counts: &[u64], n: usize, h: &Ball, d: &Ball, reach: usize, wp: u64) -> Vec<Ball> {
    let root = Ball::const_pi(wp).mul_2exp(&Integer::from(1)).sqrt(wp);
    let scale = ball(n as u64).mul(h, wp).mul(&root, wp).inv(wp);
    let weights: Vec<Ball> = (0..=reach)
        .map(|k| {
            let t = d.mul(&ball(k as u64), wp).div(h, wp);
            gauss(&t.mul(&t, wp), wp).mul(&scale, wp)
        })
        .collect();
    let bins = counts.len();
    let mut out = vec![Ball::ZERO; bins];
    for (i, &c) in counts.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let c = ball(c);
        for (j, y) in out
            .iter_mut()
            .enumerate()
            .take((i + reach).min(bins - 1) + 1)
            .skip(i.saturating_sub(reach))
        {
            y.add_assign(&c.mul(&weights[i.abs_diff(j)], wp), wp);
        }
    }
    out
}

/// `Σ exp(−(x − xᵢ)²/(2h²)) / (h n √(2π))` (:7221-7229), evaluated.
fn exact_form(ev: &mut Eval<'_>, x: &Expr, data: &[Expr], h: &Expr) -> Result<Expr, EvalError> {
    let twice = Expr::from(2) * h.pow(2);
    let terms: Vec<Expr> = data
        .iter()
        .map(|xi| Expr::call(Symbol::EXP, [-((x - xi.clone()).pow(2)) / twice.clone()]))
        .collect();
    let n = i64::try_from(data.len()).unwrap_or(i64::MAX);
    let root = Expr::call(Symbol::SQRT, [Expr::from(2) * Expr::from(Symbol::PI)]);
    let sum = Expr::call(Symbol::PLUS, terms);
    ev.eval(&(sum / (h.clone() * Expr::from(n) * root)))
}

/// giac's linear interpolation of the list (:7264-7282): `Piecewise[x <
/// g₀, 0, x < g₁, s₀(x − g₀) + y₀, …, x <= g_last, …, 0]`.
fn piecewise(
    ev: &mut Eval<'_>,
    x: &Expr,
    a: &Real,
    d: &Real,
    ys: &[Ball],
) -> Result<Expr, EvalError> {
    let wp = ev.context().working_precision();
    let at = |i: usize| a.plus(&Rational::from(i), d, wp);
    let bins = ys.len();
    let mut args = vec![
        Expr::call(Symbol::LESS, [x.clone(), at(0).expr()]),
        Expr::from(0),
    ];
    for i in 1..bins {
        let relation = if i + 1 < bins {
            Symbol::LESS
        } else {
            Symbol::LESS_EQUAL
        };
        args.push(Expr::call(relation, [x.clone(), at(i).expr()]));
        let slope = ys[i].sub(&ys[i - 1], wp).div(&d.ball, wp);
        args.push(Expr::from(slope) * (x - at(i - 1).expr()) + Expr::from(&ys[i - 1]));
    }
    args.push(Expr::from(0));
    ev.eval(&Expr::call(Symbol::PIECEWISE, args))
}
