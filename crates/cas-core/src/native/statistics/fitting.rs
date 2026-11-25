//! Distribution fitting (CAS-629): Xcas's `fitdistr` (giac
//! optimization.cc `weibull_mle` :7415, `cauchy_mle` :7433, `_fitdistr`
//! :7460; misc.cc quartiles :3108-3171).
//!
//! - **The call.** `FitDistribution[L, law]`: at least two real numbers
//!   (rationals, decimals, numeric constants) and a law's name as giac
//!   compares it (`normal`, `normald`, `NORMALD`; `poisson`, `POISSON`;
//!   `exp`, `EXP`, `exponential`, `exponentiald`; `geometric`; `gammad`,
//!   `Gamma`; `betad`, `Beta`; `cauchy`, `cauchyd`; `weibull`, `weibulld`).
//!   The answer is the law's object (the random variables' module doc), the
//!   maximum likelihood estimate as giac finds it, except that the mean is
//!   `Σx/N` where giac divides by `N − 1` (:7481); the variance is `Σ(x −
//!   m)²/(N − 1)`. Rationals stay exact; a decimal's midpoint is taken
//!   (C0064).
//! - **Closed forms**, exact on exact data, balls otherwise: `NormalPDF[m,
//!   √v]`, `PoissonPDF[m]`, `ExponentialPDF[1/m]`, `GeometricPDF[1/m]`.
//! - **Iterative fits** in certified balls from giac's starting points:
//!   gamma `ln a − ψ(a) = ln m − Σln x/N` from `m²/v` (a step past 0 halves
//!   `a`), `GammaPDF[a, a/m]`; beta `ψ(a) − ψ(a + b) = Σln x/N` and `ψ(b) −
//!   ψ(a + b) = Σln(1 − x)/N` from giac's moment start; Cauchy giac's Newton
//!   system on the score equations from the median and half the quartiles'
//!   gap (giac's lower values `v[⌈N/4⌉]`, `v[⌈N/2⌉]`, `v[⌈3N/4⌉]`), the
//!   scale's sign dropped (the equations are even in it); Weibull giac's
//!   `k_init` (Newton on `Γ(1 + 2/k)/Γ(1 + 1/k) − v/m² − 1`, unsquared as
//!   giac has it, from `max(1, 1/v)`; that start itself when the solve
//!   fails), then Newton on the profile score `Σxᵏ ln x − Σxᵏ (1/k + Σln
//!   x/N)`, and `λ = (Σxᵏ/N)^(1/k)` at the root. Newton runs in balls at the
//!   working precision plus 32 bits, each iterate an exact midpoint (C0064),
//!   the inverse Jacobian in `f64`, until a step falls below the working
//!   precision; a Krawczyk step on a box about the iterate then certifies
//!   the root (C0061), whose ball is the answer (giac stops at a step of
//!   10⁻⁵).
//! - **Errors** leave the call with giac's message as the warning, in
//!   giac's order: no second argument ("Bad Argument Type"), more than two
//!   ("Expected a sequence of two arguments"), data no list ("Expected a
//!   list"), fewer than two ("At least two samples are required"), a datum
//!   no real number ("Expected a real constant"); then the law's checks:
//!   Poisson data not integers `≥ 0` ("Expected a nonnegative integer"),
//!   geometric data not integers `≥ 1` ("Expected a positive integer"),
//!   exponential, gamma or Weibull data `≤ 0` ("Expected a positive
//!   number", giac's `is_positive(-x)`), Weibull data of variance 0
//!   ("Expected a nonzero variance", checked first), beta data outside `[0,
//!   1]` ("Expected a number in [0,1]"), and a 0 or a 1 among them, a
//!   variance 0 or no root ("Unable to compute parameters for Beta
//!   distribution"; giac answers its start for the 0, `betad(-0.0,-0.0)`);
//!   any other law, the binomial, uniform, lognormal and negative binomial
//!   among them ("Unrecognized probability distribution").
//! - **Limits.** The iterative fits share one budget of `Limits::terms`
//!   evaluations of a datum's terms (a score over all data costs `N`, one
//!   over giac's sums 1; C0019, C0062). Past it, or with no root found (a
//!   step leaving the domain, 100 steps without converging, no
//!   certificate), the call stays, without a warning but for beta's.

#![allow(
    clippy::many_single_char_names,
    reason = "a, b, k, m, v, x: giac's names for the parameters and the moments"
)]

use std::cmp::Ordering;

use cas_num::ball::{Ball, Float, Mag};
use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::numeric::approximation;
use crate::native::random::variates::Law;
use crate::native::statistics::random_variables::{Name, name, object};
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `gentypeerr`.
const TYPE: &str = "Bad Argument Type";
const TWO: &str = "Expected a sequence of two arguments";
const LIST: &str = "Expected a list";
const SAMPLES: &str = "At least two samples are required";
const REAL: &str = "Expected a real constant";
const NONNEGATIVE_INTEGER: &str = "Expected a nonnegative integer";
const POSITIVE: &str = "Expected a positive number";
const POSITIVE_INTEGER: &str = "Expected a positive integer";
const UNIT: &str = "Expected a number in [0,1]";
const BETA: &str = "Unable to compute parameters for Beta distribution";
const VARIANCE: &str = "Expected a nonzero variance";
const UNRECOGNIZED: &str = "Unrecognized probability distribution";

/// Guard bits over the working precision for the fits' balls.
const GUARD: u64 = 32;
/// Newton steps of one solve at most.
const STEPS: usize = 100;
/// Times a Krawczyk box is widened (by 2⁸) before the fit gives up.
const WIDENINGS: usize = 8;

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

fn ball(k: u64) -> Ball {
    Ball::from_integer(&Integer::from(k))
}

/// The laws giac fits.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Fit {
    Normal,
    Poisson,
    Exponential,
    Geometric,
    Gamma,
    Beta,
    Cauchy,
    Weibull,
}

/// The law `f` names, as giac's `_fitdistr` compares it (:7488-7534): the
/// builtins that are its names' aliases, and the capitalised names (C0008).
fn fit_of(f: &Expr) -> Option<Fit> {
    let s = f.as_symbol()?;
    match s.name() {
        "NORMALD" => return Some(Fit::Normal),
        "POISSON" => return Some(Fit::Poisson),
        "EXP" => return Some(Fit::Exponential),
        _ => {}
    }
    if s == Symbol::RANDOM_NORMAL {
        return None;
    }
    Some(match name(f)? {
        Name::Law(Law::Normal) => Fit::Normal,
        Name::Law(Law::Poisson) => Fit::Poisson,
        Name::Law(Law::Exponential) => Fit::Exponential,
        Name::Law(Law::Geometric) => Fit::Geometric,
        Name::Law(Law::Gamma) => Fit::Gamma,
        Name::Law(Law::Beta) => Fit::Beta,
        Name::Law(Law::Cauchy) => Fit::Cauchy,
        Name::Law(Law::Weibull) => Fit::Weibull,
        _ => return None,
    })
}

/// A real number: exact when a rational, and its ball at the fit's
/// precision.
#[derive(Clone, Debug)]
struct Real {
    exact: Option<Rational>,
    ball: Ball,
}

impl Real {
    fn rational(q: Rational, p: u64) -> Real {
        Real {
            ball: Ball::from_rational(&q, p),
            exact: Some(q),
        }
    }

    fn ball(ball: Ball) -> Real {
        Real { exact: None, ball }
    }

    fn expr(&self) -> Expr {
        match &self.exact {
            Some(q) => Expr::from(q),
            None => Expr::from(&self.ball),
        }
    }

    /// Its order against `k`, exactly or by its midpoint (C0028).
    fn against(&self, k: i64) -> Ordering {
        match &self.exact {
            Some(q) => q.cmp(&Rational::from(Integer::from(k))),
            #[allow(clippy::cast_precision_loss, reason = "0 and 1")]
            None => self.ball.mid().to_f64().total_cmp(&(k as f64)),
        }
    }

    fn is_zero(&self) -> bool {
        match &self.exact {
            Some(q) => q.sign() == Ordering::Equal,
            None => self.ball.is_exact() && self.ball.mid().is_zero(),
        }
    }

    /// `1/self`.
    fn inverse(&self, p: u64) -> Real {
        match &self.exact {
            Some(q) => Real::rational(&Rational::from(Integer::ONE) / q, p),
            None => Real::ball(self.ball.inv(p)),
        }
    }
}

/// A datum as giac's `is_real_number` reads it: a rational exactly, a
/// numeric expression by its approximation's midpoint; `None` otherwise.
fn real(ev: &Eval<'_>, x: &Expr, p: u64) -> Option<Real> {
    if let Some(q) = x.to_rational() {
        return Some(Real::rational(q, p));
    }
    if !x.is_numeric() {
        return None;
    }
    match approximation(ev.context(), x, p).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(Real::ball(Ball::exact(b.mid().clone()))),
        _ => None,
    }
}

/// The mean `Σx/N` and the variance `Σ(x − m)²/(N − 1)`, exact on exact
/// data.
fn moments(xs: &[Real], p: u64) -> (Real, Real) {
    let n = xs.len() as u64;
    let exact: Option<Vec<&Rational>> = xs.iter().map(|x| x.exact.as_ref()).collect();
    if let Some(q) = exact {
        let zero = Rational::from(Integer::ZERO);
        let m = &q.iter().fold(zero.clone(), |s, q| &s + *q) / &Rational::from(Integer::from(n));
        let v = q.iter().fold(zero, |s, q| {
            let d = *q - &m;
            &s + &(&d * &d)
        });
        let v = &v / &Rational::from(Integer::from(n - 1));
        return (Real::rational(m, p), Real::rational(v, p));
    }
    let mut s = Ball::ZERO;
    for x in xs {
        s = s.add(&x.ball, p);
    }
    let m = s.div(&ball(n), p);
    let mut v = Ball::ZERO;
    for x in xs {
        let d = x.ball.sub(&m, p);
        v = v.add(&d.mul(&d, p), p);
    }
    (Real::ball(m), Real::ball(v.div(&ball(n - 1), p)))
}

/// A score's values and its Jacobian, at a point or over a box.
type Score = (Vec<Ball>, Vec<Vec<Ball>>);

/// A score function of one or two unknowns, `None` off its domain.
type ScoreFn<'a> = dyn Fn(&[Ball]) -> Option<Score> + 'a;

fn finite(s: Score) -> Option<Score> {
    (s.0.iter().all(Ball::is_finite) && s.1.iter().flatten().all(Ball::is_finite)).then_some(s)
}

/// The inverse of the midpoints of `j` (one or two unknowns), in `f64`, as
/// exact balls.
fn inverse(j: &[Vec<Ball>]) -> Option<Vec<Vec<Ball>>> {
    let f = |b: &Ball| b.mid().to_f64();
    let y: Vec<Vec<f64>> = match j {
        [r] => vec![vec![1.0 / f(&r[0])]],
        [r, s] => {
            let (a, b, c, d) = (f(&r[0]), f(&r[1]), f(&s[0]), f(&s[1]));
            let det = a * d - b * c;
            vec![vec![d / det, -b / det], vec![-c / det, a / det]]
        }
        _ => return None,
    };
    y.iter()
        .map(|r| {
            r.iter()
                .map(|&v| Float::from_f64(v).map(Ball::exact))
                .collect()
        })
        .collect()
}

/// `−Y F`.
fn step(y: &[Vec<Ball>], f: &[Ball], p: u64) -> Vec<Ball> {
    y.iter()
        .map(|row| {
            row.iter()
                .zip(f)
                .fold(Ball::ZERO, |s, (a, b)| s.sub(&a.mul(b, p), p))
        })
        .collect()
}

/// The fits' shared work and precisions.
struct Run {
    /// Evaluations of a datum's terms left.
    left: usize,
    /// The balls' precision.
    p: u64,
    /// The working precision.
    wp: u64,
}

impl Run {
    fn charge(&mut self, cost: usize) -> Option<()> {
        self.left = self.left.checked_sub(cost)?;
        Some(())
    }

    /// Newton's method from `x`, each score costing `cost`: each iterate an
    /// exact midpoint (C0064); `positive` the unknowns that must stay above
    /// 0, `halve` halving one that a step takes past it instead (the gamma
    /// fit). The iterate once every step is below the working precision;
    /// `None` off the domain, past the budget or [`STEPS`].
    fn newton(
        &mut self,
        score: &ScoreFn<'_>,
        cost: usize,
        mut x: Vec<Float>,
        positive: &[bool],
        halve: bool,
    ) -> Option<Vec<Float>> {
        let p = self.p;
        let scale = Integer::from(self.wp + 4);
        let one = Float::from_f64(1.0)?;
        let above = |x: &[Float]| {
            x.iter()
                .zip(positive)
                .all(|(xi, &pos)| !pos || xi.sign() == Ordering::Greater)
        };
        if !above(&x) {
            return None;
        }
        for _ in 0..STEPS {
            self.charge(cost)?;
            let at: Vec<Ball> = x.iter().cloned().map(Ball::exact).collect();
            let (f, j) = finite(score(&at)?)?;
            let d = step(&inverse(&j)?, &f, p);
            let mut done = true;
            let mut next = Vec::with_capacity(x.len());
            for ((xi, di), &pos) in x.iter().zip(&d).zip(positive) {
                let di = di.mid();
                let mut ni = Ball::exact(xi.clone())
                    .add(&Ball::exact(di.clone()), p)
                    .mid()
                    .clone();
                if halve && pos && ni.sign() != Ordering::Greater {
                    ni = xi.mul_2exp(&Integer::from(-1));
                    done = false;
                }
                let size = if xi.is_zero() { &one } else { xi };
                if !di.is_zero() && di.mul_2exp(&scale).cmp_abs(size) == Ordering::Greater {
                    done = false;
                }
                next.push(ni);
            }
            if !above(&next) {
                return None;
            }
            if done {
                return Some(next);
            }
            x = next;
        }
        None
    }

    /// A Krawczyk step (C0061) on the box about `x` of relative radius
    /// `2^(8 − wp)`, widened by 2⁸ up to [`WIDENINGS`] times: `K = x − Y
    /// F(x) + (I − Y J(X))(X − x)` inside `X` certifies a root in `K ∩ X`.
    #[allow(clippy::needless_range_loop, reason = "matrix indices")]
    fn certify(&mut self, score: &ScoreFn<'_>, cost: usize, x: &[Float]) -> Option<Vec<Ball>> {
        let p = self.p;
        let n = x.len();
        self.charge(cost)?;
        let at: Vec<Ball> = x.iter().cloned().map(Ball::exact).collect();
        let (f, j) = finite(score(&at)?)?;
        let y = inverse(&j)?;
        let shift = step(&y, &f, p);
        let one = Float::from_f64(1.0)?;
        let mut e = 8 - i64::try_from(self.wp).ok()?;
        for _ in 0..=WIDENINGS {
            let radii: Vec<Mag> = x
                .iter()
                .map(|xi| {
                    let size = if xi.is_zero() { one.clone() } else { xi.abs() };
                    Mag::from_float(&size.mul_2exp(&Integer::from(e)))
                })
                .collect();
            let boxes: Vec<Ball> = x
                .iter()
                .zip(&radii)
                .map(|(xi, r)| Ball::exact(xi.clone()).add_error(r))
                .collect();
            self.charge(cost)?;
            let (_, jx) = finite(score(&boxes)?)?;
            let mut k = Vec::with_capacity(n);
            for i in 0..n {
                let mut ki = at[i].add(&shift[i], p);
                for c in 0..n {
                    let mut mic = if i == c { ball(1) } else { Ball::ZERO };
                    for l in 0..n {
                        mic = mic.sub(&y[i][l].mul(&jx[l][c], p), p);
                    }
                    ki = ki.add(&mic.mul(&Ball::ZERO.add_error(&radii[c]), p), p);
                }
                k.push(ki);
            }
            if boxes.iter().zip(&k).all(|(b, ki)| b.contains_ball(ki)) {
                return boxes
                    .iter()
                    .zip(&k)
                    .map(|(b, ki)| b.intersection(ki, p))
                    .collect();
            }
            e += 8;
        }
        None
    }
}

/// `Σ f(x)/N`.
fn average(xs: &[Real], f: impl Fn(&Ball) -> Ball, p: u64) -> Ball {
    let mut s = Ball::ZERO;
    for x in xs {
        s = s.add(&f(&x.ball), p);
    }
    s.div(&ball(xs.len() as u64), p)
}

/// The gamma fit (:7505-7514): `(a, a/m)`.
fn gamma(run: &mut Run, xs: &[Real], m: &Real, v: &Real) -> Option<(Ball, Ball)> {
    let p = run.p;
    let c = m.ball.log(p).sub(&average(xs, |x| x.log(p), p), p);
    let score = |a: &[Ball]| {
        let a = &a[0];
        Some((
            vec![a.log(p).sub(&a.digamma(p), p).sub(&c, p)],
            vec![vec![a.inv(p).sub(&a.polygamma(1, p), p)]],
        ))
    };
    let start = m.ball.mul(&m.ball, p).div(&v.ball, p);
    if !start.is_finite() {
        return None;
    }
    let x = run.newton(&score, 1, vec![start.mid().clone()], &[true], true)?;
    let [a] = &run.certify(&score, 1, &x)?[..] else {
        return None;
    };
    Some((a.clone(), a.div(&m.ball, p)))
}

/// The beta fit (:7515-7529), data inside `(0, 1)`.
fn beta(run: &mut Run, xs: &[Real], m: &Real, v: &Real) -> Option<(Ball, Ball)> {
    let p = run.p;
    let one = ball(1);
    let s1 = average(xs, |x| x.log(p), p);
    let s2 = average(xs, |x| one.sub(x, p).log(p), p);
    let score = |x: &[Ball]| {
        let (a, b) = (&x[0], &x[1]);
        let ab = a.add(b, p);
        let (dab, t) = (ab.digamma(p), ab.polygamma(1, p));
        Some((
            vec![
                a.digamma(p).sub(&dab, p).sub(&s1, p),
                b.digamma(p).sub(&dab, p).sub(&s2, p),
            ],
            vec![
                vec![a.polygamma(1, p).sub(&t, p), t.neg()],
                vec![t.neg(), b.polygamma(1, p).sub(&t, p)],
            ],
        ))
    };
    let (m, v) = (&m.ball, &v.ball);
    let fac = v.add(&m.mul(m, p), p).sub(m, p).div(v, p);
    let a = m.neg().mul(&fac, p);
    let b = m.sub(&one, p).mul(&fac, p);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let start = vec![a.mid().clone(), b.mid().clone()];
    let x = run.newton(&score, 1, start, &[true, true], false)?;
    let [a, b] = &run.certify(&score, 1, &x)?[..] else {
        return None;
    };
    Some((a.clone(), b.clone()))
}

/// giac's lower quantile `v[⌈N k/4⌉]` of the sorted data (misc.cc:3108-3171).
fn quarter(sorted: &[Real], k: usize) -> &Real {
    &sorted[(k * sorted.len()).div_ceil(4) - 1]
}

/// The Cauchy fit (`cauchy_mle`, :7433-7457), from the median and half
/// the quartiles' gap.
fn cauchy(run: &mut Run, xs: &[Real]) -> Option<(Ball, Ball)> {
    let p = run.p;
    let mut sorted = xs.to_vec();
    if sorted.iter().all(|x| x.exact.is_some()) {
        sorted.sort_by(|a, b| a.exact.cmp(&b.exact));
    } else {
        sorted.sort_by(|a, b| a.ball.mid().to_f64().total_cmp(&b.ball.mid().to_f64()));
    }
    let two = ball(2);
    let half = ball(xs.len() as u64).div(&two, p);
    let score = |x: &[Ball]| {
        let (x0, g) = (&x[0], &x[1]);
        let gsq = g.mul(g, p);
        let (mut f, mut gg) = (Ball::ZERO, half.neg());
        let (mut fx0, mut fg, mut ggam) = (Ball::ZERO, Ball::ZERO, Ball::ZERO);
        for xi in xs {
            let d = xi.ball.sub(x0, p);
            let dsq = d.mul(&d, p);
            let den = gsq.add(&dsq, p);
            let densq = den.mul(&den, p);
            f = f.add(&d.div(&den, p), p);
            gg = gg.add(&gsq.div(&den, p), p);
            fx0 = fx0.add(&dsq.sub(&gsq, p).div(&densq, p), p);
            fg = fg.add(&d.div(&densq, p), p);
            ggam = ggam.add(&dsq.div(&densq, p), p);
        }
        let gx0 = two.mul(&gsq, p).mul(&fg, p);
        let fgam = two.mul(g, p).mul(&fg, p).neg();
        let ggam = two.mul(g, p).mul(&ggam, p);
        Some((vec![f, gg], vec![vec![fx0, fgam], vec![gx0, ggam]]))
    };
    let x0 = quarter(&sorted, 2).ball.mid().clone();
    let g = quarter(&sorted, 3)
        .ball
        .sub(&quarter(&sorted, 1).ball, p)
        .div(&two, p);
    let cost = xs.len();
    let x = run.newton(
        &score,
        cost,
        vec![x0, g.mid().clone()],
        &[false, false],
        false,
    )?;
    let [x0, g] = &run.certify(&score, cost, &x)?[..] else {
        return None;
    };
    Some((x0.clone(), g.abs()))
}

/// The Weibull fit (:7534-7543 and `weibull_mle`): giac's `k_init`, then
/// the profile score's root `k` and `λ = (Σxᵏ/N)^(1/k)`; data `> 0`.
fn weibull(run: &mut Run, xs: &[Real], m: &Real, v: &Real) -> Option<(Ball, Ball)> {
    let p = run.p;
    let one = ball(1);
    let two = ball(2);
    let n = ball(xs.len() as u64);
    let logs: Vec<Ball> = xs.iter().map(|x| x.ball.log(p)).collect();
    let mean_log = average(xs, |x| x.log(p), p);
    let rhs = v.ball.div(&m.ball.mul(&m.ball, p), p).add(&one, p);
    let moment = |k: &[Ball]| {
        let k = &k[0];
        let t1 = one.add(&k.inv(p), p);
        let t2 = one.add(&two.div(k, p), p);
        let r = t2.gamma(p).div(&t1.gamma(p), p);
        let dr = r
            .mul(&t1.digamma(p).sub(&two.mul(&t2.digamma(p), p), p), p)
            .div(&k.mul(k, p), p);
        Some((vec![r.sub(&rhs, p)], vec![vec![dr]]))
    };
    let sums = |k: &Ball| {
        let (mut s0, mut s1, mut s2) = (Ball::ZERO, Ball::ZERO, Ball::ZERO);
        for l in &logs {
            let t = k.mul(l, p).exp(p);
            let tl = t.mul(l, p);
            s2 = s2.add(&tl.mul(l, p), p);
            s1 = s1.add(&tl, p);
            s0 = s0.add(&t, p);
        }
        (s0, s1, s2)
    };
    let profile = |k: &[Ball]| {
        let k = &k[0];
        let (s0, s1, s2) = sums(k);
        let ik = k.inv(p);
        let w = ik.add(&mean_log, p);
        Some((
            vec![s1.sub(&s0.mul(&w, p), p)],
            vec![vec![
                s2.sub(&s1.mul(&w, p), p)
                    .add(&s0.mul(&ik.mul(&ik, p), p), p),
            ]],
        ))
    };
    let inv = v.ball.inv(p);
    let start = if inv.is_finite() && inv.mid().to_f64() > 1.0 {
        inv.mid().clone()
    } else {
        Float::from_f64(1.0)?
    };
    let k0 = run
        .newton(&moment, 1, vec![start.clone()], &[true], false)
        .unwrap_or_else(|| vec![start]);
    let cost = xs.len();
    let x = run.newton(&profile, cost, k0, &[true], false)?;
    let [k] = &run.certify(&profile, cost, &x)?[..] else {
        return None;
    };
    run.charge(cost)?;
    let (s0, _, _) = sums(k);
    let lambda = s0.div(&n, p).log(p).div(k, p).exp(p);
    lambda.is_finite().then(|| (k.clone(), lambda))
}

/// `FitDistribution[L, law]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's checks and laws in its order, one arm each"
)]
pub(in crate::native) fn fit_distribution(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (data, law) = match e.args() {
        [d, l] => (d, l),
        [_, _, _, ..] => return refuse(ev, e, TWO),
        _ => return refuse(ev, e, TYPE),
    };
    if !is_list(data) {
        return refuse(ev, e, LIST);
    }
    let given = data.args();
    if given.len() < 2 {
        return refuse(ev, e, SAMPLES);
    }
    let wp = ev.context().working_precision();
    let p = wp + GUARD;
    let Some(xs) = given
        .iter()
        .map(|x| real(ev, x, p))
        .collect::<Option<Vec<_>>>()
    else {
        return refuse(ev, e, REAL);
    };
    let (m, v) = moments(&xs, p);
    let Some(fit) = fit_of(law) else {
        return refuse(ev, e, UNRECOGNIZED);
    };
    let integers = |least: i64| {
        given
            .iter()
            .all(|x| x.to_integer().is_some_and(|k| k >= Integer::from(least)))
    };
    let positive = || xs.iter().all(|x| x.against(0) == Ordering::Greater);
    let mut run = Run {
        left: ev.context().limits().terms,
        p,
        wp,
    };
    let pair = |r: Option<(Ball, Ball)>| r.map(|(a, b)| vec![Expr::from(&a), Expr::from(&b)]);
    let (law, parameters) = match fit {
        Fit::Normal => {
            let s = match &v.exact {
                Some(q) => ev.eval(&Expr::call(Symbol::SQRT, [Expr::from(q)]))?,
                None => Expr::from(&v.ball.sqrt(p)),
            };
            (Law::Normal, vec![m.expr(), s])
        }
        Fit::Poisson => {
            if !integers(0) {
                return refuse(ev, e, NONNEGATIVE_INTEGER);
            }
            (Law::Poisson, vec![m.expr()])
        }
        Fit::Exponential => {
            if !positive() {
                return refuse(ev, e, POSITIVE);
            }
            (Law::Exponential, vec![m.inverse(p).expr()])
        }
        Fit::Geometric => {
            if !integers(1) {
                return refuse(ev, e, POSITIVE_INTEGER);
            }
            (Law::Geometric, vec![m.inverse(p).expr()])
        }
        Fit::Gamma => {
            if !positive() {
                return refuse(ev, e, POSITIVE);
            }
            let Some(ps) = pair(gamma(&mut run, &xs, &m, &v)) else {
                return Ok(None);
            };
            (Law::Gamma, ps)
        }
        Fit::Beta => {
            if xs
                .iter()
                .any(|x| x.against(0) == Ordering::Less || x.against(1) == Ordering::Greater)
            {
                return refuse(ev, e, UNIT);
            }
            let inside = xs
                .iter()
                .all(|x| x.against(0) == Ordering::Greater && x.against(1) == Ordering::Less);
            let found = if inside && !v.is_zero() {
                pair(beta(&mut run, &xs, &m, &v))
            } else {
                None
            };
            let Some(ps) = found else {
                return refuse(ev, e, BETA);
            };
            (Law::Beta, ps)
        }
        Fit::Cauchy => {
            let Some(ps) = pair(cauchy(&mut run, &xs)) else {
                return Ok(None);
            };
            (Law::Cauchy, ps)
        }
        Fit::Weibull => {
            if v.is_zero() {
                return refuse(ev, e, VARIANCE);
            }
            if !positive() {
                return refuse(ev, e, POSITIVE);
            }
            let Some(ps) = pair(weibull(&mut run, &xs, &m, &v)) else {
                return Ok(None);
            };
            (Law::Weibull, ps)
        }
    };
    Ok(Some(object(law, parameters)))
}
