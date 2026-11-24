//! giac's `simulated_annealing` (CAS-610), after optimization.cc `_siman`
//! and GSL 2.8 siman/siman.c `gsl_siman_solve`: simulated annealing (S.
//! Kirkpatrick, C. D. Gelatt and M. P. Vecchi, *Optimization by simulated
//! annealing*, Science 220 (1983) 671-680). A port of GPL code in our own
//! words.
//!
//! - **`SimulatedAnnealing[x0, cost, dist, step, {k, t_initial, mu_t,
//!   t_min}(, {n_tries, iters_fixed_T}(, maxstep))]`**: `cost`, `dist` and
//!   `step` functions with named parameters (`Function[x, body]` or
//!   `Function[{x, …}, body]`), of 1, 2 and 1 parameters, `step` of 2 when
//!   `maxstep` is given; the cooling schedule four real numbers;
//!   `n_tries` and `iters_fixed_T` integers at least 1 (`10` and `100`);
//!   `maxstep` a real number above 0. Anything else leaves the call.
//! - Before any work, with a warning and the call left: `k > 0`,
//!   `mu_t > 1`, `0 < t_min < t_initial`, and the steps,
//!   `iters_fixed_T` times the number of temperatures (at most
//!   `⌈ln(t_initial/t_min)/ln mu_t⌉ + 1`), at most the limit `terms`
//!   (C0019).
//! - GSL's loop: `E = cost(x0)`, the best `x0`; at each temperature `T`
//!   (from `t_initial`), `iters_fixed_T` times: `new = step(x)` (or
//!   `step(x, maxstep)`), `cost(new)` at most the best's makes `new` the
//!   best, below `E` makes it `x`, else it becomes `x` when a uniform draw
//!   is below `exp(−(cost(new) − E)/(k·T))`; then `T` is divided by `mu_t`,
//!   and the loop ends once `T < t_min`. The answer: the best point, as
//!   `step` gave it, but its decimals made exact at their midpoints (GSL's
//!   doubles; radii would compound from step to step). `dist` is never applied (GSL's neither); `n_tries` is
//!   unused (GSL's too).
//! - `cost`'s value, numericised (`N`, giac's `evalf_double`), must be a
//!   finite real number (a decimal by its midpoint), and `step`'s value not
//!   `Indeterminate` or `ComplexInfinity`; else a warning and the call
//!   stays. The host's interrupt is polled at each step.
//!
//! Departures from giac:
//!
//! - The draws, ours and those of the user's functions (`Random`), come
//!   from the session's generator, so `SeedRandom[s]` makes a run repeat:
//!   not GSL's default generator, nor giac's `rand`.
//! - giac's errors in `cost` or `step` leave the call here.
//! - The schedule's checks are ours: GSL loops forever when `mu_t ≤ 1`.

use cas_num::Rational;
use cas_num::ball::{Ball, ComplexBall, Float, Round};
use cas_num::random::Generator;

use super::find_minimum::{interrupted, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// giac's defaults: `n_tries` (unused) and `iters_fixed_T`.
const ITERS_FIXED_T: usize = 100;

const K: &str = "simulated_annealing: k should be positive";
const MU_T: &str = "simulated_annealing: mu_t should be greater than 1";
const T_MIN: &str = "simulated_annealing: t_min should be positive and below t_initial";
const LIMITS: &str = "simulated_annealing: past the limits";
const COST: &str = "simulated_annealing: the cost should be a real number";
const STEP: &str = "simulated_annealing: the step gave no point";

fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// A finite real number as a double: a rational or a finite decimal's
/// midpoint.
fn real(e: &Expr) -> Option<f64> {
    let v = match e.view() {
        View::Ball(b) if b.is_finite() => b.mid().to_f64(),
        View::Ball(_) => return None,
        _ => to_f64(&e.to_rational()?),
    };
    v.is_finite().then_some(v)
}

/// The number of named parameters of `Function[x, body]` or
/// `Function[{x, …}, body]`.
fn arity(f: &Expr) -> Option<usize> {
    if f.head().and_then(Expr::as_symbol) != Some(Symbol::FUNCTION) {
        return None;
    }
    let [params, _] = f.args() else {
        return None;
    };
    match params.view() {
        View::Symbol(_) => Some(1),
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::LIST) => args
            .iter()
            .all(|a| a.as_symbol().is_some())
            .then_some(args.len()),
        _ => None,
    }
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// An integer at least 1, saturated to `usize::MAX`.
fn count(e: &Expr) -> Option<usize> {
    let n = e.to_integer()?;
    if !n.is_positive() {
        return None;
    }
    Some(
        n.to_i64()
            .and_then(|k| usize::try_from(k).ok())
            .unwrap_or(usize::MAX),
    )
}

/// `cost(x)` numericised as a double, or `None` with the warning.
fn energy(ev: &mut Eval<'_>, e: &Expr, cost: &Expr, x: &Expr) -> Result<Option<f64>, EvalError> {
    let value = ev.eval(&Expr::call(
        Symbol::N,
        [Expr::apply(cost, std::slice::from_ref(x))],
    ))?;
    let v = real(&value);
    if v.is_none() {
        warn(ev, e, COST);
    }
    Ok(v)
}

/// `e` with each ball made its exact midpoint: GSL's point, a double. A
/// ball's radius would compound through `step` (in `Random[x − m, x + m]`
/// both ends carry `x`'s radius, so the draw's is about `r(1 + 2u)`), and
/// `N[cost[x]]` would warn of it at every step (`InputLimited`); once it
/// passes the cost's scale, the cost's midpoint is no longer its value.
fn points(e: &Expr) -> Expr {
    match e.view() {
        View::Ball(b) => Expr::from(Ball::exact(b.mid().clone())),
        View::ComplexBall(z) => Expr::from(ComplexBall::new(
            Ball::exact(z.re().mid().clone()),
            Ball::exact(z.im().mid().clone()),
        )),
        View::Apply { head, args } => {
            Expr::apply(head, &args.iter().map(points).collect::<Vec<_>>())
        }
        _ => e.clone(),
    }
}

/// GSL's `gsl_rng_uniform`: a double in `[0, 1)` from 53 drawn bits.
fn uniform(ev: &mut Eval<'_>) -> f64 {
    let bits = ev.draw(Generator::next_u64) >> 11;
    #[allow(clippy::cast_precision_loss, reason = "53 bits are exact")]
    let u = bits as f64;
    u * f64::powi(2.0, -53)
}

/// `SimulatedAnnealing[x0, cost, dist, step, schedule(, iterations(,
/// maxstep))]`.
pub(in crate::native) fn simulated_annealing(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x0, cost, dist, step, schedule, rest @ ..] = e.args() else {
        return Ok(None);
    };
    if rest.len() > 2 || arity(cost) != Some(1) || arity(dist) != Some(2) {
        return Ok(None);
    }
    let maxstep = rest.get(1);
    if arity(step) != Some(if maxstep.is_some() { 2 } else { 1 }) {
        return Ok(None);
    }
    if !is_list(schedule) {
        return Ok(None);
    }
    let Some([k, t_initial, mu_t, t_min]) = schedule
        .args()
        .iter()
        .map(real)
        .collect::<Option<Vec<_>>>()
        .and_then(|v| <[f64; 4]>::try_from(v).ok())
    else {
        return Ok(None);
    };
    let iters = match rest.first() {
        None => ITERS_FIXED_T,
        Some(l) => match l.args() {
            // n_tries is parsed and unused, as GSL's.
            [n_tries, iters] if is_list(l) && count(n_tries).is_some() => match count(iters) {
                Some(n) => n,
                None => return Ok(None),
            },
            _ => return Ok(None),
        },
    };
    if let Some(m) = maxstep {
        if !real(m).is_some_and(|m| m > 0.0) {
            return Ok(None);
        }
    }
    // The schedule's checks and the work bounded, before any (C0019).
    let checks = [
        (k > 0.0, K),
        (mu_t > 1.0, MU_T),
        (t_min > 0.0 && t_min < t_initial, T_MIN),
    ];
    if let Some((_, message)) = checks.iter().find(|(ok, _)| !ok) {
        warn(ev, e, message);
        return Ok(None);
    }
    let temperatures = ((t_initial / t_min).ln() / mu_t.ln()).ceil() + 1.0;
    #[allow(clippy::cast_precision_loss, reason = "a bound in doubles")]
    let steps = temperatures * iters as f64;
    #[allow(clippy::cast_precision_loss, reason = "a bound in doubles")]
    let terms = ev.context().limits().terms as f64;
    if !steps.is_finite() || steps > terms {
        warn(ev, e, LIMITS);
        return Ok(None);
    }
    // GSL's loop (siman.c 53-179).
    let Some(mut energy_x) = energy(ev, e, cost, x0)? else {
        return Ok(None);
    };
    let mut x = x0.clone();
    let mut best = x0.clone();
    let mut energy_best = energy_x;
    let mut temperature = t_initial;
    let factor = 1.0 / mu_t;
    loop {
        for _ in 0..iters {
            interrupted(ev)?;
            let call = match maxstep {
                Some(size) => Expr::apply(step, &[x.clone(), size.clone()]),
                None => Expr::apply(step, std::slice::from_ref(&x)),
            };
            let new = points(&ev.eval(&call)?);
            if matches!(
                new.as_symbol(),
                Some(Symbol::INDETERMINATE | Symbol::COMPLEX_INFINITY)
            ) {
                warn(ev, e, STEP);
                return Ok(None);
            }
            let Some(energy_new) = energy(ev, e, cost, &new)? else {
                return Ok(None);
            };
            if energy_new <= energy_best {
                best = new.clone();
                energy_best = energy_new;
            }
            if energy_new < energy_x
                || uniform(ev) < (-(energy_new - energy_x) / (k * temperature)).exp()
            {
                x = new;
                energy_x = energy_new;
            }
        }
        temperature *= factor;
        if temperature < t_min {
            break;
        }
    }
    Ok(Some(best))
}
