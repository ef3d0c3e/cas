//! `Newton` and `NSolve` (CAS-591), Xcas's `newton` and `nSolve`: a root
//! of one equation by giac's damped Newton iteration with restarts.
//!
//! - `Newton[f, x]`, `Newton[f, x == a]`, `Newton[f, x, a, p, eps]` ports
//!   giac's `newton` (solve.cc `newton`, `_newton`): an equation is
//!   `lhs − rhs`; the start `a`, a range's midpoint, or restart points;
//!   `p` (an integer, giac's 20 by default) the steps of each loop, `eps`
//!   (another number) the second loop's tolerance. A first loop of damped
//!   steps `a + λ d`, `d = −f(a)/f'(a)`, `λ` halved when `|f|` grows, else
//!   times 6/5 up to 1, until `|d|` is below `10⁻¹⁰` of `max(1, |b|)`; a
//!   second of plain steps until below `eps`. A loop that runs out
//!   restarts with 5 steps more, up to 20 times; past the fifth restart
//!   the answer is giac's `undef`, even when the iteration converged.
//! - The arithmetic is exact rationals rounded to the working precision
//!   and 32 bits more, `f` and its derivative by `D` their balls'
//!   midpoints there; giac's decisions compare those values.
//! - A converged iterate is a ball certified by Krawczyk's test
//!   ([`systems::certified_root`]); when it is not, the call stays with a
//!   warning, as `FSolve`'s. With `p` given, a first loop that runs out
//!   answers giac's iterate as a ball at the working precision: no root is
//!   claimed.
//! - `NSolve[eqn, x]`, `NSolve[eqn, x == a]` (giac's ti89.cc `_nSolve`):
//!   `Newton` from 0 or `a`, the first loop to `10⁻⁵`; `x == a .. b` is
//!   `FSolve`'s, a list.
//! - Departures: giac restarts from random points, here the fixed
//!   `(−1)^⌊(j−1)/2⌋ 3j/4` of restart `j`; giac's `p` equal to 20 behaves
//!   as no `p`, here any `p` given stops. `NSolve` of a system (a list of
//!   equations, a list of unknowns) is `FSolve`'s without a start, its
//!   first real solution, where giac's depends on its random restarts.
//! - Complex `f` or start, parameters, a `p` past [`MAX_ITERATIONS`] or
//!   whose product with the size of `f` passes the term limit leave the
//!   call.

use cas_num::{Integer, Rational};

use super::fsolve::{self, at, numeric, point, range, rounded, systems};
use super::sides;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `NEWTON_DEFAULT_ITERATION`: the steps of each loop.
const DEFAULT_ITERATIONS: usize = 20;

/// giac's `NEWTON_MAX_RANDOM_RESTART`.
const RESTARTS: usize = 20;

/// The most steps of each loop a call may ask for.
const MAX_ITERATIONS: usize = 1000;

/// giac's tolerances and its iteration count.
struct Settings {
    niter: usize,
    /// Whether a first loop that runs out answers its iterate.
    out: bool,
    eps1: Rational,
    eps2: Rational,
}

/// How the iteration ended.
enum Outcome {
    Converged(Rational),
    Stopped(Rational),
    Undefined,
}

/// One start's attempt.
enum Attempt {
    Converged(Rational),
    Stopped(Rational),
    Failed,
}

/// A step `d`, or giac's `inv(0)` times a nonzero `f`, or an undefined
/// one.
enum Step {
    Finite(Rational),
    Infinite,
    Undefined,
}

/// `Newton[f, x]`, `Newton[f, x == a]` or `Newton[f, x, a, p, eps]`.
pub(in crate::native) fn newton(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (eqn, x, guess, rest) = match e.args() {
        [eqn, spec] => match sides(spec) {
            Some((x, a)) => (eqn, x, Some(a), &[][..]),
            None => (eqn, spec, None, &[][..]),
        },
        [eqn, x, a, rest @ ..] if rest.len() <= 2 => (eqn, x, Some(a), rest),
        _ => return Ok(None),
    };
    let mut settings = Settings {
        niter: DEFAULT_ITERATIONS,
        out: false,
        eps1: power_of_ten(10),
        eps2: power_of_ten(12),
    };
    for arg in rest {
        if let Some(p) = arg.to_integer() {
            match usize::try_from(p) {
                Ok(p) if (1..=MAX_ITERATIONS).contains(&p) => {
                    settings.niter = p;
                    settings.out = true;
                }
                _ => return Ok(None),
            }
        } else {
            match point(ev.context(), arg) {
                Some(eps) if eps > Rational::ZERO => settings.eps2 = eps,
                _ => return Ok(None),
            }
        }
    }
    solved(ev, e, "Newton", eqn, x, guess, settings)
}

/// `NSolve[eqn, x]`, `NSolve[eqn, x == a]`, `NSolve[eqn, x == a .. b]`, or
/// of a system.
pub(in crate::native) fn nsolve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [eqn, var] = e.args() else {
        return Ok(None);
    };
    if head(eqn) == Some(Symbol::LIST) {
        if head(var) != Some(Symbol::LIST)
            || var.args().len() != eqn.args().len()
            || var.args().is_empty()
            || var.args().len() > systems::MAX_UNKNOWNS
        {
            return Ok(None);
        }
        let call = Expr::call(Symbol::FSOLVE, [eqn.clone(), var.clone()]);
        let Some(found) = systems::system(ev, &call)? else {
            return Ok(None);
        };
        // A list of tuples, or Newton's one tuple of numbers.
        let first = match found.args().first() {
            Some(t) if head(t) == Some(Symbol::LIST) => t.clone(),
            Some(_) => found,
            None => return Ok(None),
        };
        ev.step("NSolve", e, &first);
        return Ok(Some(first));
    }
    let (x, guess) = match sides(var) {
        Some((_, a)) if range(a).is_some() => {
            let call = Expr::call(Symbol::FSOLVE, [eqn.clone(), var.clone()]);
            let found = fsolve::fsolve(ev, &call)?;
            if let Some(r) = &found {
                ev.step("NSolve", e, r);
            }
            return Ok(found);
        }
        Some((x, a)) => (x, a.clone()),
        None => (var, Expr::from(Integer::from(0))),
    };
    let settings = Settings {
        niter: DEFAULT_ITERATIONS,
        out: false,
        eps1: power_of_ten(5),
        eps2: power_of_ten(12),
    };
    solved(ev, e, "NSolve", eqn, x, Some(&guess), settings)
}

/// The answer of giac's `newton` to `eqn` in `x` from `guess`.
fn solved(
    ev: &mut Eval<'_>,
    e: &Expr,
    name: &'static str,
    eqn: &Expr,
    x: &Expr,
    guess: Option<&Expr>,
    settings: Settings,
) -> Native {
    if x.as_symbol().is_none() || has_i(eqn) || guess.is_some_and(has_i) {
        return Ok(None);
    }
    let size = usize::try_from(eqn.size()).unwrap_or(usize::MAX);
    if settings.niter.saturating_mul(size) > ev.context().limits().terms {
        return Ok(None);
    }
    let start = match guess {
        None => None,
        Some(g) => {
            let ctx = ev.context();
            let q = match range(g) {
                Some((a, b)) => match (point(ctx, a), point(ctx, b)) {
                    (Some(a), Some(b)) => &(&a + &b) / &Rational::from(Integer::from(2)),
                    _ => return Ok(None),
                },
                None => match point(ctx, g) {
                    Some(q) => q,
                    None => return Ok(None),
                },
            };
            Some(q)
        }
    };
    let f = match sides(eqn) {
        Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
        None => eqn.clone(),
    };
    if has_i(&f) || !numeric(ev.context(), &f, x) {
        return Ok(None);
    }
    let df = ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?;
    let bits = ev.context().working_precision();
    let outcome = iterate(ev.context(), &f, &df, x, start, settings, bits);
    let result = match outcome {
        Outcome::Undefined => Expr::from(Symbol::INDETERMINATE),
        Outcome::Stopped(q) => {
            match super::super::approximation(ev.context(), &Expr::from(q), bits) {
                Ok(b) => b,
                Err(_) => return Ok(None),
            }
        }
        Outcome::Converged(q) => {
            let root =
                systems::certified_root(ev, vec![f], std::slice::from_ref(x), vec![q.clone()])?;
            let Some(r) = root.and_then(|mut r| r.pop()) else {
                ev.warn(Warning::Native {
                    message: "Newton's iteration found no root it could certify".into(),
                    expr: Expr::from(q),
                });
                return Ok(None);
            };
            r
        }
    };
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// giac's `newton` loop (solve.cc), restarts at [`start`].
fn iterate(
    ctx: &Context,
    f: &Expr,
    df: &Expr,
    x: &Expr,
    mut guess: Option<Rational>,
    mut settings: Settings,
    bits: u64,
) -> Outcome {
    let prec = bits.saturating_add(32);
    let mut j = 1;
    let mut root = None;
    while j <= RESTARTS {
        let a = guess.take().unwrap_or_else(|| start(j));
        match attempt(ctx, f, df, x, a, j, &settings, prec) {
            Attempt::Stopped(a) => return Outcome::Stopped(a),
            Attempt::Converged(a) => {
                root = Some(a);
                break;
            }
            Attempt::Failed => {}
        }
        j += 1;
        settings.niter += 5;
    }
    match root {
        Some(a) if j <= 5 => Outcome::Converged(a),
        _ => Outcome::Undefined,
    }
}

/// One start's two loops; a value or derivative that is no real number
/// where giac would throw fails it.
#[allow(clippy::too_many_arguments, reason = "the iteration's state")]
fn attempt(
    ctx: &Context,
    f: &Expr,
    df: &Expr,
    x: &Expr,
    mut a: Rational,
    j: usize,
    settings: &Settings,
    prec: u64,
) -> Attempt {
    let value = |q: &Rational| at(ctx, f, x, q, prec);
    let one = Rational::from(Integer::from(1));
    let scale = |q: &Rational| q.abs().max(one.clone());
    let niter = settings.niter;
    let Some(mut fa) = value(&a) else {
        return Attempt::Failed;
    };
    let mut lambda = one.clone();
    let mut k = 0;
    let mut broke = false;
    while k < niter {
        let step = match at(ctx, df, x, &a, prec) {
            Some(s) if !s.is_zero() => Step::Finite(rounded(&-(&fa / &s), prec)),
            Some(_) if !fa.is_zero() => Step::Infinite,
            _ => Step::Undefined,
        };
        let d = match step {
            Step::Finite(d) if !(k == 0 && d.is_zero() && &fa * &fa > settings.eps2) => d,
            // giac's infinite `b` ends the loop as run out.
            Step::Infinite => break,
            // giac's restart inside the loop, which takes a step.
            _ => {
                a = start(j);
                let Some(v) = value(&a) else {
                    return Attempt::Failed;
                };
                fa = v;
                k += 1;
                continue;
            }
        };
        let b = rounded(&(&a + &(&lambda * &d)), prec);
        if d.abs() < &settings.eps1 * &scale(&b) {
            a = b;
            broke = true;
            break;
        }
        match value(&b) {
            Some(fb) if fb.abs() <= fa.abs() => {
                lambda = if lambda > Rational::new(Integer::from(9), Integer::from(10)) {
                    one.clone()
                } else {
                    rounded(
                        &(&lambda * &Rational::new(Integer::from(6), Integer::from(5))),
                        prec,
                    )
                };
                a = b;
                fa = fb;
            }
            _ => lambda = &lambda / &Rational::from(Integer::from(2)),
        }
        k += 1;
    }
    if !broke {
        return if settings.out {
            Attempt::Stopped(a)
        } else {
            Attempt::Failed
        };
    }
    let mut d = Rational::ZERO;
    let mut k = 0;
    while k < niter {
        let (Some(v), Some(s)) = (value(&a), at(ctx, df, x, &a, prec)) else {
            return Attempt::Failed;
        };
        if s.is_zero() {
            return Attempt::Failed;
        }
        d = rounded(&-(&v / &s), prec);
        a = rounded(&(&a + &d), prec);
        if d.abs() < &settings.eps2 * &scale(&a) {
            break;
        }
        k += 1;
    }
    if k != niter || d.abs() < &settings.eps1 * &scale(&a) {
        Attempt::Converged(a)
    } else {
        Attempt::Failed
    }
}

/// Restart `j`'s point, `(−1)^⌊(j−1)/2⌋ 3j/4`, in place of giac's random
/// one.
fn start(j: usize) -> Rational {
    let sign = if (j - 1) / 2 % 2 == 0 { 1 } else { -1 };
    let j = i64::try_from(j).unwrap_or(i64::MAX);
    Rational::new(Integer::from(sign * 3 * j), Integer::from(4))
}

/// `10^−n`.
fn power_of_ten(n: i64) -> Rational {
    Rational::from(Integer::from(10)).pow(-n)
}

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

/// Whether `e` has a complex number or `I`.
fn has_i(e: &Expr) -> bool {
    e.distinct().any(|z| {
        matches!(z.view(), View::Complex(_) | View::ComplexBall(_))
            || z.as_symbol() == Some(Symbol::I)
    })
}
