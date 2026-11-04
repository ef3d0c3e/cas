//! Iterative linear solvers (CAS-DOC-51 §5, CAS-499): Xcas's
//! `jacobi_linsolve` and `gauss_seidel_linsolve`.
//!
//! - **`JacobiLinearSolve[A, b, ε, k]`**: `x_{m+1} = D⁻¹ (b − (L + U) x_m)`
//!   from `x₀ = 0`, `k` steps at most (1000 by default), until `‖x_{m+1} −
//!   x_m‖∞ < ε`.
//! - **`GaussSeidelLinearSolve[ω, A, b, ε, k]`**: the forward sweep, each
//!   entry from the ones just found, relaxed by `ω` (1 by default, giac's
//!   leading argument), `0 < ω < 2`.
//!
//! The iterations run on the midpoints of balls at the working precision,
//! for real entries. The answer is certified by Varah's bound: for `A`
//! strictly diagonally dominant by rows, `δ = minᵢ (|aᵢᵢ| − Σⱼ≠ᵢ |aᵢⱼ|) >
//! 0` and `‖x − x̃‖∞ ≤ ‖b − A x̃‖∞ / δ`, the residual and `δ` in balls, so
//! that each entry's ball holds the solution of every system in the balls
//! of `A` and `b`. A matrix not shown dominant, or an iteration that does
//! not settle, leaves the call with a warning; giac answers either way, and
//! its `jacobi_linsolve([[4,1],[1,3]],[1,2],1e-12)` is `[0.1333,0.4667]`,
//! where the solution is `[1/11, 7/11]`.

use std::cmp::Ordering;

use cas_num::ball::{Ball, Mag};

use super::{rows, warn};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::approximation;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The default bound on the number of steps, and the largest one taken.
const STEPS: u64 = 1000;
const MAX_STEPS: u64 = 1_000_000;

/// The real ball of `x` at `bits` bits; `None` for another value.
fn ball(ev: &mut Eval<'_>, x: &Expr, bits: u64) -> Result<Option<Ball>, EvalError> {
    if let Some(q) = x.to_rational() {
        return Ok(Some(Ball::from_rational(&q, bits)));
    }
    let value = approximation(ev.context(), x, bits)?;
    Ok(match value.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    })
}

/// A lower bound of `|x|` over the ball.
fn lower(b: &Ball) -> Mag {
    Mag::from_float_lower(b.mid()).sub_lower(b.rad())
}

/// An upper bound of `|x|` over the ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// The system: `A`'s rows and `b`, as balls.
struct System {
    a: Vec<Vec<Ball>>,
    b: Vec<Ball>,
    bits: u64,
}

impl System {
    fn new(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<System>, EvalError> {
        let bits = ev.context().working_precision().max(64);
        let Some(m) = rows(a) else {
            return Ok(None);
        };
        let n = m.len();
        if m[0].len() != n || b.head() != Some(&Expr::from(Symbol::LIST)) || b.args().len() != n {
            return Ok(None);
        }
        let mut rows = Vec::with_capacity(n);
        for row in &m {
            let mut entries = Vec::with_capacity(n);
            for x in row {
                let Some(x) = ball(ev, x, bits)? else {
                    return Ok(None);
                };
                entries.push(x);
            }
            rows.push(entries);
        }
        let mut rhs = Vec::with_capacity(n);
        for x in b.args() {
            let Some(x) = ball(ev, x, bits)? else {
                return Ok(None);
            };
            rhs.push(x);
        }
        Ok(Some(System {
            a: rows,
            b: rhs,
            bits,
        }))
    }

    /// `bᵢ − Σⱼ≠ᵢ aᵢⱼ xⱼ`.
    fn off_diagonal(&self, i: usize, x: &[Ball]) -> Ball {
        let mut s = self.b[i].clone();
        for (j, (a, y)) in self.a[i].iter().zip(x).enumerate() {
            if j != i {
                s.sub_assign(&a.mul(y, self.bits), self.bits);
            }
        }
        s
    }

    /// Varah's lower bound `δ` of `‖A⁻¹‖∞⁻¹`; `None` unless `A` is shown
    /// strictly diagonally dominant.
    fn dominance(&self) -> Option<Mag> {
        let mut least: Option<Mag> = None;
        for (i, row) in self.a.iter().enumerate() {
            let mut d = row[i].abs();
            for (j, a) in row.iter().enumerate() {
                if j != i {
                    d.sub_assign(&a.abs(), self.bits);
                }
            }
            if d.sign() != Some(Ordering::Greater) {
                return None;
            }
            let d = lower(&d);
            least = Some(match least {
                Some(m) if m <= d => m,
                _ => d,
            });
        }
        least
    }

    /// `x̃` with each entry widened by Varah's bound; `None` without one.
    fn certified(&self, x: &[Ball]) -> Option<Vec<Ball>> {
        let delta = self.dominance()?;
        let mut residual = Mag::from_integer(&cas_num::Integer::ZERO);
        for i in 0..x.len() {
            let mut r = self.off_diagonal(i, x);
            r.sub_assign(&self.a[i][i].mul(&x[i], self.bits), self.bits);
            let r = upper(&r);
            if r > residual {
                residual = r;
            }
        }
        let error = residual.div(&delta);
        Some(x.iter().map(|y| y.add_error(&error)).collect())
    }
}

/// `ε` and the bound on the steps from the trailing arguments.
fn tolerance(ev: &mut Eval<'_>, rest: &[Expr]) -> Result<Option<(Ball, u64)>, EvalError> {
    let (eps, steps) = match rest {
        [eps] => (eps, STEPS),
        [eps, k] => match k.as_i64().and_then(|k| u64::try_from(k).ok()) {
            Some(k) if (1..=MAX_STEPS).contains(&k) => (eps, k),
            _ => return Ok(None),
        },
        _ => return Ok(None),
    };
    let Some(eps) = ball(ev, eps, 64)? else {
        return Ok(None);
    };
    Ok((eps.sign() == Some(Ordering::Greater)).then_some((eps, steps)))
}

/// Iterates `step` from 0 until two iterates are within `ε`, then
/// certifies the last one.
fn solve(
    ev: &mut Eval<'_>,
    e: &Expr,
    system: &System,
    eps: &Ball,
    steps: u64,
    step: &dyn Fn(&System, &mut Vec<Ball>),
) -> Native {
    let n = system.b.len();
    let mut x = vec![Ball::from_integer(&cas_num::Integer::ZERO); n];
    let mut settled = false;
    for _ in 0..steps {
        let old = x.clone();
        step(system, &mut x);
        for y in &mut x {
            *y = Ball::exact(y.mid().clone());
        }
        let mut change = Mag::from_integer(&cas_num::Integer::ZERO);
        for (y, z) in x.iter().zip(&old) {
            let d = upper(&y.sub(z, system.bits));
            if d > change {
                change = d;
            }
        }
        if !x.iter().all(Ball::is_finite) {
            break;
        }
        let change = Ball::new(
            change.to_float().expect("finite"),
            Mag::from_integer(&cas_num::Integer::ZERO),
        );
        if change.cmp(eps) == Some(Ordering::Less) {
            settled = true;
            break;
        }
    }
    if !settled {
        warn(ev, "the iteration does not settle within its steps", e);
        return Ok(None);
    }
    let Some(x) = system.certified(&x) else {
        warn(
            ev,
            "the matrix is not shown strictly diagonally dominant: no error bound",
            e,
        );
        return Ok(None);
    };
    Ok(Some(Expr::call(
        Symbol::LIST,
        x.into_iter().map(Expr::from).collect::<Vec<_>>(),
    )))
}

pub(in crate::native) fn jacobi_linear_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b, rest @ ..] = e.args() else {
        return Ok(None);
    };
    let Some((eps, steps)) = tolerance(ev, rest)? else {
        return Ok(None);
    };
    let Some(system) = System::new(ev, a, b)? else {
        return Ok(None);
    };
    let jacobi = |s: &System, x: &mut Vec<Ball>| {
        let old = x.clone();
        for (i, y) in x.iter_mut().enumerate() {
            *y = s.off_diagonal(i, &old).div(&s.a[i][i], s.bits);
        }
    };
    let result = solve(ev, e, &system, &eps, steps, &jacobi)?;
    if let Some(r) = &result {
        ev.step("jacobi_linsolve", e, r);
    }
    Ok(result)
}

pub(in crate::native) fn gauss_seidel_linear_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    // giac's leading relaxation factor, when the first argument is not the
    // matrix.
    let (omega, args) = match args {
        [w, rest @ ..] if rows(w).is_none() && rest.len() >= 3 => (Some(w), rest),
        _ => (None, args),
    };
    let [a, b, rest @ ..] = args else {
        return Ok(None);
    };
    let Some((eps, steps)) = tolerance(ev, rest)? else {
        return Ok(None);
    };
    let Some(system) = System::new(ev, a, b)? else {
        return Ok(None);
    };
    let omega = match omega {
        None => Ball::from_integer(&cas_num::Integer::from(1)),
        Some(w) => match ball(ev, w, system.bits)? {
            Some(w)
                if w.sign() == Some(Ordering::Greater)
                    && w.cmp(&Ball::from_integer(&cas_num::Integer::from(2)))
                        == Some(Ordering::Less) =>
            {
                w
            }
            _ => return Ok(None),
        },
    };
    let one = Ball::from_integer(&cas_num::Integer::from(1));
    let keep = one.sub(&omega, system.bits);
    let sweep = move |s: &System, x: &mut Vec<Ball>| {
        for i in 0..x.len() {
            let gs = s.off_diagonal(i, x).div(&s.a[i][i], s.bits);
            x[i] = keep.mul(&x[i], s.bits).add(&omega.mul(&gs, s.bits), s.bits);
        }
    };
    let result = solve(ev, e, &system, &eps, steps, &sweep)?;
    if let Some(r) = &result {
        ev.step("gauss_seidel_linsolve", e, r);
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
