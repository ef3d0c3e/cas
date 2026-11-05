//! `FSolve` of systems (CAS-518), Xcas's `fsolve([eqns], [vars], …)`:
//! square systems of `n` equations in `n` unknowns.
//!
//! - `FSolve[{eqns}, {vars}, {x0, …}]`: Newton's iteration from the start
//!   point, the Jacobian by `D` and each step by Gaussian elimination at
//!   the working precision and 32 bits more; the limit as a list of balls,
//!   certified by Krawczyk's test, giac's `[1.93185165258,0.517638090205]`.
//!   giac's multivariate method names (`newton_solver`, `dnewton_solver`,
//!   `newtonj_solver`, `hybrid_solver`, `hybrids_solver`, `hybridj_solver`,
//!   `hybridsj_solver`) may follow; each is this iteration.
//! - Without a start point, a polynomial system's real solutions by
//!   `Solve`, approximated: a list of tuples, as giac's
//!   `fsolve([x^2+y^2=1,x=y],[x,y])`. Another system, or one `Solve`
//!   leaves, by Newton's iteration from the origin: one tuple.
//! - Krawczyk's operator over a box `X` about the limit `q`, `C` the
//!   inverse of the Jacobian at `q`, is `q − C F(q) + (I − C J(X))(X − q)`;
//!   inside `X` it proves one root there, and holds it. The box's radius
//!   starts a few bits past `q`'s and widens sixteenfold up to four times;
//!   no inclusion leaves the call with a warning.
//! - Parameters, a system that is not square, or unknowns that are not
//!   distinct symbols leave the call.

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};
use cas_poly::Cancel;

use super::super::super::approximation;
use super::super::sides;
use super::{MAX_NEWTON, list, point, rounded, settled};
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{Bounds, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// giac's multivariate root finders' names: all Newton's iteration here.
const METHODS: [&str; 7] = [
    "newton_solver",
    "dnewton_solver",
    "newtonj_solver",
    "hybrid_solver",
    "hybrids_solver",
    "hybridj_solver",
    "hybridsj_solver",
];

/// The most unknowns.
const MAX_UNKNOWNS: usize = 16;

/// A square system: its functions `F_i`, unknowns and Jacobian.
struct System {
    fs: Vec<Expr>,
    vars: Vec<Expr>,
    jacobian: Vec<Vec<Expr>>,
}

/// `FSolve[{eqns}, {vars}]` or `FSolve[{eqns}, {vars}, {x0, …}]`, with a
/// method last.
pub(in crate::native) fn system(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = match e.args() {
        [rest @ .., last]
            if rest.len() >= 2 && METHODS.iter().any(|m| *last == Expr::symbol(m)) =>
        {
            rest
        }
        args => args,
    };
    let (eqns, vars, guess) = match args {
        [eqns, vars] => (eqns, vars, None),
        [eqns, vars, guess] => (eqns, vars, Some(guess)),
        _ => return Ok(None),
    };
    let (Some(eqns), Some(vars)) = (listed(eqns), listed(vars)) else {
        return Ok(None);
    };
    let n = vars.len();
    if n == 0 || n > MAX_UNKNOWNS || eqns.len() != n || !distinct_symbols(vars) {
        return Ok(None);
    }
    let mut fs = Vec::with_capacity(n);
    for eqn in eqns {
        fs.push(match sides(eqn) {
            Some((l, r)) => ev.eval(&(l.clone() - r.clone()))?,
            None => eqn.clone(),
        });
    }
    if !numeric(ev.context(), &fs, vars) {
        return Ok(None);
    }
    let start = match guess {
        Some(g) => {
            let Some(g) = listed(g).filter(|g| g.len() == n) else {
                return Ok(None);
            };
            let mut q = Vec::with_capacity(n);
            for x in g {
                let Some(x) = point(ev.context(), x) else {
                    return Ok(None);
                };
                q.push(x);
            }
            Some(q)
        }
        None => None,
    };
    if start.is_none()
        && polynomial(ev.context(), &fs, vars)
        && let Some(found) = solved(ev, eqns, vars)?
    {
        ev.step("FSolve", e, &found);
        return Ok(Some(found));
    }
    let mut jacobian = Vec::with_capacity(n);
    for f in &fs {
        let mut row = Vec::with_capacity(n);
        for x in vars {
            row.push(ev.eval(&Expr::call(Symbol::D, [f.clone(), x.clone()]))?);
        }
        jacobian.push(row);
    }
    let system = System {
        fs,
        vars: vars.to_vec(),
        jacobian,
    };
    let bits = ev.context().working_precision();
    let q = start.unwrap_or_else(|| vec![Rational::ZERO; n]);
    let ctx = ev.context();
    let root = newton(ctx, &system, q, bits).and_then(|q| krawczyk(ctx, &system, &q, bits));
    let Some(root) = root else {
        ev.warn(Warning::Native {
            message: "Newton's iteration found no root it could certify".into(),
            expr: guess
                .cloned()
                .unwrap_or_else(|| list(vec![Expr::from(0); n])),
        });
        return Ok(None);
    };
    let result = list(root);
    ev.step("FSolve", e, &result);
    Ok(Some(result))
}

fn listed(e: &Expr) -> Option<&[Expr]> {
    (e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)).then(|| e.args())
}

fn distinct_symbols(vars: &[Expr]) -> bool {
    vars.iter()
        .enumerate()
        .all(|(i, v)| v.as_symbol().is_some() && !vars[..i].contains(v))
}

/// The substitution of `values` for `vars`.
fn at(vars: &[Expr], values: impl IntoIterator<Item = Expr>) -> Vec<(Expr, Expr)> {
    vars.iter().cloned().zip(values).collect()
}

/// Whether every `F_i` is a number at a point: no parameters.
fn numeric(ctx: &Context, fs: &[Expr], vars: &[Expr]) -> bool {
    let point = (0..vars.len()).map(|i| {
        let i = i64::try_from(i).unwrap_or(0);
        Expr::from(Rational::new(
            Integer::from(2 * i + 1),
            Integer::from(3 * i + 7),
        ))
    });
    let rules = at(vars, point);
    fs.iter()
        .all(|f| real(ctx, &f.substitute(&rules), 64).is_some())
}

/// `e`'s value as a finite real ball at `bits`.
fn real(ctx: &Context, e: &Expr, bits: u64) -> Option<Ball> {
    match approximation(ctx, e, bits).ok()?.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        View::ComplexBall(z) if z.is_finite() && z.im().contains_rational(&Rational::ZERO) => {
            Some(z.re().clone())
        }
        _ => None,
    }
}

/// Whether the `F_i` are polynomials in the unknowns alone.
fn polynomial(ctx: &Context, fs: &[Expr], vars: &[Expr]) -> bool {
    let cancel: &dyn Cancel = ctx.interrupt_flag();
    let bounds = Bounds {
        terms: ctx.limits().terms,
        bits: ctx.limits().number_bits,
    };
    Polynomials::new(fs, bounds, cancel)
        .is_ok_and(|p| p.generators().iter().all(|g| vars.contains(g)))
}

/// A polynomial system's real solutions by `Solve`, approximated; `None`
/// when `Solve` gives no list of tuples.
fn solved(ev: &mut Eval<'_>, eqns: &[Expr], vars: &[Expr]) -> Native {
    let exact = ev.eval(&Expr::call(
        Symbol::SOLVE,
        [list(eqns.to_vec()), list(vars.to_vec())],
    ))?;
    let Some(tuples) = listed(&exact) else {
        return Ok(None);
    };
    let bits = ev.context().working_precision();
    let mut found = Vec::new();
    for tuple in tuples {
        let Some(parts) = listed(tuple).filter(|p| p.len() == vars.len()) else {
            return Ok(None);
        };
        let mut values = Vec::with_capacity(parts.len());
        let mut complex = false;
        for part in parts {
            let Ok(v) = approximation(ev.context(), part, bits) else {
                return Ok(None);
            };
            match v.view() {
                View::Ball(b) if b.is_finite() => values.push(v.clone()),
                View::ComplexBall(_) => complex = true,
                _ => return Ok(None),
            }
        }
        if !complex {
            found.push(list(values));
        }
    }
    Ok(Some(list(found)))
}

/// The midpoints of the `F_i` at `q`, and of the Jacobian.
fn midpoints(
    ctx: &Context,
    system: &System,
    q: &[Rational],
    prec: u64,
) -> Option<(Vec<Rational>, Vec<Vec<Rational>>)> {
    let rules = at(&system.vars, q.iter().cloned().map(Expr::from));
    let mid = |e: &Expr| real(ctx, &e.substitute(&rules), prec).map(|b| b.mid().to_rational());
    let fv = system.fs.iter().map(mid).collect::<Option<Vec<_>>>()?;
    let jv = system
        .jacobian
        .iter()
        .map(|row| row.iter().map(mid).collect::<Option<Vec<_>>>())
        .collect::<Option<Vec<_>>>()?;
    Some((fv, jv))
}

/// Newton's iteration from `q`: its limit, uncertified.
fn newton(
    ctx: &Context,
    system: &System,
    mut q: Vec<Rational>,
    bits: u64,
) -> Option<Vec<Rational>> {
    let prec = bits.saturating_add(32);
    for _ in 0..MAX_NEWTON {
        let (fv, jv) = midpoints(ctx, system, &q, prec)?;
        if fv.iter().all(Rational::is_zero) {
            return Some(q);
        }
        let step = solve_linear(jv, fv)?;
        let mut done = true;
        for (qi, si) in q.iter_mut().zip(&step) {
            let next = rounded(&(&*qi - si), prec);
            done &= settled(si, &next, bits);
            *qi = next;
        }
        if done {
            return Some(q);
        }
    }
    None
}

/// The solution of `a s = b` by Gaussian elimination with partial
/// pivoting; `None` for a singular `a`.
fn solve_linear(mut a: Vec<Vec<Rational>>, mut b: Vec<Rational>) -> Option<Vec<Rational>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().cmp(&a[j][col].abs()))?;
        if a[pivot][col].is_zero() {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = &a[row][col] / &a[col][col];
            if factor.is_zero() {
                continue;
            }
            let pivot_row = a[col].clone();
            for (x, p) in a[row][col..].iter_mut().zip(&pivot_row[col..]) {
                *x = &*x - &(&factor * p);
            }
            let t = &factor * &b[col];
            b[row] = &b[row] - &t;
        }
    }
    let mut s = vec![Rational::ZERO; n];
    for row in (0..n).rev() {
        let mut acc = b[row].clone();
        for k in row + 1..n {
            acc = &acc - &(&a[row][k] * &s[k]);
        }
        s[row] = &acc / &a[row][row];
    }
    Some(s)
}

/// The inverse of `a`, rounded to `prec` bits; `None` for a singular one.
fn inverse(a: &[Vec<Rational>], prec: u64) -> Option<Vec<Vec<Rational>>> {
    let n = a.len();
    let mut columns = Vec::with_capacity(n);
    for j in 0..n {
        let e = (0..n)
            .map(|i| {
                if i == j {
                    Rational::from(Integer::from(1))
                } else {
                    Rational::ZERO
                }
            })
            .collect();
        columns.push(solve_linear(a.to_vec(), e)?);
    }
    Some(
        (0..n)
            .map(|i| (0..n).map(|j| rounded(&columns[j][i], prec)).collect())
            .collect(),
    )
}

/// The ball `Ball[mid, rad]` as an expression.
fn ball(mid: &Rational, rad: &Rational) -> Expr {
    Expr::call(
        Symbol::BALL,
        [Expr::from(mid.clone()), Expr::from(rad.clone())],
    )
}

/// Krawczyk's test about `q`: balls that each hold one root's part.
fn krawczyk(ctx: &Context, system: &System, q: &[Rational], bits: u64) -> Option<Vec<Expr>> {
    let prec = bits.saturating_add(32);
    let n = q.len();
    let (_, jv) = midpoints(ctx, system, q, prec)?;
    let c = inverse(&jv, prec)?;
    let exact = at(&system.vars, q.iter().cloned().map(Expr::from));
    let fq = system
        .fs
        .iter()
        .map(|f| real(ctx, &f.substitute(&exact), prec).map(Expr::from))
        .collect::<Option<Vec<_>>>()?;
    let exponent = i64::try_from(bits).unwrap_or(i64::MAX);
    let two = Rational::from(Integer::from(2));
    let scale = q
        .iter()
        .map(Rational::abs)
        .max()
        .unwrap_or(Rational::ZERO)
        .max(two.pow(-exponent));
    let mut r = &scale * &two.pow(4 - exponent);
    'widen: for _ in 0..4 {
        let boxed = at(&system.vars, q.iter().map(|qi| ball(qi, &r)));
        let jx = system
            .jacobian
            .iter()
            .map(|row| {
                row.iter()
                    .map(|d| real(ctx, &d.substitute(&boxed), prec).map(Expr::from))
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()?;
        let deviation = ball(&Rational::ZERO, &r);
        let mut found = Vec::with_capacity(n);
        for i in 0..n {
            let mut k = Expr::from(q[i].clone());
            for (cik, fk) in c[i].iter().zip(&fq) {
                k -= Expr::from(cik.clone()) * fk.clone();
            }
            for j in 0..n {
                let mut m = Expr::from(i64::from(i == j));
                for (cik, row) in c[i].iter().zip(&jx) {
                    m -= Expr::from(cik.clone()) * row[j].clone();
                }
                k += m * deviation.clone();
            }
            let Some(b) = real(ctx, &k, prec) else {
                r = &r * &Rational::from(Integer::from(16));
                continue 'widen;
            };
            let Some(rad) = b.rad().to_float().map(|f| f.to_rational()) else {
                r = &r * &Rational::from(Integer::from(16));
                continue 'widen;
            };
            if &(&b.mid().to_rational() - &q[i]).abs() + &rad >= r {
                r = &r * &Rational::from(Integer::from(16));
                continue 'widen;
            }
            found.push(Expr::from(b.round_to_accuracy(bits)));
        }
        return Some(found);
    }
    None
}
