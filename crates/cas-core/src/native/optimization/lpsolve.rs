//! Linear programs (CAS-DOC-38 §4, CAS-417): `LPSolve`, Xcas's `lpsolve`.
//!
//! `LPSolve[c, constraints, bounds…, options…]`:
//!
//! - `c` and the constraints' sides are linear in the variables with
//!   rational coefficients; the constraints are `≤`, `≥`, `=`, and `<`, `>`
//!   taken as their closures, in a list; the variables are the symbols of
//!   `c` and the constraints in the order they first appear. Anything else
//!   leaves the call.
//! - The bounds are `x = a..b`, ends rational or infinite. The options are
//!   `maximize` (or `lp_maximize`) and giac's `assume = lp_nonnegative`,
//!   every variable at least 0. A variable is free otherwise.
//! - Each variable becomes nonnegative: `x = a + x′` for a lower end `a`,
//!   with `x′ ≤ b − a` for an upper one, `x = b − x′` for an upper end
//!   alone, `x = x⁺ − x⁻` for none. The program goes to [`simplex`], which
//!   minimizes, `−c` for a maximum.
//! - The answer is giac's `{optimum, {x = …, …}}`. An unbounded program is
//!   `{-Infinity, {}}`, or `{Infinity, {}}` for a maximum, and an
//!   infeasible one `{}`, each with a warning, as giac's.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of linear programs: c·x, rows a·x ≤ b, sides l and r"
)]

use cas_num::Rational;

use super::simplex::{Outcome, Relation, Row, minimize};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// The symbols of `e` that are not builtins, in the order they first
/// appear, depth first and left to right, after those in `found`.
fn symbols_in_order(e: &Expr, found: &mut Vec<Expr>) {
    if let Some(s) = e.as_symbol() {
        if !s.is_builtin() && !found.contains(e) {
            found.push(e.clone());
        }
        return;
    }
    for a in e.args() {
        symbols_in_order(a, found);
    }
}

/// `e` as `Σ a_j x_j + k`, rational coefficients, evaluated; `None` when
/// it is not linear so.
fn linear(
    ev: &mut Eval<'_>,
    e: &Expr,
    xs: &[Expr],
) -> Result<Option<(Vec<Rational>, Rational)>, EvalError> {
    let mut coefficients = Vec::with_capacity(xs.len());
    for x in xs {
        let a = ev.eval(&Expr::call(Symbol::D, [e.clone(), x.clone()]))?;
        let Some(a) = a.to_rational() else {
            return Ok(None);
        };
        coefficients.push(a);
    }
    let zero: Vec<(Expr, Expr)> = xs.iter().map(|x| (x.clone(), Expr::from(0))).collect();
    let Some(k) = ev.eval(&e.substitute(&zero))?.to_rational() else {
        return Ok(None);
    };
    let mut rebuilt = Expr::from(k.clone());
    for (a, x) in coefficients.iter().zip(xs) {
        rebuilt += Expr::from(a.clone()) * x.clone();
    }
    let gap = ev.eval(&Expr::call(Symbol::EXPAND, [e.clone() - rebuilt]))?;
    Ok((gap == Expr::from(0)).then_some((coefficients, k)))
}

/// A rational end of a range, `None` for an infinite one; `Err(())` for
/// anything else.
fn end(e: &Expr) -> Result<Option<Rational>, ()> {
    if let Some(q) = e.to_rational() {
        return Ok(Some(q));
    }
    let infinite = e.may_contain(Symbol::INFINITY) || e.may_contain(Symbol::DIRECTED_INFINITY);
    if infinite { Ok(None) } else { Err(()) }
}

/// A program: the objective's coefficients and constant, the rows, the
/// variables' bounds, and whether it maximizes.
pub(super) struct Program {
    pub(super) xs: Vec<Expr>,
    pub(super) objective: Vec<Rational>,
    pub(super) constant: Rational,
    pub(super) rows: Vec<Row>,
    pub(super) bounds: Vec<(Option<Rational>, Option<Rational>)>,
    pub(super) maximum: bool,
    pub(super) integer: Vec<bool>,
}

/// The program of a call; `None` when it is not one.
pub(super) fn program(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<Program>, EvalError> {
    let [objective, rest @ ..] = e.args() else {
        return Ok(None);
    };
    let (constraints, rest) = match rest {
        [c, more @ ..] if head(c) == Some(Symbol::LIST) => (c.args(), more),
        _ => (&[][..], rest),
    };
    let mut maximum = false;
    let mut nonnegative = false;
    let mut kind: Option<Expr> = None;
    let mut ranges: Vec<(&Expr, &Expr, &Expr)> = Vec::new();
    for option in rest {
        if option.as_symbol() == Some(Symbol::MAXIMIZE) || *option == Expr::symbol("lp_maximize") {
            maximum = true;
            continue;
        }
        let (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) = (head(option), option.args()) else {
            return Ok(None);
        };
        if l.as_symbol() == Some(Symbol::ASSUME) {
            if *r == Expr::symbol("lp_nonnegative") {
                nonnegative = true;
            } else {
                kind = Some(r.clone());
            }
            continue;
        }
        let (Some(Symbol::RANGE), [a, b]) = (head(r), r.args()) else {
            return Ok(None);
        };
        ranges.push((l, a, b));
    }
    let mut xs = Vec::new();
    symbols_in_order(objective, &mut xs);
    for c in constraints {
        symbols_in_order(c, &mut xs);
    }
    for (x, _, _) in &ranges {
        if x.as_symbol().is_none_or(Symbol::is_builtin) {
            return Ok(None);
        }
        symbols_in_order(x, &mut xs);
    }
    let Some((objective, constant)) = linear(ev, objective, &xs)? else {
        return Ok(None);
    };
    let mut rows = Vec::new();
    for c in constraints {
        let (Some(h), [l, r]) = (head(c), c.args()) else {
            return Ok(None);
        };
        let relation = match h {
            Symbol::LESS_EQUAL | Symbol::LESS => Relation::AtMost,
            Symbol::GREATER_EQUAL | Symbol::GREATER => Relation::AtLeast,
            Symbol::EQUATION | Symbol::EQUAL => Relation::Equal,
            _ => return Ok(None),
        };
        let Some((a, k)) = linear(ev, &(l.clone() - r.clone()), &xs)? else {
            return Ok(None);
        };
        rows.push((a, relation, -k));
    }
    let start = if nonnegative {
        Some(Rational::ZERO)
    } else {
        None
    };
    let mut bounds = vec![(start, None); xs.len()];
    for (x, a, b) in ranges {
        let j = xs.iter().position(|v| v == x).expect("collected");
        let (Ok(low), Ok(high)) = (end(a), end(b)) else {
            return Ok(None);
        };
        bounds[j] = (low, high);
    }
    let integer = match kind {
        None => vec![false; xs.len()],
        Some(k) => match super::integer::kind(&k, &mut bounds) {
            Some(flags) => flags,
            None => return Ok(None),
        },
    };
    Ok(Some(Program {
        xs,
        objective,
        constant,
        rows,
        bounds,
        maximum,
        integer,
    }))
}

/// How a variable goes to nonnegative ones: `x = a + x′`, `x = b − x′`, or
/// `x = x⁺ − x⁻`, by the column of `x′` or `x⁺`.
#[derive(Clone, Copy)]
enum Shift {
    From(usize),
    Below(usize),
    Free(usize),
}

/// The continuous relaxation's best point and value, in the program's
/// variables.
pub(super) fn relaxation(p: &Program, bounds: &[(Option<Rational>, Option<Rational>)]) -> Outcome {
    let mut shifts = Vec::with_capacity(p.xs.len());
    let mut columns = 0;
    for (low, high) in bounds {
        shifts.push(match (low, high) {
            (Some(_), _) => Shift::From(columns),
            (None, Some(_)) => Shift::Below(columns),
            (None, None) => Shift::Free(columns),
        });
        columns += if matches!((low, high), (None, None)) {
            2
        } else {
            1
        };
    }
    // a·x in the new columns, and the constant it adds.
    let transform = |a: &[Rational]| -> (Vec<Rational>, Rational) {
        let mut row = vec![Rational::ZERO; columns];
        let mut constant = Rational::ZERO;
        for ((coefficient, shift), (low, high)) in a.iter().zip(&shifts).zip(bounds) {
            match *shift {
                Shift::From(c) => {
                    row[c] = coefficient.clone();
                    constant = &constant + &(coefficient * low.as_ref().expect("a lower end"));
                }
                Shift::Below(c) => {
                    row[c] = -coefficient.clone();
                    constant = &constant + &(coefficient * high.as_ref().expect("an upper end"));
                }
                Shift::Free(c) => {
                    row[c] = coefficient.clone();
                    row[c + 1] = -coefficient.clone();
                }
            }
        }
        (row, constant)
    };
    let sign = if p.maximum {
        -Rational::ONE
    } else {
        Rational::ONE
    };
    let objective: Vec<Rational> = p.objective.iter().map(|a| a * &sign).collect();
    let (cost, _) = transform(&objective);
    let mut rows: Vec<Row> = Vec::new();
    for (a, relation, b) in &p.rows {
        let (row, constant) = transform(a);
        rows.push((row, *relation, b - &constant));
    }
    for (shift, (low, high)) in shifts.iter().zip(bounds) {
        if let (Shift::From(c), Some(low), Some(high)) = (*shift, low, high) {
            let mut row = vec![Rational::ZERO; columns];
            row[c] = Rational::ONE;
            rows.push((row, Relation::AtMost, high - low));
        }
    }
    match minimize(&cost, &rows) {
        Outcome::Optimal(y, _) => {
            let x: Vec<Rational> = shifts
                .iter()
                .zip(bounds)
                .map(|(shift, (low, high))| match *shift {
                    Shift::From(c) => low.as_ref().expect("a lower end") + &y[c],
                    Shift::Below(c) => high.as_ref().expect("an upper end") - &y[c],
                    Shift::Free(c) => &y[c] - &y[c + 1],
                })
                .collect();
            let value = p
                .objective
                .iter()
                .zip(&x)
                .fold(p.constant.clone(), |s, (a, v)| &s + &(a * v));
            Outcome::Optimal(x, value)
        }
        other => other,
    }
}

/// `LPSolve[…]`.
pub(in crate::native) fn lp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some(p) = program(ev, e)? else {
        return Ok(None);
    };
    let outcome = if p.integer.iter().any(|&b| b) {
        super::integer::branch_and_bound(&p)
    } else {
        relaxation(&p, &p.bounds)
    };
    let result = match outcome {
        Outcome::Optimal(x, value) => list(vec![
            Expr::from(value),
            list(
                p.xs.iter()
                    .zip(x)
                    .map(|(v, q)| Expr::call(Symbol::EQUATION, [v.clone(), Expr::from(q)]))
                    .collect(),
            ),
        ]),
        Outcome::Unbounded => {
            warn(ev, e, "the program is unbounded");
            let infinity = Expr::from(Symbol::INFINITY);
            let value = if p.maximum { infinity } else { -infinity };
            list(vec![value, list(Vec::new())])
        }
        Outcome::Infeasible => {
            warn(ev, e, "the program has no feasible point");
            list(Vec::new())
        }
        Outcome::GaveUp => {
            warn(ev, e, "the program takes too many steps");
            return Ok(None);
        }
    };
    ev.step("lpsolve", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
