//! `Ball[c, r]`, the language's `ball(c, r)` (CAS-DOC-5 §6.5): the ball of
//! midpoint `c` and radius `r`.

use std::cmp::Ordering;

use cas_num::ball::{Ball, Float, Mag, Round};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

/// An upper bound of `r >= 0` as a radius; `None` for a negative or
/// symbolic `r`.
fn radius(r: &Expr) -> Option<Mag> {
    if let View::Ball(b) = r.view() {
        return (b.sign() != Some(Ordering::Less)).then(|| Mag::from_float(b.mid()).add(b.rad()));
    }
    let q = r.to_rational()?;
    if q.is_negative() {
        return None;
    }
    let (up, _) = Float::from_rational(&q, 30, Round::Up);
    Some(Mag::from_float(&up))
}

/// `Ball[c, r]`: an exact real `c` rounded to the working precision, its
/// radius grown by `r`, which may be exact or a ball (its upper bound); a
/// ball `c` with `r` more radius. It threads over a list `c`, with a
/// scalar `r` for every element or a list `r` of the same length paired with
/// it, other shapes being [`EvalError::Shape`] (CAS-DOC-5 §6.5). Symbols, a
/// negative `r` and a complex `c` leave it unevaluated.
pub(super) fn ball(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [c, r] = e.args() else {
        return Ok(None);
    };
    if is_list(c) {
        let items: Vec<Expr> = if is_list(r) {
            if r.args().len() != c.args().len() {
                return Err(EvalError::Shape);
            }
            c.args()
                .iter()
                .zip(r.args())
                .map(|(c, r)| Expr::apply(&Symbol::BALL.into(), &[c.clone(), r.clone()]))
                .collect()
        } else {
            c.args()
                .iter()
                .map(|c| Expr::apply(&Symbol::BALL.into(), &[c.clone(), r.clone()]))
                .collect()
        };
        return Ok(Some(Expr::apply(&Symbol::LIST.into(), &items)));
    }
    if is_list(r) {
        return Err(EvalError::Shape);
    }
    let Some(r) = radius(r) else {
        return Ok(None);
    };
    let center = match c.view() {
        View::Ball(b) => b.clone(),
        _ => match c.to_rational() {
            Some(q) => Ball::from_rational(&q, ev.context().working_precision()),
            None => return Ok(None),
        },
    };
    Ok(Some(Expr::from(center.add_error(&r))))
}

#[cfg(test)]
mod tests;
