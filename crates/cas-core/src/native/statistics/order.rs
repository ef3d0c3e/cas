//! Order statistics (CAS-DOC-22 §4, CAS-295): a sample sorted by the signs
//! of the differences of its values, as `Max` and `Min` compare them, and
//! the values at which its cumulative frequency reaches a part of the
//! total. A pair of values that the comparisons do not tell apart, and that
//! are not one expression, leaves the call, and so does a frequency that is
//! not a non-negative real number.

use std::cmp::Ordering;

use super::super::rounding::real_sign;
use super::Sample;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;

type Native = Result<Option<Expr>, EvalError>;

/// Whether `p` is a real number in `[0, 1]`.
pub(super) fn is_fraction(ev: &mut Eval<'_>, p: &Expr) -> Result<bool, EvalError> {
    let rest = ev.eval(&(1 - p))?;
    let ctx = ev.context();
    let at_least_0 =
        |x: &Expr| matches!(real_sign(ctx, x), Some(Ordering::Greater | Ordering::Equal));
    Ok(at_least_0(p) && at_least_0(&rest))
}

/// How two values compare, by the sign of their difference; `None` when it
/// is not decided.
fn compare(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<Ordering>, EvalError> {
    if a == b {
        return Ok(Some(Ordering::Equal));
    }
    if let (Some(p), Some(q)) = (a.to_rational(), b.to_rational()) {
        return Ok(Some(p.cmp(&q)));
    }
    let difference = ev.eval(&(a - b))?;
    Ok(real_sign(ev.context(), &difference))
}

/// The pairs of a value and its frequency sorted by value, stably: a merge
/// sort, whose comparisons may fail; `None` when one does.
fn merge_sort(
    ev: &mut Eval<'_>,
    mut items: Vec<(Expr, Expr)>,
) -> Result<Option<Vec<(Expr, Expr)>>, EvalError> {
    if items.len() < 2 {
        return Ok(Some(items));
    }
    let right = items.split_off(items.len() / 2);
    let Some(left) = merge_sort(ev, items)? else {
        return Ok(None);
    };
    let Some(right) = merge_sort(ev, right)? else {
        return Ok(None);
    };
    let mut out = Vec::with_capacity(left.len() + right.len());
    let mut right = right.into_iter().peekable();
    for item in left {
        // The right's smaller values first; equal ones after the left's.
        while let Some((x, _)) = right.peek() {
            match compare(ev, x, &item.0)? {
                Some(Ordering::Less) => out.push(right.next().expect("peeked")),
                Some(_) => break,
                None => return Ok(None),
            }
        }
        out.push(item);
    }
    out.extend(right);
    Ok(Some(out))
}

/// A sample's values of positive frequency, ascending, with their
/// cumulative frequencies, the last of which is the total.
pub(super) struct Sorted {
    values: Vec<Expr>,
    cumulative: Vec<Expr>,
}

impl Sorted {
    /// `s` sorted; `None` when a frequency is not a non-negative real
    /// number, all are 0, or two values are not told apart.
    pub(super) fn of(ev: &mut Eval<'_>, s: &Sample) -> Result<Option<Sorted>, EvalError> {
        let mut items = Vec::with_capacity(s.values.len());
        for (i, x) in s.values.iter().enumerate() {
            let f = s
                .frequencies
                .as_ref()
                .map_or_else(|| Expr::from(1), |f| f[i].clone());
            match real_sign(ev.context(), &f) {
                Some(Ordering::Greater) => items.push((x.clone(), f)),
                Some(Ordering::Equal) => {}
                _ => return Ok(None),
            }
        }
        if items.is_empty() {
            return Ok(None);
        }
        let Some(items) = merge_sort(ev, items)? else {
            return Ok(None);
        };
        let mut values = Vec::with_capacity(items.len());
        let mut cumulative = Vec::with_capacity(items.len());
        let mut sum = Expr::from(0);
        for (x, f) in items {
            sum = ev.eval(&(sum + f))?;
            cumulative.push(sum.clone());
            values.push(x);
        }
        Ok(Some(Sorted { values, cumulative }))
    }

    /// The least and the greatest value.
    pub(super) fn ends(&self) -> (Expr, Expr) {
        let last = self.values.len() - 1;
        (self.values[0].clone(), self.values[last].clone())
    }

    fn total(&self) -> &Expr {
        &self.cumulative[self.cumulative.len() - 1]
    }

    /// The first value whose cumulative frequency reaches `target`, by its
    /// index, and how it compares to it; `None` when a comparison is not
    /// decided.
    fn reach(
        &self,
        ev: &mut Eval<'_>,
        target: &Expr,
    ) -> Result<Option<(usize, Ordering)>, EvalError> {
        for (k, c) in self.cumulative.iter().enumerate() {
            let difference = ev.eval(&(c - target))?;
            match real_sign(ev.context(), &difference) {
                Some(Ordering::Less) => {}
                Some(o) => return Ok(Some((k, o))),
                None => return Ok(None),
            }
        }
        Ok(None)
    }

    /// The value at which the cumulative frequency passes half the total,
    /// or the half-sum of the two around it when it reaches half exactly.
    pub(super) fn median(&self, ev: &mut Eval<'_>) -> Native {
        let half = ev.eval(&(self.total() / 2))?;
        let Some((k, o)) = self.reach(ev, &half)? else {
            return Ok(None);
        };
        if o == Ordering::Equal && k + 1 < self.values.len() {
            let sum = &self.values[k] + &self.values[k + 1];
            return Ok(Some(ev.eval(&(sum / 2))?));
        }
        Ok(Some(self.values[k].clone()))
    }

    /// The smallest value whose cumulative frequency is at least `p` times
    /// the total, for `0 ≤ p ≤ 1`.
    pub(super) fn quantile(&self, ev: &mut Eval<'_>, p: &Expr) -> Native {
        let target = ev.eval(&(p * self.total()))?;
        Ok(self
            .reach(ev, &target)?
            .map(|(k, _)| self.values[k].clone()))
    }
}
