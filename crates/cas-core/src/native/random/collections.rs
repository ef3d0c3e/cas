//! Random vectors, matrices and polynomials, permutations and samples
//! (CAS-DOC-33 §5, §6, CAS-373): Xcas's `randvector` (`ranv`),
//! `randmatrix` (`ranm`), `randpoly`, `randperm` and `shuffle`, and
//! `sample`.
//!
//! - `RandomVector[n, …]`, `RandomMatrix[n, p, …]` (`n × n` for
//!   `RandomMatrix[n]`) and `RandomPolynomial[n, x, …]`, `x^n + Σ_{k<n}
//!   c_k x^k` with giac's leading coefficient 1, in `x` by default and
//!   with the size and the variable in either order, as giac takes them.
//!   What follows the sizes says how each entry is drawn: nothing, an
//!   integer in `[−99, 99]`, giac's default; an integer `m > 0`, one in
//!   `[0, m)`; a range `a .. b`, a real in `[a, b)`; the name of a law's
//!   density and its parameters (`normald, 0, 1`), a variate.
//! - `RandomPermutation[n]` (`randperm`, `shuffle`), a permutation of `0 …
//!   n − 1`, as lists are indexed from 0; `RandomPermutation[l]`, the
//!   elements of `l` in a random order: Fisher–Yates.
//! - `RandomSample[l, k]` (`sample`): `k` elements of `l` at distinct
//!   places, in a random order, `Random[k, l]`.
//!
//! Sizes are natural numbers, and the entries drawn at most
//! `Limits::terms`; past it, or with an entry that does not draw, the call
//! stays.

use cas_num::{Integer, Natural};

use super::variates::{Law, variate};
use super::{between, count, distinct, is_list, list, refuse, sample, uniform};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// How each entry is drawn.
enum Entry<'e> {
    /// An integer between the bounds, both included.
    Integer(Integer, Integer),
    /// A real between the ends.
    Real(&'e Expr, &'e Expr),
    /// A variate of a law with its parameters.
    Variate(Law, &'e [Expr]),
}

/// What follows the sizes, read; `None` when no form fits.
fn entry(spec: &[Expr]) -> Option<Entry<'_>> {
    match spec {
        [] => Some(Entry::Integer(Integer::from(-99), Integer::from(99))),
        [range] if range.head() == Some(&Expr::from(Symbol::RANGE)) => match range.args() {
            [a, b] => Some(Entry::Real(a, b)),
            _ => None,
        },
        [m] if m.is_number() => {
            let m = m.to_integer().filter(Integer::is_positive)?;
            Some(Entry::Integer(Integer::ZERO, &m - &Integer::ONE))
        }
        [law, parameters @ ..] => {
            let law = Law::of_density(law.as_symbol()?)?;
            Some(Entry::Variate(law, parameters))
        }
    }
}

/// `count` entries drawn as `entry` says; `None` when one does not draw.
fn entries(
    ev: &mut Eval<'_>,
    entry: &Entry<'_>,
    count: u64,
) -> Result<Option<Vec<Expr>>, EvalError> {
    let mut out = Vec::new();
    for _ in 0..count {
        let drawn = match entry {
            Entry::Integer(low, high) => Some(Expr::from(between(ev, low, high))),
            Entry::Real(a, b) => uniform(ev, a, b)?,
            Entry::Variate(law, parameters) => variate(ev, *law, parameters)?,
        };
        let Some(drawn) = drawn else {
            return Ok(None);
        };
        out.push(drawn);
    }
    Ok(Some(out))
}

/// A size: a natural number within `Limits::terms`; `Err` after the
/// call's warning when it is a number that does not fit, `Ok(None)` for a
/// symbol.
fn size(ev: &mut Eval<'_>, e: &Expr, n: &Expr) -> Result<Option<u64>, Native> {
    match count(ev, n) {
        Ok(found) => Ok(found),
        Err(why) => Err(refuse(ev, e, why)),
    }
}

/// `RandomVector[n, …]`.
pub(in crate::native) fn random_vector(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((n, spec)) = e.args().split_first() else {
        return Ok(None);
    };
    let n = match size(ev, e, n) {
        Ok(Some(n)) => n,
        Ok(None) => return Ok(None),
        Err(refused) => return refused,
    };
    let Some(entry) = entry(spec) else {
        return Ok(None);
    };
    Ok(entries(ev, &entry, n)?.map(list))
}

/// `RandomMatrix[n]`, `RandomMatrix[n, p, …]`.
pub(in crate::native) fn random_matrix(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (n, p, spec) = match e.args() {
        [n] => (n, n, &[][..]),
        [n, p, spec @ ..] => (n, p, spec),
        [] => return Ok(None),
    };
    let (n, p) = match (size(ev, e, n), size(ev, e, p)) {
        (Ok(Some(n)), Ok(Some(p))) => (n, p),
        (Err(refused), _) | (_, Err(refused)) => return refused,
        _ => return Ok(None),
    };
    let most = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    if n.checked_mul(p).is_none_or(|all| all > most) {
        return refuse(ev, e, "too many draws");
    }
    let Some(entry) = entry(spec) else {
        return Ok(None);
    };
    let mut rows = Vec::new();
    for _ in 0..n {
        let Some(row) = entries(ev, &entry, p)? else {
            return Ok(None);
        };
        rows.push(list(row));
    }
    Ok(Some(list(rows)))
}

/// `RandomPolynomial[n]`, `RandomPolynomial[n, x, …]`,
/// `RandomPolynomial[x, n, …]`.
pub(in crate::native) fn random_polynomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let x = Expr::from(Symbol::new("x"));
    let (n, variable, spec) = match e.args() {
        [n] => (n, &x, &[][..]),
        [first, second, spec @ ..] if first.as_symbol().is_some() && second.is_number() => {
            (second, first, spec)
        }
        [n, variable, spec @ ..] => (n, variable, spec),
        [] => return Ok(None),
    };
    if variable.as_symbol().is_none() {
        return Ok(None);
    }
    let n = match size(ev, e, n) {
        Ok(Some(n)) => n,
        Ok(None) => return Ok(None),
        Err(refused) => return refused,
    };
    let Some(entry) = entry(spec) else {
        return Ok(None);
    };
    // The coefficients of x^(n−1) down to x^0, as giac draws them.
    let Some(coefficients) = entries(ev, &entry, n)? else {
        return Ok(None);
    };
    let mut terms = vec![variable.clone().pow(Expr::from(Integer::from(n)))];
    for (c, k) in coefficients.into_iter().zip((0..n).rev()) {
        terms.push(c * variable.clone().pow(Expr::from(Integer::from(k))));
    }
    ev.eval(&Expr::call(Symbol::PLUS, terms)).map(Some)
}

/// `RandomPermutation[n]`, `RandomPermutation[l]`.
pub(in crate::native) fn random_permutation(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if is_list(arg) {
        let items = arg.args();
        let n = items.len() as u64;
        let order = distinct(ev, n, &Natural::from(n));
        let shuffled = order.iter().map(|i| {
            let i = usize::try_from(i).expect("below the length");
            items[i].clone()
        });
        return Ok(Some(list(shuffled.collect())));
    }
    let n = match size(ev, e, arg) {
        Ok(Some(n)) => n,
        Ok(None) => return Ok(None),
        Err(refused) => return refused,
    };
    let order = distinct(ev, n, &Natural::from(n));
    Ok(Some(list(
        order
            .into_iter()
            .map(|i| Expr::from(Integer::from(i)))
            .collect(),
    )))
}

/// `RandomSample[l, k]`.
pub(in crate::native) fn random_sample(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, k] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return Ok(None);
    }
    sample(ev, e, k, l.args())
}
