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
//! - Over random variables (CAS-624, [`random_variables`]):
//!   `RandomSample[X, n]` `n` values of the object `X` and
//!   `RandomSample[X]` one; `RandomSample[e, n]` `n` values of an
//!   expression holding objects, each occurrence drawn afresh in each
//!   (giac vecteur.cc `find_randvars`), `n` times the occurrences at most
//!   `Limits::terms`; an expression without one stays.
//!   `RandomVector[n, X]` draws `X`; `RandomVector[n, multinomial, P, V]`
//!   `n` elements of `V`, and `RandomVector[n, multinomial, P]` the counts
//!   of each index of `P` in `n` draws, as giac answers.
//!
//! [`random_variables`]: crate::native::statistics::random_variables
//!
//! Sizes are natural numbers, and the entries drawn at most
//! `Limits::terms`; past it, or with an entry that does not draw, the call
//! stays.

use cas_num::{Integer, Natural};
use foldhash::fast::FixedState;
use hashbrown::HashMap;

use super::variates::{Law, variate};
use super::{between, count, distinct, is_list, list, refuse, sample, uniform};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::statistics::random_variables::{
    Variable, draw, is_multinomial, multinomial, variable,
};
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
    /// A value of a random variable.
    Variable(Variable<'e>),
}

/// What follows the sizes, read, decimal probabilities summed at `wp`
/// bits; `None` when no form fits.
fn entry(spec: &[Expr], wp: u64) -> Option<Entry<'_>> {
    match spec {
        [] => Some(Entry::Integer(Integer::from(-99), Integer::from(99))),
        [x] if variable(x).is_some() => variable(x).map(Entry::Variable),
        [m, p, v] if is_multinomial(m) && is_list(p) && is_list(v) => {
            multinomial(p.args(), Some(v.args()), wp).map(Entry::Variable)
        }
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
            Entry::Variable(x) => draw(ev, x)?,
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
    let wp = ev.context().working_precision();
    if let [m, p] = spec
        && is_multinomial(m)
        && is_list(p)
    {
        return counts(ev, n, p.args(), wp);
    }
    let Some(entry) = entry(spec, wp) else {
        return Ok(None);
    };
    Ok(entries(ev, &entry, n)?.map(list))
}

/// `RandomVector[n, multinomial, P]`: how many of `n` draws gave each
/// index of `P`.
fn counts(ev: &mut Eval<'_>, n: u64, p: &[Expr], wp: u64) -> Native {
    let Some(x) = multinomial(p, None, wp) else {
        return Ok(None);
    };
    let mut counts = vec![0u64; p.len()];
    for _ in 0..n {
        let Some(i) = draw(ev, &x)?.and_then(|i| i.as_i64()) else {
            return Ok(None);
        };
        counts[usize::try_from(i).expect("an index")] += 1;
    }
    Ok(Some(list(
        counts
            .into_iter()
            .map(|c| Expr::from(Integer::from(c)))
            .collect(),
    )))
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
    let Some(entry) = entry(spec, ev.context().working_precision()) else {
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
    let Some(entry) = entry(spec, ev.context().working_precision()) else {
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

/// `RandomSample[l, k]`, `RandomSample[X]`, `RandomSample[e, n]`.
pub(in crate::native) fn random_sample(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [x] => match variable(x) {
            Some(x) => draw(ev, &x),
            None => Ok(None),
        },
        [l, k] if is_list(l) => sample(ev, e, k, l.args()),
        [x, n] => values(ev, e, x, n),
        _ => Ok(None),
    }
}

/// The deepest expression whose random variables `RandomSample` finds.
const DEPTH: u16 = 1000;

/// How many objects `e` holds, counted per occurrence (saturating), once
/// per distinct subexpression.
fn occurrences(e: &Expr, seen: &mut HashMap<Expr, u64, FixedState>) -> u64 {
    if variable(e).is_some() {
        return 1;
    }
    if e.is_atom() {
        return 0;
    }
    if let Some(&n) = seen.get(e) {
        return n;
    }
    let n = e
        .args()
        .iter()
        .fold(0u64, |n, a| n.saturating_add(occurrences(a, seen)));
    seen.insert(e.clone(), n);
    n
}

/// `e` with each occurrence of an object replaced by a draw of its own;
/// `None` when one does not draw.
fn filled(
    ev: &mut Eval<'_>,
    e: &Expr,
    seen: &HashMap<Expr, u64, FixedState>,
) -> Result<Option<Expr>, EvalError> {
    if let Some(x) = variable(e) {
        return draw(ev, &x);
    }
    let Some(head) = e.head() else {
        return Ok(Some(e.clone()));
    };
    if seen.get(e).is_none_or(|&n| n == 0) {
        return Ok(Some(e.clone()));
    }
    let mut args = Vec::with_capacity(e.args().len());
    for a in e.args() {
        let Some(a) = filled(ev, a, seen)? else {
            return Ok(None);
        };
        args.push(a);
    }
    Ok(Some(Expr::apply(head, &args)))
}

/// `RandomSample[X, n]` and `RandomSample[e, n]`.
fn values(ev: &mut Eval<'_>, e: &Expr, x: &Expr, n: &Expr) -> Native {
    if x.depth() > DEPTH {
        return Ok(None);
    }
    let mut seen = HashMap::default();
    let each = occurrences(x, &mut seen);
    if each == 0 {
        return Ok(None);
    }
    let n = match size(ev, e, n) {
        Ok(Some(n)) => n,
        Ok(None) => return Ok(None),
        Err(refused) => return refused,
    };
    let most = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    if n.checked_mul(each).is_none_or(|all| all > most) {
        return refuse(ev, e, "too many draws");
    }
    let mut out = Vec::new();
    for _ in 0..n {
        let Some(drawn) = filled(ev, x, &seen)? else {
            return Ok(None);
        };
        out.push(ev.eval(&drawn)?);
    }
    Ok(Some(list(out)))
}
