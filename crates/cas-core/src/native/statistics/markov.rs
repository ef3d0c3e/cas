//! Markov chains (CAS-628): Xcas's `markov` and `randmarkov` (giac
//! misc.cc `is_stochastic` :7247, `proba2adjacence` :7082, `strongconnect`
//! :7108 and `tarjan` :7151, `classify_scc` :7163, `_markov` :7279,
//! `_randmarkov` :7411; vecteur.cc `dichotomy` :11196).
//!
//! - **The chain.** A square matrix of real numbers, exact rationals or
//!   decimals (anything else stays), none negative, whose rows sum to 1, or
//!   else whose columns do (then it is transposed), as giac's
//!   `is_stochastic`: exactly for a row of rationals, within giac's
//!   `epsilon` 10⁻¹² in `f64` for a row holding a decimal. Otherwise the
//!   call stays with giac's "Not a stochastic matrix!".
//! - **The classes.** State `i` leads to `j` when the entry `(i, j)` is not
//!   0. Tarjan's search from state 0, successors in ascending order, gives
//!   the strongly connected components in giac's order, each listing its
//!   states as they leave the stack (`[[2,1,0]]`); a class no edge leaves
//!   is recurrent, the others transient (`[]` when there are none).
//! - **markov** answers `Sequence[recurrent classes, invariant vectors,
//!   transient classes, absorption]`. The invariant vector of a recurrent
//!   class `C` solves `π (P_C − I) = 0, Σ π = 1` (one equation replaced by
//!   the sum), zero-padded to every state. The absorption matrix has one
//!   row per state and one column per recurrent class: a recurrent state's
//!   row is its class's unit vector, a transient state's row is
//!   `B = (I − Q)⁻¹ R` over the transient states (`Q` among them, `R`
//!   from them into each class); one recurrent class gives all ones (giac
//!   :7349). giac iterates the chain where it finds no kernel and so prints
//!   floats (`1.0`, `3812798742493/7625597484987`); we solve, exactly on
//!   exact input (R67), in balls when an entry is a decimal; and we answer
//!   where giac's iteration never ends (a transient state feeding a
//!   periodic class). `markov(M, eps)` accepts any real number `eps` and
//!   leaves it unused; another `eps` stays with giac's Bad Argument Value.
//! - **randmarkov** (`RandomMarkovChain`): `randmarkov(M, i0, n)` the path
//!   of `n + 1` states from `i0`, each drawn from its row; `randmarkov(M,
//!   [i…], n)` one path row per start; `randmarkov(v, nt)` a random
//!   stochastic matrix in giac's block layout, recurrent classes of the
//!   sizes in `v` (or of the one size `v`) then `nt` transient states, the
//!   entries of each row uniform draws in `[0, 1)` divided by their sum,
//!   decimals. A row of rationals draws exactly (integer masses below their
//!   common denominator), a row holding a decimal in `f64` by cumulative
//!   search, clamped to its last state of positive mass. giac's errors stay
//!   with its message: Bad Argument Value, Invalid dimension (a start
//!   outside the states, a size not positive), "Not a stochastic matrix!";
//!   a fourth (plot) argument stays without one.
//! - **Limits** (`Limits::terms`): `markov` takes at most `terms` entries,
//!   a path at most `terms` states in all, `randmarkov(v, nt)` a matrix of
//!   at most `terms` entries; past them the call stays with "too many
//!   states" or "too many draws", before any work.

#![allow(
    clippy::many_single_char_names,
    reason = "p, e, t, b, r, n, k: the chain, the call, the transient states, B, R and sizes"
)]

use cas_num::ball::{Ball, Float};
use cas_num::{Integer, Natural, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::matrix::inverse_of;
use crate::native::numeric::approximation;
use crate::native::random::unit;
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `gensizeerr`.
const BAD: &str = "Bad Argument Value";
/// giac's `gendimerr`.
const DIMENSION: &str = "Invalid dimension";
/// giac's `is_stochastic` failing.
const STOCHASTIC: &str = "Not a stochastic matrix!";
/// giac's `epsilon` (10⁻¹²): how far from 1 a decimal row may sum. giac
/// keeps fewer bits of its doubles, so its own cut-off lies within a few
/// units of the last place around it.
const EPSILON: f64 = 1e-12;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn state(i: usize) -> Expr {
    Expr::from(Integer::from(i as u64))
}

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// The rows of a square matrix, `None` for another expression.
fn square(m: &Expr) -> Option<Vec<&[Expr]>> {
    if !is_list(m) || m.args().is_empty() {
        return None;
    }
    let n = m.args().len();
    m.args()
        .iter()
        .map(|r| (is_list(r) && r.args().len() == n).then(|| r.args()))
        .collect()
}

/// An entry's value in `f64`: a rational's, a finite real ball's midpoint;
/// `None` for another expression.
fn value(x: &Expr) -> Option<f64> {
    if let Some(q) = x.to_rational() {
        return Some(Ball::from_rational(&q, 53).mid().to_f64());
    }
    match x.view() {
        View::Ball(b) if b.is_finite() => Some(b.mid().to_f64()),
        _ => None,
    }
}

/// Whether an entry is 0: an exact 0, or a ball whose midpoint is 0.
fn is_zero(x: &Expr) -> bool {
    match x.view() {
        View::Ball(b) => b.mid().is_zero(),
        _ => x.to_rational().is_some_and(|q| q.is_zero()),
    }
}

/// Whether each row sums to 1: exactly for rationals, within [`EPSILON`]
/// in `f64` for a row holding a decimal.
fn rows_sum_to_one(rows: &[Vec<Expr>]) -> bool {
    rows.iter().all(|r| {
        let exact: Option<Vec<Rational>> = r.iter().map(Expr::to_rational).collect();
        match exact {
            Some(q) => q.iter().fold(Rational::from(0), |s, q| &s + q) == Rational::from(1),
            None => (r.iter().filter_map(value).sum::<f64>() - 1.0).abs() <= EPSILON,
        }
    })
}

/// The rows of `m` as a chain's, rows summing to 1, and whether an entry
/// is a decimal; `None` when `m` is no stochastic matrix (giac's
/// `is_stochastic`).
fn stochastic(rows: &[&[Expr]]) -> Option<(Vec<Vec<Expr>>, bool)> {
    let mut decimal = false;
    for x in rows.iter().flat_map(|r| r.iter()) {
        let v = value(x)?;
        if v < 0.0 || x.to_rational().is_some_and(|q| q.is_negative()) {
            return None;
        }
        decimal |= x.to_rational().is_none();
    }
    let rows: Vec<Vec<Expr>> = rows.iter().map(|r| r.to_vec()).collect();
    if rows_sum_to_one(&rows) {
        return Some((rows, decimal));
    }
    let n = rows.len();
    let columns: Vec<Vec<Expr>> = (0..n)
        .map(|j| (0..n).map(|i| rows[i][j].clone()).collect())
        .collect();
    rows_sum_to_one(&columns).then_some((columns, decimal))
}

/// The strongly connected components of the graph `next` (each state's
/// successors, ascending), by Tarjan's search from state 0 with an explicit
/// stack, in giac's order (`strongconnect`): components as they complete,
/// each its states as they leave the stack.
fn components(next: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = next.len();
    let mut search = Search {
        index: vec![Search::NEW; n],
        low: vec![0; n],
        on_stack: vec![false; n],
        stack: Vec::new(),
        calls: Vec::new(),
        count: 0,
    };
    let mut found = Vec::new();
    for root in 0..n {
        if search.index[root] != Search::NEW {
            continue;
        }
        search.enter(root);
        while let Some(&(v, tried)) = search.calls.last() {
            if let Some(&w) = next[v].get(tried) {
                if let Some(top) = search.calls.last_mut() {
                    top.1 += 1;
                }
                if search.index[w] == Search::NEW {
                    search.enter(w);
                } else if search.on_stack[w] {
                    search.low[v] = search.low[v].min(search.index[w]);
                }
                continue;
            }
            search.calls.pop();
            if search.low[v] == search.index[v] {
                let mut component = Vec::new();
                while let Some(s) = search.stack.pop() {
                    search.on_stack[s] = false;
                    component.push(s);
                    if s == v {
                        break;
                    }
                }
                found.push(component);
            }
            if let Some(&(u, _)) = search.calls.last() {
                search.low[u] = search.low[u].min(search.low[v]);
            }
        }
    }
    found
}

/// Tarjan's search state: each state's index and low link, the component
/// stack, and the calls in progress (a state and the successors it tried).
struct Search {
    index: Vec<usize>,
    low: Vec<usize>,
    on_stack: Vec<bool>,
    stack: Vec<usize>,
    calls: Vec<(usize, usize)>,
    count: usize,
}

impl Search {
    /// The index of a state not yet visited.
    const NEW: usize = usize::MAX;

    /// Visits `v`: giac's `strongconnect` entered.
    fn enter(&mut self, v: usize) {
        self.index[v] = self.count;
        self.low[v] = self.count;
        self.count += 1;
        self.stack.push(v);
        self.on_stack[v] = true;
        self.calls.push((v, 0));
    }
}

/// The recurrent and the transient classes of `next`, each in Tarjan's
/// order (giac's `classify_scc`): a class is recurrent when no edge leaves
/// it.
fn classes(next: &[Vec<usize>]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let found = components(next);
    let mut class = vec![0; next.len()];
    for (c, states) in found.iter().enumerate() {
        for &s in states {
            class[s] = c;
        }
    }
    let mut recurrent = Vec::new();
    let mut transient = Vec::new();
    for (c, states) in found.into_iter().enumerate() {
        if states
            .iter()
            .all(|&s| next[s].iter().all(|&w| class[w] == c))
        {
            recurrent.push(states);
        } else {
            transient.push(states);
        }
    }
    (recurrent, transient)
}

/// Whether `eps` is a real number, as giac's `evalf_double` reads it.
fn real(ev: &Eval<'_>, eps: &Expr) -> bool {
    if eps.to_rational().is_some() {
        return true;
    }
    let wp = ev.context().working_precision();
    approximation(ev.context(), eps, wp)
        .is_ok_and(|v| matches!(v.view(), View::Ball(b) if b.is_finite()))
}

/// The square matrix of `entries` (row by row) as a list of lists.
fn matrix(entries: Vec<Vec<Expr>>) -> Expr {
    list(entries.into_iter().map(list).collect())
}

/// `MarkovChain[m]`, `MarkovChain[m, eps]`: Xcas's `markov`.
pub(in crate::native) fn markov(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let m = match e.args() {
        [m] => m,
        [m, eps] => {
            if !real(ev, eps) {
                return refuse(ev, e, BAD);
            }
            m
        }
        _ => return Ok(None),
    };
    let Some(rows) = square(m) else {
        return refuse(ev, e, STOCHASTIC);
    };
    let n = rows.len();
    if n.checked_mul(n)
        .is_none_or(|all| all > ev.context().limits().terms)
    {
        return refuse(ev, e, "too many states");
    }
    let Some((mut p, decimal)) = stochastic(&rows) else {
        return refuse(ev, e, STOCHASTIC);
    };
    let next: Vec<Vec<usize>> = p
        .iter()
        .map(|r| (0..n).filter(|&j| !is_zero(&r[j])).collect())
        .collect();
    let (recurrent, transient) = classes(&next);
    // Every entry a ball when one is a decimal, so the solves are numeric.
    let wp = ev.context().working_precision();
    let one = if decimal {
        for x in p.iter_mut().flatten() {
            if let Some(q) = x.to_rational() {
                *x = Expr::from(Ball::from_rational(&q, wp));
            }
        }
        Expr::from(Ball::from_rational(&Rational::from(1), wp))
    } else {
        Expr::from(1)
    };
    let mut invariant = Vec::new();
    for class in &recurrent {
        let Some(w) = invariant_of(ev, e, &p, class, &one)? else {
            return Ok(None);
        };
        invariant.push(w);
    }
    let r = recurrent.len();
    let absorption = if r == 1 {
        vec![list(vec![Expr::from(1)]); n]
    } else {
        let mut t: Vec<usize> = transient.iter().flatten().copied().collect();
        t.sort_unstable();
        let b = absorbed(ev, e, &p, &t, &recurrent)?;
        let Some(b) = b else {
            return Ok(None);
        };
        (0..n)
            .map(|s| {
                if let Some(c) = recurrent.iter().position(|c| c.contains(&s)) {
                    list((0..r).map(|j| Expr::from(i64::from(j == c))).collect())
                } else {
                    let i = t.binary_search(&s).expect("a transient state");
                    list(b[i].clone())
                }
            })
            .collect()
    };
    let states = |classes: &[Vec<usize>]| {
        list(
            classes
                .iter()
                .map(|c| list(c.iter().map(|&s| state(s)).collect()))
                .collect(),
        )
    };
    Ok(Some(Expr::call(
        Symbol::SEQUENCE,
        [
            states(&recurrent),
            list(invariant),
            states(&transient),
            list(absorption),
        ],
    )))
}

/// The invariant vector of the recurrent `class` of the chain `p` (its
/// entries balls or rationals, as `one`): `π (P_C − I) = 0`, `Σ π = 1`,
/// solved as the last column of the inverse of `(P_C − I)ᵀ` with its last
/// row the ones, zero-padded to every state; `None` after `inverse_of`'s
/// warning.
fn invariant_of(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &[Vec<Expr>],
    class: &[usize],
    one: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let mut v = class.to_vec();
    v.sort_unstable();
    let k = v.len();
    let mut a = Vec::with_capacity(k);
    for (i, &si) in v.iter().enumerate().take(k - 1) {
        let mut row = Vec::with_capacity(k);
        for (j, &sj) in v.iter().enumerate() {
            let x = p[sj][si].clone();
            row.push(if i == j {
                ev.eval(&(x - one.clone()))?
            } else {
                x
            });
        }
        a.push(row);
    }
    a.push(vec![one.clone(); k]);
    let Some(inverse) = inverse_of(ev, e, &matrix(a))? else {
        return Ok(None);
    };
    let mut w = vec![Expr::from(0); p.len()];
    for (row, &s) in inverse.args().iter().zip(&v) {
        w[s] = row.args()[k - 1].clone();
    }
    Ok(Some(list(w)))
}

/// `B = (I − Q)⁻¹ R` over the transient states `t` (ascending), one column
/// per recurrent class; `None` after `inverse_of`'s warning.
fn absorbed(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &[Vec<Expr>],
    t: &[usize],
    recurrent: &[Vec<usize>],
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    if t.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let mut a = Vec::with_capacity(t.len());
    for &i in t {
        let mut row = Vec::with_capacity(t.len());
        for &j in t {
            row.push(ev.eval(&(Expr::from(i64::from(i == j)) - p[i][j].clone()))?);
        }
        a.push(row);
    }
    let Some(inverse) = inverse_of(ev, e, &matrix(a))? else {
        return Ok(None);
    };
    // R by columns: from each transient state into each class.
    let mut into = Vec::with_capacity(recurrent.len());
    for class in recurrent {
        let mut column = Vec::with_capacity(t.len());
        for &i in t {
            let terms: Vec<Expr> = class.iter().map(|&j| p[i][j].clone()).collect();
            column.push(ev.eval(&Expr::call(Symbol::PLUS, terms))?);
        }
        into.push(column);
    }
    let mut b = Vec::with_capacity(t.len());
    for row_of_inverse in inverse.args() {
        let mut row = Vec::with_capacity(recurrent.len());
        for column in &into {
            let terms: Vec<Expr> = row_of_inverse
                .args()
                .iter()
                .zip(column)
                .map(|(x, y)| x.clone() * y.clone())
                .collect();
            row.push(ev.eval(&Expr::call(Symbol::PLUS, terms))?);
        }
        b.push(row);
    }
    Ok(Some(b))
}

/// One row of a chain, ready to draw the next state from.
enum Row {
    /// Rationals: cumulative integer masses over their common denominator.
    Exact {
        cumulative: Vec<Integer>,
        total: Natural,
    },
    /// A decimal among them: cumulative sums in `f64`, and the last state
    /// of positive mass.
    Decimal { cumulative: Vec<f64>, last: usize },
}

impl Row {
    fn new(r: &[Expr]) -> Row {
        let exact: Option<Vec<Rational>> = r.iter().map(Expr::to_rational).collect();
        if let Some(q) = exact {
            let d = q.iter().fold(Integer::ONE, |d, q| d.lcm(q.denom()));
            let mut acc = Integer::ZERO;
            let cumulative: Vec<Integer> = q
                .iter()
                .map(|q| {
                    acc = &acc + &(q.numer() * &(&d / q.denom()));
                    acc.clone()
                })
                .collect();
            let total = Natural::try_from(&acc).expect("a positive mass");
            return Row::Exact { cumulative, total };
        }
        let masses: Vec<f64> = r.iter().map(|x| value(x).unwrap_or(0.0)).collect();
        let last = masses.iter().rposition(|&m| m > 0.0).unwrap_or(0);
        let mut acc = 0.0;
        let cumulative = masses
            .iter()
            .map(|m| {
                acc += m;
                acc
            })
            .collect();
        Row::Decimal { cumulative, last }
    }

    /// The next state, drawn.
    fn draw(&self, ev: &mut Eval<'_>) -> usize {
        match self {
            Row::Exact { cumulative, total } => {
                let r = Integer::from(ev.draw(|g| g.below(total)));
                cumulative.partition_point(|c| *c <= r)
            }
            Row::Decimal { cumulative, last } => {
                let total = cumulative.last().copied().unwrap_or(0.0);
                let target = unit(ev).mid().to_f64() * total;
                cumulative.partition_point(|&c| c <= target).min(*last)
            }
        }
    }
}

/// A natural number that fits a `usize`; `None` for anything else.
fn natural(x: &Expr) -> Option<usize> {
    x.to_integer()
        .filter(|k| !k.is_negative())
        .and_then(|k| usize::try_from(&k).ok())
}

/// `RandomMarkovChain[m, i0, n]`, `RandomMarkovChain[m, {i…}, n]`,
/// `RandomMarkovChain[v, nt]`: Xcas's `randmarkov`.
pub(in crate::native) fn random_markov(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [v, nt] => random_chain(ev, e, v, nt),
        [m, starts, n] => path(ev, e, m, starts, n),
        _ => Ok(None),
    }
}

/// `randmarkov(m, i0, n)` and `randmarkov(m, [i…], n)`.
fn path(ev: &mut Eval<'_>, e: &Expr, m: &Expr, starts: &Expr, n: &Expr) -> Native {
    let one = !is_list(starts);
    let given = if one {
        std::slice::from_ref(starts)
    } else {
        starts.args()
    };
    if given.iter().any(|s| s.to_integer().is_none()) {
        return refuse(ev, e, BAD);
    }
    let Some(steps) = n
        .to_integer()
        .filter(|k| !k.is_negative())
        .map(|k| usize::try_from(&k).ok())
    else {
        return refuse(ev, e, BAD);
    };
    let Some((p, _)) = square(m).and_then(|rows| stochastic(&rows)) else {
        return refuse(ev, e, STOCHASTIC);
    };
    let size = p.len();
    let Some(from) = given
        .iter()
        .map(|s| natural(s).filter(|&i| i < size))
        .collect::<Option<Vec<usize>>>()
    else {
        return refuse(ev, e, DIMENSION);
    };
    let terms = ev.context().limits().terms;
    if steps
        .and_then(|k| k.checked_add(1))
        .and_then(|k| k.checked_mul(from.len()))
        .is_none_or(|all| all > terms)
    {
        return refuse(ev, e, "too many draws");
    }
    let steps = steps.expect("within the limit");
    let rows: Vec<Row> = p.iter().map(|r| Row::new(r)).collect();
    ev.volatile();
    let mut paths = Vec::with_capacity(from.len());
    for &start in &from {
        let mut i = start;
        let mut states = Vec::with_capacity(steps + 1);
        states.push(state(i));
        for _ in 0..steps {
            i = rows[i].draw(ev);
            states.push(state(i));
        }
        paths.push(list(states));
    }
    Ok(Some(if one {
        paths.pop().expect("one path")
    } else {
        list(paths)
    }))
}

/// `randmarkov(v, nt)`: recurrent classes of the sizes in `v`, then `nt`
/// transient states (giac misc.cc:7431-7475).
fn random_chain(ev: &mut Eval<'_>, e: &Expr, v: &Expr, nt: &Expr) -> Native {
    if nt.to_integer().is_none_or(|k| k.is_negative()) {
        return refuse(ev, e, BAD);
    }
    let sizes = if is_list(v) {
        v.args()
    } else {
        std::slice::from_ref(v)
    };
    if sizes.iter().any(|s| s.to_integer().is_none()) {
        return refuse(ev, e, BAD);
    }
    if sizes.iter().any(|s| {
        s.to_integer()
            .is_some_and(|k| k.sign() != std::cmp::Ordering::Greater)
    }) {
        return refuse(ev, e, DIMENSION);
    }
    let terms = ev.context().limits().terms;
    let total = sizes
        .iter()
        .map(natural)
        .chain([natural(nt)])
        .try_fold(0usize, |s, k| s.checked_add(k?));
    let Some(total) = total.filter(|&t| t.checked_mul(t).is_some_and(|all| all <= terms)) else {
        return refuse(ev, e, "too many states");
    };
    let recurrent: Vec<usize> = sizes.iter().filter_map(natural).collect();
    let mut rows = Vec::with_capacity(total);
    let mut start = 0;
    for &k in &recurrent {
        for _ in 0..k {
            rows.push(drawn(ev, total, start..start + k));
        }
        start += k;
    }
    while rows.len() < total {
        rows.push(drawn(ev, total, 0..total));
    }
    ev.volatile();
    Ok(Some(list(rows)))
}

/// A row of `size` decimals: uniform draws in `[0, 1)` over `block`, 0
/// elsewhere, divided by their sum (equal shares if every draw is 0).
fn drawn(ev: &mut Eval<'_>, size: usize, block: std::ops::Range<usize>) -> Expr {
    let mut row = vec![0.0; size];
    for x in &mut row[block.clone()] {
        *x = unit(ev).mid().to_f64();
    }
    let sum: f64 = row.iter().sum();
    if sum == 0.0 {
        for x in &mut row[block.clone()] {
            *x = 1.0;
        }
    }
    let sum: f64 = row.iter().sum();
    list(
        row.into_iter()
            .map(|x| {
                let f = Float::from_f64(x / sum).expect("a finite share");
                Expr::from(Ball::exact(f))
            })
            .collect(),
    )
}
