//! Xcas's list commands (CAS-555 to CAS-560), giac's answers: positions
//! from 0, counts, differences, the small builders, maps, folds, products,
//! insertions, substitutions and heaps. giac's errors leave the call, with
//! the error as a warning.
//!
//! - `FindPositions[x, l]` (`find`, CAS-555): the positions of `x` in `l`,
//!   from 0 (misc.cc `_find`); the arguments exchanged when only the first
//!   is a list, as giac's. `FirstPosition[l, x]` (`index`, prog.cc
//!   `_index`) is the first of them, and leaves the call when there is none
//!   (giac's Bad Argument Value). Not `MemberIndex`, which counts from 1.
//!   An element is `x` when it is the same expression, or when both are
//!   numbers that `Equal` makes equal: giac finds `1` at `1.0`. The string
//!   forms stay (phase 1.6), as does `find`'s start position.
//! - `CountEqual[x, l]`, `CountLess[x, l]`, `CountGreater[x, l]`
//!   (`count_eq`, `count_inf`, `count_sup`, maple.cc `count_eq0` to
//!   `count_eq2`): the elements of `l` at every depth equal to, less than,
//!   greater than `x`, an `l` that is not a list counting itself; an order
//!   counts only where `Less` gives `True`, as giac's `is_strictly_greater`
//!   counts only what it decides (`count_inf(x,[1,2])` is 0). With `row`
//!   (`MatrixRow`, C0008) or `col` (`MatrixColumn`), the counts of each row
//!   or column of a matrix, in a list; another option is ignored.
//! - `Differences[l]` (`deltalist`, ti89.cc `_deltalist`): successive
//!   differences, row by row for a matrix (`[[1,2],[4,8]]` gives
//!   `[[3,6]]`); fewer than two elements is giac's Invalid dimension.
//! - `AccumulateHeadTail[l, p, q]` (`accumulate_head_tail`, maple.cc): the
//!   sum of the first `p` elements, the middle ones, the sum of the last
//!   `q`; `p, q ≥ 1` integers with `p + q ≤ length`.
//! - `Lengths[l]` (`sizes`, permu.cc `_sizes`): the lengths of the
//!   sublists.
//! - `MakeVector[a, …]` (`makevector`): the list of the arguments, one list
//!   itself; more than [`MAX_ELEMENTS`] arguments leave the call (C0019).
//!   `MakeSequence[l]` (`makesuite`): the sequence of `l`'s elements, a
//!   non-list `x` the list `[x]`. `NewList[n]` (`newList`): `|n|` zeros, at
//!   most [`MAX_ELEMENTS`] (giac's `LIST_SIZE_LIMIT`, larger). `Enumerate[l,
//!   s]` (`enumerate`, misc.cc): `[[s + k, l_k]]`, `s` 0 by default.
//!   `Copy[x]` (`copy`) is `x`.
//! - `ApplyEach[f, l]` (`apply`, CAS-558, prog.cc `_apply`): `Map`'s list
//!   of `f(x)`, giac's order of arguments; with `matrix` (`Matrix`, C0008)
//!   `f` of every entry of a matrix, giac's Bad Argument Type on another
//!   list. A vector's square is giac's dot product and ours stays a power
//!   (`apply(x->x^2,[[1,2],[3,4]])` is giac's `[5,25]`); several lists stay.
//! - `FoldLeft[f, init, a, …]`, `FoldRight[f, init, a, …]` (`foldl`,
//!   `foldr`, graphtheory.cc): `f(f(init, a), b)…` and `f(a, f(b, …
//!   init))`, built for the evaluator, a symbol `f` staying a call; fewer
//!   than three arguments are giac's errors.
//! - `Multiply[a, …]` (`multiply`, giac's `_prod`): `*` from the left, two
//!   vectors giving their dot product over the shorter
//!   (`multiply([1,2],[3,4,5])` is 11), vectors and matrices the products
//!   of `*` ([`linear::dot`], no conjugate); one list multiplies its
//!   elements. Dimensions that do not fit are giac's Invalid dimension.
//! - `InsertAt[l, k, x]` (`insert`, misc.cc `_insert`): `x` before the
//!   position `k` from 0, `0 ≤ k ≤ length`.
//! - `SubstituteOperands[l, k = v, …]` (`subsop`, prog.cc `subsop`): the
//!   element at `k` from 0, negative from the end, or at the nested
//!   position `[i, j, …]`, replaced by `v`, deleted for `NULL`; each rule
//!   acts on the original positions, so `subsop([0,1,2,3],0=NULL,1=NULL)`
//!   is `[2,3]`. `subsop(l, k, v)` is one rule; a rule that is not an
//!   equation with an integer position is skipped, as giac's. A first
//!   argument that is not a list stays (giac's operands of a call differ).
//! - `Heapify[l, f]`, `HeapPush[l, x, f]`, `HeapPop[l, f]` (`heapify`,
//!   `heappush`, `heappop`, CAS-560, prog.cc): libstdc++'s max-heap
//!   (`stl_heap.h`: `make_heap`, `push_heap`, `pop_heap`, `__adjust_heap`)
//!   step for step, `a` before `b` when `Less[a, b]`, or the optional
//!   `f(a, b)`, is `True` (giac's `gen_sort`); a comparison that is neither
//!   `True` nor `False` leaves the call. giac changes its variable in place
//!   and ours returns the list: `heapify` the heap (giac prints 1),
//!   `heappush` the heap with `x` pushed (giac's answer), `heappop` the heap
//!   without its top (what giac leaves in the variable; it returns the
//!   element past the end, 0 or a crash). giac's `heapify` and `heappop`
//!   never read `f`. `heappop([])` is giac's Invalid dimension.

use super::{MAX_ELEMENTS, is_list, list};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::linear::{self, Shape, shape};
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

const BAD_VALUE: &str = "bad argument value";
const BAD_TYPE: &str = "bad argument type";
const SIZE_LIMIT: &str = "calculation size limit exceeded";
const INVALID_DIMENSION: &str = "invalid dimension";

/// giac's error for `e`: the call stays, with `why` as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

fn integer(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("a length in memory"))
}

fn is_string(e: &Expr) -> bool {
    e.as_str().is_some()
}

/// Whether `ev` decides `test` is `True`.
fn holds(ev: &mut Eval<'_>, test: Symbol, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    let value = ev.eval(&Expr::call(test, [a.clone(), b.clone()]))?;
    Ok(value.as_symbol() == Some(Symbol::TRUE))
}

/// giac's `==`: the same expression, or numbers `Equal` makes equal.
fn same(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    if a == b {
        return Ok(true);
    }
    if a.is_number() && b.is_number() {
        return holds(ev, Symbol::EQUAL, a, b);
    }
    Ok(false)
}

/// The positions of `x` in the list `l`, from 0.
fn positions(ev: &mut Eval<'_>, x: &Expr, l: &Expr) -> Result<Vec<usize>, EvalError> {
    let mut out = Vec::new();
    for (i, y) in l.args().iter().enumerate() {
        if same(ev, x, y)? {
            out.push(i);
        }
    }
    Ok(out)
}

/// `x` and the list to search, giac's exchange included; `None` when
/// neither is a list.
fn searched<'e>(x: &'e Expr, l: &'e Expr) -> Option<(&'e Expr, &'e Expr)> {
    if is_list(l) {
        Some((x, l))
    } else if is_list(x) {
        Some((l, x))
    } else {
        None
    }
}

/// `FindPositions[x, l]`.
pub(in crate::native) fn find_positions(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, l] = e.args() else {
        return Ok(None);
    };
    if is_string(x) || is_string(l) {
        return Ok(None);
    }
    let Some((x, l)) = searched(x, l) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let found = positions(ev, x, l)?;
    Ok(Some(list(found.into_iter().map(integer).collect())))
}

/// `FirstPosition[l, x]`.
pub(in crate::native) fn first_position(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, x] = e.args() else {
        return refuse(ev, e, BAD_VALUE);
    };
    if is_string(x) || is_string(l) {
        return Ok(None);
    }
    let Some((x, l)) = searched(x, l) else {
        return refuse(ev, e, BAD_VALUE);
    };
    match positions(ev, x, l)?.first() {
        Some(&i) => Ok(Some(integer(i))),
        None => refuse(ev, e, BAD_VALUE),
    }
}

/// What `count_eq`, `count_inf` and `count_sup` count.
#[derive(Clone, Copy)]
enum Relation {
    Equal,
    Less,
    Greater,
}

/// The elements of `v` at every depth in `relation` to `x`; `v` itself
/// when it is not a list.
fn counted(ev: &mut Eval<'_>, relation: Relation, x: &Expr, v: &Expr) -> Result<usize, EvalError> {
    if is_list(v) {
        let mut total = 0;
        for y in v.args() {
            total += counted(ev, relation, x, y)?;
        }
        return Ok(total);
    }
    let yes = match relation {
        Relation::Equal => same(ev, v, x)?,
        Relation::Less => holds(ev, Symbol::LESS, v, x)?,
        Relation::Greater => holds(ev, Symbol::LESS, x, v)?,
    };
    Ok(usize::from(yes))
}

/// giac's `ckmatrix`: a non-empty list of lists of one non-zero length.
fn is_matrix(m: &Expr) -> bool {
    let rows = m.args();
    is_list(m)
        && rows
            .first()
            .is_some_and(|r| is_list(r) && !r.args().is_empty())
        && rows
            .iter()
            .all(|r| is_list(r) && r.args().len() == rows[0].args().len())
}

fn count(ev: &mut Eval<'_>, e: &Expr, relation: Relation) -> Native {
    let (x, v, option) = match e.args() {
        [x, v] => (x, v, None),
        [x, v, o, ..] => (x, v, o.as_symbol()),
        _ => return refuse(ev, e, BAD_VALUE),
    };
    let rows: Vec<Expr> = match option {
        Some(Symbol::MATRIX_ROW | Symbol::MATRIX_COLUMN) if !is_matrix(v) => {
            return refuse(ev, e, BAD_TYPE);
        }
        Some(Symbol::MATRIX_ROW) => v.args().to_vec(),
        Some(Symbol::MATRIX_COLUMN) => (0..v.args()[0].args().len())
            .map(|j| list(v.args().iter().map(|r| r.args()[j].clone()).collect()))
            .collect(),
        _ => return Ok(Some(integer(counted(ev, relation, x, v)?))),
    };
    let mut out = Vec::with_capacity(rows.len());
    for r in &rows {
        out.push(integer(counted(ev, relation, x, r)?));
    }
    Ok(Some(list(out)))
}

/// `CountEqual[x, l]`.
pub(in crate::native) fn count_equal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Equal)
}

/// `CountLess[x, l]`.
pub(in crate::native) fn count_less(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Less)
}

/// `CountGreater[x, l]`.
pub(in crate::native) fn count_greater(ev: &mut Eval<'_>, e: &Expr) -> Native {
    count(ev, e, Relation::Greater)
}

/// `Differences[l]`.
pub(in crate::native) fn differences(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    if l.args().len() < 2 {
        return refuse(ev, e, INVALID_DIMENSION);
    }
    let mut out = Vec::with_capacity(l.args().len() - 1);
    for pair in l.args().windows(2) {
        out.push(ev.eval(&(pair[1].clone() - pair[0].clone()))?);
    }
    Ok(Some(list(out)))
}

/// `AccumulateHeadTail[l, p, q]`.
pub(in crate::native) fn accumulate_head_tail(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, p, q] = e.args() else {
        return Ok(None);
    };
    let length = |n: &Expr| n.to_integer().and_then(|n| n.to_i64());
    let (Some(p), Some(q)) = (length(p), length(q)) else {
        return refuse(ev, e, BAD_VALUE);
    };
    let w = l.args();
    let s = i64::try_from(w.len()).expect("a length in memory");
    if !is_list(l) || p < 1 || q < 1 || p.saturating_add(q) > s {
        return refuse(ev, e, BAD_VALUE);
    }
    let (p, t) = (
        usize::try_from(p).expect("in range"),
        usize::try_from(s - q).expect("in range"),
    );
    let mut out = Vec::with_capacity(t - p + 2);
    out.push(ev.eval(&Expr::call(Symbol::PLUS, w[..p].to_vec()))?);
    out.extend_from_slice(&w[p..t]);
    out.push(ev.eval(&Expr::call(Symbol::PLUS, w[t..].to_vec()))?);
    Ok(Some(list(out)))
}

/// `Lengths[l]`.
pub(in crate::native) fn lengths(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l] = e.args() else {
        return Ok(None);
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_TYPE);
    }
    if !l.args().iter().all(is_list) {
        return refuse(ev, e, BAD_VALUE);
    }
    Ok(Some(list(
        l.args().iter().map(|r| integer(r.args().len())).collect(),
    )))
}

/// `MakeVector[a, …]`.
pub(in crate::native) fn make_vector(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) => Ok(Some(l.clone())),
        args if i64::try_from(args.len()).is_ok_and(|n| n <= MAX_ELEMENTS) => {
            Ok(Some(list(args.to_vec())))
        }
        _ => refuse(ev, e, SIZE_LIMIT),
    }
}

/// `MakeSequence[l]`.
pub(in crate::native) fn make_sequence(_: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(Some(match e.args() {
        [l] if is_list(l) => Expr::call(Symbol::SEQUENCE, l.args().to_vec()),
        [x] => list(vec![x.clone()]),
        args => Expr::call(Symbol::SEQUENCE, args.to_vec()),
    }))
}

/// `NewList[n]`.
pub(in crate::native) fn new_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n] = e.args() else {
        return refuse(ev, e, BAD_VALUE);
    };
    let Some(n) = n.to_integer() else {
        return refuse(ev, e, BAD_VALUE);
    };
    match n.to_i64().and_then(i64::checked_abs) {
        Some(n) if n <= MAX_ELEMENTS => Ok(Some(list(vec![
            Expr::from(0);
            usize::try_from(n).expect("in range")
        ]))),
        _ => refuse(ev, e, SIZE_LIMIT),
    }
}

/// `Enumerate[l, s]`.
pub(in crate::native) fn enumerate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, start) = match e.args() {
        [] => return Ok(Some(list(Vec::new()))),
        [l] => (l, Expr::from(0)),
        [l, s] => (l, s.clone()),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_TYPE);
    }
    let mut out = Vec::with_capacity(l.args().len());
    for (k, x) in l.args().iter().enumerate() {
        let index = ev.eval(&(start.clone() + integer(k)))?;
        out.push(list(vec![index, x.clone()]));
    }
    Ok(Some(list(out)))
}

/// `Copy[x]`.
pub(in crate::native) fn copy(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    Ok(Some(x.clone()))
}

/// `ApplyEach[f, l]`.
pub(in crate::native) fn apply_each(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let each = |f: &Expr, l: &Expr| {
        list(
            l.args()
                .iter()
                .map(|x| Expr::apply(f, std::slice::from_ref(x)))
                .collect(),
        )
    };
    match e.args() {
        [f, l] if is_list(l) => Ok(Some(each(f, l))),
        [f, m, option] if is_matrix_option(option) => {
            if is_list(m) && m.args().is_empty() {
                Ok(Some(m.clone()))
            } else if is_matrix(m) {
                Ok(Some(list(m.args().iter().map(|r| each(f, r)).collect())))
            } else {
                refuse(ev, e, BAD_TYPE)
            }
        }
        _ => Ok(None),
    }
}

/// Whether `option` is giac's `matrix` (`Matrix`, C0008).
fn is_matrix_option(option: &Expr) -> bool {
    option
        .as_symbol()
        .is_some_and(|s| s == Symbol::MATRIX || s.name() == "matrix")
}

/// `f` folded over `e`'s arguments `f, init, a, …` from the left or the
/// right, unevaluated.
fn fold(ev: &mut Eval<'_>, e: &Expr, left: bool) -> Native {
    let why = if e.args().len() == 1 {
        BAD_TYPE
    } else {
        BAD_VALUE
    };
    let [f, init, items @ ..] = e.args() else {
        return refuse(ev, e, why);
    };
    if items.is_empty() {
        return refuse(ev, e, why);
    }
    let folded = if left {
        items
            .iter()
            .fold(init.clone(), |acc, x| Expr::apply(f, &[acc, x.clone()]))
    } else {
        items
            .iter()
            .rev()
            .fold(init.clone(), |acc, x| Expr::apply(f, &[x.clone(), acc]))
    };
    Ok(Some(folded))
}

/// `FoldLeft[f, init, a, …]`.
pub(in crate::native) fn fold_left(ev: &mut Eval<'_>, e: &Expr) -> Native {
    fold(ev, e, true)
}

/// `FoldRight[f, init, a, …]`.
pub(in crate::native) fn fold_right(ev: &mut Eval<'_>, e: &Expr) -> Native {
    fold(ev, e, false)
}

/// giac's `a*b`, evaluated: two vectors give their dot product over the
/// shorter, vectors and matrices their product ([`linear::dot`]), anything
/// else `Times`; `None` for dimensions that do not fit.
fn product(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<Option<Expr>, EvalError> {
    let unevaluated = match (shape(a), shape(b)) {
        (Some(Shape::Vector(n)), Some(Shape::Vector(k))) => {
            let n = n.min(k);
            linear::dot(&list(a.args()[..n].to_vec()), &list(b.args()[..n].to_vec()))
        }
        (Some(_), Some(_)) => linear::dot(a, b),
        _ => Ok(Expr::call(Symbol::TIMES, [a.clone(), b.clone()])),
    };
    match unevaluated {
        Ok(p) => Ok(Some(ev.eval(&p)?)),
        Err(EvalError::Shape) => Ok(None),
        Err(err) => Err(err),
    }
}

/// `Multiply[a, …]`.
pub(in crate::native) fn multiply(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let factors = match e.args() {
        [l] if is_list(l) => l.args(),
        args => args,
    };
    let Some((first, rest)) = factors.split_first() else {
        return Ok(Some(Expr::from(1)));
    };
    let mut out = first.clone();
    for x in rest {
        match product(ev, &out, x)? {
            Some(p) => out = p,
            None => return refuse(ev, e, INVALID_DIMENSION),
        }
    }
    Ok(Some(out))
}

/// `InsertAt[l, k, x]`.
pub(in crate::native) fn insert_at(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, k, x] = e.args() else {
        return refuse(ev, e, BAD_VALUE);
    };
    if is_string(l) {
        return Ok(None);
    }
    let Some(k) = k.as_i64() else {
        return refuse(ev, e, BAD_VALUE);
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    let Some(k) = usize::try_from(k).ok().filter(|&k| k <= l.args().len()) else {
        return refuse(ev, e, INVALID_DIMENSION);
    };
    let mut out = l.args().to_vec();
    out.insert(k, x.clone());
    Ok(Some(list(out)))
}

/// Whether `v` deletes in a `subsop` rule: `NULL` (the language's user
/// symbol `$NULL`), giac's empty sequence.
fn is_null(v: &Expr) -> bool {
    v.as_symbol()
        .is_some_and(|s| s == Symbol::NULL || s.name() == "$NULL")
        || v.head().and_then(Expr::as_symbol) == Some(Symbol::SEQUENCE) && v.args().is_empty()
}

/// The two sides of an equation `a = b`.
fn sides(rule: &Expr) -> Option<(&Expr, &Expr)> {
    match (rule.head().and_then(Expr::as_symbol), rule.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [a, b]) => Some((a, b)),
        _ => None,
    }
}

/// The position `i` of a list of `len` elements, a negative one from the end.
fn position(i: &Expr, len: usize) -> Option<usize> {
    let i = i.as_i64()?;
    let len = i64::try_from(len).ok()?;
    let i = if i < 0 { i + len } else { i };
    (0..len)
        .contains(&i)
        .then(|| usize::try_from(i).expect("in range"))
}

/// The elements `items` with the substitutions `(position, value)` of
/// `rules` made at their original positions, `None` for those deleted
/// (prog.cc `subsop`); `Err(())` for giac's Invalid dimension. A rule that
/// is not an equation with an integer or a list of integers on the left
/// is skipped; a nested position inside an element that is not a list
/// leaves it.
fn substituted(items: &[Expr], rules: &[(Expr, Expr)]) -> Result<Vec<Option<Expr>>, ()> {
    let mut out: Vec<Option<Expr>> = items.iter().cloned().map(Some).collect();
    for (at, value) in rules {
        let at = if is_list(at) {
            match at.args() {
                [] => return Err(()),
                [i] => i,
                [i, rest @ ..] => {
                    let Some(i) = position(i, out.len()) else {
                        return Err(());
                    };
                    let Some(element) = &out[i] else {
                        continue;
                    };
                    if is_list(element) {
                        let inner =
                            substituted(element.args(), &[(list(rest.to_vec()), value.clone())])?;
                        out[i] = Some(list(inner.into_iter().flatten().collect()));
                    }
                    continue;
                }
            }
        } else {
            at
        };
        if at.as_i64().is_none() {
            continue;
        }
        let Some(i) = position(at, out.len()) else {
            return Err(());
        };
        out[i] = (!is_null(value)).then(|| value.clone());
    }
    Ok(out)
}

/// `SubstituteOperands[l, k = v, …]`.
pub(in crate::native) fn substitute_operands(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, rules @ ..] = e.args() else {
        return refuse(ev, e, INVALID_DIMENSION);
    };
    if rules.is_empty() {
        return refuse(ev, e, INVALID_DIMENSION);
    }
    if !is_list(l) {
        return Ok(None);
    }
    let rules: Vec<(Expr, Expr)> = match rules {
        [k, v] if sides(k).is_none() => vec![(k.clone(), v.clone())],
        _ => rules
            .iter()
            .filter_map(|r| sides(r).map(|(a, b)| (a.clone(), b.clone())))
            .collect(),
    };
    match substituted(l.args(), &rules) {
        Ok(out) => Ok(Some(list(out.into_iter().flatten().collect()))),
        Err(()) => refuse(ev, e, INVALID_DIMENSION),
    }
}

/// Why a heap operation stopped: a comparison neither `True` nor `False`,
/// or an error.
enum Halt {
    Undecided,
    Failed(EvalError),
}

impl From<EvalError> for Halt {
    fn from(err: EvalError) -> Self {
        Self::Failed(err)
    }
}

/// Whether `a` comes before `b`: `Less[a, b]`, or `f(a, b)`, evaluated.
fn before(ev: &mut Eval<'_>, f: Option<&Expr>, a: &Expr, b: &Expr) -> Result<bool, Halt> {
    let test = match f {
        Some(f) => Expr::apply(f, &[a.clone(), b.clone()]),
        None => Expr::call(Symbol::LESS, [a.clone(), b.clone()]),
    };
    match ev.eval(&test)?.as_symbol() {
        Some(Symbol::TRUE) => Ok(true),
        Some(Symbol::FALSE) => Ok(false),
        _ => Err(Halt::Undecided),
    }
}

type Before<'f> = dyn FnMut(&Expr, &Expr) -> Result<bool, Halt> + 'f;

/// libstdc++'s `__push_heap`: `value` up from `hole` to at most `top`.
fn sift_up(
    v: &mut [Expr],
    mut hole: usize,
    top: usize,
    value: Expr,
    less: &mut Before<'_>,
) -> Result<(), Halt> {
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&v[parent], &value)? {
            break;
        }
        v[hole] = v[parent].clone();
        hole = parent;
    }
    v[hole] = value;
    Ok(())
}

/// libstdc++'s `__adjust_heap` over the first `len` elements.
fn adjust(
    v: &mut [Expr],
    mut hole: usize,
    len: usize,
    value: Expr,
    less: &mut Before<'_>,
) -> Result<(), Halt> {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[second], &v[second - 1])? {
            second -= 1;
        }
        v[hole] = v[second].clone();
        hole = second;
    }
    if len % 2 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[hole] = v[second - 1].clone();
        hole = second - 1;
    }
    sift_up(v, hole, top, value, less)
}

/// libstdc++'s `make_heap`.
fn make_heap(v: &mut [Expr], less: &mut Before<'_>) -> Result<(), Halt> {
    let len = v.len();
    if len < 2 {
        return Ok(());
    }
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[parent].clone();
        adjust(v, parent, len, value, less)?;
        if parent == 0 {
            return Ok(());
        }
        parent -= 1;
    }
}

/// libstdc++'s `pop_heap`, then the last element dropped.
fn pop_heap(v: &mut Vec<Expr>, less: &mut Before<'_>) -> Result<(), Halt> {
    if v.len() > 1 {
        let last = v.len() - 1;
        let value = v[last].clone();
        v[last] = v[0].clone();
        adjust(v, 0, last, value, less)?;
    }
    v.pop();
    Ok(())
}

/// libstdc++'s `push_heap` of `x` appended.
fn push_heap(v: &mut Vec<Expr>, x: Expr, less: &mut Before<'_>) -> Result<(), Halt> {
    let hole = v.len();
    v.push(x.clone());
    sift_up(v, hole, 0, x, less)
}

/// `operation` on `l`'s elements ordered by `f`: the list, `None` when a
/// comparison is undecided.
fn heap(
    ev: &mut Eval<'_>,
    l: &Expr,
    f: Option<&Expr>,
    operation: impl FnOnce(&mut Vec<Expr>, &mut Before<'_>) -> Result<(), Halt>,
) -> Native {
    let mut v = l.args().to_vec();
    let mut less = |a: &Expr, b: &Expr| before(ev, f, a, b);
    match operation(&mut v, &mut less) {
        Ok(()) => Ok(Some(list(v))),
        Err(Halt::Undecided) => Ok(None),
        Err(Halt::Failed(err)) => Err(err),
    }
}

/// `Heapify[l, f]`.
pub(in crate::native) fn heapify(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, f) = match e.args() {
        [l] => (l, None),
        [l, f] => (l, Some(f)),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    heap(ev, l, f, |v, less| make_heap(v, less))
}

/// `HeapPush[l, x, f]`.
pub(in crate::native) fn heap_push(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, x, f) = match e.args() {
        [l, x] => (l, x, None),
        [l, x, f] => (l, x, Some(f)),
        [_] | [] => return refuse(ev, e, BAD_VALUE),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    heap(ev, l, f, |v, less| push_heap(v, x.clone(), less))
}

/// `HeapPop[l, f]`.
pub(in crate::native) fn heap_pop(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (l, f) = match e.args() {
        [l] => (l, None),
        [l, f] => (l, Some(f)),
        _ => return Ok(None),
    };
    if !is_list(l) {
        return refuse(ev, e, BAD_VALUE);
    }
    if l.args().is_empty() {
        return refuse(ev, e, INVALID_DIMENSION);
    }
    heap(ev, l, f, |v, less| pop_heap(v, less))
}
