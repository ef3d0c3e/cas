//! Xcas's string commands (CAS-611, CAS-612, CAS-613): `Char`, `Asc`,
//! `Ord`, `InString`, `Contains`, `Mid`, `Left`, `Right`, `Levenshtein`,
//! `Hamming`, `Split`, `StringJoin`, `Trim` and `Strip`, ported from giac
//! 2.0.0.21: ti89.cc `_inString`, `_left`, `_right`, `_mid`, `_ord`;
//! prog.cc `_contains`, `_char`, `_asc`; optimization.cc `_levenshtein`;
//! permu.cc `_split`, `_join`; signalprocessing.cc `_trim`; misc.cc
//! `_strip`. giac has no `hamming` command (an aide entry only): ours counts
//! the positions that differ. Strings in other natives (CAS-612): `Length`,
//! `Join` (concat), `At`, the comparisons and `Sort` take them in their own
//! modules.
//!
//! `split(e, [x, y])` is giac's: the factors of `FactorList[e]` free of `x`
//! multiply to the `y` part, the others free of `y` to the `x` part, and
//! `[0]` when a factor has both. `trim` takes giac's whitespace (space and
//! codes 9 to 13) off a string, and the zeros off the ends of a list: exact
//! 0 or a decimal within giac's `epsilon` 1e-12, or with a threshold `t > 0`
//! the real entries below `t` in absolute value (another entry met makes it
//! stay). `strip` takes spaces, or the given characters.
//!
//! Positions are 0-based and count Unicode scalar values. A call giac
//! refuses, or whose shape is not one of these, stays unevaluated.
//!
//! Departures from giac:
//! - giac counts bytes of UTF-8: `asc("é")` is `[195, 169]` and `ord("é")`
//!   `-61` there, `[233]` and `233` here; `char(233)` is `"é"` here, an
//!   invalid byte there. A code that is not a Unicode scalar value
//!   (negative, a surrogate, past `0x10FFFF`) stays.
//! - `asc(s, i)` out of range stays (giac reads past the string); `ord("")`
//!   and `char(65.5)` stay (giac's errors).
//! - Counts and positions are integers only: giac's `left("abc", 2.0)` is
//!   `"ab"`, ours stays.
//! - `left`/`right` with one argument take an equation (`Equation` or
//!   `Equal`) and a range `a..b`; giac's other one-argument forms (bits of
//!   an integer, a list, the first operand of an expression) stay.
//! - `levenshtein` of a string and a list stays with giac's message as a
//!   warning; past the `terms` limit (`n·m` cells) it stays.
//! - `hamming` answers equal-length strings or lists, where giac stays.
//! - `split(s, "")` stays (giac's "Invalid dimension"); `join` of anything
//!   but strings stays (giac prints them); a `trim` threshold that is not a
//!   positive real stays (giac's "Invalid threshold specification").
//! - Results past the `terms` limit (pieces, characters) stay.

use std::cmp::Ordering;

use cas_num::Rational;
use cas_num::ball::Ball;

use super::expressions::EPSILON;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// giac's message for `levenshtein` of a string and a list.
const MIXED: &str = "Both arguments must be either lists or strings";

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn int(n: usize) -> Expr {
    Expr::from(i64::try_from(n).expect("a length in memory"))
}

/// The characters of a string atom.
fn chars(e: &Expr) -> Option<Vec<char>> {
    e.as_str().map(|s| s.chars().collect())
}

fn string(cs: &[char]) -> Expr {
    Expr::string(&cs.iter().collect::<String>())
}

/// The two arguments of `e` when its head is one of `heads`.
fn pair(e: &Expr, heads: &[Symbol]) -> Option<(Expr, Expr)> {
    let head = e.head()?.as_symbol()?;
    match e.args() {
        [a, b] if heads.contains(&head) => Some((a.clone(), b.clone())),
        _ => None,
    }
}

/// `Char[n]`, `Char[List[n…]]`: the string of these code points.
pub(super) fn char(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    let codes = if is_list(a) {
        a.args()
    } else {
        std::slice::from_ref(a)
    };
    if codes.len() > ev.context().limits().terms {
        return Ok(None);
    }
    let mut s = String::with_capacity(codes.len());
    for c in codes {
        let Some(c) = c
            .as_i64()
            .and_then(|c| u32::try_from(c).ok())
            .and_then(char::from_u32)
        else {
            return Ok(None);
        };
        s.push(c);
    }
    Ok(Some(Expr::string(&s)))
}

/// `Asc[s]`: the code points of `s`; `Asc[s, i]`: the one at `i`.
pub(super) fn asc(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (s, i) = match e.args() {
        [s] => (s, None),
        [s, i] => (s, Some(i)),
        _ => return Ok(None),
    };
    let Some(cs) = chars(s) else {
        return Ok(None);
    };
    let code = |c: char| Expr::from(i64::from(u32::from(c)));
    let Some(i) = i else {
        return Ok(Some(list(cs.into_iter().map(code).collect())));
    };
    let Some(c) = i
        .as_i64()
        .and_then(|i| usize::try_from(i).ok())
        .and_then(|i| cs.get(i))
    else {
        return Ok(None);
    };
    Ok(Some(code(*c)))
}

/// The code of the first character of a string, through lists.
fn ord_of(x: &Expr) -> Option<Expr> {
    if is_list(x) {
        return x.args().iter().map(ord_of).collect::<Option<_>>().map(list);
    }
    let c = x.as_str()?.chars().next()?;
    Some(Expr::from(i64::from(u32::from(c))))
}

/// `Ord[s]`: the code of the first character of `s`; a list maps.
pub(super) fn ord(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a] = e.args() else {
        return Ok(None);
    };
    Ok(ord_of(a))
}

/// The first position of `sub` in `s` from `start`, in characters.
fn find(s: &[char], sub: &[char], start: usize) -> Option<usize> {
    if sub.len() > s.len() {
        return None;
    }
    (start..=s.len() - sub.len()).find(|&k| s[k..k + sub.len()] == *sub)
}

/// `InString[s, sub]`, `InString[s, sub, start]`: the position of `sub` in
/// `s` from `start`, or -1.
pub(super) fn in_string(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (s, sub, start) = match e.args() {
        [s, sub] => (s, sub, Some(0)),
        [s, sub, start] => (s, sub, start.as_i64()),
        _ => return Ok(None),
    };
    let (Some(s), Some(sub), Some(start)) = (chars(s), chars(sub), start) else {
        return Ok(None);
    };
    // A negative start wraps to giac's npos, past the end.
    let found = usize::try_from(start)
        .ok()
        .filter(|&start| start <= s.len())
        .and_then(|start| find(&s, &sub, start));
    Ok(Some(found.map_or(Expr::from(-1), int)))
}

/// A real number: exact, or a ball's midpoint.
enum Real {
    Exact(Rational),
    Float(f64),
}

fn real(e: &Expr) -> Option<Real> {
    match e.view() {
        View::Ball(b) => Some(Real::Float(b.mid().to_f64())),
        _ => e.to_rational().map(Real::Exact),
    }
}

impl Real {
    fn to_f64(&self) -> f64 {
        match self {
            Real::Exact(q) => Ball::from_rational(q, 64).mid().to_f64(),
            Real::Float(x) => *x,
        }
    }

    fn compare(&self, other: &Real) -> Ordering {
        match (self, other) {
            (Real::Exact(a), Real::Exact(b)) => a.cmp(b),
            _ => self.to_f64().total_cmp(&other.to_f64()),
        }
    }
}

/// `Contains[a, b]`: the position of `b` in `a` from 1, or 0. Strings
/// search a substring; a list whose only element is `Range[lo, hi]` of
/// reals is giac's interval, 1 when the real `b` is in it; another list
/// the first element structurally equal to `b`.
pub(super) fn contains(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    if let (Some(s), Some(sub)) = (chars(a), chars(b)) {
        return Ok(Some(
            find(&s, &sub, 0).map_or(Expr::from(0), |k| int(k + 1)),
        ));
    }
    if !is_list(a) {
        return Ok(None);
    }
    if let [range] = a.args()
        && let Some((lo, hi)) = pair(range, &[Symbol::RANGE])
        && let (Some(lo), Some(hi)) = (real(&lo), real(&hi))
    {
        if is_list(b) || b.head() == Some(&Expr::from(Symbol::RANGE)) {
            return Ok(None);
        }
        let inside = real(b).is_some_and(|x| lo.compare(&x).is_le() && x.compare(&hi).is_le());
        return Ok(Some(Expr::from(i64::from(inside))));
    }
    let k = a.args().iter().position(|x| x == b);
    Ok(Some(k.map_or(Expr::from(0), |k| int(k + 1))))
}

/// A string's characters or a list's elements, and how to rebuild them.
enum Seq {
    Chars(Vec<char>),
    Items(Vec<Expr>),
}

impl Seq {
    fn of(e: &Expr) -> Option<Seq> {
        if let Some(cs) = chars(e) {
            Some(Seq::Chars(cs))
        } else if is_list(e) {
            Some(Seq::Items(e.args().to_vec()))
        } else {
            None
        }
    }

    fn len(&self) -> usize {
        match self {
            Seq::Chars(cs) => cs.len(),
            Seq::Items(xs) => xs.len(),
        }
    }

    /// The elements in `from..to`, of the same kind.
    fn slice(&self, from: usize, to: usize) -> Expr {
        match self {
            Seq::Chars(cs) => string(&cs[from..to]),
            Seq::Items(xs) => list(xs[from..to].to_vec()),
        }
    }
}

/// `n` clamped to `0..=len`.
fn clamp(n: i64, len: usize) -> usize {
    usize::try_from(n.max(0)).map_or(len, |n| n.min(len))
}

/// `Mid[x, pos]`, `Mid[x, pos, n]`: `n` elements of a string or list from
/// `pos`, or all the rest; empty from a position out of range.
pub(super) fn mid(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, pos, n) = match e.args() {
        [x, pos] => (x, pos.as_i64(), None),
        [x, pos, n] => match n.as_i64() {
            Some(n) => (x, pos.as_i64(), Some(n)),
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    let (Some(pos), Some(seq)) = (pos, Seq::of(x)) else {
        return Ok(None);
    };
    let len = seq.len();
    let Some(pos) = usize::try_from(pos).ok().filter(|&pos| pos < len) else {
        return Ok(Some(seq.slice(0, 0)));
    };
    let to = n.map_or(len, |n| pos.saturating_add(clamp(n, len)).min(len));
    Ok(Some(seq.slice(pos, to)))
}

/// `Left[x, n]`, `Right[x, n]` of a string or a list; `Left[x]`, `Right[x]`
/// of an equation or a range: its sides.
fn side(e: &Expr, left: bool) -> Native {
    let (x, n) = match e.args() {
        [x] => {
            let heads = [Symbol::EQUATION, Symbol::EQUAL, Symbol::RANGE];
            return Ok(pair(x, &heads).map(|(l, r)| if left { l } else { r }));
        }
        [x, n] => (x, n),
        _ => return Ok(None),
    };
    let (Some(n), Some(seq)) = (n.as_i64(), Seq::of(x)) else {
        return Ok(None);
    };
    let len = seq.len();
    Ok(Some(if left {
        // giac's substr(0, n): a negative n wraps to npos on strings.
        let n = if n < 0 && matches!(seq, Seq::Chars(_)) {
            len
        } else {
            clamp(n, len)
        };
        seq.slice(0, n)
    } else {
        seq.slice(len - clamp(n, len), len)
    }))
}

/// `Left[x, n]`: the first `n` elements; `Left[a = b]` is `a`.
pub(super) fn left(_: &mut Eval<'_>, e: &Expr) -> Native {
    side(e, true)
}

/// `Right[x, n]`: the last `n` elements; `Right[a = b]` is `b`.
pub(super) fn right(_: &mut Eval<'_>, e: &Expr) -> Native {
    side(e, false)
}

/// The edit distance of `a` and `b`, two rows of the table.
fn edit_distance<T: PartialEq>(a: &[T], b: &[T]) -> usize {
    let (a, b) = if a.len() > b.len() { (b, a) } else { (a, b) };
    let mut v0: Vec<usize> = (0..=b.len()).collect();
    let mut v1 = vec![0; b.len() + 1];
    for (i, x) in a.iter().enumerate() {
        v1[0] = i + 1;
        for (j, y) in b.iter().enumerate() {
            v1[j + 1] = (v0[j + 1] + 1)
                .min(v1[j] + 1)
                .min(v0[j] + usize::from(x != y));
        }
        std::mem::swap(&mut v0, &mut v1);
    }
    v0[b.len()]
}

/// `Levenshtein[a, b]`: the edit distance of two strings or two lists.
pub(super) fn levenshtein(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let distance = match (Seq::of(a), Seq::of(b)) {
        (Some(Seq::Chars(a)), Some(Seq::Chars(b))) => {
            fits(ev, a.len(), b.len()).then(|| edit_distance(&a, &b))
        }
        (Some(Seq::Items(a)), Some(Seq::Items(b))) => {
            fits(ev, a.len(), b.len()).then(|| edit_distance(&a, &b))
        }
        (Some(_), Some(_)) => {
            ev.warn(Warning::Native {
                message: MIXED.into(),
                expr: e.clone(),
            });
            None
        }
        _ => None,
    };
    Ok(distance.map(int))
}

/// Whether an `n` by `m` table is within the `terms` limit.
fn fits(ev: &Eval<'_>, n: usize, m: usize) -> bool {
    n.checked_mul(m)
        .is_some_and(|cells| cells <= ev.context().limits().terms)
}

/// `Hamming[a, b]`: the positions where two strings or two lists of the
/// same length differ.
pub(super) fn hamming(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    let differ = match (Seq::of(a), Seq::of(b)) {
        (Some(Seq::Chars(a)), Some(Seq::Chars(b))) if a.len() == b.len() => {
            a.iter().zip(&b).filter(|(x, y)| x != y).count()
        }
        (Some(Seq::Items(a)), Some(Seq::Items(b))) if a.len() == b.len() => {
            a.iter().zip(&b).filter(|(x, y)| x != y).count()
        }
        _ => return Ok(None),
    };
    Ok(Some(int(differ)))
}

/// Whether `n` elements are within the `terms` limit.
fn within(ev: &Eval<'_>, n: usize) -> bool {
    n <= ev.context().limits().terms
}

/// `Split[s, sep]`: the pieces of `s` between the occurrences of `sep`,
/// a last empty piece dropped; `Split[e, List[x, y]]`: `e` as a product
/// `List[f(x), g(y)]` of its factors, or `List[0]` when it is not one.
pub(super) fn split(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [s, sep] = e.args() else {
        return Ok(None);
    };
    if let [x, y] = sep.args()
        && is_list(sep)
    {
        return split_variables(ev, s, x, y);
    }
    let (Some(s), Some(sep)) = (s.as_str(), sep.as_str()) else {
        return Ok(None);
    };
    if sep.is_empty() || !within(ev, s.matches(sep).count().saturating_add(1)) {
        return Ok(None);
    }
    let mut pieces: Vec<&str> = s.split(sep).collect();
    if pieces.last() == Some(&"") {
        pieces.pop();
    }
    Ok(Some(list(pieces.into_iter().map(Expr::string).collect())))
}

/// giac's `split(e, [x, y])`: the factors of `FactorList[e]` free of `x`
/// make the `y` part, the others free of `y` the `x` part.
fn split_variables(ev: &mut Eval<'_>, e: &Expr, x: &Expr, y: &Expr) -> Native {
    if e.as_str().is_some() || is_list(e) {
        return Ok(None);
    }
    let factors = ev.eval(&Expr::call(Symbol::FACTOR_LIST, [e.clone()]))?;
    if !is_list(&factors) || factors.args().len() % 2 != 0 {
        return Ok(None);
    }
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for pair in factors.args().chunks_exact(2) {
        let power = Expr::call(Symbol::POWER, [pair[0].clone(), pair[1].clone()]);
        if pair[0].free_of(x) {
            ys.push(power);
        } else if pair[0].free_of(y) {
            xs.push(power);
        } else {
            return Ok(Some(list(vec![Expr::from(0)])));
        }
    }
    let xs = ev.eval(&Expr::call(Symbol::TIMES, xs))?;
    let ys = ev.eval(&Expr::call(Symbol::TIMES, ys))?;
    Ok(Some(list(vec![xs, ys])))
}

/// `StringJoin[sep, List[s…]]`: the strings `s…` with `sep` between them.
pub(super) fn string_join(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [sep, items] = e.args() else {
        return Ok(None);
    };
    let (Some(sep), true) = (sep.as_str(), is_list(items)) else {
        return Ok(None);
    };
    let Some(items) = items
        .args()
        .iter()
        .map(Expr::as_str)
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    let gaps = items.len().saturating_sub(1);
    let length = items
        .iter()
        .try_fold(0usize, |n, s| n.checked_add(s.chars().count()))
        .and_then(|n| gaps.checked_mul(sep.chars().count())?.checked_add(n));
    if length.is_none_or(|n| !within(ev, n)) {
        return Ok(None);
    }
    Ok(Some(Expr::string(&items.join(sep))))
}

/// giac's whitespace: space, tab, line feed, vertical tab, form feed and
/// carriage return.
fn is_whitespace(c: char) -> bool {
    c == ' ' || ('\t'..='\r').contains(&c)
}

/// Whether `x` is zero for `trim`: below the threshold `t` in absolute
/// value, `None` when `x` is not a real number; without one exactly 0, or a
/// decimal within [`EPSILON`].
fn trimmed(x: &Expr, t: Option<&Real>) -> Option<bool> {
    let Some(t) = t else {
        return Some(match x.view() {
            View::Ball(b) => b.mid().to_f64().abs() <= EPSILON,
            _ => x.to_rational().is_some_and(|q| q.is_zero()),
        });
    };
    let size = match real(x)? {
        Real::Exact(q) => Real::Exact(q.abs()),
        Real::Float(f) => Real::Float(f.abs()),
    };
    Some(size.compare(t).is_lt())
}

/// `Trim[s]`: `s` without whitespace at its ends; `Trim[v]`, `Trim[v, t]`:
/// the list `v` without its zeros at its ends, below `t` when given. A last
/// option `Left` or `Right` trims that end only.
pub(super) fn trim(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, options @ ..] = e.args() else {
        return Ok(None);
    };
    let (options, left, right) = match options.split_last() {
        Some((side, rest)) if side.as_symbol() == Some(Symbol::LEFT) => (rest, true, false),
        Some((side, rest)) if side.as_symbol() == Some(Symbol::RIGHT) => (rest, false, true),
        _ => (options, true, true),
    };
    if let Some(s) = a.as_str() {
        if !options.is_empty() {
            return Ok(None);
        }
        let s = if left {
            s.trim_start_matches(is_whitespace)
        } else {
            s
        };
        let s = if right {
            s.trim_end_matches(is_whitespace)
        } else {
            s
        };
        return Ok(Some(Expr::string(s)));
    }
    if !is_list(a) {
        return Ok(None);
    }
    let threshold = match options {
        [] => None,
        [t] => match real(t) {
            Some(Real::Exact(q)) if !q.is_negative() && !q.is_zero() => Some(Real::Exact(q)),
            Some(Real::Float(f)) if f > 0.0 => Some(Real::Float(f)),
            _ => return Ok(None),
        },
        _ => return Ok(None),
    };
    let items = a.args();
    let (mut i, mut j) = (0, items.len());
    if left {
        while i < j {
            let Some(zero) = trimmed(&items[i], threshold.as_ref()) else {
                return Ok(None);
            };
            if !zero {
                break;
            }
            i += 1;
        }
    }
    if right {
        while j > i {
            let Some(zero) = trimmed(&items[j - 1], threshold.as_ref()) else {
                return Ok(None);
            };
            if !zero {
                break;
            }
            j -= 1;
        }
    }
    Ok(Some(list(items[i..j].to_vec())))
}

/// `Strip[s]`, `Strip[s, chars]`: `s` without the characters of `chars`, a
/// space by default, at its ends.
pub(super) fn strip(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (s, set) = match e.args() {
        [s] => (s.as_str(), Some(" ")),
        [s, set] => (s.as_str(), set.as_str()),
        _ => return Ok(None),
    };
    let (Some(s), Some(set)) = (s, set) else {
        return Ok(None);
    };
    Ok(Some(Expr::string(s.trim_matches(|c| set.contains(c)))))
}

#[cfg(test)]
mod tests;
