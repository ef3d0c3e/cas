//! giac's list, type, format, date and polynomial leftovers of the gap
//! ledger (CAS-644, CAS-645, CAS-647), under giac's names, ported from giac
//! 2.0.0.21: plot.cc `_back`; misc.cc `_range`, `_suppress`, `_dayofweek`,
//! `_truncate`, `_roots` (`_froot`, `addfactors`); ti89.cc `_getType`,
//! `_format` (gen.cc `print_DOUBLE_`), `_shift` (`shiftrotate`); usual.cc
//! `_rdiv`; ifactor.cc `_divis` (`divis`). Where giac errs, the call stays
//! with giac's message as a warning:
//!
//! - `LastElement[l]` (`back`, CAS-644): the last element of a list, the
//!   last character of a string, the last argument of several; an empty
//!   list or string and anything else is itself;
//! - `ArithmeticRange[n]`, `[a, b]`, `[a, b, p]` (`range`, CAS-644): `[0,
//!   …, n − 1]`, or `a, a + p, …` before `b`, `p` 1 by default, with `a ≤
//!   b, p > 0` or `a ≥ b, p < 0`; integers when all are, else decimals
//!   (giac's doubles) and `⌈(b − a)/p⌉` of them;
//! - `ShiftList[l, n]` (`shift`, CAS-644): `n > 0` drops `n` elements in
//!   front and pads `n` zeros at the end, `n < 0` the other way, `n` −1
//!   when left out or not an integer, at most the length. A string pads
//!   spaces, keeping, for `n < 0`, its first `|n|` characters at the end as
//!   giac's; an integer shifts its bits, `⌊a/2ⁿ⌋` or `a·2^|n|`;
//!   anything else is itself;
//! - `DeleteAt[l, i]` (`suppress`, CAS-644): `l` without the element at the
//!   0-based `i`, unchanged out of range; a list of indices removes each, a
//!   range `a..b` or `DeleteAt[l, a, b]` the span (out of range: giac's
//!   "Invalid dimension"). A string without its `i`th character, `""` out
//!   of range as giac's, or without the characters of a string `i`;
//! - `ElementAt[l, j]`, `[m, [j, k]]` (`at`, CAS-644): `At[l, j]`,
//!   `At[m, j, k]`, 0-based; an integer index out of range stays with
//!   giac's "Index outside range";
//! - `Divide[a, b]` (`rdiv`, CAS-644): `a/b`;
//! - `GetType[e]` (`getType`, CAS-644): giac's type name, `"STR"`, `"NUM"`
//!   (reals, `true`, `false`), `"LIST"`, `"MAT"`, `"VAR"` (names, `pi`,
//!   `euler_gamma`, `undef`, `infinity`, `inf`), `"FUNC"` (a function),
//!   `"EXPR"` (complex numbers, `e`, other expressions) or `"OTHER"` (a
//!   builtin's name);
//! - `FormatNumber[x, f]` (`format`, CAS-644): the real `x` as C's `printf`
//!   prints its double with `%.nf` (`f`), `%.ne` (`e`) or `%.ng` (any other
//!   letter), `f` a letter and one or two digits (17 at most), else
//!   `"invalid format"`; giac's double printing first clears the five low
//!   bits of the mantissa. An integer beyond 32 bits keeps its digits, as
//!   giac's: cut to `d.ddde+k` for `e` and `s`. Hexadecimal `a` and a
//!   non-real `x` stay;
//! - `DayOfWeek[d, m, y]` (`dayofweek`, CAS-647): 0 for Sunday to 6, the
//!   proleptic Gregorian calendar, giac's formula (`29, 2, 2001` is March
//!   1); a day outside 1..31 or a month outside 1..12 stays;
//! - `RootsWithMultiplicity[p, x]` (`roots`, CAS-645): `[[r, m], …]`, the
//!   complex roots of the numerator of `p` in `x` (`x` by default), each
//!   irreducible factor of `FactorList` solved by `CSolve`; giac solves its
//!   square-free factors, so the order departs;
//! - `TruncateDegree[p, n]`, `[p, vars, n]` (`truncate`, CAS-645): the
//!   terms of the expanded numerator whose degree in `vars` (`x` by
//!   default, other factors degree 0) is at most `n`, 5 by default, over
//!   the denominator; a denominator of positive degree in the first
//!   variable, or `n` no integer, stays with "Bad Argument Value";
//! - `PolynomialDivisors[p]` (`divis`, CAS-645): the divisors of `p` from
//!   `FactorList`, giac's enumeration: `[1]`, then for each factor `d` of
//!   multiplicity `k` the list so far times `d`, …, `dᵏ` appended; a list
//!   maps;
//! - `RandomChoice[l]` (`choice`, CAS-646, misc.cc `_choice`): `Random[l]`,
//!   an element of the non-empty list `l`, the same draw from the session
//!   generator; anything else stays with "Bad Argument Value".
//!
//! giac's `subs` is `subst` (`ReplaceAll`); `subs(x=2, e)`, giac's error,
//! stays. `convert(e, list)` is `e` (`Convert`).

use cas_num::ball::{Ball, Float};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `gensizeerr` message.
const BAD: &str = "Bad Argument Value";
/// giac's `gendimerr` message.
const DIMENSION: &str = "Invalid dimension";

/// The call stays, with giac's error as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// Whether `x` is a list.
fn is_list(x: &Expr) -> bool {
    x.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// An integer of giac's 32-bit `_INT_`.
fn int32(x: &Expr) -> Option<i64> {
    x.as_i64().filter(|n| i32::try_from(*n).is_ok())
}

/// Whether `count` elements fit the `terms` limit.
fn fits(ev: &Eval<'_>, count: impl TryInto<usize>) -> bool {
    count
        .try_into()
        .is_ok_and(|c| c <= ev.context().limits().terms)
}

/// The double of a real number, exact or a ball.
fn double(x: &Expr) -> Option<f64> {
    match x.view() {
        View::Ball(b) => Some(b.mid().to_f64()),
        _ => x
            .to_rational()
            .map(|q| Ball::from_rational(&q, 53).mid().to_f64()),
    }
}

/// giac's roundoff (gen.cc `print_DOUBLE_`): the five low bits of the
/// mantissa cleared.
fn roundoff(x: f64) -> f64 {
    f64::from_bits(x.to_bits() & !0x1f)
}

/// A decimal of the double `x`.
fn decimal(x: f64) -> Expr {
    Float::from_f64(x).map_or_else(|| Expr::from(0), |f| Expr::from(Ball::exact(f)))
}

/// `LastElement[l]`, giac's `back`.
pub(super) fn last_element(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [] => Ok(None),
        [a] => {
            if let Some(s) = a.as_str() {
                return Ok(Some(
                    s.chars()
                        .last()
                        .map_or_else(|| a.clone(), |c| Expr::string(&c.to_string())),
                ));
            }
            if is_list(a) {
                return Ok(Some(a.args().last().unwrap_or(a).clone()));
            }
            Ok(Some(a.clone()))
        }
        [.., last] => Ok(Some(last.clone())),
    }
}

/// `ArithmeticRange[n]`, `[a, b]`, `[a, b, p]`, giac's `range`.
pub(super) fn arithmetic_range(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (from, to, by) = match e.args() {
        [n] => {
            let Some(n) = n.as_i64().filter(|n| *n >= 0) else {
                return refuse(ev, e, BAD);
            };
            if !fits(ev, n) {
                return Ok(None);
            }
            return Ok(Some(list((0..n).map(Expr::from).collect())));
        }
        [from, to] => (from, to, Expr::from(1)),
        [from, to, by] => (from, to, by.clone()),
        _ => return Ok(None),
    };
    if let (Some(start), Some(end), Some(step)) = (int32(from), int32(to), int32(&by)) {
        if (start <= end && step > 0) || (start >= end && step < 0) {
            // ⌈(end − start)/step⌉, the two of one sign; 32-bit inputs.
            let count = (end - start + step - step.signum()) / step;
            if !fits(ev, count) {
                return Ok(None);
            }
            return Ok(Some(list(
                (0..count).map(|i| Expr::from(start + i * step)).collect(),
            )));
        }
        return refuse(ev, e, BAD);
    }
    let (Some(start), Some(end), Some(step)) = (double(from), double(to), double(&by)) else {
        return refuse(ev, e, BAD);
    };
    // giac's doubles come with their low bits cleared: its `range(0,1,1/3)`
    // has 4 elements, `range(0,1,0.1)` 11.
    let (start, end, step) = (roundoff(start), roundoff(end), roundoff(step));
    if !((start <= end && step > 0.0) || (start >= end && step < 0.0)) {
        return refuse(ev, e, BAD);
    }
    // giac's doubles: `start + i·step` for `i < ⌈(end − start)/step⌉`.
    let count = ((end - start) / step).ceil();
    let terms = ev.context().limits().terms;
    let mut items = Vec::new();
    let mut i = 0.0;
    while i < count {
        if items.len() >= terms {
            return Ok(None);
        }
        items.push(decimal(start + i * step));
        i += 1.0;
    }
    Ok(Some(list(items)))
}

/// `ShiftList[l, n]`, giac's `shift`.
pub(super) fn shift_list(_: &mut Eval<'_>, e: &Expr) -> Native {
    let (a, n) = match e.args() {
        [a] => (a, -1),
        [a, n] => (a, int32(n).unwrap_or(-1)),
        _ => return Ok(None),
    };
    // giac's `shiftrotate(g, true)`: right by `n`, left for `n < 0`.
    let right = n >= 0;
    // `n` is a 32-bit integer.
    let n = n.abs();
    if let Some(k) = int32(a) {
        if right {
            return Ok(Some(Expr::from(k >> n.min(63))));
        }
        if n > 64 {
            return Ok(None);
        }
        return Ok(Some(a.clone() * Expr::from(2).pow(Expr::from(n))));
    }
    let n = usize::try_from(n).unwrap_or(usize::MAX);
    if is_list(a) {
        let items = a.args();
        let n = n.min(items.len());
        let zeros = std::iter::repeat_n(Expr::from(0), n);
        let shifted: Vec<Expr> = if right {
            items[n..].iter().cloned().chain(zeros).collect()
        } else {
            zeros
                .chain(items[..items.len() - n].iter().cloned())
                .collect()
        };
        return Ok(Some(list(shifted)));
    }
    if let Some(s) = a.as_str() {
        let chars: Vec<char> = s.chars().collect();
        let n = n.min(chars.len());
        let spaces = |k: usize| " ".repeat(k);
        let shifted = if right {
            chars[n..].iter().collect::<String>() + &spaces(n)
        } else {
            spaces(chars.len() - n) + &chars[..n].iter().collect::<String>()
        };
        return Ok(Some(Expr::string(&shifted)));
    }
    Ok(Some(a.clone()))
}

/// `l` without its elements or characters `i..=j`, giac's three-argument
/// `suppress`.
fn delete_span(ev: &mut Eval<'_>, e: &Expr, l: &Expr, i: i64, j: i64) -> Native {
    let (Ok(i), Ok(j)) = (usize::try_from(i), usize::try_from(j)) else {
        return refuse(ev, e, DIMENSION);
    };
    if i > j {
        return refuse(ev, e, DIMENSION);
    }
    if is_list(l) {
        let items = l.args();
        if j >= items.len() {
            return refuse(ev, e, DIMENSION);
        }
        let kept = items[..i].iter().chain(&items[j + 1..]).cloned().collect();
        return Ok(Some(list(kept)));
    }
    if let Some(s) = l.as_str() {
        let chars: Vec<char> = s.chars().collect();
        if j >= chars.len() {
            return refuse(ev, e, DIMENSION);
        }
        let kept: String = chars[..i].iter().chain(&chars[j + 1..]).collect();
        return Ok(Some(Expr::string(&kept)));
    }
    refuse(ev, e, BAD)
}

/// `l` without its element or character `i`, or a string without the
/// characters of the string `i`; `None` for giac's error.
fn delete_one(l: &Expr, i: &Expr) -> Option<Expr> {
    if let Some(s) = l.as_str() {
        let chars: Vec<char> = s.chars().collect();
        let kept: String = if let Some(remove) = i.as_str() {
            chars.iter().filter(|c| !remove.contains(**c)).collect()
        } else {
            match int32(i).and_then(|k| usize::try_from(k).ok()) {
                Some(k) if k < chars.len() => chars[..k].iter().chain(&chars[k + 1..]).collect(),
                // giac's empty result.
                _ => String::new(),
            }
        };
        return Some(Expr::string(&kept));
    }
    let k = int32(i)?;
    if !is_list(l) {
        return None;
    }
    let kept = l
        .args()
        .iter()
        .enumerate()
        .filter(|(j, _)| i64::try_from(*j).ok() != Some(k))
        .map(|(_, x)| x.clone())
        .collect();
    Some(list(kept))
}

/// `DeleteAt[l, i]`, `DeleteAt[l, i, j]`, giac's `suppress`.
pub(super) fn delete_at(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l, i, j] => match (int32(i), int32(j)) {
            (Some(i), Some(j)) => delete_span(ev, e, l, i, j),
            _ => refuse(ev, e, BAD),
        },
        [l, i] => {
            if i.head().and_then(Expr::as_symbol) == Some(Symbol::RANGE)
                && let [a, b] = i.args()
                && let (Some(a), Some(b)) = (int32(a), int32(b))
            {
                return delete_span(ev, e, l, a, b);
            }
            if is_list(i) {
                // giac sorts the indices down and removes each.
                let Some(mut ks) = i.args().iter().map(int32).collect::<Option<Vec<_>>>() else {
                    return refuse(ev, e, BAD);
                };
                ks.sort_unstable_by(|a, b| b.cmp(a));
                let mut l = l.clone();
                for k in ks {
                    let Some(next) = delete_one(&l, &Expr::from(k)) else {
                        return refuse(ev, e, BAD);
                    };
                    l = next;
                }
                return Ok(Some(l));
            }
            match delete_one(l, i) {
                Some(found) => Ok(Some(found)),
                None => refuse(ev, e, BAD),
            }
        }
        _ => Ok(None),
    }
}

/// `ElementAt[l, j]`, `ElementAt[m, [j, k]]`, giac's `at`.
pub(super) fn element_at(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [v, index] = e.args() else {
        return Ok(None);
    };
    let indices = if is_list(index) {
        index.args().to_vec()
    } else {
        vec![index.clone()]
    };
    if indices.is_empty() {
        return Ok(None);
    }
    // Integer indices out of range are giac's error, before `At`'s.
    let mut level = Some(v.clone());
    for k in &indices {
        let Some(at) = level.take() else {
            break;
        };
        let len = if is_list(&at) {
            at.args().len()
        } else if let Some(s) = at.as_str() {
            s.chars().count()
        } else {
            break;
        };
        let Some(k) = k.as_i64() else {
            break;
        };
        let len = i64::try_from(len).unwrap_or(i64::MAX);
        if k >= len || k < -len {
            return refuse(ev, e, "Index outside range");
        }
        if is_list(&at) {
            let k = if k < 0 { k + len } else { k };
            level = usize::try_from(k).ok().map(|k| at.args()[k].clone());
        }
    }
    let mut args = vec![v.clone()];
    args.extend(indices);
    Ok(Some(Expr::call(Symbol::AT, args)))
}

/// `Divide[a, b]`, giac's `rdiv`.
pub(super) fn divide(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    Ok(Some(a.clone() / b.clone()))
}

/// Whether `x` is a matrix as giac's `ckmatrix`: a list of lists of one
/// nonzero length.
fn is_matrix(x: &Expr) -> bool {
    let rows = x.args();
    is_list(x)
        && !rows.is_empty()
        && rows.iter().all(is_list)
        && !rows[0].args().is_empty()
        && rows.iter().all(|r| r.args().len() == rows[0].args().len())
}

/// `GetType[e]`, giac's `getType`.
pub(super) fn get_type(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let name = match x.view() {
        View::String(_) => "STR",
        View::Int(_) | View::Integer(_) | View::Rational(_) | View::Ball(_) => "NUM",
        View::Complex(_) | View::ComplexBall(_) => "EXPR",
        View::Symbol(s) => match s {
            Symbol::TRUE | Symbol::FALSE => "NUM",
            Symbol::E | Symbol::I => "EXPR",
            Symbol::PI
            | Symbol::EULER_GAMMA
            | Symbol::INDETERMINATE
            | Symbol::COMPLEX_INFINITY
            | Symbol::INFINITY => "VAR",
            s if !s.is_builtin() => "VAR",
            _ => "OTHER",
        },
        View::Apply { head, .. } => match head.as_symbol() {
            Some(Symbol::LIST) if is_matrix(x) => "MAT",
            Some(Symbol::LIST) => "LIST",
            Some(Symbol::FUNCTION) => "FUNC",
            _ => "EXPR",
        },
    };
    Ok(Some(Expr::string(name)))
}

/// C's `%.{p}e`: the exponent signed, two digits at least.
fn printf_e(x: f64, p: usize, upper: bool) -> String {
    let s = format!("{x:.p$e}");
    let (mantissa, exponent) = s.split_once('e').unwrap_or((&s, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let sign = if exponent < 0 { '-' } else { '+' };
    let e = if upper { 'E' } else { 'e' };
    format!("{mantissa}{e}{sign}{:02}", exponent.unsigned_abs())
}

/// `s` without the trailing zeros of its fraction, and its point then.
fn trim_fraction(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

/// C's `%.{p}g`.
fn printf_g(x: f64, p: usize, upper: bool) -> String {
    let p = p.max(1);
    let probe = format!("{x:.prec$e}", prec = p - 1);
    let exponent: i64 = probe
        .split_once('e')
        .and_then(|(_, k)| k.parse().ok())
        .unwrap_or(0);
    let precision = i64::try_from(p).unwrap_or(i64::MAX);
    if exponent >= -4 && exponent < precision {
        let decimals = usize::try_from(precision - 1 - exponent).unwrap_or(0);
        trim_fraction(&format!("{x:.decimals$}"))
    } else {
        let s = printf_e(x, p - 1, upper);
        let at = s.find(['e', 'E']).unwrap_or(s.len());
        trim_fraction(&s[..at]) + &s[at..]
    }
}

/// giac's `format` of an integer beyond 32 bits: its digits, cut to
/// `d.ddde+k` for a format other than `f`.
fn format_integer(digits: &str, form: &str) -> String {
    let mut chars = form.chars();
    let Some(ch) = chars.next() else {
        return digits.to_string();
    };
    if ch.eq_ignore_ascii_case(&'f') || form.chars().count() < 2 {
        return digits.to_string();
    }
    // C's `atol`, then giac's `unsigned int` arithmetic.
    let rest: String = chars.take_while(char::is_ascii_digit).collect();
    let n = rest.parse::<u32>().unwrap_or(0);
    let mut kept = n.wrapping_sub(1);
    if ch.eq_ignore_ascii_case(&'e') {
        kept = kept.wrapping_add(1);
    }
    let kept = kept.max(2);
    let len = digits.len();
    if (kept.wrapping_add(1) as usize) < len {
        let end = len.min(1 + kept as usize);
        return format!("{}.{}e+{}", &digits[..1], &digits[1..end], len - 1);
    }
    digits.to_string()
}

/// `FormatNumber[x, f]`, giac's `format`.
pub(super) fn format_number(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, form] = e.args() else {
        return Ok(None);
    };
    let Some(form) = form.as_str() else {
        return Ok(None);
    };
    if let Some(n) = x.to_integer()
        && n.to_i64().is_none_or(|k| i32::try_from(k).is_err())
    {
        return Ok(Some(Expr::string(&format_integer(&n.to_string(), form))));
    }
    let chars: Vec<char> = form.chars().collect();
    let Some(&first) = chars.first() else {
        return Ok(None);
    };
    let ch = match first.to_ascii_lowercase() {
        'f' | 'e' | 'g' | 'a' => first,
        _ => 'g',
    };
    if chars.len() < 2 || chars.len() > 3 || !chars[1].is_ascii_digit() {
        return Ok(Some(Expr::string("invalid format")));
    }
    if ch.eq_ignore_ascii_case(&'a') || chars.get(2).is_some_and(|c| !c.is_ascii_digit()) {
        return Ok(None);
    }
    let digits: String = chars[1..].iter().collect();
    let p = digits.parse::<usize>().unwrap_or(0).min(17);
    let value = if x.to_rational().is_some() || matches!(x.view(), View::Ball(_)) {
        x.clone()
    } else {
        ev.eval(&Expr::call(Symbol::N, [x.clone()]))?
    };
    let Some(d) = double(&value).filter(|d| d.is_finite()) else {
        return Ok(None);
    };
    let d = roundoff(d);
    let upper = ch.is_ascii_uppercase();
    let text = match ch.to_ascii_lowercase() {
        'f' => format!("{d:.p$}"),
        'e' => printf_e(d, p, upper),
        _ => printf_g(d, p, upper),
    };
    Ok(Some(Expr::string(&text)))
}

/// `DayOfWeek[d, m, y]`, giac's `dayofweek`.
pub(super) fn day_of_week(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [day, month, year] = e.args() else {
        return refuse(ev, e, BAD);
    };
    let (Some(day), Some(month), Some(year)) = (int32(day), int32(month), int32(year)) else {
        return refuse(ev, e, BAD);
    };
    if !(1..=31).contains(&day) || !(1..=12).contains(&month) {
        return refuse(ev, e, BAD);
    }
    // giac's formula, C's truncating division.
    let before = if month < 3 { year - 1 } else { year };
    let n = (23 * month) / 9 + day + 4 + year + before / 4 - before / 100 + before / 400;
    let weekday = if month < 3 { n % 7 } else { (n - 2) % 7 };
    Ok(Some(Expr::from(weekday)))
}

/// The variable of `args[n]`, `x` when left out; `None` for a builtin or
/// another argument.
fn variable(args: &[Expr], n: usize) -> Option<Expr> {
    match args.get(n) {
        None => Some(Expr::symbol("x")),
        Some(v) if v.as_symbol().is_some_and(|s| !s.is_builtin()) => Some(v.clone()),
        Some(_) => None,
    }
}

/// The flat `[f₁, k₁, f₂, k₂, …]` of `FactorList[p]`, `None` otherwise.
fn factor_pairs(ev: &mut Eval<'_>, p: &Expr) -> Result<Option<Vec<(Expr, i64)>>, EvalError> {
    let factors = ev.eval(&Expr::call(Symbol::FACTOR_LIST, [p.clone()]))?;
    if !is_list(&factors) || factors.args().len() % 2 != 0 {
        return Ok(None);
    }
    Ok(factors
        .args()
        .chunks(2)
        .map(|pair| Some((pair[0].clone(), pair[1].as_i64()?)))
        .collect())
}

/// `RootsWithMultiplicity[p, x]`, giac's `roots`.
pub(super) fn roots_with_multiplicity(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if args.is_empty() || args.len() > 2 {
        return Ok(None);
    }
    let Some(x) = variable(args, 1) else {
        return refuse(ev, e, BAD);
    };
    let Some(pairs) = factor_pairs(ev, &args[0])? else {
        return Ok(None);
    };
    let mut roots = Vec::new();
    for (f, k) in pairs {
        if k <= 0 || !f.contains(&x) {
            continue;
        }
        let equation = Expr::call(Symbol::EQUATION, [f, Expr::from(0)]);
        let found = ev.eval(&Expr::call(Symbol::CSOLVE, [equation, x.clone()]))?;
        if !is_list(&found) {
            return Ok(None);
        }
        for r in found.args() {
            roots.push(list(vec![r.clone(), Expr::from(k)]));
        }
    }
    Ok(Some(list(roots)))
}

/// The degree in `vars` of the term `t`: a factor `v` or `v^k`, `k` a
/// positive integer, counts; any other factor is degree 0, giac's other
/// variables.
fn degree(t: &Expr, vars: &[Expr]) -> i64 {
    let factors = if t.head().and_then(Expr::as_symbol) == Some(Symbol::TIMES) {
        t.args()
    } else {
        std::slice::from_ref(t)
    };
    factors
        .iter()
        .map(|f| {
            if vars.contains(f) {
                return 1;
            }
            match (f.head().and_then(Expr::as_symbol), f.args()) {
                (Some(Symbol::POWER), [b, k]) if vars.contains(b) => {
                    k.as_i64().filter(|k| *k > 0).unwrap_or(0)
                }
                _ => 0,
            }
        })
        .fold(0, i64::saturating_add)
}

/// The terms of `e` once expanded.
fn terms(ev: &mut Eval<'_>, e: &Expr) -> Result<Vec<Expr>, EvalError> {
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [e.clone()]))?;
    Ok(
        if expanded.head().and_then(Expr::as_symbol) == Some(Symbol::PLUS) {
            expanded.args().to_vec()
        } else {
            vec![expanded]
        },
    )
}

/// `TruncateDegree[p, n]`, `TruncateDegree[p, vars, n]`, giac's
/// `truncate`.
pub(super) fn truncate_degree(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, vars, n) = match e.args() {
        [p] => (p, vec![Expr::symbol("x")], Expr::from(5)),
        [p, n] => (p, vec![Expr::symbol("x")], n.clone()),
        [p, v, n] => {
            let vars = if is_list(v) {
                v.args().to_vec()
            } else {
                vec![v.clone()]
            };
            (p, vars, n.clone())
        }
        _ => return Ok(None),
    };
    let Some(n) = int32(&n) else {
        return refuse(ev, e, BAD);
    };
    if vars.is_empty() {
        return Ok(None);
    }
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [p.clone()]))?;
    let num = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let den = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together]))?;
    // giac refuses a denominator of positive degree in the first variable.
    if terms(ev, &den)?.iter().any(|t| degree(t, &vars[..1]) > 0) {
        return refuse(ev, e, BAD);
    }
    let kept: Vec<Expr> = terms(ev, &num)?
        .into_iter()
        .filter(|t| degree(t, &vars) <= n)
        .collect();
    Ok(Some(kept.into_iter().sum::<Expr>() / den))
}

/// `PolynomialDivisors[p]`, giac's `divis`.
pub(super) fn polynomial_divisors(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    if is_list(p) {
        let each = p
            .args()
            .iter()
            .map(|q| Expr::call(Symbol::POLYNOMIAL_DIVISORS, [q.clone()]))
            .collect();
        return Ok(Some(list(each)));
    }
    let Some(pairs) = factor_pairs(ev, p)? else {
        return Ok(None);
    };
    let count = pairs.iter().try_fold(1u128, |acc, (_, k)| {
        acc.checked_mul(u128::try_from(k.saturating_add(1)).unwrap_or(1))
    });
    if !count.is_some_and(|c| fits(ev, c)) {
        return Ok(None);
    }
    let mut divisors = vec![Expr::from(1)];
    for (d, k) in pairs {
        let mut more = Vec::new();
        for j in 1..=k {
            let dj = if j == 1 {
                d.clone()
            } else {
                d.pow(Expr::from(j))
            };
            more.extend(divisors.iter().map(|l| l.clone() * dj.clone()));
        }
        divisors.extend(more);
    }
    Ok(Some(list(divisors)))
}

/// `RandomChoice[l]` (`choice`): `Random[l]` on a non-empty list.
pub(super) fn random_choice(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [l] if is_list(l) && !l.args().is_empty() => super::random::element(ev, l.args()),
        _ => refuse(ev, e, BAD),
    }
}

#[cfg(test)]
mod tests;
