//! Xcas's `make_symbol` and `symbol_array` (CAS-623): `MakeSymbol` and
//! `SymbolArray`, ported from giac 2.0.0.21 ti89.cc `make_symbol` and
//! optimization.cc `_symbol_array`.
//!
//! - `MakeSymbol[s]`: the symbol named `s` when `s` is an identifier, an
//!   ASCII letter or `_` first, then ASCII letters, digits and `_`; another
//!   string stays with a warning (giac: `make_symbol("x y")` undef,
//!   `make_symbol("1a")` `a`, `make_symbol("_x1")` `1_x1`, a unit).
//! - A builtin's spelling: cas-core sees only canonical names
//!   (`Symbol::new`), not the language's (`Names` lives in cas-lang). A
//!   canonical builtin name (`"Plus"`, `"Sin"`) stays with a warning, so
//!   that no call answers a builtin head; a language spelling such as
//!   `"sin"` answers a plain symbol named `sin`, unrelated to `Sin`, as
//!   giac's `make_symbol("sin")` answers the quoted name `'sin'`.
//! - `SymbolArray[s, d₁, …, dₖ]`: nested lists, `d₁` by … by `dₖ`, of the
//!   symbols `s` with the indices from 0: each `%` in `s` replaced by the
//!   next index (`"a%%"` gives `a00`, …; their count must be `k`), and with
//!   no `%`, `s` then the indices joined by `_`, after a `_` for `k > 1`
//!   (`symbol_array("x", 5)` is `[x0, …, x4]`, `symbol_array("a", 2, 3)`
//!   `[[a_0_0, …], …]`). The template is checked as giac's (letters,
//!   digits, `_`, `%`, no digit first) and its first character must be a
//!   letter or `_`, so that every name is an identifier (giac's
//!   `symbol_array("%a", 2)` is `[0a, 1a]`; ours stays). Dimensions are
//!   positive integers, as giac's ("Expected a sequence of positive
//!   integers"); their product is checked against `Limits::terms` before
//!   anything is built (C0019). A name that is a builtin's (`Log%`, 2)
//!   leaves the call. giac's trailing `purge` option is not taken.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

fn warn(ev: &mut Eval<'_>, message: &'static str, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// Whether `s` is an identifier: an ASCII letter or `_` first, then ASCII
/// letters, digits and `_`.
fn is_identifier(s: &str) -> bool {
    let mut chars = s.chars();
    chars
        .next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// `MakeSymbol[s]`.
pub(super) fn make_symbol(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [s] = e.args() else {
        return Ok(None);
    };
    let Some(name) = s.as_str() else {
        return Ok(None);
    };
    if !is_identifier(name) {
        warn(ev, "the string is not an identifier", e);
        return Ok(None);
    }
    let symbol = Symbol::new(name);
    if symbol.is_builtin() {
        warn(ev, "the string names a builtin", e);
        return Ok(None);
    }
    Ok(Some(Expr::from(symbol)))
}

/// The template with its `%` positions, a suffix added as giac does when
/// there is none; `None` for a template giac or the identifier rule
/// refuses.
fn template(s: &str, k: usize) -> Option<(String, Vec<usize>)> {
    let first = s.chars().next()?;
    if !(first.is_ascii_alphabetic() || first == '_')
        || !s
            .chars()
            .all(|c| c == '%' || c == '_' || c.is_ascii_alphanumeric())
    {
        return None;
    }
    let mut holes: Vec<usize> = s
        .char_indices()
        .filter(|&(_, c)| c == '%')
        .map(|(i, _)| i)
        .collect();
    if !holes.is_empty() {
        return (holes.len() == k).then(|| (s.to_owned(), holes));
    }
    let mut t = s.to_owned();
    if k > 1 {
        t.push('_');
    }
    for i in 0..k {
        if i > 0 {
            t.push('_');
        }
        holes.push(t.len());
        t.push('%');
    }
    Some((t, holes))
}

/// The name at `indices`: the template's holes, from the last, replaced by
/// the indices.
fn name(template: &str, holes: &[usize], indices: &[usize]) -> String {
    let mut name = template.to_owned();
    for (&at, i) in holes.iter().zip(indices).rev() {
        name.replace_range(at..=at, &i.to_string());
    }
    name
}

/// The nested lists from dimension `depth` on, `indices` filled before it;
/// `None` when a name is a builtin's.
fn array(
    template: &str,
    holes: &[usize],
    dims: &[usize],
    indices: &mut Vec<usize>,
) -> Option<Expr> {
    let depth = indices.len();
    let mut items = Vec::with_capacity(dims[depth]);
    for i in 0..dims[depth] {
        indices.push(i);
        let item = if depth + 1 == dims.len() {
            let symbol = Symbol::new(&name(template, holes, indices));
            (!symbol.is_builtin()).then(|| Expr::from(symbol))
        } else {
            array(template, holes, dims, indices)
        };
        indices.pop();
        items.push(item?);
    }
    Some(Expr::call(Symbol::LIST, items))
}

/// `SymbolArray[s, d₁, …, dₖ]`.
pub(super) fn symbol_array(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [s, dims @ ..] = e.args() else {
        return Ok(None);
    };
    let Some(s) = s.as_str() else {
        return Ok(None);
    };
    if dims.is_empty() {
        return Ok(None);
    }
    let Some(dims) = dims
        .iter()
        .map(|d| {
            d.as_i64()
                .filter(|&d| d > 0)
                .and_then(|d| usize::try_from(d).ok())
        })
        .collect::<Option<Vec<usize>>>()
    else {
        return Ok(None);
    };
    let limit = ev.context().limits().terms;
    let fits = dims
        .iter()
        .try_fold(1usize, |n, &d| n.checked_mul(d).filter(|&n| n <= limit))
        .is_some();
    if !fits {
        warn(ev, "the array exceeds the terms limit", e);
        return Ok(None);
    }
    let Some((t, holes)) = template(s, dims.len()) else {
        warn(ev, "an invalid name template", e);
        return Ok(None);
    };
    let Some(result) = array(&t, &holes, &dims, &mut Vec::with_capacity(dims.len())) else {
        warn(ev, "a name is a builtin's", e);
        return Ok(None);
    };
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
