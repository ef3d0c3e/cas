//! Expressions as text and text as expressions (CAS-614): `ToString`
//! (giac's `string`, `str`), `ToExpression` (`expr`, `execute`) and `Cat`
//! (`cat`), and the printing of non-strings that `StringJoin` (`join`) and
//! `Join` (`concat`) share with `Cat`. Ported from giac 2.0.0.21: prog.cc
//! `_string` (:10930) and `_expr` (:10892-10924), maple.cc `_cat`
//! (:946-958), permu.cc `_join` (:1803-1824). They print and read through
//! the context's [`Language`](crate::Language) hook in its input form; with
//! no hook, what needs one stays.
//!
//! `ToString[e]` is `e` printed. `ToExpression[s]` reads the string `s` as
//! one expression and evaluates it; a blank `s` is `Indeterminate` (giac's
//! `undef`). `ToExpression[f, args…]`, `f` not a string, applies `f` to
//! `args` when `f` is a builtin, a `Function`, or a symbol with a value or
//! rules of its own, and stays otherwise, as giac's `expr(f,x,y)` does.
//! `Cat[a…]` concatenates its arguments, strings as they are and the others
//! printed; `Cat[]` is `""`.
//!
//! Departures from giac:
//! - The text is this CAS's input form, not giac's: `string(x+1)` is
//!   `"x + 1"` here, `"x+1"` there, so `size(string(x+1))` is 5, not 3.
//! - `expr` reads with this CAS's grammar: `expr("0123")` is 123 (giac's
//!   octal 83); `expr("a:=3")` (giac assigns, 3) and `expr("1;2")` (giac
//!   `1,2`) stay with the reader's message as a warning, as does a parse
//!   error (`expr("x+")`, giac's error).
//! - Text longer than the `terms` limit (characters) stays: a string to
//!   read, with a warning; a print, a `cat`, a `join`.

use std::borrow::Cow;

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::language::Form;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// `ToString[e]`: `e` in the input form.
pub(super) fn to_string(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let context = ev.context();
    let Some(language) = context.language() else {
        return Ok(None);
    };
    match language.print(x, Form::Input, context, context.limits().terms) {
        Ok(text) => Ok(Some(Expr::string(&text))),
        Err(message) => {
            warn(ev, message, e);
            Ok(None)
        }
    }
}

/// `ToExpression[s]`: the string `s` read and evaluated;
/// `ToExpression[f, args…]`: `f[args…]` evaluated, for a function `f`.
pub(super) fn to_expression(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if let [s] = e.args()
        && let Some(text) = s.as_str()
    {
        return read(ev, text, e);
    }
    match e.args() {
        [f, args @ ..] if !args.is_empty() && f.as_str().is_none() && applies(ev, f) => {
            Ok(Some(ev.eval(&Expr::apply(f, args))?))
        }
        _ => Ok(None),
    }
}

/// The text `s` read through the hook and evaluated.
fn read(ev: &mut Eval<'_>, s: &str, e: &Expr) -> Native {
    if s.trim().is_empty() {
        return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
    }
    let context = ev.context();
    let Some(language) = context.language() else {
        return Ok(None);
    };
    if s.chars().count() > context.limits().terms {
        warn(ev, "the text is longer than the terms limit", e);
        return Ok(None);
    }
    match language.read(s, context) {
        Ok(x) => Ok(Some(ev.eval(&x)?)),
        Err(message) => {
            warn(ev, message, e);
            Ok(None)
        }
    }
}

/// Whether `ToExpression` applies `f`: a builtin, a `Function`, or a
/// symbol with a value or rules of its own.
fn applies(ev: &Eval<'_>, f: &Expr) -> bool {
    if f.head() == Some(&Expr::from(Symbol::FUNCTION)) {
        return true;
    }
    f.as_symbol().is_some_and(|s| {
        s.is_builtin() || ev.context().value(s).is_some() || ev.context().rules().has_down_rules(s)
    })
}

/// `Cat[a…]`: the arguments as text, one after the other.
pub(super) fn cat(ev: &mut Eval<'_>, e: &Expr) -> Native {
    Ok(concatenated(ev, e.args(), "").map(|text| Expr::string(&text)))
}

/// `items` as text with `sep` between them: strings as they are, the others
/// in the input form. `None` when a non-string meets no hook or cannot be
/// printed, or when the text would be longer than the `terms` limit in
/// characters (each print capped at what is left).
pub(super) fn concatenated(ev: &Eval<'_>, items: &[Expr], sep: &str) -> Option<String> {
    let context = ev.context();
    let terms = context.limits().terms;
    let sep_chars = sep.chars().count();
    let mut text = String::new();
    let mut length = 0usize;
    for (k, item) in items.iter().enumerate() {
        if k > 0 {
            length = length.checked_add(sep_chars).filter(|&n| n <= terms)?;
            text.push_str(sep);
        }
        let printed;
        let piece = if let Some(s) = item.as_str() {
            s
        } else {
            printed = context
                .language()?
                .print(item, Form::Input, context, terms - length)
                .ok()?;
            printed.as_str()
        };
        length = length
            .checked_add(piece.chars().count())
            .filter(|&n| n <= terms)?;
        text.push_str(piece);
    }
    Some(text)
}

fn warn(ev: &mut Eval<'_>, message: impl Into<Cow<'static, str>>, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

#[cfg(test)]
mod tests;
