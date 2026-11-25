//! `Type` (`type`, CAS-621), Xcas's, and its answers `DOM_INT`, `DOM_RAT`,
//! `DOM_FLOAT`, `DOM_COMPLEX`, `DOM_IDENT`, `DOM_LIST`, `DOM_SYMBOLIC`,
//! `DOM_STRING` and `DOM_FUNC`, symbols without a value (giac prog.cc
//! `_type`, which answers an integer printed as its type's name). Their
//! names here are `DomInt` … `DomFunc` (a builtin's name is an ASCII
//! capital, then letters and digits); the language spells them `DOM_INT`
//! … `DOM_FUNC` (`names/en.txt`).
//!
//! The kind of the evaluated argument: an integer of any size `DOM_INT`
//! (giac's small and big integers both `_ZINT`), a rational `DOM_RAT`, a
//! decimal ball `DOM_FLOAT`, a non-real number `DOM_COMPLEX`, a list
//! `DOM_LIST`, a string `DOM_STRING`, a `Function` or a builtin with a native
//! (`sin`, `sqrt`) `DOM_FUNC`, any other symbol `DOM_IDENT` (`pi`, `undef`,
//! `infinity` as giac; `true` too, where giac's is an integer), anything else
//! `DOM_SYMBOLIC`. Another number of arguments stays.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// `Type[x]`.
pub(super) fn type_of(_ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let domain = match x.view() {
        View::Int(_) | View::Integer(_) => Symbol::DOM_INT,
        View::Rational(_) => Symbol::DOM_RAT,
        View::Ball(_) => Symbol::DOM_FLOAT,
        View::Complex(_) | View::ComplexBall(_) => Symbol::DOM_COMPLEX,
        View::String(_) => Symbol::DOM_STRING,
        View::Symbol(s) if super::get(s).is_some() => Symbol::DOM_FUNC,
        View::Symbol(_) => Symbol::DOM_IDENT,
        View::Apply { head, .. } => match head.as_symbol() {
            Some(Symbol::LIST) => Symbol::DOM_LIST,
            Some(Symbol::FUNCTION) => Symbol::DOM_FUNC,
            _ => Symbol::DOM_SYMBOLIC,
        },
    };
    Ok(Some(Expr::from(domain)))
}
