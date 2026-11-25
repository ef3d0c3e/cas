//! Tests of the language hook: printing, reading, what reading refuses, and
//! the round trip through a context.

use cas_core::{Context, Expr, Form, Language, Symbol};

use super::Hook;
use crate::names::Names;
use crate::print::{Balls, print_capped};
use crate::{PrintErrorKind, print_with};

fn context() -> Context {
    let mut context = Context::new();
    context.set_language(Some(Hook::english()));
    context
}

#[test]
fn prints_the_input_form() {
    let context = context();
    let hook = context.language().unwrap();
    let x = Expr::from(Symbol::new("x"));
    let e = Expr::apply(&Symbol::PLUS.into(), &[x, 1.into()]);
    assert_eq!(
        hook.print(&e, Form::Input, &context, usize::MAX).unwrap(),
        "x + 1"
    );
    assert!(hook.print(&e, Form::Input, &context, 4).is_err());
}

#[test]
fn reads_unevaluated() {
    let context = context();
    let e = Hook::english().read("1+1", &context).unwrap();
    assert_eq!(e, Expr::apply(&Symbol::PLUS.into(), &[1.into(), 1.into()]));
}

#[test]
fn refuses_what_is_not_one_expression() {
    let context = context();
    let hook = Hook::new(Names::english().clone());
    for text in ["a:=3", "assume(x>0)", "purge(x)", "x+", "1;2", ""] {
        assert!(hook.read(text, &context).is_err(), "{text}");
    }
}

#[test]
fn reads_back_what_it_prints() {
    let context = context();
    let hook = Hook::english();
    for text in [
        "x^2 + 2*x + 1",
        "[1, x, [2]]",
        "3/4",
        "\"h\\u{e9}\\\"llo\"",
        "0.1 + 0.2",
        "f(x, -y)",
    ] {
        let e = context.eval(&hook.read(text, &context).unwrap()).unwrap();
        let printed = hook.print(&e, Form::Input, &context, usize::MAX).unwrap();
        let again = context
            .eval(&hook.read(&printed, &context).unwrap())
            .unwrap();
        let equal = Expr::apply(&Symbol::EQUAL.into(), &[again.clone(), e.clone()]);
        assert!(
            again == e || context.eval(&equal).unwrap().as_symbol() == Some(Symbol::TRUE),
            "{text} printed {printed}"
        );
    }
}

#[test]
fn a_capped_print_refuses_longer_text() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let e = Expr::apply(
        &Symbol::LIST.into(),
        &(0..100).map(Expr::from).collect::<Vec<_>>(),
    );
    let text = print_with(&e, names, balls).unwrap();
    assert_eq!(print_capped(&e, names, balls, text.len()).unwrap(), text);
    let error = print_capped(&e, names, balls, text.len() - 1).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooLong);
    assert_eq!(error.part, e);
    let error = print_capped(&e, names, balls, 50).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooLong);
}
