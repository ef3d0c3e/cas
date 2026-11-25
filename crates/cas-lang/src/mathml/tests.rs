//! Tests of giac's `mathml` port: the 301 rows of `giac_mathml_doc.tsv`
//! (1.6/T5.c3), each giac's printed answer read back into giac's shape by
//! `markup`'s test parser and written as giac's `mathml()` did; the inline
//! form; the limits.

use cas_core::{Expr, Symbol};

use super::{document, mathml};
use crate::markup::tests::parse;
use crate::markup::tree::Gen;
use crate::names::Names;
use crate::print::{Balls, MAX_PRINT_DEPTH, PrintErrorKind, power_order};

/// giac's parser builds a quoted `a*b*c` as `prod(prod(a, b), c)`, which
/// `aplatir_fois_plus` then multiplies out (C0069; the prototype's `nest`).
fn nest(g: Gen) -> Gen {
    match g {
        Gen::Vect(items, sub) => Gen::Vect(items.into_iter().map(nest).collect(), sub),
        Gen::Symb(op, arg) => match nest(*arg) {
            Gen::Vect(items, _) if op == "prod" && items.len() > 2 => {
                let mut items = items.into_iter();
                let first = items.next().expect("three factors");
                items.fold(first, |acc, x| Gen::symb_n("prod", vec![acc, x]))
            }
            a => Gen::symb(&op, a),
        },
        g => g,
    }
}

#[test]
fn giacs_answers_are_written_as_giacs_mathml_writes_them() {
    let mut failures = Vec::new();
    let mut rows = 0;
    for line in include_str!("../markup/giac_mathml_doc.tsv").lines() {
        if line.starts_with('#') {
            continue;
        }
        let [input, mode, shape, want] = line.split('\t').collect::<Vec<_>>()[..] else {
            panic!("a row of four columns: {line}");
        };
        rows += 1;
        // giac's doubles printed without a point (the prototype's check.py).
        let g = match input {
            "1e10" => Gen::Dbl {
                text: "10000000000".to_owned(),
                neg: false,
            },
            "12345678901.0" => Gen::Dbl {
                text: "12345678901".to_owned(),
                neg: false,
            },
            _ if mode == "quote" => nest(parse(shape, true)),
            _ => parse(shape, false),
        };
        let got = document(&g, usize::MAX)
            .unwrap_or_else(|_| panic!("{input}: no cap"))
            // giac's table shows a newline as ' | ' (C0047).
            .replace('\n', " | ");
        if got != want {
            failures.push(format!("{input} ({shape}): {got:?}, giac {want:?}"));
        }
    }
    assert_eq!(rows, 301);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn a_list_with_a_machine_integer_is_the_formula_inline() {
    let got = document(&parse("[x,-1]", false), usize::MAX).unwrap_or_else(|_| panic!("no cap"));
    assert_eq!(
        got,
        "\n<math display=\"inline\" xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\n<mi>x</mi>\n\n</math>\n"
    );
}

#[test]
fn too_deep_is_refused_before_any_work() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let list = |e: Expr| Expr::apply(&Symbol::LIST.into(), &[e]);
    let mut e = Expr::from(1);
    // An atom is 1 deep.
    for _ in 1..MAX_PRINT_DEPTH {
        e = list(e);
    }
    assert!(mathml(&e, names, balls, usize::MAX).is_ok());
    let error = mathml(&list(e), names, balls, usize::MAX).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooDeep);
}

#[test]
fn text_past_the_cap_is_refused() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let e = Expr::apply(
        &Symbol::LIST.into(),
        &(0..100).map(Expr::from).collect::<Vec<_>>(),
    );
    let text = mathml(&e, names, balls, usize::MAX).unwrap();
    assert_eq!(mathml(&e, names, balls, text.len()).unwrap(), text);
    let error = mathml(&e, names, balls, text.len() - 1).unwrap_err();
    assert_eq!((error.kind, error.part), (PrintErrorKind::TooLong, e));
    // A long atom deep inside is refused as its text is made.
    let long = Expr::string(&"a".repeat(1000));
    let e = Expr::apply(&Symbol::PLUS.into(), &[long, Expr::from(Symbol::new("x"))]);
    let error = mathml(&e, names, balls, 100).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooLong);
}

/// `increasing_power` off (CAS-618): the sum by decreasing degree, as
/// giac's `mathml(expand((x+1)^3))` with the flag off (1.6/T6.c3 probe).
#[test]
fn a_sum_of_one_symbol_by_decreasing_degree() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let x = Expr::from(Symbol::new("x"));
    let apply = |head: Symbol, args: &[Expr]| Expr::apply(&head.into(), args);
    let cubic = apply(
        Symbol::PLUS,
        &[
            Expr::from(1),
            apply(Symbol::TIMES, &[Expr::from(3), x.clone()]),
            apply(Symbol::TIMES, &[Expr::from(3), x.pow(2)]),
            x.pow(3),
        ],
    );
    let text = mathml(&power_order(&cubic, false), names, balls, usize::MAX).unwrap();
    let giac = "<msup><mrow> <mi>x</mi></mrow><mrow><mn>3</mn></mrow></msup><mo>+</mo>\
                <mn>3</mn><msup><mrow> <mi>x</mi></mrow><mrow><mn>2</mn></mrow></msup>\
                <mo>+</mo><mn>3</mn><mi>x</mi><mo>+</mo><mn>1</mn>";
    assert!(text.contains(giac), "{text}");
    let text = mathml(&power_order(&cubic, true), names, balls, usize::MAX).unwrap();
    assert!(
        text.contains("<mn>1</mn><mo>+</mo><mn>3</mn><mi>x</mi>"),
        "{text}"
    );
}
