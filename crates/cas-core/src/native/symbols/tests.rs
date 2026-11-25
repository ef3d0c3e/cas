//! Tests of `MakeSymbol` and `SymbolArray` (CAS-623).

use crate::{Context, Expr, Symbol};

/// `(input, expected)`: the input evaluates to the expected expression.
#[test]
fn symbols() {
    let ctx = Context::new();
    let s = |text: &str| Expr::string(text);
    let int = |n: i64| Expr::from(n);
    let call = |head: Symbol, args: Vec<Expr>| Expr::call(head, args);
    let list = |items: &[&str]| {
        Expr::call(
            Symbol::LIST,
            items.iter().map(|n| Expr::from(Symbol::new(n))),
        )
    };
    let cases = [
        (
            call(Symbol::MAKE_SYMBOL, vec![s("abc")]),
            Expr::from(Symbol::new("abc")),
        ),
        (
            call(Symbol::MAKE_SYMBOL, vec![s("_x1")]),
            Expr::from(Symbol::new("_x1")),
        ),
        (
            call(Symbol::SYMBOL_ARRAY, vec![s("x"), int(5)]),
            list(&["x0", "x1", "x2", "x3", "x4"]),
        ),
        (
            call(Symbol::SYMBOL_ARRAY, vec![s("a"), int(2), int(3)]),
            Expr::call(
                Symbol::LIST,
                [
                    list(&["a_0_0", "a_0_1", "a_0_2"]),
                    list(&["a_1_0", "a_1_1", "a_1_2"]),
                ],
            ),
        ),
        (
            call(Symbol::SYMBOL_ARRAY, vec![s("a%%"), int(2), int(3)]),
            Expr::call(
                Symbol::LIST,
                [list(&["a00", "a01", "a02"]), list(&["a10", "a11", "a12"])],
            ),
        ),
        (
            call(Symbol::SYMBOL_ARRAY, vec![s("b%_c%"), int(1), int(2)]),
            Expr::call(Symbol::LIST, [list(&["b0_c0", "b0_c1"])]),
        ),
    ];
    for (input, expected) in cases {
        assert_eq!(ctx.eval(&input).unwrap(), expected, "{input:?}");
    }
}

/// What stays: strings that are not identifiers, builtins' names, bad
/// templates and dimensions, and arrays past `Limits::terms`.
#[test]
fn stays() {
    let ctx = Context::new();
    let s = |text: &str| Expr::string(text);
    let int = |n: i64| Expr::from(n);
    let inputs = [
        Expr::call(Symbol::MAKE_SYMBOL, [s("x y")]),
        Expr::call(Symbol::MAKE_SYMBOL, [s("1a")]),
        Expr::call(Symbol::MAKE_SYMBOL, [s("")]),
        Expr::call(Symbol::MAKE_SYMBOL, [s("Plus")]),
        Expr::call(Symbol::MAKE_SYMBOL, [int(1)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("x")]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("x"), int(0)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("x"), int(-1)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("%a"), int(2)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("a%"), int(2), int(3)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("a b"), int(2)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("Log%"), int(3)]),
        Expr::call(Symbol::SYMBOL_ARRAY, [s("x"), int(1001), int(1000)]),
    ];
    for input in inputs {
        assert_eq!(ctx.eval(&input).unwrap(), input, "{input:?}");
    }
}
