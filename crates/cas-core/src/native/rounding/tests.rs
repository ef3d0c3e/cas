//! Tests of `Floor`, `Ceiling`, `Round`, `Sign`, `Max` and `Min` (CAS-241).

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(ctx: &Context, cases: &[(&str, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn rounding() {
    check(
        &Context::new(),
        &[
            ("Floor[7/2]", "3"),
            ("Ceiling[7/2]", "4"),
            ("Ceiling[-7/2]", "-3"),
            ("Floor[-1/3]", "-1"),
            ("Round[7/2]", "4"),
            ("Round[5/2]", "3"),
            ("Round[-5/2]", "-2"),
            ("Round[-7/2]", "-3"),
            ("Round[-1/2]", "0"),
            ("Floor[Power[2, 1/2]]", "1"),
            ("Ceiling[Pi]", "4"),
            (
                "Floor[Plus[1000000000000000000000000000000, 1/2]]",
                "1000000000000000000000000000000",
            ),
            ("Floor[Infinity]", "Infinity"),
            ("Floor[x]", "Floor[x]"),
            ("Floor[Plus[3, x]]", "Plus[3, Floor[x]]"),
            ("Floor[List[1/2, 3/2]]", "List[0, 1]"),
        ],
    );
    let integer = Context::new().assuming(read("Element[n, Integers]"));
    check(
        &integer,
        &[("Floor[n]", "n"), ("Round[Plus[n, 1/2]]", "Plus[1, n]")],
    );
}

#[test]
fn signs() {
    check(
        &Context::new(),
        &[
            ("Sign[-3]", "-1"),
            ("Sign[0]", "0"),
            ("Sign[5/2]", "1"),
            ("Sign[Plus[Pi, -3]]", "1"),
            ("Sign[Plus[3, Times[4, I]]]", "Plus[3/5, Times[4/5, I]]"),
            ("Sign[x]", "Sign[x]"),
            ("Sign[Times[-2, x]]", "Times[-1, Sign[x]]"),
        ],
    );
    let positive = Context::new().assuming(read("Greater[x, 0]"));
    check(&positive, &[("Sign[x]", "1"), ("Sign[Times[-1, x]]", "-1")]);
}

#[test]
fn extrema() {
    check(
        &Context::new(),
        &[
            ("Max[1, 2, 3]", "3"),
            ("Min[1, 2, 3]", "1"),
            ("Max[2, Power[5, 1/2]]", "Power[5, 1/2]"),
            ("Max[List[1, 5, 3]]", "5"),
            ("Min[-1/2, -1/3]", "-1/2"),
            ("Min[Power[2, 1/2], 3/2]", "Power[2, 1/2]"),
            ("Max[Pi, 22/7]", "22/7"),
            ("Max[1]", "1"),
            ("Max[x, 1]", "Max[1, x]"),
            ("Max[-1, x, 3]", "Max[3, x]"),
            ("Max[]", "Times[-1, Infinity]"),
            ("Min[]", "Infinity"),
        ],
    );
    let large = Context::new().assuming(read("Greater[x, 2]"));
    check(&large, &[("Max[x, 1]", "x"), ("Min[x, 1]", "1")]);
}

/// Values whose balls are too long to write exactly, `2` to a huge power
/// (the `lang_session` fuzz target): they stay rather than panic.
#[test]
fn huge_balls() {
    let ctx = Context::new();
    for input in [
        "Floor[Power[2, 3777777777777]]",
        "Sign[Plus[-7, Power[2, 3777777777777]]]",
        "Max[7, Power[2, 3777777777777]]",
        "Greater[Power[2, 3777777777777], 7]",
        "Round[Power[Ball[3, 0], 377777777]]",
    ] {
        assert!(ctx.eval(&read(input)).is_ok(), "{input}");
    }
}
