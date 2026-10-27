//! Tests of `Re`, `Im`, `Conjugate` and `Arg` (CAS-240).

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
fn numbers() {
    check(
        &Context::new(),
        &[
            ("Re[Plus[3, Times[4, I]]]", "3"),
            ("Im[Plus[3, Times[4, I]]]", "4"),
            ("Conjugate[Plus[3, Times[4, I]]]", "Plus[3, Times[-4, I]]"),
            ("Arg[Plus[1, I]]", "Times[1/4, Pi]"),
            ("Arg[-1]", "Pi"),
            ("Arg[I]", "Times[1/2, Pi]"),
            ("Arg[Times[-2, I]]", "Times[-1/2, Pi]"),
            ("Arg[0]", "0"),
            ("Arg[Plus[1, Times[I, Power[3, 1/2]]]]", "Times[1/3, Pi]"),
            ("Arg[Plus[-3, Times[4, I]]]", "Plus[Pi, ArcTan[-4/3]]"),
            ("Re[Power[-2, 1/2]]", "0"),
            ("Im[Power[-2, 1/2]]", "Power[2, 1/2]"),
        ],
    );
}

/// In the Real domain a symbol is real.
#[test]
fn expressions() {
    check(
        &Context::new(),
        &[
            ("Re[x]", "x"),
            ("Im[x]", "0"),
            ("Conjugate[x]", "x"),
            ("Re[Plus[x, Times[I, y]]]", "x"),
            ("Im[Plus[x, Times[I, y]]]", "y"),
            (
                "Conjugate[Plus[x, Times[I, y]]]",
                "Plus[x, Times[-1, I, y]]",
            ),
            ("Re[Exp[Times[I, x]]]", "Cos[x]"),
            ("Im[Exp[Times[I, x]]]", "Sin[x]"),
            (
                "Re[Power[Plus[x, Times[I, y]], -1]]",
                "Times[x, Power[Plus[Power[x, 2], Power[y, 2]], -1]]",
            ),
            ("Im[Power[x, 2]]", "0"),
            ("Re[Sin[x]]", "Sin[x]"),
            // Not negative: 0, at 0 too.
            ("Arg[Power[x, 2]]", "0"),
            ("Arg[Sin[x]]", "Arg[Sin[x]]"),
            // An unknown function stays.
            ("Re[f[x]]", "Re[f[x]]"),
        ],
    );
    let positive = Context::new().assuming(read("Greater[x, 0]"));
    check(
        &positive,
        &[
            ("Arg[x]", "0"),
            ("Arg[Times[-1, x]]", "Pi"),
            ("Arg[Plus[x, Times[I, x]]]", "Times[1/4, Pi]"),
        ],
    );
}
