//! Tests of the bits, digits and falling factorials (CAS-464): giac's
//! answers, negative operands, and what stays.

use crate::{Context, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
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
fn bitwise() {
    check(&[
        // giac's 8, 14, 6 and 0; with -12, giac's (and Python's) 0, -2, -2.
        ("BitAnd[12, 10]", "8"),
        ("BitOr[12, 10]", "14"),
        ("BitXor[12, 10]", "6"),
        ("BitAnd[12, 10, 6]", "0"),
        ("BitAnd[-12, 10]", "0"),
        ("BitOr[-12, 10]", "-2"),
        ("BitXor[-12, 10]", "-2"),
        // Python's -12 & -10 = -12, -12 | -10 = -10, -12 ^ -10 = 2.
        ("BitAnd[-12, -10]", "-12"),
        ("BitOr[-12, -10]", "-10"),
        ("BitXor[-12, -10]", "2"),
        // Past a word: 2^70 + 5 and 2^70 + 3.
        (
            "BitAnd[Plus[Power[2, 70], 5], Plus[Power[2, 70], 3]]",
            "Plus[Power[2, 70], 1]",
        ),
        // giac's 2 and 4; a negative stays.
        ("HammingDistance[5, 3]", "2"),
        ("HammingDistance[15, 0]", "4"),
        ("HammingDistance[-1, 0]", "HammingDistance[-1, 0]"),
        ("BitAnd[x, 1]", "BitAnd[x, 1]"),
    ]);
}

#[test]
fn digits_and_falling_factorials() {
    check(&[
        // giac's [15,15], [-15,-15], [0], [0,1,0,1] and 255, 5.
        ("Convert[255, base, 16]", "List[15, 15]"),
        ("Convert[-255, base, 16]", "List[-15, -15]"),
        ("Convert[0, base, 2]", "List[0]"),
        ("Convert[10, base, 2]", "List[0, 1, 0, 1]"),
        ("Convert[List[15, 15], base, 16]", "255"),
        ("Convert[List[1, 0, 1], base, 2]", "5"),
        ("Convert[10, base, 1]", "Convert[10, base, 1]"),
        // giac's 720, 20, 0 and n!/((n-3)!) as the product.
        ("FactorialPower[10, 3]", "720"),
        ("FactorialPower[5, 2]", "20"),
        ("FactorialPower[3, 5]", "0"),
        ("FactorialPower[n, 0]", "1"),
        ("FactorialPower[n, 3]", "Times[n, Plus[n, -1], Plus[n, -2]]"),
        ("FactorialPower[1/2, 2]", "-1/4"),
    ]);
}
