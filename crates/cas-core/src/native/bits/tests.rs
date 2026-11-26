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

#[test]
fn convert_options() {
    check(&[
        // giac's matrix[[1,2],[3,4]], matrix[1,2,3], x+1 and 5 for array (CAS-639).
        (
            "Convert[List[List[1, 2], List[3, 4]], array]",
            "List[List[1, 2], List[3, 4]]",
        ),
        ("Convert[List[1, 2, 3], array]", "List[1, 2, 3]"),
        ("Convert[Plus[x, 1], array]", "Plus[x, 1]"),
        ("Convert[5, array]", "5"),
        // giac's `expln` is its `trig2exp` (CAS-640).
        ("Convert[Sin[x], expln]", "TrigToExp[Sin[x]]"),
        // giac's 1/2/(x-1)-1/2/(x+1) for partfrac, parfrac and fullparfrac.
        (
            "Convert[Power[Plus[Power[x, 2], -1], -1], parfrac]",
            "Apart[Power[Plus[Power[x, 2], -1], -1]]",
        ),
        (
            "Convert[Power[Plus[Power[x, 2], -1], -1], fullparfrac]",
            "Plus[Times[1/2, Power[Plus[x, -1], -1]], Times[-1/2, Power[Plus[x, 1], -1]]]",
        ),
        // `partfrac` reads as `Apart`.
        (
            "Convert[Power[Plus[Power[x, 2], -1], -1], Apart]",
            "Apart[Power[Plus[Power[x, 2], -1], -1]]",
        ),
        // giac's [3,2], [0,3]; dfc's [3,7,16]; the third argument ignored.
        ("Convert[7/2, confrac]", "List[3, 2]"),
        ("Convert[1/3, confrac]", "List[0, 3]"),
        ("Convert[355/113, confrac, L]", "List[3, 7, 16]"),
        // giac's x-1/6*x^3+1/120*x^5, and x^2+1 itself.
        (
            "Convert[Series[Sin[x], List[x, 0, 5]], polynom]",
            "Plus[x, Times[-1/6, Power[x, 3]], Times[1/120, Power[x, 5]]]",
        ),
        (
            "Convert[Plus[Power[x, 2], 1], polynom]",
            "Plus[Power[x, 2], 1]",
        ),
        ("Convert[x, other]", "Convert[x, other]"),
        // giac's interval as a ball (CAS-592): `N`, with digits too.
        ("Convert[Power[2, 1/2], Interval]", "N[Power[2, 1/2]]"),
        (
            "Convert[Power[2, 1/2], Interval, 20]",
            "N[Power[2, 1/2], 20]",
        ),
        ("Convert[List[1/3, Pi], Interval]", "N[List[1/3, Pi]]"),
        ("Convert[Plus[x, 1], Interval]", "N[Plus[x, 1]]"),
        ("Convert[2, Interval, 0]", "Convert[2, Interval, 0]"),
    ]);
}
