//! Tests of the discrete Fourier transform (CAS-465): giac's answers within
//! 10⁻¹², lengths 1 to 8, the inverse, modulo a prime, and convolutions.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`: lists of one length whose entries are within
/// 10⁻¹² of the exact ones.
fn check_near(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        // The exact values as balls: `N` of a difference near 0 has no
        // relative accuracy to reach.
        let expected = ctx
            .eval(&Expr::call(Symbol::N, [read(expected), Expr::from(40)]))
            .unwrap();
        let near = found.args().len() == expected.args().len()
            && found.args().iter().zip(expected.args()).all(|(f, x)| {
                let gap = Expr::call(
                    Symbol::ABS,
                    [Expr::call(
                        Symbol::PLUS,
                        [f.clone(), Expr::from(-1) * x.clone()],
                    )],
                );
                let small = Expr::call(Symbol::LESS, [gap, read("Power[10, -12]")]);
                ctx.eval(&small).ok().and_then(|t| t.as_symbol()) == Some(Symbol::TRUE)
            });
        if !near || found.head() != Some(&Expr::from(Symbol::LIST)) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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
fn numeric() {
    check_near(&[
        // giac's [10.0,-2.0+2.0*i,-2.0,-2.0-2.0*i] and back.
        (
            "DiscreteFourierTransform[List[1, 2, 3, 4]]",
            "List[10, Complex[-2, 2], -2, Complex[-2, -2]]",
        ),
        (
            "InverseDiscreteFourierTransform[List[10, Complex[-2, 2], -2, Complex[-2, -2]]]",
            "List[1, 2, 3, 4]",
        ),
        // giac's [6.0,-1.5+0.866025403784*i,-1.5-0.866025403784*i]: −3/2 ± (√3/2) i.
        (
            "DiscreteFourierTransform[List[1, 2, 3]]",
            "List[6, Plus[-3/2, Times[1/2, Power[3, 1/2], I]], \
             Plus[-3/2, Times[-1/2, Power[3, 1/2], I]]]",
        ),
        // giac's [1.5,-0.5]; one entry; a constant's transform.
        (
            "InverseDiscreteFourierTransform[List[1, 2]]",
            "List[3/2, -1/2]",
        ),
        ("DiscreteFourierTransform[List[5]]", "List[5]"),
        (
            "DiscreteFourierTransform[List[1, 1, 1, 1, 1]]",
            "List[5, 0, 0, 0, 0]",
        ),
        // A delta at 1 of length 8: ω^k = e^(−πik/4).
        (
            "DiscreteFourierTransform[List[0, 1, 0, 0, 0, 0, 0, 0]]",
            "Table[Power[E, Times[-1/4, Pi, I, k]], List[k, 0, 7]]",
        ),
    ]);
    check(&[(
        "DiscreteFourierTransform[List[1, x]]",
        "DiscreteFourierTransform[List[1, x]]",
    )]);
}

#[test]
fn modular_and_convolutions() {
    check(&[
        // giac's fft([1,2,3,4],2,5): [0,2,3,4], and back.
        (
            "DiscreteFourierTransform[List[1, 2, 3, 4], 2, 5]",
            "List[0, 2, 3, 4]",
        ),
        (
            "InverseDiscreteFourierTransform[List[0, 2, 3, 4], 2, 5]",
            "List[1, 2, 3, 4]",
        ),
        // 4 has order 2 modulo 5, not 4.
        (
            "DiscreteFourierTransform[List[1, 2, 3, 4], 4, 5]",
            "DiscreteFourierTransform[List[1, 2, 3, 4], 4, 5]",
        ),
        // giac's [1.0,1.0,1.0,-3.0], exactly; symbols.
        (
            "Convolution[List[1, 2, 3], List[1, -1]]",
            "List[1, 1, 1, -3]",
        ),
        (
            "Convolution[List[x, 1], List[1, y]]",
            "List[x, Plus[1, Times[x, y]], y]",
        ),
    ]);
}
