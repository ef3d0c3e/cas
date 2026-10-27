//! Tests of `Sum` and `Product` (CAS-224): giac 1.9's finite sums and
//! products, Karr's convention, linearity and polynomial antidifferences.

use crate::{Context, Expr, Symbol};

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

/// `(input, expected)`: the same polynomial, their difference expanded 0.
fn check_equal(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let difference = Expr::call(Symbol::EXPAND, [read(input) - read(expected)]);
        let found = ctx.eval(&difference).unwrap();
        if found != Expr::from(0) {
            wrong.push(format!("{input}: differs by {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn finite_bounds() {
    check(&[
        ("Sum[k, k, 1, 10]", "55"),
        ("Sum[Power[k, 2], k, 1, 4]", "30"),
        ("Sum[k, Equation[k, Range[1, 10]]]", "55"),
        // Karr's convention, as giac's sum(k,k,5,1) is -9.
        ("Sum[k, k, 5, 1]", "-9"),
        ("Sum[k, k, 5, 4]", "0"),
        ("Product[k, k, 1, 5]", "120"),
        ("Product[k, k, 5, 1]", "1/24"),
        ("Product[k, k, 5, 4]", "1"),
        (
            "Sum[Power[x, k], k, 0, 3]",
            "Plus[1, x, Power[x, 2], Power[x, 3]]",
        ),
    ]);
}

#[test]
fn polynomial_sums() {
    check(&[
        ("Sum[k, k, 1, n]", "Times[1/2, n, Plus[1, n]]"),
        (
            "Sum[Power[k, 2], k, 1, n]",
            "Times[1/6, n, Plus[1, n], Plus[1, Times[2, n]]]",
        ),
        (
            "Sum[Power[k, 3], k, 1, n]",
            "Times[1/4, Power[n, 2], Power[Plus[1, n], 2]]",
        ),
        ("Sum[1, k, 1, n]", "n"),
        ("Sum[c, k, 1, n]", "Times[c, n]"),
        ("Sum[Times[2, k], k, 1, n]", "Times[n, Plus[1, n]]"),
        (
            "Sum[k, Equation[k, Range[1, n]]]",
            "Times[1/2, n, Plus[1, n]]",
        ),
    ]);
    check_equal(&[
        (
            "Sum[Power[k, 2], k, m, n]",
            "Plus[Times[1/6, n, Plus[n, 1], Plus[Times[2, n], 1]], \
             Times[-1/6, Plus[m, -1], m, Plus[Times[2, m], -1]]]",
        ),
        (
            "Sum[Plus[Times[a, k], b], k, 1, n]",
            "Plus[Times[1/2, a, n, Plus[n, 1]], Times[b, n]]",
        ),
        (
            "Sum[Power[k, 4], k, 1, n]",
            "Times[1/30, n, Plus[n, 1], Plus[Times[2, n], 1], Plus[Times[3, Power[n, 2]], Times[3, n], -1]]",
        ),
        ("Sum[k, k, 0, n]", "Times[1/2, n, Plus[n, 1]]"),
    ]);
}

/// What has no closed form here stays, alone or under one `Sum`.
#[test]
fn what_stays() {
    check(&[
        ("Sum[Factorial[k], k, 0, n]", "Sum[Factorial[k], k, 0, n]"),
        (
            "Sum[Plus[k, Factorial[k]], k, 0, n]",
            "Plus[Times[1/2, n, Plus[1, n]], Sum[Factorial[k], k, 0, n]]",
        ),
        ("Product[k, k, 1, n]", "Product[k, k, 1, n]"),
        ("Product[c, k, 1, n]", "Power[c, n]"),
    ]);
}

#[test]
fn binomials() {
    check(&[
        ("Binomial[5, 2]", "10"),
        ("Binomial[5, 7]", "0"),
        ("Binomial[5, -1]", "0"),
        ("Binomial[-3, 2]", "6"),
        ("Binomial[1/2, 2]", "-1/8"),
        ("Binomial[n, 2]", "Times[1/2, n, Plus[-1, n]]"),
        ("Binomial[n, k]", "Binomial[n, k]"),
    ]);
}

/// The closed form of `Sum[t, k, a, n]` at `n = 7` is the sum of its eight
/// terms or so, written out.
fn check_closed_forms(cases: &[(&str, i64)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for &(t, a) in cases {
        let closed = ctx.eval(&read(&format!("Sum[{t}, k, {a}, n]"))).unwrap();
        if closed.head() == Some(&Expr::from(Symbol::SUM)) {
            wrong.push(format!("{t}: no closed form"));
            continue;
        }
        let at = ctx
            .eval(&closed.substitute(&[(read("n"), Expr::from(7))]))
            .unwrap();
        let terms = ctx.eval(&read(&format!("Sum[{t}, k, {a}, 7]"))).unwrap();
        let difference = ctx
            .eval(&Expr::call(Symbol::TOGETHER, [at - terms]))
            .unwrap();
        if difference != Expr::from(0) {
            wrong.push(format!("{t}: {closed:?} is off by {difference:?} at 7"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Hypergeometric terms by Gosper's algorithm (CAS-225).
#[test]
fn hypergeometric_sums() {
    check(&[(
        "Sum[Power[2, k], k, 0, n]",
        "Plus[-1, Power[2, Plus[1, n]]]",
    )]);
    check_closed_forms(&[
        ("Power[2, k]", 0),
        ("Power[x, k]", 0),
        ("Times[k, Power[2, k]]", 0),
        ("Times[Power[k, 2], Power[2, k]]", 0),
        ("Power[Times[k, Plus[k, 1]], -1]", 1),
        ("Times[k, Power[Factorial[Plus[k, 1]], -1]]", 0),
        (
            "Power[Times[Plus[Times[2, k], -1], Plus[Times[2, k], 1]], -1]",
            1,
        ),
        ("Power[Times[k, Plus[k, 1], Plus[k, 2]], -1]", 1),
        ("Power[Plus[Power[k, 2], -1], -1]", 2),
        ("Times[Plus[k, 1], Power[2, k]]", 0),
        ("Times[k, Factorial[k]]", 0),
    ]);
    // No hypergeometric antidifference: it stays.
    check(&[
        ("Sum[Power[k, -1], k, 1, n]", "Sum[Power[k, -1], k, 1, n]"),
        ("Sum[Factorial[k], k, 0, n]", "Sum[Factorial[k], k, 0, n]"),
    ]);
}
