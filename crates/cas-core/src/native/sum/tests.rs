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

/// `(input, expected)` in full form, both evaluated: the same expression,
/// or the same number by `N`.
fn check_values(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let n = |e: &Expr| ctx.eval(&Expr::call(Symbol::N, [e.clone()])).unwrap();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found == expected {
            continue;
        }
        let equal = Expr::call(Symbol::EQUAL, [n(&found), n(&expected)]);
        if ctx.eval(&equal).unwrap().as_symbol() != Some(Symbol::TRUE) {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input evaluates to a `Sum`: it stays.
fn check_stays(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let found = ctx.eval(&read(input)).unwrap();
        if found.head() != Some(&Expr::from(Symbol::SUM)) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Infinite sums (CAS-226): limits of antidifferences, known series,
/// divergence, and what stays.
#[test]
fn infinite_sums() {
    check_values(&[
        // Gosper's antidifferences, to the limit.
        ("Sum[Power[2, Times[-1, k]], k, 0, Infinity]", "2"),
        ("Sum[Power[Times[k, Plus[k, 1]], -1], k, 1, Infinity]", "1"),
        (
            "Sum[Power[Plus[Power[k, 2], -1], -1], k, 2, Infinity]",
            "3/4",
        ),
        ("Sum[Times[k, Power[2, Times[-1, k]]], k, 1, Infinity]", "2"),
        (
            "Sum[Times[Power[k, 2], Power[2, Times[-1, k]]], k, 0, Infinity]",
            "6",
        ),
        (
            "Sum[Power[Times[k, Plus[k, 2]], -1], k, 1, Infinity]",
            "3/4",
        ),
        // Geometric, exponential and logarithmic series.
        ("Sum[Power[-1/2, k], k, 0, Infinity]", "2/3"),
        ("Sum[Power[Factorial[k], -1], k, 0, Infinity]", "E"),
        (
            "Sum[Power[Factorial[k], -1], k, 2, Infinity]",
            "Plus[-2, E]",
        ),
        (
            "Sum[Times[Power[3, k], Power[Factorial[k], -1]], k, 0, Infinity]",
            "Power[E, 3]",
        ),
        (
            "Sum[Times[Power[-1, k], Power[k, -1]], k, 1, Infinity]",
            "Times[-1, Log[2]]",
        ),
        (
            "Sum[Times[Power[2, Times[-1, k]], Power[k, -1]], k, 1, Infinity]",
            "Log[2]",
        ),
        // ζ(2m), and its alternating series.
        (
            "Sum[Power[k, -2], k, 1, Infinity]",
            "Times[1/6, Power[Pi, 2]]",
        ),
        (
            "Sum[Power[k, -4], k, 1, Infinity]",
            "Times[1/90, Power[Pi, 4]]",
        ),
        (
            "Sum[Power[k, -6], k, 1, Infinity]",
            "Times[1/945, Power[Pi, 6]]",
        ),
        (
            "Sum[Power[k, -2], k, 2, Infinity]",
            "Plus[-1, Times[1/6, Power[Pi, 2]]]",
        ),
        (
            "Sum[Times[Power[-1, k], Power[k, -2]], k, 1, Infinity]",
            "Times[-1/12, Power[Pi, 2]]",
        ),
        // Term by term, and from -Infinity.
        (
            "Sum[Plus[Power[2, Times[-1, k]], Power[k, -2]], k, 1, Infinity]",
            "Plus[1, Times[1/6, Power[Pi, 2]]]",
        ),
        (
            "Sum[Power[k, -2], k, Times[-1, Infinity], -1]",
            "Times[1/6, Power[Pi, 2]]",
        ),
        // Divergent.
        ("Sum[1, k, 1, Infinity]", "Infinity"),
        ("Sum[-1, k, 1, Infinity]", "Times[-1, Infinity]"),
        ("Sum[k, k, 1, Infinity]", "Infinity"),
        ("Sum[Power[k, -1], k, 1, Infinity]", "Infinity"),
        ("Sum[Power[2, k], k, 0, Infinity]", "Infinity"),
        (
            "Sum[Plus[Power[k, -1], Power[k, -2]], k, 1, Infinity]",
            "Infinity",
        ),
    ]);
    // A symbolic start: 2^(1 − m), 1/16 at m = 5.
    let ctx = Context::new();
    let found = ctx
        .eval(&read("Sum[Power[1/2, k], k, m, Infinity]"))
        .unwrap();
    let at = ctx
        .eval(&found.substitute(&[(read("m"), Expr::from(5))]))
        .unwrap();
    assert_eq!(at, read("1/16"), "{found:?}");
    check_stays(&[
        // Odd zeta values, as giac's Zeta(3).
        "Sum[Power[k, -3], k, 1, Infinity]",
        // No limit: giac's undef.
        "Sum[Power[-1, k], k, 0, Infinity]",
        // A pole at k = 0.
        "Sum[Power[k, -2], k, 0, Infinity]",
        // |x| < 1 unknown, and the sign of c.
        "Sum[Power[x, k], k, 0, Infinity]",
        "Sum[c, k, 1, Infinity]",
    ]);
}
