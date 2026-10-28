//! Random numbers against giac (CAS-DOC-33 §7, CAS-PLAN-27): giac's
//! generator is not this CAS's, so no draw compares. Each case is a
//! property both answers have, which holds the conventions to giac's: the
//! ranges, the sizes, permutations from 0, monic polynomials.

use cas_core::{Context, Expr, Integer, Rational, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// A number as a rational: itself, or a ball's midpoint.
fn value(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Ball(b) => b.mid().checked_to_rational(),
        _ => e.to_rational(),
    }
}

fn q(n: i64) -> Rational {
    Rational::from(Integer::from(n))
}

/// An integer in `[low, high]`.
fn integer_in(low: i64, high: i64) -> impl Fn(&Expr) -> bool {
    move |e| value(e).is_some_and(|v| v.is_integer() && q(low) <= v && v <= q(high))
}

/// A real in `[low, high]`.
fn real_in(low: i64, high: i64) -> impl Fn(&Expr) -> bool {
    move |e| value(e).is_some_and(|v| q(low) <= v && v <= q(high))
}

/// A real number.
fn real(e: &Expr) -> bool {
    value(e).is_some()
}

/// A list of `n` items that each pass `item`, distinct where `distinct`.
fn list_of(n: usize, distinct: bool, item: impl Fn(&Expr) -> bool) -> impl Fn(&Expr) -> bool {
    move |e| {
        let items = e.args();
        e.head() == Some(&Expr::from(Symbol::LIST))
            && items.len() == n
            && items.iter().all(&item)
            && (!distinct
                || items
                    .iter()
                    .enumerate()
                    .all(|(i, v)| !items[..i].contains(v)))
    }
}

/// The items of `expected`, in some order.
fn permutation_of(expected: &str) -> impl Fn(&Expr) -> bool {
    let mut want = ours(expected).expect("readable").args().to_vec();
    want.sort();
    move |e| {
        let mut found = e.args().to_vec();
        found.sort();
        e.head() == Some(&Expr::from(Symbol::LIST)) && found == want
    }
}

/// One of the items of `expected`.
fn one_of(expected: &str) -> impl Fn(&Expr) -> bool {
    let items = ours(expected).expect("readable").args().to_vec();
    move |e| items.contains(e)
}

/// A monic polynomial in `x` of degree `n`, its other coefficients
/// passing `coefficient`.
fn monic(n: i64, coefficient: impl Fn(&Expr) -> bool) -> impl Fn(&Expr) -> bool {
    move |p| {
        let ctx = Context::new();
        let x = Expr::symbol("x");
        let at = |head: Symbol, args: Vec<Expr>| ctx.eval(&Expr::call(head, args)).ok();
        at(Symbol::EXPONENT, vec![p.clone(), x.clone()]) == Some(Expr::from(n))
            && at(
                Symbol::COEFFICIENT,
                vec![p.clone(), x.clone(), Expr::from(n)],
            ) == Some(Expr::from(1))
            && (0..n).all(|k| {
                at(
                    Symbol::COEFFICIENT,
                    vec![p.clone(), x.clone(), Expr::from(k)],
                )
                .is_some_and(|c| coefficient(&c))
            })
    }
}

/// A property of an answer.
type Property<'p> = &'p dyn Fn(&Expr) -> bool;

/// giac's answers and this CAS's, each with its property.
fn check(cases: &[(&str, Property<'_>)]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, property) in cases {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if property(&theirs) && property(&ours) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `rand`, `randint` and `srand` (CAS-371).
#[test]
fn uniform() {
    check(&[
        ("rand()", &real_in(0, 1)),
        ("rand(10)", &integer_in(0, 9)),
        ("rand(-5)", &integer_in(-4, 0)),
        ("rand(2,5)", &real_in(2, 5)),
        ("rand(5,2)", &real_in(2, 5)),
        ("rand(3,1,49)", &list_of(3, true, integer_in(1, 49))),
        ("rand([a,b,c])", &one_of("[a,b,c]")),
        ("rand(2,[a,b,c,d])", &list_of(2, true, one_of("[a,b,c,d]"))),
        ("randint(6)", &integer_in(1, 6)),
        ("randint(-2,2)", &integer_in(-2, 2)),
        ("randint(5,2)", &integer_in(2, 5)),
        ("srand(5)", &integer_in(5, 5)),
    ]);
}

/// The variates (CAS-372).
#[test]
fn variates() {
    check(&[
        ("randnorm(0,1)", &real),
        ("randnorm(1,2)", &real),
        ("randexp(2)", &real_in(0, 1000)),
        ("randbinomial(10,1/2)", &integer_in(0, 10)),
        ("randpoisson(3)", &integer_in(0, 1000)),
        ("randgeometric(1/3)", &integer_in(1, 1000)),
        ("randchisquare(3)", &real_in(0, 1000)),
        ("randstudent(3)", &real),
        ("randfisher(2,3)", &real_in(0, 100_000)),
    ]);
}

/// Vectors, matrices, polynomials, permutations and samples (CAS-373).
#[test]
fn collections() {
    let row = |n: usize, low: i64, high: i64| list_of(n, false, integer_in(low, high));
    check(&[
        ("randvector(5)", &row(5, -99, 99)),
        ("randvector(4,6)", &row(4, 0, 5)),
        ("randvector(3,5..7)", &list_of(3, false, real_in(5, 7))),
        ("randvector(3,normald,0,1)", &list_of(3, false, real)),
        (
            "randvector(2,binomial,10,1/2)",
            &list_of(2, false, integer_in(0, 10)),
        ),
        ("randmatrix(2,3)", &list_of(2, false, row(3, -99, 99))),
        ("randmatrix(2,2,6)", &list_of(2, false, row(2, 0, 5))),
        ("randpoly(3,x)", &monic(3, integer_in(-99, 99))),
        ("randpoly(2,x,5)", &monic(2, integer_in(0, 4))),
        ("randperm(5)", &permutation_of("[0,1,2,3,4]")),
        ("randperm([a,b,c])", &permutation_of("[a,b,c]")),
        ("shuffle([a,b,c])", &permutation_of("[a,b,c]")),
        (
            "sample([a,b,c,d],2)",
            &list_of(2, true, one_of("[a,b,c,d]")),
        ),
    ]);
}
