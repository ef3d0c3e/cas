//! Tests of products of lists: tables checked against giac 1.9, and
//! properties on random matrices.

use proptest::prelude::*;

use crate::{Context, EvalError, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

const A: &str = "List[List[1, 2], List[3, 4]]";
const B: &str = "List[List[5, 6], List[7, 8]]";

fn times(factors: &[&str]) -> String {
    format!("Times[{}]", factors.join(", "))
}

#[test]
fn products() {
    check(&[
        // As in Xcas: the dot product, matrix times vector, a vector as a
        // row times a matrix, and products of matrices.
        ("Times[List[1, 2], List[3, 4]]", "11"),
        (&times(&[A, "List[5, 6]"]), "List[17, 39]"),
        (&times(&["List[5, 6]", A]), "List[23, 34]"),
        (&times(&[A, B]), "List[List[19, 22], List[43, 50]]"),
        // `A B` is not `B A`.
        (&times(&[B, A]), "List[List[23, 34], List[31, 46]]"),
        (
            "Times[List[List[1, 2, 3], List[4, 5, 6]], List[List[1, 2], List[3, 4], List[5, 6]]]",
            "List[List[22, 28], List[49, 64]]",
        ),
        (
            "Times[List[List[1, 2]], List[List[3], List[4]]]",
            "List[List[11]]",
        ),
        (
            "Times[List[List[3], List[4]], List[List[1, 2]]]",
            "List[List[3, 6], List[4, 8]]",
        ),
        // From the left: `(u·v) w`.
        ("Times[List[1, 2], List[3, 4], List[5, 6]]", "List[55, 66]"),
        // Scalars scale every entry, symbols too.
        ("Times[2, List[1, 2]]", "List[2, 4]"),
        ("Times[x, List[1, 2]]", "List[x, Times[2, x]]"),
        ("Times[List[1, 2], x]", "List[x, Times[2, x]]"),
        (
            "Times[2, x, A]".replace('A', A).as_str(),
            "List[List[Times[2, x], Times[4, x]], List[Times[6, x], Times[8, x]]]",
        ),
        // Symbolic entries.
        (
            "Times[List[a, b], List[c, d]]",
            "Plus[Times[a, c], Times[b, d]]",
        ),
        ("Times[List[], List[]]", "0"),
        // Lists that are neither vectors nor matrices stay.
        (
            "Times[List[1, List[2]], List[3]]",
            "Times[List[1, List[2]], List[3]]",
        ),
    ]);
}

#[test]
fn powers_and_sums() {
    check(&[
        (&format!("Power[{A}, 2]"), "List[List[7, 10], List[15, 22]]"),
        (
            &format!("Power[{A}, 3]"),
            "List[List[37, 54], List[81, 118]]",
        ),
        (&format!("Power[{A}, 1]"), A),
        (&format!("Power[{A}, 0]"), "List[List[1, 0], List[0, 1]]"),
        // The inverse belongs to Linear algebra; a vector has no power.
        (&format!("Power[{A}, -1]"), &format!("Power[{A}, -1]")),
        ("Power[List[1, 2], 2]", "Power[List[1, 2], 2]"),
        ("Power[List[1, 2], 0]", "Power[List[1, 2], 0]"),
        (
            "Power[List[List[a, 0], List[0, b]], 2]",
            "List[List[Power[a, 2], 0], List[0, Power[b, 2]]]",
        ),
        // Fibonacci by squaring.
        (
            "Power[List[List[1, 1], List[1, 0]], 30]",
            "List[List[1346269, 832040], List[832040, 514229]]",
        ),
        // Sums are entrywise.
        ("Plus[List[1, 2], List[3, 4]]", "List[4, 6]"),
        (&format!("Plus[{A}, {B}]"), "List[List[6, 8], List[10, 12]]"),
        // `.*` multiplies entries.
        ("ElementwiseTimes[List[1, 2], List[3, 4]]", "List[3, 8]"),
        (
            &format!("ElementwiseTimes[{A}, {B}]"),
            "List[List[5, 12], List[21, 32]]",
        ),
        ("ElementwiseTimes[2, 3]", "6"),
        ("ElementwiseTimes[x, y]", "Times[x, y]"),
    ]);
}

/// `./` and `.^`, as in Xcas on lists of one shape.
#[test]
fn elementwise_quotients_and_powers() {
    check(&[
        (
            "ElementwiseDivide[List[1, 2, 3], List[2, 4, 6]]",
            "List[1/2, 1/2, 1/2]",
        ),
        (
            "ElementwisePower[List[1, 2, 3], List[2, 2, 2]]",
            "List[1, 4, 9]",
        ),
        (
            &format!("ElementwiseDivide[{A}, {B}]"),
            "List[List[1/5, 1/3], List[3/7, 1/2]]",
        ),
        (
            &format!("ElementwisePower[{A}, {A}]"),
            "List[List[1, 4], List[27, 256]]",
        ),
        ("ElementwiseDivide[x, y]", "Times[x, Power[y, -1]]"),
        ("ElementwisePower[x, 2]", "Power[x, 2]"),
        ("ElementwiseDivide[6, 3]", "2"),
        ("ElementwiseDivide[1, 0]", "ComplexInfinity"),
        ("ElementwiseDivide[1, 2, 3]", "ElementwiseDivide[1, 2, 3]"),
    ]);
}

#[test]
fn dimensions_that_do_not_fit() {
    let ctx = Context::new();
    for input in [
        "Times[List[1, 2], List[3, 4, 5]]".to_owned(),
        format!("Times[{A}, List[1, 2, 3]]"),
        "Times[List[List[1, 2]], List[List[1, 2]]]".to_owned(),
        "Times[List[1, 2, 3], List[List[1, 2], List[3, 4]]]".to_owned(),
        "Power[List[List[1, 2, 3], List[4, 5, 6]], 2]".to_owned(),
        // No broadcasting (O8).
        "Plus[List[1, 2], 1]".to_owned(),
        "ElementwiseTimes[List[1, 2], 3]".to_owned(),
        "ElementwiseTimes[List[1, 2], List[1, 2, 3]]".to_owned(),
        // `./` and `.^` too, where Xcas spreads a scalar.
        "ElementwiseDivide[List[1, 2], 2]".to_owned(),
        "ElementwisePower[List[1, 2, 3], 2]".to_owned(),
        "ElementwiseDivide[List[1, 2], List[1, 2, 3]]".to_owned(),
    ] {
        assert_eq!(ctx.eval(&read(&input)), Err(EvalError::Shape), "{input}");
    }
}

/// The entries are built unevaluated, for the loop, and the scalar reaches
/// each entry of a matrix: `Times[c, row]` would give the same result after
/// one more pass through `Times` per row.
#[test]
fn entries_are_built_unevaluated() {
    assert_eq!(
        super::scale(&Expr::from(2), &read(A)),
        read("List[List[Times[2, 1], Times[2, 2]], List[Times[2, 3], Times[2, 4]]]")
    );
    assert_eq!(
        super::dot(&read("List[a, b]"), &read(A)),
        Ok(read(
            "List[Plus[Times[a, 1], Times[b, 3]], Plus[Times[a, 2], Times[b, 4]]]"
        ))
    );
}

/// A matrix with small integer entries, `rows × columns`.
fn matrix(rows: usize, columns: usize) -> impl Strategy<Value = Expr> {
    prop::collection::vec(prop::collection::vec(-4i64..=4, columns), rows).prop_map(|m| {
        let rows: Vec<Expr> = m
            .into_iter()
            .map(|r| {
                let r: Vec<Expr> = r.into_iter().map(Expr::from).collect();
                Expr::apply(&Symbol::LIST.into(), &r)
            })
            .collect();
        Expr::apply(&Symbol::LIST.into(), &rows)
    })
}

fn product(a: &Expr, b: &Expr) -> Expr {
    Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    /// `(A B) C = A (B C)` for matrices of fitting sizes.
    #[test]
    fn products_of_matrices_are_associative(
        (a, b, c) in (1usize..4, 1usize..4, 1usize..4, 1usize..4)
            .prop_flat_map(|(m, n, p, q)| (matrix(m, n), matrix(n, p), matrix(p, q)))
    ) {
        let ctx = Context::new();
        let left = ctx.eval(&product(&product(&a, &b), &c)).unwrap();
        let right = ctx.eval(&product(&a, &product(&b, &c))).unwrap();
        prop_assert_eq!(left, right);
    }

    /// A power by squaring is the product of as many factors, from the left.
    #[test]
    fn powers_are_repeated_products(m in (1usize..4).prop_flat_map(|n| matrix(n, n)), k in 1i64..7) {
        let ctx = Context::new();
        let power = ctx.eval(&Expr::apply(&Symbol::POWER.into(), &[m.clone(), k.into()])).unwrap();
        let mut repeated = m.clone();
        for _ in 1..k {
            repeated = ctx.eval(&product(&repeated, &m)).unwrap();
        }
        prop_assert_eq!(power, repeated);
    }
}
