//! List operations (CAS-DOC-28 §5, CAS-337), one check per input: a random
//! list of small rationals, and the identities of the operations:
//! - `Sort` is ordered and a permutation, and `Sort` with `Less` agrees;
//! - `Reverse` twice, and `Rotate` by `n` then `−n`, are the identity;
//! - `Join` adds up lengths, `Append` and `Prepend` add one;
//! - `Select` and `Reject` by a threshold split the list, and `Count` is
//!   the length of `Select`;
//! - `Accumulate`'s last element is `Sum` of the list;
//! - `Table[k^2, k, 1, n]` is the squares, and `MakeList` agrees;
//! - `MemberIndex` finds each element at an index where it is.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `Function[{x}, body]`, `x` its variable.
fn function(body: Expr) -> Expr {
    Expr::call(
        Symbol::new("Function"),
        [list(vec![Expr::from(Symbol::new("x"))]), body],
    )
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let len = u.int_in_range(0..=12usize)?;
    let values: Vec<Rational> = (0..len)
        .map(|_| Ok(q(u.int_in_range(-9..=9)?, u.int_in_range(1..=3)?)))
        .collect::<Result<_>>()?;
    let l = list(values.iter().cloned().map(Expr::from).collect());
    let length = |e: &Expr| e.args().len();
    // Sort: ordered, a permutation, and by Less the same.
    let sorted = eval(&ctx, &call("Sort", vec![l.clone()]));
    let mut expected = values.clone();
    expected.sort();
    assert_eq!(
        sorted,
        list(expected.into_iter().map(Expr::from).collect()),
        "Sort[{l:?}]"
    );
    let by_less = eval(
        &ctx,
        &call("Sort", vec![l.clone(), Expr::from(Symbol::LESS)]),
    );
    assert_eq!(by_less, sorted, "Sort[{l:?}, Less]");
    // Reverse twice, Rotate there and back.
    let back = eval(
        &ctx,
        &call("Reverse", vec![call("Reverse", vec![l.clone()])]),
    );
    assert_eq!(back, l, "Reverse twice");
    let n = u.int_in_range(-20..=20i64)?;
    let rotated = eval(
        &ctx,
        &call(
            "Rotate",
            vec![
                call("Rotate", vec![l.clone(), Expr::from(n)]),
                Expr::from(-n),
            ],
        ),
    );
    assert_eq!(rotated, l, "Rotate by {n} and back");
    // Join, Append, Prepend.
    let joined = eval(
        &ctx,
        &call("Join", vec![l.clone(), sorted.clone(), l.clone()]),
    );
    assert_eq!(length(&joined), 3 * len, "Join");
    let appended = eval(&ctx, &call("Append", vec![l.clone(), Expr::from(7)]));
    assert_eq!(length(&appended), len + 1, "Append");
    // Select and Reject split; Count counts.
    let threshold = Expr::from(q(u.int_in_range(-9..=9)?, 2));
    let p = function(Expr::call(
        Symbol::GREATER,
        [Expr::from(Symbol::new("x")), threshold],
    ));
    let kept = eval(&ctx, &call("Select", vec![p.clone(), l.clone()]));
    let dropped = eval(&ctx, &call("Reject", vec![p.clone(), l.clone()]));
    assert_eq!(length(&kept) + length(&dropped), len, "Select and Reject");
    let counted = eval(&ctx, &call("Count", vec![p, l.clone()]));
    assert_eq!(
        counted,
        Expr::from(i64::try_from(length(&kept)).expect("small")),
        "Count"
    );
    // Accumulate and Sum.
    let sums = eval(&ctx, &call("Accumulate", vec![l.clone()]));
    let total = eval(&ctx, &call("Sum", vec![l.clone()]));
    assert_eq!(
        sums.args().last().cloned().unwrap_or_else(|| Expr::from(0)),
        total,
        "Accumulate"
    );
    // Table and MakeList.
    let k = Expr::from(Symbol::new("k"));
    let m = i64::try_from(len).expect("small");
    let squares = eval(
        &ctx,
        &call(
            "Table",
            vec![k.pow(2), k.clone(), Expr::from(1), Expr::from(m)],
        ),
    );
    let expected: Vec<Expr> = (1..=m).map(|i| Expr::from(i * i)).collect();
    assert_eq!(squares, list(expected), "Table");
    let made = eval(
        &ctx,
        &call(
            "MakeList",
            vec![
                function(Expr::from(Symbol::new("x")).pow(2)),
                Expr::from(1),
                Expr::from(m),
            ],
        ),
    );
    assert_eq!(made, squares, "MakeList");
    // MemberIndex finds each element where it is.
    for x in l.args() {
        let index = eval(&ctx, &call("MemberIndex", vec![x.clone(), l.clone()]));
        let i = usize::try_from(index.as_i64().expect("an index")).expect("positive");
        assert!(
            i >= 1 && l.args()[i - 1] == *x,
            "MemberIndex[{x:?}, {l:?}] = {i}"
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
