//! Tests of `simplex_reduce` (CAS-606): giac's answers to its probe rows,
//! both forms, an unbounded program, decimals, the negative right side
//! that leaves the call, and shapes that leave it.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// giac's text in full form: its lists, its `+infinity`, its decimals as
/// hex floats (C0047).
fn full(giac: &str) -> String {
    giac.replace('[', "List[")
        .replace("+infinity", "Infinity")
        .replace("1.5", "0x1.8p+0")
        .replace("0.75", "0x1.8p-1")
        .replace("0.25", "0x1p-2")
        .replace("0.5", "0x1p-1")
        .replace(',', ", ")
}

/// `simplex_reduce(args)` bracketed, ours, and giac's answer bracketed.
fn call(args: &str) -> Expr {
    read(&format!("List[SimplexReduce[{}]]", full(args)))
}

/// giac's probe rows (probe3.out 1-8) and its answers.
#[test]
fn probes() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (args, answer) in [
        (
            "[[3,2,2],[1,1,1]],[3,4],[1,2,3]",
            "9/2,[0,0,3/2,0,5/2],[[3/2,1,1,1/2,0,3/2],[-1/2,0,0,-1/2,1,5/2],[7/2,1,0,3/2,0,9/2]]",
        ),
        (
            "[[3,2,2,1,0,3],[1,1,1,0,1,4],[-1,-2,-3,0,0,0]]",
            "9/2,[0,0,3/2,0,5/2],[[3/2,1,1,1/2,0,3/2],[-1/2,0,0,-1/2,1,5/2],[7/2,1,0,3/2,0,9/2]]",
        ),
        (
            "[[-3,2],[1,1]],[3,4],[1,2]",
            "7,[1,3,0,0],[[0,1,1/5,3/5,3],[1,0,-1/5,2/5,1],[0,0,1/5,8/5,7]]",
        ),
        (
            "[[-3,2,1,0,3],[1,1,0,1,4],[-1,-2,0,0,0]]",
            "7,[1,3,0,0],[[0,1,1/5,3/5,3],[1,0,-1/5,2/5,1],[0,0,1/5,8/5,7]]",
        ),
        (
            "[[2,1,1,1,0,0,2],[1,2,3,0,1,0,5],[2,2,1,0,0,1,6],[-3,-1,-3,1,-1,2,0]]",
            "-1,[1,0,0,0,4,4],[[1,1/2,1/2,1/2,0,0,1],[0,3/2,5/2,-1/2,1,0,4],[0,1,0,-1,0,1,4],\
             [0,0,1,4,0,0,-1]]",
        ),
        (
            "[[1,1],[1,-1]],[2,1],[1,1]",
            "2,[3/2,1/2,0,0],[[0,1,1/2,-1/2,1/2],[1,0,1/2,1/2,3/2],[0,0,1,0,2]]",
        ),
        ("[[1,-1]],[1],[1,1]", "+infinity,[],[[1,-1,1,1],[0,-2,1,1]]"),
        (
            "[[1.5,2],[1,1]],[3,4],[1,2]",
            "3,[0,3/2,0,5/2],[[0.75,1,1/2,0,3/2],[0.25,0,-1/2,1,5/2],[0.5,0,1,0,3]]",
        ),
    ] {
        let found = ctx.eval(&call(args)).unwrap();
        let expected = ctx.eval(&read(&format!("List[{}]", full(answer)))).unwrap();
        if found != expected {
            wrong.push(format!("{args}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Inputs that leave the call: a negative right side (probe3.out 9, where
/// giac pivots from an infeasible basis), giac's dimension error, ragged
/// and symbolic tableaux, other arities.
#[test]
fn left() {
    let ctx = Context::new();
    for args in [
        "[[3,2,2],[1,1,1]],[-3,4],[1,2,3]",
        "[[1,2,3],[4,5,6],[7,8,9]]",
        "[[1,2,3],[4,5]]",
        "[[1,a,3],[4,5,6]]",
        "[[1,2],[3,4]],[1,2]",
        "[]",
    ] {
        let found = ctx.eval(&call(args)).unwrap();
        let inner = &found.args()[0];
        assert_eq!(
            inner.head().and_then(Expr::as_symbol),
            Some(Symbol::SIMPLEX_REDUCE),
            "{args}: {found:?}"
        );
    }
}
