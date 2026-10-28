//! List operations against giac (CAS-DOC-28 §5, CAS-PLAN-22): the answers
//! compared as expressions after evaluation. Not `seq(f,k=a..b)`, which
//! giac answers with a sequence where this CAS answers a list, nor
//! `count(2,l)`, which giac takes for the constant function 2.

use cas_core::Context;
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

#[test]
fn lists() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "seq(k^2,k,1,5)",
        "seq(k,k,1,10,3)",
        "seq(j,j,5,1,-1)",
        "makelist(x->x^2,1,9,2)",
        "sort([3,1,2])",
        "sort([3,1/2,sqrt(2)])",
        "sort([3,1,2],(a,b)->a>b)",
        "sort([b,a,c])",
        "revlist([1,2,3])",
        "rotate([1,2,3,4],2)",
        "rotate([1,2,3])",
        "concat([1,2],[3,4])",
        "concat([1,2],3)",
        "append([1,2],3)",
        "prepend([1,2],0)",
        "head([1,2,3])",
        "tail([1,2,3])",
        "flatten([[1,2],[3,[4,5]]])",
        "dim([[1,2,3],[4,5,6]])",
        "dim([1,2,3])",
        "select(x->x>2,[1,2,3,4])",
        "remove(x->x>2,[1,2,3,4])",
        "remove(2,[1,2,2,3])",
        "count(x->x>2,[1,2,3,4])",
        "member(3,[1,2,3])",
        "member(5,[1,2,3])",
        "zip((x,y)->x*y,[1,2],[3,4])",
        "cumsum([1,2,3])",
        "cumsum([a,b,c])",
        "sum([1,2,3,a])",
        "product([1,2,3,4])",
    ] {
        let theirs = giac.value(input).map(|e| ctx.eval(&e));
        match (theirs, ours(input)) {
            (Ok(Ok(theirs)), Ok(ours)) if theirs == ours => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
