//! Xcas's list commands against giac (CAS-555 to CAS-557): `find` and
//! `index` from 0, `count_eq`, `count_inf` and `count_sup` with `row` and
//! `col`, `deltalist`, `accumulate_head_tail`, `sizes`, `makevector`,
//! `makesuite`, `newList`, `enumerate` and `copy`, giac's answers read back
//! and compared with ours (`Giac::agrees`); and the departures, each with
//! giac's answer as it prints it: its errors, where ours leaves the call
//! with a warning, and the forms ours leaves.

use cas_core::{Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| format!("{e}")))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `find` and `index`, positions from 0 (CAS-555).
#[test]
fn positions() {
    check(&[
        "find(1,[3,x,1,2,1,3])",
        "find(2,[0,1,3,2,4,2,5])",
        "find(5,[1,2])",
        "find(1,[])",
        "find(1,[[1,2],1])",
        "find([1,2,1],1)",
        "find(1,[1.0,2,1])",
        "index([3,x,1,2,1,3],1)",
        "index([0,1,3,2,4,2,5],2)",
        "index(1,[3,1])",
        "index([3,x,1],x)",
    ]);
}

/// `count_eq`, `count_inf` and `count_sup`, flattened or by rows and
/// columns (CAS-555).
#[test]
fn counts() {
    check(&[
        "count_eq(1,[-2,1,0,1,2,-3])",
        "count_eq(4,[[3,4],[4,1]])",
        "count_eq(4,[[3,4],[4,1]],row)",
        "count_eq(4,[[3,4],[4,1]],col)",
        "count_eq(1,[[1,1],[1,2]],col)",
        "count_eq(4,[[3,4],[4,1]],matrix)",
        "count_eq(x,[x,y,x])",
        "count_eq(1,[1.0,1])",
        "count_eq(1,[[1,2],[1]])",
        "count_eq(1,1)",
        "count_eq(1,[])",
        "count_inf(1,[-2,-1,0,1,2,3])",
        "count_inf(4,[[3,5],[4,1]])",
        "count_inf(4,[[3,5],[4,1]],row)",
        "count_inf(4,[[3,5],[4,1]],col)",
        "count_inf(2,[[1,3],[0,5]],col)",
        "count_inf(x,[1,2])",
        "count_sup(1,[-2,-1,0,1,2,3])",
        "count_sup(3,[[3,5],[4,1]])",
        "count_sup(3,[[3,5],[4,1]],row)",
        "count_sup(3,[[3,5],[4,1]],col)",
        "count_sup(x,[1,2])",
        "count_sup(0,[sqrt(2),-pi,x])",
        "count_sup(0,[[1,2],[3]])",
    ]);
}

/// `deltalist`, `accumulate_head_tail` and `sizes` (CAS-556).
#[test]
fn differences_and_sums() {
    check(&[
        "deltalist([1,4,8,9])",
        "deltalist([1,8,4,9])",
        "deltalist([[1,2],[4,8]])",
        "deltalist([x,y,z])",
        "accumulate_head_tail([0,1,2,3,4,5,6,7,8,9],3,2)",
        "accumulate_head_tail([1,2,3],1,1)",
        "accumulate_head_tail([a,b,c,d],2,1)",
        "sizes([[1,2,3],[1,2],[1]])",
        "sizes([])",
        "sizes([[],[a]])",
    ]);
}

/// `makevector`, `makesuite` (its sequence bracketed to be read),
/// `newList`, `enumerate` and `copy` (CAS-557).
#[test]
fn builders() {
    check(&[
        "makevector(1,2,3)",
        "makevector(1)",
        "makevector()",
        "makevector([1,2])",
        "makevector(x,[1])",
        "[makesuite([1,2,3])]",
        "[makesuite(1,2)]",
        "makesuite(1)",
        "newList(4)",
        "newList(0)",
        "newList(-1)",
        "newList(-3)",
        "enumerate([a,b,c])",
        "enumerate([a,b,c],1)",
        "enumerate([a,b],-1)",
        "enumerate([a,b],x)",
        "enumerate([])",
        "enumerate()",
        "copy([1,2,3])",
        "copy(x+1)",
        "copy([[1,2],x])",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
const DEPARTURES: [(&str, &str, &str); 22] = [
    // giac's errors; ours leaves the call, with the error as a warning.
    (
        "index([1,2],5)",
        "GIAC_ERROR: index([1,2],5) \n Error: Bad Argument Value",
        "index([1,2],5)",
    ),
    (
        "index([1,2,1],1,1)",
        "GIAC_ERROR: index([1,2,1],1,1) \n Error: Bad Argument Value",
        "index([1,2,1],1,1)",
    ),
    (
        "find(x,y)",
        "GIAC_ERROR: find(x,y) \n Error: Bad Argument Value",
        "find(x,y)",
    ),
    (
        "count_eq(1)",
        "GIAC_ERROR: count_eq(1) \n Error: Bad Argument Value",
        "count_eq(1)",
    ),
    (
        "count_eq(1,[1,2],row)",
        "GIAC_ERROR: Error: Bad Argument Type",
        "count_eq(1,[1,2],row)",
    ),
    (
        "deltalist([1])",
        "GIAC_ERROR: deltalist([1]) \n Error: Invalid dimension",
        "deltalist([1])",
    ),
    (
        "deltalist([])",
        "GIAC_ERROR: deltalist([]) \n Error: Invalid dimension",
        "deltalist([])",
    ),
    (
        "deltalist(5)",
        "GIAC_ERROR: deltalist(5) \n Error: Bad Argument Value",
        "deltalist(5)",
    ),
    (
        "accumulate_head_tail([1,2,3],2,2)",
        "GIAC_ERROR: Error: Bad Argument Value",
        "accumulate_head_tail([1,2,3],2,2)",
    ),
    (
        "accumulate_head_tail([1,2,3,4],0,0)",
        "GIAC_ERROR: Error: Bad Argument Value",
        "accumulate_head_tail([1,2,3,4],0,0)",
    ),
    (
        "accumulate_head_tail([1,2,3],1,1.0)",
        "GIAC_ERROR: Error: Bad Argument Value",
        "accumulate_head_tail([1,2,3],1,1.0)",
    ),
    (
        "accumulate_head_tail(5,1,1)",
        "GIAC_ERROR: Error: Bad Argument Value",
        "accumulate_head_tail(5,1,1)",
    ),
    (
        "sizes([1,[1,2]])",
        "GIAC_ERROR: sizes([1,[1,2]]) \n Error: Bad Argument Value",
        "sizes([1,[1,2]])",
    ),
    (
        "sizes(5)",
        "GIAC_ERROR: sizes(5) \n Error: Bad Argument Type",
        "sizes(5)",
    ),
    (
        "newList(x)",
        "GIAC_ERROR: newList(x) \n Error: Bad Argument Value",
        "newList(x)",
    ),
    (
        "newList(3/2)",
        "GIAC_ERROR: newList(3/2) \n Error: calculation size limit exceeded",
        "newList(3/2)",
    ),
    (
        "newList(2,3)",
        "GIAC_ERROR: newList(2,3) \n Error: calculation size limit exceeded",
        "newList(2,3)",
    ),
    (
        "enumerate(x)",
        "GIAC_ERROR: enumerate(x) \n Error: Bad Argument Type",
        "enumerate(x)",
    ),
    // The string forms of find and index wait for phase 1.6.
    (
        "find(\"ab\",\"abracadabrant\")",
        "[0,7]",
        "find(\"ab\",\"abracadabrant\")",
    ),
    (
        "index(\"abracadabrant\",\"c\")",
        "4",
        "index(\"abracadabrant\",\"c\")",
    ),
    // find's start position, and enumerate of more than two arguments,
    // stay.
    ("find(1,[1,2,1],1)", "[2]", "find(1,[1,2,1],1)"),
    (
        "enumerate([a,b],1,2)",
        "[[0,[a,b]],[1,1],[2,2]]",
        "enumerate([a,b],1,2)",
    ),
];

/// The departures of [`DEPARTURES`], then `newList` beyond this CAS's
/// 65536 elements.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, answer, decided) in DEPARTURES {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac's LIST_SIZE_LIMIT is larger than lists.rs's MAX_ELEMENTS: ours
    // leaves the call.
    let (input, size) = ("size(newList(65537))", "65537");
    let printed = giac.eval(input).unwrap();
    if printed != size {
        wrong.push(format!("{input}: giac prints {printed}, recorded {size}"));
    }
    let stays = ours("newList(65537)")
        .is_ok_and(|e| e.head().and_then(Expr::as_symbol) == Some(Symbol::NEW_LIST));
    if !stays {
        wrong.push("newList(65537): ours does not stay".to_owned());
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
