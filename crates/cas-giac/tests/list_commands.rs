//! Xcas's list commands against giac (CAS-555 to CAS-561): `find` and
//! `index` from 0, `count_eq`, `count_inf` and `count_sup` with `row` and
//! `col`, `deltalist`, `accumulate_head_tail`, `sizes`, `makevector`,
//! `makesuite`, `newList`, `enumerate`, `copy`, `apply`, `foldl`, `foldr`,
//! `multiply`, `insert`, `subsop`, `heappush`, `nextperm`, `prevperm` and
//! `groupermu`, giac's answers read back
//! and compared with ours (`Giac::agrees`); and the departures, each with
//! giac's answer as it prints it: its errors, where ours leaves the call
//! with a warning, and the forms ours leaves; and random small integer and
//! rational lists a run to each family but the departures.

use cas_core::{Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

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

/// `apply`, with and without `matrix`, `foldl` and `foldr` (CAS-558).
#[test]
fn maps_and_folds() {
    check(&[
        "apply(x->x^3,[1,2,3])",
        "apply(x->x+1,[[1,2,3],[1,2,3]],matrix)",
        "apply(x->x+1,[[1,2],[3,4]],matrix)",
        "apply(f,[[1,2],[3,4]],matrix)",
        "apply(f,[[1,2]],matrix)",
        "apply(x->x+1,[],matrix)",
        "apply(f,[1,2])",
        "apply(x->x^2,[1,2])",
        "apply(x->[x,x],[1,2])",
        "apply(sin,[0,pi])",
        "foldl(F,init,a,b,c)",
        "foldr(F,init,a,b,c)",
        "foldl(F,init,a)",
        "foldr(F,init,a)",
        "foldr(F,init,a,b,c,d)",
        "foldl(F,init,[1,2])",
        "foldl((x,y)->x+y,0,1,2,3)",
        "foldl((x,y)->x*y,1,2,3,4)",
        "foldr((x,y)->x-y,0,1,2,3)",
    ]);
}

/// `multiply`, giac's `*` (CAS-558): dot products over the shorter vector,
/// matrix products, scalars.
#[test]
fn products() {
    check(&[
        "multiply(41,-4)",
        "multiply([4,1],[-4,2])",
        "multiply([1,2],[3,4,5])",
        "multiply([1,2,3],[3,4])",
        "multiply([1,i],[1,i])",
        "multiply([x,y],[a,b])",
        "multiply([[4,1],[-4,1]],[[4,1],[-4,1]])",
        "multiply([[1,2]],[[3],[4]])",
        "multiply([1,2],[[1,0],[0,1]])",
        "multiply([[4,1],[-4,1]],[1,2])",
        "multiply([1,2],[3,4],[5,6])",
        "multiply([[1,2],[3,4]],[[1,2],[3,4]],[1,1])",
        "multiply(2,[[1,2],[3,4]])",
        "multiply([1,2],3)",
        "multiply([1,2],x)",
        "multiply(x,y)",
        "multiply(2,3,4)",
        "multiply([1,2])",
        "multiply(3)",
        "multiply()",
    ]);
}

/// `insert` and `subsop` from 0 (CAS-559); giac's quoted `'1=NULL'` is
/// written `1=NULL`, which this language reads.
#[test]
fn insertions_and_substitutions() {
    check(&[
        "insert([0,1,2,3],2,5)",
        "insert([0,1,2,3],0,5)",
        "insert([0,1,2,3],4,5)",
        "insert([],0,5)",
        "insert([0,1,2,3],2,[5,6])",
        "insert([[0,1],[2,3]],1,[9,9])",
        "subsop([0,1,2,3,4],2=sqrt(2)+1)",
        "subsop([[1,2],[3,4]],[1,1]=5)",
        "subsop([[1,2],[3,4]],[1,-1]=9)",
        "subsop([[1,2],[3,4]],[-1,0]=9)",
        "subsop([[1,2],[3,4]],[0]=5)",
        "subsop([[1,2],[3,4]],1=[10,8])",
        "subsop([[1,2],[3,4]],[1,1]=5,0=7)",
        "subsop([[1,2],[3,4]],1=7,[1,1]=5)",
        "subsop([0,1,2,3],1=NULL)",
        "subsop([0,1,2,3],0=NULL,1=NULL)",
        "subsop([0,1],0=NULL,1=NULL)",
        "subsop([0,1,2,3],1=9,1=NULL)",
        "subsop([0,1,2,3],1=NULL,1=9)",
        "subsop([0,1,2,3],1=9,2=8)",
        "subsop([0,1,2,3],1=9,1=8)",
        "subsop([0,1,2,3],-1=9)",
        "subsop([0,1,2,3],-4=9)",
        "subsop([0,1,2,3],[1]=9)",
        "subsop([0,1,2,3],1,9)",
        "subsop([0,1,2,3],1=[])",
        "subsop([0,1,2,3],x=9)",
        "subsop([0,1,2,3],a)",
        "subsop([0,1,2,3],[1,0]=9)",
        "subsop([[1,2],[3,4]],[1,1,0]=9)",
    ]);
}

/// `heappush`, libstdc++'s `push_heap`, by `<` or a comparator (CAS-560);
/// `heapify` and `heappop` are departures.
#[test]
fn heaps() {
    check(&[
        "heappush([9,6,5,1,1,4,2,3],7)",
        "heappush([9,6,5,1,1,4,2,3],10)",
        "heappush([1],7)",
        "heappush([],7)",
        "heappush([3,1],2,(a,b)->a>b)",
        "heappush([1,3,2],0,(a,b)->a>b)",
        "heappush([9,6,5,1,1,4,2,3],7,(a,b)->a<b)",
        "heappush([1/2,1/3],2/3)",
        "heappush([2.5,1],3)",
    ]);
}

/// `nextperm`, `prevperm` and `groupermu`, giac's order (CAS-561).
#[test]
fn permutations() {
    check(&[
        "nextperm([0,2,1,3])",
        "nextperm([0,3,2,1])",
        "nextperm([0,1,2])",
        "nextperm([2,1,0])",
        "nextperm([])",
        "nextperm([0,1,2,3])",
        "prevperm([0,1,3,2])",
        "prevperm([0,1,2,3])",
        "prevperm([2,1,0])",
        "prevperm([1,0,2])",
        "prevperm([3,2,1,0])",
        "groupermu([1,0,2],[0,2,1])",
        "groupermu([[1,0,2],[0,2,1]])",
        "groupermu([1,2,0],[3,1,2,0])",
        "groupermu([0,1],[1,0])",
        "groupermu([1,0],[0,1,2])",
        "groupermu([0,1,2],[0,1,2])",
        "groupermu([],[])",
        "groupermu([0],[])",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
const DEPARTURES: [(&str, &str, &str); 65] = [
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
    // newList of a non-integer: giac reads `g.val` before its type test
    // (ti89.cc `_newList`), stray bits whose error differs between runs
    // (CI 45113559662: newList(x) gave both errors): ours-only, in cas-core.
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
    // apply: giac's x^2 of a vector is its dot product with itself; apply of
    // several lists, or of a non-list, stays; matrix on a vector is giac's
    // extra argument, ours leaves the call (giac's error elsewhere).
    ("apply(x->x^2,[[1,2],[3,4]])", "[5,25]", "[[1,2]^2,[3,4]^2]"),
    (
        "apply(f,[1,2],matrix)",
        "[f(1,'matrix'),f(2,'matrix')]",
        "apply(f,[1,2],matrix)",
    ),
    (
        "apply(x->x+1,[[1,2],[3]],matrix)",
        "GIAC_ERROR: Error: Bad Argument Type",
        "apply(x->x+1,[[1,2],[3]],matrix)",
    ),
    ("apply(f,x)", "[]", "apply(f,x)"),
    (
        "apply(f,[1,2],[3,4])",
        "[f([1,3]),f([2,4])]",
        "apply(f,[1,2],[3,4])",
    ),
    // foldl and foldr of fewer than three arguments: giac's errors.
    (
        "foldl(F)",
        "GIAC_ERROR: foldl(F) \n Error: Bad Argument Type",
        "foldl(F)",
    ),
    (
        "foldl(F,init)",
        "GIAC_ERROR: Wrong number of arguments. Error: Bad Argument Value",
        "foldl(F,init)",
    ),
    (
        "foldr(F,init)",
        "GIAC_ERROR: Wrong number of arguments. Error: Bad Argument Value",
        "foldr(F,init)",
    ),
    // multiply: a vector longer than the matrix's rows stays, as do
    // dimensions that do not fit (giac's error).
    (
        "multiply([1,2,3],[[1,0],[0,1]])",
        "[1,2]",
        "multiply([1,2,3],[[1,0],[0,1]])",
    ),
    (
        "multiply([[1,2],[3,4]],[[1,2,3]])",
        "GIAC_ERROR:  Error: Invalid dimension",
        "multiply([[1,2],[3,4]],[[1,2,3]])",
    ),
    // insert: giac's errors; a decimal position, and strings (phase 1.6),
    // stay.
    (
        "insert([0,1,2,3],5,5)",
        "GIAC_ERROR: insert([0,1,2,3],5,5) \n Error: Invalid dimension",
        "insert([0,1,2,3],5,5)",
    ),
    (
        "insert([0,1,2,3],-1,5)",
        "GIAC_ERROR: insert([0,1,2,3],-1,5) \n Error: Invalid dimension",
        "insert([0,1,2,3],-1,5)",
    ),
    (
        "insert(x,0,1)",
        "GIAC_ERROR: insert(x,0,1) \n Error: Bad Argument Value",
        "insert(x,0,1)",
    ),
    (
        "insert([0,1],x,5)",
        "GIAC_ERROR: insert([0,1],x,5) \n Error: Bad Argument Value",
        "insert([0,1],x,5)",
    ),
    (
        "insert([0,1,2,3],2,5,6)",
        "GIAC_ERROR: insert([0,1,2,3],2,5,6) \n Error: Bad Argument Value",
        "insert([0,1,2,3],2,5,6)",
    ),
    (
        "insert([0,1,2,3],2.0,5)",
        "[0,1,5,2,3]",
        "insert([0,1,2,3],2.0,5)",
    ),
    (
        "insert(\"0123\",2,\"5\")",
        "\"01523\"",
        "insert(\"0123\",2,\"5\")",
    ),
    // subsop: giac's errors; one argument, a non-list and strings stay;
    // a nested NULL deletes, where giac's unquoted NULL leaves the list.
    (
        "subsop([0,1,2,3],4=9)",
        "GIAC_ERROR: subsop([0,1,2,3],4=9) \n Error: Invalid dimension",
        "subsop([0,1,2,3],4=9)",
    ),
    (
        "subsop([0,1,2,3],-5=9)",
        "GIAC_ERROR: subsop([0,1,2,3],-5=9) \n Error: Invalid dimension",
        "subsop([0,1,2,3],-5=9)",
    ),
    (
        "subsop([[1,2],[3,4]],[1,5]=9)",
        "GIAC_ERROR: subsop([[1,2],[3,4]],[1,5]=9) \n Error: Invalid dimension",
        "subsop([[1,2],[3,4]],[1,5]=9)",
    ),
    (
        "subsop([0,1,2,3],[]=9)",
        "GIAC_ERROR: subsop([0,1,2,3],[]=9) \n Error: Invalid dimension",
        "subsop([0,1,2,3],[]=9)",
    ),
    ("subsop([0,1,2,3])", "0", "subsop([0,1,2,3])"),
    ("subsop(5,1=2)", "5", "subsop(5,1=2)"),
    ("subsop(f(a,b),1=c)", "f(c)", "subsop(f(a,b),1=c)"),
    (
        "subsop(\"abc\",1=\"x\")",
        "\"abc\"",
        "subsop(\"abc\",1=\"x\")",
    ),
    (
        "subsop([[1,2],[3,4]],[1,1]=NULL)",
        "[[1,2],[3,4]]",
        "[[1,2],[3]]",
    ),
    ("subsop([[1,2],[3,4]],[1]=NULL)", "[[1,2],[3,4]]", "[[1,2]]"),
    // Heaps (CAS-560): giac works in place and prints 1 for heapify; ours
    // returns the heap. giac orders symbols where `Less` does not decide:
    // ours stays. heappop of a non-empty list reads past the end in giac
    // (undefined: 0 in a fresh process, SIGSEGV in this shared one, CI
    // 45112693111): no giac row, ours in cas-core's unit tests.
    ("heapify([3,1,4,1,5,9,2,6])", "1", "[9,6,4,1,5,3,2,1]"),
    ("heapify([])", "1", "[]"),
    ("heappush([5,3],x)", "[x,3,5]", "heappush([5,3],x)"),
    ("heapify([x,y])", "1", "heapify([x,y])"),
    // giac's errors on heaps and permutations (CAS-560, CAS-561); ours
    // leaves the call, with the error as a warning.
    (
        "heappush(5,1)",
        "GIAC_ERROR: heappush(5,1) \n Error: Bad Argument Value",
        "heappush(5,1)",
    ),
    (
        "heappush([1])",
        "GIAC_ERROR: heappush([1]) \n Error: Bad Argument Value",
        "heappush([1])",
    ),
    (
        "heapify(5)",
        "GIAC_ERROR: heapify(5) \n Error: Bad Argument Value",
        "heapify(5)",
    ),
    (
        "heappop(5)",
        "GIAC_ERROR: heappop(5) \n Error: Bad Argument Value",
        "heappop(5)",
    ),
    (
        "heappop([])",
        "GIAC_ERROR: heappop([]) \n Error: Invalid dimension",
        "heappop([])",
    ),
    (
        "nextperm([1,2,3])",
        "GIAC_ERROR: nextperm([1,2,3]) \n Error: Bad Argument Value",
        "nextperm([1,2,3])",
    ),
    (
        "nextperm([1,1,0])",
        "GIAC_ERROR: nextperm([1,1,0]) \n Error: Bad Argument Value",
        "nextperm([1,1,0])",
    ),
    (
        "nextperm(5)",
        "GIAC_ERROR: nextperm(5) \n Error: Bad Argument Type",
        "nextperm(5)",
    ),
    (
        "prevperm(x)",
        "GIAC_ERROR: prevperm(x) \n Error: Bad Argument Type",
        "prevperm(x)",
    ),
    (
        "groupermu([1,1],[0,1])",
        "GIAC_ERROR: groupermu([1,1],[0,1]) \n Error: Bad Argument Value",
        "groupermu([1,1],[0,1])",
    ),
    (
        "groupermu([1,0,2,3])",
        "GIAC_ERROR: groupermu([1,0,2,3]) \n Error: Bad Argument Type",
        "groupermu([1,0,2,3])",
    ),
    (
        "groupermu([1,2,0])",
        "GIAC_ERROR: groupermu([1,2,0]) \n Error: Bad Argument Type",
        "groupermu([1,2,0])",
    ),
    (
        "groupermu([1,0,2],[0,2,1],[1,2,0])",
        "GIAC_ERROR: groupermu([1,0,2],[0,2,1],[1,2,0]) \n Error: Bad Argument Type",
        "groupermu([1,0,2],[0,2,1],[1,2,0])",
    ),
    (
        "groupermu(1,[0])",
        "GIAC_ERROR: groupermu(1,[0]) \n Error: Bad Argument Type",
        "groupermu(1,[0])",
    ),
    (
        "groupermu([1,0],[0,2])",
        "GIAC_ERROR: groupermu([1,0],[0,2]) \n Error: Bad Argument Value",
        "groupermu([1,0],[0,2])",
    ),
];

/// The departures of [`DEPARTURES`], then `newList` beyond this CAS's
/// 65536 elements and `apply` adding 1 to the rows of a matrix.
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
    // giac adds 1 to a vector as a polynomial's coefficients; ours refuses a
    // list plus a number: an error.
    let (input, answer) = ("apply(x->x+1,[[1,2],[3,4]])", "[poly1[1,3],poly1[3,5]]");
    let printed = giac.eval(input).unwrap();
    if printed != answer {
        wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
    }
    if ours(input).is_ok() {
        wrong.push(format!("{input}: ours is not an error"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A small integer or rational, as text: `k` or `k/d`.
fn number() -> impl Strategy<Value = String> {
    (-5i64..=5, 1i64..=4).prop_map(|(k, d)| {
        if d == 1 || k % d == 0 {
            (k / d).to_string()
        } else {
            format!("{k}/{d}")
        }
    })
}

/// Up to `max` numbers.
fn numbers(min: usize, max: usize) -> impl Strategy<Value = Vec<String>> {
    prop::collection::vec(number(), min..=max)
}

/// A `rows` by `cols` matrix of numbers, rows and columns 1 to 4.
fn matrix() -> impl Strategy<Value = Vec<Vec<String>>> {
    (1usize..=4, 1usize..=4)
        .prop_flat_map(|(r, c)| prop::collection::vec(prop::collection::vec(number(), c), r))
}

/// `[a,b,…]`.
fn list<T: ToString>(items: &[T]) -> String {
    let items: Vec<String> = items.iter().map(T::to_string).collect();
    format!("[{}]", items.join(","))
}

/// `[[…],…]`.
fn rows(m: &[Vec<String>]) -> String {
    let rows: Vec<String> = m.iter().map(|r| list(r)).collect();
    format!("[{}]", rows.join(","))
}

/// giac's answer to `input`, read back, is ours; else the failure names
/// the input, giac's answer and ours.
fn agree(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    giac.agrees(input)
        .map_err(|e| TestCaseError::fail(e.to_string()))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `find` and `index` of an element of the list (`index` of a missing
    /// one is giac's error, a departure).
    #[test]
    fn random_positions(l in numbers(1, 6), k in 0usize..6) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let v = &l[k % l.len()];
        agree(giac, &format!("find({v},{})", list(&l)))?;
        agree(giac, &format!("index({},{v})", list(&l)))?;
    }

    /// `count_eq`, `count_inf` and `count_sup` of a list, and of a matrix
    /// flattened, by `row` or by `col` (`row` of a vector is giac's error).
    #[test]
    fn random_counts(
        v in number(),
        l in numbers(0, 6),
        m in matrix(),
        how in 0usize..3,
        name in prop::sample::select(vec!["count_eq", "count_inf", "count_sup"]),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("{name}({v},{})", list(&l)))?;
        let by = ["", ",row", ",col"][how];
        agree(giac, &format!("{name}({v},{}{by})", rows(&m)))?;
    }

    /// `deltalist` of two or more, `accumulate_head_tail` with head and
    /// tail at least 1 and fitting the list (else giac's errors,
    /// departures), and `sizes` of lists.
    #[test]
    fn random_differences_and_sums(
        l in numbers(2, 6),
        h in 1usize..=5,
        t in 1usize..=5,
        ls in prop::collection::vec(numbers(0, 3), 0..=4),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("deltalist({})", list(&l)))?;
        if h + t <= l.len() {
            agree(giac, &format!("accumulate_head_tail({},{h},{t})", list(&l)))?;
        }
        let ls: Vec<String> = ls.iter().map(|l| list(l)).collect();
        agree(giac, &format!("sizes({})", list(&ls)))?;
    }

    /// `makevector`, `newList`, `enumerate` from an integer start and
    /// `copy`.
    #[test]
    fn random_builders(l in numbers(0, 6), n in -2i64..=6, k in -3i64..=3) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("makevector({})", l.join(",")))?;
        agree(giac, &format!("newList({n})"))?;
        agree(giac, &format!("enumerate({},{k})", list(&l)))?;
        agree(giac, &format!("copy({})", list(&l)))?;
    }

    /// `apply` of a power or a shift to a list, and to a matrix with
    /// `matrix`; `foldl` and `foldr` of a difference.
    #[test]
    fn random_maps_and_folds(
        l in numbers(0, 6),
        m in matrix(),
        // From 1: giac's 0^0 is undef, a power's convention, not apply's.
        e in 1u32..=3,
        c in number(),
        init in number(),
        args in numbers(1, 5),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("apply(x->x^{e},{})", list(&l)))?;
        agree(giac, &format!("apply(x->x+{c},{},matrix)", rows(&m)))?;
        let args = args.join(",");
        agree(giac, &format!("foldl((x,y)->x-y,{init},{args})"))?;
        agree(giac, &format!("foldr((x,y)->x-y,{init},{args})"))?;
    }

    /// `multiply` of vectors (over the shorter), of matrices that fit,
    /// of a matrix and a vector of its columns' or rows' length (a vector
    /// longer than the rows is a departure), and of a scalar.
    #[test]
    fn random_products(
        u in numbers(1, 4),
        w in numbers(1, 4),
        (a, b) in (1usize..=3, 1usize..=3, 1usize..=3).prop_flat_map(|(r, k, c)| (
            prop::collection::vec(prop::collection::vec(number(), k), r),
            prop::collection::vec(prop::collection::vec(number(), c), k),
        )),
        s in number(),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("multiply({},{})", list(&u), list(&w)))?;
        agree(giac, &format!("multiply({},{})", rows(&a), rows(&b)))?;
        let column: Vec<String> = b.iter().map(|r| r[0].clone()).collect();
        agree(giac, &format!("multiply({},{})", rows(&a), list(&column)))?;
        agree(giac, &format!("multiply({},{})", list(&a[0]), rows(&b)))?;
        agree(giac, &format!("multiply({s},{})", rows(&a)))?;
    }

    /// `insert` at a position within the list, `subsop` of an index from
    /// `-len` to `len-1` by a value or `NULL` (outside: giac's errors,
    /// departures).
    #[test]
    fn random_insertions_and_substitutions(
        l in numbers(1, 6),
        p in 0usize..=6,
        k in -6i64..6,
        v in number(),
        delete in any::<bool>(),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let len = i64::try_from(l.len()).expect("a short list");
        let p = p % (l.len() + 1);
        agree(giac, &format!("insert({},{p},{v})", list(&l)))?;
        let k = k.rem_euclid(2 * len) - len;
        let v = if delete { "NULL" } else { v.as_str() };
        agree(giac, &format!("subsop({},{k}={v})", list(&l)))?;
    }

    /// `heappush` onto any list of numbers (libstdc++'s `push_heap` sifts
    /// up whether or not it is a heap), by `<` or by `>`.
    #[test]
    fn random_heappush(l in numbers(0, 6), v in number(), greater in any::<bool>()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let by = if greater { ",(a,b)->a>b" } else { "" };
        agree(giac, &format!("heappush({},{v}{by})", list(&l)))?;
    }

    /// `nextperm` and `prevperm` of a permutation of `0..n`, and
    /// `groupermu` of two of the same length (other lists: giac's errors,
    /// departures).
    #[test]
    fn random_permutations(
        (p, q) in (1usize..=6).prop_flat_map(|n| {
            let perm = Just((0..n).collect::<Vec<_>>()).prop_shuffle();
            (perm.clone(), perm)
        }),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("nextperm({})", list(&p)))?;
        agree(giac, &format!("prevperm({})", list(&p)))?;
        agree(giac, &format!("groupermu({},{})", list(&p), list(&q)))?;
    }
}
