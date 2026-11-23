//! Xcas's list commands against giac (CAS-555 to CAS-561): `find` and
//! `index` from 0, `count_eq`, `count_inf` and `count_sup` with `row` and
//! `col`, `deltalist`, `accumulate_head_tail`, `sizes`, `makevector`,
//! `makesuite`, `newList`, `enumerate`, `copy`, `apply`, `foldl`, `foldr`,
//! `multiply`, `insert`, `subsop`, `heappush`, `nextperm`, `prevperm` and
//! `groupermu`, giac's answers read back
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
const DEPARTURES: [(&str, &str, &str); 69] = [
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
    // Heaps (CAS-560): giac works in place and prints 1 for heapify, the
    // element past the end for heappop; ours returns the heap. giac orders
    // symbols where `Less` does not decide: ours stays.
    ("heapify([3,1,4,1,5,9,2,6])", "1", "[9,6,4,1,5,3,2,1]"),
    ("heapify([])", "1", "[]"),
    ("heappop([5,3,4,1])", "0", "[4,3,1]"),
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
