//! Xcas's `type` and its `DOM_` answers (CAS-621), `iPart`, `fPart`, `frac`,
//! `trunc`, `SortA`, `SortD`, `sortperm` and `RandSeed` (CAS-622) against
//! giac: giac's answers read back are this CAS's value of the same input,
//! exact answers against giac's decimals
//! within their error (C0006, C0032); departures, giac's answer as caseval
//! prints it, whitespace collapsed (C0047), and this CAS's own.

use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Each input: giac's answer, read back, is this CAS's value.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `type`, compared through `==` since giac prints its answer as the
/// type's name (probe1, probe2, probe6).
#[test]
fn types() {
    check(&[
        "type(2)==DOM_INT",
        "type(2^100)==DOM_INT",
        "type(1/2)==DOM_RAT",
        "type(2.0)==DOM_FLOAT",
        "type(1+i)==DOM_COMPLEX",
        "type(1.0+i)==DOM_COMPLEX",
        "type(x)==DOM_IDENT",
        "type(pi)==DOM_IDENT",
        "type(infinity)==DOM_IDENT",
        "type(undef)==DOM_IDENT",
        "type(\"a\")==DOM_STRING",
        "type([1,2])==DOM_LIST",
        "type(sin)==DOM_FUNC",
        "type(x->x^2)==DOM_FUNC",
        "type(sqrt(2))==DOM_SYMBOLIC",
        "type(x=1)==DOM_SYMBOLIC",
    ]);
}

/// `iPart`, `fPart` and `frac`; a sequence bracketed (C0023; probe2,
/// probe6).
#[test]
fn parts() {
    check(&[
        "iPart(4.3)",
        "iPart(sqrt(2))",
        "iPart(-4.3)",
        "iPart(7/2)",
        "iPart(-7/2)",
        "[iPart(4.3,sqrt(2))]",
        "iPart([1.5,-2.5])",
        "fPart(1/2)",
        "fPart(-1/2)",
        "fPart(1.2)",
        "fPart(-1.2)",
        "fPart(7/2)",
        "fPart(-2.5)",
        "fPart(sqrt(2))",
        "fPart([3.4,sqrt(2)])",
        "frac(7/3)",
    ]);
}

/// `trunc`, to an integer or to `n` decimals (probe2, probe6).
#[test]
fn truncation() {
    check(&[
        "trunc(4.3)",
        "trunc(-4.7)",
        "trunc(7/2)",
        "trunc(-7/2)",
        "trunc(sqrt(2),3)",
        "trunc(pi,4)",
        "trunc(1/3,2)",
        "trunc(2/3,1)",
        "trunc(12345.678,-2)",
        "trunc([4.3333,sqrt(2)])",
        "trunc([4.3333,sqrt(2)],2)",
        "trunc(sqrt(2)+i*sqrt(5),4)",
    ]);
}

/// `SortA`, `sorta`, `SortD`, `sortd` and `sortperm` (CAS-622): lists,
/// sequences, matrices and sequences of rows by their first row, ties
/// (`SortD` reverses them, as giac), symbols and strings; `RandSeed`, `srand`
/// (c2probe1, c2probe2).
#[test]
fn sorting() {
    check(&[
        "SortA([3,4,2])",
        "SortA(3,4,2)",
        "SortA([[3,4,2],[6,4,5]])",
        "SortA([3,4,2],[6,4,5])",
        "SortD([3,4,2])",
        "SortD([[3,4,2],[6,4,5]])",
        "sorta([3,1,2],[6,5,4])",
        "sortd(3,4,2)",
        "SortA([x,1,y])",
        "SortD([x,1,y])",
        "SortD([\"a\",\"c\",\"b\"])",
        "SortD([1/2,0.3,1])",
        "SortA([2,1,2,1])",
        "SortA([[2,1,2],[3,2,1]])",
        "SortD([[2,1,2],[3,2,1]])",
        "SortA([[x,3,1],[1,2,3]])",
        "SortA([[3,1],[6,5],[9,8]])",
        "SortA([5])",
        "sortperm([30,25,40,10,20])",
        "sortperm([3,1,2,1])",
        "sortperm([\"b\",\"a\",\"c\"])",
        "sortperm([5])",
        "sortperm([\"ab\",\"cd\",\"ad\",\"bc\",\"de\"],[3,4,1,0,2])",
        "sortperm([3,1,2],[2,0,1])",
        "RandSeed(321)",
        "RandSeed(5)==srand(5)",
    ]);
}

/// Departures, giac's answer as caseval prints it (probe1, probe2, probe6)
/// and ours: giac's types are integers, printed as the type's name, and
/// `true` is one; giac writes `fPart` and `trunc` of a symbol with `floor`,
/// where ours stays; giac's `trunc(x, n)` floors.
#[test]
fn departures() {
    let cases = [
        ("DOM_INT+0", "2", "DOM_INT"),
        ("type(2)", "integer", "DOM_INT"),
        ("type(true)", "integer", "DOM_IDENT"),
        ("fPart(x)", "x-sign(x)*floor(x*sign(x))", "fPart(x)"),
        ("trunc(x)", "floor(x)", "trunc(x)"),
        ("trunc(-sqrt(2),3)", "-1.415", "-707/500"),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, mine) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours(input);
        if found.as_ref().ok() != ours(mine).ok().as_ref() {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
