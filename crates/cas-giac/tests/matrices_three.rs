//! giac's matrix constructors, slices and row operations (CAS-531,
//! CAS-532, CAS-533) against giac 2.0: its answers read back are this
//! CAS's, and the departures, each with giac's answer as it prints it and
//! this CAS's as decided.

use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac's answers to `inputs`, read back, are this CAS's.
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

#[test]
fn constructors() {
    check(&[
        "matrix(3,2,(j,k)->j+k)",
        "matrix(2,3)",
        "matrix(2,2,7)",
        "matrix(2,3,[1,2,3,4,5,6])",
        "matrix(2,2,[1,2,3])",
        "matrix(2,2,[1,2,3,4,5])",
        "matrix(2,2,(j,k)->j*k+1)",
        "makemat((j,k)->j+k,3,2)",
        "makemat(2,2,3)",
        "list2mat([1,8,4,9],1)",
        "list2mat([1,8,4,9,5],2)",
        "list2mat([],2)",
        "mat2list([[1,8],[4,9]])",
        "mat2list([1,2,3])",
        "mat2list([[1,[2,3]],[4,5]])",
        "blockmatrix(2,3,[idn(2),idn(2),idn(2),idn(2),idn(2),idn(2)])",
        "blockmatrix(2,1,[[[1,2]],[[3,4]]])",
        "blockmatrix(1,2,[[[1],[2]],[[3],[4]]])",
        "JordanBlock(7,3)",
        "JordanBlock(a,1)",
        "BlockDiagonal([[1,2],[3,4]])",
        "BlockDiagonal([1,2,3])",
        "BlockDiagonal([[[1]],[[2,3],[4,5]]])",
    ]);
}

#[test]
fn slices() {
    check(&[
        "rowdim([[1,2,3],[4,5,6]])",
        "nrows([[1,2,3],[4,5,6]])",
        "coldim([[1,2,3],[4,5,6]])",
        "row([[1,2],[3,4],[5,6]],0)",
        "row([[1,2],[3,4],[5,6]],-1)",
        "row([1,2,3],0)",
        "col([[1,2,3],[4,5,6],[7,8,9]],1)",
        "col([[1,2],[3,4]],-1)",
        "subMat([[1,2],[3,4],[5,6]],1,0,2,1)",
        "subMat([[1,2,3],[4,5,6],[7,8,9]],1,1)",
        "subMat([[1,2,3],[4,5,6],[7,8,9]])",
        "delrows([[1,2],[3,4],[5,6]],1..2)",
        "delrows([[1,2,3],[4,5,6],[7,8,9]],1)",
        "delrows([[1,2],[3,4]],0..1)",
        "delcols([[1,2,3],[4,5,6],[7,8,9]],1)",
        "delcols([[1,2,3],[4,5,6],[7,8,9]],0..1)",
    ]);
}

#[test]
fn row_operations() {
    check(&[
        "mRow(12,[[1,2],[3,4],[5,6]],0)",
        "mRow(2,[[1,2],[3,4],[5,6]],0)",
        "mRowAdd(12,[[1,2],[3,4],[5,6]],0,2)",
        "mRowAdd(x,[[1,2],[3,4]],1,0)",
        "scaleadd([[1,2],[3,4],[5,6]],12,0,2)",
        "SCALEADD([[1,2],[3,4],[5,6]],12,0,2)",
        "scale([[1,2],[3,4],[5,6]],12,0)",
        "SCALE([[1,2],[3,4],[5,6]],12,0)",
        "rowAdd([[1,2],[3,4],[5,6]],1,2)",
        "rowAdd([[1,2],[3,4]],0,0)",
        "rowSwap([[1,2],[3,4],[5,6]],1,2)",
        "rowSwap([[1,2],[3,4]],0,0)",
        "swaprow([[1,2],[3,4],[5,6]],0,2)",
        "colSwap([[1,2],[3,4],[5,6]],0,1)",
        "swapcol([[1,2],[3,4],[5,6]],0,1)",
    ]);
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it (an error by its message), and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // A list of rows or columns for a range, where giac's is a sequence.
        (
            "row([[1,2],[3,4],[5,6]],1..2)",
            "[3,4],[5,6]",
            "[[3,4],[5,6]]",
        ),
        (
            "col([[1,2,3],[4,5,6],[7,8,9]],0..1)",
            "[1,4,7],[2,5,8]",
            "[[1,4,7],[2,5,8]]",
        ),
        // At least one row and one column, where giac takes 0 for 1.
        ("makemat(1,0,2)", "[[1,1]]", "makemat(1,0,2)"),
        (
            "list2mat([1,2,3],0)",
            "[[1],[2],[3]]",
            "list2mat([1,2,3],0)",
        ),
        // Reversed bounds stay.
        (
            "subMat([[1,2,3],[4,5,6],[7,8,9]],2,2,0,0)",
            "[[9]]",
            "subMat([[1,2,3],[4,5,6],[7,8,9]],2,2,0,0)",
        ),
        // A negative row counts from the end, as every index does.
        (
            "mRow(2,[[1,2],[3,4]],-1)",
            "Error: Invalid dimension",
            "[[1,2],[6,8]]",
        ),
        // Calls that do not apply stay, where giac errs.
        (
            "rowDim([1,2,3])",
            "Error: Bad Argument Value",
            "rowDim([1,2,3])",
        ),
        (
            "JordanBlock(2,0)",
            "Error: Invalid dimension",
            "JordanBlock(2,0)",
        ),
        (
            "BlockDiagonal([[[1,2]],[[3],[4]]])",
            "Error: Bad Argument Type",
            "BlockDiagonal([[[1,2]],[[3],[4]]])",
        ),
        (
            "blockmatrix(2,2,[1,2,3,4])",
            "Error: Bad Argument Type",
            "blockmatrix(2,2,[1,2,3,4])",
        ),
    ];
    let normal = |s: &str| s.replace("matrix[", "[").replace(' ', "");
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        let recorded = if answer.starts_with("Error") {
            printed.contains(answer)
        } else {
            normal(&printed) == normal(answer)
        };
        if !recorded {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
