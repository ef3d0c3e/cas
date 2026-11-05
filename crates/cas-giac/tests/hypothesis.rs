//! The hypothesis tests against giac (CAS-DOC-52, CAS-509), where the
//! conventions agree: giac's answers read back, entry by entry, its
//! decimals within 10⁻¹¹. giac takes its alternatives as quoted operators
//! (`'<'`) and this language as strings (`"<"`), so each side may have its
//! own input; giac's `poly1` laws are read through `[op(…)]`.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The leaves of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// Each `(giac's input, ours)`: the same values, decimals within 10⁻¹¹.
fn check(pairs: &[(&str, &str)]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (theirs_input, our_input) in pairs {
        let (Ok(theirs), Ok(mine)) = (giac.value(theirs_input), ours(our_input)) else {
            wrong.push(format!("{theirs_input} / {our_input}: not read"));
            continue;
        };
        let (t, o) = (leaves(&theirs), leaves(&mine));
        let same = t.len() == o.len()
            && t.iter().zip(&o).all(|(a, b)| {
                if a == b {
                    return true;
                }
                let numeric = |e: &Expr| {
                    ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
                        .unwrap()
                };
                let gap = Expr::call(Symbol::ABS, [numeric(a) - numeric(b)]);
                let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-11)]);
                ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
            });
        if !same {
            wrong.push(format!("{our_input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input the same on both sides.
fn same(inputs: &[&str]) {
    let pairs: Vec<(&str, &str)> = inputs.iter().map(|s| (*s, *s)).collect();
    check(&pairs);
}

/// `normalt` and `studentt` (CAS-503).
#[test]
fn location() {
    check(&[
        (
            "normalt([10,30],0.5,0.02,'!=',0.1)",
            "normalt([10,30],0.5,0.02,\"!=\",0.1)",
        ),
        (
            "normalt([0.48,50],0.5,0.1,'<')",
            "normalt([0.48,50],0.5,0.1,\"<\")",
        ),
        (
            "normalt([10,30],0.5,0.02,'>')",
            "normalt([10,30],0.5,0.02,\">\")",
        ),
        (
            "studentt([0.48,20],0.5,0.1,'<')",
            "studentt([0.48,20],0.5,0.1,\"<\")",
        ),
        (
            "studentt([1,2,3,4,5],2,'>')",
            "studentt([1,2,3,4,5],2,\">\")",
        ),
        (
            "studentt([1,2,3,4,5],0,'>')",
            "studentt([1,2,3,4,5],0,\">\")",
        ),
    ]);
}

/// `chisquaret` (CAS-504).
#[test]
fn chi_square() {
    same(&[
        "chisquaret([57,30],[0.6,0.4])",
        "chisquaret([57,30],[0.6,0.4],0.01)",
        "chisquaret([1,1,1,5,2,4])",
        "chisquaret([57,54])",
        "chisquaret([17,15,12,15],[15,13,13,14])",
        "chisquaret([11,16,17,22,14,10],[1/6,1/6,1/6,1/6,1/6,1/6])",
        "chisquaret([1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,0,0],[0.4,0.6])",
        "chisquaret([[17,15,12,15],[15,13,13,14],[1,2,3,4]])",
        "chisquaret([57,30],[3,2])",
    ]);
}

/// `kolmogorovd` and `kolmogorovt` (CAS-505, CAS-506); giac's
/// `kolmogorovd` loses its digits below about 10⁻¹⁶.
#[test]
fn kolmogorov() {
    same(&[
        "kolmogorovd(1.0)",
        "kolmogorovd(0.5)",
        "kolmogorovd(0.9999)",
        "kolmogorovd(1.0001)",
        "kolmogorovd(3.0)",
        "kolmogorovt([1,2,3],[1.5,2.5,8])",
        "kolmogorovt([1,2,3],[3/2,5/2,8])",
        "kolmogorovt([0.1,0.5,0.9],normald(0,1))",
        "kolmogorovt([0.1,0.5,0.9],uniformd(0,2))",
        "kolmogorovt([0.3,1.2,2.5,0.7],exponentiald(1))",
    ]);
}

/// `wilcoxonp`, `wilcoxons` and `wilcoxont` (CAS-507, CAS-508), on data
/// without ties or zero differences, and two-sided or one-sample where the
/// alternatives agree.
#[test]
fn wilcoxon() {
    check(&[
        ("[op(wilcoxonp(4))]", "wilcoxonp(4)"),
        ("[op(wilcoxonp(7))]", "wilcoxonp(7)"),
        ("[op(wilcoxonp(3,2))]", "wilcoxonp(3,2)"),
        ("[op(wilcoxonp(4,5))]", "wilcoxonp(4,5)"),
        (
            "wilcoxont([1,2,3,4,6,7,8,9,11,13],2.5,'>')",
            "wilcoxont([1,2,3,4,6,7,8,9,11,13],2.5,\">\")",
        ),
        (
            "wilcoxont([1,2,3,4,6,7,8,9,11,13],8.5,'<')",
            "wilcoxont([1,2,3,4,6,7,8,9,11,13],8.5,\"<\")",
        ),
    ]);
    same(&[
        "wilcoxons([1,3,4,20],[2,6,10])",
        "wilcoxons([5,6,7],[1,2,3,4])",
        "wilcoxons([1,2,9,7],5.25)",
        "wilcoxont([10,30,20,40],[1,2,3,4])",
        "wilcoxont([1,3,4,20],[2,6,10])",
        "wilcoxont([1,3,4,5,7,8,8,12,15,17],[2,6,10,11,13,14,15,18,19,20])",
        "wilcoxont([1,3,4,5,7,8,8,12,15,17],[2,6,10,11,13,14,15,18,19,20],0.01)",
        "wilcoxont([1,2,3,4,6,7,8,9,11,13],2.5)",
    ]);
}
