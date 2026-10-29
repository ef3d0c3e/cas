//! Integer matrices against giac (CAS-DOC-42, CAS-PLAN-36): giac's
//! Hermite and Smith forms, which are unique, are this CAS's entry by
//! entry; its `lll`'s reduced basis spans this CAS's lattice. giac answers
//! sequences, read in brackets.

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

const MATRICES: [&str; 5] = [
    "[[9,-36,30],[-36,192,-180],[30,-180,180]]",
    "[[2,4],[1,3]]",
    "[[2,4],[6,8]]",
    "[[1,2,3],[4,5,6],[7,8,10]]",
    "[[6,4,2],[3,9,12],[1,5,7]]",
];

/// giac's `[command(A)]`'s part at `index` and this CAS's `command(A)`'s.
fn check_part(command: &str, index: usize) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for a in MATRICES {
        let input = format!("{command}({a})");
        let theirs = giac.value(&format!("[{input}]"));
        let ours = ours(&input);
        match (theirs, ours) {
            (Ok(theirs), Ok(ours))
                if theirs.args().len() > index
                    && ours.args().len() > index
                    && theirs.args()[index] == ours.args()[index] => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Hermite normal forms (CAS-445).
#[test]
fn hermite() {
    check_part("ihermite", 1);
}

/// Smith normal forms (CAS-446).
#[test]
fn smith() {
    check_part("ismith", 1);
}

/// The Hermite form of the rows of `m`, by this CAS.
fn lattice(ctx: &Context, m: &Expr) -> Expr {
    let h = ctx
        .eval(&Expr::call(Symbol::HERMITE_DECOMPOSITION, [m.clone()]))
        .unwrap();
    h.args()[1].clone()
}

/// LLL (CAS-447): giac's reduced basis and this CAS's span one lattice.
#[test]
fn lll() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for b in [
        "[[1,1,1],[-1,0,2],[3,5,6]]",
        "[[1,0,0,12345],[0,1,0,23456],[0,0,1,34567]]",
        "[[105,821,404,328],[881,667,644,927],[181,483,87,500],[893,834,732,441]]",
    ] {
        let input = format!("lll({b})");
        let theirs = giac.value(&format!("[{input}]"));
        let ours = ours(&input);
        match (theirs, ours) {
            (Ok(theirs), Ok(ours))
                if !theirs.args().is_empty()
                    && lattice(&ctx, &theirs.args()[0]) == lattice(&ctx, &ours) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
