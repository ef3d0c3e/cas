//! Set operations against giac (CAS-DOC-32 §4, CAS-PLAN-26): giac's
//! `set[…]` answers, read back as `{…}`, are ours exactly. giac writes
//! sets as `set[…]` and this CAS as `{…}`, so each case has both inputs.
//! Not `is_included`, which giac 2.0.0 answers by its arguments' sizes:
//! `is_included([3],[1,2])` is 1 there, and `False` here.

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
fn set_operations() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, theirs_input) in [
        ("{1,2} union {2,3}", "set[1,2] union set[2,3]"),
        ("[1,2] union [2,3]", "[1,2] union [2,3]"),
        ("{1,2} intersect {2,3}", "set[1,2] intersect set[2,3]"),
        ("[1,2] intersect [2,3]", "[1,2] intersect [2,3]"),
        ("{1,2,3} minus {2}", "set[1,2,3] minus set[2]"),
        ("[1,2,3] minus [2]", "[1,2,3] minus [2]"),
        ("{1,2} union [5]", "set[1,2] union [5]"),
        ("{x,y} union {y,z}", "set[x,y] union set[y,z]"),
        ("{} union {1}", "set[] union set[1]"),
        ("size({1,2,2})", "size(set[1,2,2])"),
    ] {
        match (giac.value(theirs_input), ours(input)) {
            (Ok(theirs), Ok(ours)) if theirs == ours => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
