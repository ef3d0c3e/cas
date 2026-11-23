//! Xcas's calculus helpers against giac (CAS-574, CAS-577, CAS-578):
//! `ibpdv` and `ibpu` (one step, the bounds, the finish), `preval`, and
//! `function_diff` (`fonction_derivee`) applied, giac's answers read back
//! and compared with ours (`Giac::agrees`); and the departures, each with
//! giac's answer as libgiac prints it.

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

/// `ibpdv` and `ibpu` (CAS-574): giac's examples.
#[test]
fn integration_by_parts() {
    check(&[
        "ibpdv(ln(x),x)",
        "ibpdv(ln(x),x,x,1,3)",
        "ibpdv(x*ln(x),x^2/2)",
        "ibpdv([x*ln(x),-1],0)",
        "ibpdv(ibpdv(ln(x),x,x,2,3),0,x,2,3)",
        "ibpdv(ln(y),y,y)",
        "ibpdv(x*exp(x),exp(x))",
        "ibpu(ln(x),ln(x))",
        "ibpu(ln(x),ln(x),x,1,3)",
        "ibpu(x*ln(x),ln(x))",
        "ibpu([x*ln(x),-1],0)",
        "ibpu(ibpu(ln(x),ln(x),x,2,3),0,x,2,3)",
        "ibpu(x^2,1)",
        "ibpu(x*exp(x),x)",
        "ibpu(exp(x)*sin(x),sin(x))",
    ]);
}

/// `preval` (CAS-578).
#[test]
fn values_at_bounds() {
    check(&[
        "preval(x^2-2,2,3)",
        "preval(y^2-2,2,3,y)",
        "preval(int(y,y),0,1,y)",
        "preval(x^2,x=1,2)",
        "preval(sin(x)/x,0,pi)",
        "preval(1/x,-1,1)",
    ]);
}

/// `function_diff` (CAS-577), applied: giac's functions print with a
/// hidden parameter `` ` x` `` that the language cannot read.
#[test]
fn derivatives_of_functions() {
    check(&[
        "function_diff(sin)(x)",
        "fonction_derivee(sin)(x)",
        "function_diff(sin+cos)(t)",
        "function_diff(2*exp)(t)",
        "function_diff(ln)(x)",
        "function_diff(ln)([x,y])",
        "function_diff(sin)([x,y])",
        "function_diff(unapply(y^3,y))(2)",
        "function_diff(f)",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac's `int` takes `x` when no variable is given; ours needs one.
        ("preval(int(x),0,1)", "1/2", "int(1)-int(0)"),
        // giac's `id` is the identity; ours has none, so `id` is an
        // unknown function.
        (
            "function_diff(sin+id)",
            "(` x`)->cos(` x`)+1",
            "function_diff(sin+id)",
        ),
        // A function of one parameter: ours ignores the arguments past it,
        // giac multiplies.
        ("function_diff(ln)(x,y)", "1/x/y", "1/x"),
        // giac's hidden parameter; ours is `x`.
        ("function_diff(sin)", "(` x`)->cos(` x`)", "(x)->cos(x)"),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        match (ours(input), ours(decided)) {
            (Ok(found), Ok(decided)) if found == decided => {}
            (found, _) => wrong.push(format!("{input}: ours {found:?}, not {decided}")),
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac's `@` composition and `@@` powers: the language has neither.
    for (input, answer) in [
        (
            "function_diff(sq@sin+id)",
            "(` x`)->2*cos(` x`)*sin(` x`)+1",
        ),
        ("(function_diff @@3)(ln)('x')", "2/x^3"),
    ] {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input).is_ok() {
            wrong.push(format!("{input}: ours reads it now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
