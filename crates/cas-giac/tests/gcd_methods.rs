//! The gcd algorithms by name against giac (CAS-526): `heugcd`, `modgcd`,
//! `psrgcd`, `ezgcd` and `icontent`, giac's answers read back; the
//! departures, each as decided: `psrgcd`'s leading coefficient positive,
//! and `ezgcd` where giac gives up or refuses one variable.

use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac 2.0's answers, read back, are this CAS's.
#[test]
fn agreement() {
    let Some(giac) = giac() else {
        return;
    };
    let inputs = [
        "heugcd(x^2-1,x^2+2*x+1)",
        "modgcd(x^2-1,x^2+2*x+1)",
        "heugcd(6*x^2-6,4*x+4)",
        "modgcd(x*y-y,x^2*y-y)",
        "heugcd(x*y-y,x^2*y-y)",
        "psrgcd(2*x^2-2,4*x+4)",
        "psrgcd(x^3-1,x^2-1)",
        "psrgcd(x^2*y-y,x*y-y)",
        "ezgcd((x+y)*(x-y+1),(x+y)*(x+2*y))",
        "ezgcd(2*x*y-2,4*x*y+4*y)",
        "heugcd(12,18)",
        "icontent(6*x^2+4*x+2)",
        "icontent(6*x*y+4*y)",
        "icontent(3*x/2+3)",
        "icontent(-6*x-4)",
        "icontent(12)",
    ];
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it prints
/// it (an error as its message), and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // The leading coefficient positive, where psrgcd keeps the sign of
        // its remainders.
        ("psrgcd(x^2-1,x^2+2*x+1)", "-x-1", "x+1"),
        // EZ-GCD gives up where giac's does not, or Brown's in one variable.
        (
            "ezgcd(x^2*y^2-1,x*y+1)",
            "GCD not successful Error: Bad Argument Value",
            "x*y+1",
        ),
        (
            "ezgcd(x^2*y-y,x*y^2-y^2)",
            "GCD not successful Error: Bad Argument Value",
            "x*y-y",
        ),
        (
            "ezgcd(x^2-1,x+1)",
            "Args must be multivariate polynomials Error: Bad Argument Value",
            "x+1",
        ),
    ];
    // giac prints an error as a string, quoted.
    let unquoted = |s: &str| s.trim_matches('"').to_owned();
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if unquoted(&printed) != answer {
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
