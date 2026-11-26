//! The gcd algorithms by name against giac (CAS-526): `heugcd`, `modgcd`,
//! `psrgcd`, `ezgcd` and `icontent`, giac's answers read back; the
//! departures, each as decided: `psrgcd`'s leading coefficient positive,
//! and `ezgcd` where giac gives up or refuses one variable; and 128 random
//! inputs a run to each.

use std::collections::BTreeMap;

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

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
        // giac's EZ-GCD gives up where ours answers; in one variable giac
        // rejects the call, ours runs Brown's algorithm.
        (
            "ezgcd(x^2*y^2-1,x*y+1)",
            "GIAC_ERROR: GCD not successful Error: Bad Argument Value",
            "x*y+1",
        ),
        (
            "ezgcd(x^2*y-y,x*y^2-y^2)",
            "GIAC_ERROR: GCD not successful Error: Bad Argument Value",
            "x*y-y",
        ),
        // giac's modgcd in one variable drops the integer content of the
        // gcd (its own gcd keeps it); ours is gcd's answer. T7's sweep
        // found it: shrunk input modgcd(((-2))*((-2)*x),((-2))*((-2)*x)).
        ("modgcd(4*x,4*x)", "x", "4*x"),
        ("modgcd(4*x,6*x)", "x", "2*x"),
        // In y alone too; T9: shrunk input modgcd(((-2))*(2*y),((-1))*(2*y)).
        ("modgcd(6*y,4*y)", "y", "2*y"),
        (
            "ezgcd(x^2-1,x+1)",
            "GIAC_ERROR: Args must be multivariate polynomials Error: Bad Argument Value",
            "x+1",
        ),
    ];
    // caseval prints an error as `GIAC_ERROR: <message>`; quotes dropped.
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

const VARIABLES: [&str; 2] = ["x", "y"];

/// The terms of a polynomial in `x, y`: exponents to coefficient.
type Terms = BTreeMap<Vec<u32>, i64>;

/// Up to `len` terms with exponents below 3 and coefficients in `-5..=5`,
/// none 0, so a nonzero polynomial.
fn nonzero(len: usize) -> impl Strategy<Value = Terms> {
    prop::collection::btree_map(
        prop::collection::vec(0u32..3, 2),
        prop_oneof![-5i64..=-1, 1i64..=5],
        1..=len,
    )
}

/// The polynomial as Xcas text, in parentheses.
fn text(terms: &Terms) -> String {
    let parts: Vec<String> = terms
        .iter()
        .filter(|&(_, &c)| c != 0)
        .map(|(exps, c)| {
            let mut factors = vec![if *c < 0 {
                format!("({c})")
            } else {
                c.to_string()
            }];
            for (v, &k) in VARIABLES.iter().zip(exps) {
                match k {
                    0 => {}
                    1 => factors.push((*v).to_owned()),
                    _ => factors.push(format!("{v}^{k}")),
                }
            }
            factors.join("*")
        })
        .collect();
    if parts.is_empty() {
        "(0)".to_owned()
    } else {
        format!("({})", parts.join("+"))
    }
}

/// The terms in `x` alone: the other exponents dropped, like terms added,
/// those that cancel left out.
fn in_x(terms: Terms) -> Terms {
    let mut only = Terms::new();
    for (exps, c) in terms {
        let mut e = vec![0; exps.len()];
        e[0] = exps[0];
        *only.entry(e).or_default() += c;
    }
    only.retain(|_, c| *c != 0);
    only
}

/// `name(a*g,b*g)` in this CAS and in giac, the same up to sign, which giac
/// takes from the inputs. `None` for a draw that cancels to 0.
fn gcd_input(name: &str, univariate: bool, a: Terms, b: Terms, g: Terms) -> Option<String> {
    let strip = |t: Terms| if univariate { in_x(t) } else { t };
    let (a, b, g) = (strip(a), strip(b), strip(g));
    // giac returns a nonzero input as it is when the other is 0.
    if [&a, &b, &g].iter().any(|t| t.is_empty()) {
        return None;
    }
    let (a, b, g) = (text(&a), text(&b), text(&g));
    Some(format!("{name}({a}*{g},{b}*{g})"))
}

/// Whether this CAS's answer to `input` is giac's, or giac's negated.
fn up_to_sign(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    let fail = |e: cas_giac::OracleError| TestCaseError::fail(format!("{input}: {e}"));
    let ctx = Context::new();
    let ours = ours(input).map_err(fail)?;
    let theirs = giac.value(input).map_err(fail)?;
    let expand = |e: Expr| {
        ctx.eval(&Expr::call(Symbol::EXPAND, [e]))
            .map_err(|e| TestCaseError::fail(format!("{input}: {e}")))
    };
    let (theirs, negated) = (expand(theirs.clone())?, expand(-theirs)?);
    prop_assert!(
        ours == theirs || ours == negated,
        "{}: giac {}, ours {}",
        input,
        theirs,
        ours
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `heugcd` of two products with a common factor, in one or two
    /// variables: giac's answer up to sign.
    #[test]
    fn random_heugcd(univariate in any::<bool>(), a in nonzero(3), b in nonzero(3), g in nonzero(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = gcd_input("heugcd", univariate, a, b, g) else {
            return Ok(());
        };
        up_to_sign(giac, &input)?;
    }

    /// `modgcd`, as [`random_heugcd`]. Excluded, the departure: one
    /// variable with a gcd whose integer content is not 1, which giac's
    /// modgcd drops (`modgcd(4*x,6*x)` is `x`, `modgcd(6*y,4*y)` is `y`);
    /// those draws are skipped.
    #[test]
    fn random_modgcd(univariate in any::<bool>(), a in nonzero(3), b in nonzero(3), g in nonzero(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = gcd_input("modgcd", univariate, a, b, g) else {
            return Ok(());
        };
        // One variable is x alone or y alone: modgcd(6*y,4*y) is y too (T9).
        if univariate || !input.contains('x') || !input.contains('y') {
            let content = format!("icontent({})", input.replacen("modgcd", "gcd", 1));
            let printed = giac.eval(&content).map_err(|e| TestCaseError::fail(e.to_string()))?;
            if printed.trim() != "1" {
                return Ok(());
            }
        }
        up_to_sign(giac, &input)?;
    }

    /// `psrgcd`, as [`random_heugcd`]. Up to sign covers the departure:
    /// ours has a positive leading coefficient where giac keeps the sign
    /// of its remainders.
    #[test]
    fn random_psrgcd(univariate in any::<bool>(), a in nonzero(3), b in nonzero(3), g in nonzero(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = gcd_input("psrgcd", univariate, a, b, g) else {
            return Ok(());
        };
        up_to_sign(giac, &input)?;
    }

    /// `ezgcd` in two variables, up to sign. Excluded, the departures: one
    /// variable, which giac refuses ("Args must be multivariate
    /// polynomials", also when a draw has no `y`), and the inputs where
    /// giac gives up ("GCD not successful"): draws giac answers with an
    /// error are skipped.
    #[test]
    fn random_ezgcd(a in nonzero(3), b in nonzero(3), g in nonzero(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = gcd_input("ezgcd", false, a, b, g) else {
            return Ok(());
        };
        let printed = giac.eval(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
        if printed.trim_matches('"').starts_with("GIAC_ERROR") {
            return Ok(());
        }
        up_to_sign(giac, &input)?;
    }

    /// `icontent` of a multiple of a product, in two variables: giac's
    /// answer.
    #[test]
    fn random_icontent(k in 1i64..=6, a in nonzero(3), g in nonzero(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let input = format!("icontent({k}*{}*{})", text(&a), text(&g));
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }
}
