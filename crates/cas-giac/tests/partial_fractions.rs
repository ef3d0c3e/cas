//! Partial fractions and fraction forms against giac (CAS-529, CAS-540,
//! CAS-544): `cpartfrac`, `convert` with `partfrac`, `parfrac`,
//! `fullparfrac`, `confrac` and `polynom`, `propfrac`, `f2nd`, `fxnd` and
//! `ratnormal`, giac's answers read back and compared after `Expand` and
//! `Together` of the difference, `i` and square roots included. The
//! departures, each with giac's answer and ours as decided. And 128 random
//! fractions a run to `cpartfrac`, `propfrac` and `f2nd`.

use cas_core::{Context, Expr, Symbol};
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

/// Whether `a − b` is 0 once expanded and put together, lists element by
/// element.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| same(ctx, x, y));
    }
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// giac 2.0's answers, read back, are this CAS's.
#[test]
fn agreement() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "cpartfrac(1/(x^2+1))",
        "cpartfrac(1/(x^2+2))",
        "cpartfrac(x/(4-x^2))",
        "cpartfrac(1/(x^2+x+1))",
        "cpartfrac(1/((x^2+1)^2))",
        "cpartfrac(1/(x^3-1))",
        "cpartfrac(x^3/(x^2+1))",
        "cpartfrac(1/(x^2-2))",
        "cpartfrac(1/(y^2+1),y)",
        "convert(1/(x^2+1),partfrac)",
        "convert(1/(x^2-1),partfrac)",
        "convert(1/(x^2-1),parfrac)",
        "convert(1/(x^2-1),fullparfrac)",
        "convert(x^3/(x^2-2),partfrac)",
        "convert(7/2,confrac)",
        "convert(1/3,confrac)",
        "convert(series(sin(x),x=0,5),polynom)",
        "convert(series(exp(x),x=0,3),polynom)",
        "convert(x^2+1,polynom)",
        "propfrac(x^3/(x^2+1))",
        "propfrac((x^2+2*x+3)/(x+1))",
        "propfrac(x/(x^2+1))",
        "propfrac(x^2+1)",
        "propfrac((x^3+y)/(x+y))",
        "propFrac(x^3/(x^2+1))",
        "propfrac(28/12)",
        "propfrac(y^2/(y+1),y)",
        "f2nd(x/(x^2+x))",
        "f2nd((x^2-1)/(x-1))",
        "f2nd(3/6)",
        "f2nd(x+1/x)",
        "fxnd(x/(x^2+x))",
        "fxnd(7)",
        "ratnormal((x^2-1)/(x-1))",
        "ratnormal(1/x+1/y)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: the input, giac's answer
/// as it prints it, ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac splits x^2±i in sqrt(2)(1±i)/2; ours keeps the quadratics
        // over ℚ(i) whole.
        (
            "cpartfrac(1/(x^4+1))",
            "1/((2+2*i)*sqrt(2))/(x+(1-i)*sqrt(2)/2)-1/((2+2*i)*sqrt(2))/(x+(-1+i)*sqrt(2)/2)-1/((-2+2*i)*sqrt(2))/(x+(1+i)*sqrt(2)/2)+1/((-2+2*i)*sqrt(2))/(x+(-1-i)*sqrt(2)/2)",
            "-i/2/(x^2-i)+i/2/(x^2+i)",
        ),
        // giac's convert drops the last quotient; ours is dfc's.
        ("convert(355/113,confrac)", "[3,7,15,1]", "[3,7,16]"),
        // giac loses the sign; ours is dfc's, the floor first.
        ("convert(-7/3,confrac)", "[2,3]", "[-3,1,2]"),
        // giac's 15 quotients from a float; ours dfc's exact period.
        (
            "convert(sqrt(2),confrac)",
            "[1,2,2,2,2,2,2,2,2,2,2,2,2,2,2]",
            "[1,[2]]",
        ),
        // giac's quotients of a float; ours leaves dfc(pi), which needs a
        // count.
        ("convert(pi,confrac)", "[3,7,15,1,292,1,1,1,2]", "dfc(pi)"),
    ];
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        match (ours(input), ours(decided)) {
            (Ok(found), Ok(decided)) if found == decided || same(&ctx, &found, &decided) => {}
            (found, _) => wrong.push(format!("{input}: ours {found:?}, not {decided}")),
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac's third argument names a variable it stores the list in; ours
    // takes and ignores it.
    let input = "convert(sqrt(2),confrac,'L'),L";
    let printed = giac.eval(input).unwrap();
    if printed != "[1,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[1,2,2,2,2,2,2,2,2,2,2,2,2,2,2]" {
        wrong.push(format!("{input}: giac prints {printed}"));
    }
    if ours("convert(sqrt(2),confrac,L)").ok() != ours("[1,[2]]").ok() {
        wrong.push("convert(sqrt(2),confrac,L): ours is not [1,[2]]".to_owned());
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A numerator of degree below 5, coefficients in `-5..=5` from `x^0` up.
fn numerator() -> impl Strategy<Value = Vec<i64>> {
    prop::collection::vec(-5i64..=5, 1..=5)
}

/// A denominator: one to three linear factors `a*x-b`, `a` in `1..=2`, `b`
/// in `-3..=3`, and maybe a quadratic `x^2+c`, `c` in `1..=5`, roots
/// `±i*sqrt(c)`. Never a quartic like `x^4+1`: there giac splits `x^2±i`
/// and ours keeps it whole (README departure).
fn denominator() -> impl Strategy<Value = (Vec<(i64, i64)>, Option<i64>)> {
    (
        prop::collection::vec((1i64..=2, -3i64..=3), 1..=3),
        prop::option::of(1i64..=5),
    )
}

/// `name(p/q)` as Xcas text; `None` for a zero numerator.
fn fraction(name: &str, p: &[i64], q: &(Vec<(i64, i64)>, Option<i64>)) -> Option<String> {
    let terms: Vec<String> = p
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c != 0)
        .map(|(k, c)| format!("({c})*x^{k}"))
        .collect();
    if terms.is_empty() {
        return None;
    }
    let mut factors: Vec<String> = q.0.iter().map(|(a, b)| format!("({a}*x-({b}))")).collect();
    if let Some(c) = q.1 {
        factors.push(format!("(x^2+{c})"));
    }
    Some(format!(
        "{name}(({})/({}))",
        terms.join("+"),
        factors.join("*")
    ))
}

/// Whether giac's answer to `input` and ours differ by 0 ([`same`]).
fn check(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(ours)) if same(&Context::new(), &ours, &theirs) => Ok(()),
        (theirs, ours) => Err(TestCaseError::fail(format!(
            "{input}: ours {ours:?}, giac's {theirs:?}"
        ))),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `cpartfrac` of `p/q`, [`denominator`]'s `q`: real linear factors and
    /// maybe `x^2+c`, split over `ℂ`. Excluded: quartics like `x^4+1`
    /// (README departure, giac splits `x^2±i` in square roots of `i`).
    #[test]
    fn random_cpartfrac(p in numerator(), q in denominator()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = fraction("cpartfrac", &p, &q) else {
            return Ok(());
        };
        check(giac, &input)?;
    }

    /// `propfrac` of `p/q`, as [`random_cpartfrac`].
    #[test]
    fn random_propfrac(p in numerator(), q in denominator()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = fraction("propfrac", &p, &q) else {
            return Ok(());
        };
        check(giac, &input)?;
    }

    /// `f2nd` of `p/q`, as [`random_cpartfrac`]: numerator and denominator,
    /// each compared.
    #[test]
    fn random_f2nd(p in numerator(), q in denominator()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let Some(input) = fraction("f2nd", &p, &q) else {
            return Ok(());
        };
        check(giac, &input)?;
    }
}
