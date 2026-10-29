//! Modular arithmetic against giac (CAS-DOC-43, CAS-PLAN-37): giac's `%`
//! objects, read back, are this CAS's `Modular` values (its distributed
//! coefficients taken in by one value): residues, inverses, sums, products
//! and powers, expressions reduced, derivatives and substitutions, and
//! over a prime field gcds and divisions. giac's `factors` lists are this
//! CAS's `Factor` factor by factor, each read as a modular value, in any
//! order. A result 0 giac writes without its modulus.

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

fn modular(a: Expr, p: i64) -> Expr {
    Expr::call(Symbol::MODULAR, [a, Expr::from(p)])
}

/// Residues and their arithmetic (CAS-450).
#[test]
fn residues() {
    check(&[
        "3 % 7",
        "5 % 7",
        "-3 % 7",
        "1/2 % 7",
        "-7 % 7",
        "2 % 6",
        "(10^20+3) % 97",
        "(3 % 7)^-1",
        "(3 % 7)*5",
        "(3 % 7)+(6 % 7)",
        "(7 % 10)*(9 % 10)",
        "(6 % 9)*(3 % 9)",
        "(3 % 101)^1007",
        "(5 % 12)^-1",
        "(2 % 1000003)^-5",
    ]);
}

/// Expressions reduced, derivatives and substitutions (CAS-450, CAS-452).
#[test]
fn expressions() {
    check(&[
        "((x+1)^5) % 5",
        "((x+y)^3) % 3",
        "(x^2+3*x+2) % 2",
        "(7*x^2+3*x+10) % 5",
        "(x^2+1/2*x) % 7",
        "((x+1) % 3)^2",
        "(x % 7)*(x+3)",
        "expand(((x+2)^4) % 3)",
        "diff((x^3+2*x) % 5,x)",
        "subst((x^2+1) % 7,x=3)",
    ]);
}

/// Over a prime field: gcds, quotients and remainders (CAS-453).
#[test]
fn gcds_and_divisions() {
    check(&[
        "gcd((x^2-1) % 5,(x^2+3*x+2) % 5)",
        "gcd((x^4-1) % 7,(x^6-1) % 7)",
        "gcd((x^3+x+1) % 2,(x^2+x+1) % 2)",
        "quorem((x^3+1) % 5,(x^2+1) % 5)",
        "quo((x^5+3*x+1) % 11,(2*x^2+1) % 11)",
        "rem((x^5+3*x+1) % 11,(2*x^2+1) % 11)",
    ]);
}

/// Results 0, which giac writes without the modulus: this CAS's
/// `Modular[0, p]`.
#[test]
fn zeros() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, p) in [
        ("(3 % 7)+1/2", 7),
        ("rem((x^101-x) % 101,(x^2+1) % 101)", 101),
    ] {
        let theirs = giac.value(input);
        let ours = ours(input);
        match (&theirs, &ours) {
            (Ok(t), Ok(o)) if *t == Expr::from(0) && *o == modular(Expr::from(0), p) => {}
            _ => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The factors with their multiplicities, each as a modular value,
/// sorted; the constant 1 left out.
fn pairs(ctx: &Context, items: Vec<(Expr, Expr)>, p: i64) -> Vec<(Expr, Expr)> {
    let one = modular(Expr::from(1), p);
    let mut pairs: Vec<(Expr, Expr)> = items
        .into_iter()
        .map(|(base, k)| (ctx.eval(&modular(base, p)).unwrap(), k))
        .filter(|(base, _)| *base != one)
        .collect();
    pairs.sort_by_key(|(base, k)| format!("{base:?} {k:?}"));
    pairs
}

/// `Factor`'s product as factors and multiplicities.
fn items(product: &Expr) -> Vec<(Expr, Expr)> {
    let factors = if product.head() == Some(&Expr::from(Symbol::TIMES)) {
        product.args().to_vec()
    } else {
        vec![product.clone()]
    };
    factors
        .into_iter()
        .map(|f| match (f.head(), f.args()) {
            (Some(h), [base, k]) if *h == Expr::from(Symbol::POWER) => (base.clone(), k.clone()),
            _ => (f, Expr::from(1)),
        })
        .collect()
}

/// Factorizations over `𝔽ₚ` (CAS-453): giac's `factors`.
#[test]
fn factorizations() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (f, p) in [
        ("x^4+1", 3),
        ("x^8-x", 2),
        ("x^9-x", 3),
        ("x^16+x", 2),
        ("x^25-x", 5),
        ("x^6+2*x^3+x+5", 7),
        ("2*x^2+2", 3),
        ("3*x^6+x^4+4*x+1", 7),
        ("x^5+x+1", 101),
        ("(x^2+1)^3*(x+2)", 5),
        ("x^3+x^2+x+1", 2),
        ("x^4+x+1", 2),
        ("x^10+x^3+1", 2),
        ("x^7+3*x+2", 11),
        ("x^12-1", 13),
        ("x^4+2", 5),
    ] {
        let theirs = giac.value(&format!("factors(({f}) % {p})")).map(|list| {
            let items = list
                .args()
                .chunks(2)
                .map(|c| (c[0].clone(), c[1].clone()))
                .collect();
            pairs(&ctx, items, p)
        });
        let ours = ours(&format!("factor(({f}) % {p})")).map(|e| pairs(&ctx, items(&e), p));
        match (&theirs, &ours) {
            (Ok(t), Ok(o)) if t == o => {}
            _ => wrong.push(format!("{f} mod {p}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
