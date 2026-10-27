use crate::context::Context;
use crate::expr::Expr;

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Primes (CAS-267).
#[test]
fn primes() {
    check(&[
        ("PrimeQ[97]", "True"),
        ("PrimeQ[91]", "False"),
        ("PrimeQ[1]", "False"),
        ("PrimeQ[0]", "False"),
        ("PrimeQ[-7]", "False"),
        ("PrimeQ[2]", "True"),
        ("PrimeQ[Plus[Power[2, 61], -1]]", "True"),
        ("PrimeQ[Plus[Power[2, 64], 13]]", "True"),
        ("PrimeQ[Plus[Power[2, 64], 1]]", "False"),
        ("PrimeQ[2/3]", "False"),
        ("PrimeQ[x]", "PrimeQ[x]"),
        ("PrimeQ[List[2, 3, 4]]", "List[True, True, False]"),
        ("NextPrime[100]", "101"),
        ("NextPrime[-5]", "2"),
        ("NextPrime[1]", "2"),
        ("NextPrime[Power[2, 64]]", "18446744073709551629"),
        ("NextPrime[x]", "NextPrime[x]"),
        ("PreviousPrime[100]", "97"),
        ("PreviousPrime[3]", "2"),
        (
            "PreviousPrime[Plus[Power[2, 64], 13]]",
            "18446744073709551557",
        ),
        ("PreviousPrime[2]", "PreviousPrime[2]"),
        ("PreviousPrime[-4]", "PreviousPrime[-4]"),
    ]);
}
