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

/// Factorisation and divisors (CAS-268).
#[test]
fn factorisation() {
    check(&[
        (
            "IntegerFactor[360]",
            "Hold[Times[Power[2, 3], Power[3, 2], 5]]",
        ),
        ("IntegerFactor[-12]", "Hold[Times[-1, Power[2, 2], 3]]"),
        ("IntegerFactor[8]", "Hold[Power[2, 3]]"),
        ("IntegerFactor[-8]", "Hold[Times[-1, Power[2, 3]]]"),
        ("IntegerFactor[7]", "7"),
        ("IntegerFactor[-7]", "-7"),
        ("IntegerFactor[0]", "0"),
        ("IntegerFactor[1]", "1"),
        ("IntegerFactor[-1]", "-1"),
        (
            "IntegerFactor[Plus[Power[2, 64], 1]]",
            "Hold[Times[274177, 67280421310721]]",
        ),
        ("IntegerFactorList[360]", "List[2, 3, 3, 2, 5, 1]"),
        ("IntegerFactorList[-12]", "List[-1, 1, 2, 2, 3, 1]"),
        ("IntegerFactorList[1]", "List[]"),
        ("IntegerFactorList[0]", "IntegerFactorList[0]"),
        ("Divisors[12]", "List[1, 2, 3, 4, 6, 12]"),
        ("Divisors[-12]", "List[1, 2, 3, 4, 6, 12]"),
        ("Divisors[1]", "List[1]"),
        ("Divisors[0]", "Divisors[0]"),
        ("EulerPhi[36]", "12"),
        ("EulerPhi[1]", "1"),
        ("EulerPhi[0]", "0"),
        ("EulerPhi[97]", "96"),
        ("EulerPhi[-5]", "EulerPhi[-5]"),
    ]);
    // Large and smooth: trial division and perfect powers at once.
    let ctx = Context::new();
    assert_eq!(
        ctx.eval(&read("EulerPhi[Power[2, 1000]]")).unwrap(),
        ctx.eval(&read("Power[2, 999]")).unwrap()
    );
}
