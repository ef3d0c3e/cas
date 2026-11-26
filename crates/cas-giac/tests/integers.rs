//! Xcas's integer functions against giac (CAS-548 to CAS-553): `chrem`,
//! `lgcd`, `fracmod` and `iratrecon`, `icomp`, `bernoulli_mod` and
//! `is_prime`, giac's answers read back and compared with ours
//! (`Giac::agrees`: the same expression, or `==`); and the departures, each
//! with giac's answer as it prints it: its errors, where ours leaves the
//! call with a warning, its Horner form, its sign of a polynomial gcd,
//! and the inputs that crash libgiac or whose certificate rests on giac's
//! random ECM, checked on our side only; and random inputs a run to each
//! function but `iratrecon`, `is_prime` without its certificate.

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

/// `chrem` of integers, and of polynomials coefficient by coefficient
/// modulo the product of the moduli (CAS-548).
#[test]
fn chinese_remainders() {
    check(&[
        "chrem([2,3],[5,7])",
        "chrem([2,3],[7,5])",
        "chrem([2,4,6],[3,5,7])",
        "chrem([2,4,6,7],[3,5,7,11])",
        "chrem([2,3,1],[5,7,11])",
        "chrem([1,3],[2,4])",
        "chrem([-2,3],[5,7])",
        "chrem([-1000,3],[5,7])",
        "chrem([3,-1000],[5,7])",
        "chrem([2,3],[2^40,7])",
        "chrem([[1,2],[3,4]])",
        "chrem([2],[5])",
        "chrem([-12],[5])",
        "chrem([x],[5])",
        "chrem([2*x+1,4*x+2,6*x-1,x+1],[3,5,7,11])",
        "chrem([x,1],[3,5])",
        "chrem([3*x+1,1],[3,5])",
        "chrem([7*x,8],[3,5])",
        "chrem([x,0],[3,5])",
        "chrem([x,x],[3,5])",
        "chrem([x*y,1],[3,5])",
        "chrem([x+1,2,3],[3,5,7])",
        "chrem([x+1,x,3],[3,5,7])",
        "chrem([x+1,x+4],[3,6])",
        "chrem([x+1,x+3],[4,6])",
        "chrem([x+5,x+11],[6,15])",
        "chrem([x+7,x+1],[10,4])",
        "chrem([x-1,x+5],[4,10])",
        "chrem([x+1,x+2],[2^40,3])",
        "chrem([x+2^39,x],[2^40,3])",
        "chrem([x-2^39,x],[2^40,3])",
        "chrem([x+2^35,x+2],[2^40,6])",
        "chrem([2147483647*x,1],[3,5])",
        "chrem([x,1],[2147483647,5])",
        "chrem([x+1,2],[2147483647,2147483629])",
    ]);
}

/// `lgcd` of lists and sequences (CAS-549).
#[test]
fn list_gcds() {
    check(&[
        "lgcd([])",
        "lgcd()",
        "lgcd([4])",
        "lgcd([-4])",
        "lgcd([0,0])",
        "lgcd([-4,6])",
        "lgcd(-4,-6)",
        "lgcd(18,12,30)",
        "lgcd([18,12,30])",
        "lgcd([45,75,20,15])",
        "lgcd([1/2,1/3])",
        "lgcd([2*x+2,4*x+4])",
        "lgcd([-2*x-2,4*x+4])",
        "lgcd([6,x])",
        "lgcd([x])",
        "lgcd(x,y)",
        "lgcd([x^2-1,x^2-2*x+1])",
        "lgcd([x^2-2*x+1,x^3-1,x-1])",
        "lgcd(x^2-2*x+1,x^3-1,x-1)",
        "lgcd(6)",
    ]);
}

/// `fracmod` and its alias `iratrecon` (CAS-550).
#[test]
fn fractions_modulo() {
    check(&[
        "fracmod(13,121)",
        "fracmod(41,121)",
        "fracmod(3,7)",
        "fracmod(10,7)",
        "fracmod(-3,7)",
        "fracmod(0,7)",
        "fracmod(5,1)",
        "fracmod([13,41],121)",
        "fracmod([13,[41]],121)",
        "fracmod(123456789,10^20+39)",
        "iratrecon(41,121)",
        "iratrecon(13,121)",
        "fracmod(13)",
        "fracmod(13,121,3)",
    ]);
}

/// `icomp` in NEXCOM's order, with `zeros` (CAS-553).
#[test]
fn compositions() {
    check(&[
        "icomp(4,1)",
        "icomp(3,3)",
        "icomp(4,2)",
        "icomp(4,2,zeros=true)",
        "icomp(5,2,zeros=false)",
        "icomp(3,3,zeros=false)",
        "icomp(6,3,zeros=false)",
        "icomp(3,2,zeros=1)",
    ]);
}

/// `bernoulli_mod`, giac's signed values (CAS-552).
#[test]
fn bernoulli_numbers_modulo() {
    check(&[
        "bernoulli_mod(10,13)",
        "bernoulli_mod(14,13)",
        "bernoulli_mod(1234,1009)",
        "bernoulli_mod(20,10007)",
        "bernoulli_mod(2,7)",
        "bernoulli_mod(8,17)",
        "bernoulli_mod(2,5)",
        "bernoulli_mod(1234,10^6+3)",
        "bernoulli_mod(50,1000003)",
    ]);
}

/// `is_prime`'s 0, 1 and 2 and giac's own certificates, `[[q, a, f], …]`
/// in the order of giac's Pollard rho (CAS-551).
#[test]
fn primes() {
    check(&[
        "is_prime(2147483647)",
        "is_prime(2147483659)",
        "is_prime(4294967311)",
        "is_prime(-2147483659)",
        "is_prime(-9856989898997)",
        "is_prime(0)",
        "is_prime(1)",
        "is_prime(2)",
        "is_prime(-7)",
        "is_prime(1999)",
        "is_prime(9856989898997)",
        "is_prime(9856989898997789789)",
        "is_prime(4294967297)",
        "is_prime(1000000007)",
        "is_prime(1000000007,0)",
        "is_prime(1000000007,2)",
        "is_prime(9856989898997,2)",
        "is_prime(2^64+13)",
        "is_prime(2^89-1)",
        "is_prime(2^127-1)",
        "is_prime(10^30+57)",
        "is_prime(10^30+57,0)",
        "is_prime([7,8,9])",
        "is_prime([2^89-1,4])",
        "is_prime(2,1)",
        "is_prime(3,1)",
        "is_prime(4,1)",
        "is_prime(1,1)",
        "is_prime(9,1)",
        "is_prime(-7,1)",
        "is_prime(1999,1)",
        "is_prime(1000000007,1)",
        "is_prime(9856989898997,1)",
        "is_prime(9856989898997789789,1)",
        "is_prime(2147483659,1)",
        "is_prime(4294967311,1)",
        "is_prime(4294967297,1)",
        "is_prime(18446744073709551629,1)",
        "is_prime(2^31-1,1)",
        "is_prime(2^61-1,1)",
        "is_prime(2^89-1,1)",
        "is_prime(2^127-1,1)",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
const DEPARTURES: [(&str, &str, &str); 30] = [
    // giac's errors; ours leaves the call, with the error as a warning.
    (
        "chrem([2,3],[4,6])",
        "GIAC_ERROR: No Integer Solution Error: Bad Argument Value",
        "chrem([2,3],[4,6])",
    ),
    (
        "chrem([x^2+1,2*x],[3,6])",
        "GIAC_ERROR: No Integer Solution Error: Bad Argument Value",
        "chrem([x^2+1,2*x],[3,6])",
    ),
    (
        "chrem([2,3,4],[5,7])",
        "GIAC_ERROR: chrem([2,3,4],[5,7]) \n Error: Bad Argument Value",
        "chrem([2,3,4],[5,7])",
    ),
    (
        "chrem([],[])",
        "GIAC_ERROR: chrem([],[]) \n Error: Bad Argument Value",
        "chrem([],[])",
    ),
    (
        "chrem(2,5)",
        "GIAC_ERROR: chrem(2,5) \n Error: Bad Argument Value",
        "chrem(2,5)",
    ),
    (
        "chrem([1/2,1],[3,5])",
        "GIAC_ERROR: Non integer argument Error: Bad Argument Type",
        "chrem([1/2,1],[3,5])",
    ),
    (
        "chrem([1,2],[1,5])",
        "GIAC_ERROR: Bad mod value Error: Bad Argument Type",
        "chrem([1,2],[1,5])",
    ),
    (
        "fracmod(2,4)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(2,4)",
    ),
    (
        "fracmod(5,8)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(5,8)",
    ),
    (
        "fracmod(13,x)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(13,x)",
    ),
    (
        "fracmod(1/2,7)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(1/2,7)",
    ),
    (
        "fracmod(1,0)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(1,0)",
    ),
    (
        "icomp(3,4)",
        "GIAC_ERROR: icomp(3,4) \n Error: Bad Argument Value",
        "icomp(3,4)",
    ),
    (
        "icomp(0,1)",
        "GIAC_ERROR: icomp(0,1) \n Error: Bad Argument Value",
        "icomp(0,1)",
    ),
    (
        "icomp(4,0)",
        "GIAC_ERROR: icomp(4,0) \n Error: Bad Argument Value",
        "icomp(4,0)",
    ),
    (
        "icomp(3,2,zeros=x)",
        "GIAC_ERROR: icomp(3,2,'zeros'=x) \n Error: Bad Argument Value",
        "icomp(3,2,zeros=x)",
    ),
    (
        "icomp(5)",
        "GIAC_ERROR: icomp(5) \n Error: Bad Argument Type",
        "icomp(5)",
    ),
    (
        "bernoulli_mod(4,5)",
        "GIAC_ERROR: Not invertible Error: Bad Argument Value",
        "bernoulli_mod(4,5)",
    ),
    (
        "bernoulli_mod(2,3)",
        "GIAC_ERROR: Not invertible Error: Bad Argument Value",
        "bernoulli_mod(2,3)",
    ),
    (
        "bernoulli_mod(3,13)",
        "GIAC_ERROR: bernoulli_mod(3,13) \n Error: Bad Argument Type",
        "bernoulli_mod(3,13)",
    ),
    (
        "bernoulli_mod(10,12)",
        "GIAC_ERROR: bernoulli_mod(10,12) \n Error: Bad Argument Type",
        "bernoulli_mod(10,12)",
    ),
    (
        "bernoulli_mod(0,7)",
        "GIAC_ERROR: bernoulli_mod(0,7) \n Error: Bad Argument Type",
        "bernoulli_mod(0,7)",
    ),
    (
        "bernoulli_mod(4,9)",
        "GIAC_ERROR: bernoulli_mod(4,9) \n Error: Bad Argument Type",
        "bernoulli_mod(4,9)",
    ),
    (
        "is_prime(x)",
        "GIAC_ERROR: is_prime(x) \n Error: Bad Argument Type",
        "is_prime(x)",
    ),
    (
        "is_prime(7/2)",
        "GIAC_ERROR: is_prime(7/2) \n Error: Bad Argument Type",
        "is_prime(7/2)",
    ),
    (
        "is_prime(7,x)",
        "GIAC_ERROR: is_prime(7,x) \n Error: Bad Argument Type",
        "is_prime(7,x)",
    ),
    (
        "is_prime(x,1)",
        "GIAC_ERROR: is_probab_prime_p Error: Bad Argument Type",
        "is_prime(x,1)",
    ),
    (
        "is_prime(7/2,1)",
        "GIAC_ERROR: is_probab_prime_p Error: Bad Argument Type",
        "is_prime(7/2,1)",
    ),
    (
        "is_prime([7,8],1)",
        "GIAC_ERROR: is_probab_prime_p Error: Bad Argument Type",
        "is_prime([7,8],1)",
    ),
    // GCD makes the leading coefficient positive; giac keeps the sign
    // of a single polynomial.
    ("lgcd([-2*x-2])", "-2*x-2", "2*x+2"),
];

/// The departures of [`DEPARTURES`], then giac's Horner form and the
/// inputs that crash libgiac.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, answer, decided) in DEPARTURES {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac prints the polynomial in Horner form; expanded, it is ours.
    let (input, answer) = ("chrem([x^2+1,2*x],[3,5])", "[(-5*x+12)*x-5,15]");
    let printed = giac.eval(input).unwrap();
    if printed != answer {
        wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
    }
    if ours(input) != ours("[-5*x^2+12*x-5,15]")
        || ours(input) != ours(&format!("expand({answer})"))
    {
        wrong.push(format!("{input}: ours is not giac's expanded"));
    }
    // libgiac crashes on polynomial moduli, on p = 2, on a negative
    // certificate above 2^31, and in PARI on primes its n - 1 test cannot
    // prove, so giac is not asked: ours folds through chinrem, leaves
    // bernoulli_mod(2,2), answers 0, and leaves the unproven primes (n - 1 =
    // 192 p q, p and q primes near 2^100).
    for (input, decided) in [
        ("chrem([x+1,x],[x^2+1,x+2])", "chinrem([x+1,x^2+1],[x,x+2])"),
        ("chrem([x+1,x],[x^2+1,x-1])", "chinrem([x+1,x^2+1],[x,x-1])"),
        ("bernoulli_mod(2,2)", "bernoulli_mod(2,2)"),
        ("is_prime(-2147483659,1)", "0"),
        (
            "is_prime(308532104498006741218424255238752343301423685716839787438005313)",
            "is_prime(308532104498006741218424255238752343301423685716839787438005313)",
        ),
        ("is_prime(2^607-1)", "is_prime(2^607-1)"),
        // giac's Pollard rho gives up on 290240017 * 454197539 and its ECM
        // (random curves) finds either factor first: 454197539 locally,
        // 290240017 on CI; ours checked alone.
        (
            "is_prime(10^30+57,1)",
            "[[2,5,1],[3,2,1],[79043,2,1],[3998741,2,1],[290240017,2,1]]",
        ),
    ] {
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answer to `input`, read back, is ours; else a failure naming the
/// input and both answers.
fn agree(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    giac.agrees(input)
        .map_err(|e| TestCaseError::fail(e.to_string()))
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

/// A list of integers in the language.
fn list(values: &[i64]) -> String {
    let items: Vec<String> = values.iter().map(i64::to_string).collect();
    format!("[{}]", items.join(","))
}

/// The odd primes from 5 to 113 and two larger ones.
const PRIMES: [i64; 30] = [
    5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 1009, 10007,
];

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `chrem` of random residues and moduli `2..=30`, coprime or not.
    #[test]
    fn random_chrem(pairs in prop::collection::vec((-30i64..=30, 2i64..=30), 1..=3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        // Residues that disagree modulo a common factor: giac's "No
        // Integer Solution" (departure row chrem([2,3],[4,6])).
        for (i, &(a, m)) in pairs.iter().enumerate() {
            for &(b, n) in &pairs[i + 1..] {
                if (a - b) % gcd(m, n) != 0 {
                    return Ok(());
                }
            }
        }
        let (residues, moduli): (Vec<i64>, Vec<i64>) = pairs.into_iter().unzip();
        agree(giac, &format!("chrem({},{})", list(&residues), list(&moduli)))?;
    }

    /// `lgcd` of a random list of integers.
    #[test]
    fn random_lgcd(values in prop::collection::vec(-60i64..=60, 1..=4)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("lgcd({})", list(&values)))?;
    }

    /// `fracmod` of a random integer modulo a prime (a composite modulus
    /// may reconstruct a denominator it shares: departure rows).
    #[test]
    fn random_fracmod(a in -300i64..=300, p in prop::sample::select(&PRIMES[..])) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("fracmod({a},{p})"))?;
    }

    /// `bernoulli_mod` of an even index modulo a prime from 5 (p = 2
    /// crashes libgiac, p = 3 never inverts).
    #[test]
    fn random_bernoulli_mod(half in 1i64..=150, p in prop::sample::select(&PRIMES[..])) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let n = 2 * half;
        // p − 1 dividing n: giac's "Not invertible" (departure row
        // bernoulli_mod(4,5)).
        if n % (p - 1) == 0 {
            return Ok(());
        }
        agree(giac, &format!("bernoulli_mod({n},{p})"))?;
    }

    /// `icomp(n, k)`, with and without `zeros`, `k ≤ n` (giac errs on
    /// `k > n`: departure row icomp(3,4)).
    #[test]
    fn random_icomp(
        (n, k) in (1u32..=7).prop_flat_map(|n| (Just(n), 1..=n)),
        zeros in prop::sample::select(&["", ",zeros=true", ",zeros=false"][..]),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("icomp({n},{k}{zeros})"))?;
    }

    /// `is_prime` of a random integer below 2^64, small or large, signed;
    /// never the certificate form `is_prime(n,1)`, whose PARI and ECM
    /// paths crash libgiac or draw at random (C0021, departures).
    #[test]
    fn random_is_prime(
        n in prop_oneof![
            (-1000i128..=100_000),
            any::<u64>().prop_map(i128::from),
            any::<u32>().prop_map(i128::from),
        ],
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("is_prime({n})"))?;
    }
}
