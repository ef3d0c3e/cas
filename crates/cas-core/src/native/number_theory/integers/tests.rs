//! Tests of Xcas's integer functions (CAS-548 to CAS-553): giac's answers
//! from the phase 1.3 probes, and the calls its errors leave.

use super::{bernoulli_mod, certificate_factors, compositions, ichinrem, invmod, pollard};
use crate::{Context, Expr};
use cas_num::{Integer, Natural};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input in full form stays as it is written (giac's errors).
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let found = ctx.eval(&read(input)).unwrap();
        if found != read(input) {
            wrong.push(format!("{input}: {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn chinese_remainders_of_integers() {
    check(&[
        (
            "ChineseRemainderLists[List[2, 3], List[5, 7]]",
            "List[17, 35]",
        ),
        (
            "ChineseRemainderLists[List[2, 4, 6], List[3, 5, 7]]",
            "List[104, 105]",
        ),
        (
            "ChineseRemainderLists[List[1, 3], List[2, 4]]",
            "List[3, 4]",
        ),
        (
            "ChineseRemainderLists[List[-2, 3], List[5, 7]]",
            "List[3, 35]",
        ),
        (
            "ChineseRemainderLists[List[-1000, 3], List[5, 7]]",
            "List[10, 35]",
        ),
        (
            "ChineseRemainderLists[List[List[1, 2], List[3, 4]]]",
            "List[10, 12]",
        ),
        // One pair is itself.
        ("ChineseRemainderLists[List[-12], List[5]]", "List[-12, 5]"),
        ("ChineseRemainderLists[List[2], List[1]]", "List[2, 1]"),
    ]);
    stays(&[
        "ChineseRemainderLists[List[2, 3], List[4, 6]]",
        "ChineseRemainderLists[List[2, 3, 4], List[5, 7]]",
        "ChineseRemainderLists[List[], List[]]",
        "ChineseRemainderLists[2, 5]",
        "ChineseRemainderLists[List[1, 2], List[1, 5]]",
        "ChineseRemainderLists[List[1/2, 1], List[3, 5]]",
    ]);
}

#[test]
fn chinese_remainders_of_polynomials() {
    check(&[
        (
            "ChineseRemainderLists[List[Plus[Power[x, 2], 1], Times[2, x]], List[3, 5]]",
            "List[Plus[Times[-5, Power[x, 2]], Times[12, x], -5], 15]",
        ),
        (
            "ChineseRemainderLists[List[Plus[Times[2, x], 1], Plus[Times[4, x], 2], \
             Plus[Times[6, x], -1], Plus[x, 1]], List[3, 5, 7, 11]]",
            "List[Plus[Times[419, x], 727], 1155]",
        ),
        (
            "ChineseRemainderLists[List[x, 1], List[3, 5]]",
            "List[Plus[Times[-5, x], 6], 15]",
        ),
        // The leading coefficient vanishes, and the integer is reduced.
        (
            "ChineseRemainderLists[List[Plus[Times[3, x], 1], 1], List[3, 5]]",
            "List[1, 15]",
        ),
        (
            "ChineseRemainderLists[List[Times[x, y], 1], List[3, 5]]",
            "List[Plus[Times[-5, x, y], 6], 15]",
        ),
        // Moduli not coprime: modulo the product.
        (
            "ChineseRemainderLists[List[Plus[x, 1], Plus[x, 4]], List[3, 6]]",
            "List[Plus[x, 4], 18]",
        ),
        (
            "ChineseRemainderLists[List[Plus[x, 5], Plus[x, 11]], List[6, 15]]",
            "List[Plus[x, -19], 90]",
        ),
        (
            "ChineseRemainderLists[List[Plus[x, 1], 2, 3], List[3, 5, 7]]",
            "List[Plus[Times[70, x], -53], 105]",
        ),
        (
            "ChineseRemainderLists[List[Times[2147483647, x], 1], List[3, 5]]",
            "List[Plus[Times[2147483635, x], 6], 15]",
        ),
        (
            "ChineseRemainderLists[List[Plus[x, 1], 2], List[2147483647, 2147483629]]",
            "List[Plus[Times[256204776296123963, x], -256204776296123961], 4611685975477714963]",
        ),
        (
            "ChineseRemainderLists[List[Plus[x, Power[2, 39]], x], List[Power[2, 40], 3]]",
            "List[Plus[x, 1649267441664], 3298534883328]",
        ),
    ]);
    stays(&["ChineseRemainderLists[List[Plus[1, Power[x, 2]], Times[2, x]], List[3, 6]]"]);
}

#[test]
fn coefficient_step() {
    let i = Integer::from;
    assert_eq!(ichinrem(&i(1), &i(0), &i(3), &i(5)), Some(i(-5)));
    assert_eq!(ichinrem(&i(0), &i(2), &i(3), &i(5)), Some(i(12)));
    assert_eq!(ichinrem(&i(-1), &i(5), &i(4), &i(10)), Some(i(15)));
    assert_eq!(ichinrem(&i(1), &i(0), &i(3), &i(6)), None);
    assert_eq!(invmod(3, 5), Some(2));
    assert_eq!(invmod(0, 5), None);
}

#[test]
fn list_gcds() {
    check(&[
        ("ListGCD[List[]]", "1"),
        ("ListGCD[]", "1"),
        ("ListGCD[List[4]]", "4"),
        ("ListGCD[List[-4]]", "4"),
        ("ListGCD[List[0, 0]]", "0"),
        ("ListGCD[List[-4, 6]]", "2"),
        ("ListGCD[-4, -6]", "2"),
        ("ListGCD[List[45, 75, 20, 15]]", "5"),
        ("ListGCD[List[1/2, 1/3]]", "1/6"),
        (
            "ListGCD[List[Plus[Times[2, x], 2], Plus[Times[4, x], 4]]]",
            "Plus[Times[2, x], 2]",
        ),
        ("ListGCD[List[6, x]]", "1"),
        ("ListGCD[List[x]]", "x"),
    ]);
    stays(&["ListGCD[6]"]);
}

#[test]
fn fractions_modulo() {
    check(&[
        ("FractionMod[13, 121]", "-4/9"),
        ("FractionMod[41, 121]", "2/3"),
        ("FractionMod[10, 7]", "-1/2"),
        ("FractionMod[-3, 7]", "1/2"),
        ("FractionMod[0, 7]", "0"),
        ("FractionMod[0, x]", "0"),
        (
            "FractionMod[List[13, List[41]], 121]",
            "List[-4/9, List[2/3]]",
        ),
        (
            "FractionMod[123456789, Plus[Power[10, 20], 39]]",
            "123456789",
        ),
    ]);
    stays(&[
        "FractionMod[2, 4]",
        "FractionMod[5, 8]",
        "FractionMod[13, x]",
        "FractionMod[1/2, 7]",
        "FractionMod[1, 0]",
        "FractionMod[1, -7]",
        "FractionMod[13]",
    ]);
}

#[test]
fn compositions_in_nexcom_order() {
    assert_eq!(compositions(4, 1), [[4]]);
    assert_eq!(
        compositions(3, 3),
        [
            [3, 0, 0],
            [2, 1, 0],
            [1, 2, 0],
            [0, 3, 0],
            [2, 0, 1],
            [1, 1, 1],
            [0, 2, 1],
            [1, 0, 2],
            [0, 1, 2],
            [0, 0, 3]
        ]
    );
    check(&[
        (
            "IntegerCompositions[4, 2]",
            "List[List[4, 0], List[3, 1], List[2, 2], List[1, 3], List[0, 4]]",
        ),
        (
            "IntegerCompositions[5, 2, Equation[Solve, False]]",
            "List[List[4, 1], List[3, 2], List[2, 3], List[1, 4]]",
        ),
        (
            "IntegerCompositions[3, 3, Equal[Solve, False]]",
            "List[List[1, 1, 1]]",
        ),
        (
            "IntegerCompositions[3, 2, Equation[Solve, 1]]",
            "List[List[3, 0], List[2, 1], List[1, 2], List[0, 3]]",
        ),
    ]);
    stays(&[
        "IntegerCompositions[3, 4]",
        "IntegerCompositions[0, 1]",
        "IntegerCompositions[4, 0]",
        "IntegerCompositions[5]",
        "IntegerCompositions[3, 2, Equation[Solve, x]]",
        // C(10⁴+9, 9) entries pass the limit on terms.
        "IntegerCompositions[10000, 10]",
    ]);
}

#[test]
fn bernoulli_numbers_modulo() {
    // giac's values, both summations and the reduction of k.
    for (k, p, b) in [
        (10, 13, -8),
        (14, 13, 12),
        (1234, 1009, -483),
        (20, 10007, 3504),
        (2, 7, -1),
        (8, 17, -4),
        (2, 5, 1),
        (1234, 1_000_003, 351_633),
        (50, 1_000_003, -739_424),
    ] {
        assert_eq!(bernoulli_mod(k, p), Some(b), "bernoulli_mod({k}, {p})");
    }
    check(&[("BernoulliMod[10, 13]", "-8")]);
    stays(&[
        "BernoulliMod[4, 5]",
        "BernoulliMod[2, 3]",
        "BernoulliMod[2, 2]",
        "BernoulliMod[3, 13]",
        "BernoulliMod[10, 12]",
        "BernoulliMod[0, 7]",
        "BernoulliMod[4, -5]",
        "BernoulliMod[4]",
    ]);
}

/// giac's 0, 1, 2 and certificates (CAS-551); a prime the limits do not
/// prove stays.
#[test]
fn primes() {
    check(&[
        ("IsPrime[1999]", "2"),
        ("IsPrime[2147483647]", "2"),
        ("IsPrime[2147483659]", "1"),
        ("IsPrime[Plus[Power[2, 89], -1]]", "1"),
        ("IsPrime[4294967297]", "0"),
        ("IsPrime[-7]", "0"),
        ("IsPrime[1]", "0"),
        ("IsPrime[List[7, 8, 9]]", "List[2, 0, 0]"),
        ("IsPrime[1000000007, 2]", "2"),
        ("IsPrime[2, 1]", "List[]"),
        ("IsPrime[4, 1]", "0"),
        ("IsPrime[-7, 1]", "0"),
        (
            "IsPrime[1000000007, 1]",
            "List[List[2, 5, 1], List[500000003, 2, 1]]",
        ),
        (
            "IsPrime[9856989898997, 1]",
            "List[List[2, 2, 1], List[101, 2, 1], List[24398489849, 2, 0]]",
        ),
        (
            "IsPrime[Plus[Power[2, 127], -1], 1]",
            "List[List[2, 3, 1], List[3, 5, 1], List[7, 3, 1], List[19, 3, 1], List[43, 3, 1], \
             List[73, 3, 1], List[127, 2, 1], List[337, 3, 1], List[5419, 3, 1], List[92737, 3, 1]]",
        ),
    ]);
    stays(&[
        "IsPrime[x]",
        "IsPrime[x, 1]",
        "IsPrime[7, x]",
        "IsPrime[List[7, 8], 1]",
        "IsPrime[7, 1, 2]",
        // 2^607 - 1: n - 1 is not factored far enough within the limits.
        "IsPrime[531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127]",
    ]);
}

/// giac's rho: the factor it finds first decides the certificate's order;
/// it gives up on `290240017 · 454197539` within 15000 iterations.
#[test]
fn giac_pollard() {
    let n = Natural::from(131_826_301_440_718_163u64);
    assert_eq!(pollard(&n, 1, 1 << 22), None);
    let n = Natural::from(101u64 * 103);
    let b = pollard(&n, 1, 1 << 22).unwrap();
    assert!(b == Natural::from(101u64) || b == Natural::from(103u64));
    let o = Natural::from(1_000_000_006u64);
    let found = certificate_factors(&o).unwrap();
    assert_eq!(found, [Natural::from(2u64), Natural::from(500_000_003u64)]);
}
