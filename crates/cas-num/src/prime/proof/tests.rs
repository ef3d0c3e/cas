//! Primality proofs on known primes above 2^64, on composites (strong
//! pseudoprimes among them), and an `n − 1` the limits cannot factor.

use super::Proof;
use crate::Natural;

/// The limits of cas-core's factorisations up to 128 bits.
const RHO: u64 = 1 << 22;
const B1: u64 = 11_000;
const CURVES: u32 = 40;

fn natural(s: &str) -> Natural {
    s.parse().unwrap()
}

fn mersenne(p: u64) -> Natural {
    Natural::from(2u64).pow(p) - 1
}

fn prove(n: &Natural) -> Proof {
    n.prove_prime(RHO, B1, CURVES)
}

#[test]
fn primes_above_2_64() {
    for n in [
        mersenne(89),
        mersenne(107),
        mersenne(127),
        natural("1000000000000000000000000000057"),
        natural("18446744073709551629"),
    ] {
        assert_eq!(prove(&n), Proof::Proven, "{n}");
    }
}

#[test]
fn small_values() {
    for (n, proof) in [
        (0u64, Proof::Composite),
        (1, Proof::Composite),
        (2, Proof::Proven),
        (1999, Proof::Proven),
        (2_147_483_647, Proof::Proven),
        (4_294_967_297, Proof::Composite),
        (3_215_031_751, Proof::Composite),
        (3_825_123_056_546_413_051, Proof::Composite),
        (9_856_989_898_997_789_789, Proof::Proven),
    ] {
        assert_eq!(prove(&Natural::from(n)), proof, "{n}");
    }
}

/// Strong pseudoprimes to the first primes as bases, a Fermat number and a
/// product of two primes near 2^40.
#[test]
fn composites() {
    for n in [
        natural("318665857834031151167461"),
        natural("3317044064679887385961981"),
        Natural::from(2u64).pow(128) + 1,
        natural("1099511627791") * &natural("1099512676421"),
    ] {
        assert_eq!(prove(&n), Proof::Composite, "{n}");
    }
}

/// `n − 1 = 12 p q`, `p` and `q` primes near 2^40: rho finds them with the
/// default effort, not with one iteration and no curves.
#[test]
fn unknown_under_tight_limits() {
    let n = natural("14507123671541913288192133");
    assert_eq!(n.prove_prime(1, 0, 0), Proof::Unknown);
    assert_eq!(prove(&n), Proof::Proven);
}
