//! Miller–Rabin, the strong Lucas test, Baillie–PSW and `next_prime`.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, d, s, q, u, v, k)"
)]

use super::Primality;
use super::small::{SIEVE_LIMIT, TRIAL_LIMIT, groups, is_small_prime, primes};
use crate::modular::{ModDouble, ModRing, ModWord, Residue};
use crate::{Integer, Natural};

#[cfg(test)]
mod tests;

/// Miller–Rabin bases that decide every `n < 2^64` (Jim Sinclair, 2011), a
/// base divisible by `n` being skipped.
const U64_BASES: [u64; 7] = [2, 325, 9375, 28178, 450_775, 9_780_504, 1_795_265_022];

impl Natural {
    /// Whether `self` is prime, labelled by how certain the answer is.
    ///
    /// Values below 2^16 are looked up. Larger ones are trial-divided by the
    /// primes below 2^12; below 2^64 the Miller–Rabin test with a fixed set
    /// of seven bases then decides, so the answer is [`Primality::Prime`] or
    /// [`Primality::NotPrime`]. Above, the Baillie–PSW test (a base-2
    /// Miller–Rabin test and a strong Lucas test with Selfridge's
    /// parameters) answers [`Primality::ProbablePrime`] or
    /// [`Primality::NotPrime`].
    #[must_use]
    pub fn primality(&self) -> Primality {
        let proven = |p| {
            if p {
                Primality::Prime
            } else {
                Primality::NotPrime
            }
        };
        match *self.limbs() {
            [] => Primality::NotPrime,
            [n] => proven(is_prime_u64(n)),
            _ => {
                if has_small_factor(self) {
                    Primality::NotPrime
                } else if baillie_psw(self) {
                    Primality::ProbablePrime
                } else {
                    Primality::NotPrime
                }
            }
        }
    }

    /// The strong probable-prime (Miller–Rabin) test of `self` to `base`:
    /// with `self - 1 = d 2^s`, `d` odd, whether `base^d ≡ 1` or
    /// `base^(d 2^r) ≡ -1 (mod self)` for some `r < s`. Every odd prime
    /// passes for every base; a base divisible by `self` passes trivially.
    /// Even values pass only if they are 2, and 0 and 1 never pass.
    #[must_use]
    pub fn is_strong_probable_prime(&self, base: &Natural) -> bool {
        if !self.bit(0) || self.is_one() {
            return self.limbs() == [2];
        }
        let ring = ModRing::new(self.clone());
        let a = ring.elem(base);
        ring.is_zero(&a) || strong_probable_prime(&ring, self, &a)
    }

    /// The smallest prime greater than `self`. Above 2^64 it is a probable
    /// prime in the sense of [`Natural::primality`].
    #[must_use]
    pub fn next_prime(&self) -> Natural {
        if let Some(p) = u64::try_from(self).ok().and_then(next_prime_u64) {
            return Natural::from(p);
        }
        // Past the last prime below 2^64, so every candidate is above it.
        let start = self + 1u64;
        let mut c = if start.bit(0) { start } else { start + 1u64 };
        while !c.primality().is_probable_prime() {
            c += 2u64;
        }
        c
    }
}

/// The smallest prime above `n`, if it fits a word.
fn next_prime_u64(n: u64) -> Option<u64> {
    if n < 2 {
        return Some(2);
    }
    // The smallest odd value above n.
    let mut c = n.checked_add(1)? | 1;
    while !is_prime_u64(c) {
        c = c.checked_add(2)?;
    }
    Some(c)
}

/// Deterministic primality of a word.
pub(crate) fn is_prime_u64(n: u64) -> bool {
    if n < u64::from(SIEVE_LIMIT) {
        return is_small_prime(u32::try_from(n).expect("below 2^16"));
    }
    // A few divisions reject most composites before the powers.
    if primes()[..16].iter().any(|&p| n % u64::from(p) == 0) {
        return false;
    }
    let w = ModWord::new(n);
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    U64_BASES.iter().all(|&b| {
        let a = w.reduce(b);
        if a == 0 {
            return true;
        }
        let mut y = w.pow(a, d);
        if y == 1 || y == n - 1 {
            return true;
        }
        for _ in 1..s {
            y = w.mul(y, y);
            if y == n - 1 {
                return true;
            }
        }
        false
    })
}

/// Whether a prime below [`TRIAL_LIMIT`] divides `n` (which exceeds them).
fn has_small_factor(n: &Natural) -> bool {
    let primes = primes();
    let trial = primes.partition_point(|&p| p < TRIAL_LIMIT);
    groups()
        .iter()
        .take_while(|(_, range)| range.start < trial)
        .any(|(product, range)| {
            let r = n.rem_limb(*product);
            primes[range.clone()].iter().any(|&p| r % u64::from(p) == 0)
        })
}

/// The ring operations the probable-prime tests use, so that a two-limb
/// modulus can run them on `u128` ([`ModDouble`]) and any other on
/// [`ModRing`].
trait TestRing {
    type Elem: Clone + PartialEq;
    fn small(&self, a: i64) -> Self::Elem;
    fn one(&self) -> Self::Elem;
    fn is_zero(&self, x: &Self::Elem) -> bool;
    fn add(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn sub(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn neg(&self, x: &Self::Elem) -> Self::Elem;
    fn half(&self, x: &Self::Elem) -> Self::Elem;
    fn mul(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn square(&self, x: &Self::Elem) -> Self::Elem;
    fn pow(&self, x: &Self::Elem, e: &Natural) -> Self::Elem;
}

impl TestRing for ModRing {
    type Elem = Residue;
    fn small(&self, a: i64) -> Residue {
        self.elem_integer(&Integer::from(a))
    }
    fn one(&self) -> Residue {
        ModRing::one(self)
    }
    fn is_zero(&self, x: &Residue) -> bool {
        ModRing::is_zero(self, x)
    }
    fn add(&self, x: &Residue, y: &Residue) -> Residue {
        ModRing::add(self, x, y)
    }
    fn sub(&self, x: &Residue, y: &Residue) -> Residue {
        ModRing::sub(self, x, y)
    }
    fn neg(&self, x: &Residue) -> Residue {
        ModRing::neg(self, x)
    }
    fn half(&self, x: &Residue) -> Residue {
        ModRing::half(self, x)
    }
    fn mul(&self, x: &Residue, y: &Residue) -> Residue {
        ModRing::mul(self, x, y)
    }
    fn square(&self, x: &Residue) -> Residue {
        ModRing::square(self, x)
    }
    fn pow(&self, x: &Residue, e: &Natural) -> Residue {
        ModRing::pow(self, x, e)
    }
}

impl TestRing for ModDouble {
    type Elem = u128;
    fn small(&self, a: i64) -> u128 {
        // |a| is far below n > 2^64.
        let m = self.to_mont(u128::from(a.unsigned_abs()));
        if a < 0 { ModDouble::neg(*self, m) } else { m }
    }
    fn one(&self) -> u128 {
        ModDouble::one(*self)
    }
    fn is_zero(&self, x: &u128) -> bool {
        *x == 0
    }
    fn add(&self, x: &u128, y: &u128) -> u128 {
        ModDouble::add(*self, *x, *y)
    }
    fn sub(&self, x: &u128, y: &u128) -> u128 {
        ModDouble::sub(*self, *x, *y)
    }
    fn neg(&self, x: &u128) -> u128 {
        ModDouble::neg(*self, *x)
    }
    fn half(&self, x: &u128) -> u128 {
        ModDouble::half(*self, *x)
    }
    fn mul(&self, x: &u128, y: &u128) -> u128 {
        ModDouble::mul(*self, *x, *y)
    }
    fn square(&self, x: &u128) -> u128 {
        ModDouble::mul(*self, *x, *x)
    }
    fn pow(&self, x: &u128, e: &Natural) -> u128 {
        self.pow_mont(*x, e)
    }
}

/// Baillie–PSW for odd `n` above 2^64 without small factors, on `u128` for
/// two limbs.
fn baillie_psw(n: &Natural) -> bool {
    if let [lo, hi] = *n.limbs() {
        let ring = ModDouble::new(u128::from(lo) | u128::from(hi) << 64);
        return baillie_psw_in(&ring, n);
    }
    baillie_psw_in(&ModRing::new(n.clone()), n)
}

fn baillie_psw_in<R: TestRing>(ring: &R, n: &Natural) -> bool {
    let two = ring.small(2);
    strong_probable_prime(ring, n, &two) && !n.is_perfect_square() && strong_lucas(ring, n)
}

/// Miller–Rabin to base `a` (an element of `ring`, modulo odd `n > 1`).
fn strong_probable_prime<R: TestRing>(ring: &R, n: &Natural, a: &R::Elem) -> bool {
    let n1 = n - 1u64;
    let s = n1.trailing_zeros().expect("n > 1");
    let minus_one = ring.neg(&ring.one());
    let mut y = ring.pow(a, &(&n1 >> s));
    if y == ring.one() || y == minus_one {
        return true;
    }
    for _ in 1..s {
        y = ring.square(&y);
        if y == minus_one {
            return true;
        }
        if y == ring.one() {
            return false;
        }
    }
    false
}

/// The strong Lucas probable-prime test with Selfridge's parameters
/// (method A): `D` is the first of 5, -7, 9, -11, … with Jacobi symbol
/// `(D/n) = -1`, `P = 1`, `Q = (1 - D)/4`. With `n + 1 = d 2^s`, `d` odd,
/// `n` passes if `U_d ≡ 0` or `V_(d 2^r) ≡ 0 (mod n)` for some `r < s`.
///
/// `n` is odd, above the trial-division bound and not a square (else no
/// such `D` exists).
fn strong_lucas<R: TestRing>(ring: &R, n: &Natural) -> bool {
    let n_int = Integer::from(n);
    let mut d = 5i64;
    loop {
        match Integer::from(d).jacobi(&n_int) {
            -1 => break,
            // n > |D|, so a common factor proves n composite.
            0 => return false,
            _ => d = if d > 0 { -d - 2 } else { -d + 2 },
        }
    }
    let q = (1 - d) / 4;
    let dd = ring.small(d);
    let qq = ring.small(q);
    let n1 = n + 1u64;
    let s = n1.trailing_zeros().expect("n + 1 > 0");
    let k = &n1 >> s;
    // (U_j, V_j, Q^j) from j = 1, doubling and incrementing along the bits
    // of k: U_2j = U_j V_j, V_2j = V_j^2 - 2 Q^j, and with P = 1,
    // U_(j+1) = (U_j + V_j)/2, V_(j+1) = (D U_j + V_j)/2.
    let mut u = ring.one();
    let mut v = ring.one();
    let mut qk = qq.clone();
    for i in (0..k.bit_length() - 1).rev() {
        u = ring.mul(&u, &v);
        v = ring.sub(&ring.square(&v), &ring.add(&qk, &qk));
        qk = ring.square(&qk);
        if k.bit(i) {
            let u1 = ring.half(&ring.add(&u, &v));
            v = ring.half(&ring.add(&ring.mul(&dd, &u), &v));
            u = u1;
            qk = ring.mul(&qk, &qq);
        }
    }
    if ring.is_zero(&u) || ring.is_zero(&v) {
        return true;
    }
    for _ in 1..s {
        v = ring.sub(&ring.square(&v), &ring.add(&qk, &qk));
        if ring.is_zero(&v) {
            return true;
        }
        qk = ring.square(&qk);
    }
    false
}
