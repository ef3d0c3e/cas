//! Miller–Rabin, the strong Lucas test, Baillie–PSW, `next_prime` and
//! `previous_prime`.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, d, s, q, u, v, k)"
)]

use super::Primality;
use super::small::{SIEVE_LIMIT, TRIAL_LIMIT, groups, is_small_prime, primes};
use crate::modular::{Fixed, ModDouble, ModFixed, ModRing, ModWord, Residue, RingWork};
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
        ring.is_zero(&a) || strong_probable_prime(&RingWork::new(&ring), self, &a)
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

    /// The greatest prime below `self`, `None` below 3. Above 2^64 it is a
    /// probable prime in the sense of [`Natural::primality`], and so is every
    /// candidate it passes over composite.
    #[must_use]
    pub fn previous_prime(&self) -> Option<Natural> {
        if let Ok(n) = u64::try_from(self) {
            return previous_prime_u64(n).map(Natural::from);
        }
        // The greatest odd value below self; primality decides below 2^64
        // too, should the search cross it.
        let below = self - 1u64;
        let mut c = if below.bit(0) { below } else { below - 1u64 };
        while !c.primality().is_probable_prime() {
            c -= 2u64;
        }
        Some(c)
    }
}

/// The greatest prime below `n`, `None` for `n ≤ 2`.
fn previous_prime_u64(n: u64) -> Option<u64> {
    match n {
        0..=2 => None,
        3 => Some(2),
        _ => {
            // The greatest odd value below n, at least 3, which is prime.
            let mut c = (n - 1) | 1;
            if c >= n {
                c -= 2;
            }
            while !is_prime_u64(c) {
                c -= 2;
            }
            Some(c)
        }
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
/// modulus can run them on `u128` ([`ModDouble`]), one of three to eight
/// limbs on stack arrays ([`ModFixed`]) and any other on [`ModRing`] with
/// reused buffers ([`RingWork`]). Each operation overwrites its first
/// argument, so the loops allocate nothing per step.
trait TestRing {
    type Elem: Clone + PartialEq;
    fn small(&self, a: i64) -> Self::Elem;
    fn elem(&self, a: &Natural) -> Self::Elem;
    fn one(&self) -> Self::Elem;
    fn is_zero(&self, x: &Self::Elem) -> bool;
    fn add(&self, r: &mut Self::Elem, x: &Self::Elem, y: &Self::Elem);
    fn sub(&self, r: &mut Self::Elem, x: &Self::Elem, y: &Self::Elem);
    fn neg(&self, r: &mut Self::Elem, x: &Self::Elem);
    fn mul(&self, r: &mut Self::Elem, x: &Self::Elem, y: &Self::Elem);
    fn square(&self, r: &mut Self::Elem, x: &Self::Elem);
    fn pow(&self, x: &Self::Elem, e: &Natural) -> Self::Elem;
    /// One doubling step of the Lucas sequence `U` for `P = 1`, `Q = q`:
    /// `(u, u1) = (U_j, U_(j+1))` becomes `(U_(2j), U_(2j+1))`, or
    /// `(U_(2j+1), U_(2j+2))` if `up`; `t` is scratch.
    fn lucas_double(
        &self,
        u: &mut Self::Elem,
        u1: &mut Self::Elem,
        q: i64,
        up: bool,
        t: &mut [Self::Elem; 4],
    ) {
        lucas_double(self, u, u1, q, up, t);
    }
}

impl TestRing for RingWork<'_> {
    type Elem = Residue;
    fn small(&self, a: i64) -> Residue {
        self.ring().elem_integer(&Integer::from(a))
    }
    fn elem(&self, a: &Natural) -> Residue {
        self.ring().elem(a)
    }
    fn one(&self) -> Residue {
        self.ring().one()
    }
    fn is_zero(&self, x: &Residue) -> bool {
        self.ring().is_zero(x)
    }
    fn add(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        RingWork::add(self, r, x, y);
    }
    fn sub(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        RingWork::sub(self, r, x, y);
    }
    fn neg(&self, r: &mut Residue, x: &Residue) {
        RingWork::neg(self, r, x);
    }
    fn mul(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        RingWork::mul(self, r, x, y);
    }
    fn square(&self, r: &mut Residue, x: &Residue) {
        RingWork::square(self, r, x);
    }
    fn pow(&self, x: &Residue, e: &Natural) -> Residue {
        self.ring().pow(x, e)
    }
    fn lucas_double(
        &self,
        u: &mut Residue,
        u1: &mut Residue,
        q: i64,
        up: bool,
        t: &mut [Residue; 4],
    ) {
        if q.unsigned_abs() < 1 << 8 {
            RingWork::lucas_double(self, u, u1, q, up);
        } else {
            lucas_double(self, u, u1, q, up, t);
        }
    }
}

impl TestRing for ModDouble {
    type Elem = u128;
    fn small(&self, a: i64) -> u128 {
        // |a| is far below n > 2^64.
        let m = self.to_mont(u128::from(a.unsigned_abs()));
        if a < 0 { ModDouble::neg(*self, m) } else { m }
    }
    fn elem(&self, a: &Natural) -> u128 {
        let limbs = a.limbs();
        let a = u128::from(limbs.first().copied().unwrap_or(0))
            | u128::from(limbs.get(1).copied().unwrap_or(0)) << 64;
        debug_assert!(limbs.len() <= 2);
        self.to_mont(a)
    }
    fn one(&self) -> u128 {
        ModDouble::one(*self)
    }
    fn is_zero(&self, x: &u128) -> bool {
        *x == 0
    }
    fn add(&self, r: &mut u128, x: &u128, y: &u128) {
        *r = ModDouble::add(*self, *x, *y);
    }
    fn sub(&self, r: &mut u128, x: &u128, y: &u128) {
        *r = ModDouble::sub(*self, *x, *y);
    }
    fn neg(&self, r: &mut u128, x: &u128) {
        *r = ModDouble::neg(*self, *x);
    }
    fn mul(&self, r: &mut u128, x: &u128, y: &u128) {
        *r = ModDouble::mul(*self, *x, *y);
    }
    fn square(&self, r: &mut u128, x: &u128) {
        *r = ModDouble::mul(*self, *x, *x);
    }
    fn pow(&self, x: &u128, e: &Natural) -> u128 {
        self.pow_mont(*x, e)
    }
}

impl<const K: usize> TestRing for ModFixed<K> {
    type Elem = Fixed<K>;
    fn small(&self, a: i64) -> Fixed<K> {
        let m = self.elem(&Natural::from(a.unsigned_abs()));
        if a < 0 { ModFixed::neg(self, m) } else { m }
    }
    fn elem(&self, a: &Natural) -> Fixed<K> {
        ModFixed::elem(self, a)
    }
    fn one(&self) -> Fixed<K> {
        ModFixed::one(self)
    }
    fn is_zero(&self, x: &Fixed<K>) -> bool {
        *x == Fixed::default()
    }
    fn add(&self, r: &mut Fixed<K>, x: &Fixed<K>, y: &Fixed<K>) {
        *r = ModFixed::add(self, *x, *y);
    }
    fn sub(&self, r: &mut Fixed<K>, x: &Fixed<K>, y: &Fixed<K>) {
        *r = ModFixed::sub(self, *x, *y);
    }
    fn neg(&self, r: &mut Fixed<K>, x: &Fixed<K>) {
        *r = ModFixed::neg(self, *x);
    }
    fn mul(&self, r: &mut Fixed<K>, x: &Fixed<K>, y: &Fixed<K>) {
        *r = ModFixed::mul(self, *x, *y);
    }
    fn square(&self, r: &mut Fixed<K>, x: &Fixed<K>) {
        *r = ModFixed::square(self, *x);
    }
    fn pow(&self, x: &Fixed<K>, e: &Natural) -> Fixed<K> {
        ModFixed::pow(self, *x, e)
    }
}

/// Baillie–PSW for odd `n` above 2^64 without small factors, on `u128` for
/// two limbs and stack arrays up to eight.
fn baillie_psw(n: &Natural) -> bool {
    if let [lo, hi] = *n.limbs() {
        let ring = ModDouble::new(u128::from(lo) | u128::from(hi) << 64);
        return baillie_psw_in(&ring, n);
    }
    match n.limb_len() {
        3 => baillie_psw_in(&ModFixed::<3>::new(n), n),
        4 => baillie_psw_in(&ModFixed::<4>::new(n), n),
        5 => baillie_psw_in(&ModFixed::<5>::new(n), n),
        6 => baillie_psw_in(&ModFixed::<6>::new(n), n),
        7 => baillie_psw_in(&ModFixed::<7>::new(n), n),
        8 => baillie_psw_in(&ModFixed::<8>::new(n), n),
        _ => baillie_psw_in(&RingWork::new(&ModRing::new(n.clone())), n),
    }
}

fn baillie_psw_in<R: TestRing>(ring: &R, n: &Natural) -> bool {
    let two = ring.small(2);
    strong_probable_prime(ring, n, &two) && !n.is_perfect_square() && strong_lucas(ring, n)
}

/// Miller–Rabin to base `a` (an element of `ring`, modulo odd `n > 1`).
fn strong_probable_prime<R: TestRing>(ring: &R, n: &Natural, a: &R::Elem) -> bool {
    let n1 = n - 1u64;
    let s = n1.trailing_zeros().expect("n > 1");
    let one = ring.one();
    let mut minus_one = one.clone();
    ring.neg(&mut minus_one, &one);
    let mut y = ring.pow(a, &(&n1 >> s));
    if y == one || y == minus_one {
        return true;
    }
    let mut t = y.clone();
    for _ in 1..s {
        ring.square(&mut t, &y);
        std::mem::swap(&mut y, &mut t);
        if y == minus_one {
            return true;
        }
        if y == one {
            return false;
        }
    }
    false
}

/// `r = c x`, by doubling and adding for a small `|c|` (as Selfridge's `Q`
/// almost always is), else by a product; `t` is scratch.
fn times<R: TestRing + ?Sized>(ring: &R, r: &mut R::Elem, x: &R::Elem, c: i64, t: &mut R::Elem) {
    let m = c.unsigned_abs();
    if m == 0 || m >= 1 << 8 {
        ring.mul(r, x, &ring.small(c));
        return;
    }
    // From the top bit of m, which gives x itself.
    r.clone_from(x);
    for i in (0..u64::BITS - 1 - m.leading_zeros()).rev() {
        ring.add(t, r, r);
        if m >> i & 1 == 1 {
            ring.add(r, t, x);
        } else {
            std::mem::swap(r, t);
        }
    }
    if c < 0 {
        ring.neg(t, r);
        std::mem::swap(r, t);
    }
}

/// [`TestRing::lucas_double`] by three squares, each reduced: with
/// `P = 1`, `Q U_(j-1) = U_j - U_(j+1)`, so for `A = U_j^2`, `B =
/// U_(j+1)^2` and `C = (U_j - U_(j+1))^2 = Q^2 U_(j-1)^2`,
/// `U_(2j+1) = U_(j+1)^2 - Q U_j^2 = B - Q A`, `U_(2j) = U_(2j+1) +
/// Q U_(2j-1) = B - C` (from `U_(2j-1) = U_j^2 - Q U_(j-1)^2`), and
/// `U_(2j+2) = U_(2j+1) - Q U_(2j)`.
fn lucas_double<R: TestRing + ?Sized>(
    ring: &R,
    u: &mut R::Elem,
    u1: &mut R::Elem,
    q: i64,
    up: bool,
    [a, b, c, t]: &mut [R::Elem; 4],
) {
    ring.square(a, u);
    ring.square(b, u1);
    ring.sub(t, u, u1);
    ring.square(c, t);
    ring.sub(u, b, c);
    times(ring, t, a, q, c);
    ring.sub(u1, b, t);
    if up {
        times(ring, t, u, q, c);
        ring.sub(u, u1, t);
        std::mem::swap(u, u1);
    }
}

/// `4^(-1) mod n` for odd `n`: `(n + 1) / 4` or `(3n + 1) / 4`, whichever
/// is whole.
fn quarter(n: &Natural) -> Natural {
    let m = if n.bit(1) { n + 1u64 } else { n * 3u64 + 1u64 };
    m >> 2
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
    let n1 = n + 1u64;
    let s = n1.trailing_zeros().expect("n + 1 > 0");
    let k = &n1 >> s;
    // The pair (U_j, U_(j+1)) from j = 1, doubling along the bits of k.
    let mut u = ring.one();
    let mut u1 = ring.one();
    let mut tmp = [u.clone(), u.clone(), u.clone(), u.clone()];
    for i in (0..k.bit_length() - 1).rev() {
        ring.lucas_double(&mut u, &mut u1, q, k.bit(i), &mut tmp);
    }
    let [a, b, c, t] = &mut tmp;
    if ring.is_zero(&u) {
        return true;
    }
    // V_k = 2 U_(k+1) - P U_k.
    let mut v = u.clone();
    ring.add(t, &u1, &u1);
    ring.sub(&mut v, t, &u);
    if ring.is_zero(&v) {
        return true;
    }
    if s == 1 {
        return false;
    }
    // Q^k = (V_k^2 - D U_k^2) / 4, for V_(2m) = V_m^2 - 2 Q^m.
    ring.square(a, &v);
    ring.square(b, &u);
    times(ring, c, b, d, t);
    ring.sub(b, a, c);
    let mut qk = u1;
    ring.mul(&mut qk, b, &ring.elem(&quarter(n)));
    for _ in 1..s {
        ring.square(a, &v);
        ring.add(b, &qk, &qk);
        ring.sub(&mut v, a, b);
        if ring.is_zero(&v) {
            return true;
        }
        ring.square(a, &qk);
        std::mem::swap(&mut qk, a);
    }
    false
}
