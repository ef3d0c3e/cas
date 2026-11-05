//! Lenstra's elliptic curve method (CAS-510).
//!
//! Each curve is a Montgomery curve `B y² = x³ + A x² + x` with Suyama's
//! parametrisation, whose group order modulo a prime is divisible by 12;
//! points are kept as `X:Z` and multiplied by the Montgomery ladder. Stage 1
//! multiplies the starting point by every prime power up to `B1`; stage 2
//! looks for one more prime in `(B1, 100 B1]` by baby steps `j Q` and giant
//! steps `k D Q`, `D = 210`, accumulating `X_k Z_j − X_j Z_k` for every
//! prime `k D ± j`. A prime factor `p` is found when the curve's order
//! modulo `p` is that smooth.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, x, z, u, v, g, k, j)"
)]

use super::sieve::sieve;
use crate::Natural;
use crate::modular::{ModRing, Residue};

#[cfg(test)]
mod tests;

/// The giant step of stage 2.
const D: u64 = 210;

/// Stage 2's bound, as a multiple of `B1`.
const STAGE_TWO: u64 = 100;

/// The least `σ` of Suyama's parametrisation; 0, 1, 3 and 5 give singular
/// or degenerate curves.
const FIRST_SIGMA: u64 = 6;

impl Natural {
    /// A non-trivial factor of `self` found by the elliptic curve method
    /// with stage-1 bound `b1` (stage 2 to `100 b1`) on at most `curves`
    /// curves, the `i`-th being Suyama's curve `σ = 6 + seed + i`; `None`
    /// when none of them finds one, always for a prime, 0 and 1.
    ///
    /// An even `self` above 2 gives 2. A factor `p` of an odd `self` is
    /// found by a curve when the curve's group order modulo `p` is a
    /// product of primes up to `b1` and one up to `100 b1`: `b1 = 2000`
    /// finds factors of about 15 digits in some 25 curves, `11 000` of 20
    /// digits in some 90, `50 000` of 25 digits in some 300.
    #[must_use]
    pub fn ecm(&self, b1: u64, curves: u32, seed: u64) -> Option<Natural> {
        if *self <= Natural::from(3u64) {
            return None;
        }
        if !self.bit(0) {
            return Some(Natural::from(2u64));
        }
        if self.primality().is_probable_prime() {
            return None;
        }
        // B1 ≥ D keeps the giant steps from k₀ ≥ 1.
        let b1 = b1.max(D);
        let plan = Plan::new(b1);
        let ring = ModRing::new(self.clone());
        (0..u64::from(curves)).find_map(|i| {
            let sigma = FIRST_SIGMA.saturating_add(seed).saturating_add(i);
            curve(&ring, self, sigma, &plan)
        })
    }
}

/// What every curve of one call shares: the prime powers of stage 1 and
/// the giant and baby steps of stage 2.
struct Plan {
    /// The largest power up to `B1` of every prime up to `B1`.
    powers: Vec<u64>,
    /// The first giant step `k₀ ≥ 1`, with `k₀ D − D/2 ≤ B1`.
    first: u64,
    /// For each prime of `(B1, B2]`, increasing: its giant step `k`, as an
    /// offset from `first`, and its baby step's index, `p = k D ± j`.
    pairs: Vec<(u64, usize)>,
    /// The baby steps `j`: odd, below `D/2`, coprime to `D`.
    babies: Vec<u64>,
}

impl Plan {
    fn new(b1: u64) -> Plan {
        let mut powers = Vec::new();
        sieve(b1, |p| {
            let mut q = p;
            while let Some(next) = q.checked_mul(p).filter(|&x| x <= b1) {
                q = next;
            }
            powers.push(q);
            true
        });
        let babies: Vec<u64> = (1..D / 2)
            .step_by(2)
            .filter(|j| crate::kernel::gcd_1(*j, D) == 1)
            .collect();
        let mut index = [usize::MAX; (D / 2) as usize];
        for (i, &j) in babies.iter().enumerate() {
            index[usize::try_from(j).expect("small")] = i;
        }
        let b2 = b1.saturating_mul(STAGE_TWO);
        let first = (b1 + D / 2) / D;
        let mut pairs = Vec::new();
        sieve(b2, |p| {
            if p > b1 {
                let k = (p + D / 2) / D;
                let j = p.abs_diff(k * D);
                pairs.push((k - first, index[usize::try_from(j).expect("small")]));
            }
            true
        });
        Plan {
            powers,
            first,
            pairs,
            babies,
        }
    }
}

/// A point `X:Z` of a Montgomery curve.
#[derive(Clone)]
struct Point {
    x: Residue,
    z: Residue,
}

/// The arithmetic of one curve modulo `n`, through `a24 = (A + 2)/4`.
struct Curve<'a> {
    ring: &'a ModRing,
    a24: Residue,
}

impl Curve<'_> {
    /// `2 P`.
    fn double(&self, p: &Point) -> Point {
        let r = self.ring;
        let s = r.square(&r.add(&p.x, &p.z));
        let d = r.square(&r.sub(&p.x, &p.z));
        let t = r.sub(&s, &d);
        Point {
            x: r.mul(&s, &d),
            z: r.mul(&t, &r.add(&d, &r.mul(&self.a24, &t))),
        }
    }

    /// `P + Q`, given `P − Q`.
    fn add(&self, p: &Point, q: &Point, difference: &Point) -> Point {
        let r = self.ring;
        let u = r.mul(&r.sub(&p.x, &p.z), &r.add(&q.x, &q.z));
        let v = r.mul(&r.add(&p.x, &p.z), &r.sub(&q.x, &q.z));
        Point {
            x: r.mul(&difference.z, &r.square(&r.add(&u, &v))),
            z: r.mul(&difference.x, &r.square(&r.sub(&u, &v))),
        }
    }

    /// `k P` for `k ≥ 1`, by the Montgomery ladder.
    fn multiply(&self, k: u64, p: &Point) -> Point {
        debug_assert!(k >= 1);
        let (mut low, mut high) = (p.clone(), self.double(p));
        for bit in (0..k.ilog2()).rev() {
            if k >> bit & 1 == 1 {
                low = self.add(&high, &low, p);
                high = self.double(&high);
            } else {
                high = self.add(&low, &high, p);
                low = self.double(&low);
            }
        }
        low
    }
}

/// The small natural `x` in the ring.
fn small(ring: &ModRing, x: u64) -> Residue {
    ring.elem(&Natural::from(x))
}

/// `gcd(x, n)` when it is a proper factor.
fn proper(ring: &ModRing, n: &Natural, x: &Residue) -> Option<Natural> {
    // n is odd, so the Montgomery factor does not change the GCD.
    let g = ring.to_natural(x).gcd(n);
    (!g.is_one() && g != *n).then_some(g)
}

/// Runs the curve of Suyama's `σ` on the odd `n`.
fn curve(ring: &ModRing, n: &Natural, sigma: u64, plan: &Plan) -> Option<Natural> {
    // u = σ² − 5, v = 4σ, the point (u³ : v³), and
    // a24 = (v − u)³ (3u + v) / (16 u³ v).
    let s = small(ring, sigma);
    let u = ring.sub(&ring.square(&s), &small(ring, 5));
    let v = ring.mul(&s, &small(ring, 4));
    let u3 = ring.mul(&ring.square(&u), &u);
    let v3 = ring.mul(&ring.square(&v), &v);
    let w = ring.sub(&v, &u);
    let numerator = ring.mul(
        &ring.mul(&ring.square(&w), &w),
        &ring.add(&ring.mul(&u, &small(ring, 3)), &v),
    );
    let denominator = ring.mul(&ring.mul(&u3, &v), &small(ring, 16));
    let Some(inverse) = ring.inv(&denominator) else {
        return proper(ring, n, &denominator);
    };
    let curve = Curve {
        ring,
        a24: ring.mul(&numerator, &inverse),
    };
    let mut q = Point { x: u3, z: v3 };
    for &power in &plan.powers {
        q = curve.multiply(power, &q);
    }
    if let found @ Some(_) = proper(ring, n, &q.z) {
        return found;
    }
    if ring.is_zero(&q.z) {
        // The order modulo every prime factor is smooth: nothing to split.
        return None;
    }
    stage_two(&curve, n, &q, plan)
}

/// Stage 2 from the stage-1 point `q`.
fn stage_two(curve: &Curve<'_>, n: &Natural, q: &Point, plan: &Plan) -> Option<Natural> {
    let ring = curve.ring;
    // The odd multiples j Q up to D/2, kept for the coprime j.
    let twice = curve.double(q);
    let mut babies = Vec::with_capacity(plan.babies.len());
    let (mut previous, mut current) = (q.clone(), q.clone());
    let mut j = 1;
    let mut wanted = plan.babies.iter().peekable();
    while let Some(&&next) = wanted.peek() {
        if j == next {
            babies.push(current.clone());
            wanted.next();
        }
        let step = if j == 1 {
            curve.add(&twice, &current, &current)
        } else {
            curve.add(&current, &twice, &previous)
        };
        previous = std::mem::replace(&mut current, step);
        j += 2;
    }
    let giant = curve.multiply(D, q);
    let mut at = curve.multiply(plan.first * D, q);
    // (k₀ − 1) D Q; for k₀ = 1 the point at infinity, which the first
    // giant step does not use.
    let mut before = if plan.first == 1 {
        at.clone()
    } else {
        curve.multiply((plan.first - 1) * D, q)
    };
    let mut k = 0;
    let mut product = ring.one();
    for &(target, index) in &plan.pairs {
        while k < target {
            let next = if plan.first == 1 && k == 0 {
                curve.double(&at)
            } else {
                curve.add(&at, &giant, &before)
            };
            before = std::mem::replace(&mut at, next);
            k += 1;
        }
        let b = &babies[index];
        let term = ring.sub(&ring.mul(&at.x, &b.z), &ring.mul(&b.x, &at.z));
        product = ring.mul(&product, &term);
    }
    proper(ring, n, &product)
}
