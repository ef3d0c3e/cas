//! Multiplication by number-theoretic transforms modulo three primes below
//! 2^62, recombined by the Chinese remainder theorem. Transform lengths are
//! `2^k` or `3 2^k`: the primes have roots of order `3 2^53`.
//!
//! Each limb is one coefficient, so the product's coefficients are below
//! `min(an, bn) 2^128`, far under the primes' product (above 2^183): their
//! residues modulo each prime, from a cyclic convolution, determine them
//! exactly. Just above a power of two, the transform is half the covering
//! length instead, and the coefficients that wrap around are recomputed
//! from the operands' top limbs (see [`convolve_mod`]); with the lengths
//! `3 2^k` in between, this flattens most of the power-of-two staircase in
//! cost. The transforms are radix 2, after one radix-3 layer for the lengths
//! `3 2^k`. The forward one decimates in frequency (natural order in,
//! bit-reversed order out within each third) and the inverse one is its exact
//! mirror, decimating in time back to natural order, so no permutation pass is
//! needed. Residues stay
//! lazily in `[0, 2p)` (Harvey's butterflies), twiddle products use Shoup's
//! precomputed quotients, and the pointwise products use Montgomery reduction,
//! whose `2^-64` factor is folded into the `1 / N` scaling. Transforms longer
//! than [`NTT_BLOCK`] recurse depth-first, so each half is finished while it is
//! still in cache.
//!
//! Everything is portable `u64` / `u128` arithmetic. The default targets have
//! no vector 64x64->128 multiply, so the butterflies stay scalar; vector or
//! assembly butterflies belong to the pluggable kernels of CAS-30.
#![allow(
    clippy::cast_possible_truncation,
    reason = "splitting double limbs into limbs"
)]
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, b, f, m, n, p, r)"
)]

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, PoisonError, RwLock};

use super::div::{div_2by1, reciprocal_2by1};
use crate::Limb;

/// `(p, g)`: primes `p = c 2^53 + 1` in `(2^61, 2^62)` with `3 | c`, and an
/// element `g` of order exactly `3 2^53` modulo `p` (`x^(c / 3)` for the
/// smallest `x` that gives that order).
const PRIMES: [(Limb, Limb); 3] = [
    (0x3ea0_0000_0000_0001, 0x26bc_11b4_7562_4031),
    (0x3ae0_0000_0000_0001, 0x0238_519f_4524_a89f),
    (0x3960_0000_0000_0001, 0x368c_a63f_7d75_28a0),
];

/// `log2` of the power-of-two part of the order of the roots in [`PRIMES`]:
/// the longest transforms are `2^53` and `3 2^52`.
const ROOT_ORDER_LOG2: u32 = 53;

/// `(m, three)` for a transform length `n = m` or `n = 3 m`, `m` a power of
/// two.
fn radix_parts(n: usize) -> (usize, bool) {
    if n % 3 == 0 {
        (n / 3, true)
    } else {
        debug_assert!(n.is_power_of_two());
        (n, false)
    }
}

/// An element of order exactly `n` (`2^k` or `3 2^k`) modulo `m`, from the
/// prime's root `root` of order `3 2^53`.
fn root_of_order(m: Modulus, root: Limb, n: usize) -> Limb {
    let order = 3u64 << ROOT_ORDER_LOG2;
    debug_assert_eq!(order % n as u64, 0);
    m.pow(root, order / n as u64)
}

/// Transforms of at most this many coefficients run layer by layer; longer
/// ones split into halves first.
pub(crate) const NTT_BLOCK: usize = crate::tuning::NTT_BLOCK;

/// Arithmetic modulo one of [`PRIMES`].
#[derive(Clone, Copy, Debug)]
struct Modulus {
    p: Limb,
    two_p: Limb,
    /// `-p^-1 mod 2^64`, for Montgomery reduction.
    neg_inv: Limb,
    /// `reciprocal_2by1(4 p)`, for Shoup quotients.
    v: Limb,
}

impl Modulus {
    fn new(p: Limb) -> Self {
        debug_assert!(p >> 61 == 1 && p & 1 == 1);
        // Newton's iteration for p^-1 mod 2^64; p p = 1 mod 8 gives 3 correct
        // bits and each step doubles them.
        let mut inv = p;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(p.wrapping_mul(inv)));
        }
        debug_assert_eq!(p.wrapping_mul(inv), 1);
        Self {
            p,
            two_p: 2 * p,
            neg_inv: inv.wrapping_neg(),
            v: reciprocal_2by1(p << 2),
        }
    }

    /// Shoup's quotient `floor(w 2^64 / p)` for `w < p`.
    fn shoup(self, w: Limb) -> Limb {
        debug_assert!(w < self.p);
        // w 2^66 / 4p, with 4p normalised and w 2^2 < 4p.
        div_2by1(w << 2, 0, self.p << 2, self.v).0
    }

    /// `x w mod p` in `[0, 2p)`, for any `x`, `w < p` and `ws = shoup(w)`.
    #[inline]
    fn mul_shoup(self, x: Limb, w: Limb, ws: Limb) -> Limb {
        let q = ((u128::from(x) * u128::from(ws)) >> Limb::BITS) as Limb;
        x.wrapping_mul(w).wrapping_sub(q.wrapping_mul(self.p))
    }

    /// `t 2^-64 mod p` in `[0, 2p)`, for `t < 2^64 p`.
    #[inline]
    fn redc(self, t: u128) -> Limb {
        let m = (t as Limb).wrapping_mul(self.neg_inv);
        ((t + u128::from(m) * u128::from(self.p)) >> Limb::BITS) as Limb
    }

    /// `[0, 4p)` to `[0, 2p)`: below `2p` the subtraction wraps above `x`.
    #[inline]
    fn reduce_2p(self, x: Limb) -> Limb {
        x.min(x.wrapping_sub(self.two_p))
    }

    /// `[0, 2p)` to `[0, p)`.
    #[inline]
    fn reduce_p(self, x: Limb) -> Limb {
        x.min(x.wrapping_sub(self.p))
    }

    /// Any limb to `[0, 2p)`; `p > 2^61`, so a limb is below `8p`.
    #[inline]
    fn reduce_limb(self, x: Limb) -> Limb {
        self.reduce_2p(x.min(x.wrapping_sub(2 * self.two_p)))
    }

    /// `x y mod p` for `x, y < p`, by division; for setup only.
    fn mul(self, x: Limb, y: Limb) -> Limb {
        (u128::from(x) * u128::from(y) % u128::from(self.p)) as Limb
    }

    fn pow(self, mut x: Limb, mut e: u64) -> Limb {
        let mut r = 1;
        while e != 0 {
            if e & 1 == 1 {
                r = self.mul(r, x);
            }
            x = self.mul(x, x);
            e >>= 1;
        }
        r
    }

    fn inv(self, x: Limb) -> Limb {
        debug_assert!(x % self.p != 0);
        self.pow(x, self.p - 2)
    }
}

/// Fills `tw` (`2n` limbs, `n` a power of two) with `(w, shoup(w))` pairs:
/// pair `h + j`, for `h` a power of two below `n` and `j < h`, holds `w_2h^j`,
/// where `w_2h = omega^(n / 2h)` and `omega` has order `n`.
fn twiddles(m: Modulus, tw: &mut [Limb], omega: Limb) {
    let n = tw.len() / 2;
    let mut h = n / 2;
    if h == 0 {
        return;
    }
    let omega_s = m.shoup(omega);
    let mut w = 1;
    for pair in tw[2 * h..].chunks_exact_mut(2) {
        pair[0] = w;
        pair[1] = m.shoup(w);
        w = m.reduce_p(m.mul_shoup(w, omega, omega_s));
    }
    // w_h^j = w_2h^(2j): each level is every other pair of the one above.
    while h > 1 {
        h /= 2;
        for j in 0..h {
            let from = 4 * (h + j);
            tw.copy_within(from..from + 2, 2 * (h + j));
        }
    }
}

/// One decimation-in-frequency layer on `a` of length `2h`, with `tw` the `h`
/// twiddle pairs of the level: `(x, y) -> (x + y, (x - y) w^j)`.
#[inline]
fn dif_layer(m: Modulus, a: &mut [Limb], tw: &[Limb]) {
    let (lo, hi) = a.split_at_mut(a.len() / 2);
    for ((x, y), w) in lo.iter_mut().zip(hi).zip(tw.chunks_exact(2)) {
        let (u, v) = (*x, *y);
        *x = m.reduce_2p(u + v);
        *y = m.mul_shoup(u + m.two_p - v, w[0], w[1]);
    }
}

/// The inverse of [`dif_layer`] up to a factor 2:
/// `(x, y) -> (x + w^-j y, x - w^-j y)`. Since `w^h = -1`,
/// `w^-j = -w^(h - j)` for `0 < j < h`, so the twiddles run backwards.
#[inline]
fn dit_layer(m: Modulus, a: &mut [Limb], tw: &[Limb]) {
    let (lo, hi) = a.split_at_mut(a.len() / 2);
    let (Some((x0, lo)), Some((y0, hi))) = (lo.split_first_mut(), hi.split_first_mut()) else {
        return;
    };
    let (u, t) = (*x0, *y0);
    *x0 = m.reduce_2p(u + t);
    *y0 = m.reduce_2p(u + m.two_p - t);
    for ((x, y), w) in lo.iter_mut().zip(hi).zip(tw[2..].chunks_exact(2).rev()) {
        let (u, t) = (*x, m.mul_shoup(*y, w[0], w[1]));
        *x = m.reduce_2p(u + m.two_p - t);
        *y = m.reduce_2p(u + t);
    }
}

/// The layer of [`dif_layer`] and [`dit_layer`] with `h = 1`, where the only
/// twiddle is 1: `(x, y) -> (x + y, x - y)` on each pair.
fn layer_1(m: Modulus, a: &mut [Limb]) {
    for pair in a.chunks_exact_mut(2) {
        let (u, v) = (pair[0], pair[1]);
        pair[0] = m.reduce_2p(u + v);
        pair[1] = m.reduce_2p(u + m.two_p - v);
    }
}

/// Forward transform of `a` (a power-of-two length, residues in `[0, 2p)`),
/// output in bit-reversed order, residues in `[0, 2p)`.
fn forward(m: Modulus, a: &mut [Limb], tw: &[Limb]) {
    let n = a.len();
    if n <= NTT_BLOCK {
        let mut h = n / 2;
        while h >= 2 {
            for block in a.chunks_exact_mut(2 * h) {
                dif_layer(m, block, &tw[2 * h..4 * h]);
            }
            h /= 2;
        }
        layer_1(m, a);
        return;
    }
    dif_layer(m, a, &tw[n..2 * n]);
    let (lo, hi) = a.split_at_mut(n / 2);
    forward(m, lo, tw);
    forward(m, hi, tw);
}

/// Inverse of [`forward`] times `a.len()`: bit-reversed order in, natural
/// order out, residues in `[0, 2p)`.
fn inverse(m: Modulus, a: &mut [Limb], tw: &[Limb]) {
    let n = a.len();
    if n <= NTT_BLOCK {
        layer_1(m, a);
        let mut h = 2;
        while h < n {
            for block in a.chunks_exact_mut(2 * h) {
                dit_layer(m, block, &tw[2 * h..4 * h]);
            }
            h *= 2;
        }
        return;
    }
    let (lo, hi) = a.split_at_mut(n / 2);
    inverse(m, lo, tw);
    inverse(m, hi, tw);
    dit_layer(m, a, &tw[n..2 * n]);
}

/// Limbs of the twiddle tables for transforms up to length `n`: the radix-2
/// pairs for its power-of-two part `m`, and for `n = 3 m` the radix-3 table
/// (see [`radix3_twiddles`]).
fn tables_len(n: usize) -> usize {
    let (m, three) = radix_parts(n);
    2 * m + if three { 2 + 8 * m } else { 0 }
}

/// The twiddle tables of a transform of length `n` in `tw`
/// (`tables_len(n)` limbs); returns the radix-2 part and the radix-3 part,
/// empty for a power of two. The radix-2 table serves every power of two up to
/// `m` as a prefix.
fn tables(m: Modulus, root: Limb, n: usize, tw: &mut [Limb]) -> (&[Limb], &[Limb]) {
    let (half, three) = radix_parts(n);
    let (tw2, tw3) = tw.split_at_mut(2 * half);
    twiddles(m, tw2, root_of_order(m, root, half));
    if three {
        radix3_twiddles(m, tw3, root_of_order(m, root, n));
    }
    (tw2, if three { tw3 } else { &mut [] })
}

/// Transforms of at most this many coefficients take their twiddle tables
/// from a cache shared by every product ([`with_tables`]): computing them
/// costs a fifth of a forward transform for powers of two and over a third
/// for `3 2^k`, per prime and product. Longer transforms compute them in
/// scratch, which bounds what the cache keeps to about 18 MiB.
const TABLE_CACHE_MAX: usize = 1 << 17;

/// Per prime: the radix-2 table for the longest power-of-two part cached so
/// far. Shorter ones are its prefixes: pair `h + j` holds `w_2h^j` with
/// `w_2h` the root of order `2h`, whatever the table's length.
static RADIX2: [RwLock<Option<Arc<Vec<Limb>>>>; 3] = [const { RwLock::new(None) }; 3];

/// Per prime: the radix-3 tables by transform length.
static RADIX3: [Mutex<BTreeMap<usize, Arc<Vec<Limb>>>>; 3] =
    [const { Mutex::new(BTreeMap::new()) }; 3];

/// Calls `f` with the twiddle tables of a transform of length `n` modulo
/// prime number `prime` of [`PRIMES`] (as [`tables`] returns them): from
/// the cache up to [`TABLE_CACHE_MAX`], else computed in `tw`
/// (`tables_len(n)` limbs).
fn with_tables<R>(
    prime: usize,
    n: usize,
    tw: &mut [Limb],
    f: impl FnOnce((&[Limb], &[Limb])) -> R,
) -> R {
    let (p, root) = PRIMES[prime];
    let m = Modulus::new(p);
    if n > TABLE_CACHE_MAX {
        return f(tables(m, root, n, tw));
    }
    let (half, three) = radix_parts(n);
    let r2 = radix2_table(prime, m, root, half);
    let r3 = three.then(|| radix3_table(prime, m, root, n));
    f((
        &r2[..2 * half],
        r3.as_deref().map_or(&[][..], Vec::as_slice),
    ))
}

/// The cached radix-2 table of prime `prime`, covering power-of-two parts
/// up to at least `half`: extended (recomputed at `half`) if shorter. A
/// poisoned lock still holds a whole table, as tables are swapped in whole.
fn radix2_table(prime: usize, m: Modulus, root: Limb, half: usize) -> Arc<Vec<Limb>> {
    let lock = &RADIX2[prime];
    let covers = |t: &Option<Arc<Vec<Limb>>>| t.as_ref().filter(|t| t.len() >= 2 * half).cloned();
    if let Some(t) = covers(&lock.read().unwrap_or_else(PoisonError::into_inner)) {
        return t;
    }
    let mut guard = lock.write().unwrap_or_else(PoisonError::into_inner);
    if let Some(t) = covers(&guard) {
        return t;
    }
    let mut t = vec![0; 2 * half];
    twiddles(m, &mut t, root_of_order(m, root, half));
    let t = Arc::new(t);
    *guard = Some(Arc::clone(&t));
    t
}

/// The cached radix-3 table of prime `prime` for the length `n = 3 m`.
fn radix3_table(prime: usize, m: Modulus, root: Limb, n: usize) -> Arc<Vec<Limb>> {
    let mut map = RADIX3[prime].lock().unwrap_or_else(PoisonError::into_inner);
    Arc::clone(map.entry(n).or_insert_with(|| {
        let mut t = vec![0; 2 + 8 * (n / 3)];
        radix3_twiddles(m, &mut t, root_of_order(m, root, n));
        Arc::new(t)
    }))
}

/// The radix-3 table for `n = 3 m` from `w` of order `n` (`2 + 8 m` limbs):
/// `(ω, shoup(ω))` for the cube root of unity `ω = w^m`, then for each `j < m`
/// the pairs of `w^j` and `w^2j`, then of `w^-j` and `w^-2j`.
fn radix3_twiddles(m: Modulus, tw: &mut [Limb], w: Limb) {
    let half = (tw.len() - 2) / 8;
    let omega = m.pow(w, half as u64);
    tw[0] = omega;
    tw[1] = m.shoup(omega);
    let (fwd, inv) = tw[2..].split_at_mut(4 * half);
    for (w, quad) in [(w, fwd), (m.inv(w), inv)] {
        let ws = m.shoup(w);
        let mut x = 1;
        for q in quad.chunks_exact_mut(4) {
            let x2 = m.mul(x, x);
            q.copy_from_slice(&[x, m.shoup(x), x2, m.shoup(x2)]);
            x = m.reduce_p(m.mul_shoup(x, w, ws));
        }
    }
}

/// The radix-3 decimation-in-frequency layer on `a` (`3 m` residues in
/// `[0, 2p)`): `(x0, x1, x2) -> (x0 + x1 + x2, (x0 + ω x1 + ω^2 x2) w^j,
/// (x0 + ω^2 x1 + ω x2) w^2j)` on the thirds, with `ω^2 = -1 - ω` so that
/// one product by `ω` serves both: `x0 - x2 + ω (x1 - x2)` and
/// `x0 - x1 - ω (x1 - x2)`.
fn radix3_forward(m: Modulus, a: &mut [Limb], tw3: &[Limb]) {
    let third = a.len() / 3;
    let (omega, omega_s) = (tw3[0], tw3[1]);
    let (x0, rest) = a.split_at_mut(third);
    let (x1, x2) = rest.split_at_mut(third);
    let fwd = tw3[2..].chunks_exact(4);
    for (((x0, x1), x2), q) in x0.iter_mut().zip(x1).zip(x2).zip(fwd) {
        let (u, v, t) = (*x0, *x1, *x2);
        *x0 = m.reduce_2p(m.reduce_2p(u + v) + t);
        let d = m.mul_shoup(v + m.two_p - t, omega, omega_s);
        let s = m.reduce_2p(u + m.two_p - t);
        let r = m.reduce_2p(u + m.two_p - v);
        *x1 = m.mul_shoup(s + d, q[0], q[1]);
        *x2 = m.mul_shoup(r + m.two_p - d, q[2], q[3]);
    }
}

/// The inverse of [`radix3_forward`] up to a factor 3: with `z1 = y1 w^-j` and
/// `z2 = y2 w^-2j`, `(y0, y1, y2) -> (y0 + z1 + z2, y0 - z1 + ω (z2 - z1),
/// y0 - z2 - ω (z2 - z1))`.
fn radix3_inverse(m: Modulus, a: &mut [Limb], tw3: &[Limb]) {
    let third = a.len() / 3;
    let (omega, omega_s) = (tw3[0], tw3[1]);
    let (y0, rest) = a.split_at_mut(third);
    let (y1, y2) = rest.split_at_mut(third);
    let inv = tw3[2 + 4 * third..].chunks_exact(4);
    for (((y0, y1), y2), q) in y0.iter_mut().zip(y1).zip(y2).zip(inv) {
        let u = *y0;
        let z1 = m.mul_shoup(*y1, q[0], q[1]);
        let z2 = m.mul_shoup(*y2, q[2], q[3]);
        let d = m.mul_shoup(z2 + m.two_p - z1, omega, omega_s);
        *y0 = m.reduce_2p(m.reduce_2p(u + z1) + z2);
        *y1 = m.reduce_2p(m.reduce_2p(u + m.two_p - z1) + d);
        *y2 = m.reduce_2p(m.reduce_2p(u + m.two_p - z2) + m.two_p - d);
    }
}

/// [`forward`] for a length `m` or `3 m`, with the tables of [`tables`].
fn forward_any(m: Modulus, a: &mut [Limb], tw2: &[Limb], tw3: &[Limb]) {
    if a.len().is_power_of_two() {
        forward(m, a, tw2);
    } else {
        radix3_forward(m, a, tw3);
        for third in a.chunks_exact_mut(a.len() / 3) {
            forward(m, third, tw2);
        }
    }
}

/// [`inverse`] for a length `m` or `3 m`: the inverse of [`forward_any`] times
/// `a.len()`.
fn inverse_any(m: Modulus, a: &mut [Limb], tw2: &[Limb], tw3: &[Limb]) {
    if a.len().is_power_of_two() {
        inverse(m, a, tw2);
    } else {
        for third in a.chunks_exact_mut(a.len() / 3) {
            inverse(m, third, tw2);
        }
        radix3_inverse(m, a, tw3);
    }
}

/// `f = x mod p`, zero-padded, residues in `[0, 2p)`.
fn load(m: Modulus, f: &mut [Limb], x: &[Limb]) {
    let (head, tail) = f.split_at_mut(x.len());
    for (f, &x) in head.iter_mut().zip(x) {
        *f = m.reduce_limb(x);
    }
    tail.fill(0);
}

/// A product of `an` and `bn` limbs wraps around a transform shorter than it
/// when its top coefficients, recomputed by a product of at most this fraction
/// of the transform's power-of-two part, are all that wrap.
pub(crate) const WRAP_FRACTION: usize = crate::tuning::NTT_WRAP_FRACTION;

/// The transform length for operands of `an, bn >= 1` limbs, by
/// [`plan`].
fn transform_len(an: usize, bn: usize, radix3: bool) -> usize {
    plan(an, bn, radix3).0
}

/// `(n, cost)`: the transform length for operands of `an, bn >= 1` limbs, a
/// power of two or (with `radix3`) three times one, and its estimated cost.
///
/// The candidates are the shortest such lengths covering `len = an + bn - 1`
/// and the next shorter ones, which hold both operands and leave `e = len - n`
/// coefficients that wrap and are recomputed (see [`convolve_mod`]); `2e - 1`
/// at most the power-of-two part of `n` over [`WRAP_FRACTION`], so their
/// transform is a power of two using a prefix of the same twiddles. A length
/// costs `n` times its radix-2 layers, plus 2 for a radix-3 layer and 3 for
/// loading, pointwise products and recombination, plus the cost of the
/// wrapped coefficients; the cheapest wins.
fn plan(an: usize, bn: usize, radix3: bool) -> (usize, u64) {
    let len = an + bn - 1;
    let pow2 = len.next_power_of_two();
    let mut three = 3;
    while three < len {
        three *= 2;
    }
    let candidates = [pow2, pow2 / 2, three, three / 2];
    let mut best = (pow2, u64::MAX);
    for (i, &n) in candidates.iter().enumerate() {
        let (half, is_three) = radix_parts(n.max(1));
        if n == 0 || (i >= 2 && (!radix3 || n % 3 != 0)) || n < an.max(bn) {
            continue;
        }
        let e = len.saturating_sub(n);
        if e > 0 && 2 * e - 1 > half / WRAP_FRACTION {
            continue;
        }
        let layers = u64::from(half.trailing_zeros()) + if is_three { 2 } else { 0 } + 3;
        let mut cost = n as u64 * layers;
        if e > 0 {
            cost += plan(e, e, false).1;
        }
        if cost < best.1 {
            best = (n, cost);
        }
    }
    best
}

/// The cost of [`mul`] for operands of `an, bn >= 1` limbs in [`plan`]'s
/// model.
pub(crate) fn cost(an: usize, bn: usize) -> u64 {
    plan(an, bn, true).1
}

/// The cost of a cyclic product of transform length `n` in [`plan`]'s model.
pub(crate) fn cyclic_cost(n: usize) -> u64 {
    let (half, three) = radix_parts(n);
    n as u64 * (u64::from(half.trailing_zeros()) + if three { 2 } else { 0 } + 3)
}

/// How close the transform for operands of `an, bn >= 1` limbs comes to the
/// ideal cost `len (log2 len + 3)` of [`plan`]'s model, in `(0, 1]`: low just
/// above the steps of the staircase, where the transform is much longer than
/// the product or the wrapped part is large.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
pub(crate) fn efficiency(an: usize, bn: usize) -> f64 {
    let len = (an + bn - 1) as f64;
    let ideal = len * (len.log2() + 3.0);
    (ideal / plan(an, bn, true).1 as f64).min(1.0)
}

/// Limbs of the residue buffer for operands of `an` and `bn` limbs: the
/// product's coefficients, and room for the transform.
fn buffer_len(an: usize, bn: usize, radix3: bool) -> usize {
    transform_len(an, bn, radix3).max(an + bn - 1)
}

/// Work limbs [`convolve_mod`] needs beyond its buffer and the twiddles; the
/// wrapped coefficients take power-of-two transforms.
fn work_len(an: usize, bn: usize, square: bool, radix3: bool) -> usize {
    let n = transform_len(an, bn, radix3);
    let second = if square { 0 } else { n };
    let e = (an + bn - 1).saturating_sub(n);
    if e == 0 {
        second
    } else {
        second.max(buffer_len(e, e, false) + work_len(e, e, square, false))
    }
}

fn scratch_len(an: usize, bn: usize, square: bool) -> usize {
    let n = transform_len(an, bn, true);
    3 * buffer_len(an, bn, true) + tables_len(n) + work_len(an, bn, square, true)
}

/// Scratch limbs [`mul`] needs for operands of `an, bn >= 1` limbs.
pub(crate) fn mul_scratch_len(an: usize, bn: usize) -> usize {
    scratch_len(an, bn, false)
}

/// Scratch limbs [`sqr`] needs for an operand of `n >= 1` limbs.
pub(crate) fn sqr_scratch_len(n: usize) -> usize {
    scratch_len(n, n, true)
}

/// `r = a * b` with `r.len() == a.len() + b.len()`, both operands non-empty
/// and `scratch.len() >= mul_scratch_len(a.len(), b.len())`.
pub(crate) fn mul(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    convolve(r, a, Some(b), scratch);
}

/// `r = a^2` with `r.len() == 2 * a.len()`, `a` non-empty and
/// `scratch.len() >= sqr_scratch_len(a.len())`.
pub(crate) fn sqr(r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    convolve(r, a, None, scratch);
}

/// `r = a * b`, or `a^2` without `b`.
fn convolve(r: &mut [Limb], a: &[Limb], b: Option<&[Limb]>, scratch: &mut [Limb]) {
    let (an, bn) = (a.len(), b.map_or(a.len(), <[Limb]>::len));
    assert!(an != 0 && bn != 0, "ntt: empty operand");
    debug_assert_eq!(r.len(), an + bn);
    let n = transform_len(an, bn, true);
    assert!(
        radix_parts(n).0.trailing_zeros() <= ROOT_ORDER_LOG2,
        "ntt: operands too long"
    );
    let stride = buffer_len(an, bn, true);
    let (residues, rest) = scratch.split_at_mut(3 * stride);
    let (tw, work) = rest.split_at_mut(tables_len(n));
    for (prime, f) in residues.chunks_exact_mut(stride).enumerate() {
        let m = Modulus::new(PRIMES[prime].0);
        // The smaller transforms of the wrapped coefficients use a prefix.
        with_tables(prime, n, tw, |tw| convolve_mod(m, f, a, b, tw, work, true));
    }
    let len = an + bn - 1;
    let (f0, f12) = residues.split_at(stride);
    let (f1, f2) = f12.split_at(stride);
    let carry = recombine(&mut r[..len], &f0[..len], &f1[..len], &f2[..len]);
    r[len] = carry as Limb;
    debug_assert_eq!(carry >> Limb::BITS, 0, "ntt: product overflows");
}

/// Whether `n` is a cyclic transform length: `2^k` or `3 2^k`.
pub(crate) fn is_transform_len(n: usize) -> bool {
    n.is_power_of_two() || (n % 3 == 0 && (n / 3).is_power_of_two())
}

/// Scratch limbs [`mul_cyclic`] needs for a transform of `n` limbs.
pub(crate) fn cyclic_scratch_len(n: usize) -> usize {
    4 * n + tables_len(n)
}

/// `r = a * b mod (B^n - 1)` for `n = r.len()` a transform length
/// ([`is_transform_len`]) and operands of `1..=n` limbs, with
/// `scratch.len() >= cyclic_scratch_len(n)`. The result is a residue in
/// `[0, B^n)`: `B^n - 1` may stand for zero.
///
/// `B^n = 1` modulo `B^n - 1`, so this is the cyclic convolution of the limbs,
/// with no top coefficients to recompute. Each coefficient still sums at
/// most `min(an, bn)` limb products, within the bound the primes cover.
pub(crate) fn mul_cyclic(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    check_cyclic(n, a.len());
    assert!((1..=n).contains(&b.len()), "ntt: cyclic length mismatch");
    let (residues, rest) = scratch.split_at_mut(3 * n);
    let (tw, g) = rest.split_at_mut(tables_len(n));
    for (prime, f) in residues.chunks_exact_mut(n).enumerate() {
        let m = Modulus::new(PRIMES[prime].0);
        with_tables(prime, n, tw, |tw| {
            cyclic_mod(m, f, a, Some(b), tw, &mut g[..n]);
        });
    }
    recombine_cyclic(r, residues);
}

/// Panics unless `n` is a transform length within the roots' order and
/// `an` in `1..=n`.
fn check_cyclic(n: usize, an: usize) {
    assert!(
        is_transform_len(n) && (1..=n).contains(&an),
        "ntt: cyclic length mismatch"
    );
    assert!(
        radix_parts(n).0.trailing_zeros() <= ROOT_ORDER_LOG2,
        "ntt: operands too long"
    );
}

/// `r` from the residues of a cyclic product (`3 r.len()` limbs, one block
/// per prime), the carry out of the top wrapped around.
fn recombine_cyclic(r: &mut [Limb], residues: &[Limb]) {
    let n = r.len();
    let (f0, f12) = residues.split_at(n);
    let (f1, f2) = f12.split_at(n);
    let carry = recombine(r, f0, &f1[..n], &f2[..n]);
    // The carry is worth carry B^n = carry.
    add_wrapping(r, 0, carry as Limb);
    add_wrapping(r, 1 % n, (carry >> Limb::BITS) as Limb);
}

/// Limbs of an operand prepared for [`mul_cyclic_prepared`] at the cyclic
/// length `n`: its transform modulo each prime.
pub(crate) fn prepared_len(n: usize) -> usize {
    3 * n
}

/// Scratch limbs [`prepare_cyclic`] and [`mul_cyclic_prepared`] need at the
/// cyclic length `n`.
pub(crate) fn cyclic_prepared_scratch_len(n: usize) -> usize {
    3 * n + tables_len(n)
}

/// Prepares `b` (`1..=n` limbs) for products modulo `B^n - 1` by
/// [`mul_cyclic_prepared`], into `prep` (`prepared_len(n)` limbs): each
/// product then transforms only its other operand and the result.
pub(crate) fn prepare_cyclic(prep: &mut [Limb], n: usize, b: &[Limb], scratch: &mut [Limb]) {
    check_cyclic(n, b.len());
    assert_eq!(prep.len(), prepared_len(n), "ntt: prepared length mismatch");
    let tw = &mut scratch[..tables_len(n)];
    for (prime, g) in prep.chunks_exact_mut(n).enumerate() {
        let m = Modulus::new(PRIMES[prime].0);
        load(m, g, b);
        with_tables(prime, n, tw, |(tw2, tw3)| forward_any(m, g, tw2, tw3));
    }
}

/// [`mul_cyclic`] for `n = r.len()` and the operand `b` prepared at `n` in
/// `prep` by [`prepare_cyclic`], with `scratch.len() >=
/// cyclic_prepared_scratch_len(n)`.
pub(crate) fn mul_cyclic_prepared(r: &mut [Limb], a: &[Limb], prep: &[Limb], scratch: &mut [Limb]) {
    let n = r.len();
    check_cyclic(n, a.len());
    assert_eq!(prep.len(), prepared_len(n), "ntt: prepared length mismatch");
    let (residues, rest) = scratch.split_at_mut(3 * n);
    let tw = &mut rest[..tables_len(n)];
    for (prime, (g, f)) in prep
        .chunks_exact(n)
        .zip(residues.chunks_exact_mut(n))
        .enumerate()
    {
        let m = Modulus::new(PRIMES[prime].0);
        load(m, f, a);
        with_tables(prime, n, tw, |tw| cyclic_with(m, f, g, tw));
    }
    recombine_cyclic(r, residues);
}

/// `r += x B^at` modulo `B^n - 1` (`n = r.len()`, `at < n`): a carry out of
/// the top is worth `B^n = 1`, so it comes back in at the bottom.
fn add_wrapping(r: &mut [Limb], at: usize, x: Limb) {
    let mut carry = super::add_1_assign(&mut r[at..], x);
    while carry != 0 {
        carry = super::add_1_assign(r, carry);
    }
}

/// `f[..len] = a * b mod p` (or `a^2`), the linear convolution of
/// `len = an + bn - 1` coefficients, as residues in `[0, p)`; `f` has
/// `buffer_len(an, bn, radix3)` limbs and `tw` covers the transform length
/// (the radix-3 part only with `radix3`).
///
/// When the transform length `n` is below `len`, the cyclic convolution adds
/// each coefficient `c_(n + k)`, `k < e = len - n`, onto `c_k`. Those `c_(n + k)`
/// are the top `e` coefficients of the product of the operands' top `e`
/// limbs (every pair of indices summing to at least `n` lies there), which a
/// recursive call computes so they can be subtracted.
fn convolve_mod(
    m: Modulus,
    f: &mut [Limb],
    a: &[Limb],
    b: Option<&[Limb]>,
    tw: (&[Limb], &[Limb]),
    work: &mut [Limb],
    radix3: bool,
) {
    let (an, bn) = (a.len(), b.map_or(a.len(), <[Limb]>::len));
    let len = an + bn - 1;
    let n = transform_len(an, bn, radix3);
    let (f, top) = f.split_at_mut(n);
    cyclic_mod(m, f, a, b, tw, &mut work[..if b.is_some() { n } else { 0 }]);

    let e = len.saturating_sub(n);
    if e == 0 {
        return;
    }
    let (t, work) = work.split_at_mut(buffer_len(e, e, false));
    let (a_top, b_top) = (&a[an - e..], b.map(|b| &b[bn - e..]));
    convolve_mod(m, t, a_top, b_top, (tw.0, &[]), work, false);
    let wrapped = &t[e - 1..2 * e - 1];
    for ((x, top), &w) in f.iter_mut().zip(&mut top[..e]).zip(wrapped) {
        *x = m.reduce_p(*x + m.p - w);
        *top = w;
    }
}

/// `f = a * b mod (x^n - 1, p)` (or `a^2`) for `n = f.len()`, the cyclic
/// convolution of operands of at most `n` limbs, as residues in `[0, p)`;
/// `g` holds `b`'s transform (`n` limbs, empty when squaring).
fn cyclic_mod(
    m: Modulus,
    f: &mut [Limb],
    a: &[Limb],
    b: Option<&[Limb]>,
    (tw2, tw3): (&[Limb], &[Limb]),
    g: &mut [Limb],
) {
    load(m, f, a);
    if let Some(b) = b {
        load(m, g, b);
        forward_any(m, g, tw2, tw3);
        cyclic_with(m, f, g, (tw2, tw3));
    } else {
        forward_any(m, f, tw2, tw3);
        let (k, ks) = pointwise_scale(m, f.len());
        for x in f.iter_mut() {
            *x = m.mul_shoup(m.redc(u128::from(*x) * u128::from(*x)), k, ks);
        }
        finish_cyclic(m, f, (tw2, tw3));
    }
}

/// `f = f * g mod (x^n - 1, p)` for `f` loaded (residues in `[0, 2p)`) and
/// `g` already transformed, `n = f.len()`; as residues in `[0, p)`.
fn cyclic_with(m: Modulus, f: &mut [Limb], g: &[Limb], (tw2, tw3): (&[Limb], &[Limb])) {
    forward_any(m, f, tw2, tw3);
    let (k, ks) = pointwise_scale(m, f.len());
    for (x, &y) in f.iter_mut().zip(g) {
        *x = m.mul_shoup(m.redc(u128::from(*x) * u128::from(y)), k, ks);
    }
    finish_cyclic(m, f, (tw2, tw3));
}

/// `(k, shoup(k))` for the factor of each pointwise product of a transform of
/// length `n`: Montgomery reduction leaves `2^-64` and the inverse transform
/// `n`, so `k = 2^64 / n`: `2^(64 - log m)`, over 3 for `n = 3 m`.
fn pointwise_scale(m: Modulus, n: usize) -> (Limb, Limb) {
    let (half, three) = radix_parts(n);
    let mut k = m.pow(2, u64::from(Limb::BITS - half.trailing_zeros()));
    if three {
        k = m.mul(k, m.inv(3));
    }
    (k, m.shoup(k))
}

/// The inverse transform of the pointwise products in `f`, reduced to
/// `[0, p)`.
fn finish_cyclic(m: Modulus, f: &mut [Limb], (tw2, tw3): (&[Limb], &[Limb])) {
    inverse_any(m, f, tw2, tw3);
    for x in f.iter_mut() {
        *x = m.reduce_p(*x);
    }
}

/// `r + carry B^r.len() = sum_k x_k B^k`, where `x_k` is the coefficient with
/// residues `(f0[k], f1[k], f2[k])` modulo the three primes (Garner's
/// algorithm), one per limb of `r`; returns the carry.
fn recombine(r: &mut [Limb], f0: &[Limb], f1: &[Limb], f2: &[Limb]) -> u128 {
    let [m0, m1, m2] = PRIMES.map(|(p, _)| Modulus::new(p));
    // The primes are within a factor 2 of each other, so reduce_p reduces a
    // residue modulo one prime to a residue modulo another.
    let c01 = m1.inv(m1.reduce_p(m0.p));
    let c02 = m2.inv(m2.reduce_p(m0.p));
    let c12 = m2.inv(m2.reduce_p(m1.p));
    let (c01s, c02s, c12s) = (m1.shoup(c01), m2.shoup(c02), m2.shoup(c12));
    let p01 = u128::from(m0.p) * u128::from(m1.p);
    let (p01_lo, p01_hi) = (u128::from(p01 as Limb), p01 >> Limb::BITS);

    let mut carry: u128 = 0;
    for (((out, &x0), &x1), &x2) in r.iter_mut().zip(f0).zip(f1).zip(f2) {
        // x = x0 + p0 y1 + p0 p1 y2 with y1 < p1 and y2 < p2.
        let y1 = m1.reduce_p(m1.mul_shoup(x1 + m1.p - m1.reduce_p(x0), c01, c01s));
        let t = m2.reduce_p(m2.mul_shoup(x2 + m2.p - m2.reduce_p(x0), c02, c02s));
        let y2 = m2.reduce_p(m2.mul_shoup(t + m2.p - m2.reduce_p(y1), c12, c12s));
        let s = u128::from(x0) + u128::from(m0.p) * u128::from(y1);
        let lo = u128::from(y2) * p01_lo;
        let hi = u128::from(y2) * p01_hi;
        let low = u128::from(s as Limb) + u128::from(lo as Limb) + u128::from(carry as Limb);
        *out = low as Limb;
        carry = (low >> Limb::BITS)
            + (s >> Limb::BITS)
            + (lo >> Limb::BITS)
            + (carry >> Limb::BITS)
            + hi;
    }
    carry
}

#[cfg(test)]
mod tests;
