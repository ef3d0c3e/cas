//! Multiplication by number-theoretic transforms modulo three primes below
//! 2^62, recombined by the Chinese remainder theorem.
//!
//! Each limb is one coefficient, so the product's coefficients are below
//! `min(an, bn) 2^128`, far under the primes' product (above 2^183): their
//! residues modulo each prime, from a cyclic convolution, determine them
//! exactly. Just above a power of two, the transform is half the covering
//! length instead, and the coefficients that wrap around are recomputed
//! from the operands' top limbs (see [`convolve_mod`]); this flattens most of
//! the power-of-two staircase in cost. The transforms are radix 2. The forward
//! one decimates in frequency (natural order in, bit-reversed order out) and
//! the inverse one is its exact mirror, decimating in time from bit-reversed
//! order back to natural order, so no permutation pass is needed. Residues stay
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

use super::div::{div_2by1, reciprocal_2by1};
use crate::Limb;

/// `(p, g)`: primes `p = c 2^54 + 1` in `(2^61, 2^62)`, and an element `g` of
/// order exactly `2^54` modulo `p` (a quadratic non-residue raised to `c`).
const PRIMES: [(Limb, Limb); 3] = [
    (0x3a00_0000_0000_0001, 0x0fc2_62fd_2103_9b12),
    (0x2c40_0000_0000_0001, 0x2a67_1570_b78e_304e),
    (0x28c0_0000_0000_0001, 0x0127_0e41_54da_e99b),
];

/// `log2` of the order of the roots in [`PRIMES`]: the longest transform.
const ROOT_ORDER_LOG2: u32 = 54;

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

/// `f = x mod p`, zero-padded, residues in `[0, 2p)`.
fn load(m: Modulus, f: &mut [Limb], x: &[Limb]) {
    let (head, tail) = f.split_at_mut(x.len());
    for (f, &x) in head.iter_mut().zip(x) {
        *f = m.reduce_limb(x);
    }
    tail.fill(0);
}

/// A product of `an` and `bn` limbs wraps around a transform of half the
/// covering power of two when its top coefficients, recomputed by a product of
/// at most this fraction of the transform's length, are all that wrap.
pub(crate) const WRAP_FRACTION: usize = crate::tuning::NTT_WRAP_FRACTION;

/// The transform length for operands of `an, bn >= 1` limbs: the power of two
/// `n >= len = an + bn - 1`, or half of it when the `e = len - n / 2`
/// coefficients that wrap are cheap to recompute (see [`convolve_mod`]).
fn transform_len(an: usize, bn: usize) -> usize {
    let len = an + bn - 1;
    let n = len.next_power_of_two();
    let half = n / 2;
    if half > 0 && an.max(bn) <= half && 2 * (len - half) - 1 <= half / WRAP_FRACTION {
        half
    } else {
        n
    }
}

/// Limbs of the residue buffer for operands of `an` and `bn` limbs: the
/// product's coefficients, and room for the transform.
fn buffer_len(an: usize, bn: usize) -> usize {
    transform_len(an, bn).max(an + bn - 1)
}

/// Work limbs [`convolve_mod`] needs beyond its buffer and the twiddles.
fn work_len(an: usize, bn: usize, square: bool) -> usize {
    let n = transform_len(an, bn);
    let second = if square { 0 } else { n };
    let e = (an + bn - 1).saturating_sub(n);
    if e == 0 {
        second
    } else {
        second.max(buffer_len(e, e) + work_len(e, e, square))
    }
}

fn scratch_len(an: usize, bn: usize, square: bool) -> usize {
    3 * buffer_len(an, bn) + 2 * transform_len(an, bn) + work_len(an, bn, square)
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
    let n = transform_len(an, bn);
    let log_n = n.trailing_zeros();
    assert!(log_n <= ROOT_ORDER_LOG2, "ntt: operands too long");
    let stride = buffer_len(an, bn);
    let (residues, rest) = scratch.split_at_mut(3 * stride);
    let (tw, work) = rest.split_at_mut(2 * n);
    for (&(p, root), f) in PRIMES.iter().zip(residues.chunks_exact_mut(stride)) {
        let m = Modulus::new(p);
        // The smaller transforms of the wrapped coefficients use a prefix.
        twiddles(m, tw, m.pow(root, 1 << (ROOT_ORDER_LOG2 - log_n)));
        convolve_mod(m, f, a, b, tw, work);
    }
    let len = an + bn - 1;
    let (f0, f12) = residues.split_at(stride);
    let (f1, f2) = f12.split_at(stride);
    recombine(r, &f0[..len], &f1[..len], &f2[..len]);
}

/// `f[..len] = a * b mod p` (or `a^2`), the linear convolution of
/// `len = an + bn - 1` coefficients, as residues in `[0, p)`; `f` has
/// `buffer_len(an, bn)` limbs and `tw` covers the transform length.
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
    tw: &[Limb],
    work: &mut [Limb],
) {
    let (an, bn) = (a.len(), b.map_or(a.len(), <[Limb]>::len));
    let len = an + bn - 1;
    let n = transform_len(an, bn);
    let log_n = n.trailing_zeros();
    let (f, top) = f.split_at_mut(n);
    load(m, f, a);
    forward(m, f, tw);
    // Montgomery reduction leaves 2^-64 and the inverse transform n, so each
    // pointwise product is scaled by 2^64 / n = 2^(64 - log n).
    let k = m.pow(2, u64::from(Limb::BITS - log_n));
    let ks = m.shoup(k);
    if let Some(b) = b {
        let g = &mut work[..n];
        load(m, g, b);
        forward(m, g, tw);
        for (x, &y) in f.iter_mut().zip(&*g) {
            *x = m.mul_shoup(m.redc(u128::from(*x) * u128::from(y)), k, ks);
        }
    } else {
        for x in f.iter_mut() {
            *x = m.mul_shoup(m.redc(u128::from(*x) * u128::from(*x)), k, ks);
        }
    }
    inverse(m, f, tw);
    for x in f.iter_mut() {
        *x = m.reduce_p(*x);
    }

    let e = len.saturating_sub(n);
    if e == 0 {
        return;
    }
    let (t, work) = work.split_at_mut(buffer_len(e, e));
    convolve_mod(m, t, &a[an - e..], b.map(|b| &b[bn - e..]), tw, work);
    let wrapped = &t[e - 1..2 * e - 1];
    for ((x, top), &w) in f.iter_mut().zip(&mut top[..e]).zip(wrapped) {
        *x = m.reduce_p(*x + m.p - w);
        *top = w;
    }
}

/// `r = sum_k x_k B^k`, where `x_k` is the coefficient with residues
/// `(f0[k], f1[k], f2[k])` modulo the three primes (Garner's algorithm); the
/// last limb of `r` only takes the carry.
fn recombine(r: &mut [Limb], f0: &[Limb], f1: &[Limb], f2: &[Limb]) {
    let [m0, m1, m2] = PRIMES.map(|(p, _)| Modulus::new(p));
    // The primes are within a factor 2 of each other, so reduce_p reduces a
    // residue modulo one prime to a residue modulo another.
    let c01 = m1.inv(m1.reduce_p(m0.p));
    let c02 = m2.inv(m2.reduce_p(m0.p));
    let c12 = m2.inv(m2.reduce_p(m1.p));
    let (c01s, c02s, c12s) = (m1.shoup(c01), m2.shoup(c02), m2.shoup(c12));
    let p01 = u128::from(m0.p) * u128::from(m1.p);
    let (p01_lo, p01_hi) = (u128::from(p01 as Limb), p01 >> Limb::BITS);

    let coefficients = r.len() - 1;
    let (r, top) = r.split_at_mut(coefficients);
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
    top[0] = carry as Limb;
    debug_assert_eq!(carry >> Limb::BITS, 0, "ntt: product overflows");
}

#[cfg(test)]
mod tests;
