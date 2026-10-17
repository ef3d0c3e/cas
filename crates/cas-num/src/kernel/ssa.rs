//! Schönhage–Strassen multiplication, kept beside [`super::ntt`] so the two
//! can be compared (CAS-31).
//!
//! A product of `L` limbs is computed modulo `F_N = 2^(64 N) + 1` for some
//! `N >= L`, which leaves it unchanged. With `K = 2^k` pieces of `m = N / K`
//! limbs, `a b mod F_N` is the negacyclic convolution of the pieces evaluated
//! at `2^(64 m)`. The convolution runs in the ring `R = Z / (2^n' + 1)`, with
//! `n' = 64 t` bits large enough to hold its coefficients (`n' >= 2 64 m + k +
//! 1`, signed) and `K | n'`. There `θ = 2^(n' / K)` has order `2K`, so
//! weighting piece `i` by `θ^i` turns the negacyclic convolution into a cyclic
//! one, transformed with `ω = θ^2`. Every root is a power of two: the
//! transforms only shift, add and subtract. The `K` pointwise products are
//! again products modulo a Fermat number `2^n' + 1`, which recurse into this
//! algorithm when large ([`FERMAT_THRESHOLD`]) and otherwise use [`mul`](super::mul)
//! and one subtraction.
//!
//! Residues are `t + 1` limbs, normalised to `[0, 2^n']`: the top limb is 0,
//! or 1 with every other limb 0. The transforms are radix 2: decimation in
//! frequency forward, its mirror in time backward, so no permutation pass is
//! needed, and both recurse depth-first.
#![allow(clippy::cast_possible_truncation, reason = "bit offsets within a limb")]
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, b, k, m, n, r, t, w)"
)]

use super::{
    add_1_assign, add_assign, add_n_assign, shl, sub_1_assign, sub_assign, sub_n, sub_n_assign,
};
use crate::Limb;

/// Modular products of at least this many limbs (`t`) recurse into
/// Schönhage–Strassen; smaller ones multiply and subtract. Must stay below
/// the NTT threshold, so the comparison never runs the NTT inside.
pub(crate) const FERMAT_THRESHOLD: usize = 384;

const _: () = assert!(FERMAT_THRESHOLD < super::mul::MUL_NTT_THRESHOLD);

/// `(L, k)`: products of at least `L` limbs split into `2^k` pieces. Tuned
/// with `tune_ssa`.
const MUL_K: [(usize, u32); 8] = [
    (0, 4),
    (512, 5),
    (1024, 6),
    (2048, 7),
    (6144, 8),
    (12288, 9),
    (98304, 10),
    (393_216, 11),
];

/// `(t, k)`: modular products of at least `t` limbs split into `2^k` pieces.
const FERMAT_K: [(usize, u32); 5] = [(0, 5), (768, 6), (2048, 7), (6144, 8), (16384, 9)];

fn lookup(table: &[(usize, u32)], n: usize) -> u32 {
    table
        .iter()
        .rev()
        .find(|&&(from, _)| n >= from)
        .map_or(table[0].1, |&(_, k)| k)
}

/// The split of a modular product of `t` limbs, if it recurses.
fn inner_k(t: usize) -> Option<u32> {
    if t < FERMAT_THRESHOLD {
        return None;
    }
    let k = lookup(&FERMAT_K, t).min(t.trailing_zeros());
    (k >= 4).then_some(k)
}

/// Limbs `t` of the ring for `2^k` pieces of `m` limbs: `64 t >= 128 m + k +
/// 1`, a multiple of `2^k` bits, and (when the pointwise products recurse) a
/// multiple of their own piece count.
fn residue_len(m: usize, k: u32) -> usize {
    let t = 2 * m + 1;
    let mut align = (1usize << k).div_ceil(64);
    if t >= FERMAT_THRESHOLD {
        align = align.max(1 << lookup(&FERMAT_K, t));
    }
    t.next_multiple_of(align)
}

/// `x += F` modulo `B^(t + 1)`: brings a value in `[-F, 0)` that wrapped
/// around back into `[0, F)`.
fn add_f_wrapping(x: &mut [Limb]) {
    let t = x.len() - 1;
    add_1_assign(x, 1);
    x[t] = x[t].wrapping_add(1);
}

/// Normalises `x`, whose top limb may be any small value.
fn norm(x: &mut [Limb]) {
    let t = x.len() - 1;
    let h = std::mem::take(&mut x[t]);
    // x = lo + h 2^n' = lo - h.
    if sub_1_assign(&mut x[..t], h) != 0 {
        // lo - h wrapped to lo - h + 2^n', one below lo - h + F.
        x[t] = add_1_assign(&mut x[..t], 1);
    }
}

/// `x = x + y mod F`.
fn add_mod(x: &mut [Limb], y: &[Limb]) {
    add_n_assign(x, y);
    norm(x);
}

/// `x = x - y mod F`.
fn sub_mod(x: &mut [Limb], y: &[Limb]) {
    if sub_n_assign(x, y) != 0 {
        add_f_wrapping(x);
    }
}

/// `x = -x mod F`.
fn neg_mod(x: &mut [Limb]) {
    if x.iter().all(|&l| l == 0) {
        return;
    }
    // !x + 2 + 2^n' = F - x modulo B^(t + 1), and F - x lies in [1, 2^n'].
    for l in x.iter_mut() {
        *l = !*l;
    }
    add_f_wrapping(x);
    add_1_assign(x, 1);
}

/// `r = a 2^s mod F` for `s < 2 n'`, with `tmp` of at least `t + 1` limbs.
/// With `a 2^s = L + H 2^n'` split at bit `n'`, `r = L - H`, or `H - L` when
/// `2^n' = -1` flips the sign: one shift pass and one subtraction.
fn shift_mod(r: &mut [Limb], a: &[Limb], s: usize, tmp: &mut [Limb]) {
    let t = a.len() - 1;
    let n = 64 * t;
    debug_assert!(s < 2 * n);
    let (s, negate) = if s >= n { (s - n, true) } else { (s, false) };
    let (q, b) = (s / 64, (s % 64) as u32);
    if a[t] != 0 {
        // a = 2^n' = -1, so a 2^s = -2^s.
        r.fill(0);
        r[q] = 1 << b;
        if !negate {
            neg_mod(r);
        }
        return;
    }
    // L is a[..t - q] shifted to limb q; H is a[t - q..t] shifted, plus the
    // bits shifted out of L.
    let (l, h) = if negate {
        (&mut *tmp, &mut *r)
    } else {
        (&mut *r, &mut *tmp)
    };
    let out = shl(&mut l[q..t], &a[..t - q], b);
    l[..q].fill(0);
    h[q] = shl(&mut h[..q], &a[t - q..t], b);
    h[0] |= out;
    let borrow = if negate {
        r[q + 1..].fill(0);
        sub_assign(r, &tmp[..t])
    } else {
        r[t] = 0;
        sub_assign(r, &tmp[..=q])
    };
    if borrow != 0 {
        add_f_wrapping(r);
    }
}

/// `d = a - b` and `a = a + b` modulo F.
fn add_sub(a: &mut [Limb], b: &[Limb], d: &mut [Limb]) {
    if sub_n(d, a, b) != 0 {
        add_f_wrapping(d);
    }
    add_mod(a, b);
}

/// Forward transform of the residues in `x` (`w = t + 1` limbs each, a
/// power-of-two count) for the root `2^step` of that order; bit-reversed
/// output. `tmp` has `2 w` limbs.
fn forward(x: &mut [Limb], w: usize, step: usize, tmp: &mut [Limb]) {
    let h = x.len() / w / 2;
    if h == 0 {
        return;
    }
    let (lo, hi) = x.split_at_mut(h * w);
    let (d, t2) = tmp.split_at_mut(w);
    for (j, (a, b)) in lo
        .chunks_exact_mut(w)
        .zip(hi.chunks_exact_mut(w))
        .enumerate()
    {
        // (a, b) -> (a + b, (a - b) ω^j).
        add_sub(a, b, d);
        shift_mod(b, d, j * step, t2);
    }
    forward(lo, w, 2 * step, tmp);
    forward(hi, w, 2 * step, tmp);
}

/// Inverse of [`forward`] times the residue count: bit-reversed input,
/// natural output.
fn inverse(x: &mut [Limb], w: usize, step: usize, tmp: &mut [Limb]) {
    let h = x.len() / w / 2;
    if h == 0 {
        return;
    }
    let (lo, hi) = x.split_at_mut(h * w);
    inverse(lo, w, 2 * step, tmp);
    inverse(hi, w, 2 * step, tmp);
    let two_n = 2 * 64 * (w - 1);
    let (u, t2) = tmp.split_at_mut(w);
    for (j, (a, b)) in lo
        .chunks_exact_mut(w)
        .zip(hi.chunks_exact_mut(w))
        .enumerate()
    {
        // (a, b) -> (a + b ω^-j, a - b ω^-j), with ω^-j = 2^(2n' - j step).
        let s = if j == 0 { 0 } else { two_n - j * step };
        shift_mod(u, b, s, t2);
        add_sub(a, u, b);
    }
}

/// `r = a b mod F` for normalised residues of `t + 1` limbs (`a^2` without
/// `b`); `scratch` has `pointwise_scratch_len(t, b.is_none())` limbs.
fn mul_mod(r: &mut [Limb], a: &[Limb], b: Option<&[Limb]>, scratch: &mut [Limb]) {
    let t = a.len() - 1;
    let other = b.unwrap_or(a);
    // 2^n' = -1.
    if a[t] != 0 || other[t] != 0 {
        r.copy_from_slice(if a[t] != 0 { other } else { a });
        neg_mod(r);
        return;
    }
    if let Some(k) = inner_k(t) {
        fermat_mul(r, &a[..t], b.map(|b| &b[..t]), k, scratch);
        return;
    }
    let (p, scratch) = scratch.split_at_mut(2 * t);
    match b {
        Some(b) => super::mul(p, &a[..t], &b[..t], scratch),
        None => super::sqr(p, &a[..t], scratch),
    }
    // p = lo + hi 2^n' = lo - hi.
    r[..t].copy_from_slice(&p[..t]);
    r[t] = 0;
    if sub_n_assign(&mut r[..t], &p[t..]) != 0 {
        r[t] = add_1_assign(&mut r[..t], 1);
    }
}

fn pointwise_scratch_len(t: usize, square: bool) -> usize {
    match inner_k(t) {
        Some(k) => fermat_scratch_len(t, k, square),
        None if square => 2 * t + super::sqr_scratch_len(t),
        None => 2 * t + super::mul_scratch_len(t, t),
    }
}

/// Scratch limbs [`fermat_mul`] needs for `n` limbs and `2^k` pieces.
fn fermat_scratch_len(n: usize, k: u32, square: bool) -> usize {
    let pieces = 1 << k;
    let t = residue_len(n / pieces, k);
    let w = t + 1;
    let second = if square { 0 } else { pieces * w };
    pieces * w + second.max(n + t + 2) + 3 * w + pointwise_scratch_len(t, square)
}

/// `out = a b mod 2^(64 n) + 1` (`a^2` without `b`), `n = out.len() - 1` a
/// multiple of `2^k`, operands of at most `n` limbs; normalised output.
fn fermat_mul(out: &mut [Limb], a: &[Limb], b: Option<&[Limb]>, k: u32, scratch: &mut [Limb]) {
    let n = out.len() - 1;
    let pieces = 1usize << k;
    debug_assert!(n % pieces == 0 && a.len() <= n && b.is_none_or(|b| b.len() <= n));
    let m = n / pieces;
    let t = residue_len(m, k);
    let (w, bits) = (t + 1, 64 * t);
    let (unit, step) = (bits / pieces, 2 * bits / pieces);
    let square = b.is_none();
    let second = if square { 0 } else { pieces * w };
    let (fa, rest) = scratch.split_at_mut(pieces * w);
    let (area, rest) = rest.split_at_mut(second.max(n + t + 2));
    let (tmp, rest) = rest.split_at_mut(2 * w);
    let (prod, rest) = rest.split_at_mut(w);

    // Piece i weighted by θ^i = 2^(i n' / K).
    let decompose = |f: &mut [Limb], x: &[Limb], tmp: &mut [Limb]| {
        let (src, t2) = tmp.split_at_mut(w);
        for (i, piece) in f.chunks_exact_mut(w).enumerate() {
            let (lo, hi) = ((i * m).min(x.len()), ((i + 1) * m).min(x.len()));
            src.fill(0);
            src[..hi - lo].copy_from_slice(&x[lo..hi]);
            shift_mod(piece, src, i * unit, t2);
        }
    };
    decompose(fa, a, tmp);
    forward(fa, w, step, tmp);
    if let Some(b) = b {
        let fb = &mut area[..pieces * w];
        decompose(fb, b, tmp);
        forward(fb, w, step, tmp);
        for (x, y) in fa.chunks_exact_mut(w).zip(fb.chunks_exact(w)) {
            mul_mod(prod, x, Some(y), rest);
            x.copy_from_slice(prod);
        }
    } else {
        for x in fa.chunks_exact_mut(w) {
            mul_mod(prod, x, None, rest);
            x.copy_from_slice(prod);
        }
    }
    inverse(fa, w, step, tmp);

    // Coefficient i is 2^-k θ^-i times residue i, signed: |c_i| < 2^(n' - 1).
    // They accumulate in two's complement, at limb offsets i m.
    let acc = &mut area[..n + t + 2];
    acc.fill(0);
    for (i, x) in fa.chunks_exact(w).enumerate() {
        shift_mod(prod, x, 2 * bits - k as usize - i * unit, tmp);
        if prod[t] != 0 || prod[t - 1] >> (Limb::BITS - 1) == 1 {
            neg_mod(prod);
            sub_assign(&mut acc[i * m..], prod);
        } else {
            add_assign(&mut acc[i * m..], prod);
        }
    }

    fold(out, acc);
}

/// `out = lo - hi mod 2^(64 n) + 1`, normalised, for `acc = lo + hi 2^(64 n)`
/// with `n = out.len() - 1`, `lo` of `n` limbs and `hi` a two's complement
/// value of at most `n + 1` limbs whose magnitude is below `2^(64 n)`.
fn fold(out: &mut [Limb], acc: &mut [Limb]) {
    let n = out.len() - 1;
    let (lo, hi) = acc.split_at_mut(n);
    out[..n].copy_from_slice(lo);
    out[n] = 0;
    if hi[hi.len() - 1] >> (Limb::BITS - 1) == 0 {
        if sub_assign(out, hi) != 0 {
            add_f_wrapping(out);
        }
    } else {
        for l in hi.iter_mut() {
            *l = !*l;
        }
        add_1_assign(hi, 1);
        add_assign(out, hi);
        norm(out);
    }
}

/// `(k, n)` for a product of `l` limbs: `2^k` pieces, modulo `2^(64 n) + 1`.
fn plan(l: usize) -> (u32, usize) {
    let k = lookup(&MUL_K, l);
    (k, l.div_ceil(1 << k) << k)
}

fn scratch_len(l: usize, square: bool) -> usize {
    let (k, n) = plan(l);
    n + 1 + fermat_scratch_len(n, k, square)
}

/// Scratch limbs [`mul`] needs for operands of `an` and `bn` limbs.
pub(crate) fn mul_scratch_len(an: usize, bn: usize) -> usize {
    scratch_len(an + bn, false)
}

/// Scratch limbs [`sqr`] needs for an operand of `n` limbs.
pub(crate) fn sqr_scratch_len(n: usize) -> usize {
    scratch_len(2 * n, true)
}

/// `r = a * b` with `r.len() == a.len() + b.len()` and
/// `scratch.len() >= mul_scratch_len(a.len(), b.len())`.
pub(crate) fn mul(r: &mut [Limb], a: &[Limb], b: &[Limb], scratch: &mut [Limb]) {
    product(r, a, Some(b), scratch);
}

/// `r = a^2` with `r.len() == 2 * a.len()` and
/// `scratch.len() >= sqr_scratch_len(a.len())`.
pub(crate) fn sqr(r: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    product(r, a, None, scratch);
}

fn product(r: &mut [Limb], a: &[Limb], b: Option<&[Limb]>, scratch: &mut [Limb]) {
    let l = r.len();
    debug_assert_eq!(l, a.len() + b.map_or(a.len(), <[Limb]>::len));
    let (k, n) = plan(l);
    // The product is below 2^(64 n), so reducing it modulo 2^(64 n) + 1
    // leaves it unchanged.
    let (out, rest) = scratch.split_at_mut(n + 1);
    fermat_mul(out, a, b, k, rest);
    r.copy_from_slice(&out[..l]);
    debug_assert!(out[l..].iter().all(|&x| x == 0));
}

#[cfg(test)]
mod tests;
