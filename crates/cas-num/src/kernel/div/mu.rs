//! Division by an approximate reciprocal: Newton's iteration for the
//! reciprocal ([`invert`]), then Barrett ("mu") quotient estimates block by
//! block ([`mu_div_qr`]), as in GMP's `mpn_mu_div_qr` and Brent and
//! Zimmermann's *Modern Computer Arithmetic* (§3.4.1).
//!
//! Both steps know the product they need up to a small window, so they take
//! it modulo `B^m - 1` ([`mul_wrap`]: a cyclic transform of `m` coefficients
//! instead of about twice as many), with `m` a power of two, and recover the
//! exact value from one or two low limbs computed directly.

use super::{bz_div_qr, bz_div_qr_scratch_len, reciprocal_2by1};
use crate::Limb;
use crate::kernel::mul::{fold, mul_wrap, mul_wrap_scratch_len};
use crate::kernel::{
    add_1_assign, add_assign, add_n, cmp_n, mul, mul_scratch_len, sub_1_assign, sub_assign,
    sub_n_assign,
};

#[cfg(test)]
use super::hits;

/// Smallest divisor and quotient length at which division uses Newton's
/// reciprocal instead of divide and conquer. From the tuning configuration.
pub(crate) const MU_DIV_THRESHOLD: usize = crate::tuning::MU_DIV_THRESHOLD;

/// Smallest length at which [`invert`] takes a Newton step instead of an
/// exact division (at least 3, which the step needs).
pub(crate) const INV_NEWTON_THRESHOLD: usize = crate::tuning::INV_NEWTON_THRESHOLD;

/// Two's-complement `t += v` modulo `B^t.len()` (`t.len() >= 2`).
fn add_signed(t: &mut [Limb], v: i128) {
    let m = v.unsigned_abs();
    #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    let limbs = [m as Limb, (m >> Limb::BITS) as Limb];
    if v >= 0 {
        add_assign(t, &limbs);
    } else {
        sub_assign(t, &limbs);
    }
}

/// Whether the two's-complement `t` is zero or negative.
fn is_nonpositive(t: &[Limb]) -> bool {
    t[t.len() - 1] >> (Limb::BITS - 1) == 1 || t.iter().all(|&x| x == 0)
}

/// `x = I` with `X = B^n + I` such that `A X < B^(2n) < A (X + 2)`, for the
/// normalised `n`-limb `A = a` and `x.len() == n`: `X` is
/// `floor((B^(2n) - 1) / A)` or one less. Needs `invert_scratch_len(n)`.
///
/// Below [`INV_NEWTON_THRESHOLD`], `X` is that quotient exactly. Above, with
/// `h = floor(n / 2) + 1` and `l = n - h < h`, let `Xh` be the reciprocal of
/// `A`'s top `h` limbs `Ah`, by recursion, and decrease it until
/// `E = B^(n + h) - A Xh > 0`. Then `X = Xh B^l + floor(floor(E / B^l) Xh /
/// B^(2h - l))` satisfies the bound:
///
/// - `A Xh >= Ah Xh B^l > (B^(2h) - 2 Ah) B^l >= B^(n + h) - 2A`, so `E < 2A`
///   at first, and a decrease only happens while `E <= 0`, after which
///   `E <= A`. Also `A Xh < B^(n + h) + 2 B^n`, since `Xh < 2 B^h`, so each of
///   the at most four decreases keeps `Xh >= B^h`.
/// - Exactly, `B^(2n) / A = Xh B^l + E Xh / B^(2h) + E^2 / (A B^(2h))`. The
///   computed `X` drops `(E mod B^l) Xh / B^(2h) < 2 B^(l - h)` from the middle
///   term, less than 1 from the floor, and the last term, in `(0, 4 B^(l - h))`
///   since `0 < E < 2A`. With `l < h`, `B^(2n) / A - X` lies in `(0, 1 + 6 /
///   B)`: `X < B^(2n) / A < X + 2`.
///
/// `E` is in `(-2 B^n, 2 B^n)` before the decreases. It comes from `A Xh`
/// modulo `B^m - 1` with `m >= n - 1` a power of two, and `A Xh mod B^2`: the
/// two fix `A Xh` modulo `B^2 (B^m - 1) > 4 B^n`.
pub(crate) fn invert(x: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    debug_assert!(n >= 1 && x.len() == n && a[n - 1] >> (Limb::BITS - 1) == 1);
    if n < INV_NEWTON_THRESHOLD {
        invert_exact(x, a, scratch);
        return;
    }
    let (h, l, m, w) = newton_lens(n);
    let (x_lo, x_hi) = x.split_at_mut(l);
    invert(x_hi, &a[l..], scratch);
    let (xh, rest) = scratch.split_at_mut(h + 1);
    let (e, rest) = rest.split_at_mut(w);
    xh[..h].copy_from_slice(x_hi);
    xh[h] = 1;

    // e[..m] = (B^(n + h) - A Xh) mod (B^m - 1), in [0, B^m - 1).
    {
        let (p, rest) = rest.split_at_mut(m);
        mul_wrap(p, a, xh, rest);
        e[..m].fill(0);
        e[(n + h) % m] = 1;
        if sub_n_assign(&mut e[..m], p) != 0 {
            // The borrow added B^m; the residue needs B^m - 1.
            sub_1_assign(&mut e[..m], 1);
        }
    }
    // E = e + j (B^m - 1) with |j| <= 2 B + 3, and E mod B^2 is known.
    let low = (u128::from(a[0]) * u128::from(xh[0])).wrapping_add(
        u128::from(
            a[0].wrapping_mul(xh[1])
                .wrapping_add(a[1].wrapping_mul(xh[0])),
        ) << Limb::BITS,
    );
    let residue = u128::from(e[0]) | u128::from(e[1]) << Limb::BITS;
    #[allow(clippy::cast_possible_wrap, reason = "j is small, of either sign")]
    let j = residue.wrapping_sub(low.wrapping_neg()) as i128;
    debug_assert!(
        j.unsigned_abs() <= 2 * (1 << Limb::BITS) + 3,
        "invert: j = {j}"
    );
    e[m..].fill(0);
    add_signed(&mut e[m..], j);
    add_signed(e, -j);

    while is_nonpositive(e) {
        #[cfg(test)]
        hits::INV_DECREMENT.with(|c| c.set(c.get() + 1));
        let borrow = sub_1_assign(xh, 1);
        debug_assert!(borrow == 0 && xh[h] == 1, "invert: Xh fell below B^h");
        add_assign(e, a);
    }
    debug_assert!(
        e[n + 1..].iter().all(|&x| x == 0) && e[n] <= 1,
        "invert: E >= 2A"
    );

    // I = (Xh - B^h) B^l + floor(floor(E / B^l) Xh / B^(2h - l)).
    let (u, rest) = rest.split_at_mut(2 * h + 2);
    mul(u, &e[l..=n], xh, rest);
    debug_assert!(u[2 * h + 1] == 0 && u[2 * h] < 4);
    x_lo.copy_from_slice(&u[2 * h - l..2 * h]);
    x_hi.copy_from_slice(&xh[..h]);
    let carry = add_1_assign(x_hi, u[2 * h]);
    debug_assert_eq!(carry, 0, "invert: X >= 2 B^n");
}

/// `(h, l, m, w)` of a Newton step on `n` limbs: the recursive length, the
/// limbs it leaves out, the modulus length and the limbs of `E`.
fn newton_lens(n: usize) -> (usize, usize, usize, usize) {
    let h = n / 2 + 1;
    let m = (n - 1).next_power_of_two();
    (h, n - h, m, m.max(n + 1) + 2)
}

/// Scratch limbs [`invert`] needs for `n` limbs.
pub(crate) fn invert_scratch_len(n: usize) -> usize {
    if n < INV_NEWTON_THRESHOLD {
        return if n == 1 {
            0
        } else {
            2 * n + bz_div_qr_scratch_len(2 * n, n)
        };
    }
    let (h, _, m, w) = newton_lens(n);
    let residual = m + mul_wrap_scratch_len(m, n, h + 1);
    let correction = 2 * h + 2 + mul_scratch_len(h + 1, h + 1);
    invert_scratch_len(h).max(h + 1 + w + residual.max(correction))
}

/// [`invert`] by dividing `B^(2n) - 1` by `a`.
fn invert_exact(x: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    if n == 1 {
        x[0] = reciprocal_2by1(a[0]);
        return;
    }
    let (num, rest) = scratch.split_at_mut(2 * n);
    num.fill(Limb::MAX);
    let top = bz_div_qr(x, num, a, rest);
    debug_assert_eq!(top, 1);
}

/// The quotient block length for `qn` quotient limbs and a `dn`-limb divisor
/// (GMP's choice): blocks of about `dn` limbs when the quotient is longer,
/// else two halves, or one block for a quotient of at most `dn / 3` limbs.
fn block_len(qn: usize, dn: usize) -> usize {
    if qn > dn {
        qn.div_ceil(qn.div_ceil(dn))
    } else if 3 * qn > dn {
        qn.div_ceil(2)
    } else {
        qn
    }
}

/// Division of `n` by the normalised `d` (`d.len() >= 2`, `n.len() > d.len()`)
/// with `q.len() == n.len() - d.len()`, by quotient blocks of `k` limbs (see
/// [`block_len`]) from the top. Leaves the remainder in `n[..d.len()]` and
/// returns the top quotient limb, as [`super::div_qr_normalized`].
///
/// One reciprocal serves every block: `X = B^k + I`, from [`invert`] on
/// `D' = floor(D B^(k + 1 - dn)) + 1`, `k + 1` limbs, with its low limb
/// dropped (`I = 0` if `D'` overflows to `B^(k + 1)`). Since
/// `D' B^(dn - k - 1) > D`, `X <= B^(dn + k) / D`, and the truncations of `D`
/// and of the reciprocal cost `X > B^(dn + k) / D - 1 - 6 / B`.
///
/// A block `W` of `dn + k` limbs, below `D B^k`, with top `k` limbs `Wt`, has
/// the estimate `Q^ = floor(Wt X / B^k) = Wt + floor(Wt I / B^k)`. Then `Q^ <=
/// Wt B^dn / D <= W / D`, so `Q^` never exceeds the quotient `Q`, and `Q - Q^
/// < (W mod B^dn) / D + (1 + 6 / B) + 1 < 4.1`. The last, shorter, block uses
/// the reciprocal's top limbs, whose truncation adds 1: at most 5 corrections.
pub(crate) fn mu_div_qr(q: &mut [Limb], n: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) -> Limb {
    let dn = d.len();
    let qn = n.len() - dn;
    debug_assert!(dn >= 2 && qn >= 1 && q.len() == qn);
    let qh = cmp_n(&n[qn..], d).is_ge();
    if qh {
        sub_n_assign(&mut n[qn..], d);
    }
    let k = block_len(qn, dn);
    let (inv, scratch) = scratch.split_at_mut(k);
    block_inverse(inv, d, scratch);
    let m = dn.next_power_of_two();
    let mut top = qn;
    while top > 0 {
        let kb = k.min(top);
        let lo = top - kb;
        mu_block(
            &mut q[lo..top],
            &mut n[lo..top + dn],
            d,
            &inv[k - kb..],
            m,
            scratch,
        );
        top = lo;
    }
    Limb::from(qh)
}

/// Scratch limbs [`mu_div_qr`] needs for `nn` dividend and `dn` divisor limbs.
pub(crate) fn mu_div_qr_scratch_len(nn: usize, dn: usize) -> usize {
    let qn = nn - dn;
    let k = block_len(qn, dn);
    let m = dn.next_power_of_two();
    let last = qn % k;
    let blocks = mu_block_scratch_len(k, dn, m).max(if last == 0 {
        0
    } else {
        mu_block_scratch_len(last, dn, m)
    });
    k + (2 * k + 2 + invert_scratch_len(k + 1)).max(blocks)
}

/// `inv = I`, the `k = inv.len()` limbs of the reciprocal [`mu_div_qr`] uses.
pub(crate) fn block_inverse(inv: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (k, dn) = (inv.len(), d.len());
    let (t, rest) = scratch.split_at_mut(k + 1);
    if dn > k {
        t.copy_from_slice(&d[dn - k - 1..]);
    } else {
        t[0] = 0;
        t[1..].copy_from_slice(d);
    }
    if add_1_assign(t, 1) != 0 {
        inv.fill(0);
        return;
    }
    let (y, rest) = rest.split_at_mut(k + 1);
    invert(y, t, rest);
    inv.copy_from_slice(&y[1..]);
}

/// One block of [`mu_div_qr`]: the `kb = q.len()` quotient limbs of the window
/// `w` (`d.len() + kb` limbs, below `d B^kb`), with `inv` the reciprocal's top
/// `kb` limbs and `m` the power of two at least `d.len()`. Leaves the remainder
/// in `w[..d.len()]`; the limbs above it are unspecified.
///
/// `R = W - Q^ D` is in `(0, 6D)`: positive since `X D < B^(dn + k)` strictly
/// (equality would take `D = B^dn / 2` and `X = 2 B^k`, which `k` limbs of `I`
/// cannot hold), so `Q^ D < Wt B^dn <= W`. It comes from `Q^ D` modulo
/// `B^m - 1` and the low limb `R mod B`: `R = r + j (B^m - 1)` for the residue
/// `r` and `0 <= j < 6`, and `j = r - R (mod B)`. As `R > 0`, either form of
/// a zero residue, `0` or `B^m - 1`, gives `R` this way.
fn mu_block(
    q: &mut [Limb],
    w: &mut [Limb],
    d: &[Limb],
    inv: &[Limb],
    m: usize,
    scratch: &mut [Limb],
) {
    let (kb, dn) = (q.len(), d.len());
    let (t, rest) = scratch.split_at_mut(2 * kb);
    mul(t, &w[dn..], inv, rest);
    let carry = add_n(q, &t[kb..], &w[dn..]);
    debug_assert_eq!(carry, 0, "mu_div_qr: estimate above B^kb");

    let (p, rest) = scratch.split_at_mut(m);
    let (r, rest) = rest.split_at_mut(m);
    mul_wrap(p, q, d, rest);
    fold(r, w);
    if sub_n_assign(r, p) != 0 {
        // As in `invert`: B^m - 1, not B^m.
        sub_1_assign(r, 1);
    }
    let j = r[0].wrapping_sub(w[0].wrapping_sub(q[0].wrapping_mul(d[0])));
    debug_assert!(j < 6, "mu_div_qr: j = {j}");
    #[cfg(test)]
    if j != 0 {
        hits::MU_WRAPPED.with(|c| c.set(c.get() + 1));
    }
    let mut top = if m == dn {
        j - sub_1_assign(r, j)
    } else {
        debug_assert!(j == 0 && r[dn + 1..].iter().all(|&x| x == 0));
        r[dn]
    };
    w[..dn].copy_from_slice(&r[..dn]);

    let mut corrections = 0;
    while top != 0 || cmp_n(&w[..dn], d).is_ge() {
        #[cfg(test)]
        hits::MU_CORRECTION.with(|c| c.set(c.get() + 1));
        top -= sub_n_assign(&mut w[..dn], d);
        let carry = add_1_assign(q, 1);
        debug_assert_eq!(carry, 0);
        corrections += 1;
    }
    debug_assert!(corrections <= 5, "mu_div_qr: {corrections} corrections");
}

fn mu_block_scratch_len(kb: usize, dn: usize, m: usize) -> usize {
    (2 * kb + mul_scratch_len(kb, kb)).max(2 * m + mul_wrap_scratch_len(m, kb, dn))
}
