//! Approximate reciprocal, quotient, inverse square root and square root with
//! proven error bounds, for ball midpoints (CAS-DOC-5 §5.2, §10): a ball needs
//! a value near the true one and a bound on the distance, not the correctly
//! rounded value, so these skip the final correction step of an exact
//! algorithm.
//!
//! Every result `r` is a **lower bound**: `r <= floor(t) <= r + e` for the true
//! value `t` and the error `e` given by the matching `*_error` function, in
//! units of the result's last limb. Where an exact algorithm is cheaper
//! (below the Newton thresholds), the dispatching functions use it and the
//! error is 0.
//!
//! The reciprocal and the inverse square root are Newton's iterations
//! (\[BZ10\], chapter 3). The inverse square root's residual is wrapped
//! modulo `B^m - 1` where the cyclic transform applies: the wrap-around
//! trick, as in GMP's `mpn_invertappr`. The exact square root is
//! Zimmermann's Karatsuba square root \[Zim99\]. Algorithms and proofs are
//! documented in `doc/algorithms/kernels.md`; bracketed keys refer to
//! `doc/algorithms/references.md`.

use super::mu::{INV_NEWTON_THRESHOLD, MU_DIV_THRESHOLD, invert, invert_scratch_len};
use super::{div_rem, div_rem_scratch_len};
use crate::kernel::mul::mulhi::{mulhi, mulhi_error, mulhi_scratch_len};
use crate::kernel::mul::mullo::{mullo, mullo_scratch_len};
use crate::kernel::mul::{mul_wrap, mul_wrap_scratch_len, wraps_cyclically};
use crate::kernel::{
    add_1_assign, add_assign, add_n_assign, mul, mul_scratch_len, shl, shr, sqr, sqr_scratch_len,
    sub_1_assign, sub_assign, sub_n_assign,
};
use crate::{Limb, Natural};

/// Smallest length at which [`rsqrt_approx`] takes a Newton step instead of
/// computing `floor(sqrt(B^(3n) / A))` exactly (at least 3, which the step
/// needs: its recursive length `n / 2 + 1` must be below `n`).
const RSQRT_NEWTON_THRESHOLD: usize = 3;

/// Smallest length at which [`sqrt_approx`] goes through [`rsqrt_approx`]
/// instead of the exact integer square root: never. Measured (`tune_approx`,
/// CAS-90), the Newton path takes 1.4–2.1× the exact Karatsuba root from 64
/// to 8192 limbs (1.59× at 8192), and is within 10% only around 16–32
/// limbs. Its residual `A Yh^2` is wrapped where the cyclic transform
/// applies (from about 1400 limbs in the default tuning), which saved only
/// 3–5% with power-of-two moduli; CAS-98 takes `3 2^k` moduli too.
const SQRT_NEWTON_THRESHOLD: Option<usize> = None;

fn sqrt_by_newton(n: usize) -> bool {
    SQRT_NEWTON_THRESHOLD.is_some_and(|t| n >= t)
}

/// `x = X` for the normalised `n`-limb `D = d` (top bit set) and
/// `x.len() == n + 1`, with `X <= floor(B^(2n) / D) <= X + recip_approx_error(n)`
/// and `B^n <= X < 2 B^n`. Needs `recip_approx_scratch_len(n)`.
///
/// This is division's Newton reciprocal (`mu::invert`): `D X < B^(2n) < D (X +
/// 2)`, and `X` is `floor((B^(2n) - 1) / D)` exactly below
/// `INV_NEWTON_THRESHOLD` limbs, or one less above. That floor is one below
/// `floor(B^(2n) / D)` when `D` is a power of two.
pub fn recip_approx(x: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = d.len();
    assert!(n >= 1 && x.len() == n + 1, "recip_approx: length mismatch");
    assert!(
        d[n - 1] >> (Limb::BITS - 1) == 1,
        "recip_approx: divisor not normalised"
    );
    assert!(
        scratch.len() >= recip_approx_scratch_len(n),
        "recip_approx: scratch too short"
    );
    invert(&mut x[..n], d, scratch);
    x[n] = 1;
}

/// Scratch limbs [`recip_approx`] needs for `n` limbs.
#[must_use]
pub fn recip_approx_scratch_len(n: usize) -> usize {
    invert_scratch_len(n)
}

/// The error bound of [`recip_approx`] for `n` limbs: 1 for the exact
/// quotient of `B^(2n) - 1`, 2 for Newton's.
#[must_use]
pub fn recip_approx_error(n: usize) -> usize {
    1 + usize::from(n >= INV_NEWTON_THRESHOLD)
}

/// `q = Q` with `Q <= floor(A B^n / D) <= Q + div_approx_error(n)`, for
/// `n`-limb `A = a` and normalised `D = d`, and `q.len() == n + 1`
/// (`Q < 2 B^n`). Needs `div_approx_scratch_len(n)`.
///
/// Below `MU_DIV_THRESHOLD` limbs the quotient is exact; above, see
/// `div_approx_newton`.
pub fn div_approx(q: &mut [Limb], a: &[Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = d.len();
    assert!(
        n >= 1 && a.len() == n && q.len() == n + 1,
        "div_approx: length mismatch"
    );
    assert!(
        d[n - 1] >> (Limb::BITS - 1) == 1,
        "div_approx: divisor not normalised"
    );
    assert!(
        scratch.len() >= div_approx_scratch_len(n),
        "div_approx: scratch too short"
    );
    if n >= MU_DIV_THRESHOLD {
        div_approx_newton(q, a, d, scratch);
        return;
    }
    let (num, rest) = scratch.split_at_mut(2 * n);
    let (r, rest) = rest.split_at_mut(n);
    num[..n].fill(0);
    num[n..].copy_from_slice(a);
    div_rem(q, r, num, d, rest);
}

/// Scratch limbs [`div_approx`] needs for `n` limbs.
#[must_use]
pub fn div_approx_scratch_len(n: usize) -> usize {
    if n >= MU_DIV_THRESHOLD {
        div_approx_newton_scratch_len(n)
    } else {
        3 * n + div_rem_scratch_len(2 * n, n)
    }
}

/// The error bound of [`div_approx`] for `n` limbs.
#[must_use]
pub fn div_approx_error(n: usize) -> usize {
    if n >= MU_DIV_THRESHOLD {
        div_approx_newton_error(n)
    } else {
        0
    }
}

/// [`div_approx`] by the reciprocal: `Q = A + mulhi(A, I)` for
/// `X = B^n + I` from [`recip_approx`].
///
/// `A B^n / D - A X / B^n = A (B^(2n) / D - X) / B^n < 2`, since
/// `B^(2n) / D < X + 2` and `A < B^n`; `A X / B^n = A + A I / B^n`, and
/// `mulhi` loses below `mulhi_error(n) + 1` of the last term. So
/// `A B^n / D - Q < 3 + mulhi_error(n)`, and `Q` never exceeds it.
pub(crate) fn div_approx_newton(q: &mut [Limb], a: &[Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = d.len();
    let (x, rest) = scratch.split_at_mut(n);
    invert(x, d, rest);
    let (m, rest) = rest.split_at_mut(n);
    mulhi(m, a, x, rest);
    q[..n].copy_from_slice(a);
    q[n] = add_assign(&mut q[..n], m);
}

pub(crate) fn div_approx_newton_scratch_len(n: usize) -> usize {
    n + invert_scratch_len(n).max(n + mulhi_scratch_len(n))
}

pub(crate) fn div_approx_newton_error(n: usize) -> usize {
    2 + mulhi_error(n)
}

/// `y = Y` with `Y <= floor(sqrt(B^(3n) / A)) <= Y + 1`, for the `n`-limb
/// `A = a` normalised for square roots (`A >= B^n / 4`: one of the top two
/// bits set) and `y.len() == n + 1` (`B^n <= Y <= 2 B^n`). In other words
/// `Y / B^n` is `1 / sqrt(A / B^n)` from below, within 2 ulps. Needs
/// `rsqrt_approx_scratch_len(n)`.
///
/// Below `RSQRT_NEWTON_THRESHOLD` limbs, `Y = floor(sqrt(floor(B^(3n) / A)))`
/// exactly (the inner floor doesn't change the root, as squares are
/// integers), through `Natural`, which allocates. Above, with `h = n / 2 + 1`
/// and `l = n - h`, let `Yh` be the result for `A`'s top `h` limbs `Ah`, and
/// write `a = A / B^n`, `y = Yh / B^h = (1 - ε) / sqrt(a)`. One Newton step
/// `y + y (1 - a y^2) / 2` is, exactly, `(1 - ε^2 (3 - ε) / 2) / sqrt(a)`, a
/// lower bound. Its computation:
///
/// - `R = B^(n + 2h) - A Yh^2` exactly, so `1 - a y^2 = R / B^(n + 2h)`, and
///   the step is `Y = Yh B^l + Yh R / (2 B^(3h))`. `R` comes from the full
///   product, or from a wrapped one (`rsqrt_residual_wrapped`) where the
///   cyclic transform applies.
/// - `R` is truncated to `floor(R / B^k)` with `k = 2h - 1`, and the
///   correction to `floor(Yh floor(R / B^k) / (2 B^(h + 1)))`: floors also when
///   `R < 0`, so `Y` stays below the exact step. As `Yh <= 2 B^h`, the first
///   truncation loses at most `1 / B`, the second below 1.
/// - `ε` is small: `Yh` is within 2 ulps of `1 / sqrt(Ah / B^h)` from below by
///   recursion, which is within a factor `1 + 2 B^-h` of `1 / sqrt(a)` from
///   above (`Ah B^l <= A < (Ah + 1) B^l` and `Ah >= B^h / 4`). So
///   `|ε| <= 2 B^-h`, and `|R| <= (2|ε| + ε^2) B^(n + 2h) < B^(n + h + 1)`.
///   The step itself is then within `12.4 B^(n - 2h) <= 12.4 / B` of
///   `B^n / sqrt(a)` (`2h >= n + 1`).
///
/// In all `sqrt(B^(3n) / A) - Y < 1 + 14 / B`, so `floor` of it is at most
/// `Y + 1`. A result below `B^n` (when `A` is near `B^n`) is raised to `B^n`,
/// which the true value exceeds.
pub fn rsqrt_approx(y: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    assert!(n >= 1 && y.len() == n + 1, "rsqrt_approx: length mismatch");
    assert!(
        a[n - 1] >> (Limb::BITS - 2) != 0,
        "rsqrt_approx: not normalised for square roots"
    );
    assert!(
        scratch.len() >= rsqrt_approx_scratch_len(n),
        "rsqrt_approx: scratch too short"
    );
    rsqrt_unchecked(y, a, scratch);
}

/// `(h, l)` of a Newton step on `n` limbs.
fn rsqrt_lens(n: usize) -> (usize, usize) {
    let h = n / 2 + 1;
    (h, n - h)
}

/// Scratch limbs [`rsqrt_approx`] needs for `n` limbs.
#[must_use]
pub fn rsqrt_approx_scratch_len(n: usize) -> usize {
    if n < RSQRT_NEWTON_THRESHOLD {
        return 0;
    }
    let (h, l) = rsqrt_lens(n);
    let residual = match rsqrt_wrap_lens(n) {
        Some((m, w)) => m + mul_wrap_scratch_len(m, n, 2 * h + 2).max(2 * w + mullo_scratch_len(w)),
        None => mul_scratch_len(n, 2 * h + 2),
    };
    let step = (2 * h + 2)
        + (n + 2 * h + 2)
        + sqr_scratch_len(h + 1)
            .max(residual)
            .max((l + 2) + (n + 3) + (l + 2) + mul_scratch_len(h + 1, l + 2));
    h + 1 + rsqrt_approx_scratch_len(h).max(step)
}

/// The error bound of [`rsqrt_approx`]: 0 for the exact small case, else 1.
#[must_use]
pub fn rsqrt_approx_error(n: usize) -> usize {
    usize::from(n >= RSQRT_NEWTON_THRESHOLD)
}

fn rsqrt_unchecked(y: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    if n < RSQRT_NEWTON_THRESHOLD {
        rsqrt_exact(y, a);
        return;
    }
    let (h, l) = rsqrt_lens(n);
    let (yh, rest) = scratch.split_at_mut(h + 1);
    rsqrt_unchecked(yh, &a[l..], rest);
    let (t, rest) = rest.split_at_mut(2 * h + 2);
    let (p, rest) = rest.split_at_mut(n + 2 * h + 2);
    sqr(t, yh, rest);
    let negative = match rsqrt_wrap_lens(n) {
        Some((m, w)) => rsqrt_residual_wrapped(p, a, t, h, m, w, rest),
        None => rsqrt_residual_full(p, a, t, h, rest),
    };
    debug_assert!(
        p[n + h + 1..].iter().all(|&x| x == 0),
        "rsqrt_approx: residual too large"
    );

    // floor(R / B^k), as a magnitude rounded up when R < 0.
    let k = 2 * h - 1;
    let (mt, rest) = rest.split_at_mut(l + 2);
    mt.copy_from_slice(&p[k..=n + h]);
    if negative && p[..k].iter().any(|&x| x != 0) {
        let carry = add_1_assign(mt, 1);
        debug_assert_eq!(carry, 0);
    }
    // The correction Yh floor(R / B^k) / (2 B^(h + 1)), floored.
    let (u, rest) = rest.split_at_mut(n + 3);
    let (c, rest) = rest.split_at_mut(l + 2);
    mul(u, yh, mt, rest);
    let out = shr(c, &u[h + 1..], 1);
    if negative && (out != 0 || u[..=h].iter().any(|&x| x != 0)) {
        let carry = add_1_assign(c, 1);
        debug_assert_eq!(carry, 0);
    }

    y[..l].fill(0);
    y[l..].copy_from_slice(yh);
    if negative {
        let borrow = sub_assign(y, c);
        debug_assert_eq!(borrow, 0);
    } else {
        let carry = add_assign(y, c);
        debug_assert_eq!(carry, 0);
    }
    if y[n] == 0 {
        y.fill(0);
        y[n] = 1;
    }
    debug_assert!(y[n] <= 2, "rsqrt_approx: above 2 B^n");
}

/// `(m, w)` of the wrapped residual for `n` limbs, when the cyclic transform
/// applies: `m` the largest transform length (`2^k` or `3 2^k`) up to the
/// span `n + h + 2`, and `w = n + h + 3 - m` low limbs (at least 1, at most
/// `m`, `n` and `2h + 2`). Consecutive transform lengths are at most 1.5
/// apart, so `w` is at most a third of the span, plus one.
fn rsqrt_wrap_lens(n: usize) -> Option<(usize, usize)> {
    if n < RSQRT_NEWTON_THRESHOLD {
        return None;
    }
    let (h, _) = rsqrt_lens(n);
    let span = n + h + 2;
    let pow2 = 1usize << (usize::BITS - 1 - span.leading_zeros());
    let three = 3 * (pow2 / 2);
    let m = if pow2 >= 2 && three <= span {
        three
    } else {
        pow2
    };
    let w = span + 1 - m;
    (wraps_cyclically(m) && w <= n.min(2 * h + 2)).then_some((m, w))
}

/// `|R|` for `R = B^(n + 2h) - A T` into `p` (all of it, `n + 2h + 2` limbs),
/// from the full product; whether `R <= 0`.
fn rsqrt_residual_full(
    p: &mut [Limb],
    a: &[Limb],
    t: &[Limb],
    h: usize,
    scratch: &mut [Limb],
) -> bool {
    let n = a.len();
    mul(p, a, t, scratch);
    let top = n + 2 * h;
    let negative = p[top..].iter().any(|&x| x != 0);
    if negative {
        let borrow = sub_1_assign(&mut p[top..], 1);
        debug_assert_eq!(borrow, 0);
    } else {
        // B^top - P = (B^top - 1 - P) + 1, with P > 0.
        for x in &mut p[..top] {
            *x = !*x;
        }
        let carry = add_1_assign(&mut p[..=top], 1);
        debug_assert_eq!(carry, 0);
    }
    negative
}

/// [`rsqrt_residual_full`] from `A T mod (B^m - 1)` (a cyclic transform of
/// `m` coefficients, where the full product has `n + 2h + 2`; any `m >= w`
/// gives the same result, the transform only applies at its lengths) and
/// `A T mod B^w` (a short product of `w` limbs), for `(m, w)` from
/// [`rsqrt_wrap_lens`].
///
/// `|R| < B^(n + h + 1)` (see [`rsqrt_approx`]). With `e = R mod (B^m - 1)`
/// in `[0, B^m - 1)`, `R = e + j (B^m - 1)` for an integer `j` with
/// `|j| < B^(n + h + 1 - m) + 1 <= 2 B^(w - 2) + 1` (`j` in `{-1, 0}` at
/// `w = 1`). Modulo `B^w` (`m >= w`), `R = e - j`, and `R = -A T` as
/// `w <= n + 2h`: so `j = e + A T mod B^w`, taken in `(-B^w / 2, B^w / 2)`.
/// Then `R = e - j + j B^m` in two's complement on `m + w = n + h + 3`
/// limbs, which holds `R` with room for its sign.
fn rsqrt_residual_wrapped(
    p: &mut [Limb],
    a: &[Limb],
    t: &[Limb],
    h: usize,
    m: usize,
    w: usize,
    scratch: &mut [Limb],
) -> bool {
    let n = a.len();
    let (prod, rest) = scratch.split_at_mut(m);
    mul_wrap(prod, a, t, rest);
    // e = (B^((n + 2h) mod m) - W) mod (B^m - 1).
    p.fill(0);
    p[(n + 2 * h) % m] = 1;
    if sub_n_assign(&mut p[..m], prod) != 0 {
        // The borrow added B^m; the residue needs B^m - 1.
        sub_1_assign(&mut p[..m], 1);
    }
    // j = e + A T mod B^w, its sign from the top bit.
    let (j, rest) = rest.split_at_mut(w);
    let (lo, rest) = rest.split_at_mut(w);
    mullo(lo, &a[..w], &t[..w], rest);
    j.copy_from_slice(&p[..w]);
    add_n_assign(j, lo);
    let j_negative = j[w - 1] >> (Limb::BITS - 1) == 1;
    if j_negative {
        for x in j.iter_mut() {
            *x = !*x;
        }
        add_1_assign(j, 1);
    }
    // R = e - j + j B^m modulo B^(m + w), carries out of the top dropped.
    let r = &mut p[..m + w];
    if j_negative {
        add_assign(r, j);
        sub_n_assign(&mut r[m..], j);
    } else {
        sub_assign(r, j);
        add_n_assign(&mut r[m..], j);
    }
    let negative = r[m + w - 1] >> (Limb::BITS - 1) == 1;
    if negative {
        for x in r.iter_mut() {
            *x = !*x;
        }
        add_1_assign(r, 1);
    }
    // R = 0 counts as R <= 0, as in the full product.
    negative || r.iter().all(|&x| x == 0)
}

/// `floor(sqrt(floor(B^(3n) / A)))`, through `Natural`.
fn rsqrt_exact(y: &mut [Limb], a: &[Limb]) {
    let n = a.len() as u64;
    let num = Natural::from_limbs(&[1]) << (3 * n * u64::from(Limb::BITS));
    let (q, _) = num.div_rem(&Natural::from_limbs(a));
    let s = q.isqrt();
    y.fill(0);
    y[..s.limbs().len()].copy_from_slice(s.limbs());
}

/// `s = S` with `S <= floor(sqrt(A B^n)) <= S + sqrt_approx_error(n)`, for the
/// `n`-limb `A = a` normalised for square roots (`A >= B^n / 4`) and
/// `s.len() == n`. Needs `sqrt_approx_scratch_len(n)`.
///
/// Below `SQRT_NEWTON_THRESHOLD` limbs the root is exact (through `Natural`,
/// which allocates); above, see `sqrt_approx_newton`.
pub fn sqrt_approx(s: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    assert!(n >= 1 && s.len() == n, "sqrt_approx: length mismatch");
    assert!(
        a[n - 1] >> (Limb::BITS - 2) != 0,
        "sqrt_approx: not normalised for square roots"
    );
    assert!(
        scratch.len() >= sqrt_approx_scratch_len(n),
        "sqrt_approx: scratch too short"
    );
    if sqrt_by_newton(n) {
        sqrt_approx_newton(s, a, scratch);
        return;
    }
    let shift = n as u64 * u64::from(Limb::BITS);
    let r = (Natural::from_limbs(a) << shift).isqrt();
    s.fill(0);
    s[..r.limbs().len()].copy_from_slice(r.limbs());
}

/// Scratch limbs [`sqrt_approx`] needs for `n` limbs.
#[must_use]
pub fn sqrt_approx_scratch_len(n: usize) -> usize {
    if sqrt_by_newton(n) {
        sqrt_approx_newton_scratch_len(n)
    } else {
        0
    }
}

/// The error bound of [`sqrt_approx`] for `n` limbs.
#[must_use]
pub fn sqrt_approx_error(n: usize) -> usize {
    if sqrt_by_newton(n) {
        sqrt_approx_newton_error(n)
    } else {
        0
    }
}

/// [`sqrt_approx`] through the inverse square root: `sqrt(A B^n) = A Y* /
/// B^n` for `Y* = sqrt(B^(3n) / A)`, so with `Y = B^n + J` from
/// [`rsqrt_approx`], `S = A + mulhi(A, J)` (or `2A` when `Y = 2 B^n`).
///
/// `A (Y* - Y) / B^n < 2 A / B^n < 2`, and `mulhi` loses below
/// `mulhi_error(n) + 1`, so `sqrt(A B^n) - S < 3 + mulhi_error(n)`, from
/// above.
pub(crate) fn sqrt_approx_newton(s: &mut [Limb], a: &[Limb], scratch: &mut [Limb]) {
    let n = a.len();
    let (y, rest) = scratch.split_at_mut(n + 1);
    rsqrt_unchecked(y, a, rest);
    if y[n] == 2 {
        let out = shl(s, a, 1);
        debug_assert_eq!(out, 0, "sqrt_approx: Y = 2 B^n needs A = B^n / 4");
        return;
    }
    let (m, rest) = rest.split_at_mut(n);
    mulhi(m, a, &y[..n], rest);
    s.copy_from_slice(a);
    let carry = add_assign(s, m);
    debug_assert_eq!(carry, 0, "sqrt_approx: above the true root");
}

pub(crate) fn sqrt_approx_newton_scratch_len(n: usize) -> usize {
    n + 1 + rsqrt_approx_scratch_len(n).max(n + mulhi_scratch_len(n))
}

pub(crate) fn sqrt_approx_newton_error(n: usize) -> usize {
    2 + mulhi_error(n)
}

#[cfg(test)]
mod tests;
