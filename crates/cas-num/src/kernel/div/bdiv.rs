//! Hensel (2-adic) division, from the low end (GMP's `bdiv` family): each
//! quotient limb comes from the lowest remaining limb times `d^-1 mod B`,
//! with no estimate to correct.
//!
//! [`bdiv_q`] is the 2-adic quotient `Q = A D^-1 mod B^qn` of any `A` by an
//! odd `D`, which needs only the low `qn` limbs of both; when `D` divides `A`
//! and `qn = an - dn + 1`, it is the true quotient, as `A / D < B^qn`.
//! Below [`BDIV_MU_THRESHOLD`] divisor limbs it is schoolbook; from it, it
//! works in blocks, each a short product by the 2-adic inverse of `D`'s low
//! limbs ([`binvert`], Newton's iteration) and a product to update the limbs
//! above.
//!
//! [`divexact`] also proves that the division is exact. Schoolbook Hensel
//! division leaves the remainder above `Q` for free; from
//! [`DIVEXACT_DIV_THRESHOLD`] on, where the check would cost a product more,
//! division with remainder is cheaper.

use super::super::mul::limb_inverse;
use super::super::mul::mullo::{mullo, mullo_scratch_len};
use super::{div_rem, div_rem_scratch_len};
use crate::Limb;
use crate::kernel::{add_1_assign, addmul_1, mul, mul_scratch_len, sub_assign, submul_1};

/// Smallest divisor length at which [`bdiv_q`] works in blocks with the
/// 2-adic inverse instead of limb by limb. From the tuning configuration.
pub(crate) const BDIV_MU_THRESHOLD: usize = crate::tuning::BDIV_MU_THRESHOLD;

/// Smallest length at which [`binvert`] takes a Newton step instead of a
/// schoolbook Hensel division of 1 (at least 2, which the step needs).
pub(crate) const BINV_NEWTON_THRESHOLD: usize = crate::tuning::BINV_NEWTON_THRESHOLD;

/// Smallest divisor length at which [`divexact`] divides with remainder
/// instead of by schoolbook Hensel division.
pub(crate) const DIVEXACT_DIV_THRESHOLD: usize = crate::tuning::DIVEXACT_DIV_THRESHOLD;

/// Scratch limbs [`bdiv_q`] needs for `qn` quotient limbs and `dn` divisor
/// limbs.
#[must_use]
pub fn bdiv_q_scratch_len(qn: usize, dn: usize) -> usize {
    let dn = dn.min(qn);
    if dn < BDIV_MU_THRESHOLD {
        qn
    } else {
        qn + mu_bdiv_q_scratch_len(qn, dn)
    }
}

/// `q = a d^-1 mod B^qn` for `qn = q.len() <= a.len()` and an odd `d`, with
/// `scratch.len() >= bdiv_q_scratch_len(qn, d.len())`. Uses only the low `qn`
/// limbs of `a` and `d`. When `d` divides `a` and `qn = a.len() - d.len() +
/// 1` (with `d`'s top limb non-zero), `q = a / d`.
pub fn bdiv_q(q: &mut [Limb], a: &[Limb], d: &[Limb], scratch: &mut [Limb]) {
    let qn = q.len();
    assert!(qn <= a.len(), "bdiv_q: quotient longer than dividend");
    assert!(
        d.first().is_some_and(|&d0| d0 & 1 == 1),
        "bdiv_q: even divisor"
    );
    assert!(
        scratch.len() >= bdiv_q_scratch_len(qn, d.len()),
        "bdiv_q: scratch too short"
    );
    let d = &d[..d.len().min(qn)];
    let (t, rest) = scratch.split_at_mut(qn);
    t.copy_from_slice(&a[..qn]);
    if d.len() < BDIV_MU_THRESHOLD {
        sb_bdiv_q(q, t, d);
    } else {
        mu_bdiv_q(q, t, d, rest);
    }
}

/// Scratch limbs [`divexact`] needs for `an` and `dn` limbs.
#[must_use]
pub fn divexact_scratch_len(an: usize, dn: usize) -> usize {
    if dn < DIVEXACT_DIV_THRESHOLD {
        an
    } else {
        dn + div_rem_scratch_len(an, dn)
    }
}

/// `q = a / d` for an odd `d` with a non-zero top limb, `a.len() >= d.len()`
/// and `q.len() == a.len() - d.len() + 1`, with `scratch.len() >=
/// divexact_scratch_len(a.len(), d.len())`. Returns whether `d` divides `a`;
/// when it does not, `q` holds an unspecified value.
pub fn divexact(q: &mut [Limb], a: &[Limb], d: &[Limb], scratch: &mut [Limb]) -> bool {
    let (an, dn) = (a.len(), d.len());
    assert!(
        dn >= 1 && an >= dn,
        "divexact: divisor longer than dividend"
    );
    assert_eq!(q.len(), an - dn + 1, "divexact: quotient length");
    assert!(d[0] & 1 == 1, "divexact: even divisor");
    assert!(d[dn - 1] != 0, "divexact: divisor not normalised");
    assert!(
        scratch.len() >= divexact_scratch_len(an, dn),
        "divexact: scratch too short"
    );
    if dn < DIVEXACT_DIV_THRESHOLD {
        let t = &mut scratch[..an];
        t.copy_from_slice(a);
        sb_bdiv_qr(q, t, d)
    } else {
        let (r, rest) = scratch.split_at_mut(dn);
        div_rem(q, r, a, d, rest);
        r.iter().all(|&x| x == 0)
    }
}

/// `-q` modulo `B^q.len()`.
fn negate(q: &mut [Limb]) {
    for x in q.iter_mut() {
        *x = !*x;
    }
    add_1_assign(q, 1);
}

/// [`bdiv_q`] by schoolbook on `t = a mod B^qn` (clobbered), `qn = q.len() =
/// t.len() >= d.len()`: as [`sb_bdiv_qr`], with each row cut at `B^qn`.
fn sb_bdiv_q(q: &mut [Limb], t: &mut [Limb], d: &[Limb]) {
    let (qn, dn) = (q.len(), d.len());
    let ninv = limb_inverse(d[0]).wrapping_neg();
    let mut cy = false;
    for i in 0..qn {
        let qi = t[i].wrapping_mul(ninv);
        q[i] = qi;
        let len = dn.min(qn - i);
        let hi = addmul_1(&mut t[i..i + len], &d[..len], qi);
        debug_assert_eq!(t[i], 0);
        if i + len < qn {
            let (s, c) = hi.overflowing_add(Limb::from(cy));
            let (s, c2) = t[i + len].overflowing_add(s);
            t[i + len] = s;
            cy = c | c2;
        }
    }
    negate(q);
}

/// [`divexact`] by schoolbook Hensel division, on `t = A` (clobbered).
///
/// Adds `q'_i D B^i` with `q'_i = t_i (-d^-1) mod B`, which clears limb `i`
/// (GMP's sign: the carries are added one limb up, at most one bit at a time,
/// instead of a borrow rippling). Then `A + Q' D = B^qn S` with `S` the
/// limbs left above `qn` and the top carry, and `Q' = -A / D mod B^qn`. If `D`
/// divides `A`, `Q' = B^qn - Q` and `S = D`, or `A = Q' = S = 0`; conversely
/// `S = D` gives `A = (B^qn - Q') D`, and `S = 0` gives `A = Q' = 0`.
fn sb_bdiv_qr(q: &mut [Limb], t: &mut [Limb], d: &[Limb]) -> bool {
    let (qn, dn) = (q.len(), d.len());
    let ninv = limb_inverse(d[0]).wrapping_neg();
    let mut cy = false;
    let mut top = (0, false);
    for i in 0..qn {
        let qi = t[i].wrapping_mul(ninv);
        q[i] = qi;
        let hi = addmul_1(&mut t[i..i + dn], d, qi);
        debug_assert_eq!(t[i], 0);
        let (s, c) = hi.overflowing_add(Limb::from(cy));
        if i + 1 < qn {
            let (s, c2) = t[i + dn].overflowing_add(s);
            t[i + dn] = s;
            cy = c | c2;
        } else {
            top = (s, c);
        }
    }
    // S = t[qn..] + top B^(dn - 1), with top = (limb, carry out of it).
    let s = &t[qn..];
    let exact = if top.1 {
        false
    } else if top.0 == 0 && s.iter().all(|&x| x == 0) {
        true
    } else {
        top.0 == d[dn - 1] && s == &d[..dn - 1]
    };
    negate(q);
    exact
}

/// The quotient block length of [`mu_bdiv_q`] (GMP's choice): blocks of at
/// most `dn` limbs, equal but for the last, when the quotient is longer;
/// else two halves, whose inverse costs less than one of the whole.
fn bdiv_block_len(qn: usize, dn: usize) -> usize {
    if qn > dn {
        qn.div_ceil(qn.div_ceil(dn))
    } else {
        qn.div_ceil(2)
    }
}

/// `(block start, block length, divisor limbs used)` for each block of
/// [`mu_bdiv_q`], and the block length.
fn bdiv_blocks(qn: usize, dn: usize) -> (usize, impl Iterator<Item = (usize, usize, usize)>) {
    let b0 = bdiv_block_len(qn, dn);
    let blocks = (0..qn)
        .step_by(b0)
        .map(move |i| (i, b0.min(qn - i), dn.min(qn - i)));
    (b0, blocks)
}

/// Scratch limbs [`mu_bdiv_q`] needs for `qn` quotient limbs and `dn <= qn`
/// divisor limbs.
fn mu_bdiv_q_scratch_len(qn: usize, dn: usize) -> usize {
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let work = blocks
        .map(|(_, b, dl)| mullo_scratch_len(b).max(mul_scratch_len(dl, b)))
        .max()
        .unwrap_or(0);
    b0 + (b0 + dn) + work.max(binvert_scratch_len(b0))
}

/// `q = t d^-1 mod B^qn` for `qn = q.len() = t.len() >= dn = d.len()`, in
/// blocks: each block of `b` quotient limbs is `t`'s next `b` limbs times `d^-1
/// mod B^b`, a short product, and its product with `d` is subtracted from the
/// limbs of `t` above it, modulo `B^qn`. Clobbers `t`.
fn mu_bdiv_q(q: &mut [Limb], t: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (qn, dn) = (q.len(), d.len());
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let (inv, rest) = scratch.split_at_mut(b0);
    let (p, rest) = rest.split_at_mut(b0 + dn);
    binvert(inv, d, rest);
    for (i, b, dl) in blocks {
        let hi = i + b;
        mullo(&mut q[i..hi], &t[i..hi], &inv[..b], rest);
        if hi == qn {
            break;
        }
        let p = &mut p[..b + dl];
        mul(p, &d[..dl], &q[i..hi], rest);
        debug_assert_eq!(&p[..b], &t[i..hi]);
        let len = (qn - hi).min(dl);
        sub_assign(&mut t[hi..], &p[b..b + len]);
    }
}

/// Scratch limbs [`binvert`] needs for `n` limbs.
pub(crate) fn binvert_scratch_len(n: usize) -> usize {
    if n < BINV_NEWTON_THRESHOLD {
        return n;
    }
    let h = n - n / 2;
    let step = n + h + mul_scratch_len(n, h).max((n - h) + mullo_scratch_len(n - h));
    binvert_scratch_len(h).max(step)
}

/// `x = d^-1 mod B^n` for `n = x.len() <= d.len()` and an odd `d`, with
/// `scratch.len() >= binvert_scratch_len(n)`.
///
/// A Newton step from `X = d^-1 mod B^h`, `h = ceil(n / 2)`: `d X = 1 + E B^h`,
/// and `X (2 - d X) = X - X E B^h` is the inverse modulo `B^2h`, so the new
/// limbs are `-(X E) mod B^(n - h)`, a short product.
pub(crate) fn binvert(x: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = x.len();
    debug_assert!(n <= d.len() && d[0] & 1 == 1);
    if n < BINV_NEWTON_THRESHOLD {
        // Hensel division of 1 by d, modulo B^n.
        let inv = limb_inverse(d[0]);
        let t = &mut scratch[..n];
        t.fill(0);
        if let Some(t0) = t.first_mut() {
            *t0 = 1;
        }
        for i in 0..n {
            let xi = t[i].wrapping_mul(inv);
            x[i] = xi;
            submul_1(&mut t[i..], &d[..n - i], xi);
        }
        return;
    }
    let h = n - n / 2;
    binvert(&mut x[..h], d, scratch);
    let (p, rest) = scratch.split_at_mut(n + h);
    mul(p, &d[..n], &x[..h], rest);
    debug_assert!(p[0] == 1 && p[1..h].iter().all(|&l| l == 0));
    let (c, rest) = rest.split_at_mut(n - h);
    let (x_lo, x_hi) = x.split_at_mut(h);
    mullo(c, &x_lo[..n - h], &p[h..n], rest);
    for (y, &v) in x_hi.iter_mut().zip(c.iter()) {
        *y = !v;
    }
    add_1_assign(x_hi, 1);
}
