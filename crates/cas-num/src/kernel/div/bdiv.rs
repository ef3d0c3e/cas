//! Hensel (2-adic) division, from the low end (GMP's `bdiv` family): each
//! quotient limb comes from the lowest remaining limb times `d^-1 mod B`,
//! with no estimate to correct.
//!
//! [`bdiv_q`] is the 2-adic quotient `Q = A D^-1 mod B^qn` of any `A` by an
//! odd `D`, which needs only the low `qn` limbs of both; when `D` divides `A`
//! and `qn = an - dn + 1`, it is the true quotient, as `A / D < B^qn`.
//! Below [`BDIV_Q_DC_THRESHOLD`] divisor limbs it is schoolbook; up to
//! [`BDIV_MU_THRESHOLD`], divide and conquer, whose halves are Hensel
//! divisions with remainder and a short product; from there, it works in
//! blocks, each a short product by the 2-adic inverse of `D`'s low limbs
//! ([`binvert`], Newton's iteration) and a product to update the limbs above.
//!
//! [`divexact`] also proves that the division is exact. Hensel division with
//! remainder (schoolbook, then divide and conquer from
//! [`BDIV_QR_DC_THRESHOLD`]) leaves the remainder above `Q` for free; from
//! [`DIVEXACT_DIV_THRESHOLD`] on, division with remainder is cheaper.

use super::super::mul::mullo::{cost as mullo_cost, mullo, mullo_scratch_len};
use super::super::mul::{
    addmul_2_carries, limb_inverse, mul_add, mul_wrap, mul_wrap_scratch_len, wrap_len,
    wraps_cyclically,
};
use super::{div_rem, div_rem_scratch_len};
use crate::Limb;
use crate::kernel::{
    add_1_assign, add_assign, addmul_1, mul, mul_scratch_len, ntt, sub_1_assign, sub_assign, sub_n,
};

/// Smallest divisor length at which [`bdiv_q`] divides by divide and conquer
/// instead of schoolbook (when below [`BDIV_MU_THRESHOLD`]). From the tuning
/// configuration, as the others.
pub(crate) const BDIV_Q_DC_THRESHOLD: usize = crate::tuning::BDIV_Q_DC_THRESHOLD;

/// Smallest length at which Hensel division with remainder ([`divexact`]'s,
/// and the halves of [`bdiv_q`]'s divide and conquer) splits in two instead
/// of schoolbook.
pub(crate) const BDIV_QR_DC_THRESHOLD: usize = crate::tuning::BDIV_QR_DC_THRESHOLD;

/// Smallest divisor length at which [`bdiv_q`] works in blocks with the
/// 2-adic inverse.
pub(crate) const BDIV_MU_THRESHOLD: usize = crate::tuning::BDIV_MU_THRESHOLD;

/// Smallest length at which [`binvert`] takes a Newton step instead of a
/// Hensel division of 1 (at least 2, which the step needs).
pub(crate) const BINV_NEWTON_THRESHOLD: usize = crate::tuning::BINV_NEWTON_THRESHOLD;

/// Smallest divisor length at which [`divexact`] divides with remainder
/// instead of by Hensel division.
pub(crate) const DIVEXACT_DIV_THRESHOLD: usize = crate::tuning::DIVEXACT_DIV_THRESHOLD;

/// How [`bdiv_q`] divides by a divisor of `dn` limbs (after cutting it to
/// the quotient's length): the one place that orders the thresholds, so the
/// scratch length and the division agree whatever their values.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum BdivQ {
    Schoolbook,
    DivideAndConquer,
    Blocks,
}

fn bdiv_q_algorithm(dn: usize) -> BdivQ {
    if dn >= BDIV_MU_THRESHOLD {
        BdivQ::Blocks
    } else if dn >= BDIV_Q_DC_THRESHOLD {
        BdivQ::DivideAndConquer
    } else {
        BdivQ::Schoolbook
    }
}

/// Scratch limbs [`bdiv_q`] needs for `qn` quotient limbs and `dn` divisor
/// limbs.
#[must_use]
pub fn bdiv_q_scratch_len(qn: usize, dn: usize) -> usize {
    let dn = dn.min(qn);
    match bdiv_q_algorithm(dn) {
        BdivQ::Schoolbook => 0,
        BdivQ::DivideAndConquer => dc_bdiv_q_scratch_len(qn, dn),
        BdivQ::Blocks => qn + mu_bdiv_q_scratch_len(qn, dn),
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
    if qn == 0 {
        return;
    }
    let d = &d[..d.len().min(qn)];
    match bdiv_q_algorithm(d.len()) {
        BdivQ::Schoolbook => {
            q.copy_from_slice(&a[..qn]);
            sb_bdiv_q(q, d);
            negate(q);
        }
        BdivQ::DivideAndConquer => {
            q.copy_from_slice(&a[..qn]);
            dc_bdiv_q(q, d, scratch);
            negate(q);
        }
        BdivQ::Blocks => {
            let (t, rest) = scratch.split_at_mut(qn);
            t.copy_from_slice(&a[..qn]);
            mu_bdiv_q(q, t, d, rest);
        }
    }
}

/// Scratch limbs [`divexact`] needs for `an` and `dn` limbs.
#[must_use]
pub fn divexact_scratch_len(an: usize, dn: usize) -> usize {
    if dn < DIVEXACT_DIV_THRESHOLD {
        an + 1 + bdiv_qr_scratch_len(an - dn + 1, dn)
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
        let qn = q.len();
        let (t, rest) = scratch.split_at_mut(an + 1);
        t[..an].copy_from_slice(a);
        t[an] = 0;
        let c = bdiv_qr(t, d, qn, rest);
        q.copy_from_slice(&t[..qn]);
        negate(q);
        exact_remainder(&t[qn..], c, d)
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

/// The quotient limbs `(q0, q1)` of Hensel division that clear `t0` and the
/// limb `t1` above it when `q0 d + q1 d B` is added, for `ninv = -d^-1 mod B`
/// and `d`'s low limbs `d0` and `d1`. `q0`'s row leaves `t1 + hi(q0 d0) +
/// lo(q0 d1)` at limb 1, plus the carry out of limb 0, which is one unless
/// `t0` (and so `q0`) is zero.
#[inline]
fn hensel_pair(t0: Limb, t1: Limb, d0: Limb, d1: Limb, ninv: Limb) -> (Limb, Limb) {
    let q0 = t0.wrapping_mul(ninv);
    let (_, hi) = mul_add(q0, d0, 0);
    let u = t1
        .wrapping_add(hi)
        .wrapping_add(q0.wrapping_mul(d1))
        .wrapping_add(Limb::from(t0 != 0));
    (q0, u.wrapping_mul(ninv))
}

/// Adds the two limbs `(c0, c1)` and the pending carry `cy` at `t[k]` and
/// `t[k + 1]`, dropping what falls at or above `t.len()`, and returns the
/// carry out of `t[k + 1]`, pending at `t[k + 2]`. `t[k] + c0 + cy` and
/// `t[k + 1] + c1` plus a carry each carry at most one.
fn add_pair(t: &mut [Limb], k: usize, (c0, c1): (Limb, Limb), cy: bool) -> bool {
    let Some(x) = t.get_mut(k) else {
        return false;
    };
    let (s, a) = x.overflowing_add(c0);
    let (s, b) = s.overflowing_add(Limb::from(cy));
    *x = s;
    let Some(x) = t.get_mut(k + 1) else {
        return false;
    };
    let (s, a2) = x.overflowing_add(c1);
    let (s, b2) = s.overflowing_add(Limb::from(a | b));
    *x = s;
    a2 | b2
}

/// The rows of [`sb_bdiv_q`] and [`sb_bdiv_qr`] two at a time, for the first
/// `rows` (even) quotient limbs and `dn = d.len() >= 2`: one pass of
/// [`addmul_2_carries`] adds both rows, so each limb of `t` is loaded and
/// stored once per two quotient limbs, and the two limbs it clears take the
/// quotient. Rows are cut at `t.len()`. The carry of the pass is pending one
/// limb above it in `cy`, as for single rows: at `i + dn` before the pair at
/// `i`, which `dn >= 2` keeps above the two limbs the pair's quotient reads.
/// Returns `rows`.
fn sb_bdiv_pairs(t: &mut [Limb], d: &[Limb], ninv: Limb, rows: usize, cy: &mut bool) -> usize {
    debug_assert!(rows & 1 == 0 && rows <= t.len() && d.len() >= 2);
    let (tn, dn) = (t.len(), d.len());
    for i in (0..rows).step_by(2) {
        let (q0, q1) = hensel_pair(t[i], t[i + 1], d[0], d[1], ninv);
        let len = dn.min(tn - i);
        let c = addmul_2_carries(&mut t[i..i + len], &d[..len], q0, q1, 0);
        debug_assert!(t[i] == 0 && t[i + 1] == 0);
        (t[i], t[i + 1]) = (q0, q1);
        *cy = add_pair(t, i + len, c, *cy);
    }
    rows
}

/// The negated 2-adic quotient `-T d^-1 mod B^qn` by schoolbook, in place on
/// `t` of `qn >= d.len()` limbs, which holds `T` on entry: as
/// [`sb_bdiv_qr`], with each row cut at `B^qn`.
fn sb_bdiv_q(t: &mut [Limb], d: &[Limb]) {
    let (qn, dn) = (t.len(), d.len());
    let ninv = limb_inverse(d[0]).wrapping_neg();
    let mut cy = false;
    let start = if dn < 2 {
        0
    } else {
        sb_bdiv_pairs(t, d, ninv, qn & !1, &mut cy)
    };
    for i in start..qn {
        let qi = t[i].wrapping_mul(ninv);
        let len = dn.min(qn - i);
        let hi = addmul_1(&mut t[i..i + len], &d[..len], qi);
        debug_assert_eq!(t[i], 0);
        t[i] = qi;
        if i + len < qn {
            let (s, c) = hi.overflowing_add(Limb::from(cy));
            let (s, c2) = t[i + len].overflowing_add(s);
            t[i + len] = s;
            cy = c | c2;
        }
    }
}

/// Schoolbook Hensel division with remainder, in place on `t` of `qn + dn`
/// limbs (`dn = d.len()`), which holds `T` on entry: leaves the 2-adic
/// quotient negated, `Q' = -T / D mod B^qn`, in `t[..qn]` and `S` in
/// `t[qn..]`, and returns `c` with `T + Q' D = B^qn (S + c B^dn)` (GMP's
/// `sbpi1_bdiv_qr`).
///
/// Adds `q'_i D B^i` with `q'_i = t_i (-d^-1) mod B`, which clears limb `i`
/// for `q'_i` to take (GMP's sign: the carries are added one limb up, at
/// most one bit at a time, instead of a borrow rippling); two rows at a time
/// ([`sb_bdiv_pairs`]) but for a one-limb `d` or the last of an odd count.
fn sb_bdiv_qr(t: &mut [Limb], d: &[Limb], qn: usize) -> bool {
    let dn = d.len();
    debug_assert_eq!(t.len(), qn + dn);
    let ninv = limb_inverse(d[0]).wrapping_neg();
    let mut cy = false;
    let start = if dn < 2 {
        0
    } else {
        sb_bdiv_pairs(t, d, ninv, qn & !1, &mut cy)
    };
    for i in start..qn {
        let qi = t[i].wrapping_mul(ninv);
        let hi = addmul_1(&mut t[i..i + dn], d, qi);
        debug_assert_eq!(t[i], 0);
        t[i] = qi;
        let (s, c) = hi.overflowing_add(Limb::from(cy));
        let (s, c2) = t[i + dn].overflowing_add(s);
        t[i + dn] = s;
        cy = c | c2;
    }
    cy
}

/// Whether `D` divides `A`, from the remainder of a Hensel division of `A`
/// by `D` to `qn` limbs with `A < B^qn D` ([`sb_bdiv_qr`]'s `S` and carry
/// `c`). If `D` divides `A`, `Q' = B^qn - A / D` (as `A / D < B^qn`) and `S =
/// D`, or `A = Q' = S = 0`; conversely `S = D` gives `A = (B^qn - Q') D`, and
/// `S = 0` gives `A = Q' = 0`.
fn exact_remainder(s: &[Limb], c: bool, d: &[Limb]) -> bool {
    !c && (s == d || s.iter().all(|&x| x == 0))
}

/// Scratch limbs [`bdiv_qr_n`] needs for `n` limbs.
fn bdiv_qr_n_scratch_len(n: usize) -> usize {
    if n < BDIV_QR_DC_THRESHOLD {
        return 0;
    }
    let (lo, hi) = (n / 2, n - n / 2);
    (n + mul_scratch_len(hi, lo)).max(bdiv_qr_n_scratch_len(hi))
}

/// [`sb_bdiv_qr`] for `qn = dn = n` (`t` of `2n` limbs), by divide and
/// conquer from [`BDIV_QR_DC_THRESHOLD`] (GMP's `dc_bdiv_qr_n`): the low
/// `lo = floor(n / 2)` quotient limbs from `t`'s low `2 lo` limbs and `d`'s
/// low `lo`, then the rest of `d` times them added at `B^lo`; the same for
/// the high `hi = n - lo` limbs, from `t[lo..n + hi]` and `d`'s low `hi`,
/// whose product with `d`'s high `lo` lands at `B^n`. Each half's carry is
/// added with its product, which it cannot overflow, as the product is at
/// most `(B^lo - 1) (B^hi - 1) < B^n - B^hi`.
fn bdiv_qr_n(t: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) -> bool {
    let n = d.len();
    debug_assert_eq!(t.len(), 2 * n);
    if n < BDIV_QR_DC_THRESHOLD {
        return sb_bdiv_qr(t, d, n);
    }
    let (lo, hi) = (n / 2, n - n / 2);
    let c = bdiv_qr_n(&mut t[..2 * lo], &d[..lo], scratch);
    let (p, rest) = scratch.split_at_mut(n);
    mul(p, &d[lo..], &t[..lo], rest);
    let out = add_1_assign(&mut p[lo..], Limb::from(c));
    debug_assert_eq!(out, 0);
    let mut top = add_assign(&mut t[lo..], p);
    let c = bdiv_qr_n(&mut t[lo..n + hi], &d[..hi], scratch);
    let (p, rest) = scratch.split_at_mut(n);
    mul(p, &t[lo..n], &d[hi..], rest);
    let out = add_1_assign(&mut p[hi..], Limb::from(c));
    debug_assert_eq!(out, 0);
    top += add_assign(&mut t[n..], p);
    // T + Q' D < 2 B^2n, so at most one carry leaves the top in all.
    debug_assert!(top <= 1);
    top != 0
}

/// The first quotient block of [`bdiv_qr`] and [`dc_bdiv_q`]: `qn` itself up
/// to `dn`, else the part of `qn` above a multiple of `dn`, in `1..=dn`
/// (GMP does the typically smaller block first).
fn first_block(qn: usize, dn: usize) -> usize {
    if qn > dn { (qn - 1) % dn + 1 } else { qn }
}

/// Scratch limbs [`bdiv_qr`] needs for `qn` quotient and `dn` divisor limbs.
fn bdiv_qr_scratch_len(qn: usize, dn: usize) -> usize {
    if dn < BDIV_QR_DC_THRESHOLD {
        return 0;
    }
    let k = first_block(qn, dn);
    let update = if k < dn {
        dn + mul_scratch_len(k, dn - k)
    } else {
        0
    };
    let blocks = if qn > k { bdiv_qr_n_scratch_len(dn) } else { 0 };
    update.max(bdiv_qr_n_scratch_len(k)).max(blocks)
}

/// [`sb_bdiv_qr`] for any `qn` and `dn` (`t` of `qn + dn` limbs), from
/// [`BDIV_QR_DC_THRESHOLD`] divisor limbs in blocks of [`bdiv_qr_n`]
/// (GMP's `dc_bdiv_qr`): a first block of `k` ([`first_block`]) with `d`'s
/// low `k` limbs, whose product with the rest of `d` (with the block's carry)
/// is added at `B^k`, then blocks of `dn`, each block's carry pending at the
/// limb above the next block's dividend.
fn bdiv_qr(t: &mut [Limb], d: &[Limb], qn: usize, scratch: &mut [Limb]) -> bool {
    let dn = d.len();
    debug_assert_eq!(t.len(), qn + dn);
    if dn < BDIV_QR_DC_THRESHOLD {
        return sb_bdiv_qr(t, d, qn);
    }
    let k = first_block(qn, dn);
    let mut c = bdiv_qr_n(&mut t[..2 * k], &d[..k], scratch);
    let mut top = 0;
    if k < dn {
        let (p, rest) = scratch.split_at_mut(dn);
        mul(p, &t[..k], &d[k..], rest);
        let out = add_1_assign(&mut p[k..], Limb::from(c));
        debug_assert_eq!(out, 0);
        top = add_assign(&mut t[k..], p);
        c = false;
    }
    for i in (k..qn).step_by(dn) {
        top += add_1_assign(&mut t[i + dn..], Limb::from(c));
        c = bdiv_qr_n(&mut t[i..i + 2 * dn], d, scratch);
    }
    // As for `bdiv_qr_n`, at most one carry leaves the top.
    debug_assert!(top + Limb::from(c) <= 1);
    top != 0 || c
}

/// Scratch limbs [`bdiv_q_n`] needs for `n` limbs.
fn bdiv_q_n_scratch_len(n: usize) -> usize {
    if n < BDIV_Q_DC_THRESHOLD {
        return 0;
    }
    let (lo, hi) = (n / 2, n - n / 2);
    bdiv_qr_n_scratch_len(lo)
        .max(lo + mullo_scratch_len(lo))
        .max(bdiv_q_n_scratch_len(hi))
}

/// [`sb_bdiv_q`] for `qn = dn = n`, in place, by divide and conquer from
/// [`BDIV_Q_DC_THRESHOLD`] (GMP's `dc_bdiv_q_n`): the low `lo = floor(n /
/// 2)` limbs by [`bdiv_qr_n`], which clears `t`'s low `lo` limbs and leaves
/// its remainder above them; then `d`'s other `hi = n - lo` limbs times them,
/// at `B^lo` and cut at `B^n`, a short product (and for odd `n` one row more,
/// as `hi = lo + 1`); then the high `hi` limbs from `t[lo..]` and `d`'s low
/// `hi` limbs, the same way.
fn bdiv_q_n(t: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = d.len();
    debug_assert_eq!(t.len(), n);
    if n < BDIV_Q_DC_THRESHOLD {
        sb_bdiv_q(t, d);
        return;
    }
    let (lo, hi) = (n / 2, n - n / 2);
    let c = bdiv_qr_n(&mut t[..2 * lo], &d[..lo], scratch);
    let (p, rest) = scratch.split_at_mut(lo);
    let (q, above) = t.split_at_mut(lo);
    mullo(p, q, &d[hi..], rest);
    add_assign(&mut above[hi - lo..], p);
    if lo < hi {
        // d[lo] times the low quotient at B^lo, its carry and the remainder's
        // at B^(2 lo) = B^(n - 1).
        let out = addmul_1(&mut above[..lo], q, d[lo]);
        above[lo] = above[lo].wrapping_add(out).wrapping_add(Limb::from(c));
    }
    bdiv_q_n(above, &d[..hi], scratch);
}

/// Scratch limbs [`dc_bdiv_q`] needs for `qn >= dn` limbs.
fn dc_bdiv_q_scratch_len(qn: usize, dn: usize) -> usize {
    if qn == dn {
        return bdiv_q_n_scratch_len(dn);
    }
    let k = first_block(qn, dn);
    let update = if k < dn {
        dn + mul_scratch_len(k, dn - k)
    } else {
        0
    };
    let blocks = if qn - k > dn {
        bdiv_qr_n_scratch_len(dn)
    } else {
        0
    };
    update
        .max(bdiv_qr_n_scratch_len(k))
        .max(blocks)
        .max(bdiv_q_n_scratch_len(dn))
}

/// [`sb_bdiv_q`] for `qn = t.len() >= dn = d.len()`, in place, by divide and
/// conquer (GMP's `dc_bdiv_q`): as [`bdiv_qr`] for all but the last `dn`
/// quotient limbs, cut at `B^qn`, then [`bdiv_q_n`] for those.
fn dc_bdiv_q(t: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (qn, dn) = (t.len(), d.len());
    debug_assert!(qn >= dn);
    let mut i = 0;
    if qn > dn {
        let k = first_block(qn, dn);
        let mut c = bdiv_qr_n(&mut t[..2 * k], &d[..k], scratch);
        if k < dn {
            let (p, rest) = scratch.split_at_mut(dn);
            mul(p, &t[..k], &d[k..], rest);
            let out = add_1_assign(&mut p[k..], Limb::from(c));
            debug_assert_eq!(out, 0);
            add_assign(&mut t[k..], p);
            c = false;
        }
        i = k;
        while qn - i > dn {
            add_1_assign(&mut t[i + dn..], Limb::from(c));
            c = bdiv_qr_n(&mut t[i..i + 2 * dn], d, scratch);
            i += dn;
        }
        // The last carry is pending at B^qn, beyond the cut.
    }
    bdiv_q_n(&mut t[i..], d, scratch);
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

/// The modulus length `m` at which [`mul_known_low`] wraps `x y` for `n` limbs
/// of `x` and `y` together, `k` of them known, if it does.
fn known_low_wrap(n: usize, k: usize) -> Option<usize> {
    let m = wrap_len(n - k);
    (m < n && wraps_cyclically(m)).then_some(m)
}

/// Scratch limbs [`mul_known_low`] needs for `xn` and `yn` limbs, `k` known.
fn mul_known_low_scratch_len(xn: usize, yn: usize, k: usize) -> usize {
    known_low_wrap(xn + yn, k).map_or_else(
        || mul_scratch_len(xn, yn),
        |m| mul_wrap_scratch_len(m, xn, yn),
    )
}

/// The limbs of `x y` from `k` on into `p[k..]` (`p.len() = x.len() +
/// y.len()`, the low `k` limbs unspecified), given those low `k = low.len()
/// <= min(x.len(), y.len())` limbs: by a cyclic transform modulo `B^m - 1`
/// with `m = wrap_len(p.len() - k)` when that wraps (as GMP's `mu_bdiv_q`
/// and `binvert` use `mulmod_bnm1`, see [`unwrap_known_low`]), else the full
/// product.
fn mul_known_low(p: &mut [Limb], x: &[Limb], y: &[Limb], low: &[Limb], scratch: &mut [Limb]) {
    let (n, k) = (p.len(), low.len());
    debug_assert!(n == x.len() + y.len() && k <= x.len().min(y.len()));
    let Some(m) = known_low_wrap(n, k) else {
        mul(p, x, y, scratch);
        debug_assert_eq!(&p[..k], low);
        return;
    };
    mul_wrap(&mut p[..m], x, y, scratch);
    unwrap_known_low(p, m, low);
}

/// The limbs of a product `P` of `p.len()` limbs from `k` on, into `p[k..]`,
/// from its residue modulo `B^m - 1` in `p[..m]` (for a non-zero multiple
/// either one, `0` or `B^m - 1`; for zero, `0`) and its low `k = low.len()`
/// limbs, with `p.len() - m <= k` and `P`'s factors at least as long.
///
/// With `P = L + B^m H`, `L < B^m` and `H < B^w`, `w = p.len() - m <= k`,
/// the residue `r` is `L + H` or, when that reaches `B^m - 1`, that less
/// `B^m - 1`. Subtracting `L`'s known low `w` limbs from `r`'s gives `H`, and
/// the borrow, taken from limb `w` up, gives `L`'s other limbs: in the second
/// case `r < B^w` (`L`'s limbs from `w` on are all ones), and the borrow
/// turns `r`'s zero limbs back into ones and takes the one back from `H`.
/// `H = B^w - 1` would need `P >= B^(m + w) - B^w`, which is above the
/// largest product `(B^xn - 1) (B^yn - 1)` as `w <= min(xn, yn)`.
fn unwrap_known_low(p: &mut [Limb], m: usize, low: &[Limb]) {
    let (w, k) = (p.len() - m, low.len());
    let (r, h) = p.split_at_mut(m);
    let borrow = sub_n(h, &r[..w], &low[..w]);
    let out = sub_1_assign(&mut p[w..], borrow);
    debug_assert!(out == 0 && p[w..k] == low[w..]);
}

/// The transform length at which [`mu_bdiv_q`] prepares the inverse for its
/// blocks' short products, if it does: `wrap_len(2 b0)`, which holds each
/// whole product, when one forward transform and then two per block cost
/// less than a short product per block in the NTT's cost model.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
pub(crate) fn bdiv_prepared_inverse(qn: usize, dn: usize) -> Option<usize> {
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let len = wrap_len(2 * b0);
    let (uses, fresh) = blocks.fold((0_usize, 0.0), |(uses, cost), (_, b, _)| {
        (uses + 1, cost + mullo_cost(b))
    });
    let prepared = ntt::cyclic_cost(len) as f64 * (1 + 2 * uses) as f64 / 3.0;
    (uses >= 2 && wraps_cyclically(len) && prepared < fresh).then_some(len)
}

/// The modulus length at which [`mu_bdiv_q`] prepares the divisor's
/// transform for its updates, if it does: `m = wrap_len(dn)`, when two
/// blocks or more of `b0` limbs update by a product with all of `d` that
/// wraps there. One forward transform then saves one per update after the
/// first.
pub(crate) fn bdiv_prepared_divisor(qn: usize, dn: usize) -> Option<usize> {
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let m = known_low_wrap(b0 + dn, b0)?;
    let updates = blocks
        .filter(|&(i, b, dl)| i + b < qn && b == b0 && dl == dn)
        .count();
    (updates >= 2).then_some(m)
}

/// Scratch limbs [`mu_bdiv_q`] needs for `qn` quotient limbs and `dn <= qn`
/// divisor limbs.
fn mu_bdiv_q_scratch_len(qn: usize, dn: usize) -> usize {
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let (prep_i, prep_d) = (bdiv_prepared_inverse(qn, dn), bdiv_prepared_divisor(qn, dn));
    let cyclic = |len: Option<usize>| len.map_or(0, ntt::cyclic_prepared_scratch_len);
    let work = blocks
        .map(|(_, b, dl)| {
            let short = match prep_i {
                Some(len) => len + ntt::cyclic_prepared_scratch_len(len),
                None => mullo_scratch_len(b),
            };
            let update = match prep_d {
                Some(m) if b == b0 && dl == dn => ntt::cyclic_prepared_scratch_len(m),
                _ => mul_known_low_scratch_len(dl, b, b),
            };
            short.max(update)
        })
        .max()
        .unwrap_or(0)
        .max(cyclic(prep_i))
        .max(cyclic(prep_d));
    let prepared = prep_i.map_or(0, ntt::prepared_len) + prep_d.map_or(0, ntt::prepared_len);
    b0 + (b0 + dn) + (prepared + work).max(binvert_scratch_len(b0))
}

/// `q = t d^-1 mod B^qn` for `qn = q.len() = t.len() >= dn = d.len()`, in
/// blocks: each block of `b` quotient limbs is `t`'s next `b` limbs times `d^-1
/// mod B^b`, a short product, and its product with `d`, whose low `b` limbs
/// are that block of `t` ([`mul_known_low`]), is subtracted from the limbs of
/// `t` above it, modulo `B^qn`. Clobbers `t`.
///
/// With the inverse's transform prepared ([`bdiv_prepared_inverse`]), each
/// short product is the low `b` limbs of a whole product by all `b0` limbs
/// of the inverse, whose limbs from `b` on only reach the product's from `b`
/// on. With the divisor's ([`bdiv_prepared_divisor`]), the updates of full
/// blocks by all of `d` take its transform.
fn mu_bdiv_q(q: &mut [Limb], t: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let (qn, dn) = (q.len(), d.len());
    let (b0, blocks) = bdiv_blocks(qn, dn);
    let (inv, rest) = scratch.split_at_mut(b0);
    let (p, rest) = rest.split_at_mut(b0 + dn);
    binvert(inv, d, rest);
    let (prep_i, prep_d) = (bdiv_prepared_inverse(qn, dn), bdiv_prepared_divisor(qn, dn));
    let (inv_t, rest) = rest.split_at_mut(prep_i.map_or(0, ntt::prepared_len));
    let (d_t, rest) = rest.split_at_mut(prep_d.map_or(0, ntt::prepared_len));
    if let Some(len) = prep_i {
        ntt::prepare_cyclic(inv_t, len, inv, rest);
    }
    if let Some(m) = prep_d {
        ntt::prepare_cyclic(d_t, m, d, rest);
    }
    for (i, b, dl) in blocks {
        let hi = i + b;
        if let Some(len) = prep_i {
            let (c, rest) = rest.split_at_mut(len);
            ntt::mul_cyclic_prepared(c, &t[i..hi], inv_t, rest);
            debug_assert!(c[b + b0..].iter().all(|&x| x == 0));
            q[i..hi].copy_from_slice(&c[..b]);
        } else {
            mullo(&mut q[i..hi], &t[i..hi], &inv[..b], rest);
        }
        if hi == qn {
            break;
        }
        let p = &mut p[..b + dl];
        match prep_d {
            Some(m) if b == b0 && dl == dn => {
                ntt::mul_cyclic_prepared(&mut p[..m], &q[i..hi], d_t, rest);
                unwrap_known_low(p, m, &t[i..hi]);
            }
            _ => mul_known_low(p, &d[..dl], &q[i..hi], &t[i..hi], rest),
        }
        let len = (qn - hi).min(dl);
        sub_assign(&mut t[hi..], &p[b..b + len]);
    }
}

/// The modulus length at which [`binvert`]'s Newton step to `n` limbs
/// prepares `X`'s transform for both its products, if it does: when `d X`
/// wraps at `m`, and the five transforms there cost less than that
/// product's three and the short product of `n - h` limbs in the NTT's cost
/// model.
#[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
pub(crate) fn binvert_shared(n: usize) -> Option<usize> {
    let h = n - n / 2;
    let m = known_low_wrap(n + h, h)?;
    let cyclic = ntt::cyclic_cost(m) as f64;
    (2.0 * cyclic < 3.0 * mullo_cost(n - h)).then_some(m)
}

/// Scratch limbs [`binvert`] needs for `n` limbs.
pub(crate) fn binvert_scratch_len(n: usize) -> usize {
    if n < BINV_NEWTON_THRESHOLD {
        return bdiv_q_n_scratch_len(n);
    }
    let h = n - n / 2;
    let products = match binvert_shared(n) {
        Some(m) => ntt::prepared_len(m) + m + ntt::cyclic_prepared_scratch_len(m),
        None => mul_known_low_scratch_len(n, h, h).max(mullo_scratch_len(n - h)),
    };
    let step = (n + h) + (n - h) + h + products;
    binvert_scratch_len(h).max(step)
}

/// `x = d^-1 mod B^n` for `n = x.len() <= d.len()` and an odd `d`, with
/// `scratch.len() >= binvert_scratch_len(n)`.
///
/// A Newton step from `X = d^-1 mod B^h`, `h = ceil(n / 2)`: `d X = 1 + E B^h`
/// (a product whose low `h` limbs are known, [`mul_known_low`]), and `X (2 -
/// d X) = X - X E B^h` is the inverse modulo `B^2h`, so the new limbs are `-(X
/// E) mod B^(n - h)`, a short product. Where it pays ([`binvert_shared`]),
/// both products share `X`'s transform: `X E` has `n <= m` limbs, so it does
/// not wrap, and its low `n - h` limbs are the short product's.
pub(crate) fn binvert(x: &mut [Limb], d: &[Limb], scratch: &mut [Limb]) {
    let n = x.len();
    debug_assert!(n <= d.len() && d[0] & 1 == 1);
    if n < BINV_NEWTON_THRESHOLD {
        // Hensel division of 1 by d, modulo B^n.
        if n == 0 {
            return;
        }
        x[0] = 1;
        x[1..].fill(0);
        bdiv_q_n(x, &d[..n], scratch);
        negate(x);
        return;
    }
    let h = n - n / 2;
    binvert(&mut x[..h], d, scratch);
    let (p, rest) = scratch.split_at_mut(n + h);
    let (c, rest) = rest.split_at_mut(n - h);
    let (one, rest) = rest.split_at_mut(h);
    one[0] = 1;
    one[1..].fill(0);
    let (x_lo, x_hi) = x.split_at_mut(h);
    if let Some(m) = binvert_shared(n) {
        let (prep, rest) = rest.split_at_mut(ntt::prepared_len(m));
        ntt::prepare_cyclic(prep, m, x_lo, rest);
        ntt::mul_cyclic_prepared(&mut p[..m], &d[..n], prep, rest);
        unwrap_known_low(p, m, one);
        let (t, rest) = rest.split_at_mut(m);
        ntt::mul_cyclic_prepared(t, &p[h..n], prep, rest);
        debug_assert!(t[n..].iter().all(|&v| v == 0));
        c.copy_from_slice(&t[..n - h]);
    } else {
        mul_known_low(p, &d[..n], x_lo, one, rest);
        mullo(c, &x_lo[..n - h], &p[h..n], rest);
    }
    for (y, &v) in x_hi.iter_mut().zip(c.iter()) {
        *y = !v;
    }
    add_1_assign(x_hi, 1);
}
