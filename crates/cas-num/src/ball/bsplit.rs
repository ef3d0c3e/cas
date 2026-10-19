//! Binary splitting for hypergeometric-type series (Haible and Papanikolaou
//! \[HP98\]; \[BZ10\] §4.9; CAS-DOC-5 §7.1).
//!
//! A [`Series`] is
//!
//! ```text
//! S = Σ_{k < n} a(k) / b(k) · p(0) p(1) … p(k) / (q(0) q(1) … q(k))
//! ```
//!
//! with integer `p, q, a, b`. Over a range `[l, r)` of terms, with the
//! products starting at `l`, binary splitting keeps `P = Π p`, `Q = Π q`,
//! `B = Π b` and `T = B Q S(l, r)`, all integers. A leaf is `P = p(k)`,
//! `Q = q(k)`, `B = b(k)`, `T = a(k) p(k)`, and two halves merge as
//!
//! ```text
//! P = Pl Pr,  Q = Ql Qr,  B = Bl Br,  T = Br Qr Tl + Bl Pl Tr,
//! ```
//!
//! since `S(l, r) = S(l, m) + (Pl / Ql) S(m, r)`. The sum is `T / (B Q)`, one
//! big division; everything before it is products of growing size, so the
//! cost is our multiplication's. The tail beyond `n` terms is each series'
//! business: callers bound it and add it to the ball's radius.
//!
//! A power of two in `q`, as in the bit-burst series of `x = p / 2^b` (CAS-95,
//! a trick of ours), is kept apart: `q(k) 2^s(k)` with `s` from
//! [`Series::shift`], and a range carries `E = Σ s` instead of multiplying it
//! into `Q`. Then `S = T / (B Q 2^E)` and the merge becomes
//!
//! ```text
//! T = 2^Er Br Qr Tl + Bl Pl Tr,   E = El + Er,
//! ```
//!
//! a shift where there was a product, and `Q` and every product with it
//! stay as short as the odd part.

use super::Ball;
use crate::Integer;

#[cfg(test)]
mod tests;

/// The integer sequences of a series; `None` stands for 1 and must be the
/// same choice for every `k`, so a factor that is always 1 costs nothing.
pub(crate) trait Series {
    fn p(&self, k: u64) -> Option<Integer>;
    fn q(&self, k: u64) -> Integer;
    fn a(&self, k: u64) -> Option<Integer>;
    fn b(&self, k: u64) -> Option<Integer>;

    /// The power of two in `q(k)`'s place, kept out of the products: the
    /// term's denominator is `q(k) 2^shift(k)`.
    fn shift(&self, _k: u64) -> u64 {
        0
    }
}

/// `P, Q, B, T` of a range; `None` is 1 (or, for `p`, not computed: the
/// rightmost ranges never need it).
pub(super) struct Split {
    p: Option<Integer>,
    q: Integer,
    b: Option<Integer>,
    t: Integer,
    e: u64,
}

/// `x 2^k`.
fn shl(x: Integer, k: u64) -> Integer {
    if k == 0 {
        return x;
    }
    let (neg, mag) = x.into_parts();
    Integer::from_parts(neg, mag << k)
}

/// `x · y`, with `None` as 1.
fn mul_opt(x: Integer, y: Option<&Integer>) -> Integer {
    match y {
        Some(y) => &x * y,
        None => x,
    }
}

/// `x · y` of two optional factors.
fn mul_opts(x: Option<Integer>, y: Option<&Integer>) -> Option<Integer> {
    match (x, y) {
        (Some(x), y) => Some(mul_opt(x, y)),
        (None, y) => y.cloned(),
    }
}

/// Terms `[lo, hi)`; `P` only when `need_p`.
pub(super) fn split<S: Series>(s: &S, lo: u64, hi: u64, need_p: bool) -> Split {
    debug_assert!(lo < hi);
    if hi - lo == 1 {
        let p = s.p(lo);
        let t = match (s.a(lo), &p) {
            (Some(a), p) => mul_opt(a, p.as_ref()),
            (None, Some(p)) => p.clone(),
            (None, None) => Integer::ONE,
        };
        return Split {
            p,
            q: s.q(lo),
            b: s.b(lo),
            t,
            e: s.shift(lo),
        };
    }
    let mid = lo + (hi - lo) / 2;
    let l = split(s, lo, mid, true);
    let r = split(s, mid, hi, need_p);
    // T = 2^Er Br Qr Tl + Bl Pl Tr.
    let mut t = shl(mul_opt(&r.q * &l.t, r.b.as_ref()), r.e);
    t += &mul_opt(mul_opt(r.t, l.p.as_ref()), l.b.as_ref());
    Split {
        p: if need_p {
            mul_opts(l.p, r.p.as_ref())
        } else {
            None
        },
        q: &l.q * &r.q,
        b: mul_opts(l.b, r.b.as_ref()),
        t,
        e: l.e + r.e,
    }
}

/// The first `n` terms as `(T, B Q, Q)`: `S = T / (B Q)`, and `Q` alone,
/// which tail bounds are often stated in. `n` must be at least 1; the
/// series has no [`Series::shift`] (see [`sum_shifted`]).
pub(crate) fn sum<S: Series>(s: &S, n: u64) -> (Integer, Integer, Integer) {
    assert!(n >= 1, "bsplit::sum: no terms");
    let Split { q, b, t, e, .. } = split(s, 0, n, false);
    debug_assert_eq!(e, 0, "bsplit::sum of a shifted series");
    let d = mul_opt(q.clone(), b.as_ref());
    (t, d, q)
}

/// The first `n` terms as `(T, B Q, E)`: `S = T / (B Q 2^E)`. `n` must be
/// at least 1.
pub(crate) fn sum_shifted<S: Series>(s: &S, n: u64) -> (Integer, Integer, u64) {
    assert!(n >= 1, "bsplit::sum_shifted: no terms");
    let Split { q, b, t, e, .. } = split(s, 0, n, false);
    (t, mul_opt(q, b.as_ref()), e)
}

/// `t / d` as a ball at `prec` bits; both are rounded to `prec` bits first,
/// so that the division costs `prec`-bit operands however long they are.
pub(crate) fn quotient(t: &Integer, d: &Integer, prec: u64) -> Ball {
    let t = Ball::from_integer(t).round_to_accuracy(prec);
    let d = Ball::from_integer(d).round_to_accuracy(prec);
    t.div(&d, prec)
}
