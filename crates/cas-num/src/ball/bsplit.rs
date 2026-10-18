//! Binary splitting for hypergeometric-type series (CAS-DOC-5 §7.1).
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
}

/// `P, Q, B, T` of a range; `None` is 1 (or, for `p`, not computed: the
/// rightmost ranges never need it).
pub(super) struct Split {
    p: Option<Integer>,
    q: Integer,
    b: Option<Integer>,
    t: Integer,
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
        };
    }
    let mid = lo + (hi - lo) / 2;
    let l = split(s, lo, mid, true);
    let r = split(s, mid, hi, need_p);
    // T = Br Qr Tl + Bl Pl Tr.
    let mut t = mul_opt(&r.q * &l.t, r.b.as_ref());
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
    }
}

/// The first `n` terms as `(T, B Q, Q)`: `S = T / (B Q)`, and `Q` alone,
/// which tail bounds are often stated in. `n` must be at least 1.
pub(crate) fn sum<S: Series>(s: &S, n: u64) -> (Integer, Integer, Integer) {
    assert!(n >= 1, "bsplit::sum: no terms");
    let Split { q, b, t, .. } = split(s, 0, n, false);
    let d = mul_opt(q.clone(), b.as_ref());
    (t, d, q)
}

/// `t / d` as a ball at `prec` bits; both are rounded to `prec` bits first,
/// so that the division costs `prec`-bit operands however long they are.
pub(crate) fn quotient(t: &Integer, d: &Integer, prec: u64) -> Ball {
    let t = Ball::from_integer(t).round_to_accuracy(prec);
    let d = Ball::from_integer(d).round_to_accuracy(prec);
    t.div(&d, prec)
}
