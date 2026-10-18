//! Greatest common divisor and least common multiple, by Lehmer's algorithm.
//!
//! Lehmer's algorithm (Knuth, TAOCP vol. 2, §4.5.2, Algorithm L) runs Euclid's
//! algorithm on the leading 128 bits of the operands for as long as the
//! steps provably leave the full values non-negative (a condition in the
//! manner of Jebelean's and of GMP's `mpn_hgcd2`), collecting them in a 2×2
//! matrix, then applies the matrix to the full values with one pass of
//! single-limb multiplications. When the leading bits cannot decide a single
//! step, a full division step is taken instead. Once both values fit a limb
//! the binary algorithm (or a word-level extended Euclid) finishes. Large
//! operands are reduced by half-GCD steps first (`hgcd`).
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithm uses Knuth's names (a, b, u, v, q) and matrix entries"
)]

use std::cmp::Ordering;

use super::Natural;
use crate::{Integer, Limb, kernel};

/// Bits of the leading parts the inner loop works on: two limbs, so that
/// each matrix advances about 64 bits with entries that fit a limb.
const LEAD_BITS: u64 = 128;

/// The matrix `((a, b), (c, d))` from Lehmer's inner loop: the pair `(x, y)`
/// becomes `(a x + b y, c x + d y)`. It is a product of Euclidean steps, so
/// it is unimodular and the entries of each row have opposite signs (or one
/// is zero); the loop takes only steps that leave both results non-negative,
/// though not necessarily in order. Every entry has a magnitude below
/// `2^64`.
#[derive(Clone, Copy, Debug)]
struct Matrix {
    a: i128,
    b: i128,
    c: i128,
    d: i128,
}

/// One reduction of the pair `(x, y)`, reported to cofactor tracking.
enum Step<'a> {
    /// `(x, y)` became `M (x, y)`.
    Matrix(&'a Matrix),
    /// `(x, y)` became `M^-1 (x, y)`.
    Hgcd(&'a HgcdMatrix),
    /// `(x, y)` became `(y, x - q y)`, with `q` the given limbs.
    Quotient(&'a [Limb]),
}

mod hgcd;
#[cfg(test)]
mod tests;

use hgcd::HgcdMatrix;

impl Natural {
    /// The greatest common divisor; `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, other: &Natural) -> Natural {
        self.gcd_from(other, hgcd::GCD_DC_THRESHOLD)
    }

    /// [`Natural::gcd`] with half-GCD steps from `dc` limbs.
    fn gcd_from(&self, other: &Natural, dc: usize) -> Natural {
        let (mut a, mut b) = ordered(self, other);
        reduce(&mut a, &mut b, dc, false, |_| {});
        match *b {
            [] => Natural::from_vec(a),
            [d] => Natural::from(kernel::gcd_1(kernel::rem_1(&a, d), d)),
            _ => unreachable!("reduce leaves at most one limb"),
        }
    }

    /// The least common multiple; zero if either operand is zero.
    #[must_use]
    pub fn lcm(&self, other: &Natural) -> Natural {
        if self.is_zero() || other.is_zero() {
            return Natural::ZERO;
        }
        self.div_exact_unchecked(&self.gcd(other)) * other
    }

    /// `(g, s)` with `g = gcd(self, other)` and `g = s self + t other` for
    /// some integer `t`. `s` is not reduced; [`Integer::gcd_ext`] normalises
    /// the cofactors.
    pub(crate) fn gcd_cofactor(&self, other: &Natural) -> (Natural, Integer) {
        self.gcd_cofactor_from(other, hgcd::GCD_DC_THRESHOLD)
    }

    /// [`Natural::gcd_cofactor`] with half-GCD steps from `dc` limbs.
    fn gcd_cofactor_from(&self, other: &Natural, dc: usize) -> (Natural, Integer) {
        let (mut a, mut b) = ordered(self, other);
        // sa and sb are the coefficients of `self` in a and b.
        let (mut sa, mut sb) = if self >= other {
            (Integer::ONE, Integer::ZERO)
        } else {
            (Integer::ZERO, Integer::ONE)
        };
        reduce(&mut a, &mut b, dc, true, |step| {
            track(&mut sa, &mut sb, &step);
        });
        if b.is_empty() {
            return (Natural::from_vec(a), sa);
        }
        if a.len() > 1 {
            let q = div_step(&mut a, &mut b);
            track(&mut sa, &mut sb, &Step::Quotient(&q));
        }
        // Both fit a limb: finish with a word-level extended Euclid, whose
        // cofactors are bounded by the operands and so fit an i128.
        let (mut x, mut y) = (single(&a), single(&b));
        let (mut p, mut q) = ((1i128, 0i128), (0i128, 1i128));
        while y != 0 {
            let k = i128::from(x / y);
            (x, y) = (y, x % y);
            let next = (p.0 - k * q.0, p.1 - k * q.1);
            (p, q) = (q, next);
        }
        let s = &sa * &Integer::from(p.0) + &sb * &Integer::from(p.1);
        (Natural::from(x), s)
    }
}

/// The limbs of the larger and the smaller operand.
fn ordered(x: &Natural, y: &Natural) -> (Vec<Limb>, Vec<Limb>) {
    let (a, b) = if x >= y { (x, y) } else { (y, x) };
    (a.limbs().to_vec(), b.limbs().to_vec())
}

/// The value of at most one limb.
fn single(x: &[Limb]) -> Limb {
    match *x {
        [] => 0,
        [x] => x,
        _ => unreachable!("more than one limb"),
    }
}

/// Applies `step` to the coefficients `(sa, sb)` of the pair it reduced.
fn track(sa: &mut Integer, sb: &mut Integer, step: &Step<'_>) {
    match step {
        Step::Matrix(m) => {
            let (a, b) = (std::mem::take(sa), std::mem::take(sb));
            *sa = &a * &Integer::from(m.a) + &b * &Integer::from(m.b);
            *sb = &a * &Integer::from(m.c) + &b * &Integer::from(m.d);
        }
        Step::Hgcd(m) => {
            // M^-1 = ±[[u11, -u01], [-u10, u00]].
            let [[u00, u01], [u10, u11]] =
                m.m.each_ref().map(|row| row.each_ref().map(Integer::from));
            let (a, b) = (std::mem::take(sa), std::mem::take(sb));
            let (x, y) = (&u11 * &a - &u01 * &b, &u00 * &b - &u10 * &a);
            (*sa, *sb) = if m.negative { (-x, -y) } else { (x, y) };
        }
        Step::Quotient(q) => {
            let next = &*sa - &(Integer::from(Natural::from_limbs(q)) * &*sb);
            *sa = std::mem::replace(sb, next);
        }
    }
}

/// Reduces normalised `a >= b` until `b` has at most one limb, keeping
/// `gcd(a, b)` and `a >= b`, by half-GCD steps while `b` has at least `dc`
/// limbs and Lehmer steps below. With `track`, reports every step to
/// `on_step`; without, half-GCD steps are not reported, which saves building
/// their matrices.
fn reduce(
    a: &mut Vec<Limb>,
    b: &mut Vec<Limb>,
    dc: usize,
    track: bool,
    mut on_step: impl FnMut(Step<'_>),
) {
    let mut buf = [Vec::new(), Vec::new()];
    while b.len() >= 2 {
        if b.len() >= dc {
            let progress = if track {
                hgcd::hgcd(a, b).map(|mut m| {
                    if kernel::cmp(a, b) == Ordering::Less {
                        std::mem::swap(a, b);
                        m.swap_columns();
                    }
                    on_step(Step::Hgcd(&m));
                })
            } else {
                hgcd::hgcd_values(a, b).then(|| {
                    if kernel::cmp(a, b) == Ordering::Less {
                        std::mem::swap(a, b);
                    }
                })
            };
            if progress.is_none() {
                let q = div_step(a, b);
                on_step(Step::Quotient(&q));
            }
        } else if let Some(mut m) = lehmer_matrix(a, b, None) {
            apply(&m, a, b, &mut buf);
            if kernel::cmp(a, b) == Ordering::Less {
                std::mem::swap(a, b);
                m = Matrix {
                    a: m.c,
                    b: m.d,
                    c: m.a,
                    d: m.b,
                };
            }
            on_step(Step::Matrix(&m));
        } else {
            let q = div_step(a, b);
            on_step(Step::Quotient(&q));
        }
        debug_assert!(kernel::is_normalized(a) && kernel::is_normalized(b));
        // Equal only when a half-GCD stopped at a gcd of at least B^s.
        debug_assert_ne!(kernel::cmp(a, b), Ordering::Less);
    }
}

/// `(a, b) = (b, a mod b)` for `a >= b > 0`; returns the quotient.
fn div_step(a: &mut Vec<Limb>, b: &mut Vec<Limb>) -> Vec<Limb> {
    let mut q = vec![0; a.len() - b.len() + 1];
    let mut r = vec![0; b.len()];
    let mut scratch = vec![0; kernel::div_rem_scratch_len(a.len(), b.len())];
    kernel::div_rem(&mut q, &mut r, a, b, &mut scratch);
    r.truncate(kernel::normalized_len(&r));
    q.truncate(kernel::normalized_len(&q));
    *a = std::mem::replace(b, r);
    q
}

/// Bits `k .. k + 128` of `x`.
fn leading(x: &[Limb], k: u64) -> u128 {
    let i = usize::try_from(k / 64).expect("bit index fits usize");
    let s = k % 64;
    let limb = |j: usize| u128::from(x.get(j).copied().unwrap_or(0));
    let v = limb(i) | limb(i + 1) << 64;
    if s == 0 {
        v
    } else {
        v >> s | limb(i + 2) << (128 - s)
    }
}

/// Lehmer's inner loop on the leading 128 bits of `a >= b`, where `a` has at
/// least two limbs: a double-limb step in the manner of GMP's `mpn_hgcd2`.
/// With a `floor` of `t`, only steps whose remainder is provably at least
/// `2^t` are taken. `None` when not even the first step is.
fn lehmer_matrix(a: &[Limb], b: &[Limb], floor: Option<u64>) -> Option<Matrix> {
    let k = kernel::bit_length(a).saturating_sub(LEAD_BITS);
    // With a = 2^k u + e_a, 0 <= e_a < 2^k (and alike for b), a remainder of
    // the leading parts r = s u + t v, with s and t of opposite signs, belongs
    // to the remainder s a + t b = 2^k r + s e_a + t e_b of the full values,
    // which is at least 2^k (r - |s|) or 2^k (r - |t|), whichever is
    // negative, and is exact when k = 0. A step is taken only while that
    // bound is at least 2^t (or non-negative): then the matrix takes (a, b)
    // to a non-negative pair, with no need for the quotients to be those of
    // the full values. Every entry is also checked to fit a limb, which with
    // the bound stops the loop near the middle of the 128 bits.
    let need: u128 = match floor {
        None => 0,
        Some(t) if t.saturating_sub(k) < LEAD_BITS => 1 << t.saturating_sub(k),
        Some(_) => return None,
    };
    let mut l = Lehmer {
        u: leading(a, k),
        v: leading(b, k),
        m: [1, 0, 0, 1],
        need,
        exact: k == 0,
    };
    debug_assert!(l.v <= l.u);
    // Two steps a round, so that each knows where its negative entry is.
    let odd = loop {
        if !l.step::<false>() {
            break false;
        }
        if !l.step::<true>() {
            break true;
        }
    };
    let [ma, mb, mc, md] = l.m.map(i128::from);
    let signed = |x: i128, negative: bool| if negative { -x } else { x };
    (mb != 0).then(|| Matrix {
        a: signed(ma, odd),
        b: signed(mb, !odd),
        c: signed(mc, !odd),
        d: signed(md, odd),
    })
}

/// The state of Lehmer's inner loop: the leading parts `u >= v` and the
/// magnitudes `[ma, mb, mc, md]` of their rows `(s, t)`, with `r = s u0 + t
/// v0` for the starting parts `u0, v0`. The signs alternate: after an even
/// number of steps the negative entry of `u`'s row is `t` and that of `v`'s
/// is `s`, and the other way round after an odd number.
struct Lehmer {
    u: u128,
    v: u128,
    m: [Limb; 4],
    /// The least admissible `r - |negative entry|`.
    need: u128,
    /// Whether the leading parts are the whole values.
    exact: bool,
}

impl Lehmer {
    /// One Euclidean step, after an odd number of them if `ODD`; whether it
    /// was admissible and taken.
    #[inline]
    fn step<const ODD: bool>(&mut self) -> bool {
        if self.v == 0 {
            return false;
        }
        let (q, w) = div(self.u, self.v);
        let [ma, mb, mc, md] = self.m;
        let (Some(c), Some(d)) = (cofactor(q, ma, mc), cofactor(q, mb, md)) else {
            return false;
        };
        // w's row comes two after u's, so its negative entry is in the same
        // place as u's.
        let negative = if self.exact {
            0
        } else if ODD {
            c
        } else {
            d
        };
        if w < self.need + u128::from(negative) {
            return false;
        }
        self.m = [mc, md, c, d];
        (self.u, self.v) = (self.v, w);
        true
    }
}

/// `p + q c` if it fits a limb.
#[inline]
fn cofactor(q: u128, p: Limb, c: Limb) -> Option<Limb> {
    let q = Limb::try_from(q).ok()?;
    let x = u128::from(q) * u128::from(c) + u128::from(p);
    Limb::try_from(x).ok()
}

/// `(u / v, u % v)` for `u >= v > 0`. Most quotients are small, and a
/// 128-bit division is far slower than a few subtractions.
#[expect(
    clippy::cast_possible_truncation,
    reason = "both values are shifted to below 2^64"
)]
fn div(u: u128, v: u128) -> (u128, u128) {
    let mut r = u - v;
    for q in 1..=4 {
        if r < v {
            return (q, r);
        }
        r -= v;
    }
    if let (Ok(u), Ok(v)) = (u64::try_from(u), u64::try_from(v)) {
        return (u128::from(u / v), u128::from(u % v));
    }
    // A 64-bit division of the leading bits of both, which never
    // overestimates; when v keeps at least 32 of those bits, q is at most
    // 2^32 and the estimate falls short of it by at most a couple.
    let shift = 64 - u.leading_zeros();
    let (u1, v1) = ((u >> shift) as u64, (v >> shift) as u64);
    if v1 >> 32 == 0 {
        return (u / v, u % v);
    }
    let mut q = u128::from(u1 / (v1 + 1));
    let mut r = u - q * v;
    while r >= v {
        r -= v;
        q += 1;
    }
    (q, r)
}

/// `(a, b) = M (a, b)`, with `buf` as working space; the old buffers of `a`
/// and `b` become `buf`, so repeated calls do not allocate.
fn apply(m: &Matrix, a: &mut Vec<Limb>, b: &mut Vec<Limb>, buf: &mut [Vec<Limb>; 2]) {
    let n = a.len();
    b.resize(n, 0);
    let [x, y] = buf;
    x.resize(n, 0);
    y.resize(n, 0);
    combine(x, y, m, a, b);
    x.truncate(kernel::normalized_len(x));
    y.truncate(kernel::normalized_len(y));
    std::mem::swap(a, x);
    std::mem::swap(b, y);
}

/// `x = m.a a + m.b b` and `y = m.c a + m.d b` in one pass over `a` and `b`,
/// where each row has entries of opposite signs (or a zero) and the results
/// are known to be non-negative and to fit `x` and `y`.
fn combine(x: &mut [Limb], y: &mut [Limb], m: &Matrix, a: &[Limb], b: &[Limb]) {
    let (x_rows, y_rows) = (Row::new(m.a, m.b), Row::new(m.c, m.d));
    let (mut cx, mut cy) = (Row::CLEAR, Row::CLEAR);
    for (((x, y), &a), &b) in x.iter_mut().zip(y.iter_mut()).zip(a).zip(b) {
        *x = x_rows.limb(a, b, &mut cx);
        *y = y_rows.limb(a, b, &mut cy);
    }
    assert!(
        cx.0 == cx.1 && cy.0 == cy.1,
        "Lehmer step produced an out-of-range value"
    );
}

/// A row of a Lehmer matrix as the magnitudes of the entry for the
/// positive operand and of the one for the negative.
#[derive(Clone, Copy)]
struct Row {
    plus: Limb,
    minus: Limb,
    /// Whether `a` is the positive operand.
    a_plus: bool,
}

impl Row {
    /// No carry and no borrow.
    const CLEAR: (Limb, Limb) = (0, 0);

    fn new(p: i128, q: i128) -> Row {
        if q <= 0 {
            Row {
                plus: magnitude(p),
                minus: magnitude(-q),
                a_plus: true,
            }
        } else {
            Row {
                plus: magnitude(q),
                minus: magnitude(-p),
                a_plus: false,
            }
        }
    }

    /// The next limb of the row applied to `(a, b)`, with the running carry
    /// of the positive product and borrow of the negative one.
    #[inline]
    #[expect(clippy::cast_possible_truncation, reason = "taking the low limb")]
    fn limb(self, a: Limb, b: Limb, (carry, borrow): &mut (Limb, Limb)) -> Limb {
        let (p, n) = if self.a_plus { (a, b) } else { (b, a) };
        let t = u128::from(self.plus) * u128::from(p) + u128::from(*carry);
        let s = u128::from(self.minus) * u128::from(n) + u128::from(*borrow);
        let (lo, under) = (t as Limb).overflowing_sub(s as Limb);
        *carry = (t >> 64) as Limb;
        *borrow = (s >> 64) as Limb + Limb::from(under);
        lo
    }
}

/// A non-negative matrix entry as a limb.
fn magnitude(x: i128) -> Limb {
    Limb::try_from(x).expect("Lehmer cofactor in 0..2^64")
}
