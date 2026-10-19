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

use hgcd::{HgcdMatrix, HgcdRow};

impl Natural {
    /// The greatest common divisor; `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, other: &Natural) -> Natural {
        self.gcd_from(other, hgcd::GCD_DC_THRESHOLD)
    }

    /// [`Natural::gcd`] with half-GCD steps from `dc` limbs.
    fn gcd_from(&self, other: &Natural, dc: usize) -> Natural {
        let (mut a, mut b) = ordered(self, other);
        halve(&mut a, &mut b, dc);
        lehmer(&mut a, &mut b, |_| {});
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
        self.div_exact(&self.gcd(other)) * other
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
        let mut s = Cofactors::new(self < other, a.len());
        let mut halved = false;
        while b.len() >= dc.max(2) {
            s.halve(&mut a, &mut b);
            halved = true;
        }
        if !halved || b.len() < 2 {
            return finish(a, b, s);
        }
        // As in GMP's `mpn_gcdext`, the Lehmer steps that remain start from
        // fresh cofactors, x of a and y of b with g = x a + y b, which stay
        // as short as a rather than growing from the length of s. The
        // coefficient of the first operand is then x sa + y sb = σ (x u -
        // y v), where x and y have opposite signs (or x is zero).
        let (a0, b0) = (Natural::from_limbs(&a), Natural::from_limbs(&b));
        let n = a.len();
        let (g, x) = finish(a, b, Cofactors::new(false, n));
        let (x_negative, x) = x.into_parts();
        let positive = !x_negative && !x.is_zero();
        // y = (g - x a0) / b0, of sign opposite to x.
        let mut p = &x * &a0;
        if positive {
            p -= &g;
        } else {
            p += &g;
        }
        let y = p.div_exact(&b0);
        let mag = &x * &Natural::from_vec(s.u) + &y * &Natural::from_vec(s.v);
        let negative = if positive { s.negative } else { !s.negative };
        (g, Integer::from_parts(negative, mag))
    }
}

/// The GCD of the normalised `a >= b`, by Lehmer steps, and the coefficient
/// of the first operand, whose coefficients in `a` and `b` are `s`.
fn finish(mut a: Vec<Limb>, mut b: Vec<Limb>, mut s: Cofactors) -> (Natural, Integer) {
    lehmer(&mut a, &mut b, |step| s.update(&step));
    if b.is_empty() {
        return (Natural::from_vec(a), s.combine(1, 0));
    }
    if a.len() > 1 {
        let q = div_step(&mut a, &mut b);
        s.update(&Step::Quotient(&q));
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
    (Natural::from(x), s.combine(p.0, p.1))
}

/// The coefficients of the first operand in the pair `(a, b)` being reduced:
/// `σ u` in `a` and `-σ v` in `b`, with `σ = -1` if `negative`, as GMP's
/// `mpn_gcdext` keeps them. Every step subtracts a multiple of one value from
/// the other or swaps them, so the coefficients keep opposite signs (or one
/// is zero), and each row of a step's matrix adds the magnitudes.
struct Cofactors {
    u: Vec<Limb>,
    v: Vec<Limb>,
    negative: bool,
}

impl Cofactors {
    /// The coefficients `(1, 0)`, or `(0, 1)` if the operands were swapped,
    /// with room for `n + 2` limbs: the coefficients stay below the larger
    /// operand, of `n` limbs, while both values are positive, and
    /// [`mul_row`] works on two limbs more.
    fn new(swapped: bool, n: usize) -> Cofactors {
        let room = || Vec::with_capacity(n + 2);
        let (mut u, mut v) = (room(), room());
        if swapped { &mut v } else { &mut u }.push(1);
        Cofactors {
            u,
            v,
            negative: swapped,
        }
    }

    /// Follows `step` on the pair.
    fn update(&mut self, step: &Step<'_>) {
        match step {
            Step::Matrix(m) => {
                // σ (m.a u - m.b v) and σ (m.c u - m.d v), where the entries
                // of each row have opposite signs and the rows opposite
                // patterns.
                let flip = m.a < 0 || m.b > 0;
                debug_assert_eq!(flip, m.c > 0 || m.d < 0);
                let n = [[entry(m.a), entry(m.c)], [entry(m.b), entry(m.d)]];
                mul_row(&mut self.u, &mut self.v, &n);
                self.negative ^= flip;
            }
            Step::Hgcd(m) => {
                // M^-1 = ±[[u11, -u01], [-u10, u00]] takes (σ u, -σ v) to
                // ±(σ (u11 u + u01 v), -σ (u10 u + u00 v)).
                let [[u00, u01], [u10, u11]] = &m.m;
                let u = Natural::from_vec(std::mem::take(&mut self.u));
                let v = Natural::from_vec(std::mem::take(&mut self.v));
                self.u = (&(u11 * &u) + &(u01 * &v)).limbs().to_vec();
                self.v = (&(u10 * &u) + &(u00 * &v)).limbs().to_vec();
                self.negative ^= m.negative;
            }
            Step::Quotient(q) => {
                // (σ u, -σ v) becomes (-σ v, σ (u + q v)).
                let w = add_product(&self.u, q, &self.v);
                self.u = std::mem::replace(&mut self.v, w);
                self.negative = !self.negative;
            }
        }
    }

    /// One half-GCD step on `(a, b)`, normalised with `a >= b`, or a
    /// division step if it can take none, followed on the coefficients. The
    /// row `(v, u)` is multiplied by the step's matrix on the right: while
    /// it is shorter than half of `a`, about the length the matrix entries
    /// reach, step by step inside the half-GCD at its top level, which
    /// costs less than building the matrix there; otherwise by the matrix.
    fn halve(&mut self, a: &mut Vec<Limb>, b: &mut Vec<Limb>) {
        let progress = if 2 * self.u.len().max(self.v.len()) < a.len() {
            let mut row = HgcdRow {
                r: [&self.v, &self.u].map(|x| Natural::from_limbs(x)),
                negative: false,
            };
            let progress = hgcd::hgcd_row(a, b, &mut row);
            let [v, u] = &row.r;
            (self.u, self.v) = (u.limbs().to_vec(), v.limbs().to_vec());
            self.negative ^= row.negative;
            progress
        } else if let Some(m) = hgcd::hgcd(a, b) {
            self.update(&Step::Hgcd(&m));
            true
        } else {
            false
        };
        if !progress {
            let q = div_step(a, b);
            self.update(&Step::Quotient(&q));
        }
        debug_assert!(kernel::is_normalized(a) && kernel::is_normalized(b));
        debug_assert_ne!(kernel::cmp(a, b), Ordering::Less);
    }

    /// The coefficient `p sa + q sb` for coefficients `sa` of `a` and `sb` of
    /// `b`, where `p` and `q` have opposite signs (or one is zero).
    fn combine(mut self, p: i128, q: i128) -> Integer {
        let flip = p < 0 || q > 0;
        mul_row(&mut self.u, &mut self.v, &[[entry(p), 0], [entry(q), 0]]);
        Integer::from_parts(self.negative ^ flip, Natural::from_vec(self.u))
    }
}

/// The magnitude of a matrix entry or word-level cofactor, below `2^64`.
#[expect(clippy::cast_possible_truncation, reason = "checked in debug builds")]
fn entry(x: i128) -> Limb {
    debug_assert!(x.unsigned_abs() <= u128::from(Limb::MAX));
    x.unsigned_abs() as Limb
}

/// `u + q v`.
fn add_product(u: &[Limb], q: &[Limb], v: &[Limb]) -> Vec<Limb> {
    if q.is_empty() || v.is_empty() {
        return u.to_vec();
    }
    let len = q.len() + v.len();
    let mut w = vec![0; len.max(u.len()) + 1];
    let mut scratch = vec![0; kernel::mul_scratch_len(q.len(), v.len())];
    kernel::mul(&mut w[..len], q, v, &mut scratch);
    let carry = kernel::add_assign(&mut w, u);
    debug_assert_eq!(carry, 0);
    w.truncate(kernel::normalized_len(&w));
    w
}

/// `(x, y) = (x n00 + y n10, x n01 + y n11)`, the row `(x, y)` times the
/// matrix `n`, in place and in one pass.
fn mul_row(x: &mut Vec<Limb>, y: &mut Vec<Limb>, n: &[[Limb; 2]; 2]) {
    let len = x.len().max(y.len());
    // Each result is below 2 B^(len + 1).
    x.resize(len + 2, 0);
    y.resize(len + 2, 0);
    let [[n00, n01], [n10, n11]] = *n;
    let (mut cx, mut cy) = (Row::CLEAR, Row::CLEAR);
    for (x, y) in x[..len].iter_mut().zip(&mut y[..len]) {
        let (a, b) = (*x, *y);
        *x = mul_add_2(n00, a, n10, b, &mut cx);
        *y = mul_add_2(n01, a, n11, b, &mut cy);
    }
    for (r, (c0, c1)) in [(x, cx), (y, cy)] {
        let (lo, carry) = c0.overflowing_add(c1);
        (r[len], r[len + 1]) = (lo, Limb::from(carry));
        r.truncate(kernel::normalized_len(r));
    }
}

/// The next limb of `a x + b y`, with the running carries of both products.
#[inline]
#[expect(clippy::cast_possible_truncation, reason = "taking the low limb")]
fn mul_add_2(a: Limb, x: Limb, b: Limb, y: Limb, (c0, c1): &mut (Limb, Limb)) -> Limb {
    let s = u128::from(a) * u128::from(x) + u128::from(*c0);
    let t = u128::from(b) * u128::from(y) + u128::from(*c1);
    let (lo, carry) = (s as Limb).overflowing_add(t as Limb);
    // s < 2^128 - 2^64 + 1, so its high limb is below 2^64 - 1 unless its
    // low limb is zero, and the carry fits.
    *c0 = (s >> 64) as Limb + Limb::from(carry);
    *c1 = (t >> 64) as Limb;
    lo
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

/// Reduces normalised `a >= b` by half-GCD steps (or a division step where
/// one takes none) while `b` has at least `dc` limbs, keeping `gcd(a, b)`
/// and `a >= b`, without building their matrices at the top level.
fn halve(a: &mut Vec<Limb>, b: &mut Vec<Limb>, dc: usize) {
    while b.len() >= dc.max(2) {
        if !hgcd::hgcd_values(a, b) {
            div_step(a, b);
        }
        debug_assert!(kernel::is_normalized(a) && kernel::is_normalized(b));
        // Equal only when a half-GCD stopped at a gcd of at least B^s.
        debug_assert_ne!(kernel::cmp(a, b), Ordering::Less);
    }
}

/// Reduces normalised `a >= b` until `b` has at most one limb, keeping
/// `gcd(a, b)` and `a >= b`, by Lehmer steps, and division steps where the
/// leading bits cannot decide one; reports every step to `on_step`.
fn lehmer(a: &mut Vec<Limb>, b: &mut Vec<Limb>, mut on_step: impl FnMut(Step<'_>)) {
    let mut buf = [Vec::new(), Vec::new()];
    while b.len() >= 2 {
        if let Some(mut m) = lehmer_matrix(a, b, None) {
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
