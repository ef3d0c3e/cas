//! Greatest common divisor and least common multiple, by Lehmer's algorithm.
//!
//! Lehmer's algorithm (Knuth, TAOCP vol. 2, §4.5.2, Algorithm L) runs Euclid's
//! algorithm on the leading 63 bits of the operands for as long as the
//! quotients provably agree with those of the full values, collecting them in a
//! 2×2 matrix, then applies the matrix to the full values with one pass of
//! single-limb multiplications. When the leading bits cannot decide a single
//! quotient, a full division step is taken instead. Once both values fit a
//! limb the binary algorithm (or a word-level extended Euclid) finishes.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithm uses Knuth's names (a, b, u, v, q) and matrix entries"
)]

use std::cmp::Ordering;

use super::Natural;
use crate::{Integer, Limb, kernel};

/// Bits of the leading parts the inner loop works on. With 63 bits every
/// cofactor and every sum in the loop fits comfortably in an `i128`, and every
/// matrix entry has a magnitude below `2^63`.
const LEAD_BITS: u64 = 63;

/// The matrix `((a, b), (c, d))` from Lehmer's inner loop: the pair `(x, y)`
/// becomes `(a x + b y, c x + d y)`. It is the product of the Euclidean steps
/// of the full values, so it is unimodular, the entries of each row have
/// opposite signs (or one is zero), and both results are non-negative.
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
        // sa and sb are the coefficients of `self` in a and b.
        let (mut sa, mut sb) = if self >= other {
            (Integer::ONE, Integer::ZERO)
        } else {
            (Integer::ZERO, Integer::ONE)
        };
        reduce(&mut a, &mut b, dc, true, |step| {
            track(&mut sa, &mut sb, &step)
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
        } else if let Some(m) = lehmer_matrix(a, b, None) {
            apply(&m, a, b, &mut buf);
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

/// Bits `k .. k + 63` of `x`.
fn leading(x: &[Limb], k: u64) -> u64 {
    let i = usize::try_from(k / 64).expect("bit index fits usize");
    let s = k % 64;
    let lo = x.get(i).map_or(0, |&l| l >> s);
    let hi = match x.get(i + 1) {
        Some(&l) if s != 0 => l << (64 - s),
        _ => 0,
    };
    lo | hi
}

/// Lehmer's inner loop on the leading bits of `a >= b`, where `a` has at
/// least two limbs. With a `floor` of `t`, only steps whose remainder is
/// provably at least `2^t` are taken. `None` when not even the first quotient
/// is decided.
fn lehmer_matrix(a: &[Limb], b: &[Limb], floor: Option<u64>) -> Option<Matrix> {
    let k = kernel::bit_length(a) - LEAD_BITS;
    // With e as below, the new remainder is mc a + md b = 2^k v + mc e_a +
    // md e_b, where mc and md have opposite signs: it is at least
    // 2^k (v + min(mc, md)), and so at least 2^t once that sum is 2^(t - k).
    let need = match floor {
        None => None,
        Some(t) if t <= k => Some(1i128),
        Some(t) if t - k < LEAD_BITS => Some(1i128 << (t - k)),
        Some(_) => return None,
    };
    let (mut u, mut v) = (i128::from(leading(a, k)), i128::from(leading(b, k)));
    debug_assert!(u < 1 << LEAD_BITS && v <= u);
    let (mut ma, mut mb, mut mc, mut md) = (1i128, 0i128, 0i128, 1i128);
    // a = u 2^k + e with 0 <= e < 2^k (and alike for b), so the true quotient
    // of the current pair lies between (u + ma) / (v + mc) and
    // (u + mb) / (v + md); when both agree, it is decided. The four sums stay
    // in 0 ..= 2^63 (Knuth), so only a zero denominator needs checking.
    loop {
        let (n0, d0, n1, d1) = (u + ma, v + mc, u + mb, v + md);
        if d0 == 0 || d1 == 0 {
            break;
        }
        let q = quotient(n0, d0);
        if q != quotient(n1, d1) {
            break;
        }
        let (c, d, w) = (ma - q * mc, mb - q * md, u - q * v);
        if need.is_some_and(|need| w + c.min(d) < need) {
            break;
        }
        (ma, mc) = (mc, c);
        (mb, md) = (md, d);
        (u, v) = (v, w);
    }
    (mb != 0).then_some(Matrix {
        a: ma,
        b: mb,
        c: mc,
        d: md,
    })
}

/// `n / d` for `0 <= n <= 2^63` and `0 < d`, by a 64-bit division, which is
/// much cheaper than a 128-bit one.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "both operands are in 0 ..= 2^63"
)]
fn quotient(n: i128, d: i128) -> i128 {
    debug_assert!((0..=1 << 63).contains(&n) && (1..=1 << 63).contains(&d));
    i128::from(n as u64 / d as u64)
}

/// `(a, b) = M (a, b)`, with `buf` as working space; the old buffers of `a`
/// and `b` become `buf`, so repeated calls do not allocate.
fn apply(m: &Matrix, a: &mut Vec<Limb>, b: &mut Vec<Limb>, buf: &mut [Vec<Limb>; 2]) {
    let n = a.len();
    b.resize(n, 0);
    let [x, y] = buf;
    x.resize(n, 0);
    y.resize(n, 0);
    combine(x, m.a, a, m.b, b);
    combine(y, m.c, a, m.d, b);
    x.truncate(kernel::normalized_len(x));
    y.truncate(kernel::normalized_len(y));
    std::mem::swap(a, x);
    std::mem::swap(b, y);
}

/// `r = p a + q b` for `p` and `q` of opposite signs (or one zero), where the
/// result is known to be non-negative and to fit `r`.
fn combine(r: &mut [Limb], p: i128, a: &[Limb], q: i128, b: &[Limb]) {
    let (plus, pm, minus, mm) = if q <= 0 { (a, p, b, -q) } else { (b, q, a, -p) };
    let carry = kernel::mul_1(r, plus, magnitude(pm));
    let borrow = kernel::submul_1(r, minus, magnitude(mm));
    assert_eq!(carry, borrow, "Lehmer step produced an out-of-range value");
}

/// A non-negative matrix entry as a limb.
fn magnitude(x: i128) -> Limb {
    Limb::try_from(x).expect("Lehmer cofactor in 0..2^64")
}
