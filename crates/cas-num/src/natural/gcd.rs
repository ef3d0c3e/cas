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
    /// `(x, y)` became `(y, x - q y)`, with `q` the given limbs.
    Quotient(&'a [Limb]),
}

#[cfg(test)]
mod tests;

impl Natural {
    /// The greatest common divisor; `gcd(0, 0) = 0`.
    #[must_use]
    pub fn gcd(&self, other: &Natural) -> Natural {
        let (mut a, mut b) = ordered(self, other);
        reduce(&mut a, &mut b, |_| {});
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
        let (mut a, mut b) = ordered(self, other);
        // sa and sb are the coefficients of `self` in a and b.
        let (mut sa, mut sb) = if self >= other {
            (Integer::ONE, Integer::ZERO)
        } else {
            (Integer::ZERO, Integer::ONE)
        };
        reduce(&mut a, &mut b, |step| track(&mut sa, &mut sb, &step));
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
        Step::Quotient(q) => {
            let next = &*sa - &(Integer::from(Natural::from_limbs(q)) * &*sb);
            *sa = std::mem::replace(sb, next);
        }
    }
}

/// Reduces normalised `a >= b` until `b` has at most one limb, keeping
/// `gcd(a, b)` and `a >= b`, and reports every step to `on_step`.
fn reduce(a: &mut Vec<Limb>, b: &mut Vec<Limb>, mut on_step: impl FnMut(Step<'_>)) {
    while b.len() >= 2 {
        if let Some(m) = lehmer_matrix(a, b) {
            apply(&m, a, b);
            on_step(Step::Matrix(&m));
        } else {
            let q = div_step(a, b);
            on_step(Step::Quotient(&q));
        }
        debug_assert!(kernel::is_normalized(a) && kernel::is_normalized(b));
        debug_assert_eq!(kernel::cmp(a, b), Ordering::Greater);
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
/// least two limbs. `None` when not even the first quotient is decided.
fn lehmer_matrix(a: &[Limb], b: &[Limb]) -> Option<Matrix> {
    let k = kernel::bit_length(a) - LEAD_BITS;
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
        let q = n0 / d0;
        if q != n1 / d1 {
            break;
        }
        (ma, mc) = (mc, ma - q * mc);
        (mb, md) = (md, mb - q * md);
        (u, v) = (v, u - q * v);
    }
    (mb != 0).then_some(Matrix {
        a: ma,
        b: mb,
        c: mc,
        d: md,
    })
}

/// `(a, b) = M (a, b)`.
fn apply(m: &Matrix, a: &mut Vec<Limb>, b: &mut Vec<Limb>) {
    let n = a.len();
    b.resize(n, 0);
    let mut x = vec![0; n];
    let mut y = vec![0; n];
    combine(&mut x, m.a, a, m.b, b);
    combine(&mut y, m.c, a, m.d, b);
    x.truncate(kernel::normalized_len(&x));
    y.truncate(kernel::normalized_len(&y));
    (*a, *b) = (x, y);
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
