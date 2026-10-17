//! Subquadratic half-GCD (Möller, "On Schönhage's algorithm and subquadratic
//! integer GCD computation", Math. Comp. 77 (2008); GMP's `mpn_hgcd`).
//!
//! [`hgcd`] takes a pair of at most `n` limbs and reduces it by unimodular
//! steps for as long as both values stay at least `B^s`, with `B = 2^64` and
//! `s = ⌊n/2⌋ + 1`. Above [`HGCD_THRESHOLD`] it recurses: a half-GCD of the
//! top half of the limbs, lifted to the full values, one division step, a
//! half-GCD of the top of what remains, lifted again, and Lehmer steps to
//! finish. Below it, Lehmer steps do all the work.
//!
//! Why lifting is exact. Let `(a, b) = M (x, y)` with `M` non-negative of
//! determinant ±1 and `x, y >= B^s`. Every entry of `M` multiplies `x` or `y`
//! in `a` or `b`, so is at most `max(a, b) / min(x, y) < B^(n - s)`, and
//! `n - s < s`. If `(a, b)` are the limbs of `(A, B)` from limb `p` up, then
//! `M^-1 (A, B) = B^p (x, y) + M^-1 (A_lo, B_lo)`, and each component of the
//! second term exceeds `-B^p` times an entry of `M`. The lifted values are
//! therefore above `B^p (B^s - B^(n - s)) >= B^(p + s - 1)`: `M` reduces the
//! full values too, and keeps them at least `B^(p + s - 1)`. Both recursive
//! calls choose `p` so that this is at least the caller's own `B^s`.
#![allow(
    clippy::many_single_char_names,
    reason = "x, y, p, q, s, m follow the notation of the module comment"
)]

use std::cmp::Ordering;

use super::{Matrix, apply, lehmer_matrix};
use crate::{Limb, Natural, kernel};

/// Smaller operand length (limbs) from which `Natural::gcd` reduces by
/// half-GCD steps instead of Lehmer steps.
pub(super) const GCD_DC_THRESHOLD: usize = crate::tuning::GCD_DC_THRESHOLD;

/// Length (limbs) from which [`hgcd`] recurses; at least 3, the smallest `n`
/// for which the first recursive call keeps the values above `B^s`.
const HGCD_THRESHOLD: usize = crate::tuning::HGCD_THRESHOLD;

/// The steps of a half-GCD: a non-negative matrix `M` of determinant ±1 with
/// `(a, b) = M (x, y)`, where `(a, b)` is the pair before the steps and
/// `(x, y)` the pair after them.
#[derive(Clone, Debug)]
pub(super) struct HgcdMatrix {
    /// `[[u00, u01], [u10, u11]]`.
    pub(super) m: [[Natural; 2]; 2],
    /// Whether the determinant is -1.
    pub(super) negative: bool,
}

impl HgcdMatrix {
    fn identity() -> HgcdMatrix {
        HgcdMatrix {
            m: [
                [Natural::one(), Natural::ZERO],
                [Natural::ZERO, Natural::one()],
            ],
            negative: false,
        }
    }

    /// `M P` with `P` the swap, for `(x, y)` that changed places.
    pub(super) fn swap_columns(&mut self) {
        for row in &mut self.m {
            row.swap(0, 1);
        }
        self.negative = !self.negative;
    }

    /// `M N` for `N = [[1, q], [0, 1]]`, the step `x -= q y`.
    fn add_column(&mut self, q: &Natural) {
        for [u0, u1] in &mut self.m {
            *u1 += &(&*u0 * q);
        }
    }

    /// `M L^-1` for the Lehmer matrix `L` that took `(x, y)` to `L (x, y)`,
    /// in place with `buf` as working space. `L^-1` is the product of
    /// `[[q, 1], [1, 0]]` over its quotients, so its entries are the
    /// magnitudes of `L`'s, rearranged.
    fn mul_lehmer(&mut self, l: &Matrix, buf: &mut [Vec<Limb>; 2]) {
        let mag = |x: i128| Limb::try_from(x.unsigned_abs()).expect("Lehmer entry below 2^64");
        let n = [[mag(l.d), mag(l.b)], [mag(l.c), mag(l.a)]];
        for [u0, u1] in &mut self.m {
            for (j, r) in buf.iter_mut().enumerate() {
                add_products(r, u0.limbs(), n[0][j], u1.limbs(), n[1][j]);
            }
            assign(u0, &buf[0]);
            assign(u1, &buf[1]);
        }
        // det L^-1 = det L, and each Lehmer step has determinant -1.
        self.negative ^= l.a * l.d - l.b * l.c < 0;
    }

    /// The product `self other`.
    fn mul(&self, other: &HgcdMatrix) -> HgcdMatrix {
        let (a, b) = (&self.m, &other.m);
        let entry = |i: usize, j: usize| &(&a[i][0] * &b[0][j]) + &(&a[i][1] * &b[1][j]);
        HgcdMatrix {
            m: [[entry(0, 0), entry(0, 1)], [entry(1, 0), entry(1, 1)]],
            negative: self.negative != other.negative,
        }
    }
}

/// Reduces `(x, y)`, normalised with `x >= y` and of at most `n` limbs, by
/// steps that keep both values at least `B^s`, `s = ⌊n/2⌋ + 1`, until no
/// further step can (`|x - y| < B^s`). Returns the matrix of the steps, or
/// `None`, with `(x, y)` unchanged, if none was possible.
pub(super) fn hgcd(x: &mut Vec<Limb>, y: &mut Vec<Limb>) -> Option<HgcdMatrix> {
    let mut m = Some(HgcdMatrix::identity());
    if run(x, y, &mut m) { m } else { None }
}

/// [`hgcd`] without the matrix, which saves building it at this level;
/// whether any step was taken.
pub(super) fn hgcd_values(x: &mut Vec<Limb>, y: &mut Vec<Limb>) -> bool {
    run(x, y, &mut None)
}

/// [`hgcd`], appending its steps to `m` if there is one; whether any step was
/// taken.
fn run(x: &mut Vec<Limb>, y: &mut Vec<Limb>, m: &mut Option<HgcdMatrix>) -> bool {
    debug_assert_ne!(kernel::cmp(x, y), Ordering::Less);
    let n = x.len();
    let s = n / 2 + 1;
    // B^s has s + 1 limbs.
    if y.len() <= s {
        return false;
    }
    let mut progress = false;
    if n >= HGCD_THRESHOLD {
        // The top n - p limbs, s' = ⌊(n - p)/2⌋ + 1: p + s' - 1 >= s for n >= 3.
        progress |= reduce_top(x, y, n / 2, m);
        // A division step, so that a large quotient cannot leave the second
        // call with nothing to do and Lehmer steps with half the work.
        order(x, y, m);
        if division_step(x, y, s, m) {
            progress = true;
            order(x, y, m);
            // With len limbs left, the top 2 (len - s) - 1 of them give
            // s' = len - s, and p + s' - 1 = s exactly.
            let len = x.len();
            if len > s + 1 {
                reduce_top(x, y, 2 * s + 1 - len, m);
            }
        }
    }
    lehmer_steps(x, y, s, m) || progress
}

/// Half-GCD of the limbs of `(x, y)`, `x >= y`, from `p` up, applied to the
/// full values and appended to `m`; whether any step was taken.
fn reduce_top(x: &mut Vec<Limb>, y: &mut Vec<Limb>, p: usize, m: &mut Option<HgcdMatrix>) -> bool {
    let top = |v: &[Limb]| v.get(p..).map_or_else(Vec::new, <[Limb]>::to_vec);
    let (mut xh, mut yh) = (top(x), top(y));
    let Some(m1) = hgcd(&mut xh, &mut yh) else {
        return false;
    };
    let low = |v: &[Limb]| Natural::from_limbs(&v[..p.min(v.len())]);
    let (xl, yl) = (low(x), low(y));
    // (x, y) = M1^-1 (X, Y) = ±(u11 X - u01 Y, u00 Y - u10 X), split at limb p.
    let [[u00, u01], [u10, u11]] = &m1.m;
    let (mut xp, mut xm) = (u11 * &xl, u01 * &yl);
    let (mut yp, mut ym) = (u00 * &yl, u10 * &xl);
    if m1.negative {
        std::mem::swap(&mut xp, &mut xm);
        std::mem::swap(&mut yp, &mut ym);
    }
    *x = lift(xh, p, &xp, &xm);
    *y = lift(yh, p, &yp, &ym);
    if let Some(m) = m {
        *m = m.mul(&m1);
    }
    true
}

/// `hi B^p + plus - minus`, known to be non-negative.
fn lift(hi: Vec<Limb>, p: usize, plus: &Natural, minus: &Natural) -> Vec<Limb> {
    let mut v = Natural::from_vec(hi) << kernel::limbs_to_bits(p);
    v += plus;
    v -= minus;
    v.limbs().to_vec()
}

/// Puts the larger value first.
fn order(x: &mut Vec<Limb>, y: &mut Vec<Limb>, m: &mut Option<HgcdMatrix>) {
    if kernel::cmp(x, y) == Ordering::Less {
        std::mem::swap(x, y);
        if let Some(m) = m {
            m.swap_columns();
        }
    }
}

/// `r = a p + b q`, normalised, for `p, q < 2^63`.
fn add_products(r: &mut Vec<Limb>, a: &[Limb], p: Limb, b: &[Limb], q: Limb) {
    let (a, p, b, q) = if a.len() >= b.len() {
        (a, p, b, q)
    } else {
        (b, q, a, p)
    };
    r.resize(a.len() + 1, 0);
    r[a.len()] = kernel::mul_1(&mut r[..a.len()], a, p);
    let carry = kernel::addmul_1(&mut r[..b.len()], b, q);
    // a p + b q < 2 B^len(a) 2^63, so the top limb absorbs the carry.
    let overflow = kernel::add_1_assign(&mut r[b.len()..], carry);
    debug_assert_eq!(overflow, 0);
    r.truncate(kernel::normalized_len(r));
}

/// Sets `x` to the normalised limbs `v`, reusing its buffer.
fn assign(x: &mut Natural, v: &[Limb]) {
    x.limbs.resize(v.len());
    x.limbs.as_mut_slice().copy_from_slice(v);
}

/// Lehmer and division steps on `(x, y)`, keeping both at least `B^s`, until
/// `|x - y| < B^s`; whether any step was taken.
fn lehmer_steps(
    x: &mut Vec<Limb>,
    y: &mut Vec<Limb>,
    s: usize,
    m: &mut Option<HgcdMatrix>,
) -> bool {
    let floor = kernel::limbs_to_bits(s);
    let (mut values, mut entries) = ([Vec::new(), Vec::new()], [Vec::new(), Vec::new()]);
    let mut progress = false;
    loop {
        order(x, y, m);
        if let Some(l) = lehmer_matrix(x, y, Some(floor)) {
            apply(&l, x, y, &mut values);
            if let Some(m) = m {
                m.mul_lehmer(&l, &mut entries);
            }
        } else if !division_step(x, y, s, m) {
            return progress;
        }
        progress = true;
    }
}

/// `x -= q y` for `x >= y >= B^s` with the largest `q` that leaves
/// `x >= B^s`; `false`, with nothing changed, if `x - y < B^s`.
fn division_step(x: &mut Vec<Limb>, y: &[Limb], s: usize, m: &mut Option<HgcdMatrix>) -> bool {
    let mut q = vec![0; x.len() - y.len() + 1];
    let mut r = vec![0; y.len() + 1];
    let mut scratch = vec![0; kernel::div_rem_scratch_len(x.len(), y.len())];
    kernel::div_rem(&mut q, &mut r[..y.len()], x, y, &mut scratch);
    q.truncate(kernel::normalized_len(&q));
    if kernel::normalized_len(&r) <= s {
        if q == [1] {
            return false;
        }
        // q >= 2, so q - 1 >= 1 and r + y >= y >= B^s.
        kernel::sub_1_assign(&mut q, 1);
        q.truncate(kernel::normalized_len(&q));
        let carry = kernel::add_assign(&mut r[..y.len()], y);
        r[y.len()] = carry;
    }
    r.truncate(kernel::normalized_len(&r));
    if let Some(m) = m {
        m.add_column(&Natural::from_vec(q));
    }
    *x = r;
    true
}
