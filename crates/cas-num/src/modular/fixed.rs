//! Montgomery arithmetic modulo an odd modulus of a fixed, small number of
//! limbs, on stack arrays.
//!
//! [`super::ModRing`] allocates a [`super::Residue`] and the product buffers
//! for every operation, which dominates the arithmetic itself at three to
//! eight limbs. [`ModFixed`] keeps its elements in `[Limb; K]` arrays and
//! its double-length products in a stack buffer, so a power or a Lucas
//! sequence runs without touching the heap. With `K` known, the loops over
//! limbs unroll and the operands stay in registers: a product interleaves
//! the multiplication with the reduction a limb at a time, and a square
//! computes each cross product once before a separate REDC.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, k, r, t, x, y)"
)]

use super::ring::limb_inverse;
use super::word::sliding_window;
use crate::kernel::{self, adc, sbb};
use crate::{Limb, Natural};

#[cfg(test)]
mod tests;

/// The largest modulus length, in limbs, that [`ModFixed`] serves.
pub(crate) const FIXED_MAX: usize = 8;

/// An element of a [`ModFixed`], in Montgomery form: `x R mod n` with
/// `R = B^K`, below `n`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Fixed<const K: usize>([Limb; K]);

impl<const K: usize> Default for Fixed<K> {
    fn default() -> Self {
        Fixed([0; K])
    }
}

/// Montgomery arithmetic modulo an odd `K`-limb `n`, `2 <= K <=
/// FIXED_MAX`. The arithmetic methods take elements of this ring (below
/// `n`); others give an unspecified result.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ModFixed<const K: usize> {
    n: [Limb; K],
    /// `-n^(-1) mod B`.
    ninv: Limb,
    /// `R mod n`, the element 1.
    one: Fixed<K>,
    /// `R^2 mod n`, for conversion into Montgomery form.
    r2: Fixed<K>,
}

/// The limbs of `a < B^K`, zero-padded.
fn load<const K: usize>(a: &Natural) -> [Limb; K] {
    let mut x = [0; K];
    x[..a.limb_len()].copy_from_slice(a.limbs());
    x
}

impl<const K: usize> ModFixed<K> {
    /// The arithmetic modulo `n`, odd and of exactly `K` limbs.
    pub(crate) fn new(n: &Natural) -> Self {
        assert!(
            (2..=FIXED_MAX).contains(&K) && n.limb_len() == K && n.bit(0),
            "ModFixed: modulus must be odd and of K limbs"
        );
        let limbs = load(n);
        let r = Natural::power_of_two(kernel::limbs_to_bits(K));
        ModFixed {
            n: limbs,
            ninv: limb_inverse(limbs[0]).wrapping_neg(),
            one: Fixed(load(&(&r % n))),
            r2: Fixed(load(&(r.square() % n))),
        }
    }

    /// The element `a mod n`. A value of up to `K` limbs needs no division:
    /// its product with `R^2 mod n` is reduced (see [`ModFixed::mul`]).
    pub(crate) fn elem(&self, a: &Natural) -> Fixed<K> {
        let x = if a.limb_len() <= K {
            load(a)
        } else {
            load(&(a % &Natural::from_limbs(&self.n)))
        };
        self.mul(Fixed(x), self.r2)
    }

    /// The value of `x`, in `0..n`.
    pub(crate) fn value(&self, x: Fixed<K>) -> Natural {
        let mut t = [0; 2 * FIXED_MAX];
        t[..K].copy_from_slice(&x.0);
        Natural::from_limbs(&self.redc(&mut t).0)
    }

    /// The element 1.
    pub(crate) fn one(&self) -> Fixed<K> {
        self.one
    }

    /// `x + y`.
    pub(crate) fn add(&self, x: Fixed<K>, y: Fixed<K>) -> Fixed<K> {
        let (r, carry) = add(&x.0, &y.0);
        self.normalize(r, Limb::from(carry))
    }

    /// `x - y`.
    pub(crate) fn sub(&self, x: Fixed<K>, y: Fixed<K>) -> Fixed<K> {
        let (r, borrow) = sub(&x.0, &y.0);
        if borrow {
            Fixed(add(&r, &self.n).0)
        } else {
            Fixed(r)
        }
    }

    /// `-x`.
    pub(crate) fn neg(&self, x: Fixed<K>) -> Fixed<K> {
        if x.0 == [0; K] {
            return x;
        }
        Fixed(sub(&self.n, &x.0).0)
    }

    /// `x y`, by Montgomery multiplication with the product and the reduction
    /// interleaved a limb at a time (coarsely integrated operand scanning):
    /// `t = (t + x_i y + q n) / B` with `q = (t + x_i y) (-n^(-1)) mod B`.
    /// `t` stays below `2n`, in `K` limbs and a top bit: with `t < 2n`,
    /// `x_i, q < B` and `y < n`, the new `t` is below `(2n + (B - 1) n +
    /// (B - 1) n) / B = 2n`. So only `y` must be reduced; `x` may be any
    /// `K`-limb value.
    pub(crate) fn mul(&self, x: Fixed<K>, y: Fixed<K>) -> Fixed<K> {
        let (x, y, n) = (&x.0, &y.0, &self.n);
        let mut t = [0; K];
        let mut top: Limb = 0;
        for &xi in x {
            let mut c = 0;
            for (tj, &yj) in t.iter_mut().zip(y) {
                (*tj, c) = mac(*tj, xi, yj, c);
            }
            let (high, over) = top.overflowing_add(c);
            let q = t[0].wrapping_mul(self.ninv);
            let (_, mut c) = mac(t[0], q, n[0], 0);
            for j in 1..K {
                (t[j - 1], c) = mac(t[j], q, n[j], c);
            }
            let carry;
            (t[K - 1], carry) = high.overflowing_add(c);
            top = Limb::from(over) + Limb::from(carry);
        }
        self.normalize(t, top)
    }

    /// `x^2`: the off-diagonal products `x_i x_j`, `i < j`, once, doubled,
    /// plus the squares `x_i^2`, then REDC a limb at a time.
    pub(crate) fn square(&self, x: Fixed<K>) -> Fixed<K> {
        let x = &x.0;
        let mut t = [0; 2 * FIXED_MAX];
        for i in 0..K - 1 {
            let mut c = 0;
            for j in i + 1..K {
                (t[i + j], c) = mac(t[i + j], x[i], x[j], c);
            }
            t[i + K] = c;
        }
        let mut high = 0;
        for l in &mut t[..2 * K] {
            (*l, high) = (*l << 1 | high, *l >> (Limb::BITS - 1));
        }
        let mut carry = false;
        for i in 0..K {
            let (lo, hi) = mac(0, x[i], x[i], 0);
            let c;
            (t[2 * i], c) = adc(t[2 * i], lo, carry);
            (t[2 * i + 1], carry) = adc(t[2 * i + 1], hi, c);
        }
        self.redc(&mut t)
    }

    /// `x^e`, with `x^0 = 1`, by sliding windows.
    pub(crate) fn pow(&self, x: Fixed<K>, e: &Natural) -> Fixed<K> {
        let bits = e.bit_length();
        if bits == 0 {
            return self.one;
        }
        // The window loop squares as `mul(r, r)`; equal operands take the
        // cheaper squaring, which gives the same value.
        sliding_window(bits, &|i| e.bit(i), x, |p, q| {
            if p == q {
                self.square(p)
            } else {
                self.mul(p, q)
            }
        })
    }

    /// `t R^(-1) mod n` for `t < n R` in `t[..2K]`, eliminating a limb per
    /// pass and keeping each pass's carry in the limb it cleared, which no
    /// later pass reads, until the two halves are added.
    fn redc(&self, t: &mut [Limb; 2 * FIXED_MAX]) -> Fixed<K> {
        let n = &self.n;
        for i in 0..K {
            let q = t[i].wrapping_mul(self.ninv);
            let mut c = 0;
            for j in 0..K {
                (t[i + j], c) = mac(t[i + j], q, n[j], c);
            }
            t[i] = c;
        }
        let mut r = [0; K];
        let mut carry = false;
        for i in 0..K {
            (r[i], carry) = adc(t[K + i], t[i], carry);
        }
        self.normalize(r, Limb::from(carry))
    }

    /// `t + top B^K - n` if that is not negative, else `t`, for a value
    /// below `2n`.
    fn normalize(&self, t: [Limb; K], top: Limb) -> Fixed<K> {
        let (d, borrow) = sub(&t, &self.n);
        if top != 0 || !borrow {
            Fixed(d)
        } else {
            Fixed(t)
        }
    }
}

/// `(low, high)` of `a + b c + d`, which fits a double limb.
#[inline]
fn mac(a: Limb, b: Limb, c: Limb, d: Limb) -> (Limb, Limb) {
    let t = u128::from(a) + u128::from(b) * u128::from(c) + u128::from(d);
    #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
    (t as Limb, (t >> Limb::BITS) as Limb)
}

/// `x + y` and the carry out.
#[inline]
fn add<const K: usize>(x: &[Limb; K], y: &[Limb; K]) -> ([Limb; K], bool) {
    let mut r = [0; K];
    let mut carry = false;
    for i in 0..K {
        (r[i], carry) = adc(x[i], y[i], carry);
    }
    (r, carry)
}

/// `x - y` and the borrow out.
#[inline]
fn sub<const K: usize>(x: &[Limb; K], y: &[Limb; K]) -> ([Limb; K], bool) {
    let mut r = [0; K];
    let mut borrow = false;
    for i in 0..K {
        (r[i], borrow) = sbb(x[i], y[i], borrow);
    }
    (r, borrow)
}
