//! The ring ℤ/nℤ for a multi-precision modulus.
//!
//! Odd moduli use Montgomery's representation: a residue `x` is stored as
//! `x R mod n` with `R = B^k` for a `k`-limb modulus, and products are
//! reduced by REDC, which divides by `R` instead of by `n`. Below
//! [`REDC_N_THRESHOLD`] limbs REDC eliminates two limbs per pass over `n`
//! (GMP's `redc_2`); from it on, it multiplies by a precomputed
//! `-n^(-1) mod R` (a short product) and by `n` with the subquadratic
//! multiplication (`redc_n`). Even moduli store
//! plain values and reduce products by division.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, k, r, t, q, u)"
)]

use std::cmp::Ordering;

use crate::kernel::mul::mullo::{mullo, mullo_scratch_len};
use crate::kernel::mul::{addmul_2_carries, mul_wrap, mul_wrap_scratch_len, wraps_cyclically};
use crate::kernel::{
    self, add_1_assign, add_n, add_n_assign, cmp_n, div_rem, div_rem_scratch_len, mul_scratch_len,
    sqr_scratch_len, sub_1_assign, sub_assign, sub_n, sub_n_assign,
};
use crate::{Integer, Limb, Natural};

#[cfg(test)]
mod tests;

/// Smallest modulus length, in limbs, at which REDC multiplies by
/// `-n^(-1) mod R` instead of eliminating one limb at a time. Set by
/// the tuning configuration (`build.rs`).
pub(crate) const REDC_N_THRESHOLD: usize = crate::tuning::REDC_N_THRESHOLD;

/// The ring ℤ/nℤ for a non-zero modulus `n` of any size.
///
/// Elements are [`Residue`]s created by this ring ([`ModRing::elem`],
/// [`ModRing::zero`], [`ModRing::one`]) and combined by its methods. A
/// residue of a different ring gives an unspecified result (see the
/// [module documentation](super)).
#[derive(Clone, Debug)]
pub struct ModRing {
    n: Natural,
    form: Form,
}

/// An element of a [`ModRing`].
///
/// The representation is internal (Montgomery form for odd moduli) but
/// unique, so equality and hashing are by value within one ring.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Residue {
    /// Exactly as many limbs as the modulus, below it.
    limbs: Box<[Limb]>,
}

#[derive(Clone, Debug)]
enum Form {
    /// Odd modulus: residues are `x R mod n`.
    Montgomery {
        /// `-n^(-1) mod B`.
        ninv: Limb,
        /// `-n^(-1) mod B^2`, for `redc_2`.
        ninv2: u128,
        /// `-n^(-1) mod R`, for `redc_n` (empty below `REDC_N_THRESHOLD`).
        ninv_full: Vec<Limb>,
        /// `R mod n`, the residue 1.
        one: Box<[Limb]>,
        /// `R^2 mod n`, for conversion into Montgomery form.
        r2: Box<[Limb]>,
    },
    /// Even modulus: residues are plain values.
    Plain,
}

/// Buffers for one product and its reduction, reused across the steps of a
/// power.
struct Scratch {
    /// The double-length product.
    t: Vec<Limb>,
    /// `redc_n`'s two products, or the division quotient.
    u: Vec<Limb>,
    /// Scratch for the multiplication and division kernels.
    work: Vec<Limb>,
}

impl ModRing {
    /// The ring ℤ/nℤ.
    ///
    /// # Panics
    ///
    /// If `n` is zero.
    #[must_use]
    pub fn new(n: Natural) -> ModRing {
        assert!(!n.is_zero(), "zero modulus");
        let k = n.limb_len();
        let form = if n.bit(0) {
            let r = Natural::power_of_two(kernel::limbs_to_bits(k));
            let ninv_full = if k >= REDC_N_THRESHOLD {
                let inv = Integer::from(&n)
                    .mod_inverse(&Integer::from(&r))
                    .expect("odd modulus is invertible modulo a power of two")
                    .unsigned_abs();
                padded(&(&r - &inv), k).into_vec()
            } else {
                Vec::new()
            };
            let inv = limb_inverse(n.limbs()[0]);
            // One more Newton step takes the inverse to 128 bits.
            let low = u128::from(n.limbs()[0])
                | u128::from(n.limbs().get(1).copied().unwrap_or(0)) << Limb::BITS;
            let inv2 =
                u128::from(inv).wrapping_mul(2u128.wrapping_sub(low.wrapping_mul(u128::from(inv))));
            debug_assert_eq!(low.wrapping_mul(inv2), 1);
            Form::Montgomery {
                ninv: inv.wrapping_neg(),
                ninv2: inv2.wrapping_neg(),
                ninv_full,
                one: padded(&(&r % &n), k),
                r2: padded(&(r.square() % &n), k),
            }
        } else {
            Form::Plain
        };
        ModRing { n, form }
    }

    /// The modulus `n`.
    #[must_use]
    pub fn modulus(&self) -> &Natural {
        &self.n
    }

    fn len(&self) -> usize {
        self.n.limb_len()
    }

    fn residue(limbs: Box<[Limb]>) -> Residue {
        Residue { limbs }
    }

    fn check(&self, x: &Residue) {
        debug_assert!(
            x.limbs.len() == self.len()
                && cmp_n(&x.limbs, &padded(&self.n, self.len())) == Ordering::Less,
            "residue of another ring"
        );
    }

    /// The residue 0.
    #[must_use]
    pub fn zero(&self) -> Residue {
        Self::residue(vec![0; self.len()].into_boxed_slice())
    }

    /// The residue 1 (which is 0 modulo 1).
    #[must_use]
    pub fn one(&self) -> Residue {
        match &self.form {
            Form::Montgomery { one, .. } => Self::residue(one.clone()),
            Form::Plain => self.elem(&Natural::one()),
        }
    }

    /// The residue of `a`.
    #[must_use]
    pub fn elem(&self, a: &Natural) -> Residue {
        let x = padded(&(a % &self.n), self.len());
        match &self.form {
            Form::Montgomery { r2, .. } => self.mul_limbs(&x, r2),
            Form::Plain => Self::residue(x),
        }
    }

    /// The residue of `a`, also for negative `a`.
    #[must_use]
    pub fn elem_integer(&self, a: &Integer) -> Residue {
        let x = self.elem(&a.unsigned_abs());
        if a.is_negative() { self.neg(&x) } else { x }
    }

    /// The value of `x`, in `0..n`.
    #[must_use]
    pub fn to_natural(&self, x: &Residue) -> Natural {
        self.check(x);
        match &self.form {
            Form::Montgomery { .. } => {
                let k = self.len();
                let mut s = self.scratch();
                s.t[..k].copy_from_slice(&x.limbs);
                s.t[k..].fill(0);
                let mut r = vec![0; k];
                self.redc(&mut r, &mut s);
                Natural::from_vec(r)
            }
            Form::Plain => Natural::from_limbs(&x.limbs),
        }
    }

    /// Whether `x` is 0.
    #[must_use]
    pub fn is_zero(&self, x: &Residue) -> bool {
        self.check(x);
        x.limbs.iter().all(|&l| l == 0)
    }

    /// `x + y`.
    #[must_use]
    pub fn add(&self, x: &Residue, y: &Residue) -> Residue {
        self.check(x);
        self.check(y);
        let mut r = vec![0; self.len()];
        let carry = add_n(&mut r, &x.limbs, &y.limbs);
        if carry != 0 || cmp_n(&r, self.n.limbs()) != Ordering::Less {
            sub_n_assign(&mut r, self.n.limbs());
        }
        Self::residue(r.into_boxed_slice())
    }

    /// `x - y`.
    #[must_use]
    pub fn sub(&self, x: &Residue, y: &Residue) -> Residue {
        self.check(x);
        self.check(y);
        let mut r = vec![0; self.len()];
        if sub_n(&mut r, &x.limbs, &y.limbs) != 0 {
            add_n_assign(&mut r, self.n.limbs());
        }
        Self::residue(r.into_boxed_slice())
    }

    /// `-x`.
    #[must_use]
    pub fn neg(&self, x: &Residue) -> Residue {
        if self.is_zero(x) {
            return x.clone();
        }
        let mut r = vec![0; self.len()];
        sub_n(&mut r, self.n.limbs(), &x.limbs);
        Self::residue(r.into_boxed_slice())
    }

    /// `x / 2`, for an odd modulus.
    pub(crate) fn half(&self, x: &Residue) -> Residue {
        self.check(x);
        debug_assert!(self.n.bit(0), "halving needs an odd modulus");
        let mut r = x.limbs.to_vec();
        // x or x + n is even; halving is linear, so it commutes with the
        // Montgomery representation.
        let carry = if r[0] & 1 == 1 {
            add_n_assign(&mut r, self.n.limbs())
        } else {
            0
        };
        kernel::shr_assign(&mut r, 1);
        *r.last_mut().expect("non-empty") |= carry << (Limb::BITS - 1);
        Self::residue(r.into_boxed_slice())
    }

    /// `x y`.
    #[must_use]
    pub fn mul(&self, x: &Residue, y: &Residue) -> Residue {
        self.check(x);
        self.check(y);
        self.mul_limbs(&x.limbs, &y.limbs)
    }

    /// `x^2`.
    #[must_use]
    pub fn square(&self, x: &Residue) -> Residue {
        self.check(x);
        let mut s = self.scratch();
        let mut r = vec![0; self.len()];
        self.sqr_into(&mut r, &x.limbs, &mut s);
        Self::residue(r.into_boxed_slice())
    }

    fn mul_limbs(&self, x: &[Limb], y: &[Limb]) -> Residue {
        let mut s = self.scratch();
        let mut r = vec![0; self.len()];
        self.mul_into(&mut r, x, y, &mut s);
        Self::residue(r.into_boxed_slice())
    }

    /// `x^e`, with `x^0 = 1`.
    ///
    /// Left-to-right sliding-window exponentiation: the odd powers
    /// `x, x^3, …, x^(2^w - 1)` are precomputed, and each window of up to `w`
    /// exponent bits ending in a one costs one multiplication.
    #[must_use]
    pub fn pow(&self, x: &Residue, e: &Natural) -> Residue {
        self.check(x);
        let bits = e.bit_length();
        if bits == 0 {
            return self.one();
        }
        let k = self.len();
        let mut s = self.scratch();
        let w = window_bits(bits);
        // table[i] = x^(2i + 1).
        let mut table = vec![x.limbs.to_vec()];
        if w > 1 {
            let mut x2 = vec![0; k];
            self.sqr_into(&mut x2, &x.limbs, &mut s);
            for i in 1..1usize << (w - 1) {
                let mut next = vec![0; k];
                self.mul_into(&mut next, &table[i - 1], &x2, &mut s);
                table.push(next);
            }
        }
        // The window of bits i - 1 down to l, at most w long and ending in a
        // one, for a one bit i - 1; returns l and the window's value.
        let window = |i: u64| {
            let mut l = i.saturating_sub(w);
            while !e.bit(l) {
                l += 1;
            }
            let v = (l..i)
                .rev()
                .fold(0usize, |v, j| (v << 1) | usize::from(e.bit(j)));
            (l, v)
        };
        // The top bit is set, so the first window starts there.
        let (mut i, v) = window(bits);
        let mut r = table[v >> 1].clone();
        let mut tmp = vec![0; k];
        // Bit i - 1 is the next one to process.
        while i > 0 {
            if !e.bit(i - 1) {
                self.sqr_into(&mut tmp, &r, &mut s);
                std::mem::swap(&mut r, &mut tmp);
                i -= 1;
                continue;
            }
            let (l, v) = window(i);
            for _ in l..i {
                self.sqr_into(&mut tmp, &r, &mut s);
                std::mem::swap(&mut r, &mut tmp);
            }
            self.mul_into(&mut tmp, &r, &table[v >> 1], &mut s);
            std::mem::swap(&mut r, &mut tmp);
            i = l;
        }
        Self::residue(r.into_boxed_slice())
    }

    /// The inverse of `x`, or `None` if it is not a unit. Modulo 1 the
    /// inverse of 0 is 0.
    #[must_use]
    pub fn inv(&self, x: &Residue) -> Option<Residue> {
        let n = Integer::from(&self.n);
        let a = Integer::from(self.to_natural(x));
        let inv = a.mod_inverse(&n)?;
        Some(self.elem(&inv.unsigned_abs()))
    }

    fn scratch(&self) -> Scratch {
        let k = self.len();
        let (u, work) = match &self.form {
            Form::Montgomery { ninv_full, .. } => {
                let product = mul_scratch_len(k, k).max(sqr_scratch_len(k));
                if ninv_full.is_empty() {
                    (0, product)
                } else {
                    (4 * k, product.max(redc_n_scratch_len(k)))
                }
            }
            Form::Plain => (
                k + 1,
                mul_scratch_len(k, k)
                    .max(sqr_scratch_len(k))
                    .max(div_rem_scratch_len(2 * k, k)),
            ),
        };
        Scratch {
            t: vec![0; 2 * k],
            u: vec![0; u],
            work: vec![0; work],
        }
    }

    /// `r = x y` in the ring's representation.
    fn mul_into(&self, r: &mut [Limb], x: &[Limb], y: &[Limb], s: &mut Scratch) {
        kernel::mul(&mut s.t, x, y, &mut s.work);
        self.reduce_product(r, s);
    }

    /// `r = x^2` in the ring's representation.
    fn sqr_into(&self, r: &mut [Limb], x: &[Limb], s: &mut Scratch) {
        kernel::sqr(&mut s.t, x, &mut s.work);
        self.reduce_product(r, s);
    }

    /// Reduces the product in `s.t` into `r`.
    fn reduce_product(&self, r: &mut [Limb], s: &mut Scratch) {
        match &self.form {
            Form::Montgomery { .. } => self.redc(r, s),
            Form::Plain => div_rem(&mut s.u, r, &s.t, self.n.limbs(), &mut s.work),
        }
    }

    /// `r = t R^(-1) mod n` for `t = s.t < n R`, in `0..n`.
    fn redc(&self, r: &mut [Limb], s: &mut Scratch) {
        let Form::Montgomery {
            ninv,
            ninv2,
            ninv_full,
            ..
        } = &self.form
        else {
            unreachable!("REDC needs an odd modulus");
        };
        let n = self.n.limbs();
        let carry = if ninv_full.is_empty() {
            redc_2(r, &mut s.t, n, *ninv, *ninv2)
        } else {
            redc_n(r, &s.t, n, ninv_full, &mut s.u, &mut s.work)
        };
        if carry != 0 || cmp_n(r, n) != Ordering::Less {
            sub_n_assign(r, n);
        }
    }
}

/// `r = (t + q n) / R` for the `q < R` that makes it exact, eliminating two
/// limbs of `t` per pass with `q`'s next two limbs `t[i..i + 2] (-n^(-1))
/// mod B^2`, and a last single limb for odd `k`. Returns the carry out of
/// `r`; `r + carry B^k < 2n`. Clobbers `t`.
///
/// A pass clears `t[i]` and `t[i + 1]`; its two carry limbs belong at
/// positions `i + k` and `i + k + 1` and are kept in the cleared limbs until
/// the end, as no later `q` depends on them.
fn redc_2(r: &mut [Limb], t: &mut [Limb], n: &[Limb], ninv: Limb, ninv2: u128) -> Limb {
    let k = n.len();
    let mut i = 0;
    while i + 1 < k {
        let low = u128::from(t[i]) | u128::from(t[i + 1]) << Limb::BITS;
        let q = low.wrapping_mul(ninv2);
        #[allow(clippy::cast_possible_truncation, reason = "splitting a double limb")]
        let (q0, q1) = (q as Limb, (q >> Limb::BITS) as Limb);
        let (c0, c1) = addmul_2_carries(&mut t[i..i + k], n, q0, q1, 0);
        debug_assert!(t[i] == 0 && t[i + 1] == 0);
        (t[i], t[i + 1]) = (c0, c1);
        i += 2;
    }
    if i < k {
        let q = t[i].wrapping_mul(ninv);
        t[i] = kernel::addmul_1(&mut t[i..i + k], n, q);
    }
    let (lo, hi) = t.split_at(k);
    add_n(r, hi, lo)
}

/// `r = (t + q n) / R` as [`redc_2`], with `q = t ninv_full mod R` a short
/// product ([`mullo`]) and `q n` from a full product, or, when a cyclic
/// transform of `m = next_pow2(k)` coefficients applies, modulo `B^m -
/// 1` ([`wrapped_high`]). `u` holds `4k` limbs, `work` at least
/// [`redc_n_scratch_len`]`(k)`.
fn redc_n(
    r: &mut [Limb],
    t: &[Limb],
    n: &[Limb],
    ninv_full: &[Limb],
    u: &mut [Limb],
    work: &mut [Limb],
) -> Limb {
    let k = n.len();
    let (q, rest) = u.split_at_mut(k);
    mullo(q, &t[..k], ninv_full, work);
    // The low halves of t and q n add up to 0 mod R: to R exactly, with a
    // carry, unless both are zero.
    let low_carry = Limb::from(t[..k].iter().any(|&l| l != 0));
    let carry = if redc_wraps(k) {
        let (w, h) = rest.split_at_mut(k.next_power_of_two());
        wrapped_high(&mut h[..k], w, q, n, &t[..k], work);
        add_n(r, &t[k..], &h[..k])
    } else {
        let qn = &mut rest[..2 * k];
        kernel::mul(qn, q, n, work);
        add_n(r, &t[k..], &qn[k..])
    };
    carry + add_1_assign(r, low_carry)
}

/// Whether [`redc_n`] takes `q n` modulo `B^m - 1`, `m = next_pow2(k)`
/// (always below `2k`).
fn redc_wraps(k: usize) -> bool {
    wraps_cyclically(k.next_power_of_two())
}

/// Scratch limbs `work` that [`redc_n`] needs for a `k`-limb modulus.
fn redc_n_scratch_len(k: usize) -> usize {
    let product = if redc_wraps(k) {
        mul_wrap_scratch_len(k.next_power_of_two(), k, k)
    } else {
        mul_scratch_len(k, k)
    };
    product.max(mullo_scratch_len(k))
}

/// `h = floor(q n / B^k)` for `k = n.len() = q.len()`, from the product
/// modulo `B^m - 1` (`m = w.len()`, `k <= m < 2k`) and the low half `L = q n
/// mod B^k`, which REDC knows: `B^k - t_lo`, or 0.
///
/// Write `H = h`: `q n = L + H B^k`, and `H B^k = H_lo B^k + H_hi B^m` for
/// the low `m - k` limbs `H_lo` and the top `2k - m` limbs `H_hi` of `H`, so
/// `q n - L = H_lo B^k + H_hi` modulo `B^m - 1`. The two parts sit side by
/// side below `B^m`, and the sum is not `B^m - 1`: that needs `H_hi = B^k -
/// 1`, but `H_hi` has fewer than `k` limbs when `m > k`, and `H < n < B^k`
/// when `m = k`. So the canonical residue holds `H` exactly. Clobbers `w`; `h`
/// is the scratch for `L`.
fn wrapped_high(
    h: &mut [Limb],
    w: &mut [Limb],
    q: &[Limb],
    n: &[Limb],
    t_lo: &[Limb],
    work: &mut [Limb],
) {
    let (k, m) = (n.len(), w.len());
    debug_assert!(k <= m && m < 2 * k);
    mul_wrap(w, q, n, work);
    // L = B^k - t_lo, or 0.
    for (l, &x) in h.iter_mut().zip(t_lo) {
        *l = !x;
    }
    add_1_assign(h, 1);
    if sub_assign(w, h) != 0 {
        // The borrow added B^m; the residue needs B^m - 1.
        sub_1_assign(w, 1);
    }
    debug_assert!(
        w[2 * k - m..k].iter().all(|&x| x == 0),
        "redc_n: H_lo and H_hi overlap"
    );
    h[..m - k].copy_from_slice(&w[k..]);
    h[m - k..].copy_from_slice(&w[..2 * k - m]);
}

/// The window width for an exponent of `bits` bits, minimising the
/// multiplications: about `bits / (w + 1)` windows plus a table of `2^(w-1)`
/// odd powers, while the squarings are `bits` for every width. Width `w + 1`
/// wins once `bits > 2^(w-1) (w + 1) (w + 2)`: from 7, 25, 81, 241, 673,
/// 1793, 4609 bits on.
pub(super) fn window_bits(bits: u64) -> u64 {
    let mut w = 1;
    while w < 16 && bits > (1 << (w - 1)) * (w + 1) * (w + 2) {
        w += 1;
    }
    w
}

/// `a^(-1) mod B` for odd `a`, by Newton's iteration: each step doubles the
/// number of correct low bits, and `a` is its own inverse modulo 8.
pub(super) fn limb_inverse(a: Limb) -> Limb {
    debug_assert_eq!(a & 1, 1);
    let mut x = a;
    for _ in 0..5 {
        x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
    }
    x
}

/// The limbs of `a < B^k`, zero-padded to `k`.
fn padded(a: &Natural, k: usize) -> Box<[Limb]> {
    let mut v = a.limbs().to_vec();
    debug_assert!(v.len() <= k);
    v.resize(k, 0);
    v.into_boxed_slice()
}
