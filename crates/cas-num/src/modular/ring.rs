//! The ring ℤ/nℤ for a multi-precision modulus.
//!
//! Odd moduli use Montgomery's representation: a residue `x` is stored as
//! `x R mod n` with `R = B^k` for a `k`-limb modulus, and products are
//! reduced by REDC, which divides by `R` instead of by `n`. Below
//! [`REDC_N_THRESHOLD`] limbs REDC eliminates two limbs per pass over `n`
//! (GMP's `redc_2`). From it on (`redc_n`), it takes `q = t (-n^(-1)) mod R`
//! by a short product and `q n` in full or modulo `B^m - 1`, from which the
//! known low half recovers the high one. Below the cyclic transform that product
//! splits `B^m - 1` ([`crate::kernel::mul::split`]); at transform sizes, `n`
//! and `-n^(-1)` are transformed once per ring, so that each product
//! transforms only `q` (or `t`'s low half) and the result. Even moduli store
//! plain values and reduce products by division.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (n, k, r, t, q, u)"
)]

use std::cell::RefCell;
use std::cmp::Ordering;

use crate::kernel::div::bdiv::{binvert, binvert_scratch_len};
use crate::kernel::div::{div_2by1, reciprocal_2by1};
use crate::kernel::mul::mullo::{cost as mullo_cost, mullo, mullo_scratch_len};
use crate::kernel::mul::split::{split_wrap_len, splits};
use crate::kernel::mul::{
    addmul_2_carries, mul_wrap, mul_wrap_scratch_len, wrap_len, wraps_cyclically,
};
use crate::kernel::{
    self, add_1_assign, add_n, add_n_assign, addmul_1, cmp_n, div_rem, div_rem_scratch_len,
    mul_scratch_len, ntt, sqr_scratch_len, sub_1_assign, sub_assign, sub_n, sub_n_assign, submul_1,
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
#[derive(PartialEq, Eq, Hash, Debug)]
pub struct Residue {
    /// Exactly as many limbs as the modulus, below it.
    limbs: Box<[Limb]>,
}

impl Clone for Residue {
    fn clone(&self) -> Self {
        Residue {
            limbs: self.limbs.clone(),
        }
    }

    /// Reuses the buffer for a residue of the same ring.
    fn clone_from(&mut self, source: &Self) {
        if self.limbs.len() == source.limbs.len() {
            self.limbs.copy_from_slice(&source.limbs);
        } else {
            *self = source.clone();
        }
    }
}

#[derive(Clone, Debug)]
enum Form {
    /// Odd modulus: residues are `x R mod n`.
    Montgomery {
        /// `-n^(-1) mod B`.
        ninv: Limb,
        /// `-n^(-1) mod B^2`, for `redc_2`.
        ninv2: u128,
        /// What `redc_n` needs, from `REDC_N_THRESHOLD` limbs on.
        redc: Option<Redc>,
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
                redc: (k >= REDC_N_THRESHOLD).then(|| Redc::new(n.limbs())),
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
        let mut r = self.zero();
        self.add_into(&mut r.limbs, x, y);
        r
    }

    /// `x - y`.
    #[must_use]
    pub fn sub(&self, x: &Residue, y: &Residue) -> Residue {
        let mut r = self.zero();
        self.sub_into(&mut r.limbs, x, y);
        r
    }

    /// `-x`.
    #[must_use]
    pub fn neg(&self, x: &Residue) -> Residue {
        let mut r = self.zero();
        self.neg_into(&mut r.limbs, x);
        r
    }

    fn add_into(&self, r: &mut [Limb], x: &Residue, y: &Residue) {
        self.check(x);
        self.check(y);
        let carry = add_n(r, &x.limbs, &y.limbs);
        if carry != 0 || cmp_n(r, self.n.limbs()) != Ordering::Less {
            sub_n_assign(r, self.n.limbs());
        }
    }

    fn sub_into(&self, r: &mut [Limb], x: &Residue, y: &Residue) {
        self.check(x);
        self.check(y);
        if sub_n(r, &x.limbs, &y.limbs) != 0 {
            add_n_assign(r, self.n.limbs());
        }
    }

    fn neg_into(&self, r: &mut [Limb], x: &Residue) {
        if self.is_zero(x) {
            r.fill(0);
        } else {
            sub_n(r, self.n.limbs(), &x.limbs);
        }
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
            Form::Montgomery { redc, .. } => {
                let product = mul_scratch_len(k, k).max(sqr_scratch_len(k));
                match redc {
                    None => (0, product),
                    Some(redc) => (4 * k, product.max(redc.scratch_len(k))),
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
        let carry = self.redc_carry(r, s);
        if carry != 0 || cmp_n(r, self.n.limbs()) != Ordering::Less {
            sub_n_assign(r, self.n.limbs());
        }
    }

    /// `r + carry B^k = (t + q n) / R` for `t = s.t` (any `2k` limbs) and
    /// the `q < R` that makes it exact: congruent to `t R^(-1)`, and below
    /// `B^k + n`. Clobbers `s.t`.
    fn redc_carry(&self, r: &mut [Limb], s: &mut Scratch) -> Limb {
        let Form::Montgomery {
            ninv, ninv2, redc, ..
        } = &self.form
        else {
            unreachable!("REDC needs an odd modulus");
        };
        let n = self.n.limbs();
        match redc {
            None => redc_2(r, &mut s.t, n, *ninv, *ninv2),
            Some(redc) => redc_n(r, &s.t, n, redc, &mut s.u, &mut s.work),
        }
    }
}

/// A [`ModRing`] with its product buffers, for loops of steps that write
/// into existing residues and so allocate nothing (the probable-prime
/// tests). Each method's output is a residue of the ring, overwritten.
pub(crate) struct RingWork<'a> {
    ring: &'a ModRing,
    s: RefCell<Scratch>,
    /// Four products of `2k + 1` limbs and a `k`-limb difference, for
    /// [`RingWork::lucas_double`].
    wide: RefCell<Vec<Limb>>,
    /// `n`'s top limb after shifting `n` left by `shift` to set its top bit,
    /// and that limb's reciprocal, for [`RingWork::reduce_top`].
    shift: u32,
    top: Limb,
    recip: Limb,
}

impl<'a> RingWork<'a> {
    pub(crate) fn new(ring: &'a ModRing) -> Self {
        let n = ring.n.limbs();
        let k = n.len();
        let shift = n[k - 1].leading_zeros();
        let below = if k >= 2 { n[k - 2] } else { 0 };
        let top = if shift == 0 {
            n[k - 1]
        } else {
            n[k - 1] << shift | below >> (Limb::BITS - shift)
        };
        RingWork {
            ring,
            s: RefCell::new(ring.scratch()),
            wide: RefCell::new(vec![0; 4 * (2 * k + 1) + k]),
            shift,
            top,
            recip: reciprocal_2by1(top),
        }
    }

    /// The ring.
    pub(crate) fn ring(&self) -> &'a ModRing {
        self.ring
    }

    /// `r = x + y`.
    pub(crate) fn add(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        self.ring.add_into(&mut r.limbs, x, y);
    }

    /// `r = x - y`.
    pub(crate) fn sub(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        self.ring.sub_into(&mut r.limbs, x, y);
    }

    /// `r = -x`.
    pub(crate) fn neg(&self, r: &mut Residue, x: &Residue) {
        self.ring.neg_into(&mut r.limbs, x);
    }

    /// `r = x y`.
    pub(crate) fn mul(&self, r: &mut Residue, x: &Residue, y: &Residue) {
        self.ring.check(x);
        self.ring.check(y);
        let s = &mut self.s.borrow_mut();
        self.ring.mul_into(&mut r.limbs, &x.limbs, &y.limbs, s);
    }

    /// `r = x^2`.
    pub(crate) fn square(&self, r: &mut Residue, x: &Residue) {
        self.ring.check(x);
        let s = &mut self.s.borrow_mut();
        self.ring.sqr_into(&mut r.limbs, &x.limbs, s);
    }

    /// One doubling step of the Lucas sequence `U` for `P = 1` and `Q = q`,
    /// `|q| < 2^8`, modulo the odd modulus: `(u, u1) = (U_j, U_(j+1))`
    /// becomes `(U_(2j), U_(2j+1))`, or `(U_(2j+1), U_(2j+2))` if `up`.
    ///
    /// With `A = U_j^2`, `B = U_(j+1)^2` and `C = (U_j - U_(j+1))^2`,
    /// `U_(2j) = B - C`, `U_(2j+1) = B - Q A` and `U_(2j+2) = U_(2j+1) -
    /// Q U_(2j)`. The three squares stay unreduced, and only the two
    /// combinations needed are reduced: two REDCs per step instead of
    /// three. Each combination gets a multiple of `n R` that keeps it
    /// non-negative, which REDC maps to a multiple of `n`.
    pub(crate) fn lucas_double(&self, u: &mut Residue, u1: &mut Residue, q: i64, up: bool) {
        self.ring.check(u);
        self.ring.check(u1);
        let n = self.ring.n.limbs();
        let k = n.len();
        let m = q.unsigned_abs();
        debug_assert!(m < 1 << 8, "lucas_double: |Q| too large");
        let s = &mut *self.s.borrow_mut();
        let wide = &mut *self.wide.borrow_mut();
        let (a, rest) = wide.split_at_mut(2 * k + 1);
        let (b, rest) = rest.split_at_mut(2 * k + 1);
        let (x, rest) = rest.split_at_mut(2 * k + 1);
        let (c, d) = rest.split_at_mut(2 * k + 1);
        kernel::sqr(&mut a[..2 * k], &u.limbs, &mut s.work);
        kernel::sqr(&mut b[..2 * k], &u1.limbs, &mut s.work);
        // |U_j - U_(j+1)|, whose square is C.
        if cmp_n(&u.limbs, &u1.limbs) == Ordering::Less {
            sub_n(d, &u1.limbs, &u.limbs);
        } else {
            sub_n(d, &u.limbs, &u1.limbs);
        }
        kernel::sqr(&mut c[..2 * k], d, &mut s.work);
        // x = B - C + n R, in (0, 2 n R).
        x[..2 * k].copy_from_slice(&b[..2 * k]);
        let carry = add_n_assign(&mut x[k..2 * k], n);
        let borrow = sub_n_assign(&mut x[..2 * k], &c[..2 * k]);
        x[2 * k] = carry - borrow;
        // b = B - Q A, plus Q n R for Q > 0: in [0, (1 + |Q|) n R).
        b[2 * k] = if q < 0 {
            addmul_1(&mut b[..2 * k], &a[..2 * k], m)
        } else {
            let carry = addmul_1(&mut b[k..2 * k], n, m);
            carry - submul_1(&mut b[..2 * k], &a[..2 * k], m)
        };
        if up {
            // c = b - Q x, plus 2 Q n R for Q > 0: in [0, (1 + 3|Q|) n R).
            c.copy_from_slice(b);
            if q < 0 {
                let carry = addmul_1(c, x, m);
                debug_assert_eq!(carry, 0);
            } else {
                c[2 * k] += addmul_1(&mut c[k..2 * k], n, 2 * m);
                let borrow = submul_1(c, x, m);
                debug_assert_eq!(borrow, 0);
            }
            self.redc_wide(&mut u.limbs, b, s);
            self.redc_wide(&mut u1.limbs, c, s);
        } else {
            self.redc_wide(&mut u.limbs, x, s);
            self.redc_wide(&mut u1.limbs, b, s);
        }
    }

    /// `r = w R^(-1) mod n` for `w` of `2k + 1` limbs, below `2^61 n R`:
    /// `(w + q n) / R` is below `2^62 n`.
    fn redc_wide(&self, r: &mut [Limb], w: &[Limb], s: &mut Scratch) {
        let k = r.len();
        s.t.copy_from_slice(&w[..2 * k]);
        let carry = self.ring.redc_carry(r, s);
        self.reduce_top(r, carry + w[2 * k]);
    }

    /// `r = (r + top B^k) mod n` for a value `V` below `2^62 n`, by one
    /// quotient estimate. With `n` shifted left to set its top bit, `N_top`
    /// its top limb, and `v` the two limbs of `V` shifted alike at the same
    /// positions, `q_est = floor(v / N_top)` is at least `q = floor(V / n)`
    /// (`v >= q N_top`), and at most one more: `v < (q + 1) (N_top + 1)`
    /// gives `q_est <= q + 1 + floor((q + 1) / N_top)`, and `q + 1 < 2^63
    /// <= N_top`. So `V - q_est n` needs `n` added back at most once.
    fn reduce_top(&self, r: &mut [Limb], mut top: Limb) {
        let n = self.ring.n.limbs();
        let k = n.len();
        if top == 0 && cmp_n(r, n) == Ordering::Less {
            return;
        }
        let r1 = r[k - 1];
        let r2 = if k >= 2 { r[k - 2] } else { 0 };
        let (v1, v0) = if self.shift == 0 {
            (top, r1)
        } else {
            let back = Limb::BITS - self.shift;
            (
                top << self.shift | r1 >> back,
                r1 << self.shift | r2 >> back,
            )
        };
        let (q, _) = div_2by1(v1, v0, self.top, self.recip);
        top = top.wrapping_sub(submul_1(r, n, q));
        if top != 0 {
            top = top.wrapping_add(add_n_assign(r, n));
            debug_assert_eq!(top, 0, "reduce_top: quotient estimate off by more than one");
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

/// What [`redc_n`] needs besides the modulus, chosen once per ring:
/// `-n^(-1) mod R` and how to take its two products.
#[derive(Clone, Debug)]
struct Redc {
    /// `-n^(-1) mod R`.
    ninv: Vec<Limb>,
    /// The modulus length `m` of `q n` modulo `B^m - 1` ([`redc_wrap_len`]),
    /// or 0 for the full product.
    wrap: usize,
    /// `n`'s transform at the cyclic length `wrap` ([`prepare_cyclic`]),
    /// when `wrap` takes the cyclic transform; else empty.
    n_prep: Vec<Limb>,
    /// A cyclic length `m >= 2k` at which `ninv` is prepared in
    /// `ninv_prep`, so that `q = t ninv mod R` is the low half of a cyclic
    /// product that does not wrap; or 0 for [`mullo`].
    ninv_len: usize,
    ninv_prep: Vec<Limb>,
}

impl Redc {
    /// The choice for an odd `k`-limb modulus `n`. `q n` is taken modulo
    /// `B^m - 1` as [`redc_wrap_len`] says, with `n` prepared when that is
    /// cyclic. `q` comes from `ninv` prepared at the shortest transform
    /// length from `2k` when that product's two transforms cost less than
    /// [`mullo`]. Each prepared product transforms only the other operand
    /// and the result, two transforms instead of three.
    #[allow(clippy::cast_precision_loss, reason = "a cost estimate")]
    fn new(n: &[Limb]) -> Redc {
        let k = n.len();
        let m = wrap_len(2 * k);
        let prepared = ntt::cyclic_cost(m) as f64 * 2.0 / 3.0;
        let ninv_len = if wraps_cyclically(m) && prepared < mullo_cost(k) {
            m
        } else {
            0
        };
        Redc::with(n, redc_wrap_len(k), ninv_len)
    }

    /// The products for `n` with the given `wrap` and `ninv_len`; `n` is
    /// prepared whenever `wrap` takes the cyclic transform.
    fn with(n: &[Limb], wrap: usize, ninv_len: usize) -> Redc {
        let ninv = neg_inverse(n);
        let prepare = |m: usize, x: &[Limb]| {
            let mut p = vec![0; ntt::prepared_len(m)];
            let mut scratch = vec![0; ntt::cyclic_prepared_scratch_len(m)];
            ntt::prepare_cyclic(&mut p, m, x, &mut scratch);
            p
        };
        let n_prep = if wrap != 0 && wraps_cyclically(wrap) {
            prepare(wrap, n)
        } else {
            Vec::new()
        };
        let ninv_prep = if ninv_len == 0 {
            Vec::new()
        } else {
            prepare(ninv_len, &ninv)
        };
        Redc {
            ninv,
            wrap,
            n_prep,
            ninv_len,
            ninv_prep,
        }
    }

    /// Scratch limbs `work` that [`redc_n`] needs for a `k`-limb modulus.
    fn scratch_len(&self, k: usize) -> usize {
        let q = if self.ninv_len == 0 {
            mullo_scratch_len(k)
        } else {
            ntt::cyclic_prepared_scratch_len(self.ninv_len)
        };
        let qn = if self.wrap == 0 {
            mul_scratch_len(k, k)
        } else if self.n_prep.is_empty() {
            mul_wrap_scratch_len(self.wrap, k, k)
        } else {
            ntt::cyclic_prepared_scratch_len(self.wrap)
        };
        q.max(qn)
    }
}

/// `r = (t + q n) / R` as [`redc_2`], with `q = t ninv mod R` a short product
/// ([`mullo`], or the low half of a prepared cyclic product) and `q n` from a
/// full product or modulo `B^m - 1` ([`wrapped_high`]), as `redc` says. `u`
/// holds `4k` limbs, `work` at least `redc.scratch_len(k)`.
fn redc_n(
    r: &mut [Limb],
    t: &[Limb],
    n: &[Limb],
    redc: &Redc,
    u: &mut [Limb],
    work: &mut [Limb],
) -> Limb {
    let k = n.len();
    let (q, rest) = u.split_at_mut(k);
    if redc.ninv_len == 0 {
        mullo(q, &t[..k], &redc.ninv, work);
    } else {
        // t_lo ninv < B^2k <= B^m does not wrap; ninv_len < 3k.
        let p = &mut rest[..redc.ninv_len];
        ntt::mul_cyclic_prepared(p, &t[..k], &redc.ninv_prep, work);
        q.copy_from_slice(&p[..k]);
    }
    // The low halves of t and q n add up to 0 mod R: to R exactly, with a
    // carry, unless both are zero.
    let low_carry = Limb::from(t[..k].iter().any(|&l| l != 0));
    let carry = if redc.wrap == 0 {
        let qn = &mut rest[..2 * k];
        kernel::mul(qn, q, n, work);
        add_n(r, &t[k..], &qn[k..])
    } else {
        let (w, h) = rest.split_at_mut(redc.wrap);
        if redc.n_prep.is_empty() {
            mul_wrap(w, q, n, work);
        } else {
            ntt::mul_cyclic_prepared(w, q, &redc.n_prep, work);
        }
        wrapped_high(&mut h[..k], w, &t[..k]);
        add_n(r, &t[k..], &h[..k])
    };
    carry + add_1_assign(r, low_carry)
}

/// The modulus length `m` with which [`redc_n`] takes `q n` modulo `B^m - 1`
/// for a `k`-limb modulus, or 0 for the full product: the shortest transform
/// length from `k` on when it takes the cyclic transform, else `k` rounded
/// up for [`splits`] when the split applies. Both are below `2k`.
fn redc_wrap_len(k: usize) -> usize {
    let m = wrap_len(k);
    if wraps_cyclically(m) {
        return m;
    }
    let m = split_wrap_len(k);
    if splits(m, k, k) { m } else { 0 }
}

/// `h = floor(q n / B^k)` for `k = h.len()`, from a residue `w` of `q n`
/// modulo `B^m - 1` (`m = w.len()`, `k <= m < 2k`) and the low half `L = q n
/// mod B^k`, which REDC knows: `B^k - t_lo`, or 0. `w` may be `B^m - 1` for
/// zero, as a cyclic transform leaves it: with `L > 0` both give `B^m - 1 -
/// L` below, and `L = 0` means `t_lo = 0`, so `q = 0`, whose product is
/// exactly zero.
///
/// Write `H = h`: `q n = L + H B^k`, and `H B^k = H_lo B^k + H_hi B^m` for
/// the low `m - k` limbs `H_lo` and the top `2k - m` limbs `H_hi` of `H`, so
/// `q n - L = H_lo B^k + H_hi` modulo `B^m - 1`. The two parts sit side by
/// side below `B^m`, and the sum is not `B^m - 1`: that needs `H_hi = B^k -
/// 1`, but `H_hi` has fewer than `k` limbs when `m > k`, and `H < n < B^k`
/// when `m = k`. So the canonical residue holds `H` exactly. Clobbers `w`; `h`
/// is the scratch for `L`.
fn wrapped_high(h: &mut [Limb], w: &mut [Limb], t_lo: &[Limb]) {
    let (k, m) = (h.len(), w.len());
    debug_assert!(k <= m && m < 2 * k);
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

/// `-n^(-1) mod B^k` for an odd `k`-limb `n`, by the 2-adic Newton iteration
/// ([`binvert`]) rather than an extended GCD.
fn neg_inverse(n: &[Limb]) -> Vec<Limb> {
    let k = n.len();
    let mut x = vec![0; k];
    let mut scratch = vec![0; binvert_scratch_len(k)];
    binvert(&mut x, n, &mut scratch);
    // B^k - x, as x is odd and so non-zero.
    for l in &mut x {
        *l = !*l;
    }
    add_1_assign(&mut x, 1);
    x
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
