//! Division and square root: correctly rounded, and approximate entries with
//! an error bound for ball midpoints (CAS-DOC-5 §3.4, §5.2).
#![allow(
    clippy::many_single_char_names,
    reason = "a, b, q, r, d, e, n, s: the operands and textbook names"
)]

use std::cmp::Ordering;

use super::{Float, LIMB_BITS, Round, with_limbs};
use crate::kernel;
use crate::{Integer, Limb, Natural};

/// The error bound, in units in the last place at the requested precision,
/// of [`Float::div_approx`], [`Float::sqrt_approx`] and
/// [`Float::rsqrt_approx`]: truncation loses below one unit, the kernels'
/// error (a few units of a limb at least 63 bits further down) below another.
pub const APPROX_ULPS: u64 = 2;

impl Float {
    /// `self / rhs` correctly rounded to `prec` bits, and `rounded.cmp(exact)`;
    /// `None` if `rhs` is zero.
    ///
    /// The mantissa of `self` gets zero limbs below it until the exact integer
    /// quotient has at least `prec + 2` bits; the remainder is the sticky
    /// bit.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn checked_div_round(
        &self,
        rhs: &Float,
        prec: u64,
        rnd: Round,
    ) -> Option<(Float, Ordering)> {
        assert!(prec >= 1, "Float::div_round: precision 0");
        if rhs.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some((Float::ZERO, Ordering::Equal));
        }
        let (a, b) = (self.mant(), rhs.mant());
        let (ka, kb) = (a.len(), b.len());
        // q >= 2^(64 (ka + t - kb) - 1) has 64 (ka + t - kb) bits or more.
        let need = usize::try_from(prec.div_ceil(LIMB_BITS)).expect("precision fits") + 1;
        let t = (need + kb).saturating_sub(ka);
        let na = ka + t;
        let nq = na - kb + 1;
        let scratch = kernel::div_rem_scratch_len(na, kb);
        let neg = self.neg != rhs.neg;
        // a / b = (A B^t / D) 2^(ea - eb - 64 (ka + t - kb)).
        let e = &self.exp - &rhs.exp - Integer::from(LIMB_BITS * (na - kb) as u64);
        Some(with_limbs(na + nq + kb + scratch, |buf| {
            let (num, rest) = buf.split_at_mut(na);
            let (q, rest) = rest.split_at_mut(nq);
            let (r, scratch) = rest.split_at_mut(kb);
            num[t..].copy_from_slice(a);
            kernel::div_rem(q, r, num, b, scratch);
            let sticky = r.iter().any(|&x| x != 0);
            Float::build(neg, q, &e, sticky, Some((prec, rnd)))
        }))
    }

    /// `self / rhs` correctly rounded to `prec` bits, and `rounded.cmp(exact)`.
    ///
    /// # Panics
    ///
    /// If `rhs` is zero or `prec` is 0.
    #[must_use]
    pub fn div_round(&self, rhs: &Float, prec: u64, rnd: Round) -> (Float, Ordering) {
        self.checked_div_round(rhs, prec, rnd)
            .expect("attempt to divide by zero")
    }

    /// `sqrt(self)` correctly rounded to `prec` bits, and `rounded.cmp(exact)`;
    /// `None` if `self` is negative.
    ///
    /// The mantissa is scaled by an even power of two to at least
    /// `2 (prec + 2)` bits, so its integer square root has `prec + 2` bits or
    /// more; the remainder is the sticky bit.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn checked_sqrt_round(&self, prec: u64, rnd: Round) -> Option<(Float, Ordering)> {
        assert!(prec >= 1, "Float::sqrt_round: precision 0");
        if self.neg {
            return None;
        }
        if self.is_zero() {
            return Some((Float::ZERO, Ordering::Equal));
        }
        // self = M 2^E with E = exp - 64 k; shift M left by s, with E - s even.
        let (m, e) = self.to_integer_exp();
        let bits = self.mant_bits_stored();
        let mut s = (2 * (prec + 2)).saturating_sub(bits);
        if is_odd(&e) != (s % 2 == 1) {
            s += 1;
        }
        let n = Natural::try_from(&m).expect("non-negative") << s;
        let (root, rem) = n.sqrt_rem();
        // sqrt(self) = sqrt(n) 2^((E - s) / 2).
        let half = (&e - &Integer::from(s)).div_exact(&Integer::from(2));
        let mut limbs = root.limbs().to_vec();
        Some(Float::build(
            false,
            &mut limbs,
            &half,
            !rem.is_zero(),
            Some((prec, rnd)),
        ))
    }

    /// `sqrt(self)` correctly rounded to `prec` bits, and `rounded.cmp(exact)`.
    ///
    /// # Panics
    ///
    /// If `self` is negative or `prec` is 0.
    #[must_use]
    pub fn sqrt_round(&self, prec: u64, rnd: Round) -> (Float, Ordering) {
        self.checked_sqrt_round(prec, rnd)
            .expect("square root of a negative number")
    }

    /// `q` near `self / rhs` with `|q| <= |self / rhs|` and
    /// `|self / rhs - q| < APPROX_ULPS` units of `2^(exp(q) - prec)`, `q` of at
    /// most `prec` bits; `None` if `rhs` is zero. For ball midpoints.
    ///
    /// Where the kernel would divide exactly (below its Newton threshold),
    /// this is [`Float::div_round`] toward zero, within one ulp. Otherwise
    /// both mantissas are padded to `n = max(ka, kb, ceil(prec / 64) + 1)`
    /// limbs and divided by `kernel::div_approx`, whose quotient `Q` is a lower
    /// bound within `e = div_approx_error(n)` of `floor(A B^n / D)`, in units
    /// of `2^(ea - eb - 64 n)`. As `Q >= B^n / 2 - e`, `exp(q) >= ea - eb - 1`,
    /// so a unit is at most `2^(prec + 1 - 64 n) <= 2^-63` ulps: the kernel's
    /// error and the floor stay below one ulp, the final truncation below
    /// another.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn div_approx(&self, rhs: &Float, prec: u64) -> Option<Float> {
        assert!(prec >= 1, "Float::div_approx: precision 0");
        if rhs.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some(Float::ZERO);
        }
        let n = approx_len(prec, self.mant().len().max(rhs.mant().len()));
        if kernel::div_approx_error(n) == 0 {
            // The kernel divides exactly at this size: the correctly rounded
            // quotient, truncated, is cheaper and within one ulp.
            return Some(self.div_round(rhs, prec, Round::Down).0);
        }
        let neg = self.neg != rhs.neg;
        let e = &self.exp - &rhs.exp - Integer::from(LIMB_BITS * n as u64);
        Some(with_limbs(
            3 * n + 1 + kernel::div_approx_scratch_len(n),
            |buf| {
                let (a, rest) = buf.split_at_mut(n);
                let (d, rest) = rest.split_at_mut(n);
                let (q, scratch) = rest.split_at_mut(n + 1);
                pad(a, self.mant());
                pad(d, rhs.mant());
                kernel::div_approx(q, a, d, scratch);
                Float::build(neg, q, &e, false, Some((prec, Round::Down))).0
            },
        ))
    }

    /// `s` near `sqrt(self)` with `s <= sqrt(self)` and `sqrt(self) - s <
    /// APPROX_ULPS` units of `2^(exp(s) - prec)`, `s` of at most `prec` bits;
    /// `None` if `self` is negative. For ball midpoints.
    ///
    /// Where the kernel's root is exact, this is [`Float::sqrt_round`] toward
    /// zero. Otherwise the mantissa, padded to `n` limbs as in [`Float::div_approx`] (one more,
    /// halved, when the exponent is odd, so that it stays even and at least
    /// `B^n / 4`), goes to `kernel::sqrt_approx`; `S >= B^n / 2 - e` gives the
    /// same accounting.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn sqrt_approx(&self, prec: u64) -> Option<Float> {
        assert!(prec >= 1, "Float::sqrt_approx: precision 0");
        if self.neg {
            return None;
        }
        if self.is_zero() {
            return Some(Float::ZERO);
        }
        let (n, half) = self.sqrt_layout(prec);
        if kernel::sqrt_approx_error(n) == 0 {
            // Exact at this size, as for division.
            return Some(self.sqrt_round(prec, Round::Down).0);
        }
        // sqrt(A / B^n 2^(2 half)) = sqrt(A B^n) / B^n 2^half.
        let e = half - Integer::from(LIMB_BITS * n as u64);
        Some(with_limbs(
            2 * n + kernel::sqrt_approx_scratch_len(n),
            |buf| {
                let (a, rest) = buf.split_at_mut(n);
                let (s, scratch) = rest.split_at_mut(n);
                self.sqrt_mantissa(a);
                kernel::sqrt_approx(s, a, scratch);
                Float::build(false, s, &e, false, Some((prec, Round::Down))).0
            },
        ))
    }

    /// `y` near `1 / sqrt(self)` with `y <= 1 / sqrt(self)` and
    /// `1 / sqrt(self) - y < APPROX_ULPS` units of `2^(exp(y) - prec)`, `y` of
    /// at most `prec` bits; `None` if `self` is zero or negative. For ball
    /// midpoints.
    ///
    /// As [`Float::sqrt_approx`], through `kernel::rsqrt_approx`:
    /// `1 / sqrt(A / B^n 2^(2 half)) = Y* / B^n 2^-half` with `Y* =
    /// sqrt(B^(3n) / A) >= B^n`, so a unit of `Y` is at most `2^-63` ulps.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn rsqrt_approx(&self, prec: u64) -> Option<Float> {
        assert!(prec >= 1, "Float::rsqrt_approx: precision 0");
        if self.neg || self.is_zero() {
            return None;
        }
        let (n, half) = self.sqrt_layout(prec);
        let e = -half - Integer::from(LIMB_BITS * n as u64);
        Some(with_limbs(
            2 * n + 1 + kernel::rsqrt_approx_scratch_len(n),
            |buf| {
                let (a, rest) = buf.split_at_mut(n);
                let (y, scratch) = rest.split_at_mut(n + 1);
                self.sqrt_mantissa(a);
                kernel::rsqrt_approx(y, a, scratch);
                Float::build(false, y, &e, false, Some((prec, Round::Down))).0
            },
        ))
    }

    /// A ball midpoint for `self / rhs`: `(q, k)` with `|q| <= |self / rhs|`
    /// and `|self / rhs| - |q| < k` units of `2^(exp(q) - prec)`, `k = 0`
    /// exactly when `q` is the exact quotient (known where the quotient is
    /// computed exactly, below the kernel's Newton threshold); `None` for a
    /// zero divisor.
    pub(crate) fn div_mid(&self, rhs: &Float, prec: u64) -> Option<(Float, u64)> {
        let n = approx_len(prec, self.mant().len().max(rhs.mant().len()));
        if kernel::div_approx_error(n) == 0 {
            let (q, o) = self.checked_div_round(rhs, prec, Round::Down)?;
            return Some((q, u64::from(o != Ordering::Equal)));
        }
        let q = self.div_approx(rhs, prec)?;
        // An exact quotient that fits `prec` bits is q, q + u or q + 2u
        // (the error is below 2 ulps); c rhs = self only for that one. The
        // product of odd mantissas of b and c bits has b + c - 1 or b + c
        // bits, which rules out most candidates without multiplying.
        let (bits_a, bits_b) = (self.bit_length(), rhs.bit_length());
        let exact = candidates(&q, prec).find(|c| {
            let bits_c = c.bit_length();
            bits_a + 1 >= bits_c + bits_b && bits_a <= bits_c + bits_b && c.mul(rhs) == *self
        });
        Some(match exact {
            Some(c) => (c, 0),
            None => (q, APPROX_ULPS),
        })
    }

    /// A ball midpoint for `sqrt(self)`, as [`Float::div_mid`]; `None` for a
    /// negative argument.
    pub(crate) fn sqrt_mid(&self, prec: u64) -> Option<(Float, u64)> {
        let (n, _) = if self.is_zero() || self.neg {
            (1, Integer::ZERO)
        } else {
            self.sqrt_layout(prec)
        };
        if kernel::sqrt_approx_error(n) == 0 {
            let (s, o) = self.checked_sqrt_round(prec, Round::Down)?;
            return Some((s, u64::from(o != Ordering::Equal)));
        }
        let s = self.sqrt_approx(prec)?;
        // As for the quotient: an exact root is s, s + u or s + 2u, and c^2
        // has 2c - 1 or 2c bits.
        let bits_a = self.bit_length();
        let exact = candidates(&s, prec).find(|c| {
            let bits_c = c.bit_length();
            bits_a + 1 >= 2 * bits_c && bits_a <= 2 * bits_c && c.mul(c) == *self
        });
        Some(match exact {
            Some(c) => (c, 0),
            None => (s, APPROX_ULPS),
        })
    }

    /// `(n, half)` with `self = A / B^n 2^(2 half)` for the mantissa `A` of
    /// [`Float::sqrt_mantissa`] (`n` limbs, `A >= B^n / 4`).
    fn sqrt_layout(&self, prec: u64) -> (usize, Integer) {
        let n = approx_len(prec, self.mant().len());
        if is_odd(&self.exp) {
            (
                n + 1,
                (&self.exp + &Integer::ONE).div_exact(&Integer::from(2)),
            )
        } else {
            (n, self.exp.div_exact(&Integer::from(2)))
        }
    }

    /// Writes the mantissa padded to `a.len()` limbs, halved (with the extra
    /// limb of [`Float::sqrt_layout`] below it) when the exponent is odd.
    fn sqrt_mantissa(&self, a: &mut [Limb]) {
        pad(a, self.mant());
        if is_odd(&self.exp) {
            let out = kernel::shr_assign(a, 1);
            debug_assert_eq!(out, 0, "the extra low limb takes the bit");
        }
    }
}

/// `q`, `q + u`, `q + 2u` for `u` one ulp of `q` at `prec` bits (away from
/// zero), those of at most `prec` bits: where an exactly representable
/// result lies when `q` is a lower bound within 2 ulps.
fn candidates(q: &Float, prec: u64) -> impl Iterator<Item = Float> + '_ {
    let unit = (!q.is_zero()).then(|| {
        let u = Float::from_integer_exp(&Integer::ONE, &(q.exponent() - &Integer::from(prec)));
        if q.is_negative() { u.neg() } else { u }
    });
    let next = move |c: &Float| unit.as_ref().map(|u| c.add(u));
    std::iter::successors(Some(q.clone()), next)
        .take(3)
        .filter(move |c| c.bit_length() <= prec)
}

/// Limbs for an approximate entry: one limb beyond the precision, and at
/// least the operands' length.
fn approx_len(prec: u64, operand: usize) -> usize {
    let need = usize::try_from(prec.div_ceil(LIMB_BITS)).expect("precision fits") + 1;
    need.max(operand)
}

fn is_odd(x: &Integer) -> bool {
    match x.to_i64() {
        Some(v) => v & 1 == 1,
        None => x.unsigned_abs().limbs()[0] & 1 == 1,
    }
}

/// Writes the top-aligned `m` at the top of `v`, zeros below.
fn pad(v: &mut [Limb], m: &[Limb]) {
    let low = v.len() - m.len();
    v[..low].fill(0);
    v[low..].copy_from_slice(m);
}
