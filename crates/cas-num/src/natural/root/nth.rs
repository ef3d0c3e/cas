//! The kth root for `k >= 3`, on limbs: one Newton step per doubling of the
//! precision, as GMP's `mpn_rootrem` takes (Brent and Zimmermann, *Modern
//! Computer Arithmetic*, §1.5.2).
//!
//! From the root `s'` of the input's top bits and its remainder, a division
//! by `k s'^(k-1)` gives the next bits of the root, at most two too large.
//! The powers `s^(k-1)` and `s^k` then check them and give the next step's
//! divisor and remainder. Without a remainder, the last step skips the
//! powers: it computes guard bits below the root, which absorb the
//! estimate's error.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithm uses the textbook names (k, s, r, p, w, u, d, q, b, m, n)"
)]

#[cfg(test)]
use std::cell::Cell;

use super::super::bits::split_bit_index;
use super::Natural;
use crate::kernel::div::{divappr_q, divappr_q_scratch_len};
use crate::kernel::mul::{
    MUL_KARATSUBA_THRESHOLD, MUL_WRAP_NTT_THRESHOLD, SQR_KARATSUBA_THRESHOLD,
};
use crate::{Limb, kernel};

// Counters for the rarely taken branches, so tests can check that their
// inputs reach them.
#[cfg(test)]
thread_local! {
    /// A Newton step's estimate was too large: its powers were computed
    /// again.
    pub(super) static STEP_DECREMENT: Cell<u64> = const { Cell::new(0) };
    /// A Newton step's estimate reached `2^b`, past its bits.
    pub(super) static STEP_CLAMP: Cell<u64> = const { Cell::new(0) };
    /// The first step's floating-point estimate was too small.
    pub(super) static FIRST_INCREMENT: Cell<u64> = const { Cell::new(0) };
    /// A step took fewer bits to fit its division to a transform length.
    pub(super) static STEP_FIT: Cell<u64> = const { Cell::new(0) };
    /// The guard bits did not decide the root, so a power did.
    pub(super) static GUARD_CHECK: Cell<u64> = const { Cell::new(0) };
}

macro_rules! hit {
    ($counter:ident) => {
        #[cfg(test)]
        $counter.with(|c| c.set(c.get() + 1));
    };
}

/// Bits of the root that the first step estimates in floating point, to
/// within about `2^-48.5` relative, so off by under one: it settles with a
/// power or two.
const BASE_ROOT_BITS: u64 = 46;

/// Limbs of [`newton_root`]'s buffers that it takes from the stack rather
/// than the heap: enough for inputs of up to 26 limbs.
const STACK_LIMBS: usize = 320;

/// Guard bits that [`root`] computes below the root, from
/// [`GUARD_MIN_BITS`] on, rather than the last step's powers.
const GUARD_BITS: u64 = 32;

/// Root lengths, in bits, from which [`root`] computes guard bits: below,
/// the last step's powers cost less than a longer root.
const GUARD_MIN_BITS: u64 = 256;

/// The bit length of `⌊a^(1/k)⌋` for a normalised non-zero `a`.
fn root_bits(a: &[Limb], k: u32) -> u64 {
    kernel::bit_length(a).div_ceil(u64::from(k))
}

/// `⌊x^(1/k)⌋` for `k >= 3`.
///
/// From [`GUARD_MIN_BITS`] on, this is the root `S` of `x 2^(k G)`, for
/// `G` guard bits, shifted right by `G`: the last step gives `S` up to 2
/// too large, which changes `⌊S / 2^G⌋` only if `S mod 2^G` is 0 or 1.
/// Only then is the root checked by a power.
pub(super) fn root(x: &Natural, k: u32) -> Natural {
    let a = x.limbs();
    if a.is_empty() {
        return Natural::ZERO;
    }
    if root_bits(a, k) < GUARD_MIN_BITS {
        return newton_root(a, k, true).0;
    }
    // A root of 256 bits needs 255 k bits, so the guard bits add under an
    // eighth to x.
    let (big, _) = newton_root((x << (u64::from(k) * GUARD_BITS)).limbs(), k, false);
    let s = &big >> GUARD_BITS;
    if big.limbs()[0] & ((1 << GUARD_BITS) - 1) >= 2 {
        return s;
    }
    hit!(GUARD_CHECK);
    if s.pow(u64::from(k)) <= *x {
        s
    } else {
        s - 1u64
    }
}

/// `(s, x - s^k)` with `s = ⌊x^(1/k)⌋`, for `k >= 3`.
pub(super) fn root_rem(x: &Natural, k: u32) -> (Natural, Natural) {
    let a = x.limbs();
    if a.is_empty() {
        return (Natural::ZERO, Natural::ZERO);
    }
    let (s, r) = newton_root(a, k, true);
    (s, r.expect("an exact last step gives the remainder"))
}

/// The root `s` of the normalised non-zero `a` and its remainder `a - s^k`,
/// for `k >= 3`, by one Newton step per doubling of the precision from the
/// root of `a`'s top bits. When `exact` is false, the last step skips the
/// powers: `s` is then up to 2 too large, and there is no remainder.
///
/// A step goes from `s'`, the root of `u'` with `m` bits, to the root
/// `s = s' 2^b + q` of `u = u' 2^(kb) + u₀` (`u₀ < 2^(kb)`). In
/// `u - s'^k 2^(kb) = (s' 2^b + q)^k - s'^k 2^(kb) + r`, the first term of
/// the binomial expansion is `k s'^(k-1) q 2^((k-1) b)`, so `q` is
/// estimated by `⌊D / (k s'^(k-1))⌋` with `D = ⌊(u - s'^k 2^(kb)) /
/// 2^((k-1) b)⌋`. The other terms and the remainder `r` only raise the
/// estimate, so it is at least `q`. With `b <= m - c` and `2^c >= 64k`,
/// they raise it by less than 2: the remainder by less than
/// `(1 + 1/s')^(k-1) < e^(1/64)`, and the rest, mostly the second term
/// `(k - 1) q² / (2 s' 2^b)`, by less than `1/64 + 1/s'`. So the estimate is
/// too large only when the remainder is that close to its largest, about
/// one step in 200. [`RootWork::step`]'s approximate division adds at most
/// one more, rarely. The powers `s^(k-1)` and `s^k` then check the
/// estimate, and are computed again when it needs decrementing.
fn newton_root(a: &[Limb], k: u32, exact: bool) -> (Natural, Option<Natural>) {
    let k64 = u64::from(k);
    let n = root_bits(a, k);
    // c = ⌈log₂ k⌉ + 6, at most 38, so a step from more than
    // BASE_ROOT_BITS appends at least one bit.
    let c = u64::from(u32::BITS - (k - 1).leading_zeros()) + 6;
    // Root lengths from n down to the first step's: at most 64, as each
    // step about halves the length.
    let mut sizes = [0; 64];
    let mut steps = 0;
    let mut m = n;
    while m > BASE_ROOT_BITS {
        sizes[steps] = m;
        steps += 1;
        m -= step_bits(m, c);
    }
    // Every power is of an estimate whose kth power is below 2a, plus a limb
    // of slack for each operand of a product. The first step's estimates
    // are near u by logarithms; the later ones exceed roots of at least
    // 2^46 by at most 2, and (1 + 2^-45)^k < 2 for k < 2^32. Short roots
    // work on the stack.
    let len = a.len() + 4;
    let need = 9 * len + len + 16;
    let mut stack = [0; STACK_LIMBS];
    let mut heap = Vec::new();
    let buf = if need <= STACK_LIMBS {
        &mut stack[..need]
    } else {
        heap.resize(need, 0);
        &mut heap[..]
    };
    let (bufs, fixed) = buf.split_at_mut(9 * len);
    let mut parts = bufs.chunks_exact_mut(len);
    let mut part = || parts.next().expect("nine buffers");
    let mut w = RootWork {
        k,
        u: part(),
        s: part(),
        r: part(),
        p: part(),
        w: part(),
        d: part(),
        q: part(),
        t: part(),
        tmp: part(),
        fixed,
        extra: Vec::new(),
        ul: 0,
        sl: 0,
        rl: 0,
        pl: 0,
        wl: 0,
    };
    w.load(a, k64 * (n - m));
    w.first();
    for &next in sizes[..steps].iter().rev() {
        w.load(a, k64 * (n - next));
        w.step(next - m);
        if next == n && !exact {
            return (Natural::from_limbs(&w.s[..w.sl]), None);
        }
        if w.settle() {
            hit!(STEP_DECREMENT);
        }
        m = next;
    }
    (
        Natural::from_limbs(&w.s[..w.sl]),
        Some(Natural::from_limbs(&w.r[..w.rl])),
    )
}

/// The bits `b` that the step to a root of `n` bits appends: at most `n - m`
/// for the root of `m` bits it starts from, less `c` (see [`newton_root`]),
/// so about half. The step's division ([`divappr_q`]) works with a divisor
/// and a quotient of up to `⌈(b + 1) / 64⌉ + 2` limbs (its guard limb
/// included). Where that would just pass a transform length of the NTT,
/// which costs a third or a half more than the length below, the step
/// takes a few limbs fewer and leaves them to the one before.
fn step_bits(n: u64, c: u64) -> u64 {
    let b = n - (n + c).div_ceil(2);
    let bits = u64::from(Limb::BITS);
    let dn = (b + 1).div_ceil(bits) + 2;
    let pow2 = 1 << (u64::BITS - 1 - dn.leading_zeros());
    // The longest transform length (2^j or 3 2^j) up to dn.
    let fit = if 3 * pow2 / 2 <= dn {
        3 * pow2 / 2
    } else {
        pow2
    };
    let threshold = u64::try_from(MUL_WRAP_NTT_THRESHOLD).unwrap_or(u64::MAX);
    if dn > fit && fit >= threshold && dn - fit <= fit / 32 {
        hit!(STEP_FIT);
        bits * (fit - 2) - 1
    } else {
        b
    }
}

/// The buffers of [`newton_root`], with the normalised lengths of the
/// values they hold between steps.
struct RootWork<'a> {
    k: u32,
    /// The input's top bits for the current step, `ul` limbs.
    u: &'a mut [Limb],
    /// The root of `u`, `sl` limbs, and its remainder `u - s^k`, `rl`.
    s: &'a mut [Limb],
    r: &'a mut [Limb],
    /// `p = s^(k-1)` and `w = k p`, `pl` and `wl` limbs.
    p: &'a mut [Limb],
    w: &'a mut [Limb],
    /// A step's dividend and quotient, `s^k`, and a power's other buffer.
    d: &'a mut [Limb],
    q: &'a mut [Limb],
    t: &'a mut [Limb],
    tmp: &'a mut [Limb],
    /// Scratch for products and divisions, and more when that is short.
    fixed: &'a mut [Limb],
    extra: Vec<Limb>,
    ul: usize,
    sl: usize,
    rl: usize,
    pl: usize,
    wl: usize,
}

impl RootWork<'_> {
    /// `u = a >> shift` (bits).
    fn load(&mut self, a: &[Limb], shift: u64) {
        let (limbs, bits) = split_bit_index(shift);
        let len = a.len() - limbs;
        kernel::shr(&mut self.u[..len], &a[limbs..], bits);
        self.ul = kernel::normalized_len(&self.u[..len]);
    }

    /// The root of `u`, which has at most [`BASE_ROOT_BITS`] bits: exact
    /// when `u` fits a limb, else from its top limb `t` in floating point.
    /// With `u` about `t 2^e` and `e = qk + r`, the root is about
    /// `t^(1/k) 2^(r/k) 2^q`, where `powf` and `exp2` are within an ulp or
    /// two and `1/k` rounds off at most `2^-53` relative (`2^-49` on the
    /// root, as `ln t < 45`). Logarithms settle candidates whose power is
    /// far from `u`, so that the exact powers stay about `u`'s size: with a
    /// root far below `k`, `(s + 1)^k` is far above `u`.
    fn first(&mut self) {
        let u = &self.u[..self.ul];
        let bits = kernel::bit_length(u);
        if bits <= u64::from(Limb::BITS) {
            self.s[0] = root_u64(u[0], self.k);
            self.sl = 1;
            self.settle();
            return;
        }
        let e = bits - u64::from(Limb::BITS);
        let k = f64::from(self.k);
        #[allow(
            clippy::cast_precision_loss,
            reason = "logarithms, off by under 2^-40 relative"
        )]
        let (log_u, log) = ((bits_at(u, e) as f64).log2() + e as f64, |s: Limb| {
            k * (s as f64).log2()
        });
        let eps = (k + log_u) * f64::powi(2.0, -40);
        let (q, r) = (e / u64::from(self.k), e % u64::from(self.k));
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "an estimate, settled exactly; q < 46"
        )]
        let mut s = ((bits_at(u, e) as f64).powf(k.recip())
            * (r as f64 / k).exp2()
            * (1u64 << q) as f64) as Limb;
        s = s.max(1);
        loop {
            if log(s) > log_u + eps {
                s -= 1;
            } else if log(s + 1) < log_u - eps {
                s += 1;
            } else {
                break;
            }
        }
        loop {
            self.s[0] = s;
            self.sl = 1;
            // Done if s was too large, or if (s + 1)^k > u: when
            // u - s^k < k s^(k-1) or by logarithms.
            if self.settle()
                || kernel::cmp(&self.r[..self.rl], &self.w[..self.wl]).is_lt()
                || log(s + 1) > log_u + eps
            {
                return;
            }
            hit!(FIRST_INCREMENT);
            s += 1;
        }
    }

    /// Appends `b` bits to the root: `s = s 2^b + q` for the estimate `q`
    /// of [`newton_root`]'s doc, through `d`, from an approximate quotient
    /// ([`divappr_q`]): the true one or one more. Only the top limbs of `w`
    /// matter, and no remainder is computed.
    fn step(&mut self, b: u64) {
        let k = u64::from(self.k);
        // D: the remainder followed by bits (k - 1) b .. k b of u. They
        // lie below the remainder's bits, so D keeps the shifted
        // remainder's length; if the remainder is zero, then D < 2^b < w
        // and the quotient is zero anyway.
        let dl = shl_bits(self.d, &self.r[..self.rl], b);
        or_field(self.d, &self.u[..self.ul], (k - 1) * b, b);
        let wl = self.wl;
        let mut ql = 0;
        if dl >= wl {
            ql = dl - wl + 1;
            let scratch = take(self.fixed, &mut self.extra, divappr_q_scratch_len(dl, wl));
            divappr_q(&mut self.q[..ql], &self.d[..dl], &self.w[..wl], scratch);
            ql = kernel::normalized_len(&self.q[..ql]);
        }
        // q < 2^b, as the root has b more bits.
        if kernel::bit_length(&self.q[..ql]) > b {
            hit!(STEP_CLAMP);
            let (limbs, bits) = split_bit_index(b);
            self.q[..limbs].fill(Limb::MAX);
            ql = limbs;
            if bits > 0 {
                self.q[limbs] = (1 << bits) - 1;
                ql += 1;
            }
        }
        // s 2^b + q, through t, which is free until the next powers.
        let sl = shl_bits(self.t, &self.s[..self.sl], b);
        for (s, &q) in self.t.iter_mut().zip(&self.q[..ql]) {
            *s |= q;
        }
        self.sl = sl.max(ql);
        std::mem::swap(&mut self.s, &mut self.t);
    }

    /// Decrements the root `s` while `s^k > u`, and returns whether it did.
    /// Leaves `p = s^(k-1)`, `w = k p` and `r = u - s^k`.
    fn settle(&mut self) -> bool {
        let mut decremented = false;
        loop {
            let e = u64::from(self.k) - 1;
            let s = &self.s[..self.sl];
            self.pl = pow(self.p, self.tmp, s, e, self.fixed, &mut self.extra);
            let tl = self.pl + self.sl;
            mul(
                &mut self.t[..tl],
                &self.p[..self.pl],
                s,
                self.fixed,
                &mut self.extra,
            );
            let tl = kernel::normalized_len(&self.t[..tl]);
            if kernel::cmp(&self.t[..tl], &self.u[..self.ul]).is_le() {
                let ul = self.ul;
                kernel::sub(&mut self.r[..ul], &self.u[..ul], &self.t[..tl]);
                self.rl = kernel::normalized_len(&self.r[..ul]);
                self.w[self.pl] =
                    kernel::mul_1(&mut self.w[..self.pl], &self.p[..self.pl], self.k.into());
                self.wl = kernel::normalized_len(&self.w[..=self.pl]);
                return decremented;
            }
            kernel::sub_1_assign(&mut self.s[..self.sl], 1);
            self.sl = kernel::normalized_len(&self.s[..self.sl]);
            decremented = true;
        }
    }
}

/// `need` limbs of scratch: `fixed` when it is long enough, else `extra`,
/// grown to fit.
fn take<'s>(fixed: &'s mut [Limb], extra: &'s mut Vec<Limb>, need: usize) -> &'s mut [Limb] {
    if need <= fixed.len() {
        &mut fixed[..need]
    } else {
        if extra.len() < need {
            extra.resize(need, 0);
        }
        &mut extra[..need]
    }
}

/// `r = a b` for `a.len() >= b.len()`: below Karatsuba's threshold, the
/// schoolbook product directly, without the dispatch.
fn mul(r: &mut [Limb], a: &[Limb], b: &[Limb], fixed: &mut [Limb], extra: &mut Vec<Limb>) {
    if b.len() < MUL_KARATSUBA_THRESHOLD {
        kernel::mul_basecase(r, a, b);
    } else {
        let scratch = take(fixed, extra, kernel::mul_scratch_len(a.len(), b.len()));
        kernel::mul(r, a, b, scratch);
    }
}

/// `dst = s^e` for a normalised non-zero `s` and `e >= 1`, by left-to-right
/// binary powering through `tmp`, both long enough for every power on the
/// way; returns the length.
fn pow(
    dst: &mut [Limb],
    tmp: &mut [Limb],
    s: &[Limb],
    e: u64,
    fixed: &mut [Limb],
    extra: &mut Vec<Limb>,
) -> usize {
    let top = u64::BITS - 1 - e.leading_zeros();
    // Start in whichever buffer leaves the last product in dst.
    let products = top + e.count_ones() - 1;
    let (mut cur, mut other) = if products % 2 == 0 {
        (dst, tmp)
    } else {
        (tmp, dst)
    };
    let sn = s.len();
    cur[..sn].copy_from_slice(s);
    let mut len = sn;
    for i in (0..top).rev() {
        if len < SQR_KARATSUBA_THRESHOLD {
            kernel::sqr_basecase(&mut other[..2 * len], &cur[..len]);
        } else {
            let scratch = take(fixed, extra, kernel::sqr_scratch_len(len));
            kernel::sqr(&mut other[..2 * len], &cur[..len], scratch);
        }
        len = kernel::normalized_len(&other[..2 * len]);
        std::mem::swap(&mut cur, &mut other);
        if (e >> i) & 1 == 1 {
            mul(&mut other[..len + sn], &cur[..len], s, fixed, extra);
            len = kernel::normalized_len(&other[..len + sn]);
            std::mem::swap(&mut cur, &mut other);
        }
    }
    len
}

/// `dst = src << shift` (bits), `dst` long enough for the result; returns
/// its length.
fn shl_bits(dst: &mut [Limb], src: &[Limb], shift: u64) -> usize {
    let (limbs, bits) = split_bit_index(shift);
    let len = src.len();
    dst[..limbs].fill(0);
    dst[limbs + len] = kernel::shl(&mut dst[limbs..limbs + len], src, bits);
    kernel::normalized_len(&dst[..=limbs + len])
}

/// ORs bits `lo .. lo + len` of `src` into the low `len` bits of `dst`,
/// which are zero.
fn or_field(dst: &mut [Limb], src: &[Limb], lo: u64, len: u64) {
    let (limbs, bits) = split_bit_index(len);
    let n = limbs + usize::from(bits > 0);
    let mut pos = lo;
    for (i, d) in dst[..n].iter_mut().enumerate() {
        let mut word = bits_at(src, pos);
        if i == limbs {
            word &= (1 << bits) - 1;
        }
        *d |= word;
        pos += u64::from(Limb::BITS);
    }
}

/// The 64 bits of `a` from bit `pos` on, zero past its end.
fn bits_at(a: &[Limb], pos: u64) -> Limb {
    let (limb, bit) = split_bit_index(pos);
    let lo = a.get(limb).map_or(0, |&x| x >> bit);
    let hi = match (bit, a.get(limb + 1)) {
        (0, _) | (_, None) => 0,
        (_, Some(&x)) => x << (Limb::BITS - bit),
    };
    lo | hi
}

/// `⌊v^(1/n)⌋` for `n >= 2`: a floating-point estimate, corrected.
pub(super) fn root_u64(v: u64, n: u32) -> u64 {
    if v < 2 || n >= u64::BITS {
        return v.min(1);
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "an estimate, corrected below"
    )]
    let mut r = (v as f64).powf(1.0 / f64::from(n)) as u64;
    let fits = |r: u64| r.checked_pow(n).is_some_and(|p| p <= v);
    while !fits(r) {
        r -= 1;
    }
    while fits(r + 1) {
        r += 1;
    }
    r
}
