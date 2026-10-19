//! Euler's constant γ by the Brent-McMillan algorithm B3 \[BM80\] with
//! binary splitting, as MPFR's `mpfr_const_euler` (CAS-104).
//!
//! With `u_k = (n^k / k!)^2`, `I_0(2n) = Σ u_k` and `S_0(2n) = Σ H_k u_k`
//! (`H_k = 1 + 1/2 + … + 1/k`), the Bessel functions give
//!
//! ```text
//! γ = S_0(2n) / I_0(2n) - K_0(2n) / I_0(2n) - log n
//!   = S_0 / I_0 - (I_0 K_0) / I_0^2 - log n,
//! ```
//!
//! and the product `I_0(2n) K_0(2n)` comes from its asymptotic expansion
//! `T = (1/4n) Σ_{k < 2n} ((2k)!)^3 / ((k!)^4 (16n)^(2k))`, within
//! `7 e^-4n / (4n)` for every integer `n >= 1` (Brent and Johansson, "A
//! bound for the error term in the Brent-McMillan algorithm", Math. Comp. 84
//! (2015), Lemma 4, \[BJ15\]). As `I_0(2n) > e^2n / √(4πn)`, that costs
//! `7π e^-8n` in γ, so `n ≈ w log(2) / 8` gives `w` bits. `S_0` and `I_0`
//! are summed to `N ≈ 4.97 n` terms, where their terms have fallen `e^-8n`
//! below `I_0` (ibid., Corollary 2 and Remark 2), and their tails are
//! bounded here. All three bounds go into the balls, so the result holds
//! whatever `n` and `N` are: they only set its accuracy.
//!
//! `S` and `I` share one splitting that carries the harmonic numbers
//! ([`split`], Haible and Papanikolaou's sums weighted by partial sums
//! \[HP98\]); `T` is a plain [`Series`]. `n` is free, so it is the first
//! integer of the form `2^a 3^b 5^c 7^d` above the bound: then `log n` is a
//! combination of four fast `atanh(1/m)` series ([`log_smooth`]), five times
//! cheaper than a general logarithm at `10^5` digits, where it would be a
//! third of the time.

use super::{GUARD, atanh_inv};
use crate::Integer;
use crate::ball::bsplit::{self, Series};
use crate::ball::{Ball, Mag};

#[cfg(test)]
mod tests;

/// `S` and `I` over the terms `k` in `[a, b)`, relative to `u_{a-1}`:
/// `P = n^(2(b-a))`, `Q = Π k^2`, `D = Π k`, `C / D = Σ 1/k`,
/// `T / Q = Σ u_k / u_{a-1}` and `V / (Q D) = Σ (H_k - H_{a-1}) u_k / u_{a-1}`.
/// `P` and `C` are `None` where not computed: the rightmost ranges never
/// need them.
struct Split {
    p: Option<Integer>,
    q: Integer,
    d: Integer,
    c: Option<Integer>,
    t: Integer,
    v: Integer,
}

/// The terms `[a, b)` (`a >= 1`) for `n2 = n^2`; `P` and `C` only when
/// `need_pc`. A leaf is `P = T = V = n^2`, `Q = k^2`, `D = k`, `C = 1`, and two
/// halves merge as
///
/// ```text
/// T = Tl Qr + Pl Tr,   V = Dr (Qr Vl + Cl Pl Tr) + Dl Pl Vr,
/// C = Cl Dr + Cr Dl,   P = Pl Pr,   Q = Ql Qr,   D = Dl Dr,
/// ```
///
/// since the right half's terms are `Pl / Ql` times its own, and its
/// harmonic numbers exceed its own by `Cl / Dl`.
fn split(n2: &Integer, a: u64, b: u64, need_pc: bool) -> Split {
    debug_assert!(1 <= a && a < b);
    if b - a == 1 {
        let k = Integer::from(a);
        return Split {
            p: Some(n2.clone()),
            q: &k * &k,
            d: k,
            c: Some(Integer::ONE),
            t: n2.clone(),
            v: n2.clone(),
        };
    }
    let m = a + (b - a) / 2;
    let l = split(n2, a, m, true);
    let r = split(n2, m, b, need_pc);
    let (Some(lp), Some(lc)) = (&l.p, &l.c) else {
        unreachable!("left halves keep P and C")
    };
    let pt = lp * &r.t;
    let mut t = &l.t * &r.q;
    t += &pt;
    let mut v = &r.q * &l.v;
    v += &(lc * &pt);
    let mut v = &v * &r.d;
    v += &(&(lp * &r.v) * &l.d);
    let (p, c) = match (need_pc, &r.p, &r.c) {
        (true, Some(rp), Some(rc)) => {
            let mut c = lc * &r.d;
            c += &(rc * &l.d);
            (Some(lp * rp), Some(c))
        }
        _ => (None, None),
    };
    Split {
        p,
        q: &l.q * &r.q,
        d: &l.d * &r.d,
        c,
        t,
        v,
    }
}

/// `T = (1/4n) Σ ((2k)!)^3 / ((k!)^4 (16n)^(2k))`: the ratio of consecutive
/// terms is `(2k-1)^3 / (32 k n^2)`, so `p(k) = (2k-1)^3`, `q(k) = k n^2`
/// with `2^5` as its shift, and `q(0) = n` with shift 2 (`p(0) = 1`).
struct Bessel(u64);

impl Series for Bessel {
    fn p(&self, k: u64) -> Option<Integer> {
        // (2k-1)^3 < 64 n^3 and k n^2 < 2 n^3 fit a u128 below n = 2^40.
        let x = 2 * u128::from(k.max(1)) - 1;
        Some(Integer::from(x * x * x))
    }

    fn q(&self, k: u64) -> Integer {
        let n = u128::from(self.0);
        Integer::from(if k == 0 { n } else { u128::from(k) * n * n })
    }

    fn a(&self, _: u64) -> Option<Integer> {
        None
    }

    fn b(&self, _: u64) -> Option<Integer> {
        None
    }

    fn shift(&self, k: u64) -> u64 {
        if k == 0 { 2 } else { 5 }
    }
}

/// Upper bounds of the tails `Σ_{k >= N} H_k u_k` of `S` and `Σ_{k >= N} u_k`
/// of `I`. From `k = N >= 4n` on, `u_{k+1} / u_k = n^2 / (k+1)^2 < 1/16` and
/// `H_{k+1} / H_k <= 1 + 1/(k+1) <= 2`, so the tails are below `8/7 H_N u_N`
/// and `16/15 u_N`. Here they are `2 H_N u_N` and `2 u_N`, with `H_N <= 1 +
/// log N < 1 + bits(N)` and `u_N = n^(2N) / (N^2 Q)` for the splitting's
/// `Q = ((N-1)!)^2`.
fn tails(n: u64, big_n: u64, q: &Integer) -> (Mag, Mag) {
    debug_assert!(big_n >= 4 * n);
    let n_sq = Mag::from_integer_lower(&Integer::from(u128::from(big_n).pow(2)));
    let fact = Mag::from_integer_lower(q).mul_lower(&n_sq);
    let u = Mag::from_integer(&Integer::from(n))
        .pow(2 * big_n)
        .div(&fact);
    let i_tail = u.mul_2exp(&Integer::ONE);
    let h = Mag::from_integer(&Integer::from(u64::from(big_n.ilog2()) + 2));
    (i_tail.mul(&h), i_tail)
}

/// `7 e^-4n / (4n)`, Brent and Johansson's bound on `|I_0(2n) K_0(2n) - T|`
/// (Lemma 4), with `e^-4n <= 2^-⌊4n λ⌋` for `λ = 1.442695040888963407`,
/// below `log2 e`.
fn bessel_error(n: u64) -> Mag {
    let four_n = 4 * u128::from(n);
    let k = four_n * 1_442_695_040_888_963_407 / 1_000_000_000_000_000_000;
    Mag::from_integer(&Integer::from(7u64))
        .mul(&Mag::power_of_two(&-Integer::from(k)))
        .div(&Mag::from_integer_lower(&Integer::from(four_n)))
}

/// The `m` of the Machin-like formulas for `log 2`, `log 3`, `log 5` and
/// `log 7`: each `atanh(1/m) = log((m + 1) / (m - 1)) / 2` is a combination
/// of the four logarithms (`126/125`, `225/224`, `2401/2400`, `4375/4374`),
/// with determinant -1.
const MACHIN: [u64; 4] = [251, 449, 4801, 8749];

/// `log p = Σ_i LOG_PRIMES[p][i] atanh(1 / MACHIN[i])` for `p = 2, 3, 5, 7`:
/// twice the inverse of that matrix.
const LOG_PRIMES: [[i64; 4]; 4] = [
    [144, 54, -38, 62],
    [228, 86, -60, 98],
    [334, 126, -88, 144],
    [404, 152, -106, 174],
];

/// The smallest `n >= lo` of the form `2^a 3^b 5^c 7^d`, and `[a, b, c, d]`.
/// Such numbers are dense enough that `n` exceeds `lo` by at most 5% from
/// `lo = 512` on, and 1.1% from `lo = 2^19`.
fn smooth_at_least(lo: u64) -> (u64, [u32; 4]) {
    let mut best = (u64::MAX, [0; 4]);
    let (mut p7, mut d) = (1u64, 0);
    loop {
        let (mut p5, mut c) = (p7, 0);
        loop {
            let (mut p3, mut b) = (p5, 0);
            loop {
                let a = lo.div_ceil(p3).next_power_of_two().ilog2();
                if p3 << a < best.0 {
                    best = (p3 << a, [a, b, c, d]);
                }
                if p3 >= lo {
                    break;
                }
                (p3, b) = (3 * p3, b + 1);
            }
            if p5 >= lo {
                break;
            }
            (p5, c) = (5 * p5, c + 1);
        }
        if p7 >= lo {
            break;
        }
        (p7, d) = (7 * p7, d + 1);
    }
    best
}

/// `log(2^a 3^b 5^c 7^d)` for the exponents `e`, at `w` bits: one
/// combination of the four `atanh(1/m)`. Its terms are at most a few times
/// the result, so `w + 8` bits in each are plenty.
fn log_smooth(e: [u32; 4], w: u64) -> Ball {
    let wp = w + 8;
    let mut sum = Ball::ZERO;
    for (i, &m) in MACHIN.iter().enumerate() {
        let k: i64 = (0..4).map(|p| i64::from(e[p]) * LOG_PRIMES[p][i]).sum();
        if k != 0 {
            let term = atanh_inv(m, wp).mul(&Ball::from_integer(&Integer::from(k)), wp);
            sum = sum.add(&term, wp);
        }
    }
    sum
}

/// Balls for the whole `S_0(2n)` and `I_0(2n)` at `wp` bits, from `N >= 4n`
/// terms and [`tails`]: `I = 1 + T / Q` and `S = V / (Q D)` over `k < N`.
fn sums(n: u64, big_n: u64, wp: u64) -> (Ball, Ball) {
    let s = split(&Integer::from(u128::from(n).pow(2)), 1, big_n, false);
    let (s_tail, i_tail) = tails(n, big_n, &s.q);
    let r = |x: &Integer| Ball::from_integer(x).round_to_accuracy(wp);
    let q = r(&s.q);
    let i = r(&(&s.q + &s.t)).div(&q, wp).add_error(&i_tail);
    let sum = r(&s.v).div(&q.mul(&r(&s.d), wp), wp).add_error(&s_tail);
    (sum, i)
}

/// The parameters `(n, N)` for `wp` bits, and the exponents of `n`:
/// `24 e^-8n <= 2^-wp` from `n >= (wp + log2 24) log(2) / 8`, with
/// `log(2) / 8 < 0.0866434`, rounded up to the next [`smooth_at_least`]; and
/// `N >= α n + 1`, `α (log α - 1) = 3` (Brent and Johansson, Remark 2).
fn parameters(wp: u64) -> (u64, [u32; 4], u64) {
    let lo = (u128::from(wp + 5) * 866_434).div_ceil(10_000_000);
    let (n, e) = smooth_at_least(u64::try_from(lo).expect("precision below 2^60"));
    let big_n = (u128::from(n) * 4_970_626).div_ceil(1_000_000) + 1;
    (n, e, u64::try_from(big_n).expect("precision below 2^60"))
}

/// γ at `prec + GUARD` bits; see the module comment.
pub(super) fn euler_uncached(prec: u64) -> Ball {
    let w = prec + GUARD;
    let (n, e, big_n) = parameters(w + EXTRA);
    euler_with(n, e, big_n, w)
}

/// Bits beyond `w` inside [`euler_with`]: `S/I - T/I^2` is `log n + γ`, and
/// subtracting `log n` cancels up to 6 bits for `n < 2^32`; the operations
/// before it cost a few more.
const EXTRA: u64 = 16;

/// γ at `w` bits from the parameter `n` (with exponents `e`, see
/// [`smooth_at_least`]) and `N >= 4n` terms of `S` and `I`. The ball holds γ
/// for any such `n` and `N`; [`parameters`] chooses them for `w` bits.
fn euler_with(n: u64, e: [u32; 4], big_n: u64, w: u64) -> Ball {
    let wp = w + EXTRA;
    let (sum, i) = sums(n, big_n, wp);
    let (t, d, shift) = bsplit::sum_shifted(&Bessel(n), 2 * n);
    let ik = bsplit::quotient(&t, &d, wp)
        .mul_2exp(&-Integer::from(shift))
        .add_error(&bessel_error(n));
    let log_n = log_smooth(e, wp);
    sum.div(&i, wp)
        .sub(&ik.div(&i.mul(&i, wp), wp), wp)
        .sub(&log_n, w)
}
