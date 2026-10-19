//! The splitting against the sums in exact rationals, the tail bounds
//! against the true tails, Brent and Johansson's bound on the Bessel product
//! against `I_0(2n) K_0(2n)` from MPFR at 8000 bits, the Machin-like
//! identities exactly, the smooth numbers against a search, and γ from too
//! small an `n` or too few terms against MPFR.

use rug::float::Constant;
use rug::ops::Pow;
use rug::{Float as MpFloat, Integer as MpInt, Rational as MpRat};

use super::{
    Bessel, LOG_PRIMES, MACHIN, bessel_error, euler_with, log_smooth, parameters, smooth_at_least,
    split, sums, tails,
};
use crate::Integer;
use crate::ball::Mag;
use crate::ball::bsplit;
use crate::testkit::ToRug;
use crate::testkit::float::{bracket, check_encloses};

fn mag_value(m: &Mag) -> MpRat {
    m.to_float().expect("finite").to_rational().to_rug()
}

/// `(u_k, H_k)` for `k < len`, `u_k = (n^k / k!)^2`, exactly.
fn terms(n: u64, len: u64) -> Vec<(MpRat, MpRat)> {
    let (mut u, mut h) = (MpRat::from(1), MpRat::new());
    (0..len)
        .map(|k| {
            if k > 0 {
                u *= MpRat::from((n * n, k * k));
                h += MpRat::from((1, k));
            }
            (u.clone(), h.clone())
        })
        .collect()
}

/// `Σ u_k` and `Σ H_k u_k` over `terms`.
fn sums_of(terms: &[(MpRat, MpRat)]) -> (MpRat, MpRat) {
    terms
        .iter()
        .fold((MpRat::new(), MpRat::new()), |(i, s), (u, h)| {
            (i + u, s + MpRat::from(h * u))
        })
}

#[test]
fn split_matches_the_sums() {
    for n in 1..=5u64 {
        let all = terms(n, 40);
        for big_n in 2..=40u64 {
            let (i_true, s_true) = sums_of(&all[..usize::try_from(big_n).unwrap()]);
            for need_pc in [false, true] {
                let s = split(&Integer::from(n * n), 1, big_n, need_pc);
                let (q, t, d, v) = (s.q.to_rug(), s.t.to_rug(), s.d.to_rug(), s.v.to_rug());
                // I = 1 + T / Q and S = V / (Q D).
                assert_eq!(MpRat::from((MpInt::from(&q + &t), q.clone())), i_true);
                assert_eq!(MpRat::from((v, MpInt::from(&q * &d))), s_true);
                if need_pc || big_n == 2 {
                    let p = s.p.expect("P kept").to_rug();
                    assert_eq!(p, MpInt::from(n).pow(2 * u32::try_from(big_n - 1).unwrap()));
                    let c = MpRat::from((s.c.expect("C kept").to_rug(), d));
                    assert_eq!(c, all[usize::try_from(big_n - 1).unwrap()].1);
                } else {
                    assert!(s.p.is_none() && s.c.is_none());
                }
            }
        }
    }
}

#[test]
fn tail_bounds() {
    for n in 1..=12u64 {
        // Past 300 more terms the rest is below 16^-300 of the first.
        let all = terms(n, 4 * n + 340);
        for big_n in 4 * n..4 * n + 40 {
            let s = split(&Integer::from(n * n), 1, big_n, false);
            let (s_tail, i_tail) = tails(n, big_n, &s.q);
            let (i_true, s_true) = sums_of(&all[usize::try_from(big_n).unwrap()..]);
            let (s_bound, i_bound) = (mag_value(&s_tail), mag_value(&i_tail));
            assert!(i_true <= i_bound, "I tail, n = {n}, N = {big_n}");
            assert!(s_true <= s_bound, "S tail, n = {n}, N = {big_n}");
            // Within the factors 2 and 2 (bits(N) + 1) / H_N of the proof.
            assert!(i_bound <= MpRat::from(&i_true * 4u32), "loose I tail");
            assert!(s_bound <= MpRat::from(&s_true * 8u32), "loose S tail");
        }
    }
}

/// Each of the balls for `S_0(2n)` and `I_0(2n)` holds the whole sum,
/// from `N = 4n` terms on, where the tails dominate their radii.
#[test]
fn sums_hold_the_series() {
    for n in 1..=12u64 {
        let all = terms(n, 4 * n + 340);
        let (i_true, s_true) = sums_of(&all);
        for big_n in [4 * n, 4 * n + 1, 4 * n + 7] {
            let (s, i) = sums(n, big_n, 300);
            for (name, ball, truth) in [("S", s, &s_true), ("I", i, &i_true)] {
                let (mid, rad) = (ball.mid().to_rational().to_rug(), mag_value(ball.rad()));
                let err = MpRat::from(&mid - truth).abs();
                assert!(
                    err <= rad,
                    "{name}, n = {n}, N = {big_n}: off by {err}, radius {rad}"
                );
            }
        }
    }
}

/// `I_0(2n) K_0(2n)` at `p` bits from `K_0(2n) = S_0(2n) - (log n + γ)
/// I_0(2n)`, the series summed until their terms fall below `2^-(p+64)`
/// of `I_0` (they shrink faster than `1/16` per term from `k = 4n` on).
fn bessel_product(n: u64, p: u32) -> MpFloat {
    let (mut u, mut h) = (MpFloat::with_val(p, 1), MpFloat::new(p));
    let (mut i0, mut s0) = (MpFloat::with_val(p, 1), MpFloat::new(p));
    let mut k = 0u64;
    while k < 4 * n || u.clone() << (p + 64) > i0 {
        k += 1;
        u *= n * n;
        u /= k * k;
        h += MpFloat::with_val(p, 1) / k;
        i0 += &u;
        s0 += MpFloat::with_val(p, &h * &u);
    }
    let log = MpFloat::with_val(p, n).ln() + MpFloat::with_val(p, Constant::Euler);
    let k0 = s0 - log * &i0;
    i0 * k0
}

#[test]
fn bessel_bound() {
    let p = 8000;
    for n in 1..=40u64 {
        let (t, d, e) = bsplit::sum_shifted(&Bessel(n), 2 * n);
        let t = MpRat::from((t.to_rug(), d.to_rug() << u32::try_from(e).unwrap()));
        let err = MpFloat::with_val(p, bessel_product(n, p) - &t).abs();
        // The product's own error is below 2^-7700 (K_0's cancellation
        // costs 6n bits), far inside the bound.
        let bound = mag_value(&bessel_error(n));
        assert!(
            err.to_rational().unwrap() <= bound,
            "n = {n}: |I0 K0 - T| = {err:.3e} above the bound"
        );
    }
}

#[test]
fn parameters_meet_the_bound() {
    for wp in [1u64, 2, 10, 64, 100, 1000, 12_345, 1 << 20, 1 << 40] {
        let (n, _, big_n) = parameters(wp);
        // 24 e^-8n <= 2^-wp, and N >= 4.9706 n + 1 >= 4n.
        #[allow(clippy::cast_precision_loss, reason = "n below 2^40")]
        let bits = 8.0 * n as f64 * std::f64::consts::LOG2_E - 24f64.log2();
        #[allow(clippy::cast_precision_loss, reason = "wp below 2^41")]
        let target = wp as f64;
        assert!(bits >= target, "wp = {wp}: n = {n} gives {bits} bits");
        assert!(
            big_n * 10_000 >= n * 49_706 + 10_000,
            "wp = {wp}: N = {big_n}"
        );
    }
}

#[test]
fn machin_identities() {
    // log p = Σ c_i atanh(1/m_i) = Σ c_i log((m_i + 1) / (m_i - 1)) / 2
    // exactly when Π ((m_i + 1) / (m_i - 1))^c_i = p^2.
    for (p, row) in [2u32, 3, 5, 7].into_iter().zip(LOG_PRIMES) {
        let mut prod = MpRat::from(1);
        for (&m, &c) in MACHIN.iter().zip(&row) {
            let r = MpRat::from((m + 1, m - 1));
            let r = if c < 0 { r.recip() } else { r };
            prod *= r.pow(u32::try_from(c.unsigned_abs()).unwrap());
        }
        assert_eq!(prod, MpRat::from(p * p), "log {p}");
    }
}

#[test]
fn smooth_numbers() {
    let smooth = |mut x: u64| {
        for p in [2, 3, 5, 7] {
            while x % p == 0 {
                x /= p;
            }
        }
        x == 1
    };
    for lo in 1..=20_000u64 {
        let (n, [a, b, c, d]) = smooth_at_least(lo);
        let want = (lo..=2 * lo).find(|&x| smooth(x)).unwrap();
        assert_eq!(n, want, "lo = {lo}");
        assert_eq!(n, 2u64.pow(a) * 3u64.pow(b) * 5u64.pow(c) * 7u64.pow(d));
    }
    // The density claimed: within 5% from 512 on, 1.1% from 2^19.
    for lo in (512..1u64 << 22).step_by(97) {
        let (n, _) = smooth_at_least(lo);
        let slack = if lo >= 1 << 19 { 1011 } else { 1050 };
        assert!(n * 1000 <= lo * slack, "lo = {lo}: n = {n}");
    }
}

#[test]
fn smooth_logarithms() {
    for (e, w) in [
        ([0, 0, 0, 0], 64),
        ([1, 0, 0, 0], 100),
        ([0, 1, 1, 1], 1000),
        ([10, 3, 0, 2], 3000),
        ([41, 0, 5, 9], 20_000),
    ] {
        let b = log_smooth(e, w);
        let n = (0..4).fold(MpInt::from(1), |acc, i| {
            acc * MpInt::from([2u32, 3, 5, 7][i]).pow(e[i])
        });
        let truth = MpFloat::with_val(u32::try_from(w).unwrap() + 64, &n).ln();
        let diff = MpFloat::with_val(64, &truth - &b.mid().to_rational().to_rug()).abs();
        let rad = mag_value(b.rad());
        assert!(
            diff.to_rational().unwrap() <= rad,
            "log of {n}: off by {diff:.3e}"
        );
        if n != 1 {
            let bits = b.rel_accuracy_bits();
            assert!(bits >= i64::try_from(w).unwrap(), "log of {n}: {bits} bits");
        }
    }
}

/// Too small an `n` or too few terms for the precision: the Bessel bound
/// (`N >= α n + 1`) or the tails (`N = 4n`) dominate the radius in turn,
/// and the ball must still hold γ.
#[test]
fn few_terms() {
    let w = 400;
    let br = bracket(u32::try_from(w).unwrap(), |p, r| {
        MpFloat::with_val_round(p, Constant::Euler, r).0
    });
    for (n, e) in [
        (1u64, [0, 0, 0, 0]),
        (2, [1, 0, 0, 0]),
        (6, [1, 1, 0, 0]),
        (10, [1, 0, 1, 0]),
        (21, [0, 1, 0, 1]),
        (35, [0, 0, 1, 1]),
    ] {
        let enough = (n * 4_970_626).div_ceil(1_000_000) + 1;
        for big_n in [4 * n, 4 * n + 3, enough] {
            let b = euler_with(n, e, big_n, w);
            let (mid, rad) = (b.mid().to_rational().to_rug(), mag_value(b.rad()));
            if let Err(err) = check_encloses(&mid, &rad, &br) {
                panic!("n = {n}, N = {big_n}: {err}");
            }
            // A meaningful containment: n = 1 alone has a radius near 1/50.
            assert!(
                rad < MpRat::from((1, 16)),
                "n = {n}, N = {big_n}: radius {rad}"
            );
        }
    }
}
