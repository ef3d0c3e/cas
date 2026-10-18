//! Fixed-point series summation (CAS-DOC-5 §7.2): values are integers `A`
//! standing for `A / 2^w`, each with an error bound `e` in units of `2^-w`
//! (`|value - A / 2^w| <= e / 2^w`), so a sum costs integer products and no
//! per-term normalisation. Error bounds are computed at run time, alongside
//! the values, by the rules written next to each step.

use crate::{Integer, Natural};

/// `x / 2^w`, truncated toward zero (error below one unit).
pub(super) fn shr(x: Integer, w: u64) -> Integer {
    let (neg, mag) = x.into_parts();
    Integer::from_parts(neg, mag >> w)
}

/// `ceil(a / b)` for `b >= 1`.
fn div_ceil(a: &Natural, b: &Natural) -> Natural {
    let (q, r) = a.div_rem(b);
    if r == Natural::ZERO {
        q
    } else {
        &q + &Natural::from(1u64)
    }
}

/// `Σ_{j < n} x^j / j!` for `x = X / 2^w` with `|x| <= 1/2`, `n >= 1`, by
/// rectangular splitting (Paterson–Stockmeyer): `(S, e)` with
/// `|Σ - S / 2^w| <= e / 2^w`. The tail beyond `n` terms is the caller's.
///
/// With blocks of `m ≈ √n` terms, `acc_i = Σ_{j >= m i} x^(j - m i) (m i)! /
/// j!` (over `j < n`) satisfies, for a block of `L` terms,
///
/// ```text
/// acc_i = [Σ_{l < L} x^l c_l + x^m acc_(i+1) / (m i + m)] / D,
/// c_l = Π_{t = l+1}^{L-1} (m i + t),  D = c_0,
/// ```
///
/// (the second term only for full blocks with a successor), since `(m i)! /
/// (m i + l)! = c_l / D` and `(m i)! / (m (i + 1))! = 1 / (D (m i + m))`.
/// So the sum takes `m` full products for the powers `x^l` and one per
/// block, the rest being products and quotients by small integers; `S =
/// acc_0`.
///
/// Errors, in units of `2^-w`, with every true value at most 2 in magnitude
/// (`|acc_i| <= Σ |x|^k <= 2`, `|x^l| <= 1`):
/// - powers `P_l = trunc(P_(l-1) X / 2^w)`: `e_l <= |x| e_(l-1) + 1`, so
///   `e_l = l - 1` from the exact `P_1 = X`;
/// - `trunc(P_m A / 2^w)`: `|x^m| e_A + |acc| e_m + e_m e_A 2^-w + 1 <= e_A +
///   2 e_m + 2` (while `e_m e_A < 2^w`), then a truncated quotient by `q`
///   adds 1 to the bound divided by `q`;
/// - the numerator `Σ P_l c_l + T` is exact integer arithmetic on inexact
///   values: `Σ e_l c_l + e_T`; its truncated quotient by `D` adds 1 to that
///   bound divided by `D`.
pub(super) fn exp_sum(x: &Integer, w: u64, n: u64) -> (Integer, Natural) {
    assert!(n >= 1, "exp_sum: no terms");
    let m = n.isqrt().max(1);
    let one = Integer::from_parts(false, Natural::from(1u64) << w);
    // Powers 0..=m and their errors.
    let mut pow = Vec::with_capacity(usize::try_from(m + 1).expect("small"));
    pow.push(one);
    pow.push(x.clone());
    for _ in 2..=m {
        let last = pow.last().expect("non-empty");
        let next = shr(last * x, w);
        pow.push(next);
    }
    let pow_err = |l: u64| Natural::from(l.saturating_sub(1));
    let blocks = n.div_ceil(m);
    let mut acc = Integer::ZERO;
    let mut acc_err = Natural::ZERO;
    for i in (0..blocks).rev() {
        let start = m * i;
        let len = m.min(n - start);
        // c_l from l = len - 1 down to 0, and the numerator with its error.
        let mut c = Natural::from(1u64);
        let mut num = Integer::ZERO;
        let mut num_err = Natural::ZERO;
        for l in (0..len).rev() {
            if l + 1 < len {
                c = &c * &Natural::from(start + l + 1);
            }
            let p = &pow[usize::try_from(l).expect("small")];
            num += &(p * &Integer::from(&c));
            num_err = &num_err + &(&pow_err(l) * &c);
        }
        let d = c;
        if i + 1 < blocks {
            // T = trunc(trunc(P_m acc / 2^w) / (m i + m)).
            let q = Natural::from(start + m);
            let prod = shr(&pow[usize::try_from(m).expect("small")] * &acc, w);
            let t = prod.div_rem_trunc(&Integer::from(&q)).0;
            let prod_err =
                &(&acc_err + &(&pow_err(m) * &Natural::from(2u64))) + &Natural::from(2u64);
            let t_err = &div_ceil(&prod_err, &q) + &Natural::from(1u64);
            num += &t;
            num_err = &num_err + &t_err;
        }
        acc = num.div_rem_trunc(&Integer::from(&d)).0;
        acc_err = &div_ceil(&num_err, &d) + &Natural::from(1u64);
    }
    (acc, acc_err)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::exp_sum;
    use crate::{Integer, Natural, Rational};

    /// `Σ_{j < n} x^j / j!` exactly.
    fn exact(x: &Rational, n: u64) -> Rational {
        let (mut s, mut t) = (Rational::from(0i64), Rational::from(1i64));
        for j in 0..n {
            if j > 0 {
                t = &(&t * x) / &Rational::from(Integer::from(j));
            }
            s = &s + &t;
        }
        s
    }

    proptest! {
        #[test]
        fn within_the_bound(w in 8u64..400, n in 1u64..80, num in any::<i64>(), neg in any::<bool>()) {
            // |x| <= 1/2: X in [-2^(w-1), 2^(w-1)].
            let half = Integer::from_parts(false, Natural::from(1u64) << (w - 1));
            let raw = Integer::from(num);
            let span = &half + &Integer::from(1i64);
            let x = raw.div_rem_euclid(&span).1;
            let x = if neg { -x } else { x };
            let (s, e) = exp_sum(&x, w, n);
            let scale = Rational::from(Integer::from_parts(false, Natural::from(1u64) << w));
            let xr = &Rational::from(x) / &scale;
            let got = &Rational::from(s) / &scale;
            let diff = &exact(&xr, n) - &got;
            let diff = if diff < Rational::from(0i64) { -diff } else { diff };
            let bound = &Rational::from(Integer::from(e.clone())) / &scale;
            prop_assert!(diff <= bound, "error {diff} above {bound}");
            // And the bound is not wildly loose: a few units per block.
            prop_assert!(Integer::from(e.clone()) <= Integer::from(4 * n + 8), "loose bound {e:?}");
        }
    }
}
