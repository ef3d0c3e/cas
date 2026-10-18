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

/// `Σ_{j < n} x^j / j!` for `x = X / 2^w` with `|x| <= 1/2`, `n >= 1`: see
/// [`hyp_sum`].
pub(super) fn exp_sum(x: &Integer, w: u64, n: u64) -> (Integer, Natural) {
    hyp_sum(x, w, n, |t| i64::try_from(t).expect("small"))
}

/// `Σ_{j < n} x^j / (q(1) q(2) … q(j))` for `x = X / 2^w` with `|x| <= 1/2`,
/// integers `|q(t)| >= 1`, `n >= 1`, by rectangular splitting
/// (Paterson–Stockmeyer): `(S, e)` with `|Σ - S / 2^w| <= e / 2^w`. The
/// tail beyond `n` terms is the caller's. `q(t) = t` is `exp`; `-(2t-1)(2t)`
/// and `-(2t)(2t+1)` are `cos` and `sin(x)/x` in `x^2`.
///
/// Write `Q(j) = q(1) … q(j)`. With blocks of `m ≈ √n` terms, `acc_i =
/// Σ_{j >= m i} x^(j - m i) Q(m i) / Q(j)` (over `j < n`) satisfies, for a
/// block of `L` terms,
///
/// ```text
/// acc_i = [Σ_{l < L} x^l c_l + x^m acc_(i+1) / q(m i + m)] / D,
/// c_l = Π_{t = l+1}^{L-1} q(m i + t),  D = c_0,
/// ```
///
/// (the second term only for full blocks with a successor), since `Q(m i) /
/// Q(m i + l) = c_l / D` and `Q(m i) / Q(m i + m) = 1 / (D q(m i + m))`.
/// So the sum takes `m` full products for the powers `x^l` and one per
/// block, the rest being products and quotients by small integers; `S =
/// acc_0`.
///
/// Errors, in units of `2^-w`, with every true value at most 2 in magnitude
/// (`|acc_i| <= Σ |x|^k <= 2` as `|Q(m i) / Q(j)| <= 1`, `|x^l| <= 1`):
/// - powers `P_l = trunc(P_(l-1) X / 2^w)`: `e_l <= |x| e_(l-1) + 1`, so
///   `e_l = l - 1` from the exact `P_1 = X`;
/// - `trunc(P_m A / 2^w)`: `|x^m| e_A + |acc| e_m + e_m e_A 2^-w + 1 <= e_A +
///   2 e_m + 2` (while `e_m e_A < 2^w`), then a truncated quotient by `q`
///   adds 1 to the bound divided by `|q|`;
/// - the numerator `Σ P_l c_l + T` is exact integer arithmetic on inexact
///   values: `Σ e_l |c_l| + e_T`; its truncated quotient by `D` adds 1 to
///   that bound divided by `|D|`.
pub(super) fn hyp_sum(x: &Integer, w: u64, n: u64, q: impl Fn(u64) -> i64) -> (Integer, Natural) {
    assert!(n >= 1, "hyp_sum: no terms");
    let m = n.isqrt().max(1);
    let pow = powers(x, w, m);
    let pow_err = |l: u64| Natural::from(l.saturating_sub(1));
    let blocks = n.div_ceil(m);
    let mut acc = Integer::ZERO;
    let mut acc_err = Natural::ZERO;
    for i in (0..blocks).rev() {
        let start = m * i;
        let len = m.min(n - start);
        // c_l from l = len - 1 down to 0, and the numerator with its error.
        let mut c = Integer::ONE;
        let mut num = Integer::ZERO;
        let mut num_err = Natural::ZERO;
        for l in (0..len).rev() {
            if l + 1 < len {
                c = &c * &Integer::from(q(start + l + 1));
            }
            let p = &pow[usize::try_from(l).expect("small")];
            num += &(p * &c);
            num_err = &num_err + &(&pow_err(l) * &c.unsigned_abs());
        }
        let d = c;
        if i + 1 < blocks {
            // T = trunc(trunc(P_m acc / 2^w) / q(m i + m)).
            let qm = Integer::from(q(start + m));
            let prod = shr(&pow[usize::try_from(m).expect("small")] * &acc, w);
            let t = prod.div_rem_trunc(&qm).0;
            let prod_err =
                &(&acc_err + &(&pow_err(m) * &Natural::from(2u64))) + &Natural::from(2u64);
            let t_err = &div_ceil(&prod_err, &qm.unsigned_abs()) + &Natural::from(1u64);
            num += &t;
            num_err = &num_err + &t_err;
        }
        acc = num.div_rem_trunc(&d).0;
        acc_err = &div_ceil(&num_err, &d.unsigned_abs()) + &Natural::from(1u64);
    }
    (acc, acc_err)
}

/// `P_0 = 2^w, P_1 = X, …, P_m`, `P_l = trunc(P_(l-1) X / 2^w)` (error
/// `l - 1` units for `|x| <= 1`).
fn powers(x: &Integer, w: u64, m: u64) -> Vec<Integer> {
    let mut pow = Vec::with_capacity(usize::try_from(m + 1).expect("small"));
    pow.push(Integer::from_parts(false, Natural::from(1u64) << w));
    pow.push(x.clone());
    for _ in 2..=m {
        let last = pow.last().expect("non-empty");
        let next = shr(last * x, w);
        pow.push(next);
    }
    pow
}

/// `Σ_{j < n} (-y)^j / (2j + 1)` for `y = Y / 2^w` with `0 <= y <= 1/4`,
/// `n >= 1` (`atan(x) / x` in `y = x^2`): `(S, e)` as for [`hyp_sum`].
///
/// The coefficients do not telescope into products, so each block is summed
/// directly, `Σ_l (-1)^j trunc(P_l / (2j + 1))` for `j = m i + l`, and
/// Horner's rule joins the blocks, `acc_i = block_i + trunc(P_m acc_(i+1) /
/// 2^w)`: still `m + n/m` full products. Errors: a term `e_l / (2j+1) + 1 <=
/// e_l + 1` (`e_l = l - 1` as before), the join `e_A + 2 e_m + 2` (`|acc| <=
/// Σ y^k <= 2`), all added.
pub(super) fn atan_sum(y: &Integer, w: u64, n: u64) -> (Integer, Natural) {
    assert!(n >= 1, "atan_sum: no terms");
    let m = n.isqrt().max(1);
    let pow = powers(y, w, m);
    let pow_err = |l: u64| Natural::from(l.saturating_sub(1));
    let blocks = n.div_ceil(m);
    let mut acc = Integer::ZERO;
    let mut acc_err = Natural::ZERO;
    for i in (0..blocks).rev() {
        let start = m * i;
        let len = m.min(n - start);
        let mut sum = Integer::ZERO;
        let mut err = Natural::ZERO;
        for l in 0..len {
            let j = start + l;
            let d = Integer::from(2 * j + 1);
            let term = pow[usize::try_from(l).expect("small")].div_rem_trunc(&d).0;
            if j % 2 == 0 {
                sum += &term;
            } else {
                sum -= &term;
            }
            err = &(&err + &pow_err(l)) + &Natural::from(1u64);
        }
        if i + 1 < blocks {
            let join = shr(&pow[usize::try_from(m).expect("small")] * &acc, w);
            sum += &join;
            err = &(&(&err + &acc_err) + &(&pow_err(m) * &Natural::from(2u64)))
                + &Natural::from(2u64);
        }
        acc = sum;
        acc_err = err;
    }
    (acc, acc_err)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::{atan_sum, exp_sum, hyp_sum};
    use crate::{Integer, Natural, Rational};

    /// `Σ_{j < n} x^j / Π_{t <= j} q(t)` exactly.
    fn exact(x: &Rational, n: u64, q: impl Fn(u64) -> i64) -> Rational {
        let (mut s, mut t) = (Rational::from(0i64), Rational::from(1i64));
        for j in 0..n {
            if j > 0 {
                t = &(&t * x) / &Rational::from(q(j));
            }
            s = &s + &t;
        }
        s
    }

    fn exact_atan(y: &Rational, n: u64) -> Rational {
        let (mut s, mut p) = (Rational::from(0i64), Rational::from(1i64));
        for j in 0..n {
            let term = &p / &Rational::from(i64::try_from(2 * j + 1).unwrap());
            s = if j % 2 == 0 { &s + &term } else { &s - &term };
            p = &p * y;
        }
        s
    }

    fn scale(w: u64) -> Rational {
        Rational::from(Integer::from_parts(false, Natural::from(1u64) << w))
    }

    /// Checks `(s, e)` against the exact value and a loose tightness bound.
    fn check(
        s: Integer,
        e: &Natural,
        want: &Rational,
        w: u64,
        loose: u64,
    ) -> Result<(), TestCaseError> {
        let got = &Rational::from(s) / &scale(w);
        let diff = want - &got;
        let diff = if diff < Rational::from(0i64) {
            -diff
        } else {
            diff
        };
        let bound = &Rational::from(Integer::from(e.clone())) / &scale(w);
        prop_assert!(diff <= bound, "error {diff} above {bound}");
        prop_assert!(
            Integer::from(e.clone()) <= Integer::from(loose),
            "loose bound {e:?}"
        );
        Ok(())
    }

    /// `X` with `|X / 2^w| <= limit / 4` (limit 2: 1/2; 1: 1/4).
    fn fixed(w: u64, num: i64, neg: bool, quarter: bool) -> Integer {
        let k = if quarter { w - 2 } else { w - 1 };
        let span = &Integer::from_parts(false, Natural::from(1u64) << k) + &Integer::ONE;
        let x = Integer::from(num).div_rem_euclid(&span).1;
        if neg { -x } else { x }
    }

    proptest! {
        #[test]
        fn hyp_within_the_bound(w in 8u64..400, n in 1u64..80, num in any::<i64>(), neg in any::<bool>(), kind in 0u8..3) {
            let x = fixed(w, num, neg, false);
            let xr = &Rational::from(x.clone()) / &scale(w);
            let (q, s, e): (fn(u64) -> i64, _, _) = match kind {
                0 => { let (s, e) = exp_sum(&x, w, n); (|t| i64::try_from(t).unwrap(), s, e) }
                1 => { let q = |t: u64| -i64::try_from((2 * t - 1) * 2 * t).unwrap(); let (s, e) = hyp_sum(&x, w, n, q); (q, s, e) }
                _ => { let q = |t: u64| -i64::try_from(2 * t * (2 * t + 1)).unwrap(); let (s, e) = hyp_sum(&x, w, n, q); (q, s, e) }
            };
            // A few units per block.
            check(s, &e, &exact(&xr, n, q), w, 4 * n + 8)?;
        }

        #[test]
        fn atan_within_the_bound(w in 8u64..400, n in 1u64..80, num in any::<i64>()) {
            let y = fixed(w, num, false, true);
            let yr = &Rational::from(y.clone()) / &scale(w);
            let (s, e) = atan_sum(&y, w, n);
            // About m^2 / 2 per block of m terms, so n √n / 2 in all.
            check(s, &e, &exact_atan(&yr, n), w, n * (n.isqrt() + 4) + 16)?;
        }
    }
}
