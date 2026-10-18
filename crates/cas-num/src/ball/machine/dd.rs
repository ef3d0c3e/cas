//! Double-double ("double-word") arithmetic for the machine-ball kernels
//! (CAS-DOC-5 §5.4): a value `hi + lo` with `|lo| <= ulp(hi) / 2`, about 106
//! bits.
//!
//! The algorithms are those of Joldes, Muller and Popescu, "Tight and
//! rigorous error bounds for basic building blocks of double-word
//! arithmetic" (ACM TOMS 44, 2017), built on the error-free two-sum and
//! two-product. With `u = 2^-53`, their proved relative error bounds are
//! `3u²/(1 - 4u)` for [`Dd::add`] (Alg. 6), `1.5u² + 4u³` for [`Dd::mul_f64`]
//! (Alg. 7), `7u²` for [`Dd::mul`] (Alg. 10) and `15u² + 56u³` for
//! [`Dd::div`] (Alg. 17): all below [`OP_ERR`] `= 2^-102 = 16u²`, the one
//! bound the kernels use.
//!
//! The two-product is exact while nothing overflows or underflows: inputs
//! here stay within `2^-600 <= |x| <= 2^600`, which the kernels guarantee
//! (smaller and larger arguments go back to `Ball`).

/// A relative error bound for every operation here, `2^-102`.
#[allow(
    dead_code,
    reason = "the unit of the kernels' error budgets; checked by the tests"
)]
pub(super) const OP_ERR: f64 = f64::from_bits((1023 - 102) << 52);

/// Knuth's two-sum: `s + e = a + b` exactly.
#[inline]
pub(super) fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

/// Dekker's fast two-sum, for `|a| >= |b|` (or `a = 0`).
#[inline]
fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

/// `p + e = a b` exactly: a fused multiply-add where the target has one,
/// Veltkamp's split and Dekker's product otherwise.
#[inline]
pub(super) fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    if cfg!(target_feature = "fma") {
        return (p, a.mul_add(b, -p));
    }
    let split = |x: f64| {
        let c = 134_217_729.0 * x;
        let h = c - (c - x);
        (h, x - h)
    };
    let ((ah, al), (bh, bl)) = (split(a), split(b));
    (p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
}

/// A double-double `hi + lo`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Dd {
    pub(super) hi: f64,
    pub(super) lo: f64,
}

impl Dd {
    pub(super) const ZERO: Dd = Dd { hi: 0.0, lo: 0.0 };
    pub(super) const ONE: Dd = Dd { hi: 1.0, lo: 0.0 };

    #[inline]
    pub(super) fn from_f64(x: f64) -> Dd {
        Dd { hi: x, lo: 0.0 }
    }

    /// The exact sum of two doubles.
    #[inline]
    pub(super) fn sum(a: f64, b: f64) -> Dd {
        let (hi, lo) = two_sum(a, b);
        Dd { hi, lo }
    }

    /// The exact product of two doubles.
    #[inline]
    pub(super) fn prod(a: f64, b: f64) -> Dd {
        let (hi, lo) = two_prod(a, b);
        Dd { hi, lo }
    }

    #[inline]
    pub(super) fn neg(self) -> Dd {
        Dd {
            hi: -self.hi,
            lo: -self.lo,
        }
    }

    /// `self + y` (JMP Alg. 6, `AccurateDWPlusDW`).
    #[inline]
    pub(super) fn add(self, y: Dd) -> Dd {
        let (sh, sl) = two_sum(self.hi, y.hi);
        let (th, tl) = two_sum(self.lo, y.lo);
        let c = sl + th;
        let (vh, vl) = fast_two_sum(sh, c);
        let w = tl + vl;
        let (hi, lo) = fast_two_sum(vh, w);
        Dd { hi, lo }
    }

    #[inline]
    pub(super) fn sub(self, y: Dd) -> Dd {
        self.add(y.neg())
    }

    /// `self × y` for a double `y` (JMP Alg. 7, `DWTimesFP1`).
    #[inline]
    pub(super) fn mul_f64(self, y: f64) -> Dd {
        let (ch, cl1) = two_prod(self.hi, y);
        let cl2 = self.lo * y;
        let (th, tl1) = fast_two_sum(ch, cl2);
        let tl2 = tl1 + cl1;
        let (hi, lo) = fast_two_sum(th, tl2);
        Dd { hi, lo }
    }

    /// `self × y` (JMP Alg. 10, `DWTimesDW1`).
    #[inline]
    pub(super) fn mul(self, y: Dd) -> Dd {
        let (ch, cl1) = two_prod(self.hi, y.hi);
        let tl1 = self.hi * y.lo;
        let tl2 = self.lo * y.hi;
        let cl2 = tl1 + tl2;
        let cl3 = cl1 + cl2;
        let (hi, lo) = fast_two_sum(ch, cl3);
        Dd { hi, lo }
    }

    /// `self / y` (JMP Alg. 17, `DWDivDW2`).
    #[inline]
    pub(super) fn div(self, y: Dd) -> Dd {
        let th = self.hi / y.hi;
        let r = y.mul_f64(th);
        let ph = self.hi - r.hi;
        let dl = self.lo - r.lo;
        let d = ph + dl;
        let tl = d / y.hi;
        let (hi, lo) = fast_two_sum(th, tl);
        Dd { hi, lo }
    }

    /// `self × 2^k`, exactly while the result stays normal.
    #[inline]
    pub(super) fn mul_pow2(self, k: i32) -> Dd {
        let s = f64::from_bits(u64::try_from(1023 + k).expect("normal scale") << 52);
        Dd {
            hi: self.hi * s,
            lo: self.lo * s,
        }
    }

    /// The nearest double (`hi + lo` rounded once).
    #[inline]
    pub(super) fn to_f64(self) -> f64 {
        self.hi + self.lo
    }
}

/// Horner's rule `Σ c_k t^k` in double-double, `c` from the constant term up.
///
/// With every operation within `OP_ERR` relative, the computed value is
/// within `γ_(2n) Σ |c_k| |t|^k` of the exact polynomial at the same `t`
/// (`γ_k = k ε / (1 - k ε)`, the classical Horner bound), `n` the degree.
pub(super) fn horner(c: &[Dd], t: Dd) -> Dd {
    let mut s = *c.last().expect("non-empty");
    for &ck in c.iter().rev().skip(1) {
        s = s.mul(t).add(ck);
    }
    s
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;
    use rug::Float as MpFloat;

    use super::{Dd, OP_ERR};

    fn exact(d: Dd) -> MpFloat {
        MpFloat::with_val(400, d.hi) + d.lo
    }

    /// A normalised double-double with |x| in roughly [2^-300, 2^300].
    fn dd() -> impl Strategy<Value = Dd> {
        (-300i32..300, 1.0f64..2.0, -1.0f64..1.0, any::<bool>()).prop_map(|(e, m, l, neg)| {
            let hi = m * 2f64.powi(e) * if neg { -1.0 } else { 1.0 };
            let lo = l * hi * 2f64.powi(-54);
            Dd::sum(hi, lo)
        })
    }

    fn rel_err(got: Dd, want: &MpFloat) -> f64 {
        let d = (exact(got) - want).abs() / want.clone().abs();
        d.to_f64()
    }

    proptest! {
        #[test]
        fn operations_within_the_bound(x in dd(), y in dd(), f in -1e300f64..1e300) {
            let (xe, ye) = (exact(x), exact(y));
            let sum = MpFloat::with_val(400, &xe + &ye);
            if !sum.is_zero() {
                prop_assert!(rel_err(x.add(y), &sum) <= OP_ERR, "add");
            }
            let prod = MpFloat::with_val(400, &xe * &ye);
            prop_assert!(rel_err(x.mul(y), &prod) <= OP_ERR, "mul");
            let quo = MpFloat::with_val(400, &xe / &ye);
            prop_assert!(rel_err(x.div(y), &quo) <= OP_ERR, "div");
            let g = if f == 0.0 { 1.0 } else { f.abs().max(1e-300) * f.signum() * 1e-200 };
            let pf = MpFloat::with_val(400, &xe * g);
            if !pf.is_zero() {
                prop_assert!(rel_err(x.mul_f64(g), &pf) <= OP_ERR, "mul_f64");
            }
            // Error-free building blocks.
            let s = Dd::sum(x.hi, y.hi);
            prop_assert_eq!(exact(s), MpFloat::with_val(400, x.hi) + y.hi);
            let p = Dd::prod(x.hi, y.hi);
            prop_assert_eq!(exact(p), MpFloat::with_val(400, x.hi) * y.hi);
        }
    }
}
