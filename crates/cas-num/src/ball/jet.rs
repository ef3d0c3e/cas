//! Jets of balls (CAS-DOC-35 §2, CAS-389): truncated Taylor series `Σ_{k<n}
//! c_k s^k` with ball coefficients, the Taylor coefficients of a function at
//! a point to the order `n`, and the elementary functions of them.
//!
//! Each function follows the recurrence its derivative gives (Moore;
//! \[Joh17\]), on the coefficients of `u = Σ u_k s^k`:
//!
//! - `v = e^u`: `v_0 = e^(u_0)` and `k v_k = Σ_{j=1}^{k} j u_j v_(k−j)`,
//!   from `v′ = u′ v`; `sinh` and `cosh` by the pair `v = (cosh u, sinh u)`,
//!   `sin` and `cos` by the same with a sign.
//! - `w = u/v`: `w_k = (u_k − Σ_{j=1}^{k} v_j w_(k−j))/v_0`.
//! - `w = √u`: `w_0 = √(u_0)` and `w_k = (u_k − Σ_{j=1}^{k−1} w_j
//!   w_(k−j))/(2 w_0)`, from `w² = u`.
//! - `v = u^a` for a real `a`: `v_0 = u_0^a` and `k u_0 v_k = Σ_{j=1}^{k}
//!   ((a + 1) j − k) u_j v_(k−j)`, from `u v′ = a u′ v`; an integer power
//!   by squaring.
//! - `log`, `atan`, `asin`, `acos`, `asinh`, `acosh` and `atanh` as the
//!   integral of their derivatives' jets, `u′/u`, `u′/(1 + u²)`, `u′/√(1 −
//!   u²)`, … from their values at `u_0`.
//! - `tan` and `tanh` as quotients.
//!
//! The identities hold of formal power series, and ball arithmetic encloses
//! each coefficient: a jet holds the Taylor coefficients of the function of
//! every series whose coefficients its argument's balls hold.

use crate::Integer;
use crate::ball::Ball;

/// A truncated Taylor series with ball coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct Jet {
    coeffs: Vec<Ball>,
}

fn small(k: usize) -> Ball {
    Ball::from_integer(&Integer::from(u64::try_from(k).expect("a small index")))
}

impl Jet {
    /// The jet of the coefficients `coeffs`, the constant first.
    ///
    /// # Panics
    ///
    /// If `coeffs` is empty.
    #[must_use]
    pub fn new(coeffs: Vec<Ball>) -> Jet {
        assert!(!coeffs.is_empty(), "Jet::new: no coefficient");
        Jet { coeffs }
    }

    /// The constant `c` to the order `n`.
    #[must_use]
    pub fn constant(c: Ball, n: usize) -> Jet {
        let mut coeffs = vec![Ball::ZERO; n.max(1)];
        coeffs[0] = c;
        Jet { coeffs }
    }

    /// `c + s`, the variable at `c`, to the order `n`.
    #[must_use]
    pub fn variable(c: Ball, n: usize) -> Jet {
        let mut jet = Jet::constant(c, n);
        if let Some(one) = jet.coeffs.get_mut(1) {
            *one = Ball::from_integer(&Integer::ONE);
        }
        jet
    }

    /// The coefficients, the constant first.
    #[must_use]
    pub fn coeffs(&self) -> &[Ball] {
        &self.coeffs
    }

    /// The order: the number of coefficients.
    #[must_use]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Never: a jet has its constant.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The common order of two jets.
    fn order(&self, other: &Jet) -> usize {
        self.len().min(other.len())
    }

    /// `self + other`.
    #[must_use]
    pub fn add(&self, other: &Jet, prec: u64) -> Jet {
        let n = self.order(other);
        Jet::new(
            (0..n)
                .map(|k| self.coeffs[k].add(&other.coeffs[k], prec))
                .collect(),
        )
    }

    /// `self − other`.
    #[must_use]
    pub fn sub(&self, other: &Jet, prec: u64) -> Jet {
        let n = self.order(other);
        Jet::new(
            (0..n)
                .map(|k| self.coeffs[k].sub(&other.coeffs[k], prec))
                .collect(),
        )
    }

    /// `−self`.
    #[must_use]
    pub fn neg(&self) -> Jet {
        Jet::new(self.coeffs.iter().map(Ball::neg).collect())
    }

    /// `c self` for a ball `c`.
    #[must_use]
    pub fn scale(&self, c: &Ball, prec: u64) -> Jet {
        Jet::new(self.coeffs.iter().map(|a| a.mul(c, prec)).collect())
    }

    /// `self + c` for a ball `c`.
    #[must_use]
    pub fn shift(&self, c: &Ball, prec: u64) -> Jet {
        let mut jet = self.clone();
        jet.coeffs[0] = jet.coeffs[0].add(c, prec);
        jet
    }

    /// `self · other`, Cauchy's product.
    #[must_use]
    pub fn mul(&self, other: &Jet, prec: u64) -> Jet {
        let n = self.order(other);
        let coeffs = (0..n)
            .map(|k| {
                (0..=k).fold(Ball::ZERO, |acc, j| {
                    acc.add(&self.coeffs[j].mul(&other.coeffs[k - j], prec), prec)
                })
            })
            .collect();
        Jet::new(coeffs)
    }

    /// `self / other`; indeterminate coefficients where `other`'s constant
    /// meets 0.
    #[must_use]
    pub fn div(&self, other: &Jet, prec: u64) -> Jet {
        let n = self.order(other);
        let mut w: Vec<Ball> = Vec::with_capacity(n);
        for k in 0..n {
            let mut acc = self.coeffs[k].clone();
            for j in 1..=k {
                acc = acc.sub(&other.coeffs[j].mul(&w[k - j], prec), prec);
            }
            w.push(acc.div(&other.coeffs[0], prec));
        }
        Jet::new(w)
    }

    /// The derivative in `s`: one coefficient shorter, at least the
    /// constant 0.
    #[must_use]
    pub fn derivative(&self, prec: u64) -> Jet {
        if self.len() == 1 {
            return Jet::constant(Ball::ZERO, 1);
        }
        Jet::new(
            (1..self.len())
                .map(|k| self.coeffs[k].mul(&small(k), prec))
                .collect(),
        )
    }

    /// The integral in `s` from `c` at 0: one coefficient longer.
    #[must_use]
    pub fn integral(&self, c: Ball, prec: u64) -> Jet {
        let mut coeffs = Vec::with_capacity(self.len() + 1);
        coeffs.push(c);
        for (k, a) in self.coeffs.iter().enumerate() {
            coeffs.push(a.div(&small(k + 1), prec));
        }
        Jet::new(coeffs)
    }

    /// `v` of `k v_k = Σ_{j=1}^{k} j u_j w_(k−j)` and `w` of the same with
    /// `v` for the pair `(v, w)` = `(cosh u, sinh u)`, or `(cos u, −…)` for
    /// `sin` and `cos` (`trig`).
    fn pair(&self, v0: Ball, w0: Ball, trig: bool, prec: u64) -> (Jet, Jet) {
        let n = self.len();
        let (mut v, mut w) = (vec![v0], vec![w0]);
        for k in 1..n {
            let (mut dv, mut dw) = (Ball::ZERO, Ball::ZERO);
            for j in 1..=k {
                let ju = self.coeffs[j].mul(&small(j), prec);
                dv = dv.add(&ju.mul(&w[k - j], prec), prec);
                dw = dw.add(&ju.mul(&v[k - j], prec), prec);
            }
            let kb = small(k);
            // (cosh u)′ = u′ sinh u, (sinh u)′ = u′ cosh u; (cos u)′ = −u′
            // sin u, (sin u)′ = u′ cos u.
            let next_v = dv.div(&kb, prec);
            v.push(if trig { next_v.neg() } else { next_v });
            w.push(dw.div(&kb, prec));
        }
        (Jet::new(v), Jet::new(w))
    }

    /// `e^self`.
    #[must_use]
    pub fn exp(&self, prec: u64) -> Jet {
        let n = self.len();
        let mut v = vec![self.coeffs[0].exp(prec)];
        for k in 1..n {
            let mut acc = Ball::ZERO;
            for j in 1..=k {
                acc = acc.add(
                    &self.coeffs[j].mul(&small(j), prec).mul(&v[k - j], prec),
                    prec,
                );
            }
            v.push(acc.div(&small(k), prec));
        }
        Jet::new(v)
    }

    /// `(sin self, cos self)`.
    #[must_use]
    pub fn sin_cos(&self, prec: u64) -> (Jet, Jet) {
        let (s, c) = self.coeffs[0].sin_cos(prec);
        let (cos, sin) = self.pair(c, s, true, prec);
        (sin, cos)
    }

    /// `(sinh self, cosh self)`.
    #[must_use]
    pub fn sinh_cosh(&self, prec: u64) -> (Jet, Jet) {
        let u0 = &self.coeffs[0];
        let (up, down) = (u0.exp(prec), u0.neg().exp(prec));
        let half = |b: Ball| b.mul_2exp(&-Integer::ONE);
        let sinh = half(up.sub(&down, prec));
        let cosh = half(up.add(&down, prec));
        let (cosh, sinh) = self.pair(cosh, sinh, false, prec);
        (sinh, cosh)
    }

    /// `tan self`.
    #[must_use]
    pub fn tan(&self, prec: u64) -> Jet {
        let (s, c) = self.sin_cos(prec);
        s.div(&c, prec)
    }

    /// `tanh self`.
    #[must_use]
    pub fn tanh(&self, prec: u64) -> Jet {
        let (s, c) = self.sinh_cosh(prec);
        s.div(&c, prec)
    }

    /// `√self`.
    #[must_use]
    pub fn sqrt(&self, prec: u64) -> Jet {
        let n = self.len();
        let w0 = self.coeffs[0].sqrt(prec);
        let twice = w0.mul_2exp(&Integer::ONE);
        let mut w = vec![w0];
        for k in 1..n {
            let mut acc = self.coeffs[k].clone();
            for j in 1..k {
                acc = acc.sub(&w[j].mul(&w[k - j], prec), prec);
            }
            w.push(acc.div(&twice, prec));
        }
        Jet::new(w)
    }

    /// `self^a` for a real ball `a`; `self`'s constant must be positive.
    #[must_use]
    pub fn pow(&self, a: &Ball, prec: u64) -> Jet {
        let n = self.len();
        let u0 = &self.coeffs[0];
        let a1 = a.add(&Ball::from_integer(&Integer::ONE), prec);
        let mut v = vec![u0.log(prec).mul(a, prec).exp(prec)];
        for k in 1..n {
            let mut acc = Ball::ZERO;
            for j in 1..=k {
                let weight = a1.mul(&small(j), prec).sub(&small(k), prec);
                acc = acc.add(
                    &weight.mul(&self.coeffs[j], prec).mul(&v[k - j], prec),
                    prec,
                );
            }
            v.push(acc.div(&u0.mul(&small(k), prec), prec));
        }
        Jet::new(v)
    }

    /// `self^m` for an integer `m`, by squaring; a negative one through
    /// the quotient.
    #[must_use]
    pub fn powi(&self, m: i64, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let mut base = self.clone();
        let mut e = m.unsigned_abs();
        let mut result = one.clone();
        while e > 0 {
            if e & 1 == 1 {
                result = result.mul(&base, prec);
            }
            e >>= 1;
            if e > 0 {
                base = base.mul(&base, prec);
            }
        }
        if m < 0 {
            one.div(&result, prec)
        } else {
            result
        }
    }

    /// The function whose value at `self`'s constant is `at0` and whose
    /// derivative is `derivative` times `self′`, integrated.
    fn integrated(&self, at0: Ball, derivative: &Jet, prec: u64) -> Jet {
        let n = self.len();
        if n == 1 {
            return Jet::constant(at0, 1);
        }
        let shorter = Jet::new(derivative.coeffs[..n - 1].to_vec());
        let du = self.derivative(prec);
        shorter.mul(&du, prec).integral(at0, prec)
    }

    /// `log self`; `self`'s constant must be positive.
    #[must_use]
    pub fn log(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let inverse = one.div(self, prec);
        self.integrated(self.coeffs[0].log(prec), &inverse, prec)
    }

    /// `1 ± self²`, the base of the inverse functions' derivatives.
    fn one_plus_square(&self, sign: i64, prec: u64) -> Jet {
        let square = self.mul(self, prec);
        let one = Ball::from_integer(&Integer::ONE);
        if sign > 0 {
            square.shift(&one, prec)
        } else {
            square.neg().shift(&one, prec)
        }
    }

    /// `atan self`.
    #[must_use]
    pub fn atan(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let d = one.div(&self.one_plus_square(1, prec), prec);
        self.integrated(self.coeffs[0].atan(prec), &d, prec)
    }

    /// `asin self`.
    #[must_use]
    pub fn asin(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let d = one.div(&self.one_plus_square(-1, prec).sqrt(prec), prec);
        self.integrated(self.coeffs[0].asin(prec), &d, prec)
    }

    /// `acos self`, `π/2 − asin self`.
    #[must_use]
    pub fn acos(&self, prec: u64) -> Jet {
        let half_pi = Ball::const_pi(prec).mul_2exp(&-Integer::ONE);
        self.asin(prec).neg().shift(&half_pi, prec)
    }

    /// `asinh self`, from `1/√(1 + u²)`.
    #[must_use]
    pub fn asinh(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let d = one.div(&self.one_plus_square(1, prec).sqrt(prec), prec);
        self.integrated(self.coeffs[0].asinh(prec), &d, prec)
    }

    /// `acosh self`, from `1/√(u² − 1)`.
    #[must_use]
    pub fn acosh(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let d = one.div(&self.one_plus_square(-1, prec).neg().sqrt(prec), prec);
        self.integrated(self.coeffs[0].acosh(prec), &d, prec)
    }

    /// `atanh self`, from `1/(1 − u²)`.
    #[must_use]
    pub fn atanh(&self, prec: u64) -> Jet {
        let one = Jet::constant(Ball::from_integer(&Integer::ONE), self.len());
        let d = one.div(&self.one_plus_square(-1, prec), prec);
        self.integrated(self.coeffs[0].atanh(prec), &d, prec)
    }
}

#[cfg(test)]
mod tests;
