//! The residues of a rational function that are quadratic irrationals
//! (CAS-DOC-9 §3, CAS-197). For an irreducible quadratic factor
//! `r(t) = t² + β t + γ` of the residue resultant, `Δ = β² − 4γ` its
//! discriminant, the residues are `α = u ± v √Δ` with `u = −β/2` and
//! `v = 1/2`, and their part of the antiderivative is `α log S + ᾱ log S̄`
//! with `S = gcd(E, P − α E′)` over ℚ(√Δ) [Tra76] and `S̄` its conjugate:
//! - `Δ > 0`, real residues: `u log|N| + v √Δ log|S / S̄|`, `N = S S̄`, the
//!   two logarithms joined as giac joins them, without the absolute values
//!   when `N` has no real root: `S` and `S̄` are then real, monic and of
//!   even degree, so positive;
//! - `Δ < 0`, complex residues `u ± i v √−Δ`: with `S = A + i B`, `A` and
//!   `B` real, `u log N + v √−Δ · i log((A + iB)/(A − iB))`, and Rioboo's
//!   conversion of the last into arctangents of polynomials, which have no
//!   poles [Bro05 ch. 2].

use cas_num::algebra::{EuclideanDomain, Field, Rationals, Ring};
use cas_num::{Integer, Natural, Rational};

use crate::dense::{Coefficients, Dense, DenseRing};

/// A number `a + b √d` of ℚ(√d), `d` known from where it is.
pub type Surd = (Rational, Rational);

/// ℚ(√d) for a rational `d` that is not a square.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct Quadratic {
    d: Rational,
}

impl Quadratic {
    pub(crate) fn new(d: Rational) -> Quadratic {
        Quadratic { d }
    }
}

impl Ring for Quadratic {
    type Elem = Surd;

    fn zero(&self) -> Surd {
        (Rational::ZERO, Rational::ZERO)
    }

    fn one(&self) -> Surd {
        (Rational::ONE, Rational::ZERO)
    }

    fn add(&self, a: &Surd, b: &Surd) -> Surd {
        (&a.0 + &b.0, &a.1 + &b.1)
    }

    fn sub(&self, a: &Surd, b: &Surd) -> Surd {
        (&a.0 - &b.0, &a.1 - &b.1)
    }

    fn neg(&self, a: &Surd) -> Surd {
        (-&a.0, -&a.1)
    }

    fn mul(&self, a: &Surd, b: &Surd) -> Surd {
        (
            &(&a.0 * &b.0) + &(&self.d * &(&a.1 * &b.1)),
            &(&a.0 * &b.1) + &(&a.1 * &b.0),
        )
    }

    fn elem_integer(&self, n: &Integer) -> Surd {
        (Rational::from(n.clone()), Rational::ZERO)
    }

    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl EuclideanDomain for Quadratic {
    fn div_rem(&self, a: &Surd, b: &Surd) -> (Surd, Surd) {
        let inverse = self.inv(b).expect("division by zero");
        (self.mul(a, &inverse), self.zero())
    }

    fn normal_form(&self, a: &Surd) -> (Surd, Surd) {
        match self.inv(a) {
            Some(inverse) => (self.one(), inverse),
            None => (self.zero(), self.one()),
        }
    }
}

impl Field for Quadratic {
    /// `1 / (a + b √d) = (a − b √d) / (a² − d b²)`, the norm 0 only for 0
    /// as `d` is not a square.
    fn inv(&self, x: &Surd) -> Option<Surd> {
        let norm = &(&x.0 * &x.0) - &(&self.d * &(&x.1 * &x.1));
        let inverse = norm.checked_recip()?;
        Some((&x.0 * &inverse, -&(&x.1 * &inverse)))
    }
}

impl Coefficients for Quadratic {}

/// A pair of conjugate residues `u ± v √d` of a rational function, `d` not
/// a square: their part of the antiderivative.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Pair {
    /// `d > 0`: `u log|N| + v √d log|S / S̄|`, `S` monic over ℚ(√d), `S̄`
    /// its conjugate, and `N = S S̄` over ℚ, with whether `N` has a real
    /// root: if not, `N`, `S` and `S̄` are positive.
    Real {
        d: Rational,
        u: Rational,
        v: Rational,
        s: Vec<Surd>,
        norm: Dense<Rational>,
        real_root: bool,
    },
    /// `d < 0`: `u log N + v √−d Σⱼ 2 atan(Qⱼ)`, `N = S S̄ > 0` over ℚ and
    /// the polynomials `Qⱼ` over ℚ(√−d), which is ℚ when `−d` is a square.
    Complex {
        d: Rational,
        u: Rational,
        v: Rational,
        norm: Dense<Rational>,
        atans: Vec<Vec<Surd>>,
    },
}

impl Pair {
    /// `N = S S̄`, the factor of the denominator at whose roots the
    /// residues are.
    #[must_use]
    pub fn norm(&self) -> &Dense<Rational> {
        match self {
            Pair::Real { norm, .. } | Pair::Complex { norm, .. } => norm,
        }
    }
}

/// The part of `∫ P / E`, `E` square-free and `deg P < deg E`, at the roots
/// of `E` whose residues are the roots of `r`, monic and irreducible of
/// degree 2 over ℚ; `None` if no root has them.
pub(super) fn pair(
    ring: &DenseRing<Rationals>,
    p: &Dense<Rational>,
    e: &Dense<Rational>,
    r: &Dense<Rational>,
) -> Option<Pair> {
    let (gamma, beta) = (&r.coeffs()[0], &r.coeffs()[1]);
    let d = &(beta * beta) - &(&Rational::from(Integer::from(4)) * gamma);
    let kx = DenseRing::new(Quadratic::new(d.clone()));
    let lift = |q: &Dense<Rational>| {
        kx.poly(
            q.coeffs()
                .iter()
                .map(|c| (c.clone(), Rational::ZERO))
                .collect(),
        )
    };
    let half = Rational::new(Integer::ONE, Integer::from(2));
    let (u, v) = (-&(beta * &half), half);
    let alpha: Surd = (u.clone(), v.clone());
    let second = kx.sub(&lift(p), &kx.scale(&lift(&ring.derivative(e)), &alpha));
    let s = EuclideanDomain::gcd(&kx, &lift(e), &second);
    if s.degree().unwrap_or(0) == 0 {
        return None;
    }
    let conjugate = kx.poly(s.coeffs().iter().map(|(a, b)| (a.clone(), -b)).collect());
    let norm = ring.poly(
        kx.mul(&s, &conjugate)
            .coeffs()
            .iter()
            .map(|(a, _)| a.clone())
            .collect(),
    );
    if d.is_positive() {
        let real_root = ring.count_real_roots(&norm, None, None) > 0;
        return Some(Pair::Real {
            d,
            u,
            v,
            s: s.coeffs().to_vec(),
            norm,
            real_root,
        });
    }
    // S = A + √d B′ = A + i √−d B′, B = √−d B′: over ℚ when √−d is
    // rational, over the field ℚ(√−d) otherwise.
    let atans = if let Some(q) = rational_sqrt(&-&d) {
        let a = ring.poly(s.coeffs().iter().map(|(a, _)| a.clone()).collect());
        let b = ring.poly(s.coeffs().iter().map(|(_, b)| b * &q).collect());
        log_to_atan(ring, &a, &b)
            .into_iter()
            .map(|p| {
                p.coeffs()
                    .iter()
                    .map(|c| (c.clone(), Rational::ZERO))
                    .collect()
            })
            .collect()
    } else {
        let rx = DenseRing::new(Quadratic::new(-&d));
        let a = rx.poly(
            s.coeffs()
                .iter()
                .map(|(a, _)| (a.clone(), Rational::ZERO))
                .collect(),
        );
        let b = rx.poly(
            s.coeffs()
                .iter()
                .map(|(_, b)| (Rational::ZERO, b.clone()))
                .collect(),
        );
        log_to_atan(&rx, &a, &b)
            .into_iter()
            .map(|p| p.coeffs().to_vec())
            .collect()
    };
    Some(Pair::Complex {
        d,
        u,
        v,
        norm,
        atans,
    })
}

/// The rational square root of `r`, if it has one.
pub(crate) fn rational_sqrt(r: &Rational) -> Option<Rational> {
    let (n, d) = (r.numer(), r.denom());
    (!r.is_negative() && n.is_perfect_square() && d.is_perfect_square())
        .then(|| Rational::new(n.isqrt(), d.isqrt()))
}

/// Rioboo's conversion [Bro05 ch. 2]: polynomials `Qⱼ` such that
/// `i log((A + iB)/(A − iB))` and `Σⱼ 2 atan(Qⱼ)` have the same derivative,
/// for coprime real `A` and `B ≠ 0`. When `B` divides `A`, `2 atan(A/B)`.
/// Otherwise, with `B D − A C = 1`, `(A + iB)(D − iC) = (A D + B C) + i`,
/// so the logarithm is that of `((A D + B C) + i)/((A D + B C) − i)`, which
/// is `2 atan(A D + B C)`, and that of `(D + iC)/(D − iC)`, of lower
/// degrees.
fn log_to_atan<F: Field + Coefficients>(
    rx: &DenseRing<F>,
    a: &Dense<F::Elem>,
    b: &Dense<F::Elem>,
) -> Vec<Dense<F::Elem>> {
    let (q, r) = rx.div_rem(a, b);
    if r.is_zero() {
        return vec![q];
    }
    if a.degree() < b.degree() {
        // (−B + iA)/(−B − iA) is −(A + iB)/(A − iB): a constant apart.
        return log_to_atan(rx, &rx.neg(b), a);
    }
    let (g, d, c) = rx.gcd_ext(b, &rx.neg(a));
    debug_assert!(rx.is_one(&g), "coprime");
    let mut atans = vec![rx.add(&rx.mul(a, &d), &rx.mul(b, &c))];
    atans.extend(log_to_atan(rx, &d, &c));
    atans
}
