//! Integration of rational functions over ℚ (CAS-DOC-9 §3, CAS-196): the
//! polynomial part; Hermite's reduction in Mack's linear form [Her1872,
//! Bro98 §1.2], which leaves a square-free denominator; and the logarithms
//! whose residues are rational, from the factors of the residue resultant
//! `R(t) = res_x(E, P − t E′)` [Tra76]: for each rational root `c` of `R`,
//! `c log gcd(E, P − c E′)`. The part whose residues are not rational comes
//! back as it is, for the logarithmic part of Lazard, Rioboo and Trager
//! (CAS-197).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithms: A/D, P/E, b u + c v = w, s and t"
)]

#[cfg(test)]
mod tests;

use cas_num::Rational;
use cas_num::algebra::{EuclideanDomain, Field, Rationals, Ring};

use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::dense::{Dense, DenseRing};
use crate::sparse::SparseRing;

/// A rational function, a numerator over a denominator.
pub type Ratio = (Dense<Rational>, Dense<Rational>);

/// A logarithm `c log S` of an antiderivative.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Logarithm {
    /// The coefficient, a residue of the integrand.
    pub coefficient: Rational,
    /// The argument, monic and square-free.
    pub argument: Dense<Rational>,
    /// Whether the argument has a real root: a real antiderivative then
    /// takes the logarithm of its absolute value, `ln|x − 1|`, and otherwise
    /// of the argument, which is positive, `ln(x² + 1)`.
    pub real_root: bool,
}

/// `∫ A / D` for a rational function over ℚ: the integrated polynomial
/// part, Hermite's rational part, the logarithms with rational residues,
/// and what is left.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RationalIntegral {
    /// The antiderivative of the polynomial part.
    pub polynomial: Dense<Rational>,
    /// Hermite's rational part, a numerator over a monic denominator, in
    /// lowest terms.
    pub rational: Ratio,
    /// The logarithms whose residues are rational.
    pub logarithms: Vec<Logarithm>,
    /// What is left to integrate, `P / E` with `deg P < deg E` and `E`
    /// square-free, whose residues are not rational; `None` when nothing
    /// is.
    pub rest: Option<Ratio>,
}

/// `a / b`, which is exact.
fn exact(ring: &DenseRing<Rationals>, a: &Dense<Rational>, b: &Dense<Rational>) -> Dense<Rational> {
    ring.div_exact(a, b).expect("an exact quotient")
}

/// The antiderivative of `p` without a constant term.
fn antiderivative(ring: &DenseRing<Rationals>, p: &Dense<Rational>) -> Dense<Rational> {
    let mut coeffs = vec![Rational::ZERO];
    coeffs.extend(p.coeffs().iter().enumerate().map(|(k, c)| {
        let k = cas_num::Integer::from(k + 1);
        c / &Rational::from(k)
    }));
    ring.poly(coeffs)
}

/// `(b, c)` with `b u + c v = w` and `deg b < deg v`, for coprime `u` and
/// `v`.
fn solve(
    ring: &DenseRing<Rationals>,
    u: &Dense<Rational>,
    v: &Dense<Rational>,
    w: &Dense<Rational>,
) -> (Dense<Rational>, Dense<Rational>) {
    let (g, s, _) = ring.gcd_ext(u, v);
    debug_assert!(ring.is_one(&g), "coprime");
    let b = ring.div_rem(&ring.mul(w, &s), v).1;
    let c = exact(ring, &ring.sub(w, &ring.mul(&b, u)), v);
    (b, c)
}

/// Hermite's reduction in Mack's linear form [Bro98 §1.2, Bro05 §2.2]:
/// `∫ A / D = G + ∫ H / D*`, `G` a fraction and `D*` the square-free part
/// of `D`. With `D⁻ = gcd(D, D′)`, each round solves
/// `B (−D* D⁻′ / D⁻) + C D⁻* = A` for `D⁻* = D⁻ / gcd(D⁻, D⁻′)`, adds
/// `B / D⁻` to `G`, and goes on with `A = C − B′ D* / D⁻*` and
/// `D⁻ = gcd(D⁻, D⁻′)`.
fn hermite(
    ring: &DenseRing<Rationals>,
    a: &Dense<Rational>,
    d: &Dense<Rational>,
    cancel: &dyn Cancel,
) -> Result<(Ratio, Ratio), Cancelled> {
    let mut a = a.clone();
    let (mut gn, mut gd) = (ring.zero(), ring.one());
    let mut dm = ring.gcd(d, &ring.derivative(d));
    let ds = exact(ring, d, &dm);
    while dm.degree().unwrap_or(0) > 0 {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let dm2 = ring.gcd(&dm, &ring.derivative(&dm));
        let dms = exact(ring, &dm, &dm2);
        let u = ring.neg(&exact(ring, &ring.mul(&ds, &ring.derivative(&dm)), &dm));
        let (b, c) = solve(ring, &u, &dms, &a);
        a = ring.sub(&c, &exact(ring, &ring.mul(&ring.derivative(&b), &ds), &dms));
        gn = ring.add(&ring.mul(&gn, &dm), &ring.mul(&b, &gd));
        gd = ring.mul(&gd, &dm);
        let common = ring.gcd(&gn, &gd);
        if common.degree().unwrap_or(0) > 0 {
            gn = exact(ring, &gn, &common);
            gd = exact(ring, &gd, &common);
        }
        dm = dm2;
    }
    // A monic denominator.
    if let Some(lead) = gd.lead().cloned() {
        let inverse = ring.base().inv(&lead).expect("non-zero");
        gn = ring.scale(&gn, &inverse);
        gd = ring.scale(&gd, &inverse);
    }
    Ok(((gn, gd), (a, ds)))
}

/// The number of real roots of a square-free `p`, by Sturm's sequence: its
/// sign changes at −∞ less those at +∞.
pub(crate) fn real_roots(ring: &DenseRing<Rationals>, p: &Dense<Rational>) -> usize {
    let mut sequence = vec![p.clone(), ring.derivative(p)];
    loop {
        let n = sequence.len();
        if sequence[n - 1].is_zero() {
            sequence.pop();
            break;
        }
        let r = ring.div_rem(&sequence[n - 2], &sequence[n - 1]).1;
        if r.is_zero() {
            break;
        }
        sequence.push(ring.neg(&r));
    }
    let changes = |signs: Vec<bool>| signs.windows(2).filter(|w| w[0] != w[1]).count();
    let negative = |q: &Dense<Rational>| q.lead().is_some_and(Rational::is_negative);
    let at_plus = changes(sequence.iter().map(negative).collect());
    let at_minus = changes(
        sequence
            .iter()
            .map(|q| negative(q) != (q.degree().unwrap_or(0) % 2 == 1))
            .collect(),
    );
    at_minus - at_plus
}

/// The logarithms of `∫ P / E`, `deg P < deg E` and `E` square-free, whose
/// residues are rational, and the part left [Tra76]. With `E_c` the product
/// of their arguments and `E_r = E / E_c`, `s E_r + t E_c = 1`, the part
/// left is `(P t mod E_r) / E_r`.
fn logarithms(
    ring: &DenseRing<Rationals>,
    p: &Dense<Rational>,
    e: &Dense<Rational>,
    cancel: &dyn Cancel,
) -> Result<(Vec<Logarithm>, Option<Ratio>), Cancelled> {
    if p.is_zero() {
        return Ok((Vec::new(), None));
    }
    // R(t) = res_x(E, P − t E′), x the variable 0 and t the variable 1.
    let two = SparseRing::new(Rationals, 2);
    let lift = |q: &Dense<Rational>| two.from_dense(q, 0);
    let second = two.sub(&lift(p), &two.mul(&two.var(1), &lift(&ring.derivative(e))));
    let r = two.resultant_with(&lift(e), &second, 0, cancel)?;
    let r = two.to_dense(&r, 1).expect("a polynomial in t");
    let factors = ring.factor_with(&r, cancel)?;
    let mut logarithms = Vec::new();
    let mut rational = ring.one();
    for (f, _) in factors.factors() {
        if f.degree() != Some(1) {
            continue;
        }
        // f = t + f₀, monic: the residue −f₀.
        let c = -&f.coeffs()[0];
        let s = ring.gcd(e, &ring.sub(p, &ring.scale(&ring.derivative(e), &c)));
        if s.degree().unwrap_or(0) == 0 {
            continue;
        }
        let real_root = real_roots(ring, &s) > 0;
        rational = ring.mul(&rational, &s);
        logarithms.push(Logarithm {
            coefficient: c,
            argument: s,
            real_root,
        });
    }
    let others = exact(ring, e, &rational);
    if others.degree().unwrap_or(0) == 0 {
        return Ok((logarithms, None));
    }
    let (_, _, t) = ring.gcd_ext(&others, &rational);
    let rest = ring.div_rem(&ring.mul(p, &t), &others).1;
    Ok((logarithms, Some((rest, others))))
}

impl DenseRing<Rationals> {
    /// `∫ a / d` for `d` not zero, as far as rational residues go
    /// (CAS-DOC-9 §3): the integrated polynomial part, Hermite's rational
    /// part, the logarithms with rational residues, and the part whose
    /// residues are not rational.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn integrate_rational(&self, a: &Dense<Rational>, d: &Dense<Rational>) -> RationalIntegral {
        uncancelled(self.integrate_rational_with(a, d, &Never))
    }

    /// As [`DenseRing::integrate_rational`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    pub fn integrate_rational_with(
        &self,
        a: &Dense<Rational>,
        d: &Dense<Rational>,
        cancel: &dyn Cancel,
    ) -> Result<RationalIntegral, Cancelled> {
        assert!(!d.is_zero(), "a zero denominator");
        let (q, r) = self.div_rem(a, d);
        let (rational, (h, e)) = hermite(self, &r, d, cancel)?;
        // Hermite's remainder can be improper.
        let (q2, p) = self.div_rem(&h, &e);
        let polynomial = antiderivative(self, &self.add(&q, &q2));
        let (logarithms, rest) = logarithms(self, &p, &e, cancel)?;
        Ok(RationalIntegral {
            polynomial,
            rational,
            logarithms,
            rest,
        })
    }
}
