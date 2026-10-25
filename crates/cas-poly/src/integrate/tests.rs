//! Tests of the integration of rational functions (CAS-196, CAS-197):
//! examples, Sturm's count of real roots, and properties. The derivative of
//! the result plus the part left is the integrand, checked by values at
//! rational points, exactly, in the field of the residues (ℚ, or ℚ(√k) for
//! a pair of quadratic residues); denominators that split over ℚ leave
//! nothing, nor do those with quadratic factors.

use cas_num::algebra::{Field, Rationals, Ring};
use cas_num::{Integer, Rational};
use proptest::prelude::*;

use super::quadratic::{Pair, Quadratic, Surd, rational_sqrt};
use super::{RationalIntegral, real_roots};
use crate::dense::{Dense, DenseRing};

fn qx() -> DenseRing<Rationals> {
    DenseRing::new(Rationals)
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn poly(coeffs: &[i64]) -> Dense<Rational> {
    qx().poly(coeffs.iter().map(|&c| q(c, 1)).collect())
}

/// `x` in a field ℚ(√k).
fn surd(x: Rational) -> Surd {
    (x, Rational::ZERO)
}

/// `p(x)` and `p′(x)` for `p` over ℚ(√k), from its coefficients.
fn values(field: &Quadratic, p: &[Surd], x: &Surd) -> (Surd, Surd) {
    let mut value = field.zero();
    let mut derivative = field.zero();
    for c in p.iter().rev() {
        derivative = field.add(&field.mul(&derivative, x), &value);
        value = field.add(&field.mul(&value, x), c);
    }
    (value, derivative)
}

fn lifted(p: &Dense<Rational>) -> Vec<Surd> {
    p.coeffs().iter().cloned().map(surd).collect()
}

/// `p′(x) / p(x)`, `None` where `p(x)` is 0.
fn log_derivative(field: &Quadratic, p: &[Surd], x: &Surd) -> Option<Surd> {
    let (value, derivative) = values(field, p, x);
    Some(field.mul(&derivative, &field.inv(&value)?))
}

/// The field of the residues of `r`, ℚ(√k), with `√k` as an element when
/// it is rational: ℚ as ℚ(√2) when every residue is rational.
fn field_of(r: &RationalIntegral) -> (Quadratic, Option<Rational>) {
    match r.pairs.first() {
        Some(Pair::Real { d, .. }) => (Quadratic::new(d.clone()), None),
        Some(Pair::Complex { d, .. }) => match rational_sqrt(&-d) {
            Some(root) => (Quadratic::new(q(2, 1)), Some(root)),
            None => (Quadratic::new(-d), None),
        },
        None => (Quadratic::new(q(2, 1)), None),
    }
}

/// The derivative of the result at `x`, plus the part left: `None` at a
/// pole of a part.
fn derivative_at(r: &RationalIntegral, x: &Rational) -> Option<Surd> {
    let (field, root) = field_of(r);
    let at = surd(x.clone());
    let mut sum = values(&field, &lifted(&r.polynomial), &at).1;
    // (n / d)′ = n′/d − n d′/d².
    let (n, d) = &r.rational;
    let (nv, nd) = values(&field, &lifted(n), &at);
    let (dv, dd) = values(&field, &lifted(d), &at);
    let inverse = field.inv(&dv)?;
    let quotient = field.sub(
        &field.mul(&nd, &inverse),
        &field.mul(&field.mul(&nv, &dd), &field.mul(&inverse, &inverse)),
    );
    sum = field.add(&sum, &quotient);
    for log in &r.logarithms {
        let term = log_derivative(&field, &lifted(&log.argument), &at)?;
        sum = field.add(&sum, &field.mul(&surd(log.coefficient.clone()), &term));
    }
    for pair in &r.pairs {
        sum = field.add(&sum, &pair_derivative(&field, root.as_ref(), pair, &at)?);
    }
    if let Some((p, e)) = &r.rest {
        let (pv, _) = values(&field, &lifted(p), &at);
        let (ev, _) = values(&field, &lifted(e), &at);
        sum = field.add(&sum, &field.mul(&pv, &field.inv(&ev)?));
    }
    Some(sum)
}

/// The derivative of a pair's part at `x`, in `field`, ℚ(√d) for a real
/// pair, ℚ(√−d) for a complex one, or ℚ with `√−d = root`.
fn pair_derivative(
    field: &Quadratic,
    root: Option<&Rational>,
    pair: &Pair,
    x: &Surd,
) -> Option<Surd> {
    match pair {
        // u N′/N + v √d (S′/S − S̄′/S̄).
        Pair::Real { u, v, s, norm, .. } => {
            let conjugate: Vec<Surd> = s.iter().map(|(a, b)| (a.clone(), -b)).collect();
            let difference = field.sub(
                &log_derivative(field, s, x)?,
                &log_derivative(field, &conjugate, x)?,
            );
            let log = log_derivative(field, &lifted(norm), x)?;
            Some(field.add(
                &field.mul(&surd(u.clone()), &log),
                &field.mul(&(Rational::ZERO, v.clone()), &difference),
            ))
        }
        // u N′/N + v √−d Σ 2 Q′/(1 + Q²).
        Pair::Complex {
            u, v, norm, atans, ..
        } => {
            let mut sum = field.zero();
            for qj in atans {
                let (value, derivative) = values(field, qj, x);
                let denominator = field.add(&field.one(), &field.mul(&value, &value));
                let term = field.mul(&derivative, &field.inv(&denominator)?);
                sum = field.add(&sum, &field.add(&term, &term));
            }
            let scale = match root {
                Some(root) => surd(v * root),
                None => (Rational::ZERO, v.clone()),
            };
            let log = log_derivative(field, &lifted(norm), x)?;
            Some(field.add(&field.mul(&surd(u.clone()), &log), &field.mul(&scale, &sum)))
        }
    }
}

/// `∫ a / d`, checked: its derivative plus the part left is `a / d` at
/// several points; the rational part in lowest terms over a monic
/// denominator, the arguments monic and square-free with their real roots
/// counted, the norms of the pairs positive where the residues are complex,
/// the part left proper.
fn check(a: &Dense<Rational>, d: &Dense<Rational>) -> Result<RationalIntegral, TestCaseError> {
    let ring = qx();
    let r = ring.integrate_rational(a, d);
    let mut points = 0;
    for x in [q(7, 3), q(-11, 5), q(13, 2), q(2, 7), q(-5, 3)] {
        let dv = ring.eval(d, &x);
        let Some(found) = derivative_at(&r, &x) else {
            continue;
        };
        if dv.is_zero() {
            continue;
        }
        let expected = surd(&ring.eval(a, &x) / &dv);
        prop_assert_eq!(&found, &expected, "{:?} / {:?} at {:?}: {:?}", a, d, x, r);
        points += 1;
    }
    prop_assert!(points >= 2, "{:?}", r);
    let (gn, gd) = &r.rational;
    prop_assert!(ring.is_one(&ring.gcd(gn, gd)) || gn.is_zero());
    prop_assert_eq!(gd.lead(), Some(&Rational::ONE));
    for log in &r.logarithms {
        prop_assert_eq!(log.argument.lead(), Some(&Rational::ONE));
        let s = &log.argument;
        prop_assert!(ring.is_one(&ring.gcd(s, &ring.derivative(s))));
        prop_assert_eq!(log.real_root, real_roots(&ring, s) > 0);
    }
    for pair in &r.pairs {
        match pair {
            Pair::Real {
                norm, real_root, ..
            } => prop_assert_eq!(*real_root, real_roots(&ring, norm) > 0),
            Pair::Complex { norm, .. } => prop_assert_eq!(real_roots(&ring, norm), 0),
        }
    }
    if let Some((p, e)) = &r.rest {
        prop_assert!(p.degree() < e.degree());
    }
    Ok(r)
}

#[test]
fn examples() {
    let ring = qx();
    // ∫ 1/x² = −1/x.
    let r = check(&poly(&[1]), &poly(&[0, 0, 1])).unwrap();
    assert_eq!(r.rational, (poly(&[-1]), poly(&[0, 1])));
    assert!(r.logarithms.is_empty() && r.pairs.is_empty() && r.rest.is_none());
    // ∫ 1/(x² (x + 1)) = −1/x − ln|x| + ln|x + 1|.
    let r = check(&poly(&[1]), &poly(&[0, 0, 1, 1])).unwrap();
    assert_eq!(r.rational, (poly(&[-1]), poly(&[0, 1])));
    let logs: Vec<(Rational, Dense<Rational>, bool)> = r
        .logarithms
        .iter()
        .map(|l| (l.coefficient.clone(), l.argument.clone(), l.real_root))
        .collect();
    assert!(logs.contains(&(q(-1, 1), poly(&[0, 1]), true)));
    assert!(logs.contains(&(q(1, 1), poly(&[1, 1]), true)));
    // ∫ 1/(x² − 1) = ln|x − 1|/2 − ln|x + 1|/2.
    let r = check(&poly(&[1]), &poly(&[-1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 2);
    assert!(r.rest.is_none());
    // ∫ 2x/(x² + 1) = ln(x² + 1): no real root, no absolute value.
    let r = check(&poly(&[0, 2]), &poly(&[1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 1);
    assert_eq!(r.logarithms[0].argument, poly(&[1, 0, 1]));
    assert!(!r.logarithms[0].real_root);
    // ∫ 1/(x² + 1): residues ±i/2, √1 rational: one arctangent.
    let r = check(&poly(&[1]), &poly(&[1, 0, 1])).unwrap();
    assert!(r.logarithms.is_empty() && r.rest.is_none());
    assert!(matches!(&r.pairs[..], [Pair::Complex { atans, .. }] if atans.len() == 1));
    // ∫ 1/(x² + x + 1): residues ∓i/√3, an arctangent over ℚ(√(4/3)).
    let r = check(&poly(&[1]), &poly(&[1, 1, 1])).unwrap();
    assert!(matches!(&r.pairs[..], [Pair::Complex { d, .. }] if *d == q(-4, 3)));
    // ∫ 1/(x² − x − 1): residues ±1/√5, real: two logarithms joined.
    let r = check(&poly(&[1]), &poly(&[-1, -1, 1])).unwrap();
    assert!(matches!(&r.pairs[..], [Pair::Real { d, .. }] if *d == q(4, 5)));
    // ∫ 1/(x² − 2): ±√2/4, (√2/4) ln|(x − √2)/(x + √2)|.
    let r = check(&poly(&[1]), &poly(&[-2, 0, 1])).unwrap();
    assert!(matches!(
        &r.pairs[..],
        [Pair::Real { d, u, s, real_root: true, .. }]
            if *d == q(1, 2) && u.is_zero() && *s == [(q(0, 1), q(-2, 1)), (q(1, 1), q(0, 1))]
    ));
    // ∫ x/(x⁴ + 1) = atan(x²)/2: residues ∓i/4, S = x² + i of degree 2.
    let r = check(&poly(&[0, 1]), &poly(&[1, 0, 0, 0, 1])).unwrap();
    assert!(matches!(
        &r.pairs[..],
        [Pair::Complex { atans, .. }] if *atans == [vec![surd(q(0, 1)), surd(q(0, 1)), surd(q(1, 1))]]
    ));
    // Bronstein's example [Bro05 §2.8]: (x⁴ − 3x² + 6)/(x⁶ − 5x⁴ + 5x² + 4),
    // R = c (t² + 1)³, S of degree 3, three arctangents.
    let r = check(&poly(&[6, 0, -3, 0, 1]), &poly(&[4, 0, 5, 0, -5, 0, 1])).unwrap();
    assert!(matches!(&r.pairs[..], [Pair::Complex { atans, .. }] if atans.len() == 3));
    assert!(r.rest.is_none());
    // ∫ (x³ + 1)/(x − 1) = x³/3 + x²/2 + x + 2 ln|x − 1|.
    let r = check(&poly(&[1, 0, 0, 1]), &poly(&[-1, 1])).unwrap();
    assert_eq!(
        r.polynomial,
        ring.poly(vec![q(0, 1), q(1, 1), q(1, 2), q(1, 3)])
    );
    // ∫ 1/(x (x² + 1)): ln|x| and −ln(x² + 1)/2, residues 1 and −1/2 twice.
    let r = check(&poly(&[1]), &poly(&[0, 1, 0, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 2);
    assert!(r.pairs.is_empty() && r.rest.is_none());
    // ∫ 1/((x − 1)(x² + 1)): ln|x − 1|/2, and over x² + 1 the residues
    // −1/4 ∓ i/4: a pair.
    let r = check(&poly(&[1]), &poly(&[-1, 1, -1, 1])).unwrap();
    assert_eq!(r.logarithms.len(), 1);
    assert_eq!(r.logarithms[0].coefficient, q(1, 2));
    assert_eq!(r.pairs.len(), 1);
    assert_eq!(r.pairs[0].norm(), &poly(&[1, 0, 1]));
    assert!(r.rest.is_none());
    // ∫ 1/(x² + 1)² = x/(2(x² + 1)) + ∫ 1/(2(x² + 1)).
    let r = check(&poly(&[1]), &poly(&[1, 0, 2, 0, 1])).unwrap();
    assert_eq!(
        r.rational,
        (ring.poly(vec![q(0, 1), q(1, 2)]), poly(&[1, 0, 1]))
    );
    assert_eq!(r.pairs.len(), 1);
    // (x² + 1)/(x⁴ + x² + 1): giac's atan((x² − 1)/(x √3)) with a jump;
    // Rioboo's arctangents of polynomials have none.
    let r = check(&poly(&[1, 0, 1]), &poly(&[1, 0, 1, 0, 1])).unwrap();
    assert!(r.rest.is_none());
    // ∫ 1/(x³ − 2): the residues of x² + 2^(1/3) x + 2^(2/3) are cubic
    // irrationals, the part stays; the cube root's residue too.
    let r = check(&poly(&[1]), &poly(&[-2, 0, 0, 1])).unwrap();
    assert!(r.rest.is_some());
}

#[test]
fn sturm_counts_real_roots() {
    let ring = qx();
    for (p, n) in [
        (poly(&[1, 0, 1]), 0),
        (poly(&[-1, 0, 1]), 2),
        (poly(&[1, 0, -1, 1]), 1),
        (poly(&[0, -1, 0, 1]), 3),
        (poly(&[1, 0, 0, 0, 1]), 0),
        (poly(&[-2, 0, 0, 0, 1]), 2),
        (poly(&[5]), 0),
        (poly(&[3, 1]), 1),
    ] {
        assert_eq!(real_roots(&ring, &p), n, "{p:?}");
    }
}

/// A product of `(x − r)^e` times at most one quadratic `x² + b x + c` with
/// a discriminant that is not a square, times a constant.
type Denominator = (Vec<(i64, u32)>, Option<(i64, i64)>, i64);

fn denominator() -> impl Strategy<Value = Denominator> {
    (
        prop::collection::vec((-6i64..7, 1u32..4), 0..4),
        prop::option::of((-5i64..6, -6i64..7)),
        prop_oneof![Just(1i64), Just(-2), Just(3)],
    )
}

fn build((roots, quadratic, k): &Denominator) -> Option<Dense<Rational>> {
    let ring = qx();
    let mut d = poly(&[*k]);
    for &(r, e) in roots {
        for _ in 0..e {
            d = ring.mul(&d, &poly(&[-r, 1]));
        }
    }
    if let Some((b, c)) = quadratic {
        // Irreducible over ℚ: a discriminant that is not a square.
        if rational_sqrt(&q(b * b - 4 * c, 1)).is_some() {
            return None;
        }
        d = ring.mul(&d, &poly(&[*c, *b, 1]));
    }
    Some(d)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    /// Any numerator over products of linear factors with multiplicities
    /// and at most one irreducible quadratic: the derivative checks, and
    /// nothing is left, the quadratic's residues rational or a pair.
    #[test]
    fn derivatives_add_up(
        a in prop::collection::vec(-9i64..10, 0..8),
        shape in denominator(),
    ) {
        let Some(d) = build(&shape) else {
            return Ok(());
        };
        let r = check(&poly(&a), &d)?;
        prop_assert!(r.rest.is_none(), "{:?}", r.rest);
    }
}
