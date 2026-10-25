//! The roots of polynomials (CAS-211): `polynomial_roots` on random
//! polynomials over ℚ and ℚ(i) of degree 1 to 12, some with repeated
//! factors. Each ball it gives holds a root: the polynomial's value on the
//! ball holds 0. With rational coefficients, roots come only for a
//! square-free polynomial, all of them; the real ones are as many as
//! Sturm's count of real roots, and each real ball's interval holds exactly
//! one root by Sturm's count.
#![no_main]

use arbitrary::Unstructured;
use cas_num::algebra::{EuclideanDomain, Rationals};
use cas_num::ball::{Ball, ComplexBall, polynomial_roots};
use cas_num::{Complex, Integer, Rational};
use cas_num_fuzz::integer;
use cas_poly::{Dense, DenseRing};
use libfuzzer_sys::fuzz_target;

/// The highest degree, each factor's at most 4.
const MAX_DEGREE: usize = 12;

/// The precision the roots are sought at, and the most they may take.
const BITS: u64 = 53;
const MAX_BITS: u64 = 1024;

/// A coefficient: mostly a small integer, sometimes a small fraction or a
/// two-limb integer.
fn coefficient(u: &mut Unstructured<'_>) -> arbitrary::Result<Rational> {
    Ok(match u.int_in_range(0u8..=7)? {
        0 => Rational::new(
            Integer::from(u.int_in_range(-20i64..=20)?),
            Integer::from(u.int_in_range(1i64..=12)?),
        ),
        1 => Rational::from(integer(u, 2)?.0),
        _ => Rational::from(Integer::from(u.int_in_range(-9i64..=9)?)),
    })
}

/// A product of up to 4 random factors of degree 1 to 4, each to a power
/// of 1 or 2, of degree [`MAX_DEGREE`] at most; coefficients with
/// imaginary parts when `complex`.
fn polynomial(u: &mut Unstructured<'_>, complex: bool) -> arbitrary::Result<Vec<Complex>> {
    let mut p = vec![Complex::from(1)];
    for _ in 0..u.int_in_range(1..=4)? {
        let degree = u.int_in_range(1usize..=4)?;
        let power = u.int_in_range(1usize..=2)?;
        if p.len() - 1 + degree * power > MAX_DEGREE {
            break;
        }
        let mut factor = Vec::with_capacity(degree + 1);
        for _ in 0..=degree {
            let im = if complex {
                coefficient(u)?
            } else {
                Rational::ZERO
            };
            factor.push(Complex::new(coefficient(u)?, im));
        }
        if factor[degree].is_zero() {
            factor[degree] = Complex::from(1);
        }
        for _ in 0..power {
            p = times(&p, &factor);
        }
    }
    Ok(p)
}

fn times(p: &[Complex], q: &[Complex]) -> Vec<Complex> {
    let mut r = vec![Complex::from(0); p.len() + q.len() - 1];
    for (i, a) in p.iter().enumerate() {
        for (j, b) in q.iter().enumerate() {
            r[i + j] = &r[i + j] + &(a * b);
        }
    }
    r
}

/// `p(z)` in ball arithmetic.
fn value(p: &[Complex], z: &ComplexBall) -> ComplexBall {
    let prec = 256;
    p.iter().rev().fold(ComplexBall::ZERO, |v, c| {
        v.mul(z, prec)
            .add(&ComplexBall::from_complex(c, prec), prec)
    })
}

/// The ends of a real ball, as rationals.
fn ends(b: &Ball) -> (Rational, Rational) {
    let mid = b.mid().to_rational();
    let rad = b.rad().to_float().expect("a finite radius").to_rational();
    (&mid - &rad, &mid + &rad)
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let Ok(complex) = u.arbitrary::<bool>() else {
        return;
    };
    let Ok(p) = polynomial(&mut u, complex) else {
        return;
    };
    if p.len() < 2 {
        return;
    }
    let roots = polynomial_roots(&p, BITS, MAX_BITS);
    if let Some(roots) = &roots {
        assert_eq!(roots.len(), p.len() - 1, "{p:?}");
        for z in roots {
            assert!(value(&p, z).may_be_zero(), "{p:?}: {z:?} holds no root");
        }
    }
    if complex {
        return;
    }
    let qx = DenseRing::new(Rationals);
    let q: Dense<Rational> = qx.poly(p.iter().map(|c| c.re().clone()).collect());
    let squarefree = EuclideanDomain::gcd(&qx, &q, &qx.derivative(&q))
        .degree()
        .is_some_and(|d| d == 0);
    let Some(roots) = roots else {
        return;
    };
    assert!(
        squarefree,
        "{p:?}: roots of a polynomial with a repeated factor"
    );
    let real: Vec<&ComplexBall> = roots.iter().filter(|z| z.is_real()).collect();
    assert_eq!(
        real.len(),
        qx.count_real_roots(&q, None, None),
        "{p:?}: {roots:?}"
    );
    for z in real {
        let (lo, hi) = ends(z.re());
        assert_eq!(
            qx.count_real_roots(&q, Some(&lo), Some(&hi)),
            1,
            "{p:?}: {z:?}"
        );
    }
});
