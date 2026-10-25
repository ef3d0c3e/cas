//! The integration of rational functions over ℚ (CAS-196, CAS-197): random
//! numerators over products of powers of random factors of degree 1 to 3.
//! The derivative of the result plus the part left is the integrand at
//! rational points, exactly: each pair's part in the field of its residues,
//! ℚ(√d) or ℚ(√−d), where it must come out rational. The rational part is
//! in lowest terms over a monic denominator; the logarithms' arguments and
//! the norms of the pairs are monic, the arguments square-free, a real
//! pair's norm `S S̄`; the part left is proper over a square-free
//! denominator, and there is none when every factor has degree 2 at most.
#![no_main]

use arbitrary::Unstructured;
use cas_num::algebra::{Rationals, Ring};
use cas_num::{Integer, Rational};
use cas_num_fuzz::integer;
use cas_poly::{Dense, DenseRing, Pair, RationalIntegral, Surd};
use libfuzzer_sys::fuzz_target;

/// At most this many factors in a denominator, each of degree 3 at most and
/// to a power 3 at most: the square-free part has degree 9 at most.
const MAX_FACTORS: usize = 3;

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

/// A polynomial of degree `degree` at most.
fn poly(
    u: &mut Unstructured<'_>,
    qx: &DenseRing<Rationals>,
    degree: usize,
) -> arbitrary::Result<Dense<Rational>> {
    let coeffs = (0..=degree)
        .map(|_| coefficient(u))
        .collect::<arbitrary::Result<_>>()?;
    Ok(qx.poly(coeffs))
}

/// `a + b √k`.
fn mul(k: &Rational, x: &Surd, y: &Surd) -> Surd {
    (
        &(&x.0 * &y.0) + &(k * &(&x.1 * &y.1)),
        &(&x.0 * &y.1) + &(&x.1 * &y.0),
    )
}

fn add(x: &Surd, y: &Surd) -> Surd {
    (&x.0 + &y.0, &x.1 + &y.1)
}

fn sub(x: &Surd, y: &Surd) -> Surd {
    (&x.0 - &y.0, &x.1 - &y.1)
}

/// `1 / x` in ℚ(√k), `k` not a square: the conjugate over the norm.
fn inv(k: &Rational, x: &Surd) -> Surd {
    let norm = &(&x.0 * &x.0) - &(k * &(&x.1 * &x.1));
    let r = norm.checked_recip().expect("a nonzero element");
    (&x.0 * &r, -&(&x.1 * &r))
}

fn rational(r: &Rational) -> Surd {
    (r.clone(), Rational::ZERO)
}

/// `p(x)` and `p′(x)` for `p` over ℚ(√k).
fn values(k: &Rational, p: &[Surd], x: &Rational) -> (Surd, Surd) {
    let at = rational(x);
    let zero = rational(&Rational::ZERO);
    let (mut value, mut derivative) = (zero.clone(), zero);
    for c in p.iter().rev() {
        derivative = add(&mul(k, &derivative, &at), &value);
        value = add(&mul(k, &value, &at), c);
    }
    (value, derivative)
}

fn lifted(p: &Dense<Rational>) -> Vec<Surd> {
    p.coeffs().iter().map(rational).collect()
}

/// `p′(x) / p(x)` over ℚ(√k), `p(x)` not 0.
fn log_derivative(k: &Rational, p: &[Surd], x: &Rational) -> Surd {
    let (value, derivative) = values(k, p, x);
    mul(k, &derivative, &inv(k, &value))
}

/// The rational square root of `r`, if it has one.
fn rational_sqrt(r: &Rational) -> Option<Rational> {
    let (n, d) = (r.numer(), r.denom());
    (!r.is_negative() && n.is_perfect_square() && d.is_perfect_square())
        .then(|| Rational::new(n.isqrt(), d.isqrt()))
}

/// The derivative of a pair's part at `x`, which is rational; checks its
/// form on the way.
fn pair_derivative(qx: &DenseRing<Rationals>, pair: &Pair, x: &Rational) -> Rational {
    let two = Rational::from(Integer::from(2));
    let (u, v, norm) = match pair {
        Pair::Real { u, v, norm, .. } | Pair::Complex { u, v, norm, .. } => (u, v, norm),
    };
    assert_eq!(*v, Rational::new(Integer::ONE, Integer::from(2)), "v = 1/2");
    assert_eq!(norm.lead(), Some(&Rational::ONE), "a monic norm");
    let log = &qx.eval(&qx.derivative(norm), x) / &qx.eval(norm, x);
    let part = match pair {
        // u N′/N + v √d (S′/S − S̄′/S̄), the difference a multiple of √d.
        Pair::Real { d, s, norm, .. } => {
            assert!(d.is_positive() && rational_sqrt(d).is_none(), "{d:?}");
            assert_eq!(s.last(), Some(&rational(&Rational::ONE)), "a monic S");
            let conjugate: Vec<Surd> = s.iter().map(|(a, b)| (a.clone(), -b)).collect();
            let product = s.iter().enumerate().fold(
                vec![rational(&Rational::ZERO); 2 * s.len() - 1],
                |mut acc, (i, a)| {
                    for (j, b) in conjugate.iter().enumerate() {
                        acc[i + j] = add(&acc[i + j], &mul(d, a, b));
                    }
                    acc
                },
            );
            assert_eq!(product, lifted(norm), "N = S S̄");
            let difference = sub(&log_derivative(d, s, x), &log_derivative(d, &conjugate, x));
            assert!(difference.0.is_zero(), "{difference:?}");
            &(v * d) * &difference.1
        }
        // u N′/N + v √−d Σ 2 Q′/(1 + Q²): over ℚ(√−d) a multiple of √−d,
        // over ℚ when √−d is rational.
        Pair::Complex { d, atans, .. } => {
            assert!(d.is_negative(), "{d:?}");
            let k = -d;
            let root = rational_sqrt(&k);
            let mut sum = rational(&Rational::ZERO);
            for q in atans {
                let (value, derivative) = values(&k, q, x);
                if root.is_some() {
                    assert!(q.iter().all(|c| c.1.is_zero()), "{q:?}");
                }
                let denominator = add(&rational(&Rational::ONE), &mul(&k, &value, &value));
                let term = mul(&k, &derivative, &inv(&k, &denominator));
                sum = add(&sum, &mul(&k, &rational(&two), &term));
            }
            match root {
                Some(root) => {
                    assert!(sum.1.is_zero(), "{sum:?}");
                    &(v * &root) * &sum.0
                }
                None => {
                    assert!(sum.0.is_zero(), "{sum:?}");
                    &(v * &k) * &sum.1
                }
            }
        }
    };
    &(u * &log) + &part
}

/// The derivative of the result at `x`, plus the part left at `x`.
fn derivative_at(qx: &DenseRing<Rationals>, r: &RationalIntegral, x: &Rational) -> Rational {
    let mut sum = qx.eval(&qx.derivative(&r.polynomial), x);
    let (n, d) = &r.rational;
    let (nv, dv) = (qx.eval(n, x), qx.eval(d, x));
    let (nd, dd) = (qx.eval(&qx.derivative(n), x), qx.eval(&qx.derivative(d), x));
    sum = &sum + &(&(&(&nd * &dv) - &(&nv * &dd)) / &(&dv * &dv));
    for log in &r.logarithms {
        let s = &log.argument;
        sum = &sum + &(&log.coefficient * &(&qx.eval(&qx.derivative(s), x) / &qx.eval(s, x)));
    }
    for pair in &r.pairs {
        sum = &sum + &pair_derivative(qx, pair, x);
    }
    if let Some((p, e)) = &r.rest {
        sum = &sum + &(&qx.eval(p, x) / &qx.eval(e, x));
    }
    sum
}

fn run(u: &mut Unstructured<'_>) -> arbitrary::Result<()> {
    let qx = DenseRing::new(Rationals);
    let mut d = qx.poly(vec![coefficient(u)?]);
    let mut small = true;
    for _ in 0..u.int_in_range(0..=MAX_FACTORS)? {
        let degree = u.int_in_range(1..=3)?;
        let factor = poly(u, &qx, degree)?;
        small &= factor.degree().unwrap_or(0) <= 2;
        for _ in 0..u.int_in_range(1..=3)? {
            d = qx.mul(&d, &factor);
        }
    }
    if d.is_zero() {
        return Ok(());
    }
    let degree = u.int_in_range(0..=d.degree().unwrap_or(0) + 2)?;
    let a = poly(u, &qx, degree)?;
    let r = qx.integrate_rational(&a, &d);
    let (gn, gd) = &r.rational;
    assert!(gn.is_zero() || qx.is_one(&qx.gcd(gn, gd)), "lowest terms");
    assert_eq!(gd.lead(), Some(&Rational::ONE), "a monic denominator");
    for log in &r.logarithms {
        let s = &log.argument;
        assert!(!log.coefficient.is_zero() && s.degree().unwrap_or(0) > 0);
        assert_eq!(s.lead(), Some(&Rational::ONE), "a monic argument");
        assert!(qx.is_one(&qx.gcd(s, &qx.derivative(s))), "square-free");
    }
    if let Some((p, e)) = &r.rest {
        assert!(!small, "a part left over factors of degree 2 at most");
        assert!(p.degree() < e.degree() && qx.is_one(&qx.gcd(e, &qx.derivative(e))));
    }
    // Three points that are not poles, among 30 distinct ones: `d` has 27
    // roots at most.
    let mut points = 0;
    for k in 1i64..=30 {
        let x = Rational::new(
            Integer::from(if k % 2 == 0 { k } else { -k }),
            Integer::from(k % 7 + 2),
        );
        let dv = qx.eval(&d, &x);
        if dv.is_zero() {
            continue;
        }
        let expected = &qx.eval(&a, &x) / &dv;
        assert_eq!(
            derivative_at(&qx, &r, &x),
            expected,
            "{a:?} / {d:?} at {x:?}: {r:?}"
        );
        points += 1;
        if points == 3 {
            break;
        }
    }
    assert_eq!(points, 3, "{d:?}");
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
