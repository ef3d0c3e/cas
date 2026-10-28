//! Tests of Gröbner bases (CAS-DOC-21 §2, §3): the orders, division, known
//! bases, and the defining properties of the results.

use cas_num::{Integer, Rational};

use super::{
    Distributed, GroebnerError, GroebnerLimits, MonomialOrder, groebner_basis, groebner_reduce,
    s_polynomial,
};
use crate::cancel::Never;

use MonomialOrder::{DegLex, DegRevLex, Lex};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// A polynomial from exponent vectors and integer coefficients.
fn poly(terms: &[(&[u64], i64)], order: MonomialOrder) -> Distributed {
    Distributed::new(
        terms.iter().map(|(m, c)| (m.to_vec(), q(*c, 1))).collect(),
        order,
    )
}

/// The reduced basis of `polys`, checked: each input reduces to 0 by it,
/// so does each S-polynomial of two of its elements, each element is monic
/// and no term of one is divisible by another's leading monomial.
fn basis(polys: &[Distributed], order: MonomialOrder) -> Vec<Distributed> {
    let g = groebner_basis(polys, order, GroebnerLimits::NONE, &Never).unwrap();
    for p in polys {
        assert!(
            groebner_reduce(p, &g, order, GroebnerLimits::NONE, &Never)
                .unwrap()
                .is_zero(),
            "{p:?} by {g:?}"
        );
    }
    for (i, a) in g.iter().enumerate() {
        assert!(a.terms()[0].1.is_one(), "{a:?} monic");
        for (j, b) in g.iter().enumerate() {
            if i < j {
                let s = s_polynomial(a, b, order);
                assert!(
                    groebner_reduce(&s, &g, order, GroebnerLimits::NONE, &Never)
                        .unwrap()
                        .is_zero(),
                    "S({a:?}, {b:?})"
                );
            }
            if i != j {
                let lead = b.leading_monomial().unwrap();
                assert!(
                    a.terms()
                        .iter()
                        .all(|(m, _)| !m.iter().zip(lead).all(|(x, y)| y <= x)),
                    "{a:?} reduced by {b:?}"
                );
            }
        }
    }
    g
}

/// The basis, each element primitive over ℤ, as `poly` writes them.
fn primitive(g: &[Distributed]) -> Vec<Distributed> {
    g.iter().map(Distributed::primitive).collect()
}

/// Cyclic-3: x + y + z, xy + yz + zx, xyz − 1 under `Lex`.
fn cyclic() -> [Distributed; 3] {
    [
        poly(&[(&[1, 0, 0], 1), (&[0, 1, 0], 1), (&[0, 0, 1], 1)], Lex),
        poly(&[(&[1, 1, 0], 1), (&[0, 1, 1], 1), (&[1, 0, 1], 1)], Lex),
        poly(&[(&[1, 1, 1], 1), (&[0, 0, 0], -1)], Lex),
    ]
}

#[test]
fn orders() {
    // x y² against y³ z⁴: lexicographically greater, of a lower degree.
    let (a, b) = ([1, 2, 0], [0, 3, 4]);
    assert!(Lex.cmp(&a, &b).is_gt());
    assert!(DegLex.cmp(&a, &b).is_lt());
    assert!(DegRevLex.cmp(&a, &b).is_lt());
    // x⁴ y⁷ z against x⁴ y² z⁴, both of degree 12 [CLO, §2.2].
    let (a, b) = ([4, 7, 1], [4, 2, 4]);
    assert!(DegLex.cmp(&a, &b).is_gt());
    assert!(DegRevLex.cmp(&a, &b).is_gt());
    // x y⁵ z² against x⁴ y z³, both of degree 8: the orders disagree.
    let (a, b) = ([1, 5, 2], [4, 1, 3]);
    assert!(DegLex.cmp(&a, &b).is_lt());
    assert!(DegRevLex.cmp(&a, &b).is_gt());
    for order in [Lex, DegLex, DegRevLex] {
        assert!(order.cmp(&a, &a).is_eq());
    }
}

/// Division by lists, the divisors' order mattering [CLO, §2.3].
#[test]
fn division() {
    let f = poly(&[(&[2, 1], 1), (&[1, 2], 1), (&[0, 2], 1)], Lex);
    let f1 = poly(&[(&[1, 1], 1), (&[0, 0], -1)], Lex);
    let f2 = poly(&[(&[0, 2], 1), (&[0, 0], -1)], Lex);
    let r = groebner_reduce(
        &f,
        &[f1.clone(), f2.clone()],
        Lex,
        GroebnerLimits::NONE,
        &Never,
    )
    .unwrap();
    assert_eq!(r, poly(&[(&[1, 0], 1), (&[0, 1], 1), (&[0, 0], 1)], Lex));
    let r = groebner_reduce(&f, &[f2, f1], Lex, GroebnerLimits::NONE, &Never).unwrap();
    assert_eq!(r, poly(&[(&[1, 0], 2), (&[0, 0], 1)], Lex));
    // Exact remainders: x³ + y² by [2y² − 1, x − y] is y/2 + 1/2, giac's.
    let order = DegRevLex;
    let p = poly(&[(&[3, 0], 1), (&[0, 2], 1)], order);
    let divisors = [
        poly(&[(&[0, 2], 2), (&[0, 0], -1)], order),
        poly(&[(&[1, 0], 1), (&[0, 1], -1)], order),
    ];
    assert_eq!(
        groebner_reduce(&p, &divisors, order, GroebnerLimits::NONE, &Never).unwrap(),
        Distributed::new(vec![(vec![0, 1], q(1, 2)), (vec![0, 0], q(1, 2))], order)
    );
}

#[test]
fn bases() {
    // x² + y² − 1, x − y: giac's [2*y^2-1, x-y] under each order, in
    // increasing order of the leading monomials, x below y² but under plex.
    for order in [Lex, DegLex, DegRevLex] {
        let polys = [
            poly(&[(&[2, 0], 1), (&[0, 2], 1), (&[0, 0], -1)], order),
            poly(&[(&[1, 0], 1), (&[0, 1], -1)], order),
        ];
        let mut expected = vec![
            poly(&[(&[0, 2], 2), (&[0, 0], -1)], order),
            poly(&[(&[1, 0], 1), (&[0, 1], -1)], order),
        ];
        if order != Lex {
            expected.reverse();
        }
        assert_eq!(primitive(&basis(&polys, order)), expected, "{order:?}");
    }
    // x³ − 2xy, x²y − 2y² + x [CLO, §2.7]: x², xy, 2y² − x under
    // revlex; y³, x − 2y² under plex, as giac's.
    let generators = |order| {
        [
            poly(&[(&[3, 0], 1), (&[1, 1], -2)], order),
            poly(&[(&[2, 1], 1), (&[0, 2], -2), (&[1, 0], 1)], order),
        ]
    };
    assert_eq!(
        primitive(&basis(&generators(DegRevLex), DegRevLex)),
        [
            poly(&[(&[0, 2], 2), (&[1, 0], -1)], DegRevLex),
            poly(&[(&[1, 1], 1)], DegRevLex),
            poly(&[(&[2, 0], 1)], DegRevLex),
        ]
    );
    assert_eq!(
        primitive(&basis(&generators(Lex), Lex)),
        [
            poly(&[(&[0, 3], 1)], Lex),
            poly(&[(&[1, 0], 1), (&[0, 2], -2)], Lex),
        ]
    );
    // The twisted cubic, x² − y and x³ − z under plex.
    let twisted = [
        poly(&[(&[2, 0, 0], 1), (&[0, 1, 0], -1)], Lex),
        poly(&[(&[3, 0, 0], 1), (&[0, 0, 1], -1)], Lex),
    ];
    assert_eq!(
        primitive(&basis(&twisted, Lex)),
        [
            poly(&[(&[0, 3, 0], 1), (&[0, 0, 2], -1)], Lex),
            poly(&[(&[1, 0, 1], 1), (&[0, 2, 0], -1)], Lex),
            poly(&[(&[1, 1, 0], 1), (&[0, 0, 1], -1)], Lex),
            poly(&[(&[2, 0, 0], 1), (&[0, 1, 0], -1)], Lex),
        ]
    );
    // Cyclic-3: z³ − 1, y² + yz + z², x + y + z, giac's.
    assert_eq!(
        primitive(&basis(&cyclic(), Lex)),
        [
            poly(&[(&[0, 0, 3], 1), (&[0, 0, 0], -1)], Lex),
            poly(&[(&[0, 2, 0], 1), (&[0, 1, 1], 1), (&[0, 0, 2], 1)], Lex),
            poly(&[(&[1, 0, 0], 1), (&[0, 1, 0], 1), (&[0, 0, 1], 1)], Lex),
        ]
    );
    // x²y − 1, xy² − x: x² − y, y² − 1, giac's.
    let polys = [
        poly(&[(&[2, 1], 1), (&[0, 0], -1)], DegRevLex),
        poly(&[(&[1, 2], 1), (&[1, 0], -1)], DegRevLex),
    ];
    assert_eq!(
        primitive(&basis(&polys, DegRevLex)),
        [
            poly(&[(&[0, 2], 1), (&[0, 0], -1)], DegRevLex),
            poly(&[(&[2, 0], 1), (&[0, 1], -1)], DegRevLex),
        ]
    );
}

/// The whole ring, the zero ideal, a rational coefficient.
#[test]
fn edges() {
    let one = poly(&[(&[0], 1)], Lex);
    let polys = [
        poly(&[(&[2], 1), (&[0], -1)], Lex),
        poly(&[(&[1], 1), (&[0], -2)], Lex),
    ];
    assert_eq!(basis(&polys, Lex), std::slice::from_ref(&one));
    assert_eq!(basis(std::slice::from_ref(&one), Lex), [one]);
    assert_eq!(basis(&[Distributed::new(Vec::new(), Lex)], Lex), []);
    assert_eq!(
        groebner_basis(&[], Lex, GroebnerLimits::NONE, &Never).unwrap(),
        []
    );
    // x − 1/2 and y: giac's [2*x-1, y].
    let polys = [
        Distributed::new(
            vec![(vec![1, 0], q(1, 1)), (vec![0, 0], q(-1, 2))],
            DegRevLex,
        ),
        poly(&[(&[0, 1], 1)], DegRevLex),
    ];
    assert_eq!(
        primitive(&basis(&polys, DegRevLex)),
        [
            poly(&[(&[0, 1], 1)], DegRevLex),
            poly(&[(&[1, 0], 2), (&[0, 0], -1)], DegRevLex),
        ]
    );
    // An exponent past the bound.
    let huge = poly(&[(&[super::MAX_EXPONENT + 1], 1)], Lex);
    assert_eq!(
        groebner_basis(&[huge], Lex, GroebnerLimits::NONE, &Never),
        Err(GroebnerError::TooLarge)
    );
}

/// Past the limits a caller sets, the work or the bits of a coefficient,
/// a computation stops (CAS-291): the fuzz target's step budget let one
/// input run for minutes.
#[test]
fn limits() {
    let little = GroebnerLimits {
        work: 8,
        coefficient_bits: u64::MAX,
    };
    assert_eq!(
        groebner_basis(&cyclic(), Lex, little, &Never),
        Err(GroebnerError::TooLarge)
    );
    // x²y + xy² + y² by xy − 1: the first step makes two coefficients.
    let f = poly(&[(&[2, 1], 1), (&[1, 2], 1), (&[0, 2], 1)], Lex);
    let divisor = poly(&[(&[1, 1], 1), (&[0, 0], -1)], Lex);
    let once = GroebnerLimits { work: 1, ..little };
    assert_eq!(
        groebner_reduce(&f, &[divisor], Lex, once, &Never),
        Err(GroebnerError::TooLarge)
    );
    // 3x − 1, 5y − x: y − 1/15 has a denominator of 4 bits.
    let polys = [
        poly(&[(&[1, 0], 3), (&[0, 0], -1)], Lex),
        poly(&[(&[0, 1], 5), (&[1, 0], -1)], Lex),
    ];
    let bits = |coefficient_bits| GroebnerLimits {
        work: u64::MAX,
        coefficient_bits,
    };
    assert_eq!(
        groebner_basis(&polys, Lex, bits(3), &Never),
        Err(GroebnerError::TooLarge)
    );
    assert_eq!(
        primitive(&groebner_basis(&polys, Lex, bits(4), &Never).unwrap()),
        [
            poly(&[(&[0, 1], 15), (&[0, 0], -1)], Lex),
            poly(&[(&[1, 0], 3), (&[0, 0], -1)], Lex),
        ]
    );
}
