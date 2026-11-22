//! Tests of normal sets, rational univariate representations and syzygies
//! (CAS-541): giac 2.0.0's answers, and the defining properties.

use cas_num::{Integer, Rational};

use super::{MAX_NORMAL_SET, Rur, div_rem, mul, normal_set, rational_univariate, syzygies};
use crate::cancel::Never;
use crate::groebner::{
    Distributed, GroebnerError, GroebnerLimits, MonomialOrder::DegRevLex, groebner_basis,
};

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

/// A polynomial from integer coefficients and exponent vectors, under
/// `DegRevLex`.
fn p(terms: &[(i64, &[u64])]) -> Distributed {
    Distributed::new(
        terms.iter().map(|(c, m)| (m.to_vec(), q(*c, 1))).collect(),
        DegRevLex,
    )
}

/// A polynomial in one variable from its integer coefficients, from
/// degree 0.
fn u(c: &[i64]) -> Vec<Rational> {
    c.iter().map(|&c| q(c, 1)).collect()
}

fn rur(polys: &[Distributed]) -> Option<Rur> {
    rational_univariate(polys, GroebnerLimits::NONE, &Never).unwrap()
}

fn expect(polys: &[Distributed], form: &[i64], f: &[Rational], f1: &[i64], g: &[&[i64]]) {
    let r = rur(polys).unwrap();
    assert_eq!(r.form, u(form), "{polys:?}");
    assert_eq!(r.minimal, f, "{polys:?}");
    assert_eq!(r.derivative, u(f1), "{polys:?}");
    let g: Vec<Vec<Rational>> = g.iter().map(|g| u(g)).collect();
    assert_eq!(r.numerators, g, "{polys:?}");
}

#[test]
fn normal_sets() {
    // x² + y² − 1, x − y: the basis x − y, y² − 1/2, the normal set 1, y.
    let g = groebner_basis(
        &[
            p(&[(1, &[2, 0]), (1, &[0, 2]), (-1, &[0, 0])]),
            p(&[(1, &[1, 0]), (-1, &[0, 1])]),
        ],
        DegRevLex,
        GroebnerLimits::NONE,
        &Never,
    )
    .unwrap();
    assert_eq!(
        normal_set(&g, DegRevLex),
        Some(vec![vec![0, 0], vec![0, 1]])
    );
    // x², xy, y²: 1, y, x in increasing order.
    let g = [p(&[(1, &[2, 0])]), p(&[(1, &[1, 1])]), p(&[(1, &[0, 2])])];
    assert_eq!(
        normal_set(&g, DegRevLex),
        Some(vec![vec![0, 0], vec![0, 1], vec![1, 0]])
    );
    assert_eq!(
        normal_set(&[p(&[(1, &[0, 0])])], DegRevLex),
        Some(Vec::new())
    );
    let circle = p(&[(1, &[2, 0]), (1, &[0, 2]), (-1, &[0, 0])]);
    assert_eq!(normal_set(&[circle], DegRevLex), None);
    assert_eq!(normal_set(&[], DegRevLex), None);
}

#[test]
fn giac_representations() {
    // x² + y² − 1, x − y.
    expect(
        &[
            p(&[(1, &[2, 0]), (1, &[0, 2]), (-1, &[0, 0])]),
            p(&[(1, &[1, 0]), (-1, &[0, 1])]),
        ],
        &[0, 1],
        &[q(1, 2), q(0, 1), q(-1, 1)],
        &[0, -2],
        &[&[-1], &[-1]],
    );
    // x² − 2, y² − 3.
    expect(
        &[
            p(&[(1, &[2, 0]), (-2, &[0, 0])]),
            p(&[(1, &[0, 2]), (-3, &[0, 0])]),
        ],
        &[1, 1],
        &u(&[-1, 0, 10, 0, -1]),
        &[0, 20, 0, -4],
        &[&[-8, 0, -8], &[12, 0, -12]],
    );
    // x², y.
    expect(
        &[p(&[(1, &[2, 0])]), p(&[(1, &[0, 1])])],
        &[1, 0],
        &u(&[0, -1]),
        &[-1],
        &[&[], &[]],
    );
    // x²(x − 1), y.
    expect(
        &[p(&[(1, &[3, 0]), (-1, &[2, 0])]), p(&[(1, &[0, 1])])],
        &[1, 0],
        &u(&[0, 1, -1]),
        &[1, -2],
        &[&[0, -1], &[]],
    );
    // x² − 1, (y − 1)².
    expect(
        &[
            p(&[(1, &[2, 0]), (-1, &[0, 0])]),
            p(&[(1, &[0, 2]), (-2, &[0, 1]), (1, &[0, 0])]),
        ],
        &[1, 0],
        &u(&[1, 0, -1]),
        &[0, -2],
        &[&[-2], &[0, -2]],
    );
    // x³ − y, y² − x.
    expect(
        &[
            p(&[(1, &[3, 0]), (-1, &[0, 1])]),
            p(&[(1, &[0, 2]), (-1, &[1, 0])]),
        ],
        &[0, 1],
        &u(&[0, 1, 0, 0, 0, 0, -1]),
        &[1, 0, 0, 0, 0, -6],
        &[&[0, 0, -5], &[0, -5]],
    );
    // x − 1, y − 2, z − 3.
    expect(
        &[
            p(&[(1, &[1, 0, 0]), (-1, &[0, 0, 0])]),
            p(&[(1, &[0, 1, 0]), (-2, &[0, 0, 0])]),
            p(&[(1, &[0, 0, 1]), (-3, &[0, 0, 0])]),
        ],
        &[0, 0, 1],
        &u(&[3, -1]),
        &[-1],
        &[&[-1], &[-2], &[-3]],
    );
    // x² − 3x + 2, y − x².
    expect(
        &[
            p(&[(1, &[2, 0]), (-3, &[1, 0]), (2, &[0, 0])]),
            p(&[(1, &[0, 1]), (-1, &[2, 0])]),
        ],
        &[0, 1],
        &u(&[-4, 5, -1]),
        &[5, -2],
        &[&[6, -3], &[8, -5]],
    );
}

#[test]
fn no_representation() {
    let circle = p(&[(1, &[2, 0]), (1, &[0, 2]), (-1, &[0, 0])]);
    assert_eq!(rur(&[circle]), None);
    // x² − y², x − y: the line x = y.
    let line = [
        p(&[(1, &[2, 0]), (-1, &[0, 2])]),
        p(&[(1, &[1, 0]), (-1, &[0, 1])]),
    ];
    assert_eq!(rur(&line), None);
    assert_eq!(rur(&[p(&[(1, &[0, 0])])]), None);
    assert_eq!(rur(&[p(&[(1, &[2]), (-2, &[0])])]), None);
    assert_eq!(rur(&[]), None);
}

#[test]
fn separating_combination() {
    // x² − 1, y² − 1: no variable separates the four points; giac picks
    // −6x − 4y, we a form of the sequence. Each generator vanishes at
    // xᵢ = gᵢ(t)/f'(t) modulo f: its numerator over f'^deg is 0 mod f.
    let polys = [
        p(&[(1, &[2, 0]), (-1, &[0, 0])]),
        p(&[(1, &[0, 2]), (-1, &[0, 0])]),
    ];
    let r = rur(&polys).unwrap();
    assert_eq!(r.minimal.len(), 5);
    for g in &polys {
        let top = g
            .terms()
            .iter()
            .map(|(m, _)| m.iter().sum::<u64>())
            .max()
            .unwrap();
        let mut sum: Vec<Rational> = Vec::new();
        for (m, c) in g.terms() {
            let mut t = vec![c.clone()];
            for (e, numerator) in m.iter().zip(&r.numerators) {
                for _ in 0..*e {
                    t = mul(&t, numerator);
                }
            }
            for _ in 0..top - m.iter().sum::<u64>() {
                t = mul(&t, &r.derivative);
            }
            let len = sum.len().max(t.len());
            sum.resize(len, Rational::ZERO);
            for (s, x) in sum.iter_mut().zip(&t) {
                *s += x;
            }
        }
        assert!(div_rem(&sum, &r.minimal).1.is_empty(), "{g:?}");
    }
}

#[test]
fn normal_set_bound() {
    // x^300 − 1, y − 1: 300 points.
    let polys = [
        p(&[(1, &[300, 0]), (-1, &[0, 0])]),
        p(&[(1, &[0, 1]), (-1, &[0, 0])]),
    ];
    const { assert!(300 > MAX_NORMAL_SET) };
    assert_eq!(
        rational_univariate(&polys, GroebnerLimits::NONE, &Never),
        Err(GroebnerError::TooLarge)
    );
}

/// The syzygies of `polys`, checked: each element linear in the `m` first
/// variables, `Σ aᵢAᵢ`, has `Σ aᵢfᵢ = 0`; the others are products of two
/// of them.
fn syzygy_basis(polys: &[Distributed]) -> Vec<Distributed> {
    let m = polys.len();
    let s = syzygies(polys, GroebnerLimits::NONE, &Never).unwrap();
    for g in &s {
        let degree = |e: &[u64]| e[..m].iter().sum::<u64>();
        let a = degree(g.leading_monomial().unwrap());
        assert!(g.terms().iter().all(|(e, _)| degree(e) == a), "{g:?}");
        if a != 1 {
            continue;
        }
        let mut terms = Vec::new();
        for (e, c) in g.terms() {
            let i = e[..m].iter().position(|&x| x == 1).unwrap();
            for (x, d) in polys[i].terms() {
                let monomial = e[m..].iter().zip(x).map(|(a, b)| a + b).collect();
                terms.push((monomial, c * d));
            }
        }
        assert!(Distributed::new(terms, DegRevLex).is_zero(), "{g:?}");
    }
    s
}

/// Whether `s` has, up to a constant factor, each of `expected`.
fn has(s: &[Distributed], expected: &[Distributed]) -> bool {
    let s: Vec<Distributed> = s.iter().map(Distributed::primitive).collect();
    expected.iter().all(|e| s.contains(&e.primitive()))
}

#[test]
fn giac_syzygies() {
    let squares = |n: usize| -> Vec<Distributed> {
        let mut out = Vec::new();
        for i in 0..2 {
            for j in i..2 {
                let mut e = vec![0; 2 + n];
                e[i] += 1;
                e[j] += 1;
                out.push(p(&[(1, &e)]));
            }
        }
        out
    };
    // x, y in [A0, A1, x, y].
    let s = syzygy_basis(&[p(&[(1, &[1, 0])]), p(&[(1, &[0, 1])])]);
    let mut e = squares(2);
    e.push(p(&[(1, &[1, 0, 0, 1]), (-1, &[0, 1, 1, 0])]));
    assert_eq!(s.len(), 4);
    assert!(has(&s, &e), "{s:?}");
    // xy, x² − y³.
    let s = syzygy_basis(&[p(&[(1, &[1, 1])]), p(&[(1, &[2, 0]), (-1, &[0, 3])])]);
    let mut e = squares(2);
    e.push(p(&[
        (1, &[1, 0, 2, 0]),
        (-1, &[1, 0, 0, 3]),
        (-1, &[0, 1, 1, 1]),
    ]));
    assert_eq!(s.len(), 4);
    assert!(has(&s, &e), "{s:?}");
    // x − 1, x² − 1 in [A0, A1, x].
    let s = syzygy_basis(&[p(&[(1, &[1]), (-1, &[0])]), p(&[(1, &[2]), (-1, &[0])])]);
    let mut e = squares(1);
    e.push(p(&[(1, &[1, 0, 1]), (1, &[1, 0, 0]), (-1, &[0, 1, 0])]));
    assert_eq!(s.len(), 4);
    assert!(has(&s, &e), "{s:?}");
    // x², xy, y² in [A0, A1, A2, x, y].
    let s = syzygy_basis(&[p(&[(1, &[2, 0])]), p(&[(1, &[1, 1])]), p(&[(1, &[0, 2])])]);
    let e = [
        p(&[(1, &[0, 1, 0, 1, 0]), (-1, &[1, 0, 0, 0, 1])]),
        p(&[(1, &[0, 0, 1, 1, 0]), (-1, &[0, 1, 0, 0, 1])]),
    ];
    assert!(has(&s, &e), "{s:?}");
}
