//! Tests of the sparse interpolation of gcd images (CAS-190), over a word
//! prime field: images from the forms of known gcds, with single and with
//! multiple scaling; wrong forms, a gcd of too large a degree and a content
//! reported.

use std::collections::BTreeSet;

use cas_num::Integer;
use cas_num::algebra::{Field, Ring, WordPrimeField};
use proptest::prelude::*;

use super::{Failure, Form, image};
use crate::cancel::Never;
use crate::gcd::modular::pgcd;
use crate::random::Random;
use crate::sparse::{Sparse, SparseRing};

/// A word prime: random points rarely collide in it.
const P: u64 = (1 << 61) - 1;

type Fp = SparseRing<WordPrimeField>;

fn ring(n: usize) -> Fp {
    SparseRing::new(WordPrimeField::new(P).unwrap(), n)
}

fn poly(ring: &Fp, terms: &[(Vec<u64>, i64)]) -> Sparse<u64> {
    ring.from_terms(
        terms
            .iter()
            .map(|(e, c)| (e.clone(), ring.base().elem_integer(&Integer::from(*c)))),
    )
}

/// `p` over its leading coefficient.
fn monic(ring: &Fp, p: &Sparse<u64>) -> Sparse<u64> {
    let inverse = ring.base().inv(&p.coeffs()[0]).unwrap();
    ring.scale(p, &inverse)
}

/// `g u` and `g v`, and their monic gcd.
fn problem(ring: &Fp, g: &Sparse<u64>, u: &Sparse<u64>, v: &Sparse<u64>) -> [Sparse<u64>; 3] {
    let (a, b) = (ring.mul(g, u), ring.mul(g, v));
    let truth = pgcd(ring, &a, &b, &Never).unwrap().unwrap();
    [a, b, truth]
}

/// The image of the gcd of `a` and `b` from `form`, made monic.
fn image_from(
    ring: &Fp,
    a: &Sparse<u64>,
    b: &Sparse<u64>,
    form: &Form,
    seed: u64,
) -> Result<Sparse<u64>, Failure> {
    let mut random = Random::new(seed);
    let got = image(ring, a, b, form, &mut random, &Never).unwrap();
    got.map(|p| monic(ring, &p))
}

/// Non-zero small coefficients.
fn coefficient() -> impl Strategy<Value = i64> {
    (-50i64..50).prop_filter("non-zero", |c| *c != 0)
}

/// Terms in `n` variables with exponents below 4.
fn terms(n: usize, count: std::ops::Range<usize>) -> impl Strategy<Value = Vec<(Vec<u64>, i64)>> {
    prop::collection::vec((prop::collection::vec(0u64..4, n), coefficient()), count)
}

/// A gcd in `n` variables whose coefficients in `x₀`, at two to four powers,
/// have two to four terms each, and two cofactors.
type Blocks = (
    usize,
    Vec<Vec<(Vec<u64>, i64)>>,
    Vec<(Vec<u64>, i64)>,
    Vec<(Vec<u64>, i64)>,
);

fn blocks() -> impl Strategy<Value = Blocks> {
    (2usize..=4).prop_flat_map(|n| {
        let monomials = prop::collection::btree_set(prop::collection::vec(0u64..4, n - 1), 2..5);
        let block = (monomials, prop::collection::vec(coefficient(), 4)).prop_map(
            |(ms, cs): (BTreeSet<Vec<u64>>, Vec<i64>)| {
                ms.into_iter().zip(cs).collect::<Vec<(Vec<u64>, i64)>>()
            },
        );
        (
            Just(n),
            prop::collection::vec(block, 2..=4),
            terms(n, 1..5),
            terms(n, 1..5),
        )
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Single scaling: the gcd has a term `x₀⁵`, alone at its power.
    #[test]
    fn single_scaling(
        (n, g, u, v) in (2usize..=4).prop_flat_map(|n| (Just(n), terms(n, 1..8), terms(n, 1..5), terms(n, 1..5))),
        seed in any::<u64>(),
    ) {
        let r = ring(n);
        let mut top = vec![0; n];
        top[0] = 5;
        let g = r.add(&poly(&r, &g), &r.monomial(1, &top));
        let (u, v) = (poly(&r, &u), poly(&r, &v));
        prop_assume!(!u.is_zero() && !v.is_zero());
        let [a, b, truth] = problem(&r, &g, &u, &v);
        let vars: Vec<usize> = (1..n).collect();
        let form = Form::of(&truth, 0, &vars);
        prop_assume!(form.single().is_some());
        prop_assert_eq!(image_from(&r, &a, &b, &form, seed), Ok(truth));
    }

    /// Multiple scaling: every coefficient of `x₀` in the gcd has two terms
    /// or more, and the gcd has no content in `x₀`.
    #[test]
    fn multiple_scaling((n, g, u, v) in blocks(), seed in any::<u64>()) {
        let r = ring(n);
        let g: Vec<(Vec<u64>, i64)> = g
            .into_iter()
            .enumerate()
            .flat_map(|(d, block)| {
                block.into_iter().map(move |(m, c)| {
                    let mut e = vec![d as u64];
                    e.extend(m);
                    (e, c)
                })
            })
            .collect();
        let (u, v) = (poly(&r, &u), poly(&r, &v));
        prop_assume!(!u.is_zero() && !v.is_zero());
        let [a, b, truth] = problem(&r, &poly(&r, &g), &u, &v);
        let vars: Vec<usize> = (1..n).collect();
        let form = Form::of(&truth, 0, &vars);
        prop_assume!(form.single().is_none());
        let content = r.content_in(&truth, 0);
        prop_assume!(content.degrees().iter().all(|&d| d == 0));
        prop_assert_eq!(image_from(&r, &a, &b, &form, seed), Ok(truth));
    }
}

/// A form missing a term is caught by the extra image, in both scaling
/// cases; one with a term too many still gives the gcd; one that lacks the
/// gcd's top power makes every point look unlucky.
#[test]
fn wrong_forms() {
    let r = ring(3);
    let (u, v) = (
        poly(
            &r,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 0], 1)],
        ),
        poly(
            &r,
            &[(vec![1, 0, 0], 1), (vec![0, 0, 1], -1), (vec![0, 0, 0], 2)],
        ),
    );
    // (y + 2z) x² + (y² + 3) x + z + 5: multiple scaling; and with x³,
    // single.
    let multiple = poly(
        &r,
        &[
            (vec![2, 1, 0], 1),
            (vec![2, 0, 1], 2),
            (vec![1, 2, 0], 1),
            (vec![1, 0, 0], 3),
            (vec![0, 0, 1], 1),
            (vec![0, 0, 0], 5),
        ],
    );
    let single = r.add(&multiple, &poly(&r, &[(vec![3, 0, 0], 7)]));
    for g in [&multiple, &single] {
        let [a, b, truth] = problem(&r, g, &u, &v);
        let form = Form::of(&truth, 0, &[1, 2]);
        assert_eq!(image_from(&r, &a, &b, &form, 1), Ok(truth.clone()));
        // Without y² in the coefficient of x.
        let mut missing = form.clone();
        let block = missing.blocks.iter_mut().find(|(d, _)| *d == 1).unwrap();
        block.1.retain(|m| m != &[2, 0]);
        assert_eq!(image_from(&r, &a, &b, &missing, 2), Err(Failure::Form));
        // With y z in the constant coefficient.
        let mut extra = form.clone();
        let block = extra.blocks.iter_mut().find(|(d, _)| *d == 0).unwrap();
        block.1.insert(0, vec![1, 1]);
        assert_eq!(image_from(&r, &a, &b, &extra, 3), Ok(truth.clone()));
        // Without the top power.
        let mut low = form.clone();
        low.blocks.remove(0);
        assert_eq!(image_from(&r, &a, &b, &low, 4), Err(Failure::Points));
    }
}

/// De Kleine, Monagan and Wittkopf's example 4 [DKMW05]: `(y + 2) x³ +
/// 12 y² + 24 y` has the content `y + 2` in `x`, which the scales absorb.
#[test]
fn a_content_leaves_the_scales_free() {
    let r = ring(2);
    let g = poly(
        &r,
        &[
            (vec![3, 1], 1),
            (vec![3, 0], 2),
            (vec![0, 2], 12),
            (vec![0, 1], 24),
        ],
    );
    let (u, v) = (
        poly(&r, &[(vec![1, 0], 1), (vec![0, 1], 1)]),
        poly(&r, &[(vec![1, 0], 1), (vec![0, 1], -1), (vec![0, 0], 3)]),
    );
    let [a, b, truth] = problem(&r, &g, &u, &v);
    let form = Form::of(&truth, 0, &[1]);
    assert_eq!(form.single(), None);
    assert_eq!(image_from(&r, &a, &b, &form, 5), Err(Failure::Content));
}

/// The main variable need not be the first: `y` in `x, y, z`.
#[test]
fn another_main_variable() {
    let r = ring(3);
    let g = poly(
        &r,
        &[
            (vec![2, 1, 0], 3),
            (vec![0, 3, 1], 1),
            (vec![1, 1, 1], -4),
            (vec![0, 0, 2], 9),
            (vec![1, 0, 0], 2),
        ],
    );
    let (u, v) = (
        poly(
            &r,
            &[(vec![0, 1, 0], 1), (vec![1, 0, 0], 1), (vec![0, 0, 0], 1)],
        ),
        poly(&r, &[(vec![0, 2, 0], 1), (vec![0, 0, 1], -1)]),
    );
    let [a, b, truth] = problem(&r, &g, &u, &v);
    let form = Form::of(&truth, 1, &[0, 2]);
    assert_eq!(image_from(&r, &a, &b, &form, 6), Ok(truth));
}
