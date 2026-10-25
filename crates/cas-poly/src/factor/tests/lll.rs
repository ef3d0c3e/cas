//! Tests of the integral LLL reduction (CAS-188).

use cas_num::{Integer, Rational};
use proptest::prelude::*;

use crate::cancel::Never;
use crate::factor::lll::Lattice;

fn dot(a: &[Rational], b: &[Rational]) -> Rational {
    a.iter()
        .zip(b)
        .fold(Rational::ZERO, |acc, (x, y)| &acc + &(x * y))
}

/// The Gram-Schmidt vectors and coefficients `μ`, over ℚ.
fn gram_schmidt(basis: &[Vec<Integer>]) -> (Vec<Vec<Rational>>, Vec<Vec<Rational>>) {
    let rows: Vec<Vec<Rational>> = basis
        .iter()
        .map(|r| r.iter().map(|x| Rational::from(x.clone())).collect())
        .collect();
    let n = rows.len();
    let mut stars: Vec<Vec<Rational>> = Vec::with_capacity(n);
    let mut mu = vec![vec![Rational::ZERO; n]; n];
    for i in 0..n {
        let mut v = rows[i].clone();
        for j in 0..i {
            mu[i][j] = &dot(&rows[i], &stars[j]) / &dot(&stars[j], &stars[j]);
            for (x, y) in v.iter_mut().zip(&stars[j]) {
                *x = &*x - &(&mu[i][j] * y);
            }
        }
        stars.push(v);
    }
    (stars, mu)
}

/// The rank over ℚ of the rows together.
fn rank(rows: &[Vec<Integer>]) -> usize {
    let mut m: Vec<Vec<Rational>> = rows
        .iter()
        .map(|r| r.iter().map(|x| Rational::from(x.clone())).collect())
        .collect();
    let width = m.first().map_or(0, Vec::len);
    let mut rank = 0;
    for col in 0..width {
        let Some(p) = (rank..m.len()).find(|&i| !m[i][col].is_zero()) else {
            continue;
        };
        m.swap(rank, p);
        let pivot = m[rank].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i != rank && !row[col].is_zero() {
                let f = &row[col] / &pivot[col];
                for (x, y) in row.iter_mut().zip(&pivot) {
                    *x = &*x - &(&f * y);
                }
            }
        }
        rank += 1;
    }
    rank
}

/// The basis, reduced, checked: size-reduced, the Lovász condition with
/// `δ = 3/4`, and the same lattice: each basis spans the other's rows with
/// integer coefficients, which for bases of one dimension means the reduced
/// rows are in the lattice of the input and the Gram determinants agree.
fn check(basis: &[Vec<Integer>]) -> Result<(), TestCaseError> {
    let reduced = Lattice::reduce(basis.to_vec(), &Never)
        .unwrap()
        .into_basis();
    prop_assert_eq!(reduced.len(), basis.len());
    let (stars, mu) = gram_schmidt(&reduced);
    let half = Rational::new(Integer::ONE, Integer::from(2));
    let three_quarters = Rational::new(Integer::from(3), Integer::from(4));
    for i in 0..reduced.len() {
        for m in mu[i].iter().take(i) {
            prop_assert!(m.abs() <= half, "size-reduced");
        }
        if i > 0 {
            let lhs = dot(&stars[i], &stars[i]);
            let rhs = &(&three_quarters - &(&mu[i][i - 1] * &mu[i][i - 1]))
                * &dot(&stars[i - 1], &stars[i - 1]);
            prop_assert!(lhs >= rhs, "Lovász's condition");
        }
    }
    // The same lattice: one rank for both together, and the same volume.
    let mut both = basis.to_vec();
    both.extend(reduced.iter().cloned());
    prop_assert_eq!(rank(&both), basis.len());
    let volume = |b: &[Vec<Integer>]| {
        let (s, _) = gram_schmidt(b);
        s.iter().fold(Rational::ONE, |acc, v| &acc * &dot(v, v))
    };
    prop_assert_eq!(volume(&reduced), volume(basis));
    Ok(())
}

#[test]
fn examples() {
    let int = |rows: &[&[i64]]| -> Vec<Vec<Integer>> {
        rows.iter()
            .map(|r| r.iter().map(|&x| Integer::from(x)).collect())
            .collect()
    };
    // A knapsack: the identity beside a column, and the modulus.
    let basis = int(&[
        &[1, 0, 0, 0, 3_000_000_001],
        &[0, 1, 0, 0, 5_000_000_002],
        &[0, 0, 1, 0, 7_000_000_003],
        &[0, 0, 0, 1, 2_000_000_000],
        &[0, 0, 0, 0, 10_000_000_000],
    ]);
    check(&basis).unwrap();
    // One vector, and a basis already reduced.
    check(&int(&[&[3, 4]])).unwrap();
    check(&int(&[&[1, 0], &[0, 1]])).unwrap();
    // A long vector first.
    check(&int(&[&[201, 37], &[1648, 297]])).unwrap();
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// Random bases of 1 to 7 independent vectors in up to 8 dimensions.
    #[test]
    fn random_bases(rows in prop::collection::vec(prop::collection::vec(-60i64..60, 8), 1..8), width in 1usize..9) {
        let basis: Vec<Vec<Integer>> = rows
            .iter()
            .map(|r| r.iter().take(width.max(rows.len())).map(|&x| Integer::from(x)).collect())
            .collect();
        prop_assume!(rank(&basis) == basis.len());
        check(&basis)?;
    }
}
