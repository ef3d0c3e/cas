//! van Hoeij's recombination of lifted factors (CAS-188) \[vH02\]: the
//! subsets of the lifts that give the true factors, found by lattice
//! reduction instead of trying subsets.
//!
//! For the lifts `Gᵢ` of `f ≡ ℓ ∏ Gᵢ (mod pᵃ)`, let `Lᵢ = (f / Gᵢ) Gᵢ'`, the
//! logarithmic derivative of `Gᵢ` times `f`. A true factor `g` of `f` is
//! `∏_(i∈S) Gᵢ` up to a constant, so `Σ_(i∈S) Lᵢ ≡ f g'/g (mod pᵃ)`, a
//! polynomial over ℤ whose coefficients are small: `f g'/g = Σ_α f/(x − α)`
//! over the roots `α` of `g`, and the coefficient `j` of each `f/(x − α)` is
//! within `C(n−1, j) M(f) ≤ C(n−1, j) ‖f‖₂` \[Mig74\]. So the 0/1 vector of
//! `S` is a short vector of a knapsack lattice: the identity, beside the
//! coefficients of the `Lᵢ` cut to their bound's digits, and `pᵃ` cut alike.
//! LLL reduction leaves the true vectors near the start of the basis, and
//! the vectors at the end whose Gram-Schmidt norm is past the true vectors'
//! cannot be needed. Columns are added one at a time until the first `r`
//! coordinates of what is left, in reduced echelon form, partition the
//! lifts: those are the factors, checked by division.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of van Hoeij's paper: f = ℓ ∏ G, L for logarithmic derivatives, p^a, r lifts"
)]

use cas_num::algebra::{Integers, Ring};
use cas_num::modular::{ModRing, Residue};
use cas_num::{Integer, Natural, Rational};

use super::lll::Lattice;
use super::{FactorError, Search};
use crate::dense::{Dense, DenseRing};

/// The number of lifts past which van Hoeij's recombination replaces
/// Zassenhaus's: below it, the subsets are few.
pub(crate) const THRESHOLD: usize = 8;

/// `p^e`.
fn power(p: u64, e: u64) -> Integer {
    Integer::from(Natural::from(p).pow(e))
}

/// The least `e` with `pᵉ ≥ x`, at least 1.
fn digits_for(x: &Natural, p: u64) -> u64 {
    let p = Natural::from(p);
    let (mut e, mut power) = (1, p.clone());
    while &power < x {
        power = &power * &p;
        e += 1;
    }
    e
}

/// `⌈‖f‖₂⌉`.
fn norm(f: &Dense<Integer>) -> Natural {
    let squares = f.coeffs().iter().fold(Natural::ZERO, |acc, c| {
        let a = c.unsigned_abs();
        &acc + &(&a * &a)
    });
    squares.isqrt() + 1u64
}

/// The bounds `n C(n−1, j) ‖f‖₂` on coefficient `j` of `f g'/g` for the
/// factors `g` of `f`.
fn bounds(f: &Dense<Integer>) -> Vec<Natural> {
    let n = f.degree().unwrap_or(0);
    let scale = &norm(f) * &Natural::from(n as u64);
    let mut binomial = Natural::one();
    (0..n)
        .map(|j| {
            let bound = &scale * &binomial;
            // C(n−1, j + 1) = C(n−1, j) (n − 1 − j) / (j + 1).
            binomial = (&binomial * &Natural::from((n - 1 - j) as u64))
                .div_exact(&Natural::from((j + 1) as u64));
            bound
        })
        .collect()
}

/// The columns in the order they are added: coefficients from both ends of
/// `f g'/g`, `n − 1, 0, n − 2, 1, …`, whose bounds are the smallest.
fn columns(n: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(n);
    let (mut low, mut high) = (0, n);
    while low < high {
        high -= 1;
        out.push(high);
        if low < high {
            out.push(low);
            low += 1;
        }
    }
    out
}

/// The digits of `p` a column needs beyond its bound for LLL to tell the
/// true vectors apart: `2^((r + 1)/2)` of LLL's approximation, a bound on
/// the true vectors, and a margin.
pub(crate) fn margin(r: usize, p: u64) -> u64 {
    let r = r as u64;
    let bits = r / 2 + 2 * u64::from(u64::BITS - (r + 1).leading_zeros()) + 24;
    digits_for(&Natural::power_of_two(bits), p)
}

/// The digits of `p` van Hoeij's recombination wants for `f` with `r`
/// lifts: those of the first four columns' bounds, and the margin.
pub(crate) fn precision(f: &Dense<Integer>, r: usize, p: u64) -> u64 {
    let n = f.degree().unwrap_or(0);
    let bounds = bounds(f);
    let first = columns(n)
        .into_iter()
        .take(4)
        .map(|j| digits_for(&bounds[j], p))
        .max()
        .unwrap_or(1);
    first + margin(r, p)
}

/// `x` modulo `m` in `(−m/2, m/2]`.
fn symmetric(x: &Integer, m: &Integer) -> Integer {
    let r = x.div_rem_euclid(m).1;
    if &r + &r > *m { r - m } else { r }
}

/// The logarithmic derivatives `Lᵢ = (f / Gᵢ) Gᵢ'` modulo `m`, their
/// coefficients in `[0, m)`.
fn derivatives(f: &Dense<Integer>, lifted: &[Dense<Integer>], m: &Integer) -> Vec<Vec<Integer>> {
    let ring = DenseRing::new(ModRing::new(m.unsigned_abs()));
    let base = ring.base();
    let residues = |p: &Dense<Integer>| -> Dense<Residue> {
        ring.poly(p.coeffs().iter().map(|c| base.elem_integer(c)).collect())
    };
    let n = f.degree().unwrap_or(0);
    let f = residues(f);
    lifted
        .iter()
        .map(|g| {
            let g = residues(g);
            let (quotient, _) = ring.pseudo_div_rem(&f, &g);
            let l = ring.mul(&quotient, &ring.derivative(&g));
            (0..n)
                .map(|j| {
                    l.coeff(j)
                        .map_or(Integer::ZERO, |c| Integer::from(base.to_natural(c)))
                })
                .collect()
        })
        .collect()
}

/// The reduced echelon form of the rows over ℚ, when it is a partition of
/// the columns into 0/1 rows: the column sets.
fn partition(rows: &[Vec<Integer>], r: usize) -> Option<Vec<Vec<usize>>> {
    let mut m: Vec<Vec<Rational>> = rows
        .iter()
        .map(|row| {
            row.iter()
                .take(r)
                .map(|x| Rational::from(x.clone()))
                .collect()
        })
        .collect();
    let mut pivot_row = 0;
    for col in 0..r {
        let Some(found) = (pivot_row..m.len()).find(|&i| !m[i][col].is_zero()) else {
            continue;
        };
        m.swap(pivot_row, found);
        let pivot = m[pivot_row][col].clone();
        for x in &mut m[pivot_row] {
            *x = &*x / &pivot;
        }
        let row = m[pivot_row].clone();
        for (i, other) in m.iter_mut().enumerate() {
            if i != pivot_row && !other[col].is_zero() {
                let factor = other[col].clone();
                for (x, y) in other.iter_mut().zip(&row) {
                    *x = &*x - &(&factor * y);
                }
            }
        }
        pivot_row += 1;
        if pivot_row == m.len() {
            break;
        }
    }
    let mut covered = vec![false; r];
    let mut groups = Vec::with_capacity(m.len());
    for row in &m {
        let mut group = Vec::new();
        for (i, x) in row.iter().enumerate() {
            if x.is_one() {
                if covered[i] {
                    return None;
                }
                covered[i] = true;
                group.push(i);
            } else if !x.is_zero() {
                return None;
            }
        }
        if group.is_empty() {
            return None;
        }
        groups.push(group);
    }
    covered.iter().all(|&c| c).then_some(groups)
}

/// The factors of `f` that the groups of lifts give, `pp(ℓ ∏ Gᵢ mods m)`,
/// when each divides `f` and they multiply to it; `None` otherwise.
fn factors_of(
    f: &Dense<Integer>,
    lifted: &[Dense<Integer>],
    groups: &[Vec<usize>],
    m: &Integer,
) -> Option<Vec<Dense<Integer>>> {
    let zx = DenseRing::new(Integers);
    let lead = f.lead()?.clone();
    let mut rest = f.clone();
    let mut out = Vec::with_capacity(groups.len());
    for group in groups {
        let mut product = zx.constant(lead.clone());
        for &i in group {
            let full = zx.mul(&product, &lifted[i]);
            product = zx.poly(full.coeffs().iter().map(|x| symmetric(x, m)).collect());
        }
        let g = zx.primitive_part(&product);
        rest = zx.div_exact(&rest, &g)?;
        out.push(g);
    }
    (rest.degree() == Some(0)).then_some(out)
}

/// The irreducible factors of `f`, square-free, primitive, with a positive
/// leading coefficient, from the lifts `lifted` modulo `pᵃ` of its local
/// factors, by van Hoeij's knapsack; `None` when the columns run out at this
/// precision. Each column is a step of `search`.
pub(crate) fn recombine(
    f: &Dense<Integer>,
    lifted: &[Dense<Integer>],
    p: u64,
    a: u64,
    search: &Search<'_>,
) -> Result<Option<Vec<Dense<Integer>>>, FactorError> {
    let r = lifted.len();
    let n = f.degree().unwrap_or(0);
    let m = power(p, a);
    let derivatives = derivatives(f, lifted, &m);
    let bounds = bounds(f);
    let margin = margin(r, p);
    // The basis: the identity, beside the columns added so far.
    let mut basis: Vec<Vec<Integer>> = (0..r)
        .map(|i| (0..r).map(|k| Integer::from(u64::from(i == k))).collect())
        .collect();
    let mut added = 0u64;
    for j in columns(n) {
        let b = digits_for(&bounds[j], p);
        if b + margin > a {
            continue;
        }
        search.step()?;
        let cut = power(p, b);
        let column: Vec<Integer> = derivatives
            .iter()
            .map(|l| (&l[j] - &symmetric(&l[j], &cut)).div_exact(&cut))
            .collect();
        let width = basis.first().map_or(r, Vec::len);
        let mut rows: Vec<Vec<Integer>> = basis
            .iter()
            .map(|row| {
                let value = row
                    .iter()
                    .take(r)
                    .zip(&column)
                    .fold(Integer::ZERO, |acc, (x, c)| &acc + &(x * c));
                let mut row = row.clone();
                row.push(value);
                row
            })
            .collect();
        let mut modulus = vec![Integer::ZERO; width];
        modulus.push(power(p, a - b));
        rows.push(modulus);
        added += 1;
        let lattice = Lattice::reduce(rows, search.cancel())?;
        // The true vectors: 0/1 on the lifts, each column within r/2 + 1.
        let half = Integer::from(r as u64 / 2 + 1);
        let bound = &Integer::from(r as u64) + &(&Integer::from(added) * &(&half * &half));
        let mut keep = lattice.basis().len();
        while keep > 1 && lattice.gram_schmidt_exceeds(keep - 1, &bound) {
            keep -= 1;
        }
        basis = lattice.into_basis();
        basis.truncate(keep);
        if let Some(groups) = partition(&basis, r)
            && let Some(factors) = factors_of(f, lifted, &groups, &m)
        {
            return Ok(Some(factors));
        }
    }
    Ok(None)
}
