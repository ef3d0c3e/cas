//! Factoring univariate polynomials over ℤ and ℚ by Zassenhaus's algorithm
//! (CAS-182) \[Zas69, vzGG13 §15.6\]. Each square-free part's local factors
//! modulo a small prime are lifted past twice Mignotte's bound \[Mig74\] and
//! recombined by subsets.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithm: f = lc(f) ∏ G modulo pᵏ, subsets of size s"
)]

use cas_num::algebra::{EuclideanDomain, Integers, Rationals, Ring, WordPrimeField};
use cas_num::{Integer, Natural, Rational};

use super::finite::{distinct_degree, split_degrees};
use super::hensel::lift;
use super::{FactorError, Factors, Search, unbounded};
use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::dense::{Dense, DenseRing};
use crate::gcd::split;

/// How many primes that keep a polynomial square-free are tried. The one
/// with the fewest local factors is lifted.
const PRIMES: usize = 3;

/// A prime and the distinct-degree factorisation of the polynomial modulo
/// it, with the number of local factors.
struct Prime {
    field: WordPrimeField,
    ddf: Vec<(Dense<u64>, usize)>,
    count: usize,
}

/// The degrees from 0 to `n` that sums of the local factors' degrees reach:
/// the degrees a factor over ℤ can have.
fn degree_sums(ddf: &[(Dense<u64>, usize)], n: usize) -> Vec<bool> {
    let mut sums = vec![false; n + 1];
    sums[0] = true;
    for (g, d) in ddf {
        for _ in 0..g.degree().unwrap_or(0) / d {
            for s in (*d..=n).rev() {
                if sums[s - d] {
                    sums[s] = true;
                }
            }
        }
    }
    sums
}

/// The prime to lift from: among the first [`PRIMES`] odd primes that do
/// not divide `lc(f)` and keep `f` square-free, the one with the fewest
/// local factors. With it come the degrees a factor of `f` can have, those
/// that the degree sums modulo every prime tried allow. `None` when they
/// prove `f` irreducible: one local factor, or no degree strictly between 0
/// and `deg f`.
fn choose_prime(
    f: &Dense<Integer>,
    cancel: &dyn Cancel,
) -> Result<Option<(Prime, Vec<bool>)>, Cancelled> {
    let n = f.degree().unwrap_or(0);
    let mut allowed = vec![true; n + 1];
    let mut best: Option<Prime> = None;
    let mut found = 0;
    let mut p = 3;
    while found < PRIMES {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let candidate = WordPrimeField::new(p);
        p += 2;
        let Some(field) = candidate else {
            continue;
        };
        let fp = DenseRing::new(field.clone());
        let image = fp.poly(f.coeffs().iter().map(|c| field.elem_integer(c)).collect());
        if image.degree() != Some(n) {
            continue;
        }
        let monic = EuclideanDomain::normal_form(&fp, &image).0;
        if fp
            .gcd_with(&monic, &fp.derivative(&monic), cancel)?
            .degree()
            != Some(0)
        {
            continue;
        }
        found += 1;
        let ddf = distinct_degree(&fp, &monic, cancel)?;
        let count = ddf.iter().map(|(g, d)| g.degree().unwrap_or(0) / d).sum();
        for (a, b) in allowed.iter_mut().zip(degree_sums(&ddf, n)) {
            *a &= b;
        }
        if count == 1 || (1..n).all(|d| !allowed[d]) {
            return Ok(None);
        }
        if best.as_ref().is_none_or(|b| count < b.count) {
            best = Some(Prime { field, ddf, count });
        }
    }
    Ok(best.map(|prime| (prime, allowed)))
}

/// The least `k` with `pᵏ > 2B`, for `B = 2^(n−1) ‖f‖₂` and `n = deg f`.
/// Let `g` be a factor of `f` of degree below `n`. By Mignotte's bound
/// \[Mig74\], the coefficients of `lc(f)/lc(g) · g` are within `B`. So it is
/// the symmetric residue modulo `pᵏ` of `lc(f) ∏ Gᵢ`, the product over the
/// lifts of its local factors. That still holds after factors of `f` are
/// divided out, since what is left has a leading coefficient that divides
/// `lc(f)`.
fn exponent(f: &Dense<Integer>, p: u64) -> u64 {
    let n = f.degree().unwrap_or(0) as u64;
    let squares = f.coeffs().iter().fold(Natural::ZERO, |acc, c| {
        let a = c.unsigned_abs();
        &acc + &(&a * &a)
    });
    // 2B, with ⌊‖f‖₂²^(1/2)⌋ + 1 above ‖f‖₂.
    let twice_bound = (squares.isqrt() + 1) << n;
    let p = Natural::from(p);
    let mut power = p.clone();
    let mut k = 1;
    while power <= twice_bound {
        power = &power * &p;
        k += 1;
    }
    k
}

/// `c` modulo `m` in `(−m/2, m/2]`.
fn symmetric(c: &Integer, m: &Integer) -> Integer {
    let r = c.div_rem_euclid(m).1;
    if &r + &r > *m { r - m } else { r }
}

/// The next subset of `subset.len()` elements of `0..r` in lexicographic
/// order; `false` after the last.
fn next_subset(subset: &mut [usize], r: usize) -> bool {
    let s = subset.len();
    for i in (0..s).rev() {
        if subset[i] < r - s + i {
            subset[i] += 1;
            for j in i + 1..s {
                subset[j] = subset[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// Zassenhaus's recombination of the lifts `lifted` modulo `m` of the
/// local factors of `f`. Subsets of the lifts are tried by increasing size.
/// A subset is a candidate when its degree is allowed and its constant
/// terms pass the test: `lc(f) ∏ Gᵢ(0)` modulo `m` must divide
/// `lc(f) f(0)`. A candidate's factor is the primitive part of `lc(f) ∏ Gᵢ`
/// modulo `m`, kept when it divides `f`; then its lifts are removed. Once no
/// subset of at most half the lifts left gives a factor, none does, and what
/// is left of `f` is irreducible.
fn recombine(
    mut f: Dense<Integer>,
    mut lifted: Vec<Dense<Integer>>,
    m: &Integer,
    allowed: &[bool],
    search: &Search<'_>,
) -> Result<Vec<Dense<Integer>>, FactorError> {
    let zx = DenseRing::new(Integers);
    let mut out = Vec::new();
    let mut s = 1;
    'sizes: while 2 * s <= lifted.len() {
        let r = lifted.len();
        let lead = f.lead().cloned().unwrap_or(Integer::ONE);
        let constant = &lead * f.coeff(0).unwrap_or(&Integer::ZERO);
        let mut subset: Vec<usize> = (0..s).collect();
        loop {
            search.step()?;
            // At half the lifts a subset and its complement give one split:
            // only the subsets with the first are tried.
            if 2 * s == r && subset[0] != 0 {
                break;
            }
            let degree: usize = subset
                .iter()
                .map(|&i| lifted[i].degree().unwrap_or(0))
                .sum();
            let c = || {
                subset.iter().fold(lead.clone(), |acc, &i| {
                    symmetric(&(&acc * lifted[i].coeff(0).unwrap_or(&Integer::ZERO)), m)
                })
            };
            if allowed[degree] && constant.is_divisible_by(&c()) {
                let mut product = zx.constant(lead.clone());
                for &i in &subset {
                    let full = zx.mul_with(&product, &lifted[i], search.cancel())?;
                    product = zx.poly(full.coeffs().iter().map(|x| symmetric(x, m)).collect());
                }
                let g = zx.primitive_part(&product);
                if let Some(q) = zx.div_exact(&f, &g) {
                    out.push(g);
                    f = q;
                    for &i in subset.iter().rev() {
                        lifted.remove(i);
                    }
                    continue 'sizes;
                }
            }
            if !next_subset(&mut subset, r) {
                break;
            }
        }
        s += 1;
    }
    if f.degree() > Some(0) {
        out.push(f);
    }
    Ok(out)
}

/// The irreducible factors of square-free, primitive `f` of positive degree
/// with a positive leading coefficient, each primitive with a positive
/// leading coefficient. A factor `x` is split off first, so that `f(0)` is
/// not 0 for the constant-term test. Each subset of the recombination is a
/// step of `search`.
pub(crate) fn factor_squarefree(
    f: &Dense<Integer>,
    search: &Search<'_>,
) -> Result<Vec<Dense<Integer>>, FactorError> {
    let cancel = search.cancel();
    let zx = DenseRing::new(Integers);
    let mut out = Vec::new();
    let mut f = f.clone();
    if f.coeff(0).is_some_and(Integer::is_zero) {
        out.push(zx.x());
        f = zx.poly(f.coeffs()[1..].to_vec());
    }
    if f.degree().unwrap_or(0) <= 1 {
        if f.degree() == Some(1) {
            out.push(f);
        }
        return Ok(out);
    }
    let Some((prime, allowed)) = choose_prime(&f, cancel)? else {
        out.push(f);
        return Ok(out);
    };
    let fp = DenseRing::new(prime.field);
    let local = split_degrees(&fp, &prime.ddf, cancel)?;
    let p = u64::try_from(&fp.characteristic()).expect("a word prime");
    let k = exponent(&f, p);
    let lifted = lift(&f, &local, &fp, k, cancel)?;
    let m = Integer::from(Natural::from(p).pow(k));
    out.extend(recombine(f, lifted, &m, &allowed, search)?);
    Ok(out)
}

/// The factors in the order of `factor` over ℤ: by degree, then by
/// coefficients from the constant term up, then by exponent.
fn sort(factors: &mut [(Dense<Integer>, u64)]) {
    factors.sort_by(|(a, e), (b, x)| {
        a.degree()
            .cmp(&b.degree())
            .then_with(|| a.coeffs().cmp(b.coeffs()))
            .then(e.cmp(x))
    });
}

impl DenseRing<Integers> {
    /// The factorisation of `p` into irreducible polynomials over ℤ,
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is the content with the sign of the
    /// leading coefficient. The `fᵢ` are distinct, irreducible and
    /// primitive, with positive leading coefficients, sorted by degree, then
    /// by coefficients from the constant term up. Zero gives `c = 0` and no
    /// factors. giac: `factor(x^6-1)` is `(x-1)*(x+1)*(x^2-x+1)*(x^2+x+1)`.
    #[must_use]
    pub fn factor(&self, p: &Dense<Integer>) -> Factors<Integer, Dense<Integer>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`, also in the recombination, which is
    /// exponential in the number of local factors in the worst case.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Dense<Integer>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Integer, Dense<Integer>>, Cancelled> {
        unbounded(self.factor_within(p, u64::MAX, cancel))
    }

    /// As `factor_with`, trying at most `bound` subsets of local factors in
    /// all.
    ///
    /// # Errors
    ///
    /// [`FactorError::Cancelled`] when `cancel` says so, and
    /// [`FactorError::Search`] past the bound.
    pub fn factor_within(
        &self,
        p: &Dense<Integer>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Integer, Dense<Integer>>, FactorError> {
        self.factor_search(p, &Search::new(cancel, bound))
    }

    /// As `factor_within`, the subsets counted by `search`.
    pub(crate) fn factor_search(
        &self,
        p: &Dense<Integer>,
        search: &Search<'_>,
    ) -> Result<Factors<Integer, Dense<Integer>>, FactorError> {
        let (constant, parts) = self.squarefree_with(p, search.cancel())?.into_parts();
        let mut factors = Vec::new();
        for (f, e) in parts {
            for g in factor_squarefree(&f, search)? {
                factors.push((g, e));
            }
        }
        sort(&mut factors);
        Ok(Factors { constant, factors })
    }
}

impl DenseRing<Rationals> {
    /// The factorisation of `p` into irreducible polynomials over ℚ,
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is the leading coefficient, and the `fᵢ`
    /// are monic, in the order of the factors of the numerator over ℤ. Zero
    /// gives `c = 0` and no factors.
    #[must_use]
    pub fn factor(&self, p: &Dense<Rational>) -> Factors<Rational, Dense<Rational>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Dense<Rational>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Rational, Dense<Rational>>, Cancelled> {
        let Some(lead) = p.lead() else {
            return Ok(Factors {
                constant: Rational::ZERO,
                factors: Vec::new(),
            });
        };
        let zx = DenseRing::new(Integers);
        let numerator = zx.poly(split::<Integers>(p.coeffs()).0);
        let factors = zx
            .factor_with(&numerator, cancel)?
            .into_parts()
            .1
            .into_iter()
            .map(|(g, e)| {
                let lc = g.lead().cloned().unwrap_or(Integer::ONE);
                let monic = g
                    .coeffs()
                    .iter()
                    .map(|c| Rational::new(c.clone(), lc.clone()))
                    .collect();
                (self.poly(monic), e)
            })
            .collect();
        Ok(Factors {
            constant: lead.clone(),
            factors,
        })
    }
}
