//! Factoring over ℚ(i) by Trager's norms (CAS-186) \[Tra76\]. A square-free
//! `f` over ℤ\[i\], primitive in a main variable `x`, is shifted to
//! `f(x − s i)` until its norm `N = f(x − s i) · f̄(x + s i)`, which has
//! integer coefficients, is square-free. Then each irreducible factor of `N`
//! over ℚ has the gcd of an irreducible factor of the shifted `f` with it,
//! and those gcds, shifted back, are the factors of `f`.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Trager's paper: f, its norm N, the shift s, the factors g and h"
)]

use std::cmp::Ordering;

use cas_num::algebra::{EuclideanDomain, GaussianIntegers, GaussianRationals, Integers, Ring};
use cas_num::{Complex, Rational};

use super::{FactorError, Factors, Search, unbounded};
use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::gcd::split_sparse;
use crate::sparse::{Sparse, SparseRing};

/// The variables that `p` has.
fn variables<E>(p: &Sparse<E>) -> Vec<usize> {
    p.degrees()
        .iter()
        .enumerate()
        .filter(|&(_, &d)| d > 0)
        .map(|(v, _)| v)
        .collect()
}

/// The shifts tried: 0, 1, −1, 2, −2, …
fn shift(k: u64) -> i64 {
    let half = i64::try_from(k.div_ceil(2)).unwrap_or(i64::MAX);
    if k % 2 == 1 { half } else { -half }
}

/// The total degree.
fn total_degree<E>(p: &Sparse<E>) -> u64 {
    (0..p.term_count())
        .map(|k| p.exponents(k).iter().sum())
        .max()
        .unwrap_or(0)
}

/// The order of the factors: by total degree, then by terms from the
/// leading one, each by exponents, then real and imaginary parts.
fn compare(a: &Sparse<Complex>, b: &Sparse<Complex>) -> Ordering {
    total_degree(a).cmp(&total_degree(b)).then_with(|| {
        a.terms()
            .zip(b.terms())
            .map(|((e, c), (f, d))| {
                e.cmp(&f)
                    .then_with(|| c.re().cmp(d.re()))
                    .then_with(|| c.im().cmp(d.im()))
            })
            .find(|o| o.is_ne())
            .unwrap_or_else(|| a.term_count().cmp(&b.term_count()))
    })
}

impl SparseRing<GaussianIntegers> {
    /// The factorisation of `p` into irreducible polynomials over ℚ(i),
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is a Gaussian integer, the content times a
    /// unit. The `fᵢ` are distinct and irreducible over ℚ(i), primitive over
    /// ℤ\[i\], with their lexicographic leading coefficients in the first
    /// quadrant, sorted by total degree, then by terms. giac:
    /// `factor(x^2+(1+2*i)*x+2*i)` is `(x+1)*(x+2*i)`, and
    /// `factor(i*x^2+i)` is `i*(x+i)*(x-i)`.
    #[must_use]
    pub fn factor(&self, p: &Sparse<Complex>) -> Factors<Complex, Sparse<Complex>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Sparse<Complex>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Complex, Sparse<Complex>>, Cancelled> {
        unbounded(self.factor_within(p, u64::MAX, cancel))
    }

    /// As `factor_with`, with at most `bound` steps of the searches in all:
    /// the shifts, and those of factoring the norms over ℤ.
    ///
    /// # Errors
    ///
    /// [`FactorError::Cancelled`] when `cancel` says so, and
    /// [`FactorError::Search`] past the bound.
    pub fn factor_within(
        &self,
        p: &Sparse<Complex>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Complex, Sparse<Complex>>, FactorError> {
        let search = Search::new(cancel, bound);
        let base = self.base();
        let (mut constant, parts) = self.squarefree_with(p, cancel)?.into_parts();
        let mut factors = Vec::new();
        for (f, e) in parts {
            let gs = self.trager(&f, &search)?;
            // f is a unit times the product of its normalised factors.
            let leads = gs
                .iter()
                .fold(base.one(), |acc, g| base.mul(&acc, &g.coeffs()[0]));
            let unit = base.div_rem(&f.coeffs()[0], &leads).0;
            constant = base.mul(&constant, &base.pow(&unit, e));
            factors.extend(gs.into_iter().map(|g| (g, e)));
        }
        factors.sort_by(|(a, e), (b, x)| compare(a, b).then(e.cmp(x)));
        Ok(Factors { constant, factors })
    }

    /// The irreducible factors over ℚ(i) of square-free `f`, primitive over
    /// ℤ\[i\], normalised. Its content in its first variable is factored in
    /// the others; the primitive part, which has no factor without that
    /// variable, by Trager's norms.
    fn trager(
        &self,
        f: &Sparse<Complex>,
        search: &Search<'_>,
    ) -> Result<Vec<Sparse<Complex>>, FactorError> {
        let cancel = search.cancel();
        let Some(&x) = variables(f).first() else {
            return Ok(Vec::new());
        };
        let (content, primitive) = self.content_split(f, x, cancel)?;
        let mut out = if self.is_one(&content) {
            Vec::new()
        } else {
            self.trager(&content, search)?
        };
        let integers = SparseRing::new(Integers, self.nvars());
        let i = Complex::new(Rational::ZERO, Rational::ONE);
        for k in 0.. {
            search.step()?;
            let s = Complex::from(shift(k));
            // f(x - s i) and its norm.
            let moved = |p: &Sparse<Complex>, s: &Complex| {
                let image = self.sub(&self.var(x), &self.constant(s * &i));
                self.substitute(p, x, &image)
            };
            let shifted = moved(&primitive, &s);
            let conjugate = self.map_into(&shifted, self, Complex::conj);
            let norm = self.mul_with(&shifted, &conjugate, cancel)?;
            debug_assert!(norm.coeffs().iter().all(Complex::is_real), "a norm");
            let norm = self.map_into(&norm, &integers, |c| c.re().numer().clone());
            if integers
                .squarefree_with(&norm, cancel)?
                .factors()
                .iter()
                .any(|&(_, e)| e > 1)
            {
                continue;
            }
            for (h, _) in integers.factor_search(&norm, search)?.factors() {
                let h = integers.map_into(h, self, |c| Complex::from(c.clone()));
                let g = self.gcd_with(&shifted, &h, cancel)?;
                if variables(&g).is_empty() {
                    continue;
                }
                out.push(self.primitive_part(&moved(&g, &-&s)));
            }
            return Ok(out);
        }
        unreachable!("a shift with a square-free norm")
    }
}

impl SparseRing<GaussianRationals> {
    /// The factorisation of `p` into irreducible polynomials over ℚ(i),
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is the lexicographic leading coefficient,
    /// and the `fᵢ` are monic in the lexicographic order, in the order of
    /// the factors of the numerator over ℤ\[i\]. Zero gives `c = 0`.
    #[must_use]
    pub fn factor(&self, p: &Sparse<Complex>) -> Factors<Complex, Sparse<Complex>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Sparse<Complex>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Complex, Sparse<Complex>>, Cancelled> {
        let Some(lead) = p.coeffs().first() else {
            return Ok(Factors {
                constant: Complex::from(0),
                factors: Vec::new(),
            });
        };
        let numerators = SparseRing::new(GaussianIntegers, self.nvars());
        let (numerator, _) = split_sparse::<GaussianIntegers>(self, &numerators, p);
        let factors = numerators
            .factor_with(&numerator, cancel)?
            .into_parts()
            .1
            .into_iter()
            .map(|(g, e)| {
                let lc = g.coeffs().first().cloned().unwrap_or(Complex::from(1));
                let monic = numerators.map_into(&g, self, |c| c / &lc);
                (monic, e)
            })
            .collect();
        Ok(Factors {
            constant: lead.clone(),
            factors,
        })
    }
}
