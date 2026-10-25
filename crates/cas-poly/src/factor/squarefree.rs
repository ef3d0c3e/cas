//! Square-free decomposition by Yun's algorithm \[Yun76\] (CAS-179).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Yun's paper: q = a b, c and d for derivatives"
)]

use cas_num::algebra::{EuclideanDomain, Ring};

use super::{CharacteristicZero, Factors};
use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::dense::{Dense, DenseRing};
use crate::gcd::{Cofactors, normalized_sparse};
use crate::sparse::{Sparse, SparseRing};

/// Yun's algorithm on `q`, primitive in the variable of `derivative`, of
/// positive degree in it. It pushes the square-free parts and their exponents
/// to `out`, leaving out parts without the variable, which are units.
///
/// With `q = ∏ fᵢ^i`, the `fᵢ` square-free and pairwise coprime, the gcd of
/// `q` and `q'` is `∏ fᵢ^(i-1)`, with cofactors `b₁ = ∏ fᵢ` and `c₁`. Round
/// `k` has `bₖ = ∏_(i≥k) fᵢ` and `cₖ - bₖ' = fₖ Σ_(i>k) (i - k) fᵢ'
/// ∏_(j>k, j≠i) fⱼ`, so `fₖ` is their gcd, with cofactors `bₖ₊₁` and `cₖ₊₁`.
/// All of this holds up to units.
fn yun<P: Ring>(
    ring: &P,
    q: &P::Elem,
    derivative: impl Fn(&P::Elem) -> P::Elem,
    is_constant: impl Fn(&P::Elem) -> bool,
    gcd: impl Fn(&P::Elem, &P::Elem) -> Result<Cofactors<P::Elem>, Cancelled>,
    out: &mut Vec<(P::Elem, u64)>,
) -> Result<(), Cancelled> {
    let (_, mut b, mut c) = gcd(q, &derivative(q))?;
    let mut k = 1;
    while !is_constant(&b) {
        let d = ring.sub(&c, &derivative(&b));
        let (f, b_next, c_next) = gcd(&b, &d)?;
        if !is_constant(&f) {
            out.push((f, k));
        }
        (b, c) = (b_next, c_next);
        k += 1;
    }
    Ok(())
}

/// The constant of `p = c ∏ fᵢ^eᵢ`, from the leading coefficients: `lead(p)
/// / ∏ lead(fᵢ)^eᵢ`, an exact division. The factors are normalised, but over
/// ℤ\[i\] a product of normalised leading coefficients may not be.
fn constant_of<R: EuclideanDomain>(
    base: &R,
    lead: &R::Elem,
    leads: impl Iterator<Item = (R::Elem, u64)>,
) -> R::Elem {
    let product = leads.fold(base.one(), |acc, (l, e)| base.mul(&acc, &base.pow(&l, e)));
    let (c, r) = base.div_rem(lead, &product);
    debug_assert!(base.is_zero(&r), "the factors' leading coefficients divide");
    c
}

impl<R: CharacteristicZero> DenseRing<R> {
    /// The square-free decomposition `p = c f₁ f₂² ⋯ fₖᵏ`. The `fᵢ` are
    /// square-free and pairwise coprime. Each is primitive, with its leading
    /// coefficient in normal form, as [`Factors`] describes. Parts of degree
    /// 0 are left out, so the exponents increase. giac:
    /// `sqrfree(12x⁴ − 24x² + 12)` is `12 (x² − 1)²`.
    #[must_use]
    pub fn squarefree(&self, p: &Dense<R::Elem>) -> Factors<R::Elem, Dense<R::Elem>> {
        uncancelled(self.squarefree_with(p, &Never))
    }

    /// As [`DenseRing::squarefree`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn squarefree_with(
        &self,
        p: &Dense<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<R::Elem, Dense<R::Elem>>, Cancelled> {
        let base = self.base();
        let Some(lead) = p.lead() else {
            return Ok(Factors {
                constant: base.zero(),
                factors: Vec::new(),
            });
        };
        let mut factors = Vec::new();
        let q = self.primitive_part(p);
        if q.degree() > Some(0) {
            yun(
                self,
                &q,
                |f| self.derivative(f),
                |f| f.degree() == Some(0),
                |a, b| self.gcd_cofactors_with(a, b, cancel),
                &mut factors,
            )?;
        }
        let leads = factors
            .iter()
            .filter_map(|(f, e)| Some((f.lead()?.clone(), *e)));
        Ok(Factors {
            constant: constant_of(base, lead, leads),
            factors,
        })
    }
}

impl<R: CharacteristicZero> SparseRing<R> {
    /// The square-free decomposition `p = c f₁ f₂² ⋯ fₖᵏ`, as
    /// [`DenseRing::squarefree`]'s. The leading coefficients are the
    /// lexicographic ones. The content in the first variable that occurs is
    /// decomposed in the others, and the primitive part by Yun's algorithm in
    /// that variable. Parts with one exponent are then multiplied together.
    #[must_use]
    pub fn squarefree(&self, p: &Sparse<R::Elem>) -> Factors<R::Elem, Sparse<R::Elem>> {
        uncancelled(self.squarefree_with(p, &Never))
    }

    /// As [`SparseRing::squarefree`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn squarefree_with(
        &self,
        p: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<R::Elem, Sparse<R::Elem>>, Cancelled> {
        let base = self.base();
        let Some(lead) = p.coeffs().first() else {
            return Ok(Factors {
                constant: base.zero(),
                factors: Vec::new(),
            });
        };
        let mut parts = Vec::new();
        self.squarefree_parts(&self.primitive_part(p), &mut parts, cancel)?;
        parts.sort_by_key(|&(_, e)| e);
        let mut factors: Vec<(Sparse<R::Elem>, u64)> = Vec::with_capacity(parts.len());
        for (f, e) in parts {
            match factors.last_mut() {
                Some((g, last)) if *last == e => {
                    *g = normalized_sparse(self, &self.mul_with(g, &f, cancel)?);
                }
                _ => factors.push((f, e)),
            }
        }
        let leads = factors.iter().map(|(f, e)| (f.coeffs()[0].clone(), *e));
        Ok(Factors {
            constant: constant_of(base, lead, leads),
            factors,
        })
    }

    /// The square-free parts of `q`, primitive with its leading coefficient
    /// in normal form, with their exponents, pushed to `out`. Those of its
    /// content in the first variable it has come first, by recursion in the
    /// remaining variables, then Yun's in that variable. A constant has
    /// none.
    fn squarefree_parts(
        &self,
        q: &Sparse<R::Elem>,
        out: &mut Vec<(Sparse<R::Elem>, u64)>,
        cancel: &dyn Cancel,
    ) -> Result<(), Cancelled> {
        let Some(var) = q.degrees().iter().position(|&d| d > 0) else {
            return Ok(());
        };
        let mut content = self.zero();
        for (_, c) in self.to_recursive(q, var) {
            content = self.gcd_with(&content, &c, cancel)?;
            if self.is_one(&content) {
                break;
            }
        }
        let primitive = if self.is_one(&content) {
            q.clone()
        } else {
            self.squarefree_parts(&content, out, cancel)?;
            self.div_exact_with(q, &content, cancel)?
                .expect("the content divides")
        };
        let constant = |f: &Sparse<R::Elem>| self.degree(f, var).unwrap_or(0) == 0;
        yun(
            self,
            &primitive,
            |f| self.derivative(f, var),
            constant,
            |a, b| self.gcd_cofactors_with(a, b, cancel),
            out,
        )
    }
}
