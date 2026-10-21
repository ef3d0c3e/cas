//! Sparse multivariate polynomials (CAS-DOC-8 §11.2): terms in decreasing
//! lexicographic order of their exponents, the first variable most
//! significant, each a non-zero coefficient and its exponents packed into
//! machine words \[MP07\]. `doc/algorithms/sparse.md` describes the
//! algorithms.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities, k for a term"
)]

mod heap;
mod layout;

#[cfg(test)]
mod tests;

use std::borrow::Cow;
use std::cmp::Ordering;

use cas_num::Integer;
use cas_num::algebra::Ring;

use crate::cancel::{Cancel, Cancelled, Never};
use crate::dense::{Coefficients, Dense, DenseRing};
use layout::Layout;

/// A polynomial in the variables `x₀ … xₙ₋₁` of its ring: its terms in
/// decreasing lexicographic order of their exponent vectors, `x₀` most
/// significant, each coefficient non-zero. The exponents are packed in the
/// narrowest layout that holds the largest of them, so that equal
/// polynomials are equal values. The operations are those of a
/// [`SparseRing`].
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Sparse<E> {
    layout: Layout,
    coeffs: Vec<E>,
    /// `layout.words()` words a term.
    exps: Vec<u64>,
}

impl<E> Sparse<E> {
    /// The number of variables of its ring.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.layout.nvars()
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The number of terms.
    #[must_use]
    pub fn term_count(&self) -> usize {
        self.coeffs.len()
    }

    /// The coefficients, in the order of the terms.
    #[must_use]
    pub fn coeffs(&self) -> &[E] {
        &self.coeffs
    }

    /// The exponents of term `k`, from the leading term.
    #[must_use]
    pub fn exponents(&self, k: usize) -> Vec<u64> {
        self.layout.unpack(self.mono(k))
    }

    /// The terms, as exponent vectors and coefficients, from the leading
    /// term.
    pub fn terms(&self) -> impl Iterator<Item = (Vec<u64>, &E)> {
        self.coeffs
            .iter()
            .enumerate()
            .map(|(k, c)| (self.exponents(k), c))
    }

    /// The packed monomial of term `k`.
    fn mono(&self, k: usize) -> &[u64] {
        let w = self.layout.words();
        &self.exps[k * w..(k + 1) * w]
    }

    /// The largest exponent of any variable.
    fn max_exponent(&self) -> u64 {
        (0..self.coeffs.len())
            .map(|k| self.layout.max_in(self.mono(k)))
            .max()
            .unwrap_or(0)
    }

    /// The same polynomial in `layout`, wide enough.
    fn repacked(&self, layout: Layout) -> Sparse<E>
    where
        E: Clone,
    {
        if layout == self.layout {
            return self.clone();
        }
        let w = layout.words();
        let mut exps = vec![0; self.coeffs.len() * w];
        for k in 0..self.coeffs.len() {
            layout.pack(&self.exponents(k), &mut exps[k * w..(k + 1) * w]);
        }
        Sparse {
            layout,
            coeffs: self.coeffs.clone(),
            exps,
        }
    }

    /// In the narrowest layout that holds its exponents.
    fn canonical(self) -> Sparse<E>
    where
        E: Clone,
    {
        let layout = Layout::holding(self.nvars(), self.max_exponent());
        if layout == self.layout {
            self
        } else {
            self.repacked(layout)
        }
    }
}

/// The ring of sparse polynomials in `nvars` variables over the ring object
/// `R`: a [`Ring`] whose elements are [`Sparse`] polynomials. Division
/// ([`SparseRing::div_rem`], [`SparseRing::div_exact`]) needs a
/// [`EuclideanDomain`](cas_num::algebra::EuclideanDomain) of coefficients.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct SparseRing<R> {
    base: R,
    nvars: usize,
}

/// Two polynomials in one layout, and the layout.
type Common<'a, E> = (Cow<'a, Sparse<E>>, Cow<'a, Sparse<E>>, Layout);

/// `a` and `b` in one layout, the wider.
fn common<'a, E: Clone>(a: &'a Sparse<E>, b: &'a Sparse<E>) -> Common<'a, E> {
    match a.layout.bits().cmp(&b.layout.bits()) {
        Ordering::Equal => (Cow::Borrowed(a), Cow::Borrowed(b), a.layout),
        Ordering::Less => (Cow::Owned(a.repacked(b.layout)), Cow::Borrowed(b), b.layout),
        Ordering::Greater => (Cow::Borrowed(a), Cow::Owned(b.repacked(a.layout)), a.layout),
    }
}

/// Terms under construction: coefficients and packed monomials.
struct Terms<E> {
    layout: Layout,
    coeffs: Vec<E>,
    exps: Vec<u64>,
}

impl<E: Clone> Terms<E> {
    fn new(layout: Layout) -> Terms<E> {
        Terms {
            layout,
            coeffs: Vec::new(),
            exps: Vec::new(),
        }
    }

    fn push(&mut self, c: E, mono: &[u64]) {
        self.coeffs.push(c);
        self.exps.extend_from_slice(mono);
    }

    /// The polynomial of terms pushed in decreasing order.
    fn done(self) -> Sparse<E> {
        Sparse {
            layout: self.layout,
            coeffs: self.coeffs,
            exps: self.exps,
        }
        .canonical()
    }
}

impl<R: Ring> SparseRing<R> {
    #[must_use]
    pub fn new(base: R, nvars: usize) -> Self {
        SparseRing { base, nvars }
    }

    /// The ring of the coefficients.
    #[must_use]
    pub fn base(&self) -> &R {
        &self.base
    }

    #[must_use]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    fn empty(&self) -> Sparse<R::Elem> {
        Sparse {
            layout: Layout::holding(self.nvars, 0),
            coeffs: Vec::new(),
            exps: Vec::new(),
        }
    }

    /// The polynomial of `terms`, exponent vectors of `nvars` exponents
    /// with their coefficients, in any order: equal monomials are added,
    /// and zero coefficients dropped.
    ///
    /// # Panics
    ///
    /// If an exponent vector has another length, or an exponent is `2^63`
    /// or more.
    pub fn from_terms(
        &self,
        terms: impl IntoIterator<Item = (Vec<u64>, R::Elem)>,
    ) -> Sparse<R::Elem> {
        let mut terms: Vec<(Vec<u64>, R::Elem)> = terms.into_iter().collect();
        for (e, _) in &terms {
            assert_eq!(
                e.len(),
                self.nvars,
                "an exponent vector of {} variables",
                e.len()
            );
        }
        // Decreasing lexicographic order, equal monomials together.
        terms.sort_by(|(a, _), (b, _)| b.cmp(a));
        let max = terms
            .iter()
            .flat_map(|(e, _)| e.iter().copied())
            .max()
            .unwrap_or(0);
        let layout = Layout::holding(self.nvars, max);
        let mut out = Terms::new(layout);
        let mut mono = vec![0; layout.words()];
        let mut iter = terms.into_iter().peekable();
        while let Some((e, mut c)) = iter.next() {
            while let Some((_, d)) = iter.next_if(|(f, _)| *f == e) {
                c = self.base.add(&c, &d);
            }
            if !self.base.is_zero(&c) {
                layout.pack(&e, &mut mono);
                out.push(c, &mono);
            }
        }
        out.done()
    }

    /// The constant `c`.
    #[must_use]
    pub fn constant(&self, c: R::Elem) -> Sparse<R::Elem> {
        self.from_terms([(vec![0; self.nvars], c)])
    }

    /// `c x^exps`.
    #[must_use]
    pub fn monomial(&self, c: R::Elem, exps: &[u64]) -> Sparse<R::Elem> {
        self.from_terms([(exps.to_vec(), c)])
    }

    /// The variable `x_var`.
    ///
    /// # Panics
    ///
    /// If `var` is not a variable of the ring.
    #[must_use]
    pub fn var(&self, var: usize) -> Sparse<R::Elem> {
        assert!(var < self.nvars, "variable {var} of {}", self.nvars);
        let mut e = vec![0; self.nvars];
        e[var] = 1;
        self.monomial(self.base.one(), &e)
    }

    /// `a + s b` by merging, `s` the identity or negation.
    fn merge(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        f: impl Fn(&R::Elem) -> R::Elem,
    ) -> Sparse<R::Elem> {
        let (a, b, layout) = common(a, b);
        let mut out = Terms::new(layout);
        let (mut i, mut j) = (0, 0);
        while i < a.coeffs.len() || j < b.coeffs.len() {
            let order = match (i < a.coeffs.len(), j < b.coeffs.len()) {
                (true, true) => a.mono(i).cmp(b.mono(j)),
                (true, false) => Ordering::Greater,
                _ => Ordering::Less,
            };
            match order {
                Ordering::Greater => {
                    out.push(a.coeffs[i].clone(), a.mono(i));
                    i += 1;
                }
                Ordering::Less => {
                    out.push(f(&b.coeffs[j]), b.mono(j));
                    j += 1;
                }
                Ordering::Equal => {
                    let c = self.base.add(&a.coeffs[i], &f(&b.coeffs[j]));
                    if !self.base.is_zero(&c) {
                        out.push(c, a.mono(i));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.done()
    }

    /// `c p`.
    #[must_use]
    pub fn scale(&self, p: &Sparse<R::Elem>, c: &R::Elem) -> Sparse<R::Elem> {
        let mut out = Terms::new(p.layout);
        for (k, a) in p.coeffs.iter().enumerate() {
            let d = self.base.mul(c, a);
            if !self.base.is_zero(&d) {
                out.push(d, p.mono(k));
            }
        }
        out.done()
    }

    /// `a b` by Monagan and Pearce's heap \[MP07\], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn mul_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Sparse<R::Elem>, Cancelled> {
        heap::mul(self, a, b, cancel)
    }

    /// `p(point)`, the coefficient ring's value at a point of `nvars`
    /// values.
    ///
    /// # Panics
    ///
    /// If `point` has another length.
    #[must_use]
    pub fn eval(&self, p: &Sparse<R::Elem>, point: &[R::Elem]) -> R::Elem {
        assert_eq!(point.len(), self.nvars, "a point of {} values", point.len());
        let mut sum = self.base.zero();
        for (k, c) in p.coeffs.iter().enumerate() {
            let mut t = c.clone();
            for (var, x) in point.iter().enumerate() {
                let e = p.layout.get(p.mono(k), var);
                if e > 0 {
                    t = self.base.mul(&t, &self.base.pow(x, e));
                }
            }
            sum = self.base.add(&sum, &t);
        }
        sum
    }

    /// `p` with `x_var` replaced by the value `x`.
    #[must_use]
    pub fn eval_var(&self, p: &Sparse<R::Elem>, var: usize, x: &R::Elem) -> Sparse<R::Elem> {
        self.from_terms(p.terms().map(|(mut e, c)| {
            let power = self.base.pow(x, e[var]);
            e[var] = 0;
            (e, self.base.mul(c, &power))
        }))
    }

    /// `p` with `x_var` replaced by `q`: Horner's rule on the powers of
    /// `x_var`, the gaps by powers of `q`.
    #[must_use]
    pub fn substitute(
        &self,
        p: &Sparse<R::Elem>,
        var: usize,
        q: &Sparse<R::Elem>,
    ) -> Sparse<R::Elem> {
        let parts = self.to_recursive(p, var);
        let mut acc = self.zero();
        let mut last: Option<u64> = None;
        for (e, coeff) in parts {
            if let Some(l) = last {
                acc = self.mul(&acc, &self.pow(q, l - e));
            }
            acc = self.add(&acc, &coeff);
            last = Some(e);
        }
        match last {
            Some(l) if l > 0 => self.mul(&acc, &self.pow(q, l)),
            _ => acc,
        }
    }

    /// The derivative in `x_var`. Taking 1 from the exponent of `x_var` in
    /// every term keeps their order.
    #[must_use]
    pub fn derivative(&self, p: &Sparse<R::Elem>, var: usize) -> Sparse<R::Elem> {
        let mut out = Terms::new(p.layout);
        let mut mono = vec![0; p.layout.words()];
        for (k, c) in p.coeffs.iter().enumerate() {
            let e = p.layout.get(p.mono(k), var);
            if e == 0 {
                continue;
            }
            let d = self.base.mul(&self.base.elem_integer(&Integer::from(e)), c);
            if self.base.is_zero(&d) {
                continue;
            }
            mono.copy_from_slice(p.mono(k));
            p.layout.set(&mut mono, var, e - 1);
            out.push(d, &mono);
        }
        out.done()
    }

    /// The degree in `x_var`; `None` for the zero polynomial.
    #[must_use]
    pub fn degree(&self, p: &Sparse<R::Elem>, var: usize) -> Option<u64> {
        (0..p.coeffs.len())
            .map(|k| p.layout.get(p.mono(k), var))
            .max()
    }

    /// The coefficient of `x_var^d`, a polynomial without `x_var`. The terms
    /// with that power keep their order once it is cleared.
    #[must_use]
    pub fn coeff(&self, p: &Sparse<R::Elem>, var: usize, d: u64) -> Sparse<R::Elem> {
        let mut out = Terms::new(p.layout);
        let mut mono = vec![0; p.layout.words()];
        for (k, c) in p.coeffs.iter().enumerate() {
            if p.layout.get(p.mono(k), var) == d {
                mono.copy_from_slice(p.mono(k));
                p.layout.set(&mut mono, var, 0);
                out.push(c.clone(), &mono);
            }
        }
        out.done()
    }

    /// The leading term in lexicographic order: its exponents and
    /// coefficient.
    #[must_use]
    pub fn lead_term<'p>(&self, p: &'p Sparse<R::Elem>) -> Option<(Vec<u64>, &'p R::Elem)> {
        p.coeffs.first().map(|c| (p.exponents(0), c))
    }

    /// The leading coefficient in `x_var`: the coefficient of its highest
    /// power, zero for zero.
    #[must_use]
    pub fn lead_coeff(&self, p: &Sparse<R::Elem>, var: usize) -> Sparse<R::Elem> {
        match self.degree(p, var) {
            Some(d) => self.coeff(p, var, d),
            None => self.zero(),
        }
    }

    /// The recursive view: `p` as a polynomial in `x_var` whose
    /// coefficients are polynomials without it, its powers decreasing.
    #[must_use]
    pub fn to_recursive(&self, p: &Sparse<R::Elem>, var: usize) -> Vec<(u64, Sparse<R::Elem>)> {
        let mut powers: Vec<u64> = (0..p.coeffs.len())
            .map(|k| p.layout.get(p.mono(k), var))
            .collect();
        powers.sort_unstable_by(|a, b| b.cmp(a));
        powers.dedup();
        powers
            .into_iter()
            .map(|d| (d, self.coeff(p, var, d)))
            .collect()
    }

    /// The polynomial of a recursive view: `Σ cᵈ x_var^d`, each `cᵈ` without
    /// `x_var`, in any order.
    ///
    /// # Panics
    ///
    /// If a coefficient has `x_var`.
    #[must_use]
    pub fn from_recursive(&self, var: usize, parts: &[(u64, Sparse<R::Elem>)]) -> Sparse<R::Elem> {
        self.from_terms(parts.iter().flat_map(|(d, coeff)| {
            coeff.terms().map(move |(mut e, c)| {
                assert_eq!(e[var], 0, "a coefficient with the main variable");
                e[var] = *d;
                (e, c.clone())
            })
        }))
    }
}

impl<R: Coefficients + Clone> SparseRing<R> {
    /// `p` as a dense polynomial in `x_var`; `None` when another variable
    /// occurs.
    #[must_use]
    pub fn to_dense(&self, p: &Sparse<R::Elem>, var: usize) -> Option<Dense<R::Elem>> {
        let Some(d) = self.degree(p, var) else {
            return Some(Dense::default());
        };
        let mut coeffs = vec![self.base.zero(); usize::try_from(d).ok()? + 1];
        for (e, c) in p.terms() {
            if e.iter().enumerate().any(|(v, &x)| v != var && x != 0) {
                return None;
            }
            coeffs[usize::try_from(e[var]).ok()?] = c.clone();
        }
        Some(DenseRing::new(self.base.clone()).poly(coeffs))
    }

    /// The dense polynomial `d` in `x_var`.
    #[must_use]
    pub fn from_dense(&self, d: &Dense<R::Elem>, var: usize) -> Sparse<R::Elem> {
        self.from_terms(d.coeffs().iter().enumerate().map(|(k, c)| {
            let mut e = vec![0; self.nvars];
            e[var] = k as u64;
            (e, c.clone())
        }))
    }
}

/// Never fails: [`Never`] is never cancelled.
fn uncancelled<T>(result: Result<T, Cancelled>) -> T {
    match result {
        Ok(value) => value,
        Err(Cancelled) => unreachable!("never cancelled"),
    }
}

impl<R: Ring> Ring for SparseRing<R> {
    type Elem = Sparse<R::Elem>;

    fn zero(&self) -> Self::Elem {
        self.empty()
    }

    fn one(&self) -> Self::Elem {
        self.constant(self.base.one())
    }

    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.merge(a, b, Clone::clone)
    }

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.merge(a, b, |c| self.base.neg(c))
    }

    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        Sparse {
            layout: a.layout,
            coeffs: a.coeffs.iter().map(|c| self.base.neg(c)).collect(),
            exps: a.exps.clone(),
        }
    }

    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        uncancelled(self.mul_with(a, b, &Never))
    }

    fn elem_integer(&self, n: &Integer) -> Self::Elem {
        self.constant(self.base.elem_integer(n))
    }

    fn characteristic(&self) -> cas_num::Natural {
        self.base.characteristic()
    }

    fn is_zero(&self, a: &Self::Elem) -> bool {
        a.is_zero()
    }
}

/// Dense polynomials over sparse ones multiply by Karatsuba's method.
impl<R: Ring> Coefficients for SparseRing<R> {}
