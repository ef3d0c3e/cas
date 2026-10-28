//! Gröbner bases (CAS-DOC-21 §2, §3, CAS-287): reduced Gröbner bases of ideals of
//! ℚ\[x₁, …, xₙ\] under the lexicographic, graded lexicographic and graded
//! reverse lexicographic orders, by Buchberger's algorithm with the
//! criteria in Gebauer and Möller's installation \[GM88\], and the
//! remainders of division by lists of polynomials.
//! `doc/algorithms/groebner.md` describes the algorithms.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities: f, g, h, a term's c and m"
)]

#[cfg(test)]
mod tests;

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use crate::cancel::{Cancel, Cancelled};

/// A monomial order, `x₁ > x₂ > … > xₙ`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum MonomialOrder {
    /// The exponents compared from `x₁`: giac's `plex`.
    Lex,
    /// The total degrees, then [`Lex`](Self::Lex): giac's `tdeg`.
    DegLex,
    /// The total degrees, then the last differing exponent, the smaller
    /// one greater: giac's `revlex`, its default.
    DegRevLex,
}

impl MonomialOrder {
    /// How the monomials of the exponent vectors `a` and `b`, of one
    /// length, compare.
    #[must_use]
    pub fn cmp(self, a: &[u64], b: &[u64]) -> Ordering {
        match self {
            MonomialOrder::Lex => a.cmp(b),
            MonomialOrder::DegLex => degree(a).cmp(&degree(b)).then_with(|| a.cmp(b)),
            MonomialOrder::DegRevLex => degree(a).cmp(&degree(b)).then_with(|| {
                a.iter()
                    .zip(b)
                    .rev()
                    .find(|(x, y)| x != y)
                    .map_or(Ordering::Equal, |(x, y)| y.cmp(x))
            }),
        }
    }
}

/// The total degree of a monomial.
fn degree(a: &[u64]) -> u128 {
    a.iter().map(|&e| u128::from(e)).sum()
}

/// The product of two monomials, saturated: an exponent past
/// [`MAX_EXPONENT`] stops a basis computation before any overflows.
fn product(a: &[u64], b: &[u64]) -> Vec<u64> {
    a.iter().zip(b).map(|(x, y)| x.saturating_add(*y)).collect()
}

/// Whether the monomial `a` divides `b`.
fn divides(a: &[u64], b: &[u64]) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}

/// `b / a`, for `a` dividing `b`.
fn quotient(b: &[u64], a: &[u64]) -> Vec<u64> {
    b.iter().zip(a).map(|(y, x)| y - x).collect()
}

fn lcm(a: &[u64], b: &[u64]) -> Vec<u64> {
    a.iter().zip(b).map(|(x, y)| *x.max(y)).collect()
}

/// Whether two monomials have no variable in common.
fn disjoint(a: &[u64], b: &[u64]) -> bool {
    a.iter().zip(b).all(|(x, y)| *x == 0 || *y == 0)
}

/// A polynomial of ℚ\[x₁, …, xₙ\] in distributed form: its terms in
/// decreasing order of their monomials under one [`MonomialOrder`], each an
/// exponent vector and a non-zero rational coefficient. The operations take
/// the order, which must be the one the terms were sorted by.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Distributed {
    terms: Vec<(Vec<u64>, Rational)>,
}

impl Distributed {
    /// The polynomial of `terms`, in any order and with like monomials
    /// added, sorted under `order`. The exponent vectors have one length.
    #[must_use]
    pub fn new(mut terms: Vec<(Vec<u64>, Rational)>, order: MonomialOrder) -> Distributed {
        terms.sort_by(|(a, _), (b, _)| order.cmp(b, a));
        let mut out: Vec<(Vec<u64>, Rational)> = Vec::with_capacity(terms.len());
        for (m, c) in terms {
            match out.last_mut() {
                Some((last, sum)) if *last == m => *sum = &*sum + &c,
                _ => out.push((m, c)),
            }
        }
        out.retain(|(_, c)| !c.is_zero());
        Distributed { terms: out }
    }

    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The terms, from the leading one.
    #[must_use]
    pub fn terms(&self) -> &[(Vec<u64>, Rational)] {
        &self.terms
    }

    /// The leading monomial; `None` for 0.
    #[must_use]
    pub fn leading_monomial(&self) -> Option<&[u64]> {
        self.terms.first().map(|(m, _)| m.as_slice())
    }

    fn lead(&self) -> Option<(&[u64], &Rational)> {
        self.terms.first().map(|(m, c)| (m.as_slice(), c))
    }

    /// Whether it is a constant other than 0.
    fn is_constant(&self) -> bool {
        matches!(self.terms.as_slice(), [(m, _)] if m.iter().all(|&e| e == 0))
    }

    /// The largest exponent.
    fn max_exponent(&self) -> u64 {
        self.terms
            .iter()
            .flat_map(|(m, _)| m.iter().copied())
            .max()
            .unwrap_or(0)
    }

    /// `self` over its leading coefficient.
    fn monic(self) -> Distributed {
        let Some((_, c)) = self.lead() else {
            return self;
        };
        if c.is_one() {
            return self;
        }
        let inverse = c.recip();
        let terms = self
            .terms
            .into_iter()
            .map(|(m, d)| (m, &d * &inverse))
            .collect();
        Distributed { terms }
    }

    /// `self` times the rational that makes its coefficients coprime
    /// integers, the leading one positive.
    #[must_use]
    pub fn primitive(&self) -> Distributed {
        let Some((_, lead)) = self.lead() else {
            return self.clone();
        };
        let denominators = self
            .terms
            .iter()
            .fold(Integer::ONE, |l, (_, c)| l.lcm(c.denom()));
        let up = Rational::from(denominators);
        let scaled: Vec<Rational> = self.terms.iter().map(|(_, c)| c * &up).collect();
        let content = scaled.iter().fold(Integer::ZERO, |g, c| g.gcd(c.numer()));
        let mut factor = Rational::new(Integer::ONE, content);
        if lead.is_negative() {
            factor = -factor;
        }
        let terms = self
            .terms
            .iter()
            .zip(&scaled)
            .map(|((m, _), c)| (m.clone(), c * &factor))
            .collect();
        Distributed { terms }
    }

    /// `self` times the monomial `m`, which keeps the order.
    fn shifted(&self, m: &[u64]) -> Distributed {
        let terms = self
            .terms
            .iter()
            .map(|(t, c)| (product(t, m), c.clone()))
            .collect();
        Distributed { terms }
    }

    /// `self − c x^m g`, merging.
    fn sub_mul(
        &self,
        c: &Rational,
        m: &[u64],
        g: &Distributed,
        order: MonomialOrder,
    ) -> Distributed {
        let other: Vec<(Vec<u64>, Rational)> = g
            .terms
            .iter()
            .map(|(t, d)| (product(t, m), -(c * d)))
            .collect();
        let mut out = Vec::with_capacity(self.terms.len() + other.len());
        let (mut i, mut j) = (0, 0);
        while i < self.terms.len() && j < other.len() {
            let ((a, x), (b, y)) = (&self.terms[i], &other[j]);
            match order.cmp(a, b) {
                Ordering::Greater => {
                    out.push((a.clone(), x.clone()));
                    i += 1;
                }
                Ordering::Less => {
                    out.push((b.clone(), y.clone()));
                    j += 1;
                }
                Ordering::Equal => {
                    let sum = x + y;
                    if !sum.is_zero() {
                        out.push((a.clone(), sum));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&self.terms[i..]);
        out.extend_from_slice(&other[j..]);
        Distributed { terms: out }
    }
}

/// The leading monomial of a polynomial known not to be 0.
fn leading(p: &Distributed) -> &[u64] {
    p.terms.first().map_or(&[], |(m, _)| m.as_slice())
}

/// The S-polynomial of `f` and `g`, both non-zero: each times the monomial
/// that lifts its leading monomial to their lcm, over its leading
/// coefficient, the second subtracted.
///
/// # Panics
///
/// If `f` or `g` is 0.
#[must_use]
pub fn s_polynomial(f: &Distributed, g: &Distributed, order: MonomialOrder) -> Distributed {
    let ((fm, fc), (gm, gc)) = (f.lead().expect("non-zero"), g.lead().expect("non-zero"));
    let l = lcm(fm, gm);
    let lifted = f.shifted(&quotient(&l, fm));
    let scaled = Distributed {
        terms: lifted
            .terms
            .into_iter()
            .map(|(m, c)| (m, &c / fc))
            .collect(),
    };
    scaled.sub_mul(&gc.recip(), &quotient(&l, gm), g, order)
}

/// The remainder of `p` divided by `divisors` under `order` (CAS-DOC-21
/// §2): each term, leading or not, reduced by the first divisor whose
/// leading monomial divides it, until none of the remainder's terms is
/// divisible by a leading monomial. The zero divisors are left out.
///
/// # Errors
///
/// [`Cancelled`] when `cancel` says so.
pub fn groebner_reduce(
    p: &Distributed,
    divisors: &[Distributed],
    order: MonomialOrder,
    cancel: &dyn Cancel,
) -> Result<Distributed, Cancelled> {
    let divisors: Vec<&Distributed> = divisors.iter().filter(|g| !g.is_zero()).collect();
    remainder(p, &divisors, order, cancel)
}

/// [`groebner_reduce`] by non-zero divisors.
fn remainder(
    p: &Distributed,
    divisors: &[&Distributed],
    order: MonomialOrder,
    cancel: &dyn Cancel,
) -> Result<Distributed, Cancelled> {
    let mut rest = p.clone();
    let mut out = Vec::new();
    while let Some((m, c)) = rest.lead() {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let divisor = divisors
            .iter()
            .find(|g| g.leading_monomial().is_some_and(|gm| divides(gm, m)));
        match divisor {
            Some(g) => {
                let (gm, gc) = g.lead().expect("non-zero");
                let (q, shift) = (c / gc, quotient(m, gm));
                rest = rest.sub_mul(&q, &shift, g, order);
            }
            None => out.push(rest.terms.remove(0)),
        }
    }
    Ok(Distributed { terms: out })
}

/// Why [`groebner_basis`] stopped.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum GroebnerError {
    /// The [`Cancel`] said so.
    Cancelled,
    /// Past the bounds of a computation: [`MAX_ELEMENTS`] polynomials,
    /// [`MAX_REDUCTIONS`] S-polynomials, or an exponent past
    /// [`MAX_EXPONENT`].
    TooLarge,
}

impl From<Cancelled> for GroebnerError {
    fn from(_: Cancelled) -> Self {
        GroebnerError::Cancelled
    }
}

/// The most polynomials a basis computation adds to its basis.
pub const MAX_ELEMENTS: usize = 1 << 12;

/// The most S-polynomials a basis computation reduces.
pub const MAX_REDUCTIONS: usize = 1 << 16;

/// The largest exponent a basis computation holds.
pub const MAX_EXPONENT: u64 = 1 << 20;

/// A critical pair: two elements of the basis and the lcm of their leading
/// monomials.
struct Pair {
    i: usize,
    j: usize,
    lcm: Vec<u64>,
}

/// The state of Buchberger's algorithm: every polynomial added, the ones
/// still in the basis, and the pairs left.
struct Buchberger {
    order: MonomialOrder,
    all: Vec<Distributed>,
    basis: Vec<usize>,
    pairs: Vec<Pair>,
}

impl Buchberger {
    fn lead(&self, k: usize) -> &[u64] {
        leading(&self.all[k])
    }

    /// The basis polynomials, the divisors of a reduction.
    fn divisors(&self) -> Vec<&Distributed> {
        self.basis.iter().map(|&k| &self.all[k]).collect()
    }

    /// Adds `h`, reduced by the basis and not 0 nor a constant, made monic,
    /// with its pairs: Gebauer and Möller's update \[GM88\], as Becker and
    /// Weispfenning give it \[BW93, UPDATE\].
    fn add(&mut self, h: Distributed) -> Result<(), GroebnerError> {
        if self.all.len() >= MAX_ELEMENTS || h.max_exponent() > MAX_EXPONENT {
            return Err(GroebnerError::TooLarge);
        }
        self.all.push(h.monic());
        let h = self.all.len() - 1;
        let lh = self.lead(h).to_vec();
        // The new pairs, those whose lcm another's divides dropped, one of
        // equal lcms kept; then the coprime ones (Buchberger's first
        // criterion), which served to drop others.
        let mut new: Vec<Pair> = self
            .basis
            .iter()
            .map(|&k| Pair {
                i: k,
                j: h,
                lcm: lcm(self.lead(k), &lh),
            })
            .collect();
        let mut kept: Vec<Pair> = Vec::with_capacity(new.len());
        while let Some(p) = new.pop() {
            let coprime = disjoint(self.lead(p.i), &lh);
            if coprime || !new.iter().chain(&kept).any(|q| divides(&q.lcm, &p.lcm)) {
                kept.push(p);
            }
        }
        kept.retain(|p| !disjoint(self.lead(p.i), &lh));
        // The old pairs whose lcm the new leading monomial divides, but
        // equals neither lcm with it.
        let pairs = std::mem::take(&mut self.pairs);
        self.pairs = pairs
            .into_iter()
            .filter(|p| {
                !divides(&lh, &p.lcm)
                    || lcm(self.lead(p.i), &lh) == p.lcm
                    || lcm(self.lead(p.j), &lh) == p.lcm
            })
            .collect();
        self.pairs.extend(kept);
        // The basis elements whose leading monomial the new one divides go.
        let basis = std::mem::take(&mut self.basis);
        self.basis = basis
            .into_iter()
            .filter(|&k| !divides(&lh, self.lead(k)))
            .collect();
        self.basis.push(h);
        Ok(())
    }
}

/// The reduced Gröbner basis of the ideal `polys` generate under `order`
/// (CAS-DOC-21 §3): each element monic, in increasing order of the leading
/// monomials; `[1]` for the whole ring, `[]` for the zero ideal. The
/// polynomials have exponent vectors of one length, sorted under `order`.
///
/// Buchberger's algorithm: the pair of the least lcm first (the normal
/// strategy), its S-polynomial reduced by the basis, added when not 0;
/// pairs dropped by Gebauer and Möller's criteria \[GM88\]. Then each
/// element reduced by the others.
///
/// # Errors
///
/// [`GroebnerError::Cancelled`] when `cancel` says so, and
/// [`GroebnerError::TooLarge`] past the bounds.
pub fn groebner_basis(
    polys: &[Distributed],
    order: MonomialOrder,
    cancel: &dyn Cancel,
) -> Result<Vec<Distributed>, GroebnerError> {
    let mut inputs: Vec<&Distributed> = polys.iter().filter(|p| !p.is_zero()).collect();
    if inputs.iter().any(|p| p.max_exponent() > MAX_EXPONENT) {
        return Err(GroebnerError::TooLarge);
    }
    let Some(first) = inputs.first() else {
        return Ok(Vec::new());
    };
    let one = Distributed::new(vec![(vec![0; leading(first).len()], Rational::ONE)], order);
    // The least leading monomials first: fewer elements replaced.
    inputs.sort_by(|a, b| order.cmp(leading(a), leading(b)));
    let mut state = Buchberger {
        order,
        all: Vec::new(),
        basis: Vec::new(),
        pairs: Vec::new(),
    };
    for p in inputs {
        let h = remainder(p, &state.divisors(), order, cancel)?;
        if h.is_constant() {
            return Ok(vec![one]);
        }
        if !h.is_zero() {
            state.add(h)?;
        }
    }
    let mut reductions = 0;
    while !state.pairs.is_empty() {
        if cancel.cancelled() {
            return Err(GroebnerError::Cancelled);
        }
        reductions += 1;
        if reductions > MAX_REDUCTIONS {
            return Err(GroebnerError::TooLarge);
        }
        let k = (0..state.pairs.len())
            .min_by(|&a, &b| order.cmp(&state.pairs[a].lcm, &state.pairs[b].lcm))
            .unwrap_or(0);
        let pair = state.pairs.swap_remove(k);
        let s = s_polynomial(&state.all[pair.i], &state.all[pair.j], state.order);
        let h = remainder(&s, &state.divisors(), order, cancel)?;
        if h.is_constant() {
            return Ok(vec![one]);
        }
        if !h.is_zero() {
            state.add(h)?;
        }
    }
    // The minimal basis, each element reduced by the others: its leading
    // monomial stays, so the order of the reductions does not matter.
    let mut basis: Vec<Distributed> = state.basis.iter().map(|&k| state.all[k].clone()).collect();
    basis.sort_by(|a, b| order.cmp(leading(a), leading(b)));
    for i in 0..basis.len() {
        let reduced = {
            let others: Vec<&Distributed> = basis
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, g)| g)
                .collect();
            remainder(&basis[i], &others, order, cancel)?
        };
        basis[i] = reduced.monic();
    }
    Ok(basis)
}
