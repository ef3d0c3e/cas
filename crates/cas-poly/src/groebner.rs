//! Gröbner bases (CAS-DOC-21 §2, §3, CAS-287): reduced Gröbner bases of ideals of
//! ℚ\[x₁, …, xₙ\] under the lexicographic, graded lexicographic and graded
//! reverse lexicographic orders, by Buchberger's algorithm with the
//! criteria in Gebauer and Möller's installation \[GM88\], and the
//! remainders of division by lists of polynomials, within the limits a
//! caller sets on their work (CAS-291).
//! `doc/algorithms/groebner.md` describes the algorithms.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of polynomial identities: f, g, h, a term's c and m"
)]

#[cfg(test)]
mod tests;

use std::cell::Cell;
use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use crate::cancel::{Cancel, Never};

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

/// A term: an exponent vector and a non-zero coefficient.
type Term = (Vec<u64>, Rational);

/// A polynomial of ℚ\[x₁, …, xₙ\] in distributed form: its terms in
/// decreasing order of their monomials under one [`MonomialOrder`], each an
/// exponent vector and a non-zero rational coefficient. The operations take
/// the order, which must be the one the terms were sorted by.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Distributed {
    terms: Vec<Term>,
}

impl Distributed {
    /// The polynomial of `terms`, in any order and with like monomials
    /// added, sorted under `order`. The exponent vectors have one length.
    #[must_use]
    pub fn new(mut terms: Vec<(Vec<u64>, Rational)>, order: MonomialOrder) -> Distributed {
        terms.sort_by(|(a, _), (b, _)| order.cmp(b, a));
        let mut out: Vec<Term> = Vec::with_capacity(terms.len());
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

    /// `self` over its leading coefficient, each coefficient made counted.
    fn monic(self, work: &Work<'_>) -> Result<Distributed, GroebnerError> {
        let Some((_, c)) = self.lead() else {
            return Ok(self);
        };
        if c.is_one() {
            return Ok(self);
        }
        let inverse = c.recip();
        let mut terms = Vec::with_capacity(self.terms.len());
        for (m, d) in self.terms {
            let scaled = &d * &inverse;
            work.made(&scaled)?;
            terms.push((m, scaled));
        }
        Ok(Distributed { terms })
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
}

/// The leading monomial of a polynomial known not to be 0.
fn leading(p: &Distributed) -> &[u64] {
    p.terms.first().map_or(&[], |(m, _)| m.as_slice())
}

/// The limits a caller sets on a computation (CAS-291): past either, it
/// stops with [`GroebnerError::TooLarge`]. A count of division steps
/// bounded no time: a step merges the whole remainder, thousands of terms
/// in a lexicographic basis, and the coefficients grow without bound.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct GroebnerLimits {
    /// The work: each coefficient a computation makes or copies, in its
    /// divisions and S-polynomials, counts 1 and 1 more for each 64 bits of
    /// its numerator and denominator; each critical pair it forms or goes
    /// through, and each term it moves to a remainder, counts 1.
    pub work: u64,
    /// The bits of a coefficient's numerator or denominator.
    pub coefficient_bits: u64,
}

impl GroebnerLimits {
    /// No limits: only the fixed bounds ([`MAX_ELEMENTS`],
    /// [`MAX_REDUCTIONS`], [`MAX_EXPONENT`]) and the [`Cancel`] stop a
    /// computation.
    pub const NONE: GroebnerLimits = GroebnerLimits {
        work: u64::MAX,
        coefficient_bits: u64::MAX,
    };
}

/// What a computation has left of its [`GroebnerLimits`], and the
/// [`Cancel`] it polls at each division step.
struct Work<'a> {
    cancel: &'a dyn Cancel,
    left: Cell<u64>,
    coefficient_bits: u64,
}

impl<'a> Work<'a> {
    fn new(limits: GroebnerLimits, cancel: &'a dyn Cancel) -> Work<'a> {
        Work {
            cancel,
            left: Cell::new(limits.work),
            coefficient_bits: limits.coefficient_bits,
        }
    }

    /// One more division step.
    fn step(&self) -> Result<(), GroebnerError> {
        if self.cancel.cancelled() {
            return Err(GroebnerError::Cancelled);
        }
        Ok(())
    }

    /// `units` more of the work.
    fn spend(&self, units: usize) -> Result<(), GroebnerError> {
        let units = u64::try_from(units).unwrap_or(u64::MAX);
        let left = self.left.get();
        if units > left {
            return Err(GroebnerError::TooLarge);
        }
        self.left.set(left - units);
        Ok(())
    }

    /// A coefficient made or copied, within the bound on its bits.
    fn made(&self, c: &Rational) -> Result<(), GroebnerError> {
        let (n, d) = (c.numer().bit_length(), c.denom().bit_length());
        if n.max(d) > self.coefficient_bits {
            return Err(GroebnerError::TooLarge);
        }
        self.spend(usize::try_from(1 + (n + d) / 64).unwrap_or(usize::MAX))
    }
}

/// `terms − c x^m g`, merged, `terms` sorted under `order`.
fn sub_mul(
    terms: &[Term],
    c: &Rational,
    m: &[u64],
    g: &Distributed,
    order: MonomialOrder,
    work: &Work<'_>,
) -> Result<Vec<Term>, GroebnerError> {
    let mut other = Vec::with_capacity(g.terms.len());
    for (t, d) in &g.terms {
        let cd = -(c * d);
        work.made(&cd)?;
        other.push((product(t, m), cd));
    }
    let mut out = Vec::with_capacity(terms.len() + other.len());
    let mut theirs = other.into_iter().peekable();
    for (a, x) in terms {
        // The greater terms of the multiple first, then a like one added.
        while let Some(term) = theirs.next_if(|(b, _)| order.cmp(a, b).is_lt()) {
            out.push(term);
        }
        if let Some((_, y)) = theirs.next_if(|(b, _)| b == a) {
            let sum = x + &y;
            work.made(&sum)?;
            if !sum.is_zero() {
                out.push((a.clone(), sum));
            }
        } else {
            work.made(x)?;
            out.push((a.clone(), x.clone()));
        }
    }
    out.extend(theirs);
    Ok(out)
}

/// [`s_polynomial`], each coefficient made counted.
fn s_poly(
    f: &Distributed,
    g: &Distributed,
    order: MonomialOrder,
    work: &Work<'_>,
) -> Result<Distributed, GroebnerError> {
    let ((fm, fc), (gm, gc)) = (f.lead().expect("non-zero"), g.lead().expect("non-zero"));
    let l = lcm(fm, gm);
    let shift = quotient(&l, fm);
    let mut scaled = Vec::with_capacity(f.terms.len());
    for (m, c) in &f.terms {
        let c = c / fc;
        work.made(&c)?;
        scaled.push((product(m, &shift), c));
    }
    let terms = sub_mul(&scaled, &gc.recip(), &quotient(&l, gm), g, order, work)?;
    Ok(Distributed { terms })
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
    s_poly(f, g, order, &Work::new(GroebnerLimits::NONE, &Never)).expect("no limits")
}

/// The remainder of `p` divided by `divisors` under `order` (CAS-DOC-21
/// §2): each term, leading or not, reduced by the first divisor whose
/// leading monomial divides it, until none of the remainder's terms is
/// divisible by a leading monomial. The zero divisors are left out.
///
/// # Errors
///
/// [`GroebnerError::Cancelled`] when `cancel` says so, and
/// [`GroebnerError::TooLarge`] past the `limits`.
pub fn groebner_reduce(
    p: &Distributed,
    divisors: &[Distributed],
    order: MonomialOrder,
    limits: GroebnerLimits,
    cancel: &dyn Cancel,
) -> Result<Distributed, GroebnerError> {
    let divisors: Vec<&Distributed> = divisors.iter().filter(|g| !g.is_zero()).collect();
    remainder(p, &divisors, order, &Work::new(limits, cancel))
}

/// [`groebner_reduce`] by non-zero divisors.
fn remainder(
    p: &Distributed,
    divisors: &[&Distributed],
    order: MonomialOrder,
    work: &Work<'_>,
) -> Result<Distributed, GroebnerError> {
    // What is left is `rest[done..]`: the terms before it went out, moved
    // rather than removed from the front one by one.
    let mut rest = p.terms.clone();
    let mut done = 0;
    let mut out = Vec::new();
    while let Some((m, c)) = rest.get(done) {
        work.step()?;
        let divisor = divisors
            .iter()
            .find(|g| g.leading_monomial().is_some_and(|gm| divides(gm, m)));
        if let Some(g) = divisor {
            let (gm, gc) = g.lead().expect("non-zero");
            let (q, shift) = (c / gc, quotient(m, gm));
            rest = sub_mul(&rest[done..], &q, &shift, g, order, work)?;
            done = 0;
        } else {
            work.spend(1)?;
            out.push(std::mem::replace(
                &mut rest[done],
                (Vec::new(), Rational::ZERO),
            ));
            done += 1;
        }
    }
    Ok(Distributed { terms: out })
}

/// Why [`groebner_basis`] or [`groebner_reduce`] stopped.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum GroebnerError {
    /// The [`Cancel`] said so.
    Cancelled,
    /// Past the bounds of a computation: [`MAX_ELEMENTS`] polynomials,
    /// [`MAX_REDUCTIONS`] S-polynomials, an exponent past [`MAX_EXPONENT`],
    /// or past its [`GroebnerLimits`].
    TooLarge,
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
    fn add(&mut self, h: Distributed, work: &Work<'_>) -> Result<(), GroebnerError> {
        if self.all.len() >= MAX_ELEMENTS || h.max_exponent() > MAX_EXPONENT {
            return Err(GroebnerError::TooLarge);
        }
        work.spend(self.basis.len() + self.pairs.len())?;
        self.all.push(h.monic(work)?);
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
/// [`GroebnerError::TooLarge`] past the bounds or the `limits`.
pub fn groebner_basis(
    polys: &[Distributed],
    order: MonomialOrder,
    limits: GroebnerLimits,
    cancel: &dyn Cancel,
) -> Result<Vec<Distributed>, GroebnerError> {
    let work = Work::new(limits, cancel);
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
        let h = remainder(p, &state.divisors(), order, &work)?;
        if h.is_constant() {
            return Ok(vec![one]);
        }
        if !h.is_zero() {
            state.add(h, &work)?;
        }
    }
    let mut reductions = 0;
    while !state.pairs.is_empty() {
        work.step()?;
        reductions += 1;
        if reductions > MAX_REDUCTIONS {
            return Err(GroebnerError::TooLarge);
        }
        work.spend(state.pairs.len())?;
        let k = (0..state.pairs.len())
            .min_by(|&a, &b| order.cmp(&state.pairs[a].lcm, &state.pairs[b].lcm))
            .unwrap_or(0);
        let pair = state.pairs.swap_remove(k);
        let s = s_poly(&state.all[pair.i], &state.all[pair.j], state.order, &work)?;
        let h = remainder(&s, &state.divisors(), order, &work)?;
        if h.is_constant() {
            return Ok(vec![one]);
        }
        if !h.is_zero() {
            state.add(h, &work)?;
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
            remainder(&basis[i], &others, order, &work)?
        };
        basis[i] = reduced.monic(&work)?;
    }
    Ok(basis)
}
