//! Zippel's sparse modular gcd \[Zip79, Zip90\] in de Kleine, Monagan and
//! Wittkopf's form, LINZIP \[DKMW05\] (CAS-152). Here, CAS-190: the image of
//! a gcd over a prime field from its assumed form, the monomials of an
//! earlier image, by univariate gcds in the main variable at the powers
//! `α, α², …` of a random point `α` for the other variables.
//!
//! The images are monic, and the gcd's at `αⁱ` is `λᵢ` times the `i`th,
//! `λᵢ` the value of its leading coefficient in the main variable there.
//! - **Single scaling**, when a coefficient of the main variable in the
//!   form is one monomial `μ`: scaled so that `μ`'s coefficient is 1, `λᵢ`
//!   is `μ(α)ⁱ` over that coefficient of the `i`th image, and each power of
//!   the main variable is a transposed Vandermonde system in the values of
//!   its monomials at `α`.
//! - **Multiple scaling** otherwise: the `λᵢ` are unknowns too, `λ₁ = 1`.
//!   Each power's first images give its coefficients in terms of the `λᵢ`,
//!   and its other images linear equations for them.
//!
//! One image more than the systems need checks them: a form with a term
//! missing makes them inconsistent.

#![cfg_attr(not(test), expect(dead_code, reason = "used by Zippel's gcd, CAS-191"))]

#[cfg(test)]
mod tests;

use std::cmp::Ordering;

use cas_num::Integer;
use cas_num::algebra::{EuclideanDomain, Field, Ring};

use super::vandermonde;
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Coefficients, DenseRing};
use crate::random::Random;
use crate::sparse::{Sparse, SparseRing};

/// The assumed form of a gcd: for each power of the main variable in it,
/// the monomials of its coefficient in the other variables.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct Form {
    /// The main variable: the images are univariate in it.
    main: usize,
    /// The other variables, evaluated at the powers of a point.
    vars: Vec<usize>,
    /// For each power of the main variable, from the largest down, the
    /// exponents of `vars` in the monomials of its coefficient.
    blocks: Vec<(u64, Vec<Vec<u64>>)>,
}

impl Form {
    /// The form of `g` in the main variable `main`, with `vars` the others
    /// that occur where its images are taken.
    pub(super) fn of<E>(g: &Sparse<E>, main: usize, vars: &[usize]) -> Form {
        let mut blocks: Vec<(u64, Vec<Vec<u64>>)> = Vec::new();
        for (e, _) in g.terms() {
            let m: Vec<u64> = vars.iter().map(|&v| e[v]).collect();
            match blocks.iter_mut().find(|(d, _)| *d == e[main]) {
                Some((_, monomials)) => monomials.push(m),
                None => blocks.push((e[main], vec![m])),
            }
        }
        blocks.sort_by_key(|(d, _)| std::cmp::Reverse(*d));
        Form {
            main,
            vars: vars.to_vec(),
            blocks,
        }
    }

    /// The block that is one monomial, for single scaling.
    fn single(&self) -> Option<usize> {
        self.blocks.iter().position(|(_, m)| m.len() == 1)
    }

    /// The images the systems need, and one more to check them: the
    /// largest block's terms; with multiple scaling, also enough equations
    /// for the scales, `N` with `Σₑ (N − tₑ) ≥ N − 1` over `s` blocks of
    /// `tₑ` terms, `N ≥ (Σ tₑ − 1) / (s − 1)` \[DKMW05\].
    fn images(&self) -> usize {
        let largest = self.blocks.iter().map(|(_, m)| m.len()).max().unwrap_or(0);
        let needed = if self.single().is_some() || self.blocks.len() < 2 {
            largest
        } else {
            let terms: usize = self.blocks.iter().map(|(_, m)| m.len()).sum();
            largest.max((terms - 1).div_ceil(self.blocks.len() - 1))
        };
        needed + 1
    }
}

/// Why a form gave no image.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Failure {
    /// The form is wrong: an image of lower degree in the main variable,
    /// a power it lacks, a missing scaling term, inconsistent systems.
    Form,
    /// The scales stay under-determined: the gcd's image has a content in
    /// the main variable, an unlucky one or its own \[DKMW05\].
    Content,
    /// Every point tried was bad or unlucky, or made two monomials' values
    /// equal: an unlucky evaluation of the variables fixed before, or a
    /// field too small.
    Points,
}

/// Points tried before the images count as unlucky.
const ATTEMPTS: usize = 4;

/// Images added to under-determined scales before they count as a content.
const MORE_IMAGES: usize = 2;

/// A random non-zero element.
fn nonzero<F: Ring>(field: &F, random: &mut Random) -> F::Elem {
    loop {
        let x = field.elem_integer(&Integer::from(random.next()));
        if !field.is_zero(&x) {
            return x;
        }
    }
}

/// The value at `alpha` of the monomial of exponents `m`.
fn monomial<F: Ring>(field: &F, alpha: &[F::Elem], m: &[u64]) -> F::Elem {
    alpha
        .iter()
        .zip(m)
        .fold(field.one(), |w, (x, &e)| field.mul(&w, &field.pow(x, e)))
}

/// A polynomial's terms ready for its values at the powers `αⁱ` of a point:
/// for each term, its power of the main variable, its value at the last
/// power, and its monomial's value at `α`, which takes it to the next.
struct Evaluation<E> {
    terms: Vec<(usize, E, E)>,
    degree: usize,
}

impl<E: Clone> Evaluation<E> {
    fn new<F: Ring<Elem = E>>(field: &F, p: &Sparse<E>, form: &Form, alpha: &[E]) -> Self {
        let mut degree = 0;
        let terms = p
            .terms()
            .map(|(e, c)| {
                debug_assert!(
                    e.iter()
                        .enumerate()
                        .all(|(v, &x)| x == 0 || v == form.main || form.vars.contains(&v)),
                    "a variable outside the form"
                );
                let k = usize::try_from(e[form.main]).expect("a degree in memory");
                degree = degree.max(k);
                let m: Vec<u64> = form.vars.iter().map(|&v| e[v]).collect();
                (k, c.clone(), monomial(field, alpha, &m))
            })
            .collect();
        Evaluation { terms, degree }
    }

    /// The coefficients of the polynomial at the next power of the point,
    /// univariate in the main variable.
    fn next<F: Ring<Elem = E>>(&mut self, field: &F) -> Vec<E> {
        let mut coeffs = vec![field.zero(); self.degree + 1];
        for (k, value, w) in &mut self.terms {
            *value = field.mul(value, w);
            coeffs[*k] = field.add(&coeffs[*k], value);
        }
        coeffs
    }
}

/// The outcome of the systems.
enum Solved<E> {
    /// The coefficients of each block.
    Coefficients(Vec<Vec<E>>),
    Inconsistent,
    Underdetermined,
    /// Two monomials of a block have one value, or one is 0.
    Singular,
}

/// Whether `Σⱼ cⱼ vⱼⁱ = rhs` for the powers `i = from + 1, …` of the
/// values `rhs`.
fn consistent<F: Field>(
    field: &F,
    v: &[F::Elem],
    c: &[F::Elem],
    rhs: &[F::Elem],
    from: usize,
) -> bool {
    let from = u64::try_from(from).expect("a count in memory");
    let mut terms: Vec<F::Elem> = c
        .iter()
        .zip(v)
        .map(|(cj, vj)| field.mul(cj, &field.pow(vj, from)))
        .collect();
    rhs.iter().all(|r| {
        let mut sum = field.zero();
        for (t, vj) in terms.iter_mut().zip(v) {
            *t = field.mul(t, vj);
            sum = field.add(&sum, t);
        }
        sum == *r
    })
}

/// Single scaling by block `s`, one monomial of value `μ = values[s][0]`:
/// `λᵢ = μⁱ / ys[i][s]`, and each block's first images give its
/// coefficients, its others check them.
fn single_scaling<F: Field>(
    field: &F,
    values: &[Vec<F::Elem>],
    ys: &[Vec<F::Elem>],
    s: usize,
) -> Solved<F::Elem> {
    let mu = &values[s][0];
    let mut power = field.one();
    let mut lambda = Vec::with_capacity(ys.len());
    for y in ys {
        power = field.mul(&power, mu);
        // The scaling term missing from an image: the form is wrong.
        let Some(inverse) = field.inv(&y[s]) else {
            return Solved::Inconsistent;
        };
        lambda.push(field.mul(&power, &inverse));
    }
    let mut coefficients = Vec::with_capacity(values.len());
    for (e, v) in values.iter().enumerate() {
        let rhs: Vec<F::Elem> = lambda
            .iter()
            .zip(ys)
            .map(|(l, y)| field.mul(l, &y[e]))
            .collect();
        let t = v.len();
        let Some(c) = vandermonde::solve(field, v, &rhs[..t]) else {
            return Solved::Singular;
        };
        if !consistent(field, v, &c, &rhs[t..], t) {
            return Solved::Inconsistent;
        }
        coefficients.push(c);
    }
    Solved::Coefficients(coefficients)
}

/// The scales `λ₁ = 1, λ₂, …, λₙ` from rows `Σₖ rₖ λₖ = 0`, by Gaussian
/// elimination; the failure when they are inconsistent or under-determined.
fn scales<F: Field>(
    field: &F,
    rows: Vec<Vec<F::Elem>>,
    n: usize,
) -> Result<Vec<F::Elem>, Solved<F::Elem>> {
    let unknowns = n - 1;
    // λ₁ = 1 goes to the right-hand side.
    let mut m: Vec<Vec<F::Elem>> = rows
        .into_iter()
        .map(|r| {
            let mut x = r[1..].to_vec();
            x.push(field.neg(&r[0]));
            x
        })
        .collect();
    let mut rank = 0;
    for col in 0..unknowns {
        let Some(p) = (rank..m.len()).find(|&i| !field.is_zero(&m[i][col])) else {
            continue;
        };
        m.swap(rank, p);
        let inverse = field.inv(&m[rank][col]).expect("a non-zero pivot");
        for x in &mut m[rank] {
            *x = field.mul(x, &inverse);
        }
        let pivot = m[rank].clone();
        for (i, row) in m.iter_mut().enumerate() {
            if i != rank && !field.is_zero(&row[col]) {
                let f = row[col].clone();
                for (x, y) in row.iter_mut().zip(&pivot) {
                    *x = field.sub(x, &field.mul(&f, y));
                }
            }
        }
        rank += 1;
    }
    if m[rank..].iter().any(|row| !field.is_zero(&row[unknowns])) {
        return Err(Solved::Inconsistent);
    }
    if rank < unknowns {
        return Err(Solved::Underdetermined);
    }
    // Full rank: the pivots are the columns in order.
    let mut lambda = vec![field.one()];
    lambda.extend(m.iter().take(unknowns).map(|row| row[unknowns].clone()));
    Ok(lambda)
}

/// Multiple scaling \[DKMW05\]: with `M` the inverse of block `e`'s first `t`
/// images' system, its coefficients are `cⱼ = Σᵢ Mⱼᵢ yᵢₑ λᵢ`; each further
/// image `k` adds the equation `Σⱼ vⱼᵏ cⱼ = yₖₑ λₖ` for the scales.
fn multiple_scaling<F: Field>(
    field: &F,
    values: &[Vec<F::Elem>],
    ys: &[Vec<F::Elem>],
) -> Solved<F::Elem> {
    let n = ys.len();
    let mut inverses = Vec::with_capacity(values.len());
    let mut rows = Vec::new();
    for (e, v) in values.iter().enumerate() {
        let t = v.len();
        let Some(m) = vandermonde::inverse(field, v) else {
            return Solved::Singular;
        };
        let start = u64::try_from(t).expect("a count in memory");
        let mut powers: Vec<F::Elem> = v.iter().map(|x| field.pow(x, start)).collect();
        for k in t..n {
            // The powers k + 1 of the values, for image k (from 0).
            for (p, x) in powers.iter_mut().zip(v) {
                *p = field.mul(p, x);
            }
            let mut row = vec![field.zero(); n];
            for (i, r) in row.iter_mut().take(t).enumerate() {
                let sum = powers.iter().zip(&m).fold(field.zero(), |s, (p, mj)| {
                    field.add(&s, &field.mul(p, &mj[i]))
                });
                *r = field.mul(&sum, &ys[i][e]);
            }
            row[k] = field.neg(&ys[k][e]);
            rows.push(row);
        }
        inverses.push(m);
    }
    let lambda = match scales(field, rows, n) {
        Ok(lambda) => lambda,
        Err(failure) => return failure,
    };
    let coefficients = inverses
        .iter()
        .enumerate()
        .map(|(e, m)| {
            m.iter()
                .map(|mj| {
                    mj.iter().enumerate().fold(field.zero(), |s, (i, x)| {
                        field.add(&s, &field.mul(&field.mul(x, &ys[i][e]), &lambda[i]))
                    })
                })
                .collect()
        })
        .collect();
    Solved::Coefficients(coefficients)
}

/// The polynomial of the form with the solved coefficients.
fn assemble<F: Ring>(
    ring: &SparseRing<F>,
    form: &Form,
    coefficients: Vec<Vec<F::Elem>>,
) -> Sparse<F::Elem> {
    let n = ring.nvars();
    ring.from_terms(
        form.blocks
            .iter()
            .zip(coefficients)
            .flat_map(|((d, monomials), cs)| {
                monomials.iter().zip(cs).map(move |(m, c)| {
                    let mut e = vec![0; n];
                    e[form.main] = *d;
                    for (&v, &x) in form.vars.iter().zip(m) {
                        e[v] = x;
                    }
                    (e, c)
                })
            }),
    )
}

/// The gcd of `a` and `b` over a prime field, up to a scalar, from its
/// assumed form \[Zip90, DKMW05\]. `a` and `b` have no variables but the
/// form's main variable and its others; `random` draws the points. A point
/// where `a` or `b` loses degree in the main variable is bad, and one whose
/// image has a larger degree than the form unlucky: another point is
/// tried.
pub(super) fn image<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    form: &Form,
    random: &mut Random,
    cancel: &dyn Cancel,
) -> Result<Result<Sparse<F::Elem>, Failure>, Cancelled> {
    let field = ring.base();
    let dense = DenseRing::new(field.clone());
    let degree = |p: &Sparse<F::Elem>| {
        ring.degree(p, form.main)
            .map(|d| usize::try_from(d).expect("a degree in memory"))
    };
    let (da, db) = (degree(a), degree(b));
    let Some(&(top, _)) = form.blocks.first() else {
        return Ok(Err(Failure::Form));
    };
    let top = usize::try_from(top).expect("a degree in memory");
    let single = form.single();
    if single.is_none() && form.blocks.len() < 2 {
        return Ok(Err(Failure::Content));
    }
    'points: for _ in 0..ATTEMPTS {
        let alpha: Vec<F::Elem> = form.vars.iter().map(|_| nonzero(field, random)).collect();
        let values: Vec<Vec<F::Elem>> = form
            .blocks
            .iter()
            .map(|(_, ms)| ms.iter().map(|m| monomial(field, &alpha, m)).collect())
            .collect();
        let (mut ea, mut eb) = (
            Evaluation::new(field, a, form, &alpha),
            Evaluation::new(field, b, form, &alpha),
        );
        let mut count = form.images();
        let mut ys: Vec<Vec<F::Elem>> = Vec::with_capacity(count);
        let mut more = 0;
        loop {
            while ys.len() < count {
                if cancel.cancelled() {
                    return Err(Cancelled);
                }
                let (ai, bi) = (dense.poly(ea.next(field)), dense.poly(eb.next(field)));
                if ai.degree() != da || bi.degree() != db {
                    continue 'points;
                }
                let g = EuclideanDomain::gcd(&dense, &ai, &bi);
                match g.degree().cmp(&Some(top)) {
                    Ordering::Greater => continue 'points,
                    Ordering::Less => return Ok(Err(Failure::Form)),
                    Ordering::Equal => {}
                }
                // A power of the main variable outside the form.
                let outside = g.coeffs().iter().enumerate().any(|(k, c)| {
                    !field.is_zero(c) && !form.blocks.iter().any(|(d, _)| *d == k as u64)
                });
                if outside {
                    return Ok(Err(Failure::Form));
                }
                ys.push(
                    form.blocks
                        .iter()
                        .map(|(d, _)| g.coeffs()[usize::try_from(*d).expect("at most top")].clone())
                        .collect(),
                );
            }
            let solved = match single {
                Some(s) => single_scaling(field, &values, &ys, s),
                None => multiple_scaling(field, &values, &ys),
            };
            match solved {
                Solved::Coefficients(c) => return Ok(Ok(assemble(ring, form, c))),
                Solved::Inconsistent => return Ok(Err(Failure::Form)),
                Solved::Singular => continue 'points,
                Solved::Underdetermined if more < MORE_IMAGES => {
                    more += 1;
                    count += 1;
                }
                Solved::Underdetermined => return Ok(Err(Failure::Content)),
            }
        }
    }
    Ok(Err(Failure::Points))
}
