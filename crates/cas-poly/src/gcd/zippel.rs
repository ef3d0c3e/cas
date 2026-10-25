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

#[cfg(test)]
mod tests;

use std::cmp::Ordering;

use cas_num::Integer;
use cas_num::algebra::{EuclideanDomain, Field, Integers, Ring};

use super::modular::{
    Primes, Runs, Step, combine, content, divide_runs, from_runs, is_constant, lift, runs,
};
use super::{Cofactors, GcdDomain, normalized_sparse, vandermonde};
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Coefficients, Dense, DenseRing};
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
    /// `tₑ` terms, `N ≥ (Σ tₑ − 1) / (s − 1)` \[DKMW05\]. Without other
    /// variables, the one image there is.
    fn images(&self) -> usize {
        if self.vars.is_empty() {
            return 1;
        }
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
    /// For `p` in the variable `main` and those of `vars`, at the powers of
    /// the point `x_vars[k] = alpha[k]`.
    fn new<F: Ring<Elem = E>>(
        field: &F,
        p: &Sparse<E>,
        main: usize,
        vars: &[usize],
        alpha: &[E],
    ) -> Self {
        let mut degree = 0;
        let terms = p
            .terms()
            .map(|(e, c)| {
                debug_assert!(
                    e.iter()
                        .enumerate()
                        .all(|(v, &x)| x == 0 || v == main || vars.contains(&v)),
                    "a variable outside the form"
                );
                let k = usize::try_from(e[main]).expect("a degree in memory");
                degree = degree.max(k);
                let m: Vec<u64> = vars.iter().map(|&v| e[v]).collect();
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
            Evaluation::new(field, a, form.main, &form.vars, &alpha),
            Evaluation::new(field, b, form.main, &form.vars, &alpha),
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

/// Restarts from a new first image before Zippel's gcd gives up.
const RESTARTS: usize = 4;

/// The seed of the points of Zippel's gcd, so that its results repeat.
const SEED: u64 = 0x5A49_5050_454C;

/// `p` at the point `x_vars[k] = point[k]`, univariate in `main`.
fn univariate<F: Field + Coefficients + Clone>(
    dense: &DenseRing<F>,
    p: &Sparse<F::Elem>,
    main: usize,
    vars: &[usize],
    point: &[F::Elem],
) -> Dense<F::Elem> {
    let field = dense.base();
    dense.poly(Evaluation::new(field, p, main, vars, point).next(field))
}

/// Whether `g` divides `a` and `b` at a random point of `vars`, univariate in
/// `main`, where `g` keeps its degree \[DKMW05\]: a test that a wrong `g`
/// fails with high probability, much cheaper than a division.
fn divides_at_a_point<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    [g, a, b]: [&Sparse<F::Elem>; 3],
    main: usize,
    vars: &[usize],
    random: &mut Random,
) -> bool {
    let dense = DenseRing::new(ring.base().clone());
    let degree = ring
        .degree(g, main)
        .map(|d| usize::try_from(d).expect("a degree in memory"));
    for _ in 0..ATTEMPTS {
        let point: Vec<F::Elem> = vars.iter().map(|_| nonzero(ring.base(), random)).collect();
        let at = |p| univariate(&dense, p, main, vars, &point);
        let g1 = at(g);
        if g1.degree() == degree {
            return dense.div_rem(&at(a), &g1).1.is_zero()
                && dense.div_rem(&at(b), &g1).1.is_zero();
        }
    }
    false
}

/// The monic gcd of non-zero `a` and `b` over a prime field by Zippel's
/// algorithm in LINZIP's form \[Zip79, DKMW05\] (CAS-191); `None` when it
/// gives up, after `RESTARTS` wrong forms or failed tests, or in a field with
/// too few points. The first variable that occurs is the main one.
///
/// In one variable, Euclid's algorithm. Otherwise `a` and `b` are taken in
/// `F[y][X]`, `y` the last variable that occurs, as Brown's algorithm takes
/// them: their contents in `F[y]` come out, and `γ`, the gcd of their leading
/// coefficients in `X`, is the leading coefficient each image is scaled to.
/// The first image, at a random point of `y`, comes from this algorithm in
/// one variable fewer and gives the form of the others, which the sparse
/// interpolation computes. Newton's interpolation in `y` combines them until
/// it stops changing or passes the degree bound; its primitive part in `y`
/// is then tested by division at a random point of the other variables.
pub(super) fn gcd<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    random: &mut Random,
    cancel: &dyn Cancel,
) -> Result<Option<Sparse<F::Elem>>, Cancelled> {
    let base = ring.base();
    let (da, db) = (a.degrees(), b.degrees());
    let occurring: Vec<usize> = (0..ring.nvars())
        .filter(|&v| da[v] > 0 || db[v] > 0)
        .collect();
    let (Some(&main), Some(&y)) = (occurring.first(), occurring.last()) else {
        return Ok(Some(ring.one()));
    };
    let dense = DenseRing::new(base.clone());
    if main == y {
        let to = |p: &Sparse<F::Elem>| ring.to_dense(p, y).expect("one variable");
        let g = EuclideanDomain::gcd(&dense, &to(a), &to(b));
        return Ok(Some(ring.from_dense(&g, y)));
    }
    let (ra, rb) = (runs(&dense, a, y), runs(&dense, b, y));
    let (ca, cb) = (content(&dense, &ra), content(&dense, &rb));
    let c = ring.from_dense(&EuclideanDomain::gcd(&dense, &ca, &cb), y);
    let (ra, rb) = (divide_runs(&dense, ra, &ca), divide_runs(&dense, rb, &cb));
    let gamma = EuclideanDomain::gcd(&dense, &ra[0].1, &rb[0].1);
    let (pa, pb) = (from_runs(ring, &ra, y), from_runs(ring, &rb, y));
    // The gcd scaled to the leading coefficient γ has at most this degree
    // in y.
    let degree = |r: &Runs<F::Elem>| r.iter().filter_map(|(_, q)| q.degree()).max().unwrap_or(0);
    let bound = gamma.degree().unwrap_or(0) + degree(&ra).min(degree(&rb));
    let others = &occurring[1..occurring.len() - 1];
    // A point of y where neither γ nor `m`, the product of the y − v so
    // far, vanishes, and γ there.
    let point = |random: &mut Random, m: &Dense<F::Elem>| {
        (0..ATTEMPTS).find_map(|_| {
            let v = nonzero(base, random);
            let at = dense.eval(&gamma, &v);
            (!base.is_zero(&at) && !base.is_zero(&dense.eval(m, &v))).then_some((v, at))
        })
    };
    let linear = |v: &F::Elem| dense.poly(vec![base.neg(v), base.one()]);
    for _ in 0..RESTARTS {
        let Some((v, at)) = point(random, &dense.one()) else {
            return Ok(None);
        };
        let (av, bv) = (ring.eval_var(&pa, y, &v), ring.eval_var(&pb, y, &v));
        let Some(first) = gcd(ring, &av, &bv, random, cancel)? else {
            return Ok(None);
        };
        if is_constant(&first) {
            return Ok(Some(normalized_sparse(ring, &c)));
        }
        let form = Form::of(&first, main, others);
        let lead = first.exponents(0);
        let mut h = ring.scale(&first, &at);
        let mut m = linear(&v);
        let mut count = 1;
        let mut unlucky = 0;
        let interpolated = loop {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            let Some((v, at)) = point(random, &m) else {
                break None;
            };
            let (av, bv) = (ring.eval_var(&pa, y, &v), ring.eval_var(&pb, y, &v));
            let g = match image(ring, &av, &bv, &form, random, cancel)? {
                Ok(g) => g,
                // An unlucky point of y, likely: another.
                Err(Failure::Points) if unlucky < ATTEMPTS => {
                    unlucky += 1;
                    continue;
                }
                Err(_) => break None,
            };
            // Its leading term is the form's, γ(v) its coefficient.
            if g.is_zero() || g.exponents(0) != lead {
                break None;
            }
            let scale = base.mul(&at, &base.inv(&g.coeffs()[0]).expect("non-zero"));
            let g = ring.scale(&g, &scale);
            // Newton's step: h + (g − h(v)) m / m(v).
            let diff = ring.sub(&g, &ring.eval_var(&h, y, &v));
            let unchanged = diff.is_zero();
            if !unchanged {
                let scale = base.inv(&dense.eval(&m, &v)).expect("a new point");
                h = ring.add(
                    &h,
                    &ring.mul(&ring.scale(&diff, &scale), &ring.from_dense(&m, y)),
                );
            }
            m = dense.mul(&m, &linear(&v));
            count += 1;
            if unchanged || count > bound {
                break Some(h);
            }
        };
        let Some(h) = interpolated else {
            continue;
        };
        let r = runs(&dense, &h, y);
        let candidate = from_runs(
            ring,
            &divide_runs(&dense, r.clone(), &content(&dense, &r)),
            y,
        );
        if divides_at_a_point(ring, [&candidate, &pa, &pb], main, &occurring[1..], random) {
            return Ok(Some(normalized_sparse(ring, &ring.mul(&c, &candidate))));
        }
    }
    Ok(None)
}

/// `a` and `b` as `cₐ pₐ` and `c_b p_b`, `cₐ` and `c_b` their contents in
/// the main variable, by gcds in fewer variables, and `c` the gcd of those.
struct Contents<E> {
    ca: Sparse<E>,
    cb: Sparse<E>,
    c: Sparse<E>,
    pa: Sparse<E>,
    pb: Sparse<E>,
}

impl<E> Contents<E> {
    fn of<R: GcdDomain<Elem = E>>(
        ring: &SparseRing<R>,
        [a, b]: [&Sparse<E>; 2],
        main: usize,
        cancel: &dyn Cancel,
    ) -> Result<Contents<E>, Cancelled> {
        let content_in = |p: &Sparse<E>| -> Result<Sparse<E>, Cancelled> {
            let mut g = ring.zero();
            for (_, c) in ring.to_recursive(p, main) {
                g = ring.gcd_with(&g, &c, cancel)?;
                if ring.is_one(&g) {
                    break;
                }
            }
            Ok(g)
        };
        let (ca, cb) = (content_in(a)?, content_in(b)?);
        let c = ring.gcd_with(&ca, &cb, cancel)?;
        let pa = exact(ring, a, &ca, cancel)?;
        let pb = exact(ring, b, &cb, cancel)?;
        Ok(Contents { ca, cb, c, pa, pb })
    }
}

/// `p / q`, which is exact.
fn exact<R: GcdDomain>(
    ring: &SparseRing<R>,
    p: &Sparse<R::Elem>,
    q: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    Ok(ring.div_exact_with(p, q, cancel)?.expect("a divisor"))
}

/// The first variable of `a` or `b`, the main one.
fn main_variable<E>(ring_vars: usize, a: &Sparse<E>, b: &Sparse<E>) -> Option<usize> {
    let (da, db) = (a.degrees(), b.degrees());
    (0..ring_vars).find(|&v| da[v] > 0 || db[v] > 0)
}

/// The gcd of non-zero `a` and `b` over a prime field by Zippel's algorithm,
/// monic, their contents in the main variable taken out first, as LINZIP
/// takes them \[DKMW05\]: a content of the gcd would leave the scales of
/// multiple scaling free. `None` when the algorithm gives up.
pub(super) fn field_gcd<F: GcdDomain + Field>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Sparse<F::Elem>>, Cancelled> {
    let Some(main) = main_variable(ring.nvars(), a, b) else {
        return Ok(Some(ring.one()));
    };
    let parts = Contents::of(ring, [a, b], main, cancel)?;
    let mut random = Random::new(SEED);
    Ok(gcd(ring, &parts.pa, &parts.pb, &mut random, cancel)?
        .map(|g| normalized_sparse(ring, &ring.mul(&parts.c, &g))))
}

/// Forms taken from a new prime before Zippel's gcd over ℤ gives up.
const FORMS: usize = 4;

/// The gcd of non-zero `a` and `b` over ℤ by Zippel's algorithm in
/// LINZIP's form \[DKMW05\] (CAS-192), with its cofactors; `None` when it
/// gives up, after `FORMS` forms or more primes than the coefficients
/// should need, and Brown's algorithm takes over.
///
/// The contents in the main variable come out first, by gcds in fewer
/// variables. With `γ` the gcd of the leading coefficients of the primitive
/// parts, the first prime's image, by Zippel's algorithm over `𝔽ₚ`, gives
/// the form of the gcd. The next primes' images come from the sparse
/// interpolation alone, in all the other variables at once. Each image is
/// scaled to `γ` modulo its prime, and they are combined by Chinese
/// remaindering as in Brown's algorithm until the combination stops
/// changing; its primitive part is then tested by division. A wrong form
/// is taken again from the next prime.
pub(super) fn integers(
    ring: &SparseRing<Integers>,
    a: &Sparse<Integer>,
    b: &Sparse<Integer>,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Sparse<Integer>>>, Cancelled> {
    let Some(main) = main_variable(ring.nvars(), a, b) else {
        let (x, y) = (&a.coeffs()[0], &b.coeffs()[0]);
        let g = x.gcd(y);
        let constant = |c: Integer| ring.constant(c);
        return Ok(Some((
            constant(g.clone()),
            constant(x.div_exact(&g)),
            constant(y.div_exact(&g)),
        )));
    };
    let Contents { ca, cb, c, pa, pb } = Contents::of(ring, [a, b], main, cancel)?;
    // a = c (cₐ / c) g qₐ, b = c (c_b / c) g q_b.
    let done = |g: &Sparse<Integer>, qa: &Sparse<Integer>, qb: &Sparse<Integer>| {
        Ok::<_, Cancelled>(Some((
            ring.mul(&c, g),
            ring.mul(&exact(ring, &ca, &c, cancel)?, qa),
            ring.mul(&exact(ring, &cb, &c, cancel)?, qb),
        )))
    };
    if is_constant(&pa) || is_constant(&pb) {
        return done(&ring.one(), &pa, &pb);
    }
    let (da, db) = (pa.degrees(), pb.degrees());
    let others: Vec<usize> = (0..ring.nvars())
        .filter(|&v| v != main && (da[v] > 0 || db[v] > 0))
        .collect();
    let gamma = pa.coeffs()[0].gcd(&pb.coeffs()[0]);
    // More primes than this, and an undetected wrong image is likelier than
    // coefficients so large.
    let bits = |p: &Sparse<Integer>| {
        p.coeffs()
            .iter()
            .map(Integer::bit_length)
            .max()
            .unwrap_or(0)
    };
    let total = |d: &[u64]| d.iter().sum::<u64>();
    let most = 64 + (bits(&pa).max(bits(&pb)) + total(&da).max(total(&db))) / 30;
    let mut random = Random::new(SEED);
    let mut form: Option<Form> = None;
    let mut forms = 0;
    let mut h: Option<(Sparse<Integer>, Integer)> = None;
    for (p, field) in Primes::new().take(usize::try_from(most).unwrap_or(usize::MAX)) {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let gp = field.elem_integer(&gamma);
        if gp == 0 {
            continue;
        }
        let fp = SparseRing::new(field.clone(), ring.nvars());
        let image_of = |q: &Sparse<Integer>| ring.map_into(q, &fp, |x| field.elem_integer(x));
        let (ap, bp) = (image_of(&pa), image_of(&pb));
        let g = if let Some(f) = &form {
            match image(&fp, &ap, &bp, f, &mut random, cancel)? {
                Ok(g) => g,
                // An unlucky prime, likely.
                Err(Failure::Points) => continue,
                // A wrong form: the next prime gives another.
                Err(_) => {
                    form = None;
                    continue;
                }
            }
        } else {
            if forms == FORMS {
                return Ok(None);
            }
            forms += 1;
            let Some(g) = gcd(&fp, &ap, &bp, &mut random, cancel)? else {
                return Ok(None);
            };
            if is_constant(&g) {
                return done(&ring.one(), &pa, &pb);
            }
            form = Some(Form::of(&g, main, &others));
            g
        };
        let Some(inverse) = field.inv(&g.coeffs()[0]) else {
            continue;
        };
        let g = fp.scale(&g, &field.mul(&gp, &inverse));
        let lead = g.exponents(0);
        h = match h {
            Some((old, m)) if old.exponents(0) == lead => {
                let step = Step::new(&field, &m, p);
                let new = combine(ring, &old, &g, &step);
                if new == old {
                    let content = ring.content(&new);
                    let candidate = ring.map_into(&new, ring, |x| x.div_exact(&content));
                    let candidate = normalized_sparse(ring, &candidate);
                    if let Some(qa) = ring.div_exact_with(&pa, &candidate, cancel)?
                        && let Some(qb) = ring.div_exact_with(&pb, &candidate, cancel)?
                    {
                        return done(&candidate, &qa, &qb);
                    }
                }
                Some((new, step.mp))
            }
            // An unlucky prime: its gcd is larger.
            Some((old, m)) if old.exponents(0) < lead => Some((old, m)),
            // The first prime, or the earlier ones were unlucky.
            _ => Some((fp.map_into(&g, ring, |&r| lift(r, p)), Integer::from(p))),
        };
    }
    Ok(None)
}
