//! Brown's modular gcd \[Bro71\]. Over ℤ, the gcds of the images modulo
//! word primes, each scaled so that its leading coefficient is the image
//! of `γ`, the gcd of the leading coefficients, combined by Chinese
//! remaindering until the result stops changing and its primitive part
//! divides both. Over a prime field, in several variables, the same with
//! the images at points of the last variable, combined by interpolation.
//! An image whose gcd has a larger leading monomial than the others comes
//! from an unlucky prime or point and is dropped; a smaller one shows the
//! earlier ones were unlucky, and the combination starts over.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Brown's paper: a, b, g, h, m, p, γ"
)]

use cas_num::Integer;
use cas_num::algebra::{EuclideanDomain, Field, Integers, Ring, WordPrimeField};

use super::subres::div_coeffs;
use super::{Cofactors, normalized_sparse};
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Coefficients, Dense, DenseRing};
use crate::sparse::{Sparse, SparseRing};

/// Word primes, decreasing from `2^62`, each with its field.
pub(crate) struct Primes {
    next: u64,
}

impl Primes {
    pub(crate) fn new() -> Primes {
        Primes {
            next: (1 << 62) - 1,
        }
    }
}

impl Iterator for Primes {
    type Item = (u64, WordPrimeField);

    fn next(&mut self) -> Option<(u64, WordPrimeField)> {
        while self.next > 2 {
            let p = self.next;
            self.next -= 2;
            if let Some(field) = WordPrimeField::new(p) {
                return Some((p, field));
            }
        }
        None
    }
}

/// `r` modulo `p` as the integer in `(-p/2, p/2)`.
fn lift(r: u64, p: u64) -> Integer {
    if r > p / 2 {
        Integer::from(r) - Integer::from(p)
    } else {
        Integer::from(r)
    }
}

/// One step of Chinese remaindering: from `x ≡ h (mod m)`, `h` in
/// `(-m/2, m/2]`, and `x ≡ g (mod p)`, the `x` in `(-mp/2, mp/2]`, as
/// `h + m ((g - h) m⁻¹ mod p)` brought into range.
struct Step<'a> {
    field: &'a WordPrimeField,
    m: &'a Integer,
    mp: Integer,
    inverse: u64,
}

impl<'a> Step<'a> {
    fn new(field: &'a WordPrimeField, m: &'a Integer, p: u64) -> Step<'a> {
        let inverse = field
            .inv(&field.elem_integer(m))
            .expect("a modulus of other primes");
        Step {
            field,
            m,
            mp: m * &Integer::from(p),
            inverse,
        }
    }

    fn combine(&self, h: &Integer, g: u64) -> Integer {
        let f = self.field;
        let t = f.mul(&f.sub(&g, &f.elem_integer(h)), &self.inverse);
        let x = h + &(self.m * &Integer::from(t));
        if &x + &x > self.mp { x - &self.mp } else { x }
    }
}

/// The gcd `c` of the contents, and `a / c`, `b / c`.
fn contents(ca: &Integer, cb: &Integer) -> (Integer, Integer, Integer) {
    let c = ca.gcd(cb);
    (c.clone(), ca.div_exact(&c), cb.div_exact(&c))
}

/// The gcd of non-zero dense `a` and `b` over ℤ, with its cofactors.
pub(crate) fn dense(
    ring: &DenseRing<Integers>,
    a: &Dense<Integer>,
    b: &Dense<Integer>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Dense<Integer>>, Cancelled> {
    dense_from(ring, a, b, Primes::new(), cancel)
}

/// [`dense`] over the primes of `primes`, which must not run out.
pub(crate) fn dense_from(
    ring: &DenseRing<Integers>,
    a: &Dense<Integer>,
    b: &Dense<Integer>,
    primes: impl Iterator<Item = (u64, WordPrimeField)>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Dense<Integer>>, Cancelled> {
    let (ca, cb) = (ring.content(a), ring.content(b));
    let (pa, pb) = (div_coeffs(ring, a, &ca), div_coeffs(ring, b, &cb));
    let (c, ka, kb) = contents(&ca, &cb);
    let done = |g: &Dense<Integer>, qa: &Dense<Integer>, qb: &Dense<Integer>| {
        (ring.scale(g, &c), ring.scale(qa, &ka), ring.scale(qb, &kb))
    };
    if pa.degree() == Some(0) || pb.degree() == Some(0) {
        return Ok(done(&ring.one(), &pa, &pb));
    }
    let gamma = pa
        .lead()
        .expect("non-zero")
        .gcd(pb.lead().expect("non-zero"));
    // The combination so far, and the product of its primes.
    let mut h: Option<(Dense<Integer>, Integer)> = None;
    for (p, field) in primes {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let gp = field.elem_integer(&gamma);
        if gp == 0 {
            continue;
        }
        let fp = DenseRing::new(field.clone());
        let image = |q: &Dense<Integer>| {
            fp.poly(q.coeffs().iter().map(|x| field.elem_integer(x)).collect())
        };
        let g = EuclideanDomain::gcd(&fp, &image(&pa), &image(&pb));
        if g.degree() == Some(0) {
            return Ok(done(&ring.one(), &pa, &pb));
        }
        let g = fp.scale(&g, &gp);
        h = match h {
            Some((old, m)) if old.degree() == g.degree() => {
                let step = Step::new(&field, &m, p);
                let new = ring.poly(
                    old.coeffs()
                        .iter()
                        .zip(g.coeffs())
                        .map(|(x, &y)| step.combine(x, y))
                        .collect(),
                );
                if new == old {
                    let candidate = div_coeffs(ring, &new, &ring.content(&new));
                    if let Some(qa) = ring.div_exact(&pa, &candidate)
                        && let Some(qb) = ring.div_exact(&pb, &candidate)
                    {
                        return Ok(done(&candidate, &qa, &qb));
                    }
                }
                Some((new, step.mp))
            }
            // An unlucky prime: its gcd is larger.
            Some((old, m)) if old.degree() < g.degree() => Some((old, m)),
            // The first prime, or the earlier ones were unlucky.
            _ => Some((
                ring.poly(g.coeffs().iter().map(|&r| lift(r, p)).collect()),
                Integer::from(p),
            )),
        };
    }
    unreachable!("enough word primes")
}

/// Whether `p` is a constant.
pub(super) fn is_constant<E>(p: &Sparse<E>) -> bool {
    p.degrees().iter().all(|&d| d == 0)
}

/// The gcd of non-zero sparse `a` and `b` over ℤ, with its cofactors: the
/// images modulo word primes by [`pgcd`].
pub(crate) fn sparse(
    ring: &SparseRing<Integers>,
    a: &Sparse<Integer>,
    b: &Sparse<Integer>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Sparse<Integer>>, Cancelled> {
    sparse_from(ring, a, b, Primes::new(), cancel)
}

/// [`sparse`] over the primes of `primes`, which must not run out.
pub(crate) fn sparse_from(
    ring: &SparseRing<Integers>,
    a: &Sparse<Integer>,
    b: &Sparse<Integer>,
    primes: impl Iterator<Item = (u64, WordPrimeField)>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Sparse<Integer>>, Cancelled> {
    let (ca, cb) = (ring.content(a), ring.content(b));
    let divide = |p: &Sparse<Integer>, c: &Integer| ring.map_into(p, ring, |x| x.div_exact(c));
    let (pa, pb) = (divide(a, &ca), divide(b, &cb));
    let (c, ka, kb) = contents(&ca, &cb);
    let done = |g: &Sparse<Integer>, qa: &Sparse<Integer>, qb: &Sparse<Integer>| {
        (ring.scale(g, &c), ring.scale(qa, &ka), ring.scale(qb, &kb))
    };
    if is_constant(&pa) || is_constant(&pb) {
        return Ok(done(&ring.one(), &pa, &pb));
    }
    let gamma = pa.coeffs()[0].gcd(&pb.coeffs()[0]);
    let mut h: Option<(Sparse<Integer>, Integer)> = None;
    for (p, field) in primes {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let gp = field.elem_integer(&gamma);
        if gp == 0 {
            continue;
        }
        let fp = SparseRing::new(field.clone(), ring.nvars());
        let image = |q: &Sparse<Integer>| ring.map_into(q, &fp, |x| field.elem_integer(x));
        let Some(g) = pgcd(&fp, &image(&pa), &image(&pb), cancel)? else {
            continue;
        };
        if is_constant(&g) {
            return Ok(done(&ring.one(), &pa, &pb));
        }
        let g = fp.scale(&g, &gp);
        let lead = g.exponents(0);
        h = match h {
            Some((old, m)) if old.exponents(0) == lead => {
                let step = Step::new(&field, &m, p);
                let new = combine(ring, &old, &g, &step);
                if new == old {
                    let candidate = divide(&new, &ring.content(&new));
                    if let Some(qa) = ring.div_exact_with(&pa, &candidate, cancel)?
                        && let Some(qb) = ring.div_exact_with(&pb, &candidate, cancel)?
                    {
                        return Ok(done(&candidate, &qa, &qb));
                    }
                }
                Some((new, step.mp))
            }
            Some((old, m)) if old.exponents(0) < lead => Some((old, m)),
            _ => Some((fp.map_into(&g, ring, |&r| lift(r, p)), Integer::from(p))),
        };
    }
    unreachable!("enough word primes")
}

/// `old` modulo `m` and `g` modulo `p` combined term by term, a missing
/// term 0.
fn combine(
    ring: &SparseRing<Integers>,
    old: &Sparse<Integer>,
    g: &Sparse<u64>,
    step: &Step<'_>,
) -> Sparse<Integer> {
    let mut terms = Vec::with_capacity(old.term_count().max(g.term_count()));
    let mut x = old.terms().peekable();
    let mut y = g.terms().peekable();
    loop {
        let term = match (x.peek(), y.peek()) {
            (None, None) => break,
            (Some((e, _)), Some((f, _))) if e == f => {
                let ((e, h), (_, r)) = (x.next().expect("peeked"), y.next().expect("peeked"));
                (e, step.combine(h, *r))
            }
            (Some((e, _)), Some((f, _))) if e < f => {
                let (f, r) = y.next().expect("peeked");
                (f, step.combine(&Integer::ZERO, *r))
            }
            (Some(_), _) => {
                let (e, h) = x.next().expect("peeked");
                (e, step.combine(h, 0))
            }
            (None, Some(_)) => {
                let (f, r) = y.next().expect("peeked");
                (f, step.combine(&Integer::ZERO, *r))
            }
        };
        terms.push(term);
    }
    ring.from_terms(terms)
}

/// The first point of evaluation: past the small integers, the likeliest
/// to be bad points (ours, CAS-137).
const FIRST_POINT: u64 = 1 << 31;

/// The monic gcd of non-zero `a` and `b` over a prime field, by Brown's
/// dense evaluation and interpolation of the last variable that occurs;
/// `None` when the field runs out of points.
pub(crate) fn pgcd<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    cancel: &dyn Cancel,
) -> Result<Option<Sparse<F::Elem>>, Cancelled> {
    pgcd_from(ring, a, b, FIRST_POINT, cancel)
}

/// The terms of `p` grouped by their exponents in the variables other than
/// `x_y`, the last that occurs, which makes each group a run of terms: the
/// exponents with `x_y`'s 0, and the polynomial in `x_y`.
pub(super) type Runs<E> = Vec<(Vec<u64>, Dense<E>)>;

pub(super) fn runs<F: Field + Coefficients + Clone>(
    dense: &DenseRing<F>,
    p: &Sparse<F::Elem>,
    y: usize,
) -> Runs<F::Elem> {
    let mut out: Vec<(Vec<u64>, Vec<F::Elem>)> = Vec::new();
    for (mut e, c) in p.terms() {
        let k = usize::try_from(e[y]).expect("a degree that fits");
        e[y] = 0;
        match out.last_mut() {
            Some((f, coeffs)) if *f == e => coeffs[k] = c.clone(),
            _ => {
                let mut coeffs = vec![dense.base().zero(); k + 1];
                coeffs[k] = c.clone();
                out.push((e, coeffs));
            }
        }
    }
    out.into_iter().map(|(e, c)| (e, dense.poly(c))).collect()
}

pub(super) fn from_runs<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    runs: &Runs<F::Elem>,
    y: usize,
) -> Sparse<F::Elem> {
    ring.from_terms(runs.iter().flat_map(|(e, q)| {
        q.coeffs().iter().enumerate().map(move |(k, c)| {
            let mut f = e.clone();
            f[y] = k as u64;
            (f, c.clone())
        })
    }))
}

/// The gcd in `F[x_y]` of the polynomials of `runs`, monic.
pub(super) fn content<F: Field + Coefficients>(
    dense: &DenseRing<F>,
    runs: &Runs<F::Elem>,
) -> Dense<F::Elem> {
    runs.iter()
        .fold(dense.zero(), |g, (_, q)| EuclideanDomain::gcd(dense, &g, q))
}

/// Each polynomial of `runs` over `c`, which divides them all.
pub(super) fn divide_runs<F: Field + Coefficients>(
    dense: &DenseRing<F>,
    runs: Runs<F::Elem>,
    c: &Dense<F::Elem>,
) -> Runs<F::Elem> {
    if dense.is_one(c) {
        return runs;
    }
    runs.into_iter()
        .map(|(e, q)| (e, dense.div_rem(&q, c).0))
        .collect()
}

/// [`pgcd`] from the points `start, start + 1, …`.
pub(crate) fn pgcd_from<F: Field + Coefficients + Clone>(
    ring: &SparseRing<F>,
    a: &Sparse<F::Elem>,
    b: &Sparse<F::Elem>,
    start: u64,
    cancel: &dyn Cancel,
) -> Result<Option<Sparse<F::Elem>>, Cancelled> {
    let base = ring.base();
    let (da, db) = (a.degrees(), b.degrees());
    let Some(y) = (0..ring.nvars()).rev().find(|&v| da[v] > 0 || db[v] > 0) else {
        return Ok(Some(ring.one()));
    };
    let dense = DenseRing::new(base.clone());
    if (0..y).all(|v| da[v] == 0 && db[v] == 0) {
        // One variable: Euclid's algorithm.
        let to = |p: &Sparse<F::Elem>| ring.to_dense(p, y).expect("one variable");
        let g = EuclideanDomain::gcd(&dense, &to(a), &to(b));
        return Ok(Some(ring.from_dense(&g, y)));
    }
    // a and b in F[x_y][X], X the other variables: their contents in F[x_y],
    // their primitive parts, and γ, the gcd of their leading coefficients.
    let (ra, rb) = (runs(&dense, a, y), runs(&dense, b, y));
    let (ca, cb) = (content(&dense, &ra), content(&dense, &rb));
    let c = ring.from_dense(&EuclideanDomain::gcd(&dense, &ca, &cb), y);
    let (ra, rb) = (divide_runs(&dense, ra, &ca), divide_runs(&dense, rb, &cb));
    let gamma = EuclideanDomain::gcd(&dense, &ra[0].1, &rb[0].1);
    let (pa, pb) = (from_runs(ring, &ra, y), from_runs(ring, &rb, y));
    // The gcd scaled to the leading coefficient γ has at most this degree
    // in x_y.
    let degree = |r: &Runs<F::Elem>| r.iter().filter_map(|(_, q)| q.degree()).max().unwrap_or(0);
    let bound = gamma.degree().unwrap_or(0) + degree(&ra).min(degree(&rb));
    let leading = |h: &Sparse<F::Elem>| {
        let mut e = h.exponents(0);
        e[y] = 0;
        e
    };
    let points = u64::try_from(&base.characteristic()).unwrap_or(u64::MAX);
    // The interpolation so far, its points' product m, their count.
    let mut h: Option<Sparse<F::Elem>> = None;
    let mut m = dense.one();
    let mut count = 0;
    for k in 0..points {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let alpha = base.elem_integer(&(Integer::from(start) + Integer::from(k)));
        let gamma_at = dense.eval(&gamma, &alpha);
        if base.is_zero(&gamma_at) {
            continue;
        }
        let (va, vb) = (ring.eval_var(&pa, y, &alpha), ring.eval_var(&pb, y, &alpha));
        let Some(g) = pgcd_from(ring, &va, &vb, start, cancel)? else {
            return Ok(None);
        };
        if is_constant(&g) {
            return Ok(Some(normalized_sparse(ring, &c)));
        }
        let g = ring.scale(&g, &gamma_at);
        let linear = dense.poly(vec![base.neg(&alpha), base.one()]);
        let unchanged = match h.take() {
            Some(old) if leading(&old) == g.exponents(0) => {
                // Newton's step: h + (g - h(α)) m / m(α).
                let diff = ring.sub(&g, &ring.eval_var(&old, y, &alpha));
                let unchanged = diff.is_zero();
                h = Some(if unchanged {
                    old
                } else {
                    let scale = base.inv(&dense.eval(&m, &alpha)).expect("a new point");
                    ring.add(
                        &old,
                        &ring.mul(&ring.scale(&diff, &scale), &ring.from_dense(&m, y)),
                    )
                });
                m = dense.mul(&m, &linear);
                count += 1;
                unchanged
            }
            // An unlucky point: its gcd is larger.
            Some(old) if leading(&old) < g.exponents(0) => {
                h = Some(old);
                continue;
            }
            // The first point, or the earlier ones were unlucky.
            _ => {
                h = Some(g);
                m = linear;
                count = 1;
                false
            }
        };
        if unchanged || count > bound {
            let interpolated = h.as_ref().expect("interpolated");
            let r = runs(&dense, interpolated, y);
            let primitive = from_runs(
                ring,
                &divide_runs(&dense, r.clone(), &content(&dense, &r)),
                y,
            );
            if ring.div_exact_with(&pa, &primitive, cancel)?.is_some()
                && ring.div_exact_with(&pb, &primitive, cancel)?.is_some()
            {
                return Ok(Some(normalized_sparse(ring, &ring.mul(&c, &primitive))));
            }
        }
    }
    Ok(None)
}
