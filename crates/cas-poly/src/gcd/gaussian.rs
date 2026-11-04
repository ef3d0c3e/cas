//! Brown's modular gcd over ℤ\[i\] \[Bro71\] (ours, CAS-494): the ℤ
//! algorithm of [`super::modular`] modulo word primes `p ≡ 1 (mod 4)`,
//! where `−1` has square roots `±s` and `a + b i ↦ a ± b s` are the two
//! ring maps onto 𝔽ₚ, one per Gaussian prime over `p`. Each image's gcd,
//! from [`pgcd`], is scaled so that its leading coefficient is the image of
//! `γ`, the gcd of the leading coefficients; the two images of a
//! coefficient `c = a + b i`, `v± = a ± b s`, give `a = (v₊ + v₋)/2` and `b =
//! (v₊ − v₋)/(2s)` modulo `p`, and the parts are combined by Chinese
//! remaindering until the result stops changing and its primitive part
//! divides both. A prime whose two images have different leading monomials
//! is passed over: one of them is unlucky.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Brown's paper: a, b, g, h, m, p, s, γ"
)]

use std::collections::BTreeMap;

use cas_num::algebra::{EuclideanDomain, Field, GaussianIntegers, Ring, WordPrimeField};
use cas_num::{Complex, Integer, Rational};

use super::Cofactors;
use super::modular::{Primes, Step, is_constant, lift, pgcd};
use super::subres::Exact;
use crate::cancel::{Cancel, Cancelled};
use crate::sparse::{Sparse, SparseRing};

/// Word primes `p ≡ 1 (mod 4)`, decreasing from `2^62`, each with its field
/// and a square root of `−1` modulo `p`: `t^((p − 1)/4)` for the first
/// small `t` that is not a square.
fn primes() -> impl Iterator<Item = (u64, WordPrimeField, u64)> {
    Primes::new()
        .filter(|(p, _)| p % 4 == 1)
        .filter_map(|(p, field)| {
            let minus_one = field.neg(&field.one());
            (2u64..200)
                .find_map(|t| {
                    let t = field.elem_integer(&Integer::from(t));
                    (field.pow(&t, (p - 1) / 2) == minus_one).then(|| field.pow(&t, (p - 1) / 4))
                })
                .map(|s| (p, field, s))
        })
}

/// The image of the Gaussian integer `c` under `i ↦ root`.
fn embed(field: &WordPrimeField, root: u64, c: &Complex) -> u64 {
    field.add(
        &field.elem_integer(c.re().numer()),
        &field.mul(&root, &field.elem_integer(c.im().numer())),
    )
}

/// The Gaussian integer `a + b i`.
fn gaussian(a: Integer, b: Integer) -> Complex {
    Complex::new(Rational::from(a), Rational::from(b))
}

/// The gcd of non-zero sparse `a` and `b` over ℤ\[i\], with its cofactors.
pub(crate) fn sparse(
    ring: &SparseRing<GaussianIntegers>,
    a: &Sparse<Complex>,
    b: &Sparse<Complex>,
    cancel: &dyn Cancel,
) -> Result<Cofactors<Sparse<Complex>>, Cancelled> {
    let base = ring.base();
    let (ca, cb) = (ring.content(a), ring.content(b));
    let divide = |p: &Sparse<Complex>, c: &Complex| ring.map_into(p, ring, |x| base.exact(x, c));
    let (pa, pb) = (divide(a, &ca), divide(b, &cb));
    let c = base.gcd(&ca, &cb);
    let (ka, kb) = (base.exact(&ca, &c), base.exact(&cb, &c));
    let done = |g: &Sparse<Complex>, qa: &Sparse<Complex>, qb: &Sparse<Complex>| {
        (ring.scale(g, &c), ring.scale(qa, &ka), ring.scale(qb, &kb))
    };
    if is_constant(&pa) || is_constant(&pb) {
        return Ok(done(&ring.one(), &pa, &pb));
    }
    let gamma = base.gcd(&pa.coeffs()[0], &pb.coeffs()[0]);
    // The combination so far, and the product of its primes.
    let mut h: Option<(Sparse<Complex>, Integer)> = None;
    for (p, field, s) in primes() {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let fp = SparseRing::new(field.clone(), ring.nvars());
        let mut images = Vec::with_capacity(2);
        for root in [s, field.neg(&s)] {
            let gp = embed(&field, root, &gamma);
            if gp == 0 {
                break;
            }
            let image = |q: &Sparse<Complex>| ring.map_into(q, &fp, |x| embed(&field, root, x));
            let Some(g) = pgcd(&fp, &image(&pa), &image(&pb), cancel)? else {
                break;
            };
            images.push(fp.scale(&g, &gp));
        }
        let [plus, minus] = images.as_slice() else {
            continue;
        };
        if is_constant(plus) || is_constant(minus) {
            return Ok(done(&ring.one(), &pa, &pb));
        }
        let lead = plus.exponents(0);
        if minus.exponents(0) != lead {
            continue;
        }
        let parts = parts(&field, s, plus, minus);
        h = match h {
            Some((old, m)) if old.exponents(0) == lead => {
                let step = Step::new(&field, &m, p);
                let new = combine(ring, &old, &parts, &step);
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
            // An unlucky prime: its gcd is larger.
            Some((old, m)) if old.exponents(0) < lead => Some((old, m)),
            // The first prime, or the earlier ones were unlucky.
            _ => Some((
                ring.from_terms(
                    parts
                        .into_iter()
                        .map(|(e, (x, y))| (e, gaussian(lift(x, p), lift(y, p)))),
                ),
                Integer::from(p),
            )),
        };
    }
    unreachable!("enough word primes")
}

/// The real and imaginary parts modulo `p` of each coefficient, from its
/// images `v₊` under `i ↦ s` and `v₋` under `i ↦ −s`, a missing term 0.
fn parts(
    field: &WordPrimeField,
    s: u64,
    plus: &Sparse<u64>,
    minus: &Sparse<u64>,
) -> BTreeMap<Vec<u64>, (u64, u64)> {
    let mut images: BTreeMap<Vec<u64>, (u64, u64)> = BTreeMap::new();
    for (e, &v) in plus.terms() {
        images.entry(e).or_default().0 = v;
    }
    for (e, &v) in minus.terms() {
        images.entry(e).or_default().1 = v;
    }
    let two = field.elem_integer(&Integer::from(2));
    let half = field.inv(&two).expect("an odd prime");
    let over = field.inv(&field.mul(&two, &s)).expect("s is a unit");
    images
        .into_iter()
        .map(|(e, (vp, vm))| {
            let re = field.mul(&field.add(&vp, &vm), &half);
            let im = field.mul(&field.sub(&vp, &vm), &over);
            (e, (re, im))
        })
        .collect()
}

/// `old` modulo `m` and the parts modulo `p` combined term by term, a
/// missing term 0.
fn combine(
    ring: &SparseRing<GaussianIntegers>,
    old: &Sparse<Complex>,
    parts: &BTreeMap<Vec<u64>, (u64, u64)>,
    step: &Step<'_>,
) -> Sparse<Complex> {
    let mut terms: BTreeMap<Vec<u64>, (Integer, Integer, u64, u64)> = BTreeMap::new();
    for (e, c) in old.terms() {
        let t = terms
            .entry(e)
            .or_insert((Integer::ZERO, Integer::ZERO, 0, 0));
        t.0 = c.re().numer().clone();
        t.1 = c.im().numer().clone();
    }
    for (e, &(x, y)) in parts {
        let t = terms
            .entry(e.clone())
            .or_insert((Integer::ZERO, Integer::ZERO, 0, 0));
        t.2 = x;
        t.3 = y;
    }
    ring.from_terms(
        terms
            .into_iter()
            .map(|(e, (hr, hi, x, y))| (e, gaussian(step.combine(&hr, x), step.combine(&hi, y)))),
    )
}
