//! EZ-GCD over ℤ (CAS-526) \[MY73, Wang80, GCL92 §7.6\]: the gcd of the
//! images at a point of the variables other than the main one, lifted back
//! to the gcd by the multivariate Hensel lifting of Wang's factoring
//! (CAS-183), with the gcd of the leading coefficients imposed on it.

use std::iter;

use cas_num::Integer;
use cas_num::algebra::{BigPrimeField, Field, Integers, Ring};

use super::Cofactors;
use super::modular::is_constant;
use super::zippel::{Contents, exact, main_variable};
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};
use crate::factor::lifting::{lift_variables, lifting_prime, symmetric};
use crate::random::Random;
use crate::sparse::{Sparse, SparseRing};

/// Rounds of points before EZ-GCD gives up.
const ATTEMPTS: usize = 16;

/// Points drawn at a range of values before it grows.
const TRIES: usize = 20;

/// The largest `k` of the combinations `a + k b` tried when neither `a`
/// nor `b` has an image cofactor coprime to the image gcd.
const MULTIPLES: i64 = 3;

/// The seed of the points.
const SEED: u64 = 0x455A_4743_4443;

/// The gcd of non-zero `a` and `b` over ℤ by EZ-GCD, with its cofactors;
/// `None` when no point within `ATTEMPTS` rounds gives it.
///
/// The contents in the main variable `x` come out first, by gcds in fewer
/// variables. Of two points keeping the leading coefficients in `x`, the
/// image gcd of lower degree is the likelier to be the image of the gcd.
/// Degree 0 means the primitive parts are coprime; the degree of `a` or of
/// `b` means it divides the other, or the point is unlucky. Otherwise `f`,
/// one of `a`, `b` and `a + k b`, has an image cofactor `h` coprime to the
/// image gcd `g`, and `γ f = G H`, `γ` the gcd of the leading coefficients,
/// lifts from `γ(α) f(α) = (γ(α) g / lc g) (lc g h)` with the leading
/// coefficients `γ` and `lc_x f`. The primitive part of `G` in `x` is the
/// gcd when it divides `a` and `b`.
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
    let (g, qa, qb) = if is_constant(&pa) || is_constant(&pb) {
        (ring.one(), pa, pb)
    } else {
        match primitive(ring, &pa, &pb, main, cancel)? {
            Some(parts) => parts,
            None => return Ok(None),
        }
    };
    // a = c (cₐ / c) g qₐ, b = c (c_b / c) g q_b.
    Ok(Some((
        ring.mul(&c, &g),
        ring.mul(&exact(ring, &ca, &c, cancel)?, &qa),
        ring.mul(&exact(ring, &cb, &c, cancel)?, &qb),
    )))
}

/// Points of the variables other than the main one, from a fixed seed, in
/// a range that grows.
struct Points {
    random: Random,
    range: u64,
}

impl Points {
    /// A point of `ys` at which `leads` do not vanish; `None` after `TRIES`
    /// draws.
    fn draw(
        &mut self,
        ring: &SparseRing<Integers>,
        ys: &[usize],
        leads: [&Sparse<Integer>; 2],
    ) -> Option<Vec<(usize, Integer)>> {
        let span = 2 * self.range + 1;
        let range = i64::try_from(self.range).unwrap_or(i64::MAX);
        for _ in 0..TRIES {
            let point: Vec<(usize, Integer)> = ys
                .iter()
                .map(|&y| {
                    let v = i64::try_from(self.random.next() % span).unwrap_or(0) - range;
                    (y, Integer::from(v))
                })
                .collect();
            if leads.iter().all(|l| !at(ring, l, &point).is_zero()) {
                return Some(point);
            }
        }
        None
    }
}

/// `p` with the variables of `point` at their values.
fn at(
    ring: &SparseRing<Integers>,
    p: &Sparse<Integer>,
    point: &[(usize, Integer)],
) -> Sparse<Integer> {
    point
        .iter()
        .fold(p.clone(), |q, (v, c)| ring.eval_var(&q, *v, c))
}

/// The images of `a` and `b` at a point, in the main variable alone: their
/// gcd `g` and its cofactors.
struct Image {
    point: Vec<(usize, Integer)>,
    g: Dense<Integer>,
    ha: Dense<Integer>,
    hb: Dense<Integer>,
}

impl Image {
    fn at(
        ring: &SparseRing<Integers>,
        dense: &DenseRing<Integers>,
        [a, b]: [&Sparse<Integer>; 2],
        x: usize,
        point: Vec<(usize, Integer)>,
        cancel: &dyn Cancel,
    ) -> Result<Image, Cancelled> {
        let univariate = |p| {
            ring.to_dense(&at(ring, p, &point), x)
                .expect("only the main variable is left")
        };
        let (g, ha, hb) = dense.gcd_cofactors_with(&univariate(a), &univariate(b), cancel)?;
        Ok(Image { point, g, ha, hb })
    }

    fn degree(&self) -> usize {
        self.g.degree().unwrap_or(0)
    }
}

/// The gcd of `a` and `b`, primitive in `x` and not constants, with its
/// cofactors; `None` when no round of points gives it.
fn primitive(
    ring: &SparseRing<Integers>,
    a: &Sparse<Integer>,
    b: &Sparse<Integer>,
    x: usize,
    cancel: &dyn Cancel,
) -> Result<Option<Cofactors<Sparse<Integer>>>, Cancelled> {
    let dense = DenseRing::new(Integers);
    let (da, db) = (a.degrees(), b.degrees());
    let ys: Vec<usize> = (0..ring.nvars())
        .filter(|&v| v != x && (da[v] > 0 || db[v] > 0))
        .collect();
    let (la, lb) = (ring.lead_coeff(a, x), ring.lead_coeff(b, x));
    let gamma = ring.gcd_with(&la, &lb, cancel)?;
    let mut points = Points {
        random: Random::new(SEED),
        range: 1,
    };
    for _ in 0..ATTEMPTS {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let mut best: Option<Image> = None;
        for _ in 0..2 {
            let Some(point) = points.draw(ring, &ys, [&la, &lb]) else {
                continue;
            };
            let image = Image::at(ring, &dense, [a, b], x, point, cancel)?;
            if best
                .as_ref()
                .is_none_or(|best| image.degree() < best.degree())
            {
                best = Some(image);
            }
        }
        points.range += 1;
        let Some(image) = best else {
            continue;
        };
        // The leading coefficients do not vanish at the point, so the image
        // gcd has at least the degree of the gcd.
        let degree = image.degree() as u64;
        if degree == 0 {
            return Ok(Some((ring.one(), a.clone(), b.clone())));
        }
        if Some(degree) == ring.degree(a, x) {
            if let Some(q) = ring.div_exact_with(b, a, cancel)? {
                return Ok(Some((a.clone(), ring.one(), q)));
            }
            continue;
        }
        if Some(degree) == ring.degree(b, x) {
            if let Some(q) = ring.div_exact_with(a, b, cancel)? {
                return Ok(Some((b.clone(), q, ring.one())));
            }
            continue;
        }
        let Some(g) = lift(ring, &dense, [a, b], x, &gamma, &image, cancel)? else {
            continue;
        };
        let Some(qa) = ring.div_exact_with(a, &g, cancel)? else {
            continue;
        };
        if let Some(qb) = ring.div_exact_with(b, &g, cancel)? {
            return Ok(Some((g, qa, qb)));
        }
    }
    Ok(None)
}

/// The candidate gcd from the image: the primitive part in `x` of `G` in
/// `γ f = G H`, lifted from the point for the first `f` of `a`, `b` and
/// `a + k b` whose image cofactor is coprime to the image gcd. `None` when
/// none lifts.
fn lift(
    ring: &SparseRing<Integers>,
    dense: &DenseRing<Integers>,
    [a, b]: [&Sparse<Integer>; 2],
    x: usize,
    gamma: &Sparse<Integer>,
    image: &Image,
    cancel: &dyn Cancel,
) -> Result<Option<Sparse<Integer>>, Cancelled> {
    let g = &image.g;
    let g_lead = g.lead().expect("a gcd of non-zero images");
    let gamma_value = at(ring, gamma, &image.point).coeffs()[0].clone();
    let candidates = iter::once((a.clone(), image.ha.clone()))
        .chain(iter::once((b.clone(), image.hb.clone())))
        .chain((1..=MULTIPLES).map(|k| {
            let k = Integer::from(k);
            (
                ring.add(a, &ring.scale(b, &k)),
                dense.add(&image.ha, &dense.scale(&image.hb, &k)),
            )
        }));
    for (f, h) in candidates {
        // f(α) = g h keeps the degree of f, and h is coprime to g.
        let Some(degree) = ring.degree(&f, x) else {
            continue;
        };
        if h.degree().map(|d| (d + image.degree()) as u64) != Some(degree)
            || dense.gcd_with(g, &h, cancel)?.degree() != Some(0)
        {
            continue;
        }
        let target = ring.mul(gamma, &f);
        let prime = lifting_prime(&target);
        let field = BigPrimeField::new(prime.clone()).expect("a prime");
        let Some(inverse) = field.inv(&field.elem_integer(g_lead)) else {
            continue;
        };
        // u = γ(α) g / lc g modulo the prime, v = lc g h: lc u = γ(α),
        // lc v = lc_x f (α), u v = γ(α) f(α).
        let scale = field.mul(&field.elem_integer(&gamma_value), &inverse);
        let u = dense.poly(
            g.coeffs()
                .iter()
                .map(|c| symmetric(&field, &field.mul(&field.elem_integer(c), &scale), &prime))
                .collect(),
        );
        let v = dense.scale(&h, g_lead);
        let vars = target.degrees();
        let point: Vec<(usize, Integer)> = image
            .point
            .iter()
            .filter(|(y, _)| vars[*y] > 0)
            .cloned()
            .collect();
        let lcs = [gamma.clone(), ring.lead_coeff(&f, x)];
        if let Some(lifted) = lift_variables(&target, x, &point, &[u, v], &lcs, &prime, cancel)? {
            return Ok(Some(ring.primitive_part_in(&lifted[0], x)));
        }
    }
    Ok(None)
}
