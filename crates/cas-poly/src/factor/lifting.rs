//! Multivariate Hensel lifting for Wang's algorithm (CAS-183) \[Wan78,
//! GCL92 §6.7–6.8\]. A factorisation `f(x, a) = ∏ uᵢ(x)` at a point `a` of
//! the variables other than the main one lifts to `f = ∏ Uᵢ`, one variable
//! at a time, with the leading coefficients of the `Uᵢ` in `x` imposed. The
//! arithmetic is modulo a prime `P` above twice a bound on the coefficients
//! of the `Uᵢ`, so their symmetric residues are the factors over ℤ.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithms: f = ∏ u, σ b = c, e for errors, d for degrees"
)]

use cas_num::algebra::{BigPrimeField, EuclideanDomain, Field, Integers, Ring};
use cas_num::modular::Residue;
use cas_num::{Integer, Natural};

use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};
use crate::sparse::{Sparse, SparseRing};

/// The rings of the lifting: sparse polynomials over `𝔽_P` in the
/// variables of `f`, and dense ones in the main variable.
struct Rings {
    field: BigPrimeField,
    sparse: SparseRing<BigPrimeField>,
    dense: DenseRing<BigPrimeField>,
    main: usize,
}

/// The univariate factors at the point and, for each, the inverse `tᵢ`
/// modulo `uᵢ` of the product `bᵢ` of the others. For `deg c < deg ∏ uⱼ`,
/// `σᵢ = c tᵢ mod uᵢ` solves `Σ σᵢ bᵢ = c`: both sides agree modulo every
/// `uᵢ`, so modulo their product, and both have a smaller degree.
struct Univariate {
    factors: Vec<Dense<Residue>>,
    inverses: Vec<Dense<Residue>>,
}

impl Univariate {
    /// `None` when the factors are not pairwise coprime modulo `P`.
    fn new(dense: &DenseRing<BigPrimeField>, factors: Vec<Dense<Residue>>) -> Option<Univariate> {
        let mut inverses = Vec::with_capacity(factors.len());
        for (i, u) in factors.iter().enumerate() {
            let b = factors
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .fold(dense.one(), |acc, (_, v)| {
                    dense.div_rem(&dense.mul(&acc, v), u).1
                });
            let (g, t, _) = EuclideanDomain::gcd_ext(dense, &b, u);
            if !dense.is_one(&g) {
                return None;
            }
            inverses.push(t);
        }
        Some(Univariate { factors, inverses })
    }

    fn solve(&self, dense: &DenseRing<BigPrimeField>, c: &Dense<Residue>) -> Vec<Dense<Residue>> {
        self.factors
            .iter()
            .zip(&self.inverses)
            .map(|(u, t)| {
                let c = dense.div_rem(c, u).1;
                dense.div_rem(&dense.mul(&c, t), u).1
            })
            .collect()
    }
}

/// The coefficient of `(x_var − a)^m` in `e`: `e⁽ᵐ⁾(a) / m!`, the `m`-th
/// derivative in `x_var` at `a`. `P` is above the degrees, so `m!` is
/// invertible.
fn taylor(rings: &Rings, e: &Sparse<Residue>, var: usize, a: &Residue, m: u64) -> Sparse<Residue> {
    let (field, ring) = (&rings.field, &rings.sparse);
    let mut d = e.clone();
    let mut factorial = field.one();
    for k in 1..=m {
        d = ring.derivative(&d, var);
        factorial = field.mul(&factorial, &field.elem_integer(&Integer::from(k)));
    }
    let inverse = field.inv(&factorial).expect("P above the degrees");
    ring.scale(&ring.eval_var(&d, var, a), &inverse)
}

/// `∏ pᵢ`.
fn product<R: Ring>(ring: &SparseRing<R>, ps: &[Sparse<R::Elem>]) -> Sparse<R::Elem> {
    ps.iter().fold(ring.one(), |acc, p| ring.mul(&acc, p))
}

/// The multivariate diophantine equation \[GCL92 Algorithm 6.2\]: the `σᵢ`
/// with `Σ σᵢ bᵢ = c` for `bᵢ = ∏_(j≠i) aⱼ` and `deg_x σᵢ < deg_x aᵢ`,
/// modulo `(y − α)^(d+1)` for each variable `y` of `points` with its value
/// `α`. The `aᵢ` and `c` have only `x` and those variables, and at the point
/// the `aᵢ` are the univariate factors of `base`. The solution at the point
/// of the last variable is corrected by the coefficients of the error's
/// expansion in powers of `y − α`, one power at a time.
fn diophantine(
    rings: &Rings,
    base: &Univariate,
    a: &[Sparse<Residue>],
    c: &Sparse<Residue>,
    points: &[(usize, Residue)],
    d: u64,
    cancel: &dyn Cancel,
) -> Result<Vec<Sparse<Residue>>, Cancelled> {
    let ring = &rings.sparse;
    let Some(((y, alpha), rest)) = points.split_last() else {
        let c = ring
            .to_dense(c, rings.main)
            .expect("only the main variable at the point");
        let sigma = base.solve(&rings.dense, &c);
        return Ok(sigma
            .iter()
            .map(|s| ring.from_dense(s, rings.main))
            .collect());
    };
    let all = product(ring, a);
    let b: Vec<Sparse<Residue>> = a
        .iter()
        .map(|ai| {
            ring.div_exact(&all, ai)
                .expect("a factor divides the product")
        })
        .collect();
    let at: Vec<Sparse<Residue>> = a.iter().map(|ai| ring.eval_var(ai, *y, alpha)).collect();
    let mut sigma = diophantine(
        rings,
        base,
        &at,
        &ring.eval_var(c, *y, alpha),
        rest,
        d,
        cancel,
    )?;
    let combination = |s: &[Sparse<Residue>]| {
        s.iter().zip(&b).fold(ring.zero(), |acc, (si, bi)| {
            ring.add(&acc, &ring.mul(si, bi))
        })
    };
    let mut e = ring.sub(c, &combination(&sigma));
    let linear = ring.sub(&ring.var(*y), &ring.constant(alpha.clone()));
    let mut power = ring.one();
    for m in 1..=d {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        if e.is_zero() {
            break;
        }
        power = ring.mul(&power, &linear);
        let cm = taylor(rings, &e, *y, alpha, m);
        if cm.is_zero() {
            continue;
        }
        let delta: Vec<Sparse<Residue>> = diophantine(rings, base, &at, &cm, rest, d, cancel)?
            .iter()
            .map(|s| ring.mul(s, &power))
            .collect();
        e = ring.sub(&e, &combination(&delta));
        for (s, ds) in sigma.iter_mut().zip(&delta) {
            *s = ring.add(s, ds);
        }
    }
    Ok(sigma)
}

/// `p` with its leading coefficient in the main variable replaced by `lead`.
fn with_lead(rings: &Rings, p: &Sparse<Residue>, lead: &Sparse<Residue>) -> Sparse<Residue> {
    let ring = &rings.sparse;
    let n = ring.degree(p, rings.main).unwrap_or(0);
    let mut exps = vec![0; ring.nvars()];
    exps[rings.main] = n;
    let power = ring.monomial(rings.field.one(), &exps);
    let old = ring.mul(&ring.coeff(p, rings.main, n), &power);
    ring.add(&ring.sub(p, &old), &ring.mul(lead, &power))
}

/// `x` modulo `P` in `(−P/2, P/2]`.
fn symmetric(field: &BigPrimeField, x: &Residue, prime: &Natural) -> Integer {
    let n = field.ring().to_natural(x);
    if &n + &n > *prime {
        Integer::from(n) - Integer::from(prime)
    } else {
        Integer::from(n)
    }
}

/// The lift of `f(x, a) = ∏ uᵢ(x)` to `f = ∏ Uᵢ` with `lc_x(Uᵢ) = lcs[i]`
/// \[GCL92 Algorithm 6.4\]. `x` is the variable `main` of `f`'s ring, and
/// `a` gives the values of the other variables that `f` has, in the order
/// in which they are lifted. The `uᵢ` are pairwise coprime, and
/// `lcs[i](a) = lc(uᵢ)`.
///
/// Variable `y_j` is lifted with the later ones still at their values:
/// - each factor takes `lcs[i]` at those values as its leading coefficient;
/// - the error `e = f − ∏ Uᵢ`, at those values, is corrected one power of
///   `y_j − a_j` at a time, up to the degree in `y_j`, by the diophantine
///   equation of the factors before the step.
///
/// The arithmetic is modulo `prime`. `None` when the `uᵢ` are not coprime
/// modulo it, or when the symmetric residues of the lifts do not multiply to
/// `f`: no factorisation with those leading coefficients lifts from the
/// point. One factor is `f` itself.
///
/// # Panics
///
/// If `prime` is not a prime above the degrees of `f`.
pub(crate) fn lift_variables(
    f: &Sparse<Integer>,
    main: usize,
    points: &[(usize, Integer)],
    factors: &[Dense<Integer>],
    lcs: &[Sparse<Integer>],
    prime: &Natural,
    cancel: &dyn Cancel,
) -> Result<Option<Vec<Sparse<Integer>>>, Cancelled> {
    if let [_] = factors {
        return Ok(Some(vec![f.clone()]));
    }
    let field = BigPrimeField::new(prime.clone()).expect("a prime");
    let nvars = f.nvars();
    let integers = SparseRing::new(Integers, nvars);
    let rings = Rings {
        sparse: SparseRing::new(field.clone(), nvars),
        dense: DenseRing::new(field.clone()),
        field,
        main,
    };
    let (field, ring) = (&rings.field, &rings.sparse);
    let reduce = |p: &Sparse<Integer>| integers.map_into(p, ring, |c| field.elem_integer(c));
    let alphas: Vec<(usize, Residue)> = points
        .iter()
        .map(|(v, a)| (*v, field.elem_integer(a)))
        .collect();
    // images[j]: f with the variables after the j-th at their values.
    let mut images = vec![reduce(f)];
    for (v, a) in alphas.iter().rev() {
        let next = ring.eval_var(images.last().expect("f"), *v, a);
        images.push(next);
    }
    images.reverse();
    let univariate: Vec<Dense<Residue>> = factors
        .iter()
        .map(|u| {
            let coeffs = u.coeffs().iter().map(|c| field.elem_integer(c)).collect();
            rings.dense.poly(coeffs)
        })
        .collect();
    let Some(base) = Univariate::new(&rings.dense, univariate.clone()) else {
        return Ok(None);
    };
    let d = points
        .iter()
        .filter_map(|&(v, _)| integers.degree(f, v))
        .max()
        .unwrap_or(0);
    let mut lifted: Vec<Sparse<Residue>> = univariate
        .iter()
        .map(|u| ring.from_dense(u, main))
        .collect();
    let leads: Vec<Sparse<Residue>> = lcs.iter().map(reduce).collect();
    for (j, (y, alpha)) in alphas.iter().enumerate() {
        let before = lifted.clone();
        for (u, lead) in lifted.iter_mut().zip(&leads) {
            let lead = alphas[j + 1..]
                .iter()
                .fold(lead.clone(), |l, (v, a)| ring.eval_var(&l, *v, a));
            *u = with_lead(&rings, u, &lead);
        }
        let target = &images[j + 1];
        let mut e = ring.sub(target, &product(ring, &lifted));
        let linear = ring.sub(&ring.var(*y), &ring.constant(alpha.clone()));
        let mut power = ring.one();
        for k in 1..=ring.degree(target, *y).unwrap_or(0) {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            if e.is_zero() {
                break;
            }
            power = ring.mul(&power, &linear);
            let c = taylor(&rings, &e, *y, alpha, k);
            if c.is_zero() {
                continue;
            }
            let delta = diophantine(&rings, &base, &before, &c, &alphas[..j], d, cancel)?;
            for (u, du) in lifted.iter_mut().zip(&delta) {
                *u = ring.add(u, &ring.mul(du, &power));
            }
            e = ring.sub(target, &product(ring, &lifted));
        }
    }
    let lifted: Vec<Sparse<Integer>> = lifted
        .iter()
        .map(|u| ring.map_into(u, &integers, |c| symmetric(field, c, prime)))
        .collect();
    Ok((product(&integers, &lifted) == *f).then_some(lifted))
}
