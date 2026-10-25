//! Factoring multivariate polynomials over ℤ and ℚ by Wang's algorithm
//! (CAS-184) \[Wan78\], as SymPy's `dmp_zz_wang` arranges it \[SymPy\]. A
//! square-free part loses its content in a main variable, which is factored
//! in the other variables. The primitive part is factored at a point of the
//! other variables, its factors take their leading coefficients from those
//! of the leading coefficient, and they lift back (CAS-183).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Wang's paper: f, its image u = δ ∏ h, the leading coefficient's factors F and their values E"
)]

use std::cmp::Ordering;
use std::collections::HashSet;

use cas_num::algebra::{Integers, Rationals, Ring};
use cas_num::{Integer, Natural, Rational};

use super::lifting::lift_variables;
use super::zassenhaus;
use super::{FactorError, Factors, Search, unbounded};
use crate::cancel::{Cancel, Cancelled, Never, uncancelled};
use crate::dense::{Dense, DenseRing};
use crate::gcd::split_sparse;
use crate::random::Random;
use crate::sparse::{Sparse, SparseRing};

/// How many good points are gathered before one is chosen. Those whose
/// image has the fewest factors are the likeliest to have no extraneous
/// ones.
const CONFIGURATIONS: usize = 3;

/// How many points are tried at a range of values before it grows.
const TRIES: usize = 20;

/// The seed of the points.
const SEED: u64 = 0x4341_5331_3834;

/// The variables that `p` has.
fn variables<E>(p: &Sparse<E>) -> Vec<usize> {
    p.degrees()
        .iter()
        .enumerate()
        .filter(|&(_, &d)| d > 0)
        .map(|(v, _)| v)
        .collect()
}

/// `p` at the point.
fn at(
    ring: &SparseRing<Integers>,
    p: &Sparse<Integer>,
    point: &[(usize, Integer)],
) -> Sparse<Integer> {
    point
        .iter()
        .fold(p.clone(), |q, (v, a)| ring.eval_var(&q, *v, a))
}

/// The value of a polynomial without variables.
fn value(p: &Sparse<Integer>) -> Integer {
    p.coeffs().first().cloned().unwrap_or(Integer::ZERO)
}

/// Wang's condition on the values `E_j = F_j(a)` of the irreducible factors
/// of the leading coefficient: each has a prime factor that divides neither
/// an earlier one nor `δ Ω`, for `δ` the content of the image and `Ω` the
/// constant of the leading coefficient. That prime tells which factor of the
/// image `F_j` belongs to.
fn distinct_primes(values: &[Integer], content: &Integer, constant: &Integer) -> bool {
    let mut seen = vec![(content * constant).abs()];
    for e in values {
        let mut q = e.abs();
        for r in seen.iter().rev() {
            let mut r = r.clone();
            while !r.is_one() {
                r = r.gcd(&q);
                q = q.div_exact(&r);
            }
            if q.is_one() {
                return false;
            }
        }
        seen.push(q);
    }
    true
}

/// A good point: the image `u = δ ∏ hᵢ` there, its content `δ` with the sign
/// that makes the leading coefficient of `∏ hᵢ` positive, the irreducible
/// factors `hᵢ`, and the values `E_j` of the leading coefficient's factors.
struct Configuration {
    point: Vec<(usize, Integer)>,
    content: Integer,
    norm: Natural,
    factors: Vec<Dense<Integer>>,
    values: Vec<Integer>,
}

/// The configuration at `point`: `None` unless the leading coefficient does
/// not vanish there, the image is square-free, and Wang's condition holds.
fn configuration(
    ring: &SparseRing<Integers>,
    f: &Sparse<Integer>,
    x: usize,
    point: Vec<(usize, Integer)>,
    leading: &(Integer, Vec<(Sparse<Integer>, u64)>),
    search: &Search<'_>,
) -> Result<Option<Configuration>, FactorError> {
    let zx = DenseRing::new(Integers);
    let (constant, lc_factors) = leading;
    let u = ring
        .to_dense(&at(ring, f, &point), x)
        .expect("only the main variable at the point");
    if u.degree().map(|d| d as u64) != ring.degree(f, x)
        || zx
            .gcd_with(&u, &zx.derivative(&u), search.cancel())?
            .degree()
            != Some(0)
    {
        return Ok(None);
    }
    let values: Vec<Integer> = lc_factors
        .iter()
        .map(|(g, _)| value(&at(ring, g, &point)))
        .collect();
    let lead = u.lead().cloned().unwrap_or(Integer::ONE);
    let mut content = zx.content(&u);
    if lead.is_negative() {
        content = -content;
    }
    if !distinct_primes(&values, &content, constant) {
        return Ok(None);
    }
    let h = zx.primitive_part(&u);
    let norm = h
        .coeffs()
        .iter()
        .map(Integer::unsigned_abs)
        .max()
        .unwrap_or(Natural::ZERO);
    let factors = zx
        .factor_search(&h, search)?
        .into_parts()
        .1
        .into_iter()
        .map(|(g, _)| g)
        .collect();
    Ok(Some(Configuration {
        point,
        content,
        norm,
        factors,
        values,
    }))
}

/// What the search for points finds.
enum Found {
    /// An image with one factor: `f` is irreducible.
    Irreducible,
    Points(Configuration),
}

/// The points of the search: values in `[−range, range]`, `range` growing
/// by 2 after [`TRIES`] points, none twice.
struct Points {
    random: Random,
    range: u64,
    tried: HashSet<Vec<i64>>,
}

impl Points {
    fn next(&mut self, n: usize) -> Option<Vec<i64>> {
        let span = 2 * self.range + 1;
        let range = i64::try_from(self.range).unwrap_or(i64::MAX);
        let values: Vec<i64> = (0..n)
            .map(|_| i64::try_from(self.random.next() % span).unwrap_or(0) - range)
            .collect();
        self.tried.insert(values.clone()).then_some(values)
    }
}

/// [`CONFIGURATIONS`] good points whose images have the fewest factors seen,
/// and of those the one whose image has the smallest coefficients. Each
/// point is a step of `search`.
fn find_points(
    ring: &SparseRing<Integers>,
    f: &Sparse<Integer>,
    x: usize,
    ys: &[usize],
    leading: &(Integer, Vec<(Sparse<Integer>, u64)>),
    points: &mut Points,
    search: &Search<'_>,
) -> Result<Found, FactorError> {
    let mut found: Vec<Configuration> = Vec::new();
    while found.len() < CONFIGURATIONS {
        for _ in 0..TRIES {
            search.step()?;
            let Some(values) = points.next(ys.len()) else {
                continue;
            };
            let point = ys
                .iter()
                .zip(&values)
                .map(|(&v, &a)| (v, Integer::from(a)))
                .collect();
            let Some(config) = configuration(ring, f, x, point, leading, search)? else {
                continue;
            };
            if config.factors.len() == 1 {
                return Ok(Found::Irreducible);
            }
            match found
                .first()
                .map(|c| config.factors.len().cmp(&c.factors.len()))
            {
                Some(Ordering::Greater) => continue,
                Some(Ordering::Less) => found.clear(),
                _ => {}
            }
            found.push(config);
            if found.len() == CONFIGURATIONS {
                break;
            }
        }
        if found.len() < CONFIGURATIONS {
            points.range += 2;
        }
    }
    let best = found
        .into_iter()
        .min_by(|a, b| a.norm.cmp(&b.norm))
        .expect("configurations");
    Ok(Found::Points(best))
}

/// The polynomial to lift, the factors at the point and their leading
/// coefficients.
type Lifting = (Sparse<Integer>, Vec<Dense<Integer>>, Vec<Sparse<Integer>>);

/// Wang's leading coefficients. For each factor `hᵢ` of the image, the
/// `F_j` whose values divide `lc(hᵢ) δ`, as often as they do, times a
/// constant that makes the product `lc(hᵢ)` at the point. The content `δ`
/// is spread over the factors, and what remains of it multiplies each
/// factor and its leading coefficient, `f` taking its `(r − 1)`-th power.
/// `None` when some `F_j` goes to no factor or a division is not exact: the
/// image has extraneous factors.
fn leading_coefficients(
    ring: &SparseRing<Integers>,
    f: &Sparse<Integer>,
    lc_factors: &[(Sparse<Integer>, u64)],
    config: &Configuration,
) -> Option<Lifting> {
    let zx = DenseRing::new(Integers);
    let mut assigned = vec![false; lc_factors.len()];
    let mut candidates = Vec::with_capacity(config.factors.len());
    for h in &config.factors {
        let mut c = ring.one();
        let mut d = h.lead().cloned().unwrap_or(Integer::ONE) * &config.content;
        for j in (0..lc_factors.len()).rev() {
            let e = &config.values[j];
            let mut k = 0;
            while d.is_divisible_by(e) {
                d = d.div_exact(e);
                k += 1;
            }
            if k > 0 {
                c = ring.mul(&c, &ring.pow(&lc_factors[j].0, k));
                assigned[j] = true;
            }
        }
        candidates.push(c);
    }
    if !assigned.iter().all(|&a| a) {
        return None;
    }
    let mut content = config.content.clone();
    let mut lcs = Vec::with_capacity(candidates.len());
    let mut factors = Vec::with_capacity(candidates.len());
    for (c, h) in candidates.iter().zip(&config.factors) {
        let d = value(&at(ring, c, &config.point));
        let lc = h.lead().cloned().unwrap_or(Integer::ONE);
        let (scale, h) = if content.is_one() {
            if !lc.is_divisible_by(&d) {
                return None;
            }
            (lc.div_exact(&d), h.clone())
        } else {
            let g = lc.gcd(&d);
            let d = d.div_exact(&g);
            if !content.is_divisible_by(&d) {
                return None;
            }
            content = content.div_exact(&d);
            (lc.div_exact(&g), zx.scale(h, &d))
        };
        lcs.push(ring.scale(c, &scale));
        factors.push(h);
    }
    if content.is_one() {
        return Some((f.clone(), factors, lcs));
    }
    let r = factors.len() as u64;
    Some((
        ring.scale(f, &content.pow(r - 1)),
        factors.iter().map(|h| zx.scale(h, &content)).collect(),
        lcs.iter().map(|l| ring.scale(l, &content)).collect(),
    ))
}

/// A prime above twice `2^(Σ deg) ‖g‖₂`, summing the degrees in every
/// variable. The coefficients of any factor `h` of `g` are within
/// `2^(Σ deg h) M(h) ≤ 2^(Σ deg g) M(g) ≤ 2^(Σ deg g) ‖g‖₂` by Mahler's and
/// Landau's inequalities, `M` the Mahler measure, which is at least 1 for a
/// non-zero integer polynomial.
fn lifting_prime(g: &Sparse<Integer>) -> Natural {
    let squares = g.coeffs().iter().fold(Natural::ZERO, |acc, c| {
        let a = c.unsigned_abs();
        &acc + &(&a * &a)
    });
    let degrees: u64 = g.degrees().iter().sum();
    ((squares.isqrt() + 1) << (degrees + 1)).next_prime()
}

/// The irreducible factors of `f`, primitive over ℤ and in `x`, square-free,
/// with a variable other than `x`. Points are tried until Wang's leading
/// coefficients lift: in theory with probability 1, as most points have
/// images without extraneous factors.
fn wang(
    ring: &SparseRing<Integers>,
    f: &Sparse<Integer>,
    x: usize,
    search: &Search<'_>,
) -> Result<Vec<Sparse<Integer>>, FactorError> {
    let ys: Vec<usize> = variables(f).into_iter().filter(|&v| v != x).collect();
    let leading = ring
        .factor_search(&ring.lead_coeff(f, x), search)?
        .into_parts();
    let mut points = Points {
        random: Random::new(SEED),
        range: 1,
        tried: HashSet::new(),
    };
    loop {
        let config = match find_points(ring, f, x, &ys, &leading, &mut points, search)? {
            Found::Irreducible => return Ok(vec![f.clone()]),
            Found::Points(config) => config,
        };
        if let Some((g, factors, lcs)) = leading_coefficients(ring, f, &leading.1, &config) {
            let prime = lifting_prime(&g);
            if let Some(lifted) = lift_variables(
                &g,
                x,
                &config.point,
                &factors,
                &lcs,
                &prime,
                search.cancel(),
            )? {
                return Ok(lifted.iter().map(|u| ring.primitive_part(u)).collect());
            }
        }
        points.range += 1;
    }
}

/// The total degree.
fn total_degree<E>(p: &Sparse<E>) -> u64 {
    (0..p.term_count())
        .map(|k| p.exponents(k).iter().sum())
        .max()
        .unwrap_or(0)
}

/// The order of `factor`'s factors: by total degree, then by terms from the
/// leading one, each by exponents and then coefficient, then by exponent.
fn sort(factors: &mut [(Sparse<Integer>, u64)]) {
    factors.sort_by(|(a, e), (b, x)| {
        total_degree(a)
            .cmp(&total_degree(b))
            .then_with(|| {
                a.terms()
                    .zip(b.terms())
                    .map(|((f, c), (g, d))| f.cmp(&g).then_with(|| c.cmp(d)))
                    .find(|o| o.is_ne())
                    .unwrap_or_else(|| a.term_count().cmp(&b.term_count()))
            })
            .then(e.cmp(x))
    });
}

impl SparseRing<Integers> {
    /// The factorisation of `p` into irreducible polynomials over ℤ,
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is the content with the sign of the
    /// lexicographic leading coefficient. The `fᵢ` are distinct, irreducible
    /// and primitive, with positive lexicographic leading coefficients. They
    /// are sorted by total degree, then by terms from the leading one. Zero
    /// gives `c = 0` and no factors. giac: `factor(x^4-y^4)` is
    /// `(x-y)*(x+y)*(x^2+y^2)`.
    #[must_use]
    pub fn factor(&self, p: &Sparse<Integer>) -> Factors<Integer, Sparse<Integer>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Sparse<Integer>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Integer, Sparse<Integer>>, Cancelled> {
        unbounded(self.factor_within(p, u64::MAX, cancel))
    }

    /// As `factor_with`, with at most `bound` steps of the searches in all:
    /// subsets of local factors in Zassenhaus's recombination, and points in
    /// Wang's algorithm.
    ///
    /// # Errors
    ///
    /// [`FactorError::Cancelled`] when `cancel` says so, and
    /// [`FactorError::Search`] past the bound.
    pub fn factor_within(
        &self,
        p: &Sparse<Integer>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Integer, Sparse<Integer>>, FactorError> {
        self.factor_search(p, &Search::new(cancel, bound))
    }

    /// As `factor_within`, the steps counted by `search`.
    pub(crate) fn factor_search(
        &self,
        p: &Sparse<Integer>,
        search: &Search<'_>,
    ) -> Result<Factors<Integer, Sparse<Integer>>, FactorError> {
        let (constant, parts) = self.squarefree_with(p, search.cancel())?.into_parts();
        let mut factors = Vec::new();
        for (f, e) in parts {
            for g in self.factor_squarefree(&f, search)? {
                factors.push((g, e));
            }
        }
        sort(&mut factors);
        Ok(Factors { constant, factors })
    }

    /// The irreducible factors of square-free, primitive `f` with a positive
    /// leading coefficient. With one variable, Zassenhaus's algorithm. With
    /// several, the main variable is one of least degree: the content in it
    /// is factored in the others, and the primitive part by Wang's algorithm.
    fn factor_squarefree(
        &self,
        f: &Sparse<Integer>,
        search: &Search<'_>,
    ) -> Result<Vec<Sparse<Integer>>, FactorError> {
        let cancel = search.cancel();
        let vars = variables(f);
        let degrees = f.degrees();
        let Some(&x) = vars.iter().min_by_key(|&&v| degrees[v]) else {
            return Ok(Vec::new());
        };
        if let [v] = vars.as_slice() {
            let dense = self.to_dense(f, *v).expect("one variable");
            return Ok(zassenhaus::factor_squarefree(&dense, search)?
                .iter()
                .map(|g| self.from_dense(g, *v))
                .collect());
        }
        let (content, primitive) = self.content_split(f, x, cancel)?;
        let mut out = if self.is_one(&content) {
            Vec::new()
        } else {
            self.factor_squarefree(&content, search)?
        };
        if variables(&primitive).len() == 1 {
            out.extend(self.factor_squarefree(&primitive, search)?);
        } else {
            out.extend(wang(self, &primitive, x, search)?);
        }
        Ok(out)
    }
}

impl SparseRing<Rationals> {
    /// The factorisation of `p` into irreducible polynomials over ℚ,
    /// `p = c f₁^e₁ ⋯ fₖ^eₖ`. `c` is the lexicographic leading coefficient,
    /// and the `fᵢ` are monic in the lexicographic order, in the order of
    /// the factors of the numerator over ℤ. Zero gives `c = 0`.
    #[must_use]
    pub fn factor(&self, p: &Sparse<Rational>) -> Factors<Rational, Sparse<Rational>> {
        uncancelled(self.factor_with(p, &Never))
    }

    /// As `factor`, polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    pub fn factor_with(
        &self,
        p: &Sparse<Rational>,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Rational, Sparse<Rational>>, Cancelled> {
        let Some(lead) = p.coeffs().first() else {
            return Ok(Factors {
                constant: Rational::ZERO,
                factors: Vec::new(),
            });
        };
        let integers = SparseRing::new(Integers, self.nvars());
        let (numerator, _) = split_sparse::<Integers>(self, &integers, p);
        let factors = integers
            .factor_with(&numerator, cancel)?
            .into_parts()
            .1
            .into_iter()
            .map(|(g, e)| {
                let lc = g.coeffs().first().cloned().unwrap_or(Integer::ONE);
                let monic = integers.map_into(&g, self, |c| Rational::new(c.clone(), lc.clone()));
                (monic, e)
            })
            .collect();
        Ok(Factors {
            constant: lead.clone(),
            factors,
        })
    }
}
