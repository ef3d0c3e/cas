//! Factoring square-free polynomials over word-size prime fields of odd
//! characteristic (CAS-180): the distinct-degree factorisation, then Cantor
//! and Zassenhaus's equal-degree splitting \[CZ81, vzGG13 §14\]. Hensel
//! lifting starts from these local factors.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithms: f = ∏ g, h = x^(p^d), a and b for splitting"
)]

use cas_num::Natural;
use cas_num::algebra::{Field, Ring, WordPrimeField};

use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};
use crate::random::Random;

/// `a b mod f`.
fn mul_mod(
    ring: &DenseRing<WordPrimeField>,
    a: &Dense<u64>,
    b: &Dense<u64>,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Dense<u64>, Cancelled> {
    let product = ring.mul_with(a, b, cancel)?;
    Ok(ring.div_rem_with(&product, f, cancel)?.1)
}

/// `a^e mod f` by binary powering, for `a` reduced modulo `f`.
pub(crate) fn pow_mod(
    ring: &DenseRing<WordPrimeField>,
    a: &Dense<u64>,
    e: &Natural,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Dense<u64>, Cancelled> {
    let mut r = ring.one();
    for i in (0..e.bit_length()).rev() {
        r = mul_mod(ring, &r, &r, f, cancel)?;
        if e.bit(i) {
            r = mul_mod(ring, &r, a, f, cancel)?;
        }
    }
    Ok(r)
}

/// The distinct-degree factorisation of monic, square-free `f`: for each
/// degree `d` that its irreducible factors have, in increasing order, the
/// product `g_d` of those of degree `d`, monic. Round `d` takes
/// `h = x^(p^d)` modulo what is left of `f`, whose irreducible factors of
/// degree dividing `d` are those of `h − x`, and those of degree `d` once
/// the smaller degrees are gone. When what is left has a degree below
/// `2 (d + 1)`, it is irreducible or 1.
pub(crate) fn distinct_degree(
    ring: &DenseRing<WordPrimeField>,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Vec<(Dense<u64>, usize)>, Cancelled> {
    let p = ring.characteristic();
    let x = ring.x();
    let mut rest = f.clone();
    let mut h = x.clone();
    let mut out = Vec::new();
    let mut d = 0;
    loop {
        let n = rest.degree().unwrap_or(0);
        if n < 2 * (d + 1) {
            if n > 0 {
                out.push((rest, n));
            }
            return Ok(out);
        }
        d += 1;
        h = pow_mod(ring, &h, &p, &rest, cancel)?;
        let (g, quotient, _) = ring.gcd_cofactors_with(&rest, &ring.sub(&h, &x), cancel)?;
        if g.degree() > Some(0) {
            rest = quotient;
            if rest.degree() > Some(0) {
                h = ring.div_rem_with(&h, &rest, cancel)?.1;
            }
            out.push((g, d));
        }
    }
}

/// A random polynomial of degree below `n`.
fn random_poly(ring: &DenseRing<WordPrimeField>, n: usize, random: &mut Random) -> Dense<u64> {
    let m = ring.base().ring();
    ring.poly((0..n).map(|_| m.reduce(random.next())).collect())
}

/// The irreducible factors of monic `g`, a product of distinct irreducible
/// polynomials of degree `d`, by Cantor and Zassenhaus's splitting
/// \[CZ81\]. For a random `a` coprime to `g`, `a^((p^d − 1)/2)` is ±1 modulo
/// each irreducible factor, each sign with probability about 1/2, so
/// `gcd(a^((p^d − 1)/2) − 1, g)` is a proper factor with probability about
/// `1 − 2^(1−r)` for `r` factors. The pieces split again until each has
/// degree `d`.
pub(crate) fn equal_degree(
    ring: &DenseRing<WordPrimeField>,
    g: &Dense<u64>,
    d: usize,
    random: &mut Random,
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<u64>>, Cancelled> {
    let p = ring.characteristic();
    assert!(p.bit(0), "an odd characteristic");
    let exponent = (p.pow(d as u64) - 1u64) >> 1;
    let mut pending = vec![g.clone()];
    let mut out = Vec::new();
    while let Some(g) = pending.pop() {
        let n = g.degree().unwrap_or(0);
        if n <= d {
            if n > 0 {
                out.push(g);
            }
            continue;
        }
        let factor = loop {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            let a = random_poly(ring, n, random);
            if a.degree().unwrap_or(0) == 0 {
                continue;
            }
            let (s, _, _) = ring.gcd_cofactors_with(&g, &a, cancel)?;
            if s.degree() > Some(0) {
                break s;
            }
            let b = pow_mod(ring, &a, &exponent, &g, cancel)?;
            let (s, _, _) = ring.gcd_cofactors_with(&g, &ring.sub(&b, &ring.one()), cancel)?;
            if s.degree() > Some(0) && s.degree() < Some(n) {
                break s;
            }
        };
        let (rest, _) = ring.div_rem_with(&g, &factor, cancel)?;
        pending.push(factor);
        pending.push(rest);
    }
    Ok(out)
}

/// The seed of the splitting polynomials.
const SEED: u64 = 0x4341_5331_3830;

/// The irreducible factors of the products of a distinct-degree
/// factorisation, in increasing order of degree, then of coefficients from
/// the constant term up, so that the result does not depend on the random
/// choices.
///
/// # Panics
///
/// If the characteristic is 2.
pub(crate) fn split_degrees(
    ring: &DenseRing<WordPrimeField>,
    ddf: &[(Dense<u64>, usize)],
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<u64>>, Cancelled> {
    let mut random = Random::new(SEED);
    let mut out = Vec::new();
    for (g, d) in ddf {
        out.extend(equal_degree(ring, g, *d, &mut random, cancel)?);
    }
    out.sort_by(|a, b| {
        a.degree()
            .cmp(&b.degree())
            .then_with(|| a.coeffs().cmp(b.coeffs()))
    });
    Ok(out)
}

/// The irreducible factors of monic `g`, a product of distinct irreducible
/// polynomials of degree `d` over `𝔽₂`, by the trace map: for a random `a`,
/// `T(a) = a + a² + a⁴ + ⋯ + a^(2^(d−1))` is 0 or 1 modulo each factor,
/// each with probability 1/2, so `gcd(T(a), g)` splits `g` as often
/// (CAS-453).
fn equal_degree_two(
    ring: &DenseRing<WordPrimeField>,
    g: &Dense<u64>,
    d: usize,
    random: &mut Random,
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<u64>>, Cancelled> {
    let mut pending = vec![g.clone()];
    let mut out = Vec::new();
    while let Some(g) = pending.pop() {
        let n = g.degree().unwrap_or(0);
        if n <= d {
            if n > 0 {
                out.push(g);
            }
            continue;
        }
        let factor = loop {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            let a = random_poly(ring, n, random);
            if a.degree().unwrap_or(0) == 0 {
                continue;
            }
            let mut power = a.clone();
            let mut trace = a;
            for _ in 1..d {
                power = mul_mod(ring, &power, &power, &g, cancel)?;
                trace = ring.add(&trace, &power);
            }
            let (s, _, _) = ring.gcd_cofactors_with(&g, &trace, cancel)?;
            if s.degree() > Some(0) && s.degree() < Some(n) {
                break s;
            }
        };
        let (rest, _) = ring.div_rem_with(&g, &factor, cancel)?;
        pending.push(factor);
        pending.push(rest);
    }
    Ok(out)
}

/// `f`'s `p`-th root for `f` a polynomial in `x^p`: in `𝔽ₚ`, `a^p = a`, so
/// the root of `Σ a_k x^(kp)` is `Σ a_k x^k`.
fn pth_root(ring: &DenseRing<WordPrimeField>, f: &Dense<u64>, p: usize) -> Dense<u64> {
    ring.poly(f.coeffs().iter().step_by(p).copied().collect())
}

/// The square-free decomposition of monic `f` over `𝔽ₚ`, `(part,
/// multiplicity)`, the parts coprime: Yun's steps for the multiplicities
/// that `p` does not divide, then the `p`-th root of what is left, which
/// has every exponent a multiple of `p` (CAS-453).
fn squarefree_parts(
    ring: &DenseRing<WordPrimeField>,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Vec<(Dense<u64>, u64)>, Cancelled> {
    let p = ring.base().ring().modulus();
    let (mut c, mut w, _) = ring.gcd_cofactors_with(f, &ring.derivative(f), cancel)?;
    let mut out = Vec::new();
    let mut i = 1;
    while w.degree() > Some(0) {
        let (y, z, _) = ring.gcd_cofactors_with(&w, &c, cancel)?;
        if z.degree() > Some(0) {
            out.push((z, i));
        }
        i += 1;
        c = ring.div_rem_with(&c, &y, cancel)?.0;
        w = y;
    }
    if c.degree() > Some(0) {
        let root = pth_root(ring, &c, usize::try_from(p).expect("a word"));
        for (part, e) in squarefree_parts(ring, &root, cancel)? {
            out.push((part, e * p));
        }
    }
    Ok(out)
}

/// A factorization over a prime field: the leading coefficient, and the
/// monic irreducible factors with their multiplicities.
pub type FactorsMod = (u64, Vec<(Dense<u64>, u64)>);

impl DenseRing<WordPrimeField> {
    /// The factorization of `f`, not 0, over the prime field: its leading
    /// coefficient, and its monic irreducible factors with their
    /// multiplicities, in increasing degree, then by coefficients
    /// (CAS-453). Square-free parts first, then the distinct-degree
    /// factorization and the equal-degree splitting, Cantor and
    /// Zassenhaus's for an odd `p`, the trace map's for `p = 2`.
    ///
    /// # Panics
    ///
    /// If `f` is zero.
    #[must_use]
    pub fn factor_mod(&self, f: &Dense<u64>) -> FactorsMod {
        crate::cancel::uncancelled(self.factor_mod_with(f, &crate::cancel::Never))
    }

    /// As [`DenseRing::factor_mod`], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `f` is zero.
    pub fn factor_mod_with(
        &self,
        f: &Dense<u64>,
        cancel: &dyn Cancel,
    ) -> Result<FactorsMod, Cancelled> {
        let lead = *f.lead().expect("a polynomial that is not zero");
        let inverse = self.base().inv(&lead).expect("a field");
        let monic = self.scale(f, &inverse);
        let two = self.characteristic() == Natural::from(2u64);
        let mut random = Random::new(SEED);
        let mut out = Vec::new();
        for (part, e) in squarefree_parts(self, &monic, cancel)? {
            let ddf = distinct_degree(self, &part, cancel)?;
            for (g, d) in ddf {
                let pieces = if two {
                    equal_degree_two(self, &g, d, &mut random, cancel)?
                } else {
                    equal_degree(self, &g, d, &mut random, cancel)?
                };
                out.extend(pieces.into_iter().map(|q| (q, e)));
            }
        }
        out.sort_by(|a, b| {
            a.0.degree()
                .cmp(&b.0.degree())
                .then_with(|| a.0.coeffs().cmp(b.0.coeffs()))
                .then_with(|| a.1.cmp(&b.1))
        });
        Ok((lead, out))
    }
}

/// The monic irreducible factors of monic, square-free `f` over a word-size
/// prime field of odd characteristic, in the order of [`split_degrees`].
///
/// # Panics
///
/// If the characteristic is 2.
#[cfg(test)]
pub(crate) fn factor_squarefree(
    ring: &DenseRing<WordPrimeField>,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<u64>>, Cancelled> {
    split_degrees(ring, &distinct_degree(ring, f, cancel)?, cancel)
}
