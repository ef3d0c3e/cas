//! Factoring square-free polynomials over word-size prime fields of odd
//! characteristic (CAS-180): the distinct-degree factorisation, then Cantor
//! and Zassenhaus's equal-degree splitting \[CZ81, vzGG13 §14\]. Hensel
//! lifting starts from these local factors.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithms: f = ∏ g, h = x^(p^d), a and b for splitting"
)]
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Zassenhaus's algorithm (CAS-182) is the first caller"
    )
)]

use cas_num::Natural;
use cas_num::algebra::{Ring, WordPrimeField};

use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};

/// The `SplitMix64` generator \[SLF14\], with a fixed seed so that
/// factorisations repeat.
pub(crate) struct Random {
    state: u64,
}

impl Random {
    pub(crate) fn new(seed: u64) -> Random {
        Random { state: seed }
    }

    pub(crate) fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

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

/// The monic irreducible factors of monic, square-free `f` over a word-size
/// prime field of odd characteristic, in increasing order of degree, then
/// of coefficients from the constant term up, so that the result does not
/// depend on the random choices.
///
/// # Panics
///
/// If the characteristic is 2.
pub(crate) fn factor_squarefree(
    ring: &DenseRing<WordPrimeField>,
    f: &Dense<u64>,
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<u64>>, Cancelled> {
    let mut random = Random::new(SEED);
    let mut out = Vec::new();
    for (g, d) in distinct_degree(ring, f, cancel)? {
        out.extend(equal_degree(ring, &g, d, &mut random, cancel)?);
    }
    out.sort_by(|a, b| {
        a.degree()
            .cmp(&b.degree())
            .then_with(|| a.coeffs().cmp(b.coeffs()))
    });
    Ok(out)
}
