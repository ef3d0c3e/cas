//! Chinese remaindering.
#![allow(
    clippy::many_single_char_names,
    reason = "the algorithms use the textbook names (a, b, m, n, g, s, k)"
)]

use super::ModWord;
use crate::{Integer, Limb, Natural};

#[cfg(test)]
mod tests;

/// The solution of `x ≡ a (mod m)`, `x ≡ b (mod n)`, as `(x, l)` with
/// `l = lcm(m, n)` and `x` in `0..l`, or `None` if the congruences are
/// incompatible (`a ≢ b` modulo `gcd(m, n)`). The moduli need not be coprime,
/// and `a` and `b` need not be reduced.
///
/// # Panics
///
/// If `m` or `n` is zero.
#[must_use]
pub fn crt(a: &Natural, m: &Natural, b: &Natural, n: &Natural) -> Option<(Natural, Natural)> {
    assert!(!m.is_zero() && !n.is_zero(), "zero modulus");
    let a = a % m;
    let (g, s, _) = Integer::from(m).gcd_ext(&Integer::from(n));
    // x = a + m k with m k ≡ b - a (mod n), so (m/g) k ≡ (b - a)/g (mod n/g),
    // where s (m/g) ≡ 1 (mod n/g).
    let diff = Integer::from(b % n) - Integer::from(&a);
    if !diff.is_divisible_by(&g) {
        return None;
    }
    let n_g = Integer::from(n).div_exact(&g);
    let k = (diff.div_exact(&g) * s).div_rem_euclid(&n_g).1;
    let k = k.unsigned_abs();
    let n_g = n_g.unsigned_abs();
    let x = a + &(m * &k);
    Some((x, m * &n_g))
}

/// Precomputed data to reconstruct values from their residues modulo
/// pairwise coprime word-size moduli `m_0, …, m_(k-1)`, as in multi-modular
/// algorithms.
///
/// Reconstruction is Garner's mixed-radix algorithm (Knuth §4.3.2): the
/// digits `v_i` with `x = v_0 + v_1 m_0 + v_2 m_0 m_1 + …` are found one
/// modulus at a time with word arithmetic, then `x` is assembled by Horner's
/// rule. Both steps take `O(k^2)` word operations.
#[derive(Clone, Debug)]
pub struct CrtBasis {
    moduli: Vec<ModWord>,
    /// `c_i = (m_0 ⋯ m_(i-1))^(-1) mod m_i`.
    inverses: Vec<Limb>,
    /// `m_0 ⋯ m_(k-1)`.
    product: Natural,
}

impl CrtBasis {
    /// The basis for `moduli`, or `None` if two of them share a factor. The
    /// empty basis has product 1 and reconstructs every value as 0.
    ///
    /// # Panics
    ///
    /// If a modulus is zero.
    #[must_use]
    pub fn new(moduli: &[Limb]) -> Option<CrtBasis> {
        let moduli: Vec<ModWord> = moduli.iter().map(|&m| ModWord::new(m)).collect();
        let mut inverses = Vec::with_capacity(moduli.len());
        for (i, w) in moduli.iter().enumerate() {
            let prefix = moduli[..i]
                .iter()
                .fold(w.reduce(1), |p, m| w.mul(p, w.reduce(m.modulus())));
            // Invertible iff m_i is coprime to every earlier modulus.
            inverses.push(w.inv(prefix)?);
        }
        let product = moduli.iter().fold(Natural::one(), |p, m| p * m.modulus());
        Some(CrtBasis {
            moduli,
            inverses,
            product,
        })
    }

    /// The number of moduli.
    #[must_use]
    pub fn len(&self) -> usize {
        self.moduli.len()
    }

    /// Whether the basis has no moduli.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.moduli.is_empty()
    }

    /// The moduli.
    pub fn moduli(&self) -> impl ExactSizeIterator<Item = Limb> + '_ {
        self.moduli.iter().map(ModWord::modulus)
    }

    /// The product `M` of the moduli.
    #[must_use]
    pub fn modulus(&self) -> &Natural {
        &self.product
    }

    /// The `x` in `0..M` with `x ≡ residues[i] (mod m_i)` for every `i`.
    /// Residues need not be reduced.
    ///
    /// # Panics
    ///
    /// If there are not as many residues as moduli.
    #[must_use]
    pub fn combine(&self, residues: &[Limb]) -> Natural {
        assert_eq!(residues.len(), self.moduli.len(), "one residue per modulus");
        let mut digits: Vec<Limb> = Vec::with_capacity(residues.len());
        for (i, (w, &r)) in self.moduli.iter().zip(residues).enumerate() {
            // The value of the digits so far, modulo m_i.
            let so_far = digits
                .iter()
                .zip(&self.moduli)
                .rev()
                .fold(0, |acc, (&v, m)| {
                    w.add(w.mul(acc, w.reduce(m.modulus())), w.reduce(v))
                });
            digits.push(w.mul(w.sub(w.reduce(r), so_far), self.inverses[i]));
        }
        digits
            .iter()
            .zip(&self.moduli)
            .rev()
            .fold(Natural::ZERO, |x, (&v, m)| x * m.modulus() + v)
    }

    /// The `x` in `(-M/2, M/2]` with `x ≡ residues[i] (mod m_i)` for every
    /// `i`, for values known to be small in absolute value.
    ///
    /// # Panics
    ///
    /// If there are not as many residues as moduli.
    #[must_use]
    pub fn combine_symmetric(&self, residues: &[Limb]) -> Integer {
        let x = self.combine(residues);
        if &x + &x > self.product {
            Integer::from(x) - Integer::from(&self.product)
        } else {
            Integer::from(x)
        }
    }
}
