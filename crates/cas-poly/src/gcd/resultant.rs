//! Resultants by remainders over a field, and over ℤ from the resultants
//! modulo enough word primes for Hadamard's bound, by Chinese
//! remaindering.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the identities: a, b, m, n, r"
)]

use cas_num::algebra::{Field, Ring};
use cas_num::modular::CrtBasis;
use cas_num::{Integer, Natural};

use super::modular::Primes;
use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Coefficients, Dense, DenseRing};

/// The resultant of non-zero `a` and `b` over a field, by Euclid's
/// algorithm on `res(a, b) = (-1)^(mn) lc(b)^(m - k) res(b, a mod b)`, `m`,
/// `n` and `k` the degrees of `a`, `b` and the remainder, and `res(a, c) =
/// c^m` for a constant `c`.
pub(crate) fn field<F: Field + Coefficients>(
    ring: &DenseRing<F>,
    a: &Dense<F::Elem>,
    b: &Dense<F::Elem>,
    cancel: &dyn Cancel,
) -> Result<F::Elem, Cancelled> {
    let base = ring.base();
    let mut acc = base.one();
    let (mut a, mut b) = (a.clone(), b.clone());
    loop {
        let (Some(m), Some(n)) = (a.degree(), b.degree()) else {
            unreachable!("non-zero polynomials");
        };
        let lead = b.lead().expect("non-zero").clone();
        if n == 0 {
            return Ok(base.mul(&acc, &base.pow(&lead, m as u64)));
        }
        let r = ring.div_rem_with(&a, &b, cancel)?.1;
        let Some(k) = r.degree() else {
            return Ok(base.zero());
        };
        if m % 2 == 1 && n % 2 == 1 {
            acc = base.neg(&acc);
        }
        acc = base.mul(&acc, &base.pow(&lead, (m - k) as u64));
        (a, b) = (b, r);
    }
}

/// `Σ cᵢ²`, the square of the Euclidean norm of the coefficients.
fn norm_squared(p: &Dense<Integer>) -> Natural {
    p.coeffs().iter().fold(Natural::ZERO, |s, c| {
        let m = c.unsigned_abs();
        &s + &(&m * &m)
    })
}

/// Each word prime of [`Primes`] has more bits than this.
const PRIME_BITS: u64 = 61;

/// The resultant of non-zero `a` and `b` over ℤ: the resultants of their
/// images modulo word primes that divide neither leading coefficient,
/// enough that their product passes twice Hadamard's bound on the
/// Sylvester determinant, `‖a‖^(deg b) ‖b‖^(deg a)`, combined by
/// [`CrtBasis`] in the symmetric range.
pub(crate) fn integers(
    a: &Dense<Integer>,
    b: &Dense<Integer>,
    cancel: &dyn Cancel,
) -> Result<Integer, Cancelled> {
    let (Some(m), Some(n)) = (a.degree(), b.degree()) else {
        unreachable!("non-zero polynomials");
    };
    let (la, lb) = (a.lead().expect("non-zero"), b.lead().expect("non-zero"));
    if n == 0 {
        return Ok(lb.pow(m as u64));
    }
    if m == 0 {
        return Ok(la.pow(n as u64));
    }
    // log₂ of the bound: (n log₂ ‖a‖² + m log₂ ‖b‖²) / 2.
    let bound = (n as u64 * norm_squared(a).bit_length() + m as u64 * norm_squared(b).bit_length())
        .div_ceil(2);
    let leads = la * lb;
    let (mut moduli, mut residues) = (Vec::new(), Vec::new());
    for (p, field) in Primes::new() {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        if field.elem_integer(&leads) == 0 {
            continue;
        }
        let fp = DenseRing::new(field.clone());
        let image = |q: &Dense<Integer>| {
            fp.poly(q.coeffs().iter().map(|x| field.elem_integer(x)).collect())
        };
        residues.push(self::field(&fp, &image(a), &image(b), cancel)?);
        moduli.push(p);
        if PRIME_BITS * moduli.len() as u64 > bound + 1 {
            break;
        }
    }
    let basis = CrtBasis::new(&moduli).expect("distinct primes");
    Ok(basis.combine_symmetric(&residues))
}
