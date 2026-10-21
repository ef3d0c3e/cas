//! Products of coefficient slices: schoolbook, Karatsuba, and the choice of
//! each coefficient ring ([`Coefficients`]).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Karatsuba's identity, r for the ring"
)]

use cas_num::algebra::{
    GaussianIntegers, GaussianRationals, Integers, PrimeField, Rationals, Ring,
};
use cas_num::ball::RealBalls;
use cas_num::modular::{ModRing, ModWord};
use cas_num::{Complex, Integer, Rational};

use super::kronecker;
use crate::cancel::{Cancel, Cancelled};

/// Below this many coefficients in the shorter factor, Karatsuba's method
/// multiplies by schoolbook.
pub(crate) const KARATSUBA_THRESHOLD: usize = 16;

/// The coefficient rings of dense polynomials: ring objects, each with the
/// product of coefficient slices that suits it. The default is Karatsuba's
/// method over schoolbook products; ℤ, ℤ\[i\] and the word-size moduli
/// multiply by Kronecker substitution into one integer product, and balls
/// by schoolbook, whose radii are tighter.
pub trait Coefficients: Ring {
    /// The coefficients of `a b` for non-empty `a` and `b`: `a.len() +
    /// b.len() - 1` of them, the last zero when the leading coefficients'
    /// product is.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    fn mul_coeffs(
        &self,
        a: &[Self::Elem],
        b: &[Self::Elem],
        cancel: &dyn Cancel,
    ) -> Result<Vec<Self::Elem>, Cancelled>
    where
        Self: Sized,
    {
        karatsuba(self, a, b, cancel)
    }

    /// The coefficients of `a²`, as [`Coefficients::mul_coeffs`].
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    fn square_coeffs(
        &self,
        a: &[Self::Elem],
        cancel: &dyn Cancel,
    ) -> Result<Vec<Self::Elem>, Cancelled>
    where
        Self: Sized,
    {
        self.mul_coeffs(a, a, cancel)
    }
}

/// `dst[i] += src[i]`. Past the end of `dst`, `src` holds zeros: the
/// leading terms that cancel in Karatsuba's middle product.
fn add_into<R: Ring>(r: &R, dst: &mut [R::Elem], src: &[R::Elem]) {
    for (d, s) in dst.iter_mut().zip(src) {
        *d = r.add(d, s);
    }
}

/// `dst[i] -= src[i]`, `dst` at least as long as `src`.
fn sub_into<R: Ring>(r: &R, dst: &mut [R::Elem], src: &[R::Elem]) {
    debug_assert!(dst.len() >= src.len());
    for (d, s) in dst.iter_mut().zip(src) {
        *d = r.sub(d, s);
    }
}

/// `a + b`, as long as the longer.
fn sum<R: Ring>(r: &R, a: &[R::Elem], b: &[R::Elem]) -> Vec<R::Elem> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut s = long.to_vec();
    add_into(r, &mut s, short);
    s
}

/// The schoolbook product: each coefficient of `a` times all of `b`.
pub(crate) fn schoolbook<R: Ring>(
    r: &R,
    a: &[R::Elem],
    b: &[R::Elem],
    cancel: &dyn Cancel,
) -> Result<Vec<R::Elem>, Cancelled> {
    let mut c = vec![r.zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        if r.is_zero(x) {
            continue;
        }
        for (cij, y) in c[i..].iter_mut().zip(b) {
            *cij = r.add(cij, &r.mul(x, y));
        }
    }
    Ok(c)
}

/// Karatsuba's product \[KO62\]: with `a = a₀ + xʰ a₁` and `b = b₀ + xʰ b₁`,
/// `a b = a₀b₀ + xʰ((a₀ + a₁)(b₀ + b₁) − a₀b₀ − a₁b₁) + x²ʰ a₁b₁`, three
/// half-size products for four. A factor more than twice as long as the
/// other is cut in blocks of the shorter's length.
pub(crate) fn karatsuba<R: Ring>(
    r: &R,
    a: &[R::Elem],
    b: &[R::Elem],
    cancel: &dyn Cancel,
) -> Result<Vec<R::Elem>, Cancelled> {
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let (n, m) = (a.len(), b.len());
    if m <= KARATSUBA_THRESHOLD {
        return schoolbook(r, a, b, cancel);
    }
    if cancel.cancelled() {
        return Err(Cancelled);
    }
    let mut c = vec![r.zero(); n + m - 1];
    if n >= 2 * m {
        for (k, block) in a.chunks(m).enumerate() {
            let part = karatsuba(r, block, b, cancel)?;
            add_into(r, &mut c[k * m..], &part);
        }
        return Ok(c);
    }
    // n < 2m, so b has at least h coefficients.
    let h = n.div_ceil(2);
    let (a0, a1) = a.split_at(h);
    let (b0, b1) = b.split_at(h);
    let z0 = karatsuba(r, a0, b0, cancel)?;
    let z2 = if b1.is_empty() {
        Vec::new()
    } else {
        karatsuba(r, a1, b1, cancel)?
    };
    let mut z1 = karatsuba(r, &sum(r, a0, a1), &sum(r, b0, b1), cancel)?;
    sub_into(r, &mut z1, &z0);
    sub_into(r, &mut z1, &z2);
    add_into(r, &mut c, &z0);
    add_into(r, &mut c[h..], &z1);
    add_into(r, &mut c[2 * h..], &z2);
    Ok(c)
}

/// Kronecker substitution \[Sch82\] over ℤ, schoolbook for short factors,
/// the shorter the larger the coefficients.
impl Coefficients for Integers {
    fn mul_coeffs(
        &self,
        a: &[Integer],
        b: &[Integer],
        cancel: &dyn Cancel,
    ) -> Result<Vec<Integer>, Cancelled> {
        let (a_bits, b_bits) = (kronecker::max_bits(a), kronecker::max_bits(b));
        if a.len().min(b.len()) <= kronecker::schoolbook_up_to(a_bits + b_bits) {
            return schoolbook(self, a, b, cancel);
        }
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        Ok(kronecker::mul_bits(a, b, a_bits, b_bits))
    }

    fn square_coeffs(&self, a: &[Integer], cancel: &dyn Cancel) -> Result<Vec<Integer>, Cancelled> {
        let bits = kronecker::max_bits(a);
        if a.len() <= kronecker::schoolbook_up_to(2 * bits) {
            return schoolbook(self, a, a, cancel);
        }
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        Ok(kronecker::square_bits(a, bits))
    }
}

/// Over ℤ\[i\], the parts are integer polynomials: with `a = a_re + i a_im`
/// and `b = b_re + i b_im`, three products by Kronecker substitution give
/// `a_re b_re − a_im b_im` and `(a_re + a_im)(b_re + b_im) − a_re b_re −
/// a_im b_im`, the real and imaginary parts.
impl Coefficients for GaussianIntegers {
    fn mul_coeffs(
        &self,
        a: &[Complex],
        b: &[Complex],
        cancel: &dyn Cancel,
    ) -> Result<Vec<Complex>, Cancelled> {
        let (ar, ai) = parts(a);
        let (br, bi) = parts(b);
        let bits =
            |re: &[Integer], im: &[Integer]| kronecker::max_bits(re).max(kronecker::max_bits(im));
        if a.len().min(b.len()) <= kronecker::schoolbook_up_to(bits(&ar, &ai) + bits(&br, &bi)) {
            return schoolbook(self, a, b, cancel);
        }
        if cancel.cancelled() {
            return Err(Cancelled);
        }
        let rr = kronecker::mul(&ar, &br);
        let ii = kronecker::mul(&ai, &bi);
        let sa: Vec<Integer> = ar.iter().zip(&ai).map(|(x, y)| x + y).collect();
        let sb: Vec<Integer> = br.iter().zip(&bi).map(|(x, y)| x + y).collect();
        let ss = kronecker::mul(&sa, &sb);
        Ok(rr
            .iter()
            .zip(&ii)
            .zip(&ss)
            .map(|((rr, ii), ss)| {
                let im = &(ss - rr) - ii;
                Complex::new(Rational::from(rr - ii), Rational::from(im))
            })
            .collect())
    }
}

/// The integer real and imaginary parts of Gaussian integers.
fn parts(z: &[Complex]) -> (Vec<Integer>, Vec<Integer>) {
    z.iter()
        .map(|c| {
            debug_assert!(c.is_gaussian_integer());
            (c.re().numer().clone(), c.im().numer().clone())
        })
        .unzip()
}

/// Kronecker substitution over the residues, then one reduction a
/// coefficient.
impl Coefficients for ModWord {
    fn mul_coeffs(&self, a: &[u64], b: &[u64], cancel: &dyn Cancel) -> Result<Vec<u64>, Cancelled> {
        word_product(self, a, b, cancel)
    }
}

/// The word-size prime fields: as [`ModWord`].
impl Coefficients for PrimeField<ModWord> {
    fn mul_coeffs(&self, a: &[u64], b: &[u64], cancel: &dyn Cancel) -> Result<Vec<u64>, Cancelled> {
        word_product(self.ring(), a, b, cancel)
    }
}

/// A product modulo a word: schoolbook for short factors, otherwise the
/// integer product by Kronecker substitution, reduced.
fn word_product(
    m: &ModWord,
    a: &[u64],
    b: &[u64],
    cancel: &dyn Cancel,
) -> Result<Vec<u64>, Cancelled> {
    if a.len().min(b.len()) <= kronecker::word_schoolbook_up_to(m) {
        return schoolbook(m, a, b, cancel);
    }
    if cancel.cancelled() {
        return Err(Cancelled);
    }
    Ok(kronecker::mul_mod(a, b, m))
}

impl Coefficients for Rationals {}

impl Coefficients for GaussianRationals {}

impl Coefficients for ModRing {}

impl Coefficients for PrimeField<ModRing> {}

/// Schoolbook: the subtractions of Karatsuba's method widen radii.
impl Coefficients for RealBalls {
    fn mul_coeffs(
        &self,
        a: &[Self::Elem],
        b: &[Self::Elem],
        cancel: &dyn Cancel,
    ) -> Result<Vec<Self::Elem>, Cancelled> {
        schoolbook(self, a, b, cancel)
    }
}
