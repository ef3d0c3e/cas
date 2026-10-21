//! Polynomials for the CAS engine: their types and algorithms over the ring
//! objects of cas-num (CAS-DOC-8). The crate knows nothing of expressions;
//! cas-core converts between the two.
//!
//! - [`Dense`] and [`DenseRing`]: dense univariate polynomials over a ring
//!   object `R`, coefficients lowest degree first, without trailing zeros.
//!   `DenseRing<R>` is a ring object itself, so polynomials over polynomials
//!   and the generic algorithms of `cas_num::algebra` work on it: a
//!   [`Ring`](cas_num::algebra::Ring) always, a
//!   [`EuclideanDomain`](cas_num::algebra::EuclideanDomain) over a field.
//! - [`Coefficients`]: the coefficient rings, each with the product of
//!   coefficient slices that suits it: Karatsuba by default, Kronecker
//!   substitution into one integer product over ℤ, ℤ\[i\] and word-size
//!   prime fields.
//! - [`RationalPoly`] and [`GaussianRationalPoly`]: ℚ\[x\] and ℚ(i)\[x\] as
//!   integer (or Gaussian integer) polynomials over one shared denominator,
//!   as FLINT's `fmpq_poly`, with their ring objects [`RationalPolys`] and
//!   [`GaussianRationalPolys`].
//! - [`Cancel`]: what long computations poll to stop early, with the
//!   never-cancelled default [`Never`].
//!
//! The algorithms, with their sources, are described in `doc/algorithms`.
//!
//! ```
//! use cas_num::algebra::{EuclideanDomain, Integers, Rationals, Ring};
//! use cas_num::{Integer, Rational};
//! use cas_poly::DenseRing;
//!
//! // (x + 1)(x - 1) = x^2 - 1 over ℤ.
//! let zx = DenseRing::new(Integers);
//! let p = zx.mul(&zx.from_i64s(&[1, 1]), &zx.from_i64s(&[-1, 1]));
//! assert_eq!(p, zx.from_i64s(&[-1, 0, 1]));
//! assert_eq!(zx.eval(&p, &Integer::from(3)), Integer::from(8));
//!
//! // x^2 - 1 = (x + 1)(x - 1) + 0 over ℚ.
//! let qx = DenseRing::new(Rationals);
//! let (q, r) = qx.div_rem(&qx.from_i64s(&[-1, 0, 1]), &qx.from_i64s(&[1, 1]));
//! assert_eq!(q, qx.from_i64s(&[-1, 1]));
//! assert!(qx.is_zero(&r));
//! ```

mod cancel;
mod dense;
mod shared;
mod sparse;

pub use cancel::{Cancel, Cancelled, Never};
pub use dense::{Coefficients, Dense, DenseRing};
pub use shared::{
    GaussianRationalPoly, GaussianRationalPolys, Numerators, RationalPoly, RationalPolys, Shared,
    SharedRing,
};
pub use sparse::{Sparse, SparseRing};
