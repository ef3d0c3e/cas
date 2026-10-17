//! Modular arithmetic: residues modulo a word-size or multi-precision
//! modulus, Chinese remaindering and rational reconstruction.
//!
//! The moduli follow the ring-object pattern: a [`ModWord`] or [`ModRing`]
//! holds the modulus and its precomputed constants, and the arithmetic is
//! done by its methods on plain element values (`u64` for [`ModWord`],
//! [`Residue`] for [`ModRing`]). Elements carry no reference to their ring;
//! passing an element of one ring to another is a caller bug that gives an
//! unspecified (but memory-safe) result.
//!
//! - [`ModWord`]: a modulus that fits a limb, reduced by a precomputed
//!   reciprocal (Möller–Granlund). Elements are `u64` values in `0..n`.
//! - [`ModRing`]: any modulus, with Montgomery multiplication for odd moduli
//!   and division for even ones.
//! - [`crt`] combines two congruences; [`CrtBasis`] reconstructs a value from
//!   its residues modulo many coprime word-size moduli.
//! - [`rational_reconstruction`] recovers a fraction `n/d` from `n d^(-1) mod m`.

mod crt;
mod reconstruct;
mod ring;
mod word;

pub use crt::{CrtBasis, crt};
pub use reconstruct::{rational_reconstruction, rational_reconstruction_with_bounds};
pub use ring::{ModRing, Residue};
pub use word::ModWord;
