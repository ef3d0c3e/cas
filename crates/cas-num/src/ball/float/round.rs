//! Rounding to a precision: one routine over a top-aligned mantissa, shared by
//! every rounded operation.

use std::cmp::Ordering;

use super::{Float, LIMB_BITS};
use crate::kernel;
use crate::{Integer, Limb};

/// A rounding direction.
///
/// MPFR's names differ: `Down`/`Up` here are toward/away from zero (MPFR's
/// `RNDZ`/`RNDA`), `Floor`/`Ceil` are toward −∞/+∞ (`RNDD`/`RNDU`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Round {
    /// Toward zero (truncation).
    Down,
    /// Away from zero.
    Up,
    /// Toward −∞.
    Floor,
    /// Toward +∞.
    Ceil,
    /// To the nearest value, ties to an even last bit.
    Nearest,
}

impl Round {
    /// Every mode.
    pub const ALL: [Round; 5] = [
        Round::Down,
        Round::Up,
        Round::Floor,
        Round::Ceil,
        Round::Nearest,
    ];

    /// Whether rounding the magnitude of a value of the given sign moves it up
    /// one unit in the last kept place, given the first dropped bit, whether
    /// any later dropped bit is set, and the last kept bit.
    #[allow(
        clippy::fn_params_excessive_bools,
        reason = "the four bits rounding decides on"
    )]
    fn increments(self, neg: bool, round_bit: bool, sticky: bool, last_bit: bool) -> bool {
        let inexact = round_bit || sticky;
        match self {
            Round::Down => false,
            Round::Up => inexact,
            Round::Floor => neg && inexact,
            Round::Ceil => !neg && inexact,
            Round::Nearest => round_bit && (sticky || last_bit),
        }
    }
}

/// The result of [`round_mant`].
pub(super) struct Rounded {
    /// Index of the lowest kept limb.
    pub(super) low: usize,
    /// Whether the increment carried out of the top: the kept limbs are then
    /// `0.1000…` and the exponent goes up by one.
    pub(super) carry: bool,
    /// `rounded.cmp(exact)` for the signed values.
    pub(super) ord: Ordering,
}

/// Rounds the top-aligned mantissa `m` (top bit set) of a value with the
/// given sign to `prec` bits in place: the bits below are cleared and the kept
/// part incremented if the mode says so. `sticky` stands for set bits below
/// `m`; when it is true, `m` must hold at least `prec + 1` bits, so that the
/// first dropped bit lies inside it.
pub(super) fn round_mant(
    m: &mut [Limb],
    sticky: bool,
    neg: bool,
    prec: u64,
    rnd: Round,
) -> Rounded {
    let n = m.len();
    debug_assert!(n > 0 && m[n - 1] >> (Limb::BITS - 1) == 1);
    let total = LIMB_BITS * n as u64;
    if prec >= total {
        debug_assert!(!sticky, "round_mant: no room for the round bit");
        return Rounded {
            low: 0,
            carry: false,
            ord: Ordering::Equal,
        };
    }
    // `prec < total` fits a `usize`.
    let kept = usize::try_from(prec.div_ceil(LIMB_BITS)).expect("below the limb count");
    let low = n - kept;
    #[allow(clippy::cast_possible_truncation, reason = "below 64")]
    let z = (LIMB_BITS * kept as u64 - prec) as u32;
    let (round_bit, mut rest) = if z > 0 {
        let below = m[low] & ((1 << z) - 1);
        let half = 1 << (z - 1);
        (below & half != 0, below & (half - 1) != 0)
    } else {
        let next = m[low - 1];
        (next >> (Limb::BITS - 1) == 1, next << 1 != 0)
    };
    let below_end = if z > 0 { low } else { low - 1 };
    rest |= sticky || m[..below_end].iter().any(|&x| x != 0);
    let last_bit = (m[low] >> z) & 1 == 1;
    if z > 0 {
        m[low] &= !((1 << z) - 1);
    }
    if !round_bit && !rest {
        return Rounded {
            low,
            carry: false,
            ord: Ordering::Equal,
        };
    }
    let up = rnd.increments(neg, round_bit, rest, last_bit);
    let mut carry = false;
    if up && kernel::add_1_assign(&mut m[low..], 1 << z) != 0 {
        m[n - 1] = 1 << (Limb::BITS - 1);
        carry = true;
    }
    // The magnitude went up or down; the signed value moves with the sign.
    let mag = if up {
        Ordering::Greater
    } else {
        Ordering::Less
    };
    Rounded {
        low,
        carry,
        ord: if neg { mag.reverse() } else { mag },
    }
}

/// [`round_mant`] for a two-limb mantissa `(hi, lo)` rounded to `prec <= 64`
/// bits, straight-line: the kept limb (top bit set), whether the exponent goes
/// up by one, and `rounded.cmp(exact)`. `sticky` stands for set bits below
/// `lo`.
#[inline]
pub(super) fn round_1(
    neg: bool,
    hi: Limb,
    lo: Limb,
    sticky: bool,
    prec: u64,
    rnd: Round,
) -> (Limb, bool, Ordering) {
    debug_assert!((1..=LIMB_BITS).contains(&prec) && hi >> (Limb::BITS - 1) == 1);
    #[allow(clippy::cast_possible_truncation, reason = "below 64")]
    let z = (LIMB_BITS - prec) as u32;
    let (round_bit, rest) = if z == 0 {
        (lo >> (Limb::BITS - 1) == 1, lo << 1 != 0 || sticky)
    } else {
        let half = 1 << (z - 1);
        (hi & half != 0, hi & (half - 1) != 0 || lo != 0 || sticky)
    };
    let mask = if z == 0 { 0 } else { (1 << z) - 1 };
    let kept = hi & !mask;
    if !round_bit && !rest {
        return (kept, false, Ordering::Equal);
    }
    let up = rnd.increments(neg, round_bit, rest, (hi >> z) & 1 == 1);
    let mag = if up {
        Ordering::Greater
    } else {
        Ordering::Less
    };
    let ord = if neg { mag.reverse() } else { mag };
    if !up {
        return (kept, false, ord);
    }
    match kept.overflowing_add(1 << z) {
        (m, false) => (m, false, ord),
        (_, true) => (1 << (Limb::BITS - 1), true, ord),
    }
}

impl Float {
    /// `self` rounded to `prec` bits in the given mode, and
    /// `rounded.cmp(self)`.
    ///
    /// # Panics
    ///
    /// If `prec` is 0.
    #[must_use]
    pub fn round(&self, prec: u64, rnd: Round) -> (Float, Ordering) {
        assert!(prec >= 1, "Float::round: precision 0");
        if self.is_zero() || self.bit_length() <= prec {
            return (self.clone(), Ordering::Equal);
        }
        let k = self.mant().len();
        super::with_limbs(k, |m| {
            m.copy_from_slice(self.mant());
            Float::finish_round(self.neg, m, false, self.exp.clone(), prec, rnd)
        })
    }

    /// Rounds the top-aligned `m` of `±0.m × 2^exp` (with set bits below it if
    /// `sticky`) and builds the result.
    pub(super) fn finish_round(
        neg: bool,
        m: &mut [Limb],
        sticky: bool,
        exp: Integer,
        prec: u64,
        rnd: Round,
    ) -> (Float, Ordering) {
        let r = round_mant(m, sticky, neg, prec, rnd);
        let exp = if r.carry { exp + Integer::ONE } else { exp };
        (Float::from_aligned(neg, &m[r.low..], exp), r.ord)
    }
}
