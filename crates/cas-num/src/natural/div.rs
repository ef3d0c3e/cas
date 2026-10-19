//! Division, exact division and divisibility.

use std::ops::{Div, DivAssign, Rem, RemAssign};

use super::Natural;
use super::limb_vec::LimbVec;
use crate::{Limb, kernel};

const DIV_BY_ZERO: &str = "attempt to divide by zero";
const NOT_EXACT: &str = "division is not exact";

/// Scratch limbs that Hensel division takes from the stack instead of the
/// heap, so that short exact divisions do not allocate for it.
const STACK_SCRATCH: usize = 32;

/// Runs `f` on `len` zeroed limbs, on the stack up to [`STACK_SCRATCH`].
fn with_scratch<R>(len: usize, f: impl FnOnce(&mut [Limb]) -> R) -> R {
    if len <= STACK_SCRATCH {
        f(&mut [0; STACK_SCRATCH][..len])
    } else {
        f(&mut vec![0; len])
    }
}

/// Runs `f` on `a` and `d` shifted right by `twos` bits, at most `d`'s
/// trailing zeros (so `d` stays non-zero), both normalised; when `twos > 0`,
/// into one buffer from [`with_scratch`].
fn without_twos<R>(a: &[Limb], d: &[Limb], twos: u64, f: impl FnOnce(&[Limb], &[Limb]) -> R) -> R {
    if twos == 0 {
        return f(a, d);
    }
    let skip = usize::try_from(twos / u64::from(Limb::BITS)).expect("a limb count");
    #[expect(clippy::cast_possible_truncation, reason = "below Limb::BITS")]
    let bits = (twos % u64::from(Limb::BITS)) as u32;
    let (a, d) = (&a[skip.min(a.len())..], &d[skip..]);
    with_scratch(a.len() + d.len(), |buf| {
        let (sa, sd) = buf.split_at_mut(a.len());
        kernel::shr(sa, a, bits);
        kernel::shr(sd, d, bits);
        let (an, dn) = (kernel::normalized_len(sa), kernel::normalized_len(sd));
        f(&sa[..an], &sd[..dn])
    })
}

impl Natural {
    /// `(self / d, self % d)`, or `None` if `d` is zero.
    #[must_use]
    pub fn checked_div_rem(&self, d: &Natural) -> Option<(Natural, Natural)> {
        let (num, den) = (self.limbs(), d.limbs());
        match den {
            [] => None,
            [den] => {
                let (q, r) = self.div_rem_limb(*den);
                Some((q, Natural::from(r)))
            }
            _ if num.len() < den.len() || self < d => Some((Natural::ZERO, self.clone())),
            _ => {
                let mut q = LimbVec::zeroed(num.len() - den.len() + 1);
                let mut r = LimbVec::zeroed(den.len());
                let mut scratch = vec![0; kernel::div_rem_scratch_len(num.len(), den.len())];
                kernel::div_rem(q.as_mut_slice(), r.as_mut_slice(), num, den, &mut scratch);
                let (mut q, mut r) = (Natural { limbs: q }, Natural { limbs: r });
                q.normalize();
                r.normalize();
                q.debug_check();
                r.debug_check();
                Some((q, r))
            }
        }
    }

    /// `(self / d, self % d)`.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem(&self, d: &Natural) -> (Natural, Natural) {
        self.checked_div_rem(d).expect(DIV_BY_ZERO)
    }

    /// `(self / d, self % d)` for a single limb `d`.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn div_rem_limb(&self, d: Limb) -> (Natural, Limb) {
        assert!(d != 0, "{DIV_BY_ZERO}");
        let a = self.limbs();
        let mut q = LimbVec::zeroed(a.len());
        let r = kernel::div_rem_1(q.as_mut_slice(), a, d);
        let mut q = Natural { limbs: q };
        q.normalize();
        q.debug_check();
        (q, r)
    }

    /// `self % d` for a single limb `d`, without computing the quotient.
    ///
    /// # Panics
    ///
    /// If `d` is zero.
    #[must_use]
    pub fn rem_limb(&self, d: Limb) -> Limb {
        assert!(d != 0, "{DIV_BY_ZERO}");
        kernel::rem_1(self.limbs(), d)
    }

    /// `self / d` where `d` is known to divide `self`.
    ///
    /// Short divisors use Hensel (2-adic) division from the low end, which
    /// needs no quotient estimates and whose remainder proves the division
    /// exact; longer ones divide with remainder, which is cheaper than Hensel
    /// division plus a check.
    ///
    /// # Panics
    ///
    /// If `d` is zero or does not divide `self`.
    #[must_use]
    pub fn div_exact(&self, d: &Natural) -> Natural {
        if let [d] = d.limbs() {
            return self.div_exact_limb(*d);
        }
        assert!(!d.is_zero(), "{DIV_BY_ZERO}");
        self.hensel_div(d).expect(NOT_EXACT)
    }

    /// `self / d` for a non-zero `d` that divides `self` by construction (a
    /// gcd, a found factor): [`Natural::div_exact`] without its check, by the
    /// 2-adic quotient alone. Checked in debug builds; otherwise an inexact
    /// division gives an unspecified value.
    pub(crate) fn div_exact_unchecked(&self, d: &Natural) -> Natural {
        let q = match d.limbs() {
            [] => panic!("{DIV_BY_ZERO}"),
            [d] => return self.div_exact_limb(*d),
            _ => self.bdiv_q(d),
        };
        debug_assert!(&q * d == *self, "{NOT_EXACT}");
        q
    }

    /// The 2-adic quotient of `self` and the multi-limb `d`, both without
    /// `d`'s factors of two, to `an - dn + 1` limbs: `self / d` when it is
    /// exact.
    fn bdiv_q(&self, d: &Natural) -> Natural {
        let twos = d.trailing_zeros().expect("non-zero divisor");
        without_twos(self.limbs(), d.limbs(), twos, |a, d| {
            let (an, dn) = (a.len(), d.len());
            if an < dn {
                return Natural::ZERO;
            }
            let mut q = LimbVec::zeroed(an - dn + 1);
            if let [d] = d {
                kernel::divexact_1(q.as_mut_slice(), a, *d);
            } else {
                with_scratch(kernel::bdiv_q_scratch_len(an - dn + 1, dn), |s| {
                    kernel::bdiv_q(q.as_mut_slice(), a, d, s);
                });
            }
            let mut q = Natural { limbs: q };
            q.normalize();
            q.debug_check();
            q
        })
    }

    /// `Some(self / d)` if the non-zero `d` divides `self`, by Hensel
    /// division of both without their common factors of two.
    fn hensel_div(&self, d: &Natural) -> Option<Natural> {
        let Some(s) = self.trailing_zeros() else {
            return Some(Natural::ZERO);
        };
        let twos = d.trailing_zeros().expect("non-zero divisor");
        if s < twos {
            return None;
        }
        without_twos(self.limbs(), d.limbs(), twos, |a, d| {
            let (an, dn) = (a.len(), d.len());
            if an < dn {
                return None;
            }
            let mut q = LimbVec::zeroed(an - dn + 1);
            let exact = if let [d] = d {
                kernel::divexact_1(q.as_mut_slice(), a, *d)
            } else {
                with_scratch(kernel::divexact_scratch_len(an, dn), |s| {
                    kernel::divexact(q.as_mut_slice(), a, d, s)
                })
            };
            exact.then(|| {
                let mut q = Natural { limbs: q };
                q.normalize();
                q.debug_check();
                q
            })
        })
    }

    /// `self / d` for a single limb `d` known to divide `self`.
    ///
    /// # Panics
    ///
    /// If `d` is zero or does not divide `self`.
    #[must_use]
    pub fn div_exact_limb(&self, d: Limb) -> Natural {
        assert!(d != 0, "{DIV_BY_ZERO}");
        let a = self.limbs();
        let mut q = LimbVec::zeroed(a.len());
        assert!(kernel::divexact_1(q.as_mut_slice(), a, d), "{NOT_EXACT}");
        let mut q = Natural { limbs: q };
        q.normalize();
        q.debug_check();
        q
    }

    /// Whether `d` divides `self`. Zero divides only zero.
    #[must_use]
    pub fn is_divisible_by(&self, d: &Natural) -> bool {
        match d.limbs() {
            [] => self.is_zero(),
            [d] => self.is_divisible_by_limb(*d),
            _ => self.hensel_div(d).is_some(),
        }
    }

    /// Whether the single limb `d` divides `self`. Zero divides only zero.
    #[must_use]
    pub fn is_divisible_by_limb(&self, d: Limb) -> bool {
        if d == 0 {
            self.is_zero()
        } else {
            kernel::rem_1(self.limbs(), d) == 0
        }
    }
}

/// Implements the owned/borrowed operator combinations of `/` or `%` from a
/// function of `(&Natural, &Natural)`.
macro_rules! div_ops {
    ($($op:ident::$f:ident, $assign:ident::$fa:ident, $pick:expr;)*) => {$(
        impl $op<&Natural> for &Natural {
            type Output = Natural;

            /// # Panics
            ///
            /// If `rhs` is zero.
            fn $f(self, rhs: &Natural) -> Natural {
                $pick(self.div_rem(rhs))
            }
        }

        impl $op<Natural> for &Natural {
            type Output = Natural;

            fn $f(self, rhs: Natural) -> Natural {
                $op::$f(self, &rhs)
            }
        }

        impl $op<&Natural> for Natural {
            type Output = Natural;

            fn $f(self, rhs: &Natural) -> Natural {
                $op::$f(&self, rhs)
            }
        }

        impl $op<Natural> for Natural {
            type Output = Natural;

            fn $f(self, rhs: Natural) -> Natural {
                $op::$f(&self, &rhs)
            }
        }

        impl $assign<&Natural> for Natural {
            fn $fa(&mut self, rhs: &Natural) {
                *self = $op::$f(&*self, rhs);
            }
        }

        impl $assign<Natural> for Natural {
            fn $fa(&mut self, rhs: Natural) {
                *self = $op::$f(&*self, &rhs);
            }
        }
    )*};
}

div_ops! {
    Div::div, DivAssign::div_assign, |(q, _)| q;
    Rem::rem, RemAssign::rem_assign, |(_, r)| r;
}

impl Div<Limb> for &Natural {
    type Output = Natural;

    /// # Panics
    ///
    /// If `rhs` is zero.
    fn div(self, rhs: Limb) -> Natural {
        self.div_rem_limb(rhs).0
    }
}

impl Div<Limb> for Natural {
    type Output = Natural;

    fn div(self, rhs: Limb) -> Natural {
        &self / rhs
    }
}

impl DivAssign<Limb> for Natural {
    fn div_assign(&mut self, rhs: Limb) {
        *self = &*self / rhs;
    }
}

/// The remainder is below `rhs`, so it is returned as a limb.
impl Rem<Limb> for &Natural {
    type Output = Limb;

    /// # Panics
    ///
    /// If `rhs` is zero.
    fn rem(self, rhs: Limb) -> Limb {
        self.rem_limb(rhs)
    }
}

impl Rem<Limb> for Natural {
    type Output = Limb;

    fn rem(self, rhs: Limb) -> Limb {
        self.rem_limb(rhs)
    }
}

impl RemAssign<Limb> for Natural {
    fn rem_assign(&mut self, rhs: Limb) {
        *self = Natural::from(self.rem_limb(rhs));
    }
}
