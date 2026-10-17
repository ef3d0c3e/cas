//! Word-level greatest common divisors.

use crate::Limb;

/// `gcd(a, b)` by the binary (Stein) algorithm; `gcd(0, b) = b`.
#[must_use]
pub fn gcd_1(mut a: Limb, mut b: Limb) -> Limb {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        // a is odd here; removing factors of two from b keeps the gcd.
        b >>= b.trailing_zeros();
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        b -= a;
        if b == 0 {
            return a << shift;
        }
    }
}
