//! Ball arithmetic against exact `Rational` evaluation on dyadic inputs
//! (CAS-DOC-5 §11): every result ball must contain the exact value.
//!
//! Skeleton until the ball types land: it decodes the inputs and the
//! working precision, computes the exact reference with our `Rational`, and
//! checks the reference against GMP. CAS-76 adds the ball side at each
//! `CAS-76` mark.
#![no_main]

use arbitrary::Unstructured;
use cas_num_fuzz::{check, dyadic, gmp_rat};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok(_prec), Ok((a, ag)), Ok((b, bg))) = (
        u.int_in_range(0u8..=3),
        u.int_in_range(1u32..=512),
        dyadic(&mut u, 8),
        dyadic(&mut u, 8),
    ) else {
        return;
    };
    let exact = match op {
        0 => Some(&a + &b),
        1 => Some(&a - &b),
        2 => Some(&a * &b),
        _ => a.checked_div(&b),
    };
    let gmp = match op {
        0 => Some(ag + bg),
        1 => Some(ag - bg),
        2 => Some(ag * bg),
        _ => (bg != 0).then(|| ag / bg),
    };
    check("exact reference", exact.as_ref().map(gmp_rat), gmp);
    // CAS-76: build balls from `a`, `b` at `_prec` (exact and with a radius),
    // apply `op`, and check that the result contains `exact` (and that
    // division by a ball containing 0 gives the indeterminate ball).
});
