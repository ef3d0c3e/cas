//! The constants and their cache against MPFR (CAS-113): from a cleared
//! cache, a sequence of requests for π, log 2, e, γ and G at fuzzer-chosen
//! precisions, which grow, reuse and round the cached balls. Each result
//! meets MPFR's bracket at `prec + 64` bits, knows at least `prec - 1` bits
//! and has a midpoint of at most `prec` bits, the public entries' contract
//! in the constants' tests.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::ball::{Ball, clear_constant_cache};
use cas_num_fuzz::ball::{Oracle, bracket, meets};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;
use rug::float::{Constant, Round};

type Ours = fn(u64) -> Ball;
type Mpfr = fn(u32, Round) -> MpFloat;

const CONSTANTS: [(&str, Ours, Mpfr); 5] = [
    ("pi", Ball::const_pi, |p, r| {
        MpFloat::with_val_round(p, Constant::Pi, r).0
    }),
    ("log2", Ball::const_log2, |p, r| {
        MpFloat::with_val_round(p, Constant::Log2, r).0
    }),
    ("e", Ball::const_e, |p, r| {
        MpFloat::with_val_round(p, MpFloat::with_val(1, 1u32).exp_ref(), r).0
    }),
    ("euler", Ball::const_euler, |p, r| {
        MpFloat::with_val_round(p, Constant::Euler, r).0
    }),
    ("catalan", Ball::const_catalan, |p, r| {
        MpFloat::with_val_round(p, Constant::Catalan, r).0
    }),
];

fuzz_target!(|data: &[u8]| {
    clear_constant_cache();
    let _ = run(&mut Unstructured::new(data));
});

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    for _ in 0..u.int_in_range(1u8..=6)? {
        let (name, f, mp) = CONSTANTS[u.choose_index(CONSTANTS.len())?];
        let prec: u64 = match u.int_in_range(0u8..=7)? {
            0 => u.int_in_range(4_097..=20_000)?,
            1 | 2 => u.int_in_range(257..=4_096)?,
            _ => u.int_in_range(1..=256)?,
        };
        let b = f(prec);
        match bracket(prec + 64, mp) {
            Oracle::Bracket(lo, hi) => {
                assert!(
                    meets(&b, &lo, &hi),
                    "{name} at {prec} bits: {b:?} misses [{lo:?}, {hi:?}]"
                );
            }
            Oracle::Nan | Oracle::Infinite => unreachable!("MPFR's constants are finite"),
        }
        let want = i64::try_from(prec).expect("small precision") - 1;
        let bits = b.rel_accuracy_bits();
        assert!(bits >= want, "{name} at {prec} bits: {bits} bits known");
        assert!(
            b.mid().bit_length() <= prec,
            "{name} at {prec} bits: midpoint of {} bits",
            b.mid().bit_length()
        );
    }
    Ok(())
}
