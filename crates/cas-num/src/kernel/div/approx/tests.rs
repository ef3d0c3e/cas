//! The approximate kernels against exact quotients and roots from GMP: every
//! result must be a lower bound within its documented error, through the
//! dispatching functions and through the Newton paths at every size.

use proptest::prelude::*;
use rug::Integer;

use super::{
    RSQRT_NEWTON_THRESHOLD, div_approx, div_approx_error, div_approx_newton,
    div_approx_newton_error, div_approx_newton_scratch_len, div_approx_scratch_len, recip_approx,
    recip_approx_error, recip_approx_scratch_len, rsqrt_approx, rsqrt_approx_error,
    rsqrt_approx_scratch_len, sqrt_approx, sqrt_approx_error, sqrt_approx_newton,
    sqrt_approx_newton_error, sqrt_approx_newton_scratch_len, sqrt_approx_scratch_len,
};
use crate::Limb;
use crate::kernel::div::mu::INV_NEWTON_THRESHOLD;
use crate::kernel::mul::mullo::MULLO_DC_THRESHOLD;
use crate::testkit::strategy::{self, len_around};
use crate::testkit::{Mag, ToRug, tested};

/// Filler for output and scratch buffers, so results can't rely on zeroed memory.
const GARBAGE: Limb = 0xdead_beef_dead_beef;

fn val(x: &[Limb]) -> Integer {
    Mag::from_limbs(x).to_rug()
}

fn pow_b(k: usize) -> Integer {
    Integer::from(1) << (64 * u32::try_from(k).unwrap())
}

/// Fails unless `0 <= exact - r <= bound`.
fn check_lower(
    what: &str,
    n: usize,
    exact: &Integer,
    r: &Integer,
    bound: usize,
) -> Result<(), TestCaseError> {
    let err = Integer::from(exact - r);
    prop_assert!(
        err >= 0,
        "{}, n = {}: above the exact value by {}",
        what,
        n,
        -err
    );
    prop_assert!(
        err <= bound,
        "{}, n = {}: error {} beyond the bound {}",
        what,
        n,
        err,
        bound
    );
    Ok(())
}

/// Lengths around the thresholds the kernels recurse through.
fn lens() -> BoxedStrategy<usize> {
    prop_oneof![
        3 => 1..=24usize,
        1 => len_around(&[
            RSQRT_NEWTON_THRESHOLD,
            2 * RSQRT_NEWTON_THRESHOLD,
            INV_NEWTON_THRESHOLD,
            2 * INV_NEWTON_THRESHOLD,
            MULLO_DC_THRESHOLD,
            2 * MULLO_DC_THRESHOLD,
        ]),
        1 => 25..=160usize,
    ]
    .boxed()
}

/// `n` limbs with the top bit set (a normalised divisor).
fn normalised(n: usize) -> BoxedStrategy<Vec<Limb>> {
    strategy::raw_limbs(n)
        .prop_map(|mut v| {
            *v.last_mut().unwrap() |= 1 << (Limb::BITS - 1);
            v
        })
        .boxed()
}

/// `n` limbs normalised for square roots: at least `B^n / 4`, with the edges
/// `B^n / 4` (the root is exactly `2 B^n`) and `B^n - 1` (the root barely
/// exceeds `B^n`) weighted in.
fn sqrt_normalised(n: usize) -> BoxedStrategy<Vec<Limb>> {
    let mut quarter = vec![0; n];
    quarter[n - 1] = 1 << (Limb::BITS - 2);
    prop_oneof![
        1 => Just(quarter),
        1 => Just(vec![Limb::MAX; n]),
        6 => strategy::raw_limbs(n).prop_map(|mut v| {
            let top = v.last_mut().unwrap();
            if *top >> (Limb::BITS - 2) == 0 {
                *top |= 1 << (Limb::BITS - 2);
            }
            v
        }),
    ]
    .boxed()
}

fn check_recip(d: &[Limb]) -> Result<(), TestCaseError> {
    let n = d.len();
    let mut x = vec![GARBAGE; n + 1];
    let mut scratch = vec![GARBAGE; recip_approx_scratch_len(n)];
    recip_approx(&mut x, d, &mut scratch);
    let x = val(&x);
    prop_assert!(x >= pow_b(n) && x < pow_b(n) * 2u32);
    check_lower(
        "recip",
        n,
        &(pow_b(2 * n) / val(d)),
        &x,
        recip_approx_error(n),
    )
}

fn check_div(a: &[Limb], d: &[Limb], newton: bool) -> Result<(), TestCaseError> {
    let n = d.len();
    let mut q = vec![GARBAGE; n + 1];
    let bound = if newton {
        let mut scratch = vec![GARBAGE; div_approx_newton_scratch_len(n)];
        div_approx_newton(&mut q, a, d, &mut scratch);
        div_approx_newton_error(n)
    } else {
        let mut scratch = vec![GARBAGE; div_approx_scratch_len(n)];
        div_approx(&mut q, a, d, &mut scratch);
        div_approx_error(n)
    };
    let exact = val(a) * pow_b(n) / val(d);
    check_lower("div", n, &exact, &val(&q), bound)
}

fn check_rsqrt(a: &[Limb]) -> Result<(), TestCaseError> {
    let n = a.len();
    let mut y = vec![GARBAGE; n + 1];
    let mut scratch = vec![GARBAGE; rsqrt_approx_scratch_len(n)];
    rsqrt_approx(&mut y, a, &mut scratch);
    let y = val(&y);
    prop_assert!(y >= pow_b(n) && y <= pow_b(n) * 2u32);
    let exact = (pow_b(3 * n) / val(a)).sqrt();
    check_lower("rsqrt", n, &exact, &y, rsqrt_approx_error(n))
}

fn check_sqrt(a: &[Limb], newton: bool) -> Result<(), TestCaseError> {
    let n = a.len();
    let mut s = vec![GARBAGE; n];
    let bound = if newton {
        let mut scratch = vec![GARBAGE; sqrt_approx_newton_scratch_len(n)];
        sqrt_approx_newton(&mut s, a, &mut scratch);
        sqrt_approx_newton_error(n)
    } else {
        let mut scratch = vec![GARBAGE; sqrt_approx_scratch_len(n)];
        sqrt_approx(&mut s, a, &mut scratch);
        sqrt_approx_error(n)
    };
    let exact = (val(a) * pow_b(n)).sqrt();
    check_lower("sqrt", n, &exact, &val(&s), bound)
}

proptest! {
    #[test]
    fn recip_within_bound(d in lens().prop_flat_map(normalised)) {
        check_recip(&d)?;
    }

    #[test]
    fn div_within_bound(
        (a, d) in lens().prop_flat_map(|n| (strategy::raw_limbs(n), normalised(n))),
    ) {
        check_div(&a, &d, false)?;
        check_div(&a, &d, true)?;
    }

    #[test]
    fn rsqrt_within_bound(a in lens().prop_flat_map(sqrt_normalised)) {
        check_rsqrt(&a)?;
    }

    #[test]
    fn sqrt_within_bound(a in lens().prop_flat_map(sqrt_normalised)) {
        check_sqrt(&a, false)?;
        check_sqrt(&a, true)?;
    }
}

/// All-ones operands and the normalisation edges, through the Newton paths,
/// at lengths where every level of the recursions changes shape.
#[test]
fn worst_cases() {
    let big = [
        MULLO_DC_THRESHOLD,
        MULLO_DC_THRESHOLD + 1,
        3 * MULLO_DC_THRESHOLD,
    ]
    .map(tested);
    for n in (1..=40).chain(big) {
        let ones = vec![Limb::MAX; n];
        let mut half = vec![0; n];
        half[n - 1] = 1 << (Limb::BITS - 1);
        let mut quarter = vec![0; n];
        quarter[n - 1] = 1 << (Limb::BITS - 2);
        check_recip(&ones).unwrap();
        check_recip(&half).unwrap();
        check_div(&ones, &half, true).unwrap();
        check_div(&ones, &ones, true).unwrap();
        check_div(&ones, &half, false).unwrap();
        for a in [&ones, &half, &quarter] {
            check_rsqrt(a).unwrap();
            check_sqrt(a, true).unwrap();
            check_sqrt(a, false).unwrap();
        }
    }
}

/// The exact dispatch and the error functions' shape.
#[test]
fn exact_cases_and_bounds() {
    assert_eq!(rsqrt_approx_error(RSQRT_NEWTON_THRESHOLD - 1), 0);
    assert_eq!(rsqrt_approx_error(RSQRT_NEWTON_THRESHOLD), 1);
    assert_eq!(recip_approx_error(INV_NEWTON_THRESHOLD - 1), 1);
    assert_eq!(recip_approx_error(INV_NEWTON_THRESHOLD), 2);
    assert_eq!(div_approx_error(1), 0);
    assert_eq!(sqrt_approx_error(1), 0);
    for n in 1..=400 {
        assert!(div_approx_newton_error(n) <= 2 + n.max(1000));
        assert!(sqrt_approx_newton_error(n) <= 2 + n.max(1000));
    }
}

/// Prints the time of the exact and the Newton square root and quotient by
/// size, to place `SQRT_NEWTON_THRESHOLD` (and check `MU_DIV_THRESHOLD` for
/// the approximate quotient). Run with
/// `cargo test -p cas-num --release -- --ignored --nocapture tune_approx`.
#[test]
#[ignore = "timing, not a test"]
fn tune_approx() {
    use std::hint::black_box;
    use std::time::Instant;

    fn time(mut f: impl FnMut()) -> f64 {
        f();
        let start = Instant::now();
        let mut reps = 0u32;
        while start.elapsed().as_millis() < 100 {
            f();
            reps += 1;
        }
        start.elapsed().as_secs_f64() * 1e6 / f64::from(reps)
    }

    println!(
        "{:>6} {:>12} {:>12} {:>8} {:>12} {:>12} {:>8}",
        "limbs", "sqrt exact", "sqrt newton", "ratio", "div exact", "div newton", "ratio"
    );
    for k in 2..=13 {
        let n = 1usize << k;
        let mut a: Vec<Limb> = (0..n as u64)
            .map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0x5555)
            .collect();
        a[n - 1] |= 1 << 63;
        let mut d = a.iter().rev().copied().collect::<Vec<_>>();
        d[n - 1] |= 1 << 63;
        let mut s = vec![0; n];
        let shift = n as u64 * 64;
        let exact = time(|| {
            black_box((crate::Natural::from_limbs(&a) << shift).isqrt());
        });
        let mut scratch = vec![0; sqrt_approx_newton_scratch_len(n)];
        let newton = time(|| sqrt_approx_newton(&mut s, black_box(&a), &mut scratch));
        let mut q = vec![0; n + 1];
        let mut num = vec![0; 2 * n];
        num[n..].copy_from_slice(&a);
        let mut r = vec![0; n];
        let mut dscratch = vec![0; crate::kernel::div_rem_scratch_len(2 * n, n)];
        let dexact =
            time(|| crate::kernel::div_rem(&mut q, &mut r, black_box(&num), &d, &mut dscratch));
        let mut nscratch = vec![0; div_approx_newton_scratch_len(n)];
        let dnewton = time(|| div_approx_newton(&mut q, black_box(&a), &d, &mut nscratch));
        println!(
            "{n:>6} {exact:>12.2} {newton:>12.2} {:>8.2} {dexact:>12.2} {dnewton:>12.2} {:>8.2}",
            newton / exact,
            dnewton / dexact
        );
    }
}
