//! `Rational` arithmetic, comparison and rounding against GMP.
#![no_main]

use arbitrary::Unstructured;
use cas_num_fuzz::{check, gmp_int, gmp_rat, rational};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok(Some((a, ag))), Ok(Some((b, bg)))) = (
        u.int_in_range(0u8..=5),
        rational(&mut u, 60),
        rational(&mut u, 60),
    ) else {
        return;
    };
    match op {
        0 => {
            check("add", gmp_rat(&(&a + &b)), ag.clone() + &bg);
            check("sub", gmp_rat(&(&a - &b)), ag - bg);
        }
        1 => check("mul", gmp_rat(&(&a * &b)), ag * bg),
        2 => check(
            "checked_div",
            a.checked_div(&b).map(|q| gmp_rat(&q)),
            (bg != 0).then(|| ag / bg),
        ),
        3 => check("cmp", a.cmp(&b), ag.cmp(&bg)),
        4 => {
            check(
                "floor",
                gmp_int(&a.floor()),
                ag.clone().floor().into_numer_denom().0,
            );
            check(
                "ceil",
                gmp_int(&a.ceil()),
                ag.clone().ceil().into_numer_denom().0,
            );
            check(
                "trunc",
                gmp_int(&a.trunc()),
                ag.clone().trunc().into_numer_denom().0,
            );
            check(
                "round",
                gmp_int(&a.round()),
                ag.round().into_numer_denom().0,
            );
        }
        _ => {
            let e = i64::from(u.int_in_range(-20i8..=20).unwrap_or(0));
            if a.is_zero() && e < 0 {
                return;
            }
            check(
                "pow",
                gmp_rat(&a.pow(e)),
                rug::ops::Pow::pow(ag, i32::try_from(e).unwrap()),
            );
        }
    }
});
