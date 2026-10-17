//! `Integer` arithmetic, rounding divisions, extended GCD and powers
//! against GMP.
#![no_main]

use arbitrary::Unstructured;
use cas_num::{Integer, Rounding};
use cas_num_fuzz::{check, gmp_int, integer};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok((a, ag)), Ok((b, bg))) = (
        u.int_in_range(0u8..=8),
        integer(&mut u, 200),
        integer(&mut u, 200),
    ) else {
        return;
    };
    match op {
        0 => {
            check("add", gmp_int(&(&a + &b)), ag.clone() + &bg);
            check("sub", gmp_int(&(&a - &b)), ag - bg);
        }
        1 => check("mul", gmp_int(&(&a * &b)), ag * bg),
        2 => {
            if b.is_zero() {
                return;
            }
            for (rounding, gmp) in [
                (Rounding::Trunc, ag.clone().div_rem(bg.clone())),
                (Rounding::Floor, ag.clone().div_rem_floor(bg.clone())),
                (Rounding::Ceil, ag.clone().div_rem_ceil(bg.clone())),
                (Rounding::Euclid, ag.clone().div_rem_euc(bg.clone())),
            ] {
                let (q, r) = a.checked_div_rem(&b, rounding).unwrap();
                check("div_rem", (gmp_int(&q), gmp_int(&r)), gmp);
            }
        }
        3 => {
            let (g, s, t) = a.gcd_ext(&b);
            let (gg, sg, tg) = ag.extended_gcd(bg, rug::Integer::new());
            check(
                "gcd_ext",
                (gmp_int(&g), gmp_int(&s), gmp_int(&t)),
                (gg, sg, tg),
            );
        }
        4 => check(
            "mod_inverse",
            a.mod_inverse(&b).map(|x| gmp_int(&x)),
            if b.is_zero() {
                None
            } else {
                ag.invert(&bg.abs()).ok()
            },
        ),
        5 => {
            let e = u64::from(u.int_in_range(0u16..=300).unwrap_or(0));
            if a.bit_length() * e > 1 << 20 {
                return;
            }
            check(
                "pow",
                gmp_int(&a.pow(e)),
                rug::ops::Pow::pow(ag, u32::try_from(e).unwrap()),
            );
        }
        6 => {
            if b.is_zero() || b.bit_length() > 64 * 16 {
                return;
            }
            let e = Integer::from(u.int_in_range(-1000i64..=100_000).unwrap_or(0));
            let eg = rug::Integer::from(e.to_i64().unwrap());
            check(
                "pow_mod",
                a.pow_mod(&e, &b).map(|x| gmp_int(&x)),
                ag.pow_mod(&eg, &bg.abs()).ok(),
            );
        }
        7 => {
            let n = u.int_in_range(1u32..=12).unwrap_or(1);
            if a.is_negative() && n % 2 == 0 {
                return;
            }
            let (r, rem) = a.root_rem(n);
            let (rg, remg) = ag.root_rem(rug::Integer::new(), n);
            check("root_rem", (gmp_int(&r), gmp_int(&rem)), (rg, remg));
        }
        _ => {
            check("cmp", a.cmp(&b), ag.cmp(&bg));
            check(
                "is_perfect_power",
                a.is_perfect_power(),
                ag.is_perfect_power(),
            );
            check(
                "is_perfect_square",
                a.is_perfect_square(),
                ag.is_perfect_square(),
            );
        }
    }
});
