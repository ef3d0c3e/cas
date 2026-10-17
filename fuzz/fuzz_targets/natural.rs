//! `Natural` arithmetic, bit operations and identities against GMP.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Natural;
use cas_num_fuzz::{check, gmp, natural};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok((a, ag)), Ok((b, bg))) = (
        u.int_in_range(0u8..=11),
        natural(&mut u, 400),
        natural(&mut u, 400),
    ) else {
        return;
    };
    match op {
        0 => check("add", gmp(&(&a + &b)), ag + bg),
        1 => check(
            "checked_sub",
            a.checked_sub(&b).map(|d| gmp(&d)),
            (ag >= bg).then(|| ag - bg),
        ),
        2 => check("mul", gmp(&(&a * &b)), ag * bg),
        3 => check("square", gmp(&a.square()), ag.square()),
        4 => {
            if b.is_zero() {
                return;
            }
            let (q, r) = a.div_rem(&b);
            let (qg, rg) = ag.div_rem(bg);
            check("div_rem", (gmp(&q), gmp(&r)), (qg, rg));
        }
        5 => {
            // (a / b) b + a % b == a, through the operators.
            if b.is_zero() {
                return;
            }
            check("division identity", &(&a / &b) * &b + &(&a % &b), a.clone());
        }
        6 => check("gcd", gmp(&a.gcd(&b)), ag.gcd(&bg)),
        7 => check("lcm", gmp(&a.lcm(&b)), ag.lcm(&bg)),
        8 => {
            let k = b.limbs().first().copied().unwrap_or(0) % 2000;
            check(
                "shl",
                gmp(&(&a << k)),
                ag.clone() << u32::try_from(k).unwrap(),
            );
            check("shr", gmp(&(&a >> k)), ag >> u32::try_from(k).unwrap());
        }
        9 => {
            check("and", gmp(&(&a & &b)), ag.clone() & &bg);
            check("or", gmp(&(&a | &b)), ag.clone() | &bg);
            check("xor", gmp(&(&a ^ &b)), ag ^ bg);
        }
        10 => {
            check("cmp", a.cmp(&b), ag.cmp(&bg));
            check(
                "bit_length",
                a.bit_length(),
                u64::from(ag.significant_bits()),
            );
            check(
                "count_ones",
                a.count_ones(),
                u64::from(ag.count_ones().unwrap()),
            );
        }
        _ => {
            // Exact division of a product.
            if b.is_zero() {
                return;
            }
            let p = &a * &b;
            check("div_exact", p.div_exact(&b), a.clone());
            check(
                "is_divisible_by",
                (&p + 1u64).is_divisible_by(&b),
                b == Natural::one(),
            );
        }
    }
});
