//! Roots, modular arithmetic, CRT, Jacobi symbols, primality, factorisation
//! and the previous prime against GMP.
#![no_main]

use arbitrary::Unstructured;
use cas_num::Natural;
use cas_num::modular::{ModRing, ModWord, crt};
use cas_num_fuzz::{check, gmp, gmp_int, integer, natural};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok((a, ag)), Ok((m, mg))) = (
        u.int_in_range(0u8..=9),
        natural(&mut u, 200),
        natural(&mut u, 160),
    ) else {
        return;
    };
    match op {
        0 => {
            let (s, r) = a.sqrt_rem();
            let (sg, rg) = ag.sqrt_rem(rug::Integer::new());
            check("sqrt_rem", (gmp(&s), gmp(&r)), (sg, rg));
        }
        1 => {
            // pow_mod through every reduction: ModWord, Montgomery on both
            // sides of the REDC threshold, and even moduli.
            if m.is_zero() {
                return;
            }
            let (e, eg) = natural(&mut u, 2).unwrap_or_default();
            check(
                "pow_mod",
                gmp(&a.pow_mod(&e, &m)),
                ag.pow_mod(&eg, &mg).unwrap(),
            );
        }
        2 => {
            if m.is_zero() {
                return;
            }
            let ring = ModRing::new(m.clone());
            let (b, bg) = natural(&mut u, 200).unwrap_or_default();
            let (x, y) = (ring.elem(&a), ring.elem(&b));
            check(
                "ring mul",
                gmp(&ring.to_natural(&ring.mul(&x, &y))),
                (ag.clone() * &bg).modulo(&mg),
            );
            check(
                "ring add",
                gmp(&ring.to_natural(&ring.add(&x, &y))),
                (ag.clone() + &bg).modulo(&mg),
            );
            check(
                "ring sub",
                gmp(&ring.to_natural(&ring.sub(&x, &y))),
                (ag.clone() - &bg).modulo(&mg),
            );
            check(
                "ring inv",
                ring.inv(&x).map(|i| gmp(&ring.to_natural(&i))),
                ag.invert(&mg).ok(),
            );
        }
        3 => {
            let (Some(&n), Ok(x)) = (m.limbs().first().filter(|&&n| n != 0), u.arbitrary::<u64>())
            else {
                return;
            };
            let w = ModWord::new(n);
            let (x, y) = (x % n, w.reduce_natural(&a));
            check(
                "word mul",
                w.mul(x, y),
                (rug::Integer::from(x) * y % n).to_u64().unwrap(),
            );
            check(
                "word inv",
                w.inv(y),
                rug::Integer::from(y)
                    .invert(&rug::Integer::from(n))
                    .ok()
                    .map(|i| i.to_u64().unwrap()),
            );
        }
        4 => {
            let (Ok((b, bg)), Ok((n, ng))) = (natural(&mut u, 20), natural(&mut u, 20)) else {
                return;
            };
            if m.is_zero() || n.is_zero() {
                return;
            }
            let g = mg.clone().gcd(&ng);
            let compatible = (ag.clone() - &bg).is_divisible(&g);
            match crt(&a, &m, &b, &n) {
                None => assert!(!compatible, "crt rejected compatible congruences"),
                Some((x, l)) => {
                    let (x, l) = (gmp(&x), gmp(&l));
                    check("crt lcm", l.clone(), mg.clone().lcm(&ng));
                    assert!(compatible && x < l, "crt solution out of range");
                    check("crt mod m", x.clone().modulo(&mg), ag.modulo(&mg));
                    check("crt mod n", x.modulo(&ng), bg.modulo(&ng));
                }
            }
        }
        5 => {
            let (Ok((x, xg)), true) = (integer(&mut u, 100), m.bit(0)) else {
                return;
            };
            let n = cas_num::Integer::from(m);
            check("jacobi", x.jacobi(&n), xg.jacobi(&gmp_int(&n)));
        }
        6 => {
            if a.limb_len() > 12 {
                return;
            }
            check(
                "primality",
                a.primality().is_probable_prime(),
                ag.is_probably_prime(30) != rug::integer::IsPrime::No,
            );
        }
        7 => {
            if a.is_zero() || a.limb_len() > 3 {
                return;
            }
            let f = a.factor_with_effort(1 << 12);
            check("factor product", f.value(), a.clone());
            for x in f.factors() {
                let prime = gmp(x.base()).is_probably_prime(30) != rug::integer::IsPrime::No;
                check("factor label", x.primality().is_probable_prime(), prime);
            }
        }
        9 => {
            // The previous prime (CAS-267): a prime below `a` with none
            // between, by GMP's next prime after it.
            if a.limb_len() > 3 {
                return;
            }
            match a.previous_prime() {
                None => assert!(a < Natural::from(3u64), "no previous prime"),
                Some(p) => {
                    let pg = gmp(&p);
                    assert!(pg < ag, "previous prime not below");
                    let prime = pg.is_probably_prime(30) != rug::integer::IsPrime::No;
                    check("previous prime", prime, true);
                    assert!(pg.next_prime() >= ag, "a prime between");
                }
            }
        }
        _ => {
            if let Some((b, k)) = a.perfect_power() {
                check("perfect power", b.pow(k), a.clone());
                assert!(ag.is_perfect_power(), "not a perfect power");
            } else {
                check(
                    "perfect power",
                    ag.is_perfect_power() && a > Natural::one(),
                    false,
                );
            }
        }
    }
});
