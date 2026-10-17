//! Radix conversion: printing against GMP, round trips, and parsing of
//! arbitrary strings against the primitive integers' rules.
#![no_main]

use arbitrary::Unstructured;
use cas_num::{Integer, Natural, Rational};
use cas_num_fuzz::{check, gmp_int, integer, rational};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let (Ok(op), Ok(radix)) = (u.int_in_range(0u8..=3), u.int_in_range(2u32..=36)) else {
        return;
    };
    match op {
        0 => {
            let Ok((a, ag)) = integer(&mut u, 300) else {
                return;
            };
            let s = a.to_string_radix(radix);
            check(
                "to_string_radix",
                s.clone(),
                ag.to_string_radix(i32::try_from(radix).unwrap()),
            );
            check(
                "round trip",
                Integer::from_str_radix(&s, radix).ok(),
                Some(a),
            );
        }
        1 => {
            // Short strings over the radix alphabet with signs: our parser
            // follows i128/u128::from_str_radix exactly where those fit.
            let Ok(bytes) = u.arbitrary::<Vec<u8>>() else {
                return;
            };
            let s: String = bytes
                .iter()
                .take(30)
                .map(|&b| {
                    char::from(b"+-0123456789abcdefghijklmnopqrstuvwxyzAZ_ "[usize::from(b) % 42])
                })
                .collect();
            // Valid: an optional sign, then one or more digits below the
            // radix (no '-' for naturals), as for the primitive integers.
            let valid = |signed: bool| {
                let body = s
                    .strip_prefix('+')
                    .or_else(|| s.strip_prefix('-').filter(|_| signed))
                    .unwrap_or(&s);
                !body.is_empty() && body.chars().all(|c| c.is_digit(radix))
            };
            let value = || {
                let body = s.strip_prefix('+').unwrap_or(&s);
                rug::Integer::from_str_radix(body, i32::try_from(radix).unwrap()).unwrap()
            };
            let ours_int = Integer::from_str_radix(&s, radix).ok().map(|x| gmp_int(&x));
            check("parse Integer", ours_int.clone(), valid(true).then(value));
            if let Ok(x) = i128::from_str_radix(&s, radix) {
                check("parse i128", ours_int, Some(rug::Integer::from(x)));
            }
            let ours_nat = Natural::from_str_radix(&s, radix)
                .ok()
                .map(|x| cas_num_fuzz::gmp(&x));
            check("parse Natural", ours_nat.clone(), valid(false).then(value));
            if let Ok(x) = u128::from_str_radix(&s, radix) {
                check("parse u128", ours_nat, Some(rug::Integer::from(x)));
            }
        }
        2 => {
            // Long digit strings: GMP parses them the same way.
            let Ok(n) = u.int_in_range(1usize..=3000) else {
                return;
            };
            let digits: String = (0..n)
                .map(|_| {
                    let d = u.int_in_range(0..=radix - 1).unwrap_or(0);
                    char::from_digit(d, radix).unwrap()
                })
                .collect();
            check(
                "parse long",
                Natural::from_str_radix(&digits, radix)
                    .ok()
                    .map(|x| cas_num_fuzz::gmp(&x)),
                rug::Integer::from_str_radix(&digits, i32::try_from(radix).unwrap()).ok(),
            );
        }
        _ => {
            let Ok(Some((a, _))) = rational(&mut u, 40) else {
                return;
            };
            let s = a.to_string_radix(radix);
            check(
                "rational round trip",
                Rational::from_str_radix(&s, radix).ok(),
                Some(a),
            );
        }
    }
});
