//! Transposed Vandermonde systems over a field, solved as Zippel does
//! \[Zip90\] (CAS-189): `Σⱼ cⱼ vⱼⁱ = yᵢ` for `i = 1 … t` in `O(t²)`
//! operations. Zippel's sparse interpolation (`zippel.rs`) solves one for
//! each power of the main variable of a gcd's images.

#![cfg_attr(
    not(test),
    expect(dead_code, reason = "used by the sparse interpolation, CAS-190")
)]

use cas_num::algebra::Field;

/// The `c` with `Σⱼ cⱼ vⱼⁱ = yᵢ` for `i = 1 … t`, `t` the length of `v` and
/// of `y`; `None` when two of the `vⱼ` coincide or one is 0, which makes the
/// system singular.
///
/// With `P(Z) = ∏ₗ (Z − vₗ)` and `Pⱼ(Z) = P(Z) / (Z − vⱼ) = Σₖ qⱼₖ Zᵏ`,
/// `Σᵢ yᵢ qⱼ,ᵢ₋₁ = Σₗ cₗ vₗ Pⱼ(vₗ) = cⱼ vⱼ Pⱼ(vⱼ)`, since `Pⱼ` vanishes at
/// the other `vₗ`: one synthetic division, one sum and one evaluation for
/// each `j`.
///
/// # Panics
///
/// If `y` and `v` differ in length.
pub(crate) fn solve<F: Field>(field: &F, v: &[F::Elem], y: &[F::Elem]) -> Option<Vec<F::Elem>> {
    let t = v.len();
    assert_eq!(y.len(), t, "as many values as points");
    // P, from its constant term up: monic of degree t.
    let mut p = vec![field.one()];
    for x in v {
        p.push(field.zero());
        for k in (1..p.len()).rev() {
            p[k] = field.sub(&p[k - 1], &field.mul(x, &p[k]));
        }
        p[0] = field.neg(&field.mul(x, &p[0]));
    }
    let mut q = vec![field.zero(); t];
    let mut c = Vec::with_capacity(t);
    for x in v {
        if field.is_zero(x) {
            return None;
        }
        // P = (Z − x) Q: q_{t−1} = p_t and q_{k−1} = p_k + x q_k.
        q[t - 1] = p[t].clone();
        for k in (1..t).rev() {
            q[k - 1] = field.add(&p[k], &field.mul(x, &q[k]));
        }
        let sum = y.iter().zip(&q).fold(field.zero(), |s, (yi, qk)| {
            field.add(&s, &field.mul(yi, qk))
        });
        let at = q
            .iter()
            .rev()
            .fold(field.zero(), |s, qk| field.add(&field.mul(&s, x), qk));
        // Pⱼ(vⱼ) is 0 when vⱼ is another point too.
        c.push(field.mul(&sum, &field.inv(&field.mul(x, &at))?));
    }
    Some(c)
}

#[cfg(test)]
mod tests {
    use cas_num::algebra::{BigPrimeField, Ring, WordPrimeField};
    use cas_num::{Integer, Natural};
    use proptest::prelude::*;

    use super::*;

    /// `yᵢ = Σⱼ cⱼ vⱼⁱ` for `i = 1 … t`.
    fn values<F: Field>(field: &F, c: &[F::Elem], v: &[F::Elem]) -> Vec<F::Elem> {
        (1..=v.len() as u64)
            .map(|i| {
                c.iter().zip(v).fold(field.zero(), |s, (cj, vj)| {
                    field.add(&s, &field.mul(cj, &field.pow(vj, i)))
                })
            })
            .collect()
    }

    proptest! {
        /// Random systems over word primes, small and large: the solution
        /// gives the values back. Points are distinct and not 0.
        #[test]
        fn solves_random_systems(
            p in prop::sample::select(vec![11u64, 65_537, (1 << 61) - 1, u64::MAX - 58]),
            pairs in prop::collection::btree_map(1u64..1 << 40, any::<u64>(), 0..24),
        ) {
            let f = WordPrimeField::new(p).unwrap();
            let mut v = Vec::new();
            let mut c = Vec::new();
            for (x, cj) in pairs {
                let x = f.elem_integer(&Integer::from(x));
                if !f.is_zero(&x) && !v.contains(&x) {
                    v.push(x);
                    c.push(f.elem_integer(&Integer::from(cj)));
                }
            }
            let y = values(&f, &c, &v);
            prop_assert_eq!(solve(&f, &v, &y), Some(c));
        }
    }

    #[test]
    fn degenerate_systems() {
        let f = WordPrimeField::new(13).unwrap();
        // No unknowns, one, and every non-zero point of 𝔽₁₃.
        assert_eq!(solve(&f, &[], &[]), Some(vec![]));
        assert_eq!(solve(&f, &[5], &[10]), Some(vec![2]));
        let v: Vec<u64> = (1..13).collect();
        let c: Vec<u64> = (0..12).map(|k| (k * 7 + 3) % 13).collect();
        assert_eq!(solve(&f, &v, &values(&f, &c, &v)), Some(c));
        // A repeated point, and 0.
        assert_eq!(solve(&f, &[3, 4, 3], &[1, 2, 3]), None);
        assert_eq!(solve(&f, &[3, 0], &[1, 2]), None);
        // Over a prime past a word.
        let big = BigPrimeField::new(Natural::power_of_two(127) - 1u64).unwrap();
        let int = |x: i64| big.elem_integer(&Integer::from(x));
        let (v, c) = (vec![int(2), int(-3), int(5)], vec![int(7), int(0), int(-1)]);
        assert_eq!(solve(&big, &v, &values(&big, &c, &v)), Some(c));
    }
}
