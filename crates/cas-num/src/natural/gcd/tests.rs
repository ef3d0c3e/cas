#![allow(
    clippy::many_single_char_names,
    reason = "a, b, g, s, t mirror the notation"
)]
//! `Natural` GCD and LCM checked against GMP, and the Lehmer machinery checked
//! on its own.

use proptest::prelude::*;
use rug::Integer;

use super::hgcd::{HgcdRow, classical, hgcd, hgcd_row, hgcd_values, strassen};
use super::{Step, lehmer};
use crate::Natural;
use crate::testkit::strategy::{self, DEFAULT_MAX_LIMBS};
use crate::testkit::{Mag, ToRug, differential};

fn pair() -> BoxedStrategy<(Mag, Mag)> {
    prop_oneof![
        4 => strategy::gcd_pair(DEFAULT_MAX_LIMBS),
        1 => strategy::unbalanced_pair(DEFAULT_MAX_LIMBS),
    ]
    .boxed()
}

/// Matrix and division steps `lehmer` takes on `(a, b)`.
fn steps(a: &Natural, b: &Natural) -> (usize, usize) {
    let (mut x, mut y) = super::ordered(a, b);
    let (mut matrices, mut divisions) = (0, 0);
    lehmer(&mut x, &mut y, |step| match step {
        Step::Matrix(_) => matrices += 1,
        Step::Quotient(_) => divisions += 1,
        Step::Hgcd(_) => unreachable!("half-GCD is off"),
    });
    (matrices, divisions)
}

/// Uniformly random limbs with a top limb of at least `2^32`.
fn dense(len: usize) -> impl Strategy<Value = Natural> {
    proptest::collection::vec(any::<u64>(), len).prop_map(move |mut limbs| {
        limbs[len - 1] |= 1 << 32;
        Natural::from_limbs(&limbs)
    })
}

fn fibonacci(n: u32) -> Natural {
    Mag::from_rug(&Integer::from(Integer::fibonacci(n))).to_natural()
}

/// Pairs for the half-GCD: dense, nearly equal, with a large common factor,
/// or of different lengths.
fn hgcd_pair() -> BoxedStrategy<(Natural, Natural)> {
    let n = 3usize..240;
    prop_oneof![
        3 => n.clone().prop_flat_map(|n| (dense(n), dense(n))),
        1 => n.clone().prop_flat_map(|n| (dense(n), dense(n / 2 + 1))),
        1 => n.clone().prop_flat_map(|n| (dense(n), any::<u64>()))
            .prop_map(|(a, d)| (&a + &Natural::from(d), a)),
        1 => (dense(25), dense(20), dense(20)).prop_map(|(c, a, b)| (&a * &c, &b * &c)),
        1 => (100u32..4000).prop_map(|k| (fibonacci(k + 1), fibonacci(k))),
    ]
    .boxed()
}

differential! {
    fn gcd((a, b) in pair()) {
        ours: {
            let (x, y) = (a.to_natural(), b.to_natural());
            let g = x.gcd(&y);
            prop_assert_eq!(&g, &y.gcd(&x));
            g
        },
        gmp: a.to_rug().gcd(&b.to_rug()),
    }

    fn lcm((a, b) in pair()) {
        ours: a.to_natural().lcm(&b.to_natural()),
        gmp: a.to_rug().lcm(&b.to_rug()),
    }

    /// `g = s a + t b` for an integer `t`, i.e. `b` divides `s a - g`.
    fn gcd_cofactor((a, b) in pair()) {
        ours: {
            let (g, s) = a.to_natural().gcd_cofactor(&b.to_natural());
            let rest = s.to_rug() * a.to_rug() - g.to_rug();
            let b = b.to_rug();
            prop_assert!(
                if b == 0 { rest == 0 } else { rest.is_divisible(&b) },
                "s = {:?} is not a Bézout coefficient", s
            );
            g
        },
        gmp: a.to_rug().gcd(&b.to_rug()),
    }
}

proptest! {
    /// `hgcd` returns `M` and `(x, y)` with `(a, b) = M (x, y)`, `det M` as
    /// recorded, `x >= y`, both values at least `B^s` and no step left
    /// (`|x - y| < B^s`), and `hgcd_values` the same values.
    #[test]
    fn hgcd_contract((a, b) in hgcd_pair()) {
        let (a, b) = if a >= b { (a, b) } else { (b, a) };
        let (mut x, mut y) = (a.limbs().to_vec(), b.limbs().to_vec());
        let (mut vx, mut vy) = (x.clone(), y.clone());
        let progress = hgcd_values(&mut vx, &mut vy);
        let n = x.len().max(y.len());
        let s = n / 2 + 1;
        let floor = Integer::from(1) << (64 * s);
        let (ra, rb) = (a.to_rug(), b.to_rug());
        let m = hgcd(&mut x, &mut y);
        prop_assert_eq!(progress, m.is_some());
        prop_assert!(vx == x && vy == y, "hgcd_values reduced differently");
        match m {
            None => {
                prop_assert!(x == a.limbs() && y == b.limbs());
                prop_assert!(ra.clone().min(rb.clone()) < floor || (ra - rb).abs() < floor);
            }
            Some(m) => {
                let [[u00, u01], [u10, u11]] = m.m.each_ref().map(|r| r.each_ref().map(ToRug::to_rug));
                let (rx, ry) = (Natural::from_limbs(&x).to_rug(), Natural::from_limbs(&y).to_rug());
                prop_assert!(rx >= ry, "hgcd left x < y");
                prop_assert_eq!(&ra, &(u00.clone() * &rx + u01.clone() * &ry));
                prop_assert_eq!(&rb, &(u10.clone() * &rx + u11.clone() * &ry));
                let det = u00 * &u11 - u01 * &u10;
                prop_assert_eq!(det, if m.negative { -1 } else { 1 });
                prop_assert!(rx >= floor && ry >= floor);
                prop_assert!((rx - ry).abs() < floor);
            }
        }
    }

    /// `hgcd_row` takes the steps `hgcd` takes and multiplies the row by
    /// their matrix, keeping the sign of its determinant.
    #[test]
    fn hgcd_row_is_row_times_matrix(
        (a, b) in hgcd_pair(),
        r in [(); 2].map(|()| proptest::collection::vec(any::<u64>(), 0..30)),
    ) {
        let (a, b) = if a >= b { (a, b) } else { (b, a) };
        let (mut x, mut y) = (a.limbs().to_vec(), b.limbs().to_vec());
        let (mut x2, mut y2) = (x.clone(), y.clone());
        let m = hgcd(&mut x, &mut y);
        let r = r.map(|l| Natural::from_limbs(&l));
        let mut row = HgcdRow { r: r.clone(), negative: false };
        prop_assert_eq!(hgcd_row(&mut x2, &mut y2, &mut row), m.is_some());
        prop_assert_eq!((&x2, &y2), (&x, &y));
        if let Some(m) = m {
            let entry = |j: usize| &(&r[0] * &m.m[0][j]) + &(&r[1] * &m.m[1][j]);
            prop_assert_eq!(&row.r, &[entry(0), entry(1)]);
            prop_assert_eq!(row.negative, m.negative);
        } else {
            prop_assert_eq!(&row.r, &r);
        }
    }

    /// On uniformly random operands the half-GCD does nearly half the work.
    #[test]
    fn hgcd_halves((a, b) in (3usize..240).prop_flat_map(|n| (dense(n), dense(n)))) {
        let (a, b) = if a >= b { (a, b) } else { (b, a) };
        let (mut x, mut y) = (a.limbs().to_vec(), b.limbs().to_vec());
        let s = a.limb_len() / 2 + 1;
        prop_assert!(hgcd(&mut x, &mut y).is_some());
        prop_assert!(x.len().max(y.len()) <= s + 2, "{} limbs left, s = {}", x.len().max(y.len()), s);
    }

    /// GCD and cofactor with half-GCD steps from two limbs on.
    #[test]
    fn half_gcd_everywhere((a, b) in prop_oneof![pair().prop_map(|(a, b)| (a.to_natural(), b.to_natural())), hgcd_pair()]) {
        let g = a.to_rug().gcd(&b.to_rug());
        prop_assert!(a.gcd_from(&b, 2).to_rug() == g);
        let (g2, s) = a.gcd_cofactor_from(&b, 2);
        prop_assert!(g2.to_rug() == g);
        let rest = s.to_rug() * a.to_rug() - g;
        let b = b.to_rug();
        prop_assert!(if b == 0 { rest == 0 } else { rest.is_divisible(&b) }, "s = {:?}", s);
    }
}

proptest! {
    /// Half-GCD steps down to `dc` limbs, whose cofactors are tracked by
    /// matrix products, then Lehmer steps on fresh cofactors, which join
    /// them at the end.
    #[test]
    fn fresh_cofactors_after_half_gcd((a, b) in hgcd_pair(), dc in 3usize..60) {
        let g = a.to_rug().gcd(&b.to_rug());
        let (g2, s) = a.gcd_cofactor_from(&b, dc);
        prop_assert!(g2.to_rug() == g);
        let rest = s.to_rug() * a.to_rug() - g;
        prop_assert!(rest.is_divisible(&b.to_rug()), "s = {:?}", s);
    }
}

/// A 2×2 matrix whose entries have independent lengths, some empty.
fn matrix() -> impl Strategy<Value = [[Natural; 2]; 2]> {
    let entry = || {
        prop_oneof![
            1 => Just(Natural::ZERO),
            4 => proptest::collection::vec(any::<u64>(), 1..80)
                .prop_map(|l| Natural::from_limbs(&l)),
        ]
    };
    [[entry(), entry()], [entry(), entry()]]
}

proptest! {
    /// Strassen's seven products give the eight-product result, which is
    /// the product.
    #[test]
    fn matrix_products(a in matrix(), b in matrix()) {
        let c = classical(&a, &b);
        prop_assert_eq!(&strassen(&a, &b), &c);
        let rug = |m: &[[Natural; 2]; 2]| m.each_ref().map(|r| r.each_ref().map(ToRug::to_rug));
        let ([[a00, a01], [a10, a11]], [[b00, b01], [b10, b11]]) = (rug(&a), rug(&b));
        let expect = [
            [a00.clone() * &b00 + a01.clone() * &b10, a00 * &b01 + a01 * &b11],
            [a10.clone() * &b00 + a11.clone() * &b10, a10 * &b01 + a11 * &b11],
        ];
        prop_assert_eq!(rug(&c), expect);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(6))]

    /// Operands of thousands of limbs: the half-GCD recursion at the default
    /// thresholds.
    #[test]
    fn huge_operands((a, b, c) in (dense(2500), dense(2000), dense(300))) {
        let (x, y) = (&a * &c, &b * &c);
        let g = x.gcd(&y);
        prop_assert!(g.to_rug() == x.to_rug().gcd(&y.to_rug()));
        let (g2, s) = x.gcd_cofactor(&y);
        prop_assert!(g2 == g);
        prop_assert!((s.to_rug() * x.to_rug() - g.to_rug()).is_divisible(&y.to_rug()));
    }
}

proptest! {
    /// On dense operands nearly every step is a Lehmer matrix, each about 64
    /// bits of the leading 128; the full division is the exception.
    #[test]
    fn lehmer_does_the_work(a in dense(24), b in dense(24)) {
        let (matrices, divisions) = steps(&a, &b);
        prop_assert!(matrices >= 8 * divisions, "{} matrices, {} divisions", matrices, divisions);
        // 24 limbs at 60 bits a matrix or more.
        prop_assert!((20..=26).contains(&matrices), "{} matrices", matrices);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Operands of hundreds of limbs, where division steps take the
    /// divide-and-conquer path.
    #[test]
    fn large_operands((a, b, c) in (dense(200), dense(150), dense(120))) {
        let (x, y) = (&a * &c, &b * &c);
        let g = x.gcd(&y);
        prop_assert!(g.to_rug() == x.to_rug().gcd(&y.to_rug()));
        let (g2, s) = x.gcd_cofactor(&y);
        prop_assert!(g2 == g);
        prop_assert!((s.to_rug() * x.to_rug() - g.to_rug()).is_divisible(&y.to_rug()));
    }
}

#[test]
fn fibonacci_worst_case_uses_matrices() {
    let (a, b) = (fibonacci(5001), fibonacci(5000));
    assert!(a.gcd(&b).is_one());
    let (matrices, divisions) = steps(&a, &b);
    // 3471 bits: at most 60 matrices means at least 57 bits a matrix even on
    // all-ones quotients.
    assert!(
        (40..=60).contains(&matrices) && divisions <= matrices / 8,
        "{matrices} matrices, {divisions} divisions"
    );
}

#[test]
fn small_cases() {
    let n = |x: u64| Natural::from(x);
    assert_eq!(n(0).gcd(&n(0)), n(0));
    assert_eq!(n(0).gcd(&n(12)), n(12));
    assert_eq!(n(12).gcd(&n(0)), n(12));
    assert_eq!(n(12).gcd(&n(18)), n(6));
    assert_eq!(n(0).lcm(&n(5)), n(0));
    assert_eq!(n(4).lcm(&n(6)), n(12));
    let big = Natural::power_of_two(200);
    assert_eq!(
        big.gcd(&Natural::power_of_two(130)),
        Natural::power_of_two(130)
    );
    assert_eq!(big.gcd(&n(0)), big);
    assert_eq!(n(0).gcd(&big), big);
}
