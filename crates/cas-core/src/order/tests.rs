#![allow(clippy::many_single_char_names, reason = "x, y, f mirror the notation")]

use std::cmp::Ordering;
use std::time::{Duration, Instant};

use proptest::prelude::*;

use cas_num::ball::{Ball, ComplexBall, Float, Mag};

use crate::{Complex, Expr, Integer, Rational, Symbol};

fn sym(name: &str) -> Expr {
    Symbol::new(name).into()
}

fn int(x: i64) -> Expr {
    x.into()
}

fn q(n: i64, d: i64) -> Expr {
    Rational::new(n.into(), d.into()).into()
}

fn big(x: i128) -> Expr {
    Integer::from(x).into()
}

/// The complex number `re + im·I` with integer parts.
fn cplx(re: i64, im: i64) -> Expr {
    Complex::new(re.into(), im.into()).into()
}

fn apply(head: &Expr, args: &[Expr]) -> Expr {
    Expr::apply(head, args)
}

fn times(args: &[Expr]) -> Expr {
    apply(&Symbol::TIMES.into(), args)
}

fn pow(b: &Expr, n: Expr) -> Expr {
    apply(&Symbol::POWER.into(), &[b.clone(), n])
}

/// Each element is less than the next.
fn assert_increasing(v: &[Expr]) {
    for w in v.windows(2) {
        assert_eq!(w[0].cmp(&w[1]), Ordering::Less, "{:?} < {:?}", w[0], w[1]);
        assert_eq!(
            w[1].cmp(&w[0]),
            Ordering::Greater,
            "{:?} > {:?}",
            w[1],
            w[0]
        );
    }
}

#[test]
fn numbers_by_value_and_first() {
    assert_increasing(&[int(i64::MIN), int(-5), int(0), int(3), int(i64::MAX)]);
    assert_increasing(&[
        int(i64::MAX),
        times(&[]),
        Expr::string(""),
        sym("a"),
        apply(&sym("a"), &[]),
    ]);
}

/// Integers inside and outside `i64` and rationals interleave by value.
#[test]
fn all_numbers_by_value() {
    let min = i128::from(i64::MIN);
    let max = i128::from(i64::MAX);
    assert_increasing(&[
        big(min * 4),
        big(min - 1),
        int(i64::MIN),
        q(-5, 2),
        int(-2),
        q(-1, 3),
        int(0),
        q(1, 1_000_000),
        q(1, 2),
        int(1),
        q(4, 3),
        int(i64::MAX),
        Expr::from(Rational::new((max * 2 + 1).into(), 2.into())),
        big(max + 1),
        big(max * 4),
        times(&[]),
    ]);
    // Complex numbers after every real one, by real part, then imaginary.
    assert_increasing(&[
        big(max * 4),
        cplx(-2, 5),
        cplx(0, -1),
        cplx(0, 1),
        Expr::from(Complex::new(q(1, 2).to_rational().unwrap(), 7.into())),
        cplx(1, -3),
        cplx(1, 2),
        times(&[]),
        sym("I"),
    ]);
    // A rational coefficient orders like an integer one.
    let x = sym("x");
    assert_increasing(&[
        times(&[int(-1), x.clone()]),
        times(&[q(1, 2), x.clone()]),
        x.clone(),
        times(&[q(3, 2), x.clone()]),
        times(&[int(2), x]),
    ]);
}

#[test]
fn kinds_then_names() {
    assert_increasing(&[
        Expr::string("b"),
        sym("a"),
        sym("A"),
        sym("ab"),
        sym("Ab"),
        sym("b"),
        sym("x10"),
        sym("x2"),
        apply(&sym("a"), &[]),
    ]);
    assert_increasing(&[Expr::string("a"), Expr::string("A"), Expr::string("b")]);
    // Non-ASCII: folded case first, then lowercase first, then code point.
    assert_increasing(&[sym("z"), sym("ß"), sym("é"), sym("É")]);
    // A name the language escaped sorts as its spelling, just before it.
    assert_increasing(&[
        sym("a"),
        sym("$A"),
        sym("A"),
        sym("ab"),
        sym("$Ab"),
        sym("Ab"),
        sym("b"),
        sym("$Delta"),
        sym("e"),
        sym("$E"),
        sym("E"),
        sym("f"),
    ]);
    // Strings keep their `$`.
    assert_increasing(&[Expr::string("$b"), Expr::string("a")]);
}

#[test]
fn applications_by_head_then_length_then_arguments() {
    let (f, g, a, b) = (sym("f"), sym("g"), sym("a"), sym("b"));
    assert_increasing(&[
        apply(&f, std::slice::from_ref(&b)),
        apply(&f, &[a.clone(), b.clone()]),
        apply(&f, &[b.clone(), a.clone()]),
        apply(&g, std::slice::from_ref(&a)),
        apply(&apply(&f, &[]), &[]),
    ]);
    // Cos[x] < f[x] < Sin[x]: case-insensitive names.
    let x = sym("x");
    assert_increasing(&[
        apply(&Symbol::COS.into(), std::slice::from_ref(&x)),
        apply(&f, std::slice::from_ref(&x)),
        apply(&Symbol::SIN.into(), std::slice::from_ref(&x)),
    ]);
}

#[test]
fn polynomial_terms() {
    let (x, y) = (sym("x"), sym("y"));
    // (1 + x + y)^2 = 1 + 2 x + x^2 + 2 y + 2 x y + y^2
    assert_increasing(&[
        int(1),
        times(&[int(2), x.clone()]),
        pow(&x, int(2)),
        times(&[int(2), y.clone()]),
        times(&[int(2), x.clone(), y.clone()]),
        pow(&y, int(2)),
    ]);
    // (x + y)^3 = x^3 + 3 x^2 y + 3 x y^2 + y^3
    assert_increasing(&[
        pow(&x, int(3)),
        times(&[int(3), pow(&x, int(2)), y.clone()]),
        times(&[int(3), x.clone(), pow(&y, int(2))]),
        pow(&y, int(3)),
    ]);
    // Coefficients, then exponents, then other functions.
    let sin_x = apply(&Symbol::SIN.into(), std::slice::from_ref(&x));
    assert_increasing(&[
        times(&[int(-1), x.clone()]),
        x.clone(),
        times(&[int(2), x.clone()]),
        pow(&x, int(2)),
        pow(&x, sym("n")),
        sin_x.clone(),
        times(&[x.clone(), sin_x]),
    ]);
    // A numeric base is a number: 2^x < y.
    assert_increasing(&[pow(&int(2), x.clone()), y.clone()]);
}

#[test]
fn equal_only_when_identical() {
    let x = sym("x");
    // Equal as products of powers, ordered structurally.
    assert_increasing(&[
        x.clone(),
        pow(&x, int(1)),
        times(std::slice::from_ref(&x)),
        times(&[int(1), x.clone()]),
    ]);
    assert_eq!(
        times(std::slice::from_ref(&x)).cmp(&times(&[x])),
        Ordering::Equal
    );
}

/// Ids follow creation order; the canonical order must not.
#[test]
fn independent_of_creation_order() {
    let zeta = sym("order_creation_zeta");
    let alpha = sym("order_creation_alpha");
    assert!(alpha < zeta);
    assert!(apply(&alpha, &[]) < apply(&zeta, &[]));
}

/// A ball far from an exact number orders by the exponents, at once: the
/// midpoints near `2^(±10^9)` that a power of a ball makes at a point of
/// `Trace::verify` took a second each to make rationals (CAS-169). Near
/// one, by value, an exact number before a ball of its value.
#[test]
fn far_balls_order_by_exponents() {
    let ball = |n: i64, e: i64| {
        let mid = Float::from(n).mul_2exp(&Integer::from(e));
        Expr::from(Ball::new(mid, Mag::ZERO))
    };
    let start = Instant::now();
    assert_increasing(&[
        ball(-3, 1_000_000_000),
        int(-5),
        q(-1, 3),
        ball(-3, -1_000_000_000),
        int(0),
        ball(3, -1_000_000_000),
        q(1, 1_000_000),
        int(7),
        big(3 << 70),
        ball(3, 70),
        ball(3, 1_000_000_000),
    ]);
    let took = start.elapsed();
    assert!(took < Duration::from_secs(5), "{took:?}");
}

/// Balls of small dyadic midpoints, so that they meet the exact numbers
/// above, exact or with a radius.
fn small_ball() -> impl Strategy<Value = Ball> {
    (-12i64..12, 0i64..3, prop::option::of(-4i64..0)).prop_map(|(n, e, r)| {
        let mid = Float::from(n).mul_2exp(&Integer::from(-e));
        Ball::new(
            mid,
            r.map_or(Mag::ZERO, |r| Mag::power_of_two(&Integer::from(r))),
        )
    })
}

/// A tree over a few names, biased towards `Times` and `Power`.
fn tree() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "A", "$A", "b", "f", "Times", "Power"]).prop_map(sym),
        (-2i64..3).prop_map(int),
        (-3i64..4, 1i64..4).prop_map(|(n, d)| q(n, d)),
        prop::sample::select(vec![-(1i128 << 64), 1 << 63, 3 << 70]).prop_map(big),
        (-2i64..3, prop_oneof![-2i64..0, 1i64..3]).prop_map(|(re, im)| cplx(re, im)),
        prop::sample::select(vec!["", "a", "B"]).prop_map(Expr::string),
        small_ball().prop_map(Expr::from),
        (small_ball(), small_ball()).prop_map(|(re, im)| Expr::from(ComplexBall::new(re, im))),
    ];
    leaf.prop_recursive(4, 24, 3, |inner| {
        let head = prop_oneof![
            Just(Expr::from(Symbol::TIMES)),
            Just(Expr::from(Symbol::POWER)),
            Just(sym("f")),
            inner.clone(),
        ];
        (head, prop::collection::vec(inner, 0..4)).prop_map(|(h, a)| Expr::apply(&h, &a))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    #[test]
    fn antisymmetric_and_equal_iff_identical(a in tree(), b in tree()) {
        prop_assert_eq!(a.cmp(&b), b.cmp(&a).reverse());
        prop_assert_eq!(a.cmp(&b) == Ordering::Equal, a == b);
        prop_assert_eq!(a.cmp(&a), Ordering::Equal);
    }

    /// Sorting checks consistency itself; every pair of the result must then
    /// agree with its positions, which fails if the order is not transitive.
    #[test]
    fn transitive(mut v in prop::collection::vec(tree(), 0..8)) {
        v.sort();
        for i in 0..v.len() {
            for j in i..v.len() {
                prop_assert_ne!(v[i].cmp(&v[j]), Ordering::Greater, "{:?} > {:?}", v[i], v[j]);
            }
        }
        for a in &v {
            for b in &v {
                for c in &v {
                    if a <= b && b <= c {
                        prop_assert!(a <= c, "{:?} <= {:?} <= {:?}", a, b, c);
                    }
                }
            }
        }
    }

    /// The coefficient decides only between terms with the same factors.
    #[test]
    fn coefficients_order_like_terms(m in tree(), c in (-3i64..4, 1i64..4), d in (-3i64..4, 1i64..4)) {
        prop_assume!(!m.is_number());
        let (tc, td) = (times(&[q(c.0, c.1), m.clone()]), times(&[q(d.0, d.1), m.clone()]));
        let value = |(n, d): (i64, i64)| Rational::new(n.into(), d.into());
        prop_assert_eq!(tc.cmp(&td), value(c).cmp(&value(d)));
    }

    /// A ball sorts by its midpoint among exact numbers, after an exact
    /// number of that value, and balls of one midpoint by radius
    /// (CAS-DOC-5 §9.1).
    #[test]
    fn balls_sort_by_midpoint(b in small_ball(), c in small_ball(), n in -12i64..12, d in 1i64..5) {
        let x = Rational::new(n.into(), d.into());
        let mid = b.mid().to_rational();
        let expected = mid.cmp(&x).then(Ordering::Greater);
        prop_assert_eq!(Expr::from(b.clone()).cmp(&Expr::from(&x)), expected);
        let expected = b.mid().cmp(c.mid()).then_with(|| b.rad().cmp(c.rad()));
        prop_assert_eq!(Expr::from(b).cmp(&Expr::from(c)), expected);
    }

    /// Numbers compare by value, whatever their kind.
    #[test]
    fn numbers_compare_by_value(a in (-5i128..5, 1i128..5, 0u32..3), b in (-5i128..5, 1i128..5, 0u32..3)) {
        let value = |(n, d, k): (i128, i128, u32)| {
            let scale = (0..k).fold(Integer::ONE, |acc, _| acc * Integer::from(i64::MAX));
            Rational::new(Integer::from(n) * scale, Integer::from(d))
        };
        let (x, y) = (value(a), value(b));
        prop_assert_eq!(Expr::from(&x).cmp(&Expr::from(&y)), x.cmp(&y));
    }
}
