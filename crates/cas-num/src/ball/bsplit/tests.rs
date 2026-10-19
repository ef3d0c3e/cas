//! The binary-splitting sum against the series summed term by term in exact
//! rationals, on random sequences with each factor present or always 1,
//! and with powers of two kept apart by [`Series::shift`].

use proptest::prelude::*;

use super::{Series, quotient, sum, sum_shifted};
use crate::{Integer, Natural, Rational};

/// Tabulated sequences; an empty table is the constant 1.
#[derive(Debug, Clone)]
struct Table {
    p: Vec<i64>,
    q: Vec<i64>,
    a: Vec<i64>,
    b: Vec<i64>,
    /// Shifts; empty is 0.
    s: Vec<u64>,
}

fn entry(v: &[i64], k: u64) -> Option<Integer> {
    (!v.is_empty()).then(|| Integer::from(v[usize::try_from(k).unwrap()]))
}

impl Series for Table {
    fn p(&self, k: u64) -> Option<Integer> {
        entry(&self.p, k)
    }

    fn q(&self, k: u64) -> Integer {
        entry(&self.q, k).unwrap()
    }

    fn a(&self, k: u64) -> Option<Integer> {
        entry(&self.a, k)
    }

    fn b(&self, k: u64) -> Option<Integer> {
        entry(&self.b, k)
    }

    fn shift(&self, k: u64) -> u64 {
        self.s
            .get(usize::try_from(k).unwrap())
            .copied()
            .unwrap_or(0)
    }
}

/// `2^k` as a rational.
fn pow2(k: u64) -> Rational {
    Rational::from(Integer::from_parts(false, Natural::from(1u64) << k))
}

fn direct(s: &Table, n: u64) -> Rational {
    let one = Rational::from(Integer::ONE);
    let (mut sum, mut prod) = (Rational::from(Integer::ZERO), one.clone());
    for k in 0..n {
        let get = |x: Option<Integer>| x.map_or_else(|| one.clone(), Rational::from);
        prod = &(&(&prod * &get(s.p(k))) / &Rational::from(s.q(k))) / &pow2(s.shift(k));
        sum = &sum + &(&(&prod * &get(s.a(k))) / &get(s.b(k)));
    }
    sum
}

fn nonzero() -> impl Strategy<Value = i64> {
    prop_oneof![
        -1000i64..=-1,
        1i64..=1000,
        Just(i64::MAX),
        Just(i64::MIN + 1)
    ]
}

fn table() -> impl Strategy<Value = (Table, u64)> {
    (1usize..40).prop_flat_map(|n| {
        let seq = |nz: bool| {
            let v = if nz {
                nonzero().boxed()
            } else {
                any::<i64>().boxed()
            };
            prop_oneof![Just(Vec::new()), proptest::collection::vec(v, n)]
        };
        (
            seq(false),
            proptest::collection::vec(nonzero(), n),
            seq(false),
            seq(true),
        )
            .prop_map(move |(p, q, a, b)| {
                (
                    Table {
                        p,
                        q,
                        a,
                        b,
                        s: Vec::new(),
                    },
                    n as u64,
                )
            })
    })
}

/// [`table`] with shifts in each term, from 0 to beyond a limb.
fn shifted_table() -> impl Strategy<Value = (Table, u64)> {
    table().prop_flat_map(|(t, n)| {
        let len = usize::try_from(n).unwrap();
        proptest::collection::vec(prop_oneof![Just(0u64), 0u64..8, 60u64..200], len)
            .prop_map(move |s| (Table { s, ..t.clone() }, n))
    })
}

proptest! {
    #[test]
    fn matches_term_by_term((s, n) in table()) {
        let (t, d, q) = sum(&s, n);
        prop_assert!(!d.is_zero());
        let ball = quotient(&t, &d, 64);
        let got = Rational::new(t, d);
        prop_assert_eq!(&got, &direct(&s, n));
        // d = B Q, with Q the product of the q(k).
        let qs = (0..n).fold(Integer::ONE, |acc, k| &acc * &s.q(k));
        prop_assert_eq!(q, qs);
        // The rounded quotient contains it.
        prop_assert!(ball.contains_rational(&got));
    }

    #[test]
    fn shifts_match_term_by_term((s, n) in shifted_table()) {
        let (t, d, e) = sum_shifted(&s, n);
        prop_assert!(!d.is_zero());
        let got = &Rational::new(t.clone(), d.clone()) / &pow2(e);
        prop_assert_eq!(&got, &direct(&s, n));
        prop_assert_eq!(e, (0..n).map(|k| s.shift(k)).sum::<u64>());
        // The scaled quotient contains it.
        let ball = quotient(&t, &d, 64).mul_2exp(&-Integer::from(e));
        prop_assert!(ball.contains_rational(&got));
    }
}

#[test]
#[should_panic(expected = "no terms")]
fn zero_terms() {
    let _ = sum(
        &Table {
            p: vec![],
            q: vec![1],
            a: vec![],
            b: vec![],
            s: vec![],
        },
        0,
    );
}
