//! The exact sign of a short sum of floats without forming the sum, so that
//! comparisons such as `|x - a|` against `r` cost the operands' mantissa
//! lengths whatever their exponents (`2^(10^30) - 1` is never built).

use std::cmp::Ordering;

use super::Float;
use crate::Integer;

/// The sign of `Σ ±t` over `terms` (`true` negates), exactly.
///
/// A non-zero term `t` lies in `[2^(lo t), 2^(hi t))` in magnitude and is a
/// multiple of `2^(lo t)`, with `hi t = exp(t)` and `lo t = exp(t) -
/// bits(t)`. Sorted by `hi`, the terms split into clusters: a term joins the
/// current one when its `hi` is at least `L - g`, `L` the lowest `lo` so
/// far and `g = bits(k) + 2` for `k` terms. A cluster's span is then at most
/// the sum of its mantissa lengths and `k g`, so its exact sum is short; if
/// it is not zero it is a non-zero multiple of `2^L`, and every later term is
/// below `2^(L - g)`, their sum below `k 2^(L - g) <= 2^(L - 2)`: the first
/// cluster with a non-zero sum decides the sign.
pub(crate) fn sum_sign(terms: &[(bool, &Float)]) -> Ordering {
    let mut ts: Vec<(Float, Integer)> = terms
        .iter()
        .filter(|(_, t)| !t.is_zero())
        .map(|&(neg, t)| {
            let lo = t.exponent() - &Integer::from(t.bit_length());
            (if neg { t.neg() } else { t.clone() }, lo)
        })
        .collect();
    ts.sort_by(|a, b| b.0.exponent().cmp(a.0.exponent()));
    let guard = Integer::from(u64::from(u64::BITS - ts.len().leading_zeros()) + 2);
    let mut i = 0;
    while i < ts.len() {
        let (mut acc, mut low) = (ts[i].0.clone(), ts[i].1.clone());
        i += 1;
        while i < ts.len() && *ts[i].0.exponent() >= &low - &guard {
            acc = acc.add(&ts[i].0);
            low = std::cmp::min(low, ts[i].1.clone());
            i += 1;
        }
        if !acc.is_zero() {
            return acc.sign();
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use proptest::prelude::*;

    use super::sum_sign;
    use crate::Integer;
    use crate::ball::Float;
    use crate::testkit::float::dyadic_with;

    proptest! {
        #[test]
        fn matches_exact_sums(ts in proptest::collection::vec((any::<bool>(), dyadic_with(2, -200..=200)), 1..6)) {
            let floats: Vec<(bool, Float)> = ts.iter().map(|(n, d)| (*n, d.to_ours())).collect();
            let refs: Vec<(bool, &Float)> = floats.iter().map(|(n, f)| (*n, f)).collect();
            let exact = floats.iter().fold(Float::ZERO, |acc, (n, f)| if *n { acc.sub(f) } else { acc.add(f) });
            prop_assert_eq!(sum_sign(&refs), exact.sign());
        }
    }

    #[test]
    fn huge_gaps_and_cancellation() {
        let p = |e: i64| Float::from(1i64).mul_2exp(&Integer::from(e));
        let big = Float::from(1i64).mul_2exp(&Integer::from(10u64).pow(30));
        let one = Float::from(1i64);
        // 2^(10^30) - 1 > 0, 1 - 2^(10^30) < 0, and exact cancellation of
        // the top cluster hands the decision to the next one.
        assert_eq!(sum_sign(&[(false, &big), (true, &one)]), Ordering::Greater);
        assert_eq!(sum_sign(&[(true, &big), (false, &one)]), Ordering::Less);
        assert_eq!(
            sum_sign(&[(false, &big), (true, &big), (true, &one)]),
            Ordering::Less
        );
        assert_eq!(sum_sign(&[(false, &big), (true, &big)]), Ordering::Equal);
        assert_eq!(
            sum_sign(&[(false, &p(-5)), (true, &p(-5)), (false, &p(-3000))]),
            Ordering::Greater
        );
        assert_eq!(sum_sign(&[]), Ordering::Equal);
    }
}
