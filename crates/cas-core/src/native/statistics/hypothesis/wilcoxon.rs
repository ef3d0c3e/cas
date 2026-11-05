//! Wilcoxon's signed-rank and Mann–Whitney's rank-sum statistics
//! (CAS-DOC-52 §5, CAS-507): `WilcoxonDistribution` (`wilcoxonp`) and
//! `WilcoxonStatistic` (`wilcoxons`).
//!
//! - `WilcoxonDistribution[n]`: `{P(W⁺ = 0), …, P(W⁺ = n(n+1)/2)}`, the
//!   subsets of `{1, …, n}` by their sums over `2ⁿ`, by the subset-sum
//!   recurrence; `n ≤ 200`.
//! - `WilcoxonDistribution[m, n]`: `{P(U = 0), …, P(U = mn)}`, the
//!   coefficients of the Gaussian binomial `[m+n, m]_q` over `C(m+n, m)`,
//!   by `∏ᵢ (1 − q^(n+i))/(1 − q^i)` with the smaller of `m`, `n` as `m`;
//!   `mn ≤ 40000`. giac answers both as a `poly1` of the coefficients.
//! - `WilcoxonStatistic[l, μ]`: `W⁺`, the sum of the ranks of `|xᵢ − μ|`
//!   over `xᵢ > μ`, the differences of 0 dropped and ties at their
//!   midranks, where giac keeps the zeros as positive and ranks ties apart
//!   (`wilcoxons([1,2,9,4],5)` is 4 there, 7/2 here).
//! - `WilcoxonStatistic[l₁, l₂]`: the sum of `l₂`'s ranks in the pooled
//!   sample, at midranks, giac's (as a decimal there).
//!
//! `WilcoxonTest` (`wilcoxont`, CAS-508), with the alternatives and level
//! of `ZTest`:
//!
//! - `WilcoxonTest[l, μ, alt, α]`: `W⁺` against the signed-rank law of the
//!   `k` nonzero differences, `P(W⁺ ≤ w)` for `"<"`, `P(W⁺ ≥ w)` for `">"`,
//!   and twice the smaller, at most 1, for `"!="`.
//! - `WilcoxonTest[l₁, l₂, alt, α]`: `U₁ = R₁ − n₁(n₁+1)/2` from `l₁`'s
//!   rank sum against Mann–Whitney's law, `P(U ≥ U₁)` for `">"` (`l₁` the
//!   larger), `P(U ≤ U₁)` for `"<"`; giac's one-sided alternatives are the
//!   other way round and strict (`wilcoxont([1,3,4,20],[2,6,10],'>')` has
//!   its p-value 11/35, `P(U ≤ 4)`; 24/35 here).
//! - With ties, the laws of tie-free data, the usual convention; the
//!   statistic may then be a half-integer.

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::{Alternative, Native, compare, decide, options, written};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::rounding::real_sign;
use crate::native::statistics::{is_list, refuse};
use crate::symbol::Symbol;

/// The largest `n` of the signed-rank law.
pub(super) const MAX_SIGNED: u64 = 200;
/// The largest `mn` of the rank-sum law.
pub(super) const MAX_PRODUCT: u64 = 40_000;

/// The numbers of subsets of `{1, …, n}` by their sums, `0 … n(n+1)/2`.
pub(super) fn signed_rank_counts(n: u64) -> Vec<Integer> {
    let top = usize::try_from(n * (n + 1) / 2).expect("small");
    let mut c = vec![Integer::ZERO; top + 1];
    c[0] = Integer::ONE;
    let mut reach = 0;
    for k in 1..=usize::try_from(n).expect("small") {
        reach += k;
        for s in (k..=reach).rev() {
            let below = c[s - k].clone();
            c[s] += &below;
        }
    }
    c
}

/// The numbers of the subsets of size `m` of `{1, …, m + n}` by their rank
/// sums less `m(m+1)/2`, `0 … mn`: the Gaussian binomial's coefficients.
pub(super) fn rank_sum_counts(m: u64, n: u64) -> Vec<Integer> {
    let (m, n) = (m.min(n), m.max(n));
    let (m, n) = (
        usize::try_from(m).expect("small"),
        usize::try_from(n).expect("small"),
    );
    // Room for the products before their divisions.
    let mut p = vec![Integer::ZERO; m * n + n + m + 1];
    p[0] = Integer::ONE;
    let mut degree = 0;
    for i in 1..=m {
        // Times 1 − q^(n+i).
        let shift = n + i;
        for k in (shift..=degree + shift).rev() {
            let low = p[k - shift].clone();
            p[k] -= &low;
        }
        degree += shift;
        // Over 1 − q^i: b[k] = a[k] + b[k − i], exact.
        for k in i..=degree {
            let low = p[k - i].clone();
            p[k] += &low;
        }
        degree -= i;
    }
    p.truncate(m * n + 1);
    p
}

/// The counts over their sum, as a list.
fn law(counts: &[Integer]) -> Expr {
    let total = counts.iter().fold(Integer::ZERO, |mut t, c| {
        t += c;
        t
    });
    Expr::call(
        Symbol::LIST,
        counts
            .iter()
            .map(|c| Expr::from(Rational::new(c.clone(), total.clone())))
            .collect::<Vec<_>>(),
    )
}

/// A natural number argument up to `max`.
fn natural(x: &Expr, max: u64) -> Option<u64> {
    x.as_i64()
        .and_then(|k| u64::try_from(k).ok())
        .filter(|&k| k <= max)
}

/// `WilcoxonDistribution[n]`, `WilcoxonDistribution[m, n]`, giac's
/// `wilcoxonp`.
pub(in crate::native) fn wilcoxon_distribution(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let result = match e.args() {
        [n] => {
            let Some(n) = natural(n, u64::MAX) else {
                return Ok(None);
            };
            if n > MAX_SIGNED {
                return refuse(ev, e, "the signed-rank law is built up to n = 200");
            }
            law(&signed_rank_counts(n))
        }
        [m, n] => {
            let (Some(m), Some(n)) = (natural(m, u64::MAX), natural(n, u64::MAX)) else {
                return Ok(None);
            };
            if m.checked_mul(n).is_none_or(|p| p > MAX_PRODUCT) {
                return refuse(ev, e, "the rank-sum law is built up to m n = 40000");
            }
            law(&rank_sum_counts(m, n))
        }
        _ => return Ok(None),
    };
    ev.step("wilcoxonp", e, &result);
    Ok(Some(result))
}

/// The midranks of the real `values`, 1 for the least; `None` when
/// comparisons do not order them.
pub(super) fn midranks(
    ev: &mut Eval<'_>,
    values: &[Expr],
) -> Result<Option<Vec<Rational>>, EvalError> {
    // Insertion sort of the indices, by decided comparisons.
    let mut order: Vec<usize> = Vec::with_capacity(values.len());
    for i in 0..values.len() {
        let mut at = order.len();
        while at > 0 {
            match compare(ev, &values[order[at - 1]], &values[i])? {
                Some(Ordering::Greater) => at -= 1,
                Some(_) => break,
                None => return Ok(None),
            }
        }
        order.insert(at, i);
    }
    let mut ranks = vec![Rational::ZERO; values.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() {
            match compare(ev, &values[order[end]], &values[order[start]])? {
                Some(Ordering::Equal) => end += 1,
                Some(_) => break,
                None => return Ok(None),
            }
        }
        // Ranks start + 1 … end: their mean.
        let mid = Rational::new(
            Integer::from(u64::try_from(start + 1 + end).expect("small")),
            Integer::from(2),
        );
        for &i in &order[start..end] {
            ranks[i] = mid.clone();
        }
        start = end;
    }
    Ok(Some(ranks))
}

/// `W⁺` of `l` about `mu`, and the number of nonzero differences; `None`
/// when they are not real or not ordered.
pub(super) fn signed_rank(
    ev: &mut Eval<'_>,
    l: &[Expr],
    mu: &Expr,
) -> Result<Option<(Rational, usize)>, EvalError> {
    let mut positive = Vec::new();
    let mut sizes = Vec::new();
    for x in l {
        let d = ev.eval(&(x.clone() - mu.clone()))?;
        match real_sign(ev.context(), &d) {
            Some(Ordering::Equal) => {}
            Some(s) => {
                positive.push(s == Ordering::Greater);
                sizes.push(ev.eval(&Expr::call(Symbol::ABS, [d]))?);
            }
            None => return Ok(None),
        }
    }
    let Some(ranks) = midranks(ev, &sizes)? else {
        return Ok(None);
    };
    let w = ranks
        .iter()
        .zip(&positive)
        .filter(|&(_, &p)| p)
        .fold(Rational::ZERO, |w, (r, _)| &w + r);
    Ok(Some((w, sizes.len())))
}

/// The rank sum of `b` in the pooled sample of `a` and `b`; `None` when
/// they are not ordered.
pub(super) fn rank_sum(
    ev: &mut Eval<'_>,
    a: &[Expr],
    b: &[Expr],
) -> Result<Option<Rational>, EvalError> {
    let pooled: Vec<Expr> = a.iter().chain(b).cloned().collect();
    let Some(ranks) = midranks(ev, &pooled)? else {
        return Ok(None);
    };
    Ok(Some(
        ranks[a.len()..].iter().fold(Rational::ZERO, |s, r| &s + r),
    ))
}

/// Whether every entry of the list `l` is a real number.
pub(super) fn real_list(ev: &Eval<'_>, l: &Expr) -> bool {
    is_list(l)
        && !l.args().is_empty()
        && l.args()
            .iter()
            .all(|x| real_sign(ev.context(), x).is_some())
}

/// `WilcoxonStatistic[l, μ]`, `WilcoxonStatistic[l₁, l₂]`, giac's
/// `wilcoxons`.
pub(in crate::native) fn wilcoxon_statistic(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, other] = e.args() else {
        return Ok(None);
    };
    if !real_list(ev, l) {
        return Ok(None);
    }
    let found = if is_list(other) {
        if !real_list(ev, other) {
            return Ok(None);
        }
        rank_sum(ev, l.args(), other.args())?
    } else if real_sign(ev.context(), other).is_some() {
        signed_rank(ev, l.args(), other)?.map(|(w, _)| w)
    } else {
        return Ok(None);
    };
    let Some(found) = found else {
        return refuse(ev, e, "the values are not ordered");
    };
    let result = Expr::from(found);
    ev.step("wilcoxons", e, &result);
    Ok(Some(result))
}

/// The p-value of `x` against the law of the `counts`: the lower tail
/// `P(X ≤ x)`, the upper `P(X ≥ x)`, or twice the smaller, at most 1.
fn tail_p(counts: &[Integer], x: &Rational, alternative: Alternative) -> Rational {
    let total = counts.iter().fold(Integer::ZERO, |mut t, c| {
        t += c;
        t
    });
    let mut lower = Integer::ZERO;
    let mut upper = Integer::ZERO;
    for (s, c) in counts.iter().enumerate() {
        let s = Rational::from(Integer::from(u64::try_from(s).expect("small")));
        if s <= *x {
            lower += c;
        }
        if s >= *x {
            upper += c;
        }
    }
    let (lower, upper) = (
        Rational::new(lower, total.clone()),
        Rational::new(upper, total),
    );
    match alternative {
        Alternative::Less => lower,
        Alternative::Greater => upper,
        Alternative::Unequal => {
            let twice = &lower.min(upper) * &Rational::from(Integer::from(2));
            twice.min(Rational::ONE)
        }
    }
}

/// `WilcoxonTest[l, μ, alt, α]`, `WilcoxonTest[l₁, l₂, alt, α]`, giac's
/// `wilcoxont`.
pub(in crate::native) fn wilcoxon_test(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [l, other, rest @ ..] = e.args() else {
        return Ok(None);
    };
    if !real_list(ev, l) {
        return Ok(None);
    }
    let two = is_list(other);
    if two && !real_list(ev, other) || !two && real_sign(ev.context(), other).is_none() {
        return Ok(None);
    }
    let Some((alternative, alpha)) = options(ev, e, rest)? else {
        return Ok(None);
    };
    let (statistic, counts, what, name) = if two {
        let (n1, n2) = (l.args().len(), other.args().len());
        let (m, n) = (
            u64::try_from(n1).expect("small"),
            u64::try_from(n2).expect("small"),
        );
        if m.checked_mul(n).is_none_or(|p| p > MAX_PRODUCT) {
            return refuse(ev, e, "the rank-sum law is built up to m n = 40000");
        }
        let Some(r1) = rank_sum(ev, other.args(), l.args())? else {
            return refuse(ev, e, "the values are not ordered");
        };
        let shift = Rational::new(Integer::from(m * (m + 1)), Integer::from(2));
        (
            &r1 - &shift,
            rank_sum_counts(m, n),
            "Mann–Whitney test",
            "U₁",
        )
    } else {
        let Some((w, k)) = signed_rank(ev, l.args(), other)? else {
            return refuse(ev, e, "the values are not ordered");
        };
        let k = u64::try_from(k).expect("small");
        if k > MAX_SIGNED {
            return refuse(ev, e, "the signed-rank law is built up to n = 200");
        }
        (w, signed_rank_counts(k), "Wilcoxon signed-rank test", "W⁺")
    };
    let p = Expr::from(tail_p(&counts, &statistic, alternative));
    let what = format!("{what}, H1: {}", written(alternative));
    decide(ev, e, &what, &[(name, Expr::from(statistic))], &p, &alpha)
}
