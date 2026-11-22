//! Real root isolation (CAS-538): the LMQ bounds on the positive roots and
//! Vincent–Akritas–Strzeboński (VAS) continued fractions, as giac's
//! `posubLMQ`, `poslbdLMQ`, `VAS_positive` and `VAS`, over ℚ in exact
//! arithmetic.
//!
//! - [`posub_lmq`]: Akritas's Local-Max-Quadratic upper bound on the
//!   positive roots. With the coefficients `cᵢ` of `xⁱ`, the leading one
//!   made positive, each negative `c_m` is paired with the positive `c_n`,
//!   `n > m`, giving the least `(2^t |c_m| / c_n)^(1/(n−m))`, `t` counting
//!   how often `c_n` has been used, from 1; the bound is the largest of
//!   these, times `65/64` and rounded up, as giac's margin and `ceil`
//!   (`ceil(65·e^max/64)`); 0 without a negative coefficient.
//! - [`poslbd_lmq`]: the same on the reversed polynomial (the constant made
//!   positive), whose bound `B` is on `1/x`, so `1/B` bounds the positive
//!   roots from below; giac rounds it to the power of two `2^(−⌈log₂ B⌉)`,
//!   `+∞` without a negative coefficient, 0 when a negative coefficient has
//!   no positive one below it.
//! - [`vas_positive`]: the positive roots of a square-free polynomial in
//!   disjoint intervals, by Möbius maps `x ↦ (a x + b)/(c x + d)` of the
//!   half line. Descartes's sign variations `v` bound the roots in it;
//!   `v = 0` drops the map, `v = 1` gives the interval between `b/d` and
//!   `a/c` (with `c = 0`, from `b` to `b + a·posubLMQ(f)`). Otherwise the
//!   map is moved by the lower bound when it is at least 1 (scaled first
//!   when above 16), then split at 1: `x ↦ x + 1` and `x ↦ 1/(x + 1)`,
//!   a root at the split point given exactly. The entries are sorted by
//!   their left end, a root before an interval starting there.
//! - [`vas`]: every real root of the square-free part, `p(−x)`'s positive
//!   roots mirrored for the negative ones (or the positive ones mirrored,
//!   when `p` is even or odd), the root 0 between them.
//!
//! The order of the work list, the rounding of the bounds and the tests on
//! them follow giac's `csturm.cc`, so the intervals are giac's: `VAS(x³ −
//! 7x + 7)` is `[−4, 0]`, `[1, 3/2]`, `[3/2, 2]`. giac compares the
//! bounds' logarithms in doubles; here `q₁^(1/k₁)` and `q₂^(1/k₂)` are
//! compared as `q₁^k₂` and `q₂^k₁`, and the roundings are exact, so a
//! value on a rounding edge may land one step from giac's.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: the map's a, b, c, d, the polynomial f, its variations v"
)]

use std::cmp::Ordering;

use cas_num::algebra::{EuclideanDomain, Rationals};
use cas_num::{Integer, Rational};

use crate::dense::{Dense, DenseRing};

/// An entry of an isolation: an exact root, or an interval `[lo, hi]`
/// holding exactly one root in its interior (an end may be another root,
/// given exactly).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Isolated {
    /// A rational root.
    Root(Rational),
    /// The interval from the first to the second.
    Interval(Rational, Rational),
}

/// giac's `poslbdLMQ`: 0, a power of two, or `+∞`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LowerBound {
    /// A negative coefficient with no positive one of lower degree.
    Zero,
    /// `2^e`.
    Power(i64),
    /// No negative coefficient.
    Infinity,
}

impl LowerBound {
    /// The bound as a rational, `None` for `+∞`.
    #[must_use]
    pub fn value(&self) -> Option<Rational> {
        match self {
            LowerBound::Zero => Some(Rational::ZERO),
            LowerBound::Power(e) => Some(two().pow(*e)),
            LowerBound::Infinity => None,
        }
    }
}

fn two() -> Rational {
    Rational::from(Integer::from(2))
}

/// `q^(1/k)`, kept as the pair.
#[derive(Clone)]
struct Radical {
    q: Rational,
    k: u32,
}

impl Radical {
    fn cmp(&self, other: &Radical) -> Ordering {
        self.q
            .pow(i64::from(other.k))
            .cmp(&other.q.pow(i64::from(self.k)))
    }

    /// `(2^t |num| / den)^(1/k)`.
    fn of(num: &Rational, den: &Rational, t: u32, k: usize) -> Radical {
        Radical {
            q: &(&num.abs() / den) * &two().pow(i64::from(t)),
            k: u32::try_from(k).expect("a degree in memory"),
        }
    }
}

/// The coefficients from the leading one down, as giac's `modpoly`.
fn descending(p: &Dense<Rational>) -> Vec<Rational> {
    p.coeffs().iter().rev().cloned().collect()
}

/// giac's `posubLMQ` of a polynomial given from its leading coefficient
/// down.
fn upper(p: &[Rational]) -> Integer {
    let k = p.len();
    if k <= 1 {
        return Integer::from(0);
    }
    let flip = !p[0].is_positive();
    // cl[i] is the coefficient of xⁱ.
    let cl: Vec<Rational> = p
        .iter()
        .rev()
        .map(|c| if flip { -c } else { c.clone() })
        .collect();
    let mut largest: Option<Radical> = None;
    let mut used = vec![1_u32; k + 1];
    for m in (1..k).rev() {
        if !cl[m - 1].is_negative() {
            continue;
        }
        let mut least: Option<Radical> = None;
        for n in (m + 1..=k).rev() {
            if !cl[n - 1].is_positive() {
                continue;
            }
            let t = Radical::of(&cl[m - 1], &cl[n - 1], used[n - 1], n - m);
            used[n - 1] += 1;
            if least
                .as_ref()
                .is_none_or(|l| l.cmp(&t) == Ordering::Greater)
            {
                least = Some(t);
            }
        }
        if let Some(least) = least
            && largest
                .as_ref()
                .is_none_or(|l| least.cmp(l) == Ordering::Greater)
        {
            largest = Some(least);
        }
    }
    let Some(Radical { q, k }) = largest else {
        return Integer::from(0);
    };
    // The least n with n ≥ 65/64 q^(1/k), that is (64 n/65)^k ≥ q; from
    // r = ⌊q^(1/k)⌋, within two steps.
    let ratio = Rational::new(Integer::from(64), Integer::from(65));
    let r = Rational::from(q.floor().nth_root(k));
    let mut n = (&r / &ratio).ceil();
    while (&Rational::from(n.clone()) * &ratio).pow(i64::from(k)) < q {
        n = &n + &Integer::from(1);
    }
    n
}

/// giac's `poslbdLMQ` of a polynomial given from its leading coefficient
/// down.
fn lower(p: &[Rational]) -> LowerBound {
    let k = p.len();
    if k <= 1 {
        return LowerBound::Zero;
    }
    let flip = p[k - 1].is_negative();
    let cl: Vec<Rational> = p
        .iter()
        .rev()
        .map(|c| if flip { -c } else { c.clone() })
        .collect();
    let mut largest: Option<Radical> = None;
    let mut used = vec![1_u32; k];
    for m in 1..k {
        if !cl[m].is_negative() {
            continue;
        }
        let mut least: Option<Radical> = None;
        for n in 0..m {
            if !cl[n].is_positive() {
                continue;
            }
            let t = Radical::of(&cl[m], &cl[n], used[n], m - n);
            used[n] += 1;
            if least
                .as_ref()
                .is_none_or(|l| l.cmp(&t) == Ordering::Greater)
            {
                least = Some(t);
            }
        }
        // No positive coefficient below: +∞ as the largest, so 2^(−∞).
        let Some(least) = least else {
            return LowerBound::Zero;
        };
        if largest
            .as_ref()
            .is_none_or(|l| least.cmp(l) == Ordering::Greater)
        {
            largest = Some(least);
        }
    }
    let Some(Radical { q, k }) = largest else {
        return LowerBound::Infinity;
    };
    // The least e with 2^e ≥ q^(1/k), that is 2^(e k) ≥ q; log₂ q is
    // within one of the bit lengths' difference.
    let k = i64::from(k);
    let bits = |n: &Integer| i64::try_from(n.bit_length()).expect("a size in memory");
    let d = bits(q.numer()) - bits(q.denom());
    let mut e = (d - 1).div_euclid(k) - 1;
    while two().pow(e * k) < q {
        e += 1;
    }
    LowerBound::Power(-e)
}

/// `p(x + c)`, by repeated synthetic division.
fn taylor(p: &mut [Rational], c: &Rational) {
    let n = p.len();
    for i in 0..n.saturating_sub(1) {
        for j in 1..n - i {
            let add = &p[j - 1] * c;
            p[j] += add;
        }
    }
}

/// `p(c x)`.
fn scale(p: &mut [Rational], c: &Rational) {
    let mut power = Rational::ONE;
    for coefficient in p.iter_mut().rev() {
        *coefficient = &*coefficient * &power;
        power = &power * c;
    }
}

/// The sign changes of the coefficients, zeros left out.
fn variations(p: &[Rational]) -> usize {
    let mut last = None;
    let mut count = 0;
    for c in p.iter().filter(|c| !c.is_zero()) {
        let positive = c.is_positive();
        if last.is_some_and(|l| l != positive) {
            count += 1;
        }
        last = Some(positive);
    }
    count
}

/// giac's `makeinterval`: the two ends in order.
fn interval(a: Rational, b: Rational) -> Isolated {
    if a > b {
        Isolated::Interval(b, a)
    } else {
        Isolated::Interval(a, b)
    }
}

/// The interval of a map with one root: from `b/d` to `a/c`, or from `b`
/// to `b + a·posubLMQ(f)` when `c = 0`.
fn one_root(a: &Rational, b: &Rational, c: &Rational, d: &Rational, f: &[Rational]) -> Isolated {
    if c.is_zero() {
        interval(b.clone(), b + &(a * &Rational::from(upper(f))))
    } else {
        interval(a / c, b / d)
    }
}

/// A map `x ↦ (a x + b)/(c x + d)` still to search, with `f` the
/// polynomial in the new variable and `v` its sign variations.
struct Map {
    a: Rational,
    b: Rational,
    c: Rational,
    d: Rational,
    f: Vec<Rational>,
    v: usize,
}

/// giac's `VAS_positive_roots(p, a, 0, 0, 1)`: the positive roots of `p`,
/// square-free and from its leading coefficient down, each mapped by `x ↦
/// a x`.
fn positive(p: Vec<Rational>, a: &Rational) -> Vec<Isolated> {
    let mut found = Vec::new();
    let v0 = variations(&p);
    if v0 == 0 {
        return found;
    }
    if v0 == 1 {
        found.push(interval(Rational::ZERO, a * &Rational::from(upper(&p))));
        return found;
    }
    let mut work = vec![Map {
        a: a.clone(),
        b: Rational::ZERO,
        c: Rational::ZERO,
        d: Rational::ONE,
        f: p,
        v: v0,
    }];
    let one = Rational::ONE;
    while let Some(Map {
        mut a,
        mut b,
        mut c,
        mut d,
        mut f,
        mut v,
    }) = work.pop()
    {
        // Steps 3 to 5: move by the lower bound when it is at least 1,
        // scaled to 1 first when above 16.
        let mut lb = lower(&f).value();
        if let Some(l) = lb
            .clone()
            .filter(|l| *l > Rational::from(Integer::from(16)))
        {
            scale(&mut f, &l);
            a = &l * &a;
            c = &l * &c;
            lb = Some(Rational::ONE);
        }
        if let Some(l) = lb.filter(|l| *l >= one) {
            taylor(&mut f, &l);
            b = &(&l * &a) + &b;
            d = &(&l * &c) + &d;
            if f.last().is_some_and(Rational::is_zero) {
                found.push(Isolated::Root(&b / &d));
                f.pop();
            }
            v = variations(&f);
            if v == 0 {
                continue;
            }
            if v == 1 {
                found.push(one_root(&a, &b, &c, &d, &f));
                continue;
            }
        }
        // Step 6: x ↦ x + 1, a root at 1 given exactly.
        let mut f1 = f.clone();
        taylor(&mut f1, &one);
        let (mut a1, mut b1, mut c1, mut d1) = (a.clone(), &a + &b, c.clone(), &c + &d);
        let mut r = 0;
        if f1.last().is_some_and(Rational::is_zero) {
            f1.pop();
            found.push(Isolated::Root(&b1 / &d1));
            r = 1;
        }
        let mut v1 = variations(&f1);
        // Budan's bound on (0, 1), computed when above 1: x ↦ 1/(x + 1).
        let mut v2 = v.saturating_sub(v1).saturating_sub(r);
        let (mut a2, mut b2, mut c2, mut d2) = (b, b1.clone(), d, d1.clone());
        let mut f2 = Vec::new();
        if v2 > 1 {
            f2 = f.iter().rev().cloned().collect();
            taylor(&mut f2, &one);
            if f2.last().is_some_and(Rational::is_zero) {
                f2.pop();
            }
            v2 = variations(&f2);
        }
        // Step 8: the branch with more variations first.
        if v1 < v2 {
            std::mem::swap(&mut a1, &mut a2);
            std::mem::swap(&mut b1, &mut b2);
            std::mem::swap(&mut c1, &mut c2);
            std::mem::swap(&mut d1, &mut d2);
            std::mem::swap(&mut f1, &mut f2);
            std::mem::swap(&mut v1, &mut v2);
        }
        // Steps 9 and 10.
        for (a, b, c, d, f, v) in [(a1, b1, c1, d1, f1, v1), (a2, b2, c2, d2, f2, v2)] {
            match v {
                0 => break,
                1 => found.push(one_root(&a, &b, &c, &d, &f)),
                _ => work.push(Map { a, b, c, d, f, v }),
            }
        }
    }
    // giac's vas_sort: by the left end, a root before an interval there.
    let key = |i: &Isolated| match i {
        Isolated::Root(r) => (r.clone(), 0),
        Isolated::Interval(lo, _) => (lo.clone(), 1),
    };
    found.sort_by_key(key);
    found
}

/// The square-free part of `p`.
fn squarefree_part(p: &Dense<Rational>) -> Dense<Rational> {
    let ring = DenseRing::new(Rationals);
    let g = EuclideanDomain::gcd(&ring, p, &ring.derivative(p));
    ring.div_rem(p, &g).0
}

/// giac's `posubLMQ`: an upper bound on the positive roots of `p`, 0
/// without a negative coefficient (the leading one made positive).
///
/// ```
/// use cas_num::Integer;
/// use cas_num::algebra::Rationals;
/// use cas_poly::DenseRing;
/// use cas_poly::real_roots::posub_lmq;
///
/// let qx = DenseRing::new(Rationals);
/// // x − 2: ⌈65·4/64⌉.
/// assert_eq!(posub_lmq(&qx.from_i64s(&[-2, 1])), Integer::from(5));
/// ```
#[must_use]
pub fn posub_lmq(p: &Dense<Rational>) -> Integer {
    upper(&descending(p))
}

/// giac's `poslbdLMQ`: a lower bound on the positive roots of `p`.
///
/// ```
/// use cas_num::algebra::Rationals;
/// use cas_poly::DenseRing;
/// use cas_poly::real_roots::{LowerBound, poslbd_lmq};
///
/// let qx = DenseRing::new(Rationals);
/// // x − 100: 2⁵; x³ + x + 1: none positive.
/// assert_eq!(poslbd_lmq(&qx.from_i64s(&[-100, 1])), LowerBound::Power(5));
/// assert_eq!(poslbd_lmq(&qx.from_i64s(&[1, 1, 0, 1])), LowerBound::Infinity);
/// ```
#[must_use]
pub fn poslbd_lmq(p: &Dense<Rational>) -> LowerBound {
    lower(&descending(p))
}

/// giac's `VAS_positive`: the positive roots of the square-free part of
/// `p`, isolated; empty for a constant or 0.
#[must_use]
pub fn vas_positive(p: &Dense<Rational>) -> Vec<Isolated> {
    if p.degree().is_none_or(|n| n == 0) {
        return Vec::new();
    }
    positive(descending(&squarefree_part(p)), &Rational::ONE)
}

/// giac's `VAS`: every real root of the square-free part of `p`,
/// isolated, ascending; empty for a constant or 0.
///
/// ```
/// use cas_num::{Integer, Rational};
/// use cas_num::algebra::Rationals;
/// use cas_poly::DenseRing;
/// use cas_poly::real_roots::{Isolated, vas};
///
/// let qx = DenseRing::new(Rationals);
/// let q = |n: i64, d: i64| Rational::new(Integer::from(n), Integer::from(d));
/// // x³ − 7x + 7.
/// let found = vas(&qx.from_i64s(&[7, -7, 0, 1]));
/// assert_eq!(
///     found,
///     [
///         Isolated::Interval(q(-4, 1), q(0, 1)),
///         Isolated::Interval(q(1, 1), q(3, 2)),
///         Isolated::Interval(q(3, 2), q(2, 1)),
///     ]
/// );
/// ```
#[must_use]
pub fn vas(p: &Dense<Rational>) -> Vec<Isolated> {
    if p.degree().is_none_or(|n| n == 0) {
        return Vec::new();
    }
    let mut v = descending(&squarefree_part(p));
    let has_zero = v.last().is_some_and(Rational::is_zero);
    if has_zero {
        v.pop();
    }
    // p(−x), times −1 for an odd degree.
    let n = v.len();
    let mut w: Vec<Rational> = v
        .iter()
        .enumerate()
        .map(|(i, c)| if (n - 1 - i) % 2 == 1 { -c } else { c.clone() })
        .collect();
    if n % 2 == 0 {
        w = w.iter().map(|c| -c).collect();
    }
    let symmetric = w == v;
    let positive_roots = positive(v, &Rational::ONE);
    let mut negative_roots = if symmetric {
        positive_roots
            .iter()
            .rev()
            .map(|i| match i {
                Isolated::Root(r) => Isolated::Root(-r),
                Isolated::Interval(lo, hi) => Isolated::Interval(-hi, -lo),
            })
            .collect()
    } else {
        positive(w, &-Rational::ONE)
    };
    if has_zero {
        negative_roots.push(Isolated::Root(Rational::ZERO));
    }
    negative_roots.extend(positive_roots);
    negative_roots
}

#[cfg(test)]
mod tests;
