//! Xcas's integer functions of CAS-DOC-55: `chrem`, `lgcd`, `fracmod`,
//! `icomp`, `bernoulli_mod` and `is_prime`, ports of giac's code where its
//! answers depend on it. giac's errors leave the call, with the error as a
//! warning.
//!
//! - `ChineseRemainderLists[{r₁, …}, {m₁, …}]` (`chrem`, CAS-548), or of one
//!   list `{{r₁, …}, {m₁, …}}`: giac's left fold (misc.cc `_chrem`, usual.cc
//!   `_ichinrem` and `ichinrem2`). One pair is itself, unchecked. A step
//!   from two integers is `ChineseRemainder`'s: moduli at least 2, the value
//!   in `[0, lcm)` modulo the lcm, none where the residues disagree. A step
//!   with a polynomial residue goes coefficient by coefficient over the
//!   monomials of the expanded residues, modulo the product `m n` (not the
//!   lcm), each coefficient giac's `ichinrem(a, b, m, n)` (gen.cc): when `a`,
//!   `b`, `m` and `n` fit 32 bits and `m`, `n` are coprime, `a + ((m⁻¹ (b −
//!   a)) rem n) m`, with C's truncated `rem` and giac's signed inverse
//!   ([`invmod`]); else `smod(A q m + a, m n)`, `A` GMP's cofactor of `m`,
//!   `q = (b − a)/gcd(m, n)`. Unnormalised: `chrem([x^2+1,2*x],[3,5])` is
//!   `[-5*x^2+12*x-5,15]` (giac prints `(-5*x+12)*x-5`). A step whose value
//!   is an integer takes it into `[0, m n)`. Moduli that are not all
//!   integers fold through `PolynomialChineseRemainder` (libgiac crashes).
//! - `ListGCD[{a₁, …}]` (`lgcd`, CAS-549), or of a sequence: giac's fold
//!   (modpoly.cc `lgcd`) `gcd(…gcd(gcd(a₁, a₁), a₂)…)` by `GCD`, stopped at
//!   1; `lgcd([])` is 1, and one argument that is not a list stays. `GCD`
//!   makes a polynomial's leading coefficient positive: `lgcd([-2*x-2])` is
//!   `2*x+2` (giac `-2*x-2`).
//! - `FractionMod[a, m]` (`fracmod`, `iratrecon`, CAS-550): cas-num's
//!   `fracmod`, giac's `in_fracmod`, Wang's loop without a denominator bound;
//!   a cofactor not prime to `m` leaves the call (giac's "Reconstructed
//!   denominator is not prime with modulo"), as do an `a` or `m` that is not
//!   an integer and a negative `m` (libgiac crashes); `a = 0` is 0, lists
//!   map.
//! - `IntegerCompositions[n, k]` (`icomp`, CAS-553), `1 ≤ k ≤ n`: the
//!   compositions of `n` into `k` parts, zeros allowed, in the order of
//!   Nijenhuis and Wilf's NEXCOM (giac's graphtheory.cc `nexcom`); the option
//!   `zeros = false` (`zeros` lowers to `Solve`) drops those with a 0 part.
//!   More entries than `Limits::terms` (`C(n+k−1, k−1) k`) leave the call.
//! - `BernoulliMod[k, p]` (`bernoulli_mod`, CAS-552): `B_k mod p` for an even
//!   `k ≥ 2` and a prime `p < 2³¹`, by Harvey's sums as giac computes them
//!   (ifactor.cc `bernoulli_mod`, without bernmm), its C integer arithmetic
//!   kept: the value is giac's, in `(−p, p)`, its sign included. `(p − 1) |
//!   k` leaves the call (giac's "Not invertible"), `p = 2` too (libgiac
//!   crashes), and so does a `p` whose `p/2` steps pass `1024
//!   Limits::terms`.
//! - `IsPrime[n]` (`is_prime`, CAS-551): giac's 2 for a prime below 2³¹, 1
//!   above, 0 for every other integer (negatives included); lists map. Every
//!   prime is proven, by cas-num's `Natural::prove_prime` (deterministic
//!   below 2⁶⁴, the `n − 1` test of Brillhart, Lehmer and Selfridge above)
//!   within `IntegerFactorList`'s limits; a prime it does not prove leaves
//!   the call, with a warning, where giac answers 1 through PARI's APRCL.
//!   `IsPrime[n, 1]`: giac's own certificate (ifactor.cc `prime_cert`, not
//!   PARI's): 0 unless `n` is prime, else `[[q, a, f], …]` over the primes
//!   `q` of `n − 1` that giac's `prime_factors` lists ([`certificate_factors`],
//!   with giac's Pollard rho ported, so the order is giac's), `a` the least
//!   base proving `q` (Pocklington), `f` 1 for `q < 2³¹`; where giac's rho
//!   gives up, its ECM (libecm) and ours may find another factor first
//!   (`10^30+57`). Any other integer `c` is `IsPrime[n]`; a non-integer `n`
//!   or `c` leaves the call (giac's Bad Argument Type).

use std::collections::HashMap;

use cas_num::modular::{crt, fracmod};
use cas_num::{Integer, Natural, Proof};

use super::{ECM_B1, effort, magnitude, refuse};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// A value, or giac's error for the call; evaluation's own errors outside.
type Step<T> = Result<Result<T, &'static str>, EvalError>;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn is_list(e: &Expr) -> bool {
    head(e) == Some(Symbol::LIST)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// `ChineseRemainderLists[{r₁, …}, {m₁, …}]`.
pub(in crate::native) fn chinese_remainder_lists(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (residues, moduli) = match e.args() {
        [r, m] => (r, m),
        [rm] if is_list(rm) && rm.args().len() == 2 => (&rm.args()[0], &rm.args()[1]),
        _ => return refuse(ev, e, "bad argument value"),
    };
    let (rs, ms) = (residues.args(), moduli.args());
    if !is_list(residues) || !is_list(moduli) || rs.len() != ms.len() || rs.is_empty() {
        return refuse(ev, e, "bad argument value");
    }
    if rs.len() == 1 {
        return Ok(Some(list(vec![rs[0].clone(), ms[0].clone()])));
    }
    let Some(integers) = ms.iter().map(Expr::to_integer).collect::<Option<Vec<_>>>() else {
        return polynomial_moduli(ev, rs, ms);
    };
    let two = Integer::from(2);
    let (mut r, mut m) = (rs[0].clone(), integers[0].clone());
    for (b, n) in rs[1..].iter().zip(&integers[1..]) {
        if let (Some(x), Some(y)) = (r.to_integer(), b.to_integer()) {
            if m < two || *n < two {
                return refuse(ev, e, "bad mod value");
            }
            let (mm, nn) = (magnitude(&m), magnitude(n));
            let x = magnitude(&x.div_rem_euclid(&m).1);
            let y = magnitude(&y.div_rem_euclid(n).1);
            let Some((v, l)) = crt(&x, &mm, &y, &nn) else {
                return refuse(ev, e, "no integer solution");
            };
            (r, m) = (Expr::from(Integer::from(v)), Integer::from(l));
            continue;
        }
        if m < two || *n < two {
            return refuse(ev, e, "bad mod value");
        }
        match polynomial_step(ev, &r, b, &m, n)? {
            Ok(v) => (r, m) = (v, m.clone() * n),
            Err(why) => return refuse(ev, e, why),
        }
    }
    Ok(Some(list(vec![r, Expr::from(m)])))
}

/// The fold through `PolynomialChineseRemainder`, moduli not all integers.
fn polynomial_moduli(ev: &mut Eval<'_>, rs: &[Expr], ms: &[Expr]) -> Native {
    let mut acc = list(vec![rs[0].clone(), ms[0].clone()]);
    for (r, m) in rs[1..].iter().zip(&ms[1..]) {
        let pair = list(vec![r.clone(), m.clone()]);
        acc = ev.eval(&Expr::call(
            Symbol::POLYNOMIAL_CHINESE_REMAINDER,
            [acc, pair],
        ))?;
        if !is_list(&acc) || acc.args().len() != 2 {
            return Ok(None);
        }
    }
    Ok(Some(acc))
}

/// One step of giac's `ichinrem2` with a polynomial residue: `a` modulo `m`
/// and `b` modulo `n` coefficient by coefficient; the value modulo `m n`,
/// or giac's error.
fn polynomial_step(ev: &mut Eval<'_>, a: &Expr, b: &Expr, m: &Integer, n: &Integer) -> Step<Expr> {
    let ta = match terms(ev, a)? {
        Ok(t) => t,
        Err(why) => return Ok(Err(why)),
    };
    let tb = match terms(ev, b)? {
        Ok(t) => t,
        Err(why) => return Ok(Err(why)),
    };
    // The monomials in order of appearance, each with its two coefficients.
    let mut at: HashMap<Expr, usize> = HashMap::new();
    let mut table: Vec<(Expr, Integer, Integer)> = Vec::new();
    for (k, c) in ta {
        at.insert(k.clone(), table.len());
        table.push((k, c, Integer::ZERO));
    }
    for (k, c) in tb {
        if let Some(&i) = at.get(&k) {
            table[i].2 = c;
        } else {
            table.push((k, Integer::ZERO, c));
        }
    }
    let mut sum = Vec::with_capacity(table.len());
    for (k, ca, cb) in table {
        let Some(c) = ichinrem(&ca, &cb, m, n) else {
            return Ok(Err("no integer solution"));
        };
        sum.push(Expr::call(Symbol::TIMES, [Expr::from(c), k]));
    }
    let value = ev.eval(&Expr::call(Symbol::PLUS, sum))?;
    Ok(Ok(match value.to_integer() {
        // `_ichinrem`'s `irem(x + m n, m n)`.
        Some(v) => Expr::from(v.div_rem_euclid(&(m.clone() * n)).1),
        None => value,
    }))
}

/// The monomials of `p` expanded, each with its integer coefficient (1 for
/// the constant term), or giac's "non integer argument".
fn terms(ev: &mut Eval<'_>, p: &Expr) -> Step<Vec<(Expr, Integer)>> {
    let expanded = ev.eval(&Expr::call(Symbol::EXPAND, [p.clone()]))?;
    let parts = if head(&expanded) == Some(Symbol::PLUS) {
        expanded.args().to_vec()
    } else {
        vec![expanded]
    };
    let mut out = Vec::with_capacity(parts.len());
    for t in parts {
        if let Some(c) = t.to_integer() {
            if !c.is_zero() {
                out.push((Expr::from(1), c));
            }
            continue;
        }
        let (c, monomial) = match (head(&t), t.args()) {
            (Some(Symbol::TIMES), [c, rest @ ..]) if c.is_number() => {
                let monomial = match rest {
                    [one] => one.clone(),
                    _ => Expr::call(Symbol::TIMES, rest.to_vec()),
                };
                (c.to_integer(), monomial)
            }
            _ if t.is_number() || t.is_numeric() => (None, t),
            _ => (Some(Integer::ONE), t),
        };
        match c {
            Some(c) if !monomial.is_numeric() => out.push((monomial, c)),
            _ => return Ok(Err("non integer argument")),
        }
    }
    Ok(Ok(out))
}

/// giac's `ichinrem(a, b, m, n)` (gen.cc): `x ≡ a (mod m)` and `x ≡ b (mod
/// n)` as giac finds it, or `None` for none.
fn ichinrem(a: &Integer, b: &Integer, m: &Integer, n: &Integer) -> Option<Integer> {
    let small = |x: &Integer| x.to_i64().and_then(|x| i32::try_from(x).ok());
    if let (Some(a), Some(b), Some(m), Some(n)) = (small(a), small(b), small(m), small(n))
        && gcd(m, n) == 1
    {
        let inverse = i64::from(invmod(m, n)?);
        let x =
            i64::from(a) + (inverse * (i64::from(b) - i64::from(a)) % i64::from(n)) * i64::from(m);
        return Some(Integer::from(x));
    }
    let (d, s, _) = m.gcd_ext(n);
    let diff = b - a;
    let q = if d.is_one() {
        diff
    } else if diff.is_divisible_by(&d) {
        diff.div_exact(&d)
    } else {
        return None;
    };
    let modulus = m * n;
    let r = (s * &q * m + a).div_rem_euclid(&modulus).1;
    // giac's smod: (−M/2, M/2].
    Some(if &r + &r > modulus { r - &modulus } else { r })
}

/// The gcd of two 32-bit integers, as giac's `gcd(int, int)` for moduli.
fn gcd(mut a: i32, mut b: i32) -> i32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

/// giac's `invmod(int, int)` (gen.cc): Euclid with C's truncated division,
/// so the inverse may be negative; `None` for "Not invertible".
fn invmod(mut a: i32, mut b: i32) -> Option<i32> {
    if a == 1 || a == -1 || a == 1 - b {
        return Some(a);
    }
    let (mut aa, mut ab) = (1i32, 0i32);
    while b != 0 {
        let (q, r) = (a / b, a % b);
        let ar = aa - q * ab;
        (a, b, aa, ab) = (b, r, ab, ar);
    }
    match a {
        1 => Some(aa),
        -1 => Some(-aa),
        _ => None,
    }
}

/// `ListGCD[{a₁, …}]` and `ListGCD[a₁, …]`.
pub(in crate::native) fn list_gcd(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let items = match e.args() {
        [l] if is_list(l) => l.args(),
        [_] => return Ok(None),
        args => args,
    };
    let Some(first) = items.first() else {
        return Ok(Some(Expr::from(1)));
    };
    let one = Expr::from(1);
    let mut g = first.clone();
    for a in items {
        g = ev.eval(&Expr::call(Symbol::GCD, [g, a.clone()]))?;
        if head(&g) == Some(Symbol::GCD) {
            return Ok(None);
        }
        if g == one {
            break;
        }
    }
    Ok(Some(g))
}

const NOT_PRIME: &str = "reconstructed denominator is not prime with modulo";

/// `FractionMod[a, m]`.
pub(in crate::native) fn fraction_mod(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, m] = e.args() else {
        return Ok(None);
    };
    if a.to_integer().is_some_and(|a| a.is_zero()) {
        return Ok(Some(Expr::from(0)));
    }
    let Some(m) = m.to_integer() else {
        return refuse(ev, e, NOT_PRIME);
    };
    if m.is_negative() {
        return refuse(ev, e, "negative modulus");
    }
    match reconstructed(a, &magnitude(&m)) {
        Some(v) => Ok(Some(v)),
        None => refuse(ev, e, NOT_PRIME),
    }
}

/// `a`'s fraction modulo `m`, lists entry by entry; `None` if one fails.
fn reconstructed(a: &Expr, m: &Natural) -> Option<Expr> {
    if is_list(a) {
        let items = a.args().iter().map(|x| reconstructed(x, m));
        return items.collect::<Option<Vec<_>>>().map(list);
    }
    fracmod(&a.to_integer()?, m).map(Expr::from)
}

/// `IntegerCompositions[n, k]`, `IntegerCompositions[n, k, zeros = b]`.
pub(in crate::native) fn integer_compositions(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (n, k, zeros) = match e.args() {
        [n, k] => (n, k, true),
        [n, k, option] => {
            let zeros = match (head(option), option.args()) {
                (Some(Symbol::EQUATION | Symbol::EQUAL), [name, value])
                    if name.as_symbol() == Some(Symbol::SOLVE) =>
                {
                    match value.as_symbol() {
                        Some(Symbol::TRUE) => Some(true),
                        Some(Symbol::FALSE) => Some(false),
                        _ => value.to_integer().map(|v| !v.is_zero()),
                    }
                }
                _ => None,
            };
            let Some(zeros) = zeros else {
                return refuse(ev, e, "bad argument value");
            };
            (n, k, zeros)
        }
        [_] => return refuse(ev, e, "bad argument type"),
        _ => return refuse(ev, e, "bad argument value"),
    };
    let (Some(n), Some(k)) = (n.as_i64(), k.as_i64()) else {
        return refuse(ev, e, "bad argument value");
    };
    if k < 1 || k > n {
        return refuse(ev, e, "bad argument value");
    }
    // C(n+k−1, k−1) compositions of k entries each, bounded first (C0019).
    let limit = u128::try_from(ev.context().limits().terms).unwrap_or(u128::MAX);
    let (n_wide, k_wide) = (u128::from(n.unsigned_abs()), u128::from(k.unsigned_abs()));
    let mut count: u128 = 1;
    for i in 1..k_wide {
        count = match count.checked_mul(n_wide - 1 + i) {
            Some(c) => c / i,
            None => return refuse(ev, e, "more entries than the limit on terms"),
        };
        if count.saturating_mul(k_wide) > limit {
            return refuse(ev, e, "more entries than the limit on terms");
        }
    }
    if count.saturating_mul(k_wide) > limit {
        return refuse(ev, e, "more entries than the limit on terms");
    }
    let k = usize::try_from(k).expect("k entries are within the limit");
    let found = compositions(n, k)
        .into_iter()
        .filter(|c| zeros || !c.contains(&0))
        .map(|c| list(c.into_iter().map(Expr::from).collect()))
        .collect();
    Ok(Some(list(found)))
}

/// The compositions of `n` into `k ≥ 1` parts in NEXCOM's order (giac's
/// `nexcom`): from `(n, 0, …, 0)` to `(0, …, 0, n)`.
fn compositions(n: i64, k: usize) -> Vec<Vec<i64>> {
    let mut r = vec![0; k];
    r[0] = n;
    let (mut h, mut t) = (0usize, n);
    let mut out = vec![r.clone()];
    while r[k - 1] != n {
        if t > 1 {
            h = 0;
        }
        h += 1;
        t = r[h - 1];
        r[h - 1] = 0;
        r[0] = t - 1;
        r[h] += 1;
        out.push(r.clone());
    }
    out
}

const NOT_INVERTIBLE: &str = "not invertible";

/// `BernoulliMod[k, p]`.
pub(in crate::native) fn bernoulli_mod_native(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [k, p] = e.args() else {
        return refuse(ev, e, "bad argument value");
    };
    let small = |x: &Expr| x.as_i64().and_then(|x| i32::try_from(x).ok());
    let (Some(k), Some(p)) = (small(k), small(p)) else {
        return refuse(ev, e, "bad argument type");
    };
    let prime = p >= 2 && magnitude(&Integer::from(p)).primality().is_probable_prime();
    if k < 2 || k % 2 != 0 || !prime {
        return refuse(ev, e, "bad argument type");
    }
    // giac recurses without end at p = 2 and fails to invert at (p − 1) | k.
    if p == 2 || (k > p - 3 && k % (p - 1) == 0) {
        return refuse(ev, e, NOT_INVERTIBLE);
    }
    let steps = u64::try_from(p / 2).expect("p is positive");
    let terms = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    if steps > terms.saturating_mul(1024) {
        return refuse(ev, e, "more steps than the limit on terms allows");
    }
    match bernoulli_mod(k, p) {
        Some(b) => Ok(Some(Expr::from(i64::from(b)))),
        None => refuse(ev, e, NOT_INVERTIBLE),
    }
}

/// C's conversion of a `long long` to an `int`: the low 32 bits.
#[allow(
    clippy::cast_possible_truncation,
    reason = "giac's implicit long long to int conversions, kept"
)]
fn wrap(x: i64) -> i32 {
    x as i32
}

/// giac's `powmod(int, unsigned long, int)` (gen.cc): `a` itself for `n =
/// 1`, not reduced.
fn powmod(a: i32, n: u64, m: i32) -> i32 {
    if n == 0 {
        return 1;
    }
    if n == 1 {
        return a;
    }
    let m = i64::from(m);
    if n == 2 {
        return wrap(i64::from(a) * i64::from(a) % m);
    }
    let (mut b, mut c, mut n) = (i64::from(a) % m, 1i64, n);
    while n > 0 {
        if n % 2 == 1 {
            c = c * b % m;
        }
        n /= 2;
        b = b * b % m;
    }
    wrap(c)
}

/// The prime factors of `n ≥ 1` ascending with their exponents, giac's
/// `ifactors`.
fn prime_factors(mut n: i32) -> Vec<(i32, u32)> {
    let mut out = Vec::new();
    let mut d = 2;
    while d <= n / d {
        let mut e = 0;
        while n % d == 0 {
            n /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

/// giac's `generator(p, v)` (ifactor.cc): the least generator of `(ℤ/p)*`.
fn generator(p: i32, factors: &[(i32, u32)]) -> i32 {
    let w: Vec<u64> = factors
        .iter()
        .map(|&(q, _)| u64::try_from((p - 1) / q).expect("positive"))
        .collect();
    for a in 2..p {
        let mut r = 0;
        for &wi in &w {
            r = powmod(a, wi, p);
            if r == 1 {
                break;
            }
        }
        if r != 1 {
            return a;
        }
    }
    0
}

/// giac's `znorder(2, p, p − 1, v)` (ifactor.cc): the order of 2 modulo `p`.
fn order_of_two(p: i32, factors: &[(i32, u32)]) -> i32 {
    let mut o = 1;
    for &(q, e) in factors {
        let mut a = powmod(2, u64::try_from((p - 1) / q.pow(e)).expect("positive"), p);
        while a != 1 {
            o *= q;
            a = powmod(a, u64::try_from(q).expect("positive"), p);
        }
    }
    o
}

/// giac's `precond_mulmod31` (ifactor.cc): `b₁ m mod p` up to a multiple of
/// `p`, by the precomputed `⌊2³¹ m / p⌋ + 1`.
fn precond_mulmod31(b1: i32, m: i32, p: i32, msurp: i32) -> i32 {
    let (b1, m, p, msurp) = (i64::from(b1), i64::from(m), i64::from(p), i64::from(msurp));
    wrap(b1 * m - ((b1 * msurp) >> 31) * p)
}

/// giac's `bernoulli_mod(int k, int p)` (ifactor.cc, Harvey's algorithm, no
/// bernmm), its `int` and `long long` arithmetic as i32 and i64; `None`
/// where giac fails to invert.
fn bernoulli_mod(k: i32, p: i32) -> Option<i32> {
    let p64 = i64::from(p);
    if k > p - 3 {
        let m = k % (p - 1);
        if m == 0 {
            return None;
        }
        let bm = bernoulli_mod(m, p)?;
        return Some(wrap(
            i64::from(bm) * i64::from(k) % p64 * i64::from(invmod(m, p)?) % p64,
        ));
    }
    let factors = prime_factors(p - 1);
    let g = generator(p, &factors);
    let r = powmod(g, u64::try_from(k - 1).expect("k ≥ 2"), p);
    let order = if p > 11 { order_of_two(p, &factors) } else { 0 };
    if order > 4 && k % order != 0 {
        // Faster summation.
        let n = if order % 2 != 0 { order } else { order / 2 };
        let m = (p - 1) / 2 / n;
        let exponent = u64::try_from((k - 1) % order).expect("positive");
        let twokm1 = i64::from(powmod(2, exponent, p));
        let (mut total, mut gi, mut gkm1i) = (0i64, 1i64, 1i64);
        let msurp = wrap((1i64 << 31) * twokm1 / p64 + 1);
        for _ in 0..m {
            let (mut s, mut pow2) = (0i64, 1i64);
            let mut gi2j = wrap(gi);
            for _ in 0..n {
                gi2j = gi2j.wrapping_shl(1).wrapping_sub(p);
                let sign = gi2j >> 31;
                s -= i64::from(1 + (sign << 1)) * pow2;
                gi2j -= sign * p;
                pow2 = i64::from(precond_mulmod31(wrap(pow2), wrap(twokm1), p, msurp));
            }
            gi = gi * i64::from(g) % p64;
            s = s % p64 * gkm1i % p64;
            total += s;
            if total >= p64 {
                total -= p64;
            }
            gkm1i = gkm1i * i64::from(r) % p64;
        }
        total = total * i64::from(k) % p64;
        let inverse = invmod(invmod(wrap(twokm1 % p64), p)? - 2, p)?;
        return Some(wrap(total * i64::from(inverse) % p64));
    }
    let u = if g % 2 != 0 {
        (g - 1) / 2
    } else {
        wrap(i64::from(g - 1) * i64::from(invmod(2, p)?) % p64)
    };
    let (mut s, mut x, mut y) = (0i32, 1i32, r);
    for _ in 1..=p / 2 {
        let q = wrap(i64::from(g) * i64::from(x) / p64);
        s = wrap((i64::from(s) + (i64::from(u) - i64::from(q)) * i64::from(y)) % p64);
        x = wrap(i64::from(g) * i64::from(x) % p64);
        y = wrap(i64::from(r) * i64::from(y) % p64);
    }
    let res = wrap(2 * i64::from(k) * i64::from(s) % p64);
    let inverse = invmod(1 - powmod(g, u64::try_from(k).expect("k ≥ 2"), p), p)?;
    Some(wrap(i64::from(res) * i64::from(inverse) % p64))
}

/// giac's primes are `int`s below 2³¹; larger ones it answers 1.
const INT_LIMIT: i64 = 1 << 31;

/// The most bases a certificate tries for one factor (giac's `RAND_MAX`).
const MAX_WITNESS: u64 = 10_000;

const UNPROVEN: &str = "primality not proven within the factorisation limits";

/// `IsPrime[n]`, `IsPrime[n, c]`.
pub(in crate::native) fn is_prime(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let answer = match e.args() {
        [n] => prime_level(n),
        [n, c] => match c.to_integer() {
            Some(c) if c.is_one() => certificate(n),
            Some(_) => prime_level(n),
            None => Err("bad argument type"),
        },
        _ => return Ok(None),
    };
    match answer {
        Ok(x) => Ok(Some(x)),
        Err(why) => refuse(ev, e, why),
    }
}

/// The proof of `m`'s primality within `IntegerFactorList`'s limits.
fn prove(m: &Natural) -> Result<bool, &'static str> {
    let (rho, curves) = effort(m);
    match m.prove_prime(rho, ECM_B1, curves) {
        Proof::Proven => Ok(true),
        Proof::Composite => Ok(false),
        Proof::Unknown => Err(UNPROVEN),
    }
}

/// giac's 0, 1 or 2 for `n`, mapped over lists.
fn prime_level(n: &Expr) -> Result<Expr, &'static str> {
    if is_list(n) {
        return n
            .args()
            .iter()
            .map(prime_level)
            .collect::<Result<_, _>>()
            .map(list);
    }
    let n = n.to_integer().ok_or("bad argument type")?;
    if n < Integer::from(2) || !prove(&magnitude(&n))? {
        return Ok(Expr::from(0));
    }
    Ok(Expr::from(if n < Integer::from(INT_LIMIT) { 2 } else { 1 }))
}

/// giac's `prime_cert` (ifactor.cc): 0 for a number that is not prime, else
/// `[[q, a, f], …]` over [`certificate_factors`] of `n − 1`, `a` the least
/// base with `a^(n−1) ≡ 1` and `gcd(a^((n−1)/q) − 1, n) = 1`, `f` 1 for
/// `q < 2³¹`.
fn certificate(n: &Expr) -> Result<Expr, &'static str> {
    let n = n.to_integer().ok_or("bad argument type")?;
    if n < Integer::from(2) || !prove(&magnitude(&n))? {
        return Ok(Expr::from(0));
    }
    let n = magnitude(&n);
    let o = &n - 1;
    let factors = certificate_factors(&o).ok_or("could not partially factor n - 1")?;
    let mut entries = Vec::with_capacity(factors.len());
    for q in factors {
        let cofactor = o.div_exact(&q);
        let a = (2..2 + MAX_WITNESS)
            .map(Natural::from)
            .find(|a| {
                a.pow_mod(&o, &n).is_one()
                    && (Integer::from(a.pow_mod(&cofactor, &n)) - Integer::ONE)
                        .unsigned_abs()
                        .gcd(&n)
                        .is_one()
            })
            .ok_or("no witness found")?;
        let small = q < Natural::from(INT_LIMIT.unsigned_abs());
        entries.push(list(vec![
            Expr::from(Integer::from(q)),
            Expr::from(Integer::from(a)),
            Expr::from(i64::from(small)),
        ]));
    }
    Ok(list(entries))
}

/// giac's `prime_factors(o, false)` (ifactor.cc): primes of `o` until the
/// part left is at most `√o`. Trial division to 10⁴, pushing each prime
/// once; then, while the part `g` left is above `√o` and not a probable
/// prime (pushed last), a divisor `b` of `g` by [`pollard`], else by the
/// elliptic curve method, its primes pushed in increasing order, `g / b`.
/// After trial division, a `g` with `g² ≤ o` is pushed whole. `None` when
/// a divisor is not found or not factored within the limits.
fn certificate_factors(o: &Natural) -> Option<Vec<Natural>> {
    let (rho, curves) = effort(o);
    let stop = o.isqrt();
    let mut g = o.clone();
    let mut found = Vec::new();
    // A composite d never divides g: its primes are already removed.
    for d in 2..10_000u64 {
        let d = Natural::from(d);
        let (mut quotient, mut rest) = g.div_rem(&d);
        if !rest.is_zero() {
            continue;
        }
        while rest.is_zero() {
            g = quotient;
            (quotient, rest) = g.div_rem(&d);
        }
        found.push(d);
        if stop >= g {
            return Some(found);
        }
        if *o >= &g * &g {
            found.push(g);
            return Some(found);
        }
    }
    loop {
        if stop >= g {
            return Some(found);
        }
        if g.primality().is_probable_prime() {
            found.push(g);
            return Some(found);
        }
        let b = match pollard(&g, 1, rho) {
            Some(b) => b,
            None => g
                .perfect_power()
                .map(|(base, _)| base)
                .or_else(|| g.ecm(ECM_B1, curves, 0))?,
        };
        let f = b.factor_with_ecm(rho, ECM_B1, curves);
        if !f.is_complete() {
            return None;
        }
        found.extend(f.factors().iter().map(|x| x.base().clone()));
        g = g.div_exact(&b);
    }
}

/// giac's `pollard(n, k)` (ifactor.cc): Pollard's rho on `x ↦ x² + k` from
/// `x = 2`, Brent's cycle detection with a GCD every 64 products and a
/// backtrack from the last batch; on a GCD of `n` it retries with `k = −1,
/// 3, 5, …, 49`. At most `max(15000, 1500 · 16^((log₁₀ n − 40)/10))`
/// iterations, kept under `cap` above 15000 (ours); `None` when they run
/// out. C's truncated `%` kept.
fn pollard(n: &Natural, k: i64, cap: u64) -> Option<Natural> {
    let maxiter = pollard_iterations(n).min(i64::try_from(cap).unwrap_or(i64::MAX).max(15_000));
    let modulus = Integer::from(n);
    let k_int = Integer::from(k);
    let step = |x: &Integer| &(&(x * x) + &k_int) % &modulus;
    let gcd = |p: &Integer| p.unsigned_abs().gcd(n);
    let (mut m, mut m1, mut a, mut a1, mut c) = (2i64, 2i64, 1i64, 1i64, 0);
    let mut g = Natural::one();
    let mut x = Integer::from(2);
    let (mut x1, mut y, mut y1, mut p) = (x.clone(), x.clone(), Integer::ZERO, Integer::ONE);
    while g.is_one() {
        a = 2 * a + 1;
        while g.is_one() && a > m {
            x = step(&x);
            m += 1;
            if m > maxiter {
                return None;
            }
            p = &(&p * &(&x1 - &x)) % &modulus;
            c += 1;
            if c == 64 {
                g = gcd(&p);
                if g.is_one() {
                    (y, y1, p, a1, m1, c) = (x.clone(), x1.clone(), Integer::ONE, a, m, 0);
                }
            }
        }
        if g.is_one() {
            x1 = x.clone();
            let j = 3 * (a + 1) / 2;
            for _ in m + 1..=j {
                x = step(&x);
            }
            m = j;
        }
    }
    (x, x1, g, a, m) = (y, y1, Natural::one(), (a1 - 1) / 2, m1);
    while g.is_one() {
        a = 2 * a + 1;
        while g.is_one() && a > m {
            x = step(&x);
            m += 1;
            if m > maxiter {
                return None;
            }
            g = gcd(&(&(&x1 - &x) % &modulus));
        }
        if g.is_one() {
            x1 = x.clone();
            let j = 3 * (a + 1) / 2;
            for _ in m + 1..=j {
                x = step(&x);
            }
            m = j;
        }
    }
    if g == *n {
        return match k {
            1 => pollard(n, -1, cap),
            -1 => pollard(n, 3, cap),
            k if k >= 50 => None,
            k => pollard(n, k + 2, cap),
        };
    }
    Some(g)
}

/// giac's `int(1500 · 16^((log₁₀ n − 40)/10))`, at least 15000. Past
/// 2³¹ (`n` above about 10⁹¹) C's conversion gives `INT_MIN` on x86, so
/// giac's bound stays 15000.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "giac computes the bound in doubles"
)]
fn pollard_iterations(n: &Natural) -> i64 {
    let shift = n.bit_length().saturating_sub(64);
    let top = u64::try_from(&(n >> shift)).expect("64 bits");
    let log10 = (top as f64).log10() + shift as f64 * std::f64::consts::LOG10_2;
    let bound = 1500.0 * 16f64.powf((log10 - 40.0) / 10.0);
    if bound >= 2_147_483_648.0 {
        return 15_000;
    }
    (bound as i64).max(15_000)
}
