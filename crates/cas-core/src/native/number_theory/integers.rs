//! Xcas's integer functions of CAS-DOC-55: `chrem`, `lgcd`, `fracmod`,
//! `icomp` and `bernoulli_mod`, ports of giac's code where its answers
//! depend on it. giac's errors leave the call, with the error as a warning.
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

use std::collections::HashMap;

use cas_num::modular::{crt, fracmod};
use cas_num::{Integer, Natural};

use super::{magnitude, refuse};
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
