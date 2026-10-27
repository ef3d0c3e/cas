//! Exponentials and logarithms in the search (CAS-DOC-14 §4, CAS-230):
//! candidates, each a rewriting of every part, sound where it applies.
//!
//! - Products of exponentials as one: `Exp[a] Exp[b]` as `Exp[a + b]`, the
//!   exponent expanded, with `E^a` and `Exp[a]^n` for an integer `n` among
//!   them; the logarithms of the sum come out in the next round.
//! - Logarithms out of exponents: `Exp[c Log[x] + a]` as `x^c Exp[a]` where
//!   `x^c` is `e^(c log x)`: `c` an integer, `x` positive, or the Complex
//!   domain, where powers are principal (in the Real domain an odd root of
//!   a negative is real, O7).
//! - Logarithms combined: the terms `c Log[a]` of a sum, `c` an integer and
//!   `a` positive, as one `Log[∏ a^c]`, and `c Log[a]` alone as `Log[a^c]`.
//! - Logarithms expanded: `Log[p x]` as `Log[p] + Log[x]` for a positive
//!   `p`; `Log[a^c]` as `c Log[a]` for a positive `a` and a real `c`, or as
//!   `c Log[Abs[a]]` for a real `a` and an even `c`; and `Log[q]` of a
//!   positive rational of at most 64 bits over its primes, `Log[8]` as
//!   `3 Log[2]`.

use cas_num::{Integer, Rational};

use crate::assume::{Domain, Set, Sign, Truth, element, sign_of};
use crate::context::Context;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The bits of a rational's parts whose logarithm is written over primes.
const MAX_BITS: u64 = 64;

/// The rewritings of `e` by the three kinds, those that change it.
pub(super) fn rewritings(ctx: &Context, e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    if ![Symbol::EXP, Symbol::E, Symbol::LOG]
        .into_iter()
        .any(|s| e.may_contain(s))
    {
        return out;
    }
    for rewrite in [exponentials, combined, expanded] {
        let r = e.map_bottom_up(|x| rewrite(ctx, x));
        if r != *e && !out.contains(&r) {
            out.push(r);
        }
    }
    out
}

/// Whether `x` is positive by the assumptions.
fn positive(ctx: &Context, x: &Expr) -> bool {
    sign_of(ctx, x) == Sign::Positive
}

/// Whether `x` is real by the domain or the assumptions.
fn real(ctx: &Context, x: &Expr) -> bool {
    element(ctx, x, Set::Reals) == Truth::True
}

/// A part with its exponentials as one, and logarithms out of exponents.
fn exponentials(ctx: &Context, e: &Expr) -> Expr {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::TIMES), args) => {
            let mut exponents = Vec::new();
            let mut rest = Vec::new();
            for a in args {
                match exponent(a) {
                    Some(u) => exponents.push(u),
                    None => rest.push(a.clone()),
                }
            }
            if exponents.len() < 2 {
                return e.clone();
            }
            // Expanded: evaluation leaves x − (x − 1), of exp(x)/exp(x − 1).
            let sum = Expr::call(Symbol::EXPAND, [Expr::call(Symbol::PLUS, exponents)]);
            rest.push(Expr::call(Symbol::EXP, [sum]));
            Expr::call(Symbol::TIMES, rest)
        }
        (Some(Symbol::EXP), [u]) => out_of_logs(ctx, u),
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => out_of_logs(ctx, u),
        _ => e.clone(),
    }
}

/// The exponent of `Exp[u]`, `E^u`, or of such a power to an integer.
fn exponent(e: &Expr) -> Option<Expr> {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::EXP), [u]) => Some(u.clone()),
        (Some(Symbol::POWER), [b, u]) if b.as_symbol() == Some(Symbol::E) => Some(u.clone()),
        (Some(Symbol::POWER), [b, n]) if n.to_integer().is_some() => {
            exponent(b).map(|u| n.clone() * u)
        }
        _ => None,
    }
}

/// `Exp[u]` with the terms `c Log[x]` of `u` out as `x^c`, where that is
/// its value.
fn out_of_logs(ctx: &Context, u: &Expr) -> Expr {
    let terms: Vec<Expr> = match u.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => u.args().to_vec(),
        _ => vec![u.clone()],
    };
    let mut kept = Vec::new();
    let mut powers = Vec::new();
    for t in terms {
        match log_term(&t) {
            Some((c, x))
                if c.to_integer().is_some()
                    || positive(ctx, &x)
                    || ctx.domain() == Domain::Complex =>
            {
                powers.push(x.pow(c));
            }
            _ => kept.push(t),
        }
    }
    if powers.is_empty() {
        return Expr::call(Symbol::EXP, [u.clone()]);
    }
    powers.push(Expr::call(Symbol::EXP, [Expr::call(Symbol::PLUS, kept)]));
    Expr::call(Symbol::TIMES, powers)
}

/// `t` as `c Log[x]`: `(c, x)`.
fn log_term(t: &Expr) -> Option<(Expr, Expr)> {
    let is_log =
        |a: &Expr| a.head().and_then(Expr::as_symbol) == Some(Symbol::LOG) && a.args().len() == 1;
    if is_log(t) {
        return Some((Expr::from(1), t.args()[0].clone()));
    }
    if t.head().and_then(Expr::as_symbol) != Some(Symbol::TIMES) {
        return None;
    }
    let args = t.args();
    let mut logs = args.iter().enumerate().filter(|(_, a)| is_log(a));
    let (i, log) = logs.next()?;
    if logs.next().is_some() {
        return None;
    }
    let mut c: Vec<Expr> = args
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, a)| a.clone())
        .collect();
    let c = if c.len() == 1 {
        c.pop().expect("one factor")
    } else {
        Expr::call(Symbol::TIMES, c)
    };
    Some((c, log.args()[0].clone()))
}

/// A part with its terms `c Log[a]`, `c` an integer and `a` positive, as one
/// logarithm.
fn combined(ctx: &Context, e: &Expr) -> Expr {
    let integer_log = |t: &Expr| -> Option<(Expr, Expr)> {
        let (c, x) = log_term(t)?;
        (c.to_integer().is_some() && positive(ctx, &x)).then_some((c, x))
    };
    match e.head().and_then(Expr::as_symbol) {
        Some(Symbol::PLUS) => {
            let mut product = Vec::new();
            let mut rest = Vec::new();
            for t in e.args() {
                match integer_log(t) {
                    Some((c, x)) => product.push(x.pow(c)),
                    None => rest.push(t.clone()),
                }
            }
            if product.len() < 2 {
                return e.clone();
            }
            rest.push(Expr::call(
                Symbol::LOG,
                [Expr::call(Symbol::TIMES, product)],
            ));
            Expr::call(Symbol::PLUS, rest)
        }
        Some(Symbol::TIMES) => match integer_log(e) {
            Some((c, x)) => Expr::call(Symbol::LOG, [x.pow(c)]),
            None => e.clone(),
        },
        _ => e.clone(),
    }
}

/// A part with its logarithm expanded where it can be.
fn expanded(ctx: &Context, e: &Expr) -> Expr {
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::LOG), [x]) => log_expanded(ctx, x).unwrap_or_else(|| e.clone()),
        _ => e.clone(),
    }
}

/// `Log[x]` expanded; `None` when it stays.
fn log_expanded(ctx: &Context, x: &Expr) -> Option<Expr> {
    if let Some(q) = x.to_rational() {
        return over_primes(&q);
    }
    match (x.head().and_then(Expr::as_symbol), x.args()) {
        (Some(Symbol::TIMES), args) => {
            let (out, rest): (Vec<&Expr>, Vec<&Expr>) = args.iter().partition(|f| positive(ctx, f));
            if out.is_empty() {
                return None;
            }
            let mut terms: Vec<Expr> = out.into_iter().map(|f| log_of(ctx, f)).collect();
            if !rest.is_empty() {
                let rest: Vec<Expr> = rest.into_iter().cloned().collect();
                terms.push(Expr::call(Symbol::LOG, [Expr::call(Symbol::TIMES, rest)]));
            }
            Some(Expr::call(Symbol::PLUS, terms))
        }
        (Some(Symbol::POWER), [b, c]) => {
            if positive(ctx, b) && real(ctx, c) {
                return Some(c.clone() * log_of(ctx, b));
            }
            let even = c
                .to_integer()
                .is_some_and(|n| n.unsigned_abs().trailing_zeros().is_none_or(|t| t > 0));
            (even && real(ctx, b)).then(|| {
                c.clone() * Expr::call(Symbol::LOG, [Expr::call(Symbol::ABS, [b.clone()])])
            })
        }
        _ => None,
    }
}

/// `Log[f]`, expanded where it can be.
fn log_of(ctx: &Context, f: &Expr) -> Expr {
    log_expanded(ctx, f).unwrap_or_else(|| Expr::call(Symbol::LOG, [f.clone()]))
}

/// `Log[q]` of a positive rational over its primes; `None` for one that is
/// not positive, 1, a prime, or past [`MAX_BITS`].
fn over_primes(q: &Rational) -> Option<Expr> {
    if !q.is_positive() {
        return None;
    }
    let mut powers: Vec<(i64, Integer)> = Vec::new();
    for (n, sign) in [(q.numer(), 1_i64), (q.denom(), -1)] {
        if n.bit_length() > MAX_BITS {
            return None;
        }
        for f in n.unsigned_abs().factor().factors() {
            let e = i64::try_from(f.exponent()).ok()?;
            powers.push((sign * e, Integer::from(f.base().clone())));
        }
    }
    // 1, or a prime: nothing to write.
    if matches!(powers.as_slice(), [] | [(1, _)]) {
        return None;
    }
    let terms: Vec<Expr> = powers
        .into_iter()
        .map(|(e, p)| Expr::from(e) * Expr::call(Symbol::LOG, [Expr::from(p)]))
        .collect();
    Some(Expr::call(Symbol::PLUS, terms))
}
