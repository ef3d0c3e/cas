//! Bits, digits and falling factorials (CAS-DOC-45 §3, CAS-464).
//!
//! - `BitAnd`, `BitOr` and `BitXor`, Xcas's `bitand`, `bitor` and
//!   `bitxor`: of two or more integers, folded, in two's complement for
//!   negatives as giac's and Python's: `bitor(-12,10)` is `-2`. A negative
//!   `x` is `¬(|x| − 1)`, so each operation is one of naturals and
//!   complements.
//! - `HammingDistance[a, b]`, Xcas's `hamdist`: the one bits of `a xor b`,
//!   for nonnegative integers; giac takes negatives as 32-bit words, which
//!   here leave the call.
//! - `Convert[n, base, b]`, Xcas's `convert(n, base, b)`, `b ≥ 2`: the
//!   digits of `|n|`, lowest first, negated for a negative `n`, giac's
//!   `convert(-255,base,16)` being `[-15,-15]` and `convert(0,base,2)`
//!   `[0]`; `Convert[{d₀, d₁, …}, base, b]` is `Σ dᵢ bⁱ`.
//! - `FactorialPower[n, k]`, Xcas's `perm(n, k)`: `n (n − 1) ⋯ (n − k + 1)`
//!   for an integer `0 ≤ k`, a product for a symbol up to [`MAX_SYMBOLIC`]
//!   factors; 0 for an integer `0 ≤ n < k`. A product of numbers estimated
//!   above `Limits::number_bits` bits stays.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of integers: a and b, naturals x and y, digits d in base b"
)]

use cas_num::{Integer, Natural};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most factors of a falling factorial of a symbol.
const MAX_SYMBOLIC: i64 = 64;

/// An integer as a natural `m` and whether it is `¬m`: `x ≥ 0` is `x`,
/// `x < 0` is `¬(|x| − 1)`.
fn complement_form(x: &Integer) -> (Natural, bool) {
    if x.is_negative() {
        (&x.unsigned_abs() - &Natural::one(), true)
    } else {
        (x.unsigned_abs(), false)
    }
}

/// `m` or `¬m` as an integer.
fn from_complement_form(m: Natural, not: bool) -> Integer {
    if not {
        -Integer::from(&m + &Natural::one())
    } else {
        Integer::from(m)
    }
}

#[derive(Clone, Copy)]
enum Op {
    And,
    Or,
    Xor,
}

/// `a op b` in two's complement.
fn bitwise(op: Op, a: &Integer, b: &Integer) -> Integer {
    let ((x, nx), (y, ny)) = (complement_form(a), complement_form(b));
    // With ¬ on the naturals: ¬x & ¬y = ¬(x | y), x & ¬y = x − (x & y),
    // x | ¬y = ¬(y − (x & y)), ¬x | ¬y = ¬(x & y), and xor moves the ¬s out.
    let (m, not) = match (op, nx, ny) {
        (Op::And, false, false) => (&x & &y, false),
        (Op::And, false, true) => (&x - &(&x & &y), false),
        (Op::And, true, false) => (&y - &(&x & &y), false),
        (Op::And, true, true) => (&x | &y, true),
        (Op::Or, false, false) => (&x | &y, false),
        (Op::Or, false, true) => (&y - &(&x & &y), true),
        (Op::Or, true, false) => (&x - &(&x & &y), true),
        (Op::Or, true, true) => (&x & &y, true),
        (Op::Xor, _, _) => (&x ^ &y, nx != ny),
    };
    from_complement_form(m, not)
}

/// The integer arguments, two at least.
fn integers(e: &Expr) -> Option<Vec<Integer>> {
    let values: Option<Vec<Integer>> = e.args().iter().map(Expr::to_integer).collect();
    values.filter(|v| v.len() >= 2)
}

fn fold(ev: &mut Eval<'_>, e: &Expr, op: Op, name: &'static str) -> Native {
    let Some(values) = integers(e) else {
        return Ok(None);
    };
    let (first, rest) = values.split_first().expect("two values");
    let result = Expr::from(
        rest.iter()
            .fold(first.clone(), |acc, x| bitwise(op, &acc, x)),
    );
    ev.step(name, e, &result);
    Ok(Some(result))
}

/// `BitAnd[a, b, …]`.
pub(super) fn bit_and(ev: &mut Eval<'_>, e: &Expr) -> Native {
    fold(ev, e, Op::And, "bitand")
}

/// `BitOr[a, b, …]`.
pub(super) fn bit_or(ev: &mut Eval<'_>, e: &Expr) -> Native {
    fold(ev, e, Op::Or, "bitor")
}

/// `BitXor[a, b, …]`.
pub(super) fn bit_xor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    fold(ev, e, Op::Xor, "bitxor")
}

/// `HammingDistance[a, b]`.
pub(super) fn hamming_distance(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some([a, b]) = integers(e).and_then(|v| <[Integer; 2]>::try_from(v).ok()) else {
        return Ok(None);
    };
    if a.is_negative() || b.is_negative() {
        return Ok(None);
    }
    let ones = (&a.unsigned_abs() ^ &b.unsigned_abs()).count_ones();
    let result = Expr::from(Integer::from(ones));
    ev.step("hamdist", e, &result);
    Ok(Some(result))
}

/// `Convert[n, base, b]` and `Convert[{d₀, d₁, …}, base, b]`.
pub(super) fn convert(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, option, b] = e.args() else {
        return Ok(None);
    };
    if option.as_symbol().map(Symbol::name) != Some("base") {
        return Ok(None);
    }
    let Some(b) = b.to_integer().filter(|b| *b >= Integer::from(2)) else {
        return Ok(None);
    };
    let result = if let Some(n) = x.to_integer() {
        let negative = n.is_negative();
        let mut m = n.abs();
        let mut digits = Vec::new();
        loop {
            let (q, r) = m.div_rem_euclid(&b);
            digits.push(Expr::from(if negative { -r } else { r }));
            m = q;
            if m.is_zero() {
                break;
            }
        }
        Expr::call(Symbol::LIST, digits)
    } else if x.head().and_then(Expr::as_symbol) == Some(Symbol::LIST) {
        let Some(digits) = x
            .args()
            .iter()
            .map(Expr::to_integer)
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(None);
        };
        let n = digits
            .iter()
            .rev()
            .fold(Integer::ZERO, |acc, d| &(&acc * &b) + d);
        Expr::from(n)
    } else {
        return Ok(None);
    };
    ev.step("convert", e, &result);
    Ok(Some(result))
}

/// `FactorialPower[n, k]`.
pub(super) fn factorial_power(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, k] = e.args() else {
        return Ok(None);
    };
    let Some(k) = k.as_i64().filter(|k| *k >= 0) else {
        return Ok(None);
    };
    if let Some(m) = n.to_integer()
        && !m.is_negative()
        && m < Integer::from(k)
    {
        return Ok(Some(Expr::from(0)));
    }
    if !n.is_number() && k > MAX_SYMBOLIC {
        return Ok(None);
    }
    if n.is_number() {
        // log₂ of each factor at most that of |n| + k.
        let bound = n
            .to_integer()
            .map_or(64, |m| (&m.abs() + &Integer::from(k)).bit_length());
        let bits = u64::try_from(k)
            .expect("nonnegative")
            .saturating_mul(bound.max(1));
        if bits > ev.context().limits().number_bits {
            return Ok(None);
        }
    }
    let factors = (0..k).map(|j| Expr::call(Symbol::PLUS, [n.clone(), Expr::from(-j)]));
    let result = Expr::call(Symbol::TIMES, factors.collect::<Vec<_>>());
    ev.step("perm", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
