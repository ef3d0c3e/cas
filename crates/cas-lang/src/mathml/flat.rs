//! giac's `aplatir_fois_plus` (gen.cc :12505), which `gen2mathml` runs on a
//! symbolic before writing it (mathml.cc :1929): sums nested on the left
//! flattened (`aplatir_plus` :12474), a product of two factors multiplied
//! out by `sym_mult` (gen.cc :7228), so a quoted `2*inv(4)*x` is `1/2*x`.
//! Run once over the whole tree: flattening a flattened tree changes
//! nothing, so this is giac's run at each symbolic.
//!
//! Of `sym_mult`, the branches a printed tree reaches: `a*inv(a)`, 0 and
//! 1, `neg` factors, fractions (`_FRAC2_SYMB` :1952), an integer over an
//! integer, negative numbers, and the merging of products. Not ported:
//! relations, units, programs, infinities, doubles times exact constants
//! (`evalf`), `rootof` and modular factors.

use cas_core::{Integer, Rational};

use super::print::{gtype, is_negative, is_one, neg_num, ty};
use crate::markup::tree::Gen;

/// `aplatir_fois_plus(g)`.
pub(super) fn aplatir(g: &Gen) -> Gen {
    match g {
        Gen::Vect(items, sub) => Gen::Vect(items.iter().map(aplatir).collect(), *sub),
        Gen::Symb(op, f) => match (op.as_str(), &**f) {
            ("plus", Gen::Vect(..)) => aplatir_plus(g),
            ("prod", Gen::Vect(items, _)) if items.len() == 2 => {
                sym_mult(&aplatir(&items[0]), &aplatir(&items[1]))
            }
            _ => Gen::symb(op, aplatir(f)),
        },
        _ => g.clone(),
    }
}

/// `aplatir_plus` (gen.cc :12474): `plus(plus(a, b), c)` is `plus(a, b, c)`.
fn aplatir_plus(g: &Gen) -> Gen {
    let mut items = g.args();
    if items.is_empty() {
        return Gen::int(0);
    }
    let mut v = Vec::new();
    loop {
        let (first, rest) = items.split_first().expect("not empty");
        v.extend(rest.iter().rev().map(aplatir));
        match first {
            Gen::Symb(op, f) if op == "plus" => match &**f {
                Gen::Vect(inner, _) if !inner.is_empty() => items = inner,
                _ => {
                    v.push(first.clone());
                    break;
                }
            },
            _ => {
                v.push(aplatir(first));
                break;
            }
        }
    }
    v.reverse();
    Gen::symb_n("plus", v)
}

/// A rational as giac's integer or fraction.
fn rat(q: &Rational) -> Gen {
    if q.denom().is_one() {
        Gen::Int(q.numer().clone())
    } else {
        Gen::Frac(
            Box::new(Gen::Int(q.numer().clone())),
            Box::new(Gen::Int(q.denom().clone())),
        )
    }
}

/// An exact number as a rational: an integer, a boolean, or a fraction of
/// integers.
fn exact(g: &Gen) -> Option<Rational> {
    match g {
        Gen::Int(n) => Some(Rational::from(n.clone())),
        Gen::Bool(b) => Some(Rational::from(Integer::from(i64::from(*b)))),
        Gen::Frac(a, b) => match (&**a, &**b) {
            (Gen::Int(a), Gen::Int(b)) => Rational::checked_new(a.clone(), b.clone()),
            _ => None,
        },
        _ => None,
    }
}

/// `inv(n)` for an integer `n`: `n`.
fn inv_int(g: &Gen) -> Option<&Integer> {
    match g {
        Gen::Symb(op, f) if op == "inv" => match &**f {
            Gen::Int(n) => Some(n),
            _ => None,
        },
        _ => None,
    }
}

/// `a / n` folded, `None` for `n = 0`.
fn over(a: &Rational, n: &Integer) -> Option<Gen> {
    a.checked_div(&Rational::from(n.clone())).map(|q| rat(&q))
}

/// `_FRAC2_SYMB` (gen.cc :1952): `p/q` as `inv(q)`, `-inv(q)` or
/// `p*inv(q)`.
fn frac2symb(num: &Gen, den: &Gen) -> Gen {
    let inv = Gen::symb("inv", den.clone());
    if is_one(num) {
        inv
    } else if super::print::equals(num, -1, 1) {
        Gen::symb("neg", inv)
    } else {
        Gen::symb_n("prod", vec![num.clone(), inv])
    }
}

/// The factors of a product, or the one.
fn factors(g: &Gen) -> Vec<Gen> {
    match g {
        Gen::Symb(op, _) if op == "prod" => g.args().to_vec(),
        _ => vec![g.clone()],
    }
}

/// `sym_mult(a, b)` (gen.cc :7228).
#[allow(clippy::many_single_char_names)]
fn sym_mult(a: &Gen, b: &Gen) -> Gen {
    let is_inv_of = |x: &Gen, y: &Gen| matches!(x, Gen::Symb(op, f) if op == "inv" && **f == *y);
    // :7313-7316.
    if is_inv_of(a, b) || is_inv_of(b, a) {
        return Gen::int(1);
    }
    // :7317-7328: a zero factor.
    for x in [a, b] {
        if matches!(x, Gen::Int(_) | Gen::Bool(_) | Gen::Dbl { .. }) && super::print::is_zero(x) {
            return x.clone();
        }
    }
    // :7377-7380.
    if is_one(a) {
        return b.clone();
    }
    if is_one(b) {
        return a.clone();
    }
    // :7432-7438.
    if let Gen::Symb(op, x) = a
        && op == "neg"
    {
        if let Gen::Symb(op, y) = b
            && op == "neg"
        {
            return sym_mult(x, y);
        }
        return Gen::symb("neg", sym_mult(x, b));
    }
    if let Gen::Symb(op, y) = b
        && op == "neg"
    {
        return Gen::symb("neg", sym_mult(a, y));
    }
    // :7453-7478: a fraction times an integer, an integer's inverse or a
    // number folds; times a symbolic, it is `p*inv(q)`.
    if let Gen::Frac(p, q) = a {
        if let (Some(x), Some(n)) = (exact(a), inv_int(b))
            && let Some(r) = over(&x, n)
        {
            return r;
        }
        if let (Some(x), Some(y)) = (exact(a), exact(b)) {
            return rat(&(&x * &y));
        }
        if !matches!(b, Gen::Symb(..) | Gen::Idnt(_)) {
            // A fraction times a double or a complex: not folded.
            return Gen::symb_n("prod", vec![a.clone(), b.clone()]);
        }
        return sym_mult(&frac2symb(p, q), b);
    }
    if let Gen::Frac(p, q) = b {
        if let (Some(y), Some(n)) = (exact(b), inv_int(a))
            && let Some(r) = over(&y, n)
        {
            return r;
        }
        if let (Some(x), Some(y)) = (exact(a), exact(b)) {
            return rat(&(&x * &y));
        }
        if !matches!(a, Gen::Symb(..) | Gen::Idnt(_)) {
            return Gen::symb_n("prod", vec![a.clone(), b.clone()]);
        }
        return sym_mult(a, &frac2symb(p, q));
    }
    // :7479-7482: an integer times an integer's inverse.
    if let (Gen::Int(n), Some(d)) = (a, inv_int(b))
        && let Some(r) = over(&Rational::from(n.clone()), d)
    {
        return r;
    }
    if let (Gen::Int(n), Some(d)) = (b, inv_int(a))
        && let Some(r) = over(&Rational::from(n.clone()), d)
    {
        return r;
    }
    // :7487-7490: a negative real number's minus taken out.
    if gtype(a) <= ty::REAL && is_negative(a) {
        return Gen::symb("neg", sym_mult(&neg_num(a), b));
    }
    if gtype(b) <= ty::REAL && is_negative(b) {
        return Gen::symb("neg", sym_mult(a, &neg_num(b)));
    }
    // :7529-7555: products merged; an `_INT_` after a symbolic goes first.
    let merged = match (a, b) {
        (Gen::Symb(..), Gen::Symb(..)) => [factors(a), factors(b)].concat(),
        (_, Gen::Symb(..)) => [vec![a.clone()], factors(b)].concat(),
        (Gen::Symb(op, _), _) if op == "prod" => [factors(a), vec![b.clone()]].concat(),
        (Gen::Symb(..), _) if gtype(b) == ty::INT => vec![b.clone(), a.clone()],
        _ => vec![a.clone(), b.clone()],
    };
    Gen::symb_n("prod", merged)
}
