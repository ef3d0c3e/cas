//! The canonical order on expressions (`impl Ord for Expr`).

use std::cmp::Ordering;

use cas_num::ball::{Ball, Float};
use cas_num::{Integer, Rational};

use crate::expr::{Expr, Repr, View};
use crate::symbol::Symbol;

/// The implicit exponent of a factor and coefficient of a product.
static ONE: Expr = Expr(Repr::Int(1));

/// The canonical order: deterministic and total, `Equal` exactly for equal
/// expressions.
///
/// Evaluation sorts the arguments of `Orderless` heads into this order, so it
/// decides how sums and products print. It depends only on structure, names
/// and values: never on symbol ids, hashes, addresses or interning order, so
/// it is the same in every run.
///
/// 1. **Numbers** come first: integers and rationals by value, then complex
///    numbers by real part and then imaginary part (`2 < I < 1 - I < 1 + I`).
/// 2. **Products and powers** — `Times[…]`, and `Power[b, n]` with two
///    arguments — compare as products of powers, as do all other expressions
///    against them (`x` is `x^1`, and `x^2` is the product of one factor).
///    That puts polynomial terms in the conventional order,
///    `1 + 2 x + x^2 + 2 y + 2 x y + y^2`:
///    - the factors that are not numbers are compared from the last one
///      backwards, each by base and then by exponent (implicitly 1), and a
///      product whose factors run out first comes first;
///    - then the numeric coefficient (implicitly 1), so `-x < x < 2 x`.
/// 3. **Everything else**, and ties left by rule 2 (`x` and `Times[x]`), by
///    kind: strings, then symbols, then applications. Strings and symbol
///    names are ordered alphabetically ignoring case, then lowercase before
///    uppercase (`a < A < b`), then by code point. A symbol name with a `$`
///    in front, how the language stores a user's capitalised name (`$E` for
///    `E`), sorts as the name without it, just before that name itself.
///    Applications are ordered
///    by head, then by argument count, then by their arguments from the
///    first.
///
/// Rule 2 reads `Times` and `Power` structurally, so it applies to
/// expressions that are not canonical as well, such as `Times[y, x]`, and the
/// order stays total. Comparison recurses on the tree, like evaluation:
/// expressions nested deeper than the stack allows cannot be compared.
///
/// ```
/// use cas_core::{Expr, Symbol};
///
/// let (x, y) = (Expr::from(Symbol::new("x")), Expr::from(Symbol::new("y")));
/// let pow = |b: &Expr, n: i64| Expr::apply(&Symbol::POWER.into(), &[b.clone(), n.into()]);
/// let times = |a: &[Expr]| Expr::apply(&Symbol::TIMES.into(), a);
///
/// let mut terms = vec![pow(&y, 2), times(&[x.clone(), y.clone()]), y.clone(), pow(&x, 2), x.clone(), 1.into()];
/// terms.sort();
/// assert_eq!(format!("{terms:?}"), "[1, x, Power[x, 2], y, Times[x, y], Power[y, 2]]");
/// ```
impl Ord for Expr {
    fn cmp(&self, other: &Expr) -> Ordering {
        canonical(self, other)
    }
}

impl PartialOrd for Expr {
    fn partial_cmp(&self, other: &Expr) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn canonical(a: &Expr, b: &Expr) -> Ordering {
    if a == b {
        return Ordering::Equal;
    }
    match (is_number(a), is_number(b)) {
        (true, true) => return cmp_numbers(a, b),
        (true, false) => return Ordering::Less,
        (false, true) => return Ordering::Greater,
        (false, false) => {}
    }
    if is_product(a) || is_product(b) {
        cmp_products(a, b).then_with(|| structural(a, b))
    } else {
        structural(a, b)
    }
}

fn is_number(e: &Expr) -> bool {
    e.is_number()
}

/// Real numbers by value, then complex numbers by real part and imaginary
/// part; only rationals need a `Rational` comparison. A ball sorts by its
/// midpoint (CAS-DOC-5 §9.1), after an exact number of the same value, then
/// by radius; a complex ball part by part.
fn cmp_numbers(a: &Expr, b: &Expr) -> Ordering {
    let integer = |v: View<'_>| match v {
        View::Int(x) => Some(Integer::from(x)),
        View::Integer(x) => Some(x.clone()),
        _ => None,
    };
    let complex = |v: View<'_>| matches!(v, View::Complex(_) | View::ComplexBall(_));
    match (a.view(), b.view()) {
        (View::Int(l), View::Int(r)) => l.cmp(&r),
        (l, r) if complex(l) && complex(r) => {
            let ((lre, lim), (rre, rim)) = (parts(l), parts(r));
            cmp_reals(&lre, &rre).then_with(|| cmp_reals(&lim, &rim))
        }
        (l, _) if complex(l) => Ordering::Greater,
        (_, r) if complex(r) => Ordering::Less,
        (View::Ball(_), _) | (_, View::Ball(_)) => cmp_reals(&real(a.view()), &real(b.view())),
        (l, r) => match (integer(l), integer(r)) {
            (Some(l), Some(r)) => l.cmp(&r),
            _ => rational(a).cmp(&rational(b)),
        },
    }
}

fn rational(e: &Expr) -> Rational {
    e.to_rational().expect("a number")
}

/// A real number as the order sees it: an exact value, or a ball.
enum Real<'a> {
    Exact(Rational),
    Ball(&'a Ball),
}

fn real(v: View<'_>) -> Real<'_> {
    match v {
        View::Ball(b) => Real::Ball(b),
        View::Int(x) => Real::Exact(Rational::from(x)),
        View::Integer(x) => Real::Exact(Rational::from(x.clone())),
        View::Rational(q) => Real::Exact(q.clone()),
        _ => unreachable!("a real number"),
    }
}

/// The real and imaginary parts of a complex number or ball.
fn parts(v: View<'_>) -> (Real<'_>, Real<'_>) {
    match v {
        View::Complex(z) => (Real::Exact(z.re().clone()), Real::Exact(z.im().clone())),
        View::ComplexBall(z) => (Real::Ball(z.re()), Real::Ball(z.im())),
        _ => unreachable!("a complex number"),
    }
}

/// By value, a ball by midpoint; an exact value before a ball of that
/// midpoint, and balls of one midpoint by radius.
fn cmp_reals(a: &Real<'_>, b: &Real<'_>) -> Ordering {
    match (a, b) {
        (Real::Exact(x), Real::Exact(y)) => x.cmp(y),
        (Real::Exact(x), Real::Ball(b)) => cmp_exact(x, b.mid()).then(Ordering::Less),
        (Real::Ball(a), Real::Exact(y)) => cmp_exact(y, a.mid()).reverse().then(Ordering::Greater),
        (Real::Ball(a), Real::Ball(b)) => a.mid().cmp(b.mid()).then_with(|| a.rad().cmp(b.rad())),
    }
}

/// `q` against the float `m`, exactly. A float too large or small for an
/// exact rational is beyond every rational there is memory for, so its sign
/// and exponent decide.
fn cmp_exact(q: &Rational, m: &Float) -> Ordering {
    if let Some(m) = m.checked_to_rational() {
        return q.cmp(&m);
    }
    let huge = !m.exponent().is_negative();
    match (m.is_negative(), huge) {
        (false, true) => Ordering::Less,
        (true, true) => Ordering::Greater,
        // |m| is below every nonzero |q|.
        (_, false) => q.cmp(&Rational::from(0)).then(if m.is_negative() {
            Ordering::Greater
        } else {
            Ordering::Less
        }),
    }
}

/// `Times[…]`, or `Power[b, n]`.
fn is_product(e: &Expr) -> bool {
    match e.view() {
        View::Apply { head, args } => {
            let h = head.as_symbol();
            h == Some(Symbol::TIMES) || (h == Some(Symbol::POWER) && args.len() == 2)
        }
        _ => false,
    }
}

/// The factors of `e` as a product, numbers included.
fn factors(e: &Expr) -> &[Expr] {
    match e.view() {
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::TIMES) => args,
        _ => std::slice::from_ref(e),
    }
}

/// `(base, exponent)` of a factor.
fn power(e: &Expr) -> (&Expr, &Expr) {
    match e.view() {
        View::Apply { head, args: [b, n] } if head.as_symbol() == Some(Symbol::POWER) => (b, n),
        _ => (e, &ONE),
    }
}

/// Rule 2. Terminates: a base compared here is a strict subexpression of a
/// product, and the other side's base is a subexpression of the other side.
fn cmp_products(a: &Expr, b: &Expr) -> Ordering {
    cmp_by(symbolic_factors(a), symbolic_factors(b), cmp_factors)
        .then_with(|| cmp_by(coefficient(a), coefficient(b), canonical))
}

/// By base, then by exponent.
fn cmp_factors(a: &Expr, b: &Expr) -> Ordering {
    let ((base_a, exp_a), (base_b, exp_b)) = (power(a), power(b));
    canonical(base_a, base_b).then_with(|| canonical(exp_a, exp_b))
}

/// The factors that are not numbers, last first.
fn symbolic_factors(e: &Expr) -> impl Iterator<Item = &Expr> {
    factors(e).iter().rev().filter(|f| !is_number(f))
}

/// The numeric factors, or the implicit 1.
fn coefficient(e: &Expr) -> impl Iterator<Item = &Expr> {
    let numbers = factors(e).iter().filter(|f| is_number(f));
    let implicit = numbers.clone().next().is_none().then_some(&ONE);
    numbers.chain(implicit)
}

/// Rule 3, for pairs that are not both numbers (rule 1 decides those).
fn structural(a: &Expr, b: &Expr) -> Ordering {
    debug_assert!(!(is_number(a) && is_number(b)));
    match (a.view(), b.view()) {
        (View::String(l), View::String(r)) => cmp_names(l, r),
        (View::Symbol(l), View::Symbol(r)) => cmp_symbol_names(l.name(), r.name()),
        (
            View::Apply {
                head: head_a,
                args: args_a,
            },
            View::Apply {
                head: head_b,
                args: args_b,
            },
        ) => canonical(head_a, head_b)
            .then(args_a.len().cmp(&args_b.len()))
            .then_with(|| cmp_by(args_a.iter(), args_b.iter(), canonical)),
        (l, r) => rank(l).cmp(&rank(r)),
    }
}

fn rank(v: View<'_>) -> u8 {
    match v {
        View::Int(_)
        | View::Integer(_)
        | View::Rational(_)
        | View::Complex(_)
        | View::Ball(_)
        | View::ComplexBall(_) => 0,
        View::String(_) => 1,
        View::Symbol(_) => 2,
        View::Apply { .. } => 3,
    }
}

/// Lexicographic by `cmp`; a prefix comes first.
fn cmp_by<'a>(
    a: impl IntoIterator<Item = &'a Expr>,
    b: impl IntoIterator<Item = &'a Expr>,
    cmp: fn(&Expr, &Expr) -> Ordering,
) -> Ordering {
    let mut b = b.into_iter();
    for l in a {
        let Some(r) = b.next() else {
            return Ordering::Greater;
        };
        match cmp(l, r) {
            Ordering::Equal => {}
            o => return o,
        }
    }
    if b.next().is_some() {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// Symbol names as [`cmp_names`], reading a name the language escaped (`$E`
/// for a user's `E`, CAS-DOC-7 §5) as it was written, so it sorts where its
/// spelling does. A tie with the unescaped name, the builtin `E`, goes by the
/// full names.
fn cmp_symbol_names(a: &str, b: &str) -> Ordering {
    cmp_names(unescaped(a), unescaped(b)).then_with(|| cmp_names(a, b))
}

fn unescaped(name: &str) -> &str {
    name.strip_prefix('$').unwrap_or(name)
}

/// Alphabetical ignoring case, then lowercase first, then by code point; the
/// last pass makes it `Equal` only for equal strings.
fn cmp_names<'a>(a: &'a str, b: &'a str) -> Ordering {
    let fold = |s: &'a str| s.chars().flat_map(char::to_lowercase);
    let case = |c: char| (c.is_uppercase(), c);
    fold(a)
        .cmp(fold(b))
        .then_with(|| a.chars().map(case).cmp(b.chars().map(case)))
}

#[cfg(test)]
mod tests;
