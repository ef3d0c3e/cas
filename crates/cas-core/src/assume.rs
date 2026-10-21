//! Assumptions and the number domain (CAS-DOC-4 §7, R9, R12, CAS-58).
//!
//! A context's assumptions are facts that predicates and native functions
//! consult during evaluation. They are kept as written, and read into
//! [`Facts`]: ranges of *subjects*, number sets of expressions, and the
//! assumptions that are neither, which hold as written.
//!
//! A subject is an expression in a form that does not depend on how a fact
//! was written: an expression is `scale·subject + shift`, with a rational
//! scale and shift and a monic subject, whose first term has coefficient 1.
//! So `x > y + 1`, `2 x - 2 y > 2` and `y - x < -1` all bound `x - y` below
//! by 1, and a query reads its expression the same way.
//!
//! What is known of a value is a [`Range`]: whether it is real, bounds with
//! exact rational ends, and whether it is not 0. It comes from the facts
//! about the subject and from the structure:
//! - numbers are points, and `Pi` and `E` lie between rationals;
//! - sums add their ranges, and a numeric coefficient scales them;
//! - products and powers combine signs;
//! - `Abs` is not negative, and `Exp` of a real is positive.
//!
//! In the Real domain (§7.1), the default, a symbol without facts is real.
//! The answers are sound and cheap, not complete: what cannot be told is
//! unknown.

use foldhash::fast::FixedState;
use hashbrown::HashMap;

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use crate::context::{Context, EvalError};
use crate::eval::{Eval, run};
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// What symbols without facts range over (CAS-DOC-4 §7.1, R12).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum Domain {
    /// The real numbers, as in Xcas: the default.
    #[default]
    Real,
    /// The complex numbers.
    Complex,
}

/// Whether a proposition holds, as far as can be told.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Truth {
    True,
    False,
    Unknown,
}

impl Truth {
    fn of(b: bool) -> Truth {
        if b { Truth::True } else { Truth::False }
    }

    fn not(self) -> Truth {
        match self {
            Truth::True => Truth::False,
            Truth::False => Truth::True,
            Truth::Unknown => Truth::Unknown,
        }
    }
}

/// What is known of the sign of a value. Only a real value has a sign other
/// than `Nonzero` and `Unknown`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Sign {
    Positive,
    Negative,
    Zero,
    NonNegative,
    NonPositive,
    Nonzero,
    Unknown,
}

/// The number sets of `Element`, from the smallest.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub(crate) enum Set {
    Integers,
    Rationals,
    Reals,
    Complexes,
}

impl Set {
    pub(crate) fn of(e: &Expr) -> Option<Set> {
        match e.as_symbol()? {
            Symbol::INTEGERS => Some(Set::Integers),
            Symbol::RATIONALS => Some(Set::Rationals),
            Symbol::REALS => Some(Set::Reals),
            Symbol::COMPLEXES => Some(Set::Complexes),
            _ => None,
        }
    }

    /// The smallest of these sets a number is in; a real ball is in the
    /// reals.
    fn of_number(e: &Expr) -> Option<Set> {
        if let Some(q) = e.to_rational() {
            return Some(if q.is_integer() {
                Set::Integers
            } else {
                Set::Rationals
            });
        }
        if matches!(e.view(), View::Ball(_)) {
            return Some(Set::Reals);
        }
        e.is_number().then_some(Set::Complexes)
    }
}

/// A bound of a range: a rational, reached or not.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Bound {
    value: Rational,
    strict: bool,
}

/// What is known of a value: whether it is real, its bounds when it is
/// (`None` is unbounded), and whether it is not 0.
#[derive(Clone, Debug)]
pub(crate) struct Range {
    real: bool,
    lower: Option<Bound>,
    upper: Option<Bound>,
    nonzero: bool,
}

/// Whether a real value can be negative, zero and positive.
#[derive(Clone, Copy)]
struct Signs {
    negative: bool,
    zero: bool,
    positive: bool,
}

impl Range {
    fn unknown() -> Range {
        Range {
            real: false,
            lower: None,
            upper: None,
            nonzero: false,
        }
    }

    fn real() -> Range {
        Range {
            real: true,
            ..Range::unknown()
        }
    }

    fn point(q: Rational) -> Range {
        let nonzero = !q.is_zero();
        let bound = Bound {
            value: q,
            strict: false,
        };
        Range {
            real: true,
            lower: Some(bound.clone()),
            upper: Some(bound),
            nonzero,
        }
    }

    fn above(value: Rational, strict: bool) -> Range {
        Range {
            lower: Some(Bound { value, strict }),
            ..Range::real()
        }
    }

    fn below(value: Rational, strict: bool) -> Range {
        Range {
            upper: Some(Bound { value, strict }),
            ..Range::real()
        }
    }

    /// What both ranges say.
    fn meet(self, other: Range) -> Range {
        let tighter = |a: Option<Bound>, b: Option<Bound>, lower: bool| match (a, b) {
            (Some(a), Some(b)) => Some(match a.value.cmp(&b.value) {
                std::cmp::Ordering::Equal => Bound {
                    strict: a.strict || b.strict,
                    value: a.value,
                },
                order => {
                    if (order == std::cmp::Ordering::Greater) == lower {
                        a
                    } else {
                        b
                    }
                }
            }),
            (a, b) => a.or(b),
        };
        Range {
            real: self.real || other.real,
            lower: tighter(self.lower, other.lower, true),
            upper: tighter(self.upper, other.upper, false),
            nonzero: self.nonzero || other.nonzero,
        }
    }

    /// The range of a sum.
    fn add(&self, other: &Range) -> Range {
        if !(self.real && other.real) {
            return Range::unknown();
        }
        let sum = |a: &Option<Bound>, b: &Option<Bound>| {
            Some(Bound {
                value: &a.as_ref()?.value + &b.as_ref()?.value,
                strict: a.as_ref()?.strict || b.as_ref()?.strict,
            })
        };
        Range {
            real: true,
            lower: sum(&self.lower, &other.lower),
            upper: sum(&self.upper, &other.upper),
            nonzero: false,
        }
    }

    /// The range of `c` times a value in this range.
    fn scale(&self, c: &Rational) -> Range {
        if c.is_zero() {
            return Range::point(Rational::from(0));
        }
        if !self.real {
            return Range {
                nonzero: self.nonzero,
                ..Range::unknown()
            };
        }
        let times = |b: &Option<Bound>| {
            b.as_ref().map(|b| Bound {
                value: &b.value * c,
                strict: b.strict,
            })
        };
        let (lower, upper) = if c.is_positive() {
            (times(&self.lower), times(&self.upper))
        } else {
            (times(&self.upper), times(&self.lower))
        };
        Range {
            real: true,
            lower,
            upper,
            nonzero: self.nonzero,
        }
    }

    /// Whether the value cannot be 0.
    fn excludes_zero(&self) -> bool {
        let positive = |b: &Option<Bound>| {
            b.as_ref()
                .is_some_and(|b| b.value.is_positive() || b.value.is_zero() && b.strict)
        };
        let negative = |b: &Option<Bound>| {
            b.as_ref()
                .is_some_and(|b| b.value.is_negative() || b.value.is_zero() && b.strict)
        };
        self.nonzero || self.real && (positive(&self.lower) || negative(&self.upper))
    }

    fn signs(&self) -> Signs {
        let zero_or_more = |b: &Option<Bound>| b.as_ref().is_some_and(|b| !b.value.is_negative());
        let zero_or_less = |b: &Option<Bound>| b.as_ref().is_some_and(|b| !b.value.is_positive());
        let excludes_zero = self.excludes_zero();
        Signs {
            negative: !zero_or_more(&self.lower),
            zero: !excludes_zero
                && !self.upper.as_ref().is_some_and(|b| b.value.is_negative())
                && !self.lower.as_ref().is_some_and(|b| b.value.is_positive()),
            positive: !zero_or_less(&self.upper),
        }
    }

    /// The real range of a value whose signs are `s`.
    fn of_signs(s: Signs) -> Range {
        let zero = Rational::from(0);
        let mut r = match (s.negative, s.positive) {
            (false, false) => Range::point(zero),
            (false, true) => Range::above(zero, !s.zero),
            (true, false) => Range::below(zero, !s.zero),
            (true, true) => Range::real(),
        };
        r.nonzero |= !s.zero;
        r
    }

    /// The sign of a value in this range.
    pub(crate) fn sign(&self) -> Sign {
        if !self.real {
            return if self.nonzero {
                Sign::Nonzero
            } else {
                Sign::Unknown
            };
        }
        let s = self.signs();
        match (s.negative, s.zero, s.positive) {
            (false, false, true) => Sign::Positive,
            (true, false, false) => Sign::Negative,
            (false, true, false) => Sign::Zero,
            (false, true, true) => Sign::NonNegative,
            (true, true, false) => Sign::NonPositive,
            (true, false, true) => Sign::Nonzero,
            _ => Sign::Unknown,
        }
    }

    /// Which of `a < b`, `a = b` and `a > b` are possible, for `a - b` in
    /// this range: `None` when the value need not be real.
    fn order(&self) -> Option<Signs> {
        self.real.then(|| self.signs())
    }
}

/// The signs of a product of reals with these signs.
fn product_signs(a: Signs, b: Signs) -> Signs {
    Signs {
        negative: a.negative && b.positive || a.positive && b.negative,
        zero: a.zero || b.zero,
        positive: a.positive && b.positive || a.negative && b.negative,
    }
}

/// `e` as `scale·subject + shift`, with the subject monic: `None` when `e` is
/// a rational. The subject of a sum is the sum of its terms that are not
/// numbers, divided by the coefficient of the first; of a product with a
/// numeric coefficient, the product without it.
fn linear(e: &Expr) -> (Rational, Rational, Option<Expr>) {
    if let Some(q) = e.to_rational() {
        return (q, Rational::from(0), None);
    }
    let one = Rational::from(1);
    if let View::Apply { head, args } = e.view() {
        match head.as_symbol() {
            Some(Symbol::PLUS) => {
                let mut shift = Rational::from(0);
                let mut terms = Vec::new();
                for a in args {
                    match a.to_rational() {
                        Some(q) => shift = &shift + &q,
                        None if a.is_number() => return (Rational::from(0), one, Some(e.clone())),
                        None => terms.push(coefficient(a)),
                    }
                }
                let Some(first) = terms.first().map(|(c, _)| c.clone()) else {
                    return (shift, Rational::from(0), None);
                };
                let subject = if let [(_, monomial)] = terms.as_slice() {
                    monomial.clone()
                } else {
                    let scaled: Vec<Expr> =
                        terms.iter().map(|(c, m)| times(&(c / &first), m)).collect();
                    Expr::apply(&Symbol::PLUS.into(), &scaled)
                };
                return (shift, first, Some(subject));
            }
            Some(Symbol::TIMES) => {
                let (c, monomial) = coefficient(e);
                return (Rational::from(0), c, Some(monomial));
            }
            _ => {}
        }
    }
    (Rational::from(0), one, Some(e.clone()))
}

/// A term as its rational coefficient and the rest: `3 x y` is `(3, x y)`.
fn coefficient(t: &Expr) -> (Rational, Expr) {
    if let View::Apply { head, args } = t.view() {
        if head.as_symbol() == Some(Symbol::TIMES) {
            if let Some(c) = args.first().and_then(Expr::to_rational) {
                let rest = &args[1..];
                let monomial = match rest {
                    [one] => one.clone(),
                    _ => Expr::apply(head, rest),
                };
                return (c, monomial);
            }
        }
    }
    (Rational::from(1), t.clone())
}

/// `c` times the monomial `m`, in canonical form.
fn times(c: &Rational, m: &Expr) -> Expr {
    if c.is_one() {
        return m.clone();
    }
    let mut factors = vec![Expr::from(c.clone())];
    match m.view() {
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::TIMES) => {
            factors.extend_from_slice(args);
        }
        _ => factors.push(m.clone()),
    }
    Expr::apply(&Symbol::TIMES.into(), &factors)
}

/// The facts of a context's assumptions.
#[derive(Clone, Default, Debug)]
pub(crate) struct Facts {
    /// The range of each subject.
    ranges: HashMap<Expr, Range, FixedState>,
    /// Values each subject is not equal to.
    excluded: HashMap<Expr, Vec<Rational>, FixedState>,
    /// The smallest known number set of each expression.
    sets: HashMap<Expr, Set, FixedState>,
    /// Assumptions about nothing above, true as written.
    holds: Vec<Expr>,
}

/// A relation between the two sides of a comparison.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Relation {
    Less,
    LessEqual,
    Equal,
    Unequal,
    GreaterEqual,
    Greater,
}

impl Relation {
    fn of(s: Symbol) -> Option<Relation> {
        Some(match s {
            Symbol::LESS => Relation::Less,
            Symbol::LESS_EQUAL => Relation::LessEqual,
            Symbol::EQUAL => Relation::Equal,
            Symbol::UNEQUAL => Relation::Unequal,
            Symbol::GREATER_EQUAL => Relation::GreaterEqual,
            Symbol::GREATER => Relation::Greater,
            _ => return None,
        })
    }

    /// The relation with its sides swapped.
    fn flip(self) -> Relation {
        match self {
            Relation::Less => Relation::Greater,
            Relation::LessEqual => Relation::GreaterEqual,
            Relation::GreaterEqual => Relation::LessEqual,
            Relation::Greater => Relation::Less,
            r => r,
        }
    }

    /// The relation that holds when this one does not, between reals.
    fn negate(self) -> Relation {
        match self {
            Relation::Less => Relation::GreaterEqual,
            Relation::LessEqual => Relation::Greater,
            Relation::Equal => Relation::Unequal,
            Relation::Unequal => Relation::Equal,
            Relation::GreaterEqual => Relation::Less,
            Relation::Greater => Relation::LessEqual,
        }
    }

    /// Whether the relation holds of `a - b`, given which of `a < b`,
    /// `a = b` and `a > b` are possible.
    fn holds(self, possible: Signs) -> Truth {
        let accepted = match self {
            Relation::Less => (true, false, false),
            Relation::LessEqual => (true, true, false),
            Relation::Equal => (false, true, false),
            Relation::Unequal => (true, false, true),
            Relation::GreaterEqual => (false, true, true),
            Relation::Greater => (false, false, true),
        };
        let cases = [
            (possible.negative, accepted.0),
            (possible.zero, accepted.1),
            (possible.positive, accepted.2),
        ];
        if cases.iter().all(|&(p, a)| !p || a) {
            Truth::True
        } else if cases.iter().all(|&(p, a)| !p || !a) {
            Truth::False
        } else {
            Truth::Unknown
        }
    }
}

impl Facts {
    /// Reads `assumption` into facts, evaluating the differences of
    /// comparisons in `ctx`, the context before it.
    pub(crate) fn learn(&mut self, ctx: &Context, assumption: &Expr) {
        let (Some(head), args) = (
            assumption.head().and_then(Expr::as_symbol),
            assumption.args(),
        ) else {
            self.holds.push(assumption.clone());
            return;
        };
        let zero = Expr::from(0);
        let sign = |relation| (assumption, relation);
        match (head, args) {
            (Symbol::AND, _) => {
                for part in args {
                    self.learn(ctx, part);
                }
            }
            (Symbol::POSITIVE, [value]) => self.relate(ctx, sign(Relation::Greater), value, &zero),
            (Symbol::NEGATIVE, [value]) => self.relate(ctx, sign(Relation::Less), value, &zero),
            (Symbol::NON_NEGATIVE, [value]) => {
                self.relate(ctx, sign(Relation::GreaterEqual), value, &zero);
            }
            (Symbol::NON_POSITIVE, [value]) => {
                self.relate(ctx, sign(Relation::LessEqual), value, &zero);
            }
            (Symbol::ELEMENT, [value, set]) if Set::of(set).is_some() => {
                let set = Set::of(set).expect("a set");
                let value = ctx.eval(value).unwrap_or_else(|_| value.clone());
                let known = self.sets.entry(value.clone()).or_insert(set);
                *known = (*known).min(set);
                if set <= Set::Reals {
                    if let (_, _, Some(subject)) = linear(&value) {
                        self.meet(subject, Range::real());
                    }
                }
            }
            (Symbol::NOT, [negated]) => {
                match (negated.head().and_then(Expr::as_symbol), negated.args()) {
                    // Between reals, which comparisons are about.
                    (Some(compare), [left, right]) if Relation::of(compare).is_some() => {
                        let relation = Relation::of(compare).expect("a relation").negate();
                        self.relate(ctx, sign(relation), left, right);
                    }
                    _ => self.holds.push(assumption.clone()),
                }
            }
            (compare, _) if Relation::of(compare).is_some() && args.len() >= 2 => {
                let relation = Relation::of(compare).expect("a relation");
                for pair in args.windows(2) {
                    self.relate(ctx, sign(relation), &pair[0], &pair[1]);
                }
            }
            _ => self.holds.push(assumption.clone()),
        }
    }

    /// Learns `left relation right`, from `assumption`.
    fn relate(
        &mut self,
        ctx: &Context,
        (assumption, relation): (&Expr, Relation),
        left: &Expr,
        right: &Expr,
    ) {
        let Ok(diff) = ctx.eval(&difference(left, right)) else {
            self.holds.push(assumption.clone());
            return;
        };
        // `scale·subject + shift relation 0`: a bound on the subject.
        let (shift, scale, Some(subject)) = linear(&diff) else {
            return;
        };
        let value = -(&shift / &scale);
        let relation = if scale.is_negative() {
            relation.flip()
        } else {
            relation
        };
        let range = match relation {
            Relation::Less => Range::below(value, true),
            Relation::LessEqual => Range::below(value, false),
            Relation::Equal => Range::point(value),
            Relation::GreaterEqual => Range::above(value, false),
            Relation::Greater => Range::above(value, true),
            Relation::Unequal => {
                self.excluded.entry(subject).or_default().push(value);
                return;
            }
        };
        self.meet(subject, range);
    }

    fn meet(&mut self, subject: Expr, range: Range) {
        match self.ranges.remove(&subject) {
            Some(known) => self.ranges.insert(subject, known.meet(range)),
            None => self.ranges.insert(subject, range),
        };
    }
}

/// The range of `e`, which is evaluated, from the facts of `ctx` and its
/// structure.
pub(crate) fn range(ctx: &Context, e: &Expr) -> Range {
    let (shift, scale, subject) = linear(e);
    match subject {
        Some(subject) => linear_range(ctx, &shift, &scale, &subject),
        None => Range::point(shift),
    }
}

/// The range of `scale·subject + shift`.
fn linear_range(ctx: &Context, shift: &Rational, scale: &Rational, subject: &Expr) -> Range {
    let (shift, scale, subject) = (shift.clone(), scale.clone(), subject.clone());
    if scale.is_zero() {
        return Range::point(shift);
    }
    let facts = ctx.facts();
    let mut r = structure(ctx, &subject);
    if let Some(known) = facts.ranges.get(&subject) {
        r = r.meet(known.clone());
    }
    // `e` is 0 where the subject is `-shift/scale`.
    let zero_at = -(&shift / &scale);
    let excluded = facts
        .excluded
        .get(&subject)
        .is_some_and(|values| values.contains(&zero_at));
    if !scale.is_one() {
        r = r.scale(&scale);
    }
    if !shift.is_zero() {
        // What is not 0 may be after the shift: only the bounds and the
        // excluded values tell.
        r = r.add(&Range::point(shift));
    }
    if excluded {
        r.nonzero = true;
    }
    r
}

/// The range of a ball: its closed interval, exactly (CAS-DOC-5 §9.3);
/// some real for one too large for exact ends.
fn ball_range(b: &Ball) -> Range {
    let ends = b
        .mid()
        .checked_to_rational()
        .zip(b.rad().to_float().and_then(|r| r.checked_to_rational()));
    let Some((mid, rad)) = ends else {
        return Range::real();
    };
    let (lo, hi) = (&mid - &rad, &mid + &rad);
    let zero = Rational::from(0);
    Range {
        real: true,
        nonzero: lo > zero || hi < zero,
        lower: Some(Bound {
            value: lo,
            strict: false,
        }),
        upper: Some(Bound {
            value: hi,
            strict: false,
        }),
    }
}

/// The range of a subject from its structure alone.
fn structure(ctx: &Context, e: &Expr) -> Range {
    match e.view() {
        View::Ball(b) => return ball_range(b),
        View::ComplexBall(z) => {
            return Range {
                nonzero: !z.may_be_zero(),
                ..Range::unknown()
            };
        }
        _ => {}
    }
    if e.is_number() {
        return match e.to_rational() {
            Some(q) => Range::point(q),
            None => Range {
                nonzero: true,
                ..Range::unknown()
            },
        };
    }
    match e.view() {
        View::Symbol(s) => symbol(ctx, s),
        View::Apply { head, args } => match (head.as_symbol(), args) {
            (Some(Symbol::PLUS), _) => {
                let mut sum = Range::point(Rational::from(0));
                for a in args {
                    sum = sum.add(&range(ctx, a));
                    if !sum.real {
                        break;
                    }
                }
                sum
            }
            (Some(Symbol::TIMES), _) => {
                let mut signs = Signs {
                    negative: false,
                    zero: false,
                    positive: true,
                };
                let (mut real, mut nonzero) = (true, true);
                for a in args {
                    let r = range(ctx, a);
                    nonzero &= r.excludes_zero();
                    if r.real {
                        signs = product_signs(signs, r.signs());
                    } else {
                        real = false;
                    }
                }
                if real {
                    Range::of_signs(signs)
                } else {
                    Range {
                        nonzero,
                        ..Range::unknown()
                    }
                }
            }
            (Some(Symbol::POWER), [b, n]) => power(ctx, b, n),
            (Some(Symbol::ABS), [x]) => {
                let mut r = Range::above(Rational::from(0), false);
                r.nonzero = range(ctx, x).excludes_zero();
                r
            }
            (Some(Symbol::EXP), [x]) if range(ctx, x).real => Range::above(Rational::from(0), true),
            _ => Range::unknown(),
        },
        _ => Range::unknown(),
    }
}

/// The range of a symbol: the builtin constants lie between rationals, and
/// in the Real domain any other symbol is real.
fn symbol(ctx: &Context, s: Symbol) -> Range {
    let between = |a: (i64, i64), b: (i64, i64)| {
        let q = |(n, d): (i64, i64)| Rational::new(Integer::from(n), Integer::from(d));
        Range::above(q(a), true).meet(Range::below(q(b), true))
    };
    match s {
        Symbol::PI => between((157, 50), (63, 20)),
        Symbol::E => between((271, 100), (68, 25)),
        _ if s.is_builtin() => Range::unknown(),
        _ => match ctx.domain() {
            Domain::Real => Range::real(),
            Domain::Complex => Range::unknown(),
        },
    }
}

/// The range of `b^n`.
fn power(ctx: &Context, b: &Expr, n: &Expr) -> Range {
    let Some(n) = n.to_rational() else {
        return Range::unknown();
    };
    let base = range(ctx, b);
    if !base.real {
        // A power of what is not 0 is not 0.
        return Range {
            nonzero: base.excludes_zero(),
            ..Range::unknown()
        };
    }
    let s = base.signs();
    // `0` to a negative power is no real number.
    if n.is_negative() && s.zero {
        return Range::unknown();
    }
    // A real root, when there is one: the Real domain takes the real odd
    // root of a negative (O7).
    let odd_root = !even(n.denom()) && ctx.domain() == Domain::Real;
    if !n.is_integer() && s.negative && !odd_root {
        return Range::unknown();
    }
    // The sign of `b^(p/q)` is that of `b^p`.
    Range::of_signs(if even(n.numer()) {
        Signs {
            negative: false,
            zero: s.zero,
            positive: s.negative || s.positive,
        }
    } else {
        s
    })
}

fn even(k: &Integer) -> bool {
    // 0 has no lowest set bit, and is even.
    k.unsigned_abs().trailing_zeros().is_none_or(|t| t > 0)
}

/// The smallest known number set of `e`, which is evaluated.
pub(crate) fn set(ctx: &Context, e: &Expr) -> Option<Set> {
    if let Some(s) = Set::of_number(e) {
        return Some(s);
    }
    let known = ctx.facts().sets.get(e).copied();
    let structural = match e.view() {
        View::Symbol(Symbol::PI | Symbol::E) => Some(Set::Reals),
        View::Symbol(s) if !s.is_builtin() => Some(match ctx.domain() {
            Domain::Real => Set::Reals,
            Domain::Complex => Set::Complexes,
        }),
        View::Apply { head, args } => match (head.as_symbol(), args) {
            (Some(Symbol::PLUS | Symbol::TIMES), _) => {
                let mut largest = Set::Integers;
                for a in args {
                    largest = largest.max(set(ctx, a)?);
                }
                Some(largest)
            }
            (Some(Symbol::POWER), [b, n]) => power_set(ctx, b, n),
            (Some(Symbol::ABS), [x]) => set(ctx, x).map(|s| s.min(Set::Reals)),
            _ => None,
        },
        _ => None,
    };
    let from_range = range(ctx, e).real.then_some(Set::Reals);
    [known, structural, from_range].into_iter().flatten().min()
}

fn power_set(ctx: &Context, b: &Expr, n: &Expr) -> Option<Set> {
    let base = set(ctx, b)?;
    let n = n.to_rational()?;
    if n.is_integer() {
        return Some(if n.is_negative() && base == Set::Integers {
            Set::Rationals
        } else {
            base
        });
    }
    Some(
        if range(ctx, b).real
            && range(
                ctx,
                &Expr::apply(&Symbol::POWER.into(), &[b.clone(), Expr::from(n)]),
            )
            .real
        {
            Set::Reals
        } else {
            Set::Complexes
        },
    )
}

/// The sign of `e`, which is evaluated, as far as `ctx` tells.
pub(crate) fn sign_of(ctx: &Context, e: &Expr) -> Sign {
    match e.view() {
        View::Symbol(Symbol::INFINITY) => return Sign::Positive,
        View::Apply { head, args: [d] } if head.as_symbol() == Some(Symbol::DIRECTED_INFINITY) => {
            return match d.to_rational() {
                Some(q) if q.is_positive() => Sign::Positive,
                Some(q) if q.is_negative() => Sign::Negative,
                _ => Sign::Unknown,
            };
        }
        _ => {}
    }
    range(ctx, e).sign()
}

/// `a - b`, unevaluated, the minus distributed over the terms of `b`, so
/// that like terms cancel when it is evaluated.
fn difference(a: &Expr, b: &Expr) -> Expr {
    let minus = |t: &Expr| Expr::apply(&Symbol::TIMES.into(), &[Expr::from(-1), t.clone()]);
    let mut terms = vec![a.clone()];
    match b.view() {
        View::Apply { head, args } if head.as_symbol() == Some(Symbol::PLUS) => {
            terms.extend(args.iter().map(minus));
        }
        _ => terms.push(minus(b)),
    }
    Expr::apply(&Symbol::PLUS.into(), &terms)
}

/// Whether `a r b` holds, for evaluated `a` and `b`: by the range of their
/// difference, evaluated if neither is a rational.
fn compare(ev: &mut Eval<'_>, a: &Expr, b: &Expr, r: Relation) -> Result<Truth, EvalError> {
    let ctx = ev.context();
    let difference = match (a.to_rational(), b.to_rational()) {
        // `a - q` and `q - b` are linear forms of `a` and `b`.
        (_, Some(q)) => match linear(a) {
            (shift, scale, Some(subject)) => linear_range(ctx, &(&shift - &q), &scale, &subject),
            (shift, _, None) => Range::point(&shift - &q),
        },
        (Some(q), _) => match linear(b) {
            (shift, scale, Some(subject)) => linear_range(ctx, &(&q - &shift), &-scale, &subject),
            (shift, _, None) => Range::point(&q - &shift),
        },
        (None, None) => {
            let d = ev.eval(&difference(a, b))?;
            range(ev.context(), &d)
        }
    };
    Ok(match difference.order() {
        Some(possible) => r.holds(possible),
        None => Truth::Unknown,
    })
}

/// Whether `a r b` holds, from the comparison natives, as [`compare`].
pub(crate) fn relation(
    ev: &mut Eval<'_>,
    head: Symbol,
    a: &Expr,
    b: &Expr,
) -> Result<Truth, EvalError> {
    match Relation::of(head) {
        Some(r) => compare(ev, a, b, r),
        None => Ok(Truth::Unknown),
    }
}

/// Whether `e` is in `s`: true when its known set is in `s`, false for a
/// number that is not.
pub(crate) fn element(ctx: &Context, e: &Expr, s: Set) -> Truth {
    if let Some(n) = Set::of_number(e) {
        return Truth::of(n <= s);
    }
    match set(ctx, e) {
        Some(known) if known <= s => Truth::True,
        _ => Truth::Unknown,
    }
}

/// Whether `prop` holds, evaluating its parts as needed.
fn decide_in(ev: &mut Eval<'_>, prop: &Expr) -> Result<Truth, EvalError> {
    let (Some(head), args) = (prop.head().and_then(Expr::as_symbol), prop.args()) else {
        return Ok(match prop.as_symbol() {
            Some(Symbol::TRUE) => Truth::True,
            Some(Symbol::FALSE) => Truth::False,
            _ => holds(ev.context(), prop),
        });
    };
    match (head, args) {
        (Symbol::NOT, [negated]) => Ok(decide_in(ev, negated)?.not()),
        (Symbol::AND | Symbol::OR, _) => {
            let stop = if head == Symbol::AND {
                Truth::False
            } else {
                Truth::True
            };
            let mut all = true;
            for part in args {
                match decide_in(ev, part)? {
                    truth if truth == stop => return Ok(stop),
                    Truth::Unknown => all = false,
                    _ => {}
                }
            }
            Ok(if all { stop.not() } else { Truth::Unknown })
        }
        (
            Symbol::POSITIVE | Symbol::NEGATIVE | Symbol::NON_NEGATIVE | Symbol::NON_POSITIVE,
            [value],
        ) => {
            let value = ev.eval(value)?;
            Ok(sign_test(
                head,
                sign_of(ev.context(), &value),
                value.is_exact_number(),
            ))
        }
        (Symbol::ELEMENT, [value, set]) => {
            let Some(set) = Set::of(set) else {
                return Ok(holds(ev.context(), prop));
            };
            let value = ev.eval(value)?;
            Ok(element(ev.context(), &value, set))
        }
        (compare_head, _) if Relation::of(compare_head).is_some() && args.len() >= 2 => {
            let relation = Relation::of(compare_head).expect("a relation");
            let mut all = true;
            for pair in args.windows(2) {
                let (left, right) = (ev.eval(&pair[0])?, ev.eval(&pair[1])?);
                match compare(ev, &left, &right, relation)? {
                    Truth::False => return Ok(Truth::False),
                    Truth::Unknown => all = false,
                    Truth::True => {}
                }
            }
            if all {
                Ok(Truth::True)
            } else {
                Ok(holds(ev.context(), prop))
            }
        }
        _ => Ok(holds(ev.context(), prop)),
    }
}

/// Whether `p` is an assumption, or its negation is.
fn holds(ctx: &Context, p: &Expr) -> Truth {
    let facts = ctx.facts();
    if facts.holds.contains(p) {
        return Truth::True;
    }
    let negation = Expr::apply(&Symbol::NOT.into(), std::slice::from_ref(p));
    if facts.holds.contains(&negation) {
        Truth::False
    } else {
        Truth::Unknown
    }
}

/// The answer of a sign predicate, `Positive` to `NonPositive`, for a value
/// with this sign; `number` when it is an exact number, which a non-real one
/// fails (a ball of unknown sign may be anything).
pub(crate) fn sign_test(head: Symbol, sign: Sign, number: bool) -> Truth {
    use Sign as S;
    let (yes, no): (&[Sign], &[Sign]) = match head {
        Symbol::POSITIVE => (&[S::Positive], &[S::Negative, S::Zero, S::NonPositive]),
        Symbol::NEGATIVE => (&[S::Negative], &[S::Positive, S::Zero, S::NonNegative]),
        Symbol::NON_NEGATIVE => (&[S::Positive, S::Zero, S::NonNegative], &[S::Negative]),
        Symbol::NON_POSITIVE => (&[S::Negative, S::Zero, S::NonPositive], &[S::Positive]),
        _ => return Truth::Unknown,
    };
    if yes.contains(&sign) {
        Truth::True
    } else if no.contains(&sign) || number && matches!(sign, S::Nonzero | S::Unknown) {
        Truth::False
    } else {
        Truth::Unknown
    }
}

impl Context {
    /// What symbols without facts range over.
    #[must_use]
    pub fn domain(&self) -> Domain {
        self.number_domain()
    }

    /// Whether the proposition `prop` holds under this context's assumptions
    /// and domain (CAS-DOC-4 §7): comparisons, `Positive`, `Negative`,
    /// `NonNegative`, `NonPositive`, `Element`, and `And`, `Or` and `Not` of
    /// them. Its parts are evaluated; `prop` itself is not.
    ///
    /// ```
    /// use cas_core::{Context, Expr, Truth};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let ctx = Context::new().assuming(read("Greater[x, 2]"));
    /// assert_eq!(ctx.decide(&read("Greater[Plus[x, -1], 0]")), Truth::True);
    /// assert_eq!(ctx.decide(&read("Less[x, 0]")), Truth::False);
    /// assert_eq!(ctx.decide(&read("Greater[x, 3]")), Truth::Unknown);
    /// ```
    #[must_use]
    pub fn decide(&self, prop: &Expr) -> Truth {
        run(self.limits(), |recursion| {
            decide_in(&mut Eval::new(self, recursion), prop).unwrap_or(Truth::Unknown)
        })
    }

    /// The sign of `e`, evaluated first, under this context's assumptions
    /// and domain.
    ///
    /// ```
    /// use cas_core::{Context, Expr, Sign};
    ///
    /// let read = |s: &str| Expr::parse_full_form(s).unwrap();
    /// let ctx = Context::new();
    /// // In the Real domain, `x` is real.
    /// assert_eq!(ctx.sign(&read("Plus[Power[x, 2], 1]")), Sign::Positive);
    /// ```
    #[must_use]
    pub fn sign(&self, e: &Expr) -> Sign {
        match self.eval(e) {
            Ok(e) => sign_of(self, &e),
            Err(_) => Sign::Unknown,
        }
    }
}

#[cfg(test)]
mod tests;
