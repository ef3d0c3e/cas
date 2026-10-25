//! Gruntz's algorithm (CAS-DOC-10 §5, CAS-206) [Gru96]: limits of exp-log
//! expressions at `x → +∞` where a series has an essential singularity, as
//! `exp(x)/x^10`.
//!
//! - The most rapidly varying subexpressions `Ω` of `f`: `x`, and the
//!   exponentials `exp h` whose `h` tends to `±∞`, compared through their
//!   logarithms, `lim log a / log b` 0 (`b` faster), infinite (`a` faster)
//!   or finite and not 0 (the same class); `mrv` of a sum, product, power by
//!   a constant, logarithm, or function of a subexpression is the fastest of
//!   its parts'. A power `b^c` with `x` in `c` is `exp(c log b)` first.
//! - With `x` in `Ω`, `f(e^x)` instead: one level up.
//! - Otherwise `ω = exp(s)` for an element `exp(h)` of `Ω`, `s = ∓h` so
//!   that `s → −∞` and `ω → 0⁺`; each `exp(h_i)` of `Ω` is
//!   `ω^{c_i} exp(h_i − c_i s)` with `c_i = lim h_i / s`; the series of `f`
//!   in `ω` (CAS-203), the slower subexpressions in its coefficients and
//!   `log ω = s`, has a leading term `c₀ ω^{e₀}`: 0 for `e₀ > 0`, an
//!   infinity of `c₀`'s sign for `e₀ < 0`, and the limit of `c₀`, slower,
//!   for `e₀ = 0`.
//! - The sign of `f` for large `x`: that of a constant, `+` for `x` and
//!   exponentials, products of signs, and otherwise that of `c₀`.
//!
//! Each nested limit takes one level of [`MAX_DEPTH`]; deeper, or where a
//! subexpression is not exp-log or a series fails, nothing is decided.

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::Value;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::native::series::{Expander, expansion, log_t};
use crate::symbol::Symbol;

type Found<T> = Result<Option<T>, EvalError>;

/// How deep the nested limits go.
const MAX_DEPTH: u32 = 12;

/// The orders tried for a leading term in `ω`.
const ORDERS: [i64; 4] = [2, 4, 8, 16];

/// The variable `ω` at `depth`, one a depth.
fn omega(depth: u32) -> Expr {
    Expr::from(Symbol::new(&format!("Gruntz$w{depth}")))
}

/// What Gruntz's algorithm works on: `x`, and the evaluator.
pub(super) struct Gruntz<'a, 'e> {
    ev: &'a mut Eval<'e>,
    x: Expr,
}

impl<'a, 'e> Gruntz<'a, 'e> {
    pub(super) fn new(ev: &'a mut Eval<'e>, x: Expr) -> Gruntz<'a, 'e> {
        Gruntz { ev, x }
    }

    /// `e` with its powers by an exponent in `x` as exponentials, and
    /// `E^h` as `Exp[h]`.
    fn exponentials(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        let x = self.x.clone();
        let rewritten =
            e.replace(
                |part| match (part.head().and_then(Expr::as_symbol), part.args()) {
                    (Some(Symbol::POWER), [b, c]) if b.as_symbol() == Some(Symbol::E) => {
                        Some(Expr::call(Symbol::EXP, [c.clone()]))
                    }
                    (Some(Symbol::POWER), [b, c]) if c.contains(&x) => Some(Expr::call(
                        Symbol::EXP,
                        [c.clone() * Expr::call(Symbol::LOG, [b.clone()])],
                    )),
                    _ => None,
                },
            );
        // Evaluated, but the exponentials kept: `Exp` evaluates to itself.
        self.ev.eval(&rewritten)
    }

    /// The limit of `e` at `x → +∞`.
    pub(super) fn limit(&mut self, e: &Expr, depth: u32) -> Found<Value> {
        if depth >= MAX_DEPTH {
            return Ok(None);
        }
        if !e.contains(&self.x) {
            return Ok(Some(Value::Finite(e.clone())));
        }
        if *e == self.x {
            return Ok(Some(Value::Infinite(true)));
        }
        let e = self.exponentials(e)?;
        let Some(omega_set) = self.mrv(&e, depth)? else {
            return Ok(None);
        };
        if omega_set.is_empty() {
            return Ok(None);
        }
        if omega_set.contains(&self.x) {
            let up = Expr::call(Symbol::EXP, [self.x.clone()]);
            let e = self.ev.eval(&e.substitute(&[(self.x.clone(), up)]))?;
            return self.limit(&e, depth + 1);
        }
        let Some((c0, e0)) = self.leading(&e, &omega_set, depth)? else {
            return Ok(None);
        };
        match e0.cmp(&Rational::ZERO) {
            Ordering::Greater => Ok(Some(Value::Finite(Expr::from(0)))),
            Ordering::Less => Ok(Some(match self.sign(&c0, depth + 1)? {
                Some(Ordering::Greater) => Value::Infinite(true),
                Some(Ordering::Less) => Value::Infinite(false),
                _ => Value::Unsigned,
            })),
            Ordering::Equal => self.limit(&c0, depth + 1),
        }
    }

    /// The most rapidly varying subexpressions of `e`; `None` where it is
    /// not exp-log or a comparison is not decided.
    fn mrv(&mut self, e: &Expr, depth: u32) -> Found<Vec<Expr>> {
        if !e.contains(&self.x) {
            return Ok(Some(Vec::new()));
        }
        if *e == self.x {
            return Ok(Some(vec![self.x.clone()]));
        }
        let Some(head) = e.head().and_then(Expr::as_symbol) else {
            return Ok(None);
        };
        match (head, e.args()) {
            (Symbol::PLUS | Symbol::TIMES, args) => {
                let mut fastest = Vec::new();
                for a in args {
                    let Some(m) = self.mrv(a, depth)? else {
                        return Ok(None);
                    };
                    let Some(merged) = self.max(fastest, m, depth)? else {
                        return Ok(None);
                    };
                    fastest = merged;
                }
                Ok(Some(fastest))
            }
            (Symbol::POWER, [b, c]) if !c.contains(&self.x) => self.mrv(b, depth),
            (Symbol::EXP, [h]) => {
                let Some(inner) = self.mrv(h, depth)? else {
                    return Ok(None);
                };
                match self.limit(h, depth + 1)? {
                    Some(Value::Infinite(_)) => self.max(vec![e.clone()], inner, depth),
                    Some(Value::Finite(_)) => Ok(Some(inner)),
                    _ => Ok(None),
                }
            }
            (
                Symbol::LOG
                | Symbol::SIN
                | Symbol::COS
                | Symbol::TAN
                | Symbol::ARC_TAN
                | Symbol::ARC_SIN
                | Symbol::ARC_COS
                | Symbol::ABS,
                [g],
            ) => self.mrv(g, depth),
            _ => Ok(None),
        }
    }

    /// The faster of two classes, or both when they are the same.
    fn max(&mut self, a: Vec<Expr>, b: Vec<Expr>, depth: u32) -> Found<Vec<Expr>> {
        let (Some(first), Some(second)) = (a.first().cloned(), b.first().cloned()) else {
            return Ok(Some(if a.is_empty() { b } else { a }));
        };
        Ok(match self.compare(&first, &second, depth)? {
            Some(Ordering::Greater) => Some(a),
            Some(Ordering::Less) => Some(b),
            Some(Ordering::Equal) => {
                let mut union = a;
                for e in b {
                    if !union.contains(&e) {
                        union.push(e);
                    }
                }
                Some(union)
            }
            None => None,
        })
    }

    /// `log f` for an element of a class: `log x`, or `h` for `exp h`.
    fn logarithm(f: &Expr) -> Expr {
        match f.args() {
            [h] if f.head() == Some(&Expr::from(Symbol::EXP)) => h.clone(),
            _ => Expr::call(Symbol::LOG, [f.clone()]),
        }
    }

    /// How `a` grows against `b`: `lim log a / log b`.
    fn compare(&mut self, a: &Expr, b: &Expr, depth: u32) -> Found<Ordering> {
        let ratio = self.ev.eval(&(Self::logarithm(a) / Self::logarithm(b)))?;
        Ok(match self.limit(&ratio, depth + 1)? {
            Some(Value::Finite(c)) if c == Expr::from(0) => Some(Ordering::Less),
            Some(Value::Finite(_)) => Some(Ordering::Equal),
            Some(Value::Infinite(_) | Value::Unsigned) => Some(Ordering::Greater),
            None => None,
        })
    }

    /// The leading term `c₀ ω^{e₀}` of `e` in the `ω` of its class `Ω`.
    fn leading(&mut self, e: &Expr, omega_set: &[Expr], depth: u32) -> Found<(Expr, Rational)> {
        // ω = exp(s), s → −∞, from the smallest element, which no other is
        // inside.
        let generator = omega_set
            .iter()
            .min_by_key(|f| f.size())
            .expect("a nonempty class");
        let first = Self::logarithm(generator);
        let s = match self.limit(&first, depth + 1)? {
            Some(Value::Infinite(true)) => self.ev.eval(&(-first))?,
            Some(Value::Infinite(false)) => first,
            _ => return Ok(None),
        };
        let w = omega(depth);
        let mut rules = Vec::with_capacity(omega_set.len());
        for f in omega_set {
            let h = Self::logarithm(f);
            let ratio = self.ev.eval(&(h.clone() / s.clone()))?;
            let Some(Value::Finite(c)) = self.limit(&ratio, depth + 1)? else {
                return Ok(None);
            };
            let rest = Expr::call(Symbol::EXP, [h - c.clone() * s.clone()]);
            rules.push((f.clone(), w.pow(c) * rest));
        }
        // Elements inside others, as `exp(-x)` in `exp(x + exp(-x))`, come
        // back in the rest: rewritten again.
        let mut in_omega = self.ev.eval(&e.substitute(&rules))?;
        for _ in 0..4 {
            if Self::slower(&in_omega, omega_set) {
                break;
            }
            in_omega = self.ev.eval(&in_omega.substitute(&rules))?;
        }
        if !Self::slower(&in_omega, omega_set) {
            return Ok(None);
        }
        for order in ORDERS {
            let order = Rational::from(Integer::from(order));
            let Some(series) = expansion(self.ev, &in_omega, &w, &order)? else {
                return Ok(None);
            };
            if let Some(c0) = series.coeffs.first() {
                // log ω = s.
                let c0 = self.ev.eval(&c0.substitute(&[(log_t(), s.clone())]))?;
                let e0 = Rational::new(
                    Integer::from(series.valuation),
                    Integer::from(series.ramification),
                );
                return Ok(Some((c0, e0)));
            }
        }
        Ok(None)
    }

    /// Whether no element of `Ω` is left in `e`.
    fn slower(e: &Expr, omega_set: &[Expr]) -> bool {
        !omega_set.iter().any(|f| e.contains(f))
    }

    /// The sign of `e` for large `x`.
    fn sign(&mut self, e: &Expr, depth: u32) -> Found<Ordering> {
        if depth >= MAX_DEPTH {
            return Ok(None);
        }
        if !e.contains(&self.x) {
            let mut expander = Expander::new(self.ev, self.x.clone());
            return Ok(expander.sign(e));
        }
        if *e == self.x {
            return Ok(Some(Ordering::Greater));
        }
        let e = self.exponentials(e)?;
        match (e.head().and_then(Expr::as_symbol), e.args()) {
            (Some(Symbol::EXP), [_]) => return Ok(Some(Ordering::Greater)),
            (Some(Symbol::TIMES), args) => {
                let mut sign = Ordering::Greater;
                for a in args {
                    let Some(s) = self.sign(a, depth + 1)? else {
                        return Ok(None);
                    };
                    sign = match (sign, s) {
                        (Ordering::Equal, _) | (_, Ordering::Equal) => Ordering::Equal,
                        (p, q) if p == q => Ordering::Greater,
                        _ => Ordering::Less,
                    };
                }
                return Ok(Some(sign));
            }
            (Some(Symbol::POWER), [b, c]) if !c.contains(&self.x) => {
                if self.sign(b, depth + 1)? == Some(Ordering::Greater) {
                    return Ok(Some(Ordering::Greater));
                }
            }
            _ => {}
        }
        let Some(omega_set) = self.mrv(&e, depth)? else {
            return Ok(None);
        };
        if omega_set.is_empty() {
            return Ok(None);
        }
        if omega_set.contains(&self.x) {
            let up = Expr::call(Symbol::EXP, [self.x.clone()]);
            let e = self.ev.eval(&e.substitute(&[(self.x.clone(), up)]))?;
            return self.sign(&e, depth + 1);
        }
        let Some((c0, _)) = self.leading(&e, &omega_set, depth)? else {
            return Ok(None);
        };
        self.sign(&c0, depth + 1)
    }
}
