//! Truncated Laurent and Puiseux series with expression coefficients
//! (CAS-DOC-10 §2, CAS-203), the engine of `Series` and `Limit`.
//!
//! A series in `t` is `Σ_k c_k t^((v + k)/d) + O(t^((v + n)/d))` for its
//! `n` coefficients: a valuation `v`, a ramification `d ≥ 1`, and
//! coefficients that are evaluated expressions free of `t`, `log t` among
//! them as the symbol [`log_t`] (§4.2). Over [Knu97 §4.7]:
//!
//! - sums and products coefficientwise, a product known as far as the
//!   least of each operand's order shifted by the other's valuation;
//! - powers by any exponent free of `t`, `c^α t^(vα) (1 + h)^α`, by J. C. P.
//!   Miller's recurrence `q_k = (1/k) Σ_{j=1}^{k} ((α + 1) j − k) h_j q_{k−j}`,
//!   which makes `vα` a ramification when it is not an integer: the
//!   reciprocal is `α = −1`;
//! - `Exp`, `Log`, `Sin` and `Cos` of `g = g₀ + h` (`h` of positive
//!   valuation) by the differential equations of `exp h`, `log(1 + h)`,
//!   `sin h` and `cos h`, then `e^{g₀} exp h`, `log c + v log t + log(1 + h)`
//!   and the addition formulas; `Tan` as `sin/cos`; `ArcTan`, `ArcSin` and
//!   `ArcCos` as `f(g₀) + ∫ f′(g) g′` [BK78], the derivative and the
//!   integral in `t`; the hyperbolic functions by `e^g` and `e^(−g)`, and
//!   their inverses as `ArcTan` (CAS-DOC-15 §2);
//! - an expression tree from its leaves, each node's series from its
//!   children's to one target order; where a product with a pole or a
//!   cancellation leaves less, the caller asks again with a larger target
//!   ([`expansion`]).
//!
//! A function of a series with a pole (`exp(1/t)`), of `t` in an unknown
//! head, or a power by `t` of a series with a pole has no series here.

#![allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    reason = "the names of series algebra: g = g₀ + h, s and c for sin h and cos h, sin0 and cos0"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::approximation;
use crate::assume::{Sign, sign_of};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Found = Result<Option<Series>, EvalError>;

/// `log t` in the coefficients of a series: a symbol free of `t`, for
/// `Limit` to tell apart (CAS-DOC-10 §4.2), and for `Series` to write back.
pub(crate) fn log_t() -> Expr {
    Expr::from(Symbol::new("Series$log"))
}

/// A truncated series in `t`: `Σ c_k t^((valuation + k)/ramification)` for
/// the known coefficients, the next term unknown.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Series {
    pub(crate) valuation: i64,
    pub(crate) coeffs: Vec<Expr>,
    pub(crate) ramification: i64,
}

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

impl Series {
    /// The exponent of the first unknown term, in units of `1/d`.
    pub(crate) fn order_units(&self) -> i64 {
        self.valuation + i64::try_from(self.coeffs.len()).expect("a length in memory")
    }

    /// The exponent of the first unknown term.
    pub(crate) fn order(&self) -> Rational {
        q(self.order_units(), self.ramification)
    }

    /// The coefficient at `t^(u/d)`, `u` in units; 0 below the valuation.
    fn at(&self, u: i64) -> Option<&Expr> {
        let k = u.checked_sub(self.valuation)?;
        usize::try_from(k).ok().and_then(|k| self.coeffs.get(k))
    }

    /// The same series in units of `1/(d m)`.
    fn refined(&self, m: i64) -> Series {
        if m == 1 {
            return self.clone();
        }
        let zero = Expr::from(0);
        let mut coeffs = Vec::with_capacity(self.coeffs.len() * usize::try_from(m).expect("small"));
        for c in &self.coeffs {
            coeffs.push(c.clone());
            for _ in 1..m {
                coeffs.push(zero.clone());
            }
        }
        Series {
            valuation: self.valuation * m,
            coeffs,
            ramification: self.ramification * m,
        }
    }

    /// The units in `1/d` of an absolute exponent, rounded up.
    fn units(&self, exponent: &Rational) -> i64 {
        let scaled = exponent * &Rational::from(Integer::from(self.ramification));
        scaled.ceil().to_i64().unwrap_or(i64::MAX / 4)
    }

    /// Kept to the exponents below `order`.
    pub(crate) fn truncated(mut self, order: &Rational) -> Series {
        let units = self.units(order);
        let keep = usize::try_from((units - self.valuation).max(0)).unwrap_or(0);
        self.coeffs.truncate(keep);
        self
    }

    /// The least ramification that writes the same exponents.
    fn reduced(mut self) -> Series {
        let mut g = gcd(self.ramification, self.valuation);
        for (k, c) in self.coeffs.iter().enumerate() {
            if *c != Expr::from(0) {
                g = gcd(g, self.valuation + i64::try_from(k).expect("small"));
            }
        }
        g = gcd(g, self.order_units());
        if g <= 1 {
            return self;
        }
        let step = usize::try_from(g).expect("small");
        self.coeffs = self.coeffs.into_iter().step_by(step).collect();
        self.valuation /= g;
        self.ramification /= g;
        self
    }
}

/// The series of the expressions of one variable `t`, their coefficients
/// evaluated by `ev`.
pub(crate) struct Expander<'a, 'e> {
    ev: &'a mut Eval<'e>,
    t: Expr,
}

impl<'a, 'e> Expander<'a, 'e> {
    pub(crate) fn new(ev: &'a mut Eval<'e>, t: Expr) -> Expander<'a, 'e> {
        Expander { ev, t }
    }

    fn eval(&mut self, e: &Expr) -> Result<Expr, EvalError> {
        self.ev.eval(e)
    }

    /// Whether `c` is 0, once evaluated or put over one denominator.
    pub(crate) fn is_zero(&mut self, c: &Expr) -> Result<bool, EvalError> {
        if *c == Expr::from(0) {
            return Ok(true);
        }
        if c.is_number() {
            return Ok(false);
        }
        Ok(self.eval(&Expr::call(Symbol::TOGETHER, [c.clone()]))? == Expr::from(0))
    }

    /// `s` without its leading zero coefficients.
    pub(crate) fn normalized(&mut self, mut s: Series) -> Result<Series, EvalError> {
        let mut skip = 0;
        while skip < s.coeffs.len() && self.is_zero(&s.coeffs[skip])? {
            skip += 1;
        }
        s.coeffs.drain(..skip);
        s.valuation += i64::try_from(skip).expect("small");
        Ok(s)
    }

    /// `c` to the order `order`.
    fn constant(c: &Expr, order: &Rational) -> Series {
        let n = usize::try_from(order.ceil().to_i64().unwrap_or(0).max(0)).unwrap_or(0);
        let mut coeffs = vec![Expr::from(0); n];
        if let Some(first) = coeffs.first_mut() {
            *first = c.clone();
        }
        Series {
            valuation: 0,
            coeffs,
            ramification: 1,
        }
    }

    /// `t` to the order `order`, one term at least.
    fn variable(order: &Rational) -> Series {
        let n = usize::try_from(order.ceil().to_i64().unwrap_or(0) - 1)
            .unwrap_or(0)
            .max(1);
        let mut coeffs = vec![Expr::from(0); n];
        coeffs[0] = Expr::from(1);
        Series {
            valuation: 1,
            coeffs,
            ramification: 1,
        }
    }

    /// `a` and `b` in the same units.
    fn aligned(a: &Series, b: &Series) -> (Series, Series) {
        let d = a.ramification / gcd(a.ramification, b.ramification) * b.ramification;
        (a.refined(d / a.ramification), b.refined(d / b.ramification))
    }

    pub(crate) fn add(&mut self, a: &Series, b: &Series) -> Result<Series, EvalError> {
        let (a, b) = Self::aligned(a, b);
        let v = a.valuation.min(b.valuation);
        let order = a.order_units().min(b.order_units());
        let mut coeffs = Vec::new();
        for u in v..order {
            let c = match (a.at(u), b.at(u)) {
                (Some(x), Some(y)) => self.eval(&(x.clone() + y.clone()))?,
                (Some(x), None) | (None, Some(x)) => x.clone(),
                (None, None) => Expr::from(0),
            };
            coeffs.push(c);
        }
        Ok(Series {
            valuation: v,
            coeffs,
            ramification: a.ramification,
        }
        .reduced())
    }

    pub(crate) fn mul(&mut self, a: &Series, b: &Series) -> Result<Series, EvalError> {
        let (a, b) = Self::aligned(a, b);
        let v = a.valuation + b.valuation;
        let order = (a.valuation + b.order_units()).min(b.valuation + a.order_units());
        let mut coeffs = Vec::new();
        for u in v..order {
            let mut terms = Vec::new();
            for (i, x) in a.coeffs.iter().enumerate() {
                let i = i64::try_from(i).expect("small");
                if let Some(y) = b.at(u - a.valuation - i) {
                    if *x != Expr::from(0) && *y != Expr::from(0) {
                        terms.push(x.clone() * y.clone());
                    }
                }
            }
            let sum: Expr = terms.into_iter().sum();
            coeffs.push(self.eval(&sum)?);
        }
        Ok(Series {
            valuation: v,
            coeffs,
            ramification: a.ramification,
        }
        .reduced())
    }

    fn scaled(&mut self, c: &Expr, a: &Series) -> Result<Series, EvalError> {
        let mut coeffs = Vec::with_capacity(a.coeffs.len());
        for x in &a.coeffs {
            coeffs.push(self.eval(&(c.clone() * x.clone()))?);
        }
        Ok(Series {
            coeffs,
            ..a.clone()
        })
    }

    /// `a^α` for `α` free of `t`; `None` when `a` has no known nonzero
    /// coefficient, or `α` is not rational and `a` has a nonzero valuation.
    pub(crate) fn power(&mut self, a: &Series, alpha: &Expr) -> Found {
        let a = self.normalized(a.clone())?;
        let Some(c) = a.coeffs.first().cloned() else {
            return Ok(None);
        };
        // The exponent of the leading term, `vα/d`, in units of `1/(d q)`.
        let (valuation, refine) = if a.valuation == 0 {
            (0, 1)
        } else {
            let Some(r) = alpha.to_rational() else {
                return Ok(None);
            };
            let (p, den) = (r.numer().to_i64(), r.denom().to_i64());
            let (Some(p), Some(den)) = (p, den) else {
                return Ok(None);
            };
            (a.valuation * p, den)
        };
        let n = a.coeffs.len();
        // h_j = a_j / c, q = (1 + h)^α.
        let inverse = self.eval(&c.pow(-1))?;
        let mut h = Vec::with_capacity(n);
        for x in &a.coeffs {
            h.push(self.eval(&(x.clone() * inverse.clone()))?);
        }
        let mut qs: Vec<Expr> = vec![Expr::from(1)];
        for k in 1..n {
            let kk = i64::try_from(k).expect("small");
            let mut terms = Vec::new();
            for j in 1..=k {
                if h[j] == Expr::from(0) {
                    continue;
                }
                let jj = i64::try_from(j).expect("small");
                let weight = (alpha.clone() + Expr::from(1)) * Expr::from(jj) - Expr::from(kk);
                terms.push(weight * h[j].clone() * qs[k - j].clone());
            }
            let sum: Expr = terms.into_iter().sum();
            qs.push(self.eval(&(sum / Expr::from(kk)))?);
        }
        let lead = self.eval(&c.pow(alpha.clone()))?;
        let body = Series {
            valuation: 0,
            coeffs: qs,
            ramification: a.ramification,
        }
        .refined(refine);
        let mut result = self.scaled(&lead, &body)?;
        result.valuation = valuation;
        Ok(Some(result.reduced()))
    }

    /// `g₀` and `h = g − g₀` for `g` of valuation 0 at least; `None` for a
    /// pole.
    fn split(&mut self, g: &Series) -> Result<Option<(Expr, Series)>, EvalError> {
        let g = self.normalized(g.clone())?;
        if g.valuation < 0 {
            return Ok(None);
        }
        let mut h = g.clone();
        let g0 = if g.valuation == 0 && !g.coeffs.is_empty() {
            let g0 = g.coeffs[0].clone();
            h.coeffs[0] = Expr::from(0);
            g0
        } else {
            Expr::from(0)
        };
        // h from valuation 0, so that h_j is the coefficient of s^j.
        if h.valuation > 0 {
            let lead = usize::try_from(h.valuation).expect("positive");
            let mut coeffs = vec![Expr::from(0); lead];
            coeffs.extend(h.coeffs);
            h = Series {
                valuation: 0,
                coeffs,
                ramification: h.ramification,
            };
        }
        Ok(Some((g0, h)))
    }

    /// The coefficients `j h_j` of `s h′(s)`.
    fn weighted(h: &Series) -> Vec<Expr> {
        h.coeffs
            .iter()
            .enumerate()
            .map(|(j, c)| Expr::from(i64::try_from(j).expect("small")) * c.clone())
            .collect()
    }

    /// `exp g`.
    pub(crate) fn exp(&mut self, g: &Series) -> Found {
        let Some((g0, h)) = self.split(g)? else {
            return Ok(None);
        };
        let w = Self::weighted(&h);
        let n = h.coeffs.len();
        let mut e: Vec<Expr> = vec![Expr::from(1); n.min(1)];
        for k in 1..n {
            let terms: Expr = (1..=k)
                .filter(|&j| h.coeffs[j] != Expr::from(0))
                .map(|j| w[j].clone() * e[k - j].clone())
                .sum();
            let kk = i64::try_from(k).expect("small");
            e.push(self.eval(&(terms / Expr::from(kk)))?);
        }
        let body = Series { coeffs: e, ..h };
        // e^{c + a log t} is e^c t^a for a rational a.
        let Some((c, a)) = self.logarithmic_part(&g0)? else {
            return Ok(None);
        };
        let lead = self.eval(&Expr::call(Symbol::EXP, [c]))?;
        let mut body = self.scaled(&lead, &body)?;
        if !a.is_zero() {
            let (Some(p), Some(den)) = (a.numer().to_i64(), a.denom().to_i64()) else {
                return Ok(None);
            };
            let m = den / gcd(den, body.ramification);
            body = body.refined(m);
            body.valuation += p * body.ramification / den;
        }
        Ok(Some(body.reduced()))
    }

    /// `g₀` as `c + a log t` for a rational `a`; `None` when `log t` is in
    /// it otherwise.
    fn logarithmic_part(&mut self, g0: &Expr) -> Result<Option<(Expr, Rational)>, EvalError> {
        let l = log_t();
        if !g0.contains(&l) {
            return Ok(Some((g0.clone(), Rational::ZERO)));
        }
        let coefficient =
            |k: i64| Expr::call(Symbol::COEFFICIENT, [g0.clone(), l.clone(), Expr::from(k)]);
        let a = self.eval(&coefficient(1))?;
        let c = self.eval(&coefficient(0))?;
        let Some(a) = a.to_rational() else {
            return Ok(None);
        };
        let rest = g0.clone() - c.clone() - Expr::from(a.clone()) * l;
        if !self.is_zero(&rest)? {
            return Ok(None);
        }
        Ok(Some((c, a)))
    }

    /// The sign of a real `c` for `t > 0`: by the assumptions, else by its
    /// ball when it is a number; `None` when neither tells.
    pub(crate) fn sign(&mut self, c: &Expr) -> Option<Ordering> {
        match sign_of(self.ev.context(), c) {
            Sign::Positive => return Some(Ordering::Greater),
            Sign::Negative => return Some(Ordering::Less),
            Sign::Zero => return Some(Ordering::Equal),
            _ => {}
        }
        if c.contains(&log_t()) {
            return None;
        }
        let value = approximation(self.ev.context(), c, 128).ok()?;
        match value.view() {
            View::Ball(b) => b.sign(),
            _ => None,
        }
    }

    /// `log g`, `log t` as [`log_t`].
    pub(crate) fn log(&mut self, g: &Series) -> Found {
        let g = self.normalized(g.clone())?;
        let Some(c) = g.coeffs.first().cloned() else {
            return Ok(None);
        };
        let inverse = self.eval(&c.pow(-1))?;
        let n = g.coeffs.len();
        let mut h = Vec::with_capacity(n);
        for x in &g.coeffs {
            h.push(self.eval(&(x.clone() * inverse.clone()))?);
        }
        h[0] = Expr::from(0);
        // L = log(1 + h): L_k = h_k − (1/k) Σ_{j=1}^{k−1} j L_j h_{k−j}.
        let mut l: Vec<Expr> = vec![Expr::from(0)];
        for k in 1..n {
            let kk = i64::try_from(k).expect("small");
            let terms: Expr = (1..k)
                .filter(|&j| h[k - j] != Expr::from(0) && l[j] != Expr::from(0))
                .map(|j| {
                    Expr::from(i64::try_from(j).expect("small")) * l[j].clone() * h[k - j].clone()
                })
                .sum();
            l.push(self.eval(&(h[k].clone() - terms / Expr::from(kk)))?);
        }
        let exponent = Expr::from(q(g.valuation, g.ramification));
        let lead = Expr::call(Symbol::LOG, [c]) + exponent * log_t();
        l[0] = self.eval(&lead)?;
        Ok(Some(
            Series {
                valuation: 0,
                coeffs: l,
                ramification: g.ramification,
            }
            .reduced(),
        ))
    }

    /// `sin g` (`sine`) or `cos g`.
    pub(crate) fn sin_cos(&mut self, g: &Series, sine: bool) -> Found {
        let Some((g0, h)) = self.split(g)? else {
            return Ok(None);
        };
        let w = Self::weighted(&h);
        let n = h.coeffs.len();
        let (mut s, mut c): (Vec<Expr>, Vec<Expr>) =
            (vec![Expr::from(0); n.min(1)], vec![Expr::from(1); n.min(1)]);
        for k in 1..n {
            let kk = Expr::from(i64::try_from(k).expect("small"));
            let (mut sk, mut ck) = (Vec::new(), Vec::new());
            for j in (1..=k).filter(|&j| h.coeffs[j] != Expr::from(0)) {
                sk.push(w[j].clone() * c[k - j].clone());
                ck.push(w[j].clone() * s[k - j].clone());
            }
            let (sk, ck): (Expr, Expr) = (sk.into_iter().sum(), ck.into_iter().sum());
            s.push(self.eval(&(sk / kk.clone()))?);
            c.push(self.eval(&(-ck / kk))?);
        }
        let (sin0, cos0) = (
            self.eval(&Expr::call(Symbol::SIN, [g0.clone()]))?,
            self.eval(&Expr::call(Symbol::COS, [g0]))?,
        );
        let (s, c) = (
            Series {
                coeffs: s,
                ..h.clone()
            },
            Series { coeffs: c, ..h },
        );
        let (a, b) = if sine {
            (self.scaled(&sin0, &c)?, self.scaled(&cos0, &s)?)
        } else {
            let minus_sin0 = self.eval(&(-sin0))?;
            (self.scaled(&cos0, &c)?, self.scaled(&minus_sin0, &s)?)
        };
        Ok(Some(self.add(&a, &b)?))
    }

    /// `d g/d t`.
    fn derivative(&mut self, g: &Series) -> Result<Series, EvalError> {
        let d = g.ramification;
        let mut coeffs = Vec::with_capacity(g.coeffs.len());
        for (k, c) in g.coeffs.iter().enumerate() {
            let u = g.valuation + i64::try_from(k).expect("small");
            coeffs.push(self.eval(&(Expr::from(q(u, d)) * c.clone()))?);
        }
        Ok(Series {
            valuation: g.valuation - d,
            coeffs,
            ramification: d,
        })
    }

    /// `∫ g dt` without a constant term; `None` for a term in `1/t`.
    fn integral(&mut self, g: &Series) -> Found {
        let d = g.ramification;
        let mut coeffs = Vec::with_capacity(g.coeffs.len());
        for (k, c) in g.coeffs.iter().enumerate() {
            let u = g.valuation + i64::try_from(k).expect("small");
            if u + d == 0 {
                if !self.is_zero(c)? {
                    return Ok(None);
                }
                coeffs.push(Expr::from(0));
                continue;
            }
            coeffs.push(self.eval(&(Expr::from(q(d, u + d)) * c.clone()))?);
        }
        Ok(Some(Series {
            valuation: g.valuation + d,
            coeffs,
            ramification: d,
        }))
    }

    /// `f(g)` as `f(g₀) + ∫ f′(g) g′` for `ArcTan`, `ArcSin`, `ArcCos` and
    /// the inverse hyperbolic functions;
    /// `ArcTan` of a series with a pole as `±π/2 − atan(1/g)`, the sign of
    /// `g`'s leading coefficient for `t > 0`.
    fn inverse(&mut self, f: Symbol, g: &Series) -> Found {
        let g = self.normalized(g.clone())?;
        if f == Symbol::ARC_TAN && g.valuation < 0 {
            let Some(c) = g.coeffs.first().cloned() else {
                return Ok(None);
            };
            let side = match self.sign(&c) {
                Some(Ordering::Greater) => 1,
                Some(Ordering::Less) => -1,
                _ => return Ok(None),
            };
            let Some(reciprocal) = self.power(&g, &Expr::from(-1))? else {
                return Ok(None);
            };
            let Some(arctangent) = self.inverse(f, &reciprocal)? else {
                return Ok(None);
            };
            let quarter = Expr::from(q(side, 2)) * Expr::from(Symbol::PI);
            let start = Self::constant(&quarter, &arctangent.order());
            let minus = self.scaled(&Expr::from(-1), &arctangent)?;
            return Ok(Some(self.add(&start, &minus)?));
        }
        let Some((g0, _)) = self.split(&g)? else {
            return Ok(None);
        };
        let g = self.normalized(g.clone())?;
        let square = self.mul(&g, &g)?;
        let one = Self::constant(&Expr::from(1), &g.order());
        let minus_square = self.scaled(&Expr::from(-1), &square)?;
        let minus_half = Expr::from(q(-1, 2));
        let derivative_of_f = match f {
            Symbol::ARC_TAN => {
                let w = self.add(&one, &square)?;
                self.power(&w, &Expr::from(-1))?
            }
            Symbol::ARC_TANH => {
                let w = self.add(&one, &minus_square)?;
                self.power(&w, &Expr::from(-1))?
            }
            Symbol::ARC_SINH => {
                let w = self.add(&one, &square)?;
                self.power(&w, &minus_half)?
            }
            Symbol::ARC_COSH => {
                // 1/(√(g − 1) √(g + 1)), the principal branch's.
                let minus_one = Self::constant(&Expr::from(-1), &g.order());
                let below = self.add(&g, &minus_one)?;
                let above = self.add(&g, &one)?;
                match (
                    self.power(&below, &minus_half)?,
                    self.power(&above, &minus_half)?,
                ) {
                    (Some(a), Some(b)) => Some(self.mul(&a, &b)?),
                    _ => None,
                }
            }
            _ => {
                let w = self.add(&one, &minus_square)?;
                let root = self.power(&w, &minus_half)?;
                match (f, root) {
                    (Symbol::ARC_COS, Some(r)) => Some(self.scaled(&Expr::from(-1), &r)?),
                    (_, r) => r,
                }
            }
        };
        let Some(derivative_of_f) = derivative_of_f else {
            return Ok(None);
        };
        let dg = self.derivative(&g)?;
        let integrand = self.mul(&derivative_of_f, &dg)?;
        let Some(body) = self.integral(&integrand)? else {
            return Ok(None);
        };
        let at = self.eval(&Expr::call(f, [g0]))?;
        let start = Self::constant(&at, &body.order());
        Ok(Some(self.add(&start, &body)?))
    }

    /// The series of `f` to the order `order`, maybe less (§2.5).
    pub(crate) fn expand(&mut self, f: &Expr, order: &Rational) -> Found {
        if !f.contains(&self.t) {
            return Ok(Some(Self::constant(f, order)));
        }
        if *f == self.t {
            return Ok(Some(Self::variable(order)));
        }
        let Some(head) = f.head().and_then(Expr::as_symbol) else {
            return Ok(None);
        };
        let args = f.args();
        match (head, args) {
            (Symbol::PLUS, _) => {
                let mut sum: Option<Series> = None;
                for a in args {
                    let Some(s) = self.expand(a, order)? else {
                        return Ok(None);
                    };
                    sum = Some(match sum {
                        Some(acc) => self.add(&acc, &s)?,
                        None => s,
                    });
                }
                Ok(sum)
            }
            (Symbol::TIMES, _) => {
                let mut product: Option<Series> = None;
                for a in args {
                    let Some(s) = self.expand(a, order)? else {
                        return Ok(None);
                    };
                    product = Some(match product {
                        Some(acc) => self.mul(&acc, &s)?,
                        None => s,
                    });
                }
                Ok(product)
            }
            (Symbol::POWER, [b, e]) if b.as_symbol() == Some(Symbol::E) => {
                let Some(s) = self.expand(e, order)? else {
                    return Ok(None);
                };
                self.exp(&s)
            }
            (Symbol::POWER, [b, e]) if !e.contains(&self.t) => {
                let Some(s) = self.expand(b, order)? else {
                    return Ok(None);
                };
                self.power(&s, e)
            }
            (Symbol::POWER, [b, e]) => {
                // b^e = exp(e log b).
                let exponent = e.clone() * Expr::call(Symbol::LOG, [b.clone()]);
                let Some(s) = self.expand(&exponent, order)? else {
                    return Ok(None);
                };
                self.exp(&s)
            }
            (
                Symbol::EXP
                | Symbol::LOG
                | Symbol::SIN
                | Symbol::COS
                | Symbol::TAN
                | Symbol::ARC_TAN
                | Symbol::ARC_SIN
                | Symbol::ARC_COS
                | Symbol::SINH
                | Symbol::COSH
                | Symbol::TANH
                | Symbol::ARC_SINH
                | Symbol::ARC_COSH
                | Symbol::ARC_TANH
                | Symbol::ABS,
                [g],
            ) => {
                let Some(s) = self.expand(g, order)? else {
                    return Ok(None);
                };
                match head {
                    Symbol::EXP => self.exp(&s),
                    Symbol::LOG => self.log(&s),
                    Symbol::SIN => self.sin_cos(&s, true),
                    Symbol::COS => self.sin_cos(&s, false),
                    Symbol::TAN => {
                        let (Some(sin), Some(cos)) =
                            (self.sin_cos(&s, true)?, self.sin_cos(&s, false)?)
                        else {
                            return Ok(None);
                        };
                        let Some(sec) = self.power(&cos, &Expr::from(-1))? else {
                            return Ok(None);
                        };
                        Ok(Some(self.mul(&sin, &sec)?))
                    }
                    Symbol::SINH | Symbol::COSH | Symbol::TANH => self.hyperbolic(head, &s),
                    Symbol::ABS => self.abs(&s),
                    _ => self.inverse(head, &s),
                }
            }
            _ => Ok(None),
        }
    }

    /// `sinh g`, `cosh g` or `tanh g` from `e^g` and `e^(−g)`.
    fn hyperbolic(&mut self, f: Symbol, g: &Series) -> Found {
        let minus = self.scaled(&Expr::from(-1), g)?;
        let (Some(plus), Some(minus)) = (self.exp(g)?, self.exp(&minus)?) else {
            return Ok(None);
        };
        let negated = self.scaled(&Expr::from(-1), &minus)?;
        let difference = self.add(&plus, &negated)?;
        let sum = self.add(&plus, &minus)?;
        let half = Expr::from(q(1, 2));
        Ok(match f {
            Symbol::SINH => Some(self.scaled(&half, &difference)?),
            Symbol::COSH => Some(self.scaled(&half, &sum)?),
            _ => match self.power(&sum, &Expr::from(-1))? {
                Some(inverse) => Some(self.mul(&difference, &inverse)?),
                None => None,
            },
        })
    }

    /// `|g|` for `t > 0` and a leading coefficient of known sign.
    fn abs(&mut self, g: &Series) -> Found {
        let g = self.normalized(g.clone())?;
        let Some(c) = g.coeffs.first().cloned() else {
            return Ok(None);
        };
        match self.sign(&c) {
            Some(Ordering::Greater) => Ok(Some(g)),
            Some(Ordering::Less) => Ok(Some(self.scaled(&Expr::from(-1), &g)?)),
            _ => Ok(None),
        }
    }
}

/// The series of `f` in `t` to the order `order`, asking the children for
/// more where a pole or a cancellation leaves less (§2.5): `None` when `f`
/// has none, or it is not reached.
pub(crate) fn expansion(ev: &mut Eval<'_>, f: &Expr, t: &Expr, order: &Rational) -> Found {
    let mut expander = Expander::new(ev, t.clone());
    for extra in [0, 2, 4, 8, 16] {
        let wanted = order + &Rational::from(Integer::from(extra));
        let Some(s) = expander.expand(f, &wanted)? else {
            return Ok(None);
        };
        let s = expander.normalized(s)?;
        if s.order() >= *order {
            return Ok(Some(s.truncated(order)));
        }
    }
    Ok(None)
}

mod builtin;
#[cfg(test)]
mod tests;

pub(crate) use builtin::{Point, fresh, series};
