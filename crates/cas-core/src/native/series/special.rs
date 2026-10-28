//! The special functions in the expander (CAS-DOC-31 §5, CAS-358), of a
//! series `g = g₀ + h`:
//!
//! - `Si`, `Erf` and `Erfc` everywhere, and `Ei` and `Ci` at `g₀ ≠ 0`, as
//!   `f(g₀) + ∫ f′(g) g′`, the inverse functions' path, `f′` being `sin g/g`,
//!   `±(2/√π) e^(−g²)`, `e^g/g` and `cos g/g`.
//! - `Ei` and `Ci` at `g₀ = 0`: `γ + log |g| + ∫ (e^g − 1)/g g′`, `|g|` by
//!   the sign of the leading coefficient, and `γ + log g + ∫ (cos g − 1)/g
//!   g′`, plus `iπ` for a negative one; `log` of a vanishing series keeps
//!   `log t`, so that `Ei(x)` at 0 is `γ + ln x + x + x²/4 + …`, where giac
//!   fails with a sign error.
//! - `PolyGamma[g, m]` and `Gamma` at a regular `g₀` by Taylor's
//!   coefficients `ψ^(m+k)(g₀)/k!`, composed with `h`, and `Γ(g₀ + h) =
//!   Γ(g₀) exp(Σ_{k≥1} ψ^(k−1)(g₀) h^k/k!)`; at a pole `g₀ = −n`, `n ≤`
//!   [`MAX_SHIFT`], by the recurrences `ψ^(m)(g) = ψ^(m)(g + n + 1) − (−1)^m
//!   m! Σ_{j≤n} (g + j)^(−(m+1))` and `Γ(g) = Γ(g + n + 1)/Π_{j≤n} (g + j)`.
//! - `LambertW` at `g₀ = 0` on branch 0 by `Σ_{n≥1} (−n)^(n−1)/n! g^n`, and
//!   at other `g₀` but the branch point by Newton's steps on series, `w ← w −
//!   (w e^w − g)/(e^w (1 + w))` from `W(g₀)`, each doubling the terms known.
//! - At a pole of `g`, its leading coefficient of known sign, `ν` the
//!   order of the pole: `Si(g) = ±π/2 + O(t^ν)` and `Ci(g) = 0 + O(t^ν)`
//!   (`iπ` below), the remainder oscillating within `C/|g|`; `Erf(g) = ±1`,
//!   `Erfc(g) = 0` or 2 and `Ei(g) = 0` toward `−∞` to every order, their
//!   remainders `e^(−g²)` and `e^g` beating every power of `t`; and
//!   `ψ^(m)(g)` toward `+∞` by its asymptotic series, `log g − 1/(2g) −
//!   Σ B_{2k}/(2k g^(2k))` and its derivatives, an O-bound at each order.
//!   `Ei` toward `+∞`, `Gamma` and `LambertW` have no series there.

use std::cmp::Ordering;

use cas_num::ball::bernoulli;
use cas_num::{Integer, Natural, Rational};

use super::{Expander, Found, Series};
use crate::context::EvalError;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The most steps of the recurrences from a pole of `Γ` or `ψ`.
const MAX_SHIFT: i64 = 64;

/// The most terms of `ψ`'s asymptotic series, as `Ball::gamma`'s table of
/// Bernoulli numbers.
const MAX_TERMS: u64 = 600;

fn factorial(n: u64) -> Expr {
    Expr::from(Integer::from(Natural::factorial(n)))
}

fn fraction(n: Integer, d: Integer) -> Expr {
    Expr::from(Rational::new(n, d))
}

impl Expander<'_, '_> {
    /// `f(args)` for the special functions, their arguments' series to
    /// `order`; `None` where they have none here.
    pub(super) fn special(&mut self, f: Symbol, args: &[Expr], order: &Rational) -> Found {
        let (g, extra) = match args {
            [g] => (g, None),
            [g, k] => (g, Some(k)),
            _ => return Ok(None),
        };
        let Some(s) = self.expand(g, order)? else {
            return Ok(None);
        };
        let s = self.normalized(s)?;
        match (f, extra) {
            (Symbol::POLY_GAMMA, None) => self.polygamma(&s, 0),
            (Symbol::POLY_GAMMA, Some(m)) => match m.as_i64().and_then(|m| u64::try_from(m).ok()) {
                Some(m) => self.polygamma(&s, m),
                None => Ok(None),
            },
            (Symbol::GAMMA, None) => self.gamma(&s),
            (Symbol::LAMBERT_W, None) => self.lambert(&s, false),
            (Symbol::LAMBERT_W, Some(k)) => match k.as_i64() {
                Some(0) => self.lambert(&s, false),
                Some(-1) => self.lambert(&s, true),
                _ => Ok(None),
            },
            (
                Symbol::EXP_INTEGRAL_EI
                | Symbol::SIN_INTEGRAL
                | Symbol::COS_INTEGRAL
                | Symbol::ERF
                | Symbol::ERFC,
                None,
            ) => self.integral_function(f, &s, order),
            _ => Ok(None),
        }
    }

    /// `c` with the order of `like`.
    fn constant_like(c: &Expr, like: &Series) -> Series {
        Self::constant(c, &like.order())
    }

    fn minus(&mut self, a: &Series, b: &Series) -> Result<Series, EvalError> {
        let negated = self.scaled(&Expr::from(-1), b)?;
        self.add(a, &negated)
    }

    /// `Σ c_k h^k` by Horner's scheme, `h` of positive valuation.
    fn compose(&mut self, coefficients: &[Expr], h: &Series) -> Result<Series, EvalError> {
        let order = h.order();
        let Some((last, rest)) = coefficients.split_last() else {
            return Ok(Self::constant(&Expr::from(0), &order));
        };
        let mut acc = Self::constant(last, &order);
        for c in rest.iter().rev() {
            acc = self.mul(&acc, h)?;
            let start = Self::constant(c, &order);
            acc = self.add(&acc, &start)?;
        }
        Ok(acc)
    }

    /// The sign of the leading coefficient of a normalized `s`.
    fn leading_sign(&mut self, s: &Series) -> Option<Ordering> {
        let c = s.coeffs.first()?.clone();
        self.sign(&c)
    }

    /// `Ei`, `Si`, `Ci`, `Erf` and `Erfc` of `s`.
    fn integral_function(&mut self, f: Symbol, s: &Series, order: &Rational) -> Found {
        if s.valuation < 0 {
            return self.at_infinity(f, s, order);
        }
        let Some((g0, _)) = self.split(s)? else {
            return Ok(None);
        };
        if matches!(f, Symbol::EXP_INTEGRAL_EI | Symbol::COS_INTEGRAL) && self.is_zero(&g0)? {
            return self.logarithmic_integral(f, s);
        }
        let derivative_of_f = match f {
            Symbol::ERF | Symbol::ERFC => {
                let square = self.mul(s, s)?;
                let minus_square = self.scaled(&Expr::from(-1), &square)?;
                let Some(gauss) = self.exp(&minus_square)? else {
                    return Ok(None);
                };
                let two = Expr::from(2) / Expr::call(Symbol::SQRT, [Expr::from(Symbol::PI)]);
                let two = if f == Symbol::ERF { two } else { -two };
                let two = self.eval(&two)?;
                self.scaled(&two, &gauss)?
            }
            _ => {
                let numerator = match f {
                    Symbol::SIN_INTEGRAL => self.sin_cos(s, true)?,
                    Symbol::COS_INTEGRAL => self.sin_cos(s, false)?,
                    _ => self.exp(s)?,
                };
                let (Some(numerator), Some(inverse)) = (numerator, self.power(s, &Expr::from(-1))?)
                else {
                    return Ok(None);
                };
                self.mul(&numerator, &inverse)?
            }
        };
        let dg = self.derivative(s)?;
        let integrand = self.mul(&derivative_of_f, &dg)?;
        let Some(body) = self.integral(&integrand)? else {
            return Ok(None);
        };
        let at = self.eval(&Expr::call(f, [g0]))?;
        let start = Self::constant_like(&at, &body);
        Ok(Some(self.add(&start, &body)?))
    }

    /// `Ei` or `Ci` of a vanishing `s`.
    fn logarithmic_integral(&mut self, f: Symbol, s: &Series) -> Found {
        let positive = match self.leading_sign(s) {
            Some(Ordering::Greater) => true,
            Some(Ordering::Less) => false,
            _ => return Ok(None),
        };
        let magnitude = if positive {
            s.clone()
        } else {
            self.scaled(&Expr::from(-1), s)?
        };
        let Some(log) = self.log(&magnitude)? else {
            return Ok(None);
        };
        let function = if f == Symbol::EXP_INTEGRAL_EI {
            self.exp(s)?
        } else {
            self.sin_cos(s, false)?
        };
        let (Some(function), Some(inverse)) = (function, self.power(s, &Expr::from(-1))?) else {
            return Ok(None);
        };
        let one = Self::constant_like(&Expr::from(1), &function);
        let less_one = self.minus(&function, &one)?;
        let quotient = self.mul(&less_one, &inverse)?;
        let dg = self.derivative(s)?;
        let integrand = self.mul(&quotient, &dg)?;
        let Some(body) = self.integral(&integrand)? else {
            return Ok(None);
        };
        let mut constant = Expr::from(Symbol::EULER_GAMMA);
        if f == Symbol::COS_INTEGRAL && !positive {
            constant += Expr::from(Symbol::I) * Expr::from(Symbol::PI);
        }
        let constant = self.eval(&constant)?;
        let start = Self::constant_like(&constant, &body);
        let sum = self.add(&start, &log)?;
        Ok(Some(self.add(&sum, &body)?))
    }

    /// `f(g)` for `g` with a pole, its leading coefficient of known sign.
    fn at_infinity(&mut self, f: Symbol, s: &Series, order: &Rational) -> Found {
        let up = match self.leading_sign(s) {
            Some(Ordering::Greater) => true,
            Some(Ordering::Less) => false,
            _ => return Ok(None),
        };
        // c + O(t^ν), ν the order of the pole.
        let short = |c: Expr| {
            let units = usize::try_from(-s.valuation).expect("a pole");
            let mut coeffs = vec![Expr::from(0); units];
            coeffs[0] = c;
            Series {
                valuation: 0,
                coeffs,
                ramification: s.ramification,
            }
        };
        let half_pi = Expr::from(Symbol::PI) / 2;
        let value = match (f, up) {
            (Symbol::SIN_INTEGRAL, true) => short(self.eval(&half_pi)?),
            (Symbol::SIN_INTEGRAL, false) => short(self.eval(&-half_pi)?),
            (Symbol::COS_INTEGRAL, true) => short(Expr::from(0)),
            (Symbol::COS_INTEGRAL, false) => {
                short(self.eval(&(Expr::from(Symbol::I) * Expr::from(Symbol::PI)))?)
            }
            (Symbol::ERF, _) => Self::constant(&Expr::from(if up { 1 } else { -1 }), order),
            (Symbol::ERFC, _) => Self::constant(&Expr::from(if up { 0 } else { 2 }), order),
            (Symbol::EXP_INTEGRAL_EI, false) => Self::constant(&Expr::from(0), order),
            _ => return Ok(None),
        };
        Ok(Some(value))
    }

    /// Whether `g₀` is a pole of `Γ` and `ψ`, an integer `−n` with `0 ≤ n ≤`
    /// [`MAX_SHIFT`]: `n`.
    fn pole(g0: &Expr) -> Option<i64> {
        let n = -g0.as_i64()?;
        (0..=MAX_SHIFT).contains(&n).then_some(n)
    }

    /// `s + k`.
    fn shifted(&mut self, s: &Series, k: i64) -> Result<Series, EvalError> {
        let c = Self::constant_like(&Expr::from(k), s);
        self.add(s, &c)
    }

    /// `ψ^(m)` of a normalized `s`.
    fn polygamma(&mut self, s: &Series, m: u64) -> Found {
        if s.valuation < 0 {
            if self.leading_sign(s) != Some(Ordering::Greater) {
                return Ok(None);
            }
            return self.polygamma_far(s, m);
        }
        let Some((g0, h)) = self.split(s)? else {
            return Ok(None);
        };
        if let Some(n) = Self::pole(&g0) {
            // ψ^(m)(g + n + 1) − (−1)^m m! Σ_{j≤n} (g + j)^(−(m+1)).
            let up = self.shifted(s, n + 1)?;
            let Some(regular) = self.polygamma(&up, m)? else {
                return Ok(None);
            };
            let power = Expr::from(-i64::try_from(m + 1).expect("small"));
            let mut sum: Option<Series> = None;
            for j in 0..=n {
                let gj = self.shifted(s, j)?;
                let Some(p) = self.power(&gj, &power)? else {
                    return Ok(None);
                };
                sum = Some(match sum {
                    Some(acc) => self.add(&acc, &p)?,
                    None => p,
                });
            }
            let sum = sum.expect("one term");
            let scale = if m % 2 == 0 {
                -factorial(m)
            } else {
                factorial(m)
            };
            let sum = self.scaled(&scale, &sum)?;
            return Ok(Some(self.add(&regular, &sum)?));
        }
        // Taylor's coefficients ψ^(m+k)(g₀)/k!.
        let terms = h.coeffs.len().max(1);
        let mut coefficients = Vec::with_capacity(terms);
        for k in 0..terms {
            let k = u64::try_from(k).expect("small");
            let order = Expr::from(i64::try_from(m + k).expect("small"));
            let value = Expr::call(Symbol::POLY_GAMMA, [g0.clone(), order]) / factorial(k);
            coefficients.push(self.eval(&value)?);
        }
        Ok(Some(self.compose(&coefficients, &h)?))
    }

    /// `ψ^(m)(g)` for `g → +∞` by its asymptotic series in `w = 1/g`.
    fn polygamma_far(&mut self, s: &Series, m: u64) -> Found {
        let Some(w) = self.power(s, &Expr::from(-1))? else {
            return Ok(None);
        };
        let w = self.normalized(w)?;
        if w.valuation <= 0 {
            return Ok(None);
        }
        // Powers of w up to the order: w^j has exponent j ν.
        let units = w.order_units();
        let most = u64::try_from(units / w.valuation + 1).unwrap_or(0);
        let k_terms = (most / 2 + 1).min(MAX_TERMS);
        let mut coefficients = vec![Expr::from(0); usize::try_from(most + m + 2).expect("small")];
        let set = |coefficients: &mut Vec<Expr>, j: u64, c: Expr| {
            if let Some(slot) = usize::try_from(j)
                .ok()
                .and_then(|j| coefficients.get_mut(j))
            {
                *slot = c;
            }
        };
        if m == 0 {
            // −w/2 − Σ B_{2k}/(2k) w^(2k), log g apart.
            set(
                &mut coefficients,
                1,
                fraction(Integer::from(-1), Integer::from(2)),
            );
            for k in 1..=k_terms {
                let Some(b) = bernoulli(2 * k) else {
                    break;
                };
                let c = -(&b / &Rational::from(Integer::from(2 * k)));
                set(&mut coefficients, 2 * k, Expr::from(c));
            }
        } else {
            // (−1)^(m+1) [(m − 1)! w^m + m!/2 w^(m+1) + Σ B_{2k} (2k + m − 1)!/(2k)! w^(2k+m)].
            let sign = if m % 2 == 0 { -1 } else { 1 };
            set(&mut coefficients, m, Expr::from(sign) * factorial(m - 1));
            set(
                &mut coefficients,
                m + 1,
                Expr::from(sign) * factorial(m) / 2,
            );
            for k in 1..=k_terms {
                let Some(b) = bernoulli(2 * k) else {
                    break;
                };
                let c = Expr::from(b) * factorial(2 * k + m - 1) / factorial(2 * k);
                set(&mut coefficients, 2 * k + m, Expr::from(sign) * c);
            }
        }
        let mut evaluated = Vec::with_capacity(coefficients.len());
        for c in coefficients {
            evaluated.push(self.eval(&c)?);
        }
        let body = self.compose(&evaluated, &w)?;
        if m != 0 {
            return Ok(Some(body));
        }
        let Some(log) = self.log(s)? else {
            return Ok(None);
        };
        Ok(Some(self.add(&log, &body)?))
    }

    /// `Γ` of a normalized `s`.
    fn gamma(&mut self, s: &Series) -> Found {
        if s.valuation < 0 {
            return Ok(None);
        }
        let Some((g0, h)) = self.split(s)? else {
            return Ok(None);
        };
        if let Some(n) = Self::pole(&g0) {
            // Γ(g + n + 1)/Π_{j≤n} (g + j).
            let up = self.shifted(s, n + 1)?;
            let Some(regular) = self.gamma(&up)? else {
                return Ok(None);
            };
            let mut product = self.shifted(s, 0)?;
            for j in 1..=n {
                let gj = self.shifted(s, j)?;
                product = self.mul(&product, &gj)?;
            }
            let Some(inverse) = self.power(&product, &Expr::from(-1))? else {
                return Ok(None);
            };
            return Ok(Some(self.mul(&regular, &inverse)?));
        }
        // Γ(g₀) exp(Σ_{k≥1} ψ^(k−1)(g₀) h^k/k!).
        let terms = h.coeffs.len().max(1);
        let mut coefficients = vec![Expr::from(0)];
        for k in 1..terms {
            let k = u64::try_from(k).expect("small");
            let order = Expr::from(i64::try_from(k - 1).expect("small"));
            let value = Expr::call(Symbol::POLY_GAMMA, [g0.clone(), order]) / factorial(k);
            coefficients.push(self.eval(&value)?);
        }
        let exponent = self.compose(&coefficients, &h)?;
        let Some(e) = self.exp(&exponent)? else {
            return Ok(None);
        };
        let lead = self.eval(&Expr::call(Symbol::GAMMA, [g0]))?;
        Ok(Some(self.scaled(&lead, &e)?))
    }

    /// `LambertW` of a normalized `s`, on branch −1 for `lower`.
    fn lambert(&mut self, s: &Series, lower: bool) -> Found {
        if s.valuation < 0 {
            return Ok(None);
        }
        let Some((g0, _)) = self.split(s)? else {
            return Ok(None);
        };
        if self.is_zero(&g0)? {
            if lower {
                return Ok(None);
            }
            // Σ_{n≥1} (−n)^(n−1)/n! g^n.
            let terms = usize::try_from(s.order_units().max(1)).unwrap_or(1) + 1;
            let mut coefficients = vec![Expr::from(0)];
            for n in 1..terms {
                let n = u64::try_from(n).expect("small");
                let base = Integer::from(-i64::try_from(n).expect("small"));
                let c = fraction(base.pow(n - 1), Integer::from(Natural::factorial(n)));
                coefficients.push(c);
            }
            return Ok(Some(self.compose(&coefficients, s)?));
        }
        let w0 = if lower {
            Expr::call(Symbol::LAMBERT_W, [g0, Expr::from(-1)])
        } else {
            Expr::call(Symbol::LAMBERT_W, [g0])
        };
        let w0 = self.eval(&w0)?;
        // Not at the branch point, where 1 + W = 0.
        let next = self.eval(&(Expr::from(1) + w0.clone()))?;
        if self.is_zero(&next)? {
            return Ok(None);
        }
        // Newton's steps, each doubling the terms known.
        let units = s.order_units().max(1);
        let steps = u64::BITS - u64::try_from(units).unwrap_or(1).leading_zeros() + 1;
        let mut w = Self::constant_like(&w0, s);
        for _ in 0..steps {
            let Some(e) = self.exp(&w)? else {
                return Ok(None);
            };
            let product = self.mul(&w, &e)?;
            let residual = self.minus(&product, s)?;
            let one = Self::constant_like(&Expr::from(1), &w);
            let one_plus = self.add(&one, &w)?;
            let slope = self.mul(&e, &one_plus)?;
            let Some(inverse) = self.power(&slope, &Expr::from(-1))? else {
                return Ok(None);
            };
            let step = self.mul(&residual, &inverse)?;
            w = self.minus(&w, &step)?;
        }
        Ok(Some(w))
    }
}

#[cfg(test)]
mod tests;
