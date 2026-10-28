//! The Bessel and Airy functions in the expander (CAS-DOC-34 §4, CAS-383),
//! of a series `g = g₀ + h`:
//!
//! - **At a regular point**, by their differential equations: `y(g₀ + h) =
//!   Σ c_k h^k` from `c_0 = y(g₀)` and `c_1 = y′(g₀)`, each next `c_k` by
//!   the recurrence the equation gives, then composed with `h`. Each `c_k`
//!   is `α_k c_0 + β_k c_1`, the recurrence run on `α` and `β`, which are
//!   put together. Bessel's `x² y″ + x y′ + (±x² − n²) y = 0` (`+` for `J`
//!   and `Y`, `−` for `I` and `K`) gives, for `g₀ ≠ 0`, `c_(k+2) = −(g₀
//!   (k + 1)(2k + 1) c_(k+1) + (k² ± g₀² − n²) c_k ± 2g₀ c_(k−1) ±
//!   c_(k−2)) / (g₀² (k + 1)(k + 2))`; Airy's `y″ = x y` gives `c_(k+2) =
//!   (g₀ c_k + c_(k−1))/((k + 1)(k + 2))`, at `g₀ = 0` too, and `Ai′` and
//!   `Bi′` are the series of `Ai` and `Bi` one term longer, differentiated.
//! - **At `g₀ = 0`**: `J_n` and `I_n` by their power series `Σ_k (∓g²/4)^k
//!   (g/2)^n / (k! (n + k)!)`; `Y_n` and `K_n` by DLMF 10.8.1 and 10.31.1,
//!   for a positive leading coefficient, `ln(g/2)` keeping `log t` as `log`
//!   of a vanishing series does (CAS-DOC-10 §4.2), and `ψ(k + 1) = H_k − γ`.
//! - **At a pole of `g`**, `ν` its order, its leading coefficient's sign
//!   known: `J_n(g)` and, toward `+∞`, `Y_n(g)` are `0 + O(t^(ν/2))`, as
//!   `|J_n(x)|` and `|Y_n(x)|` are at most `C x^(−1/2)`; `Ai(g)` and
//!   `Bi(g)` toward `−∞` are `0 + O(t^(ν/4))`; `K_n(g)`, `Ai(g)` and
//!   `Ai′(g)` toward `+∞` are 0 to every order, their remainders `e^(−g)`
//!   and `e^(−(2/3) g^(3/2))` beating every power of `t`. The others have no
//!   series there: `I_n`, `Bi` and `Bi′` grow toward `+∞`, `Ai′` and `Bi′`
//!   toward `−∞`, and `Y_n` and `K_n` are complex below 0.

use std::cmp::Ordering;

use cas_num::{Integer, Natural, Rational};

use super::{Expander, Found, Series};
use crate::context::EvalError;
use crate::expr::Expr;
use crate::symbol::Symbol;

fn factorial(n: u64) -> Expr {
    Expr::from(Integer::from(Natural::factorial(n)))
}

fn small(k: u64) -> Expr {
    Expr::from(Integer::from(k))
}

/// The kinds of Bessel functions.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    J,
    Y,
    I,
    K,
}

impl Kind {
    fn of(head: Symbol) -> Option<Kind> {
        Some(match head {
            Symbol::BESSEL_J => Kind::J,
            Symbol::BESSEL_Y => Kind::Y,
            Symbol::BESSEL_I => Kind::I,
            Symbol::BESSEL_K => Kind::K,
            _ => return None,
        })
    }

    fn head(self) -> Symbol {
        match self {
            Kind::J => Symbol::BESSEL_J,
            Kind::Y => Symbol::BESSEL_Y,
            Kind::I => Symbol::BESSEL_I,
            Kind::K => Symbol::BESSEL_K,
        }
    }

    /// `+1` for `J` and `Y`, whose equation has `+x²`; `−1` for `I` and `K`.
    fn sign(self) -> i64 {
        if matches!(self, Kind::J | Kind::Y) {
            1
        } else {
            -1
        }
    }

    /// `y′(x)` from the neighbouring orders.
    fn derivative(self, n: &Expr, x: &Expr) -> Expr {
        let at = |m: Expr| Expr::call(self.head(), [m, x.clone()]);
        let (below, above) = (at(n.clone() - 1), at(n.clone() + 1));
        match self {
            Kind::J | Kind::Y => (below - above) / 2,
            Kind::I => (below + above) / 2,
            Kind::K => -(below + above) / 2,
        }
    }
}

impl Expander<'_, '_> {
    /// `BesselJ[n, g]` and its kin, and Airy's functions; `None` where they
    /// have no series here.
    pub(super) fn bessel(&mut self, f: Symbol, args: &[Expr], order: &Rational) -> Found {
        let (n, g) = match args {
            [n, g] => (Some(n), g),
            [g] => (None, g),
            _ => return Ok(None),
        };
        let Some(s) = self.expand(g, order)? else {
            return Ok(None);
        };
        let s = self.normalized(s)?;
        match (Kind::of(f), n) {
            (Some(kind), Some(n)) => self.bessel_of(kind, n, &s, order),
            (None, None) => self.airy_of(f, &s, order),
            _ => Ok(None),
        }
    }

    /// The coefficients `c_0 … c_(terms−1)` of `y(g₀ + h)` from `c_0`,
    /// `c_1` and the recurrence `next(k, α_(k+1), α_k, α_(k−1), α_(k−2))`,
    /// run on `α` and `β` and put together.
    fn by_recurrence(
        &mut self,
        c0: &Expr,
        c1: &Expr,
        terms: usize,
        next: &dyn Fn(u64, [&Expr; 4]) -> Expr,
    ) -> Result<Vec<Expr>, EvalError> {
        let zero = Expr::from(0);
        let mut alpha = vec![Expr::from(1), Expr::from(0)];
        let mut beta = vec![Expr::from(0), Expr::from(1)];
        while alpha.len() < terms {
            let k = alpha.len() - 2;
            let pick = |v: &Vec<Expr>| -> [Expr; 4] {
                [
                    v[k + 1].clone(),
                    v[k].clone(),
                    if k >= 1 {
                        v[k - 1].clone()
                    } else {
                        zero.clone()
                    },
                    if k >= 2 {
                        v[k - 2].clone()
                    } else {
                        zero.clone()
                    },
                ]
            };
            let k64 = u64::try_from(k).expect("small");
            let [a1, a0, a_1, a_2] = pick(&alpha);
            let [b1, b0, b_1, b_2] = pick(&beta);
            let a = next(k64, [&a1, &a0, &a_1, &a_2]);
            let b = next(k64, [&b1, &b0, &b_1, &b_2]);
            alpha.push(self.eval(&Expr::call(Symbol::TOGETHER, [a]))?);
            beta.push(self.eval(&Expr::call(Symbol::TOGETHER, [b]))?);
        }
        let mut c = Vec::with_capacity(terms);
        for (a, b) in alpha.into_iter().zip(beta).take(terms) {
            c.push(self.eval(&(a * c0.clone() + b * c1.clone()))?);
        }
        Ok(c)
    }

    /// A Bessel function of a normalized `s`.
    fn bessel_of(&mut self, kind: Kind, n: &Expr, s: &Series, order: &Rational) -> Found {
        if s.valuation < 0 {
            return self.bessel_far(kind, s, order);
        }
        let Some((g0, h)) = self.split(s)? else {
            return Ok(None);
        };
        if self.is_zero(&g0)? {
            return self.bessel_at_zero(kind, n, s);
        }
        // Bessel's equation at g₀ ≠ 0.
        let c0 = self.eval(&Expr::call(kind.head(), [n.clone(), g0.clone()]))?;
        let c1 = self.eval(&kind.derivative(n, &g0))?;
        let terms = h.coeffs.len().max(2);
        let sign = Expr::from(kind.sign());
        let square = g0.clone() * g0.clone();
        let n_square = n.clone() * n.clone();
        let next = |k: u64, [c1, c0, c_1, c_2]: [&Expr; 4]| -> Expr {
            let (k1, k2) = (small(k + 1), small(k + 2));
            let middle = small(k * k) + sign.clone() * square.clone() - n_square.clone();
            -(g0.clone() * k1.clone() * small(2 * k + 1) * c1.clone()
                + middle * c0.clone()
                + sign.clone() * 2 * g0.clone() * c_1.clone()
                + sign.clone() * c_2.clone())
                / (square.clone() * k1 * k2)
        };
        let c = self.by_recurrence(&c0, &c1, terms, &next)?;
        Ok(Some(self.compose(&c, &h)?))
    }

    /// `J_n`, `I_n`, `Y_n` or `K_n` of a vanishing `s`, for an integer `n ≥
    /// 0`.
    fn bessel_at_zero(&mut self, kind: Kind, n: &Expr, s: &Series) -> Found {
        let Some(n) = n.as_i64().and_then(|n| u64::try_from(n).ok()) else {
            return Ok(None);
        };
        let half = self.scaled(&(Expr::from(1) / 2), s)?;
        let square = self.mul(&half, &half)?;
        // u = ∓(g/2)²: − for J and Y, + for I and K.
        let u = if kind.sign() > 0 {
            self.scaled(&Expr::from(-1), &square)?
        } else {
            square.clone()
        };
        let terms = usize::try_from(s.order_units() / u.valuation.max(1) + 2).unwrap_or(2);
        let n_i64 = i64::try_from(n).expect("small");
        let Some(lead) = self.power(&half, &Expr::from(n_i64))? else {
            return Ok(None);
        };
        // Σ_k u^k / (k! (n + k)!), times (g/2)^n.
        let mut plain = Vec::with_capacity(terms);
        for k in 0..terms {
            let k = u64::try_from(k).expect("small");
            plain.push(self.eval(&(Expr::from(1) / (factorial(k) * factorial(n + k))))?);
        }
        let body = self.compose(&plain, &u)?;
        let regular = self.mul(&lead, &body)?;
        if matches!(kind, Kind::J | Kind::I) {
            return Ok(Some(regular));
        }
        if self.leading_sign(s) != Some(Ordering::Greater) {
            return Ok(None);
        }
        // The finite sum Σ_{k<n} (n − k − 1)!/k! (±(g/2)²)^k (g/2)^(−n):
        // + (g/2)² for Y, − for K.
        let finite_u = if kind == Kind::Y {
            square.clone()
        } else {
            self.scaled(&Expr::from(-1), &square)?
        };
        let mut finite = Vec::with_capacity(usize::try_from(n).unwrap_or(0));
        for k in 0..n {
            finite.push(self.eval(&(factorial(n - k - 1) / factorial(k)))?);
        }
        let Some(inverse) = self.power(&half, &Expr::from(-n_i64))? else {
            return Ok(None);
        };
        let finite = if finite.is_empty() {
            Self::constant_like(&Expr::from(0), s)
        } else {
            let sum = self.compose(&finite, &finite_u)?;
            self.mul(&inverse, &sum)?
        };
        // Σ_k (ψ(k + 1) + ψ(n + k + 1)) u^k / (k! (n + k)!), times (g/2)^n.
        let mut digamma = Vec::with_capacity(terms);
        for k in 0..terms {
            let k = u64::try_from(k).expect("small");
            let psi = |m: u64| Expr::call(Symbol::POLY_GAMMA, [small(m)]);
            let c = (psi(k + 1) + psi(n + k + 1)) / (factorial(k) * factorial(n + k));
            digamma.push(self.eval(&c)?);
        }
        let digamma_body = self.compose(&digamma, &u)?;
        let digamma = self.mul(&lead, &digamma_body)?;
        let Some(log) = self.log(&half)? else {
            return Ok(None);
        };
        let logarithmic = self.mul(&log, &regular)?;
        let pi = Expr::from(Symbol::PI);
        let (a, b, c) = if kind == Kind::Y {
            // −(1/π) finite + (2/π) ln(g/2) J_n − (1/π) digamma.
            let inverse_pi = self.eval(&(Expr::from(1) / pi.clone()))?;
            (
                self.eval(&-inverse_pi.clone())?,
                self.eval(&(Expr::from(2) * inverse_pi.clone()))?,
                self.eval(&-inverse_pi)?,
            )
        } else {
            // (1/2) finite + (−1)^(n+1) ln(g/2) I_n + (−1)^n (1/2) digamma.
            let odd = n % 2 == 1;
            (
                Expr::from(1) / 2,
                Expr::from(if odd { 1 } else { -1 }),
                Expr::from(if odd { -1 } else { 1 }) / 2,
            )
        };
        let first = self.scaled(&a, &finite)?;
        let second = self.scaled(&b, &logarithmic)?;
        let third = self.scaled(&c, &digamma)?;
        let sum = self.add(&first, &second)?;
        Ok(Some(self.add(&sum, &third)?))
    }

    /// `0 + O(t^(ν/d))` for a pole of order `ν`.
    fn vanishing(s: &Series, d: i64) -> Series {
        let units = usize::try_from(-s.valuation).expect("a pole");
        Series {
            valuation: 0,
            coeffs: vec![Expr::from(0); units],
            ramification: s.ramification * d,
        }
    }

    /// A Bessel function at a pole of `s`.
    fn bessel_far(&mut self, kind: Kind, s: &Series, order: &Rational) -> Found {
        let up = match self.leading_sign(s) {
            Some(Ordering::Greater) => true,
            Some(Ordering::Less) => false,
            _ => return Ok(None),
        };
        Ok(match (kind, up) {
            (Kind::J, _) | (Kind::Y, true) => Some(Self::vanishing(s, 2)),
            (Kind::K, true) => Some(Self::constant(&Expr::from(0), order)),
            _ => None,
        })
    }

    /// `AiryAi`, `AiryBi`, `AiryAiPrime` or `AiryBiPrime` of a normalized
    /// `s`.
    fn airy_of(&mut self, f: Symbol, s: &Series, order: &Rational) -> Found {
        let (base, prime) = match f {
            Symbol::AIRY_AI => (Symbol::AIRY_AI, false),
            Symbol::AIRY_BI => (Symbol::AIRY_BI, false),
            Symbol::AIRY_AI_PRIME => (Symbol::AIRY_AI, true),
            Symbol::AIRY_BI_PRIME => (Symbol::AIRY_BI, true),
            _ => return Ok(None),
        };
        if s.valuation < 0 {
            let up = match self.leading_sign(s) {
                Some(Ordering::Greater) => true,
                Some(Ordering::Less) => false,
                _ => return Ok(None),
            };
            return Ok(match (base, prime, up) {
                (Symbol::AIRY_AI, _, true) => Some(Self::constant(&Expr::from(0), order)),
                (_, false, false) => Some(Self::vanishing(s, 4)),
                _ => None,
            });
        }
        let Some((g0, h)) = self.split(s)? else {
            return Ok(None);
        };
        let derivative = if base == Symbol::AIRY_AI {
            Symbol::AIRY_AI_PRIME
        } else {
            Symbol::AIRY_BI_PRIME
        };
        let c0 = self.eval(&Expr::call(base, [g0.clone()]))?;
        let c1 = self.eval(&Expr::call(derivative, [g0.clone()]))?;
        let terms = h.coeffs.len().max(2) + usize::from(prime);
        let next = |k: u64, [_, c0, c_1, _]: [&Expr; 4]| -> Expr {
            (g0.clone() * c0.clone() + c_1.clone()) / (small(k + 1) * small(k + 2))
        };
        let mut c = self.by_recurrence(&c0, &c1, terms, &next)?;
        if prime {
            // y′(g₀ + h) = Σ (k + 1) c_(k+1) h^k.
            let mut differentiated = Vec::with_capacity(c.len() - 1);
            for (k, ck) in c.iter().enumerate().skip(1) {
                let k = u64::try_from(k).expect("small");
                differentiated.push(self.eval(&(small(k) * ck.clone()))?);
            }
            c = differentiated;
        }
        Ok(Some(self.compose(&c, &h)?))
    }
}

#[cfg(test)]
mod tests;
