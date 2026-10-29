//! Rational functions (CAS-DOC-41 §4, CAS-441): `f = P/Q` over ℚ,
//! `deg Q > deg P`, `Q` without a real root (Sturm's count), by residues
//! of `f(z) e^(−i s z)`: closing below for `s > 0`,
//! `F₊ = −2πi Σ_(Im z < 0) Res`, and above for `s < 0`,
//! `F₋ = 2πi Σ_(Im z > 0) Res`. A simple pole's residue is
//! `P(z₀) e^(−i s z₀)/Q′(z₀)`, another's `Residue`'s; the poles are
//! `CSolve`'s. The answer is `F₊(|s|)` for an even `f`,
//! `sign(s) F₊(|s|)` for an odd one, and `θ(s) F₊ + θ(−s) F₋` otherwise,
//! each part `ComplexExpand`ed, real for a real `s`: giac's
//! `pi*exp(-abs(s))` for `1/(1 + x²)`. With `deg Q = deg P + 1` the
//! integral converges as a principal value, and `θ(0) = 1/2` gives the
//! mean of the sides at `s = 0`.

use std::cmp::Ordering;

use cas_num::algebra::Rationals;
use cas_poly::DenseRing;

use super::super::integrate::rational_function;
use super::super::rounding::real_sign;
use super::{Found, i, pi};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The transform of `f`, a rational function of `x` as above.
pub(super) fn rational(ev: &mut Eval<'_>, f: &Expr, x: &Expr, s: &Expr) -> Found {
    let e = Expr::call(Symbol::FOURIER_TRANSFORM, [f.clone(), x.clone(), s.clone()]);
    let Some((num, den)) = rational_function(ev, &e, f, x)? else {
        return Ok(None);
    };
    let qx = DenseRing::new(Rationals);
    let (Some(p_degree), Some(q_degree)) = (num.degree(), den.degree()) else {
        return Ok(None);
    };
    if q_degree <= p_degree || qx.count_real_roots(&den, None, None) > 0 {
        return Ok(None);
    }
    let simple = qx.gcd(&den, &qx.derivative(&den)).degree() == Some(0);
    let together = ev.eval(&Expr::call(Symbol::TOGETHER, [f.clone()]))?;
    let p = ev.eval(&Expr::call(Symbol::NUMERATOR, [together.clone()]))?;
    let q = ev.eval(&Expr::call(Symbol::DENOMINATOR, [together.clone()]))?;
    let roots = ev.eval(&Expr::call(
        Symbol::CSOLVE,
        [
            Expr::call(Symbol::EQUATION, [q.clone(), Expr::from(0)]),
            x.clone(),
        ],
    ))?;
    if roots.head() != Some(&Expr::from(Symbol::LIST)) || roots.args().is_empty() {
        return Ok(None);
    }
    let kernel = Expr::from(Symbol::E).pow(-(i() * s.clone() * x.clone()));
    let dq = ev.eval(&Expr::call(Symbol::D, [q, x.clone()]))?;
    let (mut below, mut above) = (Vec::new(), Vec::new());
    for z0 in roots.args() {
        let im = ev.eval(&Expr::call(Symbol::IM, [z0.clone()]))?;
        let side = real_sign(ev.context(), &im);
        let residue = if simple {
            let at = [(x.clone(), z0.clone())];
            (p.clone() * kernel.clone() / dq.clone()).substitute(&at)
        } else {
            let residue = ev.eval(&Expr::call(
                Symbol::RESIDUE,
                [
                    together.clone() * kernel.clone(),
                    Expr::call(Symbol::LIST, [x.clone(), z0.clone()]),
                ],
            ))?;
            if residue.contains(&Expr::from(Symbol::RESIDUE)) {
                return Ok(None);
            }
            residue
        };
        match side {
            Some(Ordering::Less) => below.push(residue),
            Some(Ordering::Greater) => above.push(residue),
            _ => return Ok(None),
        }
    }
    let two_pi_i = Expr::from(2) * pi() * i();
    let plus = -two_pi_i.clone() * Expr::call(Symbol::PLUS, below);
    let minus = two_pi_i * Expr::call(Symbol::PLUS, above);
    let real = |ev: &mut Eval<'_>, g: Expr| ev.eval(&Expr::call(Symbol::COMPLEX_EXPAND, [g]));
    // The symmetries of f.
    let mirrored = f.substitute(&[(x.clone(), -x.clone())]);
    let symmetric = |ev: &mut Eval<'_>, sum: Expr| -> Result<bool, crate::context::EvalError> {
        let gap = ev.eval(&Expr::call(Symbol::TOGETHER, [sum]))?;
        Ok(gap == Expr::from(0))
    };
    let absolute = Expr::call(Symbol::ABS, [s.clone()]);
    if symmetric(ev, mirrored.clone() - f.clone())? {
        let plus = real(ev, plus)?;
        return Ok(Some(plus.substitute(&[(s.clone(), absolute)])));
    }
    if symmetric(ev, mirrored + f.clone())? {
        let plus = real(ev, plus)?;
        return Ok(Some(
            Expr::call(Symbol::SIGN, [s.clone()]) * plus.substitute(&[(s.clone(), absolute)]),
        ));
    }
    let (plus, minus) = (real(ev, plus)?, real(ev, minus)?);
    let step = |u: Expr| Expr::call(Symbol::HEAVISIDE_THETA, [u]);
    Ok(Some(step(s.clone()) * plus + step(-s.clone()) * minus))
}
