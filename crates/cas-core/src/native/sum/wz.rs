//! `WZCertificate[f, s, n, k]`, Xcas's `wz_certificate` (CAS-580): the
//! Wilf–Zeilberger certificate of `Σ_k f(n, k) = s(n)`. A port of giac's
//! intg.cc `_wz_certificate`.
//!
//! - Arguments as giac's: `[F]`, `[f, s]` (`F = f/s`), `[F, n, k]` and
//!   `[f, s, n, k]`; `n` and `k` by default.
//! - giac takes `dF = F(n + 1, k) − F(n, k)`, its indefinite sum `G` in `k`
//!   (`G(k + 1) − G(k) = dF`, giac's `sum(dF, k)`), and answers
//!   `R = factor(G/F(n, k − 1))`, or 0 when the sum has no closed form.
//! - Ours keeps `R` a rational function with no factorials to cancel:
//!   `F(n + 1, k)/F(n, k) = r_n` and `F(n, k + 1)/F(n, k) = r_k` by
//!   Gosper's ratios, `dF = F (r_n − 1)`, Gosper's `G = y dF` for a
//!   rational `y` ([`multiplier`]), so `R = y (r_n − 1) r_k(k − 1)`,
//!   `Factor`ed. `r_n = 1` gives `dF = 0` and `R = 0`, as giac's
//!   `sum(0, k)`.
//! - 0, giac's answer when there is no certificate, also when `F` is no
//!   hypergeometric term in `n` or `k` here or Gosper's algorithm finds no
//!   `G`.

use super::gosper::{multiplier, ratio};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

pub(crate) fn wz_certificate(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let default = |s: &str| Expr::from(Symbol::new(s));
    let (f, n, k) = match e.args() {
        [f] => (f.clone(), default("n"), default("k")),
        [f, s] => (f.clone() / s.clone(), default("n"), default("k")),
        [f, n, k] => (f.clone(), n.clone(), k.clone()),
        [f, s, n, k] => (f.clone() / s.clone(), n.clone(), k.clone()),
        _ => return Ok(None),
    };
    if n.as_symbol().is_none() || k.as_symbol().is_none() {
        return Ok(None);
    }
    let f = ev.eval(&f)?;
    let result = match certificate(ev, e, &f, &n, &k)? {
        Some(r) => r,
        None => Expr::from(0),
    };
    ev.step("WZCertificate", e, &result);
    Ok(Some(result))
}

/// `R` of `F` in `n` and `k`, `None` when there is none here.
fn certificate(
    ev: &mut Eval<'_>,
    e: &Expr,
    f: &Expr,
    n: &Expr,
    k: &Expr,
) -> Result<Option<Expr>, EvalError> {
    let together = |ev: &mut Eval<'_>, x: Expr| ev.eval(&Expr::call(Symbol::TOGETHER, [x]));
    let Some(rn) = ratio(ev, f, n)? else {
        return Ok(None);
    };
    let rat = together(ev, rn - Expr::from(1))?;
    if rat == Expr::from(0) {
        return Ok(Some(rat));
    }
    let df = ev.eval(&(f.clone() * rat.clone()))?;
    let Some(y) = multiplier(ev, e, &df, k)? else {
        return Ok(None);
    };
    let Some(rk) = ratio(ev, f, k)? else {
        return Ok(None);
    };
    let back = rk.substitute(&[(k.clone(), k.clone() - Expr::from(1))]);
    let r = together(ev, y * rat * back)?;
    Ok(Some(ev.eval(&Expr::call(Symbol::FACTOR, [r]))?))
}
