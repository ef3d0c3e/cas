//! `SumRiemann[f, {n, k}]`, Xcas's `sum_riemann` (CAS-576): the limit as
//! `n → ∞` of `Σ_{k=1}^{n} f` read as a Riemann sum. A port of giac's
//! permu.cc `_sum_riemann` and `_split`.
//!
//! - `k = n x`: `a = n f(n, n x)` brought `Together` (giac's
//!   `recursive_normal`), its numerator and denominator `Factor`ed, each
//!   factor free of `n` to the `x` part and each free of `x` to the `n`
//!   part (giac's `_split`); a factor with both is not a Riemann sum.
//! - The answer is the limit at `n → ∞` of `∫₀¹ (x part) dx` times the
//!   `n` part. An infinite integral with an `n` part whose limit is 0, or
//!   the reverse, is not a Riemann sum either.
//! - Not a Riemann sum: giac's string `"Probably not a Riemann sum"`.

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's answer when the sum is not read as a Riemann sum.
const NOT_RIEMANN: &str = "Probably not a Riemann sum";

/// The `n` part and the `x` part of `p`'s factors; `None` when a factor
/// has both (giac's `_split`).
fn split(
    ev: &mut Eval<'_>,
    p: &Expr,
    n: &Expr,
    x: &Expr,
) -> Result<Option<(Expr, Expr)>, EvalError> {
    let factored = ev.eval(&Expr::call(Symbol::FACTOR, [p.clone()]))?;
    let factors = match factored.head().and_then(Expr::as_symbol) {
        Some(Symbol::TIMES) => factored.args().to_vec(),
        _ => vec![factored],
    };
    let (mut in_n, mut in_x) = (Vec::new(), Vec::new());
    for f in factors {
        if !f.contains(n) {
            in_x.push(f);
        } else if !f.contains(x) {
            in_n.push(f);
        } else {
            return Ok(None);
        }
    }
    Ok(Some((
        Expr::call(Symbol::TIMES, in_n),
        Expr::call(Symbol::TIMES, in_x),
    )))
}

/// Whether `e` has an infinite part (giac's `is_inf`).
fn infinite(e: &Expr) -> bool {
    [
        Symbol::INFINITY,
        Symbol::COMPLEX_INFINITY,
        Symbol::DIRECTED_INFINITY,
    ]
    .into_iter()
    .any(|s| e.contains(&Expr::from(s)))
}

pub(crate) fn sum_riemann(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, vars] = e.args() else {
        return Ok(None);
    };
    let (Some(Symbol::LIST), [n, k]) = (vars.head().and_then(Expr::as_symbol), vars.args()) else {
        return Ok(None);
    };
    let x = Expr::from(Symbol::new("SumRiemann$x"));
    if n.as_symbol().is_none() || k.as_symbol().is_none() || f.contains(&x) {
        return Ok(None);
    }
    let a = n.clone() * f.substitute(&[(k.clone(), n.clone() * x.clone())]);
    let a = ev.eval(&Expr::call(Symbol::TOGETHER, [a]))?;
    let num = ev.eval(&Expr::call(Symbol::NUMERATOR, [a.clone()]))?;
    let den = ev.eval(&Expr::call(Symbol::DENOMINATOR, [a]))?;
    let not_riemann = Ok(Some(Expr::string(NOT_RIEMANN)));
    let Some((n_num, x_num)) = split(ev, &num, n, &x)? else {
        return not_riemann;
    };
    let Some((n_den, x_den)) = split(ev, &den, n, &x)? else {
        return not_riemann;
    };
    let an = ev.eval(&(n_num / n_den))?;
    let ax = ev.eval(&(x_num / x_den))?;
    let bounds = Expr::call(Symbol::LIST, [x.clone(), Expr::from(0), Expr::from(1)]);
    let integral = ev.eval(&Expr::call(Symbol::INTEGRATE, [ax, bounds]))?;
    let infinity = Expr::from(Symbol::INFINITY);
    let limit = ev.eval(&Expr::call(
        Symbol::LIMIT,
        [an.clone(), n.clone(), infinity.clone()],
    ))?;
    let zero = Expr::from(0);
    if (infinite(&integral) && limit == zero) || (integral == zero && infinite(&limit)) {
        return not_riemann;
    }
    let result = ev.eval(&Expr::call(
        Symbol::LIMIT,
        [integral * an, n.clone(), infinity],
    ))?;
    ev.step("SumRiemann", e, &result);
    Ok(Some(result))
}
