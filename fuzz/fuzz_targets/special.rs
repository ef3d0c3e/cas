//! The calculus of the special functions (CAS-DOC-31 §8, CAS-361), one
//! check per input:
//! - an integrand that gives `Ei`, `Si`, `Ci` or `erf`: a Laurent
//!   polynomial over a power of a linear factor times an exponential, sine,
//!   cosine, `sinh` or `cosh` of a linear argument, a polynomial times a
//!   Gaussian, or `x^m/log x`; its antiderivative differentiates back to it,
//!   the difference a ball holding 0 at a point;
//! - a rational summand of terms `c/(k + s)^j`, `s ≥ 0`, put together: its
//!   sum from 1 to `n` is the direct sum at `n = 1 … 8`, by balls;
//! - a series of `Ei`, `Si`, `Ci`, `Erf`, `PolyGamma`, `Gamma` or
//!   `LambertW` to the order 4: its terms are the function `2^-6` away
//!   from the point, within `2^-27` of its size.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn small(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(-6..=6)?, u.int_in_range(1..=3)?))
}

fn nonzero(u: &mut Unstructured<'_>) -> Result<Rational> {
    let r = small(u)?;
    Ok(if r.is_zero() { Rational::ONE } else { r })
}

fn x() -> Expr {
    Expr::from(Symbol::new("x"))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// A random polynomial in `x` of degree up to `most`.
fn polynomial(u: &mut Unstructured<'_>, most: i64) -> Result<Expr> {
    let degree = u.int_in_range(0..=most)?;
    let mut terms = Vec::new();
    for k in 0..=degree {
        terms.push(Expr::from(small(u)?) * x().pow(k));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// `N[e]` at 30 digits as a ball, real or complex: whether it holds 0.
fn holds_zero(ctx: &Context, e: &Expr) -> Option<bool> {
    let value = ctx
        .eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
        .ok()?;
    match value.view() {
        View::Ball(b) => Some(b.contains_rational(&Rational::ZERO)),
        View::ComplexBall(z) => Some(z.may_be_zero()),
        _ => value.to_rational().map(|r| r.is_zero()),
    }
}

fn integrals(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let linear = |u: &mut Unstructured<'_>| -> Result<Expr> {
        Ok(Expr::from(nonzero(u)?) * x() + Expr::from(small(u)?))
    };
    let f = match u.int_in_range(0u8..=2)? {
        0 => {
            let arg = linear(u)?;
            let kernel = match u.int_in_range(0u8..=4)? {
                0 => Expr::from(Symbol::E).pow(arg),
                1 => Expr::call(Symbol::SIN, [arg]),
                2 => Expr::call(Symbol::COS, [arg]),
                3 => call("Sinh", vec![arg]),
                _ => call("Cosh", vec![arg]),
            };
            let x0 = Expr::from(small(u)?);
            let n = u.int_in_range(1..=3i64)?;
            polynomial(u, 3)? * kernel / (x() - x0).pow(n)
        }
        1 => {
            let alpha = -Expr::from(q(u.int_in_range(1..=4)?, u.int_in_range(1..=3)?));
            let q = alpha * x().pow(2) + Expr::from(small(u)?) * x() + Expr::from(small(u)?);
            polynomial(u, 4)? * Expr::from(Symbol::E).pow(q)
        }
        _ => {
            let m = small(u)?;
            if m == q(-1, 1) {
                return Ok(());
            }
            x().pow(Expr::from(m)) / Expr::call(Symbol::LOG, [x()])
        }
    };
    let f = eval(ctx, &f);
    let antiderivative = eval(ctx, &Expr::call(Symbol::INTEGRATE, [f.clone(), x()]));
    if antiderivative.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(());
    }
    let back = eval(ctx, &Expr::call(Symbol::D, [antiderivative.clone(), x()]));
    // At a point away from the poles: x = 37/11.
    let point = Expr::from(q(37, 11));
    let difference = (back - f.clone()).substitute(&[(x(), point)]);
    if let Some(zero) = holds_zero(ctx, &difference) {
        assert!(
            zero,
            "∫ {f:?} = {antiderivative:?} does not differentiate back"
        );
    }
    Ok(())
}

fn sums(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let k = Expr::from(Symbol::new("k"));
    let n = Expr::from(Symbol::new("n"));
    let mut terms = Vec::new();
    for _ in 0..u.int_in_range(1..=3u8)? {
        let c = Expr::from(nonzero(u)?);
        let s = Expr::from(q(u.int_in_range(0..=6)?, u.int_in_range(1..=4)?));
        let j = u.int_in_range(1..=3i64)?;
        terms.push(c / (k.clone() + s).pow(j));
    }
    let summand = eval(
        ctx,
        &Expr::call(Symbol::TOGETHER, [Expr::call(Symbol::PLUS, terms)]),
    );
    let total = eval(
        ctx,
        &Expr::call(
            Symbol::SUM,
            [summand.clone(), k.clone(), Expr::from(1), n.clone()],
        ),
    );
    if total.contains(&Expr::from(Symbol::SUM)) {
        return Ok(());
    }
    let mut direct = Expr::from(0);
    for m in 1..=8i64 {
        direct = eval(
            ctx,
            &(direct + summand.substitute(&[(k.clone(), Expr::from(m))])),
        );
        let at = total.substitute(&[(n.clone(), Expr::from(m))]);
        if let Some(zero) = holds_zero(ctx, &(at - direct.clone())) {
            assert!(zero, "Σ {summand:?} = {total:?} at n = {m}, not {direct:?}");
        }
    }
    Ok(())
}

fn series(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    // The function and a point where its series is known.
    let a = Expr::from(q(u.int_in_range(1..=12)?, 4));
    let (f, point) = match u.int_in_range(0u8..=7)? {
        0 => (call("ExpIntegralEi", vec![x()]), a),
        1 => (call("SinIntegral", vec![x()]), Expr::from(small(u)?)),
        2 => (call("CosIntegral", vec![x()]), a),
        3 => (call("Erf", vec![x()]), Expr::from(small(u)?)),
        4 => (
            call(
                "PolyGamma",
                vec![x(), Expr::from(u.int_in_range(0..=2i64)?)],
            ),
            a,
        ),
        5 => (call("Gamma", vec![x()]), a),
        6 => (
            call("LambertW", vec![x()]),
            Expr::from(q(u.int_in_range(0..=8)?, 4)),
        ),
        _ => (call("ExpIntegralEi", vec![x()]), Expr::from(0)),
    };
    let s = eval(
        ctx,
        &Expr::call(
            Symbol::SERIES,
            [
                f.clone(),
                Expr::call(Symbol::LIST, [x(), point.clone(), Expr::from(4)]),
            ],
        ),
    );
    if s.head() == Some(&Expr::from(Symbol::SERIES)) {
        return Ok(());
    }
    // The terms without the remainder, at x = a + 2^-6, against f there.
    let order_size = Symbol::new("OrderSize");
    let rest = |e: &Expr| e.contains(&Expr::from(order_size));
    let terms: Vec<Expr> = if s.head() == Some(&Expr::from(Symbol::PLUS)) {
        s.args().iter().filter(|t| !rest(t)).cloned().collect()
    } else if rest(&s) {
        Vec::new()
    } else {
        vec![s.clone()]
    };
    let at = eval(ctx, &(point + Expr::from(q(1, 64))));
    let p = Expr::call(Symbol::PLUS, terms).substitute(&[(x(), at.clone())]);
    let value = f.substitute(&[(x(), at)]);
    let difference = ctx
        .eval(&Expr::call(Symbol::N, [p - value.clone(), Expr::from(30)]))
        .ok();
    let size = ctx
        .eval(&Expr::call(
            Symbol::N,
            [Expr::call(Symbol::ABS, [value]), Expr::from(30)],
        ))
        .ok();
    let (Some(difference), Some(size)) = (difference, size) else {
        return Ok(());
    };
    let (View::Ball(d), View::Ball(v)) = (difference.view(), size.view()) else {
        return Ok(());
    };
    // |difference| ≤ 2^-27 (1 + |f|).
    let bound = (v.mid().to_f64().abs() + 1.0) * 2f64.powi(-27);
    let gap = d.mid().to_f64().abs();
    assert!(
        gap <= bound,
        "Series[{f:?}] at {:?}: {s:?} is {gap} from the function",
        x()
    );
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0u8..=2)? {
        0 => integrals(&ctx, u),
        1 => sums(&ctx, u),
        _ => series(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
