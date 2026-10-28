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
//!   `LambertW` to the order 4: its terms are the function up to a
//!   remainder of order 5, which halving the step a sixteenth of the way to
//!   the nearest singularity divides by 24 at least.
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
    // The function, a point where its series is known, and the distance
    // from it to the nearest singularity, at most 1/2.
    let a = q(u.int_in_range(1..=12)?, 4);
    let half = q(1, 2);
    let near = |r: Rational| if r < half { r } else { half.clone() };
    let (f, point, reach) = match u.int_in_range(0u8..=7)? {
        0 => (call("ExpIntegralEi", vec![x()]), a.clone(), near(a)),
        1 => (call("SinIntegral", vec![x()]), small(u)?, half.clone()),
        2 => (call("CosIntegral", vec![x()]), a.clone(), near(a)),
        3 => (call("Erf", vec![x()]), small(u)?, half.clone()),
        4 => (
            call(
                "PolyGamma",
                vec![x(), Expr::from(u.int_in_range(0..=2i64)?)],
            ),
            a.clone(),
            near(a),
        ),
        5 => (call("Gamma", vec![x()]), a.clone(), near(a)),
        6 => {
            // W's branch point −1/e, past −3/8.
            let b = q(u.int_in_range(0..=8)?, 4);
            let reach = near(&b + &q(3, 8));
            (call("LambertW", vec![x()]), b, reach)
        }
        _ => (
            call("ExpIntegralEi", vec![x()]),
            Rational::ZERO,
            half.clone(),
        ),
    };
    let s = eval(
        ctx,
        &Expr::call(
            Symbol::SERIES,
            [
                f.clone(),
                Expr::call(
                    Symbol::LIST,
                    [x(), Expr::from(point.clone()), Expr::from(4)],
                ),
            ],
        ),
    );
    if s.head() == Some(&Expr::from(Symbol::SERIES)) {
        return Ok(());
    }
    // The terms without the remainder.
    let order_size = Symbol::new("OrderSize");
    let rest = |e: &Expr| e.contains(&Expr::from(order_size));
    let terms: Vec<Expr> = if s.head() == Some(&Expr::from(Symbol::PLUS)) {
        s.args().iter().filter(|t| !rest(t)).cloned().collect()
    } else if rest(&s) {
        Vec::new()
    } else {
        vec![s.clone()]
    };
    let polynomial = Expr::call(Symbol::PLUS, terms);
    // |f − terms| at x = a + h, h a sixteenth of the reach, as a float.
    let gap = |h: &Rational| -> Option<f64> {
        let at = Expr::from(&point + h);
        let p = polynomial.substitute(&[(x(), at.clone())]);
        let value = f.substitute(&[(x(), at)]);
        let difference = ctx
            .eval(&Expr::call(Symbol::N, [p - value, Expr::from(30)]))
            .ok()?;
        match difference.view() {
            View::Ball(d) => Some(d.mid().to_f64().abs()),
            _ => None,
        }
    };
    let h = &reach / &q(16, 1);
    let (Some(far), Some(close)) = (gap(&h), gap(&(&h / &q(2, 1)))) else {
        return Ok(());
    };
    // The remainder of order 5: halving h divides it by about 32, where a
    // wrong term of order 4 or less divides it by 16 at most.
    if far > 1e-20 {
        assert!(
            close <= far / 24.0,
            "Series[{f:?}] at {point:?}: {s:?} is {far} and {close} from the function"
        );
    }
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
