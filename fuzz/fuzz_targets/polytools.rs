//! Polynomial tools (CAS-DOC-29 §4, CAS-342), one check per input:
//! - a rational function with a random numerator and a denominator of
//!   chosen linear and quadratic factors: `Apart` sums back to it;
//! - a random polynomial: `FromCoefficients` of `CoefficientList` is it,
//!   expanded; `CountRoots` over the line is the number of `RealRoots`; and
//!   the square-free part changes sign across each root's ball.
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

fn apart(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let num = polynomial(u, 4)?;
    let mut factors = Vec::new();
    for _ in 0..u.int_in_range(1..=3u8)? {
        let factor = if u.ratio(1, 2)? {
            x() - Expr::from(small(u)?)
        } else {
            x().pow(2) + Expr::from(small(u)?) * x() + Expr::from(small(u)?)
        };
        factors.push(factor.pow(u.int_in_range(1..=2i64)?));
    }
    let f = num / Expr::call(Symbol::TIMES, factors);
    let input = call("Apart", vec![f.clone()]);
    let found = eval(ctx, &input);
    if found.contains(&Expr::from(Symbol::new("Apart"))) {
        return Ok(());
    }
    let difference = eval(
        ctx,
        &Expr::call(
            Symbol::TOGETHER,
            [Expr::call(Symbol::EXPAND, [found.clone() - f])],
        ),
    );
    assert_eq!(difference, Expr::from(0), "{input:?} = {found:?}");
    Ok(())
}

fn roots(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let p = eval(ctx, &Expr::call(Symbol::EXPAND, [polynomial(u, 6)?]));
    let list = eval(ctx, &call("CoefficientList", vec![p.clone()]));
    let back = eval(ctx, &call("FromCoefficients", vec![list.clone()]));
    assert_eq!(back, p, "FromCoefficients[{list:?}]");
    if p.to_rational().is_some() {
        return Ok(());
    }
    let infinity = Expr::from(Symbol::INFINITY);
    let count = eval(
        ctx,
        &call(
            "CountRoots",
            vec![p.clone(), x(), -infinity.clone(), infinity],
        ),
    );
    let found = eval(ctx, &call("RealRoots", vec![p.clone()]));
    if found.contains(&Expr::from(Symbol::new("RealRoots"))) {
        return Ok(());
    }
    assert_eq!(
        count,
        Expr::from(i64::try_from(found.args().len()).expect("small")),
        "{p:?}: CountRoots {count:?}, RealRoots {found:?}"
    );
    // The square-free part changes sign across each ball.
    let derivative = eval(ctx, &Expr::call(Symbol::D, [p.clone(), x()]));
    let g = eval(ctx, &call("GCD", vec![p.clone(), derivative]));
    let free = eval(ctx, &call("PolynomialQuotient", vec![p.clone(), g, x()]));
    for pair in found.args() {
        let View::Ball(b) = pair.args()[0].view() else {
            panic!("{found:?}");
        };
        let (lo, hi) = (b.mid().to_rational().clone(), b.mid().to_rational().clone());
        let radius = b.rad().to_float().expect("finite").to_rational();
        let at = |v: Rational| {
            eval(ctx, &free.substitute(&[(x(), Expr::from(v))]))
                .to_rational()
                .expect("rational")
        };
        let (left, right) = (at(&lo - &radius), at(&hi + &radius));
        assert!(
            left.is_zero() || right.is_zero() || left.is_negative() != right.is_negative(),
            "{p:?}: no sign change across {b:?}"
        );
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    if u.ratio(1, 2)? {
        apart(&ctx, u)
    } else {
        roots(&ctx, u)
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
