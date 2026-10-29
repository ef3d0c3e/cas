//! Integration II (CAS-DOC-40 §6, CAS-437), one integrand per input from a
//! family the substitutions rationalise:
//! - rational functions of `e^(k x)`, `k` in {1, 2, −1, 1/2}, or of `sinh x`
//!   and `cosh x`;
//! - rational functions of `sin v` and `cos v`, `v = a x + b`, with `tan v`
//!   and `sec v` among them;
//! - rational functions of `x` and `(a x + b)^(1/n)`, `n` 2 or 3;
//! - `x^j q^(k/2)`, `k` odd, for a quadratic `q` whose leading coefficient
//!   or constant is a square, or whose roots are rational.
//!
//! An antiderivative that comes back whole has the integrand as its
//! derivative at three points where both are real, by overlap of balls; a
//! definite integral over an interval that comes back a number agrees with
//! `NIntegrate`'s ball. Points at poles, where the values are no finite
//! numbers, check nothing.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::ball::Ball;
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn x() -> Expr {
    Expr::from(Symbol::new("x"))
}

fn f(head: Symbol, u: Expr) -> Expr {
    Expr::call(head, [u])
}

fn integer(u: &mut Unstructured<'_>, lo: i64, hi: i64) -> Result<Expr> {
    Ok(Expr::from(u.int_in_range(lo..=hi)?))
}

/// `Σ c_i g^i`, `i ≤ 2`, small integer coefficients.
fn polynomial(u: &mut Unstructured<'_>, g: &Expr) -> Result<Expr> {
    let mut terms = Vec::new();
    for i in 0..=2 {
        terms.push(integer(u, -3, 3)? * g.clone().pow(i));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// A polynomial of degree at most 2 in `s` and `c`.
fn polynomial2(u: &mut Unstructured<'_>, s: &Expr, c: &Expr) -> Result<Expr> {
    let mut terms = Vec::new();
    for (i, j) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)] {
        terms.push(integer(u, -3, 3)? * s.clone().pow(i) * c.clone().pow(j));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// A quotient of two of `poly`'s, the denominator not 0.
fn quotient(
    u: &mut Unstructured<'_>,
    poly: &dyn Fn(&mut Unstructured<'_>) -> Result<Expr>,
) -> Result<Expr> {
    let num = poly(u)?;
    let mut den = poly(u)?;
    if den == Expr::from(0) || u.arbitrary()? {
        den = Expr::from(1) + den.pow(2);
    }
    Ok(num / den)
}

fn exponential(u: &mut Unstructured<'_>) -> Result<Expr> {
    if u.arbitrary()? {
        let (s, c) = (f(Symbol::SINH, x()), f(Symbol::COSH, x()));
        return quotient(u, &|u| polynomial2(u, &s, &c));
    }
    let k = u.choose(&[q(1, 1), q(2, 1), q(-1, 1), q(1, 2)])?.clone();
    let g = Expr::from(Symbol::E).pow(Expr::from(k) * x());
    quotient(u, &|u| polynomial(u, &g))
}

fn trigonometric(u: &mut Unstructured<'_>) -> Result<Expr> {
    let a = u.choose(&[q(1, 1), q(2, 1), q(1, 2), q(-1, 1)])?.clone();
    let v = Expr::from(a) * x() + integer(u, -1, 1)?;
    let (s, c) = (f(Symbol::SIN, v.clone()), f(Symbol::COS, v.clone()));
    let e = quotient(u, &|u| polynomial2(u, &s, &c))?;
    Ok(match u.int_in_range(0u8..=3)? {
        0 => e * f(Symbol::TAN, v),
        1 => e + f(Symbol::SEC, v),
        _ => e,
    })
}

fn root(u: &mut Unstructured<'_>) -> Result<Expr> {
    let n = u.int_in_range(2..=3i64)?;
    let a = integer(u, 1, 3)?;
    let r = (a * x() + integer(u, -2, 2)?).pow(Expr::from(q(1, n)));
    let g = x() + r;
    quotient(u, &|u| polynomial(u, &g))
}

fn euler(u: &mut Unstructured<'_>) -> Result<Expr> {
    let quadratic = match u.int_in_range(0u8..=2)? {
        // A square leading coefficient.
        0 => {
            let a = integer(u, 1, 2)?.pow(2);
            a * x().pow(2) + integer(u, -3, 3)? * x() + integer(u, -4, 4)?
        }
        // A square constant.
        1 => {
            let c = integer(u, 1, 3)?.pow(2);
            integer(u, -3, 3)? * x().pow(2) + integer(u, -3, 3)? * x() + c
        }
        // Rational roots.
        _ => {
            let a = integer(u, 1, 3)? * if u.arbitrary()? { -1 } else { 1 };
            let r = Expr::from(q(u.int_in_range(-4..=4)?, 2));
            let r2 = Expr::from(q(u.int_in_range(-4..=4)?, 2));
            a * (x() - r) * (x() - r2)
        }
    };
    let k = *u.choose(&[-3i64, -1, 1, 3])?;
    let j = u.int_in_range(-2..=2i64)?;
    Ok(x().pow(j) * quadratic.pow(Expr::from(q(k, 2))))
}

/// The value at `x = p`, a real ball; `None` for another value.
fn value(ctx: &Context, e: &Expr, p: &Expr) -> Option<Ball> {
    let at = e.substitute(&[(x(), p.clone())]);
    let v = ctx
        .eval(&Expr::call(Symbol::N, [at, Expr::from(30)]))
        .ok()?;
    match v.view() {
        View::Ball(b) if b.is_finite() => Some(b.clone()),
        _ => None,
    }
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let integrand = match u.int_in_range(0u8..=3)? {
        0 => exponential(u)?,
        1 => trigonometric(u)?,
        2 => root(u)?,
        _ => euler(u)?,
    };
    let integrand = ctx.eval(&integrand).unwrap_or_else(|e| panic!("{e}"));
    let integral = Expr::call(Symbol::INTEGRATE, [integrand.clone(), x()]);
    let found = ctx
        .eval(&integral)
        .unwrap_or_else(|e| panic!("{integral:?}: {e}"));
    if found.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(());
    }
    let derivative = ctx
        .eval(&Expr::call(Symbol::D, [found.clone(), x()]))
        .unwrap_or_else(|e| panic!("D[{found:?}]: {e}"));
    for _ in 0..3 {
        let p = Expr::from(q(u.int_in_range(-24..=24)?, 8));
        let (Some(d), Some(g)) = (value(&ctx, &derivative, &p), value(&ctx, &integrand, &p)) else {
            continue;
        };
        assert!(
            d.overlaps(&g),
            "{integral:?} is {found:?}, whose derivative at {p:?} is {d:?}, not {g:?}"
        );
    }
    // A definite integral against NIntegrate.
    let lo = Expr::from(q(u.int_in_range(-16..=16)?, 8));
    let hi = Expr::from(q(u.int_in_range(1..=16)?, 8)) + lo.clone();
    let definite = Expr::call(
        Symbol::INTEGRATE,
        [integrand.clone(), x(), lo.clone(), hi.clone()],
    );
    let exact = ctx
        .eval(&definite)
        .unwrap_or_else(|e| panic!("{definite:?}: {e}"));
    if exact.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(());
    }
    let Some(exact_ball) = value(&ctx, &exact, &Expr::from(0)) else {
        return Ok(());
    };
    let numeric = Expr::call(
        Symbol::NINTEGRATE,
        [integrand, Expr::call(Symbol::LIST, [x(), lo, hi])],
    );
    let Ok(numeric_value) = ctx.eval(&numeric) else {
        return Ok(());
    };
    if let View::Ball(b) = numeric_value.view()
        && b.is_finite()
    {
        assert!(
            exact_ball.overlaps(b),
            "{definite:?} is {exact:?}, {exact_ball:?}, but NIntegrate gives {b:?}"
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
