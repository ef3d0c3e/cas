//! The calculus of the Bessel and Airy functions (CAS-DOC-34 §7, CAS-386),
//! one check per input:
//! - a Bessel function of order `n ≤ 8` or one of Airy's, of `a x`, and
//!   its derivatives by `D`, satisfy its differential equation at a point,
//!   the residual a ball holding 0;
//! - an integrand `x^(n+1) C_n(a x)`, `x^(1−n) C_n(a x)`, `C_1(a x + b)`,
//!   `x Ai(a x)` or `Ai′(a x + b)`: its antiderivative, where found,
//!   differentiates back to it at a point;
//! - a series of the functions to order 4 at a regular point: its terms are
//!   the function up to a remainder of order 5, which halving the step
//!   divides by 24 at least;
//! - Bessel's and Airy's equations with random parameters: `DSolve`'s
//!   answer satisfies them at a point.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
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

/// A rational in `(0, 3]`, not 0.
fn positive(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(1..=12)?, u.int_in_range(1..=4)?))
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

/// A Bessel function's name, and the sign of `x²` in its equation.
fn kind(u: &mut Unstructured<'_>) -> Result<(&'static str, i64)> {
    Ok(*u.choose(&[
        ("BesselJ", 1),
        ("BesselY", 1),
        ("BesselI", -1),
        ("BesselK", -1),
    ])?)
}

fn equations(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let a = positive(u)?;
    let ax = Expr::from(a.clone()) * x();
    let point = Expr::from(q(u.int_in_range(1..=40)?, 7));
    let d = |f: &Expr, k: i64| -> Expr {
        if k == 1 {
            Expr::call(Symbol::D, [f.clone(), x()])
        } else {
            Expr::call(Symbol::D, [f.clone(), x(), x()])
        }
    };
    let residual = if u.arbitrary()? {
        // x² y″ + x y′ + (± a² x² − n²) y for y = C_n(a x).
        let (name, sign) = kind(u)?;
        let n = u.int_in_range(0..=8i64)?;
        let y = call(name, vec![Expr::from(n), ax]);
        let square = Expr::from(&a * &a) * x().pow(2);
        x().pow(2) * d(&y, 2) + x() * d(&y, 1) + (Expr::from(sign) * square - Expr::from(n * n)) * y
    } else {
        // y″ − a³ x y for y = Ai(a x) or Bi(a x), and its derivative's.
        let name = *u.choose(&["AiryAi", "AiryBi"])?;
        let y = call(name, vec![ax]);
        let cube = Expr::from(&(&a * &a) * &a);
        d(&y, 2) - cube * x() * y
    };
    let at = eval(ctx, &residual).substitute(&[(x(), point)]);
    if let Some(zero) = holds_zero(ctx, &at) {
        assert!(zero, "the equation's residual {residual:?} is {at:?}");
    }
    Ok(())
}

fn integrals(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let a = Expr::from(positive(u)?);
    let b = Expr::from(q(u.int_in_range(-6..=6)?, 2));
    let f = match u.int_in_range(0u8..=4)? {
        0 => {
            let (name, _) = kind(u)?;
            let n = u.int_in_range(0..=5i64)?;
            x().pow(n + 1) * call(name, vec![Expr::from(n), a * x()])
        }
        1 => {
            let (name, _) = kind(u)?;
            let n = u.int_in_range(1..=5i64)?;
            x().pow(1 - n) * call(name, vec![Expr::from(n), a * x()])
        }
        2 => {
            let (name, _) = kind(u)?;
            call(name, vec![Expr::from(1), a * x() + b])
        }
        3 => x() * call(*u.choose(&["AiryAi", "AiryBi"])?, vec![a * x()]),
        _ => call(
            *u.choose(&["AiryAiPrime", "AiryBiPrime"])?,
            vec![a * x() + b],
        ),
    };
    let f = eval(ctx, &f);
    let antiderivative = eval(ctx, &Expr::call(Symbol::INTEGRATE, [f.clone(), x()]));
    if antiderivative.contains(&Expr::from(Symbol::INTEGRATE)) {
        return Ok(());
    }
    let back = eval(ctx, &Expr::call(Symbol::D, [antiderivative.clone(), x()]));
    // At a point where every function is real: x = 37/11.
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

fn series(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    // A regular point and the distance to the nearest singularity, 0 for
    // the Bessel functions, at most 1/2.
    let x0 = q(u.int_in_range(1..=12)?, 4);
    let half = q(1, 2);
    let (f, point, reach) = if u.arbitrary()? {
        let (name, _) = kind(u)?;
        let n = u.int_in_range(0..=4i64)?;
        let reach = if x0 < half { x0.clone() } else { half.clone() };
        (call(name, vec![Expr::from(n), x()]), x0, reach)
    } else {
        let name = *u.choose(&["AiryAi", "AiryBi", "AiryAiPrime", "AiryBiPrime"])?;
        let point = q(u.int_in_range(-12..=12)?, 4);
        (call(name, vec![x()]), point, half.clone())
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
    if far > 1e-20 {
        assert!(
            close <= far / 24.0,
            "Series[{f:?}] at {point:?}: {s:?} is {far} and {close} from the function"
        );
    }
    Ok(())
}

fn equations_solved(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let y = Expr::from(Symbol::new("y"));
    let yx = Expr::apply(&y, &[x()]);
    let d = |k: i64| -> Expr {
        if k == 1 {
            Expr::call(Symbol::D, [yx.clone(), x()])
        } else {
            Expr::call(Symbol::D, [yx.clone(), x(), x()])
        }
    };
    let lhs = if u.arbitrary()? {
        let a = positive(u)?;
        let sign = if u.arbitrary()? { 1 } else { -1 };
        let n = u.int_in_range(0..=6i64)?;
        let square = Expr::from(sign) * Expr::from(&a * &a) * x().pow(2);
        x().pow(2) * d(2) + x() * d(1) + (square - Expr::from(n * n)) * yx.clone()
    } else {
        let k = Expr::from(q(u.int_in_range(-24..=24)?, u.int_in_range(1..=3)?));
        d(2) - k * x() * yx.clone()
    };
    let solution = eval(
        ctx,
        &Expr::call(
            Symbol::DSOLVE,
            [
                Expr::call(Symbol::EQUATION, [lhs.clone(), Expr::from(0)]),
                yx.clone(),
            ],
        ),
    );
    if solution.head() == Some(&Expr::from(Symbol::DSOLVE)) {
        return Ok(());
    }
    // The constants at numbers, and the residual at x = 9/7.
    let constants = [
        (Expr::from(Symbol::new("c_0")), Expr::from(q(7, 5))),
        (Expr::from(Symbol::new("c_1")), Expr::from(q(-2, 3))),
    ];
    let function = Expr::call(Symbol::FUNCTION, [x(), solution.clone()]);
    let residual = lhs.substitute(&[(y, function)]);
    let at = eval(ctx, &residual)
        .substitute(&constants)
        .substitute(&[(x(), Expr::from(q(9, 7)))]);
    if let Some(zero) = holds_zero(ctx, &at) {
        assert!(
            zero,
            "DSolve of {lhs:?} gave {solution:?}: the residual is {at:?}"
        );
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0u8..=3)? {
        0 => equations(&ctx, u),
        1 => integrals(&ctx, u),
        2 => series(&ctx, u),
        _ => equations_solved(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
