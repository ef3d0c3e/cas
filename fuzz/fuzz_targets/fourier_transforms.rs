//! Fourier transforms (CAS-DOC-41 §6, CAS-443), one function per input, a
//! sum of one or two terms of a family:
//! - Gaussians `c xⁿ e^(−a (x − b)²)`, `a` in {1/2, 1, 2}, `n ≤ 2`,
//!   modulated by `cos(k x)`, `sin(k x)` or `e^(i k x)`;
//! - two-sided exponentials `c e^(−a |p x + q|)`, `a |p| ≥ 1`;
//! - one-sided `c θ(x − q) e^(−a x)`, `a` in {1, 2};
//! - rational functions `P/Q`, `Q = (x² + α²)((x − r)² + β²)`, `deg P ≤ 3`.
//!
//! Where the transform comes back free of `FourierTransform`, it agrees at
//! a random `s₀` with `NIntegrate` of `f(x) cos(s₀ x)` and `−f(x) sin(s₀
//! x)` over `[−L, L]`, the tails past `L` below `10⁻²⁰` (not for the
//! rational ones, whose tails are long), and `InverseFourierTransform`
//! gives `f` back at points.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::ball::{Ball, Round};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn f1(head: Symbol, u: Expr) -> Expr {
    Expr::call(head, [u])
}

fn small(u: &mut Unstructured<'_>, lo: i64, hi: i64, d: i64) -> Result<Expr> {
    Ok(Expr::from(q(u.int_in_range(lo..=hi)?, d)))
}

/// A term of a family, and the `L` past which its tails are negligible;
/// `None` for a rational function.
fn term(u: &mut Unstructured<'_>, family: u8) -> Result<(Expr, Option<i64>)> {
    let x = sym("x");
    let c = small(u, -3, 3, 2)?;
    Ok(match family {
        0 => {
            let a = Expr::from(u.choose(&[q(1, 2), q(1, 1), q(2, 1)])?.clone());
            let b = small(u, -4, 4, 2)?;
            let n = u.int_in_range(0..=2i64)?;
            let gaussian = f1(Symbol::EXP, -(a * (x.clone() - b).pow(2)));
            let k = Expr::from(u.int_in_range(0..=2i64)?);
            let oscillation = match u.int_in_range(0u8..=3)? {
                0 => f1(Symbol::COS, k * x.clone()),
                1 => f1(Symbol::SIN, k * x.clone()),
                2 => f1(Symbol::EXP, Expr::from(Symbol::I) * k * x.clone()),
                _ => Expr::from(1),
            };
            (c * x.pow(n) * gaussian * oscillation, Some(16))
        }
        1 => {
            let a = Expr::from(u.int_in_range(1..=2i64)?);
            let p = Expr::from(*u.choose(&[1i64, 2, -1])?);
            let shift = small(u, -4, 4, 2)?;
            let abs = f1(Symbol::ABS, p * x.clone() + shift);
            (c * f1(Symbol::EXP, -(a * abs)), Some(50))
        }
        2 => {
            let a = Expr::from(u.int_in_range(1..=2i64)?);
            let cut = small(u, -4, 4, 2)?;
            let step = f1(Symbol::HEAVISIDE_THETA, x.clone() - cut);
            (c * step * f1(Symbol::EXP, -(a * x.clone())), Some(50))
        }
        _ => {
            let alpha = Expr::from(u.int_in_range(1..=3i64)?);
            let beta = Expr::from(u.int_in_range(1..=3i64)?);
            let r = small(u, -2, 2, 1)?;
            let den = (x.clone().pow(2) + alpha.pow(2)) * ((x.clone() - r).pow(2) + beta.pow(2));
            let mut num = Vec::new();
            for k in 0..=u.int_in_range(0..=3i64)? {
                num.push(small(u, -3, 3, 1)? * x.clone().pow(k));
            }
            (Expr::call(Symbol::PLUS, num) / den, None)
        }
    })
}

/// The number `e` as a complex pair of balls; `None` for another value.
fn complex(ctx: &Context, e: &Expr) -> Option<(Ball, Ball)> {
    let v = ctx
        .eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
        .ok()?;
    let (re, im) = match v.view() {
        View::Ball(b) => (b.clone(), Ball::ZERO),
        View::ComplexBall(z) => (z.re().clone(), z.im().clone()),
        _ => (Ball::from_rational(&v.to_rational()?, 128), Ball::ZERO),
    };
    (re.is_finite() && im.is_finite()).then_some((re, im))
}

/// Whether two balls overlap, or lie within `10⁻¹²` of each other: the
/// tails past `L` and `NIntegrate`'s own rounding.
fn close(a: &Ball, b: &Ball) -> bool {
    if a.overlaps(b) {
        return true;
    }
    let slack = Ball::from_rational(&q(1, 1_000_000_000_000), 128);
    let room = slack.sub(&a.sub(b, 128).abs(), 128);
    let Some(r) = room.rad().to_float() else {
        return false;
    };
    let lo = room.mid().sub_round(&r, 64, Round::Floor).0;
    !lo.is_negative() && !lo.is_zero()
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let (x, s) = (sym("x"), sym("s"));
    let family = u.int_in_range(0u8..=3)?;
    let (mut f, mut tail) = term(u, family)?;
    if u.arbitrary()? {
        let (g, other) = term(u, family)?;
        f = f + g;
        tail = tail.max(other);
    }
    let f = ctx.eval(&f).unwrap_or_else(|e| panic!("{f:?}: {e}"));
    let call = Expr::call(Symbol::FOURIER_TRANSFORM, [f.clone(), x.clone(), s.clone()]);
    let image = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    if image.contains(&Expr::from(Symbol::FOURIER_TRANSFORM))
        || image.contains(&Expr::from(Symbol::DIRAC_DELTA))
    {
        return Ok(());
    }
    let s0 = Expr::from(q(u.int_in_range(-24..=24)?, 8));
    let Some((re, im)) = complex(&ctx, &image.substitute(&[(s.clone(), s0.clone())])) else {
        return Ok(());
    };
    // Against NIntegrate where the tails are short.
    if let (Some(l), true) = (tail, family != 3) {
        let bounds = Expr::call(Symbol::LIST, [x.clone(), Expr::from(-l), Expr::from(l)]);
        let kernel = |h: Symbol| f.clone() * f1(h, s0.clone() * x.clone());
        let integrals = [
            Expr::call(Symbol::NINTEGRATE, [kernel(Symbol::COS), bounds.clone()]),
            Expr::call(Symbol::NINTEGRATE, [-kernel(Symbol::SIN), bounds]),
        ];
        let mut values = Vec::new();
        for integral in &integrals {
            let Some((value, _)) = complex(&ctx, integral) else {
                return Ok(());
            };
            values.push(value);
        }
        assert!(
            close(&re, &values[0]) && close(&im, &values[1]),
            "{call:?} is {image:?}, at {s0:?} {re:?} + i {im:?}, but NIntegrate gives {:?} + i {:?}",
            values[0],
            values[1]
        );
    }
    // The way back.
    let back = Expr::call(
        Symbol::INVERSE_FOURIER_TRANSFORM,
        [image.clone(), s, x.clone()],
    );
    let back = ctx.eval(&back).unwrap_or_else(|e| panic!("{back:?}: {e}"));
    if back.contains(&Expr::from(Symbol::INVERSE_FOURIER_TRANSFORM))
        || back.contains(&Expr::from(Symbol::FOURIER_TRANSFORM))
    {
        return Ok(());
    }
    for _ in 0..2 {
        let p = Expr::from(q(u.int_in_range(-24..=24)?, 8));
        let at = |e: &Expr| complex(&ctx, &e.substitute(&[(x.clone(), p.clone())]));
        let (Some((a_re, a_im)), Some((b_re, b_im))) = (at(&back), at(&f)) else {
            continue;
        };
        assert!(
            close(&a_re, &b_re) && close(&a_im, &b_im),
            "the inverse of {image:?}, the transform of {f:?}, is {back:?}: at {p:?} not f"
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
