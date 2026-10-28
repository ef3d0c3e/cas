//! Certified numerical integration (CAS-DOC-26 §6, CAS-326), one check per
//! input: a random `F`, real-analytic on the line, built from `x` and
//! small rationals by `+`, `×`, `exp`, `sin`, `cos`, `atan`, `1/(1 + u²)`,
//! `log(1 + u²)` and `√(1 + u²)`; a random segment `[a, b]` of rationals in
//! `[−3, 3]`; a working precision of 53, 128 or 256 bits.
//! `NIntegrate[D[F, x], x, a, b]` must be a ball that overlaps
//! `F(b) − F(a)`, or stay.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::ball::{Ball, ComplexBall};
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

/// `1 + g²`, positive on the line.
fn positive(g: Expr) -> Expr {
    Expr::from(1) + g.pow(2)
}

/// A random `F`, `depth` deep at most.
fn term(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    if depth == 0 || u.ratio(1, 4)? {
        return Ok(if u.ratio(2, 3)? {
            x()
        } else {
            Expr::from(small(u)?)
        });
    }
    let g = term(u, depth - 1)?;
    Ok(match u.int_in_range(0..=8u8)? {
        0 => g + term(u, depth - 1)?,
        1 => g * term(u, depth - 1)?,
        2 => Expr::call(Symbol::EXP, [g]),
        3 => Expr::call(Symbol::SIN, [g]),
        4 => Expr::call(Symbol::COS, [g]),
        5 => Expr::call(Symbol::ARC_TAN, [g]),
        6 => positive(g).pow(-1),
        7 => Expr::call(Symbol::LOG, [positive(g)]),
        _ => positive(g).pow(Expr::from(q(1, 2))),
    })
}

/// The complex ball of a number.
fn ball(e: &Expr) -> Option<ComplexBall> {
    match e.view() {
        View::Ball(b) => Some(ComplexBall::from_real(b.clone())),
        View::ComplexBall(z) => Some(z.clone()),
        _ => e
            .to_rational()
            .map(|r| ComplexBall::from_real(Ball::from_rational(&r, 512))),
    }
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let mut ctx = Context::new();
    ctx.set_working_precision([53, 128, 256][u.int_in_range(0..=2usize)?]);
    let big = term(u, 3)?;
    let ends = [q(u.int_in_range(-9..=9)?, 3), q(u.int_in_range(-9..=9)?, 3)];
    let Ok(f) = ctx.eval(&Expr::call(Symbol::D, [big.clone(), x()])) else {
        return Ok(());
    };
    let input = Expr::call(
        Symbol::new("NIntegrate"),
        [
            f.clone(),
            x(),
            Expr::from(ends[0].clone()),
            Expr::from(ends[1].clone()),
        ],
    );
    let Ok(found) = ctx.eval(&input) else {
        return Ok(());
    };
    let Some(found) = ball(&found) else {
        // Past the budget: the call stays.
        return Ok(());
    };
    let at = |r: &Rational| big.substitute(&[(x(), Expr::from(r.clone()))]);
    let difference = at(&ends[1]) - at(&ends[0]);
    let Ok(exact) = Context::new().eval(&Expr::call(Symbol::N, [difference, Expr::from(40)]))
    else {
        return Ok(());
    };
    let Some(exact) = ball(&exact) else {
        return Ok(());
    };
    assert!(
        found.overlaps(&exact),
        "{input:?} of F = {big:?}: {found:?}, but F(b) − F(a) is {exact:?}"
    );
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
