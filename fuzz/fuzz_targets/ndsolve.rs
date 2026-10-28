//! `NDSolve` (CAS-DOC-35 §5, CAS-394), one problem per input:
//! - a linear system `y′ = A y`, `A` of small integers, 2 × 2, or 3 × 3
//!   with a 2 × 2 block above a diagonal entry, so that its eigenvalues
//!   are quadratic or rational: `y(t1)` meets `MatrixExp[A t1] y0`;
//! - a scalar equation with a closed form, `y′ = a y + 1`, `y′ = −a t y`,
//!   `y′ = y²` before its pole, or `y′ = a cos(t) y`: `y(t1)` meets the
//!   closed form there;
//! - the answer at the working precision and at 100 bits overlap.
//!
//! A call that stays, as one that takes too many steps, checks nothing; nor
//! does a closed form whose value is no real ball.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::ball::Ball;
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn small(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(-8..=8)?, u.int_in_range(1..=4)?))
}

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// The real balls of `N[e, 40]`, entry by entry if it is a list; `None`
/// where one is no real number. A complex ball whose imaginary part meets 0
/// gives its real part.
fn balls(ctx: &Context, e: &Expr) -> Option<Vec<Ball>> {
    let v = ctx
        .eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(40)]))
        .ok()?;
    let items = if v.head() == Some(&Expr::from(Symbol::LIST)) {
        v.args().to_vec()
    } else {
        vec![v]
    };
    items
        .iter()
        .map(|x| match x.view() {
            View::Ball(b) => Some(b.clone()),
            View::ComplexBall(z) if z.im().overlaps(&Ball::ZERO) => Some(z.re().clone()),
            _ => x.to_rational().map(|r| Ball::from_rational(&r, 200)),
        })
        .collect()
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let t = sym("t");
    let t1 = q(u.int_in_range(1..=8)?, 4);
    let (call, expected) = if u.arbitrary()? {
        // y′ = A y.
        let n = u.int_in_range(2..=3usize)?;
        let ys: Vec<Expr> = ["y", "z", "w"][..n].iter().map(|s| sym(s)).collect();
        let mut rows = Vec::new();
        let mut rhs = Vec::new();
        for i in 0..n {
            let row: Vec<Expr> = (0..n)
                .map(|j| {
                    if i == 2 && j < 2 {
                        Ok(Expr::from(0))
                    } else {
                        u.int_in_range(-2..=2i64).map(Expr::from)
                    }
                })
                .collect::<Result<_>>()?;
            rhs.push(Expr::call(
                Symbol::PLUS,
                row.iter().zip(&ys).map(|(a, y)| a.clone() * y.clone()),
            ));
            rows.push(list(row));
        }
        let y0: Vec<Expr> = (0..n)
            .map(|_| small(u).map(Expr::from))
            .collect::<Result<_>>()?;
        let call = Expr::call(
            Symbol::NDSOLVE,
            [
                list(rhs),
                list(vec![t.clone(), list(ys)]),
                list(vec![Expr::from(0), list(y0.clone())]),
                Expr::from(t1.clone()),
            ],
        );
        let flow = Expr::call(
            Symbol::MATRIX_EXP,
            [Expr::call(Symbol::TIMES, [list(rows), Expr::from(t1)])],
        );
        (call, Expr::call(Symbol::TIMES, [flow, list(y0)]))
    } else {
        let y = sym("y");
        let a = small(u)?;
        let y0 = small(u)?;
        let start = Expr::from(y0.clone());
        let (ae, y0e, t1e) = (Expr::from(a.clone()), start.clone(), Expr::from(t1.clone()));
        let e = |x: Expr| Expr::from(Symbol::E).pow(x);
        let (f, closed) = match u.int_in_range(0u8..=3)? {
            0 => {
                // y′ = a y + 1: (y0 + 1/a) e^(a t) − 1/a, or y0 + t at a = 0.
                let f = ae.clone() * y.clone() + Expr::from(1);
                let closed = if a.is_zero() {
                    y0e + t1e
                } else {
                    let inverse = Expr::from(q(1, 1) / a);
                    (y0e + inverse.clone()) * e(ae * t1e) - inverse
                };
                (f, closed)
            }
            1 => {
                // y′ = −a t y: y0 e^(−a t²/2).
                let f = -(ae.clone() * t.clone() * y.clone());
                (f, y0e * e(-(ae * t1e.pow(2) / Expr::from(2))))
            }
            2 => {
                // y′ = y²: y0/(1 − y0 t), before the pole.
                if y0.is_positive() && q(1, 1) / y0.clone() <= t1 {
                    return Ok(());
                }
                let f = y.clone().pow(2);
                (f, y0e.clone() / (Expr::from(1) - y0e * t1e))
            }
            _ => {
                // y′ = a cos(t) y: y0 e^(a sin t).
                let f = ae.clone() * Expr::call(Symbol::COS, [t.clone()]) * y.clone();
                (f, y0e * e(ae * Expr::call(Symbol::SIN, [t1e])))
            }
        };
        let call = Expr::call(
            Symbol::NDSOLVE,
            [
                f,
                list(vec![t.clone(), y]),
                list(vec![Expr::from(0), start]),
                Expr::from(t1),
            ],
        );
        (call, closed)
    };
    let found = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    if found.head() == Some(&Expr::from(Symbol::NDSOLVE)) {
        return Ok(());
    }
    let answer = &found.args()[0];
    let ours = balls(&ctx, answer).unwrap_or_else(|| panic!("{call:?}: {found:?}"));
    if let Some(theirs) = balls(&ctx, &expected) {
        assert_eq!(ours.len(), theirs.len(), "{call:?}: {found:?}");
        for (a, b) in ours.iter().zip(&theirs) {
            assert!(
                a.overlaps(b),
                "{call:?} is {found:?}, not {expected:?}: {a:?} against {b:?}"
            );
        }
    }
    // At 100 bits too.
    let mut precise = ctx.clone();
    precise.set_working_precision(100);
    let again = precise
        .eval(&call)
        .unwrap_or_else(|e| panic!("{call:?}: {e}"));
    if again.head() != Some(&Expr::from(Symbol::NDSOLVE)) {
        let more = balls(&precise, &again.args()[0])
            .unwrap_or_else(|| panic!("{call:?}: {again:?} at 100 bits"));
        assert_eq!(ours.len(), more.len(), "{call:?}: {again:?} at 100 bits");
        for (a, b) in ours.iter().zip(&more) {
            assert!(a.overlaps(b), "{call:?}: {a:?} at 53 bits, {b:?} at 100");
        }
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
