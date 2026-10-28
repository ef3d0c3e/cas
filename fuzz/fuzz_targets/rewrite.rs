//! The rewriting forms (CAS-DOC-36 §9, CAS-403), one command per input: its
//! answer has its input's value at a point, by overlap of balls, the parts
//! of complex values apart.
//!
//! The symbols `x` and `y` take values in `(0, 1)`, and each family of
//! commands gets expressions built where its identities hold:
//! - exponentials and powers: hyperbolic and inverse hyperbolic functions
//!   (`acosh` above 1, `atanh` within `(−1, 1)`), exponentials of sums
//!   with logarithms and integer multiples, powers of positive bases;
//! - logarithms (formal rules): sums of rational multiples of logarithms
//!   of positive products and powers;
//! - trigonometric forms: polynomials and quotients in `sin`, `cos` and
//!   `tan` of `k x + c y`;
//! - inverse functions: `asin` and `acos` within `(−1, 1)`, `atan` of a
//!   real;
//! - complex parts: `x + i y` in powers, quotients, exponentials,
//!   logarithms and the circular and hyperbolic functions, and quotients
//!   with square roots;
//! - steps: `sign` and `θ` of linear forms.
//!
//! A value that is no finite number, at a pole or past a limit, checks
//! nothing.
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

fn y() -> Expr {
    Expr::from(Symbol::new("y"))
}

fn f(head: Symbol, u: Expr) -> Expr {
    Expr::call(head, [u])
}

/// A rational `p/q`, `|p| ≤ 6`, `1 ≤ q ≤ 4`.
fn small(u: &mut Unstructured<'_>) -> Result<Expr> {
    Ok(Expr::from(q(
        u.int_in_range(-6..=6)?,
        u.int_in_range(1..=4)?,
    )))
}

/// A rational in `(−1, 1)`.
fn unit_number(u: &mut Unstructured<'_>) -> Result<Expr> {
    let d = u.int_in_range(2..=8)?;
    Ok(Expr::from(q(u.int_in_range(1 - d..=d - 1)?, d)))
}

/// `k x + c y`, `k` an integer and `c` a small rational.
fn angle(u: &mut Unstructured<'_>) -> Result<Expr> {
    let k = Expr::from(u.int_in_range(-3..=3i64)?);
    Ok(k * x() + small(u)? * y())
}

/// `c x^i y^j`, `|c| < 1`: within `(−1, 1)` at the point.
fn unit(u: &mut Unstructured<'_>) -> Result<Expr> {
    let (i, j) = (u.int_in_range(0..=2i64)?, u.int_in_range(0..=2i64)?);
    Ok(unit_number(u)? * x().pow(i) * y().pow(j))
}

/// A positive expression: a monomial with a positive coefficient, `1 +
/// x²`, an exponential or a hyperbolic cosine.
fn positive(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    Ok(match u.int_in_range(0u8..=4)? {
        0 => {
            let c = Expr::from(q(u.int_in_range(1..=6)?, u.int_in_range(1..=4)?));
            let (i, j) = (u.int_in_range(-2..=3i64)?, u.int_in_range(-2..=3i64)?);
            c * x().pow(i) * y().pow(j)
        }
        1 => Expr::from(1) + x().pow(2) * small(u)?.pow(2),
        2 if depth > 0 => f(Symbol::EXP, real(u, depth - 1)?),
        3 if depth > 0 => f(Symbol::COSH, real(u, depth - 1)?),
        _ => x() + y(),
    })
}

/// A real expression of `x` and `y`.
fn real(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    if depth == 0 {
        return Ok(match u.int_in_range(0u8..=2)? {
            0 => x(),
            1 => y(),
            _ => small(u)?,
        });
    }
    let d = depth - 1;
    Ok(match u.int_in_range(0u8..=12)? {
        0 => real(u, d)? + real(u, d)?,
        1 => real(u, d)? * real(u, d)?,
        2 => real(u, d)?.pow(u.int_in_range(-2..=3i64)?),
        3 => f(Symbol::SIN, angle(u)?),
        4 => f(Symbol::COS, angle(u)?),
        5 => f(Symbol::SINH, real(u, d)?),
        6 => f(Symbol::TANH, real(u, d)?),
        7 => f(Symbol::EXP, real(u, d)?),
        8 => f(Symbol::LOG, positive(u, d)?),
        9 => f(Symbol::ARC_TAN, real(u, d)?),
        10 => f(Symbol::ARC_SIN, unit(u)?),
        11 => positive(u, d)?.pow(small(u)?),
        _ => small(u)?,
    })
}

/// An expression for the exponential and power rewritings.
fn exponential(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    if depth == 0 {
        return real(u, 0);
    }
    let d = depth - 1;
    Ok(match u.int_in_range(0u8..=10)? {
        0 => exponential(u, d)? + exponential(u, d)?,
        1 => exponential(u, d)? * exponential(u, d)?,
        2 => f(Symbol::SINH, real(u, d)?),
        3 => f(Symbol::COSH, real(u, d)?),
        4 => f(Symbol::TANH, real(u, d)?),
        5 => f(Symbol::ARC_SINH, real(u, d)?),
        6 => f(Symbol::ARC_COSH, Expr::from(1) + positive(u, d)?),
        7 => f(Symbol::ARC_TANH, unit(u)?),
        8 => {
            // e^(Σ c ln b + k u): logarithms and integer multiples.
            let mut exponent = small(u)? * f(Symbol::LOG, positive(u, d)?);
            exponent = exponent + Expr::from(u.int_in_range(-3..=3i64)?) * real(u, d)?;
            if u.arbitrary()? {
                f(Symbol::EXP, exponent)
            } else {
                Expr::from(Symbol::E).pow(exponent)
            }
        }
        9 => positive(u, d)?.pow(real(u, d)?),
        _ => f(Symbol::EXP, angle(u)?),
    })
}

/// A sum of rational multiples of logarithms of positive products.
fn logarithmic(u: &mut Unstructured<'_>) -> Result<Expr> {
    let n = u.int_in_range(1..=4)?;
    let mut terms = Vec::new();
    for _ in 0..n {
        let mut factors = Vec::new();
        for _ in 0..u.int_in_range(1..=3)? {
            factors.push(positive(u, 1)?.pow(small(u)?));
        }
        let argument = Expr::call(Symbol::TIMES, factors);
        terms.push(small(u)? * f(Symbol::LOG, argument));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// A polynomial or quotient in `sin`, `cos` and `tan` of angles.
fn trigonometric(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    if depth == 0 {
        let head = *u.choose(&[Symbol::SIN, Symbol::COS, Symbol::TAN])?;
        let g = f(head, angle(u)?);
        return Ok(small(u)? * g.pow(u.int_in_range(1..=3i64)?));
    }
    let d = depth - 1;
    Ok(match u.int_in_range(0u8..=3)? {
        0 => trigonometric(u, d)? + trigonometric(u, d)?,
        1 => trigonometric(u, d)? * trigonometric(u, d)?,
        2 => trigonometric(u, d)? / (Expr::from(2) + f(Symbol::COS, angle(u)?)),
        _ => trigonometric(u, d)? + small(u)?,
    })
}

/// Inverse trigonometric functions, `asin` and `acos` within `(−1, 1)`.
fn inverse(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    let leaf = |u: &mut Unstructured<'_>| -> Result<Expr> {
        Ok(match u.int_in_range(0u8..=2)? {
            0 => f(Symbol::ARC_SIN, unit(u)?),
            1 => f(Symbol::ARC_COS, unit(u)?),
            _ => f(Symbol::ARC_TAN, real(u, 1)?),
        })
    };
    if depth == 0 {
        return leaf(u);
    }
    let d = depth - 1;
    Ok(match u.int_in_range(0u8..=2)? {
        0 => inverse(u, d)? + small(u)? * inverse(u, d)?,
        1 => inverse(u, d)? * real(u, 0)?,
        _ => leaf(u)?,
    })
}

/// Expressions of `x + i y` and quotients with square roots.
fn complex(u: &mut Unstructured<'_>, depth: u32) -> Result<Expr> {
    let z = |u: &mut Unstructured<'_>| -> Result<Expr> {
        Ok(small(u)? * x() + Expr::from(Symbol::I) * small(u)? * y() + small(u)?)
    };
    if depth == 0 {
        return z(u);
    }
    let d = depth - 1;
    Ok(match u.int_in_range(0u8..=11)? {
        0 => complex(u, d)? + complex(u, d)?,
        1 => complex(u, d)? * complex(u, d)?,
        2 => complex(u, d)?.pow(u.int_in_range(-3..=3i64)?),
        3 => f(Symbol::EXP, z(u)?),
        4 => f(Symbol::LOG, z(u)?),
        5 => f(
            *u.choose(&[
                Symbol::SIN,
                Symbol::COS,
                Symbol::TAN,
                Symbol::SINH,
                Symbol::COSH,
                Symbol::TANH,
            ])?,
            z(u)?,
        ),
        6 => f(Symbol::ABS, z(u)?),
        7 => f(Symbol::CONJUGATE, z(u)?),
        8 => z(u)?.pow(small(u)? + Expr::from(Symbol::I) * small(u)?),
        9 => {
            // A square root in a sum, over or under.
            let root = positive(u, 0)?.pow(Expr::from(q(1, 2)));
            let sum = small(u)? * root + small(u)? + small(u)? * x();
            if u.arbitrary()? {
                real(u, 1)? / sum
            } else {
                sum
            }
        }
        10 => real(u, 1)? / z(u)?,
        _ => z(u)?,
    })
}

/// `sign` and `θ` of linear forms.
fn steps(u: &mut Unstructured<'_>) -> Result<Expr> {
    let mut terms = Vec::new();
    for _ in 0..u.int_in_range(1..=3)? {
        let head = *u.choose(&[Symbol::SIGN, Symbol::HEAVISIDE_THETA])?;
        let linear = small(u)? * x() + small(u)? * y() + small(u)?;
        terms.push(small(u)? * f(head, linear));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// The value at the point, as balls of its parts; `None` for no finite
/// number.
fn value(ctx: &Context, e: &Expr, point: &[(Expr, Expr)]) -> Option<(Ball, Ball)> {
    let v = ctx
        .eval(&Expr::call(
            Symbol::N,
            [e.substitute(point), Expr::from(30)],
        ))
        .ok()?;
    let (re, im) = match v.view() {
        View::Ball(b) => (b.clone(), Ball::ZERO),
        View::ComplexBall(z) => (z.re().clone(), z.im().clone()),
        View::Complex(z) => (
            Ball::from_rational(z.re(), 128),
            Ball::from_rational(z.im(), 128),
        ),
        _ => (Ball::from_rational(&v.to_rational()?, 128), Ball::ZERO),
    };
    (re.is_finite() && im.is_finite()).then_some((re, im))
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let one = |s: Symbol| -> Box<dyn Fn(Expr) -> Expr> { Box::new(move |e| f(s, e)) };
    let combine = |what: Expr| -> Box<dyn Fn(Expr) -> Expr> {
        Box::new(move |e| Expr::call(Symbol::COMBINE, [e, what.clone()]))
    };
    let (input, command) = match u.int_in_range(0u8..=5)? {
        0 => {
            let e = exponential(u, 2)?;
            let c = match u.int_in_range(0u8..=5)? {
                0 => one(Symbol::HYPERBOLIC_TO_EXP),
                1 => one(Symbol::EXP_TO_POWER),
                2 => one(Symbol::POWER_TO_EXP),
                3 => one(Symbol::EXP_EXPAND),
                4 => one(Symbol::TRANSCENDENTAL_SIMPLIFY),
                _ => combine(Expr::from(Symbol::EXP)),
            };
            (e, c)
        }
        1 => {
            let e = logarithmic(u)?;
            let c = match u.int_in_range(0u8..=2)? {
                0 => one(Symbol::LOG_EXPAND),
                1 => one(Symbol::LOG_COMBINE),
                _ => combine(Expr::from(Symbol::LOG)),
            };
            (e, c)
        }
        2 => {
            let e = trigonometric(u, 2)?;
            let c = match u.int_in_range(0u8..=12)? {
                0 => one(Symbol::HALF_TAN),
                1 => one(Symbol::TAN_TO_SIN_COS),
                2 => one(Symbol::SIN_TO_COS_TAN),
                3 => one(Symbol::COS_TO_SIN_TAN),
                4 => one(Symbol::TAN_TO_SIN_COS_2),
                5 => one(Symbol::TAN_TO_COS_SIN_2),
                6 => one(Symbol::TRIG_TO_COS),
                7 => one(Symbol::TRIG_TO_SIN),
                8 => one(Symbol::TRIG_TO_TAN),
                9 => one(Symbol::TRIG_COLLECT),
                10 => one(Symbol::TRANSCENDENTAL_SIMPLIFY),
                11 => one(Symbol::HALF_TAN_HYPERBOLIC_TO_EXP),
                _ => combine(Expr::from(Symbol::new("trig"))),
            };
            (e, c)
        }
        3 => {
            let e = inverse(u, 2)?;
            let c = one(*u.choose(&[
                Symbol::ARC_COS_TO_ARC_SIN,
                Symbol::ARC_SIN_TO_ARC_COS,
                Symbol::ARC_COS_TO_ARC_TAN,
                Symbol::ARC_SIN_TO_ARC_TAN,
                Symbol::ARC_TAN_TO_ARC_COS,
                Symbol::ARC_TAN_TO_ARC_SIN,
                Symbol::ARC_TRIG_TO_LOG,
            ])?);
            (e, c)
        }
        4 => {
            let e = complex(u, 2)?;
            let c = one(*u.choose(&[
                Symbol::COMPLEX_EXPAND,
                Symbol::MULTIPLY_CONJUGATE,
                Symbol::MULTIPLY_COMPLEX_CONJUGATE,
            ])?);
            (e, c)
        }
        _ => {
            let e = steps(u)?;
            let c = one(*u.choose(&[Symbol::SIGN_TO_HEAVISIDE, Symbol::HEAVISIDE_TO_SIGN])?);
            (e, c)
        }
    };
    let Ok(input) = ctx.eval(&input) else {
        return Ok(());
    };
    let call = command(input.clone());
    let output = ctx.eval(&call).unwrap_or_else(|e| panic!("{call:?}: {e}"));
    // x and y in (0, 1).
    let mut coordinate = || -> Result<Expr> {
        let d = u.int_in_range(2..=16)?;
        Ok(Expr::from(q(u.int_in_range(1..=d - 1)?, d)))
    };
    let point = [(x(), coordinate()?), (y(), coordinate()?)];
    let (Some((a, b)), Some((c, d))) = (value(&ctx, &input, &point), value(&ctx, &output, &point))
    else {
        return Ok(());
    };
    assert!(
        a.overlaps(&c) && b.overlaps(&d),
        "{call:?} is {output:?}: at {point:?}, {a:?} + {b:?} i against {c:?} + {d:?} i"
    );
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
