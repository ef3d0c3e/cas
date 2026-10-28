//! Vector calculus and curves (CAS-DOC-37 §6, CAS-413), one check per
//! input, on random fields of polynomials, exponentials and sines and
//! cosines of linear forms, by values at a random point, the balls of `N`:
//! - `curl grad f = 0`, `div curl F = 0`, `Laplacian f = div grad f`, and
//!   the Hessian symmetric with trace the Laplacian;
//! - `grad Potential[grad f] = grad f`, and `curl VectorPotential[curl G]
//!   = curl G`, where they do not stay;
//! - `ImplicitD[y = g(x), y, x]` is `g′`, and to order 2 `g″`; for `y² =
//!   h(x)`, `h > 0` at the point, `ImplicitD` at `y = √h` is `D[√h, x]`;
//! - the circle `{a + r cos(k t), b + r sin(k t)}` has curvature `sign(k)/r`;
//!   Frenet's `T` and `N` are orthonormal on random curves; a segment's arc
//!   length is its length.
//!
//! A value that is no finite number, at a pole, checks nothing.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::ball::Ball;
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn small(u: &mut Unstructured<'_>) -> Result<Expr> {
    Ok(Expr::from(q(
        u.int_in_range(-4..=4)?,
        u.int_in_range(1..=3)?,
    )))
}

/// A linear form in `vars` with small coefficients.
fn linear(u: &mut Unstructured<'_>, vars: &[Expr]) -> Result<Expr> {
    let mut terms = vec![small(u)?];
    for v in vars {
        terms.push(small(u)? * v.clone());
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// A field: a polynomial of degree at most 2 times, sometimes, an
/// exponential, sine or cosine of a linear form.
fn field(u: &mut Unstructured<'_>, vars: &[Expr]) -> Result<Expr> {
    let mut terms = vec![small(u)?];
    for _ in 0..u.int_in_range(1..=3)? {
        let i = u.choose_index(vars.len())?;
        let j = u.choose_index(vars.len())?;
        let monomial = match u.int_in_range(0u8..=2)? {
            0 => vars[i].clone(),
            1 => vars[i].clone() * vars[j].clone(),
            _ => vars[i].clone().pow(2),
        };
        terms.push(small(u)? * monomial);
    }
    let polynomial = Expr::call(Symbol::PLUS, terms);
    Ok(match u.int_in_range(0u8..=3)? {
        0 => polynomial * Expr::call(Symbol::EXP, [linear(u, vars)?]),
        1 => polynomial * Expr::call(Symbol::SIN, [linear(u, vars)?]),
        2 => polynomial + Expr::call(Symbol::COS, [linear(u, vars)?]),
        _ => polynomial,
    })
}

/// The real ball of `e` at the point, by `N`; `None` for no finite real.
fn value(ctx: &Context, e: &Expr, point: &[(Expr, Expr)]) -> Option<Ball> {
    let v = ctx
        .eval(&Expr::call(
            Symbol::N,
            [e.substitute(point), Expr::from(30)],
        ))
        .ok()?;
    let b = match v.view() {
        View::Ball(b) => b.clone(),
        _ => Ball::from_rational(&v.to_rational()?, 128),
    };
    b.is_finite().then_some(b)
}

/// Every entry of `a − b`, lists or not, holds 0 at the point.
fn agree(ctx: &Context, what: &str, a: &Expr, b: &Expr, point: &[(Expr, Expr)]) {
    let (xs, ys) = match (a.head(), b.head()) {
        (Some(h), Some(k)) if *h == Expr::from(Symbol::LIST) && h == k => {
            (a.args().to_vec(), b.args().to_vec())
        }
        _ => (vec![a.clone()], vec![b.clone()]),
    };
    assert_eq!(xs.len(), ys.len(), "{what}: {a:?} against {b:?}");
    for (x, y) in xs.iter().zip(&ys) {
        let (Some(p), Some(r)) = (value(ctx, x, point), value(ctx, y, point)) else {
            continue;
        };
        assert!(
            p.overlaps(&r),
            "{what}: {x:?} is {p:?} and {y:?} is {r:?} at {point:?}"
        );
    }
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(head: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(head, args)
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let (x, y, z, t) = (sym("x"), sym("y"), sym("z"), sym("t"));
    let vars = vec![x.clone(), y.clone(), z.clone()];
    let mut coordinate = || -> Result<Expr> { small(u) };
    let point = vec![
        (x.clone(), coordinate()?),
        (y.clone(), coordinate()?),
        (z.clone(), coordinate()?),
        (t.clone(), coordinate()?),
    ];
    let v = list(vars.clone());
    match u.int_in_range(0u8..=8)? {
        0 => {
            let f = field(u, &vars)?;
            let grad = call(Symbol::GRAD, vec![f, v.clone()]);
            let curl = eval(&ctx, &call(Symbol::CURL, vec![grad, v]));
            agree(
                &ctx,
                "curl grad",
                &curl,
                &list(vec![Expr::from(0); 3]),
                &point,
            );
        }
        1 => {
            let fs = (0..3)
                .map(|_| field(u, &vars))
                .collect::<Result<Vec<_>>>()?;
            let curl = call(Symbol::CURL, vec![list(fs), v.clone()]);
            let div = eval(&ctx, &call(Symbol::DIVERGENCE, vec![curl, v]));
            agree(&ctx, "div curl", &div, &Expr::from(0), &point);
        }
        2 => {
            let f = field(u, &vars)?;
            let laplacian = eval(&ctx, &call(Symbol::LAPLACIAN, vec![f.clone(), v.clone()]));
            let grad = call(Symbol::GRAD, vec![f.clone(), v.clone()]);
            let div = eval(&ctx, &call(Symbol::DIVERGENCE, vec![grad, v.clone()]));
            agree(&ctx, "Laplacian", &laplacian, &div, &point);
            let h = eval(&ctx, &call(Symbol::HESSIAN, vec![f, v]));
            let rows = h.args();
            let transposed: Vec<Expr> = (0..3)
                .map(|j| list(rows.iter().map(|r| r.args()[j].clone()).collect()))
                .collect();
            agree(&ctx, "Hessian symmetric", &h, &list(transposed), &point);
            let trace = Expr::call(Symbol::PLUS, (0..3).map(|k| rows[k].args()[k].clone()));
            agree(&ctx, "trace", &trace, &laplacian, &point);
        }
        3 => {
            let f = field(u, &vars)?;
            let grad = eval(&ctx, &call(Symbol::GRAD, vec![f, v.clone()]));
            let potential = eval(
                &ctx,
                &call(Symbol::POTENTIAL, vec![grad.clone(), v.clone()]),
            );
            if potential.head() != Some(&Expr::from(Symbol::POTENTIAL)) {
                let back = eval(&ctx, &call(Symbol::GRAD, vec![potential, v]));
                agree(&ctx, "grad Potential", &back, &grad, &point);
            }
        }
        4 => {
            let gs = (0..3)
                .map(|_| field(u, &vars))
                .collect::<Result<Vec<_>>>()?;
            let curl = eval(&ctx, &call(Symbol::CURL, vec![list(gs), v.clone()]));
            let a = eval(
                &ctx,
                &call(Symbol::VECTOR_POTENTIAL, vec![curl.clone(), v.clone()]),
            );
            if a.head() != Some(&Expr::from(Symbol::VECTOR_POTENTIAL)) {
                let back = eval(&ctx, &call(Symbol::CURL, vec![a, v]));
                agree(&ctx, "curl VectorPotential", &back, &curl, &point);
            }
        }
        5 => {
            // y = g(x), and y² = h(x) at y = √h.
            let g = field(u, &[x.clone()])?;
            let equation = call(Symbol::EQUATION, vec![y.clone(), g.clone()]);
            let first = eval(
                &ctx,
                &call(
                    Symbol::IMPLICIT_D,
                    vec![equation.clone(), y.clone(), x.clone()],
                ),
            );
            let g1 = eval(&ctx, &call(Symbol::D, vec![g.clone(), x.clone()]));
            agree(&ctx, "ImplicitD of y = g", &first, &g1, &point);
            let second = eval(
                &ctx,
                &call(
                    Symbol::IMPLICIT_D,
                    vec![equation, y.clone(), x.clone(), x.clone()],
                ),
            );
            let g2 = eval(
                &ctx,
                &call(Symbol::D, vec![g.clone(), x.clone(), x.clone()]),
            );
            agree(&ctx, "ImplicitD of order 2", &second, &g2, &point);
            let h = g.clone().pow(2) + Expr::from(1);
            let root = h.clone().pow(Expr::from(q(1, 2)));
            let equation = call(Symbol::EQUATION, vec![y.clone().pow(2), h]);
            let found = eval(
                &ctx,
                &call(Symbol::IMPLICIT_D, vec![equation, y.clone(), x.clone()]),
            );
            let at_root = found.substitute(&[(y.clone(), root.clone())]);
            let expected = eval(&ctx, &call(Symbol::D, vec![root, x.clone()]));
            agree(&ctx, "ImplicitD of y² = h", &at_root, &expected, &point);
        }
        6 => {
            // A circle of radius r, turning k times as fast.
            let r = Expr::from(q(u.int_in_range(1..=9)?, u.int_in_range(1..=4)?));
            let k = u.int_in_range(1..=3i64)? * if u.arbitrary()? { -1 } else { 1 };
            let kt = Expr::from(k) * t.clone();
            let circle = list(vec![
                small(u)? + r.clone() * Expr::call(Symbol::COS, [kt.clone()]),
                small(u)? + r.clone() * Expr::call(Symbol::SIN, [kt]),
            ]);
            let kappa = eval(&ctx, &call(Symbol::CURVATURE, vec![circle, t.clone()]));
            let expected = Expr::from(k.signum()) / r;
            agree(&ctx, "a circle's curvature", &kappa, &expected, &point);
        }
        7 => {
            // Frenet's T and N are orthonormal.
            let n = u.int_in_range(2..=3usize)?;
            let curve = list(
                (0..n)
                    .map(|_| field(u, &[t.clone()]))
                    .collect::<Result<_>>()?,
            );
            let t0 = small(u)?;
            let frame = eval(
                &ctx,
                &call(Symbol::FRENET, vec![curve.clone(), t.clone(), t0]),
            );
            if frame.head() == Some(&Expr::from(Symbol::LIST)) {
                let items = frame.args();
                let (tangent, normal) = if n == 2 {
                    (&items[2], &items[3])
                } else {
                    (&items[3], &items[4])
                };
                let dot = |a: &Expr, b: &Expr| {
                    Expr::call(
                        Symbol::PLUS,
                        a.args()
                            .iter()
                            .zip(b.args())
                            .map(|(p, q)| p.clone() * q.clone()),
                    )
                };
                agree(&ctx, "T·N", &dot(tangent, normal), &Expr::from(0), &point);
                agree(&ctx, "T·T", &dot(tangent, tangent), &Expr::from(1), &point);
                agree(&ctx, "N·N", &dot(normal, normal), &Expr::from(1), &point);
            }
        }
        _ => {
            // A segment's length.
            let (a, b, c, d) = (small(u)?, small(u)?, small(u)?, small(u)?);
            let (t0, t1) = (small(u)?, small(u)?);
            let segment = list(vec![a + b.clone() * t.clone(), c + d.clone() * t.clone()]);
            let length = eval(
                &ctx,
                &call(
                    Symbol::ARC_LENGTH,
                    vec![segment, t.clone(), t0.clone(), t1.clone()],
                ),
            );
            let expected = (b.pow(2) + d.pow(2)).pow(Expr::from(q(1, 2))) * (t1 - t0);
            agree(&ctx, "a segment's length", &length, &expected, &point);
        }
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
