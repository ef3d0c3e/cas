//! Interpolation and orthogonal polynomials (CAS-DOC-27 §6, CAS-332), one
//! check per input:
//! - random rational points, distinct and increasing: the
//!   `InterpolatingPolynomial` takes the values at them; Thiele's
//!   interpolant does where it is finite, or stays; the pieces of a spline
//!   of degree 1 to 4 meet their points, their derivatives of orders
//!   `1 … d − 1` agree at the inner points, and giac's end conditions hold;
//! - an orthogonal polynomial of a random degree up to 12 (Laguerre's with
//!   a random rational parameter) satisfies its differential equation, and
//!   its value at a random rational is the polynomial's.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn small(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(-9..=9)?, u.int_in_range(1..=4)?))
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

fn list(items: &[Rational]) -> Expr {
    Expr::call(
        Symbol::LIST,
        items.iter().cloned().map(Expr::from).collect::<Vec<_>>(),
    )
}

/// `f` at `x = v`, evaluated.
fn at(ctx: &Context, f: &Expr, v: &Rational) -> Expr {
    eval(ctx, &f.substitute(&[(x(), Expr::from(v.clone()))]))
}

/// The `j`-th derivative of `f` in `x`.
fn derivative(ctx: &Context, f: &Expr, j: usize) -> Expr {
    let mut g = f.clone();
    for _ in 0..j {
        g = eval(ctx, &Expr::call(Symbol::D, [g, x()]));
    }
    g
}

fn points(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let count = u.int_in_range(2..=6usize)?;
    let mut xs: Vec<Rational> = Vec::with_capacity(count);
    for _ in 0..count {
        let v = small(u)?;
        if !xs.contains(&v) {
            xs.push(v);
        }
    }
    if xs.len() < 2 {
        return Ok(());
    }
    xs.sort();
    let ys: Vec<Rational> = xs.iter().map(|_| small(u)).collect::<Result<_>>()?;
    let input = call("InterpolatingPolynomial", vec![list(&xs), list(&ys)]);
    let p = eval(ctx, &input);
    for (a, b) in xs.iter().zip(&ys) {
        assert_eq!(
            at(ctx, &p, a),
            Expr::from(b.clone()),
            "{input:?} = {p:?} at {a}"
        );
    }
    let input = call("ThieleInterpolation", vec![list(&xs), list(&ys), x()]);
    let r = eval(ctx, &input);
    if !r.contains(&Expr::from(Symbol::new("ThieleInterpolation"))) {
        for (a, b) in xs.iter().zip(&ys) {
            let value = at(ctx, &r, a);
            if value.to_rational().is_some() {
                assert_eq!(value, Expr::from(b.clone()), "{input:?} = {r:?} at {a}");
            }
        }
    }
    let d = u.int_in_range(1..=4usize)?;
    let input = call(
        "Spline",
        vec![
            list(&xs),
            list(&ys),
            x(),
            Expr::from(i64::try_from(d).expect("small")),
        ],
    );
    let spline = eval(ctx, &input);
    let pieces = spline.args();
    assert_eq!(pieces.len(), xs.len() - 1, "{input:?}: {spline:?}");
    let n = pieces.len();
    for (i, s) in pieces.iter().enumerate() {
        assert_eq!(
            at(ctx, s, &xs[i]),
            Expr::from(ys[i].clone()),
            "{input:?}: piece {i} at its start"
        );
        assert_eq!(
            at(ctx, s, &xs[i + 1]),
            Expr::from(ys[i + 1].clone()),
            "{input:?}: piece {i} at its end"
        );
    }
    for j in 1..d {
        for i in 0..n.saturating_sub(1) {
            let left = at(ctx, &derivative(ctx, &pieces[i], j), &xs[i + 1]);
            let right = at(ctx, &derivative(ctx, &pieces[i + 1], j), &xs[i + 1]);
            assert_eq!(left, right, "{input:?}: derivative {j} at {}", xs[i + 1]);
        }
    }
    let zero = Expr::from(0);
    for j in d.div_ceil(2)..d {
        assert_eq!(
            at(ctx, &derivative(ctx, &pieces[0], j), &xs[0]),
            zero,
            "{input:?}: derivative {j} at the start"
        );
    }
    for j in (d / 2 + 1)..d {
        assert_eq!(
            at(ctx, &derivative(ctx, &pieces[n - 1], j), &xs[n]),
            zero,
            "{input:?}: derivative {j} at the end"
        );
    }
    Ok(())
}

fn orthogonal(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(0..=12i64)?;
    let a = small(u)?;
    let family = u.int_in_range(0..=4u8)?;
    let (head, args) = match family {
        0 => ("LegendreP", vec![Expr::from(n), x()]),
        1 => ("HermiteH", vec![Expr::from(n), x()]),
        2 => ("LaguerreL", vec![Expr::from(n), x(), Expr::from(a.clone())]),
        3 => ("ChebyshevT", vec![Expr::from(n), x()]),
        _ => ("ChebyshevU", vec![Expr::from(n), x()]),
    };
    let input = call(head, args.clone());
    let p = eval(ctx, &input);
    let (d1, d2) = (derivative(ctx, &p, 1), derivative(ctx, &p, 2));
    let xx = x();
    let nn = Expr::from(n);
    let one = || Expr::from(1);
    let residual = match family {
        0 => {
            (one() - xx.pow(2)) * d2 - Expr::from(2) * xx.clone() * d1
                + nn.clone() * (nn.clone() + one()) * p.clone()
        }
        1 => d2 - Expr::from(2) * xx.clone() * d1 + Expr::from(2) * nn.clone() * p.clone(),
        2 => xx.clone() * d2 + (Expr::from(a) + one() - xx.clone()) * d1 + nn.clone() * p.clone(),
        3 => (one() - xx.pow(2)) * d2 - xx.clone() * d1 + nn.pow(2) * p.clone(),
        _ => {
            (one() - xx.pow(2)) * d2 - Expr::from(3) * xx.clone() * d1
                + nn.clone() * (nn.clone() + Expr::from(2)) * p.clone()
        }
    };
    let residual = eval(ctx, &Expr::call(Symbol::EXPAND, [residual]));
    assert_eq!(residual, Expr::from(0), "{input:?} = {p:?}");
    let v = small(u)?;
    let mut at_v = args;
    at_v[1] = Expr::from(v.clone());
    let value = eval(ctx, &call(head, at_v));
    assert_eq!(value, at(ctx, &p, &v), "{input:?} at {v}");
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    if u.ratio(1, 2)? {
        points(&ctx, u)
    } else {
        orthogonal(&ctx, u)
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
