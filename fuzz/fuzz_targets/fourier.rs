//! Residues and Fourier coefficients (CAS-DOC-25 §8, CAS-320), one check
//! per input:
//! - a rational function `p/q` with chosen roots, rational or Gaussian
//!   rational, of multiplicity 1 to 3, `q` expanded or not and
//!   `deg p < deg q`: its residues at the roots add up to `lc(p)/lc(q)`
//!   when `deg p = deg q − 1` and to 0 otherwise, and at a simple root `r`
//!   the residue is `p(r)/q′(r)`;
//! - a trigonometric polynomial `c + Σ_{k≤3} (a_k cos kωx + b_k sin kωx)`
//!   with rational coefficients, of period `T = 2π/ω`, `2π` or a rational,
//!   from a start `0`, `−T/2` or a rational: its coefficients come back for
//!   `n = 0 … 4`;
//! - a polynomial plus `Abs`, `Sign` or `HeavisideTheta` of `x − c` over
//!   `[−π, π]`: the coefficient for a symbolic `n`, at `n = 1 … 4`, is the
//!   coefficient for that `n`.
//!
//! Values are compared to `10⁻²⁰`, real and imaginary parts apart: `sin`
//! of a rational multiple of `π` stays as it is, so that equal values may
//! be written apart.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
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

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// Whether `a` and `b` agree to `10⁻²⁰`, real and imaginary parts apart.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// The residues of a rational function with chosen roots.
fn residues(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let z = sym("z");
    let mut roots: Vec<(Expr, u32)> = Vec::new();
    for _ in 0..u.int_in_range(1..=3u8)? {
        let im = if u.ratio(1, 2)? {
            small(u)?
        } else {
            Rational::ZERO
        };
        let r = eval(
            &ctx,
            &(Expr::from(small(u)?) + Expr::from(im) * Expr::from(Symbol::I)),
        );
        if roots.iter().all(|(s, _)| *s != r) {
            roots.push((r, u.int_in_range(1..=3)?));
        }
    }
    let degree: u32 = roots.iter().map(|(_, m)| m).sum();
    let lead = nonzero(u)?;
    let mut factors = vec![Expr::from(lead.clone())];
    for (r, m) in &roots {
        factors.push((z.clone() - r.clone()).pow(i64::from(*m)));
    }
    let mut den = Expr::call(Symbol::TIMES, factors);
    if u.ratio(1, 2)? {
        den = eval(&ctx, &Expr::call(Symbol::EXPAND, [den]));
    }
    let top = u.int_in_range(0..=degree - 1)?;
    let mut coefficients: Vec<Rational> = (0..top).map(|_| small(u)).collect::<Result<_>>()?;
    coefficients.push(nonzero(u)?);
    let num = Expr::call(
        Symbol::PLUS,
        coefficients
            .iter()
            .enumerate()
            .map(|(k, c)| Expr::from(c.clone()) * z.pow(i64::try_from(k).expect("small")))
            .collect::<Vec<_>>(),
    );
    let f = num.clone() / den.clone();
    let derivative = eval(&ctx, &Expr::call(Symbol::D, [den, z.clone()]));
    let mut total = Vec::new();
    for (r, m) in &roots {
        let residue = eval(
            &ctx,
            &call("Residue", vec![f.clone(), z.clone(), r.clone()]),
        );
        assert!(
            !residue.contains(&sym("Residue")),
            "{f:?} at {r:?}: {residue:?}"
        );
        if *m == 1 {
            let at = |g: &Expr| g.substitute(&[(z.clone(), r.clone())]);
            let simple = eval(&ctx, &(at(&num) / at(&derivative)));
            assert!(
                near(&ctx, &residue, &simple),
                "{f:?} at the simple root {r:?}: {residue:?}, p(r)/q'(r) is {simple:?}"
            );
        }
        total.push(residue);
    }
    let expected = if top + 1 == degree {
        Expr::from(&coefficients[coefficients.len() - 1] / &lead)
    } else {
        Expr::from(0)
    };
    let sum = Expr::call(Symbol::PLUS, total);
    assert!(
        near(&ctx, &sum, &expected),
        "{f:?}: the residues add up to {sum:?}, not {expected:?}"
    );
    Ok(())
}

/// The coefficients of a trigonometric polynomial.
fn trigonometric(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let x = sym("x");
    let pi = Expr::from(Symbol::PI);
    let period = if u.ratio(1, 2)? {
        Expr::from(2) * pi.clone()
    } else {
        Expr::from(q(u.int_in_range(1..=8)?, u.int_in_range(1..=3)?))
    };
    let start = match u.int_in_range(0..=2u8)? {
        0 => Expr::from(0),
        1 => -(period.clone() / Expr::from(2)),
        _ => Expr::from(small(u)?),
    };
    let count = u.int_in_range(1..=3usize)?;
    let c = small(u)?;
    let mut a = vec![Rational::ZERO];
    let mut b = vec![Rational::ZERO];
    let mut terms = vec![Expr::from(c.clone())];
    for k in 1..=count {
        let (ak, bk) = (small(u)?, small(u)?);
        let w =
            Expr::from(2) * pi.clone() * Expr::from(i64::try_from(k).expect("small")) * x.clone()
                / period.clone();
        terms.push(Expr::from(ak.clone()) * Expr::call(Symbol::COS, [w.clone()]));
        terms.push(Expr::from(bk.clone()) * Expr::call(Symbol::SIN, [w]));
        a.push(ak);
        b.push(bk);
    }
    let f = eval(&ctx, &Expr::call(Symbol::PLUS, terms));
    for n in 0..=4usize {
        let (an, bn) = if n == 0 {
            (c.clone(), Rational::ZERO)
        } else {
            (
                a.get(n).cloned().unwrap_or(Rational::ZERO),
                b.get(n).cloned().unwrap_or(Rational::ZERO),
            )
        };
        let cn = if n == 0 {
            Expr::from(c.clone())
        } else {
            (Expr::from(an.clone()) - Expr::from(Symbol::I) * Expr::from(bn.clone()))
                / Expr::from(2)
        };
        let index = Expr::from(i64::try_from(n).expect("small"));
        for (head, expected) in [
            ("FourierCosCoefficient", Expr::from(an)),
            ("FourierSinCoefficient", Expr::from(bn)),
            ("FourierCoefficient", cn),
        ] {
            let input = call(
                head,
                vec![
                    f.clone(),
                    x.clone(),
                    period.clone(),
                    index.clone(),
                    start.clone(),
                ],
            );
            let found = eval(&ctx, &input);
            assert!(
                !found.contains(&sym(head)) && near(&ctx, &found, &expected),
                "{input:?} is {found:?}, not {expected:?}"
            );
        }
    }
    Ok(())
}

/// The coefficients of a polynomial and a step for a symbolic `n` and
/// for each `n`.
fn piecewise(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let x = sym("x");
    let mut terms = Vec::new();
    for k in 0..=u.int_in_range(0..=2i64)? {
        terms.push(Expr::from(small(u)?) * x.pow(k));
    }
    let shift = x.clone() - Expr::from(q(u.int_in_range(-8..=8)?, 3));
    let step = match u.int_in_range(0..=2u8)? {
        0 => "Abs",
        1 => "Sign",
        _ => "HeavisideTheta",
    };
    terms.push(Expr::from(nonzero(u)?) * call(step, vec![shift]));
    let f = eval(&ctx, &Expr::call(Symbol::PLUS, terms));
    let head = match u.int_in_range(0..=2u8)? {
        0 => "FourierCosCoefficient",
        1 => "FourierSinCoefficient",
        _ => "FourierCoefficient",
    };
    let pi = Expr::from(Symbol::PI);
    let coefficient = |n: Expr| {
        call(
            head,
            vec![
                f.clone(),
                x.clone(),
                Expr::from(2) * pi.clone(),
                n,
                -pi.clone(),
            ],
        )
    };
    let n = sym("n");
    let generic = eval(&ctx, &coefficient(n.clone()));
    if generic.contains(&sym(head)) {
        return Ok(());
    }
    for k in 1..=4i64 {
        let given = eval(&ctx, &coefficient(Expr::from(k)));
        let at = eval(&ctx, &generic.substitute(&[(n.clone(), Expr::from(k))]));
        assert!(
            near(&ctx, &at, &given),
            "{f:?}: {head} for n is {generic:?}, {at:?} at {k}, but {given:?} for {k}"
        );
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    match u.int_in_range(0..=2u8)? {
        0 => residues(u),
        1 => trigonometric(u),
        _ => piecewise(u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
