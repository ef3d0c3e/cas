//! Z-transforms and recurrences (CAS-DOC-24 §5, CAS-311), one check per
//! input:
//! - a sequence built from the table, rational at the naturals (constants,
//!   `aⁿ`, `n aⁿ`, binomials, powers of `n`, Kronecker deltas), comes back
//!   from `InverseZTransform[ZTransform[f]]` at `n = 0 … 8`, exactly;
//! - a random linear recurrence of order 1 or 2 with rational
//!   coefficients, initial values and forcing: `RSolve`'s solution is the
//!   sequence the recurrence makes, computed here, at `n = 0 … 8`, to
//!   `10⁻²⁰` (its closed form may hold radicals, or cosines of complex
//!   roots).
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

/// `e` at `n = k`, evaluated.
fn at(ctx: &Context, e: &Expr, k: i64) -> Expr {
    let e = e.substitute(&[(sym("n"), Expr::from(k))]);
    ctx.eval(&e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// Whether `a` is within `10⁻²⁰` of the rational `b`.
fn near(ctx: &Context, a: &Expr, b: &Rational) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [a.clone() - Expr::from(b.clone()), Expr::from(30)],
        )],
    );
    let claim = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-20)]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// A term of the table, rational at the naturals.
fn term(u: &mut Unstructured<'_>) -> Result<Expr> {
    let n = sym("n");
    let c = Expr::from(nonzero(u)?);
    Ok(match u.int_in_range(0..=5u8)? {
        0 => c,
        1 => c * Expr::from(nonzero(u)?).pow(n),
        2 => c * n.clone() * Expr::from(nonzero(u)?).pow(n),
        3 => c * call("Binomial", vec![n, Expr::from(u.int_in_range(0..=3i64)?)]),
        4 => c * n.pow(u.int_in_range(1..=3i64)?),
        _ => c * call("KroneckerDelta", vec![n - u.int_in_range(0..=4i64)?]),
    })
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let n = sym("n");
    let z = sym("z");
    if u.ratio(1, 2)? {
        // A round trip through the transform.
        let count = u.int_in_range(1..=3usize)?;
        let terms = (0..count).map(|_| term(u)).collect::<Result<Vec<_>>>()?;
        let f = Expr::call(Symbol::PLUS, terms);
        let image = ctx
            .eval(&call("ZTransform", vec![f.clone(), n.clone(), z.clone()]))
            .expect("a transform");
        if image.contains(&sym("ZTransform")) {
            return Ok(());
        }
        let back = ctx
            .eval(&call("InverseZTransform", vec![image.clone(), z, n]))
            .expect("an inverse");
        if back.contains(&sym("InverseZTransform")) {
            return Ok(());
        }
        for k in 0..=8 {
            let (want, got) = (at(&ctx, &f, k), at(&ctx, &back, k));
            let want = want.to_rational().expect("rational at the naturals");
            assert!(
                near(&ctx, &got, &want),
                "{f:?} at {k}: {image:?} gave {got:?}"
            );
        }
        return Ok(());
    }
    // A recurrence Σ_{i≤k} a_i u(n + i) = g(n) with u(0 … k−1) given.
    let order = u.int_in_range(1..=2usize)?;
    let a: Vec<Rational> = (0..=order)
        .map(|i| if i == order { nonzero(u) } else { small(u) })
        .collect::<Result<_>>()?;
    let initial: Vec<Rational> = (0..order).map(|_| small(u)).collect::<Result<_>>()?;
    let (gc, gb) = (small(u)?, nonzero(u)?);
    let kind = u.int_in_range(0..=2u8)?;
    let forcing = |k: i64| -> Rational {
        match kind {
            0 => Rational::ZERO,
            1 => gc.clone(),
            _ => &gc * &(0..k).fold(Rational::ONE, |p, _| &p * &gb),
        }
    };
    let g = match kind {
        0 => Expr::from(0),
        1 => Expr::from(gc.clone()),
        _ => Expr::from(gc.clone()) * Expr::from(gb.clone()).pow(n.clone()),
    };
    let unknown = |shift: usize| call("u", vec![n.clone() + i64::try_from(shift).expect("small")]);
    let left = Expr::call(
        Symbol::PLUS,
        a.iter()
            .enumerate()
            .map(|(i, c)| Expr::from(c.clone()) * unknown(i))
            .collect::<Vec<_>>(),
    );
    let mut args = vec![
        Expr::call(Symbol::EQUATION, [left, g]),
        call("u", vec![n.clone()]),
    ];
    for (i, v) in initial.iter().enumerate() {
        let point = Expr::from(i64::try_from(i).expect("small"));
        args.push(Expr::call(
            Symbol::EQUATION,
            [call("u", vec![point]), Expr::from(v.clone())],
        ));
    }
    let input = call("RSolve", args);
    let found = ctx.eval(&input).expect("an answer");
    let [solution] = found.args() else {
        // Left unsolved, or no solution: nothing to check.
        return Ok(());
    };
    if found.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(());
    }
    // The sequence, by the recurrence.
    let mut values = initial.clone();
    for k in 0..=8i64 {
        let shift = usize::try_from(k).expect("small");
        if values.len() < shift + order + 1 {
            let mut sum = forcing(k);
            for (i, c) in a.iter().enumerate().take(order) {
                sum = &sum - &(c * &values[shift + i]);
            }
            values.push(&sum / &a[order]);
        }
        let got = at(&ctx, solution, k);
        assert!(
            near(&ctx, &got, &values[shift]),
            "{input:?}: {solution:?} at {k} is {got:?}, the recurrence makes {:?}",
            values[shift]
        );
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
