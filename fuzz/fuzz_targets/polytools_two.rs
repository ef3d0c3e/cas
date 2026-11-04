//! Polynomial tools II (CAS-DOC-48 §5, CAS-485), on random integer
//! polynomials up to degree 6, each command by its defining identity:
//! - `FromRoots` of `PolynomialRoots` the monic polynomial, within the
//!   balls; `HornerEval` the exact value; `TaylorShift` by `a` then `−a`
//!   the identity;
//! - `u a + v b = d` for `PolynomialExtendedGCD`, `d` dividing both;
//!   `PolynomialBezout`'s `u a + v b = c` with `deg u < deg b`;
//!   `PolynomialChineseRemainder`'s congruences for coprime moduli;
//!   `SeriesDivision` times `b` agreeing with `a` up to `xⁿ`;
//! - `ComplexFactor` multiplying back; `DigitPolynomial` at `b` giving `n`
//!   with balanced digits; `BerlekampMassey` a recurrence of the sequence,
//!   no longer than the one that made it.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(s: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(s, args)
}

fn x() -> Expr {
    Expr::symbol("x")
}

/// A random integer polynomial of degree `1..=max`, coefficients from the
/// constant up, its leading one not 0.
fn polynomial(u: &mut Unstructured<'_>, max: usize) -> Result<Vec<i64>> {
    let d = u.int_in_range(1..=max)?;
    let mut c: Vec<i64> = (0..d)
        .map(|_| u.int_in_range(-5..=5i64))
        .collect::<Result<_>>()?;
    let lead = u.int_in_range(1..=5i64)? * if u.arbitrary()? { 1 } else { -1 };
    c.push(lead);
    Ok(c)
}

fn expr(c: &[i64]) -> Expr {
    Expr::call(
        Symbol::PLUS,
        c.iter()
            .enumerate()
            .map(|(k, v)| Expr::from(*v) * x().pow(i64::try_from(k).expect("small")))
            .collect::<Vec<_>>(),
    )
}

/// Whether `e` expands to 0.
fn vanishes(ctx: &Context, e: Expr) -> bool {
    eval(ctx, &Expr::call(Symbol::EXPAND, [e])) == Expr::from(0)
}

fn divides(ctx: &Context, d: &Expr, a: &Expr) -> bool {
    let r = eval(
        ctx,
        &call(
            Symbol::POLYNOMIAL_REMAINDER,
            vec![a.clone(), d.clone(), x()],
        ),
    );
    vanishes(ctx, r)
}

fn degree(ctx: &Context, p: &Expr) -> i64 {
    eval(ctx, &call(Symbol::EXPONENT, vec![p.clone(), x()]))
        .as_i64()
        .unwrap_or(-1)
}

fn roots(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let c = polynomial(u, 5)?;
    let p = expr(&c);
    // FromRoots ∘ PolynomialRoots, within the balls.
    let found = eval(ctx, &call(Symbol::POLYNOMIAL_ROOTS, vec![p.clone()]));
    let back = eval(ctx, &call(Symbol::FROM_ROOTS, vec![found.clone()]));
    let lead = c[c.len() - 1];
    assert_eq!(
        back.args().len(),
        c.len(),
        "FromRoots of {found:?} for {p:?}"
    );
    for (k, b) in back.args().iter().enumerate() {
        let want = Expr::from(Rational::new(
            Integer::from(c[c.len() - 1 - k]),
            Integer::from(lead),
        ));
        // Real and imaginary parts apart: the modulus of a complex ball
        // about 0 may not be decided.
        let gap = b.clone() - want;
        for part in [Symbol::RE, Symbol::IM] {
            let size = Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]);
            let small = Expr::call(Symbol::LESS, [size, Expr::from(10).pow(-8)]);
            assert_eq!(
                eval(ctx, &small),
                Expr::from(Symbol::TRUE),
                "the roots {found:?} of {p:?} give {back:?}"
            );
        }
    }
    // Horner's value, and Taylor shifts back.
    let a = Rational::new(
        Integer::from(u.int_in_range(-6..=6i64)?),
        Integer::from(u.int_in_range(1..=3i64)?),
    );
    let value = c
        .iter()
        .rev()
        .fold(Rational::from(Integer::ZERO), |acc, v| {
            &(&acc * &a) + &Rational::from(Integer::from(*v))
        });
    let found = eval(
        ctx,
        &call(Symbol::HORNER_EVAL, vec![p.clone(), Expr::from(a.clone())]),
    );
    assert_eq!(found, Expr::from(value), "HornerEval[{p:?}, {a}]");
    let list = Expr::call(
        Symbol::LIST,
        c.iter().rev().map(|v| Expr::from(*v)).collect::<Vec<_>>(),
    );
    let shifted = eval(
        ctx,
        &call(
            Symbol::TAYLOR_SHIFT,
            vec![list.clone(), Expr::from(a.clone())],
        ),
    );
    let back = eval(
        ctx,
        &call(Symbol::TAYLOR_SHIFT, vec![shifted.clone(), Expr::from(-a)]),
    );
    assert_eq!(back, list, "TaylorShift of {list:?}: {shifted:?}");
    // Over ℚ(i).
    let factored = eval(ctx, &call(Symbol::COMPLEX_FACTOR, vec![p.clone()]));
    if factored.head() != Some(&Expr::from(Symbol::COMPLEX_FACTOR)) {
        assert!(
            vanishes(ctx, factored.clone() - p.clone()),
            "ComplexFactor[{p:?}] is {factored:?}"
        );
    }
    Ok(())
}

fn euclid(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let (a, b) = (expr(&polynomial(u, 5)?), expr(&polynomial(u, 4)?));
    let found = eval(
        ctx,
        &call(
            Symbol::POLYNOMIAL_EXTENDED_GCD,
            vec![a.clone(), b.clone(), x()],
        ),
    );
    let [s, t, d] = found.args() else {
        panic!("PolynomialExtendedGCD[{a:?}, {b:?}] is {found:?}");
    };
    assert!(
        vanishes(
            ctx,
            s.clone() * a.clone() + t.clone() * b.clone() - d.clone()
        ),
        "u a + v b = d for {a:?}, {b:?}: {found:?}"
    );
    assert!(
        divides(ctx, d, &a) && divides(ctx, d, &b),
        "d divides {a:?}, {b:?}: {found:?}"
    );
    // Bézout for a multiple of d.
    let k = expr(&polynomial(u, 2)?);
    let c = eval(ctx, &Expr::call(Symbol::EXPAND, [d.clone() * k]));
    let found = eval(
        ctx,
        &call(
            Symbol::POLYNOMIAL_BEZOUT,
            vec![a.clone(), b.clone(), c.clone(), x()],
        ),
    );
    let [s, t] = found.args() else {
        panic!("PolynomialBezout[{a:?}, {b:?}, {c:?}] is {found:?}");
    };
    assert!(
        vanishes(
            ctx,
            s.clone() * a.clone() + t.clone() * b.clone() - c.clone()
        ),
        "u a + v b = c for {a:?}, {b:?}, {c:?}: {found:?}"
    );
    assert!(degree(ctx, s) < degree(ctx, &b), "deg u < deg b: {found:?}");
    // Chinese remainders for coprime moduli.
    if degree(ctx, d) == 0 {
        let (r, q) = (expr(&polynomial(u, 3)?), expr(&polynomial(u, 3)?));
        let pair = |v: &Expr, m: &Expr| Expr::call(Symbol::LIST, [v.clone(), m.clone()]);
        let found = eval(
            ctx,
            &call(
                Symbol::POLYNOMIAL_CHINESE_REMAINDER,
                vec![pair(&r, &a), pair(&q, &b)],
            ),
        );
        let [cc, _] = found.args() else {
            panic!("PolynomialChineseRemainder for {a:?}, {b:?} is {found:?}");
        };
        assert!(
            divides(ctx, &a, &(cc.clone() - r.clone()))
                && divides(ctx, &b, &(cc.clone() - q.clone())),
            "the congruences of {found:?} modulo {a:?} and {b:?}"
        );
    }
    // Series division: (a/b) b = a up to xⁿ, b(0) ≠ 0.
    let mut bc = polynomial(u, 3)?;
    if bc[0] == 0 {
        bc[0] = 1;
    }
    let b = expr(&bc);
    let n = u.int_in_range(0..=6i64)?;
    let s = eval(
        ctx,
        &call(
            Symbol::SERIES_DIVISION,
            vec![a.clone(), b.clone(), Expr::from(n)],
        ),
    );
    let rest = eval(
        ctx,
        &Expr::call(Symbol::EXPAND, [s.clone() * b.clone() - a.clone()]),
    );
    for k in 0..=n {
        let coefficient = eval(
            ctx,
            &call(Symbol::COEFFICIENT, vec![rest.clone(), x(), Expr::from(k)]),
        );
        assert_eq!(
            coefficient,
            Expr::from(0),
            "SeriesDivision[{a:?}, {b:?}, {n}] is {s:?}"
        );
    }
    Ok(())
}

fn digits_and_recurrences(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let b = u.int_in_range(2..=16i64)?;
    // Base 2 writes no negative number: the call stays.
    let n = u.int_in_range(if b == 2 { 0 } else { -100_000 }..=100_000i64)?;
    let p = eval(
        ctx,
        &call(
            Symbol::DIGIT_POLYNOMIAL,
            vec![Expr::from(n), Expr::from(b), x()],
        ),
    );
    let at = eval(ctx, &p.substitute(&[(x(), Expr::from(b))]));
    assert_eq!(at, Expr::from(n), "DigitPolynomial[{n}, {b}] is {p:?}");
    let cs = eval(ctx, &call(Symbol::COEFFICIENT_LIST, vec![p.clone(), x()]));
    for d in cs.args() {
        let d = d.as_i64().unwrap_or_else(|| panic!("{p:?}"));
        assert!(
            -b < 2 * d && 2 * d <= b,
            "DigitPolynomial[{n}, {b}] is {p:?}"
        );
    }
    // A sequence of a recurrence of order L ≤ 3.
    let order = u.int_in_range(1..=3usize)?;
    let recurrence: Vec<i64> = (0..order)
        .map(|_| u.int_in_range(-3..=3i64))
        .collect::<Result<_>>()?;
    let mut s: Vec<i64> = (0..order)
        .map(|_| u.int_in_range(-5..=5i64))
        .collect::<Result<_>>()?;
    while s.len() < 2 * order + 2 {
        let k = s.len();
        let next: i64 = (1..=order).map(|i| recurrence[i - 1] * s[k - i]).sum();
        s.push(next);
    }
    let list = Expr::call(
        Symbol::LIST,
        s.iter().map(|v| Expr::from(*v)).collect::<Vec<_>>(),
    );
    let found = eval(ctx, &call(Symbol::BERLEKAMP_MASSEY, vec![list.clone()]));
    let c: Vec<Rational> = found
        .args()
        .iter()
        .map(|v| {
            v.to_rational()
                .unwrap_or_else(|| panic!("BerlekampMassey[{list:?}] is {found:?}"))
        })
        .collect();
    let l = c.len() - 1;
    assert!(
        l <= order,
        "BerlekampMassey[{list:?}] is {found:?}, longer than {order}"
    );
    for k in l..s.len() {
        let sum = (0..=l).fold(Rational::from(Integer::ZERO), |acc, i| {
            &acc + &(&c[i] * &Rational::from(Integer::from(s[k - i])))
        });
        assert!(
            sum.is_zero(),
            "BerlekampMassey[{list:?}] is {found:?}, not a recurrence at {k}"
        );
    }
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0..=2u8)? {
        0 => roots(&ctx, u),
        1 => euclid(&ctx, u),
        _ => digits_and_recurrences(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
