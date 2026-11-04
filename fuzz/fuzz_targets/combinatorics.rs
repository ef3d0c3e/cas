//! Permutations, bits and transforms (CAS-DOC-45 §5, CAS-466), one check
//! per input, against the target's own arithmetic:
//! - random permutations up to size 12: `PermutationQ`, the cycles giving
//!   the permutation back, `p ∘ p⁻¹` the identity, the order the least `k`
//!   with `pᵏ` the identity, signatures multiplicative, and the matrix of
//!   `p ∘ q` the product of `q`'s and `p`'s;
//! - `BitAnd`, `BitOr`, `BitXor` and `HammingDistance` against Rust's
//!   `i128` operations, `Convert` both ways, `FactorialPower` against the
//!   product;
//! - `DiscreteFourierTransform`: the inverse's balls holding the input,
//!   the first entry the sum, and Parseval's identity; modulo a prime with
//!   a root of a dividing order, against the direct sums and back;
//!   `Convolution` against the schoolbook sums.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn list(items: impl IntoIterator<Item = Expr>) -> Expr {
    Expr::call(Symbol::LIST, items.into_iter().collect::<Vec<_>>())
}

fn indices(p: &[usize]) -> Expr {
    list(
        p.iter()
            .map(|&i| Expr::from(i64::try_from(i).expect("small"))),
    )
}

fn integers(v: &[i64]) -> Expr {
    list(v.iter().map(|&x| Expr::from(x)))
}

fn read_indices(e: &Expr) -> Vec<usize> {
    e.args()
        .iter()
        .map(|x| {
            x.as_i64()
                .and_then(|k| usize::try_from(k).ok())
                .unwrap_or_else(|| panic!("{e:?} is no permutation"))
        })
        .collect()
}

fn call(s: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(s, args)
}

/// A random permutation of size `n`, by Fisher and Yates.
fn permutation(u: &mut Unstructured<'_>, n: usize) -> Result<Vec<usize>> {
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = u.int_in_range(0..=i)?;
        p.swap(i, j);
    }
    Ok(p)
}

fn compose(p: &[usize], q: &[usize]) -> Vec<usize> {
    q.iter().map(|&i| p[i]).collect()
}

fn matrix(p: &[usize]) -> Vec<Vec<i64>> {
    p.iter()
        .map(|&j| (0..p.len()).map(|k| i64::from(k == j)).collect())
        .collect()
}

fn product(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
    a.iter()
        .map(|row| {
            (0..b[0].len())
                .map(|j| row.iter().zip(b).map(|(x, r)| x * r[j]).sum())
                .collect()
        })
        .collect()
}

fn permutations(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(0..=12usize)?;
    let p = permutation(u, n)?;
    let q = permutation(u, n)?;
    let pe = indices(&p);
    let identity: Vec<usize> = (0..n).collect();
    let is = eval(ctx, &call(Symbol::PERMUTATION_Q, vec![pe.clone()]));
    assert_eq!(is, Expr::from(Symbol::TRUE), "PermutationQ[{pe:?}]");
    if n >= 2 {
        let mut bad = p.clone();
        bad[0] = bad[1];
        let is = eval(ctx, &call(Symbol::PERMUTATION_Q, vec![indices(&bad)]));
        assert_eq!(is, Expr::from(Symbol::FALSE), "PermutationQ[{bad:?}]");
    }
    let cycles = eval(ctx, &call(Symbol::PERMUTATION_CYCLES, vec![pe.clone()]));
    let back = eval(ctx, &call(Symbol::CYCLES_PERMUTATION, vec![cycles.clone()]));
    assert_eq!(back, pe, "{pe:?}: the cycles {cycles:?} give {back:?}");
    let inverse = eval(ctx, &call(Symbol::INVERSE_PERMUTATION, vec![pe.clone()]));
    let both = eval(
        ctx,
        &call(
            Symbol::PERMUTATION_PRODUCT,
            vec![pe.clone(), inverse.clone()],
        ),
    );
    assert_eq!(
        read_indices(&both),
        identity,
        "{pe:?} times its inverse {inverse:?}"
    );
    let pq = eval(
        ctx,
        &call(Symbol::PERMUTATION_PRODUCT, vec![pe.clone(), indices(&q)]),
    );
    assert_eq!(read_indices(&pq), compose(&p, &q), "{pe:?} times {q:?}");
    // The order.
    let order = eval(ctx, &call(Symbol::PERMUTATION_ORDER, vec![pe.clone()]))
        .as_i64()
        .unwrap_or_else(|| panic!("PermutationOrder[{pe:?}]"));
    let mut power = p.clone();
    for k in 1..order {
        assert_ne!(
            power, identity,
            "PermutationOrder[{pe:?}] is {order}, not the least: {k}"
        );
        power = compose(&power, &p);
    }
    assert_eq!(power, identity, "PermutationOrder[{pe:?}] is {order}");
    // Signatures.
    let sign = |r: &Expr| {
        eval(ctx, &call(Symbol::PERMUTATION_SIGNATURE, vec![r.clone()]))
            .as_i64()
            .unwrap_or_else(|| panic!("PermutationSignature[{r:?}]"))
    };
    assert_eq!(
        sign(&pq),
        sign(&pe) * sign(&indices(&q)),
        "signatures of {pe:?} and {q:?}"
    );
    // Matrices.
    if n > 0 {
        let m = eval(ctx, &call(Symbol::PERMUTATION_MATRIX, vec![pq]));
        let rows: Vec<Vec<i64>> = m
            .args()
            .iter()
            .map(|r| {
                r.args()
                    .iter()
                    .map(|x| x.as_i64().expect("an entry"))
                    .collect()
            })
            .collect();
        assert_eq!(
            rows,
            product(&matrix(&q), &matrix(&p)),
            "PermutationMatrix of {p:?} ∘ {q:?}"
        );
    }
    Ok(())
}

fn bits(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let (a, b) = (
        i128::from(u.arbitrary::<i64>()?),
        i128::from(u.arbitrary::<i64>()?),
    );
    let (ea, eb) = (
        Expr::from(i64::try_from(a).expect("i64")),
        Expr::from(i64::try_from(b).expect("i64")),
    );
    let found = |s: Symbol| eval(ctx, &call(s, vec![ea.clone(), eb.clone()]));
    let int = |x: i128| {
        eval(
            ctx,
            &Expr::call(Symbol::PLUS, [Expr::from(i64::try_from(x).expect("i64"))]),
        )
    };
    assert_eq!(found(Symbol::BIT_AND), int(a & b), "BitAnd[{a}, {b}]");
    assert_eq!(found(Symbol::BIT_OR), int(a | b), "BitOr[{a}, {b}]");
    assert_eq!(found(Symbol::BIT_XOR), int(a ^ b), "BitXor[{a}, {b}]");
    if a >= 0 && b >= 0 {
        let ones = i64::from((a ^ b).count_ones());
        assert_eq!(
            found(Symbol::HAMMING_DISTANCE),
            Expr::from(ones),
            "HammingDistance[{a}, {b}]"
        );
    }
    // Digits and back.
    let base = u.int_in_range(2..=36i64)?;
    let digits = eval(
        ctx,
        &call(
            Symbol::CONVERT,
            vec![ea.clone(), Expr::symbol("base"), Expr::from(base)],
        ),
    );
    for d in digits.args() {
        let d = d
            .as_i64()
            .unwrap_or_else(|| panic!("Convert[{a}, base, {base}] is {digits:?}"));
        assert!(
            d.abs() < base && (a >= 0 || d <= 0),
            "Convert[{a}, base, {base}] is {digits:?}"
        );
    }
    let back = eval(
        ctx,
        &call(
            Symbol::CONVERT,
            vec![digits.clone(), Expr::symbol("base"), Expr::from(base)],
        ),
    );
    assert_eq!(back, ea, "Convert of {digits:?} in base {base}");
    // Falling factorials.
    let n = u.int_in_range(-20..=20i64)?;
    let k = u.int_in_range(0..=10i64)?;
    let expected = (0..k).fold(1i128, |acc, j| acc * i128::from(n - j));
    let found = eval(
        ctx,
        &call(Symbol::FACTORIAL_POWER, vec![Expr::from(n), Expr::from(k)]),
    );
    assert_eq!(found, int(expected), "FactorialPower[{n}, {k}]");
    Ok(())
}

/// Whether `|a − b| < 10⁻⁸` by the CAS's balls, `b` exact: no `N`, whose
/// precision loop finds no relative accuracy in a difference near 0.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::PLUS, [a.clone(), Expr::from(-1) * b.clone()]);
    // Real and imaginary parts apart: the modulus of a complex ball about
    // 0 may not be decided.
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let size = Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]);
        let small = Expr::call(
            Symbol::LESS,
            [
                size,
                Expr::call(Symbol::POWER, [Expr::from(10), Expr::from(-8)]),
            ],
        );
        eval(ctx, &small) == Expr::from(Symbol::TRUE)
    })
}

fn transforms(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = u.int_in_range(1..=16usize)?;
    let x: Vec<i64> = (0..n)
        .map(|_| u.int_in_range(-100..=100i64))
        .collect::<Result<_>>()?;
    let xe = integers(&x);
    let big = eval(
        ctx,
        &call(Symbol::DISCRETE_FOURIER_TRANSFORM, vec![xe.clone()]),
    );
    assert_eq!(
        big.args().len(),
        n,
        "DiscreteFourierTransform[{xe:?}] is {big:?}"
    );
    let sum: i64 = x.iter().sum();
    assert!(
        near(ctx, &big.args()[0], &Expr::from(sum)),
        "DiscreteFourierTransform[{xe:?}] is {big:?}"
    );
    let back = eval(
        ctx,
        &call(
            Symbol::INVERSE_DISCRETE_FOURIER_TRANSFORM,
            vec![big.clone()],
        ),
    );
    assert_eq!(back.args().len(), n, "the inverse of {big:?} is {back:?}");
    for (f, v) in back.args().iter().zip(&x) {
        assert!(
            near(ctx, f, &Expr::from(*v)),
            "the inverse of {big:?} is {back:?}, not {xe:?}"
        );
    }
    // Parseval: Σ |X_k|² = n Σ x_j².
    let energy = Expr::call(
        Symbol::PLUS,
        big.args()
            .iter()
            .map(|z| {
                Expr::call(
                    Symbol::POWER,
                    [Expr::call(Symbol::ABS, [z.clone()]), Expr::from(2)],
                )
            })
            .collect::<Vec<_>>(),
    );
    let squares: i64 = x.iter().map(|v| v * v).sum();
    let n_i = i64::try_from(n).expect("small");
    assert!(
        near(ctx, &eval(ctx, &energy), &Expr::from(n_i * squares)),
        "Parseval for {xe:?}: {big:?}"
    );
    // Modulo a prime, with a root of an order that divides p − 1.
    let p = *u.choose(&[5i64, 13, 17, 97, 257])?;
    let divisors: Vec<i64> = (1..p).filter(|d| (p - 1) % d == 0 && *d <= 16).collect();
    let m = *u.choose(&divisors)?;
    let order = |g: i64| {
        (1..p)
            .find(|&k| (0..k).fold(1, |acc, _| acc * g % p) == 1)
            .expect("a unit")
    };
    let g = (2..p)
        .find(|&g| order(g) == p - 1)
        .expect("a primitive root");
    let w = (0..(p - 1) / m).fold(1, |acc, _| acc * g % p);
    let y: Vec<i64> = (0..m)
        .map(|_| u.int_in_range(0..=p - 1))
        .collect::<Result<_>>()?;
    // giac's convention: the transform with w⁻¹.
    let w_inverse = (1..p).find(|v| v * w % p == 1).expect("a unit");
    let powers: Vec<i64> = (0..m)
        .scan(1, |c, _| {
            let v = *c;
            *c = *c * w_inverse % p;
            Some(v)
        })
        .collect();
    let mu = usize::try_from(m).expect("small");
    let expected: Vec<i64> = (0..mu)
        .map(|k| {
            (0..mu)
                .map(|j| y[j] * powers[(j * k) % mu])
                .sum::<i64>()
                .rem_euclid(p)
        })
        .collect();
    let mod_call = call(
        Symbol::DISCRETE_FOURIER_TRANSFORM,
        vec![integers(&y), Expr::from(w), Expr::from(p)],
    );
    let found = eval(ctx, &mod_call);
    assert_eq!(found, integers(&expected), "{mod_call:?}");
    let back = eval(
        ctx,
        &call(
            Symbol::INVERSE_DISCRETE_FOURIER_TRANSFORM,
            vec![found, Expr::from(w), Expr::from(p)],
        ),
    );
    assert_eq!(back, integers(&y), "the inverse of {mod_call:?}");
    // Convolution.
    let (la, lb) = (u.int_in_range(1..=8usize)?, u.int_in_range(1..=8usize)?);
    let a: Vec<i64> = (0..la)
        .map(|_| u.int_in_range(-50..=50i64))
        .collect::<Result<_>>()?;
    let b: Vec<i64> = (0..lb)
        .map(|_| u.int_in_range(-50..=50i64))
        .collect::<Result<_>>()?;
    let mut c = vec![0i64; la + lb - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            c[i + j] += ai * bj;
        }
    }
    let conv = call(Symbol::CONVOLUTION, vec![integers(&a), integers(&b)]);
    assert_eq!(eval(ctx, &conv), integers(&c), "{conv:?}");
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    match u.int_in_range(0..=2u8)? {
        0 => permutations(&ctx, u),
        1 => bits(&ctx, u),
        _ => transforms(&ctx, u),
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
