//! Modular arithmetic (CAS-DOC-43 §6, CAS-454), one check per input:
//! - residues: `Modular[a, n]` of an integer or a rational is the
//!   symmetric residue, `Indeterminate` when the denominator has no
//!   inverse; sums, products and powers (inverses among them) of residues
//!   are the residues of the integers'; two moduli leave a sum;
//! - polynomials over `𝔽ₚ`, `p` in {2, 3, 5, 7, 11, 101}, their
//!   coefficients any representatives, against the target's own
//!   arithmetic: `Factor`'s factors are irreducible (Ben-Or's test) and
//!   distinct, and multiply back; `PolynomialQuotientRemainder` is
//!   Euclid's; `GCD` of `a h` and `b h` is their monic gcd.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use cas_num::{Integer, Rational};
use libfuzzer_sys::fuzz_target;

fn modular(a: Expr, n: impl Into<Integer>) -> Expr {
    Expr::call(Symbol::MODULAR, [a, Expr::from(n.into())])
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// `Modular[r, n]`, `r` the symmetric residue of `a` modulo `n`.
fn residue(a: &Integer, n: &Integer) -> Expr {
    let (_, r) = a.div_rem_euclid(n);
    let r = if &(&r + &r) > n { &r - n } else { r };
    modular(Expr::from(r), n.clone())
}

fn residues(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let n = if u.arbitrary()? {
        u.int_in_range(2..=1000i64)?
    } else {
        u.int_in_range(2..=i64::MAX)?
    };
    let n = Integer::from(n);
    let a = Integer::from(u.arbitrary::<i64>()?);
    let b = Integer::from(u.arbitrary::<i64>()?);
    let ma = modular(Expr::from(a.clone()), n.clone());
    let mb = modular(Expr::from(b.clone()), n.clone());
    assert_eq!(eval(ctx, &ma), residue(&a, &n), "{ma:?}");
    let sum = Expr::call(Symbol::PLUS, [ma.clone(), mb.clone()]);
    assert_eq!(eval(ctx, &sum), residue(&(&a + &b), &n), "{sum:?}");
    let product = Expr::call(Symbol::TIMES, [ma.clone(), mb]);
    assert_eq!(eval(ctx, &product), residue(&(&a * &b), &n), "{product:?}");
    if !b.is_zero() {
        let q = Rational::new(a.clone(), b.clone());
        let call = modular(Expr::from(q.clone()), n.clone());
        let expected = match q.denom().mod_inverse(&n) {
            Some(inverse) => residue(&(q.numer() * &inverse), &n),
            None => Expr::from(Symbol::INDETERMINATE),
        };
        assert_eq!(eval(ctx, &call), expected, "{call:?}");
    }
    let k = u.int_in_range(-1000..=1_000_000_000i64)?;
    let power = Expr::call(Symbol::POWER, [ma.clone(), Expr::from(k)]);
    let expected = match a.pow_mod(&Integer::from(k), &n) {
        Some(r) => residue(&r, &n),
        None => Expr::from(Symbol::INDETERMINATE),
    };
    assert_eq!(eval(ctx, &power), expected, "{power:?}");
    // Another modulus: the sum stays.
    let other = modular(Expr::from(b), &n + &Integer::ONE);
    let mixed = Expr::call(Symbol::PLUS, [ma, other]);
    let found = eval(ctx, &mixed);
    assert_eq!(
        found.head(),
        Some(&Expr::from(Symbol::PLUS)),
        "{mixed:?} is {found:?}"
    );
    Ok(())
}

/// A polynomial over `𝔽ₚ`: its coefficients in `0..p` from the constant
/// up, without trailing zeros, so that 0 is empty.
type Poly = Vec<u64>;

fn trim(mut f: Poly) -> Poly {
    while f.last() == Some(&0) {
        f.pop();
    }
    f
}

/// `a⁻¹` modulo the prime `p`, by Fermat.
fn inverse(a: u64, p: u64) -> u64 {
    let (mut result, mut base, mut e) = (1, a % p, p - 2);
    while e > 0 {
        if e & 1 == 1 {
            result = result * base % p;
        }
        base = base * base % p;
        e >>= 1;
    }
    result
}

fn mul(f: &[u64], g: &[u64], p: u64) -> Poly {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let mut h = vec![0; f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        for (j, b) in g.iter().enumerate() {
            h[i + j] = (h[i + j] + a * b) % p;
        }
    }
    trim(h)
}

fn sub(f: &[u64], g: &[u64], p: u64) -> Poly {
    let mut h = vec![0; f.len().max(g.len())];
    h[..f.len()].copy_from_slice(f);
    for (x, b) in h.iter_mut().zip(g) {
        *x = (*x + p - b) % p;
    }
    trim(h)
}

/// `(q, r)` with `f = q g + r` and `deg r < deg g`, `g` not 0.
fn div_rem(f: &[u64], g: &[u64], p: u64) -> (Poly, Poly) {
    let lead = inverse(*g.last().expect("not 0"), p);
    let mut q = vec![0; f.len().saturating_sub(g.len() - 1)];
    let mut r = f.to_vec();
    while r.len() >= g.len() {
        let shift = r.len() - g.len();
        let c = r.last().expect("not 0") * lead % p;
        q[shift] = c;
        for (x, b) in r[shift..].iter_mut().zip(g) {
            *x = (*x + p - c * b % p) % p;
        }
        r = trim(r);
    }
    (trim(q), r)
}

fn monic(f: Poly, p: u64) -> Poly {
    match f.last() {
        Some(&c) => {
            let c = inverse(c, p);
            f.iter().map(|x| x * c % p).collect()
        }
        None => f,
    }
}

/// The monic gcd.
fn gcd(f: &[u64], g: &[u64], p: u64) -> Poly {
    let (mut a, mut b) = (f.to_vec(), g.to_vec());
    while !b.is_empty() {
        let (_, r) = div_rem(&a, &b, p);
        a = b;
        b = r;
    }
    monic(a, p)
}

/// `h^e` modulo `f`.
fn power_mod(h: &[u64], mut e: u64, f: &[u64], p: u64) -> Poly {
    let mut result = vec![1];
    let mut base = div_rem(h, f, p).1;
    while e > 0 {
        if e & 1 == 1 {
            result = div_rem(&mul(&result, &base, p), f, p).1;
        }
        base = div_rem(&mul(&base, &base, p), f, p).1;
        e >>= 1;
    }
    result
}

/// Ben-Or's test: `f`, of degree `d ≥ 1`, is irreducible when
/// `gcd(f, x^(pⁱ) − x) = 1` for each `i ≤ d/2`.
fn irreducible(f: &[u64], p: u64) -> bool {
    let x = vec![0, 1];
    let mut h = x.clone();
    for _ in 0..(f.len() - 1) / 2 {
        h = power_mod(&h, p, f, p);
        if gcd(f, &sub(&h, &x, p), p).len() > 1 {
            return false;
        }
    }
    true
}

/// A polynomial of degree `d`, monic when `monic`.
fn random(u: &mut Unstructured<'_>, p: u64, d: usize, monic: bool) -> Result<Poly> {
    let mut f = (0..d)
        .map(|_| u.int_in_range(0..=p - 1))
        .collect::<Result<Poly>>()?;
    f.push(if monic { 1 } else { u.int_in_range(1..=p - 1)? });
    Ok(f)
}

/// `f` as an expression in `x`, each coefficient another representative
/// of its class, from `c − 2p` to `c + 2p`.
fn expression(u: &mut Unstructured<'_>, f: &[u64], p: u64) -> Result<Expr> {
    let x = Expr::symbol("x");
    let mut terms = Vec::with_capacity(f.len());
    for (k, c) in f.iter().enumerate() {
        let shift = Integer::from(u.int_in_range(-2..=2i64)?);
        let c = &Integer::from(*c) + &(&shift * &Integer::from(p));
        terms.push(Expr::from(c) * x.pow(i64::try_from(k).expect("a degree")));
    }
    Ok(Expr::call(Symbol::PLUS, terms))
}

/// The coefficients modulo `p` of a polynomial in `x` with integer
/// coefficients, by `CoefficientList`.
fn coefficients(ctx: &Context, f: &Expr, p: u64) -> Poly {
    let call = Expr::call(Symbol::COEFFICIENT_LIST, [f.clone(), Expr::symbol("x")]);
    let list = eval(ctx, &call);
    let modulus = Integer::from(p);
    // From the leading coefficient down, as giac's `symb2poly`.
    let coeffs = list.args().iter().rev().map(|c| {
        let c = c
            .to_integer()
            .unwrap_or_else(|| panic!("{call:?} is {list:?}"));
        u64::try_from(&c.div_rem_euclid(&modulus).1).expect("a residue")
    });
    trim(coeffs.collect())
}

/// The polynomial `f` of `Modular[f, p]`.
fn inner(ctx: &Context, e: &Expr, p: u64) -> Poly {
    match (e.head(), e.args()) {
        (Some(h), [f, n])
            if *h == Expr::from(Symbol::MODULAR) && *n == Expr::from(Integer::from(p)) =>
        {
            coefficients(ctx, f, p)
        }
        _ => panic!("{e:?} is no polynomial modulo {p}"),
    }
}

fn fields(ctx: &Context, u: &mut Unstructured<'_>) -> Result<()> {
    let p = *u.choose(&[2u64, 3, 5, 7, 11, 101])?;
    let d = u.int_in_range(1..=12usize)?;
    let f = random(u, p, d, false)?;
    let fe = expression(u, &f, p)?;
    // Factor: a constant times distinct irreducibles to powers.
    let call = Expr::call(Symbol::FACTOR, [modular(fe.clone(), p)]);
    let factored = eval(ctx, &call);
    let items = match factored.head() {
        Some(h) if *h == Expr::from(Symbol::TIMES) => factored.args().to_vec(),
        _ => vec![factored.clone()],
    };
    let mut product: Poly = vec![1];
    let mut seen: Vec<Poly> = Vec::new();
    for item in &items {
        if let Some(c) = item.to_integer() {
            let c = coefficients(ctx, &Expr::from(c), p);
            assert!(!c.is_empty(), "{call:?} is {factored:?}");
            product = mul(&product, &c, p);
            continue;
        }
        let (base, k) = match (item.head(), item.args()) {
            (Some(h), [b, k]) if *h == Expr::from(Symbol::POWER) => (b, k.as_i64()),
            _ => (item, Some(1)),
        };
        let k = k
            .and_then(|k| usize::try_from(k).ok())
            .filter(|k| (1..=d).contains(k))
            .unwrap_or_else(|| panic!("{call:?} is {factored:?}"));
        let g = coefficients(ctx, base, p);
        assert!(g.len() >= 2, "{call:?} is {factored:?}: {base:?}");
        let lead = *g.last().expect("not 0");
        let g = monic(g, p);
        assert!(irreducible(&g, p), "{call:?}: {base:?} is reducible");
        assert!(!seen.contains(&g), "{call:?}: {base:?} twice");
        for _ in 0..k {
            product = mul(&product, &[lead], p);
            product = mul(&product, &g, p);
        }
        seen.push(g);
    }
    assert_eq!(
        product, f,
        "{call:?} is {factored:?}, which does not multiply back"
    );
    // Euclid's quotient and remainder.
    let e = u.int_in_range(0..=d)?;
    let g = random(u, p, e, false)?;
    let ge = expression(u, &g, p)?;
    let x = Expr::symbol("x");
    let call = Expr::call(
        Symbol::POLYNOMIAL_QUOTIENT_REMAINDER,
        [modular(fe, p), modular(ge, p), x],
    );
    let found = eval(ctx, &call);
    let [q, r] = found.args() else {
        panic!("{call:?} is {found:?}");
    };
    assert_eq!(
        (inner(ctx, q, p), inner(ctx, r, p)),
        div_rem(&f, &g, p),
        "{call:?} is {found:?}"
    );
    // The monic gcd of a h and b h.
    let degrees = [
        u.int_in_range(1..=4usize)?,
        u.int_in_range(0..=4usize)?,
        u.int_in_range(0..=4usize)?,
    ];
    let h = random(u, p, degrees[0], true)?;
    let a = random(u, p, degrees[1], false)?;
    let b = random(u, p, degrees[2], false)?;
    let (ah, bh) = (mul(&a, &h, p), mul(&b, &h, p));
    let call = Expr::call(
        Symbol::GCD,
        [
            modular(expression(u, &ah, p)?, p),
            modular(expression(u, &bh, p)?, p),
        ],
    );
    let found = eval(ctx, &call);
    assert_eq!(
        inner(ctx, &found, p),
        gcd(&ah, &bh, p),
        "{call:?} is {found:?}"
    );
    Ok(())
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    if u.arbitrary()? {
        residues(&ctx, u)
    } else {
        fields(&ctx, u)
    }
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
