//! Diophantine equations (CAS-DOC-39 §4, CAS-425): `IntegerSolveLinear`
//! (`iabcuv`), `PrimeSumOfSquares` (`pa2b2`) and `IntegerSolve` (`isolve`).
//!
//! - `IntegerSolveLinear[a, b, c]`: `{u, v}` with `a u + b v = c`, `(c/g)`
//!   times the extended gcd's cofactors, giac's `[-2,2]` for `(21, 28,
//!   14)`; `g ∤ c` stays, with giac's "no solution".
//! - `PrimeSumOfSquares[p]`: `{a, b}`, `a² + b² = p`, `a > b > 0`, for a
//!   prime `p ≡ 1 mod 4` (and `{1, 1}` for 2), by Cornacchia's algorithm:
//!   `x₀² ≡ −1` from a quadratic nonresidue `r` as `r^((p−1)/4)`, then
//!   Euclid's algorithm on `(p, x₀)` to the first remainder below `√p`.
//! - `IntegerSolve[eqs]` and `IntegerSolve[eqs, vars]`, the variables the
//!   equations' symbols in the order they first appear when not given:
//!   - one equation `a x + b y = c`, `a b ≠ 0`: giac's `x = u + (b/g) n_0`,
//!     `y = v − (a/g) n_0`, `u` in `(−|b/g|/2, |b/g|/2]`, so its
//!     `[x=-1+5*_Z0,y=2-3*_Z0]` for `3x + 5y = 7` and `[x=1-2*_Z0,y=-_Z0]`
//!     for `x − 2y = 1`;
//!   - other linear equations and systems with integer coefficients, by
//!     column operations of determinant 1 to an echelon form `A U = H`, `x =
//!     U y`: the pivots' `y` solved in turn, each an integer or no solution,
//!     the others parameters `n_0, n_1, …`; with one parameter, its first
//!     nonzero coefficient made positive, `c`, and that variable's constant
//!     put in `(−c/2, c/2]`. giac's own answers there have redundant
//!     parameters, `isolve(2x+3y+5z=1)` three. The answer is giac's list `{x
//!     = …, …}`, `{}` for none;
//!   - `x y = n`: the pairs of divisors, both signs, giac's order, the
//!     divisors to `√|n|` and their cofactors;
//!   - `x² + y² = n`, `n ≤` [`MAX_SQUARES`]: every representation, `x` then
//!     `y` ascending;
//!   - Pell's `x² − D y² = 1`, `D > 0` not a square: `x + y √D = (x₁ + y₁
//!     √D)^n`, `n` a natural parameter `n_0`, giac's form, from the
//!     fundamental solution `(x₁, y₁)`, the convergent at the end of the
//!     period of `√D`'s continued fraction (CAS-424), where `Q` is 1, or of
//!     two periods for an odd period; past [`MAX_PELL`] quotients the
//!     equation stays;
//!   - other equations stay.

use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// The largest `n` of `x² + y² = n`, enumerated.
const MAX_SQUARES: i64 = 1 << 40;

/// The largest `n` of `x y = n`.
const MAX_PRODUCT: i64 = 1 << 40;

/// The most quotients of `√D` for Pell's equation.
const MAX_PELL: usize = 1 << 16;

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &str) {
    ev.warn(Warning::Native {
        message: message.to_owned().into(),
        expr: e.clone(),
    });
}

/// `IntegerSolveLinear[a, b, c]`.
pub(in crate::native) fn integer_solve_linear(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b, c] = e.args() else {
        return Ok(None);
    };
    let (Some(a), Some(b), Some(c)) = (a.to_integer(), b.to_integer(), c.to_integer()) else {
        return Ok(None);
    };
    let (g, s, t) = a.gcd_ext(&b);
    if g.is_zero() {
        return Ok(None);
    }
    let (k, rest) = c.div_rem_euclid(&g);
    if !rest.is_zero() {
        warn(ev, e, "no solution in integers");
        return Ok(None);
    }
    Ok(Some(list(vec![Expr::from(&s * &k), Expr::from(&t * &k)])))
}

/// `x₀` with `x₀² ≡ −1 (mod p)`, for a prime `p ≡ 1 mod 4`.
fn root_of_minus_one(p: &Integer) -> Option<Integer> {
    let exponent = (p - &Integer::ONE).div_exact(&Integer::from(4));
    let minus_one = p - &Integer::ONE;
    let mut r = Integer::from(2);
    while r < *p {
        let x = r.pow_mod(&exponent, p)?;
        if (&(&x * &x) % p) == minus_one {
            return Some(x);
        }
        r = &r + &Integer::ONE;
    }
    None
}

/// `PrimeSumOfSquares[p]`.
pub(in crate::native) fn prime_sum_of_squares(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [p] = e.args() else {
        return Ok(None);
    };
    let Some(p) = p.to_integer() else {
        return Ok(None);
    };
    if p == Integer::from(2) {
        return Ok(Some(list(vec![Expr::from(1), Expr::from(1)])));
    }
    let four = Integer::from(4);
    let prime = p > Integer::ONE && p.unsigned_abs().primality().is_probable_prime();
    if !prime || (&p % &four) != Integer::ONE {
        warn(ev, e, "not a prime 1 mod 4");
        return Ok(None);
    }
    let Some(x0) = root_of_minus_one(&p) else {
        return Ok(None);
    };
    // Euclid on (p, x₀) to the first remainder below √p.
    let (root, _) = p.sqrt_rem();
    let (mut a, mut b) = (p.clone(), x0);
    while b > root {
        let r = &a % &b;
        a = b;
        b = r;
    }
    let (other, rest) = (&p - &(&b * &b)).sqrt_rem();
    if !rest.is_zero() {
        return Ok(None);
    }
    let (big, small) = if b > other { (b, other) } else { (other, b) };
    Ok(Some(list(vec![Expr::from(big), Expr::from(small)])))
}

/// The symbols of `e` that are not builtins, in the order they first
/// appear.
fn symbols_in_order(e: &Expr, found: &mut Vec<Expr>) {
    if let Some(s) = e.as_symbol() {
        if !s.is_builtin() && !found.contains(e) {
            found.push(e.clone());
        }
        return;
    }
    for a in e.args() {
        symbols_in_order(a, found);
    }
}

/// `e` as integer coefficients of `xs` and a constant; `None` when it is
/// not linear so.
fn linear(
    ev: &mut Eval<'_>,
    e: &Expr,
    xs: &[Expr],
) -> Result<Option<(Vec<Integer>, Integer)>, EvalError> {
    let mut coefficients = Vec::with_capacity(xs.len());
    for x in xs {
        let a = ev.eval(&Expr::call(Symbol::D, [e.clone(), x.clone()]))?;
        let Some(a) = a.to_integer() else {
            return Ok(None);
        };
        coefficients.push(a);
    }
    let zero: Vec<(Expr, Expr)> = xs.iter().map(|x| (x.clone(), Expr::from(0))).collect();
    let Some(k) = ev.eval(&e.substitute(&zero))?.to_integer() else {
        return Ok(None);
    };
    let mut rebuilt = Expr::from(k.clone());
    for (a, x) in coefficients.iter().zip(xs) {
        rebuilt += Expr::from(a.clone()) * x.clone();
    }
    let gap = ev.eval(&Expr::call(Symbol::EXPAND, [e.clone() - rebuilt]))?;
    Ok((gap == Expr::from(0)).then_some((coefficients, k)))
}

/// The parameter `n_k`.
fn parameter(k: usize) -> Expr {
    Expr::symbol(&format!("n_{k}"))
}

/// The general solution of `A x = b` in integers: per variable, a constant
/// and the coefficients of the parameters; `None` for no solution.
fn linear_system(a: &[Vec<Integer>], b: &[Integer]) -> Option<Vec<(Integer, Vec<Integer>)>> {
    let (m, n) = (a.len(), a.first().map_or(0, Vec::len));
    let mut h: Vec<Vec<Integer>> = a.to_vec();
    // U's columns, as rows of `u`: u[j] is column j.
    let mut u: Vec<Vec<Integer>> = (0..n)
        .map(|j| (0..n).map(|i| Integer::from(i64::from(i == j))).collect())
        .collect();
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut r = 0;
    for i in 0..m {
        if r == n {
            break;
        }
        for j in r + 1..n {
            if h[i][j].is_zero() {
                continue;
            }
            let (x, y) = (h[i][r].clone(), h[i][j].clone());
            let (g, s, t) = x.gcd_ext(&y);
            let (xg, yg) = (x.div_exact(&g), y.div_exact(&g));
            // col_r ← s col_r + t col_j, col_j ← −(y/g) col_r + (x/g) col_j.
            for row in &mut h {
                let (cr, cj) = (row[r].clone(), row[j].clone());
                row[r] = &(&s * &cr) + &(&t * &cj);
                row[j] = &(&xg * &cj) - &(&yg * &cr);
            }
            let (ur, uj) = (u[r].clone(), u[j].clone());
            u[r] = ur
                .iter()
                .zip(&uj)
                .map(|(p, q)| &(&s * p) + &(&t * q))
                .collect();
            u[j] = ur
                .iter()
                .zip(&uj)
                .map(|(p, q)| &(&xg * q) - &(&yg * p))
                .collect();
        }
        if !h[i][r].is_zero() {
            pivots.push((i, r));
            r += 1;
        }
    }
    // H y = b, the pivots' y in turn.
    let mut y: Vec<Option<Integer>> = vec![None; n];
    let mut row = 0;
    for i in 0..m {
        let known = |y: &[Option<Integer>]| -> Integer {
            (0..n)
                .filter_map(|j| y[j].as_ref().map(|v| &h[i][j] * v))
                .fold(Integer::ZERO, |s, v| &s + &v)
        };
        match pivots.get(row) {
            Some(&(pi, pc)) if pi == i => {
                let rest = &b[i] - &known(&y);
                let (q, rem) = rest.div_rem_euclid(&h[i][pc]);
                if !rem.is_zero() {
                    return None;
                }
                y[pc] = Some(q);
                row += 1;
            }
            _ => {
                if known(&y) != b[i] {
                    return None;
                }
            }
        }
    }
    // x = U y, the free y parameters.
    let free: Vec<usize> = (0..n).filter(|&j| y[j].is_none()).collect();
    Some(
        (0..n)
            .map(|v| {
                let constant = (0..n)
                    .filter_map(|j| y[j].as_ref().map(|value| &u[j][v] * value))
                    .fold(Integer::ZERO, |s, t| &s + &t);
                let coefficients = free.iter().map(|&j| u[j][v].clone()).collect();
                (constant, coefficients)
            })
            .collect(),
    )
}

/// With one parameter, its first nonzero coefficient made positive, `c`,
/// and that variable's constant put in `(−c/2, c/2]`.
fn normalized(mut solution: Vec<(Integer, Vec<Integer>)>) -> Vec<(Integer, Vec<Integer>)> {
    if solution.first().is_none_or(|(_, c)| c.len() != 1) {
        return solution;
    }
    let Some(first) = solution.iter().position(|(_, c)| !c[0].is_zero()) else {
        return solution;
    };
    if solution[first].1[0].is_negative() {
        for (_, coefficients) in &mut solution {
            coefficients[0] = -coefficients[0].clone();
        }
    }
    let c = solution[first].1[0].clone();
    // The parameter shifted by k: constant − k c, the symmetric residue.
    let (_, mut target) = solution[first].0.div_rem_euclid(&c);
    if &target + &target > c {
        target = &target - &c;
    }
    let k = (&solution[first].0 - &target).div_exact(&c);
    for (constant, coefficients) in &mut solution {
        *constant = &*constant - &(&k * &coefficients[0]);
    }
    solution
}

/// `a x + b y = c`, `a b ≠ 0`: giac's `x = u + (b/g) n_0`, `y = v − (a/g)
/// n_0`, `u` in `(−|b/g|/2, |b/g|/2]`; `None` when `g ∤ c`.
fn two_variables(a: &Integer, b: &Integer, c: &Integer) -> Option<Vec<(Integer, Vec<Integer>)>> {
    let (g, s, _) = a.gcd_ext(b);
    let (k, rest) = c.div_rem_euclid(&g);
    if !rest.is_zero() {
        return None;
    }
    let (a1, b1) = (a.div_exact(&g), b.div_exact(&g));
    // a1 u ≡ k (mod b1): u ≡ s k, as a1 s ≡ 1.
    let size = Integer::from(b1.unsigned_abs());
    let (_, mut u) = (&s * &k).div_rem_euclid(&size);
    if &u + &u > size {
        u = &u - &size;
    }
    let v = (&k - &(&a1 * &u)).div_exact(&b1);
    Some(vec![(u, vec![b1]), (v, vec![-a1])])
}

/// The answer to a linear system.
fn linear_answer(xs: &[Expr], solution: Vec<(Integer, Vec<Integer>)>) -> Expr {
    list(
        xs.iter()
            .zip(solution)
            .map(|(x, (constant, coefficients))| {
                let mut value = Expr::from(constant);
                for (k, c) in coefficients.into_iter().enumerate() {
                    value += Expr::from(c) * parameter(k);
                }
                Expr::call(Symbol::EQUATION, [x.clone(), value])
            })
            .collect(),
    )
}

/// `{x = a, y = b}`.
fn pair(x: &Expr, y: &Expr, a: Integer, b: Integer) -> Expr {
    list(vec![
        Expr::call(Symbol::EQUATION, [x.clone(), Expr::from(a)]),
        Expr::call(Symbol::EQUATION, [y.clone(), Expr::from(b)]),
    ])
}

/// The coefficients of `f(x, y)` in `x²`, `y²`, `x y`, `x`, `y` and its
/// constant, integers; `None` when it is no quadratic so.
fn quadratic(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    y: &Expr,
) -> Result<Option<[Integer; 6]>, EvalError> {
    let d = |f: &Expr, v: &Expr| Expr::call(Symbol::D, [f.clone(), v.clone()]);
    let zero = [(x.clone(), Expr::from(0)), (y.clone(), Expr::from(0))];
    let half = Rational::new(Integer::ONE, Integer::from(2));
    let mut values = Vec::new();
    for (e, scale) in [
        (d(&d(f, x), x), half.clone()),
        (d(&d(f, y), y), half),
        (d(&d(f, x), y), Rational::ONE),
        (d(f, x), Rational::ONE),
        (d(f, y), Rational::ONE),
        (f.clone(), Rational::ONE),
    ] {
        // The derivative taken before its variables are replaced.
        let v = ev.eval(&e)?;
        let v = ev.eval(&v.substitute(&zero))?;
        let Some(q) = v.to_rational() else {
            return Ok(None);
        };
        let scaled = &q * &scale;
        if !scaled.is_integer() {
            return Ok(None);
        }
        values.push(scaled.numer().clone());
    }
    let [a, b, c, dx, dy, k]: [Integer; 6] = values.try_into().expect("six");
    let rebuilt = Expr::from(a.clone()) * x.pow(2)
        + Expr::from(b.clone()) * y.pow(2)
        + Expr::from(c.clone()) * x.clone() * y.clone()
        + Expr::from(dx.clone()) * x.clone()
        + Expr::from(dy.clone()) * y.clone()
        + Expr::from(k.clone());
    let gap = ev.eval(&Expr::call(Symbol::EXPAND, [f.clone() - rebuilt]))?;
    Ok((gap == Expr::from(0)).then_some([a, b, c, dx, dy, k]))
}

/// The fundamental solution of `x² − D y² = 1`; `None` past [`MAX_PELL`]
/// quotients.
fn pell(d: &Integer) -> Option<(Integer, Integer)> {
    // √D's quotients over two periods at most: the convergent `h/k` before
    // a `Q` of 1 has `h² − D k² = ±1`, the end of a period.
    let (s, _) = d.sqrt_rem();
    let (mut p, mut q) = (Integer::ZERO, Integer::ONE);
    let (mut h, mut h0, mut k, mut k0) = (Integer::ONE, Integer::ZERO, Integer::ZERO, Integer::ONE);
    for _ in 0..MAX_PELL {
        let a = Rational::new(&p + &s, q.clone()).floor();
        let (next_h, next_k) = (&(&a * &h) + &h0, &(&a * &k) + &k0);
        (h0, h, k0, k) = (h, next_h, k, next_k);
        let next = &(&a * &q) - &p;
        q = (d - &(&next * &next)).div_exact(&q);
        p = next;
        if q == Integer::ONE && &(&h * &h) - &(&(d * &k) * &k) == Integer::ONE {
            return Some((h, k));
        }
    }
    None
}

/// `x y = n`: the pairs of divisors, both signs; `None` past
/// [`MAX_PRODUCT`].
fn products(x: &Expr, y: &Expr, c: &Integer, k: &Integer) -> Option<Expr> {
    let (n, rest) = (-k.clone()).div_rem_euclid(c);
    if !rest.is_zero() {
        return Some(list(Vec::new()));
    }
    let size = n
        .to_i64()
        .filter(|v| (1..=MAX_PRODUCT).contains(&v.abs()))?;
    // |n|'s divisors ascending: those to √|n|, then their cofactors.
    let m = size.abs();
    let small: Vec<i64> = (1..=m.isqrt()).filter(|d| m % d == 0).collect();
    let large = small.iter().rev().filter(|&&d| d * d != m).map(|d| m / d);
    let divisors: Vec<i64> = small.iter().copied().chain(large).collect();
    let mut found = Vec::new();
    for sign in [1i64, -1] {
        for d in &divisors {
            found.push(pair(
                x,
                y,
                Integer::from(sign * d),
                Integer::from(size / (sign * d)),
            ));
        }
    }
    Some(list(found))
}

/// `x² + y² = n`: every representation; `None` past [`MAX_SQUARES`].
fn two_squares(x: &Expr, y: &Expr, k: &Integer) -> Option<Expr> {
    let n = (-k.clone())
        .to_i64()
        .filter(|v| (0..=MAX_SQUARES).contains(v))?;
    let (root, _) = Integer::from(n).sqrt_rem();
    let root = root.to_i64()?;
    let mut found = Vec::new();
    for u in -root..=root {
        let (v, rest) = Integer::from(n - u * u).sqrt_rem();
        if rest.is_zero() {
            if v.is_zero() {
                found.push(pair(x, y, Integer::from(u), v));
            } else {
                found.push(pair(x, y, Integer::from(u), -v.clone()));
                found.push(pair(x, y, Integer::from(u), v));
            }
        }
    }
    Some(list(found))
}

/// `x² − D y² = 1`: giac's general solution; `None` for a square `D`.
fn pell_answer(x: &Expr, y: &Expr, d: &Integer) -> Option<Expr> {
    let (_, rest) = d.sqrt_rem();
    if rest.is_zero() {
        return None;
    }
    let (x1, y1) = pell(d)?;
    let root = Expr::from(d.clone()).pow(Expr::from(Rational::new(Integer::ONE, Integer::from(2))));
    let plus = (Expr::from(x1.clone()) + Expr::from(y1.clone()) * root.clone()).pow(parameter(0));
    let minus = (Expr::from(x1) - Expr::from(y1) * root.clone()).pow(parameter(0));
    Some(list(vec![
        Expr::call(
            Symbol::EQUATION,
            [x.clone(), (plus.clone() + minus.clone()) / Expr::from(2)],
        ),
        Expr::call(
            Symbol::EQUATION,
            [y.clone(), (plus - minus) / (Expr::from(2) * root)],
        ),
    ]))
}

/// The equations' sides and the variables of a call.
fn problem(e: &Expr) -> Option<(Vec<Expr>, Vec<Expr>)> {
    let (equations, given) = match e.args() {
        [eqs] => (eqs, None),
        [eqs, vars] => (eqs, Some(vars)),
        _ => return None,
    };
    let items = match head(equations) {
        Some(Symbol::LIST) => equations.args().to_vec(),
        _ => vec![equations.clone()],
    };
    let sides: Vec<Expr> = items
        .iter()
        .map(|eq| match (head(eq), eq.args()) {
            (Some(Symbol::EQUATION | Symbol::EQUAL), [l, r]) => l.clone() - r.clone(),
            _ => eq.clone(),
        })
        .collect();
    let mut xs = Vec::new();
    match given {
        Some(vars) if head(vars) == Some(Symbol::LIST) => xs = vars.args().to_vec(),
        Some(var) => xs.push(var.clone()),
        None => {
            for f in &sides {
                symbols_in_order(f, &mut xs);
            }
        }
    }
    let symbols = xs
        .iter()
        .all(|x| x.as_symbol().is_some_and(|s| !s.is_builtin()));
    (!xs.is_empty() && symbols).then_some((sides, xs))
}

/// `IntegerSolve[eqs]` and `IntegerSolve[eqs, vars]`.
pub(in crate::native) fn integer_solve(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((sides, xs)) = problem(e) else {
        return Ok(None);
    };
    // Linear.
    let mut rows = Vec::new();
    for f in &sides {
        let Some((a, k)) = linear(ev, f, &xs)? else {
            break;
        };
        rows.push((a, -k));
    }
    if rows.len() == sides.len() {
        let (a, b): (Vec<Vec<Integer>>, Vec<Integer>) = rows.into_iter().unzip();
        let solution = match (a.as_slice(), b.as_slice()) {
            ([row], [c]) if row.len() == 2 && !row[0].is_zero() && !row[1].is_zero() => {
                two_variables(&row[0], &row[1], c)
            }
            _ => linear_system(&a, &b).map(normalized),
        };
        let result = solution.map_or_else(|| list(Vec::new()), |s| linear_answer(&xs, s));
        ev.step("isolve", e, &result);
        return Ok(Some(result));
    }
    let ([f], [x, y]) = (sides.as_slice(), xs.as_slice()) else {
        return Ok(None);
    };
    let Some([a, b, c, dx, dy, k]) = quadratic(ev, f, x, y)? else {
        return Ok(None);
    };
    let one = Integer::ONE;
    let plain = dx.is_zero() && dy.is_zero();
    let result = if plain && a.is_zero() && b.is_zero() && !c.is_zero() {
        products(x, y, &c, &k)
    } else if plain && a == one && b == one && c.is_zero() {
        two_squares(x, y, &k)
    } else if plain && a == one && c.is_zero() && k == -one.clone() && b.is_negative() {
        pell_answer(x, y, &-b)
    } else {
        None
    };
    let Some(result) = result else {
        return Ok(None);
    };
    ev.step("isolve", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
