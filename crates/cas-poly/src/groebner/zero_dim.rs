//! Zero-dimensional ideals and syzygies (CAS-541): the normal set of a
//! reduced basis, the rational univariate representation of the points of a
//! zero-dimensional ideal of ℚ\[x₁, …, xₙ\] \[Rou99\], and the syzygies of a
//! list of polynomials by elimination, all as giac 2.0.0 gives them.
//! `doc/algorithms/groebner.md` describes the algorithms.

#[cfg(test)]
mod tests;

use cas_num::Rational;

use super::{
    Distributed, GroebnerError, GroebnerLimits, MonomialOrder, Term, Work, basis_with, divides,
    product, remainder,
};
use crate::cancel::Cancel;

use MonomialOrder::{DegRevLex, Lex};

/// The most monomials of a normal set [`rational_univariate`] works with:
/// past it, the points of the ideal are too many, and it stops with
/// [`GroebnerError::TooLarge`].
pub const MAX_NORMAL_SET: usize = 256;

/// The normal set of the ideal a reduced Gröbner `basis` under `order`
/// generates: the monomials no leading monomial of the basis divides, a
/// basis of the quotient ring as a vector space, in increasing order under
/// `order`. `None` when the ideal is not zero-dimensional, some variable
/// having no pure power among the leading monomials (and for the zero
/// ideal); empty for `[1]`.
///
/// The whole set is built: its size is at most the product over the
/// variables of the least exponent of their pure powers, which a caller
/// bounds when the basis comes from a user.
#[must_use]
pub fn normal_set(basis: &[Distributed], order: MonomialOrder) -> Option<Vec<Vec<u64>>> {
    standard(basis, order, usize::MAX).ok().flatten()
}

/// [`normal_set`], [`GroebnerError::TooLarge`] past `cap` monomials.
fn standard(
    basis: &[Distributed],
    order: MonomialOrder,
    cap: usize,
) -> Result<Option<Vec<Vec<u64>>>, GroebnerError> {
    let leads: Vec<&[u64]> = basis
        .iter()
        .filter_map(Distributed::leading_monomial)
        .collect();
    let Some(n) = leads.first().map(|m| m.len()) else {
        return Ok(None);
    };
    if leads.iter().any(|m| m.iter().all(|&e| e == 0)) {
        return Ok(Some(Vec::new()));
    }
    let pure =
        |i: usize, m: &[u64]| m[i] > 0 && m.iter().enumerate().all(|(j, &e)| j == i || e == 0);
    if !(0..n).all(|i| leads.iter().any(|m| pure(i, m))) {
        return Ok(None);
    }
    let mut out = Vec::new();
    walk(&leads, 0, &mut vec![0; n], &mut out, cap)?;
    out.sort_by(|a, b| order.cmp(a, b));
    Ok(Some(out))
}

/// Adds to `out` the monomials of the normal set that agree with
/// `monomial` before the variable `k`: the exponent of `x_k` grows from 0
/// until a leading monomial divides, which a pure power of `x_k` ends, and
/// every multiple of a divisible monomial is divisible.
fn walk(
    leads: &[&[u64]],
    k: usize,
    monomial: &mut [u64],
    out: &mut Vec<Vec<u64>>,
    cap: usize,
) -> Result<(), GroebnerError> {
    if k == monomial.len() {
        if out.len() >= cap {
            return Err(GroebnerError::TooLarge);
        }
        out.push(monomial.to_vec());
        return Ok(());
    }
    while !leads.iter().any(|g| divides(g, monomial)) {
        walk(leads, k + 1, monomial, out, cap)?;
        monomial[k] += 1;
    }
    monomial[k] = 0;
    Ok(())
}

/// The rational univariate representation of the points of a
/// zero-dimensional ideal \[Rou99\]: a separating linear form `l`, and
/// polynomials in one variable `t`, each a vector of coefficients from
/// degree 0, without trailing zeros (0 is empty), such that the points are
/// `xᵢ = gᵢ(l)/f'(l)` at the roots `l` of `f`, each once.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Rur {
    /// The coefficient of each variable in the separating form `l`.
    pub form: Vec<Rational>,
    /// `f`: minus the monic minimal polynomial of `l` modulo the radical,
    /// so its leading coefficient is −1, as giac gives it.
    pub minimal: Vec<Rational>,
    /// `f'`.
    pub derivative: Vec<Rational>,
    /// `gᵢ` for each variable, of degree below that of `f`.
    pub numerators: Vec<Vec<Rational>>,
}

/// The rational univariate representation of the ideal `polys` generate
/// (CAS-541), the polynomials in `n` variables sorted under
/// [`MonomialOrder::DegRevLex`]; `None` when the ideal is the whole ring,
/// not zero-dimensional, or `n` is 1 (giac then gives the basis).
///
/// With `G` the reduced `DegRevLex` basis and `D` the size of its normal
/// set: the minimal polynomial `mᵢ` of each variable in ℚ\[X\]/I, from the
/// normal forms of its powers; the last variable with `deg mᵢ = D`, if any,
/// is `l`. The radical `I' = I + (sqfree(mᵢ)(xᵢ))` \[Sei74\], of `D'`
/// points; then the single variables from the last, then `x₁ + k·x₂ + … +
/// k^(n−1)·xₙ` for `k = 1, 2, …, n·D'²`, until a minimal polynomial modulo
/// `I'` has degree `D'`. Each `xᵢ = hᵢ(l)` modulo `I'`, and
/// `gᵢ = hᵢ·f' rem f`.
///
/// # Errors
///
/// [`GroebnerError::Cancelled`] when `cancel` says so, and
/// [`GroebnerError::TooLarge`] past the bounds of a basis computation, a
/// normal set past [`MAX_NORMAL_SET`], or past the `limits`.
pub fn rational_univariate(
    polys: &[Distributed],
    limits: GroebnerLimits,
    cancel: &dyn Cancel,
) -> Result<Option<Rur>, GroebnerError> {
    let work = Work::new(limits, cancel);
    let Some(n) = polys
        .iter()
        .find_map(Distributed::leading_monomial)
        .map(<[u64]>::len)
    else {
        return Ok(None);
    };
    if n < 2 {
        return Ok(None);
    }
    let Some(ideal) = Quotient::new(polys, &work)? else {
        return Ok(None);
    };
    let points = ideal.set.len();
    let mut minimals = vec![Vec::new(); n];
    let mut chosen = None;
    for i in (0..n).rev() {
        let (m, _) = ideal.minimal(&linear(&unit(n, i)), &work)?;
        if chosen.is_none() && m.len() == points + 1 {
            chosen = Some(i);
        }
        minimals[i] = m;
    }
    // The radical: each variable's minimal polynomial made squarefree
    // added: Seidenberg's lemma [Sei74], in characteristic 0.
    let mut generators = ideal.basis.clone();
    for (i, m) in minimals.iter().enumerate() {
        let terms = squarefree(m)
            .into_iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(k, c)| {
                let mut e = vec![0; n];
                e[i] = u64::try_from(k).unwrap_or(u64::MAX);
                (e, c)
            })
            .collect();
        generators.push(Distributed::new(terms, DegRevLex));
    }
    let Some(radical) = Quotient::new(&generators, &work)? else {
        return Ok(None);
    };
    let points = radical.set.len();
    let singles = chosen.into_iter().chain((0..n).rev()).map(|i| unit(n, i));
    let tries = u64::try_from(n * points * points).unwrap_or(u64::MAX);
    let combinations = (1..=tries).map(|k| {
        let k = Rational::from(k);
        let mut c = Rational::ONE;
        (0..n)
            .map(|_| {
                let this = c.clone();
                c = &c * &k;
                this
            })
            .collect()
    });
    for form in singles.chain(combinations) {
        let (minimal, rows) = radical.minimal(&linear(&form), &work)?;
        if minimal.len() != points + 1 {
            continue;
        }
        let f: Vec<Rational> = minimal.iter().map(|c| -c).collect();
        let derivative = derivative(&f);
        let mut numerators = Vec::with_capacity(n);
        for i in 0..n {
            let x = radical.normal_form(&linear(&unit(n, i)), &work)?;
            let h = express(&rows, radical.coordinates(&x), &work)?;
            numerators.push(div_rem(&mul(&h, &derivative), &f).1);
        }
        return Ok(Some(Rur {
            form,
            minimal: f,
            derivative,
            numerators,
        }));
    }
    // Unreached: of the `n·D'²` forms, at most `(n − 1)·D'(D' − 1)/2`
    // fail to separate the `D'` points.
    Ok(None)
}

/// The generators of the syzygies of `polys` (CAS-541): for `f₁, …, fₘ` in
/// `n` variables, sorted under [`MonomialOrder::DegRevLex`], the reduced
/// `DegRevLex` basis, in the `m + n` variables `A₁, …, Aₘ, x₁, …, xₙ`, of
/// the ideal of the `Σ aᵢAᵢ` with `Σ aᵢfᵢ = 0` and of every `AᵢAⱼ`: its
/// elements linear in the `Aᵢ` generate the module of syzygies, as giac's
/// `syzygy` gives them. `n` is 0 when every polynomial is 0.
///
/// By elimination \[KR00\]: one more variable `E`, the generators
/// `fᵢE + Aᵢ` and every monomial of degree 2 in `E, A₁, …, Aₘ`, their
/// `Lex` basis with `E` the greatest; its elements free of `E`, their basis
/// under `DegRevLex`.
///
/// # Errors
///
/// [`GroebnerError::Cancelled`] when `cancel` says so, and
/// [`GroebnerError::TooLarge`] past the bounds of a basis computation or
/// the `limits`.
pub fn syzygies(
    polys: &[Distributed],
    limits: GroebnerLimits,
    cancel: &dyn Cancel,
) -> Result<Vec<Distributed>, GroebnerError> {
    let work = Work::new(limits, cancel);
    let m = polys.len();
    let n = polys
        .iter()
        .find_map(Distributed::leading_monomial)
        .map_or(0, <[u64]>::len);
    let width = 1 + m + n;
    let mut generators = Vec::new();
    for (i, f) in polys.iter().enumerate() {
        let mut terms: Vec<Term> = f
            .terms()
            .iter()
            .map(|(x, c)| {
                let mut e = vec![0; width];
                e[0] = 1;
                e[1 + m..].copy_from_slice(x);
                (e, c.clone())
            })
            .collect();
        let mut a = vec![0; width];
        a[1 + i] = 1;
        terms.push((a, Rational::ONE));
        generators.push(Distributed::new(terms, Lex));
    }
    for i in 0..=m {
        for j in i..=m {
            let mut e = vec![0; width];
            e[i] += 1;
            e[j] += 1;
            generators.push(Distributed::new(vec![(e, Rational::ONE)], Lex));
        }
    }
    // Under Lex with E the greatest, an element whose leading monomial is
    // free of E is free of E.
    let free: Vec<Distributed> = basis_with(&generators, Lex, &work)?
        .iter()
        .filter(|g| g.leading_monomial().is_some_and(|l| l[0] == 0))
        .map(|g| {
            let terms = g.terms().iter().map(|(e, c)| (e[1..].to_vec(), c.clone()));
            Distributed::new(terms.collect(), DegRevLex)
        })
        .collect();
    basis_with(&free, DegRevLex, &work)
}

/// ℚ\[X\]/I: the reduced `DegRevLex` basis of `I` and its normal set.
struct Quotient {
    basis: Vec<Distributed>,
    set: Vec<Vec<u64>>,
}

/// A row of the echelon form of the normal forms of the powers of a form:
/// its coordinates over the normal set, 1 at `pivot` and 0 at the earlier
/// rows' pivots, and the polynomial in the form whose normal form it is.
struct Row {
    pivot: usize,
    coordinates: Vec<Rational>,
    polynomial: Vec<Rational>,
}

impl Quotient {
    /// The quotient by the ideal `polys` generate; `None` when it is the
    /// whole ring or not zero-dimensional.
    fn new(polys: &[Distributed], work: &Work<'_>) -> Result<Option<Quotient>, GroebnerError> {
        let basis = basis_with(polys, DegRevLex, work)?;
        match standard(&basis, DegRevLex, MAX_NORMAL_SET)? {
            Some(set) if !set.is_empty() => Ok(Some(Quotient { basis, set })),
            _ => Ok(None),
        }
    }

    fn normal_form(&self, p: &Distributed, work: &Work<'_>) -> Result<Distributed, GroebnerError> {
        let divisors: Vec<&Distributed> = self.basis.iter().collect();
        remainder(p, &divisors, DegRevLex, work)
    }

    /// The coordinates of a normal form over the normal set.
    fn coordinates(&self, p: &Distributed) -> Vec<Rational> {
        let mut v = vec![Rational::ZERO; self.set.len()];
        for (m, c) in p.terms() {
            let k = self
                .set
                .binary_search_by(|s| DegRevLex.cmp(s, m))
                .expect("a normal form's monomials are in the normal set");
            v[k] = c.clone();
        }
        v
    }

    /// The monic minimal polynomial of `form` in the quotient, and the rows
    /// of the normal forms of its powers below its degree: each power's
    /// normal form, the form times the last one's, reduced by the rows
    /// until it is 0 (incremental Gaussian elimination over ℚ).
    fn minimal(
        &self,
        form: &Distributed,
        work: &Work<'_>,
    ) -> Result<(Vec<Rational>, Vec<Row>), GroebnerError> {
        let one = vec![0; self.set[0].len()];
        let mut power = Distributed::new(vec![(one, Rational::ONE)], DegRevLex);
        let mut rows: Vec<Row> = Vec::new();
        for j in 0..=self.set.len() {
            work.step()?;
            let mut v = self.coordinates(&power);
            let mut p = vec![Rational::ZERO; j + 1];
            p[j] = Rational::ONE;
            eliminate(&rows, &mut v, &mut p, work)?;
            let Some(pivot) = v.iter().position(|c| !c.is_zero()) else {
                return Ok((trim(p), rows));
            };
            let inverse = v[pivot].recip();
            for c in v.iter_mut().chain(p.iter_mut()) {
                *c *= &inverse;
                work.made(c)?;
            }
            rows.push(Row {
                pivot,
                coordinates: v,
                polynomial: p,
            });
            power = self.normal_form(&times(form, &power, work)?, work)?;
        }
        unreachable!("D + 1 vectors in a space of dimension D are dependent")
    }
}

/// `v` and `p` less the multiples of the rows that clear `v` at their
/// pivots, in order: a later row is 0 at an earlier one's pivot.
fn eliminate(
    rows: &[Row],
    v: &mut [Rational],
    p: &mut Vec<Rational>,
    work: &Work<'_>,
) -> Result<(), GroebnerError> {
    for row in rows {
        let c = v[row.pivot].clone();
        if c.is_zero() {
            continue;
        }
        for (x, y) in v.iter_mut().zip(&row.coordinates) {
            if !y.is_zero() {
                *x -= &c * y;
                work.made(x)?;
            }
        }
        if p.len() < row.polynomial.len() {
            p.resize(row.polynomial.len(), Rational::ZERO);
        }
        for (x, y) in p.iter_mut().zip(&row.polynomial) {
            *x -= &c * y;
            work.made(x)?;
        }
    }
    Ok(())
}

/// The polynomial `h` in the form with `h(form) = v`, for the coordinates
/// `v` of a normal form in the span of the rows.
fn express(
    rows: &[Row],
    mut v: Vec<Rational>,
    work: &Work<'_>,
) -> Result<Vec<Rational>, GroebnerError> {
    let mut p = Vec::new();
    eliminate(rows, &mut v, &mut p, work)?;
    Ok(trim(p.into_iter().map(|c| -c).collect()))
}

/// The variable `xᵢ` of `n`, as the coefficients of a form.
fn unit(n: usize, i: usize) -> Vec<Rational> {
    let mut form = vec![Rational::ZERO; n];
    form[i] = Rational::ONE;
    form
}

/// The linear form of coefficients `form`.
fn linear(form: &[Rational]) -> Distributed {
    let terms = form
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, c)| {
            let mut e = vec![0; form.len()];
            e[i] = 1;
            (e, c.clone())
        })
        .collect();
    Distributed::new(terms, DegRevLex)
}

/// `a·b` under `DegRevLex`.
fn times(a: &Distributed, b: &Distributed, work: &Work<'_>) -> Result<Distributed, GroebnerError> {
    work.spend(a.terms().len() * b.terms().len())?;
    let terms = a
        .terms()
        .iter()
        .flat_map(|(m, c)| b.terms().iter().map(move |(t, d)| (product(m, t), c * d)))
        .collect();
    Ok(Distributed::new(terms, DegRevLex))
}

// Polynomials in one variable over ℚ: coefficients from degree 0, without
// trailing zeros.

fn trim(mut a: Vec<Rational>) -> Vec<Rational> {
    while a.last().is_some_and(Rational::is_zero) {
        a.pop();
    }
    a
}

fn derivative(a: &[Rational]) -> Vec<Rational> {
    let terms = a.iter().enumerate().skip(1);
    trim(terms.map(|(k, c)| c * &Rational::from(k)).collect())
}

fn mul(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![Rational::ZERO; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

/// The quotient and remainder of `a` by `b`, not 0.
fn div_rem(a: &[Rational], b: &[Rational]) -> (Vec<Rational>, Vec<Rational>) {
    let inverse = b.last().expect("a non-zero divisor").recip();
    let mut q = vec![Rational::ZERO; a.len().saturating_sub(b.len()) + 1];
    let mut r = trim(a.to_vec());
    while r.len() >= b.len() {
        let shift = r.len() - b.len();
        let c = r.last().expect("not shorter than b") * &inverse;
        for (j, d) in b.iter().enumerate() {
            r[shift + j] -= &c * d;
        }
        q[shift] = c;
        r = trim(r);
    }
    (trim(q), r)
}

fn monic(a: Vec<Rational>) -> Vec<Rational> {
    match a.last() {
        Some(c) => {
            let inverse = c.recip();
            a.iter().map(|d| d * &inverse).collect()
        }
        None => a,
    }
}

fn gcd(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let (mut a, mut b) = (trim(a.to_vec()), trim(b.to_vec()));
    while !b.is_empty() {
        let r = div_rem(&a, &b).1;
        a = std::mem::replace(&mut b, r);
    }
    monic(a)
}

/// `m / gcd(m, m')`, monic: the product of the distinct irreducible
/// factors of `m`, not 0.
fn squarefree(m: &[Rational]) -> Vec<Rational> {
    monic(div_rem(m, &gcd(m, &derivative(m))).0)
}
