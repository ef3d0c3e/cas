//! The differential field of a tower (risch.cc's `e2r`, `r2e` and `diff`,
//! :99): rational functions of the generators `lv` of the expressions, the
//! tower's kernels among them, as numerators and denominators over ℤ\[i\]
//! in cas-poly's sparse form, through cas-core's [`Polynomials`]. Elements
//! of `K[T]`, `K` the field of the lower generators and `T` the top of the
//! tower, are fractions whose denominator has no `T`; their quotients,
//! remainders, gcds and extended gcds in `T` are the pseudo-divisions and
//! cofactors of cas-poly, as giac's `TPseudoDivRem`, `Tabcuv` and `Tpartfrac`
//! (poly.h) on its own polynomials in `lv`.
//!
//! The coefficients are Gaussian integers in every case (giac's are when `i`
//! appears): the gcds, quotients and resultants are the same over any
//! extension, and the Rothstein–Trager resultant is factored over ℚ(i) or ℚ
//! as giac's flags say ([`Tower::factors`]). The ring has one variable more
//! than the generators, the `t` of the resultant, which never reaches an
//! expression.

use cas_num::algebra::{GaussianIntegers, Integers, Ring};
use cas_num::{Complex, Integer};
use cas_poly::{Cancel, Sparse, SparseRing};

use super::{Step, Stop};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::poly::{Bounds, Fraction, PolyError, Polynomials};

/// A polynomial of the tower's ring.
pub(super) type P = Sparse<Complex>;

/// `num / den`, without a common factor.
#[derive(Clone, Debug)]
pub(super) struct Frac {
    pub(super) num: P,
    pub(super) den: P,
}

/// The tower `v` (its top first, `x` last) and the generators of the
/// expressions it was made for.
pub(super) struct Tower {
    polys: Polynomials,
    pub(super) ring: SparseRing<GaussianIntegers>,
    /// The variable of each kernel.
    vars: Vec<usize>,
    /// The derivative of each kernel.
    dv: Vec<Frac>,
    /// The variable of the Rothstein–Trager resultant.
    pub(super) t: usize,
}

/// `exprs` as rational functions of their generators, `None` when they have
/// none; past the limits, [`Stop::Limit`].
fn convert(ev: &mut Eval<'_>, exprs: &[Expr]) -> Step<Option<Polynomials>> {
    let ctx = ev.context();
    let limits = ctx.limits();
    let bounds = Bounds {
        terms: limits.terms,
        bits: limits.number_bits,
    };
    let cancel: &dyn Cancel = ctx.interrupt_flag();
    match Polynomials::new(exprs, bounds, cancel) {
        Ok(p) => Ok(Some(p)),
        Err(PolyError::Interrupted) => Err(Stop::Eval(crate::context::EvalError::Interrupted)),
        Err(PolyError::TooLarge | PolyError::Bits) => Err(Stop::Limit),
        Err(_) => Ok(None),
    }
}

/// The generators of `e` (giac's `lvar`), `None` when it has no polynomial
/// form.
pub(super) fn generators(ev: &mut Eval<'_>, e: &Expr) -> Step<Option<Vec<Expr>>> {
    Ok(convert(ev, std::slice::from_ref(e))?.map(|p| p.generators().to_vec()))
}

impl Tower {
    /// The tower `v` of the expressions `exprs`, with their rational
    /// functions; `None` when a kernel is not a generator of them (giac's
    /// `lvar(e, l)` adds the others, risch.cc:750).
    pub(super) fn new(
        ev: &mut Eval<'_>,
        v: &[Expr],
        exprs: &[Expr],
    ) -> Step<Option<(Tower, Vec<Frac>)>> {
        let x = v.last().expect("a tower ends with x").clone();
        let mut list = exprs.to_vec();
        list.extend(v.iter().cloned());
        for k in v {
            list.push(super::derive(ev, k, &x)?);
        }
        let Some(polys) = convert(ev, &list)? else {
            return Ok(None);
        };
        let n = polys.generators().len();
        let ring = SparseRing::new(GaussianIntegers, n + 1);
        let Some(vars) = v
            .iter()
            .map(|k| polys.variable(k))
            .collect::<Option<Vec<usize>>>()
        else {
            return Ok(None);
        };
        let identity: Vec<usize> = (0..n).collect();
        let mut fracs: Vec<Frac> = polys
            .fractions()
            .iter()
            .map(|f| embed(&ring, f, &identity))
            .collect();
        let dv = fracs.split_off(exprs.len() + v.len());
        fracs.truncate(exprs.len());
        let tower = Tower {
            polys,
            ring,
            vars,
            dv,
            t: n,
        };
        Ok(Some((tower, fracs)))
    }

    /// The variable of the top of the tower.
    pub(super) fn top(&self) -> usize {
        self.vars[0]
    }

    /// `e` in the tower's generators; `None` when it has another.
    pub(super) fn e2r(&self, ev: &mut Eval<'_>, e: &Expr) -> Step<Option<Frac>> {
        let gens = self.polys.generators();
        let mut list = vec![e.clone()];
        list.extend(gens.iter().cloned());
        let Some(p) = convert(ev, &list)? else {
            return Ok(None);
        };
        let Some(map) = p
            .generators()
            .iter()
            .map(|g| gens.iter().position(|h| h == g))
            .collect::<Option<Vec<usize>>>()
        else {
            return Ok(None);
        };
        Ok(Some(embed(&self.ring, &p.fractions()[0], &map)))
    }

    /// The expression of `f`, evaluated.
    pub(super) fn r2e(&self, ev: &mut Eval<'_>, f: &Frac) -> Step<Expr> {
        debug_assert!(self.free_of(&f.num, self.t) && self.free_of(&f.den, self.t));
        Ok(ev.eval(&self.polys.ratio_expr(&f.num, &f.den))?)
    }

    /// The expression of the polynomial `p`, evaluated.
    pub(super) fn p2e(&self, ev: &mut Eval<'_>, p: &P) -> Step<Expr> {
        self.r2e(ev, &self.poly(p.clone()))
    }

    // ---------------------------------------------------------------- fractions

    /// `p / 1`.
    pub(super) fn poly(&self, p: P) -> Frac {
        Frac {
            num: p,
            den: self.ring.one(),
        }
    }

    pub(super) fn zero(&self) -> Frac {
        self.poly(self.ring.zero())
    }

    pub(super) fn one(&self) -> Frac {
        self.poly(self.ring.one())
    }

    /// `num / den` without their gcd; `den` non-zero.
    #[allow(
        clippy::needless_pass_by_value,
        reason = "callers build the parts for it"
    )]
    pub(super) fn frac(&self, num: P, den: P) -> Frac {
        if num.is_zero() {
            return self.zero();
        }
        let (_, num, den) = self.ring.gcd_cofactors(&num, &den);
        Frac { num, den }
    }

    pub(super) fn add(&self, a: &Frac, b: &Frac) -> Frac {
        if a.num.is_zero() {
            return b.clone();
        }
        if b.num.is_zero() {
            return a.clone();
        }
        let r = &self.ring;
        if a.den == b.den {
            return self.frac(r.add(&a.num, &b.num), a.den.clone());
        }
        let num = r.add(&r.mul(&a.num, &b.den), &r.mul(&b.num, &a.den));
        self.frac(num, r.mul(&a.den, &b.den))
    }

    pub(super) fn neg(&self, a: &Frac) -> Frac {
        Frac {
            num: self.ring.neg(&a.num),
            den: a.den.clone(),
        }
    }

    pub(super) fn sub(&self, a: &Frac, b: &Frac) -> Frac {
        self.add(a, &self.neg(b))
    }

    pub(super) fn mul(&self, a: &Frac, b: &Frac) -> Frac {
        let r = &self.ring;
        self.frac(r.mul(&a.num, &b.num), r.mul(&a.den, &b.den))
    }

    /// `a / b`, `b` non-zero.
    pub(super) fn div(&self, a: &Frac, b: &Frac) -> Frac {
        let r = &self.ring;
        self.frac(r.mul(&a.num, &b.den), r.mul(&a.den, &b.num))
    }

    /// `a · n` for an integer `n`.
    pub(super) fn scale(&self, a: &Frac, n: i64) -> Frac {
        self.mul(a, &self.poly(self.ring.constant(Complex::from(n))))
    }

    pub(super) fn pow(&self, a: &Frac, n: u64) -> Frac {
        Frac {
            num: self.ring.pow(&a.num, n),
            den: self.ring.pow(&a.den, n),
        }
    }

    /// `T^k`.
    pub(super) fn top_power(&self, k: u64) -> Frac {
        let mut e = vec![0; self.t + 1];
        e[self.top()] = k;
        self.poly(self.ring.monomial(Complex::from(1), &e))
    }

    /// Whether `p` has no `x_var`.
    pub(super) fn free_of(&self, p: &P, var: usize) -> bool {
        self.ring.degree(p, var).is_none_or(|d| d == 0)
    }

    /// Whether `f` is an integer, and which.
    pub(super) fn integer(f: &Frac) -> Option<Integer> {
        let n = constant(&f.num)?;
        let d = constant(&f.den)?;
        if !n.is_real() || !d.is_real() {
            return None;
        }
        let q = n.re() / d.re();
        q.is_integer().then(|| q.numer().clone())
    }

    // ---------------------------------------------------------------- the derivation

    /// The derivation of the tower on a polynomial: `Σ ∂p/∂vᵢ · vᵢ′`; the
    /// other generators are constants (risch.cc:99).
    pub(super) fn d_poly(&self, p: &P) -> Frac {
        let mut acc = self.zero();
        for (var, dv) in self.vars.iter().zip(&self.dv) {
            let dp = self.ring.derivative(p, *var);
            if dp.is_zero() || dv.num.is_zero() {
                continue;
            }
            let term = Frac {
                num: self.ring.mul(&dp, &dv.num),
                den: dv.den.clone(),
            };
            acc = self.add(&acc, &term);
        }
        acc
    }

    /// The derivation on a fraction, `(n′ d − n d′) / d²`.
    pub(super) fn d(&self, f: &Frac) -> Frac {
        let dn = self.d_poly(&f.num);
        if self.free_of_tower(&f.den) {
            return self.div(&dn, &self.poly(f.den.clone()));
        }
        let dd = self.d_poly(&f.den);
        let n = self.mul(&dn, &self.poly(f.den.clone()));
        let m = self.mul(&dd, &self.poly(f.num.clone()));
        let r = &self.ring;
        self.div(&self.sub(&n, &m), &self.poly(r.mul(&f.den, &f.den)))
    }

    /// Whether `p` has no variable of the tower.
    fn free_of_tower(&self, p: &P) -> bool {
        self.vars.iter().all(|&v| self.free_of(p, v))
    }

    // ---------------------------------------------------------------- K[T]

    /// The degree in `T`, 0 for zero (sympy's `degree` in the prototype).
    pub(super) fn deg(&self, f: &Frac) -> u64 {
        self.ring.degree(&f.num, self.top()).unwrap_or(0)
    }

    /// The leading coefficient in `T`, 0 for zero.
    pub(super) fn lc(&self, f: &Frac) -> Frac {
        self.frac(self.ring.lead_coeff(&f.num, self.top()), f.den.clone())
    }

    /// The coefficients in `T`, the highest power first; none for zero.
    pub(super) fn coeffs(&self, f: &Frac) -> Vec<Frac> {
        let Some(d) = self.ring.degree(&f.num, self.top()) else {
            return Vec::new();
        };
        (0..=d)
            .rev()
            .map(|k| self.frac(self.ring.coeff(&f.num, self.top(), k), f.den.clone()))
            .collect()
    }

    /// The coefficient of `T^0`.
    pub(super) fn const_coeff(&self, f: &Frac) -> Frac {
        self.frac(self.ring.coeff(&f.num, self.top(), 0), f.den.clone())
    }

    /// The valuation in `T`, 0 for zero.
    pub(super) fn val(&self, f: &Frac) -> u64 {
        self.ring.valuation(&f.num, self.top()).unwrap_or(0)
    }

    /// `f / T^k`, `T^k` dividing it.
    pub(super) fn shift_down(&self, f: &Frac, k: u64) -> Frac {
        if k == 0 {
            return f.clone();
        }
        let tk = self.top_power(k).num;
        Frac {
            num: self.ring.div_exact(&f.num, &tk).expect("T^k divides"),
            den: f.den.clone(),
        }
    }

    /// The quotient and remainder in `T` over `K`: `l^δ a.num = q b.num +
    /// r` gives `a / b = q b.den / (a.den l^δ) + r / (a.den l^δ b)`.
    pub(super) fn quo_rem(&self, a: &Frac, b: &Frac) -> (Frac, Frac) {
        let r = &self.ring;
        let top = self.top();
        let delta = match (r.degree(&a.num, top), r.degree(&b.num, top)) {
            (Some(m), Some(n)) if m >= n => m - n + 1,
            _ => 0,
        };
        let l = r.pow(&r.lead_coeff(&b.num, top), delta);
        let (q, rem) = r.pseudo_div_rem(&a.num, &b.num, top);
        let d = r.mul(&a.den, &l);
        (self.frac(r.mul(&q, &b.den), d.clone()), self.frac(rem, d))
    }

    pub(super) fn quo(&self, a: &Frac, b: &Frac) -> Frac {
        self.quo_rem(a, b).0
    }

    /// The gcd in `T` over `K`, primitive in `T` (the prototype's monic
    /// `tgcd`): zero for two zeros.
    pub(super) fn tgcd(&self, a: &Frac, b: &Frac) -> Frac {
        let g = self.ring.gcd(&a.num, &b.num);
        if g.is_zero() {
            return self.zero();
        }
        self.poly(self.ring.primitive_part_in(&g, self.top()))
    }

    /// `a / g` for `g` from [`Tower::tgcd`], which divides it.
    pub(super) fn exquo(&self, a: &Frac, g: &Frac) -> Frac {
        let num = self.ring.div_exact(&a.num, &g.num).expect("a gcd divides");
        self.frac(num, a.den.clone())
    }

    /// The formal derivative in `T`.
    pub(super) fn diff_top(&self, f: &Frac) -> Frac {
        self.frac(self.ring.derivative(&f.num, self.top()), f.den.clone())
    }

    /// The least common multiple in `T` of two polynomials, up to a factor
    /// without `T`.
    pub(super) fn lcm(&self, a: &Frac, b: &Frac) -> Frac {
        self.poly(self.ring.lcm(&a.num, &b.num))
    }

    /// `(u, v)` with `u a + v b = c`, `deg u < deg b` in `T` (giac's
    /// `Tabcuv`); `None` when the gcd of `a` and `b` does not divide `c`.
    pub(super) fn abcuv(&self, a: &Frac, b: &Frac, c: &Frac) -> Option<(Frac, Frac)> {
        let r = &self.ring;
        // `u a.num + v b.num = d c.num`, so `(u a.den) a + (v b.den) b = d c.den c`.
        let (u, v, d) = r.abcuv(&a.num, &b.num, &c.num, self.top())?;
        let dd = r.mul(&d, &c.den);
        Some((
            self.frac(r.mul(&u, &a.den), dd.clone()),
            self.frac(r.mul(&v, &b.den), dd),
        ))
    }

    /// The square-free factors in `T` of `den`, primitive in `T`, one per
    /// multiplicity; with `zero_apart`, the power of `T` that divides it is
    /// a last factor of its own (risch.cc:774-783).
    pub(super) fn sqff(&self, den: &P, zero_apart: bool) -> Vec<(P, u64)> {
        let r = &self.ring;
        let top = self.top();
        let zm = if zero_apart {
            r.valuation(den, top).unwrap_or(0)
        } else {
            0
        };
        let dens = if zm > 0 {
            self.shift_down(&self.poly(den.clone()), zm).num
        } else {
            den.clone()
        };
        let mut by_multiplicity: std::collections::BTreeMap<u64, P> =
            std::collections::BTreeMap::default();
        for (f, m) in r.squarefree(&dens).factors() {
            if self.free_of(f, top) {
                continue;
            }
            let f = r.primitive_part_in(f, top);
            let p = by_multiplicity.entry(*m).or_insert_with(|| r.one());
            *p = r.mul(p, &f);
        }
        let mut out: Vec<(P, u64)> = by_multiplicity.into_iter().map(|(m, f)| (f, m)).collect();
        if zm > 0 {
            out.push((self.top_power(1).num, zm));
        }
        out
    }

    /// The partial fractions of `num / den` over `factors` from
    /// [`Tower::sqff`] (giac's `Tpartfrac`): each factor's numerator over
    /// its power, and the polynomial part; `None` when cas-poly finds the
    /// factors no decomposition of `den` (never for [`Tower::sqff`]'s).
    pub(super) fn partfrac(
        &self,
        num: &P,
        den: &P,
        factors: &[(P, u64)],
    ) -> Option<(Vec<Frac>, Frac)> {
        let r = &self.ring;
        let pf = r.partial_fractions(num, den, factors, self.top())?;
        let pieces = pf
            .pieces
            .into_iter()
            .map(|(n, d)| self.frac(n, d))
            .collect();
        Some((pieces, self.frac(pf.quotient, pf.quotient_den)))
    }

    // ---------------------------------------------------------------- the resultant's variable

    /// `t`, the variable of the resultant.
    pub(super) fn t_var(&self) -> P {
        self.ring.var(self.t)
    }

    /// The factors of `p` that have `t`, over ℚ(i) when `complex` or a
    /// coefficient is not real, over ℚ otherwise (giac's `factor`).
    pub(super) fn factors(&self, p: &P, complex: bool) -> Vec<P> {
        let r = &self.ring;
        if p.is_zero() || self.free_of(p, self.t) {
            return Vec::new();
        }
        let real = p.coeffs().iter().all(Complex::is_real);
        let factors: Vec<P> = if complex || !real {
            r.factor(p)
                .factors()
                .iter()
                .map(|(f, _)| f.clone())
                .collect()
        } else {
            let z = SparseRing::new(Integers, self.t + 1);
            let pz = z.from_terms(p.terms().map(|(e, c)| (e, c.re().numer().clone())));
            z.factor(&pz)
                .factors()
                .iter()
                .map(|(f, _)| r.from_terms(f.terms().map(|(e, c)| (e, Complex::from(c.clone())))))
                .collect()
        };
        factors
            .into_iter()
            .filter(|f| !self.free_of(f, self.t))
            .collect()
    }

    /// The root of a factor `a₁ t + a₀` of degree 1 in `t`, `−a₀ / a₁`;
    /// `None` for another degree.
    pub(super) fn linear_root(&self, f: &P) -> Option<Frac> {
        let r = &self.ring;
        if r.degree(f, self.t) != Some(1) {
            return None;
        }
        let a1 = r.coeff(f, self.t, 1);
        let a0 = r.coeff(f, self.t, 0);
        Some(self.frac(r.neg(&a0), a1))
    }

    /// `p(t = root)` for a `p` of degree at most 1 in `t`, a fraction
    /// without `t`.
    pub(super) fn at_root(&self, p: &P, root: &Frac) -> Frac {
        let r = &self.ring;
        let p0 = r.coeff(p, self.t, 0);
        let p1 = r.coeff(p, self.t, 1);
        debug_assert!(r.degree(p, self.t).is_none_or(|d| d <= 1));
        let num = r.add(&r.mul(&p0, &root.den), &r.mul(&p1, &root.num));
        self.frac(num, root.den.clone())
    }

    /// The leading coefficient of `p` in lexicographic order with the
    /// variables `order` first, then the others in theirs.
    pub(super) fn lex_lead<'p>(&self, p: &'p P, order: &[usize]) -> Option<&'p Complex> {
        let mut rank: Vec<usize> = order.to_vec();
        rank.extend((0..=self.t).filter(|v| !order.contains(v)));
        p.terms()
            .max_by(|(a, _), (b, _)| {
                let ka: Vec<u64> = rank.iter().map(|&v| a[v]).collect();
                let kb: Vec<u64> = rank.iter().map(|&v| b[v]).collect();
                ka.cmp(&kb)
            })
            .map(|(_, c)| c)
    }

    /// The variables of the tower, the top first, then the other
    /// generators.
    pub(super) fn var_order(&self) -> Vec<usize> {
        let mut order = self.vars.clone();
        order.extend((0..self.t).filter(|v| !self.vars.contains(v)));
        order
    }

    /// The variables of the tower, the top first.
    pub(super) fn vars(&self) -> &[usize] {
        &self.vars
    }
}

/// The constant of a polynomial without variables.
fn constant(p: &P) -> Option<&Complex> {
    match p.term_count() {
        1 if p.terms().all(|(e, _)| e.iter().all(|&k| k == 0)) => p.coeffs().first(),
        _ => None,
    }
}

/// A conversion's fraction in the tower's ring, its variable `k` the
/// tower's `map[k]`.
fn embed(ring: &SparseRing<GaussianIntegers>, f: &Fraction, map: &[usize]) -> Frac {
    let n = ring.nvars();
    let remap = |e: Vec<u64>| {
        let mut out = vec![0; n];
        for (k, x) in e.into_iter().enumerate() {
            out[map[k]] = x;
        }
        out
    };
    let (num, den) = match f {
        Fraction::Integers { num, den } => (
            ring.from_terms(
                num.terms()
                    .map(|(e, c)| (remap(e), Complex::from(c.clone()))),
            ),
            ring.from_terms(
                den.terms()
                    .map(|(e, c)| (remap(e), Complex::from(c.clone()))),
            ),
        ),
        Fraction::GaussianIntegers { num, den } => (
            ring.from_terms(num.terms().map(|(e, c)| (remap(e), c.clone()))),
            ring.from_terms(den.terms().map(|(e, c)| (remap(e), c.clone()))),
        ),
    };
    Frac { num, den }
}
