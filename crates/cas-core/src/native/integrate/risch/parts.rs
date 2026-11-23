//! The integration of a rational function of the tower (risch.cc:517-876):
//! `in_risch`, its polynomial part `risch_poly_part`, Hermite reduction, the
//! Rothstein–Trager resultant and `integrate_ei`.

use cas_num::algebra::Ring;
use cas_num::{Complex, Integer};

use super::tower::{Frac, P, Tower};
use super::{Step, call, derive, has_i, kernel, kernel_arg, note, rde, together, within, zero};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// giac's `lncoeff`: a value, or `inf` where giac divides by `derive(1)`
/// (risch.cc:495, :292; T9.c1 departure 6).
#[derive(Clone, Debug)]
pub(super) enum Lnc {
    Value(Expr),
    Infinite,
}

impl Lnc {
    /// The value when it is an integer.
    pub(super) fn integer(&self) -> Option<i64> {
        match self {
            Lnc::Value(e) => e.as_i64(),
            Lnc::Infinite => None,
        }
    }
}

/// What `in_risch` gives: whether it integrated, the primitive, the
/// coefficient of the allowed logarithm and what remains.
#[derive(Clone, Debug)]
pub(super) struct Acc {
    pub(super) ok: bool,
    pub(super) prim: Expr,
    pub(super) lncoeff: Lnc,
    pub(super) remains: Expr,
}

impl Acc {
    pub(super) fn empty() -> Acc {
        Acc {
            ok: true,
            prim: zero(),
            lncoeff: Lnc::Value(zero()),
            remains: zero(),
        }
    }

    fn add_prim(&mut self, e: Expr) {
        self.prim = std::mem::replace(&mut self.prim, zero()) + e;
    }

    fn add_remains(&mut self, e: Expr) {
        self.remains = std::mem::replace(&mut self.remains, zero()) + e;
    }

    fn failed(mut self) -> Acc {
        self.ok = false;
        self
    }
}

/// `p` over `f^m`, as an expression.
fn over_power(ev: &mut Eval<'_>, c: &Tower, n: &Frac, f: &P, m: u64) -> Step<Expr> {
    let d = c.pow(&c.poly(f.clone()), m);
    c.r2e(ev, &c.div(n, &d))
}

/// The coefficients of `f` in `T`, the highest first, as expressions.
fn coeff_exprs(ev: &mut Eval<'_>, c: &Tower, f: &Frac) -> Step<Vec<Expr>> {
    let mut out = Vec::new();
    for cf in c.coeffs(f) {
        out.push(c.r2e(ev, &cf)?);
    }
    Ok(out)
}

/// giac's `in_risch` (risch.cc:743): `∫ e` in the tower `v`, the only new
/// logarithm allowed `ln(allowed)` when `allowed` is not 0.
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
pub(super) fn in_risch(ev: &mut Eval<'_>, e: &Expr, v: &[Expr], allowed: &Expr) -> Step<Acc> {
    let mut acc = Acc::empty();
    let Some((c, fr)) = Tower::new(ev, v, std::slice::from_ref(e))? else {
        acc.remains = e.clone();
        return Ok(acc.failed());
    };
    let r = &fr[0];
    let top = c.top();
    let degree = c
        .ring
        .degree(&r.num, top)
        .unwrap_or(0)
        .max(c.ring.degree(&r.den, top).unwrap_or(0));
    within(ev, v.len(), degree, 1)?;
    let restricted = *allowed != zero();
    let vden = c.sqff(&r.den, true);
    let Some((pfs, ip)) = c.partfrac(&r.num, &r.den, &vden) else {
        acc.remains = e.clone();
        return Ok(acc.failed());
    };
    let is_exp = kernel(&v[0]) == Some(true);
    let mut ipshift = c.zero();
    for (ni, (f, m)) in pfs.iter().zip(&vden) {
        let m = *m;
        if is_exp && f.term_count() == 1 {
            // The generalised polynomial part ni / T^m (risch.cc:784).
            let tmp = coeff_exprs(ev, &c, ni)?;
            let shift = -i64::try_from(m).unwrap_or(i64::MAX);
            if !risch_poly_part(ev, &tmp, shift, v, allowed, &mut acc)? && restricted {
                return Ok(acc.failed());
            }
            continue;
        }
        let mut tnum = ni.clone();
        if m > 1 {
            let Some((t, padd)) = hermite_reduce(ev, &c, ni, f, m, &c.zero(), &Expr::from(1))?
            else {
                let rest = over_power(ev, &c, ni, f, m)?;
                acc.add_remains(rest);
                acc.ok = false;
                continue;
            };
            tnum = t;
            acc.add_prim(padd);
        }
        if tnum.num.is_zero() {
            continue;
        }
        if restricted {
            // Only ln(allowed) may appear (risch.cc:811-826).
            if pfs.len() > 1 {
                return Ok(acc.failed());
            }
            let dal = derive(ev, allowed, v.last().expect("x"))?;
            let dal = together(ev, &dal)?;
            if dal == zero() {
                acc.lncoeff = Lnc::Infinite;
                continue;
            }
            let part = over_power(ev, &c, ni, f, m)?;
            let lnc = together(ev, &(allowed.clone() * part / dal))?;
            acc.lncoeff = Lnc::Value(lnc.clone());
            let x = v.last().expect("x");
            if {
                let d = derive(ev, &lnc, x)?;
                together(ev, &d)?
            } != zero()
            {
                return Ok(acc.failed());
            }
            continue;
        }
        let (pres, p1) = rothstein_trager_resultant(&c, &tnum, f);
        let mut add_prim = zero();
        let mut ok = true;
        let x = v.last().expect("x").clone();
        for fac in c.factors(&pres, true) {
            let Some(root) = c.linear_root(&fac) else {
                ok = false;
                break;
            };
            let root_e = c.r2e(ev, &root)?;
            // giac derive(root, x): the true derivative, foreign kernels too.
            if {
                let d = derive(ev, &root_e, &x)?;
                together(ev, &d)?
            } != zero()
            {
                ok = false;
                break;
            }
            let tmp1 = c.at_root(&p1, &root);
            let targ = giac_gcd_norm(&c, &tmp1.num, f);
            let targ_e = c.p2e(ev, &targ)?;
            add_prim += root_e * Expr::call(Symbol::LOG, [targ_e]);
            let tp = c.poly(targ);
            let q = c.quo(&c.d(&tp), &tp);
            ipshift = c.sub(&ipshift, &c.mul(&root, &q));
        }
        if ok {
            acc.add_prim(add_prim);
        } else {
            // risch.cc:837, fixed: the Hermite-reduced fraction remains.
            let rest = c.r2e(ev, &c.div(&tnum, &c.poly(f.clone())))?;
            acc.add_remains(rest);
            acc.ok = false;
        }
    }
    let mut tmp = coeff_exprs(ev, &c, &ip)?;
    let shift_e = c.r2e(ev, &ipshift)?;
    match tmp.last_mut() {
        None => tmp.push(shift_e),
        Some(last) => *last = ev.eval(&(last.clone() + shift_e))?,
    }
    if !risch_poly_part(ev, &tmp, 0, v, allowed, &mut acc)? && restricted {
        return Ok(acc.failed());
    }
    let ok = acc.ok || !restricted;
    acc.ok = ok;
    Ok(acc)
}

/// `Σ e[i] Xk^(s − 1 − i)`, built.
fn symb_horner(e: &[Expr], xk: &Expr) -> Expr {
    let s = e.len();
    let terms: Vec<Expr> = e
        .iter()
        .enumerate()
        .map(|(i, c)| c.clone() * xk.pow(Expr::from(i64::try_from(s - 1 - i).unwrap_or(0))))
        .collect();
    Expr::call(Symbol::PLUS, terms)
}

/// giac's `risch_poly_part` (risch.cc:662): `∫ Σ e[i] T^(s − 1 − i + shift)`
/// for the coefficients `e` in the lower tower, the highest first.
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
pub(super) fn risch_poly_part(
    ev: &mut Eval<'_>,
    e: &[Expr],
    shift: i64,
    v: &[Expr],
    allowed: &Expr,
    acc: &mut Acc,
) -> Step<bool> {
    let s = e.len();
    let x = v.last().expect("x").clone();
    if v.len() == 1 {
        let terms: Vec<Expr> = e
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let k = i64::try_from(s - i).unwrap_or(1);
                c.clone() * x.pow(Expr::from(k)) / Expr::from(k)
            })
            .collect();
        acc.add_prim(Expr::call(Symbol::PLUS, terms));
        return Ok(true);
    }
    let v1 = &v[1..];
    let xk = &v[0];
    match kernel(xk) {
        Some(false) => {
            // risch.cc:668-695: the coefficients from the top, each with the
            // log of the argument allowed, whose coefficient moves up.
            if s == 0 {
                return Ok(true);
            }
            let arg = kernel_arg(xk);
            let dx = {
                let d = derive(ev, xk, &x)?;
                together(ev, &d)?
            };
            let mut eprim = vec![zero(); s + 1];
            for j in (1..s).rev() {
                let k = i64::try_from(j + 1).unwrap_or(i64::MAX);
                let g = ev.eval(
                    &(e[s - 1 - j].clone() - Expr::from(k) * eprim[s - 1 - j].clone() * dx.clone()),
                )?;
                let r = in_risch(ev, &g, v1, &arg)?;
                eprim[s - j] = r.prim;
                if !r.ok {
                    acc.add_remains(symb_horner(e, xk));
                    return Ok(false);
                }
                if let Lnc::Value(lnc) = &r.lncoeff {
                    eprim[s - 1 - j] = eprim[s - 1 - j].clone() + lnc.clone() / Expr::from(k);
                }
            }
            let g = ev.eval(&(e[s - 1].clone() - eprim[s - 1].clone() * dx))?;
            // risch.cc:693 (kept): allowed_lnarg 0.
            let r = in_risch(ev, &g, v1, &zero())?;
            acc.add_prim(r.prim + symb_horner(&eprim, xk));
            acc.add_remains(r.remains);
            if !r.ok {
                acc.ok = false;
            }
            Ok(r.ok)
        }
        Some(true) => {
            let y = kernel_arg(xk);
            let dy = derive(ev, &y, &x)?;
            let mut eprim = vec![zero(); s];
            for j in (0..s).rev() {
                let idx = s - 1 - j;
                let pui = i64::try_from(j).unwrap_or(i64::MAX) + shift;
                if pui == 0 {
                    let r = in_risch(ev, &e[idx], v1, allowed)?;
                    acc.lncoeff = r.lncoeff;
                    acc.add_prim(r.prim);
                    acc.add_remains(r.remains);
                    if !r.ok {
                        acc.ok = false;
                        if *allowed != zero() {
                            return Ok(false);
                        }
                    }
                    continue;
                }
                let f = ev.eval(&(Expr::from(pui) * dy.clone()))?;
                if let Some(sol) = rde::risch_desolve(ev, &f, &e[idx], v1, true)? {
                    eprim[idx] = sol;
                    continue;
                }
                let coeff = e[idx].clone();
                if *allowed == zero()
                    && let Some((a, b)) = is_linear_wrt(ev, &y, &x)?
                    && rational_in_x(ev, &coeff, &x)?
                {
                    let pa = ev.eval(&(Expr::from(pui) * a))?;
                    if let Some((p, rem)) = integrate_ei(ev, &pa, &coeff, &x)? {
                        let shifted = Expr::call(Symbol::EXP, [Expr::from(pui) * b]);
                        acc.add_prim(shifted.clone() * p);
                        acc.add_remains(shifted * rem);
                        continue;
                    }
                    note(|| format!("integrate_ei: no answer for {coeff:?} exp({pa:?} x)"));
                }
                note(|| {
                    format!(
                        "risch_poly_part: y′ + {f:?} y = {coeff:?} unsolved; {coeff:?} {xk:?}^{pui} remains"
                    )
                });
                acc.add_remains(coeff * xk.pow(Expr::from(pui)));
                acc.ok = false;
                eprim[idx] = zero();
                if *allowed != zero() {
                    return Ok(false);
                }
            }
            acc.add_prim(symb_horner(&eprim, xk) * xk.pow(Expr::from(shift)));
            Ok(true)
        }
        None => Ok(false),
    }
}

/// Whether the kernels of `e` with `x` are `x` alone (giac's
/// `lvarx(e) == [x]`).
fn rational_in_x(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Step<bool> {
    let Some(gens) = super::tower::generators(ev, e)? else {
        return Ok(false);
    };
    let with_x: Vec<&Expr> = gens.iter().filter(|g| g.contains(x)).collect();
    Ok(with_x.len() == 1 && with_x[0] == x)
}

/// giac's `is_linear_wrt` (intg.cc:196): `(a, b)` with `e = a x + b`.
fn is_linear_wrt(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Step<Option<(Expr, Expr)>> {
    let a = {
        let d = derive(ev, e, x)?;
        together(ev, &d)?
    };
    if a.contains(x) {
        return Ok(None);
    }
    let b = together(ev, &(e.clone() - a.clone() * x.clone()))?;
    if b.contains(x) {
        return Ok(None);
    }
    Ok(Some((a, b)))
}

/// giac's `hermite_reduce` (risch.cc:517): `num / fact^mult · exp(a x)` as
/// `num′ / fact · exp(a x)` plus an integrated part; `None` when `fact′` and
/// `fact` are not coprime.
pub(super) fn hermite_reduce(
    ev: &mut Eval<'_>,
    c: &Tower,
    num: &Frac,
    fact: &P,
    mult: u64,
    a: &Frac,
    expax: &Expr,
) -> Step<Option<(Frac, Expr)>> {
    if mult == 1 {
        return Ok(Some((num.clone(), zero())));
    }
    let f = c.poly(fact.clone());
    let fprime = c.d_poly(fact);
    let mut num = num.clone();
    let mut mult = mult;
    let mut resexpr = c.zero();
    while mult > 1 {
        // Tegcdtoabcuv: fprime v + fact u = num, deg v < deg fact.
        let Some((v, u)) = c.abcuv(&fprime, &f, &num) else {
            return Ok(None);
        };
        mult -= 1;
        let k = i64::try_from(mult).unwrap_or(i64::MAX);
        // diff(v·exp(a x))/exp(a x) / mult (risch.cc:164).
        let dv = c.add(&c.d(&v), &c.mul(a, &v));
        let numtemp = c.div(&dv, &c.scale(&c.one(), k));
        num = c.add(&numtemp, &u);
        let part = c.div(&v, &c.scale(&c.pow(&f, mult), k));
        resexpr = c.sub(&resexpr, &part);
    }
    let r = c.r2e(ev, &resexpr)?;
    Ok(Some((num, r * expax.clone())))
}

/// giac's `rothstein_trager_resultant` (risch.cc:124): `p₁ = num − t den′`
/// (its numerator, the sign making its leading coefficient in `[T, t, …]`
/// positive, :140) and `res_T(p₁, den)`.
pub(super) fn rothstein_trager_resultant(c: &Tower, num: &Frac, den: &P) -> (P, P) {
    let r = &c.ring;
    let dp = c.d_poly(den);
    let l = r.lcm(&num.den, &dp.den);
    let a = r.mul(&num.num, &r.div_exact(&l, &num.den).expect("an lcm"));
    let b = r.mul(
        &r.mul(&c.t_var(), &dp.num),
        &r.div_exact(&l, &dp.den).expect("an lcm"),
    );
    let mut p1 = r.sub(&a, &b);
    let mut order = vec![c.top(), c.t];
    order.extend(c.var_order().into_iter().skip(1));
    if c.lex_lead(&p1, &order).is_some_and(negative_real) {
        p1 = r.neg(&p1);
    }
    let pres = if c.ring.degree(den, c.top()).is_some_and(|d| d > 0) {
        r.resultant(&p1, den, c.top())
    } else {
        r.one()
    };
    (pres, p1)
}

fn negative_real(z: &Complex) -> bool {
    z.is_real() && z.re().is_negative()
}

/// giac's `_gcd` of a log argument (risch.cc:852), as probed (T9.c1 notes):
/// monomial factors shifted out; `a` itself when it divides `b` (as is
/// when both are multivariate and `b/a` is a monomial), else the gcd;
/// primitive over ℤ for real coefficients; a gcd with complex ones made
/// monic, then free of denominators (the prototype's `ext_norm`); negated
/// when its leading coefficient, the tower's top first, is a negative real.
fn giac_gcd_norm(c: &Tower, a: &P, b: &P) -> P {
    let r = &c.ring;
    let ext = a.coeffs().iter().chain(b.coeffs()).any(|z| !z.is_real());
    let order = c.var_order();
    let flip = |p: P| {
        if c.lex_lead(&p, &order).is_some_and(negative_real) {
            r.neg(&p)
        } else {
            p
        }
    };
    // Over ℤ: the gcd of the integer parts, positive (sympy's `primitive`).
    let primitive = |p: &P| {
        let mut g = Integer::from(0_u64);
        for z in p.coeffs() {
            g = g.gcd(z.re().numer()).gcd(z.im().numer());
        }
        if g.is_zero() || g.is_one() {
            return p.clone();
        }
        let g = Complex::from(g);
        r.from_terms(p.terms().map(|(e, z)| (e, z / &g)))
    };
    // Monic over ℚ(i), denominators cleared: `p · conj(lc)`, primitive.
    let monic = |p: &P| match c.lex_lead(p, &order) {
        Some(l) if !l.is_real() => {
            let l = l.conj();
            primitive(&r.from_terms(p.terms().map(|(e, z)| (e, z * &l))))
        }
        _ => primitive(p),
    };
    let used = (0..c.t)
        .filter(|&v| !c.free_of(a, v) || !c.free_of(b, v))
        .count();
    if used == 0 {
        return r.one();
    }
    // giac's `_gcd` drops the integer content of `a` with complex
    // coefficients too, and tests `a | b` over ℚ(i): `ln(x² + 1)` has
    // `p₁(i) = −2 − 2ix` and giac's `ln(−1 − ix)` (CAS-585 T9.c6; the
    // prototype kept `a` whole, and ℤ[i] division sent it to the monic gcd
    // `x − i`). Over ℚ(i), `a | b` exactly when `a | b · N(lc a)` over ℤ[i]:
    // the Gaussian content of `a` divides the norm of its leading
    // coefficient (Gauss's lemma).
    let pa = without_monomial(r, &primitive(a));
    let pb = without_monomial(r, b);
    let scaled = match c.lex_lead(&pa, &order) {
        Some(l) if ext => r.mul(&pb, &r.constant(l * &l.conj())),
        _ => pb.clone(),
    };
    // `a = 0` when `p₁` vanishes at the root (any `f′/f`): giac's `_gcd(0, b)`
    // is `b`, the prototype's `pb.div(pa)` raises and falls back to the gcd;
    // dividing by `pa` there was T9.c3's division by zero (CAS-585).
    if !pa.is_zero()
        && let Some(q) = r.div_exact(&scaled, &pa)
    {
        if used >= 2 && q.term_count() == 1 {
            return pa;
        }
        return flip(pa);
    }
    let g = r.gcd(a, b);
    flip(if ext { monic(&g) } else { primitive(&g) })
}

/// `p` over the greatest monomial dividing it (gausspol.cc:5214).
fn without_monomial(r: &cas_poly::SparseRing<cas_num::algebra::GaussianIntegers>, p: &P) -> P {
    let n = r.nvars();
    let mut low = vec![u64::MAX; n];
    for (e, _) in p.terms() {
        for (l, x) in low.iter_mut().zip(&e) {
            *l = (*l).min(*x);
        }
    }
    if low.iter().all(|&l| l == 0 || l == u64::MAX) {
        return p.clone();
    }
    r.from_terms(p.terms().map(|(e, z)| {
        let e = e.iter().zip(&low).map(|(x, l)| x - l).collect();
        (e, z.clone())
    }))
}

/// giac's `integrate_ei` (risch.cc:582): `∫ coeff · exp(a x)` for `coeff`
/// rational in `x`, as `(prim, remains)` both times `exp(a x)` where giac
/// writes them so: the polynomial part in closed form, Hermite reduction,
/// then `Ei`, or `Ci` and `Si` for an imaginary `a` and a real root, at the
/// roots of each factor over ℚ(i).
pub(super) fn integrate_ei(
    ev: &mut Eval<'_>,
    a: &Expr,
    coeff: &Expr,
    x: &Expr,
) -> Step<Option<(Expr, Expr)>> {
    let expax = ev.eval(&Expr::call(Symbol::EXP, [a.clone() * x.clone()]))?;
    // re(a) = 0: a = ima i, ima a real number.
    let ima = match a.to_complex() {
        Some(z) if z.re().is_zero() && !z.im().is_zero() => Some(Expr::from(z.im().clone())),
        _ => None,
    };
    let v = [x.clone()];
    let Some((c, fr)) = Tower::new(ev, &v, &[coeff.clone(), a.clone()])? else {
        return Ok(None);
    };
    let Some(af) = c.e2r(ev, a)? else {
        return Ok(None);
    };
    let r = &fr[0];
    let vden = c.sqff(&r.den, false);
    let Some((pfs, ip)) = c.partfrac(&r.num, &r.den, &vden) else {
        return Ok(None);
    };
    let mut prim = zero();
    if !ip.num.is_zero() {
        let p = c.r2e(ev, &ip)?;
        prim = integrate_poly_exp(ev, &p, a, x, &expax)?;
    }
    let mut remains = zero();
    for (ni, (f, m)) in pfs.iter().zip(&vden) {
        let m = *m;
        let mut tnum = ni.clone();
        if m > 1 {
            let Some((t, padd)) = hermite_reduce(ev, &c, ni, f, m, &af, &expax)? else {
                return Ok(None);
            };
            tnum = t;
            prim += padd;
        }
        if tnum.num.is_zero() {
            continue;
        }
        let itden = c.diff_top(&c.poly(f.clone()));
        for fac in c.ring.factor(f).factors().iter().map(|(g, _)| g) {
            if c.free_of(fac, c.top()) {
                continue;
            }
            if c.ring.degree(fac, c.top()) != Some(1) {
                // risch.cc:640 (kept): the original fraction remains.
                let part = over_power(ev, &c, ni, f, m)?;
                note(|| {
                    format!("integrate_ei: a non-linear factor, {part:?} exp({a:?} x) remains")
                });
                remains += expax.clone() * part;
                break;
            }
            let a1 = c.ring.coeff(fac, c.top(), 1);
            let a0 = c.ring.coeff(fac, c.top(), 0);
            let root = c.frac(c.ring.neg(&a0), a1);
            let root_e = c.r2e(ev, &root)?;
            let rootn = at(ev, &c, &tnum, &root_e, x)?;
            let rootd = at(ev, &c, &itden, &root_e, x)?;
            let rr = together(ev, &(rootn / rootd))?;
            let shifted = x.clone() - root_e.clone();
            match &ima {
                Some(ima) if !has_i(&root_e) => {
                    let negative = ima.to_rational().is_some_and(|q| q.is_negative());
                    let aa = if negative { -ima.clone() } else { ima.clone() };
                    let s = if negative { -super::i() } else { super::i() };
                    let phase = Expr::call(Symbol::COS, [ima.clone() * root_e.clone()])
                        + super::i() * Expr::call(Symbol::SIN, [ima.clone() * root_e.clone()]);
                    let ci = Expr::call(Symbol::COS_INTEGRAL, [aa.clone() * shifted.clone()]);
                    let si = Expr::call(Symbol::SIN_INTEGRAL, [aa * shifted]);
                    prim += phase * rr * (ci + s * si);
                }
                _ => {
                    let ei = Expr::call(Symbol::EXP_INTEGRAL_EI, [a.clone() * shifted]);
                    let e = Expr::call(Symbol::EXP, [a.clone() * root_e]);
                    prim += rr * e * ei;
                }
            }
        }
    }
    Ok(Some((ev.eval(&prim)?, ev.eval(&remains)?)))
}

/// `f(x = root)`, evaluated.
fn at(ev: &mut Eval<'_>, c: &Tower, f: &Frac, root: &Expr, x: &Expr) -> Step<Expr> {
    let e = c.r2e(ev, f)?;
    Ok(ev.eval(&e.substitute(&[(x.clone(), root.clone())]))?)
}

/// `∫ p exp(a x) = exp(a x) Σ (−1)^k p^(k) / a^(k + 1)` for a polynomial
/// `p` (giac's `_integrate` in `integrate_ei`; the same value).
fn integrate_poly_exp(ev: &mut Eval<'_>, p: &Expr, a: &Expr, x: &Expr, expax: &Expr) -> Step<Expr> {
    let mut res = Vec::new();
    let mut q = p.clone();
    let mut k: i64 = 0;
    while q != zero() {
        let sign = if k % 2 == 0 { 1 } else { -1 };
        res.push(Expr::from(sign) * q.clone() / a.pow(Expr::from(k + 1)));
        q = call(ev, Symbol::D, [q, x.clone()])?;
        k += 1;
    }
    Ok(Expr::call(Symbol::PLUS, res) * expax.clone())
}
