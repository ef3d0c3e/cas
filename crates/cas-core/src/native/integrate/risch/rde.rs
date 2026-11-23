//! The Risch differential equation `y′ + f y = g` (risch.cc:187-515):
//! `risch_desolve`, which bounds the denominator and the degree of `y`, and
//! `SPDE`, Rothstein's degree reduction [Bronstein, ch. 6].

use super::parts::{Acc, in_risch, risch_poly_part, rothstein_trager_resultant};
use super::tower::{Frac, Tower};
use super::{Step, derive, kernel, kernel_arg, together, tsimplify, within, zero};
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

/// The deepest recursion of [`spde_x`], the prototype's.
const SPDE_X_DEPTH: usize = 60;

/// giac's `rischde_simplify` (risch.cc:187): `R`, `S`, `T` over their gcd;
/// `None` when `R` and `S` keep a common factor in `T`.
fn rischde_simplify(c: &Tower, r: &Frac, s: &Frac, t: &Frac) -> Option<(Frac, Frac, Frac)> {
    let g = c.tgcd(&c.tgcd(r, s), t);
    let (r, s, t) = if g.num.is_zero() {
        (r.clone(), s.clone(), t.clone())
    } else {
        (c.exquo(r, &g), c.exquo(s, &g), c.exquo(t, &g))
    };
    if c.deg(&c.tgcd(&r, &s)) > 0 {
        return None;
    }
    Some((r, s, t))
}

/// giac's `spde_x` (risch.cc:197): `N` with `S N + R N′ = T`, `R` without
/// `T`, by repeated division.
fn spde_x(c: &Tower, s: &Frac, r: &Frac, t: &Frac, depth: usize) -> Option<Frac> {
    if t.num.is_zero() {
        return Some(c.zero());
    }
    if s.num.is_zero() || c.deg(t) < c.deg(s) || depth > SPDE_X_DEPTH {
        return None;
    }
    let q = c.quo(t, s);
    let new_t = c.sub(&c.sub(t, &c.mul(s, &q)), &c.mul(r, &c.d(&q)));
    let n = spde_x(c, s, r, &new_t, depth + 1)?;
    Some(c.add(&q, &n))
}

/// giac's `SPDE` (risch.cc:238): `y` of degree at most `ydeg` with
/// `R y′ + S y = T`, as an expression; `None` when there is none.
#[allow(
    clippy::too_many_arguments,
    reason = "giac's signature, with the evaluator"
)]
fn spde(
    ev: &mut Eval<'_>,
    c: &Tower,
    r: &Frac,
    s: &Frac,
    t: &Frac,
    v: &[Expr],
    ydeg: i64,
) -> Step<Option<Expr>> {
    if t.num.is_zero() {
        return Ok(Some(zero()));
    }
    if ydeg < 0 {
        return Ok(None);
    }
    within(ev, v.len(), c.deg(t), ydeg)?;
    let Some((r, s, t)) = rischde_simplify(c, r, s, t) else {
        return Ok(None);
    };
    let rd = c.deg(&r);
    if rd > 0 {
        // Tabcuv(S, R, T): S V + R U = T, deg V < deg R.
        let Some((vv, uu)) = c.abcuv(&s, &r, &t) else {
            return Ok(None);
        };
        let new_s = c.add(&s, &c.d(&r));
        let new_t = c.sub(&uu, &c.d(&vv));
        let ydeg = ydeg - i64::try_from(rd).unwrap_or(i64::MAX);
        let Some(prim) = spde(ev, c, &r, &new_s, &new_t, v, ydeg)? else {
            return Ok(None);
        };
        let re = c.r2e(ev, &r)?;
        let ve = c.r2e(ev, &vv)?;
        return Ok(Some(ev.eval(&(re * prim + ve))?));
    }
    let z_kernel = &v[0];
    let x = v.last().expect("x");
    if z_kernel == x || c.deg(&s) > 0 {
        let Some(n) = spde_x(c, &s, &r, &t, 0) else {
            return Ok(None);
        };
        return Ok(Some(c.r2e(ev, &n)?));
    }
    let v1 = &v[1..];
    let rr = c.lc(&r);
    let mut tc = Vec::new();
    for cf in c.coeffs(&t) {
        tc.push(c.r2e(ev, &c.div(&cf, &rr))?);
    }
    if s.num.is_zero() {
        // risch.cc:292 (kept): the lower tower, allowed_lnarg 1.
        let mut acc = Acc::empty();
        let ok = risch_poly_part(ev, &tc, 0, v1, &Expr::from(1), &mut acc)?;
        return Ok(ok.then_some(acc.prim));
    }
    let tdeg = tc.len() - 1;
    let Some(is_exp) = kernel(z_kernel) else {
        return Ok(None);
    };
    if usize::try_from(ydeg).unwrap_or(usize::MAX) < tdeg {
        return Ok(None);
    }
    let z = kernel_arg(z_kernel).clone();
    let dz = derive(ev, &z, x)?;
    let b = c.r2e(ev, &c.div(&c.lc(&s), &rr))?;
    let mut previous = zero();
    let mut prim = Vec::new();
    for k in (0..=tdeg).rev() {
        let kk = Expr::from(i64::try_from(k).unwrap_or(i64::MAX));
        let sol = if is_exp {
            let f = ev.eval(&(kk.clone() * dz.clone() + b.clone()))?;
            risch_desolve(ev, &f, &tc[tdeg - k], v1, false)?
        } else {
            let tmp = tc[tdeg - k].clone()
                - (kk.clone() + Expr::from(1)) * previous.clone() * dz.clone() / z.clone();
            let tmp = ev.eval(&tmp)?;
            risch_desolve(ev, &b, &tmp, v1, false)?
        };
        let Some(sol) = sol else {
            return Ok(None);
        };
        prim.push(sol.clone() * z_kernel.pow(kk));
        previous = sol;
    }
    Ok(Some(ev.eval(&Expr::call(Symbol::PLUS, prim))?))
}

/// giac's `risch_desolve` (risch.cc:323): `y` with `y′ + f y = g` in the
/// tower `v`; `None` when there is none in it. `f_is_derivative` skips the
/// bound on the denominator from the simple poles of `f`.
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
pub(super) fn risch_desolve(
    ev: &mut Eval<'_>,
    f: &Expr,
    g: &Expr,
    v: &[Expr],
    f_is_derivative: bool,
) -> Step<Option<Expr>> {
    if *g == zero() {
        return Ok(Some(zero()));
    }
    let z_kernel = &v[0];
    let x = v.last().expect("x").clone();
    let Some((c, fr)) = Tower::new(ev, v, &[f.clone(), g.clone()])? else {
        return Ok(None);
    };
    let v1 = &v[1..];
    let (ff, gg) = (&fr[0], &fr[1]);
    let fnum = c.poly(ff.num.clone());
    let fden = c.poly(ff.den.clone());
    let gnum = c.poly(gg.num.clone());
    let gden = c.poly(gg.den.clone());
    let mut fdenred = fden.clone();
    let mut gdenred = gden.clone();
    let mut dd = c.one();
    let kind = kernel(z_kernel);
    if kind == Some(true) {
        // The power of T in the denominator (risch.cc:345-383).
        let val = |p: &Frac| i64::try_from(c.val(p)).unwrap_or(i64::MAX);
        let (fdenv, fnumv, gdenv, gnumv) = (val(&fden), val(&fnum), val(&gden), val(&gnum));
        fdenred = c.shift_down(&fden, c.val(&fden));
        gdenred = c.shift_down(&gden, c.val(&gden));
        let mut alpha: i64 = 0;
        let beta = fnumv - fdenv;
        let gamma = gnumv - gdenv;
        let dzr_e = derive(ev, kernel_arg(z_kernel), &x)?;
        let Some(dzr) = c.e2r(ev, &dzr_e)? else {
            return Ok(None);
        };
        let lnc = |c: &Tower| {
            let f0 = c.div(&c.const_coeff(&fnum), &c.const_coeff(&fdenred));
            let q = c.neg(&c.div(&f0, &dzr));
            Tower::integer(&q).and_then(|n| i64::try_from(n).ok())
        };
        if gamma < 0 {
            alpha = gamma;
            if beta <= 0 {
                if beta != 0 {
                    alpha = 0.min(gamma - beta);
                } else if let Some(n) = lnc(&c) {
                    alpha = alpha.min(n);
                }
            }
        }
        if gamma >= 0
            && beta == 0
            && let Some(n) = lnc(&c)
        {
            alpha = alpha.min(n);
        }
        dd = c.top_power(u64::try_from(-alpha).unwrap_or(0));
    }
    if !f_is_derivative {
        // The simple poles of f (risch.cc:386-430).
        let tmpy = c.diff_top(&fdenred);
        let tmpw = fdenred.clone();
        let g1 = c.tgcd(&tmpy, &tmpw);
        let (tmpy, tmpw) = if g1.num.is_zero() {
            (tmpy, tmpw)
        } else {
            (c.exquo(&tmpy, &g1), c.exquo(&tmpw, &g1))
        };
        let tmpy = c.sub(&tmpy, &c.diff_top(&tmpw));
        let f1 = c.tgcd(&tmpw, &tmpy);
        if c.deg(&f1) > 0 {
            let f1cofact = c.quo(&fdenred, &f1);
            // Tabcuv(f1cofact, f1, fnum): N2 f1cofact ≡ fnum mod f1.
            if let Some((n2, _)) = c.abcuv(&f1cofact, &f1, &fnum) {
                let (pres, p1) = rothstein_trager_resultant(&c, &n2, &f1.num);
                for fac in c.factors(&pres, false) {
                    let Some(root) = c.linear_root(&fac) else {
                        continue;
                    };
                    let Some(residue) = Tower::integer(&root).and_then(|n| i64::try_from(n).ok())
                    else {
                        continue;
                    };
                    if residue < 0 {
                        let tmp1 = c.at_root(&p1, &root);
                        let p1subst = c.tgcd(&c.poly(tmp1.num), &f1);
                        let k = u64::try_from(-residue).unwrap_or(0);
                        dd = c.mul(&dd, &c.pow(&p1subst, k));
                    }
                }
            }
        }
    }
    let cg = c.tgcd(&fdenred, &gdenred);
    let num = c.tgcd(&gdenred, &c.diff_top(&gdenred));
    let den = c.tgcd(&cg, &c.diff_top(&cg));
    dd = c.mul(&dd, &c.quo(&num, &den));
    let d_d = c.d(&dd);
    let l = c.lcm(&fden, &gden);
    let mut r = c.mul(&dd, &l);
    let mut s = c.sub(&c.mul(&c.mul(ff, &dd), &l), &c.mul(&d_d, &l));
    let mut t = c.mul(&c.mul(gg, &c.mul(&dd, &dd)), &l);
    let Some((r1, s1, t1)) = rischde_simplify(&c, &r, &s, &t) else {
        return Ok(None);
    };
    (r, s, t) = (r1, s1, t1);
    let deg = |p: &Frac| i64::try_from(c.deg(p)).unwrap_or(i64::MAX);
    let (rd, sd, td) = (deg(&r), deg(&s), deg(&t));
    let rr = c.lc(&r);
    let ss = c.lc(&s);
    let mut ydeg = td - sd;
    if ydeg == 0 && !ss.num.is_zero() {
        // A constant ratio T/S (risch.cc:453-466).
        let tt = c.lc(&t);
        let chk = c.sub(&c.mul(&s, &tt), &c.mul(&t, &ss));
        if chk.num.is_zero() {
            let ratio = c.r2e(ev, &c.div(&tt, &ss))?;
            if !ratio.contains(&x) {
                let de = c.r2e(ev, &dd)?;
                return Ok(Some(together(ev, &(ratio / de))?));
            }
        }
    }
    let mut expshift = Expr::from(1);
    if *z_kernel == x {
        ydeg = td - (rd - 1).max(sd);
        if rd - 1 == sd
            && !rr.num.is_zero()
            && let Some(n) = Tower::integer(&c.div(&ss, &rr)).and_then(|n| i64::try_from(n).ok())
            && n > ydeg
        {
            ydeg = n;
        }
    } else {
        let Some(is_exp) = kind else {
            return Ok(None);
        };
        let z = kernel_arg(z_kernel).clone();
        ydeg = td - rd.max(sd);
        if is_exp {
            if rd == sd && !ss.num.is_zero() {
                // risch.cc:486 (kept): Rr/Ss.
                let tmp = c.r2e(ev, &c.div(&rr, &ss))?;
                let r = in_risch(ev, &tmp, v1, &z)?;
                if r.ok
                    && let Some(n) = r.lncoeff.integer()
                {
                    ydeg = ydeg.max(n);
                }
            }
        } else if rd == sd + 1 && !ss.num.is_zero() {
            let tmp = c.r2e(ev, &c.div(&rr, &ss))?;
            let res = in_risch(ev, &tmp, v1, &Expr::from(1))?;
            if res.ok {
                let prim = tsimplify(ev, &Expr::call(Symbol::EXP, [res.prim]))?;
                if let Some(es) = c.e2r(ev, &prim)? {
                    expshift = prim;
                    let en = c.poly(es.num);
                    let ed = c.poly(es.den);
                    r = c.mul(&c.mul(&r, &rr), &ed);
                    // risch.cc:503-504 (kept): S with the updated R.
                    s = c.mul(&c.sub(&c.mul(&s, &rr), &c.mul(&r, &ss)), &ed);
                    t = c.mul(&c.mul(&t, &rr), &en);
                }
            }
        }
    }
    let Some(y) = spde(ev, &c, &r, &s, &t, v, ydeg)? else {
        return Ok(None);
    };
    let de = c.r2e(ev, &dd)?;
    Ok(Some(together(ev, &(y * expshift / de))?))
}
