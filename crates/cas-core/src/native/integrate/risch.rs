//! `Risch` (CAS-585), Xcas's `risch`: the Risch algorithm for elementary
//! integrands, ported from giac 2.0.0.21's risch.cc function for function,
//! after the checked sympy prototype of T9.c1 (Bronstein, Symbolic
//! Integration I, 2nd ed., Springer 2005, chapters 5 and 6: the
//! integration of transcendental elementary functions and the Risch
//! differential equation).
//!
//! - `Risch[f, x]`, and `Risch[f]` in `x` as giac's `_risch` (:1018); other
//!   arguments leave the call. The answer is `prim + Risch[rem, x]`, `rem`
//!   what the algorithm could not integrate (giac prints
//!   `integrate(rem, x)`), or the call stays when nothing integrates.
//! - [`risch`] (:951): trigonometric functions become exponentials
//!   (`TrigToExp`); [`risch_lin`] (:878) tries the whole integrand, then
//!   each term of its [`lin`] form (`coeff·exp(expo)`, lin.cc:134) and the
//!   sum of what they leave. Fractional powers of `x` leave the integrand;
//!   `atan(g)` is integrated as `i/2 ln((i + g)/(i − g))` and written back;
//!   the constant part of an exponent is split off (`exp(2x + 1) = exp(1)
//!   exp(2x)`).
//! - The tower (risch.cc:56-97, [`risch_tower`]): the kernels of the
//!   integrand that have `x`, recursively (giac's `rlvarx`), sorted by size
//!   so that the largest is the top; powers with `x` become `exp(b ln a)`;
//!   `TranscendentalSimplify` (giac's `tsimplify`) writes exponentials with
//!   as few generators as it can, and the smaller valid tower wins. Kernels
//!   of one size keep our canonical order, where giac keeps the order of
//!   its own printout (T9.c1 notes, tower ties): rows that tie differ in
//!   shape or split, departures, not bugs.
//! - [`parts::in_risch`] (:743): the integrand as a rational function of
//!   the top `T` over the lower field, its denominator's power of `T` apart
//!   (risch.cc:774); square-free factorisation and partial fractions in `T`;
//!   Hermite reduction (:517) of the multiple factors; the logarithmic part
//!   by the Rothstein–Trager resultant (:124) factored over ℚ(i), a root
//!   free of `x` giving `r ln gcd(p₁(r), f)`; the polynomial part by
//!   [`parts::risch_poly_part`] (:662), which solves the Risch differential
//!   equation `y′ + f y = g` ([`rde::risch_desolve`], :323, with
//!   [`rde::spde`], :238, Rothstein's SPDE) coefficient by coefficient in
//!   an exponential and integrates the coefficients from the top in a
//!   logarithm; `exp(a x) · rational` falls back on `Ei`, `Ci` and `Si`
//!   ([`parts::integrate_ei`], :582).
//! - Polynomials are giac's: in all of `lv`, the tower's variables and the
//!   parameters, through cas-core's `Polynomials`, over ℤ\[i\]
//!   ([`tower`]); quotients in `T` by pseudo-division.
//! - giac's flaws kept (T9.c1 notes): risch.cc:693 passes 0 for the
//!   allowed logarithm of the last coefficient of a logarithmic polynomial
//!   (the fix breaks `ln(ln(x))/x`, T9.c2); :292 passes the lower tower to
//!   `risch_poly_part` for `y′ = T`; :486 uses `Rr/Ss`; :503 the updated
//!   `R`. Fixed: :837 leaves the Hermite-reduced fraction, not the original
//!   one, when a residue is not a constant (giac counted the Hermite part
//!   twice).
//! - Departures: the ln-argument unit over ℚ(i) (giac `ln(i x − 1)`, ours
//!   `ln(x + i)`, a constant apart); resultant factors of degree 2 are not
//!   split by a square root, so `1/(x² − 2)`, `1/(x³ + 1)` and `1/(x⁴ + 1)`
//!   leave the call; `tsimplify` is `TranscendentalSimplify`, without
//!   giac's merging of logarithms.
//! - Safety: an answer is returned only when `D[prim] + rem − f` is 0 by
//!   `Together` (or `TranscendentalSimplify`); otherwise the call stays. The
//!   remainder's own `Risch` integrates nothing: a remainder that still
//!   integrates is folded in first, so evaluation never loops.
//! - Bounds (C0019): the tower's height times the degree in its top and the
//!   degree bound of the differential equation must not exceed
//!   `Limits::terms`; past it the call stays, with a warning, as
//!   `Integrate` does at a limit.
#![allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    reason = "the names of giac's code and Bronstein's book"
)]

mod parts;
mod rde;
#[cfg(test)]
mod tests;
mod tower;

use std::cell::RefCell;

use cas_num::{Integer, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// Why the algorithm stops: an error of evaluation, or the limits.
pub(super) enum Stop {
    Eval(EvalError),
    Limit,
}

impl From<EvalError> for Stop {
    fn from(e: EvalError) -> Stop {
        Stop::Eval(e)
    }
}

/// A step of the algorithm.
type Step<T> = Result<T, Stop>;

/// Passes over a remainder that still integrates before the call stays.
const ROUNDS: usize = 3;

thread_local! {
    /// The steps that left a remainder or refused an answer, collected
    /// while a [`traced`] call runs on this thread (CAS-585 T9.c6).
    static TRACE: RefCell<Option<Vec<String>>> = const { RefCell::new(None) };
}

/// Notes a step for [`traced`]; nothing is built when no traced call runs.
fn note(message: impl FnOnce() -> String) {
    TRACE.with(|t| {
        if let Some(v) = t.borrow_mut().as_mut() {
            v.push(message());
        }
    });
}

/// `Risch[f, x]`'s answer, as [`risch`] gives it, and the steps noted on
/// the way: which step left the remainder, or that the D check refused
/// the primitive (the tests name the step of a failing row).
#[cfg(test)]
pub(super) fn traced(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> (Option<Expr>, Vec<String>) {
    TRACE.with(|t| *t.borrow_mut() = Some(Vec::new()));
    let r = answer(ev, f, x);
    let mut notes = TRACE.with(|t| t.borrow_mut().take()).unwrap_or_default();
    match r {
        Ok(a) => (a, notes),
        Err(Stop::Eval(e)) => {
            notes.push(format!("error: {e:?}"));
            (None, notes)
        }
        Err(Stop::Limit) => {
            notes.push("the limits".into());
            (None, notes)
        }
    }
}

/// `Risch[f, x]`, `Risch[f]` in `x` (risch.cc:1018).
pub(in crate::native) fn risch(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, x) = match e.args() {
        [f] => (f.clone(), Expr::from(Symbol::new("x"))),
        [f, x] if x.as_symbol().is_some() => (f.clone(), x.clone()),
        _ => return Ok(None),
    };
    match answer(ev, &f, &x) {
        Ok(r) => Ok(r),
        Err(Stop::Eval(err)) => Err(err),
        Err(Stop::Limit) => {
            ev.warn(Warning::Native {
                message: "the tower and its degrees exceed the limits".into(),
                expr: e.clone(),
            });
            Ok(None)
        }
    }
}

/// `prim + Risch[rem, x]`, checked by differentiation; `None` when nothing
/// integrates or the check fails.
fn answer(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Step<Option<Expr>> {
    let (mut prim, mut rem) = risch_top(ev, f, x)?;
    if prim == zero() {
        note(|| format!("risch_top: nothing integrates, remains {rem:?}"));
        return Ok(None);
    }
    let mut stable = false;
    for _ in 0..ROUNDS {
        if rem == zero() {
            stable = true;
            break;
        }
        let (p, r) = risch_top(ev, &rem, x)?;
        if p == zero() {
            stable = true;
            break;
        }
        prim = ev.eval(&(prim + p))?;
        rem = r;
    }
    if !stable {
        note(|| format!("answer: the remainder still integrates after {ROUNDS} rounds: {rem:?}"));
        return Ok(None);
    }
    if !verified(ev, f, &prim, &rem, x)? {
        note(|| format!("guard: D[prim] + rem − f ≠ 0 for prim {prim:?}, rem {rem:?}"));
        return Ok(None);
    }
    // giac's form of `a^u` in the answer: `exp(u ln a)` (pow2expln).
    let prim = ev.eval(&exp_of_powers(&prim, x))?;
    if rem == zero() {
        return Ok(Some(prim));
    }
    Ok(Some(prim + Expr::call(Symbol::RISCH, [rem, x.clone()])))
}

/// Whether `D[prim] + rem − f` is 0.
fn verified(ev: &mut Eval<'_>, f: &Expr, prim: &Expr, rem: &Expr, x: &Expr) -> Step<bool> {
    let d = derive(ev, prim, x)?;
    let difference = ev.eval(&(d + rem.clone() - f.clone()))?;
    vanishes(ev, &difference)
}

/// Whether `d` is 0: by `Together`; after `TrigToExp` and
/// `TranscendentalSimplify`; or, the exponentials written over
/// independent generators ([`exp_generators`]), by `Together` of that
/// rational function. Each test is an identity of values, so a non-zero
/// `d` is never accepted; evaluation alone merges no `Exp[1] Exp[x − 1]`
/// (CAS-585 T9.c6).
pub(super) fn vanishes(ev: &mut Eval<'_>, d: &Expr) -> Step<bool> {
    let d = exps(ev, d)?;
    if together(ev, &d)? == zero() {
        return Ok(true);
    }
    let t0 = call(ev, Symbol::TRIG_TO_EXP, [d])?;
    let t = call(ev, Symbol::TRANSCENDENTAL_SIMPLIFY, [t0.clone()])?;
    if together(ev, &t)? == zero() {
        return Ok(true);
    }
    let g = exp_generators(ev, &t0)?;
    Ok(together(ev, &g)? == zero())
}

/// `e` with `E^u`, `E` and `a^u` (`a` a number other than `E`, `u` not a
/// number) as `Exp`, then each `Exp[u]` as `Π Z_t^(q L_t)`: `u` expanded
/// into terms `q t`, `q` rational (an imaginary `q i` is `q` times the atom
/// `i t`), `L_t` the lcm of the denominators of the `q` of the atom `t`, and
/// `Z_t` a fresh symbol standing for `Exp[t / L_t]`; evaluated.
fn exp_generators(ev: &mut Eval<'_>, e: &Expr) -> Step<Expr> {
    fn as_exp(e: &Expr) -> Expr {
        e.replace(|z| {
            if z.as_symbol() == Some(Symbol::E) {
                return Some(Expr::call(Symbol::EXP, [Expr::from(1)]));
            }
            let View::Apply { head, args: [b, u] } = z.view() else {
                return None;
            };
            if head.as_symbol() != Some(Symbol::POWER) {
                return None;
            }
            if b.as_symbol() == Some(Symbol::E) {
                return Some(Expr::call(Symbol::EXP, [as_exp(u)]));
            }
            if b.is_number() && !u.is_number() {
                let l = Expr::call(Symbol::LOG, [b.clone()]);
                return Some(Expr::call(Symbol::EXP, [as_exp(u) * l]));
            }
            None
        })
    }
    let e = ev.eval(&as_exp(e))?;
    let exps: Vec<Expr> = e
        .distinct()
        .filter(|z| z.head() == Some(&Expr::from(Symbol::EXP)) && z.args().len() == 1)
        .cloned()
        .collect();
    if exps.is_empty() {
        return Ok(e);
    }
    // Each exponential's terms `(q, t)`, and the atoms with their lcm.
    let mut split = Vec::with_capacity(exps.len());
    let mut atoms: Vec<(Expr, Integer)> = Vec::new();
    for k in &exps {
        let u = call(ev, Symbol::EXPAND, [k.args()[0].clone()])?;
        let terms: Vec<Expr> = if u.head() == Some(&Expr::from(Symbol::PLUS)) {
            u.args().to_vec()
        } else {
            vec![u]
        };
        let mut qs = Vec::with_capacity(terms.len());
        for term in &terms {
            let (q, t) = rational_part(ev, term)?;
            let den = q.denom().clone();
            match atoms.iter_mut().find(|(a, _)| *a == t) {
                Some((_, l)) => *l = l.lcm(&den),
                None => atoms.push((t.clone(), den)),
            }
            qs.push((q, t));
        }
        split.push(qs);
    }
    let rules: Vec<(Expr, Expr)> = exps
        .iter()
        .zip(&split)
        .map(|(k, qs)| {
            let factors: Vec<Expr> = qs
                .iter()
                .map(|(q, t)| {
                    let j = atoms.iter().position(|(a, _)| a == t).expect("an atom");
                    let z = Expr::from(Symbol::new(&format!("Risch$z{j}")));
                    let n = q * &Rational::from(&atoms[j].1);
                    z.pow(Expr::from(n))
                })
                .collect();
            (k.clone(), Expr::call(Symbol::TIMES, factors))
        })
        .collect();
    Ok(ev.eval(&e.substitute(&rules))?)
}

/// A term `q t` of an expanded exponent: `q` rational, `t` the rest (`1`
/// for a rational term, `i t′` for an imaginary coefficient).
fn rational_part(ev: &mut Eval<'_>, term: &Expr) -> Step<(Rational, Expr)> {
    let (c, rest) = match term.view() {
        View::Apply { head, args }
            if head.as_symbol() == Some(Symbol::TIMES) && !args.is_empty() =>
        {
            match args[0].to_complex() {
                Some(c) => (c, ev.eval(&Expr::call(Symbol::TIMES, args[1..].to_vec()))?),
                None => return Ok((Rational::from(Integer::from(1)), term.clone())),
            }
        }
        _ => match term.to_complex() {
            Some(c) => (c, Expr::from(1)),
            None => return Ok((Rational::from(Integer::from(1)), term.clone())),
        },
    };
    if c.im().is_zero() {
        return Ok((c.re().clone(), rest));
    }
    if c.re().is_zero() {
        return Ok((c.im().clone(), ev.eval(&(i() * rest))?));
    }
    Ok((Rational::from(Integer::from(1)), term.clone()))
}

// ---------------------------------------------------------------- risch.cc:951-1016

/// giac's `risch` (:951): `(prim, rem)` with `∫ e = prim + ∫ rem`.
fn risch_top(ev: &mut Eval<'_>, e_orig: &Expr, x: &Expr) -> Step<(Expr, Expr)> {
    if *e_orig == zero() {
        return Ok((zero(), zero()));
    }
    let e = call(ev, Symbol::TRIG_TO_EXP, [e_orig.clone()])?;
    let (mut res, try1strem) = risch_lin(ev, &e, x)?;
    let mut remains = zero();
    if try1strem != zero() {
        res = zero();
        let vexp = lin(ev, &e)?;
        let mut remsum = zero();
        for (coeff, expo) in &vexp {
            let term = ev.eval(&(coeff.clone() * Expr::call(Symbol::EXP, [expo.clone()])))?;
            let (r1, rem) = risch_lin(ev, &term, x)?;
            res = ev.eval(&(res + r1))?;
            remsum = ev.eval(&(remsum + rem))?;
        }
        if vexp.is_empty() {
            remains = remsum;
        } else if remsum != zero() {
            let (r2, rem) = risch_lin(ev, &remsum, x)?;
            res = ev.eval(&(res + r2))?;
            remains = rem;
        }
    }
    // invexptoexpneg (intg.cc:2437): exp(u)^-1 is exp(-u).
    let res = ev.eval(&res.replace(|z| match z.view() {
        View::Apply { head, args: [b, k] }
            if head.as_symbol() == Some(Symbol::POWER)
                && *k == Expr::from(-1)
                && b.head() == Some(&Expr::from(Symbol::EXP)) =>
        {
            Some(Expr::call(Symbol::EXP, [-b.args()[0].clone()]))
        }
        _ => None,
    }))?;
    if res == zero() {
        return Ok((zero(), e_orig.clone()));
    }
    let mut res = res;
    if !has_i(e_orig) && has_i(&remains) {
        let tmp = call(ev, Symbol::EXP_TO_TRIG, [remains.clone()])?;
        let tmp = call(ev, Symbol::COMPLEX_EXPAND, [tmp])?;
        let tmp = together(ev, &tmp)?;
        if has_i(&tmp) {
            note(|| format!("risch_top: the remainder {remains:?} is not real"));
            return Ok((zero(), e_orig.clone()));
        }
        remains = tmp;
    }
    if has_any(
        &res,
        &[Symbol::SIN_INTEGRAL, Symbol::COS_INTEGRAL, Symbol::EXP],
    ) {
        res = together(ev, &res)?;
        if !has_i(e_orig) && has_i(&res) {
            res = call(ev, Symbol::EXP_TO_TRIG, [res])?;
            res = together(ev, &res)?;
            if has_i(&res) {
                // re(halftan(res)) in giac: the real part, x real.
                let c = call(ev, Symbol::COMPLEX_EXPAND, [res.clone()])?;
                if !has_i(&c) {
                    res = together(ev, &c)?;
                }
            }
        }
    }
    Ok((res, remains))
}

/// giac's `risch_lin` (:878): `(prim, remains)` of one integrand.
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
fn risch_lin(ev: &mut Eval<'_>, e_orig: &Expr, x: &Expr) -> Step<(Expr, Expr)> {
    let Some(chk) = rlvarx(ev, e_orig, x)? else {
        note(|| format!("risch_lin: no polynomial form of {e_orig:?}"));
        return Ok((zero(), e_orig.clone()));
    };
    let mut chk1 = Vec::new();
    for k in &chk {
        if let View::Apply { head, args: [b, r] } = k.view()
            && head.as_symbol() == Some(Symbol::POWER)
            && r.to_rational().is_some_and(|q| !q.is_integer())
        {
            match b.view() {
                View::Apply {
                    head,
                    args: [bb, rr],
                } if head.as_symbol() == Some(Symbol::POWER) => {
                    chk1.push((
                        k.clone(),
                        Expr::call(Symbol::POWER, [bb.clone(), rr.clone() * r.clone()]),
                    ));
                }
                _ => {
                    note(|| format!("risch_lin: the fractional power {k:?} in {e_orig:?}"));
                    return Ok((zero(), e_orig.clone()));
                }
            }
        }
    }
    let mut e = e_orig.clone();
    if !chk1.is_empty() {
        e = ev.eval(&e.substitute(&chk1))?;
    }
    // atan(g) = i/2 ln((i + g)/(i − g)), written back after (:890-905).
    let vatan: Vec<Expr> = e
        .distinct()
        .filter(|z| z.head() == Some(&Expr::from(Symbol::ARC_TAN)) && z.args().len() == 1)
        .cloned()
        .collect();
    let mut vln = Vec::with_capacity(vatan.len());
    for ga in &vatan {
        let g = ga.args()[0].clone();
        vln.push(together(ev, &((i() + g.clone()) / (i() - g)))?);
    }
    // The constant part of each exponent apart (:906-925).
    let exps_of_e: Vec<Expr> = e
        .distinct()
        .filter(|z| z.head() == Some(&Expr::from(Symbol::EXP)))
        .cloned()
        .collect();
    let mut sub = Vec::new();
    for k in exps_of_e {
        let tmp = k.args()[0].clone();
        let Some(tmpcst) = constant_part(ev, &tmp, x)? else {
            continue;
        };
        if tmpcst == zero() {
            continue;
        }
        let rest = together(ev, &(tmp - tmpcst.clone()))?;
        let new =
            ev.eval(&(Expr::call(Symbol::EXP, [tmpcst]) * Expr::call(Symbol::EXP, [rest])))?;
        if new != k {
            sub.push((k, new));
        }
    }
    if !sub.is_empty() {
        e = ev.eval(&e.substitute(&sub))?;
    }
    let Some((e, v)) = risch_tower(ev, &e, x)? else {
        note(|| format!("risch_tower: no exp/ln tower of {e:?}"));
        return Ok((zero(), e_orig.clone()));
    };
    if v.is_empty() {
        return Ok((ev.eval(&(e * x.clone()))?, zero()));
    }
    let r = parts::in_risch(ev, &e, &v, &zero())?;
    let mut prim = ev.eval(&r.prim)?;
    let mut remains = ev.eval(&r.remains)?;
    if !r.ok || remains != zero() {
        note(|| {
            format!(
                "in_risch: {e:?} in the tower {v:?} leaves {remains:?} (ok {})",
                r.ok
            )
        });
    }
    if has_any(&prim, &[Symbol::SIN_INTEGRAL]) {
        prim = together(ev, &prim)?;
        if has_i(&prim) {
            let t = call(ev, Symbol::EXP_TO_TRIG, [prim])?;
            prim = together(ev, &t)?;
        }
    }
    if !vatan.is_empty() {
        // subst(prim, vln, vatan): ln((i + g)/(i − g)) is −2i atan(g).
        let logs: Vec<Expr> = prim
            .distinct()
            .filter(|z| z.head() == Some(&Expr::from(Symbol::LOG)) && z.args().len() == 1)
            .cloned()
            .collect();
        let mut back = Vec::new();
        for l in logs {
            for (arg, ga) in vln.iter().zip(&vatan) {
                if together(ev, &(l.args()[0].clone() - arg.clone()))? == zero() {
                    back.push((l.clone(), Expr::from(-2) * i() * ga.clone()));
                    break;
                }
            }
        }
        prim = together(ev, &prim.substitute(&back))?;
    }
    if prim == zero() {
        remains = e_orig.clone();
    }
    Ok((prim, remains))
}

/// The constant part of an exponent, its kernels with `x` at 0 (:906-925):
/// `None` when that has no value.
fn constant_part(ev: &mut Eval<'_>, u: &Expr, x: &Expr) -> Step<Option<Expr>> {
    let Some(v) = rlvarx(ev, u, x)? else {
        return Ok(None);
    };
    if v.is_empty() {
        return Ok(Some(u.clone()));
    }
    let mut kernels = v;
    kernels.reverse();
    let Some((c, fr)) = tower::Tower::new(ev, &kernels, std::slice::from_ref(u))? else {
        return Ok(None);
    };
    let r = &c.ring;
    let zero_at = |p: &tower::P| {
        let mut p = p.clone();
        for &var in c.vars() {
            p = r.eval_var(&p, var, &cas_num::Complex::from(0));
        }
        p
    };
    let num = zero_at(&fr[0].num);
    let den = zero_at(&fr[0].den);
    if den.is_zero() {
        return Ok(None);
    }
    Ok(Some(c.r2e(ev, &c.frac(num, den))?))
}

// ---------------------------------------------------------------- the tower, risch.cc:50-97

/// giac's `rlvarx` (series.cc:2175, :2214): the kernels of `e` with `x`,
/// and theirs, recursively, sorted by size (stably; ties by our canonical
/// order). `None` when `e` has no polynomial form.
fn rlvarx(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Step<Option<Vec<Expr>>> {
    let mut res = Vec::new();
    if !collect_kernels(ev, e, x, &mut res)? {
        return Ok(None);
    }
    res.sort_by(|a, b| a.size().cmp(&b.size()).then_with(|| a.cmp(b)));
    Ok(Some(res))
}

fn collect_kernels(ev: &mut Eval<'_>, e: &Expr, x: &Expr, res: &mut Vec<Expr>) -> Step<bool> {
    let Some(gens) = tower::generators(ev, e)? else {
        return Ok(false);
    };
    for k in gens {
        if !k.contains(x) || res.contains(&k) {
            continue;
        }
        res.push(k.clone());
        if let View::Apply { head, args } = k.view() {
            if head.as_symbol() == Some(Symbol::POWER) && args.len() == 2 {
                let (b, r) = (&args[0], &args[1]);
                if !collect_kernels(ev, b, x, res)? || !collect_kernels(ev, r, x, res)? {
                    return Ok(false);
                }
                if r.contains(x) {
                    let l = ev.eval(&Expr::call(Symbol::LOG, [b.clone()]))?;
                    if !collect_kernels(ev, &l, x, res)? {
                        return Ok(false);
                    }
                }
            } else {
                for a in args {
                    if !collect_kernels(ev, a, x, res)? {
                        return Ok(false);
                    }
                }
            }
        }
    }
    Ok(true)
}

/// giac's `in_risch_tower` (:56): `ln(exp(u))` is `u`, normalised; the
/// tower, the top first, when every kernel but `x` is `exp` or `ln`.
fn in_risch_tower(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Step<Option<(Expr, Vec<Expr>)>> {
    let g = ev.eval(&g.replace(|z| match z.view() {
        View::Apply { head, args: [a] }
            if head.as_symbol() == Some(Symbol::LOG)
                && a.head() == Some(&Expr::from(Symbol::EXP)) =>
        {
            Some(a.args()[0].clone())
        }
        _ => None,
    }))?;
    let g = together(ev, &g)?;
    let Some(mut v) = rlvarx(ev, &g, x)? else {
        return Ok(None);
    };
    if v.iter().any(|k| k != x && kernel(k).is_none()) {
        return Ok(None);
    }
    v.reverse();
    Ok(Some((g, v)))
}

/// giac's `risch_tower` (:77): of the integrand with powers rewritten and
/// of its `tsimplify`, the valid tower, the smaller when both are.
fn risch_tower(ev: &mut Eval<'_>, g: &Expr, x: &Expr) -> Step<Option<(Expr, Vec<Expr>)>> {
    let g = ev.eval(&pow2expln_x(g, x))?;
    let g0 = g.clone();
    let g = tsimplify(ev, &g)?;
    let first = in_risch_tower(ev, &g0, x)?;
    let second = in_risch_tower(ev, &g, x)?;
    Ok(match (first, second) {
        (None, s) => s,
        (f, None) => f,
        (Some(f), Some(s)) => {
            if f.1.len() < s.1.len() {
                Some(f)
            } else {
                Some(s)
            }
        }
    })
}

/// giac's `tsimplify` (subst.cc:2098): `atrig2ln`, then
/// `TranscendentalSimplify`.
fn tsimplify(ev: &mut Eval<'_>, e: &Expr) -> Step<Expr> {
    let e = call(ev, Symbol::ARC_TRIG_TO_LOG, [e.clone()])?;
    call(ev, Symbol::TRANSCENDENTAL_SIMPLIFY, [e])
}

/// giac's `pow2expln` in `x` (subst.cc:1787): `a^b` with `b` in `x`, or
/// `a` in `x` and `b` not an integer, is `exp(b ln a)`; built.
fn pow2expln_x(e: &Expr, x: &Expr) -> Expr {
    e.replace(|z| match z.view() {
        View::Apply { head, args: [b, r] }
            if head.as_symbol() == Some(Symbol::POWER)
                && (r.contains(x) || (r.to_integer().is_none() && b.contains(x))) =>
        {
            let g = pow2expln_x(b, x);
            let r = pow2expln_x(r, x);
            if g.head() == Some(&Expr::from(Symbol::EXP)) && g.args().len() == 1 {
                Some(Expr::call(Symbol::EXP, [g.args()[0].clone() * r]))
            } else {
                Some(Expr::call(Symbol::EXP, [r * Expr::call(Symbol::LOG, [g])]))
            }
        }
        _ => None,
    })
}

// ---------------------------------------------------------------- lin.cc

/// giac's `lin` (lin.cc:134): `e` as a sum of `coeff · exp(expo)`, the
/// exponentials of a product merged.
fn lin(ev: &mut Eval<'_>, e: &Expr) -> Step<Vec<(Expr, Expr)>> {
    let View::Apply { head, args } = e.view() else {
        return Ok(vec![(e.clone(), zero())]);
    };
    match head.as_symbol() {
        Some(Symbol::PLUS) => {
            let mut v = Vec::new();
            for a in args {
                v.extend(lin(ev, a)?);
            }
            compress(ev, v)
        }
        Some(Symbol::TIMES) => {
            let mut w = lin(ev, &args[0])?;
            for a in &args[1..] {
                let l = lin(ev, a)?;
                w = convolution(ev, &w, &l)?;
            }
            Ok(w)
        }
        Some(Symbol::POWER) if args.len() == 2 => {
            let (b, n) = (&args[0], &args[1]);
            let k = n.as_i64();
            if let Some(k) = k
                && k < 0
            {
                let pos = ev.eval(&b.pow(Expr::from(-k)))?;
                let w = lin(ev, &pos)?;
                if w.len() == 1 {
                    let (c, ex) = &w[0];
                    return Ok(vec![(ev.eval(&c.pow(-1))?, ev.eval(&-ex.clone())?)]);
                }
                return Ok(vec![(ev.eval(&unlin(&w).pow(-1))?, zero())]);
            }
            let w = lin(ev, b)?;
            if w.len() == 1 {
                let (mut c0, mut e0) = w[0].clone();
                if e0 == zero() && c0.to_integer().is_some() {
                    e0 = Expr::call(Symbol::LOG, [c0]);
                    c0 = Expr::from(1);
                }
                return Ok(vec![(
                    ev.eval(&c0.pow(n.clone()))?,
                    ev.eval(&(e0 * n.clone()))?,
                )]);
            }
            if let Some(k) = k {
                let limit = ev.context().limits().terms;
                let size = u32::try_from(k).ok().map(|k| w.len().saturating_pow(k));
                if size.is_some_and(|s| s <= limit) {
                    let mut r = vec![(Expr::from(1), zero())];
                    for _ in 0..k {
                        r = convolution(ev, &r, &w)?;
                    }
                    return compress(ev, r);
                }
            }
            Ok(vec![(ev.eval(&unlin(&w).pow(n.clone()))?, zero())])
        }
        Some(Symbol::EXP) if args.len() == 1 => Ok(vec![(Expr::from(1), args[0].clone())]),
        _ => Ok(vec![(e.clone(), zero())]),
    }
}

/// lin.cc:52: the terms with one exponent added up.
fn compress(ev: &mut Eval<'_>, v: Vec<(Expr, Expr)>) -> Step<Vec<(Expr, Expr)>> {
    let mut out: Vec<(Expr, Expr)> = Vec::new();
    'terms: for (c, ex) in v {
        for p in &mut out {
            if together(ev, &(p.1.clone() - ex.clone()))? == zero() {
                p.0 = together(ev, &(p.0.clone() + c))?;
                continue 'terms;
            }
        }
        out.push((c, ex));
    }
    Ok(out)
}

/// lin.cc:89-112: the product of two `lin` forms.
fn convolution(
    ev: &mut Eval<'_>,
    v: &[(Expr, Expr)],
    w: &[(Expr, Expr)],
) -> Step<Vec<(Expr, Expr)>> {
    let mut res = Vec::with_capacity(v.len() * w.len());
    for (c1, e1) in v {
        for (c2, e2) in w {
            let c = ev.eval(&(c1.clone() * c2.clone()))?;
            let e = together(ev, &(e1.clone() + e2.clone()))?;
            res.push((c, e));
        }
    }
    compress(ev, res)
}

/// lin.cc:71: the sum of a `lin` form, built.
fn unlin(v: &[(Expr, Expr)]) -> Expr {
    let terms: Vec<Expr> = v
        .iter()
        .filter(|(c, _)| *c != zero())
        .map(|(c, ex)| c.clone() * Expr::call(Symbol::EXP, [ex.clone()]))
        .collect();
    Expr::call(Symbol::PLUS, terms)
}

// ---------------------------------------------------------------- expressions

fn zero() -> Expr {
    Expr::from(0)
}

/// The imaginary unit.
fn i() -> Expr {
    Expr::from(Symbol::I)
}

/// `head[args…]`, evaluated.
fn call<const N: usize>(ev: &mut Eval<'_>, head: Symbol, args: [Expr; N]) -> Step<Expr> {
    Ok(ev.eval(&Expr::call(head, args))?)
}

/// `Together[e]`, giac's `ratnormal` and `normal`.
fn together(ev: &mut Eval<'_>, e: &Expr) -> Step<Expr> {
    call(ev, Symbol::TOGETHER, [e.clone()])
}

/// `D[e, x]`, its `E^u` written `Exp[u]` (C0036).
fn derive(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Step<Expr> {
    let d = call(ev, Symbol::D, [e.clone(), x.clone()])?;
    exps(ev, &d)
}

/// `e` with each `E^u` written `Exp[u]`, evaluated.
fn exps(ev: &mut Eval<'_>, e: &Expr) -> Step<Expr> {
    let e2 = e.replace(|z| match z.view() {
        View::Apply { head, args: [b, u] }
            if head.as_symbol() == Some(Symbol::POWER) && b.as_symbol() == Some(Symbol::E) =>
        {
            Some(Expr::call(Symbol::EXP, [u.clone()]))
        }
        _ => None,
    });
    if e2 == *e {
        return Ok(e2);
    }
    Ok(ev.eval(&e2)?)
}

/// Whether `e` has a complex number or `I`.
fn has_i(e: &Expr) -> bool {
    e.distinct()
        .any(|z| matches!(z.view(), View::Complex(_)) || z.as_symbol() == Some(Symbol::I))
}

/// Whether `e` has an application of one of `heads`.
fn has_any(e: &Expr, heads: &[Symbol]) -> bool {
    e.distinct().any(|z| {
        z.head()
            .and_then(Expr::as_symbol)
            .is_some_and(|h| heads.contains(&h))
    })
}

/// `Some(true)` for `exp(u)` and for `a^u`, `a` a number other than `E` and
/// `u` not a number (the exponential `exp(u ln a)`, which evaluation may
/// write back as `a^u`; CAS-585 T9.c6), `Some(false)` for `ln(u)`, `None`
/// otherwise.
fn kernel(k: &Expr) -> Option<bool> {
    match (k.head().and_then(Expr::as_symbol), k.args()) {
        (Some(Symbol::EXP), [_]) => Some(true),
        (Some(Symbol::LOG), [_]) => Some(false),
        (Some(Symbol::POWER), [a, u]) if power_kernel(a, u) => Some(true),
        _ => None,
    }
}

/// Whether `a^u` is an exponential kernel: `a` a number other than `E`, `u`
/// not a number.
fn power_kernel(a: &Expr, u: &Expr) -> bool {
    a.is_number() && a.as_symbol() != Some(Symbol::E) && !u.is_number()
}

/// The argument of a kernel: `u` of `exp(u)` or `ln(u)`, `u ln a` of `a^u`.
fn kernel_arg(k: &Expr) -> Expr {
    match k.args() {
        [a, u] => u.clone() * Expr::call(Symbol::LOG, [a.clone()]),
        args => args[0].clone(),
    }
}

/// `e` with each `a^u` (`a` a number other than `E`, `u` with `x`) written
/// `exp(u ln a)`, giac's form of an answer; built.
fn exp_of_powers(e: &Expr, x: &Expr) -> Expr {
    e.replace(|z| match z.view() {
        View::Apply { head, args: [a, u] }
            if head.as_symbol() == Some(Symbol::POWER) && power_kernel(a, u) && u.contains(x) =>
        {
            let u = exp_of_powers(u, x);
            Some(Expr::call(
                Symbol::EXP,
                [u * Expr::call(Symbol::LOG, [a.clone()])],
            ))
        }
        _ => None,
    })
}

/// Whether the polynomial-and-kernel sizes stay within `Limits::terms`
/// (C0019): the tower's height times its top degree times the degree bound.
fn within(ev: &Eval<'_>, height: usize, degree: u64, ydeg: i64) -> Step<()> {
    let limit = ev.context().limits().terms as u128;
    let d = u128::from(degree.max(1));
    let y = u128::try_from(ydeg.max(1)).unwrap_or(1);
    if (height as u128).saturating_mul(d).saturating_mul(y) > limit {
        return Err(Stop::Limit);
    }
    Ok(())
}
