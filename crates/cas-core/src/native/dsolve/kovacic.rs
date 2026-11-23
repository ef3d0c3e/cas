//! `KovacicSols` (CAS-584), Xcas's `kovacicsols`: the Liouvillian solutions
//! of `a y″ + b y′ + c y = 0` with `a`, `b`, `c` rational in `x`, by
//! Kovacic's algorithm (J. J. Kovacic, An algorithm for solving second order
//! linear homogeneous differential equations, J. Symbolic Comput. 2 (1986),
//! 3–43), ported from giac 2.0.0.21's desolve.cc: `_kovacicsols`
//! (:2441-2504) reads the arguments, `kovacicsols` (:2122-2430) runs the
//! three cases, with `expansion_coeff` (:1856), `kovacic_iscase1`,
//! `kovacic_iscase2`, `kovacic_iscase3` (:1875-1930) and `fullsimp` (:2111).
//!
//! - Arguments, giac's: `KovacicSols[eq]` in `y[x]`, `KovacicSols[eq, x, y]`,
//!   `KovacicSols[eq, y[x]]`, `KovacicSols[{a, b, c}]` in `x` and
//!   `KovacicSols[{a, b, c}, t]`. The equation is `l = r` or an expression
//!   meaning `= 0`, read as `DSolve` reads it, `y″` written
//!   `Derivative[2][y][x]` (`diff(y(x),x,x)`). Anything but a homogeneous
//!   linear equation of order 2 with coefficients rational in `x` is giac's
//!   Bad Argument Value: the call stays, with that warning.
//! - `y = z e^(−∫p/2)` with `p = b/a`, `q = c/a` turns the equation into
//!   `z″ = r z`, `r = p²/4 + p′/2 − q`. The poles of `r` are the roots of its
//!   denominator factored over Q, linear factors and quadratics (two roots
//!   by a square root, complex ones too), ordered as giac's `gen_map`: the
//!   real ones ascending, then the complex ones by real and imaginary part.
//!   A factor of degree 3 or more leaves the call standing.
//! - Case 1 (:2130-2230) answers one solution `P e^(∫ω)`, case 2
//!   (:2232-2325) the roots of a quadratic in `y′/y`, integrated, case 3
//!   (:2326-2430) the polynomial in `omega_` whose roots are `y′/y`. The
//!   answer is the list of solutions, `{}` for none, and `0` where giac's
//!   code returns `false` (:2157, :2183: the leading coefficient of a square
//!   root series is not a real number).
//! - giac's quirks kept: a coefficient of a partial fraction or series that
//!   is not a real number reads as 0 (:1862); case 1 skips a family of
//!   degree 0 (:2223, its system has no unknowns, which giac's `solve` never
//!   solves); the first family that solves wins; the root with `+√` comes
//!   first in case 2; free unknowns of a system are 0 (:2225).
//! - Departures, fixed rather than copied: case 3 with a `y′` term reads the
//!   shifted polynomial's coefficients lowest power first (giac's `coeff` is
//!   highest first, :2417-2422, which gives the reciprocal polynomial); the
//!   polynomial in `omega_` takes all `n + 1` coefficients for `n = 6, 12`
//!   (:2415 takes five); the vector of polynomials of case 3 starts afresh
//!   for each family (:2388 accumulates over families); a list whose entries
//!   hold `/` is read as written (giac quotes its arguments and answers `0`).
//! - Bound (C0019): the system for a family of degree `d` has `d` unknowns;
//!   its elimination takes about unknowns² × equations steps, which must not
//!   exceed `Limits::terms` (checked as `d³` before the system is built, and
//!   with the equations once it is), nor may the families tried, times the
//!   poles each sums over; past it the call stays, with a warning.
//!   `y″ = (x² − 2001) y` (`d = 1000`) stays.
#![allow(
    clippy::many_single_char_names,
    reason = "the names of Kovacic's paper and giac's code"
)]

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;
use crate::trace::Warning;

use super::{Ode, exponential, without_absolute_values};

type Native = Result<Option<Expr>, EvalError>;

/// What giac's `kovacicsols` (desolve.cc:2122) answers.
enum Answer {
    /// The solutions found, none or more.
    Solutions(Vec<Expr>),
    /// giac's `return false` (:2157, :2183), which prints `0`.
    False,
    /// Case 3: the polynomial in `omega_` whose roots are `y′/y`.
    Omega(Expr),
    /// Past this port or past the limits: the call stays, with the reason.
    Stays(&'static str),
}

/// A pole `c` of `r`, of order `order`.
struct Pole {
    c: Expr,
    order: usize,
    /// `r (x − c)^order`, regular at `c`.
    g: Expr,
    /// The monic quadratic factor `c` is a root of, which a solution's
    /// powers of `x − c` and of its conjugate combine into.
    quadratic: Option<Expr>,
    real: bool,
    re: Expr,
    im: Expr,
}

/// `KovacicSols[eq]`, `[eq, x, y]`, `[eq, y[x]]`, `[{a, b, c}]` or
/// `[{a, b, c}, t]`.
pub(in crate::native) fn kovacicsols(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((a, b, c, x)) = equation(ev, e)? else {
        return refuse(ev, e, "Bad Argument Value");
    };
    if a == Expr::from(0) {
        return refuse(ev, e, "Bad Argument Value");
    }
    let p = together(ev, &(b / a.clone()))?;
    let q = together(ev, &(c / a))?;
    if !rational_in(&p, &x) || !rational_in(&q, &x) {
        return refuse(ev, e, "Bad Argument Value");
    }
    let dp = derivative(ev, &p, &x)?;
    let r = together(
        ev,
        &(p.clone().pow(2) / Expr::from(4) + dp / Expr::from(2) - q),
    )?;
    match solve(ev, &r, &x, &p)? {
        Answer::Solutions(s) => {
            // E^u as giac's exp(u), which Exp[u] reads back as.
            let s: Vec<Expr> = s
                .iter()
                .map(|y| {
                    y.map_bottom_up(|part| {
                        match (part.head().and_then(Expr::as_symbol), part.args()) {
                            (Some(Symbol::POWER), [e, u]) if *e == Expr::from(Symbol::E) => {
                                Expr::call(Symbol::EXP, [u.clone()])
                            }
                            _ => part.clone(),
                        }
                    })
                })
                .collect();
            Ok(Some(ev.eval(&Expr::call(Symbol::LIST, s))?))
        }
        Answer::False => Ok(Some(Expr::from(0))),
        Answer::Omega(w) => Ok(Some(w)),
        Answer::Stays(why) => refuse(ev, e, why),
    }
}

/// The call stays, with `why` as a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

/// `(a, b, c, x)` of `a y″ + b y′ + c y = 0` from the arguments
/// (:2441-2495); `None` for giac's Bad Argument Value.
fn equation(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<(Expr, Expr, Expr, Expr)>, EvalError> {
    let is_list = |l: &Expr| l.head() == Some(&Expr::from(Symbol::LIST));
    let variable = |v: &Expr| v.as_symbol().is_some_and(|s| !s.is_builtin());
    let named = |s: &str| Expr::from(Symbol::new(s));
    let listed = |l: &Expr, x: Expr| match l.args() {
        [a, b, c] => Some((a.clone(), b.clone(), c.clone(), x)),
        _ => None,
    };
    let (eq, y, x) = match e.args() {
        [l] if is_list(l) => return Ok(listed(l, named("x"))),
        [l, t] if is_list(l) => return Ok(listed(l, t.clone()).filter(|_| variable(t))),
        [eq] => (eq, named("y"), named("x")),
        [eq, f] => match (f.head(), f.args()) {
            (Some(y), [x]) if variable(y) && variable(x) => (eq, y.clone(), x.clone()),
            _ => return Ok(None),
        },
        [eq, x, y] if !is_list(eq) && variable(x) && variable(y) => (eq, y.clone(), x.clone()),
        _ => return Ok(None),
    };
    let Some(ode) = Ode::of(ev, eq, &y, &x)? else {
        return Ok(None);
    };
    if ode.order != 2 {
        return Ok(None);
    }
    match ode.linear(ev)? {
        Some((a, f)) if f == Expr::from(0) => {
            Ok(Some((a[2].clone(), a[1].clone(), a[0].clone(), x)))
        }
        _ => Ok(None),
    }
}

/// Whether `e` is a rational function of `x`, its coefficients free of `x`
/// (giac's `rlvarx` test, :2497).
fn rational_in(e: &Expr, x: &Expr) -> bool {
    if e == x || !e.contains(x) {
        return true;
    }
    match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::PLUS | Symbol::TIMES), args) => args.iter().all(|a| rational_in(a, x)),
        (Some(Symbol::POWER), [b, n]) => n.to_integer().is_some() && rational_in(b, x),
        _ => false,
    }
}

/// giac's `kovacicsols(r, x, dy_coeff)` (:2122) for `z″ = r z`, `p` the
/// coefficient of `y′`.
fn solve(ev: &mut Eval<'_>, r: &Expr, x: &Expr, p: &Expr) -> Result<Answer, EvalError> {
    if *r == Expr::from(0) {
        // The square root of r at infinity is 0 (:2183).
        return Ok(Answer::False);
    }
    let s = call(ev, Symbol::NUMERATOR, [r.clone()])?;
    let t = call(ev, Symbol::DENOMINATOR, [r.clone()])?;
    let (Some(ds), Some(dt)) = (degree(ev, &s, x)?, degree(ev, &t, x)?) else {
        return Ok(Answer::Stays("r is not a rational function"));
    };
    let dinf = dt - ds;
    let Some(poles) = poles(ev, &s, &t, x)? else {
        return Ok(Answer::Stays(
            "a pole of r is a root of a factor of degree 3 or more",
        ));
    };
    let lc = |ev: &mut Eval<'_>, f: &Expr| -> Result<Expr, EvalError> {
        Ok(coefficients(ev, f, x)?
            .pop()
            .unwrap_or_else(|| Expr::from(0)))
    };
    let lc_ratio = {
        let (ls, lt) = (lc(ev, &s)?, lc(ev, &t)?);
        together(ev, &(ls / lt))?
    };
    // The coefficient of x^−2 at infinity (:2129).
    let laurent_infinity = if dinf == 2 {
        lc_ratio.clone()
    } else {
        Expr::from(0)
    };
    let infinity = Infinity {
        s: &s,
        t: &t,
        dinf,
        lc_ratio: &lc_ratio,
    };

    let case1 = !(dinf.rem_euclid(2) == 1 && dinf <= 2)
        && poles.iter().all(|p| p.order == 1 || p.order % 2 == 0);
    if case1 {
        match case_one(ev, r, x, p, &poles, &infinity)? {
            Answer::Solutions(s) if s.is_empty() => {}
            answer => return Ok(answer),
        }
    }
    if poles
        .iter()
        .any(|p| p.order == 2 || (p.order > 2 && p.order % 2 == 1))
    {
        match case_two(ev, r, x, p, &poles, dinf, &laurent_infinity)? {
            Answer::Solutions(s) if s.is_empty() => {}
            answer => return Ok(answer),
        }
    }
    if case_three_applies(ev, x, &poles, dinf)? {
        return case_three(ev, r, x, p, &poles, dinf, &laurent_infinity);
    }
    Ok(Answer::Solutions(Vec::new()))
}

/// What case 1 reads at infinity.
struct Infinity<'a> {
    s: &'a Expr,
    t: &'a Expr,
    dinf: i64,
    lc_ratio: &'a Expr,
}

/// One pole's or infinity's data in case 1: the square root part and
/// `α⁺`, `α⁻`.
struct Side {
    sq: Expr,
    plus: Expr,
    minus: Expr,
}

/// Case 1 (:2130-2230).
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
fn case_one(
    ev: &mut Eval<'_>,
    r: &Expr,
    x: &Expr,
    p: &Expr,
    poles: &[Pole],
    infinity: &Infinity<'_>,
) -> Result<Answer, EvalError> {
    let half = |e: Expr| e / Expr::from(2);
    let mut sides = Vec::with_capacity(poles.len());
    for pole in poles {
        let side = match pole.order {
            1 => Side {
                sq: Expr::from(0),
                plus: Expr::from(1),
                minus: Expr::from(1),
            },
            2 => {
                let b = laurent(ev, pole, x, 2)?;
                let b = real(ev, b)?;
                let root = sqrt(ev, &(Expr::from(1) + Expr::from(4) * b))?;
                Side {
                    sq: Expr::from(0),
                    plus: eval(ev, &half(Expr::from(1) + root.clone()))?,
                    minus: eval(ev, &half(Expr::from(1) - root))?,
                }
            }
            order => {
                // Even order 2ν ≥ 4 (:2147-2167).
                let nu = order / 2;
                let g = taylor(ev, &pole.g, x, &pole.c, nu - 1)?;
                let roots = sqrt_series(ev, &g)?;
                let mut cf = vec![Expr::from(0); nu + 1];
                let mut sq = Vec::new();
                for i in 2..=nu {
                    cf[i] = real(ev, roots[nu - i].clone())?;
                    sq.push(cf[i].clone() / (x.clone() - pole.c.clone()).pow(int(i)));
                }
                let a = cf[nu].clone();
                if a == Expr::from(0) {
                    return Ok(Answer::False);
                }
                let square: Vec<Expr> = (2..=nu)
                    .filter(|&i| nu + 1 >= i + 2 && nu + 1 - i <= nu)
                    .map(|i| cf[i].clone() * cf[nu + 1 - i].clone())
                    .collect();
                let square = eval(ev, &Expr::call(Symbol::PLUS, square))?;
                let b = laurent(ev, pole, x, nu + 1)?;
                let b = real(ev, b)? - real(ev, square)?;
                let ratio = b / a;
                let nu = int(nu);
                Side {
                    sq: eval(ev, &Expr::call(Symbol::PLUS, sq))?,
                    plus: eval(ev, &half(nu.clone() + ratio.clone()))?,
                    minus: eval(ev, &half(nu - ratio))?,
                }
            }
        };
        sides.push(side);
    }
    let at_infinity = match infinity.dinf {
        d if d > 2 => Side {
            sq: Expr::from(0),
            plus: Expr::from(0),
            minus: Expr::from(1),
        },
        2 => {
            let root = sqrt(
                ev,
                &(Expr::from(1) + Expr::from(4) * infinity.lc_ratio.clone()),
            )?;
            Side {
                sq: Expr::from(0),
                plus: eval(ev, &half(Expr::from(1) + root.clone()))?,
                minus: eval(ev, &half(Expr::from(1) - root))?,
            }
        }
        d => {
            // dinf even ≤ 0 (:2168-2190): the series of r at infinity in
            // u = 1/x is the quotient of the reversed coefficients.
            let nu = usize::try_from(-d / 2).unwrap_or(0);
            let sc = coefficients(ev, infinity.s, x)?;
            let tc = coefficients(ev, infinity.t, x)?;
            let (sc, tc): (Vec<Expr>, Vec<Expr>) = (
                sc.into_iter().rev().collect(),
                tc.into_iter().rev().collect(),
            );
            let series = divide_series(ev, &sc, &tc, nu + 1)?;
            let roots = sqrt_series(ev, &series)?;
            let mut cf = vec![Expr::from(0); nu + 1];
            for (i, c) in cf.iter_mut().enumerate() {
                *c = real(ev, roots[nu - i].clone())?;
            }
            let a = cf[nu].clone();
            if a == Expr::from(0) {
                return Ok(Answer::False);
            }
            // b: the coefficient of x^(ν−1) in propfrac(r) (:2184), read
            // term by term: for ν = 0 only a remainder c/x counts.
            let qr = call(
                ev,
                Symbol::POLYNOMIAL_QUOTIENT_REMAINDER,
                [infinity.s.clone(), infinity.t.clone(), x.clone()],
            )?;
            let [quotient, remainder] = qr.args() else {
                return Ok(Answer::Stays("r is not a rational function"));
            };
            let mut b = if nu >= 1 {
                let qc = coefficients(ev, quotient, x)?;
                qc.get(nu - 1).cloned().unwrap_or_else(|| Expr::from(0))
            } else {
                let term = together(ev, &(remainder.clone() * x.clone() / infinity.t.clone()))?;
                if term.contains(x) {
                    Expr::from(0)
                } else {
                    term
                }
            };
            b = real(ev, b)?;
            if nu >= 1 {
                let square: Vec<Expr> = (0..nu)
                    .map(|i| cf[i].clone() * cf[nu - 1 - i].clone())
                    .collect();
                b = eval(ev, &(b - Expr::call(Symbol::PLUS, square)))?;
            }
            let sq: Vec<Expr> = cf
                .iter()
                .enumerate()
                .map(|(i, c)| c.clone() * x.clone().pow(int(i)))
                .collect();
            let ratio = b / a;
            let nu = int(nu);
            Side {
                sq: eval(ev, &Expr::call(Symbol::PLUS, sq))?,
                plus: eval(ev, &half(ratio.clone() - nu.clone()))?,
                minus: eval(ev, &half(-ratio - nu))?,
            }
        }
    };

    // Step 2 (:2191-2212): a sign for infinity (bit 0) and each pole.
    let np = poles.len() + 1;
    let limit = ev.context().limits().terms;
    if np >= 40 || (1usize << np).saturating_mul(np) > limit {
        return Ok(Answer::Stays("too many families for the limits"));
    }
    let mut families: Vec<(usize, usize)> = Vec::new();
    for i in 0..1usize << np {
        let pick = |side: &Side, plus: bool| {
            if plus {
                side.plus.clone()
            } else {
                side.minus.clone()
            }
        };
        let mut terms = vec![pick(&at_infinity, i & 1 != 0)];
        for (j, side) in sides.iter().enumerate() {
            terms.push(-pick(side, i & (1 << (j + 1)) != 0));
        }
        let d = integer(ev, &Expr::call(Symbol::PLUS, terms))?;
        if let Some(d) = d.and_then(|d| usize::try_from(d).ok()) {
            families.push((d, i));
        }
    }

    // Step 3 (:2213-2229): P with P″ + 2ω P′ + (ω′ + ω² − r) P = 0.
    for (d, i) in families {
        if d == 0 {
            // giac's solve of no unknowns never succeeds (:2223).
            continue;
        }
        let sign = |plus: bool| if plus { Expr::from(1) } else { Expr::from(-1) };
        let mut fw = vec![sign(i & 1 != 0) * at_infinity.sq.clone()];
        for (j, (side, pole)) in sides.iter().zip(poles).enumerate() {
            let plus = i & (1 << (j + 1)) != 0;
            let alpha = if plus { &side.plus } else { &side.minus };
            fw.push(sign(plus) * side.sq.clone());
            fw.push(alpha.clone() / (x.clone() - pole.c.clone()));
        }
        let fw = together(ev, &Expr::call(Symbol::PLUS, fw))?;
        let dfw = derivative(ev, &fw, x)?;
        let c1 = together(ev, &(Expr::from(2) * fw.clone()))?;
        let c0 = together(ev, &(dfw + fw.clone().pow(2) - r.clone()))?;
        let coefficients = [c0, c1, Expr::from(1)];
        let Some(columns) = columns(ev, x, d, &coefficients)? else {
            return Ok(Answer::Stays("the degree bound exceeds the limits"));
        };
        let big_p = match kernel(ev, x, d, &columns)? {
            Kernel::Found(big_p) => big_p,
            Kernel::None => continue,
            Kernel::Past => return Ok(Answer::Stays("the degree bound exceeds the limits")),
        };
        // y = P e^(∫ω − p/2) (:2226), each pole's α as a power.
        let mut exponent = vec![sign(i & 1 != 0) * integral(ev, &at_infinity.sq, x)?];
        let mut powers: Vec<(Expr, Option<Expr>, Expr)> = Vec::new();
        for (j, (side, pole)) in sides.iter().zip(poles).enumerate() {
            let plus = i & (1 << (j + 1)) != 0;
            exponent.push(sign(plus) * integral(ev, &side.sq, x)?);
            let alpha = if plus { &side.plus } else { &side.minus };
            powers.push((pole.c.clone(), pole.quadratic.clone(), alpha.clone()));
        }
        exponent.push(-integral(ev, &(p.clone() / Expr::from(2)), x)?);
        let factors = combined_powers(ev, x, &powers)?;
        let e = exponential_of(ev, &Expr::call(Symbol::PLUS, exponent))?;
        let mut product = vec![big_p];
        product.extend(factors);
        product.push(e);
        let y = eval(ev, &Expr::call(Symbol::TIMES, product))?;
        return Ok(Answer::Solutions(vec![y]));
    }
    Ok(Answer::Solutions(Vec::new()))
}

/// `Π (x − c)^α` over the poles, a conjugate pair with the same `α` as its
/// monic quadratic to the `α`.
fn combined_powers(
    ev: &mut Eval<'_>,
    x: &Expr,
    powers: &[(Expr, Option<Expr>, Expr)],
) -> Result<Vec<Expr>, EvalError> {
    let mut done = vec![false; powers.len()];
    let mut out = Vec::new();
    for k in 0..powers.len() {
        if done[k] {
            continue;
        }
        done[k] = true;
        let (c, quadratic, alpha) = &powers[k];
        if let Some(f) = quadratic {
            let mate = (k + 1..powers.len())
                .find(|&m| !done[m] && powers[m].1.as_ref() == Some(f) && powers[m].2 == *alpha);
            if let Some(m) = mate {
                done[m] = true;
                out.push(f.clone().pow(alpha.clone()));
                continue;
            }
        }
        out.push(eval(ev, &(x.clone() - c.clone()))?.pow(alpha.clone()));
    }
    Ok(out)
}

/// Case 2 (:2232-2325).
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
fn case_two(
    ev: &mut Eval<'_>,
    r: &Expr,
    x: &Expr,
    p: &Expr,
    poles: &[Pole],
    dinf: i64,
    laurent_infinity: &Expr,
) -> Result<Answer, EvalError> {
    let mut sets = Vec::with_capacity(poles.len());
    for pole in poles {
        let set = match pole.order {
            1 => vec![4],
            2 => {
                let b = laurent(ev, pole, x, 2)?;
                let b = real(ev, b)?;
                exponents(ev, &b, 2, 2, &[-1, 0, 1])?
            }
            order => vec![i64::try_from(order).unwrap_or(i64::MAX)],
        };
        sets.push(set);
    }
    let set_infinity = match dinf {
        d if d > 2 => vec![0, 2, 4],
        2 => exponents(ev, laurent_infinity, 2, 2, &[-1, 0, 1])?,
        d => vec![d],
    };
    let Some(families) = families(ev, &sets, set_infinity.len())? else {
        return Ok(Answer::Stays("too many families for the limits"));
    };
    let mut fam: Vec<(usize, Expr)> = Vec::new();
    for &ei in &set_infinity {
        if families.is_empty() && ei % 2 == 0 && ei >= 0 {
            fam.push((usize::try_from(ei / 2).unwrap_or(0), Expr::from(0)));
        }
        for fm in &families {
            let sum: i64 = fm.iter().sum();
            if (ei - sum) % 2 == 0 && ei >= sum {
                let th = theta(ev, x, poles, fm)?;
                let th = together(ev, &(th / Expr::from(2)))?;
                fam.push((usize::try_from((ei - sum) / 2).unwrap_or(0), th));
            }
        }
    }
    let dr = derivative(ev, r, x)?;
    for (deg, th) in fam {
        let dth = derivative(ev, &th, x)?;
        let second = derivative(ev, &dth, x)?;
        let three = Expr::from(3);
        let four = Expr::from(4);
        let c2 = together(ev, &(three.clone() * th.clone()))?;
        let c1 = together(
            ev,
            &(three.clone() * th.clone().pow(2) + three.clone() * dth.clone()
                - four.clone() * r.clone()),
        )?;
        let c0 = together(
            ev,
            &(second + three * th.clone() * dth + th.clone().pow(3)
                - four * r.clone() * th.clone()
                - Expr::from(2) * dr.clone()),
        )?;
        let big_p = if deg == 0 {
            if c0 != Expr::from(0) {
                continue;
            }
            Expr::from(1)
        } else {
            let Some(columns) = columns(ev, x, deg, &[c0, c1, c2, Expr::from(1)])? else {
                return Ok(Answer::Stays("the degree bound exceeds the limits"));
            };
            match kernel(ev, x, deg, &columns)? {
                Kernel::Found(big_p) => big_p,
                Kernel::None => continue,
                Kernel::Past => {
                    return Ok(Answer::Stays("the degree bound exceeds the limits"));
                }
            }
        };
        // ω² − φ ω + φ′/2 + φ²/2 − r = 0, φ = θ + P′/P (:2316-2321).
        let dp = derivative(ev, &big_p, x)?;
        let ph = together(ev, &(th + dp / big_p))?;
        let dph = derivative(ev, &ph, x)?;
        let delta = together(
            ev,
            &(Expr::from(4) * r.clone() - ph.clone().pow(2) - Expr::from(2) * dph),
        )?;
        let rational = integral(ev, &((ph - p.clone()) / Expr::from(2)), x)?;
        let mut solutions = Vec::new();
        if delta == Expr::from(0) {
            solutions.push(exponential_of(ev, &rational)?);
        } else {
            let root = square_root(ev, &delta)?;
            let radical = integral(ev, &(root / Expr::from(2)), x)?;
            for sign in [1, -1] {
                let e = rational.clone() + Expr::from(sign) * radical.clone();
                solutions.push(exponential_of(ev, &e)?);
            }
        }
        return Ok(Answer::Solutions(solutions));
    }
    Ok(Answer::Solutions(Vec::new()))
}

/// The integers `base + k · step · √(1 + 4b)` for `k` in `ks`, in order,
/// repeats kept (:2245, :2341).
fn exponents(
    ev: &mut Eval<'_>,
    b: &Expr,
    base: i64,
    step: i64,
    ks: &[i64],
) -> Result<Vec<i64>, EvalError> {
    let root = sqrt(ev, &(Expr::from(1) + Expr::from(4) * b.clone()))?;
    let mut out = Vec::new();
    for &k in ks {
        let e = Expr::from(base) + Expr::from(k * step) * root.clone();
        if let Some(v) = integer(ev, &e)? {
            out.push(v);
        }
    }
    Ok(out)
}

/// The families of exponents, one from each pole's set, the first pole's
/// slowest (`build_E_families`); `None` past the limits.
fn families(
    ev: &Eval<'_>,
    sets: &[Vec<i64>],
    at_infinity: usize,
) -> Result<Option<Vec<Vec<i64>>>, EvalError> {
    if sets.is_empty() {
        return Ok(Some(Vec::new()));
    }
    let limit = ev.context().limits().terms;
    let mut count = at_infinity.max(1).saturating_mul(sets.len());
    for set in sets {
        count = count.saturating_mul(set.len());
        if count > limit {
            return Ok(None);
        }
    }
    let mut out: Vec<Vec<i64>> = vec![Vec::new()];
    for set in sets {
        let mut next = Vec::with_capacity(out.len() * set.len());
        for f in &out {
            for &e in set {
                let mut g = f.clone();
                g.push(e);
                next.push(g);
            }
        }
        out = next;
    }
    Ok(Some(out))
}

/// `Σ eₖ/(x − cₖ)` over the poles.
fn theta(ev: &mut Eval<'_>, x: &Expr, poles: &[Pole], fm: &[i64]) -> Result<Expr, EvalError> {
    let terms: Vec<Expr> = fm
        .iter()
        .zip(poles)
        .map(|(&e, pole)| Expr::from(e) / (x.clone() - pole.c.clone()))
        .collect();
    together(ev, &Expr::call(Symbol::PLUS, terms))
}

/// giac's `kovacic_iscase3` (:1897-1930).
fn case_three_applies(
    ev: &mut Eval<'_>,
    x: &Expr,
    poles: &[Pole],
    dinf: i64,
) -> Result<bool, EvalError> {
    if dinf < 2 || poles.iter().any(|p| p.order > 2) {
        return Ok(false);
    }
    let mut beta = Vec::new();
    let mut g = Vec::new();
    for pole in poles {
        if pole.order == 2 {
            let a = laurent(ev, pole, x, 2)?;
            let a = real(ev, a)?;
            g.push(a.clone());
            let root = sqrt(ev, &(Expr::from(1) + Expr::from(4) * a))?;
            if root.to_rational().is_none() {
                return Ok(false);
            }
        }
        let b = laurent(ev, pole, x, 1)?;
        let b = real(ev, b)?;
        beta.push(b.clone());
        g.push(b * pole.c.clone());
    }
    let sum = together(ev, &Expr::call(Symbol::PLUS, beta))?;
    let g = eval(ev, &Expr::call(Symbol::PLUS, g))?;
    let root = sqrt(ev, &(Expr::from(1) + Expr::from(4) * g))?;
    Ok(sum == Expr::from(0) && root.to_rational().is_some())
}

/// Case 3 (:2326-2430), for n = 4, 6, 12.
#[allow(
    clippy::too_many_lines,
    reason = "one step of giac's code after another"
)]
fn case_three(
    ev: &mut Eval<'_>,
    r: &Expr,
    x: &Expr,
    p: &Expr,
    poles: &[Pole],
    dinf: i64,
    laurent_infinity: &Expr,
) -> Result<Answer, EvalError> {
    let s: Vec<Expr> = poles
        .iter()
        .map(|pole| x.clone() - pole.c.clone())
        .collect();
    let s = call(ev, Symbol::EXPAND, [Expr::call(Symbol::TIMES, s)])?;
    let ds = derivative(ev, &s, x)?;
    let s2r = together(ev, &(s.clone().pow(2) * r.clone()))?;
    for n in [4i64, 6, 12] {
        let half = n / 2;
        let ks: Vec<i64> = (-half..=half).collect();
        let mut sets = Vec::with_capacity(poles.len());
        for pole in poles {
            let set = if pole.order == 1 {
                vec![12]
            } else {
                let a = laurent(ev, pole, x, 2)?;
                let a = real(ev, a)?;
                exponents(ev, &a, 6, 12 / n, &ks)?
            };
            sets.push(set);
        }
        let b = if dinf > 2 {
            Expr::from(0)
        } else {
            laurent_infinity.clone()
        };
        let set_infinity = exponents(ev, &b, 6, 12 / n, &ks)?;
        let Some(families) = families(ev, &sets, set_infinity.len())? else {
            return Ok(Answer::Stays("too many families for the limits"));
        };
        let mut fam: Vec<(usize, Expr)> = Vec::new();
        for &ei in &set_infinity {
            if families.is_empty() && (ei * n) % 12 == 0 && ei >= 0 {
                fam.push((usize::try_from(ei * n / 12).unwrap_or(0), Expr::from(0)));
            }
            for fm in &families {
                let sum: i64 = fm.iter().sum();
                let e = n * (ei - sum);
                if e % 12 == 0 && e >= 0 {
                    let th = theta(ev, x, poles, fm)?;
                    let th = together(ev, &(Expr::from(n) * th / Expr::from(12)))?;
                    fam.push((usize::try_from(e / 12).unwrap_or(0), th));
                }
            }
        }
        let n = usize::try_from(n).unwrap_or(4);
        for (deg, th) in fam {
            let limit = ev.context().limits().terms;
            if deg.saturating_pow(3).saturating_mul(n) > limit {
                return Ok(Answer::Stays("the degree bound exceeds the limits"));
            }
            // P[n+1] = −x^i for each power, then the recurrence (:2389-2399);
            // P starts afresh for each family (departure from :2388).
            let mut chains = Vec::with_capacity(deg + 1);
            for i in 0..=deg {
                let chain = recurrence(ev, x, n, i, &s, &ds, &s2r, &th)?;
                chains.push(chain);
            }
            let big_p = if deg == 0 {
                if chains[0][0] != Expr::from(0) {
                    continue;
                }
                chains.pop().unwrap_or_default()
            } else {
                let mut columns = Vec::with_capacity(deg + 1);
                for chain in &chains {
                    columns.push(coefficients(ev, &chain[0], x)?);
                }
                let v = match solve_monic(ev, deg, &columns)? {
                    System::Solved(v) => v,
                    System::Inconsistent => continue,
                    System::Past => {
                        return Ok(Answer::Stays("the degree bound exceeds the limits"));
                    }
                };
                let mut combined = Vec::with_capacity(n + 2);
                for k in 0..n + 2 {
                    let mut terms: Vec<Expr> = v
                        .iter()
                        .zip(&chains)
                        .map(|(vi, chain)| vi.clone() * chain[k].clone())
                        .collect();
                    terms.push(chains[deg][k].clone());
                    combined.push(together(ev, &Expr::call(Symbol::PLUS, terms))?);
                }
                combined
            };
            // ω's polynomial: Σ ωⁱ Sⁱ P[i+1]/(n−i)! (:2408-2415), all n+1
            // coefficients.
            let mut ac = Vec::with_capacity(n + 1);
            for i in 0..=n {
                let term = s.clone().pow(int(i)) * big_p[i + 1].clone() / int(factorial(n - i));
                ac.push(together(ev, &term)?);
            }
            let mut ac = strip_gcd(ev, ac)?;
            if *p != Expr::from(0) {
                // ω → ω + p/2, read lowest power first (departure from
                // :2417-2422).
                let shift = p.clone() / Expr::from(2);
                let mut shifted = Vec::with_capacity(n + 1);
                for k in 0..=n {
                    let terms: Vec<Expr> = (k..=n)
                        .map(|i| {
                            int(binomial(i, k)) * ac[i].clone() * shift.clone().pow(int(i - k))
                        })
                        .collect();
                    shifted.push(together(ev, &Expr::call(Symbol::PLUS, terms))?);
                }
                ac = strip_gcd(ev, shifted)?;
            }
            let w = Expr::from(Symbol::new("omega_"));
            let terms: Vec<Expr> = ac
                .into_iter()
                .enumerate()
                .map(|(i, c)| w.clone().pow(int(i)) * c)
                .collect();
            return Ok(Answer::Omega(eval(ev, &Expr::call(Symbol::PLUS, terms))?));
        }
    }
    Ok(Answer::Solutions(Vec::new()))
}

/// The vector `P[0..=n+1]` of case 3 for `P[n+1] = −x^i` (:2389-2399).
#[allow(
    clippy::too_many_arguments,
    reason = "the quantities of giac's recurrence"
)]
fn recurrence(
    ev: &mut Eval<'_>,
    x: &Expr,
    n: usize,
    i: usize,
    s: &Expr,
    ds: &Expr,
    s2r: &Expr,
    th: &Expr,
) -> Result<Vec<Expr>, EvalError> {
    let mut chain = vec![Expr::from(0); n + 2];
    chain[n + 1] = -x.clone().pow(int(i));
    let d = derivative(ev, &chain[n + 1], x)?;
    chain[n] = together(ev, &(-s.clone() * (d + th.clone() * chain[n + 1].clone())))?;
    for k in (0..n).rev() {
        let d = derivative(ev, &chain[k + 1], x)?;
        let m = int(n - k);
        let w = int((n - k) * (k + 1));
        let e = -s.clone() * d + (m * ds.clone() - s.clone() * th.clone()) * chain[k + 1].clone()
            - w * s2r.clone() * chain[k + 2].clone();
        let e = together(ev, &e)?;
        chain[k] = call(ev, Symbol::EXPAND, [e])?;
    }
    Ok(chain)
}

/// giac's `strip_gcd`: the entries times the gcd of the denominators over
/// the gcd of the nonzero numerators, factored.
fn strip_gcd(ev: &mut Eval<'_>, v: Vec<Expr>) -> Result<Vec<Expr>, EvalError> {
    let mut g1: Option<Expr> = None;
    let mut g2: Option<Expr> = None;
    for e in &v {
        let num = call(ev, Symbol::NUMERATOR, [e.clone()])?;
        let den = call(ev, Symbol::DENOMINATOR, [e.clone()])?;
        if num != Expr::from(0) {
            g1 = Some(match g1 {
                Some(g) => call(ev, Symbol::GCD, [g, num])?,
                None => num,
            });
        }
        g2 = Some(match g2 {
            Some(g) => call(ev, Symbol::GCD, [g, den])?,
            None => den,
        });
    }
    let (g1, g2) = (
        g1.unwrap_or_else(|| Expr::from(1)),
        g2.unwrap_or_else(|| Expr::from(1)),
    );
    let mut out = Vec::with_capacity(v.len());
    for e in v {
        let e = together(ev, &(g2.clone() / g1.clone() * e))?;
        out.push(call(ev, Symbol::FACTOR, [e])?);
    }
    Ok(out)
}

/// The columns `L(xⁱ)`, `i ≤ d`, of `L = Σ cₖ Dᵏ` over a common
/// denominator, as coefficient lists (lowest power first); `None` past the
/// limits.
fn columns(
    ev: &mut Eval<'_>,
    x: &Expr,
    d: usize,
    c: &[Expr],
) -> Result<Option<Vec<Vec<Expr>>>, EvalError> {
    if d.saturating_pow(3) > ev.context().limits().terms {
        return Ok(None);
    }
    let mut den = Expr::from(1);
    for ck in c {
        let dk = call(ev, Symbol::DENOMINATOR, [ck.clone()])?;
        den = call(ev, Symbol::LCM, [den, dk])?;
    }
    let mut out = Vec::with_capacity(d + 1);
    for i in 0..=d {
        let mut terms = Vec::with_capacity(c.len());
        let mut falling = 1usize;
        for (k, ck) in c.iter().enumerate() {
            if k > i {
                break;
            }
            terms.push(int(falling) * ck.clone() * x.clone().pow(int(i - k)));
            falling = falling.saturating_mul(i - k);
        }
        let l = together(ev, &(den.clone() * Expr::call(Symbol::PLUS, terms)))?;
        let l = call(ev, Symbol::EXPAND, [l])?;
        out.push(coefficients(ev, &l, x)?);
    }
    Ok(Some(out))
}

/// What [`kernel`] finds.
enum Kernel {
    Found(Expr),
    None,
    Past,
}

/// The monic `P` of degree `d` the columns send to 0 (:2214-2226).
fn kernel(
    ev: &mut Eval<'_>,
    x: &Expr,
    d: usize,
    columns: &[Vec<Expr>],
) -> Result<Kernel, EvalError> {
    let v = match solve_monic(ev, d, columns)? {
        System::Solved(v) => v,
        System::Inconsistent => return Ok(Kernel::None),
        System::Past => return Ok(Kernel::Past),
    };
    let mut terms: Vec<Expr> = v
        .into_iter()
        .enumerate()
        .map(|(i, vi)| vi * x.clone().pow(int(i)))
        .collect();
    terms.push(x.clone().pow(int(d)));
    Ok(Kernel::Found(eval(ev, &Expr::call(Symbol::PLUS, terms))?))
}

/// What [`solve_monic`] finds.
enum System {
    Solved(Vec<Expr>),
    Inconsistent,
    Past,
}

/// `v` with `Σ_{i<d} vᵢ colᵢ + col_d = 0`, free unknowns 0 (:2225).
fn solve_monic(ev: &mut Eval<'_>, d: usize, columns: &[Vec<Expr>]) -> Result<System, EvalError> {
    let m = columns.iter().map(Vec::len).max().unwrap_or(0);
    if d.saturating_mul(d).saturating_mul(m) > ev.context().limits().terms {
        return Ok(System::Past);
    }
    let at = |col: &Vec<Expr>, k: usize| col.get(k).cloned().unwrap_or_else(|| Expr::from(0));
    let mut a: Vec<Vec<Expr>> = (0..m)
        .map(|k| (0..d).map(|i| at(&columns[i], k)).collect())
        .collect();
    let mut b: Vec<Expr> = (0..m).map(|k| -at(&columns[d], k)).collect();
    let mut pivots = Vec::new();
    let mut row = 0;
    for col in 0..d {
        let mut found = None;
        for (k, a_k) in a.iter().enumerate().skip(row) {
            if !is_zero(ev, &a_k[col])? {
                found = Some(k);
                break;
            }
        }
        let Some(k) = found else {
            continue;
        };
        a.swap(row, k);
        b.swap(row, k);
        for k in 0..m {
            if k == row || is_zero(ev, &a[k][col])? {
                continue;
            }
            let f = eval(ev, &(a[k][col].clone() / a[row][col].clone()))?;
            let pivot = a[row].clone();
            for (entry, p) in a[k].iter_mut().zip(&pivot).skip(col) {
                if *p == Expr::from(0) {
                    continue;
                }
                *entry = eval(ev, &(entry.clone() - f.clone() * p.clone()))?;
            }
            b[k] = eval(ev, &(b[k].clone() - f * b[row].clone()))?;
        }
        pivots.push(col);
        row += 1;
    }
    for bk in &b[row..] {
        if !is_zero(ev, bk)? {
            return Ok(System::Inconsistent);
        }
    }
    let mut v = vec![Expr::from(0); d];
    for (k, &col) in pivots.iter().enumerate() {
        v[col] = eval(ev, &(b[k].clone() / a[k][col].clone()))?;
    }
    Ok(System::Solved(v))
}

/// The poles of `s/t` in giac's order (:2124-2128); `None` when a factor
/// of `t` has degree 3 or more, or coefficients that are not rational.
fn poles(ev: &mut Eval<'_>, s: &Expr, t: &Expr, x: &Expr) -> Result<Option<Vec<Pole>>, EvalError> {
    let list = call(ev, Symbol::FACTOR_LIST, [t.clone()])?;
    let pairs: Vec<(Expr, usize)> = list
        .args()
        .chunks(2)
        .filter_map(|c| match c {
            [f, m] => Some((f.clone(), usize::try_from(m.as_i64()?).ok()?)),
            _ => None,
        })
        .collect();
    let mut found = Vec::new();
    for (k, (f, m)) in pairs.iter().enumerate() {
        let cl = coefficients(ev, f, x)?;
        if cl.len() <= 1 {
            continue;
        }
        if cl.iter().any(|c| c.to_rational().is_none()) || cl.len() > 3 {
            return Ok(None);
        }
        let rest: Vec<Expr> = pairs
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != k)
            .map(|(_, (g, mg))| g.clone().pow(int(*mg)))
            .collect();
        let rest = Expr::call(Symbol::TIMES, rest);
        let mi = int(*m);
        if let [c0, c1] = cl.as_slice() {
            let c = eval(ev, &(-c0.clone() / c1.clone()))?;
            let g = together(ev, &(s.clone() / (rest * c1.clone().pow(mi))))?;
            found.push(Pole {
                re: c.clone(),
                im: Expr::from(0),
                c,
                order: *m,
                g,
                quadratic: None,
                real: true,
            });
            continue;
        }
        let [a0, a1, a2] = cl.as_slice() else {
            return Ok(None);
        };
        let disc = eval(
            ev,
            &(a1.clone().pow(2) - Expr::from(4) * a2.clone() * a0.clone()),
        )?;
        let real = disc.to_rational().is_some_and(|d| d.is_positive());
        let root = sqrt(ev, &disc)?;
        let two_a2 = Expr::from(2) * a2.clone();
        let roots = [
            eval(ev, &((-a1.clone() - root.clone()) / two_a2.clone()))?,
            eval(ev, &((-a1.clone() + root) / two_a2.clone()))?,
        ];
        let monic = call(ev, Symbol::EXPAND, [f.clone() / a2.clone()])?;
        let re = eval(ev, &(-a1.clone() / two_a2))?;
        for (j, c) in roots.iter().enumerate() {
            let other = roots[1 - j].clone();
            let g = together(
                ev,
                &(s.clone()
                    / (rest.clone()
                        * a2.clone().pow(mi.clone())
                        * (x.clone() - other).pow(mi.clone()))),
            )?;
            let im = if real {
                Expr::from(0)
            } else {
                call(ev, Symbol::IM, [c.clone()])?
            };
            found.push(Pole {
                c: c.clone(),
                order: *m,
                g,
                quadratic: Some(monic.clone()),
                real,
                re: re.clone(),
                im,
            });
        }
    }
    // Insertion into giac's order: the real poles ascending, then the
    // complex ones by real and imaginary part.
    let mut sorted: Vec<Pole> = Vec::with_capacity(found.len());
    for pole in found {
        let mut at = sorted.len();
        for (k, other) in sorted.iter().enumerate() {
            if before(ev, &pole, other)? {
                at = k;
                break;
            }
        }
        sorted.insert(at, pole);
    }
    Ok(Some(sorted))
}

/// Whether `p` comes before `q` in giac's order of poles.
fn before(ev: &mut Eval<'_>, p: &Pole, q: &Pole) -> Result<bool, EvalError> {
    match (p.real, q.real) {
        (true, false) => Ok(true),
        (false, true) => Ok(false),
        (true, true) => less(ev, &p.c, &q.c),
        (false, false) => {
            if p.re == q.re {
                less(ev, &p.im, &q.im)
            } else {
                less(ev, &p.re, &q.re)
            }
        }
    }
}

fn less(ev: &mut Eval<'_>, a: &Expr, b: &Expr) -> Result<bool, EvalError> {
    let v = call(ev, Symbol::LESS, [a.clone(), b.clone()])?;
    Ok(v.as_symbol() == Some(Symbol::TRUE))
}

/// The coefficient of `(x − c)^−k` in `r` at the pole (`expansion_coeff`
/// of the partial fractions, :1856), before the test for a real number.
fn laurent(ev: &mut Eval<'_>, pole: &Pole, x: &Expr, k: usize) -> Result<Expr, EvalError> {
    if k > pole.order {
        return Ok(Expr::from(0));
    }
    let j = pole.order - k;
    let mut coefficients = taylor(ev, &pole.g, x, &pole.c, j + 1)?;
    Ok(coefficients.pop().unwrap_or_else(|| Expr::from(0)))
}

/// The first `count` Taylor coefficients of `g` at `c`.
fn taylor(
    ev: &mut Eval<'_>,
    g: &Expr,
    x: &Expr,
    c: &Expr,
    count: usize,
) -> Result<Vec<Expr>, EvalError> {
    let mut out = Vec::with_capacity(count);
    let mut h = g.clone();
    let mut factorial = Expr::from(1);
    for k in 0..count {
        if k > 0 {
            h = derivative(ev, &h, x)?;
            factorial = eval(ev, &(factorial * int(k)))?;
        }
        let at = h.substitute(&[(x.clone(), c.clone())]);
        let v = together(ev, &(at / factorial.clone()))?;
        out.push(simplified(ev, v)?);
    }
    Ok(out)
}

/// The coefficients of `S/T` as a power series, `count` of them.
fn divide_series(
    ev: &mut Eval<'_>,
    s: &[Expr],
    t: &[Expr],
    count: usize,
) -> Result<Vec<Expr>, EvalError> {
    let zero = Expr::from(0);
    let mut out: Vec<Expr> = Vec::with_capacity(count);
    for k in 0..count {
        let mut terms = vec![s.get(k).cloned().unwrap_or_else(|| zero.clone())];
        for i in 1..=k {
            if let Some(ti) = t.get(i) {
                terms.push(-ti.clone() * out[k - i].clone());
            }
        }
        let v = eval(ev, &(Expr::call(Symbol::PLUS, terms) / t[0].clone()))?;
        out.push(v);
    }
    Ok(out)
}

/// The Taylor coefficients of `√g` from those of `g`, `g₀ ≠ 0`.
fn sqrt_series(ev: &mut Eval<'_>, g: &[Expr]) -> Result<Vec<Expr>, EvalError> {
    let Some(g0) = g.first() else {
        return Ok(Vec::new());
    };
    let s0 = sqrt(ev, g0)?;
    // q = √(g/g₀): q₀ = 1, 2 qₖ = hₖ − Σ qᵢ qₖ₋ᵢ.
    let mut q = vec![Expr::from(1)];
    for k in 1..g.len() {
        let mut terms = vec![g[k].clone() / g0.clone()];
        for i in 1..k {
            terms.push(-q[i].clone() * q[k - i].clone());
        }
        let v = together(ev, &(Expr::call(Symbol::PLUS, terms) / Expr::from(2)))?;
        q.push(simplified(ev, v)?);
    }
    let mut out = Vec::with_capacity(q.len());
    for qk in q {
        out.push(eval(ev, &(s0.clone() * qk))?);
    }
    Ok(out)
}

/// `v` when it is a real number, else 0 (`expansion_coeff`, :1862).
fn real(ev: &mut Eval<'_>, v: Expr) -> Result<Expr, EvalError> {
    if v.to_rational().is_some() {
        return Ok(v);
    }
    let n = call(ev, Symbol::N, [v.clone()])?;
    Ok(if real_number(&n) { v } else { Expr::from(0) })
}

/// Whether a number atom is real.
fn real_number(n: &Expr) -> bool {
    use crate::expr::View;
    matches!(
        n.view(),
        View::Int(_) | View::Integer(_) | View::Rational(_) | View::Ball(_)
    )
}

/// `√(N/D)` as `√c · Π fᵢ^⌊mᵢ/2⌋ · √(Π fᵢ^(mᵢ mod 2)) / D` from the
/// factors of `N D`.
fn square_root(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let num = call(ev, Symbol::NUMERATOR, [e.clone()])?;
    let den = call(ev, Symbol::DENOMINATOR, [e.clone()])?;
    let product = call(ev, Symbol::EXPAND, [num * den.clone()])?;
    let list = call(ev, Symbol::FACTOR_LIST, [product])?;
    let mut outside = vec![];
    let mut inside = vec![];
    let mut constant = Expr::from(1);
    for c in list.args().chunks(2) {
        let [f, m] = c else {
            continue;
        };
        let Some(m) = m.as_i64() else {
            continue;
        };
        if f.is_number() {
            constant = Expr::call(Symbol::TIMES, [constant, f.clone().pow(m)]);
            continue;
        }
        outside.push(f.clone().pow(m / 2));
        if m % 2 == 1 {
            inside.push(f.clone());
        }
    }
    let inside = call(ev, Symbol::EXPAND, [Expr::call(Symbol::TIMES, inside)])?;
    let root = sqrt(ev, &constant)? * Expr::call(Symbol::TIMES, outside) * sqrt(ev, &inside)? / den;
    eval(ev, &root)
}

/// `e^(∫…)` (`fullsimp`, :2111): a term `c Log[u]` of `e` as `u^c`.
fn exponential_of(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    let e = call(ev, Symbol::EXPAND, [without_absolute_values(e)])?;
    exponential(ev, &e)
}

/// `∫ f dx`, 0 for 0.
fn integral(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let f = eval(ev, f)?;
    if f == Expr::from(0) {
        return Ok(f);
    }
    let i = call(ev, Symbol::INTEGRATE, [f, x.clone()])?;
    Ok(without_absolute_values(&i))
}

/// The integer value of `e`, if it has one.
fn integer(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<i64>, EvalError> {
    let v = eval(ev, e)?;
    if let Some(n) = v.as_i64() {
        return Ok(Some(n));
    }
    if v.is_number() {
        return Ok(None);
    }
    let v = call(ev, Symbol::EXPAND, [v])?;
    let v = simplified(ev, v)?;
    Ok(v.as_i64())
}

fn is_zero(ev: &mut Eval<'_>, e: &Expr) -> Result<bool, EvalError> {
    if e.is_number() {
        return Ok(*e == Expr::from(0));
    }
    let v = call(ev, Symbol::SIMPLIFY, [e.clone()])?;
    Ok(v == Expr::from(0))
}

/// `e`, simplified when it is a number not yet in canonical form.
fn simplified(ev: &mut Eval<'_>, e: Expr) -> Result<Expr, EvalError> {
    if e.is_number() || e.to_rational().is_some() {
        return Ok(e);
    }
    call(ev, Symbol::SIMPLIFY, [e])
}

/// The coefficients of `f` in `x`, lowest power first.
fn coefficients(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<Vec<Expr>, EvalError> {
    if *f == Expr::from(0) {
        return Ok(Vec::new());
    }
    let list = call(ev, Symbol::COEFFICIENT_LIST, [f.clone(), x.clone()])?;
    if list.head() != Some(&Expr::from(Symbol::LIST)) {
        return Ok(vec![f.clone()]);
    }
    Ok(list.args().iter().rev().cloned().collect())
}

fn degree(ev: &mut Eval<'_>, f: &Expr, x: &Expr) -> Result<Option<i64>, EvalError> {
    Ok(call(ev, Symbol::EXPONENT, [f.clone(), x.clone()])?.as_i64())
}

fn sqrt(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    call(ev, Symbol::SQRT, [e.clone()])
}

fn derivative(ev: &mut Eval<'_>, e: &Expr, x: &Expr) -> Result<Expr, EvalError> {
    let d = call(ev, Symbol::D, [e.clone(), x.clone()])?;
    together(ev, &d)
}

fn together(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    call(ev, Symbol::TOGETHER, [e.clone()])
}

fn eval(ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
    ev.eval(e)
}

fn call<const N: usize>(
    ev: &mut Eval<'_>,
    head: Symbol,
    args: [Expr; N],
) -> Result<Expr, EvalError> {
    ev.eval(&Expr::call(head, args))
}

/// `n` as an integer expression.
fn int(n: usize) -> Expr {
    Expr::from(i64::try_from(n).unwrap_or(i64::MAX))
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

fn binomial(n: usize, k: usize) -> usize {
    (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
}

#[cfg(test)]
mod tests;
