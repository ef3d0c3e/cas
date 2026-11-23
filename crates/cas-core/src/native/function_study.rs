//! Function study (CAS-581, CAS-583): Xcas's `domain` and `convex`, ports
//! of giac's (derive.cc `domain`; optimization.cc `_convex` with vecteur.cc
//! `ldl`, the Bunch–Kaufman factorisation).
//!
//! `FunctionDomain[f, x, mode]` (`domain`; `x` the symbol `x` and `mode` 0
//! by default, `domain(f, 1)` mode 1 in `x`): each function of `x` in `f`,
//! its argument's first, gives a condition (mode 0) or excluded values
//! (mode 1, giac's singular values), giac's table:
//!
//! - a negative integer power `u^(−n)`: the zeros of `u` are excluded (a
//!   product's inverse factors together, giac's one `inv` kernel, their
//!   zeros in `Solve`'s order); a power `u^q`, `q` not an integer (any
//!   root, `x^(1/3)` as giac's, not O7's): `u >= 0`, `u > 0` when `q`
//!   holds `x`; with `q < 0`, the zeros of `u` excluded too; in mode 1 the
//!   zeros of `u`;
//! - `log(u)`, `Ci(u)`: `u > 0` (mode 1: the zeros of `u`); `Ei(u)`: the
//!   zeros of `u`; `acosh(u)`: `u >= 1` (mode 1: the zeros of `u² − 1`);
//!   `asin(u)`, `acos(u)`: `u² <= 1` (mode 1: the zeros of `u² − 1`);
//!   `atanh(u)`, as giac's logarithm `ln((1+u)/(1−u))/2`: the zero of
//!   `1 − u` excluded, `(1+u)/(1−u) > 0` (mode 1: the zeros of `1 + u`);
//!   `tan(u)`: the zeros of `cos(u)`, those in `[0, π/|c|]` when `f` is a
//!   function of `tan(c x)` alone, giac's assumption on a periodic `f`
//!   (derive.cc `domain_auto_assume`; other periodic `f` are not
//!   restricted);
//! - everything else (`sin`, `exp`, `abs`, user functions, …): nothing.
//!
//! The conditions are solved together by `Solve` (`x` when none); the
//! excluded values, each once, join them as giac's: none, `Solve`'s answer
//! (one condition alone); the whole line, `x ≠ e₁ ∧ …`; else the list of
//! `Solve`'s answers, each `∧ x ≠ e` for every `e` it does not exclude.
//! Mode 1 is the list of excluded values. A `Solve` that stays leaves the
//! call. [`restrictions`] gives the conditions and excluded values.
//!
//! `Convex[f, x | {x₁, …}]` and `Convex[f, vars, simplify]` (`convex`):
//! the conditions under which `f` is convex, `True` or `False` when
//! decided. One variable: the sign of `f''`. Several: the Hessian's
//! symbolic Bunch–Kaufman `L D Lᵀ` ([`ldl`], giac's pivoting: a rational
//! entry by its value, `α = (1 + √17)/8` exactly, a symbol's absolute value
//! positive, otherwise the smallest expression), then the sign of each
//! 1×1 block of `D`, and of each 2×2 block's diagonal and determinant. A
//! sign is giac's `determine_sign`: `f` normalised (`Together`, `Simplify`
//! with `simplify`), factored, and its factors' signs: rationals and known
//! constants by value; exponentials, even powers and sums or products of
//! nonnegative terms positive; `u^(−n)` as `u^n`, `u^(p/q)` as
//! `u^⌊|p/q|⌋`; any other a condition `u >= 0` (−1 times: `0 >= u`). With
//! several variables, conditions known nonnegative go, one known negative
//! gives `False`, equal ones once; with `simplify` and one symbol left,
//! `Solve` of them. giac decides a rational Hessian with LAPACK in
//! doubles; here exactly, so a Hessian a rounding makes indefinite is
//! convex here. Variables `u(t)` (giac's with `u'`) leave the call.

#![allow(
    clippy::many_single_char_names,
    reason = "giac's names: a, k, r, p, s for the factorisation, f, x, u, q, b"
)]

use std::cmp::Ordering;

use cas_num::{Integer, Rational};

use super::rounding::real_sign;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The most variables `Convex` factors a Hessian of: each step is
/// symbolic, the factorisation `n³` of them.
const MAX_VARIABLES: usize = 10;

/// `FunctionDomain[f]`, `[f, mode]`, `[f, x]`, `[f, x, mode]`.
pub(super) fn function_domain(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, x, mode) = match e.args() {
        [f] => (f, Expr::symbol("x"), Some(0)),
        [f, k] if k.as_i64().is_some() => (f, Expr::symbol("x"), k.as_i64()),
        [f, x] => (f, x.clone(), Some(0)),
        [f, x, k] => (f, x.clone(), k.as_i64()),
        _ => return Ok(None),
    };
    let Some(mode) = mode.filter(|m| (0..=1).contains(m)) else {
        return Ok(None);
    };
    if x.as_symbol().is_none() {
        return Ok(None);
    }
    domain(ev, f, &x, mode == 1)
}

/// giac's `domain(f, x, mode)`: the domain of `f` in the symbol `x` as
/// giac's conditions, or with `singular` (mode 1) the list of its excluded
/// values; `None` when a `Solve` stays.
pub(crate) fn domain(ev: &mut Eval<'_>, f: &Expr, x: &Expr, singular: bool) -> Native {
    let Some((conditions, excluded)) = restrictions(ev, f, x, singular)? else {
        return Ok(None);
    };
    if singular {
        return Ok(Some(Expr::call(Symbol::LIST, excluded)));
    }
    let mut res = if conditions.is_empty() {
        vec![x.clone()]
    } else {
        let solved = ev.eval(&Expr::call(
            Symbol::SOLVE,
            [Expr::call(Symbol::LIST, conditions), x.clone()],
        ))?;
        if solved.head() != Some(&Expr::from(Symbol::LIST)) {
            return Ok(None);
        }
        solved.args().to_vec()
    };
    if excluded.is_empty() {
        return Ok(Some(match <[Expr; 1]>::try_from(res) {
            Ok([one]) => one,
            Err(res) => Expr::call(Symbol::LIST, res),
        }));
    }
    let unequal: Vec<Expr> = excluded
        .iter()
        .map(|v| Expr::call(Symbol::UNEQUAL, [x.clone(), v.clone()]))
        .collect();
    if res.len() == 1 && res[0] == *x {
        return Ok(Some(match <[Expr; 1]>::try_from(unequal) {
            Ok([one]) => one,
            Err(all) => Expr::call(Symbol::AND, all),
        }));
    }
    // Each excluded value the answer holds is taken out of it.
    for r in &mut res {
        for (v, u) in excluded.iter().zip(&unequal) {
            let at = ev.eval(&r.substitute(&[(x.clone(), v.clone())]))?;
            if at.as_symbol() == Some(Symbol::FALSE)
                || at.to_rational().is_some_and(|q| q.is_zero())
            {
                continue;
            }
            let mut all = if r.head() == Some(&Expr::from(Symbol::AND)) {
                r.args().to_vec()
            } else {
                vec![r.clone()]
            };
            all.push(u.clone());
            *r = Expr::call(Symbol::AND, all);
        }
    }
    Ok(Some(Expr::call(Symbol::LIST, res)))
}

/// The conditions and the excluded values of [`restrictions`].
pub(crate) type Restrictions = (Vec<Expr>, Vec<Expr>);

/// giac's conditions on `x` (mode 0) and excluded values for `f` (the
/// excluded values only with `singular`), each once, in the order found;
/// `None` when a `Solve` of zeros stays.
pub(crate) fn restrictions(
    ev: &mut Eval<'_>,
    f: &Expr,
    x: &Expr,
    singular: bool,
) -> Result<Option<Restrictions>, EvalError> {
    let mut walk = Walk {
        x,
        singular,
        period: tan_period(f, x),
        conditions: Vec::new(),
        excluded: Vec::new(),
    };
    if !walk.visit(ev, f)? {
        return Ok(None);
    }
    Ok(Some((walk.conditions, walk.excluded)))
}

/// `π/|c|` when every `x` in `f` is in `tan(c x)` for one rational `c`:
/// `f` periodic, giac's assumption `0 <= x <= π/|c|`.
fn tan_period(f: &Expr, x: &Expr) -> Option<Expr> {
    let mut found: Option<&Expr> = None;
    for s in f.distinct() {
        if let (Some(Symbol::TAN), [u]) = (s.head().and_then(Expr::as_symbol), s.args())
            && u.contains(x)
        {
            match found {
                Some(t) if t != s => return None,
                _ => found = Some(s),
            }
        }
    }
    let t = found?;
    let c = match t.args() {
        [u] if u == x => Rational::from(1),
        [u] => match (u.head().and_then(Expr::as_symbol), u.args()) {
            (Some(Symbol::TIMES), [c, v]) if v == x => c.to_rational()?,
            _ => return None,
        },
        _ => return None,
    };
    if c.is_zero() || f.substitute(&[(t.clone(), Expr::from(0))]).contains(x) {
        return None;
    }
    Some(Expr::from(Symbol::PI) / Expr::from(c.abs()))
}

/// The walk of `f` for [`restrictions`].
struct Walk<'x> {
    x: &'x Expr,
    singular: bool,
    period: Option<Expr>,
    conditions: Vec<Expr>,
    excluded: Vec<Expr>,
}

impl Walk<'_> {
    fn condition(&mut self, op: Symbol, l: Expr, r: Expr) {
        if self.singular {
            return;
        }
        let c = Expr::call(op, [l, r]);
        if !self.conditions.contains(&c) {
            self.conditions.push(c);
        }
    }

    /// Excludes the real zeros of `u` in `x`; `false` when `Solve` stays.
    fn exclude_zeros(&mut self, ev: &mut Eval<'_>, u: Expr) -> Result<bool, EvalError> {
        self.exclude(ev, u, false)
    }

    /// [`Self::exclude_zeros`], kept within the period when `periodic`.
    fn exclude(&mut self, ev: &mut Eval<'_>, u: Expr, periodic: bool) -> Result<bool, EvalError> {
        let equation = Expr::call(Symbol::EQUATION, [u, Expr::from(0)]);
        let roots = ev.eval(&Expr::call(Symbol::SOLVE, [equation, self.x.clone()]))?;
        if roots.head() != Some(&Expr::from(Symbol::LIST)) {
            return Ok(false);
        }
        for r in roots.args() {
            if periodic && let Some(period) = &self.period {
                let beyond = ev.eval(&(period.clone() - r.clone()))?;
                let ctx = ev.context();
                if real_sign(ctx, r) == Some(Ordering::Less)
                    || real_sign(ctx, &beyond) == Some(Ordering::Less)
                {
                    continue;
                }
            }
            if !self.excluded.contains(r) {
                self.excluded.push(r.clone());
            }
        }
        Ok(true)
    }

    /// The restrictions of `e`'s functions of `x`, arguments first;
    /// `false` when a `Solve` stays.
    fn visit(&mut self, ev: &mut Eval<'_>, e: &Expr) -> Result<bool, EvalError> {
        let x = self.x;
        if e == x || !e.contains(x) {
            return Ok(true);
        }
        let head = e.head().and_then(Expr::as_symbol);
        match (head, e.args()) {
            (Some(Symbol::TIMES), factors) => {
                // The inverse factors are one denominator, as giac's `inv`.
                let mut denominator = Vec::new();
                for g in factors {
                    match inverse(g, x) {
                        Some(power) => {
                            if !self.visit(ev, &g.args()[0])? {
                                return Ok(false);
                            }
                            denominator.push(power);
                        }
                        None => {
                            if !self.visit(ev, g)? {
                                return Ok(false);
                            }
                        }
                    }
                }
                if !denominator.is_empty() {
                    return self.exclude_zeros(ev, Expr::call(Symbol::TIMES, denominator));
                }
                Ok(true)
            }
            (Some(Symbol::POWER), [b, q]) => self.power(ev, b, q),
            (_, args) => {
                for a in args {
                    if !self.visit(ev, a)? {
                        return Ok(false);
                    }
                }
                match head {
                    Some(f) => self.function(ev, f, args),
                    None => Ok(true),
                }
            }
        }
    }

    /// giac's rules for `b^q`.
    fn power(&mut self, ev: &mut Eval<'_>, b: &Expr, q: &Expr) -> Result<bool, EvalError> {
        if !self.visit(ev, b)? || !self.visit(ev, q)? {
            return Ok(false);
        }
        if !b.contains(self.x) {
            return Ok(true);
        }
        let zero = || Expr::from(0);
        match q.to_rational() {
            Some(q) if q.is_integer() => {
                if q.is_negative() {
                    return self.exclude_zeros(ev, b.clone());
                }
                Ok(true)
            }
            Some(q) => {
                self.condition(Symbol::GREATER_EQUAL, b.clone(), zero());
                if q.is_negative() || self.singular {
                    return self.exclude_zeros(ev, b.clone());
                }
                Ok(true)
            }
            None => {
                if self.singular {
                    return self.exclude_zeros(ev, b.clone());
                }
                let positive =
                    !q.contains(self.x) && real_sign(ev.context(), q) == Some(Ordering::Greater);
                let op = if positive {
                    Symbol::GREATER_EQUAL
                } else {
                    Symbol::GREATER
                };
                self.condition(op, b.clone(), zero());
                Ok(true)
            }
        }
    }

    /// giac's rules for the function `f` of `args`.
    fn function(&mut self, ev: &mut Eval<'_>, f: Symbol, args: &[Expr]) -> Result<bool, EvalError> {
        let Some(u) = args.last().filter(|u| u.contains(self.x)) else {
            return Ok(true);
        };
        let u = u.clone();
        let zero = || Expr::from(0);
        let one = || Expr::from(1);
        match f {
            Symbol::EXP_INTEGRAL_EI => self.exclude_zeros(ev, u),
            Symbol::LOG | Symbol::COS_INTEGRAL => {
                if self.singular {
                    return self.exclude_zeros(ev, u);
                }
                self.condition(Symbol::GREATER, u, zero());
                Ok(true)
            }
            Symbol::ARC_COSH => {
                if self.singular {
                    return self.exclude_zeros(ev, u.pow(2) - one());
                }
                self.condition(Symbol::GREATER_EQUAL, u, one());
                Ok(true)
            }
            Symbol::ARC_SIN | Symbol::ARC_COS => {
                if self.singular {
                    return self.exclude_zeros(ev, u.pow(2) - one());
                }
                self.condition(Symbol::LESS_EQUAL, u.pow(2), one());
                Ok(true)
            }
            Symbol::ARC_TANH => {
                // giac's ln((1+u)/(1-u))/2: its inverse, then its logarithm.
                if !self.exclude_zeros(ev, one() - u.clone())? {
                    return Ok(false);
                }
                if self.singular {
                    return self.exclude_zeros(ev, one() + u);
                }
                let ratio = (one() + u.clone()) / (one() - u);
                self.condition(Symbol::GREATER, ratio, zero());
                Ok(true)
            }
            Symbol::TAN => self.exclude(ev, Expr::call(Symbol::COS, [u]), true),
            _ => Ok(true),
        }
    }
}

/// `b^n` for a factor `g = b^(−n)` of a product, `n` a positive integer
/// and `b` in `x`.
fn inverse(g: &Expr, x: &Expr) -> Option<Expr> {
    match (g.head().and_then(Expr::as_symbol), g.args()) {
        (Some(Symbol::POWER), [b, q]) if b.contains(x) => {
            let q = q.to_rational()?;
            (q.is_integer() && q.is_negative()).then(|| b.pow(Expr::from(-q)))
        }
        _ => None,
    }
}

/// `Convex[f, vars]` and `Convex[f, vars, simplify]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's _convex in order: arguments, one variable, the blocks, the conditions"
)]
pub(super) fn convex(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (f, vars, simplify) = match e.args() {
        [f, v] => (f, v, false),
        [f, v, o] if o.as_symbol() == Some(Symbol::SIMPLIFY) => (f, v, true),
        _ => return Ok(None),
    };
    let given = if vars.head() == Some(&Expr::from(Symbol::LIST)) {
        vars.args()
    } else {
        std::slice::from_ref(vars)
    };
    let mut xs: Vec<Expr> = Vec::new();
    for x in given {
        if x.as_symbol().is_none() {
            return Ok(None);
        }
        if !xs.contains(x) {
            xs.push(x.clone());
        }
    }
    let n = xs.len();
    if n == 0 || n > MAX_VARIABLES {
        return Ok(None);
    }
    let h = ev.eval(&Expr::call(
        Symbol::HESSIAN,
        [f.clone(), Expr::call(Symbol::LIST, xs)],
    ))?;
    let Some(mut a) = matrix(&h, n) else {
        return Ok(None);
    };
    let signs = Signs { simplify };
    if n == 1 {
        let sign = signs.determine(ev, &a[0][0])?;
        match sign {
            (1, ref atoms) if atoms.is_empty() => return Ok(Some(Expr::from(Symbol::TRUE))),
            (-1, ref atoms) if atoms.is_empty() => return Ok(Some(Expr::from(Symbol::FALSE))),
            _ => {}
        }
        let s = strip(ev, sign)?;
        let s = signs.simp(ev, &s)?;
        let condition = if negative_coefficient(&s) {
            Expr::call(Symbol::GREATER_EQUAL, [Expr::from(0), ev.eval(&-s)?])
        } else {
            Expr::call(Symbol::GREATER_EQUAL, [s, Expr::from(0)])
        };
        return ev.eval(&condition).map(Some);
    }
    let Some(sub) = ldl(ev, &mut a)? else {
        return Ok(None);
    };
    let mut cond = Vec::new();
    let mut i = 0;
    while i < n {
        if i + 1 < n && !sub[i + 1].to_rational().is_some_and(|q| q.is_zero()) {
            let (a11, a22, a12) = (a[i][i].clone(), a[i + 1][i + 1].clone(), sub[i + 1].clone());
            cond.push(signs.determine(ev, &a11)?);
            cond.push(signs.determine(ev, &a22)?);
            let det = ev.eval(&(a11 * a22 - a12.pow(2)))?;
            cond.push(signs.determine(ev, &det)?);
            i += 2;
        } else {
            cond.push(signs.determine(ev, &a[i][i])?);
            i += 1;
        }
    }
    let mut res: Vec<Expr> = Vec::new();
    for (c, atoms) in cond {
        if atoms.is_empty() && c == -1 {
            return Ok(Some(Expr::from(Symbol::FALSE)));
        }
        if (atoms.is_empty() && c == 1) || c == 0 {
            continue;
        }
        res.push(strip(ev, (c, atoms))?);
    }
    // Equal conditions once.
    let mut i = 0;
    while i < res.len() {
        let mut j = res.len();
        while j > i + 1 {
            j -= 1;
            let difference = signs.simp(ev, &(res[j].clone() - res[i].clone()))?;
            if difference.to_rational().is_some_and(|q| q.is_zero()) {
                res.remove(j);
            }
        }
        i += 1;
    }
    let mut out = Vec::new();
    for r in res {
        let r = signs.simp(ev, &r)?;
        if nonnegative(ev, &r)?.is_some() {
            continue;
        }
        let minus = ev.eval(&-r.clone())?;
        if nonnegative(ev, &minus)? == Some(true) {
            return Ok(Some(Expr::from(Symbol::FALSE)));
        }
        out.push(Expr::call(Symbol::GREATER_EQUAL, [r, Expr::from(0)]));
    }
    if out.is_empty() {
        return Ok(Some(Expr::from(Symbol::TRUE)));
    }
    if simplify {
        let mut names: Vec<Symbol> = Vec::new();
        for c in &out {
            for s in c.symbols() {
                if !s.is_builtin() && !names.contains(&s) {
                    names.push(s);
                }
            }
        }
        if let [s] = names.as_slice() {
            let solved = ev.eval(&Expr::call(
                Symbol::SOLVE,
                [Expr::call(Symbol::LIST, out.clone()), Expr::from(*s)],
            ))?;
            if solved.head() == Some(&Expr::from(Symbol::LIST)) && !solved.args().is_empty() {
                return Ok(Some(solved));
            }
        }
    }
    ev.eval(&Expr::call(Symbol::LIST, out)).map(Some)
}

/// The rows of an `n × n` matrix `h`.
fn matrix(h: &Expr, n: usize) -> Option<Vec<Vec<Expr>>> {
    let list = Expr::from(Symbol::LIST);
    if h.head() != Some(&list) || h.args().len() != n {
        return None;
    }
    h.args()
        .iter()
        .map(|row| {
            (row.head() == Some(&list) && row.args().len() == n).then(|| row.args().to_vec())
        })
        .collect()
}

/// Whether `s` is a product with a negative rational coefficient, or a
/// negative rational: giac's `neg`.
fn negative_coefficient(s: &Expr) -> bool {
    if let Some(q) = s.to_rational() {
        return q.is_negative();
    }
    matches!(
        (s.head().and_then(Expr::as_symbol), s.args().first().and_then(Expr::to_rational)),
        (Some(Symbol::TIMES), Some(q)) if q.is_negative()
    )
}

/// A sign as giac's `_sign` leaves it: −1, 0 or 1 times the product of the
/// signs of the factors it cannot decide.
type Sign = (i8, Vec<Expr>);

/// giac's `determine_sign` with `ratnormal`, or `simplify`.
struct Signs {
    simplify: bool,
}

impl Signs {
    fn simp(&self, ev: &mut Eval<'_>, e: &Expr) -> Result<Expr, EvalError> {
        let f = if self.simplify {
            Symbol::SIMPLIFY
        } else {
            Symbol::TOGETHER
        };
        ev.eval(&Expr::call(f, [e.clone()]))
    }

    /// The sign of `e`, normalised and factored.
    fn determine(&self, ev: &mut Eval<'_>, e: &Expr) -> Result<Sign, EvalError> {
        let s = self.simp(ev, e)?;
        let g = if s.to_rational().is_some() {
            s
        } else {
            ev.eval(&Expr::call(Symbol::FACTOR, [s]))?
        };
        sign_parts(ev, &g)
    }
}

/// The sign of a factored `g`: a product's, its factors'; `u^(−n)` as
/// `u^n`, an even power positive, `u^(p/q)` as `u^⌊|p/q|⌋`; otherwise
/// by [`nonnegative`] and `real_sign`, or `g` itself undecided.
fn sign_parts(ev: &mut Eval<'_>, g: &Expr) -> Result<Sign, EvalError> {
    if let Some(q) = g.to_rational() {
        return Ok((ordering(q.sign()), Vec::new()));
    }
    match (g.head().and_then(Expr::as_symbol), g.args()) {
        (Some(Symbol::TIMES), factors) => {
            let mut c = 1;
            let mut atoms = Vec::new();
            for f in factors {
                let (fc, fa) = sign_parts(ev, f)?;
                c *= fc;
                if c == 0 {
                    return Ok((0, Vec::new()));
                }
                atoms.extend(fa);
            }
            return Ok((c, atoms));
        }
        (Some(Symbol::POWER), [b, _]) if b.as_symbol() == Some(Symbol::E) => {
            return Ok((1, Vec::new()));
        }
        (Some(Symbol::POWER), [b, q]) => {
            if let Some(q) = q.to_rational() {
                if !q.is_integer() {
                    let m = q.abs().trunc();
                    if m.is_zero() {
                        return Ok((1, Vec::new()));
                    }
                    if m == Integer::ONE {
                        return sign_parts(ev, b);
                    }
                    return sign_parts(ev, &b.pow(Expr::from(m)));
                }
                if !q.numer().unsigned_abs().bit(0) {
                    return Ok((1, Vec::new()));
                }
                if q.is_negative() {
                    return sign_parts(ev, &b.pow(Expr::from(-q)));
                }
                if q.is_one() {
                    return sign_parts(ev, b);
                }
                if nonnegative(ev, b)? == Some(true) {
                    return Ok((1, Vec::new()));
                }
                return Ok((1, vec![g.clone()]));
            }
        }
        _ => {}
    }
    if nonnegative(ev, g)? == Some(true) {
        return Ok((1, Vec::new()));
    }
    let minus = ev.eval(&-g.clone())?;
    if nonnegative(ev, &minus)? == Some(true) {
        return Ok((-1, Vec::new()));
    }
    Ok((1, vec![g.clone()]))
}

/// −1, 0 or 1.
fn ordering(o: Ordering) -> i8 {
    match o {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// giac's `strip_sign`: the sign's coefficient times its undecided factors.
fn strip(ev: &mut Eval<'_>, (c, atoms): Sign) -> Result<Expr, EvalError> {
    let mut factors = vec![Expr::from(i64::from(c))];
    factors.extend(atoms);
    ev.eval(&Expr::call(Symbol::TIMES, factors))
}

/// Whether `e` is known positive (`Some(true)`) or nonnegative
/// (`Some(false)`) for real values of its symbols, giac's
/// `is_positive_safe`: rationals and constants by value; exponentials,
/// absolute values, even powers and roots; sums and products of such.
fn nonnegative(ev: &mut Eval<'_>, e: &Expr) -> Result<Option<bool>, EvalError> {
    if let Some(q) = e.to_rational() {
        return Ok(match q.sign() {
            Ordering::Greater => Some(true),
            Ordering::Equal => Some(false),
            Ordering::Less => None,
        });
    }
    if e.is_numeric() {
        return Ok(match real_sign(ev.context(), e) {
            Some(Ordering::Greater) => Some(true),
            Some(Ordering::Equal) => Some(false),
            _ => None,
        });
    }
    Ok(match (e.head().and_then(Expr::as_symbol), e.args()) {
        (Some(Symbol::EXP), [_]) => Some(true),
        (Some(Symbol::ABS), [_]) => Some(false),
        (Some(Symbol::POWER), [b, _]) if b.as_symbol() == Some(Symbol::E) => Some(true),
        (Some(Symbol::POWER), [b, q]) => {
            let base = nonnegative(ev, b)?;
            match q.to_rational() {
                Some(q) if q.is_integer() && !q.numer().unsigned_abs().bit(0) => {
                    Some(base == Some(true))
                }
                _ if base == Some(true) => Some(true),
                Some(q) if !q.is_integer() && !q.denom().unsigned_abs().bit(0) => Some(false),
                _ => None,
            }
        }
        (Some(Symbol::TIMES | Symbol::PLUS), terms) => {
            let mut strict = Vec::with_capacity(terms.len());
            for t in terms {
                match nonnegative(ev, t)? {
                    Some(s) => strict.push(s),
                    None => return Ok(None),
                }
            }
            if e.head().and_then(Expr::as_symbol) == Some(Symbol::TIMES) {
                Some(strict.iter().all(|s| *s))
            } else {
                Some(strict.iter().any(|s| *s))
            }
        }
        _ => None,
    })
}

/// How the pivot search sees an entry: giac's `_SYMB` or not.
enum Kind {
    Number(Rational),
    Symbol,
    Compound,
}

fn kind(e: &Expr) -> Kind {
    if let Some(q) = e.to_rational() {
        Kind::Number(q)
    } else if e.as_symbol().is_some() {
        Kind::Symbol
    } else {
        Kind::Compound
    }
}

/// A modulus in the search: a rational's, or symbolic (`None`), which giac
/// holds positive.
type Modulus = Option<Rational>;

/// giac's `is_strictly_greater(abs(entry), m)` for an entry no compound.
fn strictly_greater(entry: &Kind, m: &Modulus) -> bool {
    match (entry, m) {
        (Kind::Number(q), Some(m)) => q.abs() > *m,
        // giac holds `abs(z) > 0` for a symbol `z`.
        (Kind::Symbol, Some(m)) => m.is_zero(),
        _ => false,
    }
}

fn strictly_positive(m: &Modulus) -> bool {
    m.as_ref().is_none_or(Rational::is_positive)
}

/// `m >= α l` for `α = (1 + √17)/8`, exactly, `m, l >= 0`.
fn at_least_alpha(m: &Rational, l: &Rational) -> bool {
    let t = &(&Rational::from(8) * m) - l;
    !t.is_negative() && &t * &t >= &(&Rational::from(17) * l) * l
}

/// Swaps rows and columns `k` and `l` of the lower triangle of `a`, and
/// in `perm` (vecteur.cc `ldl_swap_rc`).
#[allow(
    clippy::needless_range_loop,
    reason = "giac's ldl_swap_rc, its indices kept"
)]
fn swap_rc(a: &mut [Vec<Expr>], k: usize, l: usize, perm: &mut [usize]) {
    let (lo, hi) = (k.min(l), k.max(l));
    for i in lo + 1..hi {
        let t = a[hi][i].clone();
        a[hi][i] = a[i][lo].clone();
        a[i][lo] = t;
    }
    for i in 0..lo {
        let t = a[lo][i].clone();
        a[lo][i] = a[hi][i].clone();
        a[hi][i] = t;
    }
    for row in a.iter_mut().skip(hi + 1) {
        row.swap(lo, hi);
    }
    let t = a[lo][lo].clone();
    a[lo][lo] = a[hi][hi].clone();
    a[hi][hi] = t;
    perm.swap(k, l);
}

/// giac's symbolic Bunch–Kaufman `L D Lᵀ` of the symmetric `a` (vecteur.cc
/// `ldl`, its branch without LAPACK), in place on the lower triangle: `D`'s
/// diagonal on `a`'s; the returned list holds at `j + 1` the
/// sub-diagonal of a 2×2 block at `j`, 0 elsewhere. `None` for a singular
/// 2×2 pivot.
#[allow(
    clippy::too_many_lines,
    clippy::needless_range_loop,
    reason = "giac's ldl step by step, its indices kept"
)]
fn ldl(ev: &mut Eval<'_>, a: &mut [Vec<Expr>]) -> Result<Option<Vec<Expr>>, EvalError> {
    let n = a.len();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut ws1 = vec![Expr::from(0); n];
    let mut ws2 = vec![Expr::from(0); n];
    if n == 1 {
        return Ok(Some(ws1));
    }
    let mut k = 0;
    while k < n {
        let akk = kind(&a[k][k]);
        // The largest entry below the diagonal, or the smallest compound.
        let mut lambda: Modulus = Some(Rational::from(0));
        let mut symbolic = false;
        let mut smallest = u32::MAX;
        let mut r = k;
        for i in k + 1..n {
            let aik = kind(&a[i][k]);
            if let Kind::Compound = aik {
                if strictly_positive(&lambda) {
                    continue;
                }
                let size = a[i][k].size();
                if size < smallest {
                    smallest = size;
                    r = i;
                    symbolic = true;
                }
            } else if strictly_greater(&aik, &lambda) {
                r = i;
                lambda = match aik {
                    Kind::Number(q) => Some(q.abs()),
                    _ => None,
                };
            }
        }
        if lambda.as_ref().is_some_and(Rational::is_zero) {
            if !symbolic {
                // No pivoting.
                if k + 1 < n {
                    ws1[k + 1] = Expr::from(0);
                }
                k += 1;
                continue;
            }
            lambda = None;
        }
        let akk_modulus = match &akk {
            Kind::Number(q) => Some(q.abs()),
            _ => None,
        };
        let (s, p) = match (&akk_modulus, &lambda) {
            (Some(m), Some(l)) if at_least_alpha(m, l) => {
                r = k;
                (1, k)
            }
            _ => {
                // The largest entry of column r off the diagonal.
                let mut sigma: Modulus = Some(Rational::from(0));
                let mut symbolic = false;
                let mut smallest = u32::MAX;
                let mut p = k;
                for i in k..n {
                    if i == r {
                        continue;
                    }
                    let air = if i < r { &a[r][i] } else { &a[i][r] };
                    let kair = kind(air);
                    if let Kind::Compound = kair {
                        if strictly_positive(&sigma) {
                            continue;
                        }
                        let size = air.size();
                        if size < smallest {
                            smallest = size;
                            p = i;
                            symbolic = true;
                        }
                    } else if strictly_greater(&kair, &sigma) {
                        p = i;
                        sigma = match kair {
                            Kind::Number(q) => Some(q.abs()),
                            _ => None,
                        };
                    }
                }
                if sigma.as_ref().is_some_and(Rational::is_zero) && symbolic {
                    sigma = None;
                }
                let arr = match kind(&a[r][r]) {
                    Kind::Number(q) => Some(q.abs()),
                    _ => None,
                };
                match (&akk_modulus, &lambda, &sigma, &arr) {
                    (Some(m), Some(l), Some(s), _) if at_least_alpha(&(s * m), &(l * l)) => {
                        r = k;
                        (1, p)
                    }
                    (_, _, Some(s), Some(m)) if at_least_alpha(m, s) => (1, p),
                    _ => {
                        if p > r {
                            std::mem::swap(&mut p, &mut r);
                        }
                        (2, p)
                    }
                }
            }
        };
        if s == 1 {
            if r != k {
                swap_rc(a, r, k, &mut perm);
            }
            let akk = a[k][k].clone();
            for i in k + 1..n {
                ws1[i] = ev.eval(&(a[i][k].clone() / akk.clone()))?;
                for j in k + 1..=i {
                    a[i][j] = ev.eval(&(a[i][j].clone() - ws1[i].clone() * a[j][k].clone()))?;
                }
            }
            for i in k + 1..n {
                a[i][k] = ws1[i].clone();
            }
            if k + 1 < n {
                ws1[k + 1] = Expr::from(0);
            }
            k += 1;
        } else {
            if p == k && r != k + 1 {
                swap_rc(a, r, k + 1, &mut perm);
            } else if p == k + 1 {
                swap_rc(a, r, k, &mut perm);
            } else if p > k + 1 {
                swap_rc(a, r, k, &mut perm);
                swap_rc(a, p, k + 1, &mut perm);
            }
            let (e11, e21, e22) = (
                a[k][k].clone(),
                a[k + 1][k].clone(),
                a[k + 1][k + 1].clone(),
            );
            let d = ev.eval(&(e11.clone() * e22.clone() - e21.clone() * e21.clone()))?;
            if k + 2 < n && d.to_rational().is_some_and(|q| q.is_zero()) {
                return Ok(None);
            }
            for i in k + 2..n {
                let (aik, aik1) = (a[i][k].clone(), a[i][k + 1].clone());
                ws1[i] = ev.eval(
                    &((aik.clone() * e22.clone() - aik1.clone() * e21.clone()) / d.clone()),
                )?;
                ws2[i] = ev.eval(&((aik1 * e11.clone() - aik * e21.clone()) / d.clone()))?;
                for j in k + 2..=i {
                    let update =
                        a[j][k].clone() * ws1[i].clone() + a[j][k + 1].clone() * ws2[i].clone();
                    a[i][j] = ev.eval(&(a[i][j].clone() - update))?;
                }
            }
            for i in k + 2..n {
                a[i][k] = ws1[i].clone();
                a[i][k + 1] = ws2[i].clone();
            }
            ws1[k + 1] = e21;
            a[k + 1][k] = Expr::from(0);
            if k + 2 < n {
                ws1[k + 2] = Expr::from(0);
            }
            k += 2;
        }
    }
    Ok(Some(ws1))
}

#[cfg(test)]
mod tests;
