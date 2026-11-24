//! giac's `find_minimum` (CAS-608), after optimization.cc `_find_minimum`
//! and `find_minimum`: the global minimiser of a real function of one
//! variable on `[a, b]`.
//!
//! - **`FindMinimum[f, a, b]`**, `f` an expression in exactly one symbol
//!   that is no builtin, or `Function[{x}, body]`; `a < b` real. Then up to
//!   two more arguments in any order: a positive integer is `maxiter`
//!   (100), a real in `(0, 1]` is `eps` (giac's `epsilon^(2/3)`, `10^−8`
//!   with giac's `epsilon` `10^−12`, C0027). The answer is a ball holding
//!   the minimiser, its width at most `eps` (radius 0 when exact).
//! - The candidate, giac's: Brent's local minimiser without derivatives
//!   (R. P. Brent, *Algorithms for Minimization without Derivatives*,
//!   1973, ch. 5; giac's `local_min_rc`), at most `maxiter` evaluations per
//!   interval of a first-in first-out queue holding `[a, b]` first; an
//!   argument not within `eps` of its interval's ends and better than the
//!   best splits its interval there, both halves queued; an argument at an
//!   end only updates the best. In doubles, `f` from its compiled tree at
//!   the ball midpoint ([`search`], also `frank_wolfe`'s line search).
//! - The certificate, ours: `U` the least upper bound of `f` at the
//!   candidate, `a` and `b` (balls at the working precision and
//!   [`GUARD`] bits); a cover of `[a, b]` by pieces `m ± h`, each bounded
//!   below by the larger of the natural `f(I)` and the mean value form
//!   `f(m) + f′(I)·[−h, h]` when `D[f, x]` compiles; pieces above `U`
//!   dropped, `U` lowered by the values at the centres, the others bisected.
//!   The survivors form clusters of adjacent pieces. A cluster is closed:
//!   inside `(a, b)`, `f` without `Abs`, `Sign`, `HeavisideTheta`, and
//!   `f″ = D[f, x, x]` compiling and positive on its hull `X`, by the
//!   interval Newton step `N(X) = m − f′(m)/f″(X) ⊆ X` (R. Krawczyk, 1969;
//!   R. E. Moore), which proves one critical point in `X`, a strict local
//!   minimum, iterated on `N(X) ∩ X` while it narrows; touching `a` with
//!   `f′(X) > 0`, or `b` with `f′(X) < 0`, by the ball of that end; else
//!   by the cover alone once the hull's width is at most `eps`. Closed
//!   clusters whose value ball lies above another's are dropped.
//! - Ties: several clusters left whose values cannot be separated (as
//!   `sin` on `[0, 4π]`) give the leftmost, with a warning; it is certified
//!   only as a local minimiser whose value is within the value balls of
//!   the minimum.
//! - Limits (C0019): queue pops times `maxiter` at most `terms`, pieces at
//!   most `min(terms/16, 2^16)`; past them, or no cluster closed within the
//!   pieces (no finite bound, or a width above `eps`), a warning and the
//!   call stays. The host's interrupt is polled at each round.
//!
//! Departures from giac:
//!
//! - The answer is a certified ball, not Brent's double (giac:
//!   `find_minimum(cosh(x),-1,1)` is `1.1e-7`; ours `0`).
//! - `f` must compile to ball operations (exp, log, the trigonometric and
//!   hyperbolic functions and their inverses, powers, `Abs`, `Sign`,
//!   `HeavisideTheta`): Airy, Bessel or Gamma leave the call, as does a
//!   non-real value on the cover (giac errs at a non-real candidate value).

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the algorithm: f, a, b, c, d, e, p, q, r, u, v, w, x"
)]

use std::cmp::Ordering;

use cas_num::Integer;
use cas_num::ball::{Ball, Float, Mag, Round};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::minimax::{ball, value};
use crate::native::nintegrate::analytic::Node;
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// Bits beyond the working precision for the certificate.
const GUARD: u64 = 32;
/// The most pieces of the cover.
const MAX_PIECES: usize = 1 << 16;
/// The most interval Newton steps on a cluster.
const NEWTON_STEPS: usize = 64;

/// Why [`search`] stopped without an answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Stop {
    /// More than `max_pops` intervals popped.
    Limits,
    /// `f` was not finite at an argument (NaN or ±∞).
    NotReal,
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: &'static str) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// [`EvalError::Interrupted`] once the evaluation was interrupted.
fn interrupted(ev: &Eval<'_>) -> Result<(), EvalError> {
    if ev
        .context()
        .interrupt_flag()
        .load(std::sync::atomic::Ordering::Relaxed)
    {
        Err(EvalError::Interrupted)
    } else {
        Ok(())
    }
}

/// Brent's local minimiser of `f` on `[a, b]` (Brent 1973, ch. 5, giac's
/// `local_min_rc`): at most `maxiter` evaluations, stopped once the
/// bracket is narrower than `eps` or by Brent's own tolerance; the best
/// argument and its value.
#[expect(clippy::float_cmp, reason = "Brent's tests of which iterates coincide")]
fn brent(
    f: &dyn Fn(f64) -> f64,
    mut a: f64,
    mut b: f64,
    eps: f64,
    maxiter: usize,
) -> Result<(f64, f64), Stop> {
    let finite = |y: f64| {
        if y.is_finite() {
            Ok(y)
        } else {
            Err(Stop::NotReal)
        }
    };
    let c = 0.5 * (3.0 - 5f64.sqrt());
    let root_eps = f64::EPSILON.sqrt();
    let mut x = a + c * (b - a);
    let (mut w, mut v) = (x, x);
    let (mut d, mut e) = (0.0_f64, 0.0_f64);
    let mut fx = finite(f(x))?;
    let (mut fw, mut fv) = (fx, fx);
    for _ in 1..maxiter {
        let m = f64::midpoint(a, b);
        let tol = root_eps * x.abs() + f64::EPSILON / 3.0;
        let tol2 = 2.0 * tol;
        if (x - m).abs() <= tol2 - 0.5 * (b - a) || b - a < eps {
            break;
        }
        let mut golden = true;
        if e.abs() > tol {
            // The parabola through (v, fv), (w, fw), (x, fx).
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            let r = e;
            e = d;
            if p.abs() < (0.5 * q * r).abs() && q * (a - x) < p && p < q * (b - x) {
                d = p / q;
                let u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    d = if x < m { tol } else { -tol };
                }
                golden = false;
            }
        }
        if golden {
            e = if x < m { b - x } else { a - x };
            d = c * e;
        }
        let u = if d.abs() >= tol {
            x + d
        } else if d > 0.0 {
            x + tol
        } else {
            x - tol
        };
        let fu = finite(f(u))?;
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            (v, fv, w, fw, x, fx) = (w, fw, x, fx, u, fu);
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                (v, fv, w, fw) = (w, fw, u, fu);
            } else if fu <= fv || v == x || v == w {
                (v, fv) = (u, fu);
            }
        }
    }
    Ok((x, fx))
}

/// giac's `find_minimum` helper in doubles: Brent's minimiser on `[lo, hi]`
/// (`lo < hi`), and with `global` on each interval of giac's queue, split at
/// an interior argument better than the best; the argument and its value.
/// At most `maxiter` evaluations per interval and `max_pops` intervals
/// ([`Stop::Limits`] past them); [`Stop::NotReal`] where `f` is not finite.
pub(super) fn search(
    f: &dyn Fn(f64) -> f64,
    lo: f64,
    hi: f64,
    eps: f64,
    maxiter: usize,
    global: bool,
    max_pops: usize,
) -> Result<(f64, f64), Stop> {
    let mut queue = std::collections::VecDeque::from([(lo, hi)]);
    let mut best: Option<(f64, f64)> = None;
    let mut pops = 0;
    while let Some((l, r)) = queue.pop_front() {
        pops += 1;
        if pops > max_pops {
            return Err(Stop::Limits);
        }
        let (arg, value) = brent(f, l, r, eps, maxiter)?;
        if !global {
            return Ok((arg, value));
        }
        let better = best.is_none_or(|(_, v)| value < v);
        if arg < l + eps || arg > r - eps {
            if better {
                best = Some((arg, value));
            }
            continue;
        }
        if better {
            best = Some((arg, value));
            queue.push_back((l, arg));
            queue.push_back((arg, r));
        }
    }
    best.ok_or(Stop::Limits)
}

/// `m − r` rounded down, `None` for a ball not finite.
fn lower_end(b: &Ball, prec: u64) -> Option<Float> {
    if !b.is_finite() {
        return None;
    }
    Some(
        b.mid()
            .sub_round(&b.rad().to_float()?, prec, Round::Floor)
            .0,
    )
}

/// `m + r` rounded up, `None` for a ball not finite.
fn upper_end(b: &Ball, prec: u64) -> Option<Float> {
    if !b.is_finite() {
        return None;
    }
    Some(b.mid().add_round(&b.rad().to_float()?, prec, Round::Ceil).0)
}

/// The centre `m`, the half width as `[−h, h]` and the ball `m ± h` of the
/// piece from `l` to `r`, the ends' radii included.
fn interval(l: &Ball, r: &Ball) -> (Ball, Ball, Ball) {
    let mid = l.mid().add(r.mid()).mul_2exp(&Integer::from(-1));
    let half = Mag::from_float(&r.mid().sub(l.mid()).mul_2exp(&Integer::from(-1)));
    let radius = half.add(&std::cmp::max(l.rad().clone(), r.rad().clone()));
    (
        Ball::exact(mid.clone()),
        Ball::new(Float::ZERO, radius.clone()),
        Ball::new(mid, radius),
    )
}

/// A piece of the cover: its ends and the lower bound of `f` on it.
struct Piece {
    l: Ball,
    r: Ball,
    low: Option<Float>,
}

/// `f`, its derivatives and the precision of the certificate.
struct Cover<'a> {
    f: &'a Node,
    df: Option<&'a Node>,
    ddf: Option<&'a Node>,
    prec: u64,
}

impl Cover<'_> {
    /// The piece from `l` to `r`, and the upper bound of `f` at its centre
    /// when the centre is inside it.
    fn piece(&self, l: Ball, r: Ball) -> (Piece, Option<Float>) {
        let prec = self.prec;
        let (centre, s, i) = interval(&l, &r);
        let at_centre = value(self.f, &centre, prec);
        let inside =
            l.cmp(&centre) == Some(Ordering::Less) && centre.cmp(&r) == Some(Ordering::Less);
        let up = at_centre
            .as_ref()
            .filter(|_| inside)
            .and_then(|v| upper_end(v, prec));
        let natural = value(self.f, &i, prec).and_then(|v| lower_end(&v, prec));
        let mean = match (self.df, &at_centre) {
            (Some(df), Some(fc)) => value(df, &i, prec)
                .and_then(|dfi| lower_end(&fc.add(&dfi.mul(&s, prec), prec), prec)),
            _ => None,
        };
        let low = match (natural, mean) {
            (Some(n), Some(m)) => Some(std::cmp::max(n, m)),
            (n, m) => n.or(m),
        };
        (Piece { l, r, low }, up)
    }

    /// The interval Newton steps on `f′` from `x`: the narrowed ball holding
    /// the one critical point of `f` in `x`, a strict local minimum; `None`
    /// when `N(x) ⊆ x` or `f″(x) > 0` fails.
    fn newton(&self, mut x: Ball) -> Option<Ball> {
        let (df, ddf, prec) = (self.df?, self.ddf?, self.prec);
        for step in 0..NEWTON_STEPS {
            let m = Ball::exact(x.mid().clone());
            let d2 = value(ddf, &x, prec)?;
            if d2.sign() != Some(Ordering::Greater) {
                return None;
            }
            let n = m.sub(&value(df, &m, prec)?.div(&d2, prec), prec);
            if step == 0 && !x.contains_ball(&n) {
                return None;
            }
            let next = n.intersection(&x, prec)?;
            let narrower = next.rad() < x.rad();
            x = next;
            if !narrower {
                break;
            }
        }
        Some(x)
    }

    /// The minimiser's ball and the value ball of the cluster `pieces`, or
    /// `None` while it is not closed (the module's rules).
    fn close(
        &self,
        pieces: &[Piece],
        (lo, hi): (&Ball, &Ball),
        width: &Mag,
    ) -> Option<(Ball, Ball)> {
        let prec = self.prec;
        let (first, last) = (pieces.first()?, pieces.last()?);
        let (_, _, x) = interval(&first.l, &last.r);
        let (at_a, at_b) = (first.l == *lo, last.r == *hi);
        let slope = self
            .df
            .and_then(|df| value(df, &x, prec))
            .and_then(|d| d.sign());
        if at_a && slope == Some(Ordering::Greater) {
            return Some((lo.clone(), value(self.f, lo, prec)?));
        }
        if at_b && slope == Some(Ordering::Less) {
            return Some((hi.clone(), value(self.f, hi, prec)?));
        }
        if !at_a && !at_b {
            if let Some(k) = self.newton(x.clone()) {
                if k.rad() <= width {
                    let v = value(self.f, &k, prec)?;
                    return Some((k, v));
                }
            }
        }
        if x.rad() <= width {
            return Some((x.clone(), value(self.f, &x, prec)?));
        }
        None
    }
}

/// The clusters of adjacent pieces, as index ranges.
fn clusters(live: &[Piece]) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0;
    for k in 1..=live.len() {
        if k == live.len() || live[k - 1].r != live[k].l {
            out.push(start..k);
            start = k;
        }
    }
    out
}

/// The variable and body of `f`: a `Function[{x}, body]`, or an expression
/// in one symbol that is no builtin.
fn variable(f: &Expr) -> Option<(Expr, Expr)> {
    if f.head() == Some(&Expr::from(Symbol::FUNCTION)) {
        let [params, body] = f.args() else {
            return None;
        };
        let x = if params.head() == Some(&Expr::from(Symbol::LIST)) {
            match params.args() {
                [x] => x.clone(),
                _ => return None,
            }
        } else {
            params.clone()
        };
        x.as_symbol().filter(|s| !s.is_builtin())?;
        return Some((x, body.clone()));
    }
    let symbols: Vec<Symbol> = f
        .symbols()
        .into_iter()
        .filter(|s| !s.is_builtin())
        .collect();
    match symbols.as_slice() {
        [s] => Some((Expr::from(*s), f.clone())),
        _ => None,
    }
}

/// `maxiter` and `eps` from the arguments after `b`; `None` for anything
/// else.
fn options(given: &[Expr]) -> Option<(Option<Integer>, f64)> {
    let (mut maxiter, mut eps) = (None, 1e-8);
    for o in given {
        if let Some(n) = o.to_integer() {
            if !n.is_positive() {
                return None;
            }
            maxiter = Some(n);
            continue;
        }
        let v = match o.view() {
            View::Ball(b) if b.is_finite() => b.mid().to_f64(),
            _ => Float::from_rational(&o.to_rational()?, 53, Round::Nearest)
                .0
                .to_f64(),
        };
        if v.total_cmp(&0.0) != Ordering::Greater || v.total_cmp(&1.0) == Ordering::Greater {
            return None;
        }
        eps = v;
    }
    Some((maxiter, eps))
}

/// giac's `find_minimum`: the certified ball of the minimiser; `None`
/// where the call stays.
pub(in crate::native) fn find_minimum(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [f, a, b, given @ ..] = e.args() else {
        return Ok(None);
    };
    if given.len() > 2 {
        return Ok(None);
    }
    let (Some((x, body)), Some((maxiter, eps))) = (variable(f), options(given)) else {
        return Ok(None);
    };
    if !body.contains(&x) {
        return Ok(None);
    }
    // The work bounded before any, as a product (C0019).
    let terms = ev.context().limits().terms;
    let maxiter = match maxiter {
        None => 100,
        Some(n) => match n.to_i64().and_then(|k| usize::try_from(k).ok()) {
            Some(k) if k <= terms => k,
            _ => {
                warn(ev, e, "find_minimum: past the limits");
                return Ok(None);
            }
        },
    };
    let max_pops = terms / maxiter;
    let cap = (terms / 16).min(MAX_PIECES);
    let wp = ev.context().working_precision() + GUARD;
    let (Some(lo), Some(hi)) = (ball(ev, a, wp), ball(ev, b, wp)) else {
        return Ok(None);
    };
    if lo.cmp(&hi) != Some(Ordering::Less) {
        return Ok(None);
    }
    let (af, bf) = (lo.mid().to_f64(), hi.mid().to_f64());
    if !(af.is_finite() && bf.is_finite()) || af.total_cmp(&bf) != Ordering::Less {
        return Ok(None);
    }
    let Some(tree) = Node::compile(ev.context(), &body, &x, wp) else {
        return Ok(None);
    };
    // Step 1: giac's candidate, in doubles.
    let found = {
        let flag = ev.context().interrupt_flag();
        let fx = |t: f64| {
            if flag.load(std::sync::atomic::Ordering::Relaxed) {
                return f64::NAN;
            }
            Float::from_f64(t)
                .and_then(|z| value(&tree, &Ball::exact(z), 64))
                .map_or(f64::NAN, |v| v.mid().to_f64())
        };
        search(&fx, af, bf, eps, maxiter, true, max_pops)
    };
    interrupted(ev)?;
    let candidate = match found {
        Ok((arg, _)) => Some(arg),
        Err(Stop::Limits) => {
            warn(ev, e, "find_minimum: past the limits");
            return Ok(None);
        }
        // The cover alone decides.
        Err(Stop::NotReal) => None,
    };
    // Step 2: the cover.
    let df = ev.eval(&Expr::call(Symbol::D, [body.clone(), x.clone()]))?;
    let dtree = Node::compile(ev.context(), &df, &x, wp);
    let smooth = ![Symbol::ABS, Symbol::SIGN, Symbol::HEAVISIDE_THETA]
        .iter()
        .any(|&s| body.may_contain(s));
    let second = match &dtree {
        Some(_) if smooth => {
            let ddf = ev.eval(&Expr::call(Symbol::D, [df.clone(), x.clone()]))?;
            Node::compile(ev.context(), &ddf, &x, wp)
        }
        _ => None,
    };
    let cover = Cover {
        f: &tree,
        df: dtree.as_ref(),
        ddf: second.as_ref(),
        prec: wp,
    };
    let Some((answer, tie)) = certify(ev, &cover, (&lo, &hi), candidate, eps, cap)? else {
        warn(
            ev,
            e,
            "find_minimum: no certified minimiser within the limits",
        );
        return Ok(None);
    };
    if tie {
        warn(
            ev,
            e,
            "find_minimum: the minimum ties within precision; the leftmost minimiser answered",
        );
    }
    Ok(Some(Expr::from(answer)))
}

/// Steps 2 and 3 (the module's rules): the cover of `[lo, hi]` from the
/// upper bound at `candidate` and the ends, its clusters closed; the
/// leftmost minimiser's ball and whether others tie, or `None` past `cap`
/// pieces.
fn certify(
    ev: &Eval<'_>,
    cover: &Cover<'_>,
    (lo, hi): (&Ball, &Ball),
    candidate: Option<f64>,
    eps: f64,
    cap: usize,
) -> Result<Option<(Ball, bool)>, EvalError> {
    let wp = cover.prec;
    let mut best: Option<Float> = None;
    let lower_best = |best: &mut Option<Float>, u: Option<Float>| {
        if let Some(u) = u {
            if best.as_ref().is_none_or(|b| u < *b) {
                *best = Some(u);
            }
        }
    };
    let ends = candidate
        .and_then(Float::from_f64)
        .map(Ball::exact)
        .filter(|c| lo.cmp(c) == Some(Ordering::Less) && c.cmp(hi) == Some(Ordering::Less))
        .into_iter()
        .chain([lo.clone(), hi.clone()]);
    for point in ends {
        lower_best(
            &mut best,
            value(cover.f, &point, wp).and_then(|v| upper_end(&v, wp)),
        );
    }
    let width = Float::from_f64(eps * 0.5).map_or(Mag::ZERO, |h| Mag::from_float_lower(&h));
    let (first, up) = cover.piece(lo.clone(), hi.clone());
    lower_best(&mut best, up);
    let mut live = vec![first];
    let mut pieces = 1;
    loop {
        interrupted(ev)?;
        if let Some(u) = &best {
            live.retain(|p| p.low.as_ref().is_none_or(|l| l <= u));
        }
        // Step 3: close the clusters; bisect the open ones.
        let mut closed = Vec::new();
        let mut open = Vec::new();
        for range in clusters(&live) {
            match cover.close(&live[range.clone()], (lo, hi), &width) {
                Some(answer) => closed.push(answer),
                None => open.push(range),
            }
        }
        for (_, v) in &closed {
            lower_best(&mut best, upper_end(v, wp));
        }
        if open.is_empty() {
            let Some(u) = best.clone() else {
                break;
            };
            let mut left = closed
                .into_iter()
                .filter(|(_, v)| lower_end(v, wp).is_some_and(|l| l <= u));
            let Some((answer, _)) = left.next() else {
                break;
            };
            return Ok(Some((answer, left.next().is_some())));
        }
        let split: usize = open.iter().map(ExactSizeIterator::len).sum();
        if pieces + 2 * split > cap {
            break;
        }
        pieces += 2 * split;
        let mut next = Vec::with_capacity(live.len() + split);
        let mut todo = open.into_iter().peekable();
        for (k, p) in live.into_iter().enumerate() {
            while todo.peek().is_some_and(|r| r.end <= k) {
                todo.next();
            }
            if !todo.peek().is_some_and(|r| r.contains(&k)) {
                next.push(p);
                continue;
            }
            let m = Ball::exact(p.l.mid().add(p.r.mid()).mul_2exp(&Integer::from(-1)));
            for (l, r) in [(p.l, m.clone()), (m, p.r)] {
                let (piece, up) = cover.piece(l, r);
                lower_best(&mut best, up);
                next.push(piece);
            }
        }
        live = next;
    }
    Ok(None)
}
