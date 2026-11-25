//! Tail probabilities and the distribution functions of a law by its name
//! (CAS-DOC-23, CAS-625, CAS-626): Xcas's `UTPN`, `UTPT`, `UTPC`, `UTPF`,
//! `cdf`, `icdf` and `mgf` (giac moyal.cc `_UTPN` :271, `_UTPT` :1420,
//! `_UTPC` :1616, `_UTPF` :1811, `_mgf` :3732, `_cdf` :3807, `_icdf` :3947).
//!
//! - **Upper tails.** `UpperTailNormal` (`UTPN`), `UpperTailStudent`
//!   (`UTPT`), `UpperTailChiSquare` (`UTPC`) and `UpperTailFisher` (`UTPF`)
//!   answer `1 − F`, `F` the law's cumulative builtin at the same point,
//!   evaluated: exact where it is (`UTPN(0, 1, x)` = `(1 − erf(x/√2))/2`,
//!   `UTPC(2, 6)` = `e⁻³`), balls on decimals (giac's are doubles). `UTPN(x)`
//!   is the standard law's, `UTPN(μ, v, x)` takes the variance `v` (σ =
//!   `√v`), `v = 0` staying with giac's Bad Argument Value. giac's checks on
//!   the degrees come first, each leaving the call: `UTPT`'s degree not an
//!   integer (an integral decimal is one) with Bad Argument Value, not
//!   positive with giac's "Invalid dimension"; `UTPC` answers 1 at `x = 0`
//!   and below, stays without a warning for a degree not an integer, and
//!   with Bad Argument Value for one below 1; `UTPF` stays without a warning
//!   for a degree not an integer, answers 1 at `x ≤ 0`, and stays with
//!   "Invalid dimension" for a degree outside `1 … 300`. A symbolic degree
//!   or point gives whatever the cumulative builtin gives (giac's `UTPT(3,
//!   x)` stays).
//! - **A law by its name.** `CDF` (`cdf`), `InverseCDF` (`icdf`) and
//!   `MomentGeneratingFunction` (`mgf`) take a law as giac's
//!   `is_distribution` names it (moyal.cc:3667-3708), through the symbols
//!   those names lower to (C0008): the densities' names (`normald`,
//!   `poisson`, `student`, …), `binomial` and `exp`; `normal` is `Together`,
//!   no law. Its parameters follow, exactly as many as the law has; or the
//!   law is a random variable's object (`randvar`, CAS-624), a discrete
//!   one staying with Bad Argument Value as giac's. `cdf(law, params…, x)`
//!   and `icdf(law, params…, p)` are the law's `…CDF` and `…ICDF` builtins
//!   (the negative binomial's too, where giac has none); without `x` or `p`
//!   they take the symbol `x`. `mgf(law, params…)` is giac's closed form in
//!   the symbol `t` for the normal (its second parameter σ), binomial,
//!   negative binomial, Poisson, gamma, χ², geometric, uniform and
//!   exponential laws, and `Indeterminate` (giac's `undef`) for Student's,
//!   Fisher's, Cauchy's, Weibull's and the beta law.
//! - **A list of values** (rationals or real decimals). `cdf(L, x)` is the
//!   empirical distribution function, `#{Lᵢ ≤ x}/n`, exact, and `cdf(L)`
//!   its steps, the sorted values paired with `(i + 1)/n`. Both depart from
//!   giac, whose dichotomy (vecteur.cc:11196) answers an index: `(#{Lᵢ ≤ x}
//!   − 1)/n`, `−1/n` below the least value, and `[x₍ᵢ₎, i/n]`; ours is the
//!   true distribution function, which giac's own `icdf` inverts. A point
//!   that is not a number stays with Bad Argument Value, as giac's; giac's
//!   matrix of frequency classes is not ported and stays. `icdf(L, p)` is
//!   the `⌈p n⌉`-th least value for `0 < p < 1`, as it is in `L` (giac
//!   prints a decimal), and `Indeterminate` for any other `p` (giac's
//!   `undef`). Values are compared exactly when all are exact, else by
//!   their decimals' `f64` midpoints (C0028).

use std::cmp::Ordering;

use cas_num::ball::Ball;
use cas_num::{Integer, Rational};

use super::refuse;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::random::variates::Law;
use crate::native::statistics::random_variables::{Variable, arity, variable};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `gensizeerr`.
const BAD: &str = "Bad Argument Value";
/// giac's `gendimerr`.
const DIMENSION: &str = "Invalid dimension";
/// The most degrees of freedom of giac's `UTPF`.
const MOST_DEGREES: i64 = 300;

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

/// `1 − F[args…]`, expanded (as `F(b) − F(a)` is) and evaluated.
fn upper(ev: &mut Eval<'_>, cumulative: Symbol, args: Vec<Expr>) -> Native {
    let tail = Expr::call(Symbol::EXPAND, [1 - Expr::call(cumulative, args)]);
    Ok(Some(ev.eval(&tail)?))
}

/// A real number, to compare.
#[derive(Clone, Debug)]
enum Value {
    Exact(Rational),
    Decimal(f64),
}

impl Value {
    /// A rational, or a finite real ball's midpoint; `None` for anything
    /// else.
    fn of(e: &Expr) -> Option<Value> {
        if let Some(q) = e.to_rational() {
            return Some(Value::Exact(q));
        }
        match e.view() {
            View::Ball(b) if b.is_finite() => Some(Value::Decimal(b.mid().to_f64())),
            _ => None,
        }
    }

    fn float(&self) -> f64 {
        match self {
            Value::Exact(q) => Ball::from_rational(q, 64).mid().to_f64(),
            Value::Decimal(f) => *f,
        }
    }

    /// Exactly between rationals, else by `f64` (C0028).
    fn compare(&self, other: &Value) -> Ordering {
        match (self, other) {
            (Value::Exact(a), Value::Exact(b)) => a.cmp(b),
            _ => self.float().total_cmp(&other.float()),
        }
    }

    /// How it compares with 0.
    fn sign(&self) -> Ordering {
        self.compare(&Value::Exact(Rational::ZERO))
    }
}

/// A degree of freedom as giac's `is_integral` reads it.
enum Degree {
    /// An integer, or a decimal of integral value.
    Integer(Integer),
    /// A number that is not an integer.
    Fractional,
    /// Not a number.
    Symbolic,
}

fn degree(n: &Expr) -> Degree {
    if let Some(q) = n.to_rational() {
        return if q.is_integer() {
            Degree::Integer(q.numer().clone())
        } else {
            Degree::Fractional
        };
    }
    match n.view() {
        View::Ball(b) if b.is_finite() => match b.mid().to_rational_within(64) {
            Some(q) if q.is_integer() => Degree::Integer(q.numer().clone()),
            _ => Degree::Fractional,
        },
        _ if n.is_number() => Degree::Fractional,
        _ => Degree::Symbolic,
    }
}

/// `UpperTailNormal[x]`, `UpperTailNormal[μ, v, x]`: `P(X > x)` for the
/// normal law of mean `μ` and variance `v`, the standard one by default.
pub(in crate::native) fn upper_tail_normal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [x] => upper(ev, Symbol::NORMAL_CDF, vec![x.clone()]),
        [mu, v, x] => {
            if Value::of(v).is_some_and(|v| v.sign() == Ordering::Equal) {
                return refuse(ev, e, BAD);
            }
            let sigma = Expr::call(Symbol::SQRT, [v.clone()]);
            upper(ev, Symbol::NORMAL_CDF, vec![mu.clone(), sigma, x.clone()])
        }
        _ => refuse(ev, e, BAD),
    }
}

/// `UpperTailStudent[n, x]`: `P(T > x)` for Student's law with `n`
/// degrees of freedom.
pub(in crate::native) fn upper_tail_student(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, x] = e.args() else {
        return refuse(ev, e, BAD);
    };
    let n = match degree(n) {
        Degree::Fractional => return refuse(ev, e, BAD),
        Degree::Integer(k) if k < Integer::ONE => return refuse(ev, e, DIMENSION),
        Degree::Integer(k) => Expr::from(k),
        Degree::Symbolic => n.clone(),
    };
    upper(ev, Symbol::STUDENT_T_CDF, vec![n, x.clone()])
}

/// `UpperTailChiSquare[n, x]`: `P(X > x)` for the χ² law with `n` degrees
/// of freedom.
pub(in crate::native) fn upper_tail_chi_square(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [n, x] = e.args() else {
        return refuse(ev, e, BAD);
    };
    let place = Value::of(x).map(|x| x.sign());
    if place == Some(Ordering::Equal) {
        return Ok(Some(Expr::from(1)));
    }
    let degree = degree(n);
    if matches!(degree, Degree::Fractional) {
        return Ok(None);
    }
    if place == Some(Ordering::Less) {
        return Ok(Some(Expr::from(1)));
    }
    let n = match degree {
        Degree::Integer(k) if k < Integer::ONE => return refuse(ev, e, BAD),
        Degree::Integer(k) => Expr::from(k),
        Degree::Fractional | Degree::Symbolic => n.clone(),
    };
    upper(ev, Symbol::CHI_SQUARE_CDF, vec![n, x.clone()])
}

/// `UpperTailFisher[m, n, x]`: `P(F > x)` for Fisher's law with `m` and
/// `n` degrees of freedom.
pub(in crate::native) fn upper_tail_fisher(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [m, n, x] = e.args() else {
        return refuse(ev, e, BAD);
    };
    let degrees = [degree(m), degree(n)];
    if degrees.iter().any(|d| matches!(d, Degree::Fractional)) {
        return Ok(None);
    }
    if Value::of(x).is_some_and(|x| x.sign() != Ordering::Greater) {
        return Ok(Some(Expr::from(1)));
    }
    let range = Integer::ONE..=Integer::from(MOST_DEGREES);
    let mut args = Vec::with_capacity(3);
    for (found, given) in degrees.into_iter().zip([m, n]) {
        match found {
            Degree::Integer(count) if !range.contains(&count) => return refuse(ev, e, DIMENSION),
            Degree::Integer(count) => args.push(Expr::from(count)),
            Degree::Fractional | Degree::Symbolic => args.push(given.clone()),
        }
    }
    args.push(x.clone());
    upper(ev, Symbol::FISHER_CDF, args)
}

/// The law a name stands for, giac's `is_distribution` through our
/// symbols: a density's, `binomial`'s, and `exp` for the exponential law.
fn law_named(e: &Expr) -> Option<Law> {
    let s = e.as_symbol()?;
    if s == Symbol::EXP {
        return Some(Law::Exponential);
    }
    Law::of_density(s)
}

/// The law, its parameters and the arguments after them: `(name,
/// params…, rest…)` or `(object, rest…)`; `None` for neither, or a
/// discrete law's object.
fn law_of(args: &[Expr]) -> Option<(Law, Vec<Expr>, &[Expr])> {
    let (first, rest) = args.split_first()?;
    if let Some(law) = law_named(first) {
        let k = arity(law);
        if rest.len() < k {
            return None;
        }
        let (params, rest) = rest.split_at(k);
        return Some((law, params.to_vec(), rest));
    }
    match variable(first)? {
        Variable::Law(law, params) => Some((law, params.to_vec(), rest)),
        Variable::Discrete { .. } => None,
    }
}

/// The law's `function` at the argument after its parameters, or at the
/// symbol `x`.
fn by_law(ev: &mut Eval<'_>, e: &Expr, function: fn(Law) -> Symbol) -> Native {
    let Some((law, mut params, rest)) = law_of(e.args()) else {
        return refuse(ev, e, BAD);
    };
    match rest {
        [] => params.push(Expr::symbol("x")),
        [x] => params.push(x.clone()),
        _ => return refuse(ev, e, BAD),
    }
    Ok(Some(ev.eval(&Expr::call(function(law), params))?))
}

/// A list's values, sorted (stably), each with its expression; `None` when
/// one is not a real number, or there are none.
fn sorted(l: &[Expr]) -> Option<Vec<(Value, &Expr)>> {
    let mut values: Vec<(Value, &Expr)> = l
        .iter()
        .map(|x| Some((Value::of(x)?, x)))
        .collect::<Option<_>>()?;
    if values.is_empty() {
        return None;
    }
    if values.iter().all(|(v, _)| matches!(v, Value::Exact(_))) {
        values.sort_by(|(a, _), (b, _)| a.compare(b));
    } else {
        values.sort_by(|(a, _), (b, _)| a.float().total_cmp(&b.float()));
    }
    Some(values)
}

/// The list's values sorted, or the call stays: giac's classes (a matrix)
/// with no warning, a list past the limits or with a value that is no
/// real number with one.
fn data<'e>(
    ev: &mut Eval<'_>,
    e: &Expr,
    l: &'e [Expr],
) -> Result<Option<Vec<(Value, &'e Expr)>>, EvalError> {
    if l.first().is_some_and(is_list) {
        return Ok(None);
    }
    if l.len() > ev.context().limits().terms {
        return refuse(ev, e, "too many values").map(|_| None);
    }
    let Some(values) = sorted(l) else {
        return refuse(ev, e, BAD).map(|_| None);
    };
    Ok(Some(values))
}

fn fraction(i: usize, n: usize) -> Expr {
    Expr::from(Rational::new(Integer::from(i), Integer::from(n)))
}

/// `CDF[law, params…, x]`, `CDF[law, params…]`, `CDF[object, x]`,
/// `CDF[object]`, `CDF[L, x]`, `CDF[L]`.
pub(in crate::native) fn cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if !args.first().is_some_and(is_list) {
        return by_law(ev, e, Law::cumulative);
    }
    let l = args[0].args();
    match args {
        [_] => {
            let Some(values) = data(ev, e, l)? else {
                return Ok(None);
            };
            let n = values.len();
            let steps = values
                .iter()
                .enumerate()
                .map(|(i, (_, x))| Expr::call(Symbol::LIST, [(*x).clone(), fraction(i + 1, n)]));
            Ok(Some(Expr::call(Symbol::LIST, steps)))
        }
        [_, x] => {
            let Some(values) = data(ev, e, l)? else {
                return Ok(None);
            };
            let Some(x) = Value::of(x) else {
                return refuse(ev, e, BAD);
            };
            let at_most = values
                .iter()
                .filter(|(v, _)| v.compare(&x) != Ordering::Greater)
                .count();
            Ok(Some(fraction(at_most, values.len())))
        }
        _ => refuse(ev, e, BAD),
    }
}

/// `InverseCDF[law, params…, p]`, `InverseCDF[law, params…]`,
/// `InverseCDF[object, p]`, `InverseCDF[object]`, `InverseCDF[L, p]`.
pub(in crate::native) fn inverse_cdf(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if !args.first().is_some_and(is_list) {
        return by_law(ev, e, Law::inverse);
    }
    let [l, p] = args else {
        return refuse(ev, e, BAD);
    };
    let Some(values) = data(ev, e, l.args())? else {
        return Ok(None);
    };
    let inside = Value::of(p).filter(|p| {
        p.sign() == Ordering::Greater && p.compare(&Value::Exact(Rational::ONE)) == Ordering::Less
    });
    let Some(p) = inside else {
        return Ok(Some(Expr::from(Symbol::INDETERMINATE)));
    };
    // The least k in 1 … n with k/n ≥ p: ⌈p n⌉, exactly for an exact p.
    let n = values.len();
    let (mut low, mut high) = (1, n);
    while low < high {
        let middle = low + (high - low) / 2;
        let at = Value::Exact(Rational::new(Integer::from(middle), Integer::from(n)));
        if at.compare(&p) == Ordering::Less {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    Ok(Some(values[low - 1].1.clone()))
}

/// `MomentGeneratingFunction[law, params…]`,
/// `MomentGeneratingFunction[object]`: `E[e^(tX)]` in the symbol `t`.
pub(in crate::native) fn moment_generating_function(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((law, params, [])) = law_of(e.args()) else {
        return refuse(ev, e, BAD);
    };
    let t = Expr::symbol("t");
    let exp = |x: Expr| Expr::call(Symbol::EXP, [x]);
    let et = exp(t.clone());
    let form = match (law, &params[..]) {
        (Law::Normal, [mu, sigma]) => exp(mu * &t + sigma.pow(2) * t.pow(2) / 2),
        (Law::Binomial, [n, p]) => (1 - p + p * &et).pow(n),
        (Law::NegativeBinomial, [n, p]) => (p / (1 - (1 - p) * &et)).pow(n),
        (Law::Poisson, [lambda]) => exp(lambda * (&et - 1)),
        (Law::Gamma, [a, b]) => (1 - &t / b).pow(-a),
        (Law::ChiSquare, [k]) => (1 - 2 * &t).pow(-k / 2),
        (Law::Geometric, [p]) => p * &et / (1 - (1 - p) * &et),
        (Law::Uniform, [a, b]) => (exp(&t * b) - exp(&t * a)) / (&t * (b - a)),
        (Law::Exponential, [lambda]) => 1 / (1 - &t / lambda),
        _ => Expr::from(Symbol::INDETERMINATE),
    };
    Ok(Some(ev.eval(&form)?))
}
