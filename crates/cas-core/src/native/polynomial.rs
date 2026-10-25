//! The polynomial builtins of phase 1 (CAS-DOC-8 §10, CAS-140), Xcas's
//! `expand`, `normal`, `gcd`, `lcm`, `quo`, `rem`, `quorem`, `degree`,
//! `coeff`, `lcoeff`, `numer`, `denom`, `content`, `primpart`, `resultant`
//! and `discriminant`, and the factoring ones (CAS-185), `factor` and
//! `collect`, `factors` and `sqrfree`, over the conversion of
//! [`crate::poly`] and cas-poly's algorithms.
//!
//! - The expressions go to rational functions of their generators
//!   (CAS-DOC-8 §2), or for `Expand` to polynomials term by term
//!   ([`poly::expand`]), and the results come back through evaluation.
//! - The variable is the one given, or the one generator that contains a
//!   symbol other than a constant (§15 detail 6); with several, a warning
//!   names them and the call stays. With a variable, the input must be a
//!   polynomial in it: no other generator, and no denominator, may contain
//!   it.
//! - Arithmetic (`Expand`, `Exponent`, `Coefficient`, `LeadingCoefficient`)
//!   works with balls, each a generator that evaluation multiplies out; the
//!   exact algorithms leave the call with a warning (§9).
//! - Past `Limits::terms`, estimated before each product, the call stays,
//!   with a warning (§15 detail 3). The interrupt stops the long loops.
//! - Each call takes one step, for the trace: `Expand` and `Together` claim
//!   identities, which `Trace::verify` checks; the others are formal.

use std::borrow::Cow;

use cas_num::algebra::{EuclideanDomain, GaussianIntegers, Integers, Ring};
use cas_num::{Complex, Integer, Natural, Rational};
use cas_poly::{
    Cancel, Cancelled, CharacteristicZero, DenseRing, FactorError, Factors, GcdDomain, Sparse,
    SparseRing,
};

use crate::attributes::Attributes;
use crate::context::{Context, EvalError};
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::poly::{self, Bounds, Coefficient, Fraction, PolyError, Polynomials};
use crate::symbol::Symbol;
use crate::trace::Warning;

type Native = Result<Option<Expr>, EvalError>;

/// A coefficient ring of the conversion, ℤ or ℤ\[i\], with what the
/// builtins need of it.
pub(super) trait Numbers: GcdDomain<Elem: Coefficient> + Default {}

impl Numbers for Integers {}
impl Numbers for GaussianIntegers {}

/// Runs `$body` with `$ring`, `$num` and `$den` bound to the ring and the
/// parts of `$fraction`, over ℤ or ℤ\[i\].
macro_rules! with_fraction {
    ($fraction:expr, $nvars:expr, |$ring:ident, $num:ident, $den:ident| $body:expr) => {
        match $fraction {
            Fraction::Integers {
                num: $num,
                den: $den,
            } => {
                let $ring = SparseRing::new(Integers, $nvars);
                $body
            }
            Fraction::GaussianIntegers {
                num: $num,
                den: $den,
            } => {
                let $ring = SparseRing::new(GaussianIntegers, $nvars);
                $body
            }
        }
    };
}

/// As [`with_fraction`], for two fractions of one conversion, which share
/// their ring.
macro_rules! with_fractions {
    ($a:expr, $b:expr, $nvars:expr, |$ring:ident, $na:ident, $da:ident, $nb:ident, $db:ident| $body:expr) => {
        match ($a, $b) {
            (
                Fraction::Integers { num: $na, den: $da },
                Fraction::Integers { num: $nb, den: $db },
            ) => {
                let $ring = SparseRing::new(Integers, $nvars);
                $body
            }
            (
                Fraction::GaussianIntegers { num: $na, den: $da },
                Fraction::GaussianIntegers { num: $nb, den: $db },
            ) => {
                let $ring = SparseRing::new(GaussianIntegers, $nvars);
                $body
            }
            _ => unreachable!("one conversion, one ring"),
        }
    };
}

fn warn(ev: &mut Eval<'_>, message: impl Into<Cow<'static, str>>, e: &Expr) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// The interrupt, as cas-poly polls it.
fn interrupted(Cancelled: Cancelled) -> EvalError {
    EvalError::Interrupted
}

/// `args` as rational functions of their generators, or `None` after a
/// warning: balls are generators when `balls`, else refused.
pub(super) fn convert(
    ev: &mut Eval<'_>,
    e: &Expr,
    args: &[Expr],
    balls: bool,
) -> Result<Option<Polynomials>, EvalError> {
    let ctx = ev.context();
    let cancel: &dyn Cancel = ctx.interrupt_flag();
    let result = if balls {
        Polynomials::with_balls(args, bounds(ctx), cancel)
    } else {
        Polynomials::new(args, bounds(ctx), cancel)
    };
    refused(ev, e, result)
}

/// The bounds of a conversion in `ctx`: `Limits::terms` and
/// `Limits::number_bits` (CAS-156).
fn bounds(ctx: &Context) -> Bounds {
    let limits = ctx.limits();
    Bounds {
        terms: limits.terms,
        bits: limits.number_bits,
    }
}

/// The value of a conversion, or `None` after a warning on `e` that says
/// why there is none.
fn refused<T>(
    ev: &mut Eval<'_>,
    e: &Expr,
    result: Result<T, PolyError>,
) -> Result<Option<T>, EvalError> {
    match result {
        Ok(p) => Ok(Some(p)),
        Err(PolyError::Interrupted) => Err(EvalError::Interrupted),
        Err(PolyError::Ball(_)) => {
            warn(
                ev,
                "the coefficients must be exact, not decimals (try exact)",
                e,
            );
            Ok(None)
        }
        Err(PolyError::TooLarge) => {
            let terms = ev.context().limits().terms;
            warn(ev, format!("a polynomial could pass {terms} terms"), e);
            Ok(None)
        }
        Err(PolyError::Bits) => {
            let bits = ev.context().limits().number_bits;
            warn(ev, format!("a coefficient could pass {bits} bits"), e);
            Ok(None)
        }
        Err(PolyError::DivisionByZero) => {
            warn(ev, "a denominator is zero", e);
            Ok(None)
        }
        Err(PolyError::Infinite(_)) => {
            warn(ev, "infinities and Indeterminate are not polynomials", e);
            Ok(None)
        }
    }
}

/// Whether `e` is an infinity or `Indeterminate`, which `Expand` and
/// `Together` leave as they are.
fn infinite(e: &Expr) -> bool {
    matches!(
        e.as_symbol(),
        Some(Symbol::INFINITY | Symbol::COMPLEX_INFINITY | Symbol::INDETERMINATE)
    ) || e.head() == Some(&Expr::from(Symbol::DIRECTED_INFINITY))
}

/// `e` with `head` around the arguments of every part the conversion takes
/// as a generator, so that `Expand` and `Together` reach inside them
/// (P6): `Sin[(x + 1)^2]` becomes `Sin[Expand[(x + 1)^2]]`. Parts whose
/// head holds an argument in `ctx`, or is no symbol, stay as they are.
/// `None` when nothing changes.
fn inside(ctx: &Context, e: &Expr, head: Symbol) -> Option<Expr> {
    let View::Apply { head: h, args } = e.view() else {
        return None;
    };
    let symbol = h.as_symbol();
    let arithmetic = matches!(symbol, Some(Symbol::PLUS | Symbol::TIMES))
        || symbol == Some(Symbol::POWER) && args.len() == 2 && args[1].as_i64().is_some();
    let holds = symbol
        .is_none_or(|s| ctx.attributes(s).intersection(Attributes::HOLD_ALL) != Attributes::NONE);
    if holds && !arithmetic {
        return None;
    }
    let mut changed = false;
    let new: Vec<Expr> = args
        .iter()
        .enumerate()
        .map(|(k, a)| {
            if arithmetic {
                // The exponent of an integer power is a number.
                if symbol == Some(Symbol::POWER) && k == 1 {
                    return a.clone();
                }
                match inside(ctx, a, head) {
                    Some(b) => {
                        changed = true;
                        b
                    }
                    None => a.clone(),
                }
            } else if a.is_number() || a.as_symbol().is_some() {
                a.clone()
            } else {
                changed = true;
                Expr::call(head, [a.clone()])
            }
        })
        .collect();
    changed.then(|| Expr::call(h.clone(), new))
}

/// `arg` with [`inside`] applied and evaluated.
fn reach_inside(ev: &mut Eval<'_>, arg: &Expr, head: Symbol) -> Result<Expr, EvalError> {
    match inside(ev.context(), arg, head) {
        Some(e) => ev.eval(&e),
        None => Ok(arg.clone()),
    }
}

/// `Expand[e]`, giac's `expand` (P6): products and positive powers
/// multiplied out, each term over its own expanded denominator, inside
/// function arguments too ([`poly::expand`]).
pub(super) fn expand(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if infinite(arg) {
        return Ok(Some(arg.clone()));
    }
    let arg = reach_inside(ev, arg, Symbol::EXPAND)?;
    let ctx = ev.context();
    let expanded = poly::expand(&arg, bounds(ctx), ctx.interrupt_flag());
    let Some(result) = refused(ev, e, expanded)? else {
        return Ok(None);
    };
    ev.step("expand", e, &result);
    Ok(Some(result))
}

/// A numerator and a denominator.
type Ratio<E> = (Sparse<E>, Sparse<E>);

/// `num / den` with their gcd cancelled and the denominator's leading
/// coefficient in normal form.
fn cancelled<R: Numbers>(
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Ratio<R::Elem>, EvalError> {
    let (_, n, d) = ring
        .gcd_cofactors_with(num, den, cancel)
        .map_err(interrupted)?;
    Ok(normal_sign(ring, n, d))
}

/// `n / d` scaled by the unit that puts `d`'s leading coefficient in normal
/// form.
fn normal_sign<R: Numbers>(
    ring: &SparseRing<R>,
    n: Sparse<R::Elem>,
    d: Sparse<R::Elem>,
) -> Ratio<R::Elem> {
    let base = ring.base();
    let Some(lead) = d.coeffs().first() else {
        return (n, d);
    };
    let unit = base.normal_form(lead).1;
    if base.is_one(&unit) {
        (n, d)
    } else {
        (ring.scale(&n, &unit), ring.scale(&d, &unit))
    }
}

/// `Together[e]`, Xcas's `normal`: an expanded numerator over an expanded
/// denominator, their gcd cancelled, the denominator's leading coefficient
/// in normal form: positive over ℤ, in the first quadrant over ℤ\[i\]
/// (P7).
pub(super) fn together(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if infinite(arg) {
        return Ok(Some(arg.clone()));
    }
    let arg = reach_inside(ev, arg, Symbol::TOGETHER)?;
    let Some(p) = convert(ev, e, std::slice::from_ref(&arg), false)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::Pair).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        let (n, d) = cancelled(&ring, num, den, cancel)?;
        p.ratio_expr(&n, &d)
    });
    ev.step("normal", e, &result);
    Ok(Some(result))
}

/// The exact numbers of `args`, as Gaussian rationals, when all are.
fn numbers(args: &[Expr]) -> Option<Vec<Complex>> {
    args.iter().map(Expr::to_complex).collect()
}

/// The gcd of rationals: of their numerators, over the lcm of their
/// denominators, as giac (`gcd(1/2, 1/3)` is `1/6`).
fn rational_gcd(a: &Rational, b: &Rational) -> Rational {
    Rational::new(a.numer().gcd(b.numer()), a.denom().lcm(b.denom()))
}

/// The lcm of rationals: of their numerators, over the gcd of their
/// denominators (`lcm(1/2, 1/3)` is 1).
fn rational_lcm(a: &Rational, b: &Rational) -> Rational {
    Rational::new(a.numer().lcm(b.numer()), a.denom().gcd(b.denom()))
}

/// The polynomials of a conversion, each fraction's denominator constant;
/// `None` after a warning otherwise.
fn polynomials(ev: &mut Eval<'_>, e: &Expr, p: &Polynomials) -> Option<()> {
    let constant = p.fractions().iter().all(|f| match f {
        Fraction::Integers { den, .. } => is_constant(den),
        Fraction::GaussianIntegers { den, .. } => is_constant(den),
    });
    if constant {
        Some(())
    } else {
        warn(ev, "the arguments must be polynomials", e);
        None
    }
}

fn is_constant<E>(p: &Sparse<E>) -> bool {
    p.degrees().iter().all(|&d| d == 0)
}

/// The variables whose degrees bound the dense forms of an exact
/// algorithm (CAS-156).
#[derive(Clone, Copy)]
pub(super) enum Dense {
    /// Every variable once two of the polynomials are not constants: a gcd
    /// of two goes dense in any variable either has, even one the other
    /// lacks (`gcd(x^(10^8) - 1, y + 1)`), where a constant's gcd is its
    /// gcd with the content.
    Pair,
    /// Every variable: a resultant grows with a degree alone, as
    /// `resultant(c, b) = c^(deg b)`.
    All,
    /// Every variable but `x_k`: a content in `x_k` takes the gcds of the
    /// coefficients.
    Except(usize),
}

/// Whether the dense forms of an exact algorithm on `p` stay within
/// `Limits::terms`: no variable `rule` names may have a degree of the
/// limit or more, as gcds, resultants and divisions make arrays of a
/// degree's coefficients (CAS-156): `gcd(x^(10^8) - 1, x^(10^6) - 1)` went
/// past any memory. `None` after the warning of the limit.
pub(super) fn within_dense(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Polynomials,
    rule: Dense,
) -> Option<()> {
    let terms = u64::try_from(ev.context().limits().terms).unwrap_or(u64::MAX);
    let degrees: Vec<Vec<u64>> = p
        .fractions()
        .iter()
        .flat_map(|f| match f {
            Fraction::Integers { num, den } => [num.degrees(), den.degrees()],
            Fraction::GaussianIntegers { num, den } => [num.degrees(), den.degrees()],
        })
        .collect();
    let polynomials = degrees.iter().filter(|d| d.iter().any(|&e| e > 0)).count();
    let past = (0..p.generators().len()).any(|v| {
        let largest = degrees.iter().map(|d| d[v]).max().unwrap_or(0);
        largest >= terms
            && match rule {
                Dense::Pair => polynomials >= 2,
                Dense::All => true,
                Dense::Except(k) => v != k,
            }
    });
    if past {
        let terms = ev.context().limits().terms;
        warn(ev, format!("a polynomial could pass {terms} terms"), e);
        None
    } else {
        Some(())
    }
}

/// `GCD[a, b, …]` (§15 detail 5): of numbers, as numbers (gcd of the
/// numerators over the lcm of the denominators); of polynomials, the gcd of
/// their numerators over ℤ or ℤ\[i\], so the gcd of the contents times
/// that of the primitive parts, its leading coefficient positive.
pub(super) fn gcd(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if let Some(values) = numbers(args) {
        return Ok(
            number_fold(&values, rational_gcd, gaussian_gcd).inspect(|g| ev.step("gcd", e, g))
        );
    }
    let Some(p) = convert(ev, e, args, false)? else {
        return Ok(None);
    };
    if polynomials(ev, e, &p).is_none() || within_dense(ev, e, &p, Dense::Pair).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let n = p.generators().len();
    let result = match p.fractions()[0] {
        Fraction::Integers { .. } => fold_gcd::<Integers>(&p, n, cancel)?,
        Fraction::GaussianIntegers { .. } => fold_gcd::<GaussianIntegers>(&p, n, cancel)?,
    };
    ev.step("gcd", e, &result);
    Ok(Some(result))
}

/// The numerators of `p`'s fractions over `R`.
trait Numerator<R: Numbers> {
    fn numerators(&self) -> Vec<Ratio<R::Elem>>;
}

impl Numerator<Integers> for Polynomials {
    fn numerators(&self) -> Vec<(Sparse<Integer>, Sparse<Integer>)> {
        self.fractions()
            .iter()
            .map(|f| match f {
                Fraction::Integers { num, den } => (num.clone(), den.clone()),
                Fraction::GaussianIntegers { .. } => unreachable!("one ring"),
            })
            .collect()
    }
}

impl Numerator<GaussianIntegers> for Polynomials {
    fn numerators(&self) -> Vec<(Sparse<Complex>, Sparse<Complex>)> {
        self.fractions()
            .iter()
            .map(|f| match f {
                Fraction::GaussianIntegers { num, den } => (num.clone(), den.clone()),
                Fraction::Integers { .. } => unreachable!("one ring"),
            })
            .collect()
    }
}

fn fold_gcd<R: Numbers>(p: &Polynomials, n: usize, cancel: &dyn Cancel) -> Result<Expr, EvalError>
where
    Polynomials: Numerator<R>,
{
    let ring = SparseRing::new(R::default(), n);
    let mut g = ring.zero();
    for (num, _) in Numerator::<R>::numerators(p) {
        g = ring.gcd_with(&g, &num, cancel).map_err(interrupted)?;
    }
    Ok(p.polynomial_expr(&g))
}

/// The gcd of Gaussian rationals that are Gaussian integers, in normal
/// form; `None` otherwise.
fn gaussian_gcd(a: &Complex, b: &Complex) -> Option<Complex> {
    (a.is_gaussian_integer() && b.is_gaussian_integer()).then(|| GaussianIntegers.gcd(a, b))
}

/// The gaussian lcm, `a b / gcd(a, b)` in normal form.
fn gaussian_lcm(a: &Complex, b: &Complex) -> Option<Complex> {
    let g = gaussian_gcd(a, b)?;
    if g.is_zero() {
        return Some(Complex::from(0));
    }
    Some(GaussianIntegers.normal_form(&(&(a * b) / &g)).0)
}

/// The gcd or lcm of numbers: rationals by `real`, Gaussian integers by
/// `complex`; `None` for other mixtures.
fn number_fold(
    values: &[Complex],
    real: fn(&Rational, &Rational) -> Rational,
    complex: fn(&Complex, &Complex) -> Option<Complex>,
) -> Option<Expr> {
    let (first, rest) = values.split_first()?;
    if values.iter().all(Complex::is_real) {
        let mut acc = first.re().clone();
        for v in rest {
            acc = real(&acc, v.re());
        }
        return Some(Expr::from(acc.abs()));
    }
    let mut acc = first.clone();
    for v in rest {
        acc = complex(&acc, v)?;
    }
    Some(Expr::from(GaussianIntegers.normal_form(&acc).0))
}

/// `LCM[a, b, …]`: of numbers, as numbers; of polynomials, `a b / gcd(a,
/// b)` folded, the gcd of the numerators as in [`gcd`], expanded and
/// normalised (§15 detail 5).
pub(super) fn lcm(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let args = e.args();
    if let Some(values) = numbers(args) {
        return Ok(
            number_fold(&values, rational_lcm, gaussian_lcm).inspect(|l| ev.step("lcm", e, l))
        );
    }
    let Some(p) = convert(ev, e, args, false)? else {
        return Ok(None);
    };
    if polynomials(ev, e, &p).is_none() || within_dense(ev, e, &p, Dense::Pair).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let n = p.generators().len();
    let result = match p.fractions()[0] {
        Fraction::Integers { .. } => fold_lcm::<Integers>(&p, n, cancel)?,
        Fraction::GaussianIntegers { .. } => fold_lcm::<GaussianIntegers>(&p, n, cancel)?,
    };
    ev.step("lcm", e, &result);
    Ok(Some(result))
}

fn fold_lcm<R: Numbers>(p: &Polynomials, n: usize, cancel: &dyn Cancel) -> Result<Expr, EvalError>
where
    Polynomials: Numerator<R>,
{
    let ring = SparseRing::new(R::default(), n);
    let mut parts = Numerator::<R>::numerators(p).into_iter();
    let (mut num, mut den) = parts.next().expect("an argument");
    for (b, db) in parts {
        if num.is_zero() || b.is_zero() {
            return Ok(Expr::from(0));
        }
        // a b / gcd(a, b), over the product of the denominators.
        let (_, _, rest) = ring
            .gcd_cofactors_with(&num, &b, cancel)
            .map_err(interrupted)?;
        num = ring.mul_with(&num, &rest, cancel).map_err(interrupted)?;
        den = ring.mul(&den, &db);
    }
    let (num, den) = normal_sign(&ring, num, den);
    let (num, den) = lead_positive(&ring, num, den);
    Ok(p.ratio_expr(&num, &den))
}

/// `num / den` with `num`'s leading coefficient in normal form, for a
/// constant `den`: the normalised lcm.
fn lead_positive<R: Numbers>(
    ring: &SparseRing<R>,
    num: Sparse<R::Elem>,
    den: Sparse<R::Elem>,
) -> Ratio<R::Elem> {
    let base = ring.base();
    let Some(lead) = num.coeffs().first() else {
        return (num, den);
    };
    let unit = base.normal_form(lead).1;
    (ring.scale(&num, &unit), den)
}

/// The variable of a call.
#[derive(Clone, Copy)]
enum Variable {
    /// The generator `x_k`.
    Index(usize),
    /// None that occurs: the expressions are constants in it.
    Absent,
}

/// Whether `g` contains a symbol other than a builtin constant, so that it
/// can be the variable (§15 detail 6): `x`, `Sin[x]`, but not `Pi` or
/// `2^(1/2)`.
fn has_symbol(g: &Expr) -> bool {
    g.symbols().into_iter().any(|s| !s.is_builtin())
}

/// The variable of the call `e`: `given`, which no other generator may
/// contain, or else the one generator with a symbol (§15 detail 6). `None`
/// after a warning: the input is not a polynomial in the given variable,
/// or the variable is one of several.
fn variable(
    ev: &mut Eval<'_>,
    e: &Expr,
    p: &Polynomials,
    given: Option<&Expr>,
) -> Option<Variable> {
    let generators = p.generators();
    if let Some(x) = given {
        if x.is_number() {
            warn(ev, "the variable cannot be a number", e);
            return None;
        }
        let index = p.variable(x);
        let others = generators
            .iter()
            .enumerate()
            .filter(|&(k, g)| Some(k) != index && g.contains(x));
        if others.count() > 0 {
            warn(ev, format!("not a polynomial in {x}"), e);
            return None;
        }
        return Some(index.map_or(Variable::Absent, Variable::Index));
    }
    let candidates: Vec<usize> = (0..generators.len())
        .filter(|&k| has_symbol(&generators[k]))
        .collect();
    match candidates[..] {
        [] => Some(Variable::Absent),
        [k] => Some(Variable::Index(k)),
        _ => {
            let names: Vec<String> = candidates
                .iter()
                .rev()
                .map(|&k| generators[k].to_string())
                .collect();
            warn(
                ev,
                format!("the variable must be given: one of {}", names.join(", ")),
                e,
            );
            None
        }
    }
}

/// Whether the denominators of `p`'s fractions are free of `x_var`, as a
/// polynomial in it needs; `None` after a warning otherwise.
fn free_denominators(ev: &mut Eval<'_>, e: &Expr, p: &Polynomials, var: usize) -> Option<()> {
    let free = p.fractions().iter().all(|f| match f {
        Fraction::Integers { den, .. } => den.degrees()[var] == 0,
        Fraction::GaussianIntegers { den, .. } => den.degrees()[var] == 0,
    });
    if free {
        Some(())
    } else {
        warn(
            ev,
            format!("not a polynomial in {}", p.generators()[var]),
            e,
        );
        None
    }
}

/// The polynomial and its variable of a one-polynomial call
/// `head[p, x]` or `head[p]`, converted with balls when `balls`.
fn one_polynomial(
    ev: &mut Eval<'_>,
    e: &Expr,
    balls: bool,
) -> Result<Option<(Polynomials, Variable)>, EvalError> {
    let (p, given) = match e.args() {
        [p] => (p, None),
        [p, x] => (p, Some(x)),
        _ => return Ok(None),
    };
    let Some(polys) = convert(ev, e, std::slice::from_ref(p), balls)? else {
        return Ok(None);
    };
    let Some(var) = variable(ev, e, &polys, given) else {
        return Ok(None);
    };
    if let Variable::Index(k) = var
        && free_denominators(ev, e, &polys, k).is_none()
    {
        return Ok(None);
    }
    Ok(Some((polys, var)))
}

/// `-Infinity`, the degree of 0 (§15 detail 7).
fn minus_infinity() -> Expr {
    Expr::call(
        Symbol::TIMES,
        [Expr::from(-1), Expr::from(Symbol::INFINITY)],
    )
}

/// `Exponent[p, x]`, Xcas's `degree`: the degree of `p` in `x`, inferred
/// when not given; 0 for a nonzero constant, `-Infinity` for 0 (§15 detail
/// 7). Balls are coefficients.
pub(super) fn exponent(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, var)) = one_polynomial(ev, e, true)? else {
        return Ok(None);
    };
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        let _ = (&ring, den);
        if num.is_zero() {
            minus_infinity()
        } else {
            match var {
                Variable::Index(k) => Expr::from(Integer::from(num.degrees()[k])),
                Variable::Absent => Expr::from(0),
            }
        }
    });
    ev.step("degree", e, &result);
    Ok(Some(result))
}

/// `Coefficient[p, x, n]`, Xcas's `coeff`: the coefficient of `x^n`;
/// `Coefficient[p, n]` with `x` inferred; `Coefficient[p, x]` the list of
/// the coefficients from the highest power, as giac. Balls are
/// coefficients.
pub(super) fn coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (p, given, power) = match e.args() {
        [p, n] if n.is_exact_number() => (p, None, Some(n)),
        [p, x] => (p, Some(x), None),
        [p, x, n] => (p, Some(x), Some(n)),
        _ => return Ok(None),
    };
    let power = if let Some(n) = power {
        let Some(n) = n.as_i64() else {
            warn(ev, "the power must be an integer", e);
            return Ok(None);
        };
        Some(n)
    } else {
        None
    };
    let Some(polys) = convert(ev, e, std::slice::from_ref(p), true)? else {
        return Ok(None);
    };
    let Some(var) = variable(ev, e, &polys, given) else {
        return Ok(None);
    };
    if let Variable::Index(k) = var
        && free_denominators(ev, e, &polys, k).is_none()
    {
        return Ok(None);
    }
    // The list of coefficients has the degree's length (CAS-156).
    if let (Variable::Index(k), None) = (var, power) {
        let degree = match &polys.fractions()[0] {
            Fraction::Integers { num, .. } => num.degrees()[k],
            Fraction::GaussianIntegers { num, .. } => num.degrees()[k],
        };
        let terms = ev.context().limits().terms;
        if degree >= u64::try_from(terms).unwrap_or(u64::MAX) {
            warn(ev, format!("a polynomial could pass {terms} terms"), e);
            return Ok(None);
        }
    }
    let result = with_fraction!(
        &polys.fractions()[0],
        polys.generators().len(),
        |ring, num, den| {
            match (var, power) {
                (Variable::Index(k), Some(n)) => match u64::try_from(n) {
                    Ok(n) => polys.ratio_expr(&ring.coeff(num, k, n), den),
                    Err(_) => Expr::from(0),
                },
                (Variable::Absent, Some(n)) => {
                    if n == 0 {
                        polys.ratio_expr(num, den)
                    } else {
                        Expr::from(0)
                    }
                }
                (Variable::Index(k), None) => {
                    let degree = ring.degree(num, k).unwrap_or(0);
                    Expr::call(
                        Symbol::LIST,
                        (0..=degree)
                            .rev()
                            .map(|d| polys.ratio_expr(&ring.coeff(num, k, d), den)),
                    )
                }
                (Variable::Absent, None) => Expr::call(Symbol::LIST, [polys.ratio_expr(num, den)]),
            }
        }
    );
    ev.step("coeff", e, &result);
    Ok(Some(result))
}

/// `LeadingCoefficient[p, x]`, Xcas's `lcoeff`: the coefficient of the
/// highest power of `x`, inferred when not given; the constant itself for
/// a constant, 0 for 0. Balls are coefficients.
pub(super) fn leading_coefficient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, var)) = one_polynomial(ev, e, true)? else {
        return Ok(None);
    };
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        match var {
            Variable::Index(k) => p.ratio_expr(&ring.lead_coeff(num, k), den),
            Variable::Absent => p.ratio_expr(num, den),
        }
    });
    ev.step("lcoeff", e, &result);
    Ok(Some(result))
}

/// `Numerator[e]` and `Denominator[e]`, Xcas's `numer` and `denom`: the
/// numerator and denominator of `e` as `Together` writes it, over ℤ or
/// ℤ\[i\], the denominator's leading coefficient positive: `numer(1/2 + x)`
/// is `2 x + 1`.
fn numerator_or_denominator(ev: &mut Eval<'_>, e: &Expr, numerator: bool) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    let Some(p) = convert(ev, e, std::slice::from_ref(arg), false)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::Pair).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        let (n, d) = cancelled(&ring, num, den, cancel)?;
        p.polynomial_expr(if numerator { &n } else { &d })
    });
    ev.step(if numerator { "numer" } else { "denom" }, e, &result);
    Ok(Some(result))
}

pub(super) fn numerator(ev: &mut Eval<'_>, e: &Expr) -> Native {
    numerator_or_denominator(ev, e, true)
}

pub(super) fn denominator(ev: &mut Eval<'_>, e: &Expr) -> Native {
    numerator_or_denominator(ev, e, false)
}

/// `Content[p, x]` and `PrimitivePart[p, x]`, Xcas's `content` and
/// `primpart`: `p` as its content in `x`, the gcd of its coefficients as a
/// polynomial in `x`, normalised positive, times its primitive part, which
/// keeps the sign (`primpart(-6x^2 - 4x)` is `-3x^2 - 2x`, as giac). Over ℚ
/// the content is rational: `content(x^2/2 + x/3)` is `1/6`.
fn content_or_part(ev: &mut Eval<'_>, e: &Expr, content: bool) -> Native {
    let Some((p, var)) = one_polynomial(ev, e, false)? else {
        return Ok(None);
    };
    if let Variable::Index(k) = var
        && within_dense(ev, e, &p, Dense::Except(k)).is_none()
    {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        let c = match var {
            Variable::Index(k) => ring.content_in(num, k),
            Variable::Absent => {
                let (normal, _) = normal_lead(&ring, num);
                normal
            }
        };
        if content {
            p.ratio_expr(&c, den)
        } else if num.is_zero() {
            Expr::from(0)
        } else {
            let part = ring
                .div_exact_with(num, &c, cancel)
                .map_err(interrupted)?
                .expect("the content divides");
            p.polynomial_expr(&part)
        }
    });
    ev.step(if content { "content" } else { "primpart" }, e, &result);
    Ok(Some(result))
}

/// `p` times the unit that puts its leading coefficient in normal form, and
/// the unit.
fn normal_lead<R: Numbers>(
    ring: &SparseRing<R>,
    p: &Sparse<R::Elem>,
) -> (Sparse<R::Elem>, R::Elem) {
    let base = ring.base();
    let Some(lead) = p.coeffs().first() else {
        return (p.clone(), base.one());
    };
    let unit = base.normal_form(lead).1;
    (ring.scale(p, &unit), unit)
}

pub(super) fn content(ev: &mut Eval<'_>, e: &Expr) -> Native {
    content_or_part(ev, e, true)
}

pub(super) fn primitive_part(ev: &mut Eval<'_>, e: &Expr) -> Native {
    content_or_part(ev, e, false)
}

/// The two polynomials and their variable of a call `head[a, b, x]` or
/// `head[a, b]`, exact, their denominators free of the variable.
fn two_polynomials(
    ev: &mut Eval<'_>,
    e: &Expr,
) -> Result<Option<(Polynomials, Variable)>, EvalError> {
    let (pair, given) = match e.args() {
        [a, b] => ([a.clone(), b.clone()], None),
        [a, b, x] => ([a.clone(), b.clone()], Some(x)),
        _ => return Ok(None),
    };
    let Some(p) = convert(ev, e, &pair, false)? else {
        return Ok(None);
    };
    let Some(var) = variable(ev, e, &p, given) else {
        return Ok(None);
    };
    if let Variable::Index(k) = var
        && free_denominators(ev, e, &p, k).is_none()
    {
        return Ok(None);
    }
    Ok(Some((p, var)))
}

/// What a division gives.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Division {
    Quotient,
    Remainder,
    Both,
}

/// `PolynomialQuotient`, `PolynomialRemainder` and
/// `PolynomialQuotientRemainder[a, b, x]`, Xcas's `quo`, `rem` and
/// `quorem`: `a = q b + r` with `deg r < deg b` in `x`, over the fractions
/// of the other generators, as giac (`quo(x^2, y x + 1, x)` is `x/y -
/// 1/y^2`): from the pseudo-division `lc(b)^(δ+1) a = Q b + R`, `q = Q /
/// lc(b)^(δ+1)`, cancelled. By 0, as `Quotient` of numbers: the quotient is
/// `ComplexInfinity` and the remainder `Indeterminate`, with a warning.
#[allow(
    clippy::many_single_char_names,
    reason = "a = q b + r, as written in mathematics"
)]
fn divide(ev: &mut Eval<'_>, e: &Expr, what: Division) -> Native {
    let Some((p, var)) = two_polynomials(ev, e)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::Pair).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let n = p.generators().len();
    let parts = with_fractions!(
        &p.fractions()[0],
        &p.fractions()[1],
        n,
        |ring, na, da, nb, db| {
            if nb.is_zero() {
                None
            } else {
                Some(divided(&ring, &p, (na, da), (nb, db), var, cancel)?)
            }
        }
    );
    let Some((q, r)) = parts else {
        return Ok(Some(match what {
            Division::Quotient => {
                ev.warn(Warning::Infinite(e.clone()));
                Expr::from(Symbol::COMPLEX_INFINITY)
            }
            Division::Remainder | Division::Both => {
                ev.warn(Warning::Indeterminate(e.clone()));
                Expr::from(Symbol::INDETERMINATE)
            }
        }));
    };
    let result = match what {
        Division::Quotient => q,
        Division::Remainder => r,
        Division::Both => Expr::call(Symbol::LIST, [q, r]),
    };
    ev.step("division", e, &result);
    Ok(Some(result))
}

/// `(q, r)` with `a = q b + r`, as expressions, `b` not 0.
#[allow(
    clippy::many_single_char_names,
    reason = "a = q b + r, as written in mathematics"
)]
fn divided<R: Numbers>(
    ring: &SparseRing<R>,
    p: &Polynomials,
    (na, da): (&Sparse<R::Elem>, &Sparse<R::Elem>),
    (nb, db): (&Sparse<R::Elem>, &Sparse<R::Elem>),
    var: Variable,
    cancel: &dyn Cancel,
) -> Result<(Expr, Expr), EvalError> {
    let Variable::Index(k) = var else {
        // Constants in every variable: a / b, remainder 0.
        let num = ring.mul(na, db);
        let den = ring.mul(da, nb);
        let (n, d) = cancelled(ring, &num, &den, cancel)?;
        return Ok((p.ratio_expr(&n, &d), Expr::from(0)));
    };
    if nb.degrees()[k] == 0 {
        // A divisor constant in x_k: a / b, remainder 0, with no dense form
        // of `a` (CAS-156).
        let (n, d) = cancelled(ring, &ring.mul(na, db), &ring.mul(da, nb), cancel)?;
        return Ok((p.expanded_expr(&n, &d), Expr::from(0)));
    }
    let dense = DenseRing::new(ring.clone());
    let recursive = |x: &Sparse<R::Elem>| -> cas_poly::Dense<Sparse<R::Elem>> {
        let parts = ring.to_recursive(x, k);
        let degree = parts.first().map_or(0, |&(d, _)| d);
        let mut coeffs =
            vec![ring.zero(); usize::try_from(degree).expect("a degree in memory") + 1];
        for (d, c) in parts {
            coeffs[usize::try_from(d).expect("at most the degree")] = c;
        }
        dense.poly(coeffs)
    };
    let (a, b) = (recursive(na), recursive(nb));
    let (q, r) = dense.pseudo_div_rem(&a, &b);
    let back = |x: &cas_poly::Dense<Sparse<R::Elem>>| {
        let parts: Vec<(u64, Sparse<R::Elem>)> = x
            .coeffs()
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(d, c)| (d as u64, c.clone()))
            .collect();
        ring.from_recursive(k, &parts)
    };
    // lc(b)^(δ+1) a = Q b + R, δ = deg a - deg b, when deg a ≥ deg b.
    let lead = b.lead().expect("a divisor").clone();
    let power = match (a.degree(), b.degree()) {
        (Some(m), Some(n)) if m >= n => (m - n + 1) as u64,
        _ => 0,
    };
    let scale = ring.pow(&lead, power);
    // a = Na / da, b = Nb / db: q = Q db / (da lc^(δ+1)), r = R / (da lc^(δ+1)).
    let den = ring.mul(da, &scale);
    let (qn, qd) = cancelled(ring, &ring.mul(&back(&q), db), &den, cancel)?;
    let (rn, rd) = cancelled(ring, &back(&r), &den, cancel)?;
    Ok((p.expanded_expr(&qn, &qd), p.expanded_expr(&rn, &rd)))
}

pub(super) fn polynomial_quotient(ev: &mut Eval<'_>, e: &Expr) -> Native {
    divide(ev, e, Division::Quotient)
}

pub(super) fn polynomial_remainder(ev: &mut Eval<'_>, e: &Expr) -> Native {
    divide(ev, e, Division::Remainder)
}

pub(super) fn polynomial_quotient_remainder(ev: &mut Eval<'_>, e: &Expr) -> Native {
    divide(ev, e, Division::Both)
}

/// `Resultant[a, b, x]`: the determinant of the Sylvester matrix in `x`,
/// inferred when not given; over ℚ from the numerators, `res(n/d, n'/d') =
/// res(n, n') / (d^(deg n') d'^(deg n))`. A constant `c` gives `c^(deg b)`.
#[allow(
    clippy::many_single_char_names,
    reason = "res(n/d, n'/d'), as written in mathematics"
)]
pub(super) fn resultant(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, var)) = two_polynomials(ev, e)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::All).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let n = p.generators().len();
    let result = with_fractions!(
        &p.fractions()[0],
        &p.fractions()[1],
        n,
        |ring, na, da, nb, db| {
            match var {
                _ if na.is_zero() || nb.is_zero() => Expr::from(0),
                // Two constants: the empty determinant.
                Variable::Absent => Expr::from(1),
                Variable::Index(k) => {
                    let res = ring
                        .resultant_with(na, nb, k, cancel)
                        .map_err(interrupted)?;
                    let (m, d) = (
                        ring.degree(na, k).unwrap_or(0),
                        ring.degree(nb, k).unwrap_or(0),
                    );
                    let den = ring.mul(&ring.pow(da, d), &ring.pow(db, m));
                    let (rn, rd) = cancelled(&ring, &res, &den, cancel)?;
                    p.ratio_expr(&rn, &rd)
                }
            }
        }
    );
    ev.step("resultant", e, &result);
    Ok(Some(result))
}

/// `Discriminant[p, x]`: `(-1)^(n(n-1)/2) res(p, p') / lc(p)` in `x`,
/// inferred when not given, for a degree `n ≥ 1`, so `b^2 - 4 a c` for a
/// quadratic `a x^2 + b x + c`. Over ℚ, `disc(n / d) = disc(n) /
/// d^(2 deg n - 2)`. A constant has none: a warning, and the call stays.
pub(super) fn discriminant(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((p, var)) = one_polynomial(ev, e, false)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::All).is_none() {
        return Ok(None);
    }
    let cancel = ev.context().interrupt_flag();
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        match var {
            Variable::Absent => None,
            Variable::Index(k) => match ring
                .discriminant_with(num, k, cancel)
                .map_err(interrupted)?
            {
                None => None,
                Some(d) => {
                    let degree = ring.degree(num, k).expect("degree at least 1");
                    let den = ring.pow(den, 2 * degree - 2);
                    let (rn, rd) = cancelled(&ring, &d, &den, cancel)?;
                    Some(p.ratio_expr(&rn, &rd))
                }
            },
        }
    });
    let Some(result) = result else {
        warn(ev, "a discriminant needs a degree of at least 1", e);
        return Ok(None);
    };
    ev.step("discriminant", e, &result);
    Ok(Some(result))
}

/// The rational function of a one-argument factoring call; `None` after a
/// warning: refused by the conversion, or a degree that would go past
/// `Limits::terms` in dense form.
fn to_factor(ev: &mut Eval<'_>, e: &Expr, arg: &Expr) -> Result<Option<Polynomials>, EvalError> {
    let Some(p) = convert(ev, e, std::slice::from_ref(arg), false)? else {
        return Ok(None);
    };
    if within_dense(ev, e, &p, Dense::All).is_none() {
        return Ok(None);
    }
    Ok(Some(p))
}

/// The coefficient rings of the conversion that cas-poly factors over: ℤ
/// over ℚ, and ℤ\[i\] over ℚ(i) (CAS-186).
pub(super) trait Factoring: Numbers + CharacteristicZero {
    fn factor_within(
        ring: &SparseRing<Self>,
        p: &Sparse<Self::Elem>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Self::Elem, Sparse<Self::Elem>>, FactorError>;

    /// The real and imaginary parts, the order of `factors`' lists.
    fn parts(c: &Self::Elem) -> (Rational, Rational);

    /// Whether the coefficients are real: then `factor` splits only
    /// quadratics of a positive discriminant, as giac does without `i`.
    const REAL: bool;

    fn complex(c: &Self::Elem) -> Complex {
        let (re, im) = Self::parts(c);
        Complex::new(re, im)
    }
}

impl Factoring for Integers {
    fn factor_within(
        ring: &SparseRing<Integers>,
        p: &Sparse<Integer>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Integer, Sparse<Integer>>, FactorError> {
        ring.factor_within(p, bound, cancel)
    }

    fn parts(c: &Integer) -> (Rational, Rational) {
        (Rational::from(c.clone()), Rational::ZERO)
    }

    const REAL: bool = true;
}

impl Factoring for GaussianIntegers {
    fn factor_within(
        ring: &SparseRing<GaussianIntegers>,
        p: &Sparse<Complex>,
        bound: u64,
        cancel: &dyn Cancel,
    ) -> Result<Factors<Complex, Sparse<Complex>>, FactorError> {
        ring.factor_within(p, bound, cancel)
    }

    fn parts(c: &Complex) -> (Rational, Rational) {
        (c.re().clone(), c.im().clone())
    }

    const REAL: bool = false;
}

/// `p` times the unit that puts its leading coefficient in normal form in
/// the reverse of the conversion's order, the first generator in canonical
/// order most significant, as giac's main variable is; and the inverse of
/// the unit. So `factor(x^4-y^4)` has the factor `x - y`, as giac's does,
/// where the main variable of CAS-DOC-8 §15 detail 5 would make it
/// `y - x`.
fn readable<R: Numbers>(ring: &SparseRing<R>, p: &Sparse<R::Elem>) -> (Sparse<R::Elem>, R::Elem) {
    let base = ring.base();
    let lead = p
        .terms()
        .max_by(|(e, _), (f, _)| e.iter().rev().cmp(f.iter().rev()))
        .map(|(_, c)| c.clone());
    let Some(lead) = lead else {
        return (p.clone(), base.one());
    };
    let unit = base.normal_form(&lead).1;
    if base.is_one(&unit) {
        return (p.clone(), base.one());
    }
    (ring.scale(p, &unit), base.div_rem(&base.one(), &unit).0)
}

/// A rational function as a constant over a constant times powers of
/// polynomials, the exponents of the denominator's negative.
pub(super) struct Powers<E> {
    numerator: E,
    denominator: E,
    pub(super) factors: Vec<(Sparse<E>, i64)>,
}

impl<E: Coefficient + PartialEq> Powers<E> {
    /// The parts of a numerator's and a denominator's factorisations, each
    /// made [`readable`], the constants taking the units, and equal parts
    /// merged, those that cancel left out.
    fn new<R: Numbers<Elem = E>>(
        ring: &SparseRing<R>,
        top: Factors<E, Sparse<E>>,
        bottom: Factors<E, Sparse<E>>,
    ) -> Powers<E> {
        let base = ring.base();
        let (mut numerator, top) = top.into_parts();
        let (mut denominator, bottom) = bottom.into_parts();
        let mut factors: Vec<(Sparse<E>, i64)> = Vec::new();
        for (side, sign) in [(top, 1), (bottom, -1)] {
            for (f, e) in side {
                let (f, inverse) = readable(ring, &f);
                let power = base.pow(&inverse, e);
                if sign > 0 {
                    numerator = base.mul(&numerator, &power);
                } else {
                    denominator = base.mul(&denominator, &power);
                }
                let k = sign * i64::try_from(e).unwrap_or(i64::MAX);
                match factors.iter_mut().find(|(g, _)| *g == f) {
                    Some((_, j)) => *j += k,
                    None => factors.push((f, k)),
                }
            }
        }
        factors.retain(|&(_, k)| k != 0);
        Powers {
            numerator,
            denominator,
            factors,
        }
    }

    /// `c ∏ fᵢ^eᵢ`, built.
    fn expr(&self, p: &Polynomials) -> Expr {
        let mut parts = Vec::with_capacity(self.factors.len() + 1);
        let constant = self.numerator.over(&self.denominator);
        if constant != Expr::from(1) {
            parts.push(constant);
        }
        for (f, k) in &self.factors {
            let base = p.polynomial_expr(f);
            parts.push(if *k == 1 {
                base
            } else {
                base.pow(Expr::from(*k))
            });
        }
        match parts.len() {
            0 => Expr::from(1),
            1 => parts.pop().expect("one factor"),
            _ => Expr::call(Symbol::TIMES, parts),
        }
    }
}

/// `s` and the square-free `d` with `n = s² d`, from the factorisation of
/// `n`; a composite that resists it counts as square-free.
pub(super) fn root_part(n: &Natural) -> (Natural, Natural) {
    let (mut s, mut d) = (Natural::one(), Natural::one());
    for f in n.factor().factors() {
        let e = f.exponent();
        s = &s * &f.base().pow(e / 2);
        if e % 2 == 1 {
            d = &d * f.base();
        }
    }
    (s, d)
}

/// The linear factors of an irreducible quadratic `a x² + b x + c` over
/// square roots, when giac's `factor` splits it: a real discriminant
/// `D = b² − 4ac`, positive when the coefficients are real (CAS-187). With
/// `|D| = s² d`, `d` square-free, and `ι` the square root of the sign of
/// `D`, they are `x + (b ∓ ι s √d)/(2a)`, the fractions reduced by
/// `gcd(b, s, 2a)` over ℤ, so that `a` times them is the quadratic, written
/// as giac writes them: `factor(x^2-x-1)` is
/// `(x+(-sqrt(5)-1)/2)*(x+(sqrt(5)-1)/2)`.
#[allow(
    clippy::many_single_char_names,
    reason = "a x² + b x + c, D = s² d, as written in mathematics"
)]
fn split_quadratic(
    x: &Expr,
    a: &Complex,
    b: &Complex,
    c: &Complex,
    real: bool,
) -> Option<[Expr; 2]> {
    let four = Complex::from(4);
    let disc = &(b * b) - &(&(&four * a) * c);
    if !disc.im().is_zero() || disc.re().is_zero() || real && disc.re().is_negative() {
        return None;
    }
    let (s, d) = root_part(&disc.re().numer().unsigned_abs());
    if d.is_one() {
        return None;
    }
    let unit = if disc.re().is_negative() {
        Complex::I
    } else {
        Complex::ONE
    };
    let two_a = &Complex::from(2) * a;
    let (b, s, m) = if a.is_real() && b.is_real() {
        let (b, two_a) = (b.re().numer(), two_a.re().numer());
        let g = b.gcd(&Integer::from(&s)).gcd(two_a);
        (
            Complex::from(b.div_exact(&g)),
            Integer::from(&s).div_exact(&g),
            Complex::from(two_a.div_exact(&g)),
        )
    } else {
        (b.clone(), Integer::from(&s), two_a)
    };
    let root =
        Expr::from(Integer::from(d)).pow(Expr::from(Rational::new(Integer::ONE, Integer::from(2))));
    Some([1i64, -1].map(|sign| {
        let coefficient = &unit * &Complex::from(Integer::from(sign) * &s);
        let surd = if coefficient.is_one() {
            root.clone()
        } else {
            Expr::call(Symbol::TIMES, [Expr::from(&coefficient), root.clone()])
        };
        let term = if b.is_zero() {
            surd
        } else {
            Expr::call(Symbol::PLUS, [Expr::from(&b), surd])
        };
        let shift = if m.is_one() {
            term
        } else {
            Expr::call(Symbol::TIMES, [Expr::from(m.recip()), term])
        };
        Expr::call(Symbol::PLUS, [x.clone(), shift])
    }))
}

/// The coefficients `[c, b, a]` of `f` when it is a quadratic in one
/// generator, with that generator's index.
fn quadratic<E: Clone>(f: &Sparse<E>) -> Option<(usize, [Option<E>; 3])> {
    let degrees = f.degrees();
    let mut vars = degrees.iter().enumerate().filter(|&(_, &d)| d > 0);
    let (g, &degree) = vars.next()?;
    if degree != 2 || vars.next().is_some() {
        return None;
    }
    let mut coeffs = [None, None, None];
    for (e, c) in f.terms() {
        coeffs[usize::try_from(e[g]).ok()?] = Some(c.clone());
    }
    Some((g, coeffs))
}

/// [`Powers::expr`] with each quadratic factor in one generator split over
/// square roots where giac's `factor` splits it ([`split_quadratic`]), its
/// leading coefficient going to the constant.
fn split_expr<R: Factoring>(p: &Polynomials, powers: &Powers<R::Elem>) -> Expr {
    let mut constant = &R::complex(&powers.numerator) / &R::complex(&powers.denominator);
    let mut parts = Vec::with_capacity(powers.factors.len() + 1);
    let power = |base: Expr, k: i64| {
        if k == 1 {
            base
        } else {
            base.pow(Expr::from(k))
        }
    };
    for (f, k) in &powers.factors {
        let split = quadratic(f).and_then(|(g, [c, b, a])| {
            let zero = Complex::from(0);
            let value = |x: Option<R::Elem>| x.map_or_else(|| zero.clone(), |x| R::complex(&x));
            let a = value(a);
            let factors = split_quadratic(&p.generators()[g], &a, &value(b), &value(c), R::REAL)?;
            Some((a, factors))
        });
        match split {
            Some((a, [first, second])) => {
                constant = &constant * &a.pow(*k);
                parts.extend([power(first, *k), power(second, *k)]);
            }
            None => parts.push(power(p.polynomial_expr(f), *k)),
        }
    }
    if !constant.is_one() {
        parts.insert(0, Expr::from(constant));
    }
    match parts.len() {
        0 => Expr::from(1),
        1 => parts.pop().expect("one factor"),
        _ => Expr::call(Symbol::TIMES, parts),
    }
}

/// The factorisation of `num / den` over ℚ or ℚ(i), irreducible factors, the
/// searches of each within `Limits::terms` steps; `None` after a warning
/// past them.
pub(super) fn factored<R: Factoring>(
    ev: &mut Eval<'_>,
    e: &Expr,
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
) -> Result<Option<Powers<R::Elem>>, EvalError> {
    let terms = ev.context().limits().terms;
    let bound = u64::try_from(terms).unwrap_or(u64::MAX);
    let cancel = ev.context().interrupt_flag();
    let parts = R::factor_within(ring, num, bound, cancel)
        .and_then(|top| Ok((top, R::factor_within(ring, den, bound, cancel)?)));
    match parts {
        Ok((top, bottom)) => Ok(Some(Powers::new(ring, top, bottom))),
        Err(FactorError::Cancelled) => Err(EvalError::Interrupted),
        Err(FactorError::Search) => {
            warn(
                ev,
                format!("factoring could take more than {terms} steps of its searches"),
                e,
            );
            Ok(None)
        }
    }
}

/// `Factor[e]`, Xcas's `factor` (CAS-DOC-8 §4, CAS-185, CAS-187): the
/// numerator and the denominator factored over ℚ, or over ℚ(i) when a
/// coefficient is complex (CAS-186), the generators as variables, their
/// common factors cancelled. The result is a constant times powers of
/// irreducible, primitive polynomials, each with its leading coefficient in
/// normal form in the first generator of canonical order, as giac writes
/// them: `factor(x^4-y^4)` is `(x-y)*(x+y)*(x^2+y^2)`, `factor(sin(x)^2-1)`
/// is `(sin(x)-1)*(sin(x)+1)`, and `factor(x^2+(1+2*i)*x+2*i)` is
/// `(x+1)*(x+2*i)`. As giac's default, a quadratic factor in one generator
/// then splits over square roots ([`split_quadratic`]): `factor(x^2-2)` is
/// `(x-sqrt(2))*(x+sqrt(2))`. An infinity is given back.
pub(super) fn factor(ev: &mut Eval<'_>, e: &Expr) -> Native {
    factoring(ev, e, true)
}

/// `FactorRational[e]`, Xcas's `collect` (CAS-187): [`factor`] without the
/// square roots, over ℚ, or ℚ(i) with a complex coefficient: `collect(x^2-2)`
/// is `x^2-2`.
pub(super) fn factor_rational(ev: &mut Eval<'_>, e: &Expr) -> Native {
    factoring(ev, e, false)
}

/// [`factor`], with the square roots when `roots`.
fn factoring(ev: &mut Eval<'_>, e: &Expr, roots: bool) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if infinite(arg) {
        return Ok(Some(arg.clone()));
    }
    let Some(p) = to_factor(ev, e, arg)? else {
        return Ok(None);
    };
    let n = p.generators().len();
    let result = match &p.fractions()[0] {
        Fraction::Integers { num, den } => {
            factored(ev, e, &SparseRing::new(Integers, n), num, den)?.map(|f| {
                if roots {
                    split_expr::<Integers>(&p, &f)
                } else {
                    f.expr(&p)
                }
            })
        }
        Fraction::GaussianIntegers { num, den } => {
            factored(ev, e, &SparseRing::new(GaussianIntegers, n), num, den)?.map(|f| {
                if roots {
                    split_expr::<GaussianIntegers>(&p, &f)
                } else {
                    f.expr(&p)
                }
            })
        }
    };
    let Some(result) = result else {
        return Ok(None);
    };
    ev.step(if roots { "factor" } else { "collect" }, e, &result);
    Ok(Some(result))
}

/// A key of the order of `factors`' lists: the total degree, then the terms
/// with the first generator in canonical order most significant, each by
/// exponents, then real and imaginary parts.
type ListKey = (u64, Vec<(Vec<u64>, (Rational, Rational))>);

/// The order of `factors`' lists: by total degree, then by terms with the
/// first generator in canonical order most significant, as giac lists
/// `x-1` before `x+1` and `x-y` before `x+y`.
fn list_order<R: Factoring>(p: &Sparse<R::Elem>) -> ListKey {
    let mut terms: Vec<(Vec<u64>, (Rational, Rational))> = p
        .terms()
        .map(|(e, c)| (e.into_iter().rev().collect(), R::parts(c)))
        .collect();
    terms.sort_by(|(a, _), (b, _)| b.cmp(a));
    let degree = terms.iter().map(|(e, _)| e.iter().sum()).max().unwrap_or(0);
    (degree, terms)
}

/// The list of `FactorList`: the constant's numerator, reduced, the
/// numerator's factors, the constant's denominator, normalised, and the
/// denominator's factors, each followed by its exponent; a constant alone
/// is `[c, 1]`.
fn listed<R: Factoring>(p: &Polynomials, base: &R, mut powers: Powers<R::Elem>) -> Expr {
    powers
        .factors
        .sort_by_cached_key(|(f, k)| (*k < 0, list_order::<R>(f)));
    let g = base.gcd(&powers.numerator, &powers.denominator);
    let (n, d) = if base.is_zero(&g) {
        (powers.numerator, powers.denominator)
    } else {
        let n = base.div_rem(&powers.numerator, &g).0;
        let (d, unit) = base.normal_form(&base.div_rem(&powers.denominator, &g).0);
        (base.mul(&n, &unit), d)
    };
    let mut items = Vec::with_capacity(2 * powers.factors.len() + 4);
    if powers.factors.is_empty() {
        items.extend([n.over(&d), Expr::from(1)]);
    } else {
        let mut push = |base: Expr, k: i64| items.extend([base, Expr::from(k)]);
        if !n.is_one() {
            push(n.to_expr(), 1);
        }
        for (f, k) in powers.factors.iter().filter(|&&(_, k)| k > 0) {
            push(p.polynomial_expr(f), *k);
        }
        if !d.is_one() {
            push(d.to_expr(), -1);
        }
        for (f, k) in powers.factors.iter().filter(|&&(_, k)| k < 0) {
            push(p.polynomial_expr(f), *k);
        }
    }
    Expr::call(Symbol::LIST, items)
}

/// `FactorList[e]`, Xcas's `factors` (CAS-185): the factorisation of
/// [`factor`] as giac lists it, each factor followed by its exponent, the
/// denominator's negative: `factors(12*x^4-12)` is
/// `[12,1,x-1,1,x+1,1,x^2+1,1]`, and `factors((x^2-1)/(2*x+4))` is
/// `[x-1,1,x+1,1,2,-1,x+2,-1]`. The constant's numerator comes first and its
/// denominator before the denominator's factors, each when it is not 1; a
/// constant alone is `[c, 1]`. giac lists a constant -1 last.
pub(super) fn factor_list(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    let Some(p) = to_factor(ev, e, arg)? else {
        return Ok(None);
    };
    let n = p.generators().len();
    let result = match &p.fractions()[0] {
        Fraction::Integers { num, den } => {
            factored(ev, e, &SparseRing::new(Integers, n), num, den)?
                .map(|f| listed(&p, &Integers, f))
        }
        Fraction::GaussianIntegers { num, den } => {
            factored(ev, e, &SparseRing::new(GaussianIntegers, n), num, den)?
                .map(|f| listed(&p, &GaussianIntegers, f))
        }
    };
    let Some(result) = result else {
        return Ok(None);
    };
    ev.step("factors", e, &result);
    Ok(Some(result))
}

/// `num / den` with their gcd cancelled, as powers of square-free parts.
fn squarefree_parts<R: Numbers + CharacteristicZero>(
    ring: &SparseRing<R>,
    num: &Sparse<R::Elem>,
    den: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Powers<R::Elem>, EvalError> {
    let (n, d) = cancelled(ring, num, den, cancel)?;
    let top = ring.squarefree_with(&n, cancel).map_err(interrupted)?;
    let bottom = ring.squarefree_with(&d, cancel).map_err(interrupted)?;
    Ok(Powers::new(ring, top, bottom))
}

/// `FactorSquareFree[e]`, Xcas's `sqrfree` (CAS-185): `e` over a cancelled
/// common denominator, the numerator and the denominator as products of
/// powers of square-free, pairwise coprime parts (cas-poly's
/// `squarefree`), over ℤ or ℤ\[i\], their signs as [`factor`]'s. giac:
/// `sqrfree(x^3-x^2-x+1)` is `(x+1)*(x-1)^2`, `sqrfree((x^2-1)/(x+1)^2)` is
/// `(x-1)/(x+1)`, and `sqrfree(x^2*y^2+2*x*y^2+y^2)` is `y^2*(x+1)^2`, a
/// content's parts apart from the primitive part's. An infinity is given
/// back.
pub(super) fn factor_square_free(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [arg] = e.args() else {
        return Ok(None);
    };
    if infinite(arg) {
        return Ok(Some(arg.clone()));
    }
    let Some(p) = to_factor(ev, e, arg)? else {
        return Ok(None);
    };
    let cancel = ev.context().interrupt_flag();
    let result = with_fraction!(&p.fractions()[0], p.generators().len(), |ring, num, den| {
        squarefree_parts(&ring, num, den, cancel)?.expr(&p)
    });
    ev.step("sqrfree", e, &result);
    Ok(Some(result))
}

#[cfg(test)]
mod tests;
