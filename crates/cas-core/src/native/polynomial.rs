//! The polynomial builtins of phase 1 (CAS-DOC-8 §10, CAS-140), Xcas's
//! `expand`, `normal`, `gcd`, `lcm`, `quo`, `rem`, `quorem`, `degree`,
//! `coeff`, `lcoeff`, `numer`, `denom`, `content`, `primpart`, `resultant`
//! and `discriminant`, over the conversion of [`crate::poly`] and
//! cas-poly's algorithms.
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
use cas_num::{Complex, Integer, Rational};
use cas_poly::{Cancel, Cancelled, DenseRing, GcdDomain, Sparse, SparseRing};

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
trait Numbers: GcdDomain<Elem: Coefficient> + Default {}

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
fn convert(
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
enum Dense {
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
fn within_dense(ev: &mut Eval<'_>, e: &Expr, p: &Polynomials, rule: Dense) -> Option<()> {
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

#[cfg(test)]
mod tests;
