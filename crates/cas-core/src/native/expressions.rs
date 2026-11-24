//! Xcas's expression commands, second part (CAS-568 to CAS-571): the
//! identifiers and kernels of an expression, giac's order of complexity,
//! plain and boolean evaluation, decimals near 0, and equations turned into
//! lists and back. Each is a port of giac's own command; where giac's order
//! of operands or variables differs from this CAS's canonical order, the
//! answer keeps ours.
//!
//! - [`identifiers`]: `Identifiers[e]`, giac's `lname` and `indets`
//!   (prog.cc `_lname`): the symbols of `e` that are not builtins (giac
//!   skips `pi` and `euler_gamma`), heads of unknown functions among them,
//!   each once, depth first in order of appearance.
//! - [`kernels`]: `Kernels[e]`, giac's `lvar` (sym2poly.cc `lvar`): the
//!   kernels `e` is a rational function of. Sums, products and lists are
//!   walked, a power with an integer exponent gives its base's kernels,
//!   numbers and strings none, and anything else (a symbol but
//!   `Indeterminate`, a function, a power with another exponent) is a kernel.
//! - [`algebraic_variables`]: `AlgebraicVariables[e]`, giac's `algvar`
//!   (sym2poly.cc `alg_lvar`): the kernels as rows, the first the kernels
//!   that are not radicals, then one row per level of radicals (a power whose
//!   exponent is a fraction of denominator at most 4, giac's
//!   `MAX_ALG_EXT_ORDER_SIZE`), with the variables under it; the rows of two
//!   radicals merged as giac's `ext_glue_matrices` merges them.
//! - [`reorder`]: `Reorder[e, vars]`, giac's `reorder` (misc.cc): `e` as a
//!   rational function, giac's `e2r` then `r2e` with `vars` first, which is
//!   `Together[e]` here, as canonical order does not follow `vars`; a number
//!   stays.
//! - [`evaluate`]: `Evaluate[e]` and `Evaluate[e, n]`, giac's `eval` and
//!   `evalm`: `e`, evaluated (matrix powers evaluate already).
//! - [`evaluate_boolean`]: `EvaluateBoolean[e]`, giac's `evalb` (prog.cc):
//!   each equation and comparison `a op b` in `e` (through `And`, `Or` and
//!   `Not`) becomes `Together[a − b] op 0`, an equation `SameQ`; the result,
//!   or its `N`, decides: `True`, `False`, or a number, false when 0. Lists
//!   map; an undecided test stays.
//! - [`epsilon_to_zero`]: `EpsilonToZero[e]` and `EpsilonToZero[e, eps]`,
//!   giac's `epsilon2zero` (misc.cc): each decimal of absolute value at most
//!   `eps` (giac's `epsilon`, `10⁻¹²`) becomes 0, the real and imaginary
//!   parts of a complex decimal apart. It holds its arguments, as evaluation
//!   would merge `1e-13 + 5` into one decimal giac keeps apart: the held
//!   expression is cleared, evaluated, and cleared again.
//! - [`infinite_q`]: `InfiniteQ[e]`, giac's `isinf`: whether `e` is an
//!   infinity, signed, directed or complex.
//! - [`indeterminate_q`]: `IndeterminateQ[e]`, giac's `isnan`: whether `e`
//!   is `Indeterminate`, or a list whose first entry is (giac's `is_undef`).
//! - [`compare`]: `Compare[a, b]`, giac's `compare` (misc.cc), 1 when `a` is
//!   less complex than `b` by giac's `gen::islesscomplexthan` (gen.cc), else
//!   0. giac ranks by type first: integers of 32 bits, decimals, larger
//!   integers, complex numbers, symbols, lists, applications, fractions,
//!   strings, 0 below everything; then integers by absolute value, decimals
//!   with the larger absolute value first, complex numbers by modulus, real
//!   and imaginary parts, the larger first, symbols by name, lists by length
//!   then entries, fractions by printed text, strings by text, applications
//!   by their arguments. Between two different heads giac compares its
//!   operators' addresses, which no other program can reproduce: ours
//!   compares the heads' names (CAS-569).
//! - [`equation_to_difference`]: `EquationToDifference[e]`, giac's
//!   `equal2diff`: `a − b` for `a = b`, mapped over a list; other parts
//!   stay.
//! - [`equation_to_list`]: `EquationToList[e]`, giac's `equal2list`: the
//!   sequence `a, b` for `a = b` (spliced in a list, as giac's is), `[e, 0]`
//!   for another `e`, mapped over a list.
//! - [`list_to_equations`]: `ListToEquations[values, vars]`, giac's
//!   `list2exp` (ti89.cc, through `solvepostprocess`): `x = v` for each value
//!   of a variable, or `x = v₁ and y = v₂` for each row of values of a list of
//!   variables, joined by `Or`; one alone is itself, none is `False`.
//! - [`equations_to_list`]: `EquationsToList[e]`, giac's `exp2list`: the
//!   right sides of the equations of `e`, `Or` and `And` (and lists) giving
//!   lists of them; 0 gives `[]` and another single answer a list of one.
//!
//! Equations are `Equation` (the language's `=`) and `Equal` (`==`) alike.

use std::cmp::Ordering;

use cas_num::Rational;
use cas_num::ball::{Ball, Float, Round};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// giac's `MAX_ALG_EXT_ORDER_SIZE`: the largest root order `algvar` takes
/// for an algebraic extension.
const MAX_ROOT_ORDER: i64 = 4;

/// giac's `epsilon`, the default bound of `epsilon2zero` and of `trim`.
pub(super) const EPSILON: f64 = 1e-12;

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

fn boolean(b: bool) -> Expr {
    Expr::from(if b { Symbol::TRUE } else { Symbol::FALSE })
}

fn head_symbol(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn is_list(e: &Expr) -> bool {
    head_symbol(e) == Some(Symbol::LIST)
}

/// The two sides of an equation, `Equation` or `Equal`.
fn sides(e: &Expr) -> Option<(&Expr, &Expr)> {
    match (head_symbol(e), e.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [a, b]) => Some((a, b)),
        _ => None,
    }
}

fn push_new(out: &mut Vec<Expr>, e: &Expr) {
    if !out.contains(e) {
        out.push(e.clone());
    }
}

fn identifiers_of(e: &Expr, out: &mut Vec<Expr>) {
    match e.view() {
        View::Symbol(s) if !s.is_builtin() => push_new(out, e),
        View::Apply { head, args } => {
            identifiers_of(head, out);
            for a in args {
                identifiers_of(a, out);
            }
        }
        _ => {}
    }
}

/// `Identifiers[e, …]`: the identifiers of the arguments, in order.
pub(super) fn identifiers(_: &mut Eval<'_>, e: &Expr) -> Native {
    let mut out = Vec::new();
    for a in e.args() {
        identifiers_of(a, &mut out);
    }
    Ok(Some(list(out)))
}

fn is_integer(e: &Expr) -> bool {
    matches!(e.view(), View::Int(_) | View::Integer(_))
}

/// giac's `lvar`: the kernels of `e` appended to `out`.
fn kernels_of(e: &Expr, out: &mut Vec<Expr>) {
    match e.view() {
        View::Symbol(Symbol::INDETERMINATE) => {}
        View::Symbol(_) => push_new(out, e),
        View::Apply { head, args } => match (head.as_symbol(), args) {
            (Some(Symbol::PLUS | Symbol::TIMES | Symbol::LIST), _) => {
                for a in args {
                    kernels_of(a, out);
                }
            }
            (Some(Symbol::POWER), [base, exponent]) if is_integer(exponent) => {
                kernels_of(base, out);
            }
            _ => push_new(out, e),
        },
        _ => {}
    }
}

/// `Kernels[e]`.
pub(super) fn kernels(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let mut out = Vec::new();
    kernels_of(x, &mut out);
    Ok(Some(list(out)))
}

/// The rows of `algvar`.
type Rows = Vec<Vec<Expr>>;

/// The argument of an algebraic extension: `b` for `b^(p/q)`, `q ≤ 4`.
fn algebraic_argument(e: &Expr) -> Option<&Expr> {
    let (Some(Symbol::POWER), [base, exponent]) = (head_symbol(e), e.args()) else {
        return None;
    };
    let View::Rational(r) = exponent.view() else {
        return None;
    };
    let den = r.denom().to_i64()?;
    (den.abs() <= MAX_ROOT_ORDER).then_some(base)
}

/// giac's `equalposmat`: whether `e` is in a row of `m`.
fn in_rows(m: &Rows, e: &Expr) -> bool {
    m.iter().any(|row| row.contains(e))
}

/// giac's `addfirstrow`.
fn add_first_row(e: Expr, m: &mut Rows) {
    match m.first_mut() {
        Some(row) => row.push(e),
        None => m.push(vec![e]),
    }
}

/// giac's `ext_glue_matrices`: the rows of `a` and `b` merged, `b`'s
/// extensions kept.
fn glue(a: Rows, b: Rows) -> Rows {
    if a.len() > b.len() {
        return glue(b, a);
    }
    if a.len() == b.len() {
        for (ra, rb) in a.iter().zip(&b) {
            if ra.len() < rb.len() {
                break;
            }
            if ra.len() > rb.len() {
                return glue(b, a);
            }
        }
    }
    if a.is_empty() || b.is_empty() || a == b {
        return b;
    }
    let res: Rows = b[1..].to_vec();
    let mut out: Rows = Vec::new();
    for (k, row) in a.iter().enumerate() {
        let temp: Vec<Expr> = row.iter().filter(|x| !in_rows(&res, x)).cloned().collect();
        if k == 0 || !temp.is_empty() {
            out.push(temp);
        }
    }
    for x in &b[0] {
        if !in_rows(&out, x) {
            add_first_row(x.clone(), &mut out);
        }
    }
    out.extend(res);
    out
}

/// giac's `alg_lvar(e, m)`.
fn algebraic_rows(e: &Expr, m: &mut Rows) {
    let mut kernels = Vec::new();
    kernels_of(e, &mut kernels);
    for k in kernels {
        if let Some(argument) = algebraic_argument(&k) {
            let mut inner: Rows = vec![Vec::new()];
            algebraic_rows(argument, &mut inner);
            let mut extension: Rows = vec![Vec::new()];
            if inner.len() > 1 || (inner.len() == 1 && !inner[0].is_empty()) {
                extension.extend(inner);
            }
            *m = glue(extension, std::mem::take(m));
        } else if !in_rows(m, &k) {
            add_first_row(k, m);
        }
    }
}

/// `AlgebraicVariables[e]`.
pub(super) fn algebraic_variables(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let mut m = Rows::new();
    algebraic_rows(x, &mut m);
    Ok(Some(list(m.into_iter().map(list).collect())))
}

/// `Reorder[e, vars]`.
pub(super) fn reorder(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x, vars] = e.args() else {
        return Ok(None);
    };
    if x.is_number() {
        return Ok(Some(x.clone()));
    }
    if !is_list(vars) {
        return Ok(None);
    }
    Ok(Some(Expr::call(Symbol::TOGETHER, [x.clone()])))
}

/// `Evaluate[e]`, `Evaluate[e, n]`.
pub(super) fn evaluate(_: &mut Eval<'_>, e: &Expr) -> Native {
    match e.args() {
        [x] | [x, _] => Ok(Some(x.clone())),
        _ => Ok(None),
    }
}

/// giac's `equaltosame` and `normal` on a test: each comparison against 0.
fn normal_test(e: &Expr) -> Expr {
    let together = |a: &Expr, b: &Expr| Expr::call(Symbol::TOGETHER, [a.clone() - b.clone()]);
    match (head_symbol(e), e.args()) {
        (Some(Symbol::EQUATION | Symbol::EQUAL), [a, b]) => {
            Expr::call(Symbol::SAME_Q, [together(a, b), Expr::from(0)])
        }
        (
            Some(s @ (Symbol::LESS | Symbol::LESS_EQUAL | Symbol::GREATER | Symbol::GREATER_EQUAL)),
            [a, b],
        ) => Expr::call(s, [together(a, b), Expr::from(0)]),
        (Some(s @ (Symbol::AND | Symbol::OR | Symbol::NOT)), args) => {
            Expr::call(s, args.iter().map(normal_test))
        }
        _ => e.clone(),
    }
}

/// A decided test: `True`, `False` or a number, false when 0.
fn decided(e: &Expr) -> Option<bool> {
    match e.view() {
        View::Symbol(Symbol::TRUE) => Some(true),
        View::Symbol(Symbol::FALSE) => Some(false),
        _ => e.is_number().then(|| !is_zero_number(e)),
    }
}

/// `EvaluateBoolean[e]`.
pub(super) fn evaluate_boolean(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if is_list(x) {
        return Ok(Some(list(
            x.args()
                .iter()
                .map(|t| Expr::call(Symbol::EVALUATE_BOOLEAN, [t.clone()]))
                .collect(),
        )));
    }
    let test = ev.eval(&normal_test(x))?;
    if let Some(b) = decided(&test) {
        return Ok(Some(boolean(b)));
    }
    let numeric = ev.eval(&Expr::call(Symbol::N, [test]))?;
    Ok(decided(&numeric).map(boolean))
}

/// `e` with each decimal of absolute value at most `eps` replaced by 0.
fn clear(ev: &mut Eval<'_>, e: &Expr, eps: f64) -> Result<Expr, EvalError> {
    let small = |b: &Ball| b.mid().to_f64().abs() <= eps;
    Ok(match e.view() {
        View::Ball(b) if small(b) => Expr::from(0),
        View::ComplexBall(z) => {
            let part = |b: &Ball| {
                if small(b) {
                    Expr::from(0)
                } else {
                    Expr::from(b.clone())
                }
            };
            let (re, im) = (part(z.re()), part(z.im()));
            if re.as_i64() == Some(0) || im.as_i64() == Some(0) {
                re + im * Expr::from(Symbol::I)
            } else {
                e.clone()
            }
        }
        View::Apply { head, .. } if head.as_symbol() == Some(Symbol::DECIMAL) => {
            let ball = ev.eval(e)?;
            if ball == *e {
                ball
            } else {
                clear(ev, &ball, eps)?
            }
        }
        View::Apply { head, args } => {
            let mut cleared = Vec::with_capacity(args.len());
            for a in args {
                cleared.push(clear(ev, a, eps)?);
            }
            Expr::apply(head, &cleared)
        }
        _ => e.clone(),
    })
}

/// `EpsilonToZero[e]`, `EpsilonToZero[e, eps]`, holding its arguments.
pub(super) fn epsilon_to_zero(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let (x, eps) = match e.args() {
        [x] => (x, Some(EPSILON)),
        [x, eps] => {
            let eps = ev.eval(&Expr::call(Symbol::N, [eps.clone()]))?;
            let eps = match eps.view() {
                View::Ball(b) => Some(b.mid().to_f64()),
                _ => None,
            };
            (x, eps.filter(|eps| *eps > 0.0))
        }
        _ => return Ok(None),
    };
    let Some(eps) = eps else {
        return Ok(None);
    };
    let held = clear(ev, x, eps)?;
    let value = ev.eval(&held)?;
    Ok(Some(clear(ev, &value, eps)?))
}

/// `InfiniteQ[e]`.
pub(super) fn infinite_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let infinite = matches!(
        x.as_symbol(),
        Some(Symbol::INFINITY | Symbol::COMPLEX_INFINITY)
    ) || head_symbol(x) == Some(Symbol::DIRECTED_INFINITY);
    Ok(Some(boolean(infinite)))
}

/// `IndeterminateQ[e]`.
pub(super) fn indeterminate_q(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let first = if is_list(x) {
        x.args().first()
    } else {
        Some(x)
    };
    Ok(Some(boolean(
        first.and_then(Expr::as_symbol) == Some(Symbol::INDETERMINATE),
    )))
}

/// Whether a number is 0, a decimal by its midpoint.
fn is_zero_number(e: &Expr) -> bool {
    match e.view() {
        View::Int(n) => n == 0,
        View::Ball(b) => b.mid().is_zero(),
        View::ComplexBall(z) => z.re().mid().is_zero() && z.im().mid().is_zero(),
        _ => false,
    }
}

/// A decimal's midpoint, in doubles: exact values of huge decimals would
/// not fit.
fn decimal(e: &Expr) -> Option<f64> {
    match e.view() {
        View::Ball(b) => Some(b.mid().to_f64()),
        _ => None,
    }
}

fn double(r: &Rational) -> f64 {
    Float::from_rational(r, 64, Round::Nearest).0.to_f64()
}

/// The parts of a complex number, exact or in doubles: `(re, im)`.
enum Parts {
    Exact(Rational, Rational),
    Doubles(f64, f64),
}

fn parts(e: &Expr) -> Option<Parts> {
    match e.view() {
        View::Complex(z) => Some(Parts::Exact(z.re().clone(), z.im().clone())),
        View::ComplexBall(z) => Some(Parts::Doubles(z.re().mid().to_f64(), z.im().mid().to_f64())),
        _ => None,
    }
}

fn doubles(p: Parts) -> (f64, f64) {
    match p {
        Parts::Exact(re, im) => (double(&re), double(&im)),
        Parts::Doubles(re, im) => (re, im),
    }
}

/// giac's order of two complex numbers: the larger modulus, then real
/// part, then imaginary part first.
fn less_complex_numbers(a: Parts, b: Parts) -> bool {
    match (a, b) {
        (Parts::Exact(ra, ia), Parts::Exact(rb, ib)) => {
            let (na, nb) = (&(&ra * &ra) + &(&ia * &ia), &(&rb * &rb) + &(&ib * &ib));
            if na != nb {
                na > nb
            } else if ra != rb {
                ra > rb
            } else {
                ia >= ib
            }
        }
        (a, b) => {
            let ((ra, ia), (rb, ib)) = (doubles(a), doubles(b));
            let (na, nb) = (ra.hypot(ia), rb.hypot(ib));
            match na.total_cmp(&nb).then(ra.total_cmp(&rb)) {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => ia >= ib,
            }
        }
    }
}

/// giac's type of `e`, in giac's order: `_INT_` 0, `_DOUBLE_` 1, `_ZINT`
/// 2, `_CPLX` 4, `_IDNT` 6, `_VECT` 7, `_SYMB` 8, `_FRAC` 10, `_STRNG` 12.
fn giac_type(e: &Expr) -> u8 {
    match e.view() {
        View::Int(n) if i32::try_from(n).is_ok() => 0,
        View::Ball(_) => 1,
        View::Int(_) | View::Integer(_) => 2,
        View::Complex(_) | View::ComplexBall(_) => 4,
        View::Symbol(_) => 6,
        View::Apply { head, .. } if head.as_symbol() == Some(Symbol::LIST) => 7,
        View::Apply { .. } => 8,
        View::Rational(_) => 10,
        View::String(_) => 12,
    }
}

/// giac's `compare_VECT`.
fn less_complex_lists(v: &[Expr], w: &[Expr]) -> bool {
    if v.len() != w.len() {
        return v.len() < w.len();
    }
    v.iter()
        .zip(w)
        .find(|(a, b)| a != b)
        .is_some_and(|(a, b)| less_complex(a, b))
}

/// The arguments of an application as giac's `feuille`: the one argument,
/// or the list of them.
fn leaf(e: &Expr) -> Expr {
    match e.args() {
        [x] => x.clone(),
        args => list(args.to_vec()),
    }
}

/// The name a head compares by, ours (CAS-569).
fn head_name(e: &Expr) -> String {
    match e.head().map(Expr::view) {
        Some(View::Symbol(s)) => s.name().to_owned(),
        _ => format!("{:?}", e.head()),
    }
}

/// giac's `gen::islesscomplexthan`.
fn less_complex(a: &Expr, b: &Expr) -> bool {
    let (ta, tb) = (giac_type(a), giac_type(b));
    if ta < 6 && is_zero_number(a) {
        return false;
    }
    if tb < 6 && is_zero_number(b) {
        return true;
    }
    if ta != tb {
        return ta < tb;
    }
    if a == b {
        return false;
    }
    match ta {
        // Integers: the opposite of a positive one first, else by
        // absolute value, `_ZINT`s with ties to the first.
        0 | 2 => {
            let (Some(x), Some(y)) = (a.to_rational(), b.to_rational()) else {
                return false;
            };
            if (&x + &y).is_zero() {
                return x.is_positive();
            }
            if ta == 0 {
                x.abs() < y.abs()
            } else {
                y.abs() >= x.abs()
            }
        }
        // Decimals: the larger absolute value first.
        1 => {
            let (Some(x), Some(y)) = (decimal(a), decimal(b)) else {
                return false;
            };
            if x.total_cmp(&-y) == Ordering::Equal {
                return x > 0.0;
            }
            x.abs() >= y.abs()
        }
        4 => match (parts(a), parts(b)) {
            (Some(x), Some(y)) => less_complex_numbers(x, y),
            _ => false,
        },
        6 => a.as_symbol().map(Symbol::name) < b.as_symbol().map(Symbol::name),
        7 => less_complex_lists(a.args(), b.args()),
        8 => {
            if a.head() == b.head() {
                less_complex(&leaf(a), &leaf(b))
            } else {
                head_name(a) < head_name(b)
            }
        }
        10 => {
            let text = |e: &Expr| {
                e.to_rational()
                    .map(|r| format!("{}/{}", r.numer(), r.denom()))
            };
            text(a) < text(b)
        }
        12 => a.as_str() < b.as_str(),
        _ => false,
    }
}

/// `Compare[a, b]`: 1 or 0.
pub(super) fn compare(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [a, b] = e.args() else {
        return Ok(None);
    };
    Ok(Some(Expr::from(i64::from(less_complex(a, b)))))
}

/// `EquationToDifference[e]`.
pub(super) fn equation_to_difference(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let one = |g: &Expr| sides(g).map_or_else(|| g.clone(), |(a, b)| a.clone() - b.clone());
    Ok(Some(if is_list(x) {
        list(x.args().iter().map(one).collect())
    } else {
        one(x)
    }))
}

/// `EquationToList[e]`.
pub(super) fn equation_to_list(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    let one = |g: &Expr| match sides(g) {
        Some((a, b)) => Expr::call(Symbol::SEQUENCE, [a.clone(), b.clone()]),
        None => list(vec![g.clone(), Expr::from(0)]),
    };
    Ok(Some(if is_list(x) {
        list(x.args().iter().map(one).collect())
    } else {
        one(x)
    }))
}

/// `ListToEquations[values, vars]`.
pub(super) fn list_to_equations(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [values, vars] = e.args() else {
        return Ok(None);
    };
    if !is_list(values) {
        return Ok(None);
    }
    if values.args().is_empty() {
        return Ok(Some(boolean(false)));
    }
    let equation = |x: &Expr, v: &Expr| Expr::call(Symbol::EQUATION, [x.clone(), v.clone()]);
    let alternatives: Vec<Expr> = if is_list(vars) {
        let n = vars.args().len();
        let is_row = |r: &Expr| is_list(r) && r.args().len() == n;
        let rows = if values.args().iter().all(is_row) {
            values.args().to_vec()
        } else if is_row(values) {
            vec![values.clone()]
        } else {
            return Ok(None);
        };
        rows.iter()
            .map(|row| {
                let eqs: Vec<Expr> = vars
                    .args()
                    .iter()
                    .zip(row.args())
                    .map(|(x, v)| equation(x, v))
                    .collect();
                if let [one] = &eqs[..] {
                    one.clone()
                } else {
                    Expr::call(Symbol::AND, eqs)
                }
            })
            .collect()
    } else {
        values.args().iter().map(|v| equation(vars, v)).collect()
    };
    Ok(Some(if let [one] = &alternatives[..] {
        one.clone()
    } else {
        Expr::call(Symbol::OR, alternatives)
    }))
}

/// giac's `remove_and` and `_rightapply`: right sides, alternatives and
/// conjunctions as lists.
fn right_sides(e: &Expr) -> Expr {
    if let Some((_, b)) = sides(e) {
        return b.clone();
    }
    match head_symbol(e) {
        Some(Symbol::OR | Symbol::AND | Symbol::LIST) => {
            list(e.args().iter().map(right_sides).collect())
        }
        _ => e.clone(),
    }
}

/// `EquationsToList[e]`.
pub(super) fn equations_to_list(_: &mut Eval<'_>, e: &Expr) -> Native {
    let [x] = e.args() else {
        return Ok(None);
    };
    if x.as_i64() == Some(0) {
        return Ok(Some(list(Vec::new())));
    }
    let r = right_sides(x);
    Ok(Some(if is_list(&r) { r } else { list(vec![r]) }))
}

#[cfg(test)]
mod tests;
