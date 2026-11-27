//! Discrete wavelet transforms (CAS-634), Xcas's `dwt` and `idwt`, which
//! giac hands to GSL (signalprocessing.cc `discrete_wavelet_transform` and
//! `dwt`; GSL 2.7.1 wavelet/dwt.c). Ported here from GSL, step by step:
//!
//! - `DiscreteWaveletTransform[x, options…]` (`dwt`) and
//!   `InverseDiscreteWaveletTransform[x, options…]` (`idwt`) take a list of
//!   numbers or a matrix, up to four options.
//! - One step on `n` entries (dwt.c `dwt_step`, after Numerical Recipes'
//!   `pwt`): forward, `s_i = Σ_k h1_k a_j` and `s_(i+n/2) = Σ_k g1_k a_j`,
//!   `j = (2i + nc·n − offset + k) mod n`; inverse, `s_j += h2_k a_i + g2_k
//!   a_(i+n/2)` for the same `j`. The transform of a list
//!   (`gsl_wavelet_transform`) steps on `n, n/2, …, 2` forward and `2, 4, …,
//!   n` back; of a matrix (`gsl_wavelet2d_transform`) every row then every
//!   column forward, the columns first back; with `image`
//!   (`gsl_wavelet2d_nstransform`, the "non-standard" form) each step on the
//!   rows then the columns of the leading `i × i` block, for `i = n, n/2,
//!   …, 2`, the reverse back.
//! - The filters are GSL's (wavelet/daubechies.c, haar.c, bspline.c,
//!   `*_init`; [`filters`]), centred with `offset = nc/2`, else 0.
//! - Options, giac's loop, in order: `"daubechies"` (member 4, the
//!   default), `"haar"` (member 2) and `"bspline"` (member 103) set the
//!   wavelet; an integer sets the member, checked against the wavelet set
//!   so far (Daubechies even in 4 … 20, Haar 2, B-spline 103, 105, 202,
//!   204, 206, 208, 301, 303, 305, 307, 309), else "Wavelet not
//!   implemented"; the symbols `center` and `centre` choose the centred
//!   wavelet; `image` (which the language reads as `ColumnSpace`, C0008)
//!   the non-standard 2-D form, ignored for a list. Another string is
//!   "Unsupported wavelet type", anything else "Invalid argument k" (a bare
//!   `haar`). A member a later string leaves unimplemented
//!   (`"haar", "daubechies"` keeps member 2) is GSL's failure, "Failed to
//!   allocate memory" in giac's words. Each stays, with giac's message.
//! - Data: a list of real numbers, or a matrix (equal rows) of them; else
//!   "Data must be numeric" (a list) or "Failed to convert data to GSL
//!   matrix" (a matrix), the call staying. A list pads with zeros to the
//!   next power of 2, a matrix to `n × n`, `n` the power of 2 from the
//!   larger side, the data top left, with giac's warning "data size is not
//!   a power of 2, padding with zeros up to size N" (`N×N` for a matrix).
//!   The empty list is `{}`.
//! - Numbers: Haar on exact data, each entry `a + b√2` with `a` and `b`
//!   rational (rationals among them), is exact in the same form
//!   (`dwt([1,2,3,4],"haar")` is `[5,-2,-√2/2,-√2/2]`, giac's
//!   `[5.0,-2.0,-0.707106781187,-0.707106781187]`), so `idwt` takes it back
//!   to the data exactly. Otherwise real balls at the working precision
//!   plus a guard of `2 log₂ n + 32` bits, rounded to the working precision:
//!   GSL's decimal literals are read as exact numbers, so each ball holds
//!   the transform with GSL's filters, where giac's doubles leave rounding
//!   noise (`2.2·10⁻¹⁶` for an entry of about `10⁻³⁰`). GSL's literals are
//!   not exactly orthonormal: `idwt(dwt(x))` is `x` within about `10⁻²⁸`.
//! - Limits: the padded size (`n`, or `n²`) times the filter length within
//!   `Limits::terms`, checked before anything is built; beyond, the call
//!   stays.
//! - A data entry not structurally numeric is refused before
//!   `approximation`, which would evaluate it again (CAS-650).

#![allow(
    clippy::many_single_char_names,
    reason = "GSL's names: data a of n entries, filters h and g, steps i, j and k"
)]

use cas_num::Rational;
use cas_num::ball::Ball;

use super::approximation;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;
use crate::trace::Warning;

mod filters;
#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// GSL's wavelet families.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Daubechies,
    Haar,
    Bspline,
}

/// The B-spline members GSL implements.
const BSPLINES: [i64; 11] = [103, 105, 202, 204, 206, 208, 301, 303, 305, 307, 309];

/// The options after giac's loop.
struct Options {
    kind: Kind,
    k: i64,
    image: bool,
    centred: bool,
}

/// giac's option loop (`discrete_wavelet_transform`), the options from the
/// second argument on; `Err` with giac's message.
fn options(args: &[Expr]) -> Result<Options, String> {
    let mut o = Options {
        kind: Kind::Daubechies,
        k: 4,
        image: false,
        centred: false,
    };
    for (i, a) in args.iter().enumerate() {
        if let Some(s) = a.as_str() {
            match s {
                "daubechies" => o.kind = Kind::Daubechies,
                "haar" => (o.kind, o.k) = (Kind::Haar, 2),
                "bspline" => (o.kind, o.k) = (Kind::Bspline, 103),
                _ => return Err("Unsupported wavelet type".to_owned()),
            }
        } else if matches!(a.view(), View::Int(_) | View::Integer(_)) {
            let k = a.as_i64().unwrap_or(0);
            let implemented = match o.kind {
                Kind::Daubechies => k % 2 == 0 && (4..=20).contains(&k),
                Kind::Haar => k == 2,
                Kind::Bspline => BSPLINES.contains(&k),
            };
            if !implemented {
                return Err("Wavelet not implemented".to_owned());
            }
            o.k = k;
        } else if a.as_symbol() == Some(Symbol::COLUMN_SPACE) {
            o.image = true;
        } else if matches!(a.as_symbol().map(Symbol::name), Some("center" | "centre")) {
            o.centred = true;
        } else {
            return Err(format!("Invalid argument {}", i + 2));
        }
    }
    Ok(o)
}

/// The arithmetic of a transform: real balls, or `a + b√2` exactly.
trait Ring: Clone {
    fn zero() -> Self;
    fn add(&self, other: &Self, prec: u64) -> Self;
    fn mul(&self, other: &Self, prec: u64) -> Self;
}

impl Ring for Ball {
    fn zero() -> Ball {
        Ball::ZERO
    }

    fn add(&self, other: &Ball, prec: u64) -> Ball {
        Ball::add(self, other, prec)
    }

    fn mul(&self, other: &Ball, prec: u64) -> Ball {
        Ball::mul(self, other, prec)
    }
}

/// `a + b√2`, `a` and `b` rational: Haar's field.
#[derive(Clone)]
struct Surd {
    a: Rational,
    b: Rational,
}

impl Ring for Surd {
    fn zero() -> Surd {
        Surd {
            a: Rational::ZERO,
            b: Rational::ZERO,
        }
    }

    fn add(&self, other: &Surd, _: u64) -> Surd {
        Surd {
            a: &self.a + &other.a,
            b: &self.b + &other.b,
        }
    }

    fn mul(&self, other: &Surd, _: u64) -> Surd {
        let bb = &self.b * &other.b;
        Surd {
            a: &(&self.a * &other.a) + &(&bb + &bb),
            b: &(&self.a * &other.b) + &(&self.b * &other.a),
        }
    }
}

/// A wavelet's filters as numbers, and its offset (`gsl_wavelet`).
struct Wavelet<T> {
    offset: usize,
    h1: Vec<T>,
    g1: Vec<T>,
    h2: Vec<T>,
    g2: Vec<T>,
}

impl<T: Ring> Wavelet<T> {
    fn new(f: &filters::Filters, centred: bool, number: impl Fn(&str) -> T) -> Wavelet<T> {
        let taps = |t: &[&str]| t.iter().map(|s| number(s)).collect();
        Wavelet {
            offset: if centred { f.h1.len() >> 1 } else { 0 },
            h1: taps(f.h1),
            g1: taps(f.g1),
            h2: taps(f.h2),
            g2: taps(f.g2),
        }
    }

    /// dwt.c `dwt_step` on `a[..n]`.
    fn step(&self, a: &mut [T], n: usize, forward: bool, prec: u64) {
        let mut s = vec![T::zero(); n];
        let nmod = self.h1.len() * n - self.offset;
        let (mask, half) = (n - 1, n >> 1);
        for (ii, i) in (0..n).step_by(2).enumerate() {
            if forward {
                let (mut h, mut g) = (T::zero(), T::zero());
                for (k, (h1, g1)) in self.h1.iter().zip(&self.g1).enumerate() {
                    let x = &a[mask & (i + nmod + k)];
                    h = h.add(&h1.mul(x, prec), prec);
                    g = g.add(&g1.mul(x, prec), prec);
                }
                s[ii] = h;
                s[ii + half] = g;
            } else {
                let (ai, ai1) = (&a[ii], &a[ii + half]);
                for (k, (h2, g2)) in self.h2.iter().zip(&self.g2).enumerate() {
                    let j = mask & (i + nmod + k);
                    let t = h2.mul(ai, prec).add(&g2.mul(ai1, prec), prec);
                    s[j] = s[j].add(&t, prec);
                }
            }
        }
        for (x, y) in a.iter_mut().zip(s) {
            *x = y;
        }
    }

    /// The step sizes, `n, n/2, …, 2` forward and back up.
    fn sizes(n: usize, forward: bool) -> Vec<usize> {
        let mut sizes: Vec<usize> = std::iter::successors(Some(n), |&i| Some(i >> 1))
            .take_while(|&i| i >= 2)
            .collect();
        if !forward {
            sizes.reverse();
        }
        sizes
    }

    /// dwt.c `gsl_wavelet_transform` on a list of a power of 2 entries.
    fn transform(&self, a: &mut [T], forward: bool, prec: u64) {
        for i in Self::sizes(a.len(), forward) {
            self.step(a, i, forward, prec);
        }
    }

    /// `f` on the first `n` entries of column `j`.
    fn column(m: &mut [Vec<T>], j: usize, n: usize, f: impl FnOnce(&mut [T])) {
        let mut c: Vec<T> = m[..n].iter().map(|r| r[j].clone()).collect();
        f(&mut c);
        for (r, x) in m.iter_mut().zip(c) {
            r[j] = x;
        }
    }

    /// dwt.c `gsl_wavelet2d_transform`, on an `n × n` matrix.
    fn standard(&self, m: &mut [Vec<T>], forward: bool, prec: u64) {
        let n = m.len();
        let rows = |m: &mut [Vec<T>]| {
            for r in m.iter_mut() {
                self.transform(r, forward, prec);
            }
        };
        let columns = |m: &mut [Vec<T>]| {
            for j in 0..n {
                Self::column(m, j, n, |c| self.transform(c, forward, prec));
            }
        };
        if forward {
            rows(m);
            columns(m);
        } else {
            columns(m);
            rows(m);
        }
    }

    /// dwt.c `gsl_wavelet2d_nstransform`, on an `n × n` matrix.
    fn non_standard(&self, m: &mut [Vec<T>], forward: bool, prec: u64) {
        for i in Self::sizes(m.len(), forward) {
            let rows = |m: &mut [Vec<T>]| {
                for r in &mut m[..i] {
                    self.step(r, i, forward, prec);
                }
            };
            let columns = |m: &mut [Vec<T>]| {
                for j in 0..i {
                    Self::column(m, j, i, |c| self.step(c, i, forward, prec));
                }
            };
            if forward {
                rows(m);
                columns(m);
            } else {
                columns(m);
                rows(m);
            }
        }
    }
}

/// A list's entries, or a matrix's rows, padded with zeros.
enum Grid<T> {
    List(Vec<T>),
    Matrix(Vec<Vec<T>>),
}

fn head(e: &Expr) -> Option<Symbol> {
    e.head().and_then(Expr::as_symbol)
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// Whether `e` is a matrix: a non-empty list of non-empty lists of one
/// length (giac's `ckmatrix`).
fn is_matrix(e: &Expr) -> bool {
    head(e) == Some(Symbol::LIST)
        && !e.args().is_empty()
        && e.args().iter().all(|r| {
            head(r) == Some(Symbol::LIST)
                && !r.args().is_empty()
                && r.args().len() == e.args()[0].args().len()
        })
}

/// `data` read by `number`, padded to `n` (`n × n` for a matrix).
fn read<T: Ring>(
    data: &Expr,
    matrix: bool,
    n: usize,
    mut number: impl FnMut(&Expr) -> Option<T>,
) -> Option<Grid<T>> {
    let mut padded = |row: &[Expr]| {
        let mut v = row.iter().map(&mut number).collect::<Option<Vec<T>>>()?;
        v.resize(n, T::zero());
        Some(v)
    };
    if !matrix {
        return padded(data.args()).map(Grid::List);
    }
    let mut rows = data
        .args()
        .iter()
        .map(|r| padded(r.args()))
        .collect::<Option<Vec<_>>>()?;
    rows.resize(n, vec![T::zero(); n]);
    Some(Grid::Matrix(rows))
}

/// The transform of `grid` with `w`, as an expression by `out`.
fn run<T: Ring>(
    grid: Grid<T>,
    w: &Wavelet<T>,
    image: bool,
    forward: bool,
    prec: u64,
    out: impl Fn(T) -> Expr,
) -> Expr {
    match grid {
        Grid::List(mut a) => {
            w.transform(&mut a, forward, prec);
            list(a.into_iter().map(out).collect())
        }
        Grid::Matrix(mut m) => {
            if image {
                w.non_standard(&mut m, forward, prec);
            } else {
                w.standard(&mut m, forward, prec);
            }
            list(
                m.into_iter()
                    .map(|r| list(r.into_iter().map(&out).collect()))
                    .collect(),
            )
        }
    }
}

/// A real number as a ball at `prec` bits.
fn ball(ev: &Eval<'_>, x: &Expr, prec: u64) -> Option<Ball> {
    if let Some(q) = x.to_rational() {
        return Some(Ball::from_rational(&q, prec));
    }
    let b = match x.view() {
        View::Ball(b) => b.clone(),
        _ if !x.is_numeric() => return None,
        _ => match approximation(ev.context(), x, prec).ok()?.view() {
            View::Ball(b) => b.clone(),
            _ => return None,
        },
    };
    b.is_finite().then_some(b)
}

/// `x` as `a + b√2`: a rational, `√2` (`2^(1/2)`, `2^(−1/2)` = `√2/2`), or
/// sums and products of such, as an exact Haar transform answers.
fn surd_of(x: &Expr) -> Option<Surd> {
    if let Some(a) = x.to_rational() {
        return Some(Surd {
            a,
            b: Rational::ZERO,
        });
    }
    let root = |b: Rational| Surd {
        a: Rational::ZERO,
        b,
    };
    match (head(x)?, x.args()) {
        (Symbol::POWER, [base, exponent]) if base.as_i64() == Some(2) => {
            let half = Rational::new(1.into(), 2.into());
            match exponent.to_rational()? {
                e if e == half => Some(root(Rational::ONE)),
                e if e == -&half => Some(root(half)),
                _ => None,
            }
        }
        (Symbol::PLUS, args) => args
            .iter()
            .try_fold(Surd::zero(), |sum, y| Some(sum.add(&surd_of(y)?, 0))),
        (Symbol::TIMES, args) => args.iter().try_fold(
            Surd {
                a: Rational::ONE,
                b: Rational::ZERO,
            },
            |product, y| Some(product.mul(&surd_of(y)?, 0)),
        ),
        _ => None,
    }
}

/// `a + b√2` as an expression, to evaluate.
fn surd(s: Surd) -> Expr {
    let root = Expr::call(
        Symbol::POWER,
        [Expr::from(2), Expr::from(Rational::new(1.into(), 2.into()))],
    );
    Expr::call(
        Symbol::PLUS,
        [
            Expr::from(s.a),
            Expr::call(Symbol::TIMES, [Expr::from(s.b), root]),
        ],
    )
}

fn warn(ev: &mut Eval<'_>, e: &Expr, message: String) {
    ev.warn(Warning::Native {
        message: message.into(),
        expr: e.clone(),
    });
}

/// The transform of `data` padded to `p` (`p × p`), or giac's message for
/// data that is not numeric or a member GSL lacks.
fn compute(
    ev: &mut Eval<'_>,
    data: &Expr,
    matrix: bool,
    p: usize,
    o: &Options,
    forward: bool,
) -> Result<Result<Expr, &'static str>, EvalError> {
    let bits = ev.context().working_precision();
    let prec = bits + 2 * u64::from(usize::BITS - p.leading_zeros()) + 32;
    let exact = if o.kind == Kind::Haar {
        read(data, matrix, p, surd_of)
    } else {
        None
    };
    let balls = match exact {
        Some(_) => None,
        None => match read(data, matrix, p, |x| ball(ev, x, prec)) {
            Some(grid) => Some(grid),
            None if matrix => return Ok(Err("Failed to convert data to GSL matrix")),
            None => return Ok(Err("Data must be numeric")),
        },
    };
    let Some(f) = filters::filters(o.kind, o.k) else {
        return Ok(Err("Failed to allocate memory"));
    };
    if let Some(grid) = exact {
        // Haar's filters are all ±M_SQRT1_2 = ±√2/2.
        let half = Rational::new(1.into(), 2.into());
        let w = Wavelet::new(&f, o.centred, |s| Surd {
            a: Rational::ZERO,
            b: if s.starts_with('-') {
                -&half
            } else {
                half.clone()
            },
        });
        return ev
            .eval(&run(grid, &w, o.image, forward, prec, surd))
            .map(Ok);
    }
    let grid = balls.expect("read above");
    let w = Wavelet::new(&f, o.centred, |s| {
        Ball::from_decimal(s, prec).expect("GSL's literals")
    });
    Ok(Ok(run(grid, &w, o.image, forward, prec, |b| {
        Expr::from(b.round_to_accuracy(bits))
    })))
}

/// `DiscreteWaveletTransform[x, options…]` and
/// `InverseDiscreteWaveletTransform[x, options…]`.
pub(super) fn transform(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let forward = match head(e) {
        Some(Symbol::DISCRETE_WAVELET_TRANSFORM) => true,
        Some(Symbol::INVERSE_DISCRETE_WAVELET_TRANSFORM) => false,
        _ => return Ok(None),
    };
    let [data, rest @ ..] = e.args() else {
        return Ok(None);
    };
    if head(data) != Some(Symbol::LIST) || rest.len() > 4 {
        return Ok(None);
    }
    let o = match options(rest) {
        Ok(o) => o,
        Err(message) => {
            warn(ev, e, message);
            return Ok(None);
        }
    };
    let matrix = is_matrix(data);
    let n = if matrix {
        data.args().len().max(data.args()[0].args().len())
    } else {
        data.args().len()
    };
    if n == 0 {
        return Ok(Some(list(Vec::new())));
    }
    let nc = filters::filters(o.kind, o.k).map_or(1, |f| f.h1.len());
    let Some(p) = n.checked_next_power_of_two() else {
        return Ok(None);
    };
    let entries = if matrix { p.checked_mul(p) } else { Some(p) };
    if entries
        .and_then(|m| m.checked_mul(nc))
        .is_none_or(|m| m > ev.context().limits().terms)
    {
        return Ok(None);
    }
    if p != n {
        let size = if matrix {
            format!("{p}×{p}")
        } else {
            p.to_string()
        };
        warn(
            ev,
            e,
            format!("data size is not a power of 2, padding with zeros up to size {size}"),
        );
    }
    let result = match compute(ev, data, matrix, p, &o, forward)? {
        Ok(result) => result,
        Err(message) => {
            warn(ev, e, message.to_owned());
            return Ok(None);
        }
    };
    ev.step(if forward { "dwt" } else { "idwt" }, e, &result);
    Ok(Some(result))
}
