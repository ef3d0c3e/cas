//! B-splines and spline fits (CAS-596): Xcas's `bspline` and `fitspline`,
//! after giac's optimization.cc (`deBoor`, `bspline`, `_bspline`,
//! `subinterval_index`, `eval_basis_splines`, `_fitspline`).
//!
//! - **`BSpline[c, x, p, piecewise]`** (`_bspline`): the B-spline of degree
//!   `p` (3 when left out) with control points `c`, one row per piece, each
//!   row the piece's coordinates as polynomials in `x`.
//!   - The points: a matrix's rows; else each number `z` is `[Re z, Im z]`
//!     and another expression free of `I` is `[e, 0]` (giac's
//!     `coordonnees`, which takes identifiers real).
//!   - The knots: `x` alone, `n` equally spaced on `[0, 1]`; `x == a..b`, on
//!     `[a, b]`; `x == [t_0, …, t_{n−1}]`, those, strictly ascending; each
//!     exact, `n = |c| − p + 1 ≥ 2`. They are clamped: `t_0` and `t_{n−1}`
//!     repeated `p` more times, so the first piece starts at the first
//!     point and the last ends at the last.
//!   - Each of the `n − 1` pieces comes by de Boor's recurrence run on unit
//!     control points: the weights of the `p + 1` points that act on it, as
//!     exact polynomials in `x`; a coordinate is `Expand` of their sum.
//!   - giac's collapse: one column answers the list of pieces, and one
//!     piece of one column the polynomial itself. With `piecewise` each
//!     column is one `Piecewise`, giac's conditions: `t_0 ≤ x < t_1` first,
//!     `x ≤ t_{n−1}` last, `x < t_i` between.
//! - **`SplineFit[data, x | points, opts]`** (`_fitspline`): the least-
//!   squares spline of degree `p` through the rows `(s, y_1, …)` of a numeric
//!   matrix, each column `y_j` a series.
//!   - `n` intervals, `n = ⌊√m + 1/2⌋` for `m` rows, or an integer option,
//!     or `breakpoint = [b_0, …, b_n]` (sorted, strictly ascending, covering
//!     the data); otherwise equally spaced on `[min s, max s]`; clamped as
//!     above; `degree = p` (3 when left out).
//!   - The `n + p` basis splines at the `s_i` (giac's
//!     `eval_basis_splines`) make the matrix `A`; the control points `Q`
//!     solve `(AᵀA) Q = AᵀY` exactly, by Gauss–Jordan over the rationals
//!     (giac takes a banded Cholesky in doubles).
//!   - `derive = k` differentiates `k` times through the control points'
//!     differences `p (Q_{i+1} − Q_i)/(t_{i+p+1} − t_{i+1})`, the degree and
//!     the knots one less each time.
//!   - Points answer the spline's values by de Boor: one series a flat
//!     list, several one row per point. A variable answers `BSpline` of
//!     `Q` on the fit's own knots, `piecewise` passed on. `close` and
//!     `uniform` are taken and ignored, as giac with GSL.
//!   - Exact data answer exactly. A decimal among the data, breakpoints or
//!     points is taken as its exact rational and the answers are balls at
//!     the working precision.
//!
//! What stays: giac's errors (an argument of another form, a degree not a
//! positive integer, a wrong knot count, knots not strictly ascending,
//! `n < 1`, `n + p − 1 ≥ m`, `derive ≥ p`, breakpoints not covering the
//! data or conflicting with an integer option, a singular system), decimal
//! knots, a point outside `[min s, max s]` (giac's answer there is the
//! spline's last polynomial piece continued, not a fit), and any size past
//! the caps below. Each call that answers takes one step, its builtin's
//! name.
//!
//! Caps (C0019, C0048): [`MAX_POINTS`] control points, data rows and
//! points; [`MAX_DEGREE`]; [`MAX_COLUMNS`] coordinates or series;
//! [`MAX_SYMBOLIC_CONTROL_POINTS`] when a coordinate is not a number;
//! [`MAX_SPLINE_WORK`] for pieces `· (p + 1)⁴` (the weights) and pieces
//! `· (p + 1)² ·` columns (the sums); [`MAX_FIT_WORK`] for `m (n + p)²`,
//! `(n + p)²` · series and points `· (p + 1)²` · series; decimals of at most
//! [`MAX_DECIMAL_BITS`] bits.

use cas_num::Rational;
use cas_num::ball::Ball;

use super::interpolants::matrix_rows;
use super::{MAX_POINTS, Native, is_list};
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// The highest degree.
const MAX_DEGREE: usize = 32;

/// The most coordinates of a control point, or series of a fit.
const MAX_COLUMNS: usize = 64;

/// The most control points when a coordinate is not a number.
const MAX_SYMBOLIC_CONTROL_POINTS: usize = 32;

/// The most pieces `· (p + 1)⁴`, and pieces `· (p + 1)² ·` columns.
const MAX_SPLINE_WORK: usize = 1 << 18;

/// The most `m (n + p)²`, `(n + p)²` · series, and points `· (p + 1)²` ·
/// series of a fit.
const MAX_FIT_WORK: usize = 1 << 17;

/// The most bits of a decimal's exact numerator or denominator.
const MAX_DECIMAL_BITS: u64 = 1 << 12;

/// A polynomial in `x`, its coefficients from the constant up.
type Poly = Vec<Rational>;

fn count(n: usize) -> Rational {
    Rational::from(i64::try_from(n).expect("a size in memory"))
}

/// A real number exactly, and whether it was a decimal; `None` for another
/// expression, or a decimal past [`MAX_DECIMAL_BITS`] (C0028).
pub(super) fn real(e: &Expr) -> Option<(Rational, bool)> {
    if let Some(q) = e.to_rational() {
        return Some((q, false));
    }
    match e.view() {
        View::Ball(b) if b.is_finite() => b
            .mid()
            .to_rational_within(MAX_DECIMAL_BITS)
            .map(|q| (q, true)),
        _ => None,
    }
}

/// `n ≥ 2` equally spaced values from `a` to `b`.
fn linspace(a: &Rational, b: &Rational, n: usize) -> Vec<Rational> {
    let step = &(b - a) / &count(n - 1);
    (0..n).map(|i| a + &(&step * &count(i))).collect()
}

/// The knots clamped: the first and the last repeated `p` more times.
pub(super) fn padded(knots: &[Rational], p: usize) -> Vec<Rational> {
    let n = knots.len();
    (0..n + 2 * p)
        .map(|i| knots[i.saturating_sub(p).min(n - 1)].clone())
        .collect()
}

/// `out += (c0 + c1 x) q`.
fn add_linear_times(out: &mut Poly, q: &Poly, c0: &Rational, c1: &Rational) {
    if out.len() < q.len() + 1 {
        out.resize(q.len() + 1, Rational::from(0));
    }
    for (i, v) in q.iter().enumerate() {
        out[i] = &out[i] + &(c0 * v);
        out[i + 1] = &out[i + 1] + &(c1 * v);
    }
}

/// The weights of the control points `c_{k−p}, …, c_k` on the piece
/// `[t_k, t_{k+1}]` of the padded knots `t`, as polynomials in `x`: giac's
/// `deBoor` on unit control points.
fn weights(k: usize, t: &[Rational], p: usize) -> Vec<Poly> {
    let one = Rational::from(1);
    // d[j][l]: the weight of local control point l in d_j.
    let mut d: Vec<Vec<Poly>> = (0..=p)
        .map(|j| {
            (0..=p)
                .map(|l| {
                    if l == j {
                        vec![one.clone()]
                    } else {
                        Vec::new()
                    }
                })
                .collect()
        })
        .collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let lo = &t[j + k - p];
            let slope = &one / &(&t[j + 1 + k - r] - lo);
            // a = slope (x − lo); 1 − a = (1 + slope lo) − slope x.
            let a0 = -&(&slope * lo);
            let b0 = &one - &a0;
            let b1 = -&slope;
            d[j] = d[j]
                .iter()
                .zip(&d[j - 1])
                .map(|(u, v)| {
                    let mut w = Poly::new();
                    add_linear_times(&mut w, u, &a0, &slope);
                    add_linear_times(&mut w, v, &b0, &b1);
                    w
                })
                .collect();
        }
    }
    d.swap_remove(p)
}

/// The B-spline of degree `p` with control points `points` (rows of
/// coordinates) on the knots, in `x`, collapsed as giac's or one
/// `Piecewise` per column; `None` past the caps.
fn spline(
    ev: &mut Eval<'_>,
    x: &Expr,
    knots: &[Rational],
    points: &[Vec<Expr>],
    p: usize,
    piecewise: bool,
) -> Result<Option<Expr>, EvalError> {
    let pieces = knots.len() - 1;
    let columns = points[0].len();
    let square = (p + 1) * (p + 1);
    if pieces * square * square > MAX_SPLINE_WORK || pieces * square * columns > MAX_SPLINE_WORK {
        return Ok(None);
    }
    let t = padded(knots, p);
    let mut rows = Vec::with_capacity(pieces);
    for k in p..p + pieces {
        let w = weights(k, &t, p);
        let local = &points[k - p..=k];
        let row = (0..columns)
            .map(|c| {
                let mut terms = Vec::new();
                for (poly, point) in w.iter().zip(local) {
                    for (i, q) in poly.iter().enumerate() {
                        if !q.is_zero() {
                            let power = x.pow(i64::try_from(i).expect("a degree below MAX_DEGREE"));
                            terms.push(Expr::from(q.clone()) * point[c].clone() * power);
                        }
                    }
                }
                let sum = if terms.is_empty() {
                    Expr::from(0)
                } else {
                    Expr::call(Symbol::EXPAND, [Expr::call(Symbol::PLUS, terms)])
                };
                ev.eval(&sum)
            })
            .collect::<Result<Vec<_>, _>>()?;
        rows.push(row);
    }
    if piecewise {
        let bound = |i: usize| Expr::from(knots[i].clone());
        let mut columns_out = Vec::with_capacity(columns);
        for c in 0..columns {
            let mut args = Vec::with_capacity(2 * pieces);
            for (i, row) in rows.iter().enumerate() {
                let condition = if i == 0 {
                    Expr::call(
                        Symbol::AND,
                        [
                            Expr::call(Symbol::GREATER_EQUAL, [x.clone(), bound(0)]),
                            Expr::call(Symbol::LESS, [x.clone(), bound(1)]),
                        ],
                    )
                } else if i == pieces - 1 {
                    Expr::call(Symbol::LESS_EQUAL, [x.clone(), bound(i + 1)])
                } else {
                    Expr::call(Symbol::LESS, [x.clone(), bound(i + 1)])
                };
                args.push(condition);
                args.push(row[c].clone());
            }
            columns_out.push(ev.eval(&Expr::call(Symbol::PIECEWISE, args))?);
        }
        return Ok(Some(Expr::call(Symbol::LIST, columns_out)));
    }
    let result = if columns == 1 {
        let mut column: Vec<Expr> = rows.into_iter().map(|mut r| r.swap_remove(0)).collect();
        if column.len() == 1 {
            column.swap_remove(0)
        } else {
            Expr::call(Symbol::LIST, column)
        }
    } else {
        Expr::call(
            Symbol::LIST,
            rows.into_iter().map(|r| Expr::call(Symbol::LIST, r)),
        )
    };
    Ok(Some(result))
}

/// Whether `I` or a number not real occurs in `e`.
fn has_complex(e: &Expr) -> bool {
    matches!(e.view(), View::Complex(_) | View::ComplexBall(_))
        || e.as_symbol() == Some(Symbol::I)
        || e.args().iter().any(has_complex)
}

/// A control point's coordinates when it is not a row of a matrix.
fn coordinates(z: &Expr) -> Option<Vec<Expr>> {
    match z.view() {
        View::Complex(z) => Some(vec![Expr::from(z.re()), Expr::from(z.im())]),
        View::ComplexBall(z) => Some(vec![Expr::from(z.re()), Expr::from(z.im())]),
        _ if is_list(z) || has_complex(z) => None,
        _ => Some(vec![z.clone(), Expr::from(0)]),
    }
}

/// The rows of coordinates of a list of control points; `None` for
/// another form or past the caps.
fn control_points(c: &Expr) -> Option<Vec<Vec<Expr>>> {
    if !is_list(c) || c.args().is_empty() || c.args().len() > MAX_POINTS {
        return None;
    }
    let rows: Vec<Vec<Expr>> = match matrix_rows(c) {
        Some(rows) => rows.iter().map(|r| r.args().to_vec()).collect(),
        None => c.args().iter().map(coordinates).collect::<Option<_>>()?,
    };
    let symbolic = rows.iter().flatten().any(|v| !v.is_number());
    (rows[0].len() <= MAX_COLUMNS && !(symbolic && rows.len() > MAX_SYMBOLIC_CONTROL_POINTS))
        .then_some(rows)
}

/// The variable and the `n` knots of `x`, `x == a..b` or `x == [t…]`.
fn variable_and_knots(var: &Expr, n: usize) -> Option<(Expr, Vec<Rational>)> {
    if var.as_symbol().is_some() {
        return Some((
            var.clone(),
            linspace(&Rational::from(0), &Rational::from(1), n),
        ));
    }
    let (Some(Symbol::EQUATION | Symbol::EQUAL), [x, rh]) =
        (var.head().and_then(Expr::as_symbol), var.args())
    else {
        return None;
    };
    x.as_symbol()?;
    let knots = if rh.head().and_then(Expr::as_symbol) == Some(Symbol::RANGE) {
        let [a, b] = rh.args() else {
            return None;
        };
        let (a, b) = (a.to_rational()?, b.to_rational()?);
        if b <= a {
            return None;
        }
        linspace(&a, &b, n)
    } else if is_list(rh) && rh.args().len() == n {
        let knots = rh
            .args()
            .iter()
            .map(Expr::to_rational)
            .collect::<Option<Vec<_>>>()?;
        if knots.windows(2).any(|w| w[1] <= w[0]) {
            return None;
        }
        knots
    } else {
        return None;
    };
    Some((x.clone(), knots))
}

/// A degree option: an integer from 1 to [`MAX_DEGREE`].
fn degree(p: &Expr) -> Option<usize> {
    p.as_i64()
        .and_then(|p| usize::try_from(p).ok())
        .filter(|p| (1..=MAX_DEGREE).contains(p))
}

pub(in crate::native) fn b_spline(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let mut args = e.args();
    let piecewise = args.last().and_then(Expr::as_symbol) == Some(Symbol::PIECEWISE);
    if piecewise {
        args = &args[..args.len() - 1];
    }
    let default = Expr::symbol("x");
    let (c, var, p) = match args {
        // giac's default variable, with one argument in all.
        [c] if !piecewise => (c, &default, 3),
        [c, var] => (c, var, 3),
        [c, var, p] => match degree(p) {
            Some(p) => (c, var, p),
            None => return Ok(None),
        },
        _ => return Ok(None),
    };
    let Some(points) = control_points(c) else {
        return Ok(None);
    };
    if points.len() < p + 1 {
        return Ok(None);
    }
    let Some((x, knots)) = variable_and_knots(var, points.len() - p + 1) else {
        return Ok(None);
    };
    let Some(result) = spline(ev, &x, &knots, &points, p, piecewise)? else {
        return Ok(None);
    };
    ev.step("BSpline", e, &result);
    Ok(Some(result))
}

/// giac's `subinterval_index`: the `k` with `b_k ≤ x < b_{k+1}`, `x`
/// clamped to the knots, the last nonempty interval at the right end.
fn subinterval_index(b: &[Rational], x: &Rational) -> Option<usize> {
    let x = x.clamp(&b[0], &b[b.len() - 1]);
    let mut i = b.partition_point(|v| v <= x);
    if i == 0 {
        return None;
    }
    if i == b.len() {
        i -= 1;
    }
    while i > 1 && b[i - 1] == b[i] {
        i -= 1;
    }
    Some(i - 1)
}

/// The `n + 1` basis splines of degree `p` on the padded knots `t` at `u`:
/// giac's `eval_basis_splines`.
pub(super) fn basis(n: usize, p: usize, u: &Rational, t: &[Rational]) -> Option<Vec<Rational>> {
    let zero = Rational::from(0);
    let mut b = vec![zero; n + 1];
    let (lo, hi) = (&t[0], &t[t.len() - 1]);
    if u < lo || u > hi {
        return Some(b);
    }
    if u == lo {
        b[0] = Rational::from(1);
        return Some(b);
    }
    if u == hi {
        b[n] = Rational::from(1);
        return Some(b);
    }
    let k = subinterval_index(t, u).filter(|&k| k >= p && k <= n && k + p < t.len())?;
    b[k] = Rational::from(1);
    for d in 1..=p {
        b[k - d] = &(&(&t[k + 1] - u) / &(&t[k + 1] - &t[k - d + 1])) * &b[k - d + 1];
        for i in k - d + 1..k {
            let left = &(&(u - &t[i]) / &(&t[i + d] - &t[i])) * &b[i];
            let right = &(&(&t[i + d + 1] - u) / &(&t[i + d + 1] - &t[i + 1])) * &b[i + 1];
            b[i] = &left + &right;
        }
        b[k] = &(&(u - &t[k]) / &(&t[k + d] - &t[k])) * &b[k];
    }
    Some(b)
}

/// The solution `X` of `A X = B` over the rationals (`B` by rows), by
/// Gauss–Jordan elimination; `None` when `A` is singular.
pub(super) fn solve(
    mut a: Vec<Vec<Rational>>,
    mut b: Vec<Vec<Rational>>,
) -> Option<Vec<Vec<Rational>>> {
    let n = a.len();
    for c in 0..n {
        let pivot = (c..n).find(|&r| !a[r][c].is_zero())?;
        a.swap(c, pivot);
        b.swap(c, pivot);
        let inverse = a[c][c].recip();
        for v in a[c].iter_mut().chain(b[c].iter_mut()) {
            *v = &*v * &inverse;
        }
        for r in 0..n {
            if r == c || a[r][c].is_zero() {
                continue;
            }
            let f = a[r][c].clone();
            let (pivot_a, pivot_b) = (a[c].clone(), b[c].clone());
            for (v, w) in a[r].iter_mut().zip(&pivot_a) {
                *v = &*v - &(&f * w);
            }
            for (v, w) in b[r].iter_mut().zip(&pivot_b) {
                *v = &*v - &(&f * w);
            }
        }
    }
    Some(b)
}

/// The control points of the least-squares fit of the rows `ys` at `s` by
/// the `n + 1` basis splines of degree `p` on the padded knots `t`:
/// `(AᵀA) Q = AᵀY`; `None` when singular.
pub(super) fn least_squares(
    s: &[Rational],
    ys: &[Vec<Rational>],
    t: &[Rational],
    n: usize,
    p: usize,
) -> Option<Vec<Vec<Rational>>> {
    let zero = Rational::from(0);
    let a = s
        .iter()
        .map(|u| basis(n, p, u, t))
        .collect::<Option<Vec<_>>>()?;
    let series = ys[0].len();
    let mut ata = vec![vec![zero.clone(); n + 1]; n + 1];
    let mut aty = vec![vec![zero; series]; n + 1];
    for (row, y) in a.iter().zip(ys) {
        for i in 0..=n {
            if row[i].is_zero() {
                continue;
            }
            for j in 0..=n {
                ata[i][j] = &ata[i][j] + &(&row[i] * &row[j]);
            }
            for (k, v) in y.iter().enumerate() {
                aty[i][k] = &aty[i][k] + &(&row[i] * v);
            }
        }
    }
    solve(ata, aty)
}

/// giac's `deBoor` at a number: the spline's value at `x` on the piece `k`.
fn de_boor(k: usize, x: &Rational, t: &[Rational], q: &[Vec<Rational>], p: usize) -> Vec<Rational> {
    let one = Rational::from(1);
    let mut d: Vec<Vec<Rational>> = (0..=p).map(|j| q[j + k - p].clone()).collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let a = &(x - &t[j + k - p]) / &(&t[j + 1 + k - r] - &t[j + k - p]);
            let b = &one - &a;
            d[j] = d[j]
                .iter()
                .zip(&d[j - 1])
                .map(|(u, v)| &(&a * u) + &(&b * v))
                .collect();
        }
    }
    d.swap_remove(p)
}

/// The options of a fit.
struct Options {
    n: i64,
    p: usize,
    breakpoints: Vec<Rational>,
    derive: usize,
    piecewise: bool,
    decimal: bool,
}

/// The options after the data and the variable; `None` for another.
fn options(given: &[Expr], m: usize) -> Option<Options> {
    // giac's default: n = ⌊√m + 1/2⌋, (2n − 1)² ≤ 4m < (2n + 1)².
    let mut n = 0;
    while (2 * n + 1) * (2 * n + 1) <= 4 * m {
        n += 1;
    }
    let mut o = Options {
        n: i64::try_from(n).expect("a size in memory"),
        p: 3,
        breakpoints: Vec::new(),
        derive: 0,
        piecewise: false,
        decimal: false,
    };
    let named = |e: &Expr, name: &str| *e == Expr::symbol(name);
    for option in given {
        if let Some(k) = option.as_i64() {
            o.n = k;
            continue;
        }
        match option.as_symbol() {
            Some(Symbol::PIECEWISE) => o.piecewise = true,
            Some(Symbol::UNIFORM_PDF) => {}
            Some(_) if named(option, "close") => {}
            _ => {
                let (Some(Symbol::EQUATION | Symbol::EQUAL), [lh, rh]) =
                    (option.head().and_then(Expr::as_symbol), option.args())
                else {
                    return None;
                };
                if lh.as_symbol() == Some(Symbol::EXPONENT) {
                    o.p = degree(rh)?;
                } else if named(lh, "breakpoint") {
                    if !is_list(rh) || rh.args().is_empty() || rh.args().len() > MAX_POINTS {
                        return None;
                    }
                    o.breakpoints.clear();
                    for v in rh.args() {
                        let (q, decimal) = real(v)?;
                        o.decimal |= decimal;
                        o.breakpoints.push(q);
                    }
                    o.n = i64::try_from(o.breakpoints.len()).expect("a size in memory") - 1;
                } else if named(lh, "derive") {
                    o.derive = usize::try_from(rh.as_i64()?).ok()?;
                } else {
                    return None;
                }
            }
        }
    }
    Some(o)
}

/// The abscissas and the rows of values of a fit's data, exactly, and
/// whether a decimal was among them; `None` for another form or size.
fn data(data: &Expr) -> Option<(Vec<Rational>, Vec<Vec<Rational>>, bool)> {
    let rows = matrix_rows(data)?;
    let width = rows[0].args().len();
    if !(2..=MAX_POINTS).contains(&rows.len()) || !(2..=MAX_COLUMNS + 1).contains(&width) {
        return None;
    }
    let mut decimal = false;
    let mut s = Vec::with_capacity(rows.len());
    let mut ys = Vec::with_capacity(rows.len());
    for row in rows {
        let mut values = Vec::with_capacity(width);
        for v in row.args() {
            let (q, d) = real(v)?;
            decimal |= d;
            values.push(q);
        }
        s.push(values.remove(0));
        ys.push(values);
    }
    Some((s, ys, decimal))
}

/// The values at `points` of the spline of degree `p` with control points
/// `q` on the padded knots `t`, and whether a point was a decimal; `None`
/// for a point not a number or outside `[lo, hi]`.
fn values_at(
    points: &[Expr],
    (lo, hi): (&Rational, &Rational),
    t: &[Rational],
    q: &[Vec<Rational>],
    p: usize,
) -> Option<(Vec<Vec<Rational>>, bool)> {
    let mut decimal = false;
    let mut values = Vec::with_capacity(points.len());
    for x in points {
        let (x, d) = real(x)?;
        if x < *lo || x > *hi {
            return None;
        }
        let k = subinterval_index(t, &x).filter(|&k| k >= p && k < q.len())?;
        decimal |= d;
        values.push(de_boor(k, &x, t, q, p));
    }
    Some((values, decimal))
}

/// `v` exactly, or as a ball at `prec` bits when `decimal`.
pub(super) fn number(v: Rational, decimal: bool, prec: u64) -> Expr {
    if decimal {
        Expr::from(Ball::from_rational(&v, prec))
    } else {
        Expr::from(v)
    }
}

pub(in crate::native) fn spline_fit(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let default = Expr::symbol("x");
    let (data, var, given) = match e.args() {
        // giac's default variable, with one argument in all.
        [data] => (data, &default, &[][..]),
        [data, var, given @ ..] => (data, var, given),
        [] => return Ok(None),
    };
    let Some((s, ys, mut decimal)) = self::data(data) else {
        return Ok(None);
    };
    let m = s.len();
    let series = ys[0].len();
    let Some(o) = options(given, m) else {
        return Ok(None);
    };
    decimal |= o.decimal;
    let mut p = o.p;
    // giac's errors: n < 1, derive ≥ p, n + p − 1 ≥ m.
    let Some(n) = usize::try_from(o.n)
        .ok()
        .filter(|&n| n >= 1 && n + p - 1 < m)
    else {
        return Ok(None);
    };
    if o.derive >= p {
        return Ok(None);
    }
    let lo = s.iter().min().expect("two rows").clone();
    let hi = s.iter().max().expect("two rows").clone();
    let knots = if o.breakpoints.is_empty() {
        if lo == hi {
            return Ok(None);
        }
        linspace(&lo, &hi, n + 1)
    } else {
        let mut b = o.breakpoints;
        b.sort();
        if b.len() != n + 1 || b.windows(2).any(|w| w[1] <= w[0]) || b[0] > lo || hi > b[n] {
            return Ok(None);
        }
        b
    };
    let mut ncp = n + p - 1;
    let size = ncp + 1;
    if m * size * size > MAX_FIT_WORK || size * size * series > MAX_FIT_WORK {
        return Ok(None);
    }
    let mut t = padded(&knots, p);
    let Some(mut q) = least_squares(&s, &ys, &t, ncp, p) else {
        return Ok(None);
    };
    for _ in 0..o.derive {
        for i in 0..ncp {
            let Some(scale) = count(p).checked_div(&(&t[i + p + 1] - &t[i + 1])) else {
                return Ok(None);
            };
            q[i] = q[i + 1]
                .iter()
                .zip(&q[i])
                .map(|(a, b)| &scale * &(a - b))
                .collect();
        }
        q.pop();
        p -= 1;
        ncp -= 1;
        t = t[1..t.len() - 1].to_vec();
    }
    let prec = ev.context().working_precision();
    let result = if is_list(var) {
        let points = var.args();
        if points.len() > MAX_POINTS || points.len() * (p + 1) * (p + 1) * series > MAX_FIT_WORK {
            return Ok(None);
        }
        let Some((values, d)) = values_at(points, (&lo, &hi), &t, &q, p) else {
            return Ok(None);
        };
        decimal |= d;
        let number = |v| number(v, decimal, prec);
        let rows = values.into_iter().map(|v| {
            if series == 1 {
                number(v.into_iter().next().expect("one series"))
            } else {
                Expr::call(Symbol::LIST, v.into_iter().map(number))
            }
        });
        Expr::call(Symbol::LIST, rows)
    } else if var.as_symbol().is_some() {
        let points: Vec<Vec<Expr>> = q
            .into_iter()
            .map(|row| row.into_iter().map(|v| number(v, decimal, prec)).collect())
            .collect();
        let Some(result) = spline(ev, var, &knots, &points, p, o.piecewise)? else {
            return Ok(None);
        };
        result
    } else {
        return Ok(None);
    };
    let result = ev.eval(&result)?;
    ev.step("SplineFit", e, &result);
    Ok(Some(result))
}
