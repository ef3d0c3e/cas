//! Clustering (CAS-630): Xcas's `kmeans` and `cluster` (giac
//! optimization.cc: AS 136 `optra` :9862, `qtran` :10096, `kmns` :10292;
//! `init_cluster_centers` :10622, `kmeans` :10654; `comp_dist` :10765;
//! `cluster_crit` :10784-11091; `hclust` :11102-11443; `clustering_result` :11571, `_kmeans` :11610,
//! `_cluster` :11737), ported from crates/cas-core/proto/clustering/proto.py (R1.7-012).
//!
//! - **Data.** A list of real numbers (1-D points) or of equal-length
//!   lists of real numbers, `Matrix[l]` of such a list too; reals are
//!   rationals, decimals (their midpoints) and numeric constants, in `f64`;
//!   for `cluster` also a list of strings. Anything else leaves the call
//!   without a warning.
//! - **`KMeans[data, options…]`** (`kmeans`): Hartigan and Wong's k-means
//!   (AS 136, Burkardt's C++ as giac has it: two clusters stop after one
//!   `optra` and `qtran` pass), from centres drawn by giac's
//!   `init_cluster_centers` through the session's generator (`SeedRandom`
//!   repeats a run), up to 11 starts while a cluster comes out empty. `k`
//!   positional or `count=k`; `count=a..b`, `count_inf=`, `count_sup=` or
//!   none choose `k` by Hartigan's rule (from the least `k` up, the previous
//!   `k` once `(W_prev/W − 1)(m − k) ≤ 10`); a range of one `k` fixes it,
//!   over a given one, as giac's. `limit=` the iterations (100). `k ≤ 1` or
//!   `k ≥ m` answers `[]` with no draw (giac's fault 3); no convergence or
//!   no start answers `[]` with giac's message as the warning.
//! - **`kmeans`' `index=`** (a name, a list of them OR'd, `true` for
//!   `"calinski-harabasz"`, `false` for Hartigan's rule; the last of `k`,
//!   `count=k` and `index=` wins, as giac's): giac's `cluster_crit`, each `k`
//!   from the least up, its first start from the previous `k`'s classes,
//!   Hartigan's rough optimum capping the largest `k` at `2k`. The fourteen
//!   indices, in giac's order: Banfeld–Raftery, Davies–Bouldin, Ray–Turi,
//!   Scott–Symons (least), Calinski–Harabasz, PBM, Ratkowsky–Lance
//!   (greatest), the determinant, log determinant and log sum-of-squares
//!   ratios (least second difference), Ball–Hall, `k² det W`, `trace W`,
//!   `trace W⁻¹B` (greatest second difference); `"log-det"` names the
//!   determinant ratio, as giac's, the log one only comes with `"all"`. As
//!   giac's: the total scatter exact on exact data, a double elsewhere; the
//!   within scatter's lower triangle copied from the total scatter last
//!   computed (the identity's zeros before), the determinant `in_doublerref`'s
//!   (pivots under 1e-12 relative are zero), exact over a determinant of
//!   exactly 1.0; a division by zero `ComplexInfinity` (giac's `infinity`),
//!   `trace W⁻¹B` of a singular block of `W` too, or undefined where `B`
//!   vanishes on it. The labels of the `k` most indices suggest (the least
//!   of a tie; none: `[]`, count 0).
//! - **`Cluster[data, options…]`** (`cluster`): agglomerative clustering,
//!   single linkage by Prim's tree (the default), `type="complete"`,
//!   `"average"`, `"weighted"` and `"ward"` by the nearest-neighbour chain
//!   with Lance and Williams' updates, on squared Euclidean distances
//!   (exact on rational points, then `f64`) or `distance=f` (a builtin or a
//!   `Function`, also a builtin positionally), `N[f[p, q]]` per pair, cached
//!   as giac's; strings by giac's Levenshtein distance on their bytes
//!   (`levenshtein`, :11504-11536). The merges sorted by height; `k` given,
//!   or giac's largest gap of heights (at least `count_inf`, the least of
//!   `count=a..b`). `index=` (a name, a list of them, or `true` for the
//!   silhouette) without `k`: from `max_k` clusters down to that gap, giac's
//!   silhouette on the same distances, each cluster's mean of its points'
//!   values and the largest over clusters (`DBL_MIN` at least), from its
//!   running sums; the labels of its first best `k`. `"mcclain-rao"` and
//!   `"dunn"` compute nothing, as giac's. Departure: giac asserts and aborts
//!   when `max_k ≥ N` (the default 16 included); ours starts at `N − 1`.
//! - **Answers.** The 0-based labels, numbered by first appearance for
//!   `cluster` (`output=list`, the default), the number of clusters
//!   (`output=count`), exact integers, or the clusters by label, each the
//!   data items in order (`output=part`). `cluster`'s `output=index`: the
//!   rows `[k, value]` of the silhouette; `kmeans`': those of its first
//!   index, "not defined" where undefined (`[]` when `count=a..a` fixes
//!   `k`). `output=plot`, `tree`, `display=`, `color=` (and `labels=` for
//!   `cluster`) leave the call without a warning.
//! - **Errors** leave the call with giac's message as the warning: "No
//!   data found", kmeans' "Data must be numeric", "Invalid k-range", "Invalid
//!   output specification", "No index specified for output", "Unknown
//!   linkage type", "Expected a string", "Invalid distance function
//!   specification", "Invalid option", "Unknown index name", "Expected a
//!   string or list of strings", "Unable to transpose" (`output=index` with
//!   no silhouette values), the count parsers' texts, and "Bad distance
//!   function" for a distance not a real `≥ 0` (a symbolic point among
//!   them).
//! - **Limits** before any work (`Limits::terms`): `cluster` the pairs
//!   `m(m − 1)/2` and `m n`, with the silhouette the pairs times its steps,
//!   on strings the pairs times the longest's bytes and that length squared;
//!   `kmeans` `11 m n Σk` over the `k` tried and `m n
//!   k limit`, with an index `(m n² + k² n + n³)` per `k`; past them the call stays with "too many points". `qtran`
//!   stops after `terms/n` steps (giac's has no bound; a fault 2); the loops
//!   poll the interrupt.

#![allow(
    clippy::many_single_char_names,
    reason = "a, c, d, k, m, n: AS 136's and giac's names"
)]

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering as Atomic};

use cas_num::ball::{Ball, Float, Round};
use cas_num::{Integer, Natural, Rational};

use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::native::numeric::approximation;
use crate::native::statistics::random_variables::{equation, range};
use crate::symbol::Symbol;
use crate::trace::Warning;

#[cfg(test)]
mod tests;

type Native = Result<Option<Expr>, EvalError>;

/// giac's texts.
const NO_DATA: &str = "No data found";
const NUMERIC: &str = "Data must be numeric";
const K_RANGE: &str = "Invalid k-range";
const OUTPUT: &str = "Invalid output specification";
const NO_INDEX: &str = "No index specified for output";
const LINKAGE: &str = "Unknown linkage type";
const STRING: &str = "Expected a string";
const DISTANCE: &str = "Invalid distance function specification";
const BAD_DISTANCE: &str = "Bad distance function";
const OPTION: &str = "Invalid option";
const UNKNOWN_INDEX: &str = "Unknown index name";
const STRINGS: &str = "Expected a string or list of strings";
/// `mtran([kvals, values])` of unequal or empty lists (optimization.h:640).
const TRANSPOSE: &str = "Unable to transpose";
/// An undefined index value in `kmeans`' `output=index` (:10741).
const NOT_DEFINED: &str = "not defined";
const POSITIVE: &str = "Expected a positive integer";
const RANGE: &str = "Expected a range of positive integers";
const COUNT: &str = "Expected an integer or range";
const CENTERS: &str = "failed to generate good initial cluster centers";
const ITERATIONS: &str = "maximum number of iterations exceeded";
/// Past `Limits::terms`.
const TOO_MANY: &str = "too many points";

/// AS 136's `r8_huge` (:9822).
const HUGE: f64 = 1.0e30;
/// giac's starts of one `k`: the first and 10 more (:10675-10678).
const STARTS: usize = 11;

/// `hclust`'s index bits (optimization.h:83-86); only the silhouette is
/// computed (:11262-11443).
const SILHOUETTE: u32 = 1;

/// The call stays, with a warning.
fn refuse(ev: &mut Eval<'_>, e: &Expr, why: &'static str) -> Native {
    ev.warn(Warning::Native {
        message: why.into(),
        expr: e.clone(),
    });
    Ok(None)
}

fn is_list(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::LIST))
}

fn list(items: Vec<Expr>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

/// Whether `e` is the symbol `s` or a symbol named `name` (C0008).
fn is(e: &Expr, s: Symbol, name: &str) -> bool {
    e.as_symbol().is_some_and(|t| t == s || t.name() == name)
}

/// Whether `e` is a symbol named `name`.
fn named(e: &Expr, name: &str) -> bool {
    e.as_symbol().is_some_and(|t| t.name() == name)
}

fn interrupted(flag: &AtomicBool) -> Result<(), EvalError> {
    if flag.load(Atomic::Relaxed) {
        Err(EvalError::Interrupted)
    } else {
        Ok(())
    }
}

fn to_f64(q: &Rational) -> f64 {
    Float::from_rational(q, 53, Round::Nearest).0.to_f64()
}

/// A coordinate: exact when a rational, and its double.
#[derive(Clone, Debug)]
struct Coordinate {
    exact: Option<Rational>,
    value: f64,
}

/// A finite real number: a rational, a decimal's midpoint, a numeric
/// constant's approximation.
fn coordinate(ev: &Eval<'_>, x: &Expr) -> Option<Coordinate> {
    let (exact, value) = match x.view() {
        View::Ball(b) if b.is_finite() => (None, b.mid().to_f64()),
        View::Ball(_) => return None,
        _ => {
            if let Some(q) = x.to_rational() {
                let v = to_f64(&q);
                (Some(q), v)
            } else if x.is_numeric() {
                match approximation(ev.context(), x, 53).ok()?.view() {
                    View::Ball(b) if b.is_finite() => (None, b.mid().to_f64()),
                    _ => return None,
                }
            } else {
                return None;
            }
        }
    };
    value.is_finite().then_some(Coordinate { exact, value })
}

/// The points of `data`: a list of numbers or of equal-length lists, or
/// `Matrix[l]` of one; `None` for another shape. A point is `None` when a
/// coordinate is no real number.
struct Data<'a> {
    items: &'a [Expr],
    points: Vec<Option<Vec<Coordinate>>>,
    /// The dimension.
    n: usize,
}

fn data<'a>(ev: &Eval<'_>, d: &'a Expr) -> Option<Data<'a>> {
    let d = match (d.head().and_then(Expr::as_symbol), d.args()) {
        (Some(Symbol::MATRIX), [l]) if is_list(l) => l,
        _ => d,
    };
    if !is_list(d) {
        return None;
    }
    let items = d.args();
    let lists = items.iter().filter(|x| is_list(x)).count();
    let n = match items.first() {
        None => 0,
        Some(_) if lists == 0 => 1,
        Some(first) if lists == items.len() => {
            let n = first.args().len();
            if n == 0 || items.iter().any(|x| x.args().len() != n) {
                return None;
            }
            n
        }
        Some(_) => return None,
    };
    let points = items
        .iter()
        .map(|x| {
            if is_list(x) {
                x.args().iter().map(|c| coordinate(ev, c)).collect()
            } else {
                coordinate(ev, x).map(|c| vec![c])
            }
        })
        .collect();
    Some(Data { items, points, n })
}

/// An integer saturated to `i64`.
fn integer(e: &Expr) -> Option<i64> {
    let k = e.to_integer()?;
    Some(
        k.to_i64()
            .unwrap_or(if k.is_positive() { i64::MAX } else { i64::MIN }),
    )
}

/// A `count=` value (:11623-11638, :11750-11765).
enum Count {
    K(i64),
    Range(i64, i64),
}

fn count(rh: &Expr) -> Result<Count, &'static str> {
    if let Some(k) = integer(rh) {
        return if k <= 0 {
            Err(POSITIVE)
        } else {
            Ok(Count::K(k))
        };
    }
    if let Some((a, b)) = range(rh) {
        return match (integer(a), integer(b)) {
            (Some(a), Some(b)) if a > 0 && b >= a => Ok(Count::Range(a, b)),
            _ => Err(RANGE),
        };
    }
    Err(COUNT)
}

/// The answer asked for.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Output {
    List,
    Count,
    Index,
    /// The clusters by label, each the data items in order.
    Part,
    /// Not ported: the call stays without a warning.
    Out,
}

/// `output=` (:11694-11699, :11810-11817); `tree` only for `cluster`.
fn output(rh: &Expr, tree: bool) -> Result<Output, &'static str> {
    if is(rh, Symbol::LIST, "list") {
        Ok(Output::List)
    } else if is(rh, Symbol::COUNT, "count") {
        Ok(Output::Count)
    } else if is(rh, Symbol::FIRST_POSITION, "index") {
        Ok(Output::Index)
    } else if is(rh, Symbol::PART, "part") {
        Ok(Output::Part)
    } else if named(rh, "plot") || (tree && named(rh, "tree")) {
        Ok(Output::Out)
    } else {
        Err(OUTPUT)
    }
}

/// Options that only draw: the call stays without a warning.
fn drawing(lh: &Expr) -> bool {
    named(lh, "display") || named(lh, "color")
}

/// `count_sup=k` and `count_inf=k`: any integer.
fn bound(rh: &Expr) -> Result<i64, &'static str> {
    integer(rh).ok_or(POSITIVE)
}

/// `index=` (:11838-11853): a boolean (`true` is `on`), a name, or a
/// non-empty list of names OR'd; `names` gives a name's bits.
fn index(rh: &Expr, names: fn(&str) -> Option<u32>, on: u32) -> Result<u32, &'static str> {
    if rh.as_symbol() == Some(Symbol::TRUE) {
        return Ok(on);
    }
    if rh.as_symbol() == Some(Symbol::FALSE) {
        return Ok(0);
    }
    if let Some(name) = rh.as_str() {
        return names(name).ok_or(UNKNOWN_INDEX);
    }
    if is_list(rh) && !rh.args().is_empty() && rh.args().iter().all(|x| x.as_str().is_some()) {
        return rh.args().iter().try_fold(0, |bits, x| {
            x.as_str()
                .and_then(names)
                .map(|b| bits | b)
                .ok_or(UNKNOWN_INDEX)
        });
    }
    Err(STRINGS)
}

/// `hclust::name2index` (:11481-11487).
fn hclust_index(name: &str) -> Option<u32> {
    Some(match name {
        "silhouette" => SILHOUETTE,
        "mcclain-rao" => 2,
        "dunn" => 4,
        "all" => 7,
        _ => return None,
    })
}

/// `clustering_result` (:11571-11603): the labels, their count, or the
/// clusters by label (`output=part`) of the data `items`.
fn answer(out: Output, k: usize, labels: &[usize], items: &[Expr]) -> Expr {
    match out {
        Output::Count => Expr::from(Integer::from(k as u64)),
        Output::Part => list(
            (0..k)
                .map(|c| {
                    list(
                        labels
                            .iter()
                            .zip(items)
                            .filter(|(l, _)| **l == c)
                            .map(|(_, x)| x.clone())
                            .collect(),
                    )
                })
                .collect(),
        ),
        _ => list(
            labels
                .iter()
                .map(|&l| Expr::from(Integer::from(l as u64)))
                .collect(),
        ),
    }
}

// ---------------- AS 136 (optimization.cc:9822-10609) ----------------

/// AS 136's arrays, 1-based as in the source.
struct Kmns<'a> {
    /// Points, `a[0]` unused.
    a: &'a [Vec<f64>],
    m: usize,
    k: usize,
    c: Vec<Vec<f64>>,
    ic1: Vec<usize>,
    ic2: Vec<usize>,
    nc: Vec<usize>,
    an1: Vec<f64>,
    an2: Vec<f64>,
    ncp: Vec<i64>,
    d: Vec<f64>,
    itran: Vec<u8>,
    live: Vec<i64>,
    indx: usize,
}

fn squared(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).fold(0.0, |s, (a, b)| {
        let d = a - b;
        s + d * d
    })
}

#[allow(clippy::cast_possible_wrap, reason = "point indices, far below 2⁶³")]
fn wide(i: usize) -> i64 {
    i as i64
}

impl Kmns<'_> {
    /// Moves point `i` from `l1` to `l2`, updating the centres.
    #[allow(clippy::cast_precision_loss, reason = "cluster sizes")]
    fn transfer(&mut self, i: usize, l1: usize, l2: usize) {
        let al1 = self.nc[l1] as f64;
        let alw = al1 - 1.0;
        let al2 = self.nc[l2] as f64;
        let alt = al2 + 1.0;
        for (j, &x) in self.a[i].iter().enumerate() {
            self.c[l1][j] = (self.c[l1][j] * al1 - x) / alw;
            self.c[l2][j] = (self.c[l2][j] * al2 + x) / alt;
        }
        self.nc[l1] -= 1;
        self.nc[l2] += 1;
        self.an2[l1] = alw / al1;
        self.an1[l1] = if 1.0 < alw { alw / (alw - 1.0) } else { HUGE };
        self.an1[l2] = alt / al2;
        self.an2[l2] = alt / (alt + 1.0);
        self.ic1[i] = l2;
        self.ic2[i] = l1;
    }

    /// The optimal-transfer stage (:9862-10094).
    fn optra(&mut self) {
        let (m, k) = (self.m, self.k);
        for l in 1..=k {
            if self.itran[l] == 1 {
                self.live[l] = wide(m) + 1;
            }
        }
        for i in 1..=m {
            self.indx += 1;
            let l1 = self.ic1[i];
            let mut l2 = self.ic2[i];
            let ll = l2;
            if 1 < self.nc[l1] {
                if self.ncp[l1] != 0 {
                    self.d[i] = squared(&self.a[i], &self.c[l1]) * self.an1[l1];
                }
                let da = squared(&self.a[i], &self.c[l2]);
                let mut r2 = da * self.an2[l2];
                let wi = wide(i);
                for l in 1..=k {
                    if (wi < self.live[l1] || wi < self.live[l2]) && l != l1 && l != ll {
                        let rr = r2 / self.an2[l];
                        let dc = squared(&self.a[i], &self.c[l]);
                        if dc < rr {
                            r2 = dc * self.an2[l];
                            l2 = l;
                        }
                    }
                }
                if self.d[i] <= r2 {
                    self.ic2[i] = l2;
                } else {
                    self.indx = 0;
                    self.live[l1] = wide(m + i);
                    self.live[l2] = wide(m + i);
                    self.ncp[l1] = wi;
                    self.ncp[l2] = wi;
                    self.transfer(i, l1, l2);
                }
            }
            if self.indx == m {
                return;
            }
        }
        for l in 1..=k {
            self.itran[l] = 0;
            self.live[l] -= wide(m);
        }
    }

    /// The quick-transfer stage (:10096-10290); `false` past `steps`.
    fn qtran(&mut self, steps: i64, flag: &AtomicBool) -> Result<bool, EvalError> {
        let m = self.m;
        let mut icoun = 0;
        let mut istep: i64 = 0;
        loop {
            interrupted(flag)?;
            if istep > steps {
                return Ok(false);
            }
            for i in 1..=m {
                icoun += 1;
                istep += 1;
                let l1 = self.ic1[i];
                let l2 = self.ic2[i];
                if 1 < self.nc[l1] {
                    if istep <= self.ncp[l1] {
                        let da = squared(&self.a[i], &self.c[l1]);
                        self.d[i] = da * self.an1[l1];
                    }
                    if istep < self.ncp[l1] || istep < self.ncp[l2] {
                        let r2 = self.d[i] / self.an2[l2];
                        let dd = squared(&self.a[i], &self.c[l2]);
                        if dd < r2 {
                            icoun = 0;
                            self.indx = 0;
                            self.itran[l1] = 1;
                            self.itran[l2] = 1;
                            self.ncp[l1] = istep + wide(m);
                            self.ncp[l2] = istep + wide(m);
                            self.transfer(i, l1, l2);
                        }
                    }
                }
                if icoun == m {
                    return Ok(true);
                }
            }
        }
    }
}

/// AS 136's faults.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Fault {
    /// An empty cluster from the starting centres.
    Centers,
    /// No convergence within the iterations.
    Iterations,
    /// `k ≤ 1` or `m ≤ k`.
    Size,
}

/// A run's labels `1 … k`, each cluster's sum of squares and its centre, or
/// its fault.
type Clusters = Result<(Vec<usize>, Vec<f64>, Vec<Vec<f64>>), Fault>;

/// `kmns` (:10292-10589) on `points` from the centres `c0`: the labels
/// `1 … k` of the points, each cluster's sum of squares and its centre.
#[allow(clippy::cast_precision_loss, reason = "cluster sizes")]
#[allow(
    clippy::needless_range_loop,
    clippy::too_many_lines,
    reason = "AS 136's 1-based loops, ported literally"
)]
fn kmns(
    points: &[Vec<f64>],
    c0: &[Vec<f64>],
    k: usize,
    iterations: usize,
    steps: i64,
    flag: &AtomicBool,
) -> Result<Clusters, EvalError> {
    let m = points.len();
    if k <= 1 || m <= k {
        return Ok(Err(Fault::Size));
    }
    let n = points[0].len();
    let mut a = Vec::with_capacity(m + 1);
    a.push(Vec::new());
    a.extend(points.iter().cloned());
    let mut c = Vec::with_capacity(k + 1);
    c.push(Vec::new());
    c.extend(c0.iter().cloned());
    let mut s = Kmns {
        a: &a,
        m,
        k,
        c,
        ic1: vec![0; m + 1],
        ic2: vec![0; m + 1],
        nc: vec![0; k + 1],
        an1: vec![0.0; k + 1],
        an2: vec![0.0; k + 1],
        ncp: vec![0; k + 1],
        d: vec![0.0; m + 1],
        itran: vec![0; k + 1],
        live: vec![0; k + 1],
        indx: 0,
    };
    for i in 1..=m {
        s.ic1[i] = 1;
        s.ic2[i] = 2;
        let mut dt = [squared(&a[i], &s.c[1]), squared(&a[i], &s.c[2])];
        if dt[1] < dt[0] {
            s.ic1[i] = 2;
            s.ic2[i] = 1;
            dt.swap(0, 1);
        }
        for l in 3..=k {
            let db = squared(&a[i], &s.c[l]);
            if db < dt[1] {
                if dt[0] <= db {
                    dt[1] = db;
                    s.ic2[i] = l;
                } else {
                    dt[1] = dt[0];
                    s.ic2[i] = s.ic1[i];
                    dt[0] = db;
                    s.ic1[i] = l;
                }
            }
        }
    }
    for l in 1..=k {
        s.nc[l] = 0;
        s.c[l] = vec![0.0; n];
    }
    for i in 1..=m {
        let l = s.ic1[i];
        s.nc[l] += 1;
        for j in 0..n {
            s.c[l][j] += a[i][j];
        }
    }
    if s.nc[1..].contains(&0) {
        return Ok(Err(Fault::Centers));
    }
    for l in 1..=k {
        let aa = s.nc[l] as f64;
        for j in 0..n {
            s.c[l][j] /= aa;
        }
        s.an2[l] = aa / (aa + 1.0);
        s.an1[l] = if 1.0 < aa { aa / (aa - 1.0) } else { HUGE };
        s.itran[l] = 1;
        s.ncp[l] = -1;
    }
    s.indx = 0;
    let mut fault = Some(Fault::Iterations);
    for _ in 0..iterations {
        interrupted(flag)?;
        s.optra();
        if s.indx == m {
            fault = None;
            break;
        }
        if !s.qtran(steps, flag)? {
            break;
        }
        if k == 2 {
            fault = None;
            break;
        }
        for l in 1..=k {
            s.ncp[l] = 0;
        }
    }
    let mut wss = vec![0.0; k + 1];
    for l in 1..=k {
        s.c[l] = vec![0.0; n];
    }
    for i in 1..=m {
        let ii = s.ic1[i];
        for j in 0..n {
            s.c[ii][j] += a[i][j];
        }
    }
    for j in 0..n {
        for l in 1..=k {
            s.c[l][j] /= s.nc[l] as f64;
        }
        for i in 1..=m {
            let ii = s.ic1[i];
            let da = a[i][j] - s.c[ii][j];
            wss[ii] += da * da;
        }
    }
    if let Some(f) = fault {
        return Ok(Err(f));
    }
    Ok(Ok((
        s.ic1[1..].to_vec(),
        wss[1..].to_vec(),
        s.c.split_off(1),
    )))
}

/// A natural below `n`, `n ≥ 1`, from the session's generator.
fn below(ev: &mut Eval<'_>, n: usize) -> usize {
    let n = Natural::from(n as u64);
    let r = ev.draw(|g| g.below(&n));
    usize::try_from(&r).expect("below a usize")
}

/// `init_cluster_centers` (:10622-10646): one point of each class of
/// `ind` (all one class without), the rest distinct points of the first
/// largest class (giac's `_rand(k, list)`, a partial Fisher–Yates).
fn centers(
    ev: &mut Eval<'_>,
    x: &[Vec<f64>],
    ind: Option<&[usize]>,
    k: usize,
) -> Option<Vec<Vec<f64>>> {
    let mut partition: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..x.len() {
        partition
            .entry(ind.map_or(1, |b| b[i]))
            .or_default()
            .push(i);
    }
    let mut ch = Vec::with_capacity(k);
    let mut largest = (0, 0);
    for (&key, class) in &mut partition {
        if class.len() > largest.1 {
            largest = (key, class.len());
        }
        let pos = below(ev, class.len());
        ch.push(class.remove(pos));
    }
    if ch.len() < k {
        let need = k - ch.len();
        let pool = partition.get_mut(&largest.0)?;
        if need > pool.len() {
            return None;
        }
        for i in 0..need {
            let j = i + below(ev, pool.len() - i);
            pool.swap(i, j);
        }
        ch.extend_from_slice(&pool[..need]);
    }
    Some(ch.iter().map(|&i| x[i].clone()).collect())
}

/// One `k` (:10673-10681): starts while a cluster comes out empty, 11 at
/// most; `ind` the first start's classes.
fn run(
    ev: &mut Eval<'_>,
    x: &[Vec<f64>],
    ind: Option<&[usize]>,
    k: usize,
    limits: (usize, i64),
) -> Result<Clusters, EvalError> {
    let mut outcome = Err(Fault::Centers);
    for start in 0..STARTS {
        let Some(c) = centers(ev, x, if start == 0 { ind } else { None }, k) else {
            return Ok(Err(Fault::Centers));
        };
        let flag = ev.context().interrupt_flag();
        outcome = kmns(x, &c, k, limits.0, limits.1, flag)?;
        if outcome != Err(Fault::Centers) {
            break;
        }
    }
    Ok(outcome)
}

/// giac's `kmeans` (:10654-10763), Hartigan's rule when `k` is 0: the
/// number of clusters and the 0-based labels.
fn kmeans_core(
    ev: &mut Eval<'_>,
    x: &[Vec<f64>],
    k: usize,
    (min_k, max_k): (usize, usize),
    limits: (usize, i64),
) -> Result<Result<(usize, Vec<usize>), Fault>, EvalError> {
    let m = x.len();
    if k > 0 {
        return Ok(
            run(ev, x, None, k, limits)?.map(|(b, _, _)| (k, b.iter().map(|l| l - 1).collect()))
        );
    }
    let mut w_prev = 0.0;
    let mut best: Vec<usize> = (1..=m).collect();
    let mut b: Vec<usize> = Vec::new();
    let mut rough = false;
    let mut k = min_k;
    while k <= max_k {
        let ind = (k != min_k).then_some(&b[..]);
        let (labels, wss, _) = match run(ev, x, ind, k, limits)? {
            Ok(r) => r,
            Err(f) => return Ok(Err(f)),
        };
        b = labels;
        let w: f64 = wss.iter().fold(0.0, |s, v| s + v);
        if w > 0.0 && w_prev > 0.0 {
            #[allow(clippy::cast_precision_loss, reason = "point counts")]
            let score = (w_prev / w - 1.0) * (m - k) as f64;
            if score <= 10.0 {
                rough = true;
                k -= 1;
                break;
            }
        }
        w_prev = w;
        best.clone_from(&b);
        if w.total_cmp(&0.0) == Ordering::Equal {
            break;
        }
        k += 1;
    }
    let k = k.min(max_k);
    let labels = if rough { &best } else { &b };
    Ok(Ok((k, labels.iter().map(|l| l - 1).collect())))
}

// ---------------- cluster_crit (optimization.cc:10784-11091) ----------------

/// `kmeans`' index bits (optimization.h:89-105) in `compute_indices`' order
/// (:11040-11058), each with its direction: `-1` the least value, `1` the
/// greatest, `∓2` the least or greatest second difference (an elbow).
const KMEANS_INDICES: [(u32, i8); 14] = [
    (1, -1),   // banfeld-raftery
    (2, -1),   // davies-bouldin
    (4, -1),   // ray-turi
    (8, -1),   // scott-symons
    (16, 1),   // calinski-harabasz
    (32, 1),   // pbm
    (64, 1),   // ratkowsky-lance
    (128, -2), // det ratio
    (256, -2), // log det ratio
    (512, -2), // log ss ratio
    (2048, 2), // ball-hall
    (1024, 2), // ksq-detW
    (4096, 2), // trace-W
    (8192, 2), // trace-WiB
];

/// `index=true` on `kmeans`: Calinski and Harabasz's index (:11669).
const CALINSKI: u32 = 16;

/// giac's `epsilon` (1e-12), `in_doublerref`'s threshold.
const EPSILON: f64 = 1e-12;

/// `cluster_crit::name2index` (:11060-11076): "log-det" names the
/// determinant ratio, as giac's (its log ratio only comes with "all").
fn kmeans_index(name: &str) -> Option<u32> {
    Some(match name {
        "banfeld-raftery" => 1,
        "davies-bouldin" => 2,
        "ray-turi" => 4,
        "scott-symons" => 8,
        "calinski-harabasz" => CALINSKI,
        "pbm" => 32,
        "ratkowsky-lance" => 64,
        "det" | "log-det" => 128,
        "ksq-detW" => 1024,
        "log-ss" => 512,
        "ball-hall" => 2048,
        "trace-W" => 4096,
        "trace-WiB" => 8192,
        "all" => 16383,
        _ => return None,
    })
}

/// A number of `cluster_crit`: exact on exact data, else a double; a sum
/// or product with a double is a double, as giac's.
#[derive(Clone, Debug)]
enum Num {
    Exact(Rational),
    Real(f64),
}

impl Num {
    fn zero() -> Num {
        Num::Exact(Rational::from(Integer::ZERO))
    }

    fn value(&self) -> f64 {
        match self {
            Num::Exact(q) => to_f64(q),
            Num::Real(v) => *v,
        }
    }

    fn is_zero(&self) -> bool {
        match self {
            Num::Exact(q) => q.is_zero(),
            Num::Real(v) => *v == 0.0,
        }
    }

    fn add(&self, o: &Num) -> Num {
        match (self, o) {
            (Num::Exact(a), Num::Exact(b)) => Num::Exact(a + b),
            _ => Num::Real(self.value() + o.value()),
        }
    }

    fn sub(&self, o: &Num) -> Num {
        match (self, o) {
            (Num::Exact(a), Num::Exact(b)) => Num::Exact(a - b),
            _ => Num::Real(self.value() - o.value()),
        }
    }

    fn mul(&self, o: &Num) -> Num {
        match (self, o) {
            (Num::Exact(a), Num::Exact(b)) => Num::Exact(a * b),
            _ => Num::Real(self.value() * o.value()),
        }
    }
}

/// An index value: giac's undef, a real (an infinity giac's `infinity`), an
/// exact rational, or a value not real (the logarithm of a negative or of a
/// rational, complex or symbolic in giac).
#[derive(Clone, Debug)]
enum Value {
    Undef,
    Real(f64),
    Exact(Rational),
    NotReal,
}

impl Value {
    /// A double; NaN is undef.
    fn real(v: f64) -> Value {
        if v.is_nan() {
            Value::Undef
        } else {
            Value::Real(v)
        }
    }

    fn number(&self) -> Option<f64> {
        match self {
            Value::Real(v) => Some(*v),
            Value::Exact(q) => Some(to_f64(q)),
            Value::Undef | Value::NotReal => None,
        }
    }
}

/// giac's `is_one` on a double: a division by it leaves an exact
/// numerator exact.
fn is_one(x: f64) -> bool {
    x.total_cmp(&1.0) == Ordering::Equal
}

/// A double division as giac's: by zero an infinity, `0/0` undef.
fn ratio(a: f64, b: f64) -> Value {
    if b == 0.0 {
        if a == 0.0 {
            Value::Undef
        } else {
            Value::Real(if a > 0.0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            })
        }
    } else {
        Value::real(a / b)
    }
}

/// `in_doublerref`'s determinant (vecteur.cc:9752-9871): the column's
/// largest pivot, entries below the threshold (1e-12, then 1e-12 times the
/// pivot) zero.
#[allow(
    clippy::needless_range_loop,
    reason = "giac's indexed loops, ported literally"
)]
fn doublerref_det(a: &[Vec<f64>]) -> f64 {
    let mut m = a.to_vec();
    let n = m.len();
    let mut det = 1.0;
    let mut threshold = EPSILON;
    for c in 0..n {
        let mut pivot = m[c][c];
        if pivot.abs() < threshold {
            pivot = 0.0;
            m[c][c] = 0.0;
        }
        if c == n - 1 {
            return det * pivot;
        }
        let mut line = c;
        for l in c + 1..n {
            let mut t = m[l][c];
            if t.abs() < threshold {
                t = 0.0;
                m[l][c] = 0.0;
            }
            if t.abs() > pivot.abs() {
                pivot = t;
                line = l;
            }
        }
        if pivot == 0.0 {
            return 0.0;
        }
        threshold = (EPSILON * pivot).abs();
        if line != c {
            m.swap(c, line);
            det = -det;
        }
        for v in &mut m[c] {
            *v /= pivot;
        }
        det *= pivot;
        for l in c + 1..n {
            let f = -m[l][c];
            for j in c + 1..n {
                let t = m[c][j];
                m[l][j] += f * t;
            }
        }
    }
    det
}

/// `_det`: exact on an exact matrix, else `in_doublerref`'s on its doubles
/// (mrref, vecteur.cc:7414-7450).
#[allow(
    clippy::needless_range_loop,
    reason = "giac's indexed loops, ported literally"
)]
fn determinant(a: &[Vec<Num>]) -> Num {
    let exact: Option<Vec<Vec<Rational>>> = a
        .iter()
        .map(|r| {
            r.iter()
                .map(|v| match v {
                    Num::Exact(q) => Some(q.clone()),
                    Num::Real(_) => None,
                })
                .collect()
        })
        .collect();
    let Some(mut m) = exact else {
        let d: Vec<Vec<f64>> = a
            .iter()
            .map(|r| r.iter().map(Num::value).collect())
            .collect();
        return Num::Real(doublerref_det(&d));
    };
    let n = m.len();
    let mut det = Rational::from(Integer::from(1u64));
    for c in 0..n {
        let Some(p) = (c..n).find(|&r| !m[r][c].is_zero()) else {
            return Num::zero();
        };
        if p != c {
            m.swap(c, p);
            det = -det;
        }
        det = &det * &m[c][c];
        for r in c + 1..n {
            let f = &m[r][c] / &m[c][c];
            for j in c..n {
                let t = &f * &m[c][j];
                m[r][j] = &m[r][j] - &t;
            }
        }
    }
    Num::Exact(det)
}

/// The inverse of a double matrix by Gauss–Jordan, `None` when
/// `in_doublerref`'s determinant is zero.
#[allow(
    clippy::needless_range_loop,
    reason = "giac's indexed loops, ported literally"
)]
fn inverse(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    if doublerref_det(a) == 0.0 {
        return None;
    }
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut row = r.clone();
            row.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            row
        })
        .collect();
    for c in 0..n {
        let p = (c..n)
            .max_by(|&x, &y| m[x][c].abs().total_cmp(&m[y][c].abs()))
            .unwrap_or(c);
        m.swap(c, p);
        let pv = m[c][c];
        for v in &mut m[c] {
            *v /= pv;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                for j in 0..2 * n {
                    let t = m[c][j];
                    m[r][j] -= f * t;
                }
            }
        }
    }
    Some(m.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// `is_blockmatrix`'s split, recursively (vecteur.cc:11907, `minv`
/// :11996): the connected components of `i ~ j` when `a[i][j]` or
/// `a[j][i]` is non-zero, each increasing.
fn blocks(a: &[Vec<f64>]) -> Vec<Vec<usize>> {
    let n = a.len();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let (mut component, mut queue) = (vec![start], vec![start]);
        while let Some(i) = queue.pop() {
            for j in 0..n {
                if !seen[j] && (a[i][j] != 0.0 || a[j][i] != 0.0) {
                    seen[j] = true;
                    component.push(j);
                    queue.push(j);
                }
            }
        }
        component.sort_unstable();
        out.push(component);
    }
    out
}

/// A scatter matrix of `cluster_crit`.
type Scatter = Vec<Vec<Num>>;

/// giac's `cluster_crit`: the data's fixed parts, the clustering of the
/// current `k` (`init`), and each index's values and optimum so far.
struct Crit<'a> {
    pts: &'a [Vec<Num>],
    x: &'a [Vec<f64>],
    p: usize,
    /// The barycentre `mu()`.
    mu0: Vec<Num>,
    /// `_T`: the identity until first computed, then kept (:10880 reads it).
    t: Vec<Vec<Num>>,
    et: f64,
    k: usize,
    ci: Vec<usize>,
    size: Vec<usize>,
    mu: Vec<Vec<f64>>,
    wgss: Vec<f64>,
    wgss_sum: f64,
    have_t: bool,
    /// Each cluster's `WG(k)` and their sum, once computed for this `k`.
    wg: Option<(Vec<Scatter>, Scatter)>,
    /// `_vmap`: each index's values, by position.
    values: Vec<Vec<Value>>,
    /// `_optvals`: each index's suggested `k` and its value.
    optimal: BTreeMap<usize, (usize, f64)>,
}

#[allow(
    clippy::cast_precision_loss,
    reason = "point and cluster counts as giac's int times double"
)]
impl<'a> Crit<'a> {
    /// The constructor (:10784-10801).
    fn new(pts: &'a [Vec<Num>], x: &'a [Vec<f64>]) -> Crit<'a> {
        let (m, p) = (pts.len(), pts[0].len());
        let count = Rational::from(Integer::from(m as u64));
        let mu0: Vec<Num> = (0..p)
            .map(
                |j| match pts.iter().fold(Num::zero(), |s, q| s.add(&q[j])) {
                    Num::Exact(s) => Num::Exact(&s / &count),
                    Num::Real(s) => Num::Real(s / m as f64),
                },
            )
            .collect();
        let t = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        if i == j {
                            Num::Exact(Rational::from(Integer::from(1u64)))
                        } else {
                            Num::zero()
                        }
                    })
                    .collect()
            })
            .collect();
        let et = pts.iter().fold(0.0, |s, q| {
            let d = q.iter().zip(&mu0).fold(Num::zero(), |d, (a, b)| {
                let t = a.sub(b);
                d.add(&t.mul(&t))
            });
            s + d.value()
        });
        Crit {
            pts,
            x,
            p,
            mu0,
            t,
            et,
            k: 0,
            ci: Vec::new(),
            size: Vec::new(),
            mu: Vec::new(),
            wgss: Vec::new(),
            wgss_sum: 0.0,
            have_t: false,
            wg: None,
            values: Vec::new(),
            optimal: BTreeMap::new(),
        }
    }

    /// `init` (:10809-10831): labels `1 … k`, centres, sums of squares.
    fn init(&mut self, k: usize, labels: &[usize], centers: Vec<Vec<f64>>, wgss: Vec<f64>) {
        self.k = k;
        self.ci = labels.iter().map(|l| l - 1).collect();
        self.size = vec![0; k];
        for &c in &self.ci {
            self.size[c] += 1;
        }
        self.mu = centers;
        self.wgss_sum = wgss.iter().fold(0.0, |s, w| s + w);
        self.wgss = wgss;
        self.have_t = false;
        self.wg = None;
    }

    /// `T()` (:10832-10850), kept in `_T`.
    fn total(&mut self) -> Vec<Vec<Num>> {
        if !self.have_t {
            for i in 0..self.p {
                for j in i..self.p {
                    let v = self.pts.iter().fold(Num::zero(), |s, q| {
                        s.add(&q[i].sub(&self.mu0[i]).mul(&q[j].sub(&self.mu0[j])))
                    });
                    self.t[j][i] = v.clone();
                    self.t[i][j] = v;
                }
            }
            self.have_t = true;
        }
        self.t.clone()
    }

    /// `TSS(j)` (:10851-10861).
    fn tss(&self, j: usize) -> Num {
        self.pts.iter().fold(Num::zero(), |s, q| {
            let d = q[j].sub(&self.mu0[j]);
            s.add(&d.mul(&d))
        })
    }

    /// `WG(l)` (:10862-10884): each cluster's scatter above the diagonal;
    /// below it giac copies `_T`'s entry (the identity's 0 before any
    /// `T()`), and the sum adds those too.
    fn within(&mut self) -> &(Vec<Scatter>, Scatter) {
        if self.wg.is_none() {
            let p = self.p;
            let mut each = Vec::with_capacity(self.k);
            let mut sum: Option<Vec<Vec<Num>>> = None;
            for c in 0..self.k {
                let mut w = vec![vec![Num::zero(); p]; p];
                for i in 0..p {
                    for j in i..p {
                        let mut v = 0.0;
                        for (l, q) in self.x.iter().enumerate() {
                            if self.ci[l] == c {
                                v += (q[i] - self.mu[c][i]) * (q[j] - self.mu[c][j]);
                            }
                        }
                        w[i][j] = Num::Real(v);
                        if i != j {
                            w[j][i] = self.t[i][j].clone();
                        }
                    }
                }
                sum = Some(match sum {
                    None => w.clone(),
                    Some(s) => s
                        .iter()
                        .zip(&w)
                        .map(|(a, b)| a.iter().zip(b).map(|(u, v)| u.add(v)).collect())
                        .collect(),
                });
                each.push(w);
            }
            self.wg = Some((each, sum.unwrap_or_default()));
        }
        self.wg.as_ref().expect("computed")
    }

    fn within_det(&mut self, c: Option<usize>) -> f64 {
        let (each, sum) = self.within();
        determinant(c.map_or(sum, |c| &each[c])).value()
    }

    /// `BG()` (:10885-10900).
    #[allow(
        clippy::needless_range_loop,
        reason = "giac's indexed loops, ported literally"
    )]
    fn between(&self) -> Vec<Vec<f64>> {
        let p = self.p;
        let mut b = vec![vec![0.0; p]; p];
        for i in 0..p {
            for j in i..p {
                let mut v = 0.0;
                for c in 0..self.k {
                    v += self.size[c] as f64
                        * (self.mu[c][i] - self.mu0[i].value())
                        * (self.mu[c][j] - self.mu0[j].value());
                }
                b[i][j] = v;
                b[j][i] = v;
            }
        }
        b
    }

    /// `BGSS(j)` (:10901-10911).
    fn bgss(&self, j: usize) -> f64 {
        (0..self.k).fold(0.0, |s, c| {
            let d = self.mu[c][j] - self.mu0[j].value();
            s + self.size[c] as f64 * (d * d)
        })
    }

    fn bgss_sum(&self) -> f64 {
        (0..self.p).fold(0.0, |s, j| s + self.bgss(j))
    }

    /// `Delta(k, l)` (:10802-10808): the centres' distance.
    fn delta(&self, k: usize, l: usize) -> f64 {
        squared(&self.mu[k], &self.mu[l]).sqrt()
    }

    /// `DB(dir)` (:10912-10925): the least (`dir < 0`) or greatest distance
    /// of two centres.
    fn db(&self, least: bool) -> f64 {
        let mut v = if least {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        for k in 0..self.k {
            for l in k + 1..self.k {
                let d = self.delta(k, l);
                v = if least { v.min(d) } else { v.max(d) };
            }
        }
        v
    }

    /// The index of bit `bit` (:10926-11024).
    #[allow(
        clippy::too_many_lines,
        reason = "the fourteen indices, in giac's order"
    )]
    fn index(&mut self, bit: u32) -> Value {
        let (k, n) = (self.k as f64, self.x.len() as f64);
        match bit {
            // banfeld-raftery
            1 => {
                let mut r = 0.0;
                for c in 0..self.k {
                    if self.wgss[c] == 0.0 {
                        return Value::Undef;
                    }
                    let size = self.size[c] as f64;
                    r += size * (self.wgss[c] / size).ln();
                }
                Value::real(r)
            }
            // davies-bouldin
            2 => {
                let mut sum = 0.0;
                for c in 0..self.k {
                    let mut most = f64::NEG_INFINITY;
                    for l in (0..self.k).filter(|&l| l != c) {
                        let s =
                            self.wgss[c] / self.size[c] as f64 + self.wgss[l] / self.size[l] as f64;
                        match ratio(s, self.delta(c, l)) {
                            Value::Real(v) => most = most.max(v),
                            _ => return Value::Undef,
                        }
                    }
                    sum += most;
                }
                Value::real(sum / k)
            }
            // ray-turi
            4 => {
                let d = self.db(true);
                ratio(self.wgss_sum, n * (d * d))
            }
            // scott-symons
            8 => {
                let mut r = 0.0;
                for c in 0..self.k {
                    let d = self.within_det(Some(c));
                    let size = self.size[c] as f64;
                    if d == 0.0 || d / size < 0.0 {
                        return Value::Undef;
                    }
                    r += size * (d / size).ln();
                }
                Value::real(r)
            }
            // calinski-harabasz
            16 if self.k <= 1 => Value::Undef,
            16 => ratio((n - k) * self.bgss_sum(), (k - 1.0) * self.wgss_sum),
            // pbm
            32 => match ratio(self.et * self.db(false), k * self.wgss_sum) {
                Value::Real(v) => Value::real(v * v),
                v => v,
            },
            // ratkowsky-lance
            64 => {
                let mut r = 0.0;
                for j in 0..self.p {
                    let t = self.tss(j);
                    if t.is_zero() {
                        return Value::Undef;
                    }
                    r += self.bgss(j) / t.value();
                }
                Value::real((r / self.p as f64).sqrt())
            }
            // det ratio: `det(T)/det(WG)`, exact when `det(WG)` is 1.0
            // (giac's division by a one).
            128 => {
                let den = self.within_det(None);
                if den == 0.0 {
                    return Value::Undef;
                }
                match determinant(&self.total()) {
                    Num::Exact(q) if is_one(den) => Value::Exact(q),
                    num => Value::real(num.value() / den),
                }
            }
            // log det ratio: `T()` first.
            256 => {
                let num = determinant(&self.total());
                let den = self.within_det(None);
                if num.is_zero() || den == 0.0 {
                    return Value::Undef;
                }
                match num {
                    Num::Exact(q) if is_one(den) => {
                        if q.is_one() {
                            Value::Exact(Rational::from(Integer::ZERO))
                        } else {
                            Value::NotReal
                        }
                    }
                    num => {
                        let r = num.value() / den;
                        if r < 0.0 {
                            Value::NotReal
                        } else {
                            Value::real(n * r.ln())
                        }
                    }
                }
            }
            // log ss ratio
            512 => {
                let (b, w) = (self.bgss_sum(), self.wgss_sum);
                if b == 0.0 || w == 0.0 {
                    Value::Undef
                } else {
                    Value::real((b / w).ln())
                }
            }
            // ksq-detW
            1024 => Value::real(k * k * self.within_det(None)),
            // ball-hall
            2048 => {
                let s = (0..self.k).fold(0.0, |s, c| s + self.wgss[c] / self.size[c] as f64);
                Value::real(s / k)
            }
            // trace-W
            4096 => Value::real(self.wgss_sum),
            // trace-WiB: `minv` inverts each block; a singular block's LU
            // inverse has infinities, its trace undef where `BG` is zero on
            // it, else giac's `infinity` (probes, prototype `trace_WiB`).
            _ => {
                let w: Vec<Vec<f64>> = self
                    .within()
                    .1
                    .iter()
                    .map(|r| r.iter().map(Num::value).collect())
                    .collect();
                let b = self.between();
                let mut t = 0.0;
                for block in blocks(&w) {
                    let sub: Vec<Vec<f64>> = block
                        .iter()
                        .map(|&i| block.iter().map(|&j| w[i][j]).collect())
                        .collect();
                    let Some(inv) = inverse(&sub) else {
                        if block.iter().all(|&i| block.iter().all(|&j| b[i][j] == 0.0)) {
                            return Value::Undef;
                        }
                        t = f64::INFINITY;
                        continue;
                    };
                    for (u, &i) in block.iter().enumerate() {
                        for (v, &j) in block.iter().enumerate() {
                            t += inv[u][v] * b[j][i];
                        }
                    }
                }
                Value::real(t)
            }
        }
    }

    /// `update` (:11025-11039): a value, and the optimum of its index.
    fn update(&mut self, pos: usize, value: Value, dir: i8) {
        if self.values.len() <= pos {
            self.values.resize(pos + 1, Vec::new());
        }
        self.values[pos].push(value);
        let v = &self.values[pos];
        let elbow = dir.abs() == 2;
        let x = if elbow {
            let [a, b, c] = match v.as_slice() {
                [.., a, b, c] => [a, b, c].map(Value::number),
                _ => return,
            };
            match (a, b, c) {
                (Some(a), Some(b), Some(c)) => c - 2.0 * b + a,
                _ => return,
            }
        } else {
            match v[v.len() - 1].number() {
                Some(x) => x,
                None => return,
            }
        };
        if x.is_nan() {
            return;
        }
        let better = self
            .optimal
            .get(&pos)
            .is_none_or(|&(_, o)| (dir < 0 && o > x) || (dir > 0 && o < x));
        if better {
            self.optimal
                .insert(pos, (if elbow { self.k - 1 } else { self.k }, x));
        }
    }

    /// `compute_indices` (:11040-11058).
    fn compute(&mut self, bits: u32) {
        let mut pos = 0;
        for (bit, dir) in KMEANS_INDICES {
            if bits & bit != 0 {
                let v = self.index(bit);
                self.update(pos, v, dir);
                pos += 1;
            }
        }
    }

    /// `optimal_number_of_clusters` (:11077-11091): the `k` most indices
    /// suggest, the least of a tie; 0 for none.
    fn optimum(&self) -> usize {
        let mut votes: BTreeMap<usize, usize> = BTreeMap::new();
        for &(k, _) in self.optimal.values() {
            *votes.entry(k).or_default() += 1;
        }
        let mut best = (0, 0);
        for (k, n) in votes {
            if best.0 == 0 || n > best.1 {
                best = (k, n);
            }
        }
        best.0
    }
}

/// The answer of `kmeans`' index path: the number of clusters, its labels,
/// the `k` tried and the first index's values.
struct Indexed {
    k: usize,
    labels: Vec<usize>,
    kvals: Vec<usize>,
    values: Vec<Value>,
}

/// `kmeans`' index path (:10682-10750): each `k` from `min_k` up, its first
/// start from the previous `k`'s classes; Hartigan's rough optimum caps
/// `max_k` at `2k`; the labels of the `k` the indices of `bits` suggest
/// (none for no suggestion, giac's empty `b_map[0]`).
fn kmeans_indexed(
    ev: &mut Eval<'_>,
    (x, pts): (&[Vec<f64>], &[Vec<Num>]),
    bits: u32,
    (min_k, max_k): (usize, usize),
    limits: (usize, i64),
) -> Result<Result<Indexed, Fault>, EvalError> {
    let m = x.len();
    let mut crit = Crit::new(pts, x);
    let (mut w_prev, mut rough, mut max_k) = (0.0, false, max_k);
    let mut b: Vec<usize> = Vec::new();
    let mut runs: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut kvals = Vec::new();
    let mut k = min_k;
    while k <= max_k {
        let ind = (k != min_k).then_some(&b[..]);
        let (labels, wss, centers) = match run(ev, x, ind, k, limits)? {
            Ok(r) => r,
            Err(f) => return Ok(Err(f)),
        };
        b = labels;
        if !rough {
            let w: f64 = wss.iter().fold(0.0, |s, v| s + v);
            #[allow(clippy::cast_precision_loss, reason = "point counts")]
            if w > 0.0 && w_prev > 0.0 && (w_prev / w - 1.0) * (m - k) as f64 <= 10.0 {
                rough = true;
                max_k = max_k.min(2 * k);
            }
            w_prev = w;
        }
        interrupted(ev.context().interrupt_flag())?;
        crit.init(k, &b, centers, wss);
        crit.compute(bits);
        kvals.push(k);
        runs.insert(k, b.iter().map(|l| l - 1).collect());
        k += 1;
    }
    let k = crit.optimum();
    Ok(Ok(Indexed {
        k,
        labels: runs.remove(&k).unwrap_or_default(),
        kvals,
        values: crit.values.into_iter().next().unwrap_or_default(),
    }))
}

/// An index value as giac prints it: a double, a rational, `infinity`, or
/// "not defined" (undef; a value not real never comes first: the log det
/// ratio's bit only comes with "all").
fn index_value(v: &Value) -> Expr {
    match v {
        Value::Real(r) => Float::from_f64(*r).map_or_else(
            || Symbol::COMPLEX_INFINITY.into(),
            |f| Expr::from(Ball::exact(f)),
        ),
        Value::Exact(q) => Expr::from(q.clone()),
        Value::Undef | Value::NotReal => Expr::string(NOT_DEFINED),
    }
}

/// `KMeans[data, options…]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's checks and options in its order"
)]
pub(in crate::native) fn kmeans(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((d, options)) = e.args().split_first() else {
        return Ok(None);
    };
    let Some(data) = data(ev, d) else {
        return Ok(None);
    };
    if data.items.is_empty() {
        return refuse(ev, e, NO_DATA);
    }
    let Some(points) = data
        .points
        .iter()
        .map(|p| p.as_ref().map(|p| p.iter().map(|c| c.value).collect()))
        .collect::<Option<Vec<Vec<f64>>>>()
    else {
        return refuse(ev, e, NUMERIC);
    };
    // The points as giac's gens, for the indices: exact where rational.
    let exact: Vec<Vec<Num>> = data
        .points
        .iter()
        .flatten()
        .map(|p| {
            p.iter()
                .map(|c| c.exact.clone().map_or(Num::Real(c.value), Num::Exact))
                .collect()
        })
        .collect();
    let np = wide(points.len());
    let (mut out, mut min_k, mut max_k, mut k, mut maxiter) = (Output::List, 2, np - 1, 0, 100);
    let mut stays = false;
    for o in options {
        let step = if let Some((lh, rh)) = equation(o) {
            if is(lh, Symbol::COUNT, "count") {
                count(rh).map(|c| match c {
                    Count::K(v) => k = v,
                    Count::Range(a, b) => (min_k, max_k) = (a, b),
                })
            } else if is(lh, Symbol::COUNT_GREATER, "count_sup") {
                bound(rh).map(|v| max_k = v.min(np - 1))
            } else if is(lh, Symbol::COUNT_LESS, "count_inf") {
                bound(rh).map(|v| min_k = v.max(2))
            } else if is(lh, Symbol::LIMIT, "limit") {
                match integer(rh) {
                    Some(v) if v >= 1 => {
                        maxiter = v;
                        Ok(())
                    }
                    _ => Err(POSITIVE),
                }
            } else if is(lh, Symbol::FIRST_POSITION, "index") {
                // An index runs as `k = -bits`, `false` as Hartigan's rule.
                index(rh, kmeans_index, CALINSKI).map(|v| k = -i64::from(v))
            } else if drawing(lh) {
                stays = true;
                Ok(())
            } else if named(lh, "output") {
                output(rh, false).map(|v| out = v)
            } else {
                // giac ignores another option (:11620-11700).
                Ok(())
            }
        } else if let Some(v) = integer(o) {
            if v <= 0 {
                Err(POSITIVE)
            } else {
                k = v;
                Ok(())
            }
        } else {
            Err(OPTION)
        };
        if let Err(why) = step {
            return refuse(ev, e, why);
        }
    }
    if stays || out == Output::Out {
        return Ok(None);
    }
    if out == Output::Index && k >= 0 {
        return refuse(ev, e, NO_INDEX);
    }
    max_k = max_k.min(np - 1);
    if min_k > max_k {
        return refuse(ev, e, K_RANGE);
    }
    if min_k == max_k {
        k = max_k;
    }
    let m = points.len();
    if k > 0 && (k <= 1 || k >= np) {
        // AS 136's fault 3: giac answers [] (:10291, :11721).
        return Ok(Some(list(Vec::new())));
    }
    // 2 ≤ min_k ≤ max_k < m, or 0 < k < m: all fit a usize.
    let to = |v: i64| usize::try_from(v).unwrap_or(0);
    let (kk, lo, hi) = (to(k), to(min_k), to(max_k));
    let bits = u32::try_from(-k).unwrap_or(0);
    let (sum, top) = if kk > 0 {
        (kk as u128, kk)
    } else {
        let (l, h) = (lo as u128, hi as u128);
        ((h * (h + 1) - l * l.saturating_sub(1)) / 2, hi)
    };
    let mn = (m as u128) * (data.n as u128);
    let terms = ev.context().limits().terms as u128;
    let maxiter = to(maxiter);
    // The indices per k: the scatter matrices m n², the centres' pairs k² n,
    // determinants and inverses n³.
    let (n, h) = (data.n as u128, hi as u128);
    let crit = if bits > 0 {
        (h + 1).saturating_sub(lo as u128) * (mn * n + h * h * n + n * n * n)
    } else {
        0
    };
    if mn * sum * STARTS as u128 > terms
        || mn * (top as u128) * (maxiter as u128) > terms
        || crit > terms
    {
        return refuse(ev, e, TOO_MANY);
    }
    let steps = i64::try_from(terms / (data.n as u128)).unwrap_or(i64::MAX);
    ev.volatile();
    let outcome = if bits > 0 {
        match kmeans_indexed(ev, (&points, &exact), bits, (lo, hi), (maxiter, steps))? {
            Ok(ix) if out == Output::Index => {
                // `mtran([kvals, values])` (:10737-10745).
                return Ok(Some(list(
                    ix.kvals
                        .iter()
                        .zip(&ix.values)
                        .map(|(&k, v)| {
                            list(vec![Expr::from(Integer::from(k as u64)), index_value(v)])
                        })
                        .collect(),
                )));
            }
            r => r.map(|ix| (ix.k, ix.labels)),
        }
    } else {
        kmeans_core(ev, &points, kk, (lo, hi), (maxiter, steps))?
    };
    match outcome {
        // A k fixed with an index: giac's empty `ind` (count=a..a).
        Ok(_) if out == Output::Index => Ok(Some(list(Vec::new()))),
        Ok((k, labels)) => Ok(Some(answer(out, k, &labels, data.items))),
        Err(f) => {
            if f != Fault::Size {
                ev.warn(Warning::Native {
                    message: (if f == Fault::Centers {
                        CENTERS
                    } else {
                        ITERATIONS
                    })
                    .into(),
                    expr: e.clone(),
                });
            }
            Ok(Some(list(Vec::new())))
        }
    }
}

// ---------------- hclust (optimization.cc:11102-11443) ----------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Linkage {
    Single,
    Complete,
    Average,
    Weighted,
    Ward,
}

/// How two points' distance is taken.
enum Distance<'a> {
    /// Squared Euclidean (giac's `longueur2`).
    Squared(&'a [Option<Vec<Coordinate>>]),
    /// `N[f[p, q]]`.
    Function(&'a Expr, &'a [Expr]),
    /// Strings' Levenshtein distance (giac's for a list of strings, :11741).
    Levenshtein(&'a [Expr]),
}

/// giac's `LevenshteinDistance` (:11519-11536), on bytes as its
/// `std::string::at`.
fn levenshtein(s: &[u8], t: &[u8]) -> usize {
    let (s, t) = if s.len() > t.len() { (t, s) } else { (s, t) };
    let mut v0: Vec<usize> = (0..=t.len()).collect();
    let mut v1 = vec![0; t.len() + 1];
    for (i, &x) in s.iter().enumerate() {
        v1[0] = i + 1;
        for (j, &y) in t.iter().enumerate() {
            v1[j + 1] = (v0[j + 1] + 1)
                .min(v1[j] + 1)
                .min(v0[j] + usize::from(x != y));
        }
        std::mem::swap(&mut v0, &mut v1);
    }
    v0[t.len()]
}

/// `linkage`'s answer: the number of clusters, the labels by first
/// appearance, the `k` values of the index path and the silhouette's.
struct Linked {
    k: usize,
    labels: Vec<usize>,
    kvals: usize,
    values: Vec<(usize, f64)>,
}

/// A point's count and sum of distances: to its cluster, and to each other
/// cluster by root (giac's `slh_a`, `slh_b`).
type Sums = (usize, f64);

/// The silhouette as giac takes it (:11380-11404): each point's
/// `(b − a)/max(b, a)` (`0` alone in its cluster), each cluster's mean, the
/// largest over clusters and `DBL_MIN` at least.
#[allow(clippy::cast_precision_loss, reason = "point counts")]
fn silhouette(roots: &[usize], own: &[Sums], other: &[BTreeMap<usize, Sums>]) -> f64 {
    let mut s = vec![(0usize, 0.0f64); roots.len()];
    for (i, &r) in roots.iter().enumerate() {
        let (count, sum) = own[i];
        let ai = if count == 0 { 0.0 } else { sum / count as f64 };
        let minb = other[i]
            .values()
            .map(|&(c, d)| d / c as f64)
            .fold(f64::MAX, |m, b| if b < m { b } else { m });
        let si = if count == 0 {
            0.0
        } else {
            (minb - ai) / if minb < ai { ai } else { minb }
        };
        s[r].0 += 1;
        s[r].1 += si;
    }
    s.iter()
        .filter(|(c, _)| *c > 0)
        .map(|&(c, sum)| sum / c as f64)
        .fold(f64::MIN_POSITIVE, |m, v| if v > m { v } else { m })
}

/// A merge: its height and the two points as giac records them.
type Merge = (f64, usize, usize);

struct Hclust<'a> {
    n: usize,
    distance: Distance<'a>,
    /// The pairs' distances, `−1` not yet taken.
    cache: Vec<f64>,
}

/// The exact squared distance of rational points, else the `f64` one; `None`
/// when a point is no real point.
fn longueur2(p: Option<&Vec<Coordinate>>, q: Option<&Vec<Coordinate>>) -> Option<f64> {
    let (p, q) = (p?, q?);
    let exact: Option<Vec<(&Rational, &Rational)>> = p
        .iter()
        .zip(q)
        .map(|(a, b)| Some((a.exact.as_ref()?, b.exact.as_ref()?)))
        .collect();
    if let Some(pairs) = exact {
        let zero = Rational::from(Integer::ZERO);
        let s = pairs.iter().fold(zero, |s, (a, b)| {
            let t = *a - *b;
            &s + &(&t * &t)
        });
        let v = to_f64(&s);
        return v.is_finite().then_some(v);
    }
    let v = p.iter().zip(q).fold(0.0, |s, (a, b)| {
        let t = a.value - b.value;
        s + t * t
    });
    v.is_finite().then_some(v)
}

impl Hclust<'_> {
    fn slot(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i < j { (i, j) } else { (j, i) };
        i * (2 * self.n - i - 1) / 2 + (j - i - 1)
    }

    /// The distance, `−1` when it is no real number `≥ 0` (`comp_dist`
    /// :10765, `dist` :11120).
    fn compute(&self, ev: &mut Eval<'_>, i: usize, j: usize) -> Result<f64, EvalError> {
        let v = match &self.distance {
            Distance::Squared(points) => longueur2(points[i].as_ref(), points[j].as_ref()),
            Distance::Function(f, items) => {
                let call = Expr::apply(f, &[items[i].clone(), items[j].clone()]);
                let v = ev.eval(&Expr::call(Symbol::N, [call]))?;
                coordinate(ev, &v).map(|c| c.value)
            }
            #[allow(clippy::cast_precision_loss, reason = "string lengths")]
            Distance::Levenshtein(items) => items[i]
                .as_str()
                .zip(items[j].as_str())
                .map(|(p, q)| levenshtein(p.as_bytes(), q.as_bytes()) as f64),
        };
        Ok(v.filter(|v| *v >= 0.0).unwrap_or(-1.0))
    }

    fn dist(
        &mut self,
        ev: &mut Eval<'_>,
        i: usize,
        j: usize,
        cache: bool,
    ) -> Result<f64, EvalError> {
        if i == j {
            return Ok(0.0);
        }
        if !cache {
            return self.compute(ev, i, j);
        }
        let s = self.slot(i, j);
        if self.cache[s] < 0.0 {
            let v = self.compute(ev, i, j)?;
            if v < 0.0 {
                return Ok(-1.0);
            }
            self.cache[s] = v;
        }
        Ok(self.cache[s])
    }

    /// Lance and Williams' update (:11134-11157); `−1` fails.
    #[allow(clippy::too_many_arguments, reason = "giac's formula")]
    #[allow(clippy::cast_precision_loss, reason = "cluster sizes")]
    fn formula(
        &mut self,
        ev: &mut Eval<'_>,
        (a, b, x): (usize, usize, usize),
        (na, nb, nx): (usize, usize, usize),
        dab: f64,
        meth: Linkage,
    ) -> Result<f64, EvalError> {
        let dax = self.dist(ev, a, x, true)?;
        let dbx = self.dist(ev, b, x, true)?;
        if dax < 0.0 || dbx < 0.0 || dab < 0.0 || dax + dbx < dab {
            return Ok(-1.0);
        }
        let (na, nb, nx) = (na as f64, nb as f64, nx as f64);
        Ok(match meth {
            Linkage::Single => dax.min(dbx),
            Linkage::Complete => dax.max(dbx),
            Linkage::Average => (na * dax + nb * dbx) / (na + nb),
            Linkage::Weighted => f64::midpoint(dax, dbx),
            Linkage::Ward => ((na + nx) * dax + (nb + nx) * dbx - nx * dab) / (na + nb + nx),
        })
    }

    /// Single linkage, Prim's tree from point 0 (:11172-11209).
    fn mst(&mut self, ev: &mut Eval<'_>) -> Result<Option<Vec<Merge>>, EvalError> {
        let n = self.n;
        let mut s: Vec<usize> = (0..n).collect();
        let mut c = s[0];
        let mut dg = Vec::with_capacity(n.saturating_sub(1));
        let mut d = vec![0.0; 2 * n];
        for &v in &s {
            if v != c {
                d[v] = f64::MAX;
            }
        }
        let mut first = true;
        for _ in 1..n {
            interrupted(ev.context().interrupt_flag())?;
            let mut min_d = f64::MAX;
            let mut next = None;
            let off = if first { 0 } else { n };
            s.retain(|&v| v != c);
            for &v in &s {
                let dval = d[v + off].min(self.dist(ev, v, c, false)?);
                if dval < 0.0 {
                    return Ok(None);
                }
                if dval < min_d {
                    next = Some(v);
                    min_d = dval;
                }
                d[v + n - off] = dval;
            }
            let Some(nx) = next else {
                return Ok(None);
            };
            dg.push((d[nx + n - off], c, nx));
            c = nx;
            first = !first;
        }
        Ok(Some(dg))
    }

    /// The nearest-neighbour chain (:11210-11261).
    fn nn_chain(
        &mut self,
        ev: &mut Eval<'_>,
        meth: Linkage,
    ) -> Result<Option<Vec<Merge>>, EvalError> {
        let mut s: Vec<usize> = (0..self.n).collect();
        let mut dg = Vec::with_capacity(self.n.saturating_sub(1));
        let mut chain: Vec<usize> = Vec::new();
        let mut size = vec![1usize; self.n];
        let (mut a, mut b, mut c);
        c = 0;
        while s.len() > 1 {
            if chain.len() <= 3 {
                a = s[0];
                chain = vec![a];
                b = s[1];
            } else {
                let l = chain.len();
                a = chain[l - 4];
                b = chain[l - 3];
                chain.truncate(l - 3);
            }
            loop {
                interrupted(ev.context().interrupt_flag())?;
                let mut mind = f64::MAX;
                for &x in &s {
                    if x != a {
                        let d = self.dist(ev, x, a, true)?;
                        if d < mind || (x == b && d.total_cmp(&mind) == Ordering::Equal) {
                            c = x;
                            mind = d;
                        }
                    }
                }
                b = a;
                a = c;
                chain.push(a);
                let l = chain.len();
                if l >= 3 && a == chain[l - 3] {
                    break;
                }
                // No reciprocal pair: not a reducible distance.
                if l > s.len() + 3 {
                    return Ok(None);
                }
            }
            let dab = self.dist(ev, a, b, true)?;
            if dab < 0.0 {
                return Ok(None);
            }
            dg.push((dab, a, b));
            s.retain(|&v| v != b);
            for &x in &s {
                if x == a {
                    continue;
                }
                let f = self.formula(ev, (a, b, x), (size[a], size[b], size[x]), dab, meth)?;
                if f < 0.0 {
                    return Ok(None);
                }
                let slot = self.slot(x, a);
                self.cache[slot] = f;
            }
            size[a] += size[b];
        }
        Ok(Some(dg))
    }

    /// The first largest gap of heights (:11158-11171).
    fn k_low(&self, dg: &[Merge]) -> usize {
        let (mut dprev, mut max_step, mut k) = (0.0, -1.0, 0);
        for (i, &(d, _, _)) in dg.iter().enumerate() {
            if d - dprev > max_step {
                max_step = d - dprev;
                k = self.n - i;
            }
            dprev = d;
        }
        k
    }

    /// `linkage` (:11262-11443): `None` fails. With `k_given` 0 and index
    /// bits, giac's index path: from `max_k` clusters down, each step's `k`,
    /// and for the silhouette its value from giac's running sums, the labels
    /// those of the first best value (strictly larger replaces it).
    /// Departure: giac's first step below `max_k` asserts (:11325) and
    /// aborts libgiac; ours starts at `N − 1` clusters at most.
    #[allow(clippy::too_many_lines, reason = "giac's loop in its order")]
    fn linkage(
        &mut self,
        ev: &mut Eval<'_>,
        meth: Linkage,
        k_given: usize,
        (min_k, max_k): (usize, usize),
        indf: u32,
    ) -> Result<Option<Linked>, EvalError> {
        let dg = if meth == Linkage::Single {
            self.mst(ev)?
        } else {
            self.nn_chain(ev, meth)?
        };
        let Some(mut dg) = dg else {
            return Ok(None);
        };
        dg.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let n = self.n;
        let mut parent: Vec<usize> = (0..n).collect();
        let find = |parent: &[usize], mut i: usize| {
            while parent[i] != i {
                i = parent[i];
            }
            i
        };
        let roots = |parent: &[usize]| (0..n).map(|i| find(parent, i)).collect::<Vec<_>>();
        let indexing = k_given == 0 && indf > 0;
        let sil = indexing && indf & SILHOUETTE != 0;
        let max_k = if indexing {
            // True distances again: the chain cached Lance and Williams'.
            self.cache.fill(-1.0);
            max_k.min(n.saturating_sub(1))
        } else {
            max_k
        };
        let mut k = n;
        let target = if k_given > 0 {
            k_given
        } else {
            min_k.max(self.k_low(&dg))
        };
        let mut own: Vec<Sums> = vec![(0, 0.0); n];
        let mut other: Vec<BTreeMap<usize, Sums>> = vec![BTreeMap::new(); n];
        let (mut kvals, mut values) = (0, Vec::new());
        let (mut best, mut best_k, mut snapshot): (Option<f64>, usize, Option<Vec<usize>>) =
            (None, 0, None);
        for &(_, a, b) in &dg {
            if k == target {
                break;
            }
            k -= 1;
            let (ra, rb) = (find(&parent, a), find(&parent, b));
            if !indexing || k > max_k {
                if ra != rb {
                    parent[rb] = ra;
                }
                continue;
            }
            interrupted(ev.context().interrupt_flag())?;
            if k < max_k {
                if sil {
                    for (i, &ci) in roots(&parent).iter().enumerate() {
                        if ci == ra || ci == rb {
                            let o = if ci == ra { rb } else { ra };
                            if let Some((c, d)) = other[i].remove(&o) {
                                own[i].0 += c;
                                own[i].1 += d;
                            }
                        }
                    }
                }
                if ra != rb {
                    parent[rb] = ra;
                }
                if sil {
                    for (i, &ci) in roots(&parent).iter().enumerate() {
                        if ci != ra
                            && let Some((c, d)) = other[i].remove(&rb)
                        {
                            let e = other[i].entry(ra).or_insert((0, 0.0));
                            e.0 += c;
                            e.1 += d;
                        }
                    }
                }
            } else {
                if ra != rb {
                    parent[rb] = ra;
                }
                if sil {
                    own.fill((0, 0.0));
                    let r = roots(&parent);
                    for i in 0..n {
                        interrupted(ev.context().interrupt_flag())?;
                        for j in i + 1..n {
                            let d = self.dist(ev, i, j, true)?;
                            if r[i] == r[j] {
                                own[i].0 += 1;
                                own[i].1 += d;
                                own[j].0 += 1;
                                own[j].1 += d;
                            } else {
                                let p = other[i].entry(r[j]).or_insert((0, 0.0));
                                p.0 += 1;
                                p.1 += d;
                                let q = other[j].entry(r[i]).or_insert((0, 0.0));
                                q.0 += 1;
                                q.1 += d;
                            }
                        }
                    }
                }
            }
            kvals += 1;
            if sil {
                let r = roots(&parent);
                let v = silhouette(&r, &own, &other);
                values.push((k, v));
                if best.is_none_or(|x| x < v) {
                    best = Some(v);
                    if k != best_k {
                        snapshot = Some(r);
                        best_k = k;
                    }
                }
            }
        }
        let ind = snapshot.unwrap_or_else(|| roots(&parent));
        let mut seen: Vec<usize> = Vec::new();
        let labels = ind
            .iter()
            .map(|&r| {
                seen.iter().position(|&x| x == r).unwrap_or_else(|| {
                    seen.push(r);
                    seen.len() - 1
                })
            })
            .collect();
        Ok(Some(Linked {
            k: k.max(best_k),
            labels,
            kvals,
            values,
        }))
    }
}

/// `type=` (:11790-11803); centroid and median are `#if 0` in giac.
fn linkage(t: &str) -> Option<Linkage> {
    Some(match t {
        "single" => Linkage::Single,
        "complete" => Linkage::Complete,
        "average" => Linkage::Average,
        "weighted" => Linkage::Weighted,
        "ward" => Linkage::Ward,
        _ => return None,
    })
}

/// A distance function: a builtin symbol or a `Function`.
fn is_function(f: &Expr) -> bool {
    f.as_symbol().is_some_and(Symbol::is_builtin)
        || f.head().and_then(Expr::as_symbol) == Some(Symbol::FUNCTION)
}

/// `Cluster[data, options…]`.
#[allow(
    clippy::too_many_lines,
    reason = "giac's checks and options in its order"
)]
pub(in crate::native) fn cluster(ev: &mut Eval<'_>, e: &Expr) -> Native {
    let Some((d, options)) = e.args().split_first() else {
        return Ok(None);
    };
    let Some(data) = data(ev, d) else {
        return Ok(None);
    };
    if data.items.is_empty() {
        return refuse(ev, e, NO_DATA);
    }
    let np = wide(data.items.len());
    let (mut out, mut k, mut min_k, mut max_k) = (Output::List, 0, 0, 16);
    let mut meth = Linkage::Single;
    let mut function: Option<&Expr> = None;
    let mut indf = 0;
    let mut stays = false;
    for o in options {
        let step = if let Some((lh, rh)) = equation(o) {
            if is(lh, Symbol::COUNT, "count") {
                count(rh).map(|c| match c {
                    Count::K(v) => k = v,
                    Count::Range(a, b) => (min_k, max_k) = (a, b),
                })
            } else if is(lh, Symbol::COUNT_GREATER, "count_sup") {
                bound(rh).map(|v| max_k = v.min(np - 1))
            } else if is(lh, Symbol::COUNT_LESS, "count_inf") {
                bound(rh).map(|v| min_k = v.max(2))
            } else if is(lh, Symbol::FIRST_POSITION, "index") {
                index(rh, hclust_index, SILHOUETTE).map(|v| indf = v)
            } else if drawing(lh) || named(lh, "labels") {
                stays = true;
                Ok(())
            } else if is(lh, Symbol::TYPE, "type") {
                match rh.as_str() {
                    None => Err(STRING),
                    Some(t) => linkage(t).ok_or(LINKAGE).map(|l| meth = l),
                }
            } else if named(lh, "distance") {
                if is_function(rh) {
                    function = Some(rh);
                    Ok(())
                } else {
                    Err(DISTANCE)
                }
            } else if named(lh, "output") {
                output(rh, true).map(|v| out = v)
            } else {
                Err(OPTION)
            }
        } else if let Some(v) = integer(o) {
            if v <= 0 {
                Err(POSITIVE)
            } else {
                k = v;
                Ok(())
            }
        } else if o.as_symbol().is_some_and(Symbol::is_builtin) {
            function = Some(o);
            Ok(())
        } else {
            Err(OPTION)
        };
        if let Err(why) = step {
            return refuse(ev, e, why);
        }
    }
    if stays || out == Output::Out {
        return Ok(None);
    }
    if out == Output::Index && indf == 0 {
        return refuse(ev, e, NO_INDEX);
    }
    if max_k < min_k {
        return refuse(ev, e, K_RANGE);
    }
    if max_k == min_k {
        k = min_k;
    }
    let m = data.items.len();
    let terms = ev.context().limits().terms;
    let pairs = m.checked_mul(m - 1).map(|p| p / 2);
    if pairs.is_none_or(|p| p > terms) || m.checked_mul(data.n).is_none_or(|p| p > terms) {
        return refuse(ev, e, TOO_MANY);
    }
    let pairs = pairs.unwrap_or(0);
    let to = |v: i64| usize::try_from(v).unwrap_or(0);
    // The silhouette's steps: the pairs per `k` from `N − 1` (or `max_k`) down.
    let steps = to(max_k).min(m - 1);
    if k <= 0 && indf & SILHOUETTE != 0 && pairs.checked_mul(steps).is_none_or(|p| p > terms) {
        return refuse(ev, e, TOO_MANY);
    }
    // giac's distance for a list of strings (:11741).
    let strings = data.items.iter().all(|x| x.as_str().is_some());
    if function.is_none() && strings {
        let longest = data
            .items
            .iter()
            .filter_map(Expr::as_str)
            .map(str::len)
            .max();
        let longest = longest.unwrap_or(0);
        if pairs.checked_mul(longest).is_none_or(|p| p > terms)
            || longest.checked_mul(longest).is_none_or(|p| p > terms)
        {
            return refuse(ev, e, TOO_MANY);
        }
    }
    let mut hc = Hclust {
        n: m,
        distance: match function {
            Some(f) => Distance::Function(f, data.items),
            None if strings => Distance::Levenshtein(data.items),
            None => Distance::Squared(&data.points),
        },
        cache: vec![-1.0; pairs],
    };
    let linked = match hc.linkage(ev, meth, to(k), (to(min_k), to(max_k)), indf)? {
        Some(l) if l.k > 0 => l,
        _ => return refuse(ev, e, BAD_DISTANCE),
    };
    if out != Output::Index {
        return Ok(Some(answer(out, linked.k, &linked.labels, data.items)));
    }
    // `mtran([kvals, values])` (optimization.h:640): mcclain-rao and dunn
    // compute no value.
    if linked.values.is_empty() || linked.values.len() != linked.kvals {
        return refuse(ev, e, TRANSPOSE);
    }
    Ok(Some(list(
        linked
            .values
            .iter()
            .map(|&(k, v)| {
                let f = Float::from_f64(v).expect("a silhouette in [-1, 1] or DBL_MIN");
                list(vec![
                    Expr::from(Integer::from(k as u64)),
                    Expr::from(Ball::exact(f)),
                ])
            })
            .collect(),
    )))
}
