//! Clustering (CAS-630): Xcas's `kmeans` and `cluster` (giac
//! optimization.cc: AS 136 `optra` :9862, `qtran` :10096, `kmns` :10292;
//! `init_cluster_centers` :10622, `kmeans` :10654; `comp_dist` :10765;
//! `hclust` :11102-11443; `clustering_result` :11571, `_kmeans` :11610,
//! `_cluster` :11737), ported from crates/cas-core/proto/clustering/proto.py (R1.7-012).
//!
//! - **Data.** A list of real numbers (1-D points) or of equal-length
//!   lists of real numbers, `Matrix[l]` of such a list too; reals are
//!   rationals, decimals (their midpoints) and numeric constants, in `f64`.
//!   Anything else leaves the call without a warning.
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
//! - **`Cluster[data, options…]`** (`cluster`): agglomerative clustering,
//!   single linkage by Prim's tree (the default), `type="complete"`,
//!   `"average"`, `"weighted"` and `"ward"` by the nearest-neighbour chain
//!   with Lance and Williams' updates, on squared Euclidean distances
//!   (exact on rational points, then `f64`) or `distance=f` (a builtin or a
//!   `Function`, also a builtin positionally), `N[f[p, q]]` per pair, cached
//!   as giac's. The merges sorted by height; `k` given, or giac's largest
//!   gap of heights (at least `count_inf`, the least of `count=a..b`).
//! - **Answers.** The 0-based labels, numbered by first appearance for
//!   `cluster` (`output=list`, the default), or the number of clusters
//!   (`output=count`), exact integers. `output=plot`, `tree` and `part`,
//!   `index=`, `display=`, `color=` (and `labels=` for `cluster`) leave the
//!   call without a warning (graphics and indices, not ported).
//! - **Errors** leave the call with giac's message as the warning: "No
//!   data found", kmeans' "Data must be numeric", "Invalid k-range", "Invalid
//!   output specification", "No index specified for output", "Unknown
//!   linkage type", "Expected a string", "Invalid distance function
//!   specification", "Invalid option", the count parsers' texts, and "Bad
//!   distance function" for a distance not a real `≥ 0` (a symbolic point
//!   among them).
//! - **Limits** before any work (`Limits::terms`): `cluster` the pairs
//!   `m(m − 1)/2` and `m n`; `kmeans` `11 m n Σk` over the `k` tried and `m n
//!   k limit`; past them the call stays with "too many points". `qtran`
//!   stops after `terms/n` steps (giac's has no bound; a fault 2); the loops
//!   poll the interrupt.

#![allow(
    clippy::many_single_char_names,
    reason = "a, c, d, k, m, n: AS 136's and giac's names"
)]

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering as Atomic};

use cas_num::ball::{Float, Round};
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
    } else if is(rh, Symbol::PART, "part") || named(rh, "plot") || (tree && named(rh, "tree")) {
        Ok(Output::Out)
    } else {
        Err(OUTPUT)
    }
}

/// Options that only draw: the call stays without a warning.
fn drawing(lh: &Expr) -> bool {
    named(lh, "display") || named(lh, "color") || is(lh, Symbol::FIRST_POSITION, "index")
}

/// `count_sup=k` and `count_inf=k`: any integer.
fn bound(rh: &Expr) -> Result<i64, &'static str> {
    integer(rh).ok_or(POSITIVE)
}

fn answer(out: Output, k: usize, labels: &[usize]) -> Expr {
    match out {
        Output::Count => Expr::from(Integer::from(k as u64)),
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

/// A run's labels `1 … k` and each cluster's sum of squares, or its fault.
type Clusters = Result<(Vec<usize>, Vec<f64>), Fault>;

/// `kmns` (:10292-10589) on `points` from the centres `c0`: the labels
/// `1 … k` of the points and each cluster's sum of squares.
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
    Ok(Ok((s.ic1[1..].to_vec(), wss[1..].to_vec())))
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
            run(ev, x, None, k, limits)?.map(|(b, _)| (k, b.iter().map(|l| l - 1).collect()))
        );
    }
    let mut w_prev = 0.0;
    let mut best: Vec<usize> = (1..=m).collect();
    let mut b: Vec<usize> = Vec::new();
    let mut rough = false;
    let mut k = min_k;
    while k <= max_k {
        let ind = (k != min_k).then_some(&b[..]);
        let (labels, wss) = match run(ev, x, ind, k, limits)? {
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
    if out == Output::Index {
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
    let (sum, top) = if kk > 0 {
        (kk as u128, kk)
    } else {
        let (l, h) = (lo as u128, hi as u128);
        ((h * (h + 1) - l * l.saturating_sub(1)) / 2, hi)
    };
    let mn = (m as u128) * (data.n as u128);
    let terms = ev.context().limits().terms as u128;
    let maxiter = to(maxiter);
    if mn * sum * STARTS as u128 > terms || mn * (top as u128) * (maxiter as u128) > terms {
        return refuse(ev, e, TOO_MANY);
    }
    let steps = i64::try_from(terms / (data.n as u128)).unwrap_or(i64::MAX);
    ev.volatile();
    match kmeans_core(ev, &points, kk, (lo, hi), (maxiter, steps))? {
        Ok((k, labels)) => Ok(Some(answer(out, k, &labels))),
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

    /// `linkage` (:11262-11443), giac's path without an index: the number
    /// of clusters and the labels by first appearance; `None` fails.
    fn linkage(
        &mut self,
        ev: &mut Eval<'_>,
        meth: Linkage,
        k_given: usize,
        min_k: usize,
    ) -> Result<Option<(usize, Vec<usize>)>, EvalError> {
        let dg = if meth == Linkage::Single {
            self.mst(ev)?
        } else {
            self.nn_chain(ev, meth)?
        };
        let Some(mut dg) = dg else {
            return Ok(None);
        };
        dg.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let mut parent: Vec<usize> = (0..self.n).collect();
        let find = |parent: &[usize], mut i: usize| {
            while parent[i] != i {
                i = parent[i];
            }
            i
        };
        let mut k = self.n;
        let target = if k_given > 0 {
            k_given
        } else {
            min_k.max(self.k_low(&dg))
        };
        for &(_, a, b) in &dg {
            if k == target {
                break;
            }
            k -= 1;
            let (ra, rb) = (find(&parent, a), find(&parent, b));
            if ra != rb {
                parent[rb] = ra;
            }
        }
        let mut roots: Vec<usize> = Vec::new();
        let labels = (0..self.n)
            .map(|i| {
                let r = find(&parent, i);
                roots.iter().position(|&x| x == r).unwrap_or_else(|| {
                    roots.push(r);
                    roots.len() - 1
                })
            })
            .collect();
        Ok(Some((k, labels)))
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
    if out == Output::Index {
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
    let mut hc = Hclust {
        n: m,
        distance: match function {
            Some(f) => Distance::Function(f, data.items),
            None => Distance::Squared(&data.points),
        },
        cache: vec![-1.0; pairs.unwrap_or(0)],
    };
    let to = |v: i64| usize::try_from(v).unwrap_or(0);
    match hc.linkage(ev, meth, to(k), to(min_k))? {
        Some((k, labels)) if k > 0 => Ok(Some(answer(out, k, &labels))),
        _ => refuse(ev, e, BAD_DISTANCE),
    }
}
