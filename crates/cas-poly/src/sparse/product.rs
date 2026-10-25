//! Products of sparse polynomials (CAS-136, CAS-155), by one of two kernels,
//! accumulating in one of three ways.
//!
//! The kernels:
//! - **Arrays**, when the product's degree box is small next to its products
//!   of terms. The products whose degrees in the first variable add up to
//!   `d` add into an array over the box of the other variables, a chunk, at
//!   an index linear in their exponents, with no comparison; the chunk is
//!   read back in lexicographic order, for each `d` from the largest down.
//! - Otherwise **Johnson's heap** \[Joh74, MP07\] with lazy rows, which takes
//!   every product of the top monomial out before it puts their successors
//!   in, and chains the products of equal monomials \[MP11\].
//!
//! The accumulators \[MP11\]: products of coefficients with small integer
//! preimages (cas-num's `Ring::small_preimage`: integers, integral
//! rationals, word residues) add up in three machine words, mapped into the
//! ring once a term of the result; the others in the ring.

#![allow(
    clippy::many_single_char_names,
    reason = "a and b for the factors, i and j for their terms"
)]

use std::ops::Range;

use cas_num::algebra::Ring;
use cas_num::{Integer, Natural};

use super::heap::{Heap, Overflow};
use super::layout::Layout;
use super::{Sparse, SparseRing, Terms};
use crate::cancel::{Cancel, Cancelled};

/// A signed integer of three words, `hi 2^128 + lo`, `lo` unsigned: a sum
/// of products of two `i64`s, each at most `2^126` in absolute value, so
/// that fewer than `2^64` of them stay below `2^190`.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub(super) struct Wide {
    /// The low 128 bits as two words, so that a `Wide` takes 24 bytes.
    lo: [u64; 2],
    hi: i64,
}

impl Wide {
    /// `self += a b`.
    #[inline]
    #[allow(
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation,
        reason = "the two's complement of the product, and the words of lo"
    )]
    pub(super) fn add_product(&mut self, a: i64, b: i64) {
        let p = i128::from(a) * i128::from(b);
        let lo = u128::from(self.lo[0]) | u128::from(self.lo[1]) << 64;
        // `p as u128` is `p + 2^128` when `p` is negative.
        let (lo, carry) = lo.overflowing_add(p as u128);
        self.lo = [lo as u64, (lo >> 64) as u64];
        self.hi += i64::from(carry) - i64::from(p < 0);
    }

    fn is_zero(&self) -> bool {
        *self == Wide::default()
    }

    #[allow(
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss,
        reason = "the words of a two's complement integer"
    )]
    pub(super) fn to_integer(self) -> Integer {
        let lo = (u128::from(self.lo[0]) | u128::from(self.lo[1]) << 64) as i128;
        // An `i128` when `hi` only extends the sign of `lo`.
        if self.hi == if lo < 0 { -1 } else { 0 } {
            return Integer::from(lo);
        }
        let negative = self.hi < 0;
        let mut limbs = [self.lo[0], self.lo[1], self.hi as u64];
        if negative {
            // The magnitude: the two's complement of the three words.
            let mut carry = true;
            for limb in &mut limbs {
                (*limb, carry) = (!*limb).overflowing_add(u64::from(carry));
            }
        }
        Integer::from_parts(negative, Natural::from_limbs(&limbs))
    }
}

/// How the products of the terms of the factors add up.
trait Accumulate<R: Ring> {
    type Sum;

    fn zero(&self) -> Self::Sum;

    /// `sum += aᵢ bⱼ`.
    fn add_product(&self, sum: &mut Self::Sum, i: usize, j: usize);

    /// The coefficient of `sum`, `None` when it is zero; `sum` is left zero.
    fn take(&self, sum: &mut Self::Sum) -> Option<R::Elem>;
}

/// In the ring.
struct General<'a, R: Ring> {
    ring: &'a R,
    a: &'a [R::Elem],
    b: &'a [R::Elem],
}

impl<R: Ring> Accumulate<R> for General<'_, R> {
    type Sum = R::Elem;

    fn zero(&self) -> R::Elem {
        self.ring.zero()
    }

    #[inline]
    fn add_product(&self, sum: &mut R::Elem, i: usize, j: usize) {
        *sum = self.ring.add(sum, &self.ring.mul(&self.a[i], &self.b[j]));
    }

    fn take(&self, sum: &mut R::Elem) -> Option<R::Elem> {
        let c = std::mem::replace(sum, self.ring.zero());
        (!self.ring.is_zero(&c)).then_some(c)
    }
}

/// In words, from the small preimages of every coefficient.
struct Words<'a, R> {
    ring: &'a R,
    a: Vec<i64>,
    b: Vec<i64>,
}

impl<'a, R> Words<'a, R> {
    /// When every coefficient has a small preimage.
    fn new(ring: &'a R, a: &[Option<i64>], b: &[Option<i64>]) -> Option<Words<'a, R>> {
        Some(Words {
            ring,
            a: a.iter().copied().collect::<Option<_>>()?,
            b: b.iter().copied().collect::<Option<_>>()?,
        })
    }
}

impl<R: Ring> Accumulate<R> for Words<'_, R> {
    type Sum = Wide;

    fn zero(&self) -> Wide {
        Wide::default()
    }

    #[inline]
    fn add_product(&self, sum: &mut Wide, i: usize, j: usize) {
        sum.add_product(self.a[i], self.b[j]);
    }

    fn take(&self, sum: &mut Wide) -> Option<R::Elem> {
        let s = std::mem::take(sum);
        if s.is_zero() {
            return None;
        }
        let c = self.ring.elem_integer(&s.to_integer());
        (!self.ring.is_zero(&c)).then_some(c)
    }
}

/// Products of two small preimages in words, the others in the ring.
struct Mixed<'a, R: Ring> {
    ring: &'a R,
    a: &'a [R::Elem],
    b: &'a [R::Elem],
    small_a: Vec<Option<i64>>,
    small_b: Vec<Option<i64>>,
}

impl<R: Ring> Accumulate<R> for Mixed<'_, R> {
    type Sum = (Wide, R::Elem);

    fn zero(&self) -> (Wide, R::Elem) {
        (Wide::default(), self.ring.zero())
    }

    #[inline]
    fn add_product(&self, (words, rest): &mut (Wide, R::Elem), i: usize, j: usize) {
        match (self.small_a[i], self.small_b[j]) {
            (Some(x), Some(y)) => words.add_product(x, y),
            _ => *rest = self.ring.add(rest, &self.ring.mul(&self.a[i], &self.b[j])),
        }
    }

    fn take(&self, (words, rest): &mut (Wide, R::Elem)) -> Option<R::Elem> {
        let w = std::mem::take(words);
        let mut c = std::mem::replace(rest, self.ring.zero());
        if !w.is_zero() {
            c = self.ring.add(&c, &self.ring.elem_integer(&w.to_integer()));
        }
        (!self.ring.is_zero(&c)).then_some(c)
    }
}

/// The largest chunk of the arrays, in sums.
const MAX_CHUNK: usize = 1 << 18;

/// The arrays take a product when the sums they read back and the parts
/// they look up are at most this many a product of terms.
const DENSITY: u128 = 4;

/// How often the arrays poll `Cancel`, in products of terms.
const POLL: usize = 1 << 16;

/// The terms of `p` by their degree in `x₀`, from the largest: ranges of
/// terms, which lexicographic order keeps together.
fn parts<E>(p: &Sparse<E>) -> Vec<(u64, Range<usize>)> {
    let mut parts: Vec<(u64, Range<usize>)> = Vec::new();
    for k in 0..p.term_count() {
        let d = p.layout.get(p.mono(k), 0);
        match parts.last_mut() {
            Some((e, range)) if *e == d => range.end = k + 1,
            _ => parts.push((d, k..k + 1)),
        }
    }
    parts
}

/// The arrays' layout of a product `a b` in `n` variables: chunks over the
/// box of `x₁ … xₙ₋₁`, one for each degree in `x₀`.
struct Chunks {
    /// The stride of `xᵥ` in a chunk, for `v ≥ 1`, the first the largest so
    /// that indices order monomials as lexicographic order does;
    /// `strides[0]` is 0.
    strides: Vec<usize>,
    /// The sums of a chunk.
    size: usize,
    /// The degrees of `a` and `b` in `x₀`.
    first: (u64, u64),
    /// The largest exponent of the box.
    max: u64,
    /// The terms of `a` by their degree in `x₀`.
    parts: Vec<(u64, Range<usize>)>,
}

impl Chunks {
    /// The chunks of `a b`; `None` without variables, or when a chunk would
    /// pass `MAX_CHUNK`.
    fn of<E>(a: &Sparse<E>, b: &Sparse<E>) -> Option<Chunks> {
        let n = a.nvars();
        if n == 0 {
            return None;
        }
        let (da, db) = (a.degrees(), b.degrees());
        let mut strides = vec![0; n];
        let mut size = 1usize;
        for v in (1..n).rev() {
            strides[v] = size;
            // Exponents below 2^63: no overflow.
            let width = usize::try_from(da[v] + db[v] + 1).ok()?;
            size = size.checked_mul(width).filter(|&s| s <= MAX_CHUNK)?;
        }
        let max = da.iter().zip(&db).map(|(x, y)| x + y).max().unwrap_or(0);
        Some(Chunks {
            strides,
            size,
            first: (da[0], db[0]),
            max,
            parts: parts(a),
        })
    }

    /// Whether the arrays suit `a b` better than the heap: the sums they
    /// read back and the parts of `a` they look up, a chunk at a time, at
    /// most `DENSITY` a product of terms.
    fn suit<E>(&self, a: &Sparse<E>, b: &Sparse<E>) -> bool {
        let chunks = u128::from(self.first.0) + u128::from(self.first.1) + 1;
        let work = chunks * (self.size as u128 + self.parts.len() as u128);
        work <= DENSITY * (a.term_count() as u128 * b.term_count() as u128)
    }
}

/// An exponent or an index in a chunk.
fn to_usize(x: u64) -> usize {
    usize::try_from(x).expect("an index in memory")
}

/// `a b` through the chunks: for each degree `d` in `x₀`, from the largest
/// down, the products of the parts of `a` and `b` whose degrees add up to
/// `d` add into the chunk, which is then read back, from the largest index
/// touched to the smallest, and left zero.
fn array_product<R: Ring, A: Accumulate<R>>(
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    chunks: &Chunks,
    acc: &A,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    let n = a.nvars();
    let index = |p: &Sparse<R::Elem>, k: usize| -> usize {
        (1..n)
            .map(|v| to_usize(p.layout.get(p.mono(k), v)) * chunks.strides[v])
            .sum()
    };
    let ia: Vec<usize> = (0..a.term_count()).map(|k| index(a, k)).collect();
    let ib: Vec<usize> = (0..b.term_count()).map(|k| index(b, k)).collect();
    let mut of_b = vec![0..0; to_usize(chunks.first.1) + 1];
    for (d, range) in parts(b) {
        of_b[to_usize(d)] = range;
    }
    let mut sums: Vec<A::Sum> = (0..chunks.size).map(|_| acc.zero()).collect();
    let layout = Layout::holding(n, chunks.max);
    let mut out = Terms::new(layout);
    let mut exps = vec![0u64; n];
    let mut mono = vec![0u64; layout.words()];
    let mut since = 0;
    for d in (0..=chunks.first.0 + chunks.first.1).rev() {
        // Within a part the indices decrease with the terms, so the sums
        // touched lie between those of the parts' last terms and of their
        // first.
        let (mut low, mut high) = (usize::MAX, 0);
        for (e, ra) in &chunks.parts {
            let Some(rest) = d.checked_sub(*e) else {
                continue;
            };
            // The parts of `a` come with decreasing degrees: the rest grows.
            let Some(rb) = usize::try_from(rest).ok().and_then(|k| of_b.get(k)) else {
                break;
            };
            if rb.is_empty() {
                continue;
            }
            for i in ra.clone() {
                let x = ia[i];
                for j in rb.clone() {
                    acc.add_product(&mut sums[x + ib[j]], i, j);
                }
                since += rb.len();
                if since >= POLL {
                    since = 0;
                    if cancel.cancelled() {
                        return Err(Cancelled);
                    }
                }
            }
            low = low.min(ia[ra.end - 1] + ib[rb.end - 1]);
            high = high.max(ia[ra.start] + ib[rb.start]);
        }
        if low > high {
            continue;
        }
        since += high - low + 1;
        exps[0] = d;
        for s in (low..=high).rev() {
            let Some(c) = acc.take(&mut sums[s]) else {
                continue;
            };
            let mut rest = s;
            for (e, &stride) in exps.iter_mut().zip(&chunks.strides).skip(1) {
                *e = (rest / stride) as u64;
                rest %= stride;
            }
            layout.pack(&exps, &mut mono);
            out.push(c, &mono);
        }
    }
    Ok(out.done())
}

/// The terms of `a b` through the heap, from the exponents of `na` and `nb`
/// terms packed in `layout`. It starts row `i + 1` when it takes `(i, 0)`
/// out, so that it holds at most one product a row \[MP07\], and it takes
/// every product of the top monomial out before it puts their successors
/// in, which then chain at the top when they are the next monomial
/// \[MP11\]. `Overflow` when a product's exponent does not fit.
fn heap_product<R: Ring, A: Accumulate<R>>(
    (ea, na): (&[u64], usize),
    (eb, nb): (&[u64], usize),
    layout: Layout,
    acc: &A,
    cancel: &dyn Cancel,
) -> Result<Result<Sparse<R::Elem>, Overflow>, Cancelled> {
    let words = layout.words();
    let ma = |i: usize| &ea[i * words..(i + 1) * words];
    let mb = |j: usize| &eb[j * words..(j + 1) * words];
    let mut heap = Heap::new(words);
    let mut out: Terms<R::Elem> = Terms::new(layout);
    let mut mono = vec![0u64; words];
    let mut taken = Vec::new();
    let first = heap.slot(0, 0);
    let mut overflow = layout.add(ma(0), mb(0), heap.mono_mut(first));
    heap.insert_chained(first);
    let mut rounds = 0u32;
    while let Some(top) = heap.top() {
        rounds = rounds.wrapping_add(1);
        if rounds % 1024 == 0 && cancel.cancelled() {
            return Err(Cancelled);
        }
        mono.copy_from_slice(top);
        let mut c = acc.zero();
        while heap.top() == Some(&mono[..]) {
            let mut slot = Some(heap.pop());
            while let Some(s) = slot {
                slot = heap.take_next(s);
                let (i, j) = heap.pair(s);
                acc.add_product(&mut c, i, j);
                taken.push(s);
            }
        }
        for s in taken.drain(..) {
            let (i, j) = heap.pair(s);
            if j == 0 && i + 1 < na {
                let next = heap.slot(i + 1, 0);
                overflow |= layout.add(ma(i + 1), mb(0), heap.mono_mut(next));
                heap.insert_chained(next);
            }
            if j + 1 < nb {
                heap.set_pair(s, (i, j + 1));
                overflow |= layout.add(ma(i), mb(j + 1), heap.mono_mut(s));
                heap.insert_chained(s);
            }
        }
        if overflow {
            return Ok(Err(Overflow));
        }
        if let Some(c) = acc.take(&mut c) {
            out.push(c, &mono);
        }
    }
    Ok(Ok(out.done()))
}

/// How a product goes.
enum Kernel {
    Arrays(Chunks),
    Heap,
}

impl Kernel {
    fn of<E>(a: &Sparse<E>, b: &Sparse<E>) -> Kernel {
        match Chunks::of(a, b) {
            Some(chunks) if chunks.suit(a, b) => Kernel::Arrays(chunks),
            _ => Kernel::Heap,
        }
    }
}

/// `a b` by `kernel`, accumulating in `acc`.
fn product_by<R: Ring, A: Accumulate<R>>(
    nvars: usize,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    kernel: &Kernel,
    acc: &A,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    match kernel {
        Kernel::Arrays(chunks) => array_product(a, b, chunks, acc, cancel),
        Kernel::Heap => {
            // Fields wide enough for the sum of the largest exponents, so
            // that no product overflows; were one to, the guards would
            // tell, and the product would start over in wider fields.
            let mut max = a.max_exponent().saturating_add(b.max_exponent());
            loop {
                let layout = Layout::holding(nvars, max);
                let (ea, eb) = (a.exps_in(layout), b.exps_in(layout));
                let terms = (a.term_count(), b.term_count());
                match heap_product((&ea, terms.0), (&eb, terms.1), layout, acc, cancel)? {
                    Ok(p) => return Ok(p),
                    Err(Overflow) => {
                        max = layout.max_exponent().saturating_mul(2).saturating_add(1);
                    }
                }
            }
        }
    }
}

/// The small preimages of `p`'s coefficients.
fn preimages<R: Ring>(ring: &R, p: &Sparse<R::Elem>) -> Vec<Option<i64>> {
    p.coeffs.iter().map(|c| ring.small_preimage(c)).collect()
}

/// `a b`, by the kernel that suits it: in words when every coefficient has
/// a small preimage, partly when some have, otherwise in the ring.
pub(super) fn mul<R: Ring>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    cancel: &dyn Cancel,
) -> Result<Sparse<R::Elem>, Cancelled> {
    if a.is_zero() || b.is_zero() {
        return Ok(ring.empty());
    }
    // The heap holds a product a row of its first factor: the shorter.
    let (a, b) = if a.term_count() <= b.term_count() {
        (a, b)
    } else {
        (b, a)
    };
    let kernel = Kernel::of(a, b);
    let base = ring.base();
    let (small_a, small_b) = (preimages(base, a), preimages(base, b));
    let n = ring.nvars();
    if let Some(words) = Words::new(base, &small_a, &small_b) {
        product_by(n, a, b, &kernel, &words, cancel)
    } else if small_a.iter().any(Option::is_some) && small_b.iter().any(Option::is_some) {
        let mixed = Mixed {
            ring: base,
            a: &a.coeffs,
            b: &b.coeffs,
            small_a,
            small_b,
        };
        product_by(n, a, b, &kernel, &mixed, cancel)
    } else {
        let general = General {
            ring: base,
            a: &a.coeffs,
            b: &b.coeffs,
        };
        product_by(n, a, b, &kernel, &general, cancel)
    }
}

/// The accumulators, for the tests to choose.
#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(super) enum Accumulator {
    Words,
    Mixed,
    General,
}

/// `a b` through the arrays or the heap, accumulating as told; `None` when
/// the arrays cannot take the product, a chunk past `MAX_CHUNK`, or the
/// words its coefficients.
#[cfg(test)]
pub(super) fn mul_by<R: Ring>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    arrays: bool,
    accumulator: Accumulator,
    cancel: &dyn Cancel,
) -> Option<Result<Sparse<R::Elem>, Cancelled>> {
    if a.is_zero() || b.is_zero() {
        return Some(Ok(ring.empty()));
    }
    let (a, b) = if a.term_count() <= b.term_count() {
        (a, b)
    } else {
        (b, a)
    };
    let kernel = if arrays {
        Kernel::Arrays(Chunks::of(a, b)?)
    } else {
        Kernel::Heap
    };
    let base = ring.base();
    let (small_a, small_b) = (preimages(base, a), preimages(base, b));
    let n = ring.nvars();
    Some(match accumulator {
        Accumulator::Words => {
            let words = Words::new(base, &small_a, &small_b)?;
            product_by(n, a, b, &kernel, &words, cancel)
        }
        Accumulator::Mixed => {
            let mixed = Mixed {
                ring: base,
                a: &a.coeffs,
                b: &b.coeffs,
                small_a,
                small_b,
            };
            product_by(n, a, b, &kernel, &mixed, cancel)
        }
        Accumulator::General => {
            let general = General {
                ring: base,
                a: &a.coeffs,
                b: &b.coeffs,
            };
            product_by(n, a, b, &kernel, &general, cancel)
        }
    })
}

/// Whether `a b` goes through the arrays.
#[cfg(test)]
pub(super) fn through_arrays<E>(a: &Sparse<E>, b: &Sparse<E>) -> bool {
    let (a, b) = if a.term_count() <= b.term_count() {
        (a, b)
    } else {
        (b, a)
    };
    matches!(Kernel::of(a, b), Kernel::Arrays(_))
}

/// The terms of `a b` through the heap, from their exponents packed in
/// `layout`, accumulating in the ring: `Overflow` when a product's exponent
/// does not fit.
#[cfg(test)]
pub(super) fn heap_terms<R: Ring>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    layout: Layout,
) -> Result<Sparse<R::Elem>, Overflow> {
    let general = General {
        ring: ring.base(),
        a: &a.coeffs,
        b: &b.coeffs,
    };
    let (ea, eb) = (a.exps_in(layout), b.exps_in(layout));
    let terms = (a.term_count(), b.term_count());
    crate::cancel::uncancelled(heap_product(
        (&ea, terms.0),
        (&eb, terms.1),
        layout,
        &general,
        &crate::cancel::Never,
    ))
}
