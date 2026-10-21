//! Products and division of sparse polynomials through a heap of term
//! products, Johnson's method \[Joh74\] in Monagan and Pearce's form
//! \[MP07, MP11\]: the terms of a sum of products come out of the heap in
//! decreasing order, so the result needs no sorting, and each product's
//! monomial is stored once, in a slot reused by its successor.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of a = q b + r, i and j for term indices"
)]

use cas_num::algebra::{EuclideanDomain, Ring};

use super::layout::Layout;
use super::{Sparse, SparseRing, Terms};
use crate::cancel::{Cancel, Cancelled};

/// A max-heap of term products `(i, j)` keyed by their packed monomials,
/// which live in slots: a popped slot is refilled by the product that
/// follows it.
struct Heap {
    words: usize,
    /// Slot `s`: its monomial.
    monos: Vec<u64>,
    /// Slot `s`: its product.
    pairs: Vec<(usize, usize)>,
    /// The slots in the heap, a binary max-heap by monomial.
    order: Vec<usize>,
}

impl Heap {
    fn new(words: usize) -> Heap {
        Heap {
            words,
            monos: Vec::new(),
            pairs: Vec::new(),
            order: Vec::new(),
        }
    }

    fn mono(&self, slot: usize) -> &[u64] {
        &self.monos[slot * self.words..(slot + 1) * self.words]
    }

    fn mono_mut(&mut self, slot: usize) -> &mut [u64] {
        &mut self.monos[slot * self.words..(slot + 1) * self.words]
    }

    /// A new slot for `(i, j)`; its monomial is written by the caller.
    fn slot(&mut self, i: usize, j: usize) -> usize {
        self.pairs.push((i, j));
        self.monos.extend(std::iter::repeat_n(0, self.words));
        self.pairs.len() - 1
    }

    fn greater(&self, a: usize, b: usize) -> bool {
        self.mono(self.order[a]) > self.mono(self.order[b])
    }

    /// Puts `slot`, its monomial written, into the heap.
    fn insert(&mut self, slot: usize) {
        self.order.push(slot);
        let mut k = self.order.len() - 1;
        while k > 0 {
            let parent = (k - 1) / 2;
            if !self.greater(k, parent) {
                break;
            }
            self.order.swap(k, parent);
            k = parent;
        }
    }

    /// The largest monomial in the heap.
    fn top(&self) -> Option<&[u64]> {
        self.order.first().map(|&s| self.mono(s))
    }

    /// Takes the slot of the largest monomial out of the heap.
    fn pop(&mut self) -> usize {
        let top = self.order.swap_remove(0);
        let n = self.order.len();
        let mut k = 0;
        loop {
            let (l, r) = (2 * k + 1, 2 * k + 2);
            let mut largest = k;
            if l < n && self.greater(l, largest) {
                largest = l;
            }
            if r < n && self.greater(r, largest) {
                largest = r;
            }
            if largest == k {
                break;
            }
            self.order.swap(k, largest);
            k = largest;
        }
        top
    }
}

/// A packed exponent overflowed its field.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct Overflow;

/// The terms of `a b` in `layout`, where both are packed: the rows `a_i b`
/// merged through a heap that starts row `i + 1` when it takes `(i, 0)`
/// out, so that it holds at most one product a row. `Overflow` when a
/// product's exponent does not fit.
pub(super) fn product_terms<R: Ring>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    layout: Layout,
    cancel: &dyn Cancel,
) -> Result<Result<Sparse<R::Elem>, Overflow>, Cancelled> {
    debug_assert!(a.layout == layout && b.layout == layout);
    let base = ring.base();
    let words = layout.words();
    let mut heap = Heap::new(words);
    let mut out: Terms<R::Elem> = Terms::new(layout);
    let mut mono = vec![0u64; words];
    let first = heap.slot(0, 0);
    let mut overflow = layout.add(a.mono(0), b.mono(0), heap.mono_mut(first));
    heap.insert(first);
    let mut rounds = 0u32;
    while let Some(top) = heap.top() {
        rounds = rounds.wrapping_add(1);
        if rounds % 1024 == 0 && cancel.cancelled() {
            return Err(Cancelled);
        }
        mono.copy_from_slice(top);
        let mut c = base.zero();
        while heap.top() == Some(&mono[..]) {
            let slot = heap.pop();
            let (i, j) = heap.pairs[slot];
            c = base.add(&c, &base.mul(&a.coeffs[i], &b.coeffs[j]));
            if j == 0 && i + 1 < a.coeffs.len() {
                let next = heap.slot(i + 1, 0);
                overflow |= layout.add(a.mono(i + 1), b.mono(0), heap.mono_mut(next));
                heap.insert(next);
            }
            if j + 1 < b.coeffs.len() {
                heap.pairs[slot] = (i, j + 1);
                let (x, y) = (a.mono(i), b.mono(j + 1));
                overflow |= layout.add(x, y, heap.mono_mut(slot));
                heap.insert(slot);
            }
        }
        if overflow {
            return Ok(Err(Overflow));
        }
        if !base.is_zero(&c) {
            out.push(c, &mono);
        }
    }
    Ok(Ok(out.done()))
}

/// `a b`. Both are packed with fields wide enough for the sum of their
/// largest exponents, so no product overflows; were one to, the guards
/// would tell, and the product starts over in wider fields.
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
    let mut max = a.max_exponent().saturating_add(b.max_exponent());
    loop {
        let layout = Layout::holding(ring.nvars(), max);
        let (pa, pb) = (a.repacked(layout), b.repacked(layout));
        match product_terms(ring, &pa, &pb, layout, cancel)? {
            Ok(p) => return Ok(p),
            Err(Overflow) => max = layout.max_exponent().saturating_mul(2).saturating_add(1),
        }
    }
}

/// The quotient and the remainder of a division.
pub(super) type QuoRem<E> = (Sparse<E>, Sparse<E>);

/// A division's quotient and remainder, `None` when an exact division
/// fails, or `Overflow`.
type Divided<E> = Result<Option<QuoRem<E>>, Overflow>;

/// `(q, r)` with `a = q b + r`, as [`divide_in`] computes it. In
/// lexicographic order the exponents of later variables can grow at each
/// reduction (`x^10` over `x - y^200` leaves `y^2000`), so no bound is
/// known before: the division starts in the wider of the two layouts, and
/// when the guards report an overflow it starts over in fields twice as
/// wide.
///
/// # Panics
///
/// If `b` is zero, or an exponent reaches `2^63`.
pub(super) fn divide<R: EuclideanDomain>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    exact: bool,
    cancel: &dyn Cancel,
) -> Result<Option<QuoRem<R::Elem>>, Cancelled> {
    assert!(!b.is_zero(), "attempt to divide by zero");
    let mut bits = a.layout.bits().max(b.layout.bits());
    loop {
        let layout = Layout::new(ring.nvars(), bits);
        let (pa, pb) = (a.repacked(layout), b.repacked(layout));
        if let Ok(result) = divide_in(ring, &pa, &pb, layout, exact, cancel)? {
            return Ok(result);
        }
        assert!(bits < 64, "an exponent of 2^63 or more");
        bits = (2 * bits).min(64);
    }
}

/// `(q, r)` with `a = q b + r`, both packed in `layout`: the terms of `a −
/// q b` come out of a heap in decreasing order, `a`'s and those of `−qₖ bⱼ`
/// for `j ≥ 1` (`qₖ b₀` cancels the term `qₖ` was made from) \[MP07\]. A
/// term that the leading term of `b` divides, coefficient included, makes a
/// term of `q`; another goes to `r`, or, when `exact`, stops the division
/// with `None`. `Overflow` when a product does not fit the fields.
fn divide_in<R: EuclideanDomain>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    layout: Layout,
    exact: bool,
    cancel: &dyn Cancel,
) -> Result<Divided<R::Elem>, Cancelled> {
    let base = ring.base();
    let words = layout.words();
    let lead = &b.coeffs[0];
    let mut heap = Heap::new(words);
    let mut q: Terms<R::Elem> = Terms::new(layout);
    let mut r: Terms<R::Elem> = Terms::new(layout);
    let mut mono = vec![0u64; words];
    let mut qmono = vec![0u64; words];
    let mut i = 0;
    let mut rounds = 0u32;
    loop {
        rounds = rounds.wrapping_add(1);
        if rounds % 1024 == 0 && cancel.cancelled() {
            return Err(Cancelled);
        }
        // The largest monomial left: a's next term or the heap's top.
        match (a.coeffs.get(i), heap.top()) {
            (None, None) => break,
            (Some(_), Some(top)) if a.mono(i) < top => mono.copy_from_slice(top),
            (Some(_), _) => mono.copy_from_slice(a.mono(i)),
            (None, Some(top)) => mono.copy_from_slice(top),
        }
        let mut c = base.zero();
        if i < a.coeffs.len() && a.mono(i) == &mono[..] {
            c = a.coeffs[i].clone();
            i += 1;
        }
        while heap.top() == Some(&mono[..]) {
            let slot = heap.pop();
            let (k, j) = heap.pairs[slot];
            c = base.sub(&c, &base.mul(&q.coeffs[k], &b.coeffs[j]));
            if j + 1 < b.coeffs.len() {
                heap.pairs[slot] = (k, j + 1);
                let qk = &q.exps[k * words..(k + 1) * words];
                if layout.add(qk, b.mono(j + 1), heap.mono_mut(slot)) {
                    return Ok(Err(Overflow));
                }
                heap.insert(slot);
            }
        }
        if base.is_zero(&c) {
            continue;
        }
        let (t, rest) = if layout.divides(b.mono(0), &mono, &mut qmono) {
            base.div_rem(&c, lead)
        } else {
            (base.zero(), c.clone())
        };
        if base.is_zero(&rest) {
            q.push(t, &qmono);
            if b.coeffs.len() > 1 {
                let k = q.coeffs.len() - 1;
                let slot = heap.slot(k, 1);
                if layout.add(&qmono, b.mono(1), heap.mono_mut(slot)) {
                    return Ok(Err(Overflow));
                }
                heap.insert(slot);
            }
        } else {
            if exact {
                return Ok(Ok(None));
            }
            r.push(c, &mono);
        }
    }
    Ok(Ok(Some((q.done(), r.done()))))
}

impl<R: EuclideanDomain> SparseRing<R> {
    /// `(q, r)` with `a = q b + r`, no term of `r` a multiple of the leading
    /// term of `b`: over a field, the division with remainder by one
    /// polynomial in lexicographic order; over ℤ, a term whose coefficient
    /// `b`'s leading coefficient does not divide goes to `r`. By Monagan and
    /// Pearce's heap \[MP07\], polling `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    pub fn div_rem_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<QuoRem<R::Elem>, Cancelled> {
        Ok(divide(self, a, b, false, cancel)?.expect("a division with remainder"))
    }

    /// As [`SparseRing::div_rem_with`], never cancelled.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[must_use]
    pub fn div_rem(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> QuoRem<R::Elem> {
        match self.div_rem_with(a, b, &crate::cancel::Never) {
            Ok(qr) => qr,
            Err(Cancelled) => unreachable!("never cancelled"),
        }
    }

    /// `a / b` when `b` divides `a`; `None` otherwise, as soon as a term of
    /// the remainder appears: the test gcd algorithms make.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[must_use]
    pub fn div_exact(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> Option<Sparse<R::Elem>> {
        match divide(self, a, b, true, &crate::cancel::Never) {
            Ok(result) => result.map(|(q, _)| q),
            Err(Cancelled) => unreachable!("never cancelled"),
        }
    }
}
