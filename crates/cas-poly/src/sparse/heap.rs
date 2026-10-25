//! The heap of term products through which sparse products
//! (`product.rs`) and divisions merge their rows, Johnson's method
//! \[Joh74\] in Monagan and Pearce's form \[MP07, MP11\], and the division:
//! the terms of a sum of products come out of the heap in decreasing order,
//! so the result needs no sorting, and each product's monomial is stored
//! once, in a slot reused by its successor.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of a = q b + r, i and j for term indices"
)]

use std::cmp::Ordering;

use cas_num::algebra::EuclideanDomain;

use super::layout::Layout;
use super::{Sparse, SparseRing, Terms};
use crate::cancel::{Cancel, Cancelled};

/// No slot: the end of a chain.
const NONE: usize = usize::MAX;

/// A max-heap of term products `(i, j)` keyed by their packed monomials,
/// which live in slots: a popped slot is refilled by the product that
/// follows it. A slot put in by [`Heap::insert_chained`] can instead hang
/// in the chain of a slot of the same monomial \[MP11\].
pub(super) struct Heap {
    words: usize,
    /// Slot `s`: its monomial.
    monos: Vec<u64>,
    /// Slot `s`: its product.
    pairs: Vec<(usize, usize)>,
    /// Slot `s`: the next slot of its chain, or `NONE`.
    next: Vec<usize>,
    /// The slots in the heap, a binary max-heap by monomial.
    order: Vec<usize>,
}

impl Heap {
    pub(super) fn new(words: usize) -> Heap {
        Heap {
            words,
            monos: Vec::new(),
            pairs: Vec::new(),
            next: Vec::new(),
            order: Vec::new(),
        }
    }

    #[inline]
    fn mono(&self, slot: usize) -> &[u64] {
        &self.monos[slot * self.words..(slot + 1) * self.words]
    }

    #[inline]
    pub(super) fn mono_mut(&mut self, slot: usize) -> &mut [u64] {
        &mut self.monos[slot * self.words..(slot + 1) * self.words]
    }

    /// A new slot for `(i, j)`; its monomial is written by the caller.
    #[inline]
    pub(super) fn slot(&mut self, i: usize, j: usize) -> usize {
        self.pairs.push((i, j));
        self.next.push(NONE);
        self.monos.extend(std::iter::repeat_n(0, self.words));
        self.pairs.len() - 1
    }

    /// The product in `slot`.
    #[inline]
    pub(super) fn pair(&self, slot: usize) -> (usize, usize) {
        self.pairs[slot]
    }

    #[inline]
    pub(super) fn set_pair(&mut self, slot: usize, pair: (usize, usize)) {
        self.pairs[slot] = pair;
    }

    /// Puts `slot`, its monomial written, into the heap.
    #[inline]
    fn insert(&mut self, slot: usize) {
        let mut hole = self.order.len();
        self.order.push(slot);
        while hole > 0 {
            let parent = (hole - 1) / 2;
            if self.mono(slot) <= self.mono(self.order[parent]) {
                break;
            }
            self.order[hole] = self.order[parent];
            hole = parent;
        }
        self.order[hole] = slot;
    }

    /// Puts `slot`, its monomial written, into the heap, or into the chain
    /// of a slot of the same monomial \[MP11\]: the top's, which it is
    /// compared with first, or one it meets on its way up from the bottom.
    #[inline]
    pub(super) fn insert_chained(&mut self, slot: usize) {
        let Some(&top) = self.order.first() else {
            self.order.push(slot);
            return;
        };
        if self.mono(slot) == self.mono(top) {
            self.chain(top, slot);
            return;
        }
        // The way up ends at the first slot not below `slot`, which takes
        // it in its chain when their monomials are equal.
        let mut place = self.order.len();
        while place > 0 {
            let parent = (place - 1) / 2;
            match self.mono(slot).cmp(self.mono(self.order[parent])) {
                Ordering::Greater => place = parent,
                Ordering::Equal => {
                    self.chain(self.order[parent], slot);
                    return;
                }
                Ordering::Less => break,
            }
        }
        // The slots on the way, from the bottom up to `place`, move down.
        let mut hole = self.order.len();
        self.order.push(slot);
        while hole > place {
            let parent = (hole - 1) / 2;
            self.order[hole] = self.order[parent];
            hole = parent;
        }
        self.order[place] = slot;
    }

    /// Hangs `slot`, in no chain, in `head`'s.
    #[inline]
    fn chain(&mut self, head: usize, slot: usize) {
        debug_assert_eq!(self.next[slot], NONE);
        self.next[slot] = self.next[head];
        self.next[head] = slot;
    }

    /// The slot after `slot` in its chain, which leaves the chain.
    #[inline]
    pub(super) fn take_next(&mut self, slot: usize) -> Option<usize> {
        let next = std::mem::replace(&mut self.next[slot], NONE);
        (next != NONE).then_some(next)
    }

    /// The largest monomial in the heap.
    #[inline]
    pub(super) fn top(&self) -> Option<&[u64]> {
        self.order.first().map(|&s| self.mono(s))
    }

    /// Takes the slot of the largest monomial out of the heap, its chain
    /// with it. The hole it leaves goes down along the larger children, one
    /// comparison a level, and the last slot fills it, rising as far as it
    /// must: fewer comparisons than sifting the last slot down from the top,
    /// two a level, for monomials kept outside the heap \[MP11\].
    #[inline]
    pub(super) fn pop(&mut self) -> usize {
        let top = self.order[0];
        let last = self.order.pop().expect("a slot in the heap");
        let n = self.order.len();
        if n == 0 {
            return top;
        }
        let mut hole = 0;
        loop {
            let left = 2 * hole + 1;
            if left >= n {
                break;
            }
            let right = left + 1;
            let child = if right < n && self.mono(self.order[right]) > self.mono(self.order[left]) {
                right
            } else {
                left
            };
            self.order[hole] = self.order[child];
            hole = child;
        }
        while hole > 0 {
            let parent = (hole - 1) / 2;
            if self.mono(last) <= self.mono(self.order[parent]) {
                break;
            }
            self.order[hole] = self.order[parent];
            hole = parent;
        }
        self.order[hole] = last;
        top
    }
}

/// A packed exponent overflowed its field.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct Overflow;

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
/// An exact quotient has the degrees of `a` less those of `b` in every
/// variable, so an exact division stops at the first quotient term past
/// them (ours, CAS-137): a failed trial division then never outgrows
/// `a`'s exponents.
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
    let bound = if exact && !a.is_zero() {
        let Some(bound) = a
            .degrees()
            .iter()
            .zip(b.degrees())
            .map(|(x, y)| x.checked_sub(y))
            .collect::<Option<Vec<u64>>>()
        else {
            return Ok(None);
        };
        Some(bound)
    } else {
        None
    };
    let mut bits = a.layout.bits().max(b.layout.bits());
    loop {
        let layout = Layout::new(ring.nvars(), bits);
        let (pa, pb) = (a.repacked(layout), b.repacked(layout));
        let packed = bound.as_ref().map(|e| {
            let mut mono = vec![0; layout.words()];
            layout.pack(e, &mut mono);
            mono
        });
        if let Ok(result) = divide_in(ring, &pa, &pb, layout, packed.as_deref(), cancel)? {
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
/// term of `q`; another goes to `r`. An exact division has the `bound` of
/// its quotient's exponents, and a remainder term or a quotient term past
/// the bound stops it with `None`. `Overflow` when a product does not fit
/// the fields.
fn divide_in<R: EuclideanDomain>(
    ring: &SparseRing<R>,
    a: &Sparse<R::Elem>,
    b: &Sparse<R::Elem>,
    layout: Layout,
    bound: Option<&[u64]>,
    cancel: &dyn Cancel,
) -> Result<Divided<R::Elem>, Cancelled> {
    let exact = bound.is_some();
    let base = ring.base();
    let words = layout.words();
    let lead = &b.coeffs[0];
    let mut heap = Heap::new(words);
    let mut q: Terms<R::Elem> = Terms::new(layout);
    let mut r: Terms<R::Elem> = Terms::new(layout);
    let mut mono = vec![0u64; words];
    let mut qmono = vec![0u64; words];
    let mut spare = vec![0u64; words];
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
            if bound.is_some_and(|bound| !layout.divides(&qmono, bound, &mut spare)) {
                return Ok(Ok(None));
            }
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
    /// the remainder appears, or a term of the quotient past the degrees of
    /// `a` less those of `b`: the test gcd algorithms make. Polls `cancel`.
    ///
    /// # Errors
    ///
    /// [`Cancelled`] when `cancel` says so.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    pub fn div_exact_with(
        &self,
        a: &Sparse<R::Elem>,
        b: &Sparse<R::Elem>,
        cancel: &dyn Cancel,
    ) -> Result<Option<Sparse<R::Elem>>, Cancelled> {
        Ok(divide(self, a, b, true, cancel)?.map(|(q, _)| q))
    }

    /// As [`SparseRing::div_exact_with`], never cancelled.
    ///
    /// # Panics
    ///
    /// If `b` is zero.
    #[must_use]
    pub fn div_exact(&self, a: &Sparse<R::Elem>, b: &Sparse<R::Elem>) -> Option<Sparse<R::Elem>> {
        match self.div_exact_with(a, b, &crate::cancel::Never) {
            Ok(q) => q,
            Err(Cancelled) => unreachable!("never cancelled"),
        }
    }
}
