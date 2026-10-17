//! Limb storage with inline capacity for small values (layout E, CAS-DOC-2 §4.1).

use crate::Limb;

/// Number of limbs stored without a heap allocation (192 bits).
pub(crate) const INLINE: usize = 3;

/// A growable limb buffer: inline up to [`INLINE`] limbs, `Vec` beyond.
///
/// Values stay heap-allocated when they shrink, so that in-place loops don't
/// bounce between representations; only [`LimbVec::shrink_to_fit`] and `clone`
/// move small values back inline. Code must therefore never compare the
/// representation, only [`LimbVec::as_slice`].
#[derive(Debug)]
pub(crate) enum LimbVec {
    Inline { len: u8, buf: [Limb; INLINE] },
    Heap(Vec<Limb>),
}

impl LimbVec {
    pub(crate) const fn new() -> Self {
        LimbVec::Inline {
            len: 0,
            buf: [0; INLINE],
        }
    }

    /// `n` zero limbs.
    pub(crate) fn zeroed(n: usize) -> Self {
        if n <= INLINE {
            LimbVec::Inline {
                len: inline_len(n),
                buf: [0; INLINE],
            }
        } else {
            LimbVec::Heap(vec![0; n])
        }
    }

    pub(crate) fn from_slice(limbs: &[Limb]) -> Self {
        if limbs.len() <= INLINE {
            let mut buf = [0; INLINE];
            buf[..limbs.len()].copy_from_slice(limbs);
            LimbVec::Inline {
                len: inline_len(limbs.len()),
                buf,
            }
        } else {
            LimbVec::Heap(limbs.to_vec())
        }
    }

    /// The normalised limbs of `hi * 2^64 + lo`, inline. A fixed-size write,
    /// so the conversions from primitives avoid `from_slice`'s `memcpy` call.
    pub(crate) fn from_pair(lo: Limb, hi: Limb) -> Self {
        const { assert!(INLINE >= 2) };
        let mut buf = [0; INLINE];
        buf[0] = lo;
        buf[1] = hi;
        let len = if hi != 0 { 2 } else { u8::from(lo != 0) };
        LimbVec::Inline { len, buf }
    }

    /// Takes ownership of `limbs`, moving them inline if they fit.
    pub(crate) fn from_vec(limbs: Vec<Limb>) -> Self {
        if limbs.len() <= INLINE {
            Self::from_slice(&limbs)
        } else {
            LimbVec::Heap(limbs)
        }
    }

    pub(crate) fn as_slice(&self) -> &[Limb] {
        match self {
            LimbVec::Inline { len, buf } => &buf[..usize::from(*len)],
            LimbVec::Heap(v) => v,
        }
    }

    pub(crate) fn as_mut_slice(&mut self) -> &mut [Limb] {
        match self {
            LimbVec::Inline { len, buf } => &mut buf[..usize::from(*len)],
            LimbVec::Heap(v) => v,
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            LimbVec::Inline { len, .. } => usize::from(*len),
            LimbVec::Heap(v) => v.len(),
        }
    }

    #[cfg(test)]
    pub(crate) fn is_inline(&self) -> bool {
        matches!(self, LimbVec::Inline { .. })
    }

    pub(crate) fn push(&mut self, limb: Limb) {
        match self {
            LimbVec::Inline { len, buf } if usize::from(*len) < INLINE => {
                buf[usize::from(*len)] = limb;
                *len += 1;
            }
            LimbVec::Inline { buf, .. } => {
                let mut v = Vec::with_capacity(2 * INLINE);
                v.extend_from_slice(buf);
                v.push(limb);
                *self = LimbVec::Heap(v);
            }
            LimbVec::Heap(v) => v.push(limb),
        }
    }

    /// Grows with zero limbs or truncates to exactly `n` limbs.
    pub(crate) fn resize(&mut self, n: usize) {
        match self {
            LimbVec::Inline { len, buf } if n <= INLINE => {
                // Limbs past `len` may hold stale values from an earlier truncate.
                let old = usize::from(*len);
                if n > old {
                    buf[old..n].fill(0);
                }
                *len = inline_len(n);
            }
            LimbVec::Inline { len, buf } => {
                let mut v = Vec::with_capacity(n);
                v.extend_from_slice(&buf[..usize::from(*len)]);
                v.resize(n, 0);
                *self = LimbVec::Heap(v);
            }
            LimbVec::Heap(v) => v.resize(n, 0),
        }
    }

    pub(crate) fn truncate(&mut self, n: usize) {
        if n < self.len() {
            self.resize(n);
        }
    }

    /// Moves the limbs inline if they fit, otherwise trims spare capacity.
    pub(crate) fn shrink_to_fit(&mut self) {
        if let LimbVec::Heap(v) = self {
            if v.len() <= INLINE {
                *self = Self::from_slice(v);
            } else {
                v.shrink_to_fit();
            }
        }
    }
}

impl Clone for LimbVec {
    /// Clones into the smallest representation (inline when the value fits).
    fn clone(&self) -> Self {
        Self::from_slice(self.as_slice())
    }
}

impl Default for LimbVec {
    fn default() -> Self {
        Self::new()
    }
}

fn inline_len(n: usize) -> u8 {
    debug_assert!(n <= INLINE);
    u8::try_from(n).expect("inline length fits u8")
}
