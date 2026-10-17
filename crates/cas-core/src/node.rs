//! Interned nodes: everything that is not stored inline in an [`Expr`].

use std::cell::{Cell, RefCell};
use std::hash::{BuildHasher, Hasher};
use std::mem;
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use foldhash::fast::{FixedState, FoldHasher};

use crate::attributes::Attributes;
use crate::builtin;
use crate::expr::{Expr, Repr};
use crate::intern;
use crate::symbol::Symbol;

/// A hash-consed node. Created only by the interner, which holds a weak
/// reference to every live node; [`Drop`] removes it.
pub(crate) struct Node {
    /// Structural hash (see [`hasher`]).
    pub(crate) hash: u64,
    pub(crate) meta: Meta,
    /// `StateId` of the context state this node was last found to be a fixpoint
    /// in; 0 for never.
    pub(crate) eval_stamp: AtomicU64,
    pub(crate) kind: Kind,
}

pub(crate) enum Kind {
    String(Box<str>),
    Apply { head: Expr, args: Box<[Expr]> },
}

/// A node's contents, borrowed: what the interner looks up, so that a hit
/// allocates nothing.
#[derive(Clone, Copy)]
pub(crate) enum Key<'a> {
    String(&'a str),
    Apply { head: &'a Expr, args: &'a [Expr] },
}

/// Hash tags of the expression kinds, so that different kinds with equal
/// payloads hash apart.
pub(crate) mod tag {
    pub(crate) const SYMBOL: u8 = 0;
    pub(crate) const INT: u8 = 1;
    pub(crate) const STRING: u8 = 2;
    pub(crate) const APPLY: u8 = 3;
}

/// The structural hasher: fixed seed, no per-process randomness.
pub(crate) fn hasher() -> FoldHasher<'static> {
    const STATE: FixedState = FixedState::with_seed(0x6361_735f_636f_7265);
    STATE.build_hasher()
}

/// Context-free metadata, computed bottom-up from the children at
/// construction.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Meta {
    pub(crate) flags: Flags,
    /// 1 for atoms; saturating.
    pub(crate) depth: u16,
    /// Atoms and nodes in the tree, counting shared subterms once per
    /// occurrence; saturating.
    pub(crate) size: u32,
    /// Bloom filter of the symbols in the tree, bound or free.
    pub(crate) symbols: u64,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Flags(u16);

impl Flags {
    pub(crate) const NONE: Flags = Flags(0);
    /// Known to be a number: a number atom, a builtin constant, or a
    /// builtin `NUMERIC_FUNCTION` applied to at least one argument, all
    /// numeric. Arity is not checked. Sound but incomplete: user numeric
    /// functions are the evaluator's business.
    pub(crate) const NUMERIC: Flags = Flags(1);
    /// Contains a pattern construct ([`builtin::is_pattern_head`]), so it may
    /// match more than itself.
    pub(crate) const HAS_PATTERN: Flags = Flags(2);

    pub(crate) fn contains(self, other: Flags) -> bool {
        self.0 & other.0 == other.0
    }

    fn set(&mut self, other: Flags, on: bool) {
        if on {
            self.0 |= other.0;
        }
    }

    pub(crate) fn of_symbol(s: Symbol) -> Flags {
        let constant = s
            .builtin_attributes()
            .is_some_and(|a| a.contains(Attributes::CONSTANT));
        if constant {
            Flags::NUMERIC
        } else {
            Flags::NONE
        }
    }
}

impl Meta {
    pub(crate) fn atom(flags: Flags, symbols: u64) -> Meta {
        Meta {
            flags,
            depth: 1,
            size: 1,
            symbols,
        }
    }
}

impl Key<'_> {
    pub(crate) fn hash(self) -> u64 {
        let mut h = hasher();
        match self {
            Key::String(s) => {
                h.write_u8(tag::STRING);
                h.write_usize(s.len());
                h.write(s.as_bytes());
            }
            Key::Apply { head, args } => {
                h.write_u8(tag::APPLY);
                h.write_u64(head.structural_hash());
                h.write_usize(args.len());
                for a in args {
                    h.write_u64(a.structural_hash());
                }
            }
        }
        h.finish()
    }

    pub(crate) fn meta(self) -> Meta {
        match self {
            Key::String(_) => Meta::atom(Flags::NONE, 0),
            Key::Apply { head, args } => {
                let mut m = head.meta();
                let head_sym = head.as_symbol();
                let mut numeric = !args.is_empty()
                    && head_sym
                        .and_then(Symbol::builtin_attributes)
                        .is_some_and(|a| a.contains(Attributes::NUMERIC_FUNCTION));
                let mut pattern = m.flags.contains(Flags::HAS_PATTERN)
                    || head_sym.is_some_and(builtin::is_pattern_head);
                for a in args {
                    let am = a.meta();
                    m.depth = m.depth.max(am.depth);
                    m.size = m.size.saturating_add(am.size);
                    m.symbols |= am.symbols;
                    numeric &= am.flags.contains(Flags::NUMERIC);
                    pattern |= am.flags.contains(Flags::HAS_PATTERN);
                }
                m.flags = Flags::NONE;
                m.flags.set(Flags::NUMERIC, numeric);
                m.flags.set(Flags::HAS_PATTERN, pattern);
                m.depth = m.depth.saturating_add(1);
                m.size = m.size.saturating_add(1);
                m
            }
        }
    }

    /// Whether `kind` has exactly these contents (children by identity).
    pub(crate) fn matches(self, kind: &Kind) -> bool {
        match (self, kind) {
            (Key::String(a), Kind::String(b)) => a == &**b,
            (Key::Apply { head, args }, Kind::Apply { head: h, args: a }) => {
                head == h && args.len() == a.len() && args.iter().zip(a.iter()).all(|(x, y)| x == y)
            }
            _ => false,
        }
    }

    pub(crate) fn to_kind(self) -> Kind {
        match self {
            Key::String(s) => Kind::String(s.into()),
            Key::Apply { head, args } => Kind::Apply {
                head: head.clone(),
                args: args.into(),
            },
        }
    }
}

thread_local! {
    /// Nodes waiting to be released by the outermost [`Node::drop`] on this
    /// thread.
    static RELEASE_QUEUE: RefCell<Vec<Arc<Node>>> = const { RefCell::new(Vec::new()) };
    static RELEASING: Cell<bool> = const { Cell::new(false) };
}

impl Drop for Node {
    /// Unregisters the node, then releases its children iteratively: nested
    /// drops only queue their children, and the outermost drop on the thread
    /// drains the queue, so dropping a deep expression uses constant stack.
    fn drop(&mut self) {
        intern::remove(self);
        let Kind::Apply { head, args } = &mut self.kind else {
            return;
        };
        let head = mem::replace(head, Expr(Repr::Int(0)));
        let args = mem::take(args);
        // Only children this node holds the last reference to would recurse;
        // the rest are released here. (Losing a race to the last reference
        // costs one level of recursion, not unbounded depth.)
        let mut children = std::iter::once(head)
            .chain(args)
            .filter_map(|e| match e.0 {
                Repr::Node(n) if Arc::strong_count(&n) == 1 => Some(n),
                _ => None,
            })
            .peekable();
        if children.peek().is_none() {
            return;
        }
        // During thread-local destruction the queue may be gone; the
        // children are then dropped recursively, as plain `Arc`s would be.
        let Ok(outermost) = RELEASING.try_with(|r| !r.replace(true)) else {
            return;
        };
        // Lazily: a guard built and dropped here would end the release.
        let _reset = outermost.then(|| ResetReleasing);
        // One push per borrow: advancing `children` can drop a child that
        // became the last reference meanwhile, and its drop pushes too.
        for n in children {
            if RELEASE_QUEUE.try_with(|q| q.borrow_mut().push(n)).is_err() {
                return;
            }
        }
        if !outermost {
            return;
        }
        // Each pop is its own borrow: dropping the node pushes onto the queue.
        while let Some(n) = RELEASE_QUEUE.with(|q| q.borrow_mut().pop()) {
            drop(n);
        }
        RELEASE_QUEUE.with(|q| q.borrow_mut().shrink_to(RETAINED_QUEUE_CAPACITY));
    }
}

/// Queue capacity kept after a large release.
const RETAINED_QUEUE_CAPACITY: usize = 1024;

/// Ends the release on scope exit, unwinding included, so a panic cannot
/// leave the thread only queueing and never draining.
struct ResetReleasing;

impl Drop for ResetReleasing {
    fn drop(&mut self) {
        RELEASING.with(|r| r.set(false));
    }
}
