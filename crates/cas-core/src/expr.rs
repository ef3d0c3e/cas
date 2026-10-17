//! The expression handle.

use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::intern;
use crate::node::{Flags, Key, Kind, Meta, Node, hasher, tag};
use crate::symbol::Symbol;

/// An immutable, hash-consed expression.
///
/// Symbols and integers that fit an `i64` are stored inline; everything else
/// is a node in the global interner, shared by every expression that contains
/// it. Each value has exactly one representation, so equality is O(1):
/// inline payloads compare by value and nodes by address. Cloning is a copy
/// or a reference-count increment.
///
/// [`Expr::apply`] is purely structural: it interns exactly the head and
/// arguments it is given. The builder ([`Expr::call`], the operators in
/// [`build`](crate::build)) also flattens and sorts builtin `Flat` and
/// `Orderless` heads; arithmetic is evaluation's job.
///
/// ```
/// use cas_core::{Expr, Symbol, View};
///
/// let f = Expr::from(Symbol::new("f"));
/// let x = Expr::from(Symbol::new("x"));
/// let a = Expr::apply(&f, &[x.clone(), Expr::from(2)]);
/// let b = Expr::apply(&f, &[x.clone(), Expr::from(2)]);
/// assert_eq!(a, b); // one shared node
/// assert_ne!(a, Expr::apply(&f, &[Expr::from(2), x.clone()]));
///
/// let View::Apply { head, args } = a.view() else { unreachable!() };
/// assert_eq!(head, &f);
/// assert_eq!(args[1].as_i64(), Some(2));
/// assert_eq!(format!("{a:?}"), "f[x, 2]");
/// ```
#[derive(Clone)]
pub struct Expr(pub(crate) Repr);

#[derive(Clone)]
pub(crate) enum Repr {
    Symbol(Symbol),
    Int(i64),
    Node(Arc<Node>),
}

/// What an [`Expr`] is, borrowed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum View<'a> {
    Symbol(Symbol),
    Int(i64),
    String(&'a str),
    /// `head[args…]`; the head is any expression.
    Apply {
        head: &'a Expr,
        args: &'a [Expr],
    },
}

impl Expr {
    /// The string atom `s`.
    #[must_use]
    pub fn string(s: &str) -> Expr {
        Expr(Repr::Node(intern::intern(Key::String(s))))
    }

    /// `head[args…]`.
    #[must_use]
    pub fn apply(head: &Expr, args: &[Expr]) -> Expr {
        Expr(Repr::Node(intern::intern(Key::Apply { head, args })))
    }

    #[must_use]
    pub fn view(&self) -> View<'_> {
        match &self.0 {
            Repr::Symbol(s) => View::Symbol(*s),
            Repr::Int(x) => View::Int(*x),
            Repr::Node(n) => match &n.kind {
                Kind::String(s) => View::String(s),
                Kind::Apply { head, args } => View::Apply { head, args },
            },
        }
    }

    /// The head of an application.
    #[must_use]
    pub fn head(&self) -> Option<&Expr> {
        match self.view() {
            View::Apply { head, .. } => Some(head),
            _ => None,
        }
    }

    /// The arguments of an application; empty for atoms.
    #[must_use]
    pub fn args(&self) -> &[Expr] {
        match self.view() {
            View::Apply { args, .. } => args,
            _ => &[],
        }
    }

    #[must_use]
    pub fn is_atom(&self) -> bool {
        !matches!(self.view(), View::Apply { .. })
    }

    #[must_use]
    pub fn as_symbol(&self) -> Option<Symbol> {
        match self.0 {
            Repr::Symbol(s) => Some(s),
            _ => None,
        }
    }

    #[must_use]
    pub fn as_i64(&self) -> Option<i64> {
        match self.0 {
            Repr::Int(x) => Some(x),
            _ => None,
        }
    }

    #[must_use]
    pub fn as_str(&self) -> Option<&str> {
        match self.view() {
            View::String(s) => Some(s),
            _ => None,
        }
    }

    /// A structural hash: equal expressions have equal hashes. No per-process
    /// seed and no addresses go into it, but symbols hash by id, so it is
    /// reproducible only when symbols are created in the same order. Nothing
    /// observable (such as the canonical order) may depend on it.
    #[must_use]
    pub fn structural_hash(&self) -> u64 {
        let mut h = hasher();
        match self.0 {
            Repr::Symbol(s) => {
                h.write_u8(tag::SYMBOL);
                h.write_u32(s.id());
            }
            Repr::Int(x) => {
                h.write_u8(tag::INT);
                h.write_i64(x);
            }
            Repr::Node(ref n) => return n.hash,
        }
        h.finish()
    }

    /// Nesting depth: 1 for atoms, saturating at `u16::MAX`.
    #[must_use]
    pub fn depth(&self) -> u16 {
        self.meta().depth
    }

    /// The number of atoms and applications in the tree, counting a shared
    /// subexpression once per occurrence; saturating at `u32::MAX`.
    #[must_use]
    pub fn size(&self) -> u32 {
        self.meta().size
    }

    /// `false` if `s` certainly does not occur in the expression (as a head,
    /// an argument or bound); `true` if it may.
    #[must_use]
    pub fn may_contain(&self, s: Symbol) -> bool {
        let bit = symbol_bit(s);
        self.meta().symbols & bit == bit
    }

    /// Whether the expression is known to be a number, without evaluating:
    /// a number, a builtin constant such as `Pi`, or a builtin numeric
    /// function (`Plus`, `Sin`, …) of such. `false` means unknown: a user
    /// function can also be numeric.
    #[must_use]
    pub fn is_numeric(&self) -> bool {
        self.meta().flags.contains(Flags::NUMERIC)
    }

    /// Whether the expression contains a pattern construct (`Blank[]`,
    /// `Pattern[x, …]`, `Alternatives[…]`, …). If not, it matches only
    /// itself.
    #[must_use]
    pub fn has_pattern(&self) -> bool {
        self.meta().flags.contains(Flags::HAS_PATTERN)
    }

    pub(crate) fn meta(&self) -> Meta {
        match &self.0 {
            Repr::Symbol(s) => Meta::atom(Flags::of_symbol(*s), symbol_bit(*s)),
            Repr::Int(_) => Meta::atom(Flags::NUMERIC, 0),
            Repr::Node(n) => n.meta,
        }
    }
}

/// The symbol's bit in the Bloom filter of [`Meta::symbols`].
fn symbol_bit(s: Symbol) -> u64 {
    // Fibonacci hashing: the top 6 bits of the product are well mixed.
    1 << (s.id().wrapping_mul(0x9e37_79b9) >> 26)
}

impl From<Symbol> for Expr {
    fn from(s: Symbol) -> Expr {
        Expr(Repr::Symbol(s))
    }
}

impl From<i64> for Expr {
    fn from(x: i64) -> Expr {
        Expr(Repr::Int(x))
    }
}

impl PartialEq for Expr {
    /// Identity, which is structural equality because nodes are interned.
    fn eq(&self, other: &Expr) -> bool {
        match (&self.0, &other.0) {
            (Repr::Symbol(a), Repr::Symbol(b)) => a == b,
            (Repr::Int(a), Repr::Int(b)) => a == b,
            (Repr::Node(a), Repr::Node(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
}

impl Eq for Expr {}

impl Hash for Expr {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.structural_hash());
    }
}

/// Full form: `f[x, 1, "s"]`. Recursive; the printing ticket (CAS-23) owns the
/// real output.
impl fmt::Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.view() {
            View::Symbol(s) => write!(f, "{s:?}"),
            View::Int(x) => write!(f, "{x}"),
            View::String(s) => write!(f, "{s:?}"),
            View::Apply { head, args } => {
                write!(f, "{head:?}[")?;
                for (i, a) in args.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    write!(f, "{a:?}")?;
                }
                f.write_str("]")
            }
        }
    }
}
