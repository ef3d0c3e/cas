//! The expression handle.

use std::hash::{Hash, Hasher};
use std::sync::Arc;

use cas_num::{Complex, Integer, Rational};

use crate::intern;
use crate::node::{Flags, Key, Kind, Meta, Node, hasher, tag};
use crate::symbol::Symbol;

/// An immutable, hash-consed expression.
///
/// Symbols and integers that fit an `i64` are stored inline; everything else
/// (larger integers, rationals, complex numbers, strings, applications) is a
/// node in the global interner, shared by every expression that contains it.
/// Each value has exactly one representation (an integer is never a
/// rational, a real number never a complex one, and an `i64` never a node),
/// so equality is O(1):
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
    /// An integer that fits an `i64`.
    Int(i64),
    /// An integer outside the `i64` range.
    Integer(&'a Integer),
    /// A rational number that is not an integer.
    Rational(&'a Rational),
    /// A complex number that is not real (CAS-51).
    Complex(&'a Complex),
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
                Kind::Integer(x) => View::Integer(x),
                Kind::Rational(x) => View::Rational(x),
                Kind::Complex(z) => View::Complex(z),
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

    /// The value of an integer atom, inline or not.
    #[must_use]
    pub fn to_integer(&self) -> Option<Integer> {
        match self.view() {
            View::Int(x) => Some(Integer::from(x)),
            View::Integer(x) => Some(x.clone()),
            _ => None,
        }
    }

    /// The value of a number atom: an integer or a rational.
    ///
    /// ```
    /// use cas_core::{Expr, Integer, Rational};
    ///
    /// let half = Expr::from(Rational::new(1.into(), 2.into()));
    /// assert_eq!(half.to_string(), "1/2");
    /// assert_eq!(half.to_rational(), Some(Rational::new(2.into(), 4.into())));
    /// // An integer value is an integer atom, inline when it fits an i64.
    /// assert_eq!(Expr::from(Rational::new(4.into(), 2.into())), Expr::from(2));
    /// assert_eq!(Expr::from(2).to_rational(), Some(Rational::from(2)));
    /// ```
    #[must_use]
    pub fn to_rational(&self) -> Option<Rational> {
        match self.view() {
            View::Int(x) => Some(Rational::from(x)),
            View::Integer(x) => Some(Rational::from(x)),
            View::Rational(x) => Some(x.clone()),
            _ => None,
        }
    }

    /// The value of a number atom as a complex number: any integer,
    /// rational or complex atom.
    ///
    /// ```
    /// use cas_core::{Complex, Expr, Rational};
    ///
    /// let z = Complex::new(Rational::from(2), Rational::from(3));
    /// assert_eq!(Expr::from(z.clone()).to_complex(), Some(z));
    /// // A real value is a real atom, so each number has one form.
    /// let real = Complex::new(Rational::from(5), Rational::from(0));
    /// assert_eq!(Expr::from(real), Expr::from(5));
    /// assert_eq!(Expr::from(5).to_complex(), Some(Complex::from(5)));
    /// ```
    #[must_use]
    pub fn to_complex(&self) -> Option<Complex> {
        match self.view() {
            View::Complex(z) => Some(z.clone()),
            _ => self.to_rational().map(Complex::from),
        }
    }

    /// Whether this is a number atom: an integer, a rational or a complex
    /// number. Unlike [`Expr::is_numeric`], `Pi`, `I` and `Sin[1]` are not.
    #[must_use]
    pub fn is_number(&self) -> bool {
        matches!(
            self.view(),
            View::Int(_) | View::Integer(_) | View::Rational(_) | View::Complex(_)
        )
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
pub(crate) fn symbol_bit(s: Symbol) -> u64 {
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

impl From<Integer> for Expr {
    /// Inline if it fits an `i64`, an interned atom otherwise.
    fn from(x: Integer) -> Expr {
        Expr::from(&x)
    }
}

impl From<&Integer> for Expr {
    fn from(x: &Integer) -> Expr {
        match x.to_i64() {
            Some(x) => Expr(Repr::Int(x)),
            None => Expr(Repr::Node(intern::intern(Key::Integer(x)))),
        }
    }
}

impl From<Rational> for Expr {
    /// An integer atom if the value is an integer.
    fn from(x: Rational) -> Expr {
        Expr::from(&x)
    }
}

impl From<&Rational> for Expr {
    fn from(x: &Rational) -> Expr {
        if x.is_integer() {
            Expr::from(x.numer())
        } else {
            Expr(Repr::Node(intern::intern(Key::Rational(x))))
        }
    }
}

impl From<Complex> for Expr {
    /// A rational or integer atom if the value is real.
    fn from(z: Complex) -> Expr {
        Expr::from(&z)
    }
}

impl From<&Complex> for Expr {
    fn from(z: &Complex) -> Expr {
        if z.is_real() {
            Expr::from(z.re())
        } else {
            Expr(Repr::Node(intern::intern(Key::Complex(z))))
        }
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
