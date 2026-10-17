//! Symbol attributes.

use std::fmt;
use std::ops::{BitAnd, BitOr, BitOrAssign};

/// A set of attributes: how evaluation and matching treat a head.
///
/// Builtin symbols have fixed attributes ([`Symbol::builtin_attributes`]);
/// the attributes of other symbols belong to a `Context`.
///
/// [`Symbol::builtin_attributes`]: crate::Symbol::builtin_attributes
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Attributes(u16);

impl Attributes {
    pub const NONE: Attributes = Attributes(0);
    /// Associative: nested applications flatten, `f[a, f[b, c]]` → `f[a, b, c]`.
    pub const FLAT: Attributes = Attributes(1 << 0);
    /// Commutative: the arguments are sorted into canonical order.
    pub const ORDERLESS: Attributes = Attributes(1 << 1);
    /// In pattern matching, `f[a]` is equivalent to `a`.
    pub const ONE_IDENTITY: Attributes = Attributes(1 << 2);
    /// Threads over lists: `f[{a, b}]` → `{f[a], f[b]}`.
    pub const LISTABLE: Attributes = Attributes(1 << 3);
    /// The first argument is not evaluated.
    pub const HOLD_FIRST: Attributes = Attributes(1 << 4);
    /// The arguments after the first are not evaluated.
    pub const HOLD_REST: Attributes = Attributes(1 << 5);
    /// No argument is evaluated: [`HOLD_FIRST`](Self::HOLD_FIRST) and
    /// [`HOLD_REST`](Self::HOLD_REST).
    pub const HOLD_ALL: Attributes = Attributes(Self::HOLD_FIRST.0 | Self::HOLD_REST.0);
    /// Applied to numbers, gives a number: `f[1, 2]` is numeric.
    pub const NUMERIC_FUNCTION: Attributes = Attributes(1 << 6);
    /// A constant: its derivative with respect to any variable is zero. A
    /// builtin constant is a number.
    pub const CONSTANT: Attributes = Attributes(1 << 7);
    /// The symbol's attributes and definitions cannot be changed. Every
    /// builtin has it.
    pub const PROTECTED: Attributes = Attributes(1 << 8);

    /// Names for `Debug`, combined attributes before their parts.
    const NAMES: [(Attributes, &'static str); 10] = [
        (Self::FLAT, "Flat"),
        (Self::ORDERLESS, "Orderless"),
        (Self::ONE_IDENTITY, "OneIdentity"),
        (Self::LISTABLE, "Listable"),
        (Self::HOLD_ALL, "HoldAll"),
        (Self::HOLD_FIRST, "HoldFirst"),
        (Self::HOLD_REST, "HoldRest"),
        (Self::NUMERIC_FUNCTION, "NumericFunction"),
        (Self::CONSTANT, "Constant"),
        (Self::PROTECTED, "Protected"),
    ];

    #[must_use]
    pub const fn union(self, other: Attributes) -> Attributes {
        Attributes(self.0 | other.0)
    }

    #[must_use]
    pub const fn intersection(self, other: Attributes) -> Attributes {
        Attributes(self.0 & other.0)
    }

    #[must_use]
    pub const fn difference(self, other: Attributes) -> Attributes {
        Attributes(self.0 & !other.0)
    }

    /// Whether every attribute of `other` is in `self`.
    #[must_use]
    pub const fn contains(self, other: Attributes) -> bool {
        self.0 & other.0 == other.0
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }
}

impl BitOr for Attributes {
    type Output = Attributes;

    fn bitor(self, other: Attributes) -> Attributes {
        self.union(other)
    }
}

impl BitOrAssign for Attributes {
    fn bitor_assign(&mut self, other: Attributes) {
        *self = self.union(other);
    }
}

impl BitAnd for Attributes {
    type Output = Attributes;

    fn bitand(self, other: Attributes) -> Attributes {
        self.intersection(other)
    }
}

/// `Flat | Orderless`, `HoldAll`, or `None`.
impl fmt::Debug for Attributes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return f.write_str("None");
        }
        let mut rest = *self;
        let mut first = true;
        for (a, name) in Self::NAMES {
            if rest.contains(a) {
                if !first {
                    f.write_str(" | ")?;
                }
                f.write_str(name)?;
                first = false;
                rest = rest.difference(a);
            }
        }
        Ok(())
    }
}
