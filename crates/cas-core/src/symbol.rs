//! Symbols: permanent ids into a global, append-only name table.

use std::collections::HashMap;
use std::fmt;
use std::sync::{LazyLock, PoisonError, RwLock};

use crate::attributes::Attributes;
use crate::builtin;

/// A symbol, interned by name.
///
/// Symbols are never freed, so a `Symbol` is a plain `u32` and is stored
/// inline in an [`Expr`](crate::Expr). Two symbols are equal iff their names
/// are.
///
/// The builtins (`Symbol::PLUS`, `Symbol::LIST`, …) have fixed ids and fixed
/// attributes. Other ids depend on the order in which names are first
/// interned, so nothing that must be reproducible (such as the canonical
/// order) may depend on them.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Symbol(u32);

struct Table {
    ids: HashMap<&'static str, Symbol>,
    names: Vec<&'static str>,
}

/// Starts with the builtins, so that their ids are their indices.
static TABLE: LazyLock<RwLock<Table>> = LazyLock::new(|| {
    let names: Vec<&'static str> = builtin::BUILTINS.iter().map(|b| b.name).collect();
    let ids = (0..).zip(&names).map(|(i, &n)| (n, Symbol(i))).collect();
    RwLock::new(Table { ids, names })
});

impl Symbol {
    /// The symbol named `name`, created on first use.
    ///
    /// Symbols are never freed. Code that makes up names, such as renaming
    /// pattern variables or `Module` locals, must not create a fresh symbol
    /// per use, or the table grows without bound; it should build an
    /// expression (for example `Local[x, n]`) instead, which is freed with its
    /// uses.
    ///
    /// # Panics
    ///
    /// If more than `u32::MAX` distinct symbols are created.
    #[must_use]
    pub fn new(name: &str) -> Symbol {
        if let Some(&s) = TABLE
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .ids
            .get(name)
        {
            return s;
        }
        let mut table = TABLE.write().unwrap_or_else(PoisonError::into_inner);
        // Another thread may have added it between the two locks.
        if let Some(&s) = table.ids.get(name) {
            return s;
        }
        let id = u32::try_from(table.names.len()).expect("too many symbols");
        let name: &'static str = Box::leak(name.into());
        table.names.push(name);
        table.ids.insert(name, Symbol(id));
        Symbol(id)
    }

    /// The symbol's name.
    #[must_use]
    pub fn name(self) -> &'static str {
        if let Some(b) = builtin::get(self) {
            return b.name;
        }
        TABLE.read().unwrap_or_else(PoisonError::into_inner).names[self.0 as usize]
    }

    /// Whether this is one of the builtin symbols.
    #[must_use]
    pub fn is_builtin(self) -> bool {
        builtin::get(self).is_some()
    }

    /// The fixed attributes of a builtin; `None` for other symbols, whose
    /// attributes belong to a `Context`.
    ///
    /// ```
    /// use cas_core::{Attributes, Symbol};
    ///
    /// let plus = Symbol::new("Plus");
    /// assert_eq!(plus, Symbol::PLUS);
    /// assert!(plus.builtin_attributes().unwrap().contains(Attributes::FLAT | Attributes::ORDERLESS));
    /// assert_eq!(Symbol::new("f").builtin_attributes(), None);
    /// ```
    #[must_use]
    pub fn builtin_attributes(self) -> Option<Attributes> {
        builtin::get(self).map(|b| b.attributes)
    }

    /// The builtin symbols, in id order.
    pub fn builtins() -> impl ExactSizeIterator<Item = Symbol> {
        #[expect(clippy::cast_possible_truncation, reason = "a few dozen builtins")]
        (0..builtin::COUNT as u32).map(Symbol)
    }

    pub(crate) const fn from_id(id: u32) -> Symbol {
        Symbol(id)
    }

    /// The table index; unique per name within one process.
    pub(crate) const fn id(self) -> u32 {
        self.0
    }
}

impl fmt::Debug for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}
