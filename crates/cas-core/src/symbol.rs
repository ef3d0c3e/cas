//! Symbols: permanent ids into a global, append-only name table.

use std::collections::HashMap;
use std::fmt;
use std::sync::{LazyLock, PoisonError, RwLock};

/// A symbol, interned by name.
///
/// Symbols are never freed, so a `Symbol` is a plain `u32` and is stored
/// inline in an [`Expr`](crate::Expr). Two symbols are equal iff their names
/// are. Ids depend on the order in which names are first interned, so nothing
/// that must be reproducible (such as the canonical order) may depend on them.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Symbol(u32);

struct Table {
    ids: HashMap<&'static str, Symbol>,
    names: Vec<&'static str>,
}

static TABLE: LazyLock<RwLock<Table>> = LazyLock::new(|| {
    RwLock::new(Table {
        ids: HashMap::new(),
        names: Vec::new(),
    })
});

impl Symbol {
    /// The symbol named `name`, created on first use.
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
        TABLE.read().unwrap_or_else(PoisonError::into_inner).names[self.0 as usize]
    }

    /// The table index; unique per name within one process.
    pub(crate) fn id(self) -> u32 {
        self.0
    }
}

impl fmt::Debug for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}
