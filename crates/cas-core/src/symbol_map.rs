//! A persistent map keyed by symbol (CAS-DOC-4 §5).
//!
//! A 32-way trie on the symbol id, which is dense from 0: a lookup follows
//! one path, an update copies the nodes on its path only when they are
//! shared (`Arc::make_mut`), and a clone shares everything. So a context
//! that is not shared updates in place, and one that is pays O(log₃₂ n) per
//! change instead of a copy of the whole map.

use std::fmt;
use std::sync::Arc;

use crate::symbol::Symbol;

const BITS: u32 = 5;
const WIDTH: usize = 1 << BITS;
const MASK: u32 = (1 << BITS) - 1;

/// A persistent map from symbols to `V`.
pub(crate) struct SymbolMap<V> {
    root: Option<Arc<Node<V>>>,
    /// The bits of the id the root's children are indexed by: the trie holds
    /// ids below `WIDTH << shift`.
    shift: u32,
    len: usize,
}

enum Node<V> {
    Branch(Box<[Option<Arc<Node<V>>>; WIDTH]>),
    Leaf(V),
}

impl<V: Clone> Clone for Node<V> {
    fn clone(&self) -> Self {
        match self {
            Node::Branch(children) => Node::Branch(children.clone()),
            Node::Leaf(v) => Node::Leaf(v.clone()),
        }
    }
}

impl<V> Clone for SymbolMap<V> {
    fn clone(&self) -> Self {
        SymbolMap {
            root: self.root.clone(),
            shift: self.shift,
            len: self.len,
        }
    }
}

impl<V> Default for SymbolMap<V> {
    fn default() -> Self {
        SymbolMap {
            root: None,
            shift: 0,
            len: 0,
        }
    }
}

fn empty_branch<V>() -> Box<[Option<Arc<Node<V>>>; WIDTH]> {
    Box::new(std::array::from_fn(|_| None))
}

impl<V: Clone> SymbolMap<V> {
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// The largest id the trie holds now, plus one.
    fn capacity(&self) -> u64 {
        (WIDTH as u64) << self.shift
    }

    pub(crate) fn get(&self, s: Symbol) -> Option<&V> {
        let id = s.id();
        if u64::from(id) >= self.capacity() {
            return None;
        }
        let mut node = self.root.as_deref()?;
        let mut shift = self.shift;
        loop {
            match node {
                Node::Branch(children) => {
                    node = children[((id >> shift) & MASK) as usize].as_deref()?;
                    if shift == 0 {
                        // The level below the last branch holds the leaves.
                        return match node {
                            Node::Leaf(v) => Some(v),
                            Node::Branch(_) => None,
                        };
                    }
                    shift -= BITS;
                }
                Node::Leaf(v) => return Some(v),
            }
        }
    }

    /// The value of `s`, for a change, inserting `make()` when it has none.
    pub(crate) fn entry(&mut self, s: Symbol, make: impl FnOnce() -> V) -> &mut V {
        let id = s.id();
        while u64::from(id) >= self.capacity() {
            // Grow a level: the old root becomes the first child.
            let mut children = empty_branch();
            children[0] = self.root.take();
            self.root = Some(Arc::new(Node::Branch(children)));
            self.shift += BITS;
        }
        let mut node = self
            .root
            .get_or_insert_with(|| Arc::new(Node::Branch(empty_branch())));
        let mut shift = self.shift;
        let mut inserted = false;
        loop {
            let Node::Branch(children) = Arc::make_mut(node) else {
                unreachable!("branches down to the leaves");
            };
            let child = &mut children[((id >> shift) & MASK) as usize];
            if shift == 0 {
                let leaf = child.get_or_insert_with(|| {
                    inserted = true;
                    Arc::new(Node::Leaf(make()))
                });
                if inserted {
                    self.len += 1;
                }
                let Node::Leaf(v) = Arc::make_mut(leaf) else {
                    unreachable!("a leaf");
                };
                return v;
            }
            node = child.get_or_insert_with(|| Arc::new(Node::Branch(empty_branch())));
            shift -= BITS;
        }
    }

    pub(crate) fn insert(&mut self, s: Symbol, value: V) {
        let mut value = Some(value);
        let slot = self.entry(s, || value.take().expect("inserted once"));
        if let Some(value) = value {
            *slot = value;
        }
    }

    /// Removes `s`'s value, if any. Empty branches are kept.
    pub(crate) fn remove(&mut self, s: Symbol) -> Option<V> {
        let id = s.id();
        self.get(s)?;
        let mut node = self.root.as_mut().expect("a value");
        let mut shift = self.shift;
        loop {
            let Node::Branch(children) = Arc::make_mut(node) else {
                unreachable!("branches down to the leaves");
            };
            let child = &mut children[((id >> shift) & MASK) as usize];
            if shift == 0 {
                self.len -= 1;
                let leaf = child.take().expect("a value");
                return Some(match Arc::try_unwrap(leaf) {
                    Ok(Node::Leaf(v)) => v,
                    Ok(Node::Branch(_)) => unreachable!("a leaf"),
                    Err(shared) => match &*shared {
                        Node::Leaf(v) => v.clone(),
                        Node::Branch(_) => unreachable!("a leaf"),
                    },
                });
            }
            node = child.as_mut().expect("a value");
            shift -= BITS;
        }
    }

    /// Every entry, by increasing id.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (Symbol, &V)> {
        let mut out = Vec::with_capacity(self.len);
        if let Some(root) = &self.root {
            collect(root, 0, self.shift, &mut out);
        }
        out.into_iter()
    }
}

fn collect<'a, V>(node: &'a Node<V>, prefix: u32, shift: u32, out: &mut Vec<(Symbol, &'a V)>) {
    match node {
        Node::Leaf(v) => out.push((Symbol::from_id(prefix), v)),
        Node::Branch(children) => {
            for (i, child) in children.iter().enumerate() {
                if let Some(child) = child {
                    let i = u32::try_from(i).expect("a small index");
                    let id = prefix | (i << shift);
                    match &**child {
                        Node::Leaf(v) => out.push((Symbol::from_id(id), v)),
                        Node::Branch(_) => collect(child, id, shift - BITS, out),
                    }
                }
            }
        }
    }
}

impl<V: Clone + fmt::Debug> fmt::Debug for SymbolMap<V> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_map().entries(self.iter()).finish()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use proptest::prelude::*;

    use super::SymbolMap;
    use crate::symbol::Symbol;

    proptest! {
        /// The map agrees with a `BTreeMap` through inserts and removals,
        /// and a clone taken midway is not changed by later updates.
        #[test]
        fn behaves_like_a_map(ops in prop::collection::vec((0u32..5000, any::<bool>(), any::<u16>()), 0..200)) {
            let mut map = SymbolMap::default();
            let mut model = BTreeMap::new();
            let mut snapshot = None;
            for (i, (id, insert, value)) in ops.iter().enumerate() {
                let s = Symbol::from_id(*id);
                if *insert {
                    map.insert(s, *value);
                    model.insert(*id, *value);
                } else {
                    prop_assert_eq!(map.remove(s), model.remove(id));
                }
                if i == ops.len() / 2 {
                    snapshot = Some((map.clone(), model.clone()));
                }
                prop_assert_eq!(map.len(), model.len());
            }
            for (id, value) in &model {
                prop_assert_eq!(map.get(Symbol::from_id(*id)), Some(value));
            }
            let entries: Vec<(u32, u16)> = map.iter().map(|(s, v)| (s.id(), *v)).collect();
            let expected: Vec<(u32, u16)> = model.iter().map(|(k, v)| (*k, *v)).collect();
            prop_assert_eq!(entries, expected);
            if let Some((map, model)) = snapshot {
                for (id, value) in &model {
                    prop_assert_eq!(map.get(Symbol::from_id(*id)), Some(value));
                }
                prop_assert_eq!(map.len(), model.len());
            }
        }
    }
}
