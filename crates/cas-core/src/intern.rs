//! The global interner: 64 shards of weak references, one per hash range.
//!
//! A lookup holds its shard's lock and never drops a strong reference while
//! holding it: the last reference running [`Node`]'s drop would try to take the
//! same lock. A node removes its own entry when dropped, matched by address, so
//! a dead node is never read. Between a node's last reference going away and
//! its entry being removed, the entry cannot be upgraded; lookups skip it, so
//! the same contents may briefly have a dead and a live entry.

use std::ptr;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError, Weak};

use hashbrown::HashTable;

use crate::node::{Key, Node};

const SHARDS: usize = 64;

/// A shard grows as `HashTable` does and shrinks to twice its length when
/// below 1/8 of its capacity, so churn after a large build stays cheap.
const SHRINK_BELOW: usize = 8;
/// No shrinking below this capacity.
const MIN_CAPACITY: usize = 64;

pub(crate) struct Shard {
    /// `(hash, node)`: the hash is kept here since a dead node cannot be read.
    table: HashTable<(u64, Weak<Node>)>,
}

static TABLE: [Mutex<Shard>; SHARDS] = [const {
    Mutex::new(Shard {
        table: HashTable::new(),
    })
}; SHARDS];

/// The shard for `hash`. `HashTable` uses the low bits for the bucket and the
/// top 7 for its tag, so the shard comes from neither.
fn shard(hash: u64) -> MutexGuard<'static, Shard> {
    let i = (hash >> 40) as usize & (SHARDS - 1);
    // A panic while holding the lock cannot leave the table inconsistent in a
    // way that matters to us: at worst an entry is missing or dead.
    TABLE[i].lock().unwrap_or_else(PoisonError::into_inner)
}

/// The node with contents `key`, created if there is no live one.
pub(crate) fn intern(key: Key<'_>) -> Arc<Node> {
    let hash = key.hash();
    let meta = key.meta();
    // Upgraded candidates that did not match: dropped after the lock.
    let mut mismatched = Vec::new();
    let node = shard(hash).get_or_insert(hash, key, &mut mismatched, || {
        Arc::new(Node {
            hash,
            meta,
            eval_stamp: 0.into(),
            kind: key.to_kind(),
        })
    });
    drop(mismatched);
    node
}

/// Removes `node`'s entry; called from its drop.
#[inline]
pub(crate) fn remove(node: &Node) {
    let removed = shard(node.hash).remove(node);
    debug_assert!(removed, "a dropped node had no interner entry");
}

/// Occupancy of the global interner, for diagnostics.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub struct InternerStats {
    /// Live nodes, plus any whose last reference is being dropped right now.
    pub entries: usize,
    /// Entries the shards can hold without growing.
    pub capacity: usize,
}

/// Current occupancy of the global interner. Shards are read one at a time,
/// so under concurrent use this is not a single snapshot.
#[must_use]
pub fn interner_stats() -> InternerStats {
    let mut stats = InternerStats {
        entries: 0,
        capacity: 0,
    };
    for s in &TABLE {
        let s = s.lock().unwrap_or_else(PoisonError::into_inner);
        stats.entries += s.table.len();
        stats.capacity += s.table.capacity();
    }
    stats
}

impl Shard {
    #[cfg(test)]
    pub(crate) fn new_for_test() -> Shard {
        Shard {
            table: HashTable::new(),
        }
    }

    /// Adds an entry as is, e.g. a dead one.
    #[cfg(test)]
    pub(crate) fn insert_raw(&mut self, hash: u64, node: Weak<Node>) {
        self.table.insert_unique(hash, (hash, node), |e| e.0);
    }

    /// The addresses the entries with `hash` point at.
    #[cfg(test)]
    pub(crate) fn targets(&self, hash: u64) -> Vec<*const Node> {
        self.table
            .iter_hash(hash)
            .map(|(_, w)| w.as_ptr())
            .collect()
    }

    pub(crate) fn get_or_insert(
        &mut self,
        hash: u64,
        key: Key<'_>,
        mismatched: &mut Vec<Arc<Node>>,
        new: impl FnOnce() -> Arc<Node>,
    ) -> Arc<Node> {
        for (h, w) in self.table.iter_hash(hash) {
            if *h != hash {
                continue;
            }
            if let Some(n) = w.upgrade() {
                if key.matches(&n.kind) {
                    return n;
                }
                mismatched.push(n);
            }
        }
        let n = new();
        self.table
            .insert_unique(hash, (hash, Arc::downgrade(&n)), |e| e.0);
        n
    }

    /// Removes the entry pointing at `node`; whether there was one.
    pub(crate) fn remove(&mut self, node: &Node) -> bool {
        let Ok(entry) = self.table.find_entry(node.hash, |(h, w)| {
            *h == node.hash && ptr::eq(w.as_ptr(), node)
        }) else {
            return false;
        };
        // The weak reference dropped here is not the last: the node's own
        // allocation is still held while it is being dropped.
        entry.remove();
        let len = self.table.len();
        if self.table.capacity() > MIN_CAPACITY && len < self.table.capacity() / SHRINK_BELOW {
            self.table.shrink_to(len * 2, |e| e.0);
        }
        true
    }
}

/// Entries, live or dead, with structural hash `hash`.
#[cfg(test)]
pub(crate) fn entries_with_hash(hash: u64) -> usize {
    shard(hash)
        .table
        .iter_hash(hash)
        .filter(|(h, _)| *h == hash)
        .count()
}
