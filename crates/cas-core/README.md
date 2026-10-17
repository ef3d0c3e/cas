# cas-core

Hash-consed symbolic expressions for the CAS engine. The design is in Lific page `CAS-DOC-3`.

## Layout

| Module | Contents |
|---|---|
| `expr` | `Expr`, the 16-byte handle: a `Symbol` or an `i64` inline, or a shared interned node. Equality is O(1) (value or address) because every value has one representation. `Expr::apply`, `Expr::string`, `From<Symbol>`, `From<i64>`; `view()` and accessors; cached metadata (`depth`, `size`, `may_contain` Bloom test, `is_numeric`); a structural hash. |
| `symbol` | `Symbol`: a permanent `u32` id into a global, append-only name table. |
| `node` | Interned nodes (strings and applications) with their hash and metadata. Dropping releases children iteratively through a per-thread queue, so deep expressions use constant stack. |
| `intern` | The global interner: 64 mutex-protected shards of weak references, looked up with a borrowed key (a hit allocates nothing). A node removes its own entry when dropped, matched by address; shards shrink when sparse. `interner_stats()` reports occupancy. |

## Testing

Unit tests are in `src/tests.rs`. They cover structural equality against a plain tree model (proptest), sharing, entry removal, constant-stack drops 10⁶ deep, and concurrent construction, churn and symbol creation. The rare interleavings (full hash collisions, a dead and a live entry for the same contents) are set up deterministically on a private shard. `tests/interner.rs` runs alone in its own process and checks that the whole interner empties and shrinks.
