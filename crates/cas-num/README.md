# cas-num

Arbitrary-precision integers, rationals, modular and ball arithmetic for the CAS engine.
Correctness is the primary goal. The design is in Lific page `CAS-DOC-2`.

## Differential testing against GMP

Every operation is checked against GMP on generated inputs. The harness lives in
`src/testkit/` and is compiled only under `cfg(test)`. `rug` and `proptest` are
dev-dependencies, so they never end up in a runtime build.

### Pieces

| Item | Purpose |
|---|---|
| `testkit::Mag`, `testkit::Signed` | Normalised little-endian `u64` limb vectors (sign + magnitude). This is the neutral exchange format between our types and GMP. `Debug` prints hex with `_` at limb boundaries. |
| `ToRug` | Maps a value of ours to its GMP equivalent. It is implemented for primitives, `bool`, `Ordering`, `Option<T>`, pairs, `Mag` and `Signed`. Every cas-num type implements it as it lands. |
| `check_same(ours, gmp)` | Fails the proptest case unless the two values agree, and prints both. |
| `differential!` | Declares one `#[test]` per entry, asserting that `ours` equals `gmp`. |
| `testkit::strategy` | Adversarial generators (see below). |

Conversions go through limb slices (`Integer::from_digits` / `to_digits`), never through
strings, so the conversion is not one of the things under test. The harness tests itself
(`testkit/selftest.rs`). It checks the conversions against an independent shift-and-add
construction, and it checks that the generators actually produce each adversarial shape.

### Writing a test

```rust
use crate::testkit::{differential, strategy, ToRug};

differential! {
    /// Natural addition agrees with GMP.
    fn add(a in strategy::natural(), b in strategy::natural()) {
        ours: Natural::from_limbs(a.limbs()) + Natural::from_limbs(b.limbs()),
        gmp: a.to_rug() + b.to_rug(),
    }
}
```

An optional leading `#![proptest_config(ProptestConfig::with_cases(N))]` is forwarded to
proptest. The case count can also be changed at run time with `PROPTEST_CASES=N cargo test`.

### Strategies

Uniformly random limbs almost never hit the inputs that break arithmetic, so these
generators weight such inputs heavily:

- `limb()`: carry and sign boundary values (`0`, `1`, `u64::MAX`, `2^63 ± 1`, …),
  single bits, low masks and uniform values.
- `natural()` / `natural_up_to(max_limbs)`: zero, one, single limbs, `2^k`, `2^k ± 1` with
  `k` within one bit of a limb boundary, all-ones values, and dense, sparse or all-ones
  multi-limb values. Lengths are weighted towards 1–4 limbs.
- `natural_with_len(len)`: exactly `len` limbs.
- `len_around(&[thresholds])`: lengths `t - 1`, `t` and `t + 1`, so both sides of every
  algorithm crossover get tested.
- `unbalanced_pair(max_limbs)`: a long operand and a 1–3 limb operand, in either order.
- `signed()` / `signed_up_to(max_limbs)`: a sign plus any of the above (zero is never
  negative).
