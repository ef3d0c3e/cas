# cas-num fuzz targets

libFuzzer targets (via [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz)) that cross-check
cas-num against GMP (through `rug`). This crate is outside the main workspace because cargo-fuzz
needs a nightly toolchain and its own build profile; the library itself stays on stable.

| Target | Checks |
|---|---|
| `natural` | `+`, `checked_sub`, `*`, `square`, `div_rem`, the identity `(a / b) b + a % b = a`, `gcd`, `lcm`, shifts, bitwise ops, comparisons and bit counts, `div_exact` of products |
| `integer` | `+ -`, `*`, `checked_div_rem` with every `Rounding`, `gcd_ext` (GMP's cofactors), `mod_inverse`, `pow`, `pow_mod` (negative exponents too), `root_rem`, perfect-power predicates |
| `rational` | `+ - * /`, comparison, `floor`/`ceil`/`trunc`/`round`, `pow` with negative exponents |
| `radix` | printing in every radix against GMP and round trips (integers and rationals); parsing of arbitrary strings against the grammar the primitive integers use, and against `i128`/`u128` where those parse; long digit strings against GMP |
| `number_theory` | `sqrt_rem`, `pow_mod` through every reduction, `ModRing` and `ModWord` operations and inverses, `crt`, `jacobi`, `primality`, `factor` (product and labels), `perfect_power` |

Inputs are decoded from the fuzzer's bytes (`src/lib.rs`): limbs mix boundary values (`0`,
`u64::MAX`, `2^63`, low masks) with raw bytes, and lengths are weighted towards the algorithm
thresholds (Karatsuba, Toom-3, divide-and-conquer division, REDC).

## Running

```sh
cargo install cargo-fuzz          # once
cd fuzz
cargo fuzz run -O natural          # until stopped or a finding
cargo fuzz run -O natural -- -max_total_time=600 -rss_limit_mb=2048
cargo fuzz build -O                # all targets
```

A finding is saved under `artifacts/<target>/`; replay it with
`cargo fuzz run -O <target> artifacts/<target>/crash-…` and minimise with `cargo fuzz tmin`.
Turn a fixed finding into a regression test in cas-num's own test suite. The `corpus/`,
`artifacts/` and `target/` directories are not committed.
