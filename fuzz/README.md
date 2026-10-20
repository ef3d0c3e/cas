# cas-num fuzz targets

libFuzzer targets (via [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz)) that cross-check
cas-num against GMP (through `rug`), plus targets for the cas-core reader and evaluator and the cas-lang lexer, parser, lowering and printer. This crate is outside the
main workspace because cargo-fuzz needs a nightly toolchain and its own build profile; the
libraries themselves stay on stable.

| Target | Checks |
|---|---|
| `natural` | `+`, `checked_sub`, `*`, `square`, `div_rem`, the identity `(a / b) b + a % b = a`, `gcd`, `lcm`, shifts, bitwise ops, comparisons and bit counts, `div_exact` of products |
| `integer` | `+ -`, `*`, `checked_div_rem` with every `Rounding`, `gcd_ext` (GMP's cofactors), `mod_inverse`, `pow`, `pow_mod` (negative exponents too), `root_rem`, perfect-power predicates |
| `rational` | `+ - * /`, comparison, `floor`/`ceil`/`trunc`/`round`, `pow` with negative exponents |
| `radix` | printing in every radix against GMP and round trips (integers and rationals); parsing of arbitrary strings against the grammar the primitive integers use, and against `i128`/`u128` where those parse; long digit strings against GMP |
| `ball` | ball `+ - * /` and `sqrt` (of squares, exact answer) on dyadic midpoints with or without a radius, at fuzzer-chosen precisions: the result contains the exact `Rational` value at the midpoints (checked against GMP); division by a ball containing 0 gives the indeterminate ball |
| `ball_text` | arbitrary strings into `Ball::from_decimal` (at fuzzer-chosen precisions), `Float::from_hex` and `Mag::from_hex`: no panic; an accepted decimal with a small exponent gives a ball containing its exact value, and prints in both styles; accepted hex floats print and parse back exactly |
| `ball_fns` | The real ball functions (`exp`, `expm1`, `log`, `log1p`, `sin`, `cos`, `tan`, `sin_pi`, `cos_pi`, `atan`, `asin`, `acos`, the hyperbolic functions and inverses, `sqrt`, `root`, `pow`, `pow_ball`, `atan2`) at the midpoint, both ends and an interior point of fuzzer-chosen balls (every pair of them for two arguments): the result meets MPFR's bracket at `prec + 64` bits; a point outside the domain or at a pole gives the indeterminate ball; an exact input gives a finite result with at least `prec - 3` bits |
| `complex_fns` | The complex ball functions, `inv`, `mul`, `div`, `pow`, `pow_int`, `root`, `abs` and `arg` at the centre, the corners and an interior point: both parts meet MPC's brackets at `prec + 192` bits; a pole gives a non-finite result; an exact input gives radii within `2^(6 - prec)` of the larger part |
| `machine_fns` | The machine-ball functions and arithmetic against MPFR at 192 bits: `None` is always allowed; a `Some` meets the bracket at the midpoint, both ends and an interior point, and is `None` wherever MPFR has no finite value |
| `ball_consts` | π, log 2, e, γ and G at sequences of fuzzer-chosen precisions from a cleared cache: each meets MPFR's bracket at `prec + 64` bits, knows `prec - 1` bits and has a midpoint of at most `prec` bits |
| `number_theory` | `sqrt_rem`, `pow_mod` through every reduction, `ModRing` and `ModWord` operations and inverses, `crt`, `jacobi`, `primality`, `factor` (product and labels), `perfect_power` |
| `full_form` | The cas-core full-form reader on arbitrary text: it never panics, errors are within the text on character boundaries, and whatever reads prints in full form as text that reads back as the same expression |
| `eval` | The cas-core evaluator (CAS-60) on random expressions, `D` among their builtins with its system rules (CAS-59), and indexing, slices, sets, factorials, `./`, `.^`, `Unapply` and numbers applied as functions (CAS-123), under random rules (exact and patterned, with conditions, up- and sub-rules), attributes set before or after them, defaults, the Complex domain and assumptions, with small limits; each input runs on a thread with just the stack its recursion limit allows (`STACK_PER_LEVEL` a level). No panic (with `-a`, the debug assertions check the natives' shortcuts against their full computation), no stack overflow, no time-budget overrun; a result evaluates to itself; a trace at the Verbose level changes no result, is a consistent tree, and replays to the result; verification and the warnings' messages never panic. The limits bound the time only of programs that do not grow their work at each level (the evaluator has no budget per call, CAS-DOC-4 R8), so a right-hand side calls a user head at most once, applies nothing to parts (`Map`, `ReplaceAll`, `Function`) and uses each pattern variable at most once (`k[n___] := k[n, n, n]` triples its arguments), and a rule on a `Flat` head has no tests or conditions (a test evaluates the new group a blank binds, which matches again with a budget of its own); there is no `ReplaceRepeated` (a rule that matches what it makes grows the expression at each round); the number bits are drawn from 2⁶ to 2¹⁴ |
| `lang_lex` | The cas-lang lexer on arbitrary text (bytes read as lossy UTF-8): the tokens cover the text exactly, each lexes alone as one token of the same kind, every error token has its error, errors are ordered and on character boundaries; the values of number, string and superscript tokens and the names of identifiers never panic (number exponents are bounded) |
| `lang_parse` | The cas-lang parser on arbitrary text: the tree's text is the input, every node's children cover it exactly and in order, errors are ordered, within the text and on character boundaries, and nesting beyond `MAX_DEPTH` stops cleanly instead of overflowing the stack |
| `lang_lower` | Parsing and lowering of arbitrary text (seeded with the `lang_parse` corpus): lowering never panics, every diagnostic, fix and statement span lies within the text, and every lowered statement prints in full form |
| `lang_print` | The input form on what lowering makes of arbitrary text (seeded with the `lang_parse` corpus): every expression that prints reads back without errors or notes as itself, up to the arithmetic the text writes out (`fold_arithmetic`), and printing what it reads back is stable |

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

Long runs share the machine, so run them through the memory-capped wrapper, with libFuzzer's own
limits below the cap:

```sh
CAP_MEM=3G CAP_TIME=700 ../scripts/capped.sh cargo fuzz run -O ball -- \
    -max_total_time=600 -rss_limit_mb=2048 -malloc_limit_mb=1024
```

The ball function targets (`ball_fns`, `complex_fns`, `machine_fns`, `ball_consts`, CAS-113) and
`eval` (CAS-60) run best with debug assertions and without AddressSanitizer: cas-num and cas-core
have no `unsafe` code, ASan's quarantine more than doubles the resident memory, and the debug
assertions check the `Float` invariants of every value, and the natives' shortcuts against their
full computation:

```sh
CAP_MEM=3G CAP_TIME=3700 ../scripts/capped.sh cargo fuzz run -O -a -s none ball_fns -- \
    -max_total_time=3600 -rss_limit_mb=2048 -timeout=60
```

They fill the constant cache and the function tables once (`warm_up` in `src/ball.rs`), beyond
every precision an input asks for, so that a result never depends on earlier inputs;
`ball_consts` clears the cache before each input instead. A sample point longer than
`2 prec + 256` bits (an end across a huge exponent gap) is rounded toward the midpoint, which keeps
it in the ball and MPFR's work bounded.

A stack overflow on `eval`'s evaluation thread ends the process with SIGSEGV and no artifact,
libFuzzer's handler having no stack to run on; replay the corpus in batches
(`…/release/eval -runs=1 corpus/eval/<files>`), and the last `Running:` line before the crash
names the input.

A finding is saved under `artifacts/<target>/`; replay it with
`cargo fuzz run -O <target> artifacts/<target>/crash-…` and minimise with `cargo fuzz tmin`.
Turn a fixed finding into a regression test in cas-num's own test suite. The `corpus/`,
`artifacts/` and `target/` directories are not committed.
