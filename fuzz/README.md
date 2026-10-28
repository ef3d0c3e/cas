# cas-num fuzz targets

libFuzzer targets (via [cargo-fuzz](https://github.com/rust-fuzz/cargo-fuzz)) that cross-check
cas-num against GMP (through `rug`), plus targets for the cas-core reader and evaluator, cas-poly's products and factoring, and the cas-lang lexer, parser, lowering, printer and session. This crate is outside the
main workspace because cargo-fuzz needs a nightly toolchain and its own build profile; the
libraries themselves stay on stable.

| Target | Checks |
|---|---|
| `natural` | `+`, `checked_sub`, `*`, `square`, `div_rem`, the identity `(a / b) b + a % b = a`, `gcd`, `lcm`, shifts, bitwise ops, comparisons and bit counts, `div_exact` of products |
| `integer` | `+ -`, `*`, `checked_div_rem` with every `Rounding`, `gcd_ext` (GMP's cofactors), `mod_inverse`, `pow`, `pow_mod` (negative exponents too), `root_rem`, perfect-power predicates |
| `rational` | `+ - * /`, comparison, `floor`/`ceil`/`trunc`/`round`, `pow` with negative exponents |
| `radix` | printing in every radix against GMP and round trips (integers and rationals); parsing of arbitrary strings against the grammar the primitive integers use, and against `i128`/`u128` where those parse; long digit strings against GMP |
| `ball` | ball `+ - * /` and `sqrt` (of squares, exact answer) on dyadic midpoints with or without a radius, at fuzzer-chosen precisions: the result contains the exact `Rational` value at the midpoints (checked against GMP); division by a ball containing 0 gives the indeterminate ball |
| `ball_text` | arbitrary strings into `Ball::from_decimal` (at fuzzer-chosen precisions), `Float::from_hex` and `Mag::from_hex`: no panic; an accepted decimal with a small exponent gives a ball containing its exact value, prints in both styles, and has a literal that reads back to its midpoint (CAS-148); accepted hex floats print and parse back exactly |
| `ball_fns` | The real ball functions (`exp`, `expm1`, `log`, `log1p`, `sin`, `cos`, `tan`, `sin_pi`, `cos_pi`, `atan`, `asin`, `acos`, the hyperbolic functions and inverses, `sqrt`, `root`, `pow`, `pow_ball`, `atan2`, and the special `gamma`, `zeta`, `erf`, `erfc`, `ei` and `digamma`) at the midpoint, both ends and an interior point of fuzzer-chosen balls (every pair of them for two arguments): the result meets MPFR's bracket at `prec + 64` bits; a point outside the domain or at a pole gives the indeterminate ball; an exact input gives a finite result with at least `prec - 3` bits |
| `special_fns` | The special ball functions MPFR lacks (CAS-353): `si` and `ci` meet their series summed in MPFR far above the precision at the points of balls within 256 of 0, and agree at two precisions; `polygamma` agrees at two precisions and with its recurrence, and is indeterminate at the poles and past order 256; `lambert_w` on both branches holds `W e^W ∋ x` at every point, stays on its branch, and is indeterminate outside its domain; an exact input in the domain gives at least `prec - 3` bits |
| `complex_fns` | The complex ball functions, `inv`, `mul`, `div`, `pow`, `pow_int`, `root`, `abs` and `arg` at the centre, the corners and an interior point: both parts meet MPC's brackets at `prec + 192` bits; a pole gives a non-finite result; an exact input gives radii within `2^(6 - prec)` of the larger part |
| `machine_fns` | The machine-ball functions and arithmetic against MPFR at 192 bits: `None` is always allowed; a `Some` meets the bracket at the midpoint, both ends and an interior point, and is `None` wherever MPFR has no finite value |
| `ball_consts` | π, log 2, e, γ and G at sequences of fuzzer-chosen precisions from a cleared cache: each meets MPFR's bracket at `prec + 64` bits, knows `prec - 1` bits and has a midpoint of at most `prec` bits |
| `number_theory` | `sqrt_rem`, `pow_mod` through every reduction, `ModRing` and `ModWord` operations and inverses, `crt`, `jacobi`, `primality`, `factor` (product and labels), `perfect_power` |
| `poly_factor` | cas-poly's factoring over ℤ: univariate (CAS-182), products of powers of up to four random polynomials of degree at most 8 (total degree at most 32, so van Hoeij's recombination runs past 8 local factors, CAS-188), in two or three variables (CAS-184), products of powers of up to three random polynomials (total degree at most 12), with coefficients small or of up to three limbs, and over ℚ(i) (CAS-186), in one or two variables with small Gaussian coefficients (total degree at most 8): the factors multiply back, are sorted, distinct, primitive and positive, and each factors as itself; the factorisation is the union of the parts' factorisations, with the product of their constants |
| `poly_sparse` | cas-poly's sparse products (CAS-136, CAS-155) over ℤ, ℚ (every other term over a denominator) and word prime fields (moduli from 2 to `2^64 − 59`, residues below and above `2^63`), in one to four variables, up to 64 terms a factor: exponents in a small box, which the arrays take, now and then one up to `2^20`, which sends the product to the heap; coefficients small, at the edges of `i64`, or of up to three limbs, for the three-word, mixed and ring accumulators. The product is the naive one either way round, and divides back exactly |
| `poly_gcd` | cas-poly's gcds (CAS-137, CAS-152) of `f g` and `f h` for random sparse `f`, `g`, `h` in one to three variables, up to 8 terms of exponents at most 3, coefficients small or of two limbs, over ℤ and modulo `2^61 − 1`: every method that applies (the heuristic, Brown's, Zippel's) and the default give the gcd and cofactors of subresultants, which multiply back |
| `poly_integrate` | cas-poly's integration of rational functions over ℚ (CAS-196, CAS-197): random numerators over products of powers of up to three random factors of degree 1 to 3, coefficients small, fractions or of two limbs. At three rational points that are not poles, the derivative of the result plus the part left is the integrand, exactly, each pair of quadratic residues in its field ℚ(√d) or ℚ(√−d), where its part must come out rational; the rational part is in lowest terms over a monic denominator, the arguments and norms monic, the arguments square-free, a real pair's norm `S S̄`; the part left is proper over a square-free denominator, and there is none when every factor has degree 2 at most |
| `full_form` | The cas-core full-form reader on arbitrary text: it never panics, errors are within the text on character boundaries, and whatever reads prints in full form as text that reads back as the same expression |
| `eval` | The cas-core evaluator (CAS-60) on random expressions, `D` among their builtins with its system rules (CAS-59), and indexing, slices, sets, factorials, `./`, `.^`, `Unapply` and numbers applied as functions (CAS-123), and balls, real and complex, with `Ball` (CAS-83), and decimal literals and `exact` (CAS-148), with typed blanks on the heads of atoms (CAS-150), and `N` with the elementary functions of balls (CAS-84), and the polynomial builtins with factoring (CAS-140, CAS-185), and `Integrate` (CAS-194), and `Series` and `Limit` (CAS-204, CAS-205), under random rules (exact and patterned, with conditions, up- and sub-rules), attributes set before or after them, defaults, the Complex domain and assumptions, with small limits; each input runs on a thread with just the stack its recursion limit allows (`STACK_PER_LEVEL` a level). No panic (with `-a`, the debug assertions check the natives' shortcuts against their full computation), no stack overflow, no time-budget overrun; a result evaluates to itself; a trace at the Verbose level changes no result, is a consistent tree, and replays to the result; verification and the warnings' messages never panic. The limits bound the time only of programs that do not grow their work at each level (the evaluator has no budget per call, CAS-DOC-4 R8), so a right-hand side calls a user head at most once, applies nothing to parts (`Map`, `ReplaceAll`, `Function`) and uses each pattern variable at most once (`k[n___] := k[n, n, n]` triples its arguments), and a rule on a `Flat` head has no tests or conditions (a test evaluates the new group a blank binds, which matches again with a budget of its own); there is no `ReplaceRepeated` (a rule that matches what it makes grows the expression at each round); the number bits are drawn from 2⁶ to 2¹⁴ |
| `special` | The calculus of the special functions (CAS-361): antiderivatives of Laurent polynomials times exponentials, sines, cosines and hyperbolic functions, of polynomials times Gaussians and of `x^m/log x` differentiate back to their integrands, by balls at a point; rational sums from 1 to `n` through `PolyGamma` match direct summation for `n ≤ 8`; series of `Ei`, `Si`, `Ci`, `Erf`, `PolyGamma`, `Gamma` and `LambertW` to the order 4 match the function `2^-6` from their point |
| `sets` | Set operations (CAS-367) on random finite sets and lists of small integers and a few symbols: `Union` and `Intersection` commute and are idempotent, `|A ∪ B| = |A| + |B| − |A ∩ B|`, `A ∖ B` is included in `A` and meets `B` nowhere, De Morgan's laws hold relative to `A ∪ B`, and `IncludedQ` and `Element`, where decided, agree with the union and with each other |
| `random` | Random numbers (CAS-374) from random seeds: each draw of `Random`, `RandomInteger`, `RandomSample`, `RandomPermutation`, the variates and `RandomVector`, `RandomMatrix` and `RandomPolynomial` lies in its range or its law's support, with distinct places where they must be; the same seed draws the same again; the generator's naturals below a bound of up to 4096 bits stay below it |
| `lang_lex` | The cas-lang lexer on arbitrary text (bytes read as lossy UTF-8): the tokens cover the text exactly, each lexes alone as one token of the same kind, every error token has its error, errors are ordered and on character boundaries; the values of number, string and superscript tokens and the names of identifiers never panic (number exponents are bounded) |
| `lang_parse` | The cas-lang parser on arbitrary text: the tree's text is the input, every node's children cover it exactly and in order, errors are ordered, within the text and on character boundaries, and nesting beyond `MAX_DEPTH` stops cleanly instead of overflowing the stack |
| `lang_lower` | Parsing and lowering of arbitrary text (seeded with the `lang_parse` corpus): lowering never panics, every diagnostic, fix and statement span lies within the text, and every lowered statement prints in full form |
| `lang_session` | The cas-lang session on arbitrary text (seeded with the `lang_lower` corpus and the session scripts), with small limits and a watchdog that interrupts evaluation after a second and every 10 ms after: no panic, every output, message and fix within the text on character boundaries, and every result shown in input form parses and lowers without errors |
| `lang_print` | The input form on what lowering makes of arbitrary text (seeded with the `lang_parse` corpus): every expression that prints reads back without errors or notes as itself, up to the arithmetic the text writes out (`fold_arithmetic`), and printing what it reads back is stable |

Inputs are decoded from the fuzzer's bytes (`src/lib.rs`): limbs mix boundary values (`0`,
`u64::MAX`, `2^63`, low masks) with raw bytes, and lengths are weighted towards the algorithm
thresholds (Karatsuba, Toom-3, divide-and-conquer division, REDC).

## Running

Fuzzing runs on GitHub Actions only (`.github/workflows/fuzz.yml`) — never on the
development machine (CAS-158). Every push fuzzes every target, one job per target in
parallel, ten minutes each (CAS-159); a new push to the same ref cancels the previous run:

```sh
cargo fuzz run -O <target> corpus/<target> -- \
    -max_total_time=$FUZZ_TIME -rss_limit_mb=2048 -malloc_limit_mb=1024
```

The `asan` targets run with
`ASAN_OPTIONS=malloc_context_size=5:quarantine_size_mb=32:allocator_release_to_os_interval_ms=500`
(CAS-162): ASan keeps every distinct allocation stack for the life of the process, and the
language targets' recursion makes new ones without end; with the default 30 frames,
`lang_print`'s corpus alone took 1 GB (36 MB without ASan) and fuzzing it reached the RSS limit
within minutes.

The corpus is carried between runs with the actions cache, one cache per target, and
crash artifacts are uploaded with the failed run. A longer or narrower run is
`gh workflow run fuzz.yml -f time_per_target=<seconds> [-f targets="<names>"]`.

The ball function targets (`ball_fns`, `special_fns`, `complex_fns`, `machine_fns`, `ball_consts`, CAS-113),
`eval` (CAS-60), `poly_factor` (CAS-182), `poly_sparse` (CAS-155), `poly_gcd` (CAS-192) and `poly_integrate` (CAS-197) run with debug assertions and without AddressSanitizer (the workflow's
`assertions` job, `-O -a -s none`, with `-timeout=60`): cas-num and cas-core have no `unsafe`
code, ASan's quarantine more than doubles the resident memory, and the debug assertions check
the `Float` invariants of every value, and the natives' shortcuts against their full
computation. The other targets are the `asan` job, built `-O` with AddressSanitizer.

They fill the constant cache and the function tables once (`warm_up` in `src/ball.rs`), beyond
every precision an input asks for, so that a result never depends on earlier inputs;
`ball_consts` clears the cache before each input instead. A sample point longer than
`2 prec + 256` bits (an end across a huge exponent gap) is rounded toward the midpoint, which keeps
it in the ball and MPFR's work bounded.

A stack overflow on `eval`'s evaluation thread ends the process with SIGSEGV and no artifact,
libFuzzer's handler having no stack to run on; replay the corpus in batches
(`…/release/eval -runs=1 corpus/eval/<files>`), and the last `Running:` line before the crash
names the input.

A finding is saved under `artifacts/<target>/` and uploaded with the failed run; fetch it
with `gh run download`. Replaying that one input
(`cargo fuzz run -O <target> artifacts/<target>/crash-…`) is a debug run of a single case,
not a fuzz run, and is the one exception allowed locally. Turn a fixed finding into a
regression test in the owning crate's test suite, which CI then runs. The `corpus/`,
`artifacts/` and `target/` directories are not committed.
