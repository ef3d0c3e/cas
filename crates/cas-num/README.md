# cas-num

Arbitrary-precision integers, rationals, modular and ball arithmetic for the CAS engine.
Correctness is the primary goal. The design is in Lific page `CAS-DOC-2`.

## Layout

| Module | Contents |
|---|---|
| `kernel` | GMP `mpn`-style limb kernels on slices: add/sub (`_n`, `_1`, mixed length, in-place), compare, shifts, bitwise ops, bit counts. No allocation, no normalisation; carries and borrows are returned. |
| `kernel` (multiplication) | `mul_1`, `addmul_1`, schoolbook `mul_basecase`/`sqr_basecase` (the square doubles its off-diagonal triangle and adds the diagonal in one pass, as GMP's `mpn_sqr_diag_addlsh1`, and is unrolled up to 6 limbs), and `mul`/`sqr`, which pick schoolbook, Karatsuba, Toom-3 or a three-prime NTT (`kernel/ntt.rs`) by size and cut very unbalanced operands into balanced chunks. Products modulo `B^m - 1` (`mul_wrap`, for Newton division and Montgomery reduction) take a cyclic transform at transform lengths from `mul_wrap_ntt_threshold`; below it, from `mul_wrap_split_threshold`, they split `B^m - 1 = (B^(m/2) - 1)(B^(m/2) + 1)` into a recursive half and a half-size product, joined by the CRT (GMP's `mpn_mulmod_bnm1`). The recursive algorithms work in caller-provided scratch (`mul_scratch_len`, `sqr_scratch_len`). Thresholds come from the tuning configuration (see below). |
| `kernel` (division) | `div_rem_1`/`rem_1` (Möller–Granlund reciprocal), `divexact_1` (Hensel), and `div_rem`: Knuth's algorithm D with a 3-by-2 quotient estimate, or Burnikel–Ziegler divide and conquer from `DC_DIV_THRESHOLD` limbs, in caller-provided scratch (`div_rem_scratch_len`). Within the crate, `divappr_q` gives the quotient or one more without the remainder, as GMP's `mpn_divappr_q`: from the divisor's top limbs and a guard limb, with the last block's remainder skipped. The word-level primitives are also compiled for 8- and 16-bit words in tests, where every boundary case is reachable. |
| `kernel` (approximate) | For ball midpoints, which need a value near the true one and a bound on the distance, not correct rounding (CAS-DOC-5 §5.2). Each returns a lower bound `r` with `r <= floor(t) <= r + e`, `e` from the matching `*_error(n)`. `mulhi`: the top `n` limbs of an `n`-by-`n` product (schoolbook over the columns from `B^(n - 2)` up, two rows at a time; Mulders' divide and conquer with `mullo`'s split and thresholds, error `2 e(n1) + 5`; the full product, exact, up to 2 limbs and from `MULLO_MUL_THRESHOLD`). `recip_approx`: `floor(B^(2n) / D)` within 1 (2 with Newton), division's reciprocal. `div_approx`: `floor(A B^n / D)`, exact below `MU_DIV_THRESHOLD`, else `A + mulhi(A, I)` from the reciprocal. `rsqrt_approx`: `floor(sqrt(B^(3n) / A))` within 1, by a Newton recursion whose truncations all round down. `sqrt_approx`: `floor(sqrt(A B^n))`, exact through `Natural` (the Newton form `A + mulhi(A, J)` from the inverse root is tested but measured slower up to 8192 limbs, so not dispatched to yet). The proofs are the doc comments; the tests check each bound against GMP. |
| `kernel` (gcd) | `gcd_1`: binary GCD of two limbs. |
| `kernel` (square root) | `sqrtrem`: Zimmermann's Karatsuba square root on limbs, as GMP's `mpn_sqrtrem`. The operand is shifted once to an even length and a top limb of at least `B/4`; each level takes the root of the top half recursively, then one division by that root, normalised as it is, gives the low half. `u128::isqrt` is the base case. `sqrt` is the root alone, as GMP's `mpn_dc_sqrt`: one such step with a quotient-only division and no square, deciding from 32 or more guard bits below the root; it also says whether the operand is a square. Caller-provided scratch (`sqrtrem_scratch_len`, `sqrt_scratch_len`). |
| `Natural` | Unsigned arbitrary precision integer, limbs stored inline up to 256 bits and on the heap beyond. `std` operators in all owned/borrowed forms (`+ - * / %`, bitwise), `checked_sub*`, `square`, `div_rem`, `div_exact` (GMP's contract: the caller guarantees the division is exact; otherwise the result is unspecified, and debug builds panic), `is_divisible_by`, single-limb variants, shifts by `u64`, bit queries, primitive conversions. `gcd` and `lcm` by Lehmer's algorithm (Knuth's Algorithm L on the leading 128 bits, with a full division step when the leading bits cannot decide a quotient), finished by the binary GCD once the smaller operand fits a limb; above `gcd_dc_threshold` limbs by subquadratic half-GCD steps (Möller's algorithm, as in GMP's `mpn_hgcd`), which also serve the extended GCD. The extended GCD tracks one cofactor as GMP's `mpn_gcdext` does, as two magnitudes and a sign: one fused pass per Lehmer matrix, half-GCD matrices applied by products (or, while the cofactors are shorter than half the operands, inside the half-GCD's top level, which then builds no matrix), and the Lehmer steps after half-GCD steps on fresh cofactors joined at the end by two products. |
| `Natural` (radix) | `to_string_radix`/`from_str_radix` for radix 2–36, `Display`/`Debug` (decimal), `Binary`, `Octal`, `FromStr`. Power-of-two radices pack bits in linear time. Other radices work in chunks of the largest power of the radix that fits a limb; from `TO_STRING_DC_THRESHOLD` limbs (or `FROM_STRING_DC_THRESHOLD` digits) the number splits by divide and conquer on `base^(2^i)`. Parsing follows the primitive integers' rules: an optional sign, no whitespace, separators or prefix, and `ParseNumberError` with `ParseErrorKind::{Empty, InvalidDigit}`. |
| `Natural` (powers and roots) | `pow` (binary exponentiation of the odd part, factors of two restored by one shift) and `pow_mod` (through `ModWord` for a one-limb modulus, `ModRing` otherwise). `sqrt_rem` by the kernel's `sqrtrem` (Brent–Zimmermann Algorithm 1.12), `isqrt` by its root-only `sqrt`. `nth_root`/`root_rem` on limbs by one Newton step per doubling of the precision, as GMP's `mpn_rootrem`: a quotient-only division by `k s^(k-1)` gives the next bits of the root, and the powers `s^(k-1)` and `s^k` check them and give the next step's remainder; at large sizes `s^k` goes modulo `B^m - 1` for `m` just above the remainder's length. `nth_root` computes 32 guard bits of the root instead of the last step's powers, from roots of 256 bits. `is_perfect_square` filters by residues mod 64 and the odd primes to 47 before taking the root and remainder: what passes is almost always a square, for which the root-only `sqrt` would cost more. `perfect_power` gives `(a, k)` with the largest `k` (prime exponents in turn, filtered by trailing zeros and by `p`th-power residues mod three primes `q ≡ 1 (mod p)`); `is_perfect_power` follows GMP (0 and 1 are). |
| `Integer` | Signed integer, layout E: `Small(i64)` inline (never allocates) or a boxed sign and `Natural` magnitude, canonical (boxed iff outside `i64`), 16 bytes. Inline operands use `checked_*` with an `i128` fallback; results are demoted when they fit. `+ - * / %` in all forms (`/`, `%` truncate), `checked_div_rem` with a `Rounding` (`Trunc`, `Floor`, `Ceil`, `Euclid`) and named `div_rem_*` variants, `div_exact`, `is_divisible_by`, sign queries, conversions. `gcd`, `lcm`, `gcd_ext` (the same cofactors as GMP's `mpz_gcdext`) and `mod_inverse` (in `0..|m|`). `pow`, `pow_mod` (in `0..|m|`; a negative exponent raises the inverse), `isqrt`/`sqrt_rem` (panic on a negative), `nth_root`/`root_rem` truncated towards zero (odd degrees for negatives), `is_perfect_square`, `is_perfect_power` (GMP's convention: a negative value needs an odd exponent). Radix conversion as for `Natural`, with a leading `-`; `-0` parses as zero. |
| `Rational` | Canonical fraction of two `Integer`s: positive denominator, coprime parts, zero is `0/1`, so derived `Eq`/`Hash` are by value. `new`/`checked_new` reduce; `+ - * /` in all forms, with Henrici's cancellation (Knuth §4.5.1) so results come out reduced without full-size GCDs; `checked_div`, `recip`, `abs`, `floor`, `ceil`, `trunc`, `round` (halves away from zero), `round_ties_even`, `pow` with an `i64` exponent. Comparison decides by sign, equal denominators or a bit-length estimate before cross-multiplying. Radix conversion as `num/den` (or `num` for integers); parsing reduces, accepts a sign on the numerator only, and reports `ParseErrorKind::ZeroDenominator`. |
| `Complex` | Complex rationals ℚ(i): canonical `Rational` real and imaginary parts, so `Eq`/`Hash` are by value (a zero imaginary part is allowed; `is_real`). `+ - * /` in all forms (products skip the cross terms when a factor is real), `checked_div`, `conj`, `norm`, `checked_recip`/`recip`, `pow` with an `i64` exponent. `Display`/`FromStr` as `a + b*I` (`I`, `-1/2*I`, `1/2 - 3/4*I`). Ring objects in `algebra`: `GaussianRationals` (field) and `GaussianIntegers` (ℤ[i], Euclidean: quotient parts rounded to nearest, so `norm(r) <= norm(b)/2`; normal forms in the first quadrant). |
| `ball::Float` | Binary float `±0.m × 2^exp` (CAS-DOC-5 §3): mantissa top-aligned in limbs (top bit set, no zero low limb, inline up to 4 limbs), `Integer` exponent (never overflows), finite values only, so `Eq`/`Hash` are structural. Exact `add`, `sub`, `mul`, `neg`, `abs`, `mul_2exp`, `cmp`/`cmp_abs`; correctly rounded `round`, `add_round`, `sub_round`, `mul_round`, `div_round`, `sqrt_round` (with `checked_*` forms returning `None` on division by zero or a negative root), `from_rational` in every `Round` mode (`Down`/`Up` toward/away from zero, `Floor`, `Ceil`, `Nearest` ties to even), returning MPFR's ternary value as an `Ordering`. A rounded sum never builds more than `max(prec, mantissa) + 4` bits: an operand far below the other is replaced by a single bit that rounds alike. One-limb operands with `i64` exponents take a `u128` path, other sizes stack buffers up to 16 limbs and `i64` exponent bookkeeping, with the `Integer` path beyond. For ball midpoints, `div_approx`, `sqrt_approx` and `rsqrt_approx` return a lower bound of the magnitude within `APPROX_ULPS` (2) units in the last place, through the approximate kernels at large sizes (exact, then truncated, below their Newton thresholds). Conversions: `from_integer_exp`/`to_integer_exp`, `to_rational` (exact), `to_f64` (nearest, IEEE overflow and gradual underflow), `from_f64` (exact). |
| `ball::Mag` | Radius type (CAS-DOC-5 §4, Arb's `mag_t`): `man × 2^(exp - 30)` with a 30-bit mantissa, an `Integer` exponent, plus zero and `+∞`; 24 bytes. Every operation returns an upper bound (`*_lower`: a lower bound) within `2^-29` relative for one rounding: `add`, `sub` (of `max(a - b, 0)`), `mul`, `div`, `mul_2exp`, `pow` (within `(1 + 2^-29)^(2n)`), `exp`, `expm1`, `log1p` (fixed-point `u128` series rounded up, with explicit tail bounds; `exp` squares back at 128 bits, tight up to `2^62`, `+∞` above), conversions `from_float`/`from_integer` (and `_lower`), `to_float` (exact). Exponent bookkeeping goes through `i64` when it fits: about 6.5–8 ns per basic operation, 55–135 ns for `exp`/`log1p`. |
| `ball::Ball` | Midpoint ± radius (`Float` ± `Mag`, CAS-DOC-5 §5): every result contains the exact result at every choice of points in the inputs. Precision is passed to each operation; midpoints are rounded toward zero (or come from the approximate entries) and the error goes into the radius, so exact results stay exact (on the Newton paths too: the candidates `q`, `q + u`, `q + 2u` are checked by one exact product). `add`, `sub`, `mul`, `div`, `inv`, `sqrt`, `pow` (integer exponent), `abs`, `neg`, `mul_2exp`, `*_assign`; division by a ball containing 0 and the square root of a ball reaching below 0 give `Ball::INDETERMINATE` (`[0 ± ∞]`). Predicates: `contains` (`Float`, `contains_integer`, `contains_rational`, `contains_ball`), `overlaps`, `is_exact`, `is_finite`, `sign -> Option<Ordering>`, `cmp -> Option<Ordering>`, decided at 64 bits when possible and exactly otherwise; `union`, `intersection`, `rel_accuracy_bits`, `round_to_accuracy`; `from_integer` (exact), `from_rational` (½ ulp). `RealBalls::new(prec)` implements `Ring` only (no `EuclideanDomain`/`Field`: a ball containing 0 isn't zero). |
| `ball` constants | `Ball::const_pi`, `const_log2`, `const_e` (CAS-DOC-5 §7.1): binary splitting (`ball::bsplit`, crate-private: `P`, `Q`, `B`, `T` over `Integer`, factors that are always 1 skipped, `P` not formed on the rightmost path) of Chudnovsky's series for π, `18 atanh(1/26) − 2 atanh(1/4801) + 8 atanh(1/8749)` for log 2 and `Σ 1/k!` for e, at `prec + 32` bits, with each series' tail bounded rigorously (proofs in the doc comments, checked against the true tails) and added to the radius. One process-wide cache per constant (`RwLock`): a request at or below the cached precision rounds the cached ball (radius about ½ ulp), a higher one recomputes at `max(prec, 2·cached)`; `ball::clear_constant_cache` frees it. Against MPFR from scratch: π 1.2–1.6× at 100–1000 digits, 0.7× from 10^4 to 10^6; log 2 0.8–1.24×; e 0.24–0.9× (MPFR's `exp(1)`). |
| `ball::MachineBall` | `f64` midpoint and radius, 16 B, `Copy` (CAS-DOC-5 §5.4): `add`, `sub`, `mul`, `div`, `inv`, `sqrt`, `pow(i32)`, `neg`, `abs` with the midpoint error added exactly where an error-free transformation gives it (two-sum; two-product by FMA where the target has one, Dekker's split otherwise) and as half an ulp otherwise, and every radius step rounded up by one ulp (so exact results stay exact). Results return `None` when the midpoint overflows, underflows or comes out subnormal, or the radius overflows, and for division by a ball containing 0 or the root of one reaching below 0: the caller redoes the operation with `Ball`. A radius below the normal range is rounded up to `f64::MIN_POSITIVE`. Exact predicates `contains`, `contains_ball`, `overlaps`, `sign`, `cmp` (rounded bounds first, exact `Float`s when they cannot decide); `to_ball` (exact midpoint), `from_ball` (nearest midpoint, radius grown by the rounding), `via_ball` (a `Ball` function at 64 bits: the fallback for functions until CAS-88). About 2.4 / 5.7 / 7.4 / 8.2 ns for add / mul / div / sqrt, against 33 / 52 / 90 / 161 ns for 53-bit balls. |
| `algebra` | Ring-object traits: `Ring` (`Elem`, `zero`, `one`, `add`, `sub`, `neg`, `mul`, `elem_integer`, `characteristic`; provided `square`, `pow`, `is_zero`, `is_one`), `EuclideanDomain` (`div_rem`, `normal_form`; provided Euclidean `gcd` and `gcd_ext` with normalised results) and `Field` (`inv`, provided `div`). Instances: `Integers` (Lehmer and half-GCD), `Rationals`, `GaussianIntegers`, `GaussianRationals`, `ModWord` and `ModRing` (rings), `PrimeField<M>` (`WordPrimeField`, `BigPrimeField`; construction checks primality). |
| `modular` | Ring objects holding a modulus and its constants, with plain element values. `ModWord`: any modulus that fits a limb, elements `u64` in `0..n`, products reduced by one 2-by-1 division with a precomputed reciprocal (Möller–Granlund); `add`, `sub`, `neg`, `mul`, `pow`, `inv`, `reduce*`. `ModRing`/`Residue`: any modulus; odd moduli use Montgomery form, reduced two limbs at a time (`redc_2`) or, from `REDC_N_THRESHOLD` limbs, by `redc_n`: a short product with `-n^(-1) mod R` (from the 2-adic Newton inverse) and the product `q n` modulo `B^m - 1`. Below the cyclic transform that product splits `B^m - 1` (GMP's `mulmod_bnm1`); at transform sizes both products use transforms of `n` and `-n^(-1)` prepared once per ring, two transforms each instead of three. Even moduli store plain values and divide. `pow` uses a sliding window whose width grows with the exponent, up to 16 bits. `crt` solves two congruences with any moduli (`None` if incompatible); `CrtBasis` reconstructs from residues modulo pairwise coprime word moduli by Garner's algorithm, in `0..M` or symmetric. `rational_reconstruction` (Wang's algorithm, bounds `2 N D < m`, default `N = D`). |
| `prime` | `Natural::primality() -> Primality` (`Prime` proven, `ProbablePrime`, `NotPrime`): table lookup below 2^16; trial division by the primes below 2^12 (one remainder per group of primes whose product fits a limb); below 2^64 Miller–Rabin with Sinclair's seven bases (deterministic); above, Baillie–PSW (base-2 Miller–Rabin plus a strong Lucas test with Selfridge's parameters). `is_strong_probable_prime(base)`, `next_prime`, `Integer::jacobi`. `factor()` / `factor_with_effort(rho_iterations) -> Factorization`: trial division by the primes below 2^16, perfect powers through their base, Brent's variant of Pollard's rho (`ModWord` or `ModRing`, GCD every 128 steps); factors are sorted, merged and labelled, and composites that resist the effort stay labelled `NotPrime`. |

## Benchmarks

`benches/num.rs` (criterion, with GMP via `rug` as the baseline) times every op twice, as
`op/ours/n` and `op/gmp/n`, on the same deterministic pseudo-random `n`-limb inputs:

```sh
cargo bench -p cas-num                              # everything, sizes 1, 2, 4 … 4096 limbs
cargo bench -p cas-num -- '^mul/'                   # one group
CAS_BENCH_MAX_LIMBS=100000 cargo bench -p cas-num -- '^(mul|sqr)/'   # the full sweep
```

Groups: `add`, `mul`, `sqr`, `div_rem` (2n by n limbs), `gcd`, `pow_mod` (odd modulus, full-size
exponent, up to 64 limbs), `primality` (an n-limb prime, up to 64 limbs), `to_string`, `from_str`, and
`small`: batches of 512 sums and products of values below 2^63, against GMP and plain `i128`.
Measurements are short (20 samples, 1 s) so that the whole sweep stays in minutes. Criterion keeps
the previous run in `target/criterion` and reports changes against it.

`benches/ball.rs` does the same for binary floats and balls against MPFR (`rug::Float`), as
`op/ours/bits` and `op/mpfr/bits`: groups `add`, `mul`, `div`, `sqrt` at 64 … `CAS_BENCH_MAX_BITS`
(default 4096) bits, and `pi`, `log2`, `e`, `exp`, `log`, `sin`, `atan` at 100 …
`CAS_BENCH_MAX_DIGITS` (default 10^4) digits, with `op/approx/bits` (midpoint entries) and
`op/ball/bits` rows where they apply; the functions have MPFR rows only until their tickets land.
Filter with `cargo bench -p cas-num --bench ball -- '^mul/'`. Its `mulhi` group times `mulhi`, `mullo` and `mul`
at 1 … `CAS_BENCH_MAX_LIMBS` limbs.

## Long runs

Stress tests, long proptest runs, fuzzing and full benchmark sweeps go through
[`scripts/capped.sh`](../../scripts/capped.sh), which runs the command in a systemd scope with a
memory limit (`CAP_MEM`, default 8G, swap disabled), a timeout (`CAP_TIME`, default 3600 s) and
`CARGO_BUILD_JOBS` capped (`CAP_JOBS`, default 6). A runaway (a binary-splitting tree or a
precision loop gone wrong) is then killed instead of exhausting a shared machine:

```sh
scripts/capped.sh cargo test -p cas-num
PROPTEST_CASES=20000 CAP_MEM=4G scripts/capped.sh cargo test -p cas-num --release float
```

## Tuning configuration

The thresholds that pick an algorithm by operand size (multiplication, squaring, NTT, division,
Montgomery reduction, GCD, radix conversion) are compiled in from
[`tuning/default.toml`](tuning/default.toml), measured on the development machine. To build with
other values, point `CAS_NUM_TUNING` at another file, e.g. in `.cargo/config.toml`:

```toml
[env]
CAS_NUM_TUNING = { value = "cas-num-tuning.toml", relative = true }
```

`build.rs` rejects unknown, missing or repeated keys and any value below an algorithm's minimum,
and every accepted value gives correct results. Two stress configurations prove that:
`tuning/stress-small.toml` (every algorithm at tiny sizes) and `tuning/stress-large.toml` (only the
basecases); run the suite under each with
`CAS_NUM_TUNING=tuning/stress-small.toml cargo test -p cas-num --target-dir target/tuning-small`.
To measure the thresholds on your machine, run the auto-tuner on an idle machine:

```sh
cargo run --release -p cas-tune -- --out my-tuning.toml          # about 12 minutes
cargo run --release -p cas-tune -- --quick --out my-tuning.toml  # a coarser pass
cargo run --release -p cas-tune -- --only dc_div_threshold,redc_n_threshold --out my-tuning.toml
```

It rebuilds a small timing binary for every candidate value (inlining and code layout depend on
the value, so each is measured as compiled), times the operation the threshold controls over a
window of sizes, and keeps the candidate with the lowest mean time relative to the best at each
size. It writes the configuration and a report (`my-tuning.report.txt`) of every candidate's cost.
Timings need a quiet machine: the tuner warns when the load average shows other work and when
identical runs disagree by more than 3%, and the report gives the measured noise. The optima are
often flat (candidates within 1-2% of each other), so values that move between runs are expected.

The ignored `tune_*` tests (`tune_mul`, `tune_ntt`, `tune_div`, `tune_redc`, `tune_to_string`,
`tune_approx`) print direct algorithm-against-algorithm measurements.

## Fuzzing

libFuzzer targets live in [`/fuzz`](../../fuzz/README.md), outside the workspace (cargo-fuzz needs
nightly): `cd fuzz && cargo fuzz run -O natural`. They cover the same operations against GMP with
fuzzer-chosen inputs.

## Differential testing against GMP

Every operation is checked against GMP on generated inputs. The harness lives in
`src/testkit/` and is compiled only under `cfg(test)`. `rug` and `proptest` are
dev-dependencies, so they never end up in a runtime build.

### Pieces

| Item | Purpose |
|---|---|
| `testkit::Mag`, `testkit::Signed` | Normalised little-endian `u64` limb vectors (sign + magnitude). This is the neutral exchange format between our types and GMP. `Debug` prints hex with `_` at limb boundaries. |
| `ToRug` | Maps a value of ours to its GMP equivalent. It is implemented for primitives, `bool`, `Ordering`, `String`, `Option<T>`, pairs and triples, `Mag`, `Signed`, `Natural`, `Integer` and `Rational` (whose impls also assert the canonical representation). Every cas-num type implements it as it lands. |
| `all_forms!`, `all_forms_prim!` | Evaluate an operator in every owned/borrowed and in-place form, fail unless all agree, and yield the result. |
| `check_same(ours, gmp)` | Fails the proptest case unless the two values agree, and prints both. |
| `differential!` | Declares one `#[test]` per entry, asserting that `ours` equals `gmp`. |
| `testkit::strategy` | Adversarial generators (see below). |
| `testkit::float` | The MPFR oracle for floats and balls (CAS-DOC-5 §11): `Dyadic` (`±m·2^e`, odd mantissa, exponent as a `Signed`, so beyond MPFR's `i32` range too) with exact conversions from and to `rug::Float` and `rug::Rational`; `bracket(prec, f)`, the true value of `f` enclosed by MPFR at `prec + GUARD` (64) bits rounded toward −∞ and +∞; `check_encloses(mid, rad, bracket)`, which requires a ball to meet the bracket (the only thing a bracket can refute, so callers compute it well above the ball's accuracy); `check_rounded`, exact equality with MPFR including the ternary value; strategies `dyadic`, `dyadic_with`, `far_apart_pair` and `huge_shift` (exponents near and beyond `±2^63`, tested by scaling: `x·2^E op y·2^E` against `(x op y)·2^E`). Note that rug's `Round::Down`/`Up` are toward −∞/+∞; rug's `Zero`/`AwayZero` are our `Down`/`Up`. |

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
- `signed_near_i64_bounds()`: within ±3 of `0`, `±2^31`, `±2^32`, `±2^62`, `±2^63`, `±2^64`,
  `±2^65`, `±2^127`, `±2^128`, the promotion edges of `Integer`'s inline fast path.
- `gcd_pair(max_limbs)`: pairs for GCD tests: independent values, values with a large common
  factor, consecutive Fibonacci numbers, a value and a multiple of it, equal values, and values
  sharing their leading limbs.
- `rational()`, `rational_pair()`: canonical fractions (reduced by GMP), including integers, small
  fractions and denominators at the `i64` edges; pairs whose denominators share a factor, or whose
  numerators share a factor with the other's denominator, plus equal and opposite pairs.
- Root and power tests build their own inputs near exact powers: `s^n - 1`, `s^n`, `s^n + 1`, `s^2 + 2s` (the largest square-root remainder) and `(s + 1)^n`; `a^k` with factors of two; `2^t ± 1`.
- Modular tests build their own moduli: word moduli at the normalisation edge (`2^63 ± 1`, `u64::MAX`, powers of two), the largest 64-bit prime and small values; multi-limb moduli mostly odd, around `REDC_N_THRESHOLD`, and `2^k ± 1`; operands `0`, `1`, `m - 1`, `m + 1`. `redc_choices_match` runs all-ones and `2^k ± 1` moduli at every length where REDC changes how it takes its products (the split, prepared transforms). CRT moduli come from `gcd_pair` (shared factors) or are made pairwise coprime by dropping offenders.
- Primality tests use primes from GMP's `next_prime`, products of two primes and prime squares; factorisation tests multiply known primes of each size class (trial division, rho, BPSW) and compare the factor lists exactly, labels included.
- `integer()`: a mix of `signed()`, `signed_near_i64_bounds()` and any `i64`.
