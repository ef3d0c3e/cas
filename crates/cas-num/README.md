# cas-num

Arbitrary-precision integers, rationals, modular and ball arithmetic for the CAS engine.
Correctness is the primary goal. The design is in Lific page `CAS-DOC-2`.

## Layout

| Module | Contents |
|---|---|
| `kernel` | GMP `mpn`-style limb kernels on slices: add/sub (`_n`, `_1`, mixed length, in-place), compare, shifts, bitwise ops, bit counts. No allocation, no normalisation; carries and borrows are returned. |
| `kernel` (multiplication) | `mul_1`, `addmul_1`, schoolbook `mul_basecase`/`sqr_basecase`, and `mul`/`sqr`, which pick schoolbook, Karatsuba or Toom-3 by size and cut very unbalanced operands into balanced chunks. The recursive algorithms work in caller-provided scratch (`mul_scratch_len`, `sqr_scratch_len`). Thresholds are in `kernel/mul.rs`; retune with `cargo test --release -p cas-num tune_mul -- --ignored --nocapture`. |
| `kernel` (division) | `div_rem_1`/`rem_1` (Möller–Granlund reciprocal), `divexact_1` (Hensel), and `div_rem`: Knuth's algorithm D with a 3-by-2 quotient estimate, or Burnikel–Ziegler divide and conquer from `DC_DIV_THRESHOLD` limbs, in caller-provided scratch (`div_rem_scratch_len`). The word-level primitives are also compiled for 8- and 16-bit words in tests, where every boundary case is reachable. Retune with `cargo test --release -p cas-num tune_div -- --ignored --nocapture`. |
| `kernel` (gcd) | `gcd_1`: binary GCD of two limbs. |
| `Natural` | Unsigned arbitrary precision integer, limbs stored inline up to 192 bits and on the heap beyond. `std` operators in all owned/borrowed forms (`+ - * / %`, bitwise), `checked_sub*`, `square`, `div_rem`, `div_exact`, `is_divisible_by`, single-limb variants, shifts by `u64`, bit queries, primitive conversions. `gcd` and `lcm` by Lehmer's algorithm (Knuth's Algorithm L on the leading 63 bits, with a full division step when the leading bits cannot decide a quotient), finished by the binary GCD once the smaller operand fits a limb. |
| `Natural` (radix) | `to_string_radix`/`from_str_radix` for radix 2–36, `Display`/`Debug` (decimal), `Binary`, `Octal`, `FromStr`. Power-of-two radices pack bits in linear time. Other radices work in chunks of the largest power of the radix that fits a limb; from `TO_STRING_DC_THRESHOLD` limbs (or `FROM_STRING_DC_THRESHOLD` digits) the number splits by divide and conquer on `base^(2^i)`. Parsing follows the primitive integers' rules: an optional sign, no whitespace, separators or prefix, and `ParseNumberError` with `ParseErrorKind::{Empty, InvalidDigit}`. |
| `Natural` (powers and roots) | `pow` (binary exponentiation of the odd part, factors of two restored by one shift) and `pow_mod` (through `ModWord` for a one-limb modulus, `ModRing` otherwise). `sqrt_rem`/`isqrt` by Zimmermann's Karatsuba square root (Brent–Zimmermann Algorithm 1.12), `u128::isqrt` below 129 bits. `nth_root`/`root_rem` by Newton's iteration from the root of the top half, so each level needs a couple of steps. `is_perfect_square` filters by residues mod 64 and the odd primes to 47 before taking the root. `perfect_power` gives `(a, k)` with the largest `k` (prime exponents in turn, filtered by trailing zeros and by `p`th-power residues mod three primes `q ≡ 1 (mod p)`); `is_perfect_power` follows GMP (0 and 1 are). |
| `Integer` | Signed integer, layout E: `Small(i64)` inline (never allocates) or a boxed sign and `Natural` magnitude, canonical (boxed iff outside `i64`), 16 bytes. Inline operands use `checked_*` with an `i128` fallback; results are demoted when they fit. `+ - * / %` in all forms (`/`, `%` truncate), `checked_div_rem` with a `Rounding` (`Trunc`, `Floor`, `Ceil`, `Euclid`) and named `div_rem_*` variants, `div_exact`, `is_divisible_by`, sign queries, conversions. `gcd`, `lcm`, `gcd_ext` (the same cofactors as GMP's `mpz_gcdext`) and `mod_inverse` (in `0..|m|`). `pow`, `pow_mod` (in `0..|m|`; a negative exponent raises the inverse), `isqrt`/`sqrt_rem` (panic on a negative), `nth_root`/`root_rem` truncated towards zero (odd degrees for negatives), `is_perfect_square`, `is_perfect_power` (GMP's convention: a negative value needs an odd exponent). Radix conversion as for `Natural`, with a leading `-`; `-0` parses as zero. |
| `Rational` | Canonical fraction of two `Integer`s: positive denominator, coprime parts, zero is `0/1`, so derived `Eq`/`Hash` are by value. `new`/`checked_new` reduce; `+ - * /` in all forms, with Henrici's cancellation (Knuth §4.5.1) so results come out reduced without full-size GCDs; `checked_div`, `recip`, `abs`, `floor`, `ceil`, `trunc`, `round` (halves away from zero), `round_ties_even`, `pow` with an `i64` exponent. Comparison decides by sign, equal denominators or a bit-length estimate before cross-multiplying. Radix conversion as `num/den` (or `num` for integers); parsing reduces, accepts a sign on the numerator only, and reports `ParseErrorKind::ZeroDenominator`. |
| `algebra` | Ring-object traits: `Ring` (`Elem`, `zero`, `one`, `add`, `sub`, `neg`, `mul`, `elem_integer`, `characteristic`; provided `square`, `pow`, `is_zero`, `is_one`), `EuclideanDomain` (`div_rem`, `normal_form`; provided Euclidean `gcd` and `gcd_ext` with normalised results) and `Field` (`inv`, provided `div`). Instances: `Integers` (Lehmer GCD), `Rationals`, `ModWord` and `ModRing` (rings), `PrimeField<M>` (`WordPrimeField`, `BigPrimeField`; construction checks primality). |
| `modular` | Ring objects holding a modulus and its constants, with plain element values. `ModWord`: any modulus that fits a limb, elements `u64` in `0..n`, products reduced by one 2-by-1 division with a precomputed reciprocal (Möller–Granlund); `add`, `sub`, `neg`, `mul`, `pow`, `inv`, `reduce*`. `ModRing`/`Residue`: any modulus; odd moduli use Montgomery form, reduced one limb at a time (`redc_1`) or, from `REDC_N_THRESHOLD` limbs, by two full multiplications with `-n^(-1) mod R` (`redc_n`); even moduli store plain values and divide. `pow` uses a sliding window of up to 6 bits. `crt` solves two congruences with any moduli (`None` if incompatible); `CrtBasis` reconstructs from residues modulo pairwise coprime word moduli by Garner's algorithm, in `0..M` or symmetric. `rational_reconstruction` (Wang's algorithm, bounds `2 N D < m`, default `N = D`). Retune with `cargo test --release -p cas-num tune_redc -- --ignored --nocapture`. |
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
Measurements are short (20 samples, 1 s) so that the whole sweep stays in minutes. For tuning a
threshold, run its group alone, or use the dedicated `tune_*` tests (`tune_mul`, `tune_div`,
`tune_redc`). Criterion keeps the previous run in `target/criterion` and reports changes against it.

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
- Modular tests build their own moduli: word moduli at the normalisation edge (`2^63 ± 1`, `u64::MAX`, powers of two), the largest 64-bit prime and small values; multi-limb moduli mostly odd, around `REDC_N_THRESHOLD`, and `2^k ± 1`; operands `0`, `1`, `m - 1`, `m + 1`. CRT moduli come from `gcd_pair` (shared factors) or are made pairwise coprime by dropping offenders.
- Primality tests use primes from GMP's `next_prime`, products of two primes and prime squares; factorisation tests multiply known primes of each size class (trial division, rho, BPSW) and compare the factor lists exactly, labels included.
- `integer()`: a mix of `signed()`, `signed_near_i64_bounds()` and any `i64`.
