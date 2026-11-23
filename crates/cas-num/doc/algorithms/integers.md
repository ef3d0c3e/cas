# Exact integer and rational functions

The functions of exact integers and rationals that the CAS evaluates
exactly, beyond the arithmetic of `Natural`, `Integer` and `Rational`.

## Factorial (`Natural::factorial`)

`natural/factorial.rs` (CAS-123). `n!` exactly, for `Factorial` in cas-core.

**Base algorithm.** GMP's method for small `n` [[GMP]](references.md#gmp)
(chapter "Algorithms", subsection "Factorial"). The odd part of `n!` is the
product of the odd numbers up to `n`, times the odd part of `⌊n/2⌋!`, and so
on down:

```text
23! = (23·21·19·…·5·3)(11·9·7·5·3)(5·3) · 2^19
```

The power of two is `2^(n - s(n))`, where `s(n)` is the number of ones in the
binary form of `n` (Legendre's formula for the exponent of 2 in `n!`). It is
applied as one shift at the end.

**Tricks.**
- **Running product.** The odd products of the levels nest. The product of
  the odd numbers up to `⌊n/2^i⌋` is the one up to `⌊n/2^(i+1)⌋` times the
  odd numbers in between. So a running product, from the top level down,
  takes in each level's new odd numbers with one multiplication, and
  multiplies into the odd part with one more. The last of those
  multiplications has operands of about half the result's size each.
- **Product trees.** A level's new odd numbers are multiplied in a balanced
  tree, so that every product has operands of about the same size and runs
  through the subquadratic multiplications. Its leaves are limbs holding as
  many consecutive factors as fit in 64 bits.
- **One limb.** Up to `20!`, which fits in a limb, it is a product of `u64`.

GMP switches to Luschny's prime-swing algorithm for large `n`, through a
prime sieve. That is not done here. The odd-part method measured 1.6× GMP's
time at `n = 10^3` and `10^4`, 1.8× at `10^5` and 1.9× at `10^6`
(18.5·10^6 bits).

**Tests.** `natural/factorial/tests.rs`: against GMP (`rug`) for random `n`
up to 5000, weighted towards the one-limb results and the first levels, and
for `n = 3·10^4` and `10^5`.

## Simplest rational in an interval (`Rational::simplest_between`)

`rational/mod.rs` (CAS-148). The rational of least denominator in a closed
interval `[lo, hi]`, and of those the least in magnitude. The CAS's `exact`
uses it to turn a ball into the simplest fraction that matches the digits
it shows (CAS-DOC-5 §6.4, D4).

**Algorithm.** Continued fractions [[HW08]](references.md#hw08) (chapter X),
by the classic argument of the Stern–Brocot tree. For `0 < lo <= hi`:
- if the interval holds an integer, the least one, `⌈lo⌉`, is the answer;
- otherwise both ends have the same integer part `a`, and the answer is
  `a + 1/y`, with `y` the simplest rational in
  `[1/(hi - a), 1/(lo - a)]`.

An interval holding 0 gives 0, and a negative one the negated answer for
`[-hi, -lo]`. The terms are collected in a loop and the fraction evaluated
from the last term, so there is no recursion. Each term is a partial
quotient of both ends, so there are `O(log)` of them in their
denominators.

**Tests.** `rational/tests.rs`: a table (`[3/10, 4/10]` gives `1/3`; an
interval holding 0 gives 0; the numbers that round to `0.333333333333` at
12 digits give `1/3`, and those that round to `0.7` give `7/10`), and a
property against a search of the denominators in order.

## Primality proof (`Natural::prove_prime`)

`prime/proof.rs` (CAS-551). Whether `n` is prime, with a proof:
`Proven`, `Composite` or `Unknown`, for `is_prime` in cas-core, which
reports no prime it has not proven.

**Base algorithm.** Below 2^64 the deterministic Miller–Rabin test of
`Natural::primality` decides, as it does for every composite it finds
(trial division, Baillie–PSW). A probable prime `n` above 2^64 goes through
the `n − 1` test of Brillhart, Lehmer and Selfridge [[BLS75]](references.md#bls75):
- `n − 1` is factored by `Natural::factor_with_ecm` within the caller's
  limits (rho iterations, ECM curves and stage-1 bound; cas-core passes
  `IntegerFactorList`'s, which shrink as `(128/bits)²` above 128 bits).
  `F` is the product of the full powers of its primes `q` that are proven,
  by this test recursively for those above 2^64; the full power of each is
  divided out of `n − 1` directly, so `F` and `R = (n − 1)/F` are coprime
  even when an unsplit composite shares a prime.
- For each such `q`, Pocklington's condition: the least base `a ≥ 2` with
  `a^(n−1) ≡ 1 (mod n)` and `gcd(a^((n−1)/q) − 1, n) = 1`. Then every prime
  `p | n` has `p ≡ 1 (mod F)`. A base failing Fermat's test, or a gcd that
  is a proper factor, proves `n` composite; for a prime `n` a base fails
  with probability `1/q`, and after 1000 the answer is `Unknown`.
- `F² > n`: every prime factor of `n` exceeds `√n`, so `n` is prime
  (Pocklington's theorem, as [[BLS75]](references.md#bls75) states it).
- Otherwise Theorem 5: `F` is even and prime to `R`; with
  `R = 2Fs + r`, `1 ≤ r < 2F`, if `(F + 1)(2F² + (r − 1)F + 1) > n`,
  `n` is prime when `s = 0` or `r² − 8s` is not a square. This lowers the
  needed `F` from `√n` to about `(n/2)^(1/3)`.
- Anything else is `Unknown`, including a square `r² − 8s`, which
  [[BLS75]](references.md#bls75) makes composite: we answer only what the
  test checked.

**Limits.** An `n − 1` whose factored part stays below about `n^(1/3)`
within the limits is `Unknown`: a prime `n` with `n − 1 = 2·p·q`, `p` and
`q` beyond rho and ECM, is not proven. giac answers 1 there through PARI's
APRCL (when PARI does not fail); cas-core's `is_prime` leaves the call
with a warning. APRCL or ECPP would close the gap.

**Tests.** `prime/proof/tests.rs`: Mersenne primes `2^89 − 1`, `2^107 − 1`,
`2^127 − 1`, `10^30 + 57` and `2^64 + 13` proven; values below 2^64 (0, 1,
Fermat's `2^32 + 1`, strong pseudoprimes to the first bases, the prime
9856989898997789789) decided; strong pseudoprimes to the first twelve and thirteen prime
bases, `2^128 + 1` and a product of two primes near 2^40 composite; and
`n = 12·p·q + 1`, `p` and `q` primes near 2^40, `Unknown` with one rho
iteration and no curves, proven with the default limits.
