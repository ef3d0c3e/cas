# Exact integer functions

The functions of exact integers that the CAS evaluates exactly, beyond the
arithmetic of `Natural` and `Integer`.

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
