# Sparse multivariate polynomials

`src/sparse.rs` and `src/sparse/` (CAS-136).

## Representation

A `Sparse<E>` holds its terms in decreasing lexicographic order of their
exponent vectors, the first variable most significant, each coefficient
non-zero. The operations live on the ring object `SparseRing<R>` (`nvars`
variables over a ring object `R`), which is a `Ring`.

The coefficients are one vector. The exponents are another, of packed words,
as Monagan and Pearce do it [MP07]: dynamic arrays, without the indirection
of linked lists.

### Packed exponent vectors [MP07]

- **Fields.** Each variable gets a field of `b` bits, `b ≥ 8`, in 64-bit words,
  `⌊64 / b⌋` to a word. The first variable takes the top field of the first
  word.
- **Order.** Comparing two packed vectors word by word, as unsigned integers,
  compares them lexicographically.
- **Products.** Multiplying monomials adds their words.
- **Guard bits** [MP07]. The top bit of each field is a guard, clear in a
  packed vector, so a field holds exponents up to `2^(b−1) − 1`.
  - Two fields below `2^(b−1)` sum below `2^b`, so no carry leaves a field. A
    sum that sets a guard has overflowed.
  - For division, `(m | G) − d` is computed with the guards `G` set. A field
    of `d` above `m`'s clears its guard and borrows from no other, so `d`
    divides `m` exactly when every guard stays set, and clearing them gives
    `m − d`.
- **The canonical layout.** Every polynomial is kept in the narrowest layout
  (from 8 bits) that holds its largest exponent. Equal polynomials are then
  equal values, which derived equality and hashing need. Operands in
  different layouts are repacked to the wider before merging.

## Sums

Merging the two term lists, monomials compared as words; equal monomials
add, and zero sums drop out.

## Products: Johnson's heap [Joh74, MP07]

The products `aᵢ bⱼ` of a row `i` decrease with `j`. A binary heap merges the
rows:
- it holds one product per row, keyed by its monomial;
- the largest comes out, and the products of equal monomials add up as they
  come;
- each popped `(i, j)` is replaced by `(i, j+1)`.

The result comes out in order, and needs no sort.
- **Lazy rows** [MP07]. Row `i + 1` enters the heap when `(i, 0)` leaves it.
  The heap then holds at most one product per row of the shorter factor,
  which is taken as the first.
- **Slots.** Each product's monomial is computed when it enters the heap, and
  stored in a slot that its successor reuses. The storage is linear in the
  rows.
- **Fields.** Both factors are repacked first in fields wide enough for the
  sum of their largest exponents, so no product overflows. The guards still
  check every monomial added. Were one to overflow, the product would start
  over one bit wider; `packed_fields_at_their_boundaries` drives that path
  with fields made too narrow on purpose.
- `Cancel` is polled every 1024 monomials.

## Division [MP07, MP11]

`a = q b + r` by the leading term of `b` in lexicographic order.
- **The heap.** The terms of `a − q b` come out of a heap in decreasing
  order: `a`'s own, and those of `−qₖ bⱼ` for `j ≥ 1`. `qₖ b₀` cancels the
  term `qₖ` was made from, so it never enters. Each new quotient term enters
  as `(k, 1)`. This is the quotient heap of [MP07]; its "divisor heap", whose
  size the divisor bounds, would suit the gcd's cofactor divisions better,
  and comes later.
- **Each term.** A term whose monomial the leading monomial of `b` divides
  (by the guard test) and whose coefficient the leading coefficient divides
  exactly, in the coefficient ring, makes a quotient term. Any other goes to
  `r`.
  - Over a field this is the division with remainder by one polynomial.
  - Over ℤ it is the fraction-free division of [MP11]: `2x` divided by `3x`
    leaves `2x`.
- **Exact division** (`div_exact`, `div_exact_with`) stops at the first term
  of `r`. That is the test gcd algorithms make, and it costs nothing past the
  point where the answer is known. It also stops at the first quotient term
  past the degrees of `a` less those of `b` in some variable, which an exact
  quotient has (ours, CAS-137), so that a failed division never outgrows `a`'s
  exponents: `x^10` over `x − y^200` fails at once (gcd.md).
- **Fields** (ours, CAS-136). Unlike a product, a lexicographic division has
  no bound on its exponents known before. A reduction replaces a monomial
  `m` with `m − lm(b) + bⱼ`, which can raise the later variables at every
  step: `x^10` over `x − y^200` leaves `y^2000`, past every exponent of the
  inputs. So the division starts in the wider layout of `a` and `b`. When the
  guards report an overflow, it starts over in fields twice as wide, up to
  64 bits. The property `exponents_across_fields` found this: a first version
  sized the division by the dividend, which only holds in lexicographic order,
  not per variable.

## Other operations

- **Powers:** binary powering by heap products.
- **Evaluation:** at a point, each term's powers multiplied out.
- **Partial evaluation** (`eval_var`): rebuilds the terms, which must be
  sorted again.
- **Substitution** of a polynomial `q` for `x_var`: Horner's rule on the
  recursive view, the gaps between powers by powers of `q`.
- **The derivative** in `x_var`: takes one from that exponent in every term,
  which keeps the terms in order (the terms that differ first in `x_var` all
  lose one). Terms whose coefficient becomes zero (in characteristic `p`)
  drop out.
- **Degrees, coefficients and leading terms.**
  - The degree in a variable is a scan.
  - The coefficient of `x_var^d` keeps the terms with that power, the power
    cleared, in their order.
  - The leading term is the first.
- **The recursive view.** `p = Σ c_d x_var^d`, the coefficients without
  `x_var`, the powers decreasing. `from_recursive` puts it back.
- **Dense polynomials.** `to_dense` and `from_dense` convert with CAS-135's
  `Dense` in one variable.
- `DenseRing<SparseRing<R>>` is the dense recursive form, since `SparseRing`
  is a coefficient ring.

## Benchmarks (CAS-DOC-8 §13)

`sparse_workload`, an ignored test of `src/sparse/tests.rs`, measured with
`perf stat -e instructions:u`. The factors' own construction is subtracted.

| Product | Terms | Products of terms | Instructions | Wall time |
|---|---|---|---|---|
| Fateman [Fat03]: `f = (1 + x + y + z + t)^20`, `f (f + 1)` | 10626 × 10626 → 135751 | 113 M | 167.0 G | 9.8 s |
| Monagan–Pearce: `(1 + x + y + 2z² + 3t³ + 5u⁵)^12 (1 + u + t + 2z² + 3y³ + 5x⁵)^12` | 6188 × 6188 → 5821335 | 38.3 M | 50.3 G | 2.9 s |

That is about 1,400 instructions per product of terms. Most go to cas-num's
general integers, which allocate once a coefficient passes a word, and to
the heap. Accumulating coefficients in a few fixed words, and chaining the
heap's equal monomials [MP07], would cut both (CAS-141).

## Tests

`src/sparse/tests.rs`, against a naive reference: a map from exponent vectors
to coefficients.
- **Properties, over ℤ, ℚ and 𝔽ₚ** (word primes from 2 up):
  - sums, products both ways, scaling, powers;
  - division (the reference reduces the largest remaining term, as the heap
    does) and exact division, with `a = q b + r` checked;
  - evaluation, partial evaluation, substitution, the derivative;
  - degrees, coefficients, the leading term, the recursive view and its
    inverse.
- **Mixed layouts.** Exponents around 127 and 255, the largest of 8 and 9
  bits, so that operands come in different layouts and products cross into
  wider ones.
- **Boundaries.**
  - The guards report overflow in every field of a layout too narrow.
  - Products sized from 8 to 64 bits, repacked down when large exponents
    cancel.
  - Squares across a field; seventeen variables over several words.
  - The two division examples above.
- **The layout.** Round trips at every width, the order of packed vectors,
  the guards.
- **Sabotages** (CAS-136): fifteen changes are caught.
  - One was first missed, keeping the narrower layout in a sum; the mixed
    layouts now catch it.
  - Changing where a quotient's products start makes the division endless, and
    the test is stopped.
  - Widening a division straight to 64 bits gives the same results, more
    slowly.
