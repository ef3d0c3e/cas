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

## Products (CAS-136, CAS-155)

`src/sparse/product.rs`. A product goes through one of two kernels, and its
sums of products through one of three accumulators.

### Arrays, for products dense in their degree box

- **Chunks.** The box of the product's degrees in the variables after the
  first, `x₁ … xₙ₋₁`, is a chunk: an array of sums, one for each monomial,
  at the index `Σ eᵥ sᵥ`. The strides `sᵥ` are those of the box, `x₁`'s the
  largest, so that indices order monomials as lexicographic order does, and
  the index of a product is the sum of its factors' indices.
- **One degree of `x₀` at a time.** Lexicographic order keeps the terms of
  each degree in `x₀` together: the parts of `a` and `b`. For each degree
  `d` of the product in `x₀`, from the largest down, the parts whose degrees
  add up to `d` add their products into the chunk, with no comparison. The
  chunk is then read back from the largest index touched to the smallest,
  which gives the terms of `x₀^d` in order, and is left zero for the next
  degree. FLINT's `fmpz_mpoly_mul_array` likewise accumulates "one slice of
  the main variable" at a time [FLINT]. One chunk is all the memory, and it
  stays in the cache.
- **When.** A chunk holds at most `2^18` sums. The arrays then take a product
  when the sums read back and the parts of `a` looked up, `(deg_x₀ a +
  deg_x₀ b + 1) (chunk + parts)`, are at most 4 for each product of terms.
  - Fateman's product: 41 chunks of 41³ sums (1.6 MB of words) for 113 M
    products of terms.
  - Monagan and Pearce's: a chunk of 49³ · 73 sums, past the limit, and too
    sparse anyway (73 chunks for 38 M products).
  - In one variable, a chunk is one sum: dense products take the arrays,
    and `1 + x^1000` squared takes the heap.

### Johnson's heap [Joh74, MP07], otherwise

The products `aᵢ bⱼ` of a row `i` decrease with `j`. A binary heap merges the
rows:
- it holds one product per row, keyed by its monomial;
- the largest comes out, and the products of equal monomials add up;
- each popped `(i, j)` is replaced by `(i, j+1)`.

The result comes out in order, and needs no sort.
- **Lazy rows** [MP07]. Row `i + 1` enters the heap when `(i, 0)` leaves it.
  The heap then holds at most one product per row of the shorter factor,
  which is taken as the first.
- **Slots.** Each product's monomial is computed when it enters the heap, and
  stored in a slot that its successor reuses. The storage is linear in the
  rows, and the heap itself moves slot numbers, not monomials.
- **Chaining** [MP11 §5.2]. A product entering the heap is compared with the
  top first, then from the bottom up. When a comparison finds its monomial,
  it joins that slot's chain instead of taking a place. Every product of the
  top monomial comes out, chains and all, before their successors go in, so
  that the largest of those chain at the top.
- **Extraction** [MP11 §5.3]. The hole left at the top goes down along the
  larger children, one comparison a level. The last slot fills it, rising if
  it must. That takes fewer comparisons than sifting the last slot down from
  the top, two a level. MP11 does this for monomials kept outside the heap,
  as ours are. The division shares the heap, unchained.
- **Fields.** Both factors' exponents are repacked first (not their
  coefficients) in fields wide enough for the sum of their largest
  exponents, so no product overflows. The guards still check every monomial
  added. Were one to overflow, the product would start over one bit wider;
  `packed_fields_at_their_boundaries` drives that path with fields made too
  narrow on purpose.

`Cancel` is polled every 2^16 products of terms by the arrays, and every 1024
monomials by the heap.

### Accumulators [MP11 §5.3]

cas-num's `Ring::small_preimage(c)` is an `i64` whose image in the ring is
`c`, when the ring has one at hand: an integer that fits, an integral
rational, a word residue below `2^63`. The map from ℤ is a homomorphism, so a
sum of products of preimages, mapped into the ring once, is the sum in the
ring.
- **Words**, when every coefficient of both factors has a preimage. Each
  product of two `i64`s, at most `2^126` in absolute value, adds into three
  words, `hi 2^128 + lo`; fewer than `2^64` of them stay below `2^190`. Each
  term of the result maps its sum into the ring once, and no coefficient
  arithmetic allocates.
- **Mixed**, when some have: the products of two preimages in words, the
  others in the ring, added together at the end of each term. MP11 does the
  same with word-sized and multiprecision integers.
- **In the ring**, when none has (ℤ\[i\], ℚ(i), the big prime fields of
  Hensel lifting).

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

- **Powers:** binary powering by products.
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
giac 1.9's `time(expand(f*g))`, on the same factors expanded first, gives a
coarse comparison.

| Product | Terms | Products of terms | Instructions | Wall time | giac |
|---|---|---|---|---|---|
| Fateman [Fat03]: `f = (1 + x + y + z + t)^20`, `f (f + 1)` | 10626 × 10626 → 135751 | 113 M | 148.9 G | 9.8 s | 0.43 s |
| Monagan–Pearce: `(1 + x + y + 2z² + 3t³ + 5u⁵)^12 (1 + u + t + 2z² + 3y³ + 5x⁵)^12` | 6188 × 6188 → 5821335 | 38.3 M | 45.3 G | 2.9 s | 0.26 s |

That was about 1,300 instructions per product of terms, 23 and 11 times
giac's time, with the plain heap of CAS-136:
- **The heap** took three quarters of it (`Heap::greater` and `pop`, by
  perf), bound by its unpredictable comparisons. Marking its methods
  `#[inline]` (CAS-141), so that the products cas-core instantiates inline
  them too, cut 11% of the instructions here and 26% through cas-core's
  `expand`, but none of the time.
- **cas-num's integers** took a fifth: they allocate once a coefficient
  passes a word.

**CAS-155**, measured by the `Workload` workflow on a GitHub runner, the
best of three runs, against the commit before it (bf77326) on the same
runner; the runner has no `perf`, and the inputs' construction takes at
most 0.07 s:

| Product | Before | CAS-155 | Kernel, accumulator |
|---|---|---|---|
| Fateman | 12.62 s | 0.13 s | arrays (41 chunks of 41³ sums), words |
| Monagan–Pearce | 3.90 s | 1.56 s | chained heap, words |

- Fateman's product is 97 times faster: about one nanosecond for each
  product of terms, in a chunk that stays in the cache.
- Monagan and Pearce's is 2.5 times faster. Its 38 M products of terms go
  through the heap, and its 5.8 M terms have coefficients past a word.
- **Against giac**, coarsely: giac's times were taken on the development
  machine, where the old heap took 9.8 s and 2.9 s, three quarters of the
  runner's times. By that ratio, Fateman's product now takes about a
  quarter of giac's time, and Monagan and Pearce's about four and a half
  times giac's.
- **Next** for the heap: one-word monomials kept in the heap itself, MP11's
  "immediate monomials", so that its comparisons read no other memory.

The criterion benchmarks of `benches/poly.rs` run smaller sizes; at CAS-141,
with the plain heap, `fateman` took 0.76 ms and 58 ms at `K = 5` and `10`,
and `pearce` 0.90 ms and 109 ms at `K = 4` and `8`.

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
- **Kernels and accumulators** (CAS-155). Every kernel and accumulator that
  can take a product, forced, against the naive product:
  - over ℤ, coefficients small, at the edges of `i64` and past them, so that
    sums cross `2^127` and `2^128` and some products go to the ring;
  - over ℚ, integers and fractions together; over word prime fields,
    residues below `2^63` and above (modulo `2^64 − 59`);
  - `(1 + x + y + z)^7` squared, against the multinomial coefficients.
  - The three-word sums against the integers, at the edges of `i64` and past
    `±2^128`.
  - The kernel each product takes: arrays for Fateman's, the heap for
    Monagan and Pearce's, both in one variable; a chunk past its limit.
  - Both kernels stop when cancelled.
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
