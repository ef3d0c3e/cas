# Gröbner bases

`src/groebner.rs` (CAS-287, planned in CAS-PLAN-15; the design is
CAS-DOC-21).

## What it computes

- `groebner_basis`: the reduced Gröbner basis of the ideal a list of
  polynomials generates in ℚ\[x₁, …, xₙ\], under a monomial order: each
  element monic, in increasing order of the leading monomials; `[1]` for the
  whole ring, `[]` for the zero ideal. The reduced basis of an ideal under an
  order is unique [CLO15, §2.7, Theorem 5], so it is the answer to compare.
- `groebner_reduce`: the remainder of a polynomial divided by a list, giac's
  `greduce`; the normal form of the polynomial modulo the ideal when the list
  is a Gröbner basis.
- `s_polynomial`, for the tests and the fuzz target's checks.

## Representation

A `Distributed` polynomial is its terms in decreasing order of their
monomials under one order, each an exponent vector (`Vec<u64>`) and a
non-zero `Rational`, like monomials added on construction. The order is not
stored: every operation takes it, and it must be the one the terms were
sorted by. Multiplying by a monomial keeps the order, as every monomial order
is compatible with products [CLO15, §2.2, Definition 1].

Gröbner bases work on few terms of high multiplicity of operations, where
the packed lexicographic `Sparse` of `sparse.md` would have to be resorted
for each order: plain vectors serve, and differences are merges of the two
sorted term lists.

## Orders [CLO15, §2.2]

With `x₁ > x₂ > … > xₙ`:

- `Lex` (giac's `plex`): the exponent vectors compared from `x₁`.
- `DegLex` (`tdeg`): the total degrees, then `Lex`.
- `DegRevLex` (`revlex`, giac's default): the total degrees, then the last
  variable whose exponents differ, the smaller exponent the greater
  monomial.

## Division [CLO15, §2.3, Theorem 3]

The leading term of what is left is divided by the first divisor whose
leading monomial divides it, and the multiple subtracted; a leading term that
none divides moves to the remainder. The remainder is then reduced: none of
its terms is divisible by a leading monomial. The quotients are not kept.
The order of the divisors matters ([CLO15, §2.3, Example 4]: `x²y + xy² +
y²` by `(xy − 1, y² − 1)` leaves `x + y + 1`, by the other order `2x + 1`),
but for a Gröbner basis, where the remainder is the normal form.

## Buchberger's algorithm [Buc65]

The basis grows by the S-polynomials of its pairs, each reduced by the basis
and added when its remainder is not 0, until no pair is left [CLO15, §2.7,
Theorem 2]:

- **The normal strategy** [Buc85]: the pair of the least lcm of the leading
  monomials first.
- **Criteria** [GM88], in Becker and Weispfenning's `UPDATE` [BW93, §5.5]:
  adding `h`, of its new pairs `{g, h}` those whose lcm another new pair's
  divides go, one of equal lcms staying, and then those whose leading
  monomials are coprime (Buchberger's first criterion), which served to drop
  the others; of the old pairs `{g₁, g₂}`, those whose lcm `lm(h)` divides
  go, unless it equals the lcm of `g₁` or `g₂` with `h`; and the basis
  elements whose leading monomial `lm(h)` divides leave the basis, their
  pairs staying.
- **Each element monic**, the remainders over ℚ kept small that way.
- **The whole ring**: a constant remainder ends the computation with `[1]`.

The inputs enter in increasing order of their leading monomials, each
reduced by the basis so far, so that the basis is minimal throughout: no
element's leading monomial divides another's. At the end each element is
reduced by the others, which leaves its leading monomial, so the basis is
reduced whatever the order of these reductions [CLO15, §2.7, Proposition 6].

## Bounds

A computation polls its `Cancel` at each division step and stops at
`MAX_ELEMENTS` polynomials added, `MAX_REDUCTIONS` S-polynomials reduced or
an exponent past `MAX_EXPONENT` (products of monomials saturate rather than
overflow, so such an exponent is noticed), with `GroebnerError`: an input
whose basis runs away fails rather than exhausting the memory.

## Tests

`src/groebner/tests.rs`: the orders on [CLO15, §2.2]'s examples; division
by both orders of the divisors and an exact rational remainder, giac's
`greduce`; the bases of `x² + y² − 1, x − y` under each order, of `x³ − 2xy,
x²y − 2y² + x` [CLO15, §2.7] under `DegRevLex` and `Lex`, of the twisted
cubic, of cyclic-3, of `x²y − 1, xy² − x`, giac's answers once primitive;
the whole ring, the zero ideal, a rational coefficient, an exponent past the
bound. Each basis is checked for the defining properties: every input
reduces to 0 by it, so does every S-polynomial of two of its elements, each
element is monic and no term of one is divisible by another's leading
monomial. The fuzz target `poly_groebner` checks the same properties on
small random ideals in one to three variables under each order, within a
budget of division steps.
