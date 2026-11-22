# Gröbner bases

`src/groebner.rs` (CAS-287, planned in CAS-PLAN-15; the design is
CAS-DOC-21) and `src/groebner/zero_dim.rs` (CAS-541).

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

The caller's `GroebnerLimits` bound its time (CAS-291):

- **Work.** Each coefficient made or copied counts 1, plus 1 for each 64 bits
  of its numerator and denominator. Each critical pair formed or gone
  through counts 1, and so does each term moved to a remainder.
- **Coefficient size.** A coefficient whose numerator or denominator has
  more than `coefficient_bits` bits stops the computation.

So time follows the work, whatever the sizes. A count of division steps did
not bound it. One step merges the whole remainder with a multiple of a
divisor, which can be thousands of terms in a lexicographic basis, and the
rational coefficients grow without bound: the fuzz target ran for minutes on
one input under a budget of 200 000 steps.

The remainder keeps what is left as the tail of its term vector. The terms
that no leading monomial divides move out without the vector shifting.

## Normal sets and the rational univariate representation

`src/groebner/zero_dim.rs` (CAS-541).

- `normal_set`: the monomials of a reduced basis's quotient ring, those no
  leading monomial divides, in increasing order; `None` unless every
  variable has a pure power among the leading monomials (the finiteness
  theorem, [CLO15, §5.3]), empty for `[1]`. The monomials are walked
  variable by variable, the exponent of each growing until a leading
  monomial divides, which its pure power ends: every multiple of a
  divisible monomial is divisible, so the walk visits few more monomials
  than it keeps.
- `rational_univariate`: Rouillier's representation [Rou99] of the points
  of a zero-dimensional ideal, a separating linear form `l` and polynomials
  `f, f', gᵢ` in `t` with the points `xᵢ = gᵢ(l)/f'(l)` at the roots of
  `f`, each once. giac's `gbasis(…, rur)`.

The base algorithm works in the quotient ℚ\[X\]/I as a vector space over
the normal set: the minimal polynomial of a linear form `l` is the first
linear dependence of the normal forms of `1, l, l², …`, each the normal
form of `l` times the last, found by incremental Gaussian elimination over
ℚ, each row carrying the polynomial in `l` it is the normal form of. The
same rows give the `hᵢ` with `NF(xᵢ) = hᵢ(l)`.

Our choices, fixed by giac 2.0.0's answers (CAS-541):

1. `G` the reduced `DegRevLex` basis, `D` the size of its normal set. The
   minimal polynomial `mᵢ` of each variable, from the last to the first;
   the first of degree `D` is `l` (for `x², y`, `l = x`: `y` has degree
   1).
2. The radical `I' = I + (sqfree(mᵢ)(xᵢ))` for every `i`, `sqfree(m) =
   m / gcd(m, m')`: in characteristic 0, Seidenberg's lemma [Sei74]. Its
   normal set has `D'` monomials, one for each point.
3. The form: the variable of step 1 if any, else the single variables from
   the last, else `x₁ + k·x₂ + k²·x₃ + … + k^(n−1)·xₙ` for `k = 1, 2, …,
   n·D'²`, the first whose minimal polynomial modulo `I'` has degree `D'`.
   Two points give the same value of the form for at most `n − 1` values
   of `k`, so one of the `n·D'²` separates. giac's own combination differs
   (`−6x − 4y` for `x² − 1, y² − 1`); the representation is then another
   one of the same points.
4. `f` is minus the monic minimal polynomial of `l` modulo `I'`, as giac
   gives it, `f'` its derivative, `gᵢ = hᵢ·f' rem f`.

`None` for the whole ring, an ideal that is not zero-dimensional, and one
variable, where giac gives the basis.

Bounds: a normal set past `MAX_NORMAL_SET` (256) monomials stops the
computation with `TooLarge`, so the vectors and the elimination stay small;
the bases, normal forms and eliminations share one `GroebnerLimits` budget
and poll the `Cancel`.

## Syzygies

`syzygies` gives giac's `syzygy`: for `f₁, …, fₘ` in `n` variables, the
reduced `DegRevLex` basis, in the variables `A₁, …, Aₘ, x₁, …, xₙ`, of the
ideal of the `Σ aᵢAᵢ` with `Σ aᵢfᵢ = 0` and of every `AᵢAⱼ`. Its elements
linear in the `Aᵢ` generate the module of syzygies.

The base algorithm is elimination [KR00]: the module is the kernel of
`(a₁, …, aₘ) ↦ Σ aᵢfᵢ`, and a tag variable `E` turns it into an ideal.
With the generators `fᵢE + Aᵢ` and every monomial of degree 2 in `E, A₁,
…, Aₘ`, an element of the ideal that is free of `E` is, modulo the
products of two of `E, A`, a `Σ aᵢAᵢ` whose coefficient of `E`, `Σ aᵢfᵢ`,
vanished. Our choice: their `Lex` basis with `E` the greatest, then the
`A`, then the `x`, keeps the elements whose leading monomial is free of
`E`; the `E` exponent dropped, they are reduced again under `DegRevLex`, so
the answer is unique. Bounds: those of the two basis computations, which
share the `GroebnerLimits`.

## Ideal membership

giac's `greduce` of a polynomial by a list reduces by the list as given,
which decides membership only when the list is a Gröbner basis. Membership
is decided by reducing by the reduced basis of the ideal (`groebner_basis`,
then `groebner_reduce`): the remainder is 0 exactly for the members [CLO15,
§2.6, Corollary 2].

## Variable order for revlex

When the variables are not given, `gbasis_vars` in cas-core (CAS-541)
orders them as giac does before a `revlex` basis: a stable sort by the
product over the generators of `deg_v + 1`, `deg_v` the degree of the
generator in `v`, the ties in their order of first appearance. The order
changes the basis, not the ideal.

## Tests

`src/groebner/tests.rs`: the orders on [CLO15, §2.2]'s examples; division
by both orders of the divisors and an exact rational remainder, giac's
`greduce`; the bases of `x² + y² − 1, x − y` under each order, of `x³ − 2xy,
x²y − 2y² + x` [CLO15, §2.7] under `DegRevLex` and `Lex`, of the twisted
cubic, of cyclic-3, of `x²y − 1, xy² − x`, giac's answers once primitive;
the whole ring, the zero ideal, a rational coefficient, an exponent past the
bound; computations past each of the limits, a basis and a division past
little work, and a basis whose coefficient needs one bit more than allowed. Each basis is checked for the defining properties: every input
reduces to 0 by it, so does every S-polynomial of two of its elements, each
element is monic and no term of one is divisible by another's leading
monomial. The fuzz target `poly_groebner` checks the same properties on
small random ideals in one to three variables under each order, within
limits on the work and the coefficients.

`src/groebner/zero_dim/tests.rs`: normal sets of a reduced basis, of a
monomial ideal, of `[1]`, `None` for a curve and for the zero ideal; the
rational univariate representations giac gives for `x² + y² − 1, x − y`;
`x² − 2, y² − 3`; `x², y`; `x²(x − 1), y`; `x² − 1, (y − 1)²`; `x³ − y,
y² − x`; `x − 1, y − 2, z − 3`; `x² − 3x + 2, y − x²`; `None` for a curve,
a line, the whole ring, one variable and no polynomial; for `x² − 1, y² −
1`, where giac's form differs, the degree of `f` and each generator
vanishing at the representation modulo `f`; a normal set past the bound.
The syzygies of `x, y`; `xy, x² − y³`; `x − 1, x² − 1`; `x², xy, y²`,
giac's answers up to constant factors, each element homogeneous in the
`Aᵢ` and those linear in them syzygies.
