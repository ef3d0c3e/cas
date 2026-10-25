# Factoring

`src/factor.rs` and `src/factor/`, factoring over ℚ (CAS-142). Its steps are
planned in CAS-PLAN-1. This chapter grows with them.

## What is computed

A factorisation is a `Factors`: a constant `c` and polynomials `fᵢ` with
exponents `eᵢ ≥ 1`, where `p = c f₁^e₁ ⋯ fₖ^eₖ`.

- Each `fᵢ` has a positive degree, is primitive, and has its leading
  coefficient in the normal form of the coefficients: positive over ℤ, in the
  first quadrant over ℤ\[i\], and monic over a field. For sparse polynomials
  this is the lexicographic leading coefficient, the first variable most
  significant, as for gcds (gcd.md).
- Over ℤ, `c` is the content with the sign of the leading coefficient. Over a
  field, `c` is the leading coefficient. The zero polynomial has `c = 0` and
  no factors.

The domains are the `CharacteristicZero` ones: ℤ, ℚ, ℤ\[i\] and ℚ(i). In
characteristic `p`, a `p`-th power has derivative 0, so square-free
decomposition needs `p`-th roots as well. Factoring over ℚ does not need it.

## Square-free decomposition (CAS-179)

`squarefree` writes `p = c f₁ f₂² ⋯ fₖᵏ` with the `fᵢ` square-free and
pairwise coprime. Parts of degree 0 are left out, so the exponents increase.
giac's `sqrfree(12x⁴ − 24x² + 12)` is `12 (x² − 1)²`, which is the same
decomposition.

**Yun's algorithm [Yun76].** It runs on the primitive part `q = ∏ fᵢ^i`.

- `gcd(q, q′) = ∏ fᵢ^(i−1)`. Its cofactors are `b₁ = ∏ fᵢ` and
  `c₁ = q′ / gcd(q, q′)`.
- Round `k` has `bₖ = ∏_(i≥k) fᵢ` and
  `cₖ − bₖ′ = fₖ Σ_(i>k) (i − k) fᵢ′ ∏_(j>k, j≠i) fⱼ`. So `fₖ` is the gcd of
  `bₖ` and `cₖ − bₖ′`, and its cofactors are `bₖ₊₁` and `cₖ₊₁`.
- The rounds stop when `bₖ` is a unit.

One gcd with cofactors per round gives everything. Over ℤ the gcds of
primitive polynomials are primitive, so every division is exact. The constant
is `lead(p) / ∏ lead(fᵢ)^i`. That division is exact, and over ℤ\[i\] it also
recovers the unit: a product of first-quadrant leading coefficients need not
be in the first quadrant.

**Several variables.** The primitive part is split in the first variable `v`
that occurs, as `content_v · primitive_v`.

- The content in `v` has no `v`. It is decomposed the same way in the
  remaining variables, recursively.
- `primitive_v` has no factor without `v`. In characteristic 0, an
  irreducible factor `h` of `p` that has `v` divides `∂p/∂v` only if `h²`
  divides `p`. So Yun's algorithm with the derivative in `v` decomposes it,
  and a gcd constant in `v` is a unit.

Parts with the same exponent from different variables are then multiplied
together, so the result is the unique decomposition. For example,
`4x²y − 4y` is `4 · (x²y − y)`, where `y` and `x² − 1` share the exponent 1.

## Factors over prime fields (CAS-180)

`factor::finite` factors monic, square-free polynomials over a word-size
prime field `𝔽ₚ` of odd characteristic into their monic irreducible factors.
Factoring over ℚ needs only this: it reduces a square-free integer polynomial
modulo primes that keep it square-free, and Hensel lifting starts from the
factors there. The module is crate-internal. A public factorisation over
`𝔽ₚ` would also need square-free decomposition in characteristic `p`.

**Distinct degrees [vzGG13 §14.2].** For `d = 1, 2, …`:
- `h = x^(pᵈ)` modulo what is left of `f`, as the `p`-th power of the last
  `h` by square-and-multiply;
- `gcd(h − x, f)` is the product of the irreducible factors of degree
  dividing `d`. Smaller degrees are gone by then, so it is those of degree
  exactly `d`;
- the product is divided out, and `h` is reduced modulo what is left.

When what is left has degree below `2(d + 1)`, it is 1 or irreducible. Each
round costs `log₂ p` products modulo `f`, so small primes are cheap. This
step also gives the number of factors and their degrees without splitting
them.

**Equal degrees [CZ81; vzGG13 §14.3].** A product `g` of `r` irreducibles of
degree `d` splits by a random `a` of degree below `deg g`. If `a` is coprime
to `g`, then `b = a^((pᵈ − 1)/2)` is `±1` modulo each factor, since each
residue field has `pᵈ` elements. So `gcd(b − 1, g)` is a proper factor unless
all signs agree, which happens with probability about `2^(1−r)`. The pieces
split again, from a stack, until each has degree `d`. `gcd(a, g)` is tried
first, and catches an `a` that is not coprime.

The random coefficients come from `SplitMix64` [SLF14] with a fixed seed.
The factors are sorted by degree, then by coefficients, so the result does not
depend on the random choices anyway.

## Hensel lifting (CAS-181)

`factor::hensel::lift`, crate-internal. It starts from `f ≡ lc(f) g₁ ⋯ gᵣ
(mod p)`, where `p` does not divide `lc(f)` and the `gᵢ` are monic and
pairwise coprime modulo `p`. It returns the monic `Gᵢ ≡ gᵢ (mod p)` with
`f ≡ lc(f) G₁ ⋯ Gᵣ (mod pᵏ)`. Since the `gᵢ` are coprime, these lifts are
unique. Their coefficients are the representatives in `[0, pᵏ)`.

**The target.** Every factor is monic, so what lifts is the monic
`lc(f)⁻¹ f` modulo each power. The leading coefficient comes back when the
factors are recombined.

**The tree [vzGG13 §15.5].** The factors are the leaves of a binary tree,
split in halves. Each inner node holds:
- the product of its leaves;
- Bézout coefficients `s` and `t` of its two children `g` and `h`, with
  `s g + t h ≡ 1`, `deg s < deg h` and `deg t < deg g`. They start from the
  extended Euclidean algorithm modulo `p`.

**The Hensel step [vzGG13 Algorithm 15.10].** A node's product `f` must
become a given target. From `f ≡ g h` and `s g + t h ≡ 1` modulo `m`, it
computes the same modulo `m′`, for `m | m′ | m²`:
- `e = f − g h`; `s e = q h + r`; `g* = g + t e + q g`; `h* = h + r`;
- `b = s g* + t h* − 1`; `s b = c h* + d`; `s* = s − d`;
  `t* = t − t b − c g*`.

`h` is monic, so the divisions need no inverse. Then the children lift to
`g*` and `h*`, down to the leaves.

**The exponents.** A step at most squares the modulus. The exponents run up
`⌈k/2ʲ⌉, …, ⌈k/2⌉, k`, so the last step lands on `pᵏ` and does not overshoot
it, as FLINT does. Each step computes in `ℤ/pᵉℤ` (cas-num's `ModRing`, in
Montgomery form). Products there use Karatsuba's method.

## Univariate factoring over ℤ and ℚ (CAS-182)

`DenseRing<Integers>::factor` writes `p = c f₁^e₁ ⋯ fₖ^eₖ` with the `fᵢ`
distinct and irreducible over ℤ, primitive, with positive leading
coefficients. They are sorted by degree, then by coefficients from the
constant term up, then by exponent: giac's `factor(x⁶ − 1)` comes out as
`(x − 1)(x + 1)(x² − x + 1)(x² + x + 1)`. Over ℚ (`DenseRing<Rationals>`) the
factors are those of the numerator over ℤ, made monic, and `c` is the
leading coefficient.

The square-free decomposition (CAS-179) gives the exponents. Each part is
then factored by Zassenhaus's algorithm [Zas69; vzGG13 §15.6]:

1. **`x` first.** A factor `x` occurs at most once in a square-free part. It
   is split off, so that `f(0) ≠ 0` for the constant-term test below.
2. **The prime.** Odd primes from 3 up are tried. A prime is good when it does
   not divide `lc(f)` and keeps `f` square-free. For each of the first three
   good primes, the distinct-degree factorisation (CAS-180) gives the number
   of local factors and their degrees. The subset sums of those degrees are
   the degrees a factor over ℤ can have, and the sets of all three primes are
   intersected. One local factor, or no allowed degree strictly between 0 and
   `deg f`, proves `f` irreducible on the spot. Otherwise the prime with the
   fewest local factors is used, and its products by degree are split by
   Cantor and Zassenhaus's method.
3. **The bound.** Let `g` be a factor of `f` of degree below `n = deg f`. By
   Mignotte's bound [Mig74] the coefficients of `lc(f)/lc(g) · g` are at most
   `B = 2^(n−1) ‖f‖₂`. The local factors are lifted (CAS-181) to the least
   `pᵏ > 2B`. `lc(f)/lc(g) · g` is then the symmetric residue of `lc(f) ∏ Gᵢ`
   for the lifts `Gᵢ` of its local factors. After factors are divided out,
   what is left has a leading coefficient dividing `lc(f)`, so the bound
   still holds.
4. **Recombination.** Subsets of the lifts are tried by increasing size `s`.
   - Each subset's degree must be an allowed one.
   - The constant-term test: `lc(f) ∏ Gᵢ(0)` modulo `pᵏ` must divide
     `lc(f) f(0)`, since a true factor's constant term does. This test costs
     `s` products of integers.
   - A candidate that passes both is the primitive part of `lc(f) ∏ Gᵢ`,
     reduced symmetrically modulo `pᵏ`. It is kept when it divides `f`. Then
     `f` becomes the quotient, the subset's lifts are removed, and the search
     goes on at the same size.
   - At `s = r/2` a subset and its complement give the same split, so only
     the subsets with the first lift are tried.
   - Once no subset of at most half the lifts left gives a factor, none does,
     and what is left of `f` is irreducible.

The recombination is exponential in the number of local factors in the worst
case. Swinnerton-Dyer polynomials are that case. For `n` square roots their
Galois group is `(ℤ/2)ⁿ`, so modulo every prime they have at least `2ⁿ⁻¹`
local factors, of degree 1 or 2, and only one true factor. The one of degree
16, from √2, √3, √5 and √7, takes 127 subsets with 8 local factors. van
Hoeij's method is CAS-188. `factor_with` polls `Cancel` before each subset.

## Multivariate Hensel lifting (CAS-183)

`factor::lifting::lift_variables`, crate-internal, is the lifting half of
Wang's algorithm [Wan78], as Geddes, Czapor and Labahn give it [GCL92
§6.7–6.8].

**Input.** `f` over ℤ, a main variable `x`, and a point `a` of the other
variables `y₁, …, yₘ`, in the order in which they lift. At the point,
`f(x, a) = ∏ uᵢ(x)`, with the `uᵢ` pairwise coprime. The leading
coefficients `Lᵢ` of the factors in `x` are given, with `Lᵢ(a) = lc(uᵢ)`.
Wang's algorithm finds them (CAS-184).

**Output.** The `Uᵢ` with `f = ∏ Uᵢ`, `Uᵢ(x, a) = uᵢ` and `lc_x(Uᵢ) = Lᵢ`, or
nothing when no such factorisation exists.

**The modulus.** The arithmetic is modulo a prime `P`, chosen by the caller
above twice a bound on the coefficients of the `Uᵢ`. So the lifts' symmetric
residues are the factors over ℤ, and a final product checks them. `P` is
also above the degrees, so the factorials of Taylor coefficients are
invertible. Rather than lifting the coefficients `p`-adically too, as
[GCL92] does from a small prime, this lifts modulo one large prime, as
SymPy's `dmp_zz_wang` does [SymPy].

**Lifting `yⱼ` [GCL92 Algorithm 6.4].** The variables after `yⱼ` stay at their
values.
- Each factor takes `Lᵢ`, at those values, as its leading coefficient in `x`.
  The error `e = f − ∏ Uᵢ` then has a smaller degree in `x` than `f`.
- For `k = 1, …, deg_(yⱼ) f`, the coefficient `c` of `(yⱼ − aⱼ)ᵏ` in `e` is
  `e⁽ᵏ⁾(aⱼ)/k!`.
- The diophantine equation `Σ σᵢ ∏_(l≠i) U′ₗ = c` in the factors `U′` from
  before the step gives the corrections `Uᵢ += σᵢ (yⱼ − aⱼ)ᵏ`.

**The diophantine equation [GCL92 Algorithm 6.2].** It solves
`Σ σᵢ bᵢ = c`, with `bᵢ = ∏_(l≠i) aₗ` and `deg_x σᵢ < deg_x aᵢ`, modulo
`(y − α)^(d+1)` for each lifted variable, `d` being the largest degree of `f`
in them.
- The last variable is set to its value, and the equation is solved for the
  images, recursively.
- The error is then corrected one power of `y − α` at a time, as in the
  lifting.
- Univariate, at the point, the factors are the `uᵢ`. With `tᵢ` the inverse
  of `bᵢ` modulo `uᵢ`, computed once, `σᵢ = c tᵢ mod uᵢ`. Both sides of
  `Σ σᵢ bᵢ = c` agree modulo every `uᵢ`, so modulo their product. Both have
  a smaller degree, since `deg c < deg ∏ uᵢ`.

## Multivariate factoring over ℤ and ℚ (CAS-184)

`SparseRing<Integers>::factor` writes `p = c f₁^e₁ ⋯ fₖ^eₖ` as the dense
`factor` does. Leading coefficients are the lexicographic ones, so
positive and primitive refer to them. The factors are sorted by total
degree, then by terms from the leading one, then by exponent: giac's
`factor(x⁴ − y⁴)` is `(x − y)(x + y)(x² + y²)`. Over ℚ
(`SparseRing<Rationals>`) the factors are those of the numerator, made monic
in the lexicographic order, and `c` is the leading coefficient.

Each square-free part (CAS-179) is factored:
- **One variable.** Zassenhaus's algorithm (CAS-182).
- **Several.** The main variable `x` is one of least degree. The content in
  `x` is factored the same way in the other variables. The primitive part
  goes to Wang's algorithm [Wan78], arranged as in SymPy's `dmp_zz_wang`
  [SymPy]:

1. **The leading coefficient** `lc_x(f)` is factored first, recursively:
   `Ω ∏ F_j^(e_j)`.
2. **Points.** Values of the other variables are drawn from `[−r, r]`,
   starting with `r = 1`. `r` grows by 2 after 20 points, and no point is
   tried twice. A point is good when:
   - `lc_x(f)` does not vanish there;
   - the image `u = f(x, a)` is square-free;
   - Wang's condition holds: each value `E_j = F_j(a)` has a prime factor
     that divides no earlier `E_j`, nor `δ Ω`, `δ` the image's content.

   The image's primitive part is factored by Zassenhaus's algorithm. An
   irreducible image proves `f` irreducible: a factorisation of `f` would
   give one of the image, since the degree in `x` does not drop. Three good
   points are gathered, keeping only those whose images have the fewest
   factors. Of those, the image with the smallest coefficients is chosen.
3. **Leading coefficients.** Each factor `hᵢ` of the image takes the
   `F_j^k` whose value divides `lc(hᵢ) δ`, as often as it does. A constant
   then makes it `lc(hᵢ)` at the point. The content `δ` is spread over the
   factors. What remains of it multiplies every factor and its leading
   coefficient, and `f` takes its `(r − 1)`-th power. Wang's condition makes
   the attribution unique. An `F_j` that goes nowhere, or an inexact
   division, means extraneous factors: new points are tried.
4. **Lifting** (CAS-183) is modulo a prime above twice `2^(Σ deg) ‖f‖₂`,
   summing the degrees in every variable. That bound holds for any factor
   `h` of `f`: by Mahler's and Landau's inequalities its coefficients are
   within `2^(Σ deg h) M(h) ≤ 2^(Σ deg f) M(f) ≤ 2^(Σ deg f) ‖f‖₂`, where the
   Mahler measure `M` is at least 1 for a non-zero integer polynomial. The
   lifts' primitive parts are the factors. A lift that fails sends the
   search back to new points with a wider range.

Most points give images without extraneous factors (Hilbert's
irreducibility theorem), so the search ends with probability 1. It polls
`Cancel` at each point.

## Tests

`src/factor/tests/` has one module a step.

Square-free decomposition (`squarefree.rs`):

- **Checks of every result.** The decomposition multiplies back. Its exponents
  increase. Each factor has a positive degree, is primitive and normalised,
  and is square-free: coprime to its derivative in each variable it has, up
  to factors without that variable. The factors are pairwise coprime.
- **Random products.** `k f g² h³` over ℤ, checked and compared with Musser's
  algorithm, a reference independent of Yun's. Over ℚ the factors are the
  monic associates and the constant is the leading coefficient. Sparse cases
  have one to three variables, and in one variable match the dense result.
- **Known parts.** Exact decompositions of products of powers of known
  irreducible, pairwise coprime parts:
  - distinct linear factors and `x² + b` over ℤ;
  - linear factors with distinct Gaussian roots over ℤ\[i\] and ℚ(i);
  - distinct linear forms in one to four variables, merged by exponent. Forms
    without the first variable come from contents.
- **Examples.** giac's answers, signs and contents, and `(x + 1)^10`.

Over prime fields (`finite.rs`):

- **References.** The monic irreducibles of degrees 1 to 4 over `𝔽₃`, `𝔽₅`
  and `𝔽₇`, found by trial division by every monic polynomial of at most half
  their degree. Their counts match Gauss's formula.
- **Products of known irreducibles.** Random products of distinct
  irreducibles split back into them. The distinct-degree step gives their
  products by degree, and the equal-degree step splits each.
- **Large primes.** Random monic square-free polynomials over `𝔽ₚ` for
  `p = 65537`, `2⁶¹ − 1` and `2⁶² − 57`. The factors multiply back and pass
  Rabin's irreducibility test: `x^(pⁿ) ≡ x (mod f)` and
  `gcd(x^(p^(n/r)) − x, f) = 1` for each prime `r | n`.
- **Examples.** `x^p − x` splits into the `p` linear factors.
  `x⁹ − x` over `𝔽₃` is the product of the irreducibles of degrees 1 and 2.
  Also the generator's known outputs, and cancellation.

Hensel lifting (`hensel.rs`):

- **Random products.** Products of random polynomials, square-free modulo
  3, 11, 101 or 65537 with a leading coefficient the prime does not divide.
  Their local factors lift to exponents 1 to 11. Each lift is monic, has
  coefficients below `pᵏ`, and reduces to its factor modulo `p`. `lc(f)` times
  the product of the lifts is `f` modulo `pᵏ`. Each lift to `pᵏ` reduces to
  the lift to `pʲ` for every `j < k`, which checks uniqueness. Two groups of
  the local factors, multiplied together, lift as well.
- **An example.** `x² + 1 ≡ (x − 2)(x + 2) (mod 5)` lifts to `(x − a)(x + a)`
  with `a² ≡ −1` and `a ≡ 2`: a 5-adic square root of −1, up to `5²⁰`.

Univariate factoring (`zassenhaus.rs`):

- **Checks of every result.** The factorisation multiplies back. Its factors
  are sorted, distinct, primitive and positive, and each factors as itself.
- **Random products.** Products of powers of up to three random polynomials.
  Their factorisation is the union of the parts': the factors with their
  exponents summed, and the product of the constants. Over ℚ the factors are
  the monic associates, and the constant is the leading coefficient.
- **Known irreducibles.** Products of distinct Eisenstein polynomials factor
  into exactly them. `xⁿ − 1` for `n ≤ 30` is the product of the cyclotomic
  `Φ_d`, computed by division, for `d | n`. The Swinnerton-Dyer polynomials of
  degrees 4, 8 and 16, `x⁴ + 1` and `x² − 2` are irreducible.
- **Examples.** giac's answers, Sophie Germain's `x⁴ + 4`, coefficients of
  hundreds of bits, twelve linear factors, and cancellation.

The fuzz target `poly_factor` checks the same properties on products of
random polynomials, with coefficients of up to three limbs.

Multivariate Hensel lifting (`lifting.rs`):

- **Random factors.** Two or three random factors in two or three variables,
  with leading coefficients in the others. Where their images at a random
  point are coprime and keep their degrees, the images and leading
  coefficients lift back to exactly the factors.
- **Examples.** `x² − y²` from `y = 2`. A leading coefficient `y`. Three
  factors in three variables, lifted in both orders of the variables. A main
  variable other than the first.
- **Wrong leading coefficients.** Leading coefficients that agree with the
  images but that no factorisation has: nothing comes back.

Multivariate factoring (`wang.rs`):

- **Random products.** Products of powers of up to three random
  polynomials in one to three variables. The checks: the result multiplies
  back, and its factors are sorted, distinct, primitive and positive, and
  factor as themselves. The factorisation is the union of the parts'. In one
  variable it matches the dense factorisation. Over ℚ the factors are the
  monic associates.
- **Examples.** giac's `x⁴ − y⁴` and `x²y + xy²`. `x² − 2y²` and
  `x² + y² + 1` are irreducible. Three factors in three variables, leading
  coefficients `y` and `y²`, a content `y + 1` in the main variable, powers,
  a negative constant, and cancellation.

The fuzz target `poly_factor` also checks products in two and three
variables, of total degree at most 12.
