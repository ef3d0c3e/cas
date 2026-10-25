# Greatest common divisors, resultants and discriminants

`src/gcd.rs` and `src/gcd/` (CAS-137).

## What is computed

The coefficient domains are the `GcdDomain`s: ℤ, ℚ, ℤ\[i\], ℚ(i), and the
prime fields of word and of any size. For dense polynomials (`DenseRing`) and
sparse ones (`SparseRing`):

- **`gcd`, `gcd_cofactors`, `lcm`.** The gcd is normalised: its leading
  coefficient is in the normal form of the coefficients. For sparse
  polynomials that is the lexicographic leading coefficient, the first
  variable most significant.
  - Over ℤ and ℤ\[i\] the gcd of the contents is a factor. Over ℤ the leading
    coefficient is positive: `gcd(6x² − 6, 4x − 4) = 2x − 2`, as CAS-DOC-8
    §15 detail 5 asks. Over ℤ\[i\] it is in the first quadrant.
  - Over a field the gcd is monic. Detail 5's form over ℚ, giac's `10x − 1`
    for `gcd(x² − 1/100, x − 1/10)`, is the gcd over ℤ of the numerators,
    which the builtin (CAS-140) takes.
  - `gcd(0, b)` is `b` normalised, where giac returns `b` as it is.
  - The cofactors are `a / g` and `b / g`.
  - The lcm is `a (b / g)`, normalised as the gcd.
- **Contents.** `content` and `primitive_part`, and for sparse polynomials
  `content_in` and `primitive_part_in`, by a variable.
- **`resultant`**: the determinant of the Sylvester matrix, in a variable for
  sparse polynomials. `res(c, b) = c^(deg b)` for a constant `c`, where giac
  gives 1, and 0 when either is 0.
- **`discriminant`**: `(−1)^(n(n−1)/2) res(p, p′) / lc(p)` in degree
  `n ≥ 1`, which is 1 in degree 1. `None` for a constant or zero.
- `RationalPolys` and `GaussianRationalPolys`, ℚ\[x\] and ℚ(i)\[x\] over one
  denominator, take their gcd from the numerators. It replaces Euclid's
  algorithm over ℚ, which `EuclideanDomain` provided. They also have
  `resultant` and `discriminant`.

`gcd` chooses the algorithm. `gcd_by` runs one `GcdMethod`, for tests and
benchmarks; it returns `None` where the method does not apply or gives up.

| Domain | One variable | Several variables |
|---|---|---|
| ℤ | heuristic, then Brown's modular algorithm | the same, with images by Brown's evaluation |
| ℤ\[i\] | heuristic, then subresultants | heuristic, then subresultants |
| ℚ, ℚ(i) | the numerators over ℤ, ℤ\[i\] | the same |
| 𝔽ₚ | Euclid | Brown's evaluation, then subresultants if the field runs out of points |

`GcdMethod::Zippel` runs Zippel's sparse interpolation over prime fields of
at least `2^31` elements (CAS-191); over ℤ it comes with CAS-192.

Resultants:

| Domain | Method |
|---|---|
| ℤ | modulo primes in one variable; subresultants in several |
| 𝔽ₚ | Euclid in one variable; subresultants in several |
| ℤ\[i\] | subresultants |
| ℚ, ℚ(i) | the numerators |

## Subresultants: the reference [BT71, Bro78, Coh93]

**The remainder sequence.** For `deg a ≥ deg b`, each round:
- pseudo-divides, `lc(b)^(δ+1) a = q b + r` with `δ = deg a − deg b` (Knuth's
  Algorithm R [Knu97]);
- moves on to `(b, r / (g h^δ))`;
- sets `g = lc(b)` and `h = g^δ / h^(δ−1)`.

The divisions by `g h^δ` and `h^(δ−1)` are exact [BT71], so the coefficients
stay small without a gcd of them at every step.

**The gcd** [Coh93, §3.3.1] is the gcd of the contents times the primitive
part of the last non-zero remainder of the primitive parts.

**Several variables.**
- The sequence runs in the first variable that occurs, as a
  `DenseRing<SparseRing<R>>` whose coefficients are the polynomials in the
  others. The contents and the primitive parts of the remainders come from
  the same algorithm in fewer variables.
- Pseudo-division needs no division of coefficients, so CAS-137 relaxed
  `pseudo_div_rem` from Euclidean domains, which polynomials in further
  variables are not, to any coefficient ring.

**The resultant** [Coh93, §3.3.2] runs the same sequence, without the
contents:
- each exchange of two odd degrees changes the sign, as
  `res(a, b) = (−1)^(deg a deg b) res(b, a)`;
- when `b` reaches a non-zero constant, the resultant is
  `lc(b)^(deg a) / h^(deg a − 1)`.

Over any domain, including polynomials in further variables, this is the
reference that the faster algorithms are checked against.

## The heuristic gcd [CGG89, Par05]

1. The contents are taken out. They matter only to the gcd of the contents,
   a factor of the result.
2. The polynomials are evaluated at an integer `ξ`, and `γ`, the gcd of the
   values, is computed.
3. `γ` is read back as the polynomial of its digits in base `ξ`, each in
   `(−ξ/2, ξ/2]`.
4. If the primitive part `G` of that polynomial divides both, it is the gcd.

**Parisse's theorem** [Par05, Theorem 1] makes the check enough, in any number
of variables, once `|ξ| ≥ 2 min(|a|, |b|) + 2`, `|a|` the largest coefficient.
- In several variables the last variable that occurs is evaluated, and the
  gcd of the values comes from the same algorithm in one variable fewer, as
  Parisse writes it.
- Char, Geddes and Gonnet's own proof has a gap in several variables, which
  Parisse's note fills [Par05].
- The theorem holds over ℤ\[i\] [Par05]. There the proof needs `|α| ≤ |ξ|/2`
  for the content `α` of the digits, but a Gaussian digit's modulus is only
  below `|ξ|/√2`. So the first point over ℤ\[i\] is `4 (h + 1)`, `h` the
  smaller of the largest `|re| + |im|`. A root of a factor is then below
  `|ξ|/4`, and `|C(ξ)| > 3|ξ|/4 > |ξ|/√2` closes the proof (ours, CAS-137).

**Retries.** A failed check moves to
`ξ ← ⌊73794 ξ ⌊ξ^(1/4)⌋ / 27011⌋`, at most 6 points in all, as SymPy's
`dup_zz_heu_gcd` [SymPy]. `x + 2` and `x − 2` fail at the first point, 6:
their values 8 and 4 have the gcd 4, read back as `x − 2`. The next point,
16, gives 1.

**The size limit** (ours). The heuristic gives up when `ξ^(d+1)` would pass
`2^18` bits, `d` the degree evaluated, and the modular algorithm takes over.
In several variables the values grow by a power at each variable: for
`(1 + x + y + z)^8 + x` the third ones need between `2^16` and `2^18` bits.
The limit was chosen from the benchmarks below:
- at `2^16`, the heuristic gave that case up, and Brown's algorithm took
  three times as long;
- at `2^20`, it took twice as long as Brown's algorithm on degree 400 with
  256-bit coefficients.

## Brown's modular algorithm [Bro71]

**Over ℤ.**
1. The contents are taken out. `γ` is the gcd of the leading coefficients
   of the primitive parts.
2. For each word prime `p` from `2^62` down that does not divide `γ`:
   - the images of the primitive parts have a monic gcd modulo `p`, which is
     scaled to the leading coefficient `γ mod p`;
   - images whose gcd has the same leading monomial are combined by Chinese
     remaindering, in the symmetric range.
3. Once the combination stops changing, its primitive part is tested: if it
   divides both, it is the gcd.

**Unlucky primes.** The image of the gcd `D` divides the gcd of the images,
and keeps its leading monomial, since `lc(D)` divides `γ`. So:
- a larger leading monomial marks an unlucky prime, which is dropped;
- a smaller one shows that the combination so far came from unlucky primes,
  and it starts over.

For example, the first prime `p₀` is unlucky for `(x + p₀)(x − 1)` and
`x (x − 1)`, whose images are equal.

**Why the test is enough.** A primitive divisor of both is a divisor of `D`.
Its leading monomial is that of every image, which is at least `D`'s. So it is
`D`, up to sign.

**Over a prime field**, in several variables, Brown's dense evaluation
`pgcd`:
1. The polynomials are taken in `F[y][X]`, `y` the last variable that occurs.
   Their contents and `γ`, the gcd of their leading coefficients, are
   polynomials in `y`.
2. Points `α`, from `2^31` up (ours: past the small integers, the likeliest
   bad points), where `γ(α) ≠ 0`, give images in one variable fewer, from the
   same algorithm.
3. Each image is scaled to `γ(α)` and compared by leading monomial, as the
   primes are.
4. Images with equal leading monomials are combined by Newton's
   interpolation.
5. Once the interpolation stops changing, or passes
   `deg γ + min(deg_y a, deg_y b)` points, the degree in `y` of the scaled
   gcd, its primitive part in `X` is tested by division.

In one variable, the gcd is Euclid's. A field with too few points (`𝔽₂` when
`γ = y(y + 1)`) runs out, and `gcd` takes subresultants.

**Bad points.**
- `x + y − 1` and `x` are coprime, but at `y = 1` both are `x`.
- `(x + y − 2)(x + y)` and `x (x + y)` have the gcd `x + y`, but at `y = 2`
  both are `x (x + 2)`. That point is dropped, and `y = 1, 3` interpolate
  `x + y`.

## Zippel's sparse interpolation [Zip79, Zip90, DKMW05]

`src/gcd/zippel.rs` (CAS-152), in de Kleine, Monagan and Wittkopf's form,
LINZIP. Brown's algorithm needs `deg_y + 1` images for each variable `y`,
whatever the gcd's terms. Zippel's needs, after a first image, as many as
the gcd has terms in a coefficient of the main variable.

**The form** of a gcd is its monomials, grouped by the power of the main
variable `x`, the first that occurs: the coefficient of each power, a
polynomial in the other variables, is a block.

**An image from the form** (CAS-190), over `𝔽ₚ`, all variables but `x`
evaluated:
1. At the powers `α, α², …` of a random point `α` of the other variables,
   the monic univariate gcds `gᵢ` of `a` and `b` in `x`. Each term of `a`
   and `b` goes from one power to the next by one product, its monomial's
   value at `α`.
2. The gcd at `αⁱ` is `λᵢ gᵢ`, `λᵢ` its leading coefficient in `x` there;
   for each block, the coefficient of its power in `gᵢ` gives `Σⱼ cⱼ vⱼⁱ =
   λᵢ yᵢ`, `vⱼ` the values of its monomials at `α`.
   - **Single scaling**: when a block is one monomial `μ`, the gcd scaled
     so that `μ`'s coefficient is 1 has `λᵢ = μ(α)ⁱ / yᵢ` for that block,
     and every block is a transposed Vandermonde system [Zip90], solved in
     `O(t²)` (`src/gcd/vandermonde.rs`, CAS-189): with `P = ∏ (Z − vₗ)`
     and `Pⱼ = P / (Z − vⱼ)`, `cⱼ vⱼ Pⱼ(vⱼ) = Σᵢ yᵢ [Zⁱ⁻¹] Pⱼ`.
   - **Multiple scaling** [DKMW05], when no block is: the `λᵢ` are unknowns
     too, `λ₁ = 1`. With `M` the inverse of a block's first `t` images'
     system, its coefficients are `cⱼ = Σᵢ Mⱼᵢ yᵢ λᵢ`, and each further
     image gives a linear equation for the scales, solved by Gaussian
     elimination. That takes `N` images with `Σₑ (N − tₑ) ≥ N − 1`.
3. One image more than the systems need checks them.
- **Failures.** A point where `a` or `b` loses degree in `x` is bad, and
  one whose gcd has a larger degree unlucky: both move to a new `α`. An
  image of lower degree, a power outside the form, a missing scaling term
  or inconsistent systems mean the form is wrong. Scales that stay free
  after two more images mean a content of the gcd in `x`, which multiple
  scaling cannot separate from the scales (DKMW's example 4,
  `(y + 2) x³ + 12 y² + 24 y`).

**The gcd over `𝔽ₚ`** (LINZIP P, CAS-191):
1. In one variable, Euclid's algorithm.
2. Otherwise `a` and `b` are taken in `F[y][X]`, `y` the last variable
   that occurs, as Brown's algorithm takes them. Their contents in `F[y]`
   come out, and `γ`, the gcd of their leading coefficients in `X`, is the
   leading coefficient each image is scaled to.
3. The first image, at a random point of `y`, comes from this algorithm in
   one variable fewer, and gives the form. The others come from the form,
   each scaled to `γ(v)`.
4. Newton's interpolation in `y` combines them until it stops changing or
   passes `deg γ + min(deg_y a, deg_y b)` points. Its primitive part in `y`
   is then tested by division at a random point of all the variables but
   `x`, univariate [DKMW05]: much cheaper than a division.
5. A wrong form, unlucky points of `y` or a failed test restart from a new
   first image, up to four times; then the algorithm gives up.
- Before all that, the contents of `a` and `b` in `x` come out (by the
  domain's gcd, in fewer variables), as LINZIP requires: a content of the
  gcd in `x` would leave the scales free.
- The points are random (`SplitMix64`, a fixed seed), from fields of at
  least `2^31` elements, where they rarely collide. `GcdMethod::Zippel`
  tests the result by division and returns `None` when it gives up.

## Resultants

- **Over a field**, by Euclid's algorithm:
  `res(a, b) = (−1)^(mn) lc(b)^(m−k) res(b, a mod b)`, with `m`, `n` and `k`
  the degrees of `a`, `b` and the remainder, and `res(a, c) = c^m` for a
  constant `c`.
- **Over ℤ in one variable**, modulo word primes that divide neither
  leading coefficient, combined by `CrtBasis::combine_symmetric`.
  - A prime that divides a leading coefficient drops a degree of the image,
    and the image's resultant is then not the image of the resultant: for
    `p₀ x² + x + 1` and `x − 2` it would be 3 for `4p₀ + 3`.
  - Hadamard's inequality on the rows of the Sylvester matrix bounds
    `|res(a, b)|` by `‖a‖₂^(deg b) ‖b‖₂^(deg a)`. The primes, each above
    `2^61`, run until their product passes twice the bound.
- **Over ℚ and ℚ(i)**, from the numerators:
  `res(n/d, n′/d′) = res(n, n′) / (d^(deg n′) d′^(deg n))`.
- **Elsewhere**, by subresultants.

## Exact division

`SparseRing::div_exact`, the test of every gcd method, now stops at the first
quotient term past the degrees of `a` less those of `b` in some variable
(ours, CAS-137). A quotient has exactly those degrees, so the test answers
the same. But a failed division then never outgrows `a`'s exponents: without
the bound, `x^10` over `x − y^200` reduces down to `y^2000` before it fails.

## Benchmarks (CAS-DOC-8 §13)

`gcd_workload`, an ignored test of `src/gcd/tests.rs`, measured with
`perf stat -e instructions:u,task-clock`, the inputs' construction
subtracted:
- `dense:METHOD:DEGREE:BITS` is the gcd of `f g` and `f h` for pseudo-random
  `f`, `g`, `h` of that degree and coefficient size;
- `sparse:METHOD:N` is the same for `f = (1 + x + y + z)^N + x`, and `g`, `h`
  likewise with `y` and `z`.

giac 1.9's `time(gcd(a, b))` on the same polynomials gives a coarse
comparison.

| Workload | `gcd` | Instructions | Time | Brown's alone | giac |
|---|---|---|---|---|---|
| dense, degree 100, 64 bits | heuristic | 25.3 M | 1.1 ms | 24.9 M | 2.4 ms |
| dense, degree 400, 64 bits | heuristic | 358.9 M | 17.7 ms | 347.7 M, 16.4 ms | 10.4 ms |
| dense, degree 400, 256 bits | Brown's | 620.5 M | 23.4 ms | | 22.3 ms |
| sparse, `N = 5` | heuristic | 9.5 M | 0.5 ms | 60.6 M | 1.0 ms |
| sparse, `N = 8` | heuristic | 117.0 M | 6.4 ms | 336.8 M, 17.3 ms | 1.9 ms |
| sparse, `N = 12` | Brown's | 2.05 G | 105 ms | 1.99 G, 102 ms | 5.1 ms |

- **Dense polynomials** are within a factor 2 of giac.
- **Brown's algorithm in one variable** spends about half its time in the
  final trial division over ℤ, whose integer operations allocate, and a fifth
  in Euclid's algorithm modulo the primes.
- **Several variables.** At `N = 5` the heuristic is faster than giac. At
  `N = 8` and `N = 12` this is 3 and 20 times slower: the heuristic gives up
  at `N = 12`, and Brown's evaluation is dense in every variable. Zippel's
  sparse interpolation (CAS-152) and coefficients accumulated in place
  (CAS-155) are for later tickets.
- **Since CAS-141** the heap's methods are inlined: the sparse gcds take
  9.3 M, 113.9 M and 1.94 G instructions, in the same times. The criterion
  benchmarks of `benches/poly.rs`, `gcd_dense` at degrees 100 and 400 and
  `gcd_sparse` at `N = 5` and `8`, take 1.27 ms, 16.8 ms, 0.39 ms and 5.7
  ms.

## Tests

`src/gcd/tests.rs`:
- **Properties**, on products `f g` and `f h` of random polynomials in one to
  four variables:
  - over ℤ, ℚ, ℤ\[i\], ℚ(i) and 𝔽ₚ (`p` = 2, 3, 65537, `2^61 − 1`);
  - every method's gcd times its cofactors gives the inputs back, `f`
    divides it, it is normalised, and the cofactors are coprime;
  - coprimality is checked without a gcd: by a resultant in the first
    variable, and by the contents in it, recursively;
  - the heuristic, Brown's algorithm and subresultants agree wherever they
    apply;
  - coefficients of several limbs need several primes.
- **Dense polynomials** against Euclid over ℚ, and ℚ\[x\] over one
  denominator against coefficient-wise ℚ.
- **Resultants** against the Sylvester determinant, computed by Gaussian
  elimination over ℚ:
  - the modular, Euclid's and the subresultant resultants, over ℤ, 𝔽ₚ and ℚ;
  - coefficients of several limbs;
  - in several variables, at points where no leading coefficient vanishes;
  - 0 exactly when the gcd has the variable.
- **Targeted cases.**
  - Unlucky primes, the first and the second.
  - Bad points.
  - A field too small for its points.
  - The heuristic's first point failing, and its size limit.
  - A prime dividing a leading coefficient of a resultant.
  - Remainder sequences with gaps, Knuth's example among them, whose
    resultant 260708 is the last of its subresultants.
  - `b² − 4ac` and `−4p³ − 27q²`, and giac's resultants.
  - The normal forms of detail 5, zeros, and cancellation of every method.
- **Zippel's algorithm** (CAS-189 to CAS-191).
  - Transposed Vandermonde systems over word primes from 11 to `2^64 − 59`,
    their solutions and inverses; repeated points and 0.
  - Images from the forms of the gcds of random products, with single and
    multiple scaling; a form missing a term, one with a term too many, one
    without the top power; a content of the gcd; a main variable that is
    not the first.
  - Over `𝔽ₚ`, against subresultants wherever the field is large enough,
    in the gcd properties; non-monic gcds in three variables with single
    and multiple scaling, a content in `x`, and Fateman's gcd at the fourth
    power.

## Not yet

- Zippel's algorithm over ℤ, from the form of a first prime's image, and as
  the default for several variables (CAS-192).
- A modular gcd over ℤ\[i\].
- Modular resultants in several variables.
- The divisor heap of [MP07] for the trial divisions.
- Faster remainder sequences (Lazard's and Ducos's optimisations, the
  half-gcd).
