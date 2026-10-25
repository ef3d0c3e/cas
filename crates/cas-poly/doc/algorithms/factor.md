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

## Tests

`src/factor/tests.rs`:

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
