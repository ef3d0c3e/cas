# Integration of rational functions

`src/integrate.rs` and `src/integrate/quadratic.rs` (CAS-196 and CAS-197,
planned in CAS-PLAN-3; the design is CAS-DOC-9 §3 in Lific). cas-core's
`Integrate` turns the result into an expression.

## What is computed

For `A / D` over ℚ, `D ≠ 0`, `DenseRing<Rationals>::integrate_rational`
returns a `RationalIntegral`:

- `polynomial`: the antiderivative of the polynomial part, without a
  constant term;
- `rational`: Hermite's rational part, a numerator over a monic denominator
  in lowest terms;
- `logarithms`: the terms `c log S` whose residues `c` are rational, each
  `S` monic and square-free, with whether `S` has a real root;
- `pairs`: the parts whose residues are conjugate quadratic irrationals
  `u ± v √Δ` (CAS-197), in a real form;
- `rest`: `P / E`, `deg P < deg E` and `E` square-free, whose residues are
  algebraic of degree 3 or more, or `None`.

The antiderivative is `polynomial + rational + Σ c log S + Σ pairs + ∫ rest`.

## The polynomial part

`A = Q D + R` by long division; `Q` integrates term by term.

## Hermite's reduction [Her1872; Bro98 §1.2; Bro05 §2.2]

In Mack's linear form, which needs gcds and extended Euclidean algorithms
only, no factorisation:

1. `D⁻ = gcd(D, D′)` and `D* = D / D⁻`, the square-free part.
2. While `D⁻` is not constant: with `D⁻₂ = gcd(D⁻, D⁻′)` and
   `D⁻* = D⁻ / D⁻₂`, solve `B (−D* D⁻′ / D⁻) + C D⁻* = A` with
   `deg B < deg D⁻*` (the two are coprime), add `B / D⁻` to the rational
   part, and go on with `A = C − B′ D* / D⁻*` and `D⁻ = D⁻₂`.
3. What is left is `A / D*`, with a square-free denominator; it can be
   improper, and its polynomial part joins the first.

For `1/x²`: `D⁻ = x`, `D* = x`, `B (−1) + C x = 1` gives `B = −1`, `C = 0`:
`−1/x` and nothing left.

## Logarithms with rational residues [Tra76]

For `P / E` with `E` square-free and `deg P < deg E`, the residues are the
roots of `R(t) = res_x(E, P − t E′)`, the resultant of a polynomial in `x`
and `t` in `x` (cas-poly's sparse resultant in two variables, by
subresultants), and for a root `c`, `c log gcd(E, P − c E′)` is the part of
the antiderivative at the roots of `E` with residue `c` (the
Rothstein–Trager theorem).
- `R` is factored over ℚ (CAS-182). Each linear factor `t − c` gives its
  logarithm.
- With `E_c` the product of their arguments and of the norms of the pairs
  below, `E_r = E / E_c`, and `s E_r + t E_c = 1`, the part left is
  `(P t mod E_r) / E_r`: `P / E` is `(P s mod E_c) / E_c + (P t mod E_r) /
  E_r`, and the first part is the sum of the logarithms and the pairs. The
  residues of the second are the roots of the factors of `R` of degree 3 or
  more, which stay (CAS-DOC-9 §7).

## Pairs of quadratic residues [Tra76; Bro05 §2.8]

An irreducible factor `t² + β t + γ` of `R`, `Δ = β² − 4γ` not a square,
has the roots `α, ᾱ = u ± v √Δ`, `u = −β/2` and `v = 1/2`. Over the field
ℚ(√Δ) (`Quadratic`: pairs `a + b √Δ`, the inverse by the norm
`a² − Δ b²`, which is 0 only at 0), `S = gcd(E, P − α E′)`, monic, is the
argument of the residue `α`, and its conjugate `S̄` that of `ᾱ`. With
`N = S S̄`, over ℚ, `α log S + ᾱ log S̄ = u log N + v √Δ log(S / S̄)`.

- `Δ > 0`: `u log|N| + v √Δ log|S / S̄|`, the two logarithms joined as giac
  joins them, `(1/√5) ln|(2x − 1 − √5)/(2x − 1 + √5)|` for
  `1/(x² − x − 1)`. When `N` has no real root, neither has `S` nor `S̄`:
  real, monic and of even degree, they are positive, and the absolute
  values go, `√2 ln((x² − 2√2 x + 3)/(x² + 2√2 x + 3))` for
  `(8x² − 24)/(x⁴ − 2x² + 9)`.
- `Δ < 0`: with `√Δ = i √−Δ`, `S = A + i B` for real `A` and
  `B = √−Δ B′`, `A` and `B′` over ℚ, and `v √Δ log(S / S̄)` is
  `v √−Δ · i log((A + iB)/(A − iB))`. Rioboo's conversion writes the last
  as `Σⱼ 2 atan(Qⱼ)` for polynomials `Qⱼ`, which have no poles, so the
  antiderivative is continuous on ℝ where a single `atan(A/B)` jumps at the
  roots of `B` [Bro05 §2.8]:
  - if `B` divides `A`, `2 atan(A/B)`;
  - if `deg A < deg B`, the same for `(−B, A)`, a constant apart;
  - otherwise, with `B D − A C = 1` from the extended Euclidean algorithm,
    `(A + iB)(D − iC) = (A D + B C) + i`, so it is `2 atan(A D + B C)` and
    the same for `(D, C)`, of lower degrees.

  `N` is positive, so `u log N` takes no absolute value. The `Qⱼ` are over
  ℚ when `−Δ` is a square, `1/(x² + 1)` giving `atan(x)`; otherwise the
  conversion runs over ℚ(√−Δ), where `A` is rational and `B` a multiple of
  `√−Δ`, and each `Qⱼ` is `√−Δ` times a polynomial over ℚ: the Euclidean
  steps keep one of the pair rational and the other such a multiple.

Lazard, Rioboo and Trager's subresultants give the `S` of every residue
without an extension; with factors of degree 2 at most, one Euclidean
algorithm over ℚ(√Δ), a field of degree 2, gives `S` monic, and the real
form needs `√Δ` anyway. Bronstein's example,
`(x⁴ − 3x² + 6)/(x⁶ − 5x⁴ + 5x² + 4)`, has `R = 45796 (4t² + 1)³` and
`S = x³ + i x² − 3x − 2i`: `atan((x⁵ − 3x³ + x)/2) + atan(x³) + atan(x)`.

## Real roots

`DenseRing<Rationals>::count_real_roots` (`src/sturm.rs`, CAS-198) counts
the distinct real roots of a polynomial in a closed interval with rational
or infinite ends, by Sturm's sequence of its square-free part: `p`, `p′`,
then the negated remainders, each over the absolute value of its leading
coefficient. With `V(a)` the sign changes at `a`, zeros left out, the roots
in `(a, b]` are `V(a) − V(b)`: at a root `r`, `p′(r)` has the sign `p` takes
just after `r`, so `V(r)` is `V` just after it, and a root at `a` is counted
apart. At `±∞` the signs are those of the leading coefficients, times
`(−1)^deg` at `−∞`. A logarithm whose argument has a real root is written
with an absolute value, `ln|x − 1|`; one without, which is positive, as it
is, `ln(x² + 1)`, as giac writes them; cas-core's definite integrals look
for the poles of an integrand between the bounds.

## Tests

`src/sturm/tests.rs` counts roots on the line and between bounds, at the
bounds, repeated, and, in a proptest, those of products of linear factors
with known roots, times `x² + 1` or not, between random bounds.

`src/integrate/tests.rs` checks each result by its derivative at rational
points, exactly, in the field of its residues (ℚ(√Δ) for a real pair,
ℚ(√−Δ) or ℚ for a complex one):
- **Examples**: `1/x²`; `1/(x² (x + 1))` with its rational part and two
  logarithms; `1/(x² − 1)`; `2x/(x² + 1)` without an absolute value;
  `1/(x² + 1)`, an arctangent over ℚ; `1/(x² + x + 1)` over ℚ(√(4/3));
  `1/(x² − x − 1)`, a real pair; `1/(x² − 2)` with its `S`;
  `x/(x⁴ + 1)`, `S = x² + i`; Bronstein's example, three arctangents;
  `(x³ + 1)/(x − 1)` with its polynomial part; `1/(x (x² + 1))`, two
  logarithms; `1/((x − 1)(x² + 1))`, a logarithm and a pair;
  `1/(x² + 1)²`, Hermite's `x/(2(x² + 1))` and a pair;
  `(x² + 1)/(x⁴ + x² + 1)`, whose arctangents of polynomials have no jump
  where giac's `atan((x² − 1)/(x √3))` has one; `1/(x³ − 2)`, whose cubic
  residues stay.
- **Properties**: any numerator over products of linear factors with
  multiplicities and at most one irreducible quadratic: the derivative of
  the result is the integrand, the rational part is in lowest terms over a
  monic denominator, the arguments are monic and square-free, their real
  roots counted right, as those of the norms of the pairs, nothing is
  left.
