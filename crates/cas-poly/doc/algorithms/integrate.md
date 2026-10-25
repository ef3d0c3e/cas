# Integration of rational functions

`src/integrate.rs` (CAS-196, planned in CAS-PLAN-3; the design is CAS-DOC-9
§3 in Lific). cas-core's `Integrate` turns the result into an expression.

## What is computed

For `A / D` over ℚ, `D ≠ 0`, `DenseRing<Rationals>::integrate_rational`
returns a `RationalIntegral`:

- `polynomial`: the antiderivative of the polynomial part, without a
  constant term;
- `rational`: Hermite's rational part, a numerator over a monic denominator
  in lowest terms;
- `logarithms`: the terms `c log S` whose residues `c` are rational, each
  `S` monic and square-free, with whether `S` has a real root;
- `rest`: `P / E`, `deg P < deg E` and `E` square-free, whose residues are
  not rational, or `None`.

The antiderivative is `polynomial + rational + Σ c log S + ∫ rest`.

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
- With `E_c` the product of their arguments and `E_r = E / E_c`, and
  `s E_r + t E_c = 1`, the part left is `(P t mod E_r) / E_r`: `P / E` is
  `(P s mod E_c) / E_c + (P t mod E_r) / E_r`, and the first part is the
  sum of the logarithms. The residues of the second are the other roots of
  `R`; the logarithmic part of Lazard, Rioboo and Trager will take it
  (CAS-197).

## Real roots

`real_roots` counts the real roots of a square-free polynomial by Sturm's
sequence, `p`, `p′`, then the negated remainders: the sign changes of the
leading coefficients at `−∞` (each signed by the parity of its degree) less
those at `+∞`. A logarithm whose argument has a real root is written with an
absolute value, `ln|x − 1|`; one without, which is positive, as it is,
`ln(x² + 1)`, as giac writes them.

## Tests

`src/integrate/tests.rs`:
- **Examples**: `1/x²`; `1/(x² (x + 1))` with its rational part and two
  logarithms; `1/(x² − 1)`; `2x/(x² + 1)` without an absolute value;
  `1/(x² + 1)` left whole; `(x³ + 1)/(x − 1)` with its polynomial part;
  `1/(x (x² + 1))`, two logarithms; `1/((x − 1)(x² + 1))`, a logarithm and
  a rest; `1/(x² + 1)²`, Hermite's `x/(2(x² + 1))` and a rest.
- **Sturm's count** on polynomials with 0 to 3 real roots and constants.
- **Properties**: any numerator over products of linear factors with
  multiplicities and positive quadratics: the derivative of the result and
  the part left add up to the integrand, the rational part is in lowest
  terms over a monic denominator, the arguments are monic and square-free,
  their real roots counted right, the part left proper; without the
  quadratics, nothing is left.
