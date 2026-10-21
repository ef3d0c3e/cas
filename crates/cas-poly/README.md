# cas-poly

Polynomials for the CAS: their types and algorithms over the ring objects of
cas-num. The crate depends only on cas-num and knows nothing of expressions;
cas-core converts between the two. The design is in Lific page `CAS-DOC-8`;
phase 1 is planned in `CAS-PLAN-13`. The algorithms, with their sources, are
described in [`doc/algorithms`](doc/algorithms/README.md).

## Layout

| Module | Contents |
|---|---|
| `dense` | `Dense<E>`: a univariate polynomial as its coefficients, lowest degree first, without trailing zeros. `DenseRing<R>`: the ring object `R[x]` over a ring object of cas-num. It is a `Ring`, and a `EuclideanDomain` over a field, monic polynomials as normal forms. Polynomials over polynomials work, since `DenseRing` is itself a coefficient ring. Operations: sums, scaling, shifts, powers, products; division with remainder over a field (Knuth's Algorithm D); pseudo-division (Algorithm R) and exact division over a Euclidean domain; evaluation (Horner), derivative, content and primitive part. |
| `dense::mul` | `Coefficients`, each coefficient ring with its product of coefficient slices. The default is Karatsuba's method over schoolbook below 16 coefficients. ℤ and ℤ\[i\] use Kronecker substitution past a threshold that falls with the coefficients' size; words (`ModWord`, `WordPrimeField`) do too, reduced; balls use schoolbook. |
| `dense::kronecker` | Kronecker substitution: signed coefficients packed at `2^w` (positives and negatives apart, one subtraction), one integer product, the digits read back balanced with carries. Also the unsigned version reduced modulo a word, and the measured thresholds. |
| `dense::div` | Division over a field, pseudo-division, exact division. |
| `shared` | ℚ\[x\] and ℚ(i)\[x\] as an integer, or Gaussian integer, polynomial over one positive natural denominator, reduced, as FLINT's `fmpq_poly`. `Shared<E>` with `SharedRing<N>` over `Numerators` (ℤ, ℤ\[i\]); the aliases `RationalPoly`, `RationalPolys`, `GaussianRationalPoly`, `GaussianRationalPolys`. Division pseudo-divides the numerators, then clears the leading coefficient's power, over ℤ\[i\] through its conjugate. |
| `cancel` | `Cancel`, what long loops poll to stop early, with `Never`, the default, and an implementation for `AtomicBool`. `Cancelled` is the error of the `_with` methods. |

## Testing

Property tests (`src/dense/tests.rs`, `src/shared/tests.rs`):
- the ring laws over ℤ, 𝔽ₚ, ℚ, ℤ\[i\], ℚ(i) and ℚ\[x\];
- division (`a = q b + r` with `deg r < deg b`, `(a b) / b = a`),
  pseudo-division and exact division over each;
- Horner's rule and the derivative against products; content and primitive
  part;
- every product algorithm against schoolbook: Kronecker substitution with
  negative 256-bit coefficients and at the bounds of its digits, Karatsuba,
  the word products with zero divisors;
- the shared-denominator types against coefficient-wise arithmetic, each
  result reduced.

`products_workload` is an ignored test, a workload for instruction counts:

```sh
cargo test --release -p cas-poly --lib --no-run
CAS_POLY_WORK=kronecker:256:200:78 perf stat -e instructions:u \
    target/release/deps/cas_poly-* dense::tests::products_workload --ignored --exact
```
