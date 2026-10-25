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
| `sparse` | `Sparse<E>`: a multivariate polynomial as its terms in decreasing lexicographic order, coefficients in one vector and exponents packed into words with guard bits (`sparse::layout`), in the narrowest layout that holds them. `SparseRing<R>`: the ring object of `nvars` variables, a `Ring`. Sums merge; products and powers go through Johnson's heap with lazy rows (`sparse::heap`). Division and the exact-divisibility test go through Monagan and Pearce's heap over a Euclidean domain, widening the fields when the guards report an overflow. Evaluation, partial evaluation, substitution, the derivative, degrees, coefficients and leading terms in a variable, the recursive view, and conversions to and from `Dense`. |
| `gcd` | Greatest common divisors with cofactors, `lcm`, contents, resultants and discriminants, over the `GcdDomain`s ℤ, ℚ, ℤ\[i\], ℚ(i) and the prime fields, for `DenseRing`, `SparseRing` and the shared-denominator rings. `gcd` chooses: the heuristic gcd (`gcd::heuristic`), Brown's modular algorithm over ℤ and its dense evaluation over prime fields (`gcd::modular`), Euclid's algorithm over fields; subresultants (`gcd::subres`) are the reference and the fallback. `gcd_by` runs one `GcdMethod`. Resultants modulo primes over ℤ and by Euclid over fields (`gcd::resultant`). |
| `factor` | Factoring over ℚ (CAS-142, planned in CAS-PLAN-1). `Factors`, a constant times powers of primitive, normalised polynomials. `squarefree` for `DenseRing` and `SparseRing` over the `CharacteristicZero` domains ℤ, ℚ, ℤ\[i\] and ℚ(i), by Yun's algorithm, in several variables through the content in a variable (CAS-179). `factor::finite`, crate-internal: monic square-free polynomials over word-size prime fields of odd characteristic split into their irreducible factors, by distinct degrees and then Cantor and Zassenhaus's equal-degree splitting (CAS-180). `factor::hensel`, crate-internal: multifactor Hensel lifting of such factors to a prime power, through a tree of Bézout coefficients and quadratic steps (CAS-181). `factor` for `DenseRing<Integers>` and `DenseRing<Rationals>`: the irreducible factors over ℤ and ℚ, by Zassenhaus's algorithm on each square-free part, from the prime with the fewest local factors among three, lifted past Mignotte's bound and recombined by subsets with degree and constant-term tests (CAS-182). |
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

Sparse polynomials are checked against a naive map from exponent vectors to
coefficients for every operation, over ℤ, ℚ and 𝔽ₚ, with mixed layouts and
targeted cases at the fields' boundaries (`src/sparse/tests.rs`).

Gcds (`src/gcd/tests.rs`): every method on products with a known common
factor, in one to four variables over each domain, multiplied back from its
cofactors, normalised, with coprime cofactors, and the methods against each
other; resultants against the Sylvester determinant; unlucky primes, bad
points and the heuristic's failures as targeted cases.

Square-free decompositions (`src/factor/tests.rs`): every result multiplies
back into its input, with increasing exponents and square-free, pairwise
coprime, normalised parts; random products against Musser's algorithm; exact
decompositions of products of known coprime parts over ℤ, ℤ\[i\] and ℚ(i), in
one to four variables.
Factors over prime fields (`src/factor/tests/finite.rs`): products of
irreducibles found by brute force over 𝔽₃, 𝔽₅ and 𝔽₇ split back into them;
over large word primes the factors of random polynomials multiply back and
pass Rabin's irreducibility test.
Hensel lifts (`src/factor/tests/hensel.rs`): monic, congruent to their
factors, their product times the leading coefficient the polynomial modulo
`pᵏ`, and unique: a lift to `pᵏ` reduces to the one to `pʲ`.
Univariate factorisations (`src/factor/tests/zassenhaus.rs`): they multiply
back into distinct, primitive, positive factors that factor as themselves;
the factorisation of a product is the union of its parts'; Eisenstein
products, cyclotomic polynomials and Swinnerton-Dyer polynomials as known
answers; giac's examples. The fuzz target `poly_factor` checks the same.

`products_workload`, `sparse_workload` (`CAS_POLY_SPARSE=fateman` or
`pearce`, the products of CAS-DOC-8 §13) and `gcd_workload`
(`CAS_POLY_GCD=dense:auto:400:64` or `sparse:modular:8`, gcd.md) are ignored
tests, workloads for instruction counts:

```sh
cargo test --release -p cas-poly --lib --no-run
CAS_POLY_WORK=kronecker:256:200:78 perf stat -e instructions:u \
    target/release/deps/cas_poly-* dense::tests::products_workload --ignored --exact
```

## Benchmarks

`benches/poly.rs` (CAS-141, CAS-DOC-8 §13) times, by criterion, at sizes that
run in seconds: dense products over ℤ with 100- and 200-bit coefficients of
16 to 8192 coefficients (CAS-70's workload), Fateman's and Monagan–Pearce's
sparse products, and gcds of products with a common factor, dense and
sparse. The full sizes are the workloads above; the numbers, with giac's, are
in [`doc/algorithms`](doc/algorithms/README.md).

```sh
cargo bench -p cas-poly --bench poly -- '^dense_mul/200/'
```
