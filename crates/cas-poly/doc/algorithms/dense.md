# Dense univariate polynomials

`src/dense.rs` and `src/dense/` (CAS-135), with ℚ and ℚ(i) in
`src/shared.rs`.

## Representation

A `Dense<E>` is the vector of its coefficients, lowest degree first, without
trailing zeros. The zero polynomial has no coefficient, so two equal
polynomials are equal values, and the degree is the length minus one.

The operations live on the ring object `DenseRing<R>`, after cas-num's
pattern (`cas_num::algebra`): the ring is a value, its elements are
`Dense<R::Elem>`, and the operations are methods on `&DenseRing<R>`.
`DenseRing<R>` is a `Ring` itself, so polynomials over polynomials work, and
a `EuclideanDomain` when `R` is a field, with monic polynomials as the normal
forms. Every result goes through `DenseRing::poly`, which drops trailing
zeros: over ℤ/nℤ the product of two leading coefficients can vanish.

## Products

Each coefficient ring picks its product through the trait `Coefficients`.

| Ring | Product |
|---|---|
| ℤ, ℤ\[i\] | schoolbook for short factors, Kronecker substitution beyond |
| ℤ/nℤ and 𝔽ₚ for a word `n` | schoolbook, then Kronecker substitution modulo the word |
| balls | schoolbook |
| others (ℚ, ℚ(i), big moduli, polynomials) | Karatsuba's method over schoolbook |

### Schoolbook

Each coefficient of `a` times all of `b`, added into the result. It polls the
`Cancel` at each coefficient of `a`.

### Karatsuba's method [KO62]

With `a = a₀ + xʰ a₁` and `b = b₀ + xʰ b₁`,

```text
a b = a₀b₀ + xʰ ((a₀ + a₁)(b₀ + b₁) − a₀b₀ − a₁b₁) + x²ʰ a₁b₁,
```

three products of half the size for four. `h` is half the longer factor,
rounded up. Below 16 coefficients in the shorter factor, schoolbook.

- **Unbalanced factors.** A factor more than twice as long as the other is
  cut in blocks of the shorter's length, each multiplied by it.
- **The middle product's last coefficient.** When the longer factor has an
  odd length and `b₁` is empty, `(a₀ + a₁)(b₀ + b₁)` has one coefficient
  more than its place in the result. That coefficient is zero once `a₀b₀`
  and `a₁b₁` are subtracted, so it is dropped.
- **Balls** (CAS-135) use schoolbook: the subtractions of the middle product
  add radii that schoolbook does not have.

### Kronecker substitution [Kro82, Sch82]

A product in ℤ\[x\] becomes one product in ℤ, so that polynomials get
cas-num's integer products (Toom, NTT). If every coefficient of `c = a b`
lies in `[−2^(w−1), 2^(w−1))`, then `a(2^w) b(2^w) = c(2^w)`, and the
coefficients of `c` are the digits of that integer in base `2^w`, read with
signs.

- **The width.** The coefficients of `c` are sums of at most `min(n, m)`
  products of an `a`-bit and a `b`-bit magnitude, so below
  `2^(a + b + bits(min(n, m)))`. One more bit holds the sign:
  `w = a + b + bits(min(n, m)) + 1`.
- **Packing** (ours, CAS-135). Each coefficient is below `2^w` in magnitude,
  so the positive coefficients pack into one integer, and the negative ones
  into another, at bit offsets `k w`, without carries. One subtraction gives
  `a(2^w)`.
- **Unpacking** (ours, CAS-135). The digits of `|c(2^w)|` in base `2^w` are
  read from the low end, each with the carry of the one before. A digit `d`
  from `2^(w−1)` up stands for `d − 2^w`, and carries 1 into the next. The
  sign of `c(2^w)` then goes on every coefficient. The balanced
  representation is unique within the bounds, so these are `c`'s
  coefficients.
- **Squares** take one integer square.
- **Modulo a word** `p`, the residues are in `0..p`. So the coefficients of
  the integer product lie in `[0, min(n, m) (p − 1)²]`, the digits need no
  sign (`w = 2 bits(p) + bits(min(n, m))`), and each digit is reduced
  modulo `p`.
- **Over ℤ\[i\]** the parts are integer polynomials. With `a = a_re + i
  a_im` and `b = b_re + i b_im`, three products give the real part `a_re
  b_re − a_im b_im` and the imaginary part `(a_re + a_im)(b_re + b_im) −
  a_re b_re − a_im b_im`.
- **Later.** Harvey's multipoint variants [Har09] split the one integer
  product into two or four smaller ones, which removes most of the padding.

### Thresholds

The thresholds were measured in instructions (`perf stat -e
instructions:u`), with the ignored test `products_workload` of
`src/dense/tests.rs` as the workload.

Over ℤ, schoolbook wins while the product of two coefficients fits a word:
cas-num's integers then stay small and never allocate. Past a word,
Kronecker substitution wins almost at once. Schoolbook is used up to:

| Bits of a coefficient product | Coefficients in the shorter factor |
|---|---|
| up to 64 | 12 |
| 65 to 96 | 4 |
| more | 2 |

Modulo a word, the digits of the integer product grow with the modulus, so
the crossover does too:

| Bits of the modulus | up to 20 | 21–32 | 33–62 | 63–64 |
|---|---|---|---|---|
| schoolbook up to | 16 | 20 | 36 | 56 |

Instructions per product (CAS-135), with random signs:

| Length | Bits | Schoolbook | Karatsuba | Kronecker |
|---|---|---|---|---|
| 64 | 64 | 3.29 M | 2.48 M | 0.25 M |
| 256 | 64 | 54.2 M | 25.6 M | 1.41 M |
| 1024 | 64 | | 243 M | 8.61 M |
| 64 | 200 | 5.86 M | 4.14 M | 0.77 M |
| 256 | 200 | 92.5 M | 40.1 M | 5.03 M |
| 1024 | 200 | | 372 M | 30.3 M |

Modulo `2^64 − 59`, schoolbook takes 3.04 M at 256 coefficients and Kronecker
substitution 1.13 M; at 1024, 48.4 M and 7.34 M.

## Division

### Over a field: Knuth's Algorithm D [Knu97, §4.6.1]

Long division by the one inverse of `b`'s leading coefficient: each quotient
coefficient is the leading coefficient of the remainder times that inverse,
and `c xᵏ b` is subtracted. It gives `a = q b + r` with `deg r < deg b`, and
polls the `Cancel` at each quotient coefficient.

### Pseudo-division: Knuth's Algorithm R [Knu97, §4.6.1]

Over a ring without division, such as ℤ or polynomials in further variables
(over any coefficient ring since CAS-137, for the gcd), with `d` the leading
coefficient of `b`, it computes

```text
d^(deg a − deg b + 1) a = q b + r,   deg r < deg b.
```

Each round takes the remainder's top coefficient `c` and sets `q ← d q + c
xᵏ` and `r ← d r − c xᵏ b`. That keeps `d^s a = q b + r` and clears the top
coefficient. Exactly `deg a − deg b + 1` rounds give the power, even when the
top coefficient is already zero. When `deg a < deg b`, the result is `(0,
a)`.

### Exact division

Long division in which each quotient coefficient must be an exact quotient
in the ring, the remainder ending at zero. `None` otherwise.

The coefficient's own remainder must be checked at each step. The
subtraction never updates the top coefficient, which is assumed to vanish.
So `3x / 2x` over ℤ would otherwise pass with the quotient 1, as the test
`exact_division_checks_every_coefficient` shows.

## Evaluation, derivative, content

- **Evaluation** is Horner's rule.
- **The derivative** multiplies coefficient `k` by the image of `k` in the
  ring.
- **The content** is the gcd of the coefficients, in the ring's normal form
  (non-negative over ℤ). It stops early once the gcd is 1.
- **The primitive part** divides by the content, then multiplies by the unit
  that puts the leading coefficient in normal form, so over ℤ it has a
  positive leading coefficient: `p = ± content(p) · pp(p)`. This is the
  normalisation of CAS-DOC-8 §15, detail 5.

## ℚ and ℚ(i) over a shared denominator [FLINT]

A polynomial over ℚ is `num / den`: an integer polynomial over a positive
natural. It is reduced: `den` is coprime to every coefficient of `num`, and
zero is `0/1`. Over ℚ(i) the numerator has Gaussian integer coefficients, and
"coprime" holds against their integer parts. `Shared` and `SharedRing` are
generic over the numerator ring (`Numerators`: ℤ or ℤ\[i\]). Products then
run on the numerators by Kronecker substitution, with one gcd against the
denominator, instead of a rational sum at each step of a schoolbook product.

- **Reduction.** The gcd of the denominator and the integer parts of the
  numerator divides both, stopping as soon as it reaches 1.
- **Sums.** With `g = gcd(den(a), den(b))`, `num(a) · den(b)/g + num(b) ·
  den(a)/g` over `den(a) · den(b)/g`, then reduced.
- **Products.** The numerators' product over the denominators' product,
  reduced.
- **Division** (ours, CAS-135). It pseudo-divides the numerators. With `d`
  the leading coefficient of `num(b)` and `e = deg a − deg b + 1`, `d^e
  num(a) = Q num(b) + R`, so

  ```text
  a = (Q den(b) / (d^e den(a))) b + R / (d^e den(a)).
  ```

  `1 / d^e` is then made the inverse of a natural: over ℤ through the sign of
  `d^e`, over ℤ\[i\] through its conjugate, `1/z = z̄ / |z|²`.
- **Normal form.** Monic: the numerator divided by its leading coefficient
  (made a natural the same way), with `den / lc(num)` as the unit.

## Interruption

`Cancel` (CAS-DOC-8 §11.4) is polled at each row of a schoolbook product,
each node of Karatsuba's recursion, before a Kronecker product (the integer
product itself does not stop), and at each quotient coefficient of a
division. The `_with` methods return `Cancelled`; the `Ring` methods use
`Never`.

## Tests

`src/dense/tests.rs` and `src/shared/tests.rs`:
- **Properties.**
  - The ring laws over ℤ, 𝔽ₚ (word primes, small and large), ℚ, ℤ\[i\],
    ℚ(i) and ℚ\[x\].
  - Division: `a = q b + r` with `deg r < deg b`, `(a b) / b = a`, normal
    forms monic. Pseudo-division and exact division over ℤ and ℤ\[i\],
    exact division only claiming what is true.
  - Horner's rule and the product rule against products.
  - Content and primitive part.
  - Kronecker substitution against schoolbook with coefficients of up to 256
    bits of either sign, and at the bounds: every coefficient of the product
    the largest its width allows.
  - Karatsuba against schoolbook over ℚ and modulo `2^127 − 1`.
  - Word products, zero divisors included.
  - Every operation of the shared-denominator types against coefficient-wise
    arithmetic over `Rationals` and `GaussianRationals`, each result
    reduced.
- **Examples.** Polynomials over polynomials; balls; cancellation.
- **Sabotages** (CAS-135): seventeen changes are caught.
  - Kronecker: no carry, a digit a bit short, the negative part added, word
    digits too narrow.
  - Karatsuba: without `a₁b₁`, blocks misplaced.
  - The Gaussian imaginary part.
  - Division without the inverse; pseudo-division without scaling; exact
    division ignoring a remainder.
  - The primitive part's sign; Horner; the derivative.
  - The shared reduction, the sum, the quotient, the Gaussian conjugate.
