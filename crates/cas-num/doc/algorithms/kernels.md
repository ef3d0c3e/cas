# Approximate kernels for ball midpoints

A ball midpoint has to be near the true value, and the radius has to bound
how far it is. It does not have to be the correctly rounded value. These
limb kernels therefore skip the correction steps that exact algorithms need
at the end, and return a **lower bound with a proved error** instead:
`r <= floor(t) <= r + e` for the true value `t`, in units of the result's
last limb, with `e` given by a matching `*_error` function (CAS-DOC-5 §5.2,
CAS-72). Below the thresholds where an approximate algorithm pays, the
dispatchers call the exact one, and `e` is 0.

The float layer ([float.md](float.md#approximate-entries)) turns these into
`Float::div_approx`, `sqrt_approx` and `rsqrt_approx`, within
`APPROX_ULPS = 2` ulps.

## Short product, high half (`mulhi`)

`kernel/mul/mulhi.rs`. For `n`-limb `a` and `b` it returns a lower bound of
`floor(a b / B^n)`, the high half of the product, within
`mulhi_error(n)`.

**Base algorithm.** The schoolbook product restricted to the columns from
`B^(n-1)` up: about half the partial products. The dropped columns `k <=
n - 2` hold at most `k + 1` partial products below `B^2` each, so their sum
is below `(n - 1) B^n`. That gives the error `n - 1`, and a single limb is
exact.

**Divide and conquer.** This is Mulders' short product [[Mul00]](references.md#mul00):
- write `a = a1 B^n1 + a0` with `a1` of `n2 = n - n1` limbs, and `b` likewise;
- take the full product `a1 b1`;
- take two recursive short products of `n1` limbs for the cross terms,
  `mulhi(a1h, b0)` and `mulhi(a0, b1h)`, with `a1h` and `b1h` the top `n1`
  limbs of `a1` and `b1`;
- drop `a0 b0`, which is below `B^(2 n1) <= B^n`.

The error follows the same recursion, `e(n) = 2 e(n1) + 5`.

The cut `n1` is the one `mullo` uses for the short low product, and it
depends on the multiplication algorithm at that size: `n/2` for
schoolbook, `11n/36` for Karatsuba, `9n/40` for Toom-3 and `7n/39` from
Toom-4 on. Hanrot and Zimmermann analyse the choice of cut
[[HZ04]](references.md#hz04).

**Tricks.**
- The full product is taken, and is exact (`e = 0`), in two ranges:
  - up to 2 limbs, where it is faster than schoolbook's half (6.5 against
    7.1 ns at 2 limbs);
  - from `MULLO_MUL_THRESHOLD` (the NTT range), where a short product no
    longer pays.

  The schoolbook and divide-and-conquer thresholds are `mullo`'s, as both
  short products drop the same half.
- The schoolbook loop takes two rows at a time, as `mul_basecase` does. It
  also keeps row `j`'s partial product in column `n - 2`. That is part of a
  dropped column, so the result is still a lower bound within the same
  error.
- Since `r` never exceeds the true value, the sums inside never overflow
  `n` limbs.

**Proof.** `mulhi_error` (schoolbook, divide and conquer, full product).

**Tests.** Property tests against the exact high half at every algorithm,
all-ones operands (which maximise every dropped column) and exact cases.

## Approximate reciprocal (`recip_approx`)

`kernel/div/approx.rs`. For a normalised `n`-limb `D` (top bit set) it
returns `X` with `B^n <= X < 2 B^n` and `X <= floor(B^(2n) / D) <= X + e`.

**Algorithm.** It reuses division's Newton reciprocal (`mu::invert`, in the
integer division code): `D X < B^(2n) < D (X + 2)`. Below
`INV_NEWTON_THRESHOLD` that is the exact `floor((B^(2n) - 1) / D)`, and
`e = 1` because that floor is one below `floor(B^(2n) / D)` when `D` is a
power of two. Above the threshold it is Newton's iteration for `1/D`
[[BZ10]](references.md#bz10) (chapter 3), and `e = 2`.

## Approximate quotient (`div_approx`)

For `n`-limb `A` and normalised `D` it returns `Q < 2 B^n` with
`Q <= floor(A B^n / D) <= Q + div_approx_error(n)`.

**Algorithm.** Below `MU_DIV_THRESHOLD` the quotient is exact (the integer
division). Above, it multiplies by the reciprocal: with `X = B^n + I` from
[`recip_approx`](#approximate-reciprocal-recip_approx),

```text
Q = A + mulhi(A, I),
```

since `A X / B^n = A + A I / B^n`. The error is
`A (B^(2n) / D - X) / B^n < 2`, plus the short product's
`mulhi_error(n) + 1`. In all `A B^n / D - Q < 3 + mulhi_error(n)`, and `Q`
never exceeds the true quotient.

**Proof.** `div_approx_newton`.

## Approximate inverse square root (`rsqrt_approx`)

For an `n`-limb `A` normalised for square roots (`A >= B^n / 4`) it
returns `Y` with `B^n <= Y <= 2 B^n` and
`Y <= floor(sqrt(B^(3n) / A)) <= Y + 1`. That is, `Y / B^n` is
`1 / sqrt(A / B^n)` from below, within 2 ulps.

**Base algorithm.** Newton's iteration for the inverse square root
[[BZ10]](references.md#bz10) (chapter 3), recursive on the length:
- with `h = n/2 + 1` and `l = n - h`, take `Yh` for the top `h` limbs `Ah`;
- with `a = A / B^n` and `y = Yh / B^h = (1 - ε) / sqrt(a)`, take one step

  ```text
  y + y (1 - a y^2) / 2  =  (1 - ε^2 (3 - ε) / 2) / sqrt(a),
  ```

  which is exactly a lower bound;
- below `RSQRT_NEWTON_THRESHOLD`, take instead the exact
  `floor(sqrt(floor(B^(3n) / A)))` from the integer square root
  `Natural::sqrt_rem` [[Zim99]](references.md#zim99).

**Error.** The residual `R = B^(n + 2h) - A Yh^2` is exact, and the step is
`Y = Yh B^l + Yh R / (2 B^(3h))`. `R` is truncated to `floor(R / B^(2h-1))`,
and the correction is floored. Floors move toward −∞ when `R < 0` too, so
`Y` stays below the exact step. By recursion `|ε| <= 2 B^-h`, so
`|R| < B^(n + h + 1)`, and the step is within `12.4 / B` of
`B^n / sqrt(a)`. In all `sqrt(B^(3n) / A) - Y < 1 + 14 / B`. A result below
`B^n` is raised to `B^n`, which the true value exceeds.

**The wrapped residual (CAS-90, CAS-98).** Most of `A Yh^2` is known in
advance: it is `B^(n+2h)` up to the small `R`. So the full product is not
needed. Where the cyclic transform applies, the code computes:
- `A T mod (B^m - 1)` by a cyclic transform of length `m`, for `T = Yh^2`;
- `A T mod B^w` by a short low product of `w` limbs, with `m + w = n + h + 3`.

As `|R|` is known to be small, the two residues determine it: with
`e = R mod (B^m - 1)`, `R = e + j (B^m - 1)` for a small `j`, and `j` is
read from the residue modulo `B^w`. This is the wrap-around trick for Newton
residuals [[BZ10]](references.md#bz10) (chapter 3). GMP's `mpn_invertappr`
does the same for the reciprocal [[GMP]](references.md#gmp).

`m` is the largest transform length up to the span `n + h + 2`. Taking
lengths `3·2^k` as well as `2^k` (CAS-98) keeps consecutive lengths at most
1.5 apart, so the low part `w` stays at most a third of the span, plus
one.

**Proof.** `rsqrt_approx` (the step's error) and `rsqrt_residual_wrapped`
(recovering `R`).

## Approximate square root (`sqrt_approx`)

For the same `A` it returns `S` with `S <= floor(sqrt(A B^n)) <= S + e`.

**Algorithm.** The exact integer square root,
`Natural::sqrt_rem` [[Zim99]](references.md#zim99), at every length:
`sqrt_by_newton` is never true.

The alternative goes through the inverse square root,
`sqrt(A B^n) = A Y* / B^n` for `Y* = sqrt(B^(3n) / A)`. With `Y = B^n + J`
from [`rsqrt_approx`](#approximate-inverse-square-root-rsqrt_approx), it
takes `S = A + mulhi(A, J)`, within `3 + mulhi_error(n)`
(`sqrt_approx_newton`).

Measured (`tune_approx`, CAS-90 and CAS-98), that path took 1.4 to 2.1
times the exact Karatsuba root from 64 to 8192 limbs, and wrapping its
residual saved only 3–5%. It stays implemented and tested, so a faster
transform can switch it on by threshold.
