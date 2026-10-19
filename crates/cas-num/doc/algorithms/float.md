# Floats

`Float` is the midpoint type of balls (`ball/float/`, CAS-73, CAS-74,
CAS-91, CAS-93, CAS-99). Its operations are exact or correctly rounded in
the five modes below, as MPFR's [[Fou07]](references.md#fou07). Its
"approximate" entries return a one-sided value within a stated number of
ulps, which is all a ball midpoint needs.

## Representation

A value is `±0.m × 2^exp`, with the mantissa `m` in `[1/2, 1)`. Unlike
MPFR, whose exponent is a machine word, the exponent is an `Integer`, as in
Arb's `arf_t` [[Joh17]](references.md#joh17), so it never overflows. Only
finite values exist (CAS-DOC-5 §3.3).

- **Mantissa.** Stored top-aligned in limbs, the top limb's top bit set and
  the lowest limb non-zero. Each value has exactly one representation, so
  `Eq` and `Hash` are structural.
- **Zero.** The empty mantissa, exponent 0, no sign.
- **Length.** Exact values stay short: `3` is one limb, whatever precision
  it is later rounded to. Up to 4 limbs (256 bits) the mantissa is stored
  inline.

## Rounding

`ball/float/round.rs`. One routine, `round_mant`, rounds a top-aligned
mantissa in place to `prec` bits. Every rounded operation uses it.

- **Modes.** `Down` and `Up` round toward and away from zero (MPFR's `RNDZ`
  and `RNDA`), `Floor` and `Ceil` toward −∞ and +∞ (`RNDD`, `RNDU`), and
  `Nearest` rounds ties to even.
- **Decision.** Whether the kept part goes up one unit depends on the mode,
  the sign, the first dropped bit, a sticky bit for any later dropped bit
  (including bits the operation never formed), and the last kept bit
  [[BZ10]](references.md#bz10) (chapter 3).
- **Ternary value.** Each rounded operation returns `rounded.cmp(exact)`,
  MPFR's ternary value. Balls use it to set a radius of 0 when the rounding
  was exact.
- **Two-limb form.** A straight-line copy for a two-limb mantissa rounded to
  at most 64 bits serves the one-limb fast paths.

## Addition, subtraction, multiplication

`ball/float/arith.rs`, `ball/float/small.rs`. Each operation has an exact
form (`add`, `sub`, `mul`) and a correctly rounded one (`add_round`, …).

**General path.** Both operands are aligned in a buffer (on the stack up to
16 limbs), added or subtracted, then normalised and rounded in place by
`round_mant`. Bits shifted out of the buffer only matter through the
sticky bit.

**One limb.** Operands of one limb with `i64` exponents, at precisions up to
64 bits, take a `u128` path, rounded by the straight-line two-limb routine.

**Up to four limbs (CAS-91).** Operands of at most `N <= 4` limbs, at
precisions up to `64 N` bits with `i64` exponents, go through fixed-size
windows, so loops unroll and nothing is zeroed or allocated beyond the
window. The general path stays the reference, and tests compare both with
MPFR.

- **Sum.** `x`, the operand with the larger exponent, sits top-aligned in a
  window of `W = N + 1` limbs, and `y` is shifted right by the exponent gap
  `d`.
  - `y`'s lowest bit stays in the window while `d <= 64 (W - ky)`, so bits
    are lost only from `d > 64` on.
  - Then `y < 2^(64W - 65)` while `x >= 2^(64W - 1)`. The true value is the
    window's integer plus a fraction in `(0, 1)`, and normalising shifts by
    at most one bit either way.
  - The rounding bit lies at bit 62 of the window or above, so lost bits
    only ever count through the sticky bit. A carry out of the window
    shifts right by one bit, and that bit joins the sticky bit for the same
    reason.
- **Product.** The schoolbook product of the `N`-limb mantissas lies in
  `[B^W / 4, B^W)` and is normalised by at most one bit.

**In place (CAS-91, CAS-99).** `add_round_from` and the other `*_from`
forms write into an existing value, reusing its inline mantissa. Ball
arithmetic's in-place operations (CAS-100) use them for the midpoints.

## Division and square root, correctly rounded

`ball/float/divsqrt.rs`. Both reduce to exact integer operations, with the
remainder as the sticky bit, the usual route to correct rounding
[[BZ10]](references.md#bz10) (chapter 3):

- **Division.** The dividend's mantissa gets zero limbs below it until the
  exact integer quotient has at least `prec + 2` bits. The quotient comes
  from the exact integer division `Natural::div_rem`. It is then rounded,
  with "remainder ≠ 0" as the sticky bit.
- **Square root.** The mantissa is scaled by an even power of two to at
  least `2 (prec + 2)` bits. Its integer square root
  (`Natural::sqrt_rem`, Zimmermann's Karatsuba square root
  [[Zim99]](references.md#zim99)) then has `prec + 2` bits or more, and
  "remainder ≠ 0" is the sticky bit.

**Up to four limbs (CAS-99).** Fixed-size versions avoid allocation:
- **Division.** The schoolbook division `kernel::div::sb_div_qr` (Knuth's
  Algorithm D [[Knu97]](references.md#knu97) with the 3-by-2 reciprocal of
  [[MG11]](references.md#mg11)) runs on fixed arrays of `2N + 1` limbs.
- **Square root.** Integer Newton steps run on fixed arrays, from above
  [[BZ10]](references.md#bz10) (chapter 1):
  - the start, `r0 = (isqrt(T) + 1) 2^(t/2)` from the top 128 bits `T` of
    `X`, exceeds `sqrt(X)` by a factor below `1 + 2^-62`;
  - each step `r ← floor((r + floor(X / r)) / 2)` is the floor of an
    arithmetic mean whose geometric mean is `sqrt(X)`, so it never goes
    below `floor(sqrt(X))`, and it roughly squares the relative error;
  - after the steps `r` is within one or two of the root, and steps down
    while `r^2 > X`. If that needs more than four steps, the general path
    takes over.

**Proof.** The doc comments of `div_round`, `sqrt_round` and, for the
fixed-size forms, `small.rs`'s `quotient_n` and `root_n`. **Tests.** Bit
for bit against MPFR in all five modes, at random and edge-case operands.

## Approximate entries

`Float::div_approx`, `sqrt_approx` and `rsqrt_approx` return a value of at
most `prec` bits, on the side toward zero, within
`APPROX_ULPS = 2` ulps:

- **Small operands.** Where the kernel would divide or take the root
  exactly (below its Newton threshold), the result is the correctly
  rounded value toward zero, within one ulp.
- **Otherwise.** Both mantissas are padded to `n = max(ka, kb,
  ceil(prec/64) + 1)` limbs and passed to the kernels of
  [kernels.md](kernels.md). Their errors, a few units of a limb at least 63
  bits below the result's last place, plus the final truncation, stay
  below 2 ulps.
- **Square roots.** An odd exponent takes one more limb, halved, so that
  the exponent stays even and the mantissa at least `B^n / 4`.

Ball arithmetic adds those ulps to the radius (`Ball::with_error`, see
[balls.md](balls.md#ball-arithmetic)).

## Exact signs of short sums

`ball/float/sign.rs` (CAS-93). The ball predicates compare quantities like
`|x - a|` against `r`, whose exponents may be far apart: `2^(10^30) - 1`
must never be built. `sum_sign` finds the exact sign of `Σ ±t_i` without
forming the sum:

- **Clusters.** Sorted by top exponent, the terms are grouped: a term joins
  the current cluster while its top bit is within `g = bits(k) + 2` bits of
  the cluster's lowest bit, for `k` terms.
- **Short sums.** A cluster's exact sum spans at most its mantissa lengths
  plus `k g` bits, so it is short.
- **Decision.** If a cluster's sum is not zero, it is a non-zero multiple
  of `2^L`, while all later terms together are below `2^(L - 2)`. So the
  first cluster with a non-zero sum decides the sign.

## Huge exponents

`EXACT_BITS_MAX = 2^30` bits (CAS-93) is the largest exact result the exact
operations (`add`, `sub`, `to_rational`) build. Past it the `checked_*`
forms return `None` and the others panic, instead of exhausting memory.
Rounded operations, products, scalings and comparisons stay cheap at any
exponent.

## Conversions

- **From `Rational`.** The numerator, scaled to `prec + 2` bits or more, is
  divided by the denominator, then rounded with the remainder as the sticky
  bit.
- **To `f64`.** Round to nearest, ties to even, with IEEE overflow to ±∞ and
  gradual underflow.
- **From `f64`.** Exact.
- **Hexadecimal.** See [balls.md](balls.md#hexadecimal-floats).
