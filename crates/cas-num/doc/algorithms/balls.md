# Radii, balls and text

## Magnitudes (`Mag`)

`ball/mag.rs` (CAS-75). The radius type, after Arb's `mag_t`
[[Joh17]](references.md#joh17):

- **Representation.** A value is `man × 2^(exp - 30)`, with a 30-bit
  mantissa in `[2^29, 2^30)` and an `Integer` exponent, or it is 0, or `+∞`.
- **Direction.** Every operation returns an upper bound of the exact
  result. The `*_lower` forms return a lower bound, used for denominators.
- **Error.** One rounding is off by at most `2^-29` relative. Two 30-bit
  mantissas multiply exactly in a `u64`, then one shift rounds the product
  up; sums align, add, and round up the same way.
- **Saturation.** It only ever makes a bound looser: an upper bound that
  can't be represented is `+∞`, and a lower bound of something unknown is 0.

The transcendental bounds below are used in radius formulas and in the
series' tail bounds:

- **`exp`.**
  - Below `2^-31` it is `1 + 2x`, so the next `Mag` above 1.
  - From `2^62` it is `+∞`.
  - Otherwise `y = x 2^-s < 2^-8`, `e^y` from 14 Taylor terms in 126-bit
    fixed point (each rounded up, plus one unit for the tail), then `s <= 70`
    squarings at 128 bits, rounded up. Each squaring doubles the relative
    error and adds `2^-127`, so the final rounding to 30 bits dominates.
- **`expm1`.**
  - Below `2^-31` it is `x (1 + 2^-31)`.
  - Below `2^-8`, the series from `x` on in fixed point.
  - Above, `exp`'s 128-bit value minus 1: as `e^x >= 1 + 2^-8`, the
    subtraction magnifies the error at most 257 times.
- **`log1p`.**
  - Below `2^-31` it is `x`.
  - Otherwise `1 + x = 2^k w` with `w` in `[1, 2]`, and
    `log(1 + x) = k log 2 + 2 atanh(t)` with `t = (w - 1)/(w + 1) <= 1/3`, in
    64-bit fixed point, the tail bounded by a geometric series.
  - From `2^63`, `log x + 1/x`.
- **`pow`.** Binary exponentiation, within a factor `(1 + 2^-29)^(2n)` of the
  exact power.
- **`sqrt`.** Rounded up (or down for `sqrt_lower`).

### `SmallMag`

CAS-100. A `Copy` twin of `Mag` with an `i64` exponent. Each of its
operations is a step-for-step copy of `Mag`'s rounding, so its results are
bit-identical. Ball arithmetic computes radii in `SmallMag` when every
exponent fits an `i64` (all of them in practice) and falls back to `Mag`
otherwise. Nothing then allocates.

## Ball arithmetic

`ball/real.rs` (CAS-76). A ball `[m ± r]`, with a `Float` midpoint and a
`Mag` radius, stands for every real within `r` of `m`. This is
midpoint-radius interval arithmetic [[Moo66]](references.md#moo66)
[[Rum99]](references.md#rum99) in van der Hoeven's ball form
[[Hoe09]](references.md#hoe09), as in Arb [[Joh17]](references.md#joh17).

**The one promise.** Every result contains the exact result for every
choice of points in the inputs.

**Midpoints.** Each midpoint is computed at the requested precision by the
float operation: rounded toward zero, or one of the
[approximate entries](float.md#approximate-entries). The error that costs,
`ε`, is that many ulps of the midpoint and goes into the radius
(`Ball::with_error`). An exact result stays exact, with radius 0.

**Radii.** For inputs `[a ± r]` and `[b ± s]`:

```text
a ± b      r + s + ε
a × b      |a| s + |b| r + r s + ε
a / b      (|a| s + |b| r) / (|b| (|b| - s)) + ε     (lower bounds below)
sqrt(a)    r / (sqrt(a - r) + sqrt(a)) + ε           (lower bounds below)
a^n        by squaring and multiplying balls, at bits(n) + 8 more bits
```

**Powers.** `Ball::pow` squares and multiplies at `bits(n) + 8` bits above
`prec`, then rounds to `prec` (CAS-113). Each squaring doubles the relative
error, so at `prec` itself `x^n` keeps only about `prec - log2 n` bits.
Arb's `arb_pow_fmpz_binexp` adds `bits(n)` the same way
[[Arb]](references.md#arb).

**Indeterminate ball.** `[0 ± ∞]` is the only "no information" state. It
results from a division by a ball containing 0, the square root of a ball
reaching below 0, or a radius that overflows.

**In place (CAS-100).** `add_assign`, `mul_assign` and the others update the
midpoint in place through the floats' `*_round_from` forms, with the radius
in `SmallMag`.

### Predicates

`contains`, `contains_ball` and `overlaps` are exact set predicates.
`sign` and `cmp` are three-valued: `Some` when every point of the balls
agrees, `None` otherwise.

- **`contains`, `contains_ball`, `overlaps`.** Each compares a distance
  `|x - a|` with a radius or a sum of radii. The comparison is first tried
  at 64 bits rounded both ways. Only when that can't decide is it settled
  exactly, through `sum_sign`
  ([float.md](float.md#exact-signs-of-short-sums)), so huge exponents cost
  nothing.
- **`sign`, `cmp`.** `sign` compares `|m|` with `r` exactly. `cmp` is decided
  when the balls are disjoint, or are the same exact value.
- **`contains_rational(q)`.** `q` is bracketed by `p`-bit floats for
  `p = 64, 128, …` until both brackets lie in the ball, or one lies outside
  on its side, or `q` is exact. A dyadic `q` is reached at its bit length,
  and a non-dyadic `q` can't equal an endpoint, since endpoints are dyadic.

### Other operations

- **`union`, `intersection`.** Endpoints rounded outward.
- **`rel_accuracy_bits`.** About `exp(m) - exp(r)`.
- **`round_to_accuracy(p)`.** Rounds the midpoint to nearest at `p` bits and
  grows the radius by the rounding error. The constant cache and `N` use it.

### The ring object `RealBalls`

`ball/real/ring.rs`. Balls at a fixed precision, as a `Ring` only.

It is deliberately not a `Field` or `EuclideanDomain`: algorithms that
branch on exact zero tests (GCD, elimination, degree normalisation) would go
silently wrong on a ball that merely contains 0. `is_zero` means the exact
zero ball. A doc test checks that `RealBalls` doesn't satisfy `Field`.

## Decimal text

`ball/decimal.rs` (CAS-77).

### Printing

`Ball::to_decimal(n, style)` follows Arb's `arb_printn`
[[Arb]](references.md#arb): the digits shown are correct to within one
unit in the last place for every point of the ball. So `1 ± 10^-30` prints
as `1.000…` rather than nothing, which is what a rule printing only
certain digits would give at a digit boundary.

**Algorithm.**
1. Estimate the digits known as `log10(|m| / r)` from `rel_accuracy_bits`,
   plus two, capped at `n`.
2. For `k` digits from there down to 1:
   - with `e10 = floor(log10 |m|)`, within one, from the exponent, form one
     ball product `z = x · 10^(k - 1 - e10)`;
   - `10^s` is raised to its power as a ball, so huge exponents cost a few
     squarings rather than an exact power;
   - let `D` be `z`'s midpoint rounded to an integer;
   - if `|z - D| + rad(z) <= 1`, the `k` digits of `D` are correct, and `D`
     is printed.
3. If no `k` works, print `[± R]` with `R` bounding every point's
   magnitude.

`DecimalStyle::Radius` also prints `[digits ± R]`, where `R` bounds the
distance from the printed decimal to every point (`|z - D| + rad(z)`,
scaled back). Exact decimals print without a radius and drop trailing
zeros.

**Proof.** The doc comments of `digits` and `distance`.

**Shortest decimal (CAS-83, ours).** `Ball::to_shortest_decimal(max)` gives
the decimal with the fewest significant digits inside the ball, if it has at
most `max` of them. The CAS shows it when it has no more digits than the
Arb-style digits above (CAS-DOC-5 §6.1, revised 2026-10-20). So `0.1` read
at 53 bits shows as `0.1`, and the sum of the balls of `0.1` and `0.2` as
`0.3`.
1. The ends `m - r` and `m + r` are exact rationals. `e10 = floor(log10 |m|)`
   is computed exactly: an estimate from the bit lengths, corrected against
   powers of 10.
2. A decimal of `k` significant digits is inside when some integer `n` lies
   in `[ceil((m - r) 10^s), floor((m + r) 10^s)]`, with `s = k - 1 - e10`.
   That is monotone in `k`, since a decimal of `k` digits is one of `k + 1`
   with a trailing zero. So a binary search over `k` in `1..=max` finds the
   fewest, with `O(log max)` scalings of the ends.
3. Of the integers in that range, the one nearest `m 10^s` (ties to even) is
   the result. It always prints with a point (`2.0`, `1.0e-30`), so that it
   reads back as a decimal.

It is the interval form of shortest round-trip printing, where the interval
is the ball itself rather than a float's rounding interval.

**Bound (CAS-148).** The exact rationals of a ball of exponent `e` hold
powers of 10 of about `e` bits: printing takes 5 ms at `10^10000` and 3 s at
`10^1000000`. So past `2^16` bits of exponent, about `10^±19728`, neither
the shortest decimal nor the literal below is looked for, and the display
shows the known digits of `to_decimal`, which powers balls.

**Tests.** `ball/decimal/tests.rs`:
- the result lies in the ball; neither of its neighbours on the grid one
  digit coarser does, and a farther shorter decimal would put one of them in
  the ball;
- a decimal read at 53 bits shows with no more digits than it was typed
  with;
- tables: `0.1`, `0.1 + 0.2`, exact integers, `1/3` (too long for 15 digits),
  scientific forms, a ball around 0, a ball reaching the next power of 10,
  and a tie.

### Parsing

`Ball::from_decimal`. It accepts `[+-]digits[.digits][(e|E)[+-]digits]` and
takes the exact decimal value it spells:

- **Value.** Rounded to nearest at the requested precision, through the
  exact rational, with a radius of half an ulp, or 0 when the value is
  representable.
- **Huge exponents.** An exponent whose exact power of 10 would exceed
  `max(2^20, 4 prec)` bits is applied by ball arithmetic instead, with a
  radius a little above half an ulp.

`Ball::from_scaled(n, e, prec)` is the same from the digits `n` and the
exponent `e` of `n · 10^e`.

### Literals (CAS-148, ours)

The CAS reads a decimal literal as a ball (CAS-DOC-5 §6.4, D1), and writes a
ball back as one (§6.5).

**Precision.** `Ball::literal_bits(d, w)` is the precision a literal of `d`
significant digits is read at, with a working precision of `w` bits:
- `w`, while `d <= ⌈w log10 2⌉ + 1` (17 at 53 bits): the digits that read
  every float of `w` bits back [[Mat68]](references.md#mat68);
- beyond, `⌈d log2 10⌉ + 1`, enough for all of them, so that a long literal
  keeps its digits, as giac and Mathematica keep them.

At those bits the literal's radius, half an ulp, is under one unit of its
last digit, and a guard bit covers the radius of huge exponents. So the
shortest decimal in the ball is the literal, unless the literal ends in
zeros. `Ball::from_literal(n, e, w)` reads `n · 10^e` at that precision.

**Writing.** `Ball::to_literal(w)` is the shortest literal that
`from_literal` reads back to the ball's midpoint. For `d = 1, 2, …` digits,
it tries the two decimals of `d` digits around the midpoint, the nearer one
first, and takes the first that reads back. Both have to be tried: next to
a power of two, the rounding interval is narrower on one side.
- **Termination.** The exact expansion of the midpoint, a dyadic, reads back
  at any precision that holds it. One more digit, a trailing zero, raises
  the precision past the round-trip bound when the expansion alone does not.
- **Length.** A midpoint of at most `w` bits gets at most
  `⌈w log10 2⌉ + 1` digits. Up to that many digits, decimals are read at
  `w` bits, and decimals of that many digits are closer together than the
  floats of `w` bits [[Mat68]](references.md#mat68). That is why
  `literal_bits` keeps `w` up to the round-trip bound: with a precision that
  grew with every digit, each extra digit would also narrow the interval it
  must fall in, and nothing would bound the search.

At 53 bits, the literal of a midpoint that is an `f64` is that `f64`'s
shortest round-trip form, such as `0.30000000000000004`.

**Tests.** `ball/decimal/tests.rs`:
- the precision at the bound and past it, and no exact digits past
  `2^16` bits of exponent;
- tables, against the shortest forms of `f64`s;
- properties: a ball's literal reads back to its midpoint, at every
  precision, and within the bound when the midpoint fits it.

## Hexadecimal floats

`ball/hex.rs` (CAS-77). `Float::to_hex` and `Mag::to_hex` print exact
hexadecimal literals in C99's `%a` form [[C99]](references.md#c99) (§7.19.6.1),
such as `0x1.921fb54442d18p+1`:

- **Printing.** Canonical: a leading `1`, the fraction without trailing
  zeros, and a decimal exponent of any size.
- **Parsing.** Accepts other spellings of the same value too (`0x18p-4`,
  `0X.8P1`, no exponent).

Both directions are exact. They are `FullForm`'s spelling of a ball's
midpoint and radius.
