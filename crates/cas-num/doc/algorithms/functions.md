# Real elementary functions

`ball/fns/` (CAS-79, CAS-80, CAS-81, CAS-86, CAS-95).

## How every function is built

- **Midpoint and radius.** Each function computes its value at the exact
  midpoint `m` as a ball, at a few guard bits above the precision asked
  for. It then widens the result by a bound on the derivative over the
  input ball times its radius `ρ`. This is Arb's approach
  [[Joh17]](references.md#joh17).
- **What is proved by hand.** Only the series sums are: fixed-point sums with
  a run-time error count (`fns/fixed.rs`), and binary-splitting tails
  ([constants.md](constants.md#binary-splitting)). Argument reductions and
  reconstructions run on ball arithmetic, which tracks their errors itself.
- **Choosing a core.** After each function's own reduction (`exp` by
  `log 2`, `sin` and `cos` by `π/2`, `atan` to `[0, 1]`), the working
  precision picks the core:
  - lazy [tables](#tables) up to `TABLE_MAX_PREC = 4600` bits;
  - [bit-burst](#bit-burst) from `BIT_BURST_MIN_PREC = 10^4` bits;
  - the baselines (halvings and a series) between, for arguments the tables
    don't cover, and to build the tables.

  Both thresholds are untuned placeholders (CAS-86); measuring them is
  deferred to cas-tune.

The references for the methods are [[BZ10]](references.md#bz10) (chapter
4) and [[Mul16]](references.md#mul16).

## Fixed-point series

`fns/fixed.rs` (CAS-79). A value is an integer `A` standing for `A / 2^w`,
with an error bound `e` in units of `2^-w`. A sum then costs integer
products and no per-term normalisation. The error bounds are computed at
run time, alongside the values, by the rules written next to each step.

**`hyp_sum`.** It computes `Σ_{j < n} x^j / (q(1) … q(j))` for
`|x| <= 1/2` and small integers `q`:
- `q(t) = t` gives `exp`;
- `-(2t-1)(2t)` gives `cos`, and `-(2t)(2t+1)` gives `sin(x)/x`, both in
  `x^2`.

It uses rectangular splitting, the Paterson–Stockmeyer idea
[[PS73]](references.md#ps73) as Smith applies it to Taylor series
[[Smi89]](references.md#smi89) and as Arb's fixed-point sums do
[[Joh15]](references.md#joh15).

With blocks of `m ≈ √n` terms and `Q(j) = q(1) … q(j)`, the partial sums
`acc_i = Σ_{j >= m i} x^(j - m i) Q(m i) / Q(j)` satisfy, for a block of
`L` terms,

```text
acc_i = [Σ_{l < L} x^l c_l + x^m acc_(i+1) / q(m i + m)] / D,
c_l = Π_{t = l+1}^{L-1} q(m i + t),  D = c_0.
```

The sum costs `m` full products for the powers `x^l` and one per block. The
rest are products and quotients by small integers.

**Errors**, in units of `2^-w`, with every value at most 2 in magnitude:
- **Powers.** Each truncated power adds one unit: `e_l = l - 1`.
- **Joining blocks.** A block join adds `e_A + 2 e_m + 2`.
- **Division.** A truncated quotient by an integer adds one unit to the
  bound divided by that integer.

**`atan_sum`.** It computes `Σ (-y)^j / (2j + 1)` for `|y| <= 1/4`: `atan`
in `y = x^2`, or `atanh` with `y = -x^2`. Its coefficients don't telescope
into products, so each block is summed directly and Horner's rule joins the
blocks, still `m + n/m` full products.

**Proof.** The doc comments of `hyp_sum`, `atan_sum` and `powers`.
**Tests.** Each sum against the exact rational value, with a tightness
check on the bound.

## exp and expm1

`fns/exp.rs` (CAS-79).

**`exp(m)` for the exact midpoint.**
1. Reduce `m = k log 2 + t` with `|t| <= log 2 / 2`. `t` is a ball, from
   `log 2` at enough bits for `k`.
2. Halve: `t = x 2^h` with `|x| <= 2^-r` and `r ≈ p^(1/3)`. Each halving
   saves series terms and costs one squaring back, which balances at the
   cube root.
3. Sum the Taylor series of `exp(x)` in fixed point (`exp_sum`, through
   `hyp_sum`).
4. Square back `h` times and scale by `2^k`, as balls.

At the precisions they serve, the [tables](#tables) or the
[bit-burst](#bit-burst) replace steps 2 and 3.

**Radius.** The input radius `ρ` adds `exp(m) (e^ρ - 1)`.

**Huge arguments.** Above `HUGE_EXP` (`|m| >= 2^65536`), the reduction
would need `log 2` to that many bits. A positive argument gives the
indeterminate ball. A negative one gives `[0 ± 2^(-2^(e - 2))]`, a bound
on `exp` of every point.

**`expm1`.** It stays relative to the result:
- from `|m| >= 1/4`, `exp(m) - 1`, where the subtraction cancels at most
  3 bits;
- below that, `exp` at `-exp(m)` more bits, as `|expm1(m)| >= |m| / 2`;
- below `2^-wp`, the ball `[m, m + m^2]`, since
  `e^m - 1 - m ∈ [0, m^2]` for `|m| <= 1`.

## log and log1p

`fns/log.rs` (CAS-79). "Approximate, then certify."

1. **Approximate.** `y ≈ log(m)` comes from Newton's iteration on `exp` at
   the float level, `y ← y + m e^-y - 1`, doubling the precision each step
   from an `f64` start. This step needs no balls and no error analysis.
   Computing `log` by Newton's method on `exp` is standard
   [[BZ10]](references.md#bz10) (chapter 4).
2. **Certify.** One ball evaluation proves the result:

   ```text
   log(m) = y + log1p(u),   u = m e^-y - 1   (a ball),
   |log1p(v) - v| <= v^2    for |v| <= 1/2,
   ```

   since the series `-v^2/2 + v^3/3 - …` is at most `v^2 / (2 (1 - |v|))`.
   So `log1p([c ± ρ]) ⊆ c ± (ρ + (|c| + ρ)^2)`.

The certificate is itself a Newton step whose error `u^2` is proved, so
Newton only runs to half the precision. The cost is about one `exp` at half
precision and one at full.

**Ranges.**
- **Near 1.** Arguments within `1/4` of 1 go through `log1p` of the exact
  `m - 1`, with the step `y ← y + e^-y (x - expm1(y))`, which stays relative
  to the small result.
- **Tables.** Up to `TABLE_MAX_PREC`, both go to the [tables](#tables), with
  no Newton iteration.
- **Bit-burst.** Above, they follow `exp`'s cores, bit-burst included.

**Radius.** `r / (m - r)` for `log`, and `r / (1 + m - r)` for `log1p`. A
ball reaching 0 (or −1) is indeterminate.

## sin, cos, tan, sin_pi, cos_pi

`fns/trig.rs` (CAS-80, CAS-102).

**`sin` and `cos` at the exact midpoint.**
1. **Quadrant.** `k = round(m / (π/2))` comes from π at 64 bits beyond
   `k`'s.
2. **Reduce.** `t = m - k π/2` as a ball, with π at `p + bits(k) + extra`
   bits. While `t` has fewer than `p + 8` known bits, `extra` grows and `t`
   is recomputed, as in Ziv's strategy [[Ziv91]](references.md#ziv91).
   - `sin(10^100)` needs about 333 more bits of π.
   - A dyadic `m` can't come closer to a multiple of π/2 than π's finite
     irrationality measure allows [[Mah53]](references.md#mah53), so the
     extra stays within a multiple of `m`'s length. The loop is also
     capped.
3. **Halve and sum.** `t = x 2^h` with `|x| <= 2^-r`, `r ≈ p^(1/3)`. Then
   `sin x` and `cos x` come from their Taylor series in `y = x^2`
   (`hyp_sum`), with tails at most `2 |y|^n / (2n)!` and `2 |y|^n / (2n+1)!`.
4. **Double back.** `h` doublings `sin 2a = 2 sin a cos a` and
   `cos 2a = 1 - 2 sin^2 a`, as balls. Each doubles the absolute error, so
   the sums carry `2h` more bits. Then rotate by the quadrant `k mod 4`.

At the precisions they serve, the [tables](#tables) or the
[bit-burst](#bit-burst) replace steps 3 and 4.

**Radius.** Both functions are 1-Lipschitz, so the input radius `ρ` adds
`min(ρ, 2)`. Arguments of exponent above `2^16`, indeterminate inputs, and
radii from 2 give `[0 ± 1]`.

**`tan`.** `sin / cos` at `p + 8` bits. It is indeterminate when the
cosine's ball contains 0.

**`sin_pi`, `cos_pi`.** These need no π for the reduction:
- reduce `x` modulo 2 exactly (it is dyadic), into `[-1, 1)`;
- fold onto `[0, 1/2]` by exact differences;
- take `sin` or `cos` of `π b`, with `b` folded once more to
  `0 < b <= 1/4` through `sin(π b) = cos(π (1/2 - b))`.

So integers and half-integers give exact results, and values next to them
keep their relative accuracy. CAS-102 replaced an earlier `cos_pi` that
added `1/2` rounded at `p + 64` bits and lost all relative accuracy next to
half-integers.

## atan, asin, acos, atan2

`fns/atan.rs` (CAS-80, CAS-102).

**`atan` at the exact midpoint.**
1. **Fold.** `|m| > 1` becomes `±π/2 - atan(1/m)`.
2. **Halve.** `h` steps of `x ← x / (1 + √(1 + x^2))`, the tangent half-angle
   formula, so `atan x` is `2^h` times the result. They run as balls, until
   `|x| <= 2^-r` with `r ≈ p^(1/3) / 3`. That is fewer steps than for `exp`,
   since a step costs a root and a quotient, more than a series term.
3. **Sum.** The series `x Σ (-x^2)^j / (2j + 1)` in fixed point
   (`atan_sum`).

**Radius.** The input radius `ρ` adds `ρ / (1 + (|m| - ρ)^2)`, or `ρ` when
the ball reaches 0.

**`asin`, `acos`.** As ball compositions, with no cancellation near ±1:
`asin x = atan(x / √((1 - x)(1 + x)))`, and
`acos x = 2 atan(√((1 - x)/(1 + x)))` for `x >= 0` or `π - acos(-x)`
below. The exact values `±1` are special cases.

**`atan2(y, x)`.** The quadrant comes from the signs of the balls:
- **Ambiguous cases.** A ball across the branch cut (`x < 0`, `y ∋ 0`), or
  both balls containing 0, gives `[0 ± π]`.
- **Formula.** For `|y| < |x|`, `atan(y / x)` (`± π` on the left). Otherwise
  `±π/2 - atan(x / y)` (CAS-102), so the quotient stays at most 1. A wide
  `x` next to 0 under a larger `y` would otherwise make `y / x` a ball
  across 0.

## pow and root

`fns/pow.rs` (CAS-81).

**`x^y`.** Exact fast paths first:
- an exact integer `y` is `Ball::pow`, by squaring, for any sign of `x`;
- an exact half-integer `n + 1/2` is `sqrt(x)^(2n+1)`;
- an exact base 0 gives 0 or 1.

Otherwise `exp(y log x)` for `x > 0`, with `log x` at as many more bits as
`y log x` has above its binary point, so `exp` sees an absolute error below
its own ulp.

**`root(x, n)`.** The real `n`-th root: `exp(log(x) / n)` for `x > 0`, and
`-root(-x, n)` for odd `n` and `x < 0`.
- An odd root of a ball across 0 is `[0 ± R]`, with `R` the larger root of
  the two ends.
- An even root of a ball reaching below 0 is indeterminate.
- An exact `x` whose root is exact gives the exact ball: `c^n = x` for the
  midpoint `c` rounded to the root's length.

## Hyperbolic functions

`fns/hyp.rs` (CAS-81). These are the forms of fdlibm
[[fdlibm]](references.md#fdlibm), in which nothing cancels near 0:

```text
sinh x  = e (e + 2) / (2 (e + 1)),      e = expm1(|x|), odd
cosh x  = (E + 1/E) / 2,                E = exp(|x|)
tanh x  = e2 / (e2 + 2),                e2 = expm1(2|x|), odd
asinh x = log1p(|x| + x^2 / (1 + √(1 + x^2))), odd
acosh x = log1p(t + √(t (t + 2))),      t = x - 1
atanh x = log1p(2|x| / (1 - |x|)) / 2,  odd
```

Each is evaluated at the exact midpoint, then widened by a derivative bound
on the ball, with `z = |m| + ρ`:

| Function | Radius added |
|---|---|
| `sinh` | `ρ cosh z <= ρ (e^z + 1) / 2` |
| `cosh` | `ρ sinh z <= ρ e^z / 2` |
| `tanh`, `asinh` | `ρ` |
| `atanh` | `ρ / (1 - z)` |

`acosh`, whose derivative is unbounded at 1, is the ball composition of its
formula.

## Tables

`fns/tables.rs` (CAS-86). This is Johansson's table-based argument
reduction for medium precision [[Joh15]](references.md#joh15).

Two lookups shrink an argument from below 1 to below `2^-16`, so the series
that follows needs no halvings, and fewer terms. With
`a = trunc(|t| 2^16) = a1 2^8 + a2` and `0 <= a2 < 2^8`:

```text
exp:      e^t = e^(a1/2^8) e^(a2/2^16) e^x,          x = t - a / 2^16
sin, cos: angle addition of a1/2^8, a2/2^16 and x
log:      log(1 + u) = log(1 + a1/2^8) + log(1 + a2/2^16) + log(1 + u2),
          u1 = (2^8 u - a1) / (2^8 + a1),  u2 = (2^16 u1 - a2) / (2^16 + a2)
atan:     atan x = atan(a1/2^8) + atan(a2/2^16) + atan(x2),
          x1 = (2^8 x - a1) / (2^8 + a1 x),  x2 = (2^16 x1 - a2) / (2^16 + a2 x1)
```

For `log` and `atan`, `a2` is read off `u1` and `x1`.

**Error.** The reductions run on balls and the entries are balls, so
nothing needs proving beyond the identities: every error is in a radius.

**Building the tables.** They are built at first use by the baseline code:
- `exp`, `sin` and `cos` entries as powers (rotations) of the baseline's
  value at one step;
- `log` entries as `2 atanh(c / (2 + c))`;
- `atan` entries by the baseline `atan`.

A table records its precision. A request above it rebuilds at
`max(p, 2 · held)`, as the constants do, up to `TABLE_PREC_MAX`. The step
sizes `2^-8` and `2^-16` are untuned placeholders. At 4600 bits the four
tables take about 900 KB.

## Bit-burst

`fns/bitburst.rs` (CAS-86, CAS-95). Brent's bit-burst algorithm
[[Bre76]](references.md#bre76) [[BZ10]](references.md#bz10) (§4.9), as in
Arb's `arb_exp_arf_bb` and `arb_atan_arf_bb` [[Arb]](references.md#arb).

**Chunks.** An exact `|t| < 1` splits into chunks of doubling width:

```text
t = Σ_j p_j / 2^(b_(j+1)),   b_0 = 0, b_1 = 8, b_(j+1) = 2 b_j,
```

with `|p_j| < 2^(b_(j+1) - b_j)`. So `x_j = p_j / 2^(b_(j+1))` is below
`2^-b_j`, and its series needs about `p / b_j` terms of `b_j`-bit
integers. [Binary splitting](constants.md#binary-splitting) sums each
chunk in `O(M(p) log p)`, and there are `log2(p / 8)` chunks. Chunks stop
a few bits past the working precision, and the rest goes to the
baseline's tiny-argument bound.

**Per function.**
- **`exp`.** `exp(t) = Π exp(x_j)`, each from `Σ x^k / k!`, with tail
  `2 |x|^n / n!` for `|x| <= 1/2`.
- **`sin`, `cos`.** Angle addition over the chunks. Each chunk's series is
  alternating with decreasing terms, so each tail is at most its first term.
- **`atan`.**
  - Each step takes `atan x = atan(c) + atan((x - c) / (1 + x c))`, with `c`
    the next chunk of the current argument. The residual is a ball, so its
    rounding is in its radius.
  - `atan c` comes from `Σ (-1)^k c^(2k+1) / (2k + 1)`, with tail at most
    `|c|^(2n+1)`.
  - Arguments above 1/2 first take `atan x = π/4 - atan((1 - x)/(1 + x))`.
- **`log`.** No series of its own: its Newton iteration and certificate run
  on `exp`.

**Powers of two (CAS-95).** The powers of two of `x = p_j / 2^b` are
binary-splitting shifts, not factors of `q`. `Q` then holds only the odd
part (`k!`, `(2k+1)!` or 1), and the merges multiply integers about half
as long.

**Tests.** `fns/core_tests.rs` runs each core directly against MPFR and
against the baseline, at the precisions it serves.
