# Machine balls

`ball/machine.rs`, `ball/machine/dd.rs`, `ball/machine/kernels.rs`
(CAS-87, CAS-88, CAS-103). A `MachineBall` is an `f64` midpoint and radius,
the fast path for precisions up to 53 bits (CAS-DOC-5 §5.4).

**Fallback.** An operation returns `None` when `f64` can't keep its
promise, and the caller redoes it with `Ball` at 53 bits or more (the
CAS-87 contract):
- a midpoint that comes out infinite, subnormal or underflowed to zero;
- a radius that overflows;
- a result `Ball` would give as the indeterminate ball, which a finite
  `f64` radius can't hold.

`MachineBall::to_ball` is exact in the midpoint, since every `f64` is
dyadic.

## Arithmetic

This is midpoint-radius arithmetic in floating point, as in
[[Rum99]](references.md#rum99). The midpoint is the IEEE result, rounded to
nearest. Its rounding error goes into the radius, exactly where an
error-free transformation gives it (see [[Mul18]](references.md#mul18) for
these transformations):

| Operation | Radius |
|---|---|
| `a ± b` | `r + s + \|e\|`, with `e` from the two-sum, exact for finite operands ([[Knu97]](references.md#knu97) §4.2.2) |
| `a × b` | `\|a\| s + \|b\| r + r s + \|e\|`, with `e = a b - fl(a b)`, exact above `2^-900`, else half an ulp |
| `a / b` | `(\|a\| s + \|b\| r) / (\|b\| (\|b\| - s))`, plus half an ulp unless the residual `a - q b` is 0 |
| `sqrt(a)` | `r / (sqrt(a - r) + sqrt(a))`, plus half an ulp unless `a - s^2` is 0 |

**The residual `a b - fl(a b)`.** A fused multiply-add where the target has
one. Otherwise Dekker's two-product, which splits each factor into halves
by Veltkamp's method [[Dek71]](references.md#dek71). Rust's `mul_add` is a
software routine without the `fma` target feature.
- Where exactness isn't guaranteed, the radius takes half an ulp of the
  product instead of the exact residual:
  - a factor from `SPLIT_MAX = 2^995` up, since the split multiplies by
    `2^27 + 1` and can overflow;
  - a product below `EFT_MIN = 2^-900`, whose residual, a multiple of about
    `|p| 2^-106`, could fall below the subnormal range.
- A product that underflows to zero from non-zero factors makes the
  operation fall back.

**Rounding the radius up.** The radius is computed in round-to-nearest,
with every step pushed up by one ulp (`up`). For a correctly rounded
`fl(x)`, `x <= next_up(fl(x))` at any exponent. Lower bounds, in
denominators, step down the same way. Rump, Zimmermann, Boldo and
Melquiond study the general problem of directed bounds in round-to-nearest
[[RZBM09]](references.md#rzbm09). A radius below the normal range is
raised to `f64::MIN_POSITIVE`, which is still an upper bound.

## Double-double arithmetic

`machine/dd.rs`. A value is `hi + lo` with `|lo| <= ulp(hi) / 2`, about 106
bits.

**Algorithms.** Those of Joldes, Muller and Popescu
[[JMP17]](references.md#jmp17), built on the two-sum and the two-product
[[Dek71]](references.md#dek71) [[Knu97]](references.md#knu97). With
`u = 2^-53`, their proved relative error bounds are:

| Operation | Algorithm | Bound |
|---|---|---|
| `Dd::add` | Alg. 6 | `3u²/(1 - 4u)` |
| `Dd::mul_f64` | Alg. 7 | `1.5u² + 4u³` |
| `Dd::mul` | Alg. 10 | `7u²` |
| `Dd::div` | Alg. 17 | `15u² + 56u³` |

**`Dd::sqrt` (CAS-103, ours).** One Newton step from `s = √hi`:
- the residual `r = ((hi - p) - e) + lo` is exact, from `(p, e) =
  two_prod(s, s)`;
- then `fast_two_sum(s, r / (2s))`.

It is within `6u²`, proved next to it.

**`OP_ERR`.** Every bound above is below `OP_ERR = 2^-102 = 16u²`, the one
bound the kernels use. The two-product is exact while nothing overflows or
underflows, which the kernels guarantee by keeping inputs within
`2^-600 <= |x| <= 2^600`.

## Kernels

`machine/kernels.rs`. Each elementary function is an argument reduction
and a fixed polynomial in double-double, with an a priori relative error
budget. The `f64` midpoint then costs half an ulp plus a radius term
`KERNEL_ERR |v|`, with `KERNEL_ERR = 2^-90`. No platform libm is used.

**Constants.** Built once, at first use, from the `Ball` constants and
exact rationals, so they are correct by construction, with no hand-typed
digits:
- `log 2` and `π/2` are split Cody–Waite style
  [[CW80]](references.md#cw80):
  - `log 2 = l0 + l1 + l2 + δ`, with `l0` and `l1` of 42 bits, so `k l_i`
    is exact for `|k| < 2^11`;
  - `π/2` in four parts, with the first two of 33 bits, so `k p_i` is exact
    for `|k| < 2^20`;
- the Taylor coefficients, rounded to double-doubles;
- the table `atan(k/8)`.

**Horner's rule.** The polynomial error is bounded by `γ_n` factors
[[Hig02]](references.md#hig02). For the alternating series it is bounded
against the value, through `sinh/sin` and `cosh/cos` ratios.

| Kernel | Reduction and polynomial | Budget |
|---|---|---|
| `exp` | `m = k log 2 + t` (`k l0`, `k l1` exact, `m.hi - k l0` exact by Sterbenz's lemma), degree-24 Taylor on `\|t\| <= 0.35`, then `2^k`; for `-665 <= m <= 708` | `< 2^-94` |
| `log` | `x = f 2^e`, `f` in `[1/√2, √2)`, `log1p(f - 1)` as `2 atanh(s)`, `s = u/(2 + u)`, degree 22 in `s^2`; plus `e log 2` | `< 2^-94` |
| `sin_cos` | `m = k π/2 + t`, `\|k\| < 2^20`, `None` within `2^-40` of a non-zero multiple of π/2; degree 15 in `t^2` for `sin(t)/t` and `cos t` | `< 2^-94` |
| `atan` | `\|x\| > 1` folds to `π/2 - atan(1/x)`; `c = k/8` nearest, `y' = (y - c)/(1 + y c)`, `\|y'\| <= 1/16`, degree 14 in `y'^2`, plus the table | `< 2^-94` |
| `sinh`, `cosh` | from `E = e^\|m\|` (`(E ± 1/E)/2`); for `\|m\| < 1/2`, `sinh` by the sine series at `-m^2` | `< 2^-92.8` |
| `tanh` | `sinh / cosh`; `±1` within `2 e^-80` from `\|m\| >= 40` | `< 2^-92` |
| `asinh`, `acosh`, `atanh` | `log` where its argument is at least 1.41, `log1p` or the `atanh` series below: no cancellation in either | `< 2^-92` |
| `sin_pi`, `cos_pi` | exact reduction modulo 2 and folding, then `sin_cos` of `π b` | `< 2^-92` |
| `pow`, `root` | `exp(t)`, `t = y log x` or `log(x)/n`: `log x`'s relative error becomes `t`'s absolute one | `(1 + \|t\|) 2^-93` |
| `atan2` | the quadrant from the signs, `atan(y/x)` or `±π/2 - atan(x/y)` as for `Ball` | from `atan` |

Each budget, step by step, is in the doc comment of its kernel. Arguments
outside a kernel's range return `None` from the kernel, and the
`MachineBall` method falls back to `Ball` at 64 bits.

**A pitfall (CAS-103).** Multiplying `e^|m|` by `1/2` through
`Dd::mul_f64` splits the factor, and Veltkamp's split overflows for
`|x| > 2^996` without an FMA. A 50 000-case run found this for `sinh` and
`cosh` near `|m| = 700`. Halving is now an exact scaling (`mul_pow2`), and
the failing input is kept as a regression seed.

**Tests.**
- each kernel against MPFR at random and edge arguments, with the measured
  relative error checked below its budget;
- every `MachineBall` result checked to overlap the 64-bit `Ball` result
  and to contain MPFR's value;
- the double-double operations against MPFR at 400 bits: each operation
  within `OP_ERR`, the error-free blocks exactly.
