# Complex balls and Gaussian rationals

## Complex balls

`ball/complex.rs` and `ball/complex/` (CAS-82, CAS-102). A `ComplexBall` is
a rectangle `re + im i` of two real balls, Arb's `acb_t` representation
[[Joh17]](references.md#joh17). Every operation is a composition of real
ball operations on the parts, so it inherits their enclosure property.
The work is in choosing formulas in which nothing cancels where it need
not.

### Arithmetic

```text
(a + bi)(c + di) = (ac - bd) + (ad + bc) i                  4 products
                 = (ac - bd) + ((a+b)(c+d) - ac - bd) i     3 products, from 4096 bits
(a + bi) / (c + di) = ((ac + bd) + (bc - ad) i) / (c^2 + d^2)
```

- **Three-product form.** It trades a product for additions and a wider
  radius from the cancellation in `(a+b)(c+d) - ac - bd`, so it is used only
  from `MUL3_PREC = 4096` bits, where products dominate.
- **Real factors.** A real factor multiplies the parts directly.
- **Division.** Indeterminate when the divisor may be 0.
- **Powers.** `z^n` for an integer `n` is computed by squaring.
- **`rel_accuracy_bits`.** Relative to the larger part, as Arb's
  `acb_rel_accuracy_bits`: `max(exp(re), exp(im)) - exp(r)` for the larger
  radius `r`. A part much smaller than the other, such as the inexact 0 of
  a nearly real result, costs no accuracy. `N` measures complex results so
  (CAS-84).

### Elementary functions (CAS-82)

```text
exp(a + bi)  = e^a (cos b + i sin b)
log z        = log|z| + i arg z,   log|z| = log1p((a - 1)(a + 1) + b^2) / 2 near |z| = 1
sqrt z       = t + i b/(2t)  (a >= 0),   |b|/(2t) + i sgn(b) t  (a < 0),   t = √((|z| + |a|) / 2)
sin(a + bi)  = sin a cosh b + i cos a sinh b
cos(a + bi)  = cos a cosh b - i sin a sinh b
atan z       = (log(1 + iz) - log(1 - iz)) / (2i)
z^w          = exp(w log z)
```

- **`log`.** `log|z| = log(a^2 + b^2) / 2`, except near the unit circle,
  where it is `log1p(n - 1) / 2` with `n - 1` computed as
  `(a - 1)(a + 1) + b^2`, which does not cancel.
- **`sqrt`.** Kahan's form [[Kah87]](references.md#kah87), which never
  subtracts. A real ball takes the real square root only when all of it is
  at or above 0. Across the negative real axis the imaginary part covers
  both signs, `[0 ± t]`.
- **`atan`.** Near 0 the two logarithms are `±iz + O(z^2)` and their sum
  cancels to about `|z|` relative, so they are computed with `|log2 |z||`
  more bits.
- **`pow`.** `log z` at as many more bits as `w log z` has above its point,
  up to `2^16`, as for real balls ([functions.md](functions.md#pow-and-root)):
  `exp` turns them into relative error and an error in the angle (CAS-113).

### More functions (CAS-102)

**`tan` and the hyperbolic functions** (`complex/trig.rs`).

```text
tan(a + bi) = (sin a cos a + i sinh b cosh b) / (cos^2 a + sinh^2 b)
            = (4 e s c + i sgn(b) t (t + 2)) / (4 e c^2 + t^2),
              t = expm1(2|b|),  e = t + 1,  (s, c) = (sin a, cos a)
tanh z = -i tan(iz),   sinh z = -i sin(iz),   cosh z = cos(iz)
sin(π(a + bi)) = sin(πa) cosh(πb) + i cos(πa) sinh(πb)
```

- **Why not `sin z / cos z`.** That quotient cancels in its real part for
  large `|b|`. The form above, ours (CAS-102), has only sums of
  non-negative terms, its denominator included: `cos 2a + cosh 2b` would
  cancel next to the poles.
- **Beyond `exp`'s range.** `tan z` is `±i` within `2^(4 - 2|b|)` in each
  part.
- **`sin_pi`, `cos_pi`.** They use the real functions' exact values at
  integers and half-integers. `π b` takes as many more bits as it has
  above its point, since `cosh` and `sinh` turn its absolute error into
  relative error (CAS-113: at a fixed `p + 8`, `|π b| ≈ 551` cost 9 bits).

**Inverse functions** (`complex/inverse.rs`). `asin`, `acos`, `asinh` and
`acosh` follow Hull, Fairgrieve and Tang [[HFT97]](references.md#hft97).
With `X = |x|`, `Y = |y|` for `z = x + yi`:

```text
A = (|z + 1| + |z - 1|) / 2 >= 1,   D = √((A - X)(A + X))
asin z  = atan2(x, D) + i sgn(y) acosh A
acos z  = atan2(D, x) - i sgn(y) acosh A
acosh z = acosh A     + i sgn(y) atan2(D, x)
asinh z = asin's A, D at (|y|, |x|):  sgn(x) acosh A + i atan2(y, D)
acosh A = log1p((A - 1) + √((A - 1)(A + 1)))
```

- **Continuity.** `A` and `D` are continuous in `z`, so the only
  discontinuity across a cut is a sign. A ball across a cut gets `[0 ± v]`
  in that part, which holds both sides.
- **No cancellation.** `A - 1` and `A - X` are computed from HFT's forms,
  with `r = |z + 1|`, `s = |z - 1|`, and using
  `r - (X + 1) = Y^2 / (r + X + 1)`:

  ```text
  X < 1:  A - 1 = (Y^2/(r + X + 1) + Y^2/(s + 1 - X)) / 2
          A - X = (Y^2/(r + X + 1) + s + 1 - X) / 2
  X > 1:  A - 1 = (Y^2/(r + X + 1) + s + X - 1) / 2
          A - X = (Y^2/(r + X + 1) + Y^2/(s + X - 1)) / 2
  ```

  A ball across `X = 1` takes the sums `s ± (X - 1)` of both forms, whose
  cancellation is bounded by the ball's width.
- **Accuracy.** Near 0, `acosh A ≈ |y|` and `atan2(x, D) ≈ x` keep their
  relative accuracy.

`atanh` is direct, with `X = |x|`:

```text
Re atanh z = sgn(x) log1p(4X / ((1 - X)^2 + y^2)) / 4
Im atanh z = atan2(2y, (1 - X)(1 + X) - y^2) / 2
```

Its real part keeps its relative accuracy near 0. The quotient is `>= 0`,
so its ball is clamped to `[0, …]` when it reaches below 0.

**`expm1`, `log1p`, `root`** (`complex/near.rs`).

```text
expm1(a + bi) = (expm1(a) cos b + (cos b - 1)) + i e^a sin b,
                cos b - 1 = -sin^2 b / (1 + cos b)   while cos b > 0
log1p(a + bi) = log1p(a (2 + a) + b^2) / 2 + i atan2(b, 1 + a)
root(z, n)    = |z|^(1/n) (cos(θ/n) + i sin(θ/n)),   θ = arg z
```

- **Near 0.** Both parts of `expm1` and `log1p` stay within a few ulps of
  `|z|`: every term comes from `a` and `b` themselves, not from a rounded
  `1 + z` or `e^z`.
- **`e^a`.** It is `expm1(a) + 1` from `a = -1/2` up, and `exp(a)` below,
  where that sum would cancel.
- **Branch.** The root is the principal one: `(-8)^(1/3) = 1 + i√3`.

**Radii.** The inverse functions are ball compositions. For `tan`, `tanh`,
`sinh`, `cosh`, `sin_pi` and `cos_pi`, composing the formulas on the balls
loses the correlations between the parts: `tan(-0.3 + 4i)` came out 4000
times wider than `|tan'| ρ`. These functions therefore take the value at
the exact midpoint, widened in each part by `d √(ρ_re^2 + ρ_im^2)`
(`ComplexBall::at_midpoint`). Here `d` bounds `|f'|` on the rectangle, from
64-bit balls of the parts, and the mean value inequality for holomorphic
`f` along the segment from the midpoint gives `|f(w) - f(m)| <= d |w - m|`.
The bounds used:

```text
|tan'| = 1 / (cos^2 a + sinh^2 b)      |sinh'| = √(sinh^2 a + cos^2 b)
|cosh'| = √(sinh^2 a + sin^2 b)         |sin_pi'| = π √(cos^2(πa) + sinh^2(πb))
```

### Branch cuts

The cuts are the principal ones of C99 Annex G [[C99]](references.md#c99)
and MPC [[MPC]](references.md#mpc), after Kahan
[[Kah87]](references.md#kah87):
- `log`, `sqrt`, `pow`, `root` on the negative real axis, and `log1p` from
  −1 down;
- `atan` and `asinh` on the imaginary axis beyond `±i`;
- `asin`, `acos` and `atanh` on the real axis beyond `±1`, and `acosh`
  below 1.

**On a cut.** A point exactly on a cut takes the side that a `+0` in the
zero part selects in C99 and MPC:
- `arg` in `(-π, π]` for `log` and `sqrt`;
- the side of `Re z > 0` for `atan` and `asinh`;
- the side of `Im z > 0` for the others.

**Across a cut.** A ball straddling a cut gets an enclosure of both sides:
`arg` is `[0 ± π]`, `sqrt`'s imaginary part `[0 ± t]`, and the HFT
functions `[0 ± v]`.

**Tests.** MPC at `p + 192` bits in both directed roundings, at points
near and on every cut, near 0, and at large imaginary parts. Also
identities (`sin^2 + cos^2 ∋ 1`, `f(f^-1(z)) ∋ z`) and agreement with the
real functions on the real axis.

## Gaussian rationals

`complex/` (CAS-50). `Complex` is ℚ(i), a pair of canonical `Rational`s,
so it is canonical and `Eq`/`Hash` are by value.

- **Arithmetic.** Exact. Division is `a conj(b) / norm(b)`.
- **ℚ(i).** A field (`GaussianRationals`).
- **ℤ[i].** The values with integer parts form a Euclidean domain
  (`GaussianIntegers`) [[HW08]](references.md#hw08) (chapter XII):
  - **Division.** Each part of the exact quotient is rounded to the nearest
    integer, so each part of `a/b - q` is at most `1/2` in magnitude and
    `norm(r) <= norm(b) / 2`.
  - **GCD.** The Euclidean algorithm on that division.
  - **Normal forms.** The associates in the first quadrant (`re > 0`,
    `im >= 0`), through the unit `i^-k` for a value in quadrant `k`.
