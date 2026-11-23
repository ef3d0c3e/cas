# The Risch algorithm for transcendental elementary functions

cas-core's `Risch` (CAS-585), Xcas's `risch`: a port, function for
function, of giac 2.0.0.21's `src/risch.cc`, checked first on a sympy
prototype against libgiac (T9.c1). The theory is Bronstein's [Bro05],
chapter 5 (the transcendental Risch algorithm: the tower, Hermite reduction,
the polynomial and logarithmic parts) and chapter 6 (the Risch differential
equation). `Integrate`'s pipeline does not call it: `Risch` is a builtin of
its own.

## What is computed

`Risch[f, x]` (and `Risch[f]`, in `x`, as giac's `_risch`, risch.cc:1018)
for `f` elementary over `ℚ(i)(x)` and parameters, its kernels `exp(u)` and
`ln(u)`, trigonometric functions through exponentials and `atan` through
logarithms. The answer is `prim + Risch[rem, x]`, `∫ f = prim + ∫ rem`, where
`rem` is what the algorithm could not integrate (giac prints
`integrate(rem,x)`); when nothing integrates the call stays. `prim` is
elementary, or uses `Ei`, `Ci` and `Si` for `exp(a x) · R(x)`.

## The top level: `risch` (risch.cc:951) and `risch_lin` (:878)

`risch` rewrites trigonometric functions as exponentials (`trig2exp`, here
`TrigToExp`) and calls `risch_lin` on the whole integrand. What remains goes
through giac's `lin` (lin.cc:134): a sum of `coeff · exp(expo)`, products of
exponentials merged; `risch_lin` integrates each term, then the sum of what
the terms left. A real input whose answer keeps `i` is turned back by
`exp2trig`, and kept only if its imaginary part is 0.

`risch_lin` leaves fractional powers of `x` to the remainder, records each
`atan(g)` to integrate it as `i/2 ln((i + g)/(i − g))` and write it back as
`atan`, and splits the constant part off an exponent (`exp(2x + 1) =
exp(1) exp(2x)`, `exp(1)` a parameter), before the tower.

## The tower: `in_risch_tower` (:56), `risch_tower` (:77)

The kernels of `f` that depend on `x`, recursively (giac's `rlvarx`,
series.cc:2175), sorted by size so that the largest is the top `T`; `x` is
the bottom. Powers with `x` become `exp(b ln a)` (`pow2expln`, subst.cc:1787:
`2^x = exp(x ln 2)`), `ln(exp(u))` is `u` (`ln_exp`, :50), and giac's
`tsimplify` (here `TranscendentalSimplify`) writes the exponentials with as
few generators as it can: `exp(2x)` beside `exp(x)` is `exp(x)²`, `exp(−x)`
is `exp(x)⁻¹`. Each kernel's derivative must lie in the field below it, or
the tower is refused (`in_risch_tower`); of the tower of `f` as given and that
of `tsimplify(f)`, the smaller valid one wins. The derivation on a fraction
of the tower is `D = Σ ∂/∂vᵢ · vᵢ′` (`diff`, :99).

Kernels of one size keep our canonical order (`exp` before `ln`); giac keeps
the order of its own `ratnormal` printout, which is not portable. Rows that
tie may pick a different top and split differently: a departure.

## The integrand in `T`: `in_risch` (:743)

`f` as a rational function of `T` over the field below. The power of `T` in
its denominator (`valuation`) is split off as its own factor, since for an
exponential `T` a negative power is a Laurent part, not a pole (:774). The
rest is decomposed square-free and into partial fractions in `T`.

### Hermite reduction: `hermite_reduce` (:517)

A factor of multiplicity `m > 1`, `num/(cof · fᵐ) · exp(a x)`, is reduced
to multiplicity 1 step by step [Bro05, ch. 5]: from the extended gcd of `f′`
and `f` (`Tegcdpsr`, `Tegcdtoabcuv`) the integrated part `−v/((m − 1)
fᵐ⁻¹)`, and the rest of lower multiplicity.

### Rothstein–Trager: `rothstein_trager_resultant` (:124)

For a square-free `N/D` in `T`, the resultant in `T` of `D` and
`p₁ = N − t D′` [Bro05, ch. 5; Tra76]: its roots `t` are the residues. The
resultant is factored over `ℚ(i)`; each linear factor whose root `r` is free
of `x` gives `r ln gcd(p₁(r), D)`, the gcd normalised as giac's `_gcd`
(monomials out, sign by the leading coefficient), so the printed arguments
are giac's. In an exponential tower a logarithm of degree `d` in `T` adds
`−r d u′` to the polynomial part (`ipshift`). A factor of degree 2 or more,
or a root that depends on `x`, leaves the fraction to the remainder (the
residue criterion: the integral is not elementary in this tower).

### The polynomial part: `risch_poly_part` (:662)

- `x` on top: a polynomial integrates term by term.
- `T = ln(u)`: the coefficients from the highest down, each an `in_risch`
  one level lower, allowed one new logarithm, `ln(u)` itself, whose
  coefficient moves up one degree [Bro05, ch. 5]. A coefficient that fails
  sends the whole polynomial to the remainder.
- `T = exp(u)`: each power `Tʲ`, `j ≠ 0`, solves the Risch differential
  equation `y′ + j u′ y = cⱼ` one level lower (`risch_desolve`, below); `T⁰`
  is an `in_risch` one level lower. When the equation has no solution, `u`
  is linear in `x` and `cⱼ` rational in `x`, the term goes to `Ei`.

## The Risch differential equation: `risch_desolve` (:323), `SPDE` (:238)

`y′ + f y = g`, `y` in the field of the tower [Bro05, ch. 6]:
- the denominator `D` of `y`: for an exponential top, the negative power
  bound `α` from the valuations of `f` and `g` in `T`, `D = T^(−α)`; the
  normal part from the simple factors of `f`'s denominator, their
  Rothstein–Trager residues that are negative integers raising `D`; then
  `D · gcd(gden, gden′)/gcd(c, c′)` (`rischde_simplify`, :187, divides
  `R`, `S`, `T` by their common gcd);
- with `y = N/D` the equation becomes `R N′ + S N = T` on polynomials, and
  a bound `ydeg` on `deg N`: from the degrees of `R`, `S` and `T`, raised
  when the ratio of leading coefficients is an integer (for `x`), an
  `in_risch` of it is a logarithmic derivative (for an exponential) or a
  derivative (for a logarithm, through a change of variable `exp(prim)`);
- `SPDE`, Rothstein's degree reduction [Bro05, ch. 6]: while `deg R > 0`,
  `abcuv` writes `S V + R U = T`, `deg V < deg R`, and the equation recurses
  on `(R, S + R′, U − V′)` with `ydeg − deg R`, `N = R N₁ + V`; with `R`
  constant in `T`, `spde_x` (:197) solves by repeated pseudo-division when
  `S` has degree, `S = 0` integrates `T/R` (through `risch_poly_part` of the
  tower below, as giac does at :292), else each power `k` of `T` is an
  equation one level lower.

## `Ei`, `Ci`, `Si`: `integrate_ei` (:582)

`R(x) exp(a x)`, `R` rational, `a` free of `x`: square-free decomposition
and partial fractions of `R`; the polynomial part in closed form
`exp(a x) Σ (−1)ᵏ p⁽ᵏ⁾/aᵏ⁺¹` (giac calls `_integrate`, the same value);
Hermite reduction of the multiple factors; each linear factor `x − r` over
`ℚ(i)` gives `exp(a r) Ei(a (x − r))`, and for `a` imaginary and `r` real,
`Ci` and `Si` of `|im a| (x − r)`. A factor of higher degree stays.

## The derivative check

giac's risch.cc:693 passes 0 as the allowed logarithm of the last
coefficient in a logarithmic tower, so a nested call loses its one-log
restriction: giac then returns primitives that fail differentiation
(T9.c1's batch rows 19, 30, 116). Passing it there would break
`ln(ln(x))/x` (T9.c2), so the port keeps giac's code and checks every
answer: `D[prim] + rem − f` must be 0 by `Together`, or after `TrigToExp` and
`TranscendentalSimplify` (`verified`); otherwise the call stays. A remainder
that still integrates is integrated again, at most `ROUNDS` = 3 times, before
the check, so the remainder's own `Risch` integrates nothing and
evaluation never loops. The other flaws kept: :292 passes the lower tower to
`risch_poly_part` for `y′ = T`; :486 uses `Rr/Ss`; :503 the updated `R`.
Fixed: :837 leaves the Hermite-reduced fraction, not the original, to the
remainder (giac counted the Hermite part twice).

## Limits

The product of the tower's height, the degree in its top and the degree
bound `ydeg` of the differential equation must not exceed `Limits::terms`
(C0019, `within`); past it the call stays, with a warning, as `Integrate`
does at a limit.

## Departures from giac

- Logarithm arguments over `ℚ(i)`: giac's `ln(i x − 1)`, ours `ln(x + i)`
  where giac's gcd normalisation over Gaussian coefficients is not ported, a
  constant apart.
- Resultant factors of degree 2 are not split by a square root (giac's
  `with_sqrt` factoring): `1/(x² − 2)`, `1/(x³ + 1)`, `1/(x⁴ + 1)` stay, and
  `x³/(x² + x + 1)` keeps that part in its remainder.
- giac's last step `re(halftan(…))` is not ported: `exp(x)/(x² + 1)`
  (giac's `re`/`im` of `Ei(x ± i)` in `tan(1/2)`) and `-pi + atan(exp(x))`
  for `1/(exp(x) + exp(−x))` are ours up to that form and constant.
- `tsimplify` is `TranscendentalSimplify`, without giac's merging of
  logarithms (`decompose`, subst.cc:1522).
- Tower-size ties, above.
- Where giac's primitive is wrong (risch.cc:693), ours stays.

The oracle's departure rows, with giac's text, are in
`crates/cas-giac/tests/risch.rs` `departures`.

## Where the code lives

- `crates/cas-core/src/native/integrate/risch.rs`: `risch` (:951, :1018),
  `answer` and `verified` (the check), `risch_lin` (:878), the tower
  (`rlvarx`, `in_risch_tower`, `risch_tower`, `tsimplify`, `pow2expln_x`),
  `lin` (lin.cc:134), `within` (the limits).
- `risch/parts.rs`: `in_risch` (:743), `risch_poly_part` (:662),
  `hermite_reduce` (:517), `rothstein_trager_resultant` (:124), giac's gcd
  normalisation, `integrate_ei` (:582).
- `risch/rde.rs`: `risch_desolve` (:323), `spde` (:238), `spde_x` (:197),
  `rischde_simplify` (:187).
- `risch/tower.rs`: the tower's polynomials, in all of `lv` (the tower's
  variables and the parameters), through cas-core's `Polynomials` over
  `ℤ[i]`, and its derivation.
- cas-poly's primitives (28162ad, T9.c2): `valuation` and `pseudo_div_rem`
  (`src/sparse.rs`, [sparse.md](sparse.md)), `extended_gcd`, `abcuv` and
  `partial_fractions` in one variable (`src/gcd.rs`, [gcd.md](gcd.md)).

## Tests

`risch/tests.rs`, on T9.c1's probes (`.run/1.3/T9/probe.out`): each answer
equals giac's as a value and its derivative is the integrand; partial
answers split into the primitive and `Risch[rem, x]`, both giac's; what giac
leaves stays. `crates/cas-giac/tests/risch.rs` compares the same rows with
libgiac, and holds the departures.
