# cas-giac

libgiac, the engine of Xcas, as a test oracle (CAS-DOC-8 §12, CAS-138). Xcas is
the reference for this CAS's defaults, so tests compare with its answers. Test
support only: nothing ships with it.

The library is loaded at run time with `dlopen`, so the crate builds without
it, and checks skip when it is missing (`Giac::shared()` is `None`). The one
`unsafe` module of the workspace, `src/library.rs`, loads it and calls its two
C entry points.

## The library

- **Version.** Written against giac 1.9.0 (Debian's `libgiac0t64`
  1.9.0.93+dfsg2-5), `libgiac.so.0`.
- **Another library.** `CAS_GIAC_LIBRARY=/path/to/libgiac.so` loads another
  copy. A path that does not exist makes every check skip:

  ```sh
  CAS_GIAC_LIBRARY=/nonexistent cargo test -p cas-giac
  ```
- **Entry points.** Both are declared in giac's `gen.h` as `const char *f(const
  char *)`. They return a static buffer and share giac's global state, so
  every call takes one global lock and copies the answer.
  - `caseval` (`Giac::eval`): about 10 ms a call, whatever it computes, so
    about 100 checks a second. Loading takes about 5 ms. It suits tables and a
    few hundred random cases a run.
  - `nws_caseval` (`Giac::eval_fast`): giac's NumWorks calculator variant,
    from about 10 µs a call (0.2 ms for `factor(x^2-2)`, 1.3 ms for a
    rationalised denominator). It suits long runs.
- **State.** It persists between calls. An assignment changes later answers,
  and giac reads `x=1/2` at the top level as one. `Giac::restart` runs giac's
  `restart`. Tests that assign use names no other test uses, since tests run
  in parallel.
- **Output.** giac prints some messages on standard output (`#Primes 2`,
  `Root of …`); they are not part of the answers.

## Formatting

`nws_caseval` formats for the calculator. `eval_fast` undoes it:
- `√(…)` is `sqrt(…)`;
- an exact scalar is followed by `=` and an approximation (`sqrt(3)-sqrt(2)=0.317837245196`, `2^100` with `=1.26765060023e+30`), which is dropped. Lists and matrices have none;
- matrix rows print without commas (`[[1,2][3,4]]`); they get them back;
- a vector prints as a one-row matrix (`[[1,2]]`). It is read back as a
  vector, so the two cannot be told apart.

An equation whose right side is a decimal would lose it.

giac's own notations are rewritten where the language differs (CAS-DOC-7):
`list[…]` and the coefficient lists `poly1[…]` of `coeff(p, x)` are `[…]`,
and `set[…]` is `{…}`.

## Comparing

- `ours(text)` is the value of a text in this CAS. It is parsed and lowered
  with the English names, its decimals become balls at 53 bits as the session
  makes them (held parts too), and it is evaluated in a new context. Since the
  names follow Xcas, giac's answers read back this way.
- `Giac::value(input)` is giac's answer read back.
- `Giac::agrees(input)` compares giac's answer with this CAS's value of the
  same input. `Giac::answers(input, expected)` compares it with the value of
  an expected text.
- Two values agree when they are the same expression, or when this CAS's `==`
  finds them equal. That takes accurate overlapping balls as equal, as
  decimals from the two sides are.
- A difference is an `OracleError::Differs` that shows giac's text and both
  values in the language:

  ```text
  1/3+1/6: giac gives 1/2 (1/2 here), this CAS 1/3
  ```
- giac's `evalf` defaults to 12 digits (Xcas's `Digits`), against 15 here: ask
  for the digits on both sides.

## Tests

`tests/oracle.rs`:
- giac's answers recorded in CAS-DOC-8 (§3, §4, §7, §8: twenty of them), each
  read back and compared with the recorded text;
- agreement where this CAS computes already: arithmetic, indexing,
  derivatives, radicals, `exact`, `evalf` with digits, `irem`;
- `nws_caseval` against `caseval`;
- `restart`;
- the skip path, with a library that does not exist.

`tests/polynomials.rs`, the polynomial builtins (CAS-140):
- giac's answers to each of them, recorded as it prints them, read back and
  compared with this CAS's value of the same input;
- the departures of CAS-DOC-8 §15, each with giac's answer and this CAS's as
  decided: gcds normalised, `resultant(c, b) = c^(deg b)`, `degree(0)`, the
  inferred variable, the degree in `sin(x)`, an expanded `lcm`, one
  expanded denominator a term in `expand`, and the main variable; and giac's
  inert capitalised forms (`Gcd`, `Quo`, `Rem`, `Resultant`, `Factor`),
  which stay unevaluated user calls here (CAS-527);
- a few hundred random inputs a run (128 cases each) to `expand` (products,
  powers, a product over a polynomial, two fractions), `normal` (a sum of two
  fractions) and `gcd` (two products with a common factor). In one variable
  giac's answer is this CAS's; in two, `normal` is compared once this CAS's
  `normal` gives giac's answer the sign of §15 detail 5, and `gcd` up to the
  sign, which giac takes from its inputs. 2000 cases each agree.

`tests/gcd_methods.rs`, the gcd algorithms by name (CAS-526): giac's
`heugcd`, `modgcd`, `psrgcd`, `ezgcd` and `icontent` read back, and the
departures: `psrgcd` with a positive leading coefficient, `ezgcd` where giac
gives up or refuses one variable, with giac's message.

`tests/partial_fractions.rs`, partial fractions and fraction forms
(CAS-529, CAS-540, CAS-544): giac's `cpartfrac`, `convert` with `partfrac`,
`parfrac`, `fullparfrac`, `confrac` and `polynom`, `propfrac`, `f2nd`,
`fxnd` and `ratnormal`, read back and compared after `Expand` and
`Together` of the difference; and the departures, each with giac's answer:
`cpartfrac(1/(x^4+1))`, whose `x^2±i` giac splits and ours keeps whole;
`convert(…,confrac)` as `dfc`: `355/113` is `[3,7,16]` (giac `[3,7,15,1]`),
`-7/3` is `[-3,1,2]` (giac `[2,3]`, the sign lost), `sqrt(2)` the period
`[1,[2]]` (giac 15 quotients), `pi` stays `dfc(pi)` (giac 9 quotients); and
a third argument (giac's `'L'`, which giac assigns) taken and not assigned.

`tests/ndsolve.rs`, `NDSolve` and `BVPSolve` (CAS-DOC-35 §5, CAS-520):
- giac's `odesolve` within `10⁻⁸`, relative;
- `bvpsolve` on grids of 4 and 2 intervals against the closed forms
  (`sin x / sin 1`, `(sinh(1-x) + 2 sinh x) / sinh 1`, `x² + 16/x`), and
  against giac at a step of `1/20`, where giac's error (about `h⁴`) is
  inside the `10⁻⁴ (1 + |y|)` tolerance;
- the departure: at `n = 4` giac's step `1/2` errs by `4·10⁻²` on
  `x² + 16/x`, recorded with its values, where ours is the closed form.

`tests/area.rs`, `area` (CAS-588, CAS-589):
- giac's `area` with a count and a rule within `10⁻¹⁰`, relative: every rule
  by name and its French one (`trapeze`, `point_milieu`, `rectangle_gauche`,
  `rectangle_droit`), `rombergt` and `rombergm` at depths 0 to 3 and 5,
  Simpson's odd counts, `gauss15`, reversed bounds, `sqrt(x)` at 0;
- without a count the exact integral, compared after `Expand` and `Together`
  of the difference; no panels is giac's `undef`, ours `Indeterminate`;
- the departures, each with giac's answer as libgiac prints it: a count
  without a rule (giac's polygon branch answers an equation, ours stays
  unevaluated); `gauss15` below 4 panels on `1/x` over `1..100` (giac a pair,
  its estimate second, ours a ball around `ln(100)`); a symbolic
  integrand (giac's sum of decimals, ours exactly `a/2`, both `1.5` at
  `a = 3`); a negative count (giac steps over `[0,-1]` to `-0.5`, ours
  unevaluated); `rombergt` on a symbolic integrand (giac answers a huge
  number that changes between runs, `6.1e306` or `4.8e306`, so its text is
  not compared; ours `a/2`).

`tests/interval.rs`, `convert(…,interval)` (CAS-592): giac's intervals
`[lo..hi]`, found by their text, each widened by a unit of its last printed
digit, overlap our balls in the same places (`sqrt(2)`, with 20 digits,
`exp(1)` with 30, `1/3`, `[1/3,pi]`, `x+1`, `sqrt(2)*sin(2)+sqrt(2)*cos(2)`
with 20), and the rest, intervals and balls as symbols, is the same
expression. A departure in form only: giac prints an interval, ours a ball.

`tests/solve.rs` and `tests/nsolve.rs`, `newton` and `nSolve` (CAS-591):
- giac's `newton` (start points, step counts that stop it, a tolerance,
  `undef`) and `nSolve` (start points, a range) read back and compared by
  value, our converged roots certified balls, a stopped iterate a ball;
- the departures, each with giac's answer: giac restarts from random
  points, ours from fixed ones, so `newton(x^2-2,x)` and `nSolve(x^2=2,x)`
  are `√2` here and giac's (either of `±√2`) is checked as a root only;
  `nSolve` of a system is `fsolve`'s first real solution here
  (`[-1.93185165258,-0.517638090205]`), giac's
  `[-0.517638090205,-1.93185165258]` unseeded, `[1.93185165258,0.517638090205]`
  after `srand(42)` and `[0.517638090205,1.93185165258]` after `srand(7)`
  in the same call: each solves the equations.

`tests/real_roots.rs`, real roots (CAS-538, CAS-539): giac's `VAS`,
`VAS_positive`, `posubLMQ`, `poslbdLMQ`, `sturmseq`, `sturm` with bounds,
`rationalroot` and `crationalroot` read back and compared by value;
`complexroot` against closed forms within `10⁻¹⁰` (giac's `[a..b]`
intervals do not read back) with giac's count of roots; and the
departures, each with giac's answer: `VAS_positive` of a polynomial with a
repeated root (giac's spurious `[0,1]`), the sign of a `sturmseq` factor
when the highest multiplicity is even, and the order of rational roots
within a multiplicity.

`tests/polytools_three.rs`, polynomial tools III (CAS-540, CAS-543):
giac's `rootof`, `fcoeff`, `homogeneize`, `total_degree`, `factor_xn`,
`bezoutian` and `is_irreducible` read back and compared with ours; and the
departures, each with giac's answer: `rootof` stays a call, as there is no
algebraic-number type (R1.2-003: giac's `normal` and `evalf` of one, and
`undef` for `rootof([1,0])`); a symbolic minimal polynomial is read as its
list and one not monic left (giac errs on both); `factor_xn`,
`total_degree` and `is_irreducible` leave a non-polynomial (giac's
`1+x^-1`, `1`, an error); `factor_xn` with a second argument and
`homogeneize` of a list stay (giac echoes them); `bezoutian` with a
variable stays (giac errs); `is_irreducible(x*y)` and `(x*y+x)` are 0, where
giac's 1 is wrong; and `is_irreducible(... mod 7)` is giac's 1, as the
language does not read `mod`.

`tests/groebner_two.rs`, Gröbner tools II (CAS-541, CAS-542): giac's
`in_ideal` and `gbasis_vars` compared exactly, `syzygy`'s generators as
sets up to constants, `gbasis`'s `rur` element by element (the basis when
there is none); and the departures, each with giac's answer: `in_ideal`
by the reduced basis (giac's 0 for `x^2` in `(x+y,x-y)`), `gbasis_vars`'
ties in the evaluated order and its symbols only (`sin(x)`), `syzygy`
with variables (giac's `[1]`), giac's syzygies of `[x-1,x^2-1]` (missing
one) and `[x^2,x*y,y^2]` (not the minimal ones), another separating form
for `[x^2-1,y^2-1]`, and the calls of giac's F4 knobs `gbasis_reinject`
and `rur_gbasis` staying (`gbasis_max_pairs`, `gbasis_simult_primes` and
`rur_certify` are settings, in `tests/settings.rs`).

`tests/all_trig_solutions.rs` and `tests/with_sqrt.rs`, giac's settings
(CAS-545), one setting a file, as giac's instance is shared by the tests of
a file and giac's input sets and resets it within one call:
- `all_trig_solutions(1)`: giac's general solutions of `solve` (`sin`, `cos`,
  `tan`, `sin(2x)`, `sin(x)^2=1/4`, and `exp(x)=2` unchanged), giac's
  parameters `n_k` renamed `n_0`, `n_1`, … in their order, compared in
  giac's order after `Expand` of the difference (giac writes `nπ + π/2` as
  `-1/2*(-2*n*pi-pi)`); ours in a session with the flag set, and the
  principal solutions back with it off;
- `with_sqrt(0)` (giac refuses `with_sqrt:=0`, a reserved word): giac's
  `factor` over ℚ read back and compared by value; `factor` splits over
  square roots again with it back on;
- the departure: a setting's bare name reads its value here, 0 or 1,
  where giac gives the symbol.

`tests/increasing_power.rs`, Xcas's `increasing_power` (CAS-546, CAS-618),
giac's input setting the flag and back to its default, 0, within one call,
ours run in a session with `increasing_power := 0`:
- sums of powers of one symbol by decreasing degree in `string` (giac's text
  is ours without spaces, the term order compared, and one value read
  back), `latex`, `export_mathml` in its three forms and `mathml` (giac's
  text is ours, whitespace collapsed; the default form's annotation without
  spaces); with the flag on, increasing degree in both;
- the departures, same values: a rational coefficient (giac
  `1/2*x^2-x+1/2`, ours `x^2/2 - x + 1/2`), and sums of two symbols, which
  ours leaves in its canonical order (giac `x^2+y^2+2*x*y`, `x^2+y+1`).

`tests/settings.rs`, Xcas's settings of phase 1.6 (CAS-617, CAS-618), each
input setting and resetting within one call, ours run in a session and
compared value by value with giac's sequence:
- read by name (`approx_mode`, `complex_variables`, `angle_radian`,
  `epsilon` within giac's printed digits) and by `name()`, set by `name(v)`
  or `name := v` (`approx_mode`, `complex_mode`, `keep_algext`,
  `increasing_power`, `angle_radian(0)`, `epsilon:=1e-6`);
- giac's answers of `rur_certify(1)` (0), `gbasis_max_pairs(100)` and
  `gbasis_simult_primes(4)` (the old values, 32768 and 20; solve.cc);
- `cas_setup()` element by element, but its `Digits`, its threads and
  recursion list and its `increasing_power` (departures), and
  `cas_setup([0,0,0])` staying where giac refuses it (Invalid dimension);
- the departures, with giac's text: the bare names `complex_mode`,
  `keep_algext`, `increasing_power`, `eval_level`, `autosimplify` (giac the
  symbol, ours the value); `approx_mode:=1` (giac's library leaves
  `1/3+sqrt(2)` exact, ours a decimal); `autosimplify(factor)` (giac leaves
  `x^2-1`, ours factors); `eval_level(1)` (giac evaluates once, ours to a
  fixpoint); `autosimplify()` (giac's `regroup`, ours none, 0); and giac's
  `xcas_mode(0)`, `maple_mode(0)`, `python_compat(0)`, not builtins here
  (CAS-DOC-7's one grammar): the calls stay.

`tests/session.rs`, Xcas's session commands (CAS-619, CAS-620), each input
after giac's `restart`, ours in a new session, compared value by value:
- `VARS()` by name (`[a,b]` either way, `[a,f,zz]`, the last value only
  where our definition shows nothing), `rm_all_vars()` (`[a]`, then `a`
  bare), `rm_a_z()` (`[a,ab]` gives `[a,6]`), `restart` and `restart()` (the
  names, `[]`; `angle_radian` back to 1), and `nodisp(1+1)` (`"Done"`,
  compared raw);
- `Fill(v, A)` (CAS-639), a store: `A:=[1,2,3];Fill(7,A);A` (`[7,7,7]`
  twice), a matrix filled with `x`, `A` unassigned (`1.1`), `A:=5` (`2`), a
  ragged list `[1,[2,3]]` (`[0,0]`); `Fill(3,[[1,2],[3,4]])` errs in both
  (giac's Bad Argument Value, ours a lowering error);
- the departures, with giac's text, ours checked alone: `ans()` and
  `quest()` after `1+1` (giac's `caseval` keeps no history, `undef`; ours 2
  and `"1+1"`), `1+restart()` (giac 1, ours stays with a warning); and ours
  only: `time(1+1)` one nonnegative decimal (giac `[real, cpu]`),
  `monotonic()` a nonnegative decimal, `version()` starting `cas `.

`tests/ti_names.rs`, Xcas's `type` (CAS-621), `iPart`, `fPart`, `frac`,
`trunc`, the sorting commands and `RandSeed` (CAS-622), `make_symbol`,
`symbol_array`, `LSQ`, `LQ`, `SVD` and `SVL` (CAS-623):
- `type(e)==DOM_X` on both sides (giac prints its answer as the type's name):
  integers of any size, rationals, decimals, complex numbers exact and
  decimal, symbols (`x`, `pi`, `infinity`, `undef`), strings, lists,
  functions (`sin`, `x->x^2`) and expressions; `iPart`, `fPart` and `frac`
  of decimals, rationals, `sqrt(2)` and lists, `[iPart(4.3,sqrt(2))]`
  bracketed (C0023); `trunc` to an integer and to `n` decimals, `n < 0` by
  giac's significant digits, lists and `sqrt(2)+i*sqrt(5)`; giac's decimal
  answers (every `iPart`, `fPart` of decimals, `trunc(x, n)`) compared
  within giac's error, ours exact on exact input (`trunc(sqrt(2),3)`
  707/500 against giac's 1.414, C0006);
- the departures, with giac's text: `DOM_INT+0` (giac 2, ours `DOM_INT`),
  `type(2)` printed (giac `integer`), `type(true)` (giac `integer`, ours
  `DOM_IDENT`), `fPart(x)` (giac `x-sign(x)*floor(x*sign(x))`, ours stays),
  `trunc(x)` (giac `floor(x)`), `trunc(-sqrt(2),3)` (giac floors, -1.415;
  ours toward zero, -707/500).
- `SortA`, `sorta`, `SortD`, `sortd` and `sortperm` (CAS-622): lists,
  sequences of values, matrices and sequences of rows sorted by their first
  row, ties (`SortD` reverses the ascending order, ties too), symbols and
  strings; `sortperm(l)` and `sortperm(l, p)`; `RandSeed(321)` and
  `RandSeed(5)==srand(5)`.
- `make_symbol("abc")`, `symbol_array` with one and two dimensions and a
  `%` template, `LSQ` and `lsq` exact (vector and matrix right-hand
  sides, square and tall); `SVL` (`svl`'s alias) near giac's decimals;
  `SVD` (ours `svd`'s with `s` ascending, giac's layout) and `LQ` compared
  by their products `U diag(s) Vᴴ` and `L Q` within 10⁻⁹ (C0006), `SVD`'s
  `s` too, `LQ`'s `P` giac's identity; departures: `make_symbol("x y")`
  (giac `undef`) and `make_symbol("1a")` (giac `a`) stay,
  `symbol_array("%a",2)` (giac `[0a,1a]`) stays, `QR(A,Q,R)` and
  `LU(A,L,U,P)` store into variables in giac (ours stay),
  `SVL([[1,4],[4,1]])` giac `[5.0,3.0]`, ours `[3,5]`.

`tests/random_variables.rs`, Xcas's `randvar` objects, `sample` and
`randvector` over them, and `randmultinomial` (CAS-624):
- objects (no draws), 16 rows: the law's head and each parameter against
  giac's, `normal(ours − giac)` 0 or within giac's 12 digits (C0006, C0032):
  parameters given, from `mean=`, `stddev=`, `variance=` (`gammad(36,3)`,
  `uniformd(10-2*sqrt(3),10+2*sqrt(3))`, `geometric(1/4)`), giac's decimals
  for decimal input; invalid calls (14 rows) stay
  with giac's message as the warning (libgiac's text checked), and ours only
  `sample(randvar(poisson,3),10^7)` with "too many draws";
- draws, seeded with 42 on both sides in the same input (C0014), never giac's
  digits: supports and sizes of 21 `sample`, `randvector` and
  `randmultinomial` rows; the mean of 10,000 normal draws within 4 standard
  errors; giac's `randvector(n,multinomial,P)`, the counts of each index,
  within `4·√(n·p(1−p))` of `n·p`; ours only, `srand(7)` repeating draws;
- departures: `binomial(n,p)` from moments (giac's law; our `binomial` is the
  coefficient, ours `randvar(binomial,n,p)`), giac's `discreted(…)` alias
  tables in decimals for weights with values (ours `randvar([w…],[v…])`,
  exact), `string(randvar(normal))` (giac `"normald(0.0,1.0)"`, ours exact),
  `uniform` over a range (giac's `'uniformd',a,b` sequence in decimals, ours
  `uniformd(a,b)`), `randvar(multinomial,P)` (giac's sequence does not
  sample; ours the weights form, drawing indices), and `sample` of
  `randvar(student,3)` and `randvar(chisquare,3)` (giac's Bad Argument Value;
  ours real and positive decimals), `randvar(weibull,mean=12.5,variance=1)`
  (giac `weibulld(3.08574940721,13.9803128143)`, ours
  `weibulld(15.3481994802,12.9351890156)`; giac's denominator unsquared,
  CAS-637).

`tests/markov.rs`, Xcas's Markov chains `markov` and `randmarkov` (CAS-628),
`markov`'s four lists bracketed into one (C0023):
- exact rows (10, and three deterministic `randmarkov` paths) agree with
  giac's: classes in giac's Tarjan order, invariant vectors, transient
  classes, absorption; columns summing to 1 transposed; any real `eps`;
- giac's decimals (12 rows) within its 12 digits, or within `eps` = 10⁻³
  where giac iterates to it; ours exact where the input is, balls on decimal
  input; rows summing to 1 within giac's `epsilon` 10⁻¹² (off by 10⁻¹³
  accepted, 10⁻¹⁰ and 10⁻⁶ not, as giac's);
- invalid calls (18 rows) stay with giac's message as the warning: "Not a
  stochastic matrix!", Bad Argument Value, Invalid dimension; ours only, the
  limits ("too many states", "too many draws") and a plot argument, which
  stays without one;
- draws, seeded with 42 on both sides (C0014): paths by length, start and
  steps of positive probability, state 1's share over 10,000 steps within
  0.03 of 3/5 (giac 6043/10001); `randmarkov(v,nt)` by layout and by
  `markov`'s classes of the sizes in `v` with `nt` transient states; ours
  only, a seed repeating draws;
- departures: giac's floats where ours is exact (`markov([[1]])` giac
  `[[1.0]]`; absorption `3812798742493/7625597484987` and, with `eps`
  10⁻³, `2047/4096`, ours 1/2); a transient state feeding a periodic class
  (`markov([[1,0,0,0],[0,0,1,0],[0,1,0,0],[1/2,1/6,1/3,0]])`), where giac's
  iteration never ends, ours only; giac's cut-off for a decimal row's sum
  lies within a few units in the last place of 10⁻¹² (giac keeps fewer bits
  of its doubles: `0.5+1e-12` accepted, `0.5-1e-12` not), ours exactly
  10⁻¹² in `f64`; `randmarkov` with a plot argument (giac a graphic) stays.

`tests/kernel_density.rs`, Xcas's kernel density estimation `kernel_density`
and `kde` (CAS-629):
- the exact form (9 rows) by its values at three points within giac's 12
  digits; ours exact (no ball) on exact data with an exact bandwidth
  (`bandwidth=1/4`, `1/2`), balls for a selected one (plug-in rule, rule of
  thumb) or decimal data; a bare symbol as the variable;
- the list form (17 rows: given, plug-in and rule-of-thumb bandwidths,
  `bins`, `range=a..b`, `range=[a,b]`, a bare `a..b`, `output=list`, past
  1000 data the binned plug-in rule) entry by entry within giac's FFT noise,
  `|ours − giac| ≤ 10⁻⁹ max|entry| + 1.5·10⁻¹² |giac|` (C0084); ours certified
  balls, summed directly;
- the piecewise form (4 rows, `t=0..4` naming the variable, `interp=1`) by
  its values at seven points inside and outside the grid, within the same
  noise (C0049);
- a seeded row (C0014): a 500-draw sample's density on 50 bins has 50
  entries on both sides;
- invalid calls (18 rows) stay with giac's message as the warning, `bins=1`
  with giac's FFT's "Invalid dimension";
- departures: giac's `eval=v` (a value of the piecewise form, or its errors)
  is not ported, the call stays without a warning, as do `interp=k` for
  `k > 1` and `spline=k` (giac's spline interpolation); `bins=1` keeps giac's
  error, which comes from its FFT of one point rather than from the options
  (the grid has no spacing, `d = (b − a)/0`); `kernel_density(...)(2)`
  (applying the exact form) is not a row, its value is the exact form's at
  `x = 2`; `kernel_density(seq(irem(k*k,37),k,1,1000),bins=12)` is not a row:
  giac's FFT noise there reaches 1.3·10⁻⁹ of the largest entry (the
  prototype's port of giac's FFT 7.5·10⁻¹⁰), past the bound above.

`tests/fitdistr.rs`, Xcas's distribution fitting `fitdistr` (CAS-629):
- the Weibull, Cauchy and beta fits (6 rows: probe3's three and three of
  extra.out) within 10⁻⁸ of giac's parameters, ours certified balls (giac
  stops its iteration at a step of 10⁻⁵; probe3's Weibull `λ` differs by
  3.4·10⁻¹⁰);
- departures, giac's mean `Σx/(N − 1)` (optimization.cc:7481) and ours
  `Σx/N`, giac's answer and ours: `fitdistr([1,2,3,4,5],normal)` giac
  `normald(3.75,1.78972763291)`, ours `normald(3,sqrt(10)/2)`;
  `fitdistr([1,2,3,4,5,2,3,3],normal)` giac
  `normald(3.28571428571,1.32149751613)`, ours
  `normald(23/8,sqrt(1218)/28)`; the decimal data's normal giac
  `normald(2.28,0.947248647399)`, ours `normald(1.9,0.850881895447)`;
  `poisson` giac `3.28571428571`, ours `23/8`; `exponential` and
  `geometric` giac `0.304347826087`, ours `8/23`; `gammad` giac
  `gammad(1.8889697753,0.828495515483)`, ours
  `gammad(4.91602524944,2.58738171023)`; `fitdistr([1,2,3,4,5],weibulld)`
  giac `weibulld(2.29380667075,3.39430387666)`, ours
  `weibulld(2.29380667071,3.39429071811)` (giac's start from its mean, and
  its `λ` from the step before the last); `fitdistr([2,2,2],weibull)` giac
  `undef`, ours stays with "Expected a nonzero variance"; and
  `fitdistr([0,0.5],betad)` giac `betad(-0.0,-0.0)` (its start, ln 0 being
  −∞), ours stays with "Unable to compute parameters for Beta
  distribution";
- invalid calls (18 rows: the binomial, uniform, lognormal and negative
  binomial laws, an unknown name, giac's argument and data checks) stay
  with giac's message as the warning;
- seeded rows (C0014): `fitdistr(randvector(1000,weibulld,1/2,1),weibull)`
  within four standard errors of `(1/2, 1)` (0.05 and 0.27; giac `(0.4968,
  0.9409)`), and the normal fit of a log-normal sample's logarithms within
  four of `(5, √2)` (0.18 and 0.127; giac `(5.0334, 1.4515)`), both sides.

`tests/clustering.rs`, Xcas's `cluster` and `kmeans` (CAS-630; giac's
answers in probe4.out and probe8.out, the prototype crates/cas-core/proto/clustering/proto.py,
R1.7-012):
- `cluster` (21 rows: D6, E6, 1-D and one- and two-point data, every
  linkage, `count=`, `output=list` and `output=count`) giac's labels and
  ours as they are;
- `kmeans`' seed-independent rows (9: D6 with no `k`, 2 and `output=list`
  under no seed, `srand(42)` and `srand(1)`; four 2-D points plain and in
  `matrix(…)`, with a rational fifth; `[[1],[2],[10],[11]]`) giac's partition
  and ours up to relabelling, and `kmeans(E6,3,output=count)` = 3 exactly;
- departures, the seed-dependent rows (9: `srand(42);kmeans(D6,3)`, E6 with
  no `k`, `count=2..4`, 3, `count=3`, `output=list`, `srand(7)` and
  `srand(1)`, and `kmeans([1,2,10,11,20,21],3)`, where giac crashes on
  scalar data and is not asked): giac's and ours among the partitions the
  prototype found over 200 seeds;
- invalid calls (28 rows: giac's errors in the probes, among them the
  agreeing stays rows `kmeans(D6,2..4)` "Invalid option",
  `kmeans([[1,2],[3,4]],2)` "Invalid k-range" and `cluster(D6,type=complete)`
  "Expected a string") stay with giac's message as the warning;
- graphics (`kmeans(E6,3,output=plot)`, `cluster(D6,output=plot)` and
  `output=tree`): giac draws, ours stays without a warning;
- seeded rows (C0014): `srand(s)` then `kmeans(E6,3)` twice gives the same
  labels, s = 1, 7, 42, 2024.

`tests/sound.rs`, Xcas's sound clips (CAS-631; giac's answers in probe5,
probe9, probe10, probe11 and crates/cas-core/proto/sound/clip/dur.cas, the prototype
crates/cas-core/proto/sound/clip/proto.py): 49 value rows (properties, `channel_data` with
channels, `matrix`, `range=`, bit depths, `normalize=`, `duration=`, and
`soundsec`; several channels bracketed, C0023) giac's decimals read back
within its error; 16 of giac's errors stay with its message (C0056); 24
departures, giac's text and ours: the printed form, `string` and `about` of
a clip (ours `soundclip(44100, 16, [29204, -29204, 14602, 0])`), `type`,
`+`, `-`, `*`, `A - A` (ours 0), `abs`, `A[0]`, `S[1]`, `size`, `evalf`,
`concat`, `size=3` and the `duration=0.001` / `1/1000` rows (ours 9 samples,
giac 8: its 47-bit doubles, C0088); session rows: a clip prints as that
call, reads back and is its `string`, and `about(x)` after `assume(x>0)` is
`hold([x > 0])`. The clip edits (CAS-632; probe5, probe9 and
crates/cas-core/proto/sound/edit/probe1.cas, the prototype crates/cas-core/proto/sound/edit/proto.py, plain `f64` and
giac's 47-bit doubles alike on every row): 74 value rows of `set_channel_data`,
`stereo2mono`, `mixdown` and `splice` (offsets and fades in samples and
seconds, gain, pan, depths, silent and negative similarities, a second clip
shorter than the overlap), 34 more errors that stay, two departures (giac's
`splice` of clips of other sample rates with an exponent `p ≥ 0`, which
giac mixes and ours refuses) and `clipping`: ours warns once per call when
any sample clips, giac only when the last one does. `resample` (CAS-633;
crates/cas-core/proto/sound/edit/probe_rs.cas, the prototype crates/cas-core/proto/sound/edit/resample.py): 12 shared
errors that stay, six departures (giac's Debian build has no libsamplerate
and answers "libsamplerate is required for resampling audio"; ours answers
the clip, libsamplerate's linear and zero-order hold converters exactly at
the qualities 0 and 1, a Kaiser-windowed sinc of ours at 2 to 4) and
`resample_properties`, ours only: the lengths, a constant kept, a 440 Hz tone
kept and a 23500 Hz tone attenuated at 48000 → 44100 Hz.

`tests/transforms.rs`, Xcas's wavelet transforms `dwt` and `idwt` (CAS-634;
giac's answers in crates/cas-core/proto/wavelet/giac/p1.cas and p2.cas, the prototype
crates/cas-core/proto/wavelet/proto_dwt.py, a port of GSL's wavelet/dwt.c): 37 rows of giac's
decimals read back within its 12 digits (C0006): Daubechies, Haar and
B-spline members, `center`, inverse, padding of a list and of a matrix, the
2-D standard and `image` forms, `[]`, and `idwt(dwt(x))`; giac's rounding
noise (`2.2e-16`) against ours of about `1e-30`. Four Haar rows exact here
(`dwt([1,2,3,4],"haar")` is `[5,-2,-sqrt(2)/2,-sqrt(2)/2]`) against their
exact forms. Five of giac's errors that stay with its message ("Wavelet not
implemented", "Invalid argument 2" for a bare `haar` and for `centre=true`,
"Data must be numeric", "Unsupported wavelet type").

`tests/transforms.rs` also holds Xcas's analytic signals `hilbert` on a
list, `instfreq`, `instphase` and `emd` (CAS-635; giac's answers in
crates/cas-core/proto/wavelet/giac/p3.cas, the prototypes crates/cas-core/proto/wavelet/proto_emd.py and
proto_analytic.py, emd bit-exact with libgsl through ctypes): 25 rows within
giac's 12 digits by real and imaginary parts (C0006, C0010): `emd` bracketed
(C0023), its options, the decimal 41-sample row within giac's noise
(`5.6e-17` where giac's decimal reading differs from strtod's); `hilbert` of
exact, decimal and complex lists and `hilbert(3)`; the list forms of
`instfreq` and `instphase`. Seven expression rows by `Simplify` of the
difference (C0032; giac's `cos(x)^2+sin(x)^2` and `floor(1/2)` are ours
closed), else by values at two points. Departures: `emd([x,1,2])` stays with
"Data must be numeric", three options giac refuses stay, and
`hilbert(sin(t),t)` (giac's `-cos(t)`) stays.

`tests/tail_probabilities.rs`, Xcas's upper tails `UTPN`, `UTPT`, `UTPC` and
`UTPF` (CAS-625), and `cdf`, `icdf` and `mgf` of a law by its name, of a
`randvar` object or of a list (CAS-626):
- exact and symbolic rows (`UTPN(0,1,x)`, `cdf(normald,0,1)`, `cdf(exp,2,1)`,
  every law's `mgf` in `t` and giac's `undef` for five of them), the
  difference expanded 0 or the two within 10⁻²⁰ at a rational point; the
  calls giac leaves (`UTPC(0.5,2)`, `UTPF(2.5,2,1)`), which ours leaves;
  `icdf` of a list at `p` outside `(0, 1)`, `undef`;
- giac's decimals within its 12 digits (C0006), ours exact where the input
  is: the four tails (`UTPC(2,6)`, ours `exp(-3)`), integral decimal
  degrees, `cdf` and `icdf` of the laws and of `randvar` objects, `icdf`
  of a list (giac `1.0`, ours `1`);
- invalid calls (20 rows) stay with giac's message as the warning: Bad
  Argument Value, and "Invalid dimension" for `UTPT`'s degree below 1 and
  `UTPF`'s outside `1 … 300`;
- departures: a list's `cdf` (giac's dichotomy answers an index, `(#{Lᵢ ≤
  x} − 1)/n`, `−1/n` below the least value, steps `[x₍ᵢ₎, i/n]`; ours the
  empirical distribution function `#{Lᵢ ≤ x}/n`, steps `(i + 1)/n`, which
  giac's own `icdf` inverts), `icdf(normald,0,1,1/2)` (giac
  `-1.88884303598e-19`, ours 0), `cdf` of `gammad` and `betad` (giac keeps
  `igamma` and `Beta`, ours evaluates), of `negbinomial` (giac `undef(…)`,
  ours `negbinomial_cdf`, `negbinomial_icdf`), and `UTPT(3,x)`,
  `UTPC(2,x)` (giac stays, ours `1 − F` expanded);
- `logistic` (CAS-627), exact and symbolic as giac's, its decimals within
  12 digits, Invalid dimension for no argument or more than four;
  `logistic_regression` (CAS-627), the six entries one by one: the
  rational ones (the total, the peak production) the same value, the
  others (the two curves in `x`, sampled at four rationals; the time of
  the peak; the correlation) within giac's 12 digits; constant quotients
  `[undef,undef,infinity,infinity,undef,undef]` as giac's; a time that is
  no number stays with Bad Argument Value.

`tests/angle_radian.rs`, `tests/complex_mode.rs`,
`tests/complex_variables.rs` and `tests/epsilon.rs`, the settings that change
answers (CAS-617), one group a file, giac's flag set and reset within one
call, ours in a session:
- `angle_radian:=0` (degrees) and `:=2` (grads): numbers converted in the
  trigonometric functions and their inverses (`sin(30)`, `asin(1/2)`,
  `sec(60)`, `acot(0)`, decimals), `sin(pi)`, `sin(x+30)`,
  `asin(x)`, `atan(x)` staying, the hyperbolic functions and `exp` unchanged;
  the departure `sin(15)` in degrees (giac `(sqrt(6)-sqrt(2))/4`, ours stays:
  no closed form of `sin(pi/12)`);
- `complex_mode:=1`: `solve` as `csolve` and `factor` as `cfactor`, `[]`
  again after the reset;
- `complex_variables:=1`: `conj`, `re`, `im`, `abs` of a symbol staying; the
  departures `abs(x)^2` (giac `x^2`, ours stays) and `solve(x^2+1=0,x)`
  (giac `[]`, ours `[i,-i]`: the setting is ours' complex domain);
- `epsilon:=1e-3`: `epsilon2zero` and `trim`'s default threshold.

`tests/complex_solving.rs`, complex roots by number and real roots
(CAS-536, CAS-537):
- giac's `cfsolve` (polynomials, `exp(x)=-1`, `sin(x)=2`, a system, start
  points) read back and compared with ours as sets within `10⁻⁹ (1 + |b|)`,
  giac's 12 digits; ours against closed forms first (`iπ`,
  `π/2 ± i ln(2+√3)`, the eighth roots of unity, `i√2/2`);
- `zeros` and `cZeros`, and `surd` (`surd(-8,3)` is `-2`, `surd(x,3)` stays,
  `surd(-16,4)` by value) where giac's answer is ours;
- the departures, each with giac's answer: `zeros` keeps `Solve`'s order
  (giac's `[1,-1]`, `[3,2,1]` and its tuples reversed), `cZeros` keeps
  `CSolve`'s (giac's `[-i,i]`), `cfsolve` orders by real then imaginary
  part (giac's `x^4+1` and `sin(x)=2` orders), and `cfsolve(exp(x)=0,x)`
  is `[]` where giac gives `[undef]`.

`tests/matrices_three.rs`, matrix constructors, slices, row operations,
spaces and norms (CAS-531, CAS-532, CAS-533, CAS-534):
- giac's `matrix`, `makemat`, `list2mat`, `mat2list`, `blockmatrix`,
  `JordanBlock`, `BlockDiagonal`, `rowdim`, `coldim`, `row`, `col`,
  `subMat`, `delrows`, `delcols`, `mRow`, `mRowAdd`, `scaleadd`, `scale`,
  `rowAdd`, `rowSwap` and `colSwap` (with aliases) read back and compared;
- the departures, each with giac's answer: `row` and `col` of a range a
  list (giac's is a sequence); `makemat(1,0,2)` and `list2mat(l,0)` stay,
  as a dimension is at least 1 (giac takes 0 for 1); `subMat` with
  reversed bounds stays (giac's `[[9]]`); a negative row in `mRow` counts
  from the end, as every index does (giac errs); and where giac errs,
  `rowDim` of a vector, `JordanBlock(2,0)`, `BlockDiagonal` of non-square
  blocks and `blockmatrix` of non-matrices stay;
- `maxnorm`, `rowNorm`, `colNorm`, `syst2mat`, and `rowspace`, `colspace`
  and `ibasis` where giac's answer is in reduced row echelon form, read
  back and compared; where giac's rows are fraction-free multiples, the
  same span (one reduced row echelon form of both answers);
- the departures of these: `rownorm` and `colnorm` of a vector answer
  (giac errs), `ibasis` with no intersection is `{}` (giac errs), and
  `rowspace(A,d)` stays (giac stores the dimension in `d`).

`tests/matrices_four.rs`, minors, echelon forms and pivots (CAS-530,
CAS-528), sums of squares and rational Jordan forms (CAS-535):
- giac's `det_minor`, `adjoint_matrix`, `ref`, `rref(A,keep_pivot)`,
  `pivot`, `egvl` and `eigVl` (giac's block order) and `gauss` read back
  and compared; `rat_jordan`, a sequence, read in brackets;
- the departures, each with giac's answer: `rat_jordan`'s `P` where
  giac's kernels differ from this CAS's (three) or where `A` is already a
  rational Jordan form (`P = I` here); `gauss` of a form that is not
  quadratic stays (giac returns it); giac's inert `Det`, `Inverse`,
  `Nullspace` and `Rref` stay user calls (giac prints `det(...)` etc.).

`tests/matrices_five.rs`, reshaping, triangles, quadratic forms and changes
of basis (CAS-639):
- giac's `REDIM` (`redim`), `border`, `semi_augment`, `simult`, `newMat`,
  `upper` and `lower` (matrices and strings, ASCII letters only, as
  giac's), `trn`, `a2q`, `q2a`, `changebase`, `charpoly(A,fadeev)` and
  `pcar(A,fadeev)` (the option ignored), `convert(x,array)` (`x`) read back
  and compared; `basis` (a set, as giac's) where giac's is in reduced row
  echelon form, and the same span where its rows are fraction-free
  multiples;
- the departures, each with giac's answer: where giac errs, the call stays
  with a warning (`REDIM` with a dimension below 1 or of a scalar,
  `border` with rows that do not fit, `simult` singular, of a vector or
  not square, `upper` and `lower` of a scalar, a vector or `[]`, `a2q`
  with sizes that do not fit, `q2a` of a form not quadratic in the
  variables, `changebase` singular or of other sizes); `newMat` with a
  negative or zero dimension stays (giac takes `|r|` and `|c|`:
  `newMat(2,-1)` is `[[0],[0]]`, `newMat(2,0)` `[[],[]]`), as `makemat`'s;
  `semi_augment(5,6)` stays (giac's `0`); `trn` of a scalar stays (giac's
  `conj(tran(2+i))`);
- recorded here only: `moyal` is not built, as giac 1.9's is a stub
  (moyal.cc:41-58, `symb_moyal`) answering its own `aide` example
  unevaluated (`moyal(x^2+y^4,x^4-y^2,[x,y],5)`); `Det` and `Inverse` are
  not aliases of `det` and `inv` but user calls held as giac's inert forms
  (`Det([[1,2],[3,4]])` prints `det([[1,2],[3,4]])`, which evaluates when
  read back), the departure rows in `tests/matrices_four.rs`.

`tests/integers.rs`, Xcas's integer functions (CAS-548 to CAS-553): giac's
`chrem` of integers and of polynomials with integer moduli (coefficient by
coefficient, modulo the product), `lgcd`, `fracmod` and `iratrecon`, `icomp`
with `zeros`, `bernoulli_mod`'s signed values, and `is_prime`'s 0, 1 and 2
and giac's own certificates `is_prime(n,1)` (CAS-551), read back and
compared; and the departures, each with giac's answer: its errors (no
solution, bad arguments, a modulus below 2, a denominator not prime to the
modulus, `(p-1) | k`, a non-integer `is_prime` argument), where ours leaves
the call with a warning; `chrem`'s Horner form `(-5*x+12)*x-5`, ours
expanded; `lgcd([-2*x-2])`, `2*x+2` here as `GCD`'s sign;
and, checked on our side only, `is_prime(10^30+57,1)` (giac's rho gives
up on `290240017 * 454197539` and its random ECM finds either factor
first) and the inputs that crash libgiac: `chrem` with polynomial moduli (ours folds through `chinrem`),
`bernoulli_mod(2,2)` (stays), `is_prime(-2147483659,1)` (0), and primes
whose `n - 1` ours does not factor far enough to prove them
(`is_prime(2^607-1)`, a 62-digit `2*96*p*q+1`): giac answers 1 through
PARI's APRCL where PARI does not fail, ours leaves the call with a warning.

`tests/list_commands.rs`, Xcas's list commands (CAS-555 to CAS-561): giac's
`find` and `index` from 0 (`1` found at `1.0`), `count_eq`, `count_inf` and
`count_sup` flattened and with `row` and `col`, `deltalist` of lists and of
matrices, `accumulate_head_tail`, `sizes`, `makevector`, `makesuite`
(bracketed to be read), `newList`, `enumerate` and `copy`, read back and
compared; and the departures, each with giac's answer: its errors (bad
arguments and types, fewer than two elements for `deltalist`), where ours
leaves the call with a warning (`newList` of a non-integer is ours-only:
giac reads stray bits before its type test and its error varies); the string
forms of `find` and `index` (phase 1.6), `find`'s start position and
`enumerate` of three arguments, which stay; and `newList(65537)`, past
giac's `LIST_SIZE_LIMIT` here, which stays with a warning. Then `apply`
with and without `matrix`, `foldl` and `foldr`, `multiply` (dot products
over the shorter vector, matrix products), `insert` and `subsop` (by
position, nested position and `NULL`, giac's quoted `'1=NULL'` written
`1=NULL`); their departures: giac's errors (fewer than three arguments to
`foldl`, positions out of range, dimensions that do not fit), the forms
that stay (several lists to `apply`, a decimal position, strings, a first
argument of `subsop` that is not a list), `apply` of `x^2` to the rows of a
matrix (giac's dot products), of `x+1` (giac's `poly1[…]`, ours an error)
and a nested `NULL`, which deletes where giac's unquoted one does not.
Then `heappush`, by `<` and by a comparator, `nextperm` and `prevperm`
(`undef` past the ends) and `groupermu` in giac's order; their departures:
`heapify`, which changes giac's variable in place and prints 1 where ours
returns the heap, `heappush` of a symbol, which giac orders and ours leaves,
and giac's errors (not a list, not a permutation, `heappop([])`). giac's
`heappop` of a non-empty list reads the element past the end (undefined:
0 in a fresh process, a crash in the oracle's shared one or on
`heappop([9,6,5,1,1,4,2,3])`): no giac row, unit tests of ours (the heap
after the pop) in `cas-core`.

`tests/expressions.rs`, Xcas's expression commands (CAS-562 to CAS-567):
giac's `pow`, `comDenom` (with a variable), `getNum` and `getDenom`, `simp2`,
`regroup`, `sq`, `logb`, `log2`, `log10`, `alog10`, `normalize`, `unitV`
and `algsubs` (five of giac's elimination cases, one nested) read back and
compared; and the departures, each with giac's answer:
`powexpand`, held in ours (its value released is giac's, and where the
canonical order is giac's it prints `hold(` giac's text `)`), and `radians`
and `degrees`, exact in ours, within giac's 12 digits at `x = 1`.

`tests/expressions_two.rs`, Xcas's expression commands, second part (CAS-568
to CAS-571): giac's `lname`, `indets`, `lvar`, `algvar`, `reorder`, `part`,
`eval`, `evalm`, `evalb`, `epsilon2zero`, `isinf`, `isnan`, `compare`,
`equal2diff`, `equal2list` (bracketed), `list2exp` and `exp2list` read back
and compared, kernels, identifiers and a sum's operands as sets where giac's
order is not our canonical one; and the departures, each with giac's answer:
`part` of a sum in our order, of index 0 (the head, giac's name string), of
an atom, past the end, of a list and of `f(x,y)` (giac's `of`);
`equal2diff(equal(x,2))`, as the language has no `equal`; and `compare`
between two heads, by name in ours, by operator address in giac. With them
`linstep`, `simplifyDirac` and `periodic` (CAS-570, CAS-572, CAS-573), sums
of steps and deltas read back in our canonical order (giac's coefficients
`series(c,x,z,1)` of a split delta read back as `c`); departures: the step of
a downward quadratic, minus the box in giac, and `periodic` with a period not
known positive, a giac error, unevaluated in ours; and in print only, giac's
`series(c,x,z,1)` coefficients of `simplifyDirac(Dirac(x^2-x-6))` and
`Dirac(2*x)`, recorded as giac prints them.

`tests/calculus_tools.rs`, Xcas's calculus helpers (CAS-574, CAS-576 to
CAS-580): giac's `ibpdv` and `ibpu` (one step, with bounds, and the finish
with 0), `preval` (with a variable, `x=a`, a removable singularity, and
`int` of one argument, in `x`) and `function_diff` and `fonction_derivee`
applied to a variable, a list and a number; `revert` of polynomials, of
functions to the default order and of `x^2+x^3` (a Puiseux inverse);
`sum_riemann` on giac's examples and a sum that is not one; and
`wz_certificate` in its argument forms, read back and compared; the calculus
of variations (CAS-575): `euler_lagrange` in its argument forms and
`jacobi_equation`, compared by value (giac's `diff(f,t,2)` read as a second
derivative, its constants `K_n`, counted across its calls, renumbered from
`K_0` as ours per call; an equation `= 0` up to its sign; the Jacobi solution
up to its constant), and `conjugate_equation`; and the departures, each with
giac's answer: `function_diff(sin+id)`,
as the language has no `id`; `function_diff(ln)(x,y)`, the arguments past
the parameter ignored in ours; `function_diff(sin)`, giac's hidden parameter
`` ` x` `` against our `x`; `revert(x*exp(x))`, giac's series a term longer;
`sum_riemann(1/sqrt(n^2+k^2),[n,k])`, where ours keeps `sqrt(n^2(1+x^2))`
whole; giac's fourth `wz_certificate` example, 0 in giac and a certificate in
ours, checked by giac's `simplify` of the WZ equation; `euler_lagrange` with
two coupled unknowns, where giac's d/dt of `L_v'` leaves out `u`; and in print
only, giac's `@` and `@@` forms.

`tests/euler_lagrange_assume.rs`, giac's `euler_lagrange` example under
`assume(y>=0)` (CAS-575), alone in its binary as giac's assumptions are global,
`assume` and `purge` in one input: ours, without the assumption, compared
with giac's at a point with `y > 0`.

`tests/function_study.rs`, Xcas's `domain`, `convex`, `tabvar` and `tabsign`
(CAS-581, CAS-583, CAS-582):
`domain` in mode 0 compared as a set of points (every end of either answer,
the midpoints, a point beyond each extreme) and by its shape, a list or one
condition; in mode 1 as a set of values; `convex` condition by condition,
each comparison's difference under `Together`; and the departures, each with
giac's answer: `asinh` (giac's `[]`), giac's assumption `0 <= x <= 2π` for a
periodic `f` other than one of `tan(c x)` alone, and variables `u(t)`.
`tabvar` and `tabsign` (CAS-582): the number of rows, each row past its label
(`tail(t[k])`), cell by cell (strings equal, values whose difference
`Simplify` makes 0, `exp` as powers of `e`), and the labels `x`, `y = f` and
`tabsign`'s `f`; `undef` where giac's is. Departures: the labels with primes
(`y'=(…)` is `"y'" = …`, `(y')'` the string `"(y')'"`) and giac's parametric
`tabvar([f, g], t)`, which stays.

`tests/kovacic.rs`, Xcas's `kovacicsols` (CAS-584): giac's examples and
equations with no Liouvillian solution read back, cases 1 and 2, `[]` and
giac's `0`, the equations given as `[a,b,c]` written without `/` (giac quotes
its arguments: `y''` does not read here, `diff(y(x),x,x)` does not there);
case 2's radicals by each solution's ratio to giac's, the same at two points;
and the departures, each with giac's answer: case 3 with a `y'` term (giac's
reversed polynomial), poles `±√2` (giac's `[]`), a cubic pole (ours stays),
a `/` in the list (giac's `0`), equations in `diff(y(x),x,x)` (giac's error)
and the errors, which stay here.

`tests/risch.rs`, Xcas's `risch` (CAS-585): the probe rows giac integrates
(logarithms, exponentials, `Ei`, `Si`, `Ci`, trigonometric functions through
exponentials, parameters, `atan`, residues in ℚ(i)) read back and compared
with ours, or by value where the forms differ (`Together` of the difference,
else the same change from `x = 7/2` to `x = 9/2`: giac's `ln(i x − 1)` is
ours `ln(x + i)`, a constant apart); its partial answers by their parts, the
primitive and the remainder, giac's `integrate(rem,x)` read as an inert call
(our `Integrate` would evaluate it). Departures, each with giac's answer:
what giac leaves prints `integrate(f,x)`, which ours reads back and
evaluates (ours stays `risch(f,x)`); Rothstein–Trager factors of degree 2
split by a square root (`1/(x^2-2)`, `1/(x^3+1)`, `1/(x^4+1)`,
`x^3/(x^2+x+1)`), not ported; giac's final `re(halftan(…))` step
(`exp(x)/(x^2+1)` in `tan(1/2)`, `-pi+atan(exp(x))`), which ours skips;
tower-size ties (two seeded rows of T9.c1), where giac keeps its printout's
order and ours its canonical one, so the split differs (ours passes the
derivative check); giac's wrong answers (two seeded rows whose primitive
fails differentiation, risch.cc:693), where ours stays; and `Int`, giac's
`integrate` unevaluated, here the user's call (R1.2-001).

`tests/sampling.rs`, Xcas's `linspace`, `linear_interpolate`,
`interval2center` and `center2interval` (CAS-593): giac's answers read back
and compared element by element within `10⁻¹⁰` relatively, since ours are
exact where giac answers doubles (`19/3` for `6.33333333333`), `center2interval`
bracketed (C0023). Departures, each with giac's error as caseval prints it:
`linspace(1,2,1)` (Invalid dimension), a range past the points of
`linear_interpolate` (Values out of range) and one centre for
`center2interval` (Bad Argument Value), where ours stays.

`tests/interpolation.rs`, besides `lagrange`, `spline`, `thiele` and the
orthogonal polynomials, Xcas's `ratinterp` and `triginterp` (CAS-594):
compared by `Expand` and `Together` of the difference, decimal `triginterp`
by value at three points within `10⁻¹⁰`. Departures, each with giac's answer
as caseval prints it: `d` out of `0 … |X|` and one point for `ratinterp`,
one value for `triginterp` (giac's errors), and `d = |X|`, a repeated
abscissa and `a = b` (giac's `undef`), where ours stays.

`tests/derivatives.rs`, Xcas's `numdiff` and `nDeriv` (CAS-595): exact
`numdiff` by `Expand` and `Together` of the difference; giac's aide example
of decimals within `4·10⁻⁷` of the closed form `f″(1/π)` and `10⁻⁹` of giac;
`nDeriv` by value at three points within `10⁻¹⁰` (ours exact, `h = 1/1000`),
and `nDeriv(f(x),x,h)` symbolically, giac's `0.5` read as `1/2`. Departures,
each with giac's answer as caseval prints it: an order not a nonnegative
integer, unequal lists and two arguments (giac's errors), a repeated
abscissa (`undef`), and `nDeriv` of a variable not a symbol (`0.0`), where
ours stays.

`tests/fits.rs`, Xcas's `bspline` and `fitspline` (CAS-596): `bspline` by
`Expand` and `Together` of the difference element by element (C0032), with
`piecewise` branch by branch, conditions by their truth at points (giac's
`1>x` is our `x<1`); `fitspline` at points, ours exact equal to the closed
form and within `10⁻¹⁰` of giac's doubles; with a variable, each piece by
value at three points of its interval within `10⁻⁹`. Departures, each with
giac's answer as caseval prints it: points outside the data (`[7]` and
`[-1]`, giac's `[0.936507936508]`, its last polynomial continued), where
ours stays; `derive=1` (giac's rows of one,
`[[0.600529100529],[-0.28835978836]]`), where ours is flat,
`[227/378,-109/378]`; a variable with `breakpoint=` (giac builds on equal
intervals), where ours builds on the breakpoints and equals the points form
at a point of each piece; and too few knots (`Expected a variable`), where
ours stays. `fitpoly` (CAS-597) coefficient by coefficient within
`10⁻⁹ max(1, max |c|)` of giac's, a coefficient giac leaves out taken as 0:
functions (`exp`, `x*cos(x)`, `sin`, `x^2`, Runge's `1/(1+25*x^2)`, where
ours keeps degree 15 with odd coefficients holding 0 and giac prints 14;
`degree=` and `limit=`), and data with a variable, also against the closed
form (exact data exactly, decimals within the tolerance). Departures, each
with giac's answer as caseval prints it: no variable (giac's internal
`poly1[1,0,1.0]`, and for the aide's data
`poly1[-8.20707070707,-0.800116550117,8.89355089355,0.244755244755]`),
where ours is the coefficient list, highest first; `cos(x)` on `0..2` with
`1e-8` (giac's degree 15, its `gaussquad` of the weighted norm being off),
where ours has degree 5 near the prototype's; `abs(x)` on `-1..1`, giac's
coefficients off ours by up to `1.3·10⁻⁶` (its doubles in the degree-15
expansion), where ours is near the prototype's.

`tests/minimax.rs`, Xcas's `minimax` (CAS-599): giac's converged rows
(`exp(x)` on `0..1` degree 3, `sin(x)` on `0..1` `[2,1]`,
`abs(x)*sqrt(abs(x))` on `-2..2` degree 15, `abs(x)` on `-1..1` degree 3,
`x^3` on `-1..1` degree 2), ours' error within 2% of giac's and our `p`
within twice giac's error of giac's at 51 points (polynomials of balls by
value, never by `CoefficientList`); closed forms, `E` within `10⁻⁶ E` and
`p` within `10⁻⁸` at 51 points: `x^3` degree 2 (`3/4*x`, `1/4`), `abs(x)`
degrees 2 and 3 (`x^2+1/8`, `1/8`), `exp(x)` degree 1 (`(e-1)*x+1-E`,
`E = (1-(e-1)*(1-ln(e-1)))/2`). Departures, each with giac's answer as
caseval prints it: `abs(x)` degree 2 (giac's `0.213388376846`, its exchange
losing the alternation), and `exp(x)` degree 1 (`0.106657010751`, its stop
at `max |e| ≤ 1.02 min |e|`), where ours is the closed form;
`min(1/cosh(3*sin(x)),sin(9*x/10))` on `-3..4` degree 30 (`0.0946608433887`,
its exchange stalling), where ours is `0.0319105`, `|f − p| ≤ E (1 + 10⁻³)`
at 51 points. `tests/minimax_rand.rs`, its own file as it seeds giac:
`x*exp(-x)` on `0..10` `[5,3]` with `rand=30`, giac seeded with
`srand(42)` in the same call (`0.000287329610301`; unseeded `0.000256`),
where ours, seeded the same way, is `4.35626e-5`, `|f − p/q| ≤ E` at 51
points.

`tests/optimization.rs`, besides `extrema`, `minimize`, `maximize`,
`lpsolve` and `nlpsolve`:
- `tpsolve` (CAS-598): giac's 23 probe rows (zero supplies, balanced, unbalanced with a
  dummy row or column, forbidden routes `M`, degenerate plans), each
  bracketed, cost and plan equal; departures, where ours stays
  `TPSolve[…]`: the two rows giac finds no feasible solution for, a
  decimal cost, a negative supply and `2*M`.
- `lp_options` (CAS-600, CAS-601): giac's 26 probe rows of `lpsolve`'s options (`lp_method`,
  `lp_presolve`, `lp_verbose`, the forms of maximize, `assume`, `lp_assume`,
  `lp_integervariables`, `lp_binaryvariables`, a bare `lp_integer`, decimals); the optimum
  equal, our point checked by its constraints, integers and objective, not giac's point
  (`x+y` on `x+y>=1` has several).
- `lp_branch_and_bound` (CAS-602, CAS-603): giac's 23 probe rows of the branch and bound
  options (`lp_nodeselect`, `lp_varselect`, `lp_maxcuts`, `lp_heuristic`, `lp_gaptolerance`,
  `lp_depthlimit`, `lp_nodelimit`, `lp_iterationlimit`), `[]` where a limit is reached first;
  checked as `lp_options`. Departures: `lp_timelimit=1` (giac's 1 ms CPU clock gives `[]`,
  ours has no clock and answers the optimum) and `lp_nodeselect=foo` (giac's "Bad Argument
  Type", ours stays `LPSolve[…]`).
- `simplex_reduce` (CAS-606): giac's 8 probe rows, both forms (`A, b, c` and a tableau), an
  unbounded program (`+infinity,[]`), each bracketed, exact rows by `Together` of the
  difference, the decimal row numerically. Departure: a negative right side, from which giac
  pivots an infeasible basis (`-9/2,[0,0,-3/2,0,11/2],…`); ours warns and stays
  `SimplexReduce[…]`.

`tests/find_minimum.rs`, Xcas's `find_minimum` (CAS-608): giac's 7 probe rows that compile
(`cosh`, with `eps=1e-5, maxiter=20` too, `x^2-2*x`, `(x-1/3)^2`, `x->x^4-3*x^2+x`, `sin`,
`abs(x-1/2)`), closed form first: our ball holds it (at 60 digits) with radius at most `eps`,
giac's double within `10⁻⁶` of it (`10⁻⁴` for `eps=1e-5`). Departures: the `Airy_Ai` and
`besselJ`/`Gamma` rows (giac `-3.96344565209`, `4.82843934915`), which do not compile to ball
operations; ours stays `FindMinimum[…]`.

`tests/frank_wolfe.rs`, Xcas's `frank_wolfe` (CAS-609): giac's 3 probe rows (`(x1-2/3)^2/3+(x2-2/3)^2/5`
on the ± unit points, `(x1-1)^2+(x2-1)^2` on `[[0,0],[1,0],[0,1]]`, `x1^2+x2^2` on `[[1,1],[2,1],[1,2]]`),
by default and with `exact`, closed forms first: our point in the hull with `f − f* ≤ 10⁻³` (`f*` 1/72, 1/2,
2), giac's point within `10⁻³` by value too; `rand=5` on ours only (random starts on both sides). Departure:
`maxiter=3` on the first two rows, where giac answers `["Unable to minimize at given precision, last value",
[…]]`; ours warns and stays `FrankWolfe[…]`.

`tests/simulated_annealing.rs`, Xcas's `simulated_annealing` (CAS-610): the row
`simulated_annealing(2,x->x^2+sin(20x)/5,(x,y)->abs(x-y),(x,m)->rand(x-m,x+m),[1,0.1,1.002,1e-6],[5,20],1)`
against its closed form first, `x* = −0.0766237535726` (`f* = −0.193981966498`): giac unseeded and seeded with
`srand(42);…` and `srand(7);…` in the same input (`-0.0766196874902`, `-0.0766106601804`, `-0.0766229843721`),
ours seeded with 1, 42 and 7, each within `10⁻³` of `x*`; ours never compared with giac's. Departures: our draws
(ours and the user's `rand`) come from the session's generator, not GSL's default generator nor giac's `rand`, so
`SeedRandom` makes a run repeat; giac's errors in `cost` or `step` leave the call; our schedule checks
(`k > 0`, `mu_t > 1`, `0 < t_min < t_initial`, the steps within `terms`) warn and leave the call, where GSL
loops forever when `mu_t ≤ 1`.

`tests/nlpsolve.rs` and `tests/nlpsolve_seeded.rs` (CAS-604, CAS-605), `nlpsolve`'s options matched by
their names (`nlp_initialpoint` points and regions, `nlp_iterationlimit`, `nlp_precision`,
`nlp_tolerance`, `nlp_method`, `nlp_maximize`, `nlp_presolve`, `nlp_verbose`, `assume`,
`nlp_integervariables`, `nlp_binaryvariables`):
- `linear_programs`: 11 linear programs, which giac and ours hand to `lpsolve` and `LPSolve`,
  exact and equal (unique optima; integers, binaries, decimals, the point in sorted order).
- `from_an_initial_point` and the seeded `options` (14 rows, `srand(42)` in the same call, C0014):
  ours checked against the closed form to 10⁻¹⁰ first, giac's within its observed error (10⁻⁵,
  10⁻³ with `nlp_tolerance=1e-3`, 10⁻⁴ for the five-variable product).
- Departures, giac's text copied: an infeasible linear program (giac "The problem is
  infeasible", ours `[]`), an unbounded one (giac `+infinity`, ours `[+infinity, []]`),
  `nlp_nonnegative=true` and `nlp_verbose`, `nlp_presolve` alone (giac "Invalid argument", ours
  answers), `nlp_method=cobyla` and `nlp_method="foo"` (both refuse; ours stays `NLPSolve[…]`);
  seeded, where giac's "Failed to find a solution" and ours answers the closed form: `(x-1)^2+(y-2)^2`
  under `x+y<=2` (1/2 at (1/2, 3/2)), `x^2+y^2` under `x+y>=1`, `(x+1)^2+(y-2)^2` with
  `assume=nlp_nonnegative` (1 at (0, 2)), `nlp_method="ipt"`; and `nlp_iterationlimit=2`, giac's
  "Failed to optimize at given precision", ours stays with a warning.
- Integer and binary variables in nonlinear programs (CAS-605): `integer_programs`, giac's one
  answered row (`(x-1/3)^2`, `x=-5..5`, `assume=integer`: ours exactly `[1/9, [x=0]]`, giac's
  within 10⁻⁸); the seeded `integer_departures` (9 rows: the aide's program under `assume=integer`,
  `assume=nlp_integer`, `nlp_integervariables=[x1,x2]` and `[x1]`, the binary cubic under
  `assume=nlp_binary`, with `maximize`, and `nlp_binaryvariables`, `(x-3)^2+(y+1)^2` under
  `assume=nlp_nonnegint` and on a box with `assume=integer`), giac's "Failed to find a solution"
  copied, ours the known optimum (exact at integer points, `x2 = sqrt(5)` certified).

`tests/strings.rs`, Xcas's string commands (CAS-611, CAS-613):
- giac's answers to `char`/`chr`, `asc`, `ord`, `inString`, `contains` (with giac's interval
  `[lo..hi]`), `mid`, `left`/`lhs`, `right`/`rhs` and `levenshtein` on the probe rows, read back
  and compared with ours.
- Departures, giac's text copied: giac counts bytes of UTF-8 where ours counts Unicode scalar
  values (`asc("é")` giac `[195,169]`, ours `[233]`; `ord("é")` `-61` and `233`; `char(233)` an
  invalid byte and `"é"`); giac's errors `char(65.5)`, `ord("")` and `levenshtein("abc",[1,2])`
  stay here (the last with giac's message as a warning); giac has no `hamming` command (it
  stays), ours counts the differing positions of equal-length strings or lists (`3`, `2`).
  Not compared: non-integer counts (giac's `left("abc",2.0)` is `"ab"`, ours stays) and giac's
  other one-argument `left`/`right` forms (integer bits, lists, first operands), which stay.
- CAS-612: giac's answers to `split` (strings, and `split(e,[x,y])` of an expression through
  `factors`), `join`, `trim` (strings, lists with giac's epsilon or a threshold, `left`/`right`)
  and `strip`, and strings in `size`/`length`, `concat`, `"abc"[k]` and slices, `==`, `<`, `>`,
  `>=` and `sort` (byte order). Departures, giac's text copied: `length("héllo")` giac 6 bytes,
  ours 5 characters; `split("abc","")` ("Invalid dimension") and a `trim` threshold not a
  positive real ("Invalid threshold specification") stay; `join(",",[1,2])` (`"1,2"`) and
  `concat("a",1)` (`"a1"`) agree since CAS-614 prints non-strings; `"ab"+"cd"` (giac `"abcd"`) stays;
  `"abc"[5]` and `"abc"[-4]` are errors in both, and `"abc"[0..5]`, giac's `"abc"`, is out of
  range here, as for a list.

`tests/text.rs`, expressions as text and back (CAS-614):
- `string`/`str` and `cat`, and non-strings in `join` and `concat`: giac's text (its quotes
  removed, inner quotes undoubled) and ours are equal without spaces or read back to the same
  value; `expr`/`execute` (strings, `expr("")` as `undef`, `expr(sin,x)`, `expr(f,x,y)` staying,
  `expr(string(x^3-1))`) compared by value.
- Departures, giac's text copied: `expr("0123")` giac's octal 83, ours 123; `expr("a:=3")` (giac
  assigns, 3), `expr("1;2")` (giac `1,2`) and `expr("x+")` (giac's parse error) stay here with
  the reader's message as a warning; ours prints the input form, number first and sums spaced:
  `size(string(x+1))` giac 3, ours 5; `cat("a",x+1)` and `concat("a",x+1)` giac `"ax+1"`, ours
  `"a1 + x"`; `join(",",[1,x+1,"z"])` giac `"1,x+1,z"`, ours `"1,1 + x,z"`.

`tests/latex.rs`, `latex` and `TeX` (CAS-615): 117 rows where giac's text (quotes removed,
whitespace collapsed) is ours: the plan's probes, `cas_lang::markup`'s lowering rows that print as
giac prints its own answer, `rootof([1,0],[1,0,1,1])` (held in both), and 57 more rows of
`crates/cas-lang/src/markup/giac_latex.tsv`.
- Departures (16), ours asserted, giac's text copied: all come from giac's evaluation order, not
  the printer (the same shape prints the same). Constants first in sums: `sqrt(x+1)` giac
  `\sqrt{x+1}`, `(x^4-1)/(x^2+3)`, `2*x^2-3*x+1` (ascending degree), `1/(x+1)`, `(x+1)^2`,
  `(x+1)*(x-1)`, `-(x+1)` (ours keeps `-(1+x)`, giac `-x-1`); our product order: `sin(x)*cos(x)` (giac
  `\sin x\cdot \cos x`), `x*sqrt(2)` (giac `x \sqrt{2}`); giac rewrites what ours keeps: `x<y`
  (giac `y>x`), `x<=y` (giac `y\geq x`), `x!=y` (giac `\mathrm{true}`), `ceil(x)` (giac
  `-\left\lfloor -x\right\rfloor`), `binomial(n,k)` (giac's factorial quotient), `x^(2/3)` (giac
  `\left(x^{\frac{1}{3}}\right)^{2}`), `sum(k,k,1,n)` (giac
  `\frac{-n-1+\left(n+1\right)^{2}}{2}`, ours `\frac{1}{2} n \left(1+n\right)`, as giac prints `n*(n+1)/2`).
- Algebraic extensions and points (tex.cc's `gen2tex`, giac's `latex` on geometry) have no
  counterpart in ours.

`tests/mathml.rs`, `export_mathml` (CAS-616): `latex.rs`'s 115 inputs (not `TeX(x^2)` nor
`diff(f(x),x,x)`, which giac's MathML table lacks) in giac's three forms, `export_mathml(x,display)`,
`export_mathml(x,content)` and `export_mathml(x)`: giac's text (quotes removed, doubled quotes
undone, whitespace collapsed, C0047) is ours. In the default form the annotation
(`<annotation encoding='giac'>`) is each one's input form, so it is left out of the comparison and
ours is checked to be our `string(x)`, escaped (giac `x+1`, ours `1 + x`).
- Departures (16), ours asserted in presentation and content: `latex.rs`'s, for the same reasons
  (giac's evaluation order and rewrites; the printers agree on the same shape); giac's MathML of
  each is in `crates/cas-lang/src/markup/giac_mathml.tsv`. `x<y`: giac
  `<mi>y</mi><mo>&gt;</mo><mi>x</mi>`, ours `<mi>x</mi><mo>&lt;</mo><mi>y</mi>`.
- giac's `mathml()` (mathml.cc, a separate printer, `cas_lang::mathml`): the same 115 inputs,
  `mathml(x)`: giac's text (quotes removed, doubled quotes undone, a newline written ' | ' on both
  sides, C0047) is ours exactly. The 16 departures, ours asserted (the formula of the prototype on
  our shapes, C0068, in giac's XHTML); giac's texts are in
  `crates/cas-lang/src/markup/giac_mathml_doc.tsv`.

`src/format.rs` tests the rewriting on outputs giac gave. Nine sabotages are
caught.

`scripts/giac_gap.sh <outdir>` lists giac's commands (`aide_cas`) with and without a name in `cas-lang`.

Gap ledger (phase 1.8 T2, 2026-11-26; giac 1.9.0's `aide_cas`): every giac entry with no name in
`cas-lang` is accounted for. Regenerate the list with `bash scripts/giac_gap.sh <dir>` and
`python3 scripts/gap_describe.py missing <dir> > <dir>/missing_described.txt`. Then
`python3 scripts/gap_ledger.py <dir>/missing_described.txt` re-tallies the buckets. It fails when an entry
has no bucket, two buckets, or a named entry is gone. There are 1738 entries: 734 covered and 1004 missing,
bucketed below. Build tasks: T3 Algebra, T4 Expressions/Calculus/Integers, T5
Language/Polynomials/Statistics/Infrastructure, T6 Signal processing.
- Out of scope, 913 entries. New module, 565 (units 216, geometry 128, graph theory 221): these
  modules don't exist yet (gen 2). Plotting, 157 (turtle, graphics, plot display variants): we have no
  graphics. Syntax, 103 (operators, keywords, French control words): another syntax, not part of the language
  (CAS-DOC-7). I/O, 82 (files, TI screens and folders, GUI, debugger, sound and images): no I/O in a
  library engine. ML/music/geodesy, 6 (`neural_network`, `train`, its options, `tonnetz`,
  `geodetic2ecef`, `ecef2geodetic`): not a planned area.
- Already a departure above, 12: `QR`, `LU`, `confrac`, `Factor`, `Gcd`, `Resultant`, `Quo`, `Rem`
  (CAS-527), `xcas_mode`/`maple_mode`, `python_compat`, `gbasis_reinject`, `rur_gbasis`.
- Alias, inert or option spellings, 29, each with the name it maps to:
  - Already working: `besselY`, `besselJ` (`BesselY`, `BesselJ`); `parfrac`, `fullparfrac` and
    `polynom` (`convert` options); `base` (`convert(p,base,b)`); `mod %` (`%`, Modular); `!` postfix
    (`factorial`).
  - Inert, held as an unknown call: `Int` (`integrate`), `Rref` (`rref`), `Nullspace` (`nullspace`),
    `Fourier` (`fourier`), `Hilbert` (`hilbert(f,t)`). `Det` (`det`) and `Inverse` (`inv`) are held
    user calls here too, giac's printing `det(...)`: departure rows in `tests/matrices_four.rs` (T3).
  - Operator syntax, not built: `'` (`diff`; the lexer rejects it), `|` (`subst`; the lexer rejects
    it), `/%` and `%/` twice (`iquo`, `irem`), `div` (`iquo`), `!` prefix (`not`).
  - Not handled yet, assigned: `rdiv` (`/`), `getType` (`type`), `at` (`At`) and `subs` (`subst`) to
    T5; `maple_ifactors` (`ifactors`) to T4; `list` (`convert`) to T5; `expln` (`convert`, as trig2exp)
    to T4; `array` (`convert`, `x` itself) and `fadeev` (`charpoly`, the option ignored) built in
    T3 (`tests/matrices_five.rs`).
- Computational in an existing module, 50: the 1.8 plan said 49, and `box_constraints` moved in from
  ML. Each is marked with its task:
  - Algebra (15): `REDIM`, `border`, `semi_augment`, `simult`, `newMat`, `upper`, `lower` (matrix
    and string forms), `a2q`, `q2a`, `basis`, `changebase`, `trn` built in T3
    (`tests/matrices_five.rs`); `moyal` not built, a departure (giac's is a stub answering
    unevaluated); `Fill` built in T3 (`tests/session.rs`). `GF` is not
    built: giac's Galois-field elements need a new value type (gen 2).
  - Calculus (4): `envelope`, `avgRC`, `addtable` T4. `taux_accroissement` is not built: it is a
    French spelling (CAS-DOC-7 §5).
  - Expressions (11): `shift_phase`, `linabs`, `non_recursive_normal`, `root`, `additionally`,
    `atan2`, `has`, `lin`, `fmod`, `modf`, `period` T4.
  - Integers (2): `multinomial` T4; `order` T4 for the `(Z/nZ)^*` form (`order(3 % 7)` = 6). Its
    finite-field form is not built: it needs `GF`.
  - Numerical analysis (1): `box_constraints` T4 (`box_constraints([x,y],[[-2,2],[-5,5]])` =
    `x=-2..2,y=-5..5`).
  - Language (5): `back`, `range`, `shift`, `format`, `suppress` T5.
  - Polynomials (3): `roots`, `truncate`, `divis` T5.
  - Statistics (2): `choice` T5. `hasard` is not built: it is a French spelling, and `rand` is built.
  - Infrastructure (1): `dayofweek` T5.
  - Signal processing (6): `stft`, `istft`, `hht`, `boxcar`, `tri`, `sinc` T6.
- Option-level items under covered names: `kmeans`/`cluster` `output=part`, `index=` and clustering
  of strings T5; `hilbert(f,t)`, the transform form beside the list form, T6. `fitdistr` on
  exponential data with a 0 already agrees: giac refuses `fitdistr([0,1,2,3,4],exponential)` with
  "Expected a positive number" (caseval probe), and so does ours, because it requires data > 0
  (`fitting.rs:645,666-668`). T5 adds the oracle row.
