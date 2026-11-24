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
for `[x^2-1,y^2-1]`, and the calls of giac's F4 knobs
(`gbasis_max_pairs`, `gbasis_simult_primes`, `gbasis_reinject`,
`rur_gbasis`, `rur_certify`) staying.

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
- the departures: giac's `increasing_power` changes printing and is left for
  phase 1.6 (CAS-546); a setting's bare name reads its value here, 0 or 1,
  where giac gives the symbol.

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

`src/format.rs` tests the rewriting on outputs giac gave. Nine sabotages are
caught.

`scripts/giac_gap.sh <outdir>` lists giac's commands (`aide_cas`) with and without a name in `cas-lang`.
