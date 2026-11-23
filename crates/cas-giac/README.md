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
between two heads, by name in ours, by operator address in giac.

`src/format.rs` tests the rewriting on outputs giac gave. Nine sabotages are
caught.

`scripts/giac_gap.sh <outdir>` lists giac's commands (`aide_cas`) with and without a name in `cas-lang`.
