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

`src/format.rs` tests the rewriting on outputs giac gave. Nine sabotages are
caught.

`scripts/giac_gap.sh <outdir>` lists giac's commands (`aide_cas`) with and without a name in `cas-lang`.
