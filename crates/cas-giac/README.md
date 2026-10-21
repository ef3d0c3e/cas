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
`list[…]` is `[…]`, and `set[…]` is `{…}`.

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

`src/format.rs` tests the rewriting on outputs giac gave. Nine sabotages are
caught.
