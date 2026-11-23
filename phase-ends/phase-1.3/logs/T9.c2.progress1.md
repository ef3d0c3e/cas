# TASK_PROGRESS — T9.c2 (CAS-585 risch native), coder handoff 1

## Done so far
- Prototype fixes tested (.run/1.3/T9/risch.py, env FIX693/FIX837, defaults now 0/1): :693 REJECTED (breaks probe
  `ln(ln(x))/x`, batch 47 loses its partial answer, 19/30/89/116 become leave-all); :837 KEPT (probe unchanged 74+4+9,
  0 diff failures; `1/(x^3+x+1)^2` d=FAIL → ok). Evidence .run/1.3/T9c2/probe_fix{00,01,11}.txt, rows_fix*.txt.
- cas-poly primitives committed 28162ad: `SparseRing::valuation`, `pseudo_div_rem` (sparse.rs), `extended_gcd` →
  `Bezout`, `abcuv` → `Abcuv`, `partial_fractions` → `PartialFractions` (gcd.rs, exported in lib.rs); tests in
  sparse/tests.rs, gcd/tests.rs; doc paragraphs sparse.md, gcd.md. Local gate green.
- Log phase-ends/current/logs/T9.c2.md (timeline so far).

## In flight (uncommitted, saved in TASK_PROGRESS.patch; files also on disk, untracked)
- crates/cas-core/src/native/integrate/risch.rs: module doc, `risch` native (args, Stop/Step, `answer` with fixpoint
  ROUNDS=3 on the remainder, `verified` D-check), `risch_top` (:951), `risch_lin` (:878), `constant_part`, `rlvarx`
  (Polynomials generators recursively, sort by size then Expr Ord), `in_risch_tower`, `risch_tower`, `tsimplify`
  (ArcTrigToLog + TranscendentalSimplify), `pow2expln_x`, `lin`/`compress`/`convolution`/`unlin`, helpers
  (`derive` rewrites E^u → Exp[u], `has_i`, `kernel`, `within` bound).
- risch/tower.rs: `Tower` (Polynomials over [exprs, v, D v], ring GaussianIntegers with nvars = gens + 1, the last
  var `t` for Rothstein–Trager), `Frac`, K[T] ops (quo_rem by pseudo-division, tgcd primitive in T, abcuv, sqff with
  zero_apart, partfrac via cas-poly, factors Trager/Wang, linear_root, at_root, lex_lead).
- risch/parts.rs: `Lnc`, `Acc` (needs an `Acc::empty()` constructor: rde.rs calls `Acc::empty()`, parts.rs has
  private `Acc::new()` — rename/make pub(super)), `in_risch` (:743, fix :837 applied), `risch_poly_part` (:662),
  `rational_in_x`, `is_linear_wrt`, `hermite_reduce` (:517), `rothstein_trager_resultant` (:124, must be pub(super):
  rde.rs imports it), `giac_gcd_norm`, `without_monomial`, `integrate_ei` (:582), `integrate_poly_exp`.
- risch/rde.rs: `rischde_simplify`, `spde_x`, `spde`, `risch_desolve` (complete port, never compiled).
- NOT YET: never compiled. Not wired: no `mod risch;` in integrate.rs, no `pub(super) use`, no builtin, no NATIVES
  row, no en.txt line, no tests.rs.

## Hypotheses rejected (evidence)
- fix :693 literally (allowed_lnarg) and with lncoeff propagation: both fail probe `ln(ln(x))/x` (prim 0): the inner
  in_risch(1/x,[x],allowed=ln x) refuses ln x as it is not ln(ln x).
- Representation over Z only: rejected; GaussianIntegers always (gcd/division extension-invariant; factoring follows
  giac's complex flag via `Tower::factors(complex)`).

## Current hypothesis
- The port compiles after routine fixes; the repro crate (.run/1.3/T9c2/repro, `cargo build --release` with
  CARGO_TARGET_DIR=target/repro) evaluates giac-syntax lines on stdin (and `!gens <FullForm> ;; <FullForm>` prints
  Polynomials generators) for tracing rows before writing tests.

## Next 5 steps
1. Wire: integrate.rs `mod risch;` + `pub(super) use risch::risch;` (template bab9f39); builtin.rs
   `RISCH = "Risch" [PROTECTED];` near INTEGRATE; native.rs NATIVES `(Symbol::RISCH, integrate::risch)`; en.txt
   `Risch: risch` in [builtins] (after KovacicSols, line ~118). Create risch/tests.rs (can start empty `use`-free).
2. `cargo check -p cas-core` and fix compile errors (Acc::empty, pub(super) rothstein_trager_resultant, unused
   imports, `Complex / &Complex`, `ring.content` normal form, clippy pedantic casts in `within`).
3. Rebuild repro and run the 13 required probe rows + the 74 others (.run/1.3/T9/probe.out answers) through
   `risch(...)`; debug until most match; note unreachable rows (sqrt-split quadratics 1/(x^2-2), 1/(x^3+1),
   1/(x^4+1), x^3/(x^2+x+1); exp(x)/(x^2+1) re/im form).
4. Write risch/tests.rs: each reached row checks Equal with giac's answer (Together/TranscendentalSimplify of the
   difference; constant-difference rows only by D) AND D(answer) == f; Ok(None) rows: exp(x^2), 1/ln(x), sqrt(x),
   x^2*exp(-x^2), (x^2+1)/(x^3+x+1); `Risch[rem]` of a partial answer stays.
5. Full gate (fmt, check, stable clippy, fuzz check), commit `T9.c2: CAS-585: risch native …` with all risch files
   + builtin.rs, native.rs, integrate.rs, en.txt; finish the log; delete TASK_PROGRESS.patch once committed.

## Gotchas
- `D[Ei[u]]` gives `E^u`, not `Exp[u]` (C0036): `derive` and `verified` rewrite E^u → Exp[u].
- Polynomials turns `Exp[-x]` into generator `Exp[x]` with exponent −1, and `Exp[x/2]`,`Exp[x]` into `Exp[x/2]`:
  kernels must come from Polynomials generators (rlvarx does), else `Tower::new` returns None.
- `together(...)` is not a giac name; giac syntax `normal`/`ratnormal` → Together. `cargo fmt --all -- file`
  formats the whole workspace.
- integrate_ei keeps risch.cc's original-fraction flaw (:640 in the prototype line 1037); the D-check catches it.
- Files over 20k chars (sparse.rs, gcd.rs) can't be read/sed whole (hook): use Read with offset or python edits.

## State to carry verbatim
- Commit 28162ad (cas-poly primitives). Repro binary: target/repro/release/repro. Scratch .run/1.3/T9c2/.
- Prototype defaults: FIX693 '0', FIX837 '1' (.run/1.3/T9/risch.py:26-27).
