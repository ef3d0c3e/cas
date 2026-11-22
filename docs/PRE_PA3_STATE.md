# Pre-PA3 state

As of 2026-11-22 (PA3 installed; last product commit 2026-11-05, `866047b`, pushed: `origin/master` = `master`).

## What the project is
cas: a computer algebra system engine in Rust; the programming API is first class, syntaxes are frontends (`README.md`).
Crates: `cas-num` (integers, rationals, modular, ball arithmetic), `cas-poly` (dense/sparse polynomials, gcd, factor,
Gröbner), `cas-core` (hash-consed expressions, evaluation loop, builtins), `cas-lang` (lexer, parser, names, printer),
`cas-giac` (libgiac loaded at run time as test oracle, test support only), `cas-tune` (threshold tuning).
Xcas/giac is the reference for defaults; most builtins are named after giac commands (`crates/cas-giac/README.md`).
576 commits, one per Lific ticket: subject `CAS-N: …`.

## How work was organised before PA3
- Tracker: Lific, remote only (lific.pundalik.org, project `CAS`; CLI `--backend http`). Design docs are Lific pages
  `CAS-DOC-n`; each area is an umbrella issue plus a plan `CAS-PLAN-n` whose steps are tickets.
- Areas done, newest last: … Statistics II (CAS-490, DOC-50), Matrix decompositions III (CAS-495, DOC-51),
  hypothesis tests (DOC-52, PLAN-46, CAS-502..509), Piecewise functions (DOC-53, PLAN-47, CAS-513..515).
  No `CAS-PLAN-48` or `CAS-DOC-54` exists.
- Lific modules: Algebra, Balls, Calculus, Evaluator, Expressions, Infrastructure, Integers, Language, Polynomials,
  Statistics.

## In flight
- A run of standalone giac-parity tickets (no plan), all committed on master 2026-11-05 but not closed in Lific:
  CAS-517 FSolve method options, CAS-518 FSolve of systems (Newton + Krawczyk), CAS-519 NLPSolve (Lific: backlog),
  CAS-520 BVPSolve, CAS-521 window functions, CAS-522 lowpass/highpass, CAS-523 cross/auto correlation,
  CAS-524 threshold (Lific: active). CAS-516 (Simplify shared subexpressions, a fuzz timeout fix) is done.
- **Master's CI is red** at `866047b` (runs 42057178857 CI, 42057178575 Fuzz):
  - CI job `oracle`: `boundary_values` (`crates/cas-giac/tests/ndsolve.rs:91`) fails: our
    `bvpsolve((32+2x^3-y*diff(y(x),x))/8,[x=1..3,y],[17,43/3],4)`: interior values differ (ours y(2) = 12, giac's
    ≈ 12.02; y(3/2) 12.92 vs 12.96); endpoints agree, ours exact, giac's balls — from CAS-520.
  - Fuzz target `eval`: deadly signal, `attempt to divide by zero` at `crates/cas-num/src/rational/mod.rs:52:41`.
  - fmt, clippy, test and fuzz-check pass. CI was also red at CAS-522 and CAS-523.
- Before closing CAS-517..524: green CI on master, then status `done` in Lific.

## What is next
- Fix the two red jobs above (rule: fix what CI or the fuzzer finds right away; a fuzz finding gets a regression test
  in the owning crate).
- Then the next area, in Claude's recommended order, without asking the developer which; a brand new Lific module
  needs the developer's word (`HOW_WE_WORK.md` ## Developer).
- Backlog: CAS-30 (SIMD/assembly multiplication kernels) deferred by the developer (rule `rules/L1.md`).
- PA3 has no `GENERATION_PLAN.md` or `PROJECT_CONTEXT.md` yet; the roadmap so far lives in Lific.

## Decided; do not undo
- Tests, fuzzing and benchmarks run only on GitHub Actions (`.github/workflows/ci.yml`, `fuzz.yml`, `workload.yml`),
  never on the dev machine; local `cargo check`/`cargo build` and clippy are fine (`AGENTS.md`, `README.md`). Fuzz
  jobs stay within `-rss_limit_mb=2048 -malloc_limit_mb=1024` (CI memory).
- Work on master, no feature branches; a red master is tolerated until the first release (developer, 2026-10-25).
  Under PA3 the developer pushes; agents never do.
- `unsafe_code = "deny"` workspace-wide; the one exemption is `crates/cas-giac/src/library.rs` (`Cargo.toml:26`).
- CAS features and giac parity before kernel micro-optimisation (`rules/L1.md`).
- Every algorithm change updates its chapter in `crates/cas-num/doc/algorithms/` or `crates/cas-poly/doc/algorithms/`
  in the same commit; our own tricks are marked with their ticket.
- Numeric answers are certified balls (e.g. FSolve roots carry a sign change across the result; NLPSolve certifies
  KKT points); departures from giac are recorded in the design docs (`crates/cas-giac/README.md` ## Tests).
- The local Lific server never runs (`lific start`/`init`/`service install` are forbidden).

## Where details live
Card `HOW_WE_WORK.md` (tools, gate, conventions, living docs); techniques `cookbook/INDEX.md` (CI A/B C0001, fuzz
diagnosis C0002, stable clippy C0003, gh quota C0004); old memories `.claude-state/memory/gen3.2.md`; document list
`docs/PRE_PA3_INDEX.md`.
