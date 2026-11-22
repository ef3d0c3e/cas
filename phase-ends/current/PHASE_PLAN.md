# Phase 1.1 — Green master        (implements GENERATION_PLAN.md phase 1.1)
Milestone: Master's CI (fmt, clippy, test, fuzz-check, oracle) and Fuzz (every target) pass on the developer's push of two fixes, the eval target's divide-by-zero under a regression test in its owning crate and `boundary_values` checking BVPSolve against the closed form, and CAS-517..524 are done in Lific; `gh run list --branch master --limit 2` shows both runs green.
Approved: 2026-11-22   Planner: claude-opus-5-5/medium   Plan-hash: 62a88689aaecfc60e651c2146962869a600f42e5df41534a2e7a4df9f4145eb1

## Context
- Master red at `866047b` (docs/PRE_PA3_STATE.md:27-32): CI run 42057178857 job `oracle` fails `boundary_values`; Fuzz run 42057178575 target `eval` dies on `attempt to divide by zero` at crates/cas-num/src/rational/mod.rs:52:41. fmt, clippy, test, fuzz-check pass. Both came with the uncommitted-to-Lific tickets CAS-517..524 (CAS-519 backlog, the rest active).
- Crash artifact `fuzz-artifacts-eval-42057178575` (327 bytes) is live until 2027-02-03 (`gh api repos/ef3d0c3e/cas/actions/runs/42057178575/artifacts`, checked 2026-11-22). Download with `gh run download 42057178575 -n fuzz-artifacts-eval-42057178575 -D .run/1.1/crash`.
- Divide-by-zero: `Rational::new` panics; `Rational::checked_new` and `checked_recip` exist (R1.1-005). The callers not yet audited and the eval BUILTINS entries of CAS-517..524 are in R0-003. The crash may predate CAS-517..524: new names remap the cached corpus (C0002), so the input may reach an older builtin.
- boundary_values: our BVPSolve gives y(2) = 12 exactly for `bvpsolve((32+2x^3-y*diff(y(x),x))/8,[x=1..3,y],[17,43/3],4)`, which is the closed form y = x^2 + 16/x; giac's bvpsolve at n = 4 gives about 12.02, outside the test's 10^-4 (1 + |b|) tolerance (R0-001, R1.1-005). giac is the one that is wrong here; CAS-520 asked for closed forms.
- Verification: tests and fuzzing run only on GitHub Actions. Locally: `cargo check`, `cargo build`, clippy, and exactly one fuzz-input replay (fuzz/README.md:98-102) under `scripts/capped.sh`. Agents never push; the developer pushes at T5.
- Lific is remote: CLI always `--backend http`, never a local server. Commit lines carry the ticket after the task id (`T1.c1: CAS-525: …`).
- Modules confirmed at generation approval: Numerical analysis and Signal processing join the card's list (T3).
- Scratch: `.run/1.1/` only. The generation planner's scratch `.run/logs/gen1-gap/` holds `giac_gap.sh` and its outputs (R0-004's inputs).

## Rationale
- Two fixes, then one push stop, as the developer kept (no area pulled into 1.1). Each fix is its own task: different crates, different evidence, independent.
- T1 replays the input locally rather than through a C0002 diagnostic branch: a branch needs a push, and the one local replay is allowed. Decoding the input by hand (C0002) is the fallback if the replay does not reproduce.
- T2 does not loosen the tolerance: a looser bound would hide giac's error, not fix the test. Closed forms on small grids, giac only where its own error is inside 10^-4 (large n), and the departure recorded in the existing pattern.
- The local gate (T4) is a task of its own so the push stop gets one green, committed head.
- Lific tickets close only after green CI (T7), as the card's conventions say.
- No `effort: high`: every done-when here is arbitrated by CI or by a command.

## Interfaces
- `Rational::new(num: Integer, den: Integer) -> Rational` crates/cas-num/src/rational/mod.rs:51; panics at :52 `.expect("attempt to divide by zero")` over `Rational::checked_new(num, den) -> Option<Rational>` :57. `checked_recip(&self) -> Option<Rational>` :139; `recip` panics :154. Doc at :22 names a `checked_div` that does not exist (R1.1-005).
- `boundary_values` crates/cas-giac/tests/ndsolve.rs:91: inputs at :98-101 (n = 4, 4, 4, 2, last with `output=list`); compares `giac.value(input)` with `ours(input)` via `close_points` :68; tolerance `scale / 10_000` :78 (10^-4 (1 + |b|)); the "about 10^-5" claim :66-67 describes giac's error.
- `pub(in crate::native) fn bvp_solve(ev: &mut Eval<'_>, e: &Expr) -> Native` crates/cas-core/src/native/ndsolve/bvp.rs:133; module doc :19-23 (output `{{x_k, y_k}, …}` at n + 1 even points; "about 10^-5" claim about giac).
- Departures pattern: `fn departures()` crates/cas-giac/tests/polynomials.rs:274; cases `(input, giac_printed, ours_decided)` each with a `//` reason line, checked with `giac.eval` and `giac.agrees`. crates/cas-giac/README.md `## Tests` :75-98, one bullet per test file.
- Regression-test patterns: crates/cas-core/src/eval/tests.rs:425-428 (fuzz finding named in the doc comment with its ticket, CAS-149); crates/cas-num/src/integer/tests.rs:318-320 (`#[should_panic(expected = "attempt to divide by zero")]`).
- Replay: `cargo fuzz run -O eval <input>` (fuzz/README.md:100), run as `scripts/capped.sh cargo fuzz run -O eval <input>` from `fuzz/`; `const BUILTINS: [&str; 429]` fuzz/fuzz_targets/eval.rs:37.
- `.run/logs/gen1-gap/giac_gap.sh`: `set -eu`, run from repo root; hardcoded `out=.run/logs/gen1-gap` :5; reads `/usr/share/giac/aide_cas` :6,:13 and `crates/cas-lang/names/en.txt` `[builtins]`..`[types]` :7; writes giac.txt, ours.txt, missing.txt, present.txt, entries_covered.txt, entries_missing.txt.

## Cookbook
- C0002 — diagnosing slow, hanging or OOM fuzz inputs (decode the input; no diag branch without a push) (T1)
- C0003 — CI's stable clippy locally (T4)
- C0004 — gh API quota: one watcher per run at `--interval 60` (T6)

## Research
- R0-001 — boundary_values oracle failure: BVPSolve vs giac bvpsolve (T2)
- R0-003 — the divide-by-zero at rational/mod.rs:52, its callers, the eval BUILTINS of CAS-517..524 (T1)
- R1.1-005 — interface shapes for this phase: Rational constructors, boundary_values, bvp_solve, departures and regression patterns, replay command, giac_gap.sh (T1, T2, T3)

## Developer decides
- (none)

## Triage
- (none)

## Tasks

- T1 | done | title: fix the eval fuzz divide-by-zero at its source, with a regression test | files: crates/cas-num/src/rational/mod.rs, the offending caller (R0-003), the owning crate's tests | done-when: the crash input from run 42057178575 is replayed once locally and no longer panics, the caller that passed a zero denominator is fixed (returns an error/undefined instead of calling `Rational::new` with 0), a regression test with that input's expression sits in the owning crate's suite naming the fuzz finding and its ticket, `cargo check --all-targets` is green, and a Lific ticket for the fix exists in the owning module | verify: cd fuzz && ../scripts/capped.sh cargo fuzz run -O eval ../.run/1.1/crash/<file>
  Research: R0-003 (callers, BUILTINS), R1.1-005 (Rational API, test patterns). Cookbook: C0002. Download the artifact into `.run/1.1/crash/`. Only this one replay may run locally; if it does not reproduce, decode the input (C0002) and reason from R0-003's callers. Fix at the source, not by making `Rational::new` lenient. The ticket goes in the module of the builtin at fault (Numerical analysis or Signal processing for CAS-517..524 builtins, else the existing module); commit lines `T1.c<k>: CAS-<n>: …`.
- T2 | done | title: boundary_values checks BVPSolve against the closed form (CAS-520) | files: crates/cas-giac/tests/ndsolve.rs, crates/cas-core/src/native/ndsolve/bvp.rs, crates/cas-giac/README.md | done-when: `boundary_values` compares small-grid inputs (n = 4, 2) with closed-form values, giac is compared only on inputs where its own error is inside 10^-4 (large n), the n = 4 disagreement is recorded as a departure in the polynomials.rs:274 pattern and in crates/cas-giac/README.md ## Tests, the "about 10^-5" claims at ndsolve.rs:66-67 and bvp.rs:19-23 are corrected, the tolerance is unchanged, and `cargo check --all-targets -p cas-giac -p cas-core` is green
  Research: R0-001, R1.1-005. Commit lines `T2.c<k>: CAS-520: …`. Choose the large n by giac's documented/observed error from R0-001, not by trial runs (no local tests). If BVPSolve's algorithm doc changes, update its chapter in the same commit.
- T3 | done | title: baseline chores: giac_gap.sh into scripts/, modules on the card | files: scripts/giac_gap.sh, crates/cas-giac/README.md, HOW_WE_WORK.md | done-when: `scripts/giac_gap.sh <outdir>` exists, takes its output folder as its argument and writes the same six files as the scratch copy (one run into `.run/1.1/gap/` shows them), crates/cas-giac/README.md has one line naming it, and the card's Developer module list reads the ten modules plus Numerical analysis and Signal processing
  Research: R1.1-005 (script's hardcoded paths). The card is edited by hand at the one line (HOW_WE_WORK.md ## Developer); keep it under its cap.
- T4 | next | title: local gate green | deps: T1, T2, T3 | done-when: `cargo fmt --all --check` and `cargo clippy --workspace --all-targets -- -D warnings` pass with the stable toolchain's bin dir first on PATH (C0003), `cargo check --all-targets` in `fuzz/` passes on nightly, and the tree is committed clean | verify: bash tools/run.sh gate11 -- cargo clippy --workspace --all-targets -- -D warnings
  Cookbook: C0003. No research needed (commands only).
- T5 | title: push stop: developer pushes master | deps: T4 | review: yes | coder: none | done-when: the developer says master is pushed and `git status -sb` shows master level with origin/master
  No research needed. Present the commits since `866047b` in one short list with their tickets.
- T6 | title: CI check on the pushed head | deps: T5 | coder: none | done-when: `gh run list --branch master --limit 2` shows CI and Fuzz green on the pushed head; a red job is fixed in this phase (fuzz finding: regression test in its owning crate) and goes back through T5 as T5.1 | verify: gh run list --branch master --limit 2
  Cookbook: C0004 (one `gh run watch --interval 60` per run; logs via `gh api repos/ef3d0c3e/cas/actions/jobs/<id>/logs`, grep'd). Fuzz takes about 45 minutes. No research needed.
- T7 | title: close CAS-517..524 and the fix's ticket in Lific | deps: T6 | coder: none | done-when: `lific issue list --project CAS --backend http --json` shows CAS-517..524 and T1's ticket at status `done`
  No research needed. CAS-519 sits in backlog; set it `done` like the rest.
- T8 | title: phase cleanup | deps: T7 | done-when: what the project still needs from .run/1.1/ is in its proper place in the repo
  Also note in the summary that `.run/logs/gen1-gap/` is regenerable by `scripts/giac_gap.sh` (T3); leave it for later phases unless the developer says otherwise. No research needed.

## Risks
- The replay does not reproduce (different toolchain or corpus): detect at T1's replay; fall back to decoding the input (C0002) and R0-003's callers; if still unknown, `gh workflow run fuzz.yml -f targets=eval` after the push stop.
- Fuzz finds a second, older crash once eval stops dying early (remapped corpus, C0002): detect at T6; fix in-phase with another push stop (T5.1).
- giac's error at large n is not inside 10^-4 for some input: detect from R0-001's numbers at T2; drop giac comparison for that input and record it as a departure.
- The developer is away at T5: the loop waits there; nothing else is queued.

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-22 router: T1 next -> done
- 2026-11-22 router: T2 next -> done
- 2026-11-22 router: T3 next -> done
