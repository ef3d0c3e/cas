# How we work — cas

<!-- The card: standing facts for every agent. Never appended to; the cap is
     card.max_chars in .claude/pa.json (default 10000), enforced by the archive.
     Each role prints its slice with `tools/card.py slice <role>`.
     Written at install (interview) and edited in the same task as whatever changed it (H7). -->

## Project <!-- roles: claude expert coder planner review auditor curator -->
cas — a computer algebra system engine in Rust (crates `cas-num`, `cas-poly`, `cas-core`, `cas-lang`; `cas-giac` the
test oracle; `cas-tune`); Xcas/giac is the reference for defaults. Work is tracked in Lific (project `CAS`). Constitution `PROJECT_CONTEXT.md`; roadmap `GENERATION_PLAN.md`;
rules `rules/INDEX.md`; techniques `cookbook/INDEX.md`; ops detail `docs/ops/INDEX.md`.

## Developer <!-- roles: claude planner review auditor curator -->
Who: ef3d0c3e. Experience: advanced. Domain: cas.
Preferences: recommendations, not questions. Plain-English recaps at phase end. A notification
whenever anything waits on them. Autonomy: full inside an approved plan; stop only at the two gates. Notification channel: toast.
They ratify rules at the next plan approval.
- Pushing: agents push master themselves (`git push origin master`) when CI needs to run; no push stop for the
  developer. The remote is for CI, nothing released to preserve. Still no force-push, no history rewrite.
- Area order: when an area (Lific `CAS-PLAN-n`) closes, Claude plans the next one in its own recommended order, no
  which-area question; a brand new Lific module needs their word (existing: Algebra, Balls, Calculus, Evaluator,
  Expressions, Infrastructure, Integers, Language, Polynomials, Statistics, Numerical analysis,
  Signal processing). The plan gate still applies.
- Keep the dev machine free of test/fuzz load; be careful with CI's memory limits.

## Tools <!-- roles: expert coder planner review auditor curator -->
PY = /usr/bin/python3

| Name | Command | Purpose |
|---|---|---|
| launch | `PY tools/launch.py --state-json` | the state detector: the mod's read-only report; never typed by the developer |
| plan_edit | `PY tools/plan_edit.py` | the only writer of `PHASE_PLAN.md` (status, Changes, add/reopen) |
| status | `PY tools/status.py` | `.run/status.json` for the statusline; INBOX consume; waiting flags |
| task_log | `PY tools/task_log.py` | lint and finish a task summary |
| research_add | `PY tools/research_add.py` | allocate a report id, write the index line |
| rules_add | `PY tools/rules_add.py` | add, supersede, promote, retire a rule |
| skill_add | `PY tools/skill_add.py` | turn a `workflow:` gotcha into `.claude/skills/<name>/SKILL.md` |
| cookbook_add | `bash tools/cookbook_add.sh` | add a cookbook entry and its index line |
| phaseend_index | `PY tools/phaseend_index.py` | assemble and archive a PhaseEnd; `history` rebuilds `docs/PHASE_HISTORY.md` |
| genend_index | `PY tools/genend_index.py` | assemble a GenerationEnd |
| commit_task | `bash tools/commit_task.sh` | the only commit path; explicit paths, no trailers, never pushes |
| run | `bash tools/run.sh` | any command that may print >40 lines |
| curate | `PY tools/curate.py` | `memory --route\|--reindex`, `rules --reseed\|--demote`, `cookbook --demote`, `claude-md\|card --refresh`, `living-docs --entry`, `how-we-work --retire`; `--dry-run` first |
| lific | `lific issue list --project CAS --backend http --json` | the tracker, remote only (lific.pundalik.org; MCP via `.mcp.json`); CLI always `--backend http`; never `lific start`/`init`/`service install` |
| gh | `gh run view --job <id> --log`, `gh api repos/ef3d0c3e/cas/actions/jobs/<id>/logs` | CI results; logs grep'd, never dumped; quota trap C0004 |

## Skills <!-- roles: expert planner -->
<!-- one line per captured workflow; the SKILL.md is the canonical text -->
- (none yet)

## Paths <!-- roles: claude expert coder planner review auditor curator -->
- Oracles (the ground truth X3 names): libgiac through `crates/cas-giac` (tests compare with Xcas's answers; CI job `oracle`); giac lookups outside tests: `PY scripts/caseval.py <inputs.cas> [timeout_s]` (libgiac caseval, one child per line, C0018), `PY scripts/aide.py <name>... [--cas]` (aide_cas entries and examples), `bash scripts/giac_gap.sh <outdir>` (missing names), `PY scripts/gap_describe.py missing|all|lines` (the missing entries described)
- Generated (never hand-edited, H1): `Cargo.lock`, `target/`, `crates/*/proptest-regressions/` (proptest seeds)
- Hand-edited: `crates/*/src`, `crates/*/doc/algorithms/`, `crates/cas-lang/proto/` (Python printer ports, R74), `crates/cas-core/proto/` (1.7 statistics and signal prototypes, giac probe inputs; README), `crates/cas-lang/names/en.txt`, `fuzz/`, `.github/workflows/`
- Scratch: all of a phase's work in `.run/` goes under `.run/<phase>/` (gitignored; never the system temp), nowhere
  else in `.run/`. It is a temp folder only, never a home for long-term content: what the project needs lives in its
  proper place in the repo. The phase's cleanup task checks it and
  the close deletes it; a tool deletes its own temp folders when it finishes.

## Build / run / test <!-- roles: claude expert coder planner review auditor curator -->
- Build: `cargo build`, `cargo check` (local is fine)
- Run: library crates; `cargo run -p cas-tune` writes cas-num's tuning configuration
- Our answers are seen only on CI: no local scratch/example programs or binaries that print them (developer, 1.4
  close); Python prototypes and giac probes (`scripts/caseval.py`) are fine.
- Test / the gate: CI only, never `cargo test`/`cargo fuzz`/benchmarks locally (`AGENTS.md`). Green means every job of
  `.github/workflows/ci.yml` (fmt, clippy on stable, test, fuzz-check, oracle) and `fuzz.yml` passes on the pushed commit
- Modes: `gh workflow run fuzz.yml -f time_per_target=<s> [-f targets="<names>"]`; perf A/B `workload.yml` (C0001)
Long gates go through `tools/run.sh` with a raised timeout, never a poll loop.

## Conventions & house style <!-- roles: claude expert coder planner review auditor curator -->
- Work on master, no feature branches; a red CI on master is tolerated until the first release, but fix what CI or the
  fuzzer finds right away. Commit subject names the ticket (`CAS-N: …`); close the ticket in Lific (status `done`).
- An algorithm change updates its chapter in `crates/*/doc/algorithms/` in the same commit (base algorithm, tricks,
  error bound, references; our own tricks marked with their ticket).
- A fuzz finding gets a regression test in the owning crate's suite. `unsafe_code = "deny"`, see rule L1.
<!-- project-specific only; the general house style is in the project-architect skill -->

## Environment <!-- roles: claude expert coder planner review auditor curator -->
- OS / shells: Linux / bash
- Python: /usr/bin/python3 (3.14.10)
- Pins: local rustc nightly (bin dir on PATH, no rustup proxies); CI uses stable. Stable clippy locally: C0003
- Harness gotchas that bite here: GitHub API quota 5000/h (C0004); `gh run view --log` sometimes empty, use `gh api`

## Docs map <!-- roles: expert coder planner review auditor curator -->
- `docs/ops/INDEX.md` — setup, env, pins, per-topic ops notes
- `cookbook/INDEX.md` — techniques, grep by tag
- `rules/INDEX.md` — full rule texts
- `phase-ends/TASK_INDEX.md`, `phase-ends/RESEARCH_INDEX.md` — what was done and what was learned
- `docs/research-archive/` — migrated legacy reports
- `docs/retired/` — everything moved out of the load order

## Living docs <!-- roles: claude expert planner curator -->
- `crates/cas-num/doc/algorithms/README.md`: cas-num ball-arithmetic algorithms, error bounds, references; update in the same commit as the algorithm
- `crates/cas-poly/doc/algorithms/README.md`: cas-poly algorithms (dense, sparse, gcd, factor, integrate, groebner); update with the algorithm
- `crates/<crate>/README.md`: each crate's layout, testing, benchmarks; update when modules or test setup change
- `fuzz/README.md`: the libFuzzer targets and running them; update when a target is added or changed
- `AGENTS.md`: Lific usage and the CI-only test policy; update when the tracker or CI workflows change
