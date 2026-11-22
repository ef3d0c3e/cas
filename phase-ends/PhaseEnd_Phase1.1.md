# PhaseEnd — Phase 1.1: Green master (implements Gen 1.1)
Approved: 2026-11-22 | Closed: 2026-11-22 | Planner: claude-opus-5-5/medium | Tasks: 8 done, 0 superseded

## Milestone
Milestone: Master's CI (fmt, clippy, test, fuzz-check, oracle) and Fuzz (every target) pass on the developer's push of two fixes, the eval target's divide-by-zero under a regression test in its owning crate and `boundary_values` checking BVPSolve against the closed form, and CAS-517..524 are done in Lific; `gh run list --branch master --limit 2` shows both runs green.
Confirmed by the developer at the close.

## Tasks
- T1 — fix the eval fuzz divide-by-zero at its source, with a regression test | Done: `FromContinuedFraction` with a zero convergent denominator gives `ComplexInfinity` instead of panicking; regression test in cas-core; CAS-525 open in Lific. Summary: The nightly fuzzer crashed because turning the continued fraction {-2, -1, 0, 1} back into a number divides by zero (it is -2 + 1/0). That case now gives ComplexInfinity, the same as 1/0 anywhere else, with the usual warning. The crash input replayed locally no longer fails, and a test with the found expression guards it on CI. | commit: dee4f4b | tasks/T1.md | logs/T1.md
- T2 — boundary_values checks BVPSolve against the closed form (CAS-520) | Done: `boundary_values` checks n = 4, 2 against closed forms, giac only at h = 1/20; the n = 4 giac error is a recorded departure. Summary: The oracle test failed because giac's own boundary-value solver is off by 0.04 on a coarse grid, while ours is exact. The test now checks our small-grid answers against the exact solutions, compares with giac only on fine grids where giac is accurate, and records the coarse-grid disagreement as an intended departure. The docs' claim that giac is good to about 10^-5 was corrected. | commit: bfd3e78 | tasks/T2.md | logs/T2.md
- T3 — baseline chores: giac_gap.sh into scripts/, modules on the card | Done: `scripts/giac_gap.sh <outdir>` writes the six gap files (identical to the scratch outputs); cas-giac README names it; card's module list has 12 modules. Summary: The planning-time script that compares giac's command list with our names now lives in the repo and takes its output folder as an argument. The cas-giac README points to it, and the card's list of Lific modules now includes Numerical analysis and Signal processing, confirmed at generation approval. | commit: - | tasks/T3.md | logs/T3.md
- T4 — local gate green | Done: fmt and stable clippy (`-D warnings`, all targets) pass on the workspace; `cargo check --all-targets` in `fuzz/` passes on nightly; tracked tree clean. Summary: The local checks that mirror CI's formatting, lint and fuzz-build jobs all pass on the commit after T1 to T3, with no code changes needed. The tree is ready for the developer to push in T5. | commit: - | tasks/T4.md | logs/T4.md
- T5 — push stop: developer pushes master | Done: 21 local commits on master (866047b..8ee13e9) are listed for the developer to push; push pending Summary: Agents never push, so the phase stops here for the developer to push master. Pushing starts CI and the fuzz workflow on the fixed commits, which is the real check of T1 (divide-by-zero fix) and T2 (boundary_values fix). Master is 21 commits ahead of origin/master; nothing else is pending except an untracked session.tex that is not ours. | commit: - | tasks/T5.md | logs/T5.md
- T6 — CI check on the pushed head | Done: master head fa9ae22 is green on GitHub Actions: CI run 45047911216 (fmt, clippy, test, fuzz-check, oracle) and Fuzz run 45047911210 (all 55 targets). Summary: After the developer pushed the phase's fixes, both GitHub workflows ran on the pushed commit and passed. The oracle boundary_values failure and the eval divide-by-zero fuzz crash that kept master red are gone. No fix and no T5.1 round were needed. | commit: - | tasks/T6.md | logs/T6.md
- T7 — close CAS-517..524 and the fix's ticket in Lific | Done: CAS-517..524 and CAS-525 (T1's fix ticket) are at status `done` in Lific. Summary: The eight tickets that came with the area work (CAS-519 from backlog, the rest from active) and the ticket for T1's divide-by-zero fix were closed in the tracker. Checked by listing the project's done issues. | commit: - | tasks/T7.md | logs/T7.md
- T8 — phase cleanup | Done: every item in `.run/1.1/` checked; each is already in the repo (test, cookbook entry, script) or regenerable, so nothing moved; `.run/1.1/` is safe to delete at the close. Summary: The phase's scratch folder was audited before the close deletes it. The crash input is covered by a regression test, the two cookbook drafts are banked as entries, the test drafts are committed, and the decode crate, the cargo-fuzz install and the gap lists can all be rebuilt from documented recipes. Nothing needed moving. | commit: - | tasks/T8.md | logs/T8.md

## Decisions that still bind
-

## Rules proposed
- T1 — none

## Cookbook entries added
C0005 | Locating an eval fuzz panic and replaying it locally | fuzz,eval,panic,replay,cargo-fuzz,backtrace | 2026-11-22 | 1.1/T1 | T1 CAS-525
C0006 | Numeric oracle tests: closed form first, oracle only inside its own error | giac,oracle,numeric,tolerance,bvpsolve,departures | 2026-11-22 | 1.1/T2 | CAS-520 boundary_values

## Research
R0-001 | planning the fix of master's red CI job `oracle` (CAS-520 BVPSolve), test `boundary_values` | boundary_values oracle failure: BVPSolve vs giac bvpsolve | giac, oracle, bvpsolve, BVPSolve, ndsolve, shooting, tolerance | retriever-code | 2026-11-22 | 57 lines
R0-002 | which programming constructs the CAS surface language supports, and what is recorded as planned but not done | CAS surface-language programming constructs: what exists, what is planned | cas-lang, cas-core, programs, control-flow, strings, session, planned | retriever-code | 2026-11-22 | 60 lines
R0-003 | eval fuzz red on master (CI 42057178575); locate the panic origin, reachable builtins, regression test homes | Divide-by-zero panic at cas-num rational/mod.rs:52 and the CAS-517..524 builtins | cas-num, rational, fuzz, eval, CAS-517, CAS-524, divide-by-zero | retriever-code | 2026-11-22 | 68 lines
R0-004 | group giac 1.9.0 help entries missing from the gen1 names table (entries_missing.txt) into candidate CAS areas, plus leftovers that fit existing modules and near-name flags | Giac gap areas: the 1077 uncovered non-unit entries, grouped | giac, gen1, gap, areas, names, parity | retriever-digest | 2026-11-22 | 71 lines
R1.1-005 | exact signatures and shapes for the phase plan's Interfaces section (repo /home/user/code/cas) | Interface shapes for phase 1.1 (fuzz divide-by-zero, oracle boundary_values) | interfaces, rational, ndsolve, bvp, giac, fuzz, regression-tests, scripts | retriever-code | 2026-11-22 | 40 lines

## Audit
- seed: median 22k, max 23k, n=8 (phase 1.1)
- CLAUDE.md: 1109 bytes
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes
- .claude/skills/pa-planning/SKILL.md: 8419 bytes
- .claude/skills/project-architect/SKILL.md: 13748 bytes
- .claude-state/memory/MEMORY.md: 53 bytes
- HOW_WE_WORK.md: 7545 bytes
- cookbook/INDEX.md: 1154 bytes
- rules/INDEX.md: 3096 bytes
### Carry audit — phase 1.1 (2026-11-22T16:41:44Z → open UTC, 1 sessions, 136 requests)

| file (read) | chars | n |
|---|---|---|
| sparse.rs | 45.2k | 2 |
| print.rs | 36.7k | 5 |
| builtin.rs | 31.8k | 2 |
| integrate.rs | 31.5k | 4 |
| heap.rs | 23.9k | 2 |
| tests.rs | 20.8k | 2 |
| tests.rs | 20.8k | 1 |
| modular.rs | 19.0k | 1 |
| series.rs | 18.5k | 4 |
| layout.rs | 18.0k | 2 |
| file (write) | chars | n |
|---|---|---|
| doc31.md | 20.0k | 1 |
| doc34.md | 12.2k | 1 |
| doc15_new.md | 12.2k | 1 |
| doc11-v2.md | 11.3k | 1 |
| product.rs | 10.4k | 1 |
| doc10.md | 10.4k | 1 |
| doc14.md | 9.9k | 1 |
| doc16.md | 9.7k | 1 |
| result kind | chars | n |
|---|---|---|
| Read | 798.8k | 171 |
| bash other | 156.2k | 123 |
| Grep | 98.7k | 48 |
| tools/plan_edit.py | 50.0k | 11 |
| other | 12.8k | 66 |
| Write | 11.0k | 63 |
| run.sh | 7.8k | 17 |
| tools/commit_task.sh | 4.4k | 5 |
| tools/task_log.py | 3.3k | 3 |
| Edit | 2.9k | 18 |
- whole reads over 20.0k: 3
- seed floor: retriever-code 5,647 (n 4, prev —) · retriever-digest 6,304 (n 1, prev —)
- outline credit: 0 outlines · 0 followed by a ranged read · 0 chars credited
- warm pings: 0 pings over 0 runs, 0 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.0000 + relay $0.0000 + re-arm $0.0000 = $0.0000 vs rewrites replaced $0.0000, cheaper than one rewrite: n/a, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 153, relay 0, re-arm 3, other 38, relay cost $0.0000
- retriever re-asks: 0 of 0 retriever briefs repeat a lookup of the same run
- noise: 6 lines 233 chars — usage: : 5 lines/177, No such file or directory: 1 lines/56
- price: read $0.20/Mtok · 1h write $8.00/Mtok · output $20.00/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 8439 chars, 136 requests
- flag: 2 whole-plan reads by pa-expert
- flag: 2 tool-source reads by main, retriever-digest
  /home/user/code/cas/.claude/settings.local.json  61 chars  main
  /home/user/.claude/projects/-home-user-code-cas/1b3aa32d-02f2-407f-9868-a62b20c7002c/tool-results/toolu_01TWmht7psJbdZCJHXw9nCri.txt  133 chars  retriever-digest
- flag: 3 whole reads over 20.0k by main
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-core/src/builtin.rs  role: main  chars: 23.1k

## Agent runs
- T8 | pa-claude | partner | claude-opus-5-5 | medium | ctx 91k | $1.682 | saved $8.77 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 69k | $1.156 | saved - | completed | parent 1b3aa32d
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 62k | $0.999 | saved - | completed | parent 1b3aa32d
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 37k | $0.345 | saved - | completed | parent 1b3aa32d
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 34k | $0.265 | saved - | completed | parent 1b3aa32d
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 32k | $0.236 | saved - | completed | parent 1b3aa32d
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 36k | $0.338 | saved - | completed | parent 1b3aa32d
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 34k | $0.305 | saved - | completed | parent 1b3aa32d
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 42k | $0.377 | saved - | completed | parent 1b3aa32d

## Discussions
- (none)

## Deferred
- from T8: the close may delete `.run/1.1/` whole. (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-22 router: T1 next -> done
- 2026-11-22 router: T2 next -> done
- 2026-11-22 router: T3 next -> done
- 2026-11-22 router: T4 next -> done
- 2026-11-22 router: T5 next -> done
- 2026-11-22 router: T6 next -> done
- 2026-11-22 router: T7 next -> done
- 2026-11-22 router: T8 next -> done

## Plain-English Recap
Phase 1.1 Green master: master's CI and Fuzz are green again on fa9ae22 (CI run 45047911216, all five jobs; Fuzz run 45047911210, all 55 targets). The eval fuzz crash came from FromContinuedFraction (CAS-424), not CAS-517..524: a zero last denominator (-2 + 1/0) reached Rational::new; it now gives ComplexInfinity, under a regression test in cas-core (CAS-525). boundary_values now checks BVPSolve against closed forms on small grids and giac only at step 1/20; giac's n=4 error (0.044) is a recorded departure, and the accuracy comments were corrected (CAS-520). giac_gap.sh moved to scripts/ with an output-folder argument; the card lists Numerical analysis and Signal processing. CAS-517..525 are done in Lific. Decided during the phase: agents push master themselves when CI needs to run (card 705b703, after the developer removed the never-push fail-safe from CLAUDE.md). Open: .claude/settings.json has a trailing comma from the developer's edit (his to fix); the generation plan's standing constraints still describe a push stop, so each later phase brief says agents push; lific issue list caps at 500, so ticket checks use issue get; checkout@v4 Node 20 deprecation warning.
