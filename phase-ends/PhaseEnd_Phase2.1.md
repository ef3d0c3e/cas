# PhaseEnd — Phase 2.1: Engine carry-overs (I3) (implements Gen 2.1)
Approved: 2026-11-27 | Closed: 2026-11-28 | Planner: claude-opus-5-5/medium | Tasks: 13 done, 0 superseded

## Milestone
Milestone: The eight natives of CAS-650 no longer re-evaluate an argument under fresh limits, markov and lpsolve have random oracle sweeps, the five rational Hilbert rows agree with giac, and the needless rewrite at crates/cas-giac/tests/euler_lagrange_assume.rs:45 is gone; CI and Fuzz pass on the pushed head and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — Lific tickets and design notes for the phase, Fill backlog ticket in Language (I3) | Done: CAS-650 carries the guard decision and the eight sites (comment, status todo); CAS-651..655 filed for T3-T6; CAS-656 (Language, backlog) records Fill inside begin...end. Summary: Every change this phase makes now has its ticket in the module that already owns the code. The eight-site guard decision is on the original ticket, and Fill inside begin...end is in the Language backlog for when TI-style programs are wanted. No design page needed an edit. | commit: - | tasks/T1.md | logs/T1.md
- T2 — is_numeric guard before approximation in the eight CAS-650 natives | Done: the eight CAS-650 helpers refuse a non-numeric evaluated argument before `approximation` (call stays); nesting tests + numeric rows per entry native; docs name CAS-650; CI and Fuzz green on 9317acb. Summary: Eight numeric commands used to re-evaluate a symbolic argument from scratch, so nesting them inside their own arguments cost time doubling per level. Each now checks that the argument is a number (including forms like pi/4 or sqrt(2)) before approximating it, and otherwise leaves the call unevaluated at once. Each command has a test nesting it 30 deep that must finish in under two seconds, and a test showing a numeric non-rational argument still gives an answer. | commit: 9317acb | tasks/T2.md | logs/T2.md
- T3 — classify the remaining approximation call sites and guard those that can nest | Done: CI 46019201247 and Fuzz 46019201213 green on 17b70e4; all 51 non-T2 `approximation(` sites in native/ carry one verdict; the 17 symbolic-native sites are guarded with `is_numeric` (CAS-650), 12 nesting tests added; the numeric natives' 12 sites go to T10/T11. Summary: Many functions double their work for every level of nesting when a held call of themselves sits inside their argument, because they re-evaluate it from scratch to get a number. This task checked every remaining place that does this, found 17 in symbolic functions (Limit, Series, Sum, Solve, the matrix decompositions, Integrate bounds), and made each refuse non-numeric input there instead, with a 30-level nesting test each. The numeric solvers' sites were split off to two later tasks because a plain guard would refuse legitimate inputs there. | commit: 17b70e4 | tasks/T3.md | logs/T3.md
- T4 — random oracle sweeps for markov (exact stochastic matrices) and lpsolve (objective value) (I3) | Done: `random_markov_chains` (markov.rs:569) and `random_lp_objectives` (optimization.rs:906), 128 cases each through `agree`, both pass on CI 46037332695; README bullets at :378 and :1082. Summary: Two new random sweeps compare our Markov chain analysis and our linear program optimum with giac's on 128 fresh random inputs each CI run. Inputs where giac is known to answer with decimals, loop forever, error, or give a wrong optimum are skipped by looking at the drawn input, with the reasons written beside the tests. Both sweeps pass on CI. The run is still red because of an older, unrelated derivative sweep that drew a new failing input. | commit: 0bc06dc | tasks/T4.md | logs/T4.md
- T5 — residue rule for rational Hilbert transforms, five departures to agreeing rows (I3) | Done: `hilbert(f,t)` of a proper rational function with exact non-real poles answers exactly by residues; the five rational rows agree with giac on CI; its rows green on CI 46047327203 (hilbert_departures, expression_rows), Fuzz 46047327115 green. Summary: The Hilbert transform of rational functions such as 1/(1+t^2) now gives giac's closed form instead of staying unevaluated: the poles are found exactly, split by half-plane, and the transform is built from their principal parts. The five rational rows moved from departures to agreeing rows and pass against giac on CI. CI also showed that a random draw in the Sturm-sequence oracle test found a genuine, unrelated disagreement with giac for -x^3; that is now its own task. | commit: e0f2def | tasks/T5.md | logs/T5.md
- T6 — drop the needless ,t,2) rewrite in euler_lagrange_assume.rs (I3) | Done: `under_an_assumption` reads giac's answer as-is (no `,t,2)`→`,t,t)` rewrite); the oracle job on 1e131c4 is green. Summary: The euler_lagrange oracle test used to rewrite giac's second-derivative notation before parsing it. Since our parser now reads giac's `diff(y(t),t,2)` form, that rewrite was dropped; the test still passes against giac on CI. | commit: 1e131c4 | tasks/T6.md | logs/T6.md
- T7 — CI and Fuzz green on the pushed head, the phase's tickets done | Done: master 6b780d8 pushed; CI 46057029170 (all five ci.yml jobs) and Fuzz 46057029281 green; CAS-650..655, 657, 658, 660 done, CAS-656 backlog. Summary: The phase's code was pushed and the full CI and fuzz runs checked on it. One oracle sweep failed: a random polynomial where giac wrongly says "irreducible" and our engine correctly says it factors. Probing giac showed it misses any factor in y alone, so the sweep now skips exactly that class, the failing input is pinned as a recorded giac error, and the second run was fully green. All phase tickets except the Fill backlog ticket are closed. | commit: 6b780d8 | tasks/T7.md | logs/T7.md
- T8 — phase cleanup | Done: every file in .run/2.1/ checked; nothing there is needed outside it, so nothing moved; the close may delete it. Summary: The phase's scratch folder was sorted file by file. Everything in it is either already recorded in the repo (the site list in the T3 summary, the cookbook drafts as entries), regenerable (giac probe batches), or a log with no lasting use. The Hilbert prototype is cited by no test or document, so it is not kept. | commit: - | tasks/T8.md | logs/T8.md
- T9 — tool candidate audit_tool_source_prefix | Done: the phase audit's tool-source check matches `tools/` and `.claude/` as root-level path segments, so `polytools/` project code is no longer flagged; promoted to `~/.claude/pa3` (3.33.1 install). Summary: The audit wrongly counted reads of our own polynomial code (a folder named `polytools`) as reads of the framework's tool source, because it matched the text `tools/` anywhere in a path. It now counts only the project's top-level `tools/` and `.claude/` folders, the install folder, the framework source, and whole `tools`/`.claude` folders outside the project. The fix lives in the PA3 package source and is installed on this machine. | commit: - | tasks/T9.md | logs/T9.md
- T10 — guard the numeric natives' endpoint, guess and constant-subtree sites with is_numeric (CAS-650 class) | Done: FSolve, CFSolve, NLPSolve, NIntegrate and NDSolve refuse a non-numeric guess/start/bound/variable-free subtree before `approximation`; five nesting tests; CI and Fuzz green on 7b96d8f. Summary: Five numeric commands turned an argument into a number by re-evaluating it from scratch, so a held call of the same command nested inside cost twice as much at every level. Each now checks the argument is a number first and otherwise leaves the call unevaluated. A 30-deep nesting test for each proves it stays fast. | commit: 7b96d8f | tasks/T10.md | logs/T10.md
- T11 — guard the numeric natives' sample-point sites by refusing a function that holds a call of the native's own head (CAS-650 class) | Done: FSolve (scalar, systems), CFSolve, NLPSolve and definite Integrate refuse a function holding a nested own-head call (depth >= 2) before any sample-point `approximation`; depth 1 still answers; 10 tests green on CI 46069824120 (test) and Fuzz 46069824257. Summary: Five numeric commands evaluate the user's function at many sample points, and a call of the same command inside it reran its whole search at every point, so deep nesting grew exponentially. Each command now leaves the call unevaluated when its own command sits nested two deep inside the function, which a quick structural check detects. One level of nesting, a legitimate use, still answers. Each command has a 30-deep nesting test and a one-level test. | commit: 99abffd | tasks/T11.md | logs/T11.md
- T12 — fix or record the random_sturmseq disagreement on sturmseq(-x^3) | Done: sturmseq places a negative leading coefficient's sign as giac does beside a power of x (CAS-661); 13 oracle rows and a unit test pin it. CI green except the pre-existing giac 2.1 mathml row. Summary: giac's random sweep found that sturmseq(-x^3) differs from giac. Probing giac showed a regular convention: it takes the power of x out before splitting the polynomial into square-free factors, so the minus sign never goes on a lone x factor (except for -x times even powers). We now do the same, so this is a fix, not a recorded departure. | commit: d0e6e3a | tasks/T12.md | logs/T12.md
- T13 — pin CI's oracle giac to 2.0.0.21+dfsg2-3+b1, the version the dev machine's probes use (C0122) | Done: CI's oracle job installs libgiac0t64 2.0.0.21+dfsg2-3+b1 from snapshot.debian.org (sha1-checked, version-asserted); CI 46082504855 green incl. mathml_documents and random_n_deriv. Summary: The test oracle on CI used whatever giac Debian's rolling release carried, and a giac upgrade there broke two tests with no change on our side. CI now installs the exact giac build the developer's machine uses for probes, fetched from Debian's permanent archive and checked by hash. If that version can't be installed, the job fails right away. | commit: 5990fb7 | tasks/T13.md | logs/T13.md

## Decisions that still bind
- T13 — CI's oracle giac moves only together with the local probe package; change GIAC_VERSION/URL/SHA1 in ci.yml and the README bullet in one commit (environment fact, in the card's Paths line).
-

## Rules proposed
- (none)

## Cookbook entries added
C0096 | CoefficientList is highest degree first; blind tests read it backwards | coefficientlist,polynomial,test,r74,order | 2026-11-27 | 2.1/T2 | CI 45995577911, T2.c9 9317acb
C0118 | Nesting tests for an approximation guard: keep the inner call held at every level | native,approximation,nesting,test,held,guard | 2026-11-27 | 2.1/T3 | T3.c3 3838e8c CAS-650
C0119 | Pre-CI sweep batch: check giac against an exact reference; lpsolve maximize and markov exactness classes | giac,oracle,sweep,proptest,lpsolve,markov,departure,reference | 2026-11-27 | 2.1/T4 | T4.c1 CI 46030056576 random_lp_objectives
C0120 | Exact poles and residue rules for rational functions (RootsWithMultiplicity, half-plane split) | residue,rational,poles,hilbert,radicals,clippy | 2026-11-27 | 2.1/T5 | T5.c1 CAS-654
C0121 | Oracle sweep where giac is wrong: probe the class with caseval, then skip on an exact predicate | proptest,oracle,sweep,giac,departure,caseval | 2026-11-27 | 2.1/T7 | CI 46055615575 random_is_irreducible
C0122 | Oracle red with no related change: check the giac version (debian:sid, unpinned) | giac,oracle,ci,debian,version,pin | 2026-11-27 | 2.1/T11 | CI 46069824120, T11.c1 99abffd
C0123 | giac sqff takes x^v out first: where the sign of a negative lead goes | giac,sqff,sturmseq,sign,caseval,factor | 2026-11-28 | 2.1/T12 | CAS-661 CI 46047327203
C0124 | Pin a Debian package version in CI from snapshot.debian.org by sha1 | debian,apt,pin,snapshot,ci,giac,oracle | 2026-11-28 | 2.1/T13 | CI 46082504855, 5990fb7

## Research
R0-001 | generation 2 frontends, cas-lang analysis phase and language-server phases (Interfaces) | cas-lang surface for a language server: what exists, what is missing | cas-lang, lsp, parser, lexer, lower, session, names, print, gaps | retriever-code | 2026-11-27 | 122 lines
R0-002 | gen 2 plan (frontends): CI jobs for a new language-server crate and a desktop app; the plotting layer phase | CI jobs, workspace shape and plotting surface for the gen 2 frontends | ci, github-actions, fuzz, workspace, cargo, lints, plotting, giac, gap-ledger | retriever-code | 2026-11-27 | 104 lines
R0-003 | gen 2 plan, scope of the plotting-layer phase and the carry-over phase (question a and b from the brief) | Gen 2 scope: plotting layer (R0-004) and 1.8 carry-over details | gen2, plotting, giac-gap, carry-over, CAS-650, markov, hilbert, fill, euler-lagrange, random-oracle | retriever-digest | 2026-11-27 | 60 lines
R0-004 | Gen 2 plotting-layer phase and hover-doc source. Question (a) plot command return values and defaults; (b) robust plotting methods; (c) giac and aide_cas license. | giac/Xcas plot commands, plotting methods, and giac/aide_cas license | giac, xcas, plotting, license, GPL-3, aide_cas, hover-docs, tupper, interval-arithmetic | retriever-web | 2026-11-27 | 64 lines
R0-005 | gen-2 plan, toolkit decision at plan approval; CI test proposal | GUI toolkit choice for the CAS desktop app (Linux first, GPL-3.0-only, CI on GitHub Actions) | gui, tauri, codemirror, monaco, egui, iced, slint, gtk, qt, dioxus, ci, webdriver, licensing | retriever-web | 2026-11-27 | 64 lines
R0-006 | generation 2 plan, design and CI tests of the language-server phases (CAS server over cas-lang and cas-core) | Rust LSP frameworks, highlighting formats, and LSP result display for the CAS language server (gen 2) | lsp, rust, tower-lsp, lsp-server, async-lsp, lsp-types, semantic-tokens, textmate, tree-sitter, zed, neovim, helix, lean, testing | retriever-web | 2026-11-27 | 73 lines
R2.1-007 | plan | CAS-650: evaluation limits, nested evaluation and the eight approximation sites | cas-650,eval,limits,approximation,native,nesting,budget,fuzz | retriever-code | 2026-11-27 | 27 lines
R2.1-008 | phase 2.1 plan, `## Interfaces` (sweeps for markov and lpsolve, Hilbert residue rule, euler_lagrange cleanup, Fill backlog ticket) | Phase 2.1 interfaces: sweeps, markov, lpsolve, Hilbert, euler_lagrange, lower.rs | cas-giac, oracle, proptest, markov, lpsolve, hilbert, residue, euler-lagrange, lower | retriever-code | 2026-11-27 | 77 lines
R2.1-009 | nesting regression tests per entry native (crates/cas-core/src/native) | Guarded-argument call paths for eight native helpers | native, nesting, guard, helper-callers, dispatch | retriever-code | 2026-11-27 | 44 lines
R2.1-010 | classify 7 solve/ call sites of numeric::approximation for 2^depth blow-up under nested Solve/NSolve | approximation() call sites: can a held same-head call reach them (CAS-650 / C0117) | CAS-650, C0117, approximation, solve, nsolve, held-call, is_numeric | retriever-code | 2026-11-27 | 58 lines
R2.1-011 | CAS-650 / C0117, classify approximation call sites in crates/cas-core/src/native/solve/ for same-head 2^depth blow-up | CAS-650 / C0117: same-head nesting exposure at numeric::approximation sites (FSolve family) | CAS-650, C0117, FSolve, CFSolve, Newton, approximation, nesting | retriever-code | 2026-11-27 | 43 lines

## Audit
- seed: median 25k, max 26k, n=13 (phase 2.1)
- previous phase 1.8: median 24k, growth 5.3%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 8436 bytes
- cookbook/INDEX.md: 22281 bytes
- rules/INDEX.md: 3674 bytes
### Carry audit — phase 2.1 (2026-11-27T11:52:45Z → open UTC, 1 sessions, 942 requests)

| file (read) | chars | n |
|---|---|---|
| decompose.rs | 28.9k | 11 |
| solve.rs | 26.8k | 16 |
| fsolve.rs | 21.7k | 4 |
| GENERATION_PLAN.md | 19.4k | 1 |
| polynomial.rs | 18.1k | 6 |
| numeric.rs | 17.0k | 9 |
| transcendental.rs | 15.5k | 12 |
| echelon.rs | 14.4k | 6 |
| matrix.rs | 11.6k | 8 |
| fuzz.yml | 11.6k | 1 |
| file (write) | chars | n |
|---|---|---|
| transform.rs | 5.6k | 4 |
| retriever-code-approximation-held-head-sites.md | 5.2k | 1 |
| T3.md | 5.2k | 1 |
| retriever-code-approx-sites-same-head-blowup.md | 4.6k | 1 |
| optimization.rs | 4.1k | 3 |
| T4.c1.md | 4.0k | 1 |
| retriever-code-nesting-guard-natives.md | 4.0k | 1 |
| T11.c1.md | 3.6k | 1 |
| result kind | chars | n |
|---|---|---|
| bash other | 883.9k | 519 |
| Grep | 472.4k | 168 |
| Read | 458.9k | 236 |
| tools/plan_edit.py | 136.9k | 29 |
| Agent | 20.5k | 20 |
| run.sh | 19.1k | 50 |
| tools/commit_task.sh | 10.3k | 25 |
| tools/task_log.py | 9.7k | 8 |
| other | 9.5k | 112 |
| tools/cookbook_add.sh | 9.4k | 19 |
- whole reads over 20.0k: 0
- seed floor: retriever-code 5,830 (n 13, prev —) · retriever-digest 5,767 (n 2, prev —) · retriever-web 4,397 (n 3, prev —)
- outline credit: 0 outlines · 0 followed by a ranged read · 0 chars credited
- warm pings: 63 pings over 8 runs, 12 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $1.1043 + relay $0.0000 + re-arm $0.0000 = $1.1043 vs rewrites replaced $6.9294, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 484, relay 0, re-arm 27, other 0, relay cost $0.0000
- retriever re-asks: 4 of 10 retriever briefs repeat a lookup of the same run
- noise: 37 lines 2.0k chars — No such file or directory: 20 lines/1.2k, usage: : 16 lines/722, warning: in the working copy of: 1 lines/66
- price: read $0.19/Mtok · 1h write $7.58/Mtok · output $18.96/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 2187 chars, 942 requests
- candidate: audit_tool_source_prefix (ratified)
- flag: 2 whole-plan reads by pa-expert

## Agent runs
- T13 | pa-claude | partner | claude-opus-5-5 | medium | ctx 152k | $3.010 | saved $11.37 | running | parent d0fd883b
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 50k | $0.651 | saved $0.28 | completed | parent d0fd883b
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 85k | $1.649 | saved $2.29 | completed | parent d0fd883b
- T2 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 80k | $0.030 | saved - | completed | parent ab5a714a
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 102k | $1.896 | saved - | completed | parent ab5a714a
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 117k | $2.415 | saved $2.98 | completed | parent d0fd883b
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 38k | $0.014 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 66k | $0.027 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 48k | $0.019 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 49k | $0.016 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 30k | $0.011 | saved - | completed | parent aac814ba
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 67k | $1.095 | saved - | completed | parent aac814ba
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 67k | $1.002 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 41k | $0.013 | saved - | completed | parent a0e4b3b8
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 53k | $0.730 | saved - | completed | parent aac814ba
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 41k | $0.449 | saved - | completed | parent aac814ba
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 19k | $0.005 | saved - | completed | parent a3d2d052
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 28k | $0.012 | saved - | completed | parent a3d2d052
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 78k | $1.282 | saved $0.08 | completed | parent d0fd883b
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 88k | $2.288 | saved - | completed | parent a8667f98
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 74k | $1.192 | saved $0.12 | completed | parent d0fd883b
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 79k | $1.546 | saved - | completed | parent a65cf55e
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 43k | $0.493 | saved - | completed | parent d0fd883b
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 66k | $1.204 | saved - | completed | parent d0fd883b
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 45k | $0.412 | saved $0.06 | completed | parent d0fd883b
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 54k | $0.755 | saved $0.08 | completed | parent d0fd883b
- T10 | pa-expert | expert | claude-opus-5-5 | medium | ctx 56k | $0.695 | saved $0.18 | completed | parent d0fd883b
- T10 | pa-coder | coder | claude-opus-5-5 | medium | ctx 70k | $0.939 | saved - | completed | parent aa66cb7f
- T11 | pa-expert | expert | claude-opus-5-5 | medium | ctx 85k | $1.566 | saved $0.50 | completed | parent d0fd883b
- T11 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 12k | $0.004 | saved - | completed | parent a600b177
- T11 | pa-coder | coder | claude-opus-5-5 | medium | ctx 107k | $2.497 | saved - | completed | parent a600b177
- T12 | pa-expert | expert | claude-opus-5-5 | medium | ctx 89k | $1.816 | saved - | completed | parent d0fd883b
- T13 | pa-expert | expert | claude-opus-5-5 | medium | ctx 55k | $0.712 | saved - | completed | parent d0fd883b

## Discussions
- I3 | Gen 2 carry-overs from 1.8 -> Gen 2: (a) CAS-650 backlog: eight functions re-evaluate an argument under fresh limits lik | deferred:gen2
- I4 | giac 2.1 reference move -> Gen 2 later phase or developer's call: debian:sid moved libgiac0t64 to 2.1.0.1+dfsg2-1 (CI 46 | deferred:2.2
- I5 | -> 2.2: 2.1 closed with three loose ends: random_n_deriv (derivatives.rs:273) fails intermittently in its near check (de | deferred:2.2

## Deferred
- from T13: CI 46082504855 is green on 5990fb7. Fuzz run 46082504901 was still queued when this task ended, so check it before T7's gate (R117). (unconsumed)
- I3 (deferred:gen2): Gen 2 carry-overs from 1.8 -> Gen 2: (a) CAS-650 backlog: eight functions re-evaluate an argument under fresh limits lik
- I4 (deferred:2.2): giac 2.1 reference move -> Gen 2 later phase or developer's call: debian:sid moved libgiac0t64 to 2.1.0.1+dfsg2-1 (CI 46
- I5 (deferred:2.2): -> 2.2: 2.1 closed with three loose ends: random_n_deriv (derivatives.rs:273) fails intermittently in its near check (de

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-27 developer: added T9 — tool candidate audit_tool_source_prefix: audit_tool_source_prefix
- 2026-11-27 router: T1 next -> done
- 2026-11-27 router: T2 next -> done
- 2026-11-27 planner: added T10 — guard the numeric natives' endpoint, guess and constant-subtree sites with is_numeric (CAS-650 class)
- 2026-11-27 planner: added T11 — guard the numeric natives' sample-point sites by refusing a function that holds a call of the native's own head (CAS-650 class)
- 2026-11-27 router: T3 split at its own trigger (29 sites, ~15 families): T3 keeps the full classification in tasks/T3.md and guards the symbolic natives (Limit, Series, Sum, Solve/CSolve, JordanMatrix, the iterative solvers, PositiveDefiniteQ, LU/QR/Cholesky, definite Integrate: 17 sites); the numeric endpoint/guess sites and the sample-point sites go to two new tasks after it.
- 2026-11-27 router: T3 next -> done
- 2026-11-27 router: T4 next -> done
- 2026-11-27 planner: added T12 — fix or record the random_sturmseq disagreement on sturmseq(-x^3) found on T5's head
- 2026-11-27 router: T5 done-when clarified: green CI except random_sturmseq, which is a separate finding (sturmseq(-x^3), CI 46047327203) moved to its own new task; T5 counts done when its own rows pass on CI.
- 2026-11-27 router: T5 next -> done
- 2026-11-27 router: T6 next -> done
- 2026-11-27 router: T7 next -> done
- 2026-11-27 router: T8 next -> done
- 2026-11-27 router: T9 next -> done
- 2026-11-27 router: T10 next -> done
- 2026-11-27 router: T11 done-when clarified: the sample-point guard refuses a user function that holds an own-head call whose arguments hold another own-head call (nesting depth >= 2); depth 1, such as FSolve[x - FSolve[y^2 - x, y], x], keeps answering. Tests: a 30-level nesting test per site (C0118) and a depth-1 held-call test.
- 2026-11-27 planner: added T13 — pin CI's oracle giac to 2.0.0.21+dfsg2-3+b1, the version the dev machine's probes use (C0122)
- 2026-11-27 router: T11 done-when clarified: green except the oracle rows broken by debian:sid's giac 2.1.0.1 upgrade (mathml.rs:377, derivatives.rs:227; CI 46069824120, C0122), which the new pin task fixes; T11 counts done with its own tests green. Moving the reference to giac 2.1 is the developer's call, parked in the inbox.
- 2026-11-27 router: T11 next -> done
- 2026-11-28 router: T12 next -> done
- 2026-11-28 router: T13 next -> done

## Plain-English Recap
Phase 2.1 Engine carry-overs: closed on the developer's yes. Final gate on 5990fb7: CI 46082504855 green, Fuzz 46082504901 completed green. CAS-650 class removed: T3 classified all 59 approximation( sites in native/ (8 CAS-650 + 51); 30 sites got an is_numeric guard (T2 8, T3 17 symbolic, T10 5 endpoint/guess/constant), T11 guards 7 sample-point sites with holds_nested (refuse own-head nesting depth >= 2; depth 1 answers); every guard has a 30-level nest test under 2 s (C0118). Sweeps random_markov_chains and random_lp_objectives (128 cases; giac's wrong 3-variable maximize is a departure row, C0119). Hilbert residue rule (CAS-654): the five rational rows agree with giac (C0120). Fixes from sweeps: sturmseq sign placement (CAS-661, C0123), random_is_irreducible skip by content in y (CAS-660, C0121). euler_lagrange_assume.rs:45 rewrite gone (CAS-659). CI's oracle giac pinned to 2.0.0.21+dfsg2-3+b1 from snapshot.debian.org with SHA-1 and dpkg-query checks (C0122, C0124), after sid's 2.1.0.1 broke two rows. Fill inside begin…end is backlog CAS-656. audit_tool_source_prefix fixed in PA3's audit.py (pa3-src 96a421a). Rules M13, M14 added. Open, in the inbox for 2.2: random_n_deriv's near check fails intermittently (two inputs recorded, unpinned); CAS-661 still to set done in Lific; Integrate's DiracDelta/piecewise paths unguarded (-> 2.6); moving to giac 2.1 is the developer's call; CAS-DOC-57/58/59 drafts; commit_task.sh still sweeps unnamed files.
