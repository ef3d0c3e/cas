# PhaseEnd — Phase 1.2: Algebra and polynomials parity (implements Gen 1.2)
Approved: 2026-11-22 | Closed: 2026-11-23 | Planner: claude-opus-5-5/medium | Tasks: 12 done, 0 superseded

## Milestone
Milestone: The 89 giac Algebra and Polynomials entries this CAS lacked (R1.2-006) are accounted for: each command, option or setting this plan builds runs under its giac names with a giac oracle test or a recorded departure, its doc chapter and, for a new builtin, a fuzz entry, and the inert forms, giac's Gröbner engine knobs and increasing_power are recorded departures; CI and Fuzz pass on the pushed master head, CAS-DOC-54 records the area and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-54 and one ticket per family | Done: Lific page CAS-DOC-54 "Algebra and polynomials parity" (Scope, §1-§8 for T2-T9, §9 Oracle and tests) and 21 tickets CAS-526..546, one per family of R1.2-006 (F2 and F19 split by module), all status todo. Summary: The design page lists all 89 missing giac entries by family, each marked built, alias, option, setting or departure, with one section per build task giving the names, what they build on, the algorithm and the oracle plan. Each family has a ticket in its Lific module, so every build task commits against a ticket. giac 2.0.0 was probed for each family's example; the answers are quoted on the page. | commit: - | tasks/T1.md | logs/T1.md
- T2 — gcd by named algorithm (heugcd, modgcd, psrgcd, ezgcd) and icontent; polynomial inert forms as departures | Done: heugcd, modgcd, psrgcd, ezgcd, icontent built (new EZ-GCD method in cas-poly); Gcd, Factor, Quo, Rem, Resultant are departure rows; pushed, CI green on 6038ef7. Summary: The four giac gcd commands now run a chosen algorithm: the heuristic, Brown's modular, subresultants, and a new EZ-GCD that evaluates all variables but one, takes a one-variable gcd and lifts it back by Hensel lifting shared with multivariate factoring. icontent gives the gcd of the integer coefficients. Where giac fails (its EZ-GCD gives up or refuses one-variable input) or signs differently (psrgcd), ours answers and the difference is recorded. Capitalised inert forms stay unevaluated, as recorded departures. | commit: 6038ef7 | tasks/T2.md | logs/T2.md
- T3 — partial fractions over C, proper fractions and convert's options | Done: cpartfrac (ComplexApart), propfrac/propFrac (ProperFraction), f2nd/fxnd (NumeratorDenominator), ratnormal alias, convert's partfrac/parfrac/fullparfrac/confrac/polynom; pushed, CI green (run 45060121388). Summary: Partial fractions now also work over the complex numbers, splitting the denominator into linear factors with i and square roots as giac does. convert learned the giac options for partial fractions, continued fractions and turning a series into its polynomial; propfrac, f2nd/fxnd and ratnormal were added. Each has giac comparisons or recorded departures, and CI passed on the pushed head. | commit: bdbedf4 | tasks/T3.md | logs/T3.md
- T4 — real root isolation (VAS), root bounds, Sturm sequences, rational and complex root lists | Done: VAS, VAS_positive, posubLMQ, poslbdLMQ, sturmseq/sturm, rationalroot, crationalroot, complexroot built, matching giac's examples; oracle file real_roots.rs green on CI; CAS-538, CAS-539 closed. Summary: Eight giac root names now work: real root isolation by continued fractions with giac's own bounds, giac's Sturm sequences and counts, rational and Gaussian-rational root lists, and complex root enclosures. The algorithms were ported from giac's source, so the intervals and sequences match giac exactly rather than only being valid. A new chapter documents Sturm sequences, VAS and the LMQ bounds. | commit: 382157d | tasks/T4.md | logs/T4.md
- T5 — complex numeric solving, zero lists, surd, and the all_trig_solutions and with_sqrt settings | Done: cfsolve (certified complex roots), zeros/cZeros, surd, and the all_trig_solutions and with_sqrt settings are built, oracle-tested and green on CI; the nlpsolve oracle flake is fixed. Summary: Equations can now be solved numerically over the complex numbers with each root proved, polynomials directly and other equations by Newton's method checked by Krawczyk's test. zeros and cZeros list the solutions, surd gives real odd roots, and two giac settings work: all_trig_solutions gives trig equations' general solutions with an integer parameter, and with_sqrt(0) keeps factor over the rationals. A pre-existing oracle test that failed by chance (giac's nlpsolve depends on random draws from earlier tests) now seeds giac first. | commit: 40cc289 | tasks/T5.md | logs/T5.md
- T6 — polynomial tools: rootof, fcoeff, homogeneize, total_degree, factor_xn, bezoutian, is_irreducible | Done: giac's rootof, fcoeff, homogeneize, total_degree, factor_xn, bezoutian, is_irreducible run under their names, with unit tests, oracle rows and departure rows; CI green on the pushed head. Summary: Seven polynomial tools from giac now exist in the engine, each checked against giac's own answers. rootof stays a symbolic call as in giac, since the engine has no algebraic numbers; where giac computes with it (normal, evalf) that is recorded as a departure. Where giac answers wrongly or the engine cannot read the input, the difference is written down with its reason. | commit: 8ea8c21 | tasks/T6.md | logs/T6.md
- T7 — Gröbner: in_ideal, syzygy, gbasis_vars and gbasis's rur option; giac's engine knobs as departures | Done: in_ideal, syzygy, gbasis_vars and gbasis(…, rur) built with unit, oracle and departure tests; the five engine knobs are departure rows; CI green on 395ec60. Summary: The Gröbner tools now test ideal membership, compute the syzygies of a list of polynomials, choose a variable order as giac does, and give the rational univariate representation of a zero-dimensional system: one univariate polynomial whose roots parametrise every solution. The representation reproduces giac's answer whenever giac's separating form is a single variable or a plain sum. Where giac is wrong or arbitrary (syzygy with explicit variables, its incomplete syzygy lists, its random separating forms), ours gives the correct answer and the difference is a recorded departure. giac's tuning knobs for its modular F4 engine do not apply to our Buchberger engine and stay as calls. | commit: 395ec60 | tasks/T7.md | logs/T7.md
- T8 — matrix constructors, shape, slicing, row operations, norms and subspaces | Done: the 28 giac matrix names of F5-F8 (with all aliases) evaluate; oracle file matrices_three.rs; CI green on 63bce19; CAS-531..534 closed. Summary: Xcas's matrix builders, shape queries, row/column slicing and deletion, elementary row operations, max/row/column norms, row and column spaces, the intersection of two spans and the coefficient matrix of a linear system now work under their giac names. Indices are 0-based through the shared index rule, so negatives count from the end. Where giac's answers are only one of many equal bases, ours give the reduced echelon basis and the oracle checks the spans agree. | commit: - | tasks/T8.md | logs/T8.md
- T9 — elimination and canonical forms: det_minor, adjoint_matrix, ref, pivot, keep_pivot, egvl, gauss, rat_jordan; matrix inert forms as departures | Done: det_minor, adjoint_matrix, ref, pivot, rref keep_pivot, egvl/eigVl, gauss, rat_jordan built with oracle rows; Det, Inverse, Nullspace, Rref departure rows; CI green on 425337b (run 45085825886). Summary: Eight giac matrix commands now work: determinants by minors, the adjugate of xI−A as Faddeev–Le Verrier gives it, row echelon forms and single pivot steps exactly as giac chooses its pivots, the Jordan matrix in giac's block order, quadratic forms written as sums of squares, and the rational Jordan form. giac's pivot rule was ported from its source and checked with a prototype before the Rust. The capitalised inert forms stay unevaluated by design and are recorded as departures. | commit: 425337b | tasks/T9.md | logs/T9.md
- T10 — push master for the phase's final CI and Fuzz runs | Done: master pushed at 1c42f01 (level with origin); CI run 45087081264 and Fuzz run 45087081385 started on it. Summary: The local checks (formatting, stable clippy, fuzz crate build) passed on the committed head, and master was pushed so GitHub runs the phase's final CI and Fuzz. Both runs are queued on the same commit; T11 reads their results. | commit: 1c42f01 | tasks/T10.md | logs/T10.md
- T11 — CI and Fuzz check on the pushed head; tickets done and CAS-DOC-54 in step | Done: master 5680a2c pushed with CI 45091838361 and Fuzz 45091838357 green; CAS-526..547 all done; CAS-DOC-54 §9 lists every README departure by oracle file. Summary: The first Fuzz run of the phase to finish (on T10's push) found a timeout in the rewrite target: raising a polynomial to a huge power passed both size limits separately while the result would have needed about 800 MB. A coder bounded the product of the two limits, added a regression test, and pushed; CI and the 40-minute Fuzz run are now green. All of the phase's tickets are closed, and the design page now records every place where our answers knowingly differ from giac's, matching the oracle README. | commit: - | tasks/T11.md | logs/T11.md
- T12 — phase cleanup | Done: .run/1.2/ checked item by item; nothing in it is needed outside it, so nothing moved; the close may delete it. Summary: The phase's scratch folder was checked for anything the project still needs. Everything there is either regenerable (gap lists), already recorded (giac answers in the tests and task logs, coder briefs in the logs), already in its real home (the design pages in Lific, the fuzz finding as a regression test), or throwaway prototypes. Nothing had to move into the repo. | commit: - | tasks/T12.md | logs/T12.md

## Decisions that still bind
-

## Rules proposed
- T4 — (none)

## Cookbook entries added
C0007 | giac error answers: caseval prints GIAC_ERROR: <msg>, the CLI a quoted string | giac,oracle,departures,errors,caseval | 2026-11-22 | 1.2/T2 | CAS-526 ezgcd departures
C0008 | Option names that are builtin aliases lower to the builtin's symbol | giac,options,aliases,names,lowering,convert | 2026-11-22 | 1.2/T3 | CAS-529 convert options
C0009 | Matching giac's exact output: port its source, prove with a Python prototype | giac,oracle,source,port,prototype,debian | 2026-11-22 | 1.2/T4 | CAS-538/539 VAS, sturmseq
C0010 | Numeric tests on complex balls: compare Re and Im, never Abs | ball,complex,abs,indeterminate,oracle,numeric,tests | 2026-11-22 | 1.2/T4 | CAS-538 complexroot CI 45065408627
C0011 | Adding a session setting on the Digits precedent | settings,session,digits,lower,oracle | 2026-11-22 | 1.2/T5 | T5.c1 all_trig_solutions, with_sqrt
C0012 | Certified complex roots by real Krawczyk on the 2n real system | krawczyk,complex,roots,newton,cfsolve,certified | 2026-11-22 | 1.2/T5 | T5.c2 cfsolve
C0013 | Comparing a complex ball with a close value by parts, never by Abs | ball,complex,oracle,comparison,numeric | 2026-11-22 | 1.2/T5 | T5.c3 complex_solving tests
C0014 | giac answers that depend on random state: seed inside the same input | giac,oracle,flaky,random,nlpsolve,ci | 2026-11-22 | 1.2/T5 | CI 45070842403 nonlinear_programs
C0015 | giac CLI and libgiac answers differ: copy departure rows from the oracle's own output | giac,oracle,departures,caseval,libgiac,bezoutian | 2026-11-22 | 1.2/T6 | T6.c1 rootof, bezoutian rows
C0016 | Natives taking equations match Equation and Equal | equation,native,lowering,cas-core | 2026-11-23 | 1.2/T8 | syst2mat oracle failure
C0017 | giac row reductions: RREF_GUESS, islesscomplexthan pivots, _ref in ti89.cc | giac,rref,ref,pivot,keep_pivot,bareiss,egvl,source | 2026-11-23 | 1.2/T9 | T9.c1 ref, keep_pivot
C0018 | Exact oracle text from libgiac via Python ctypes caseval (no cargo test) | giac,oracle,caseval,ctypes,libgiac,departures | 2026-11-23 | 1.2/T9 | T9.c2 rat_jordan departure rows
C0019 | Fuzz timeout with no limit hit: bound the product of limits | fuzz,timeout,limits,polynomial,power | 2026-11-23 | 1.2/T11 | CAS-547

## Research
R1.2-001 | phase 1.2 `## Interfaces` (about 100 giac Algebra and Polynomials commands), settlement of giac's capitalised inert forms against our names | Builtin and alias touch points, re-verified, with the five open points | cas-lang, names, en.txt, builtin.rs, lowering, cas-giac, cas-poly, solve, touch-points | retriever-code | 2026-11-22 | 91 lines
R1.2-002 | phase 1.2 (Algebra and polynomials parity): for each giac name, the existing internal builtin, its implementation, and what is missing. | giac leftovers for phase 1.2: what exists to build on | giac, phase-1.2, gcd, partfrac, det, solve, modular, continued-fraction | retriever-code | 2026-11-22 | 66 lines
R1.2-003 | phase 1.2 plan, giac leftovers (gbasis/rur, polytools, Hensel/EZ-GCD, cpartfrac) | Phase 1.2 leftovers: what the CAS already has (Groebner, poly tools, Hensel/gcd, partial fractions) | groebner, gbasis, rur, polytools, factor, hensel, gcd, apart, cfactor | retriever-code | 2026-11-22 | 94 lines
R1.2-004 | giac matrix leftovers in phase 1.2 (index convention, matrix helpers, tests, docs, fuzz, per-name build-on/new) | Matrix leftovers for phase 1.2: what exists to build on | matrix, index, giac, phase-1.2, interfaces | retriever-code | 2026-11-22 | 133 lines
R1.2-005 | phase 1.2 leftovers from giac: VAS/VAS_positive, posubLMQ/poslbdLMQ, sturm/sturmseq, proot/realroot/complexroot/cfsolve/cZeros, all_trig_solutions, with_sqrt, surd, Digits precedent | Roots, solving and settings leftovers for phase 1.2: what the CAS has | roots, sturm, polynomial-roots, cas-num, fsolve, krawczyk, solve, trig, dsolve-constants, digits, settings, surd | retriever-code | 2026-11-22 | 135 lines
R1.2-006 | list the giac help entries R0-004 puts in Algebra/Polynomials (and near-name flags), plus missed ones, with kind and family, for the phase 1.2 "Algebra and polynomials parity" plan | Algebra and Polynomials leftovers: giac entries, kinds, families (phase 1.2 planning) | giac, gap, algebra, polynomials, parity, phase-1.2 | retriever-digest | 2026-11-22 | 156 lines
R1.2-007 | giac sturmseq remainder scaling, VAS/VAS_positive, posubLMQ/poslbdLMQ source reading | giac sturmseq / VAS / posubLMQ / poslbdLMQ implementation (partial, low confidence) | giac, sturm, vas, lmq, root-isolation, partial | retriever-web | 2026-11-22 | 41 lines

## Audit
- seed: median 23k, max 23k, n=12 (phase 1.2)
- previous phase 1.1: median 22k, growth 3.9%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 7545 bytes
- cookbook/INDEX.md: 3363 bytes
- rules/INDEX.md: 3096 bytes
### Carry audit — phase 1.2 (2026-11-22T18:37:14Z → open UTC, 1 sessions, 1358 requests)

| file (read) | chars | n |
|---|---|---|
| gcd.rs | 76.9k | 20 |
| groebner.rs | 69.1k | 6 |
| groebner.rs | 59.6k | 7 |
| matrix.rs | 51.5k | 12 |
| sparse.rs | 45.2k | 2 |
| csturm.cc | 44.9k | 7 |
| transcendental.rs | 43.6k | 4 |
| R1.2-006.md | 39.5k | 2 |
| polytools.rs | 38.9k | 6 |
| decompose.rs | 37.3k | 8 |
| file (write) | chars | n |
|---|---|---|
| doc31.md | 20.0k | 1 |
| echelon.rs | 19.9k | 6 |
| roots.rs | 18.1k | 1 |
| zero_dim.rs | 17.6k | 1 |
| page.md | 17.4k | 1 |
| real_roots.rs | 17.3k | 4 |
| decompose.rs | 14.3k | 5 |
| doc34.md | 12.2k | 1 |
| result kind | chars | n |
|---|---|---|
| Read | 1.7M | 361 |
| bash other | 1.7M | 993 |
| Grep | 349.2k | 145 |
| tools/plan_edit.py | 229.8k | 17 |
| run.sh | 81.4k | 103 |
| other | 67.0k | 127 |
| Agent | 18.8k | 18 |
| Write | 18.4k | 112 |
| tools/outline.py | 17.5k | 10 |
| tools/card.py | 13.2k | 3 |
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,257 (n 10, prev 6,257) · retriever-digest 6,312 (n 2, prev 6,312) · retriever-web 4,600 (n 2, prev 4,600)
- outline credit: 12 outlines · 3 followed by a ranged read · 89.5k chars credited
- warm pings: 8 pings over 2 runs, 2 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.1712 + relay $0.0000 + re-arm $0.0000 = $0.1712 vs rewrites replaced $1.2998, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 654, relay 0, re-arm 15, other 34, relay cost $0.0000
- retriever re-asks: 0 of 2 retriever briefs repeat a lookup of the same run
- noise: 18 lines 1.2k chars — No such file or directory: 14 lines/956, usage: : 4 lines/212
- price: read $0.20/Mtok · 1h write $7.97/Mtok · output $19.94/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 3093 chars, 1358 requests (prev 15978 chars, 136 requests, growth -80.6%)
- flag: 2 whole-plan reads by pa-expert
- flag: 8 tool-source reads by main, pa-coder, retriever-code, retriever-digest
  /home/user/code/cas/.claude/settings.local.json  61 chars  main
  /home/user/.claude/projects/-home-user-code-cas/1b3aa32d-02f2-407f-9868-a62b20c7002c/tool-results/toolu_01TWmht7psJbdZCJHXw9nCri.txt  133 chars  retriever-digest
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  5.6k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/tests.rs  2.1k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/roots.rs  4.9k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  1.4k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  1.4k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  646 chars  retriever-code
- flag: 3 whole reads over 20.0k by main
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-core/src/builtin.rs  role: main  chars: 23.1k

## Agent runs
- T12 | pa-claude | partner | claude-opus-5-5 | medium | ctx 119k | $2.638 | saved $30.58 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 92k | $1.424 | saved - | completed | parent 1b3aa32d
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 105k | $2.080 | saved $8.49 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 109k | $1.900 | saved $0.18 | completed | parent a71f1fc4
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 68k | $1.071 | saved - | completed | parent a71f1fc4
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 90k | $1.647 | saved $2.48 | completed | parent 1b3aa32d
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 96k | $1.739 | saved - | completed | parent a7ef80be
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 132k | $2.990 | saved $24.87 | completed | parent 1b3aa32d
- T4 | retriever-web | retriever | claude-haiku-5-5 | medium | ctx 27k | $0.010 | saved - | completed | parent ab3e777c
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 274k | $7.787 | saved - | completed | parent ab3e777c
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 66k | $0.946 | saved - | completed | parent ab3e777c
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 144k | $3.893 | saved $11.82 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 113k | $2.154 | saved - | completed | parent a826498e
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 143k | $2.944 | saved - | completed | parent a826498e
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 66k | $0.931 | saved - | completed | parent a826498e
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 89k | $1.650 | saved $0.51 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 113k | $2.270 | saved - | completed | parent a59179b8
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 138k | $2.676 | saved $12.48 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 95k | $1.704 | saved - | completed | parent a8f23a7b
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 94k | $1.807 | saved - | completed | parent a8f23a7b
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 120k | $2.526 | saved $9.55 | completed | parent 1b3aa32d
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 113k | $2.228 | saved - | completed | parent a62e0344
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 87k | $1.463 | saved - | completed | parent a62e0344
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 149k | $3.491 | saved $6.01 | completed | parent 1b3aa32d
- T9 | retriever-web | retriever | claude-haiku-5-5 | medium | ctx 30k | $0.012 | saved - | completed | parent a69d9d2f
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 166k | $4.013 | saved - | completed | parent a69d9d2f
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 157k | $3.433 | saved - | completed | parent a69d9d2f
- T10 | pa-expert | expert | claude-opus-5-5 | medium | ctx 38k | $0.318 | saved - | completed | parent 1b3aa32d
- T11 | pa-expert | expert | claude-opus-5-5 | medium | ctx 81k | $1.336 | saved $0.34 | completed | parent 1b3aa32d
- T11 | pa-coder | coder | claude-opus-5-5 | medium | ctx 81k | $1.542 | saved - | completed | parent a62b99ae
- T12 | pa-expert | expert | claude-opus-5-5 | medium | ctx 43k | $0.418 | saved - | completed | parent 1b3aa32d

## Discussions
- (none)

## Deferred
- from T12: the close deletes .run/1.2/; leave .run/logs/gen1-gap/ (GENERATION_PLAN.md:57). (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-22 router: T1 next -> done
- 2026-11-22 router: T2 next -> done
- 2026-11-22 router: T3 next -> done
- 2026-11-22 router: T4 next -> done
- 2026-11-22 router: T5 next -> done
- 2026-11-22 router: T6 next -> done
- 2026-11-22 router: T7 next -> done
- 2026-11-23 router: T8 next -> done
- 2026-11-23 router: T9 next -> done
- 2026-11-23 router: T10 next -> done
- 2026-11-23 router: T11 next -> done
- 2026-11-23 router: T12 next -> done

## Plain-English Recap
Phase 1.2 Algebra and polynomials parity: CAS-DOC-54 sorts 89 giac entries (66 built, 1 alias, 5 options, 2 settings, 15 recorded departures); tickets CAS-526..547 done. Built: named gcds (heugcd, modgcd, psrgcd, new EZ-GCD) and icontent; real root isolation (VAS), LMQ bounds, Sturm sequences, rational/complex root lists (ported from giac source, checked by Python prototype); cfsolve with certified complex roots, zeros, cZeros, surd, all_trig_solutions and with_sqrt settings; cpartfrac, propfrac, f2nd, convert options; rootof, fcoeff, homogeneize, total_degree, factor_xn, bezoutian, is_irreducible; in_ideal, syzygy, gbasis_vars, gbasis rur; 28 matrix names plus det_minor, adjoint_matrix, ref, pivot, keep_pivot, egvl, gauss, rat_jordan. Final gate green on 5680a2c (CI 45091838361, Fuzz 45091838357) after a fuzz timeout on (1+x^2/4)^46656 led to a combined size limit on polynomial powers (CAS-547). The flaky nlpsolve oracle row is fixed by seeding giac. Cookbook C0007-C0019 added. Open: the developer's ruling on local throwaway example programs; a cas-num fix for Abs of a complex ball containing 0; six log/doc commits unpushed until the next phase's first push.
