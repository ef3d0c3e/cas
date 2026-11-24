# PhaseEnd — Phase 1.4: Numerical analysis II (implements Gen 1.4)
Approved: 2026-11-23 | Closed: 2026-11-24 | Planner: claude-opus-5-5/medium | Tasks: 10 done, 0 superseded

## Milestone
Milestone: giac's quadrature rules (area with its rectangle, trapezoid, simpson, rombergt, rombergm and gauss15 methods), fsolve's twelve named solvers, newton, nSolve and convert's interval option, the sampling and interpolation commands (linspace, linear_interpolate, interval2center, center2interval, ratinterp, triginterp), the numerical derivatives nDeriv and numdiff, bspline, fitspline, fitpoly, tpsolve and minimax run under their giac names (the 41 names of ## Context), with certified balls where giac gives floats, each with an oracle test against giac or a closed form, or a recorded departure, its doc and, for a new command, a fuzz entry; CAS-DOC-56 records the area and CAS-517..520; CI and Fuzz pass on the pushed master head and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-56 and one ticket per family, recording CAS-517..520 | Done: Lific page CAS-DOC-56 "Numerical analysis II" exists (Scope, Recorded CAS-517..520, §1-§6 for T2-T7, §7 Oracle and tests); tickets CAS-588..599 exist in Numerical analysis (module 71), status todo. Summary: The phase's design page is in the tracker. It lists the 41 names this phase builds, marked built, option or alias, and the five entries left out with their reasons. It records the four earlier numerical tickets, which had no design page, and gives one design section per build task. Each of the twelve families has its own ticket, so commits can name the ticket. | commit: - | tasks/T1.md | logs/T1.md
- T2 — quadrature: area's numeric form with the rectangle, trapezoid, Simpson, Romberg and gauss15 rules | Done: `Area` (giac `area`) answers each numeric rule as a ball holding the rule's exact sum, gauss15 as NIntegrate's certified integral, and without n and method the exact integral; CI green on pushed head 0f15de6; CAS-588, CAS-589 closed. Summary: giac's `area(f,x=a..b,n,method)` now works here with all eight methods and their French spellings. Each rule's answer is a certified interval that contains the exact value of the rule's sum, so it agrees with giac's floats to their own accuracy. If the integrand has a parameter, the answer is the rule's exact sum (for example `a/2`). Where giac answers nonsense (pairs, garbage numbers, polygon output), the tests record the difference. | commit: 0f15de6 | tasks/T2.md | logs/T2.md
- T3 — solvers: fsolve's twelve named solvers, newton, nSolve and convert's interval option | Done: the 12 *_solver names and `interval` are en.txt option builtins matched by symbol; `newton` and `nSolve` are natives (giac's damped Newton, certified roots); convert(e,interval[,d]) gives balls; CI green on the pushed head. Summary: fsolve's method names used to work only by accident, as plain symbols; they are now real names of the language, so they cannot be shadowed silently. newton and nSolve port giac's Newton iteration, including its damping and restarts, checked first with a Python prototype against giac's answers. A converged answer is a ball proved to contain a root; an iteration cut short by the user's step count returns the iterate with no root claim. convert(…,interval) returns our certified balls where giac prints intervals. | commit: - | tasks/T3.md | logs/T3.md
- T4 — sampling, interpolation and numerical derivatives | Done: linspace, linear_interpolate, interval2center, center2interval, ratinterp, triginterp, numdiff, nDeriv run under their giac names, exact on exact input; CI green on 1449537 (run 45267273813); CAS-593..595 done in Lific. Summary: Eight giac commands for sampling, interpolation and numerical derivatives now work: evenly spaced points, sampling a polygonal line, the conversions between intervals and their centres, giac's rational (Floater–Hormann) and trigonometric interpolants, and finite-difference derivatives. Each answers exactly where giac prints decimals for exact data, and every giac example row is compared with libgiac on CI, with giac's error cases recorded as departures. | commit: 1449537 | tasks/T4.md | logs/T4.md
- T5 — B-splines and least-squares fits: bspline, fitspline, fitpoly | Done: bspline, fitspline and fitpoly exist as BSpline, SplineFit, PolynomialFit with unit tests, oracle rows and departures in crates/cas-giac/tests/fits.rs; CI green on 7e8c987; CAS-596, CAS-597 done. Summary: The engine now builds B-spline curves exactly from control points, fits a least-squares spline to data, and picks and fits a least-squares polynomial to data or to a function, all following giac's rules. Where giac computes in doubles on exact data, ours answers exact fractions; a fit of a function answers certified decimals. Each algorithm was first prototyped in Python against giac's answers, so the Rust matched on its first CI run. | commit: 7e8c987 | tasks/T5.md | logs/T5.md
- T6 — transportation problems: tpsolve | Done: `tpsolve` (TPSolve) answers giac's cost,plan in exact rationals; pushed; CI green on 403b003 (run 45277569828); CAS-598 done. Summary: Xcas's transportation-problem solver is now in the engine. It is a step-by-step port of giac's method: a northwest-corner first plan, then MODI improvement, with giac's exact tie-breaking, so it prints the same optimal plan giac prints. Before any Rust, a Python prototype matched giac on 165 problems, including ties, degenerate plans, forbidden routes, unbalanced totals and zero supplies. The cost is checked against the linear-program solver. | commit: 403b003 | tasks/T6.md | logs/T6.md
- T7 — minimax approximation by Remez's algorithm: minimax | Done: `MinimaxApproximation` (giac `minimax`) answers [error ball, p] or [error ball, p/q] by Remez with a certified error; CI green on pushed head 1ff9731; CAS-599 done. Summary: The engine now computes the best uniform polynomial or rational approximation of a function on an interval. It runs Remez's exchange until the error levels out to about one part in a million, where giac stops at 2%. The coefficients are exact numbers, and the error is an interval proved to contain the true maximum error of that answer. On several inputs giac's answer is not the best, and ours is visibly better (for example an error of 0.032 where giac gives 0.095); those cases are recorded as departures. | commit: 1ff9731 | tasks/T7.md | logs/T7.md
- T8 — push master for the phase's final CI and Fuzz runs | Done: master pushed at 96e80ae, level with origin/master; CI 45285948110 and Fuzz 45285948112 run on that head. Summary: The phase's finished work was checked locally (formatting, stable lints, the fuzz crate builds) and pushed to GitHub, which started the final CI and Fuzz runs. Those runs are what T9 reads; nothing may be pushed until T9 has read them, since a new push cancels Fuzz. | commit: 96e80ae | tasks/T8.md | logs/T8.md
- T9 — CI and Fuzz check on the pushed head; tickets done and CAS-DOC-56 in step | Done: CI 45317637464 and Fuzz 45317637578 green on pushed head 1a429d0; CAS-588..599 all done; CAS-DOC-56 §7 lists the README's departures. Summary: The phase's final check took five CI and Fuzz rounds. Each Fuzz run found one new problem in an older target: a false alarm in the eval harness, a memory blow-up in the evaluator's warning cache, a slow polynomial gcd, and a slow Gaussian-integer division. Each was fixed at its cause with a regression test, and the last round is all green. The tracker's twelve tickets are closed, and the design page's departure list now matches what the tests actually record. | commit: - | tasks/T9.md | logs/T9.md
- T10 — phase cleanup | Done: .run/1.4/ checked; nothing in it is still needed by the project, so the close can delete it whole. Summary: The phase's scratch folder was checked item by item. Every lesson drafted there is already in the cookbook, every code change and fuzz finding is committed with its regression test, the design page content is in the tracker, and the rest (giac probes, gap lists, giac source, ticket copies, prototypes, repro crates) can be regenerated or was a one-off check. Nothing had to be moved. | commit: - | tasks/T10.md | logs/T10.md

## Decisions that still bind
- T3 — giac's random Newton restarts become a fixed start sequence start_j = (-1)^floor((j-1)/2)·j·3/4 (j = 1..20; giac's `j > 5` → undef kept); rows whose giac answer depends on random restarts compare ours with a closed form and giac's only as a root.
- T9 — `Trace` has a public `truncated: bool`; past `MAX_STEPS` (1<<20) steps the trace is a non-replayable prefix. Route: cas-core README trace row (done) and CAS-DOC trace contract at phase end.
-

## Rules proposed
- (none)

## Cookbook entries added
C0044 | giac garbage answers change between runs: departure rows check ours only | giac,oracle,departures,nondeterministic,garbage | 2026-11-23 | 1.4/T2 | CAS-589 rombergt on a*x
C0045 | Oracle rows comparing our balls with giac's [lo..hi] intervals | giac,interval,ball,oracle,convert,overlap | 2026-11-23 | 1.4/T3 | T3.c2 CAS-592 convert interval
C0046 | Porting a float iteration to rounded Rationals: round every repeated factor | port,giac,newton,rational,rounding,iteration | 2026-11-23 | 1.4/T3 | T3.c3 CAS-591 newton
C0048 | Recurrences over Expr on symbolic inputs: cap the point count, not only the work | fuzz,timeout,limits,interpolation,symbolic | 2026-11-23 | 1.4/T4 | T4.c3 CAS-595
C0049 | giac piecewise answers in oracle rows: branch pairs, conditions by sampling | giac,oracle,piecewise,conditions,tests | 2026-11-23 | 1.4/T5 | T5.c1 CAS-596 bspline
C0050 | giac symbolic linsolve: free parameters make comparisons skip, port the skip | giac,port,linsolve,prototype,tpsolve | 2026-11-23 | 1.4/T6 | CAS-598 tpsolve zero supply
C0051 | Remez exchange: giac's one-peak-per-zero-interval step fails; use a full exchange | remez,minimax,approximation,giac,port | 2026-11-23 | 1.4/T7 | CAS-599
C0052 | Eval fuzz plain vs traced disagree at a limit: cache order, not a bug | fuzz,eval,harness,recursion-limit,cache,stamps,triage | 2026-11-24 | 1.4/T9 | Fuzz 45285948112, T9.c1 2a70d96
C0053 | Per-cache-hit side data (warnings, trace steps) needs its own count bound; OOM triage via addr2line | fuzz,oom,cache,memo,warnings,trace,limits,addr2line | 2026-11-24 | 1.4/T9 | Fuzz 45293145351, T9.c2
C0054 | Subresultant gcd: gamma/lc-scale the last remainder before its primitive part | gcd,subresultant,cas-poly,fuzz,timeout,content,performance | 2026-11-24 | 1.4/T9 | Fuzz 45304609232, T9.c5 8a82c52
C0055 | Euclid over Z[i] (and other field-backed rings): round the quotient in the ring from leading bits | gcd,gaussian,euclid,cas-num,fuzz,timeout,performance,rounding | 2026-11-24 | 1.4/T9 | Fuzz 45311134759, T9.c6 8c9ca3e

## Research
R1.4-001 | phase 1.4 plan (quadrature rules through area, fsolve solver names, newton, nSolve, convert(...,interval)) | FSolve, NIntegrate and numeric-native plumbing: what the CAS has (phase 1.4 input) | fsolve, nintegrate, quadrature, gauss-legendre, krawczyk, convert, ball-api, en.txt | retriever-code | 2026-11-23 | 61 lines
R1.4-002 | phase 1.4 plan (interpolation, splines, least-squares fitting, minimax, numerical derivatives) | giac interpolation, splines, regression, LP: what exists for phase 1.4 | interpolation, spline, thiele, pade, regression, chebyshev, lpsolve, dft, phase-1.4 | retriever-code | 2026-11-23 | 66 lines
R1.4-003 | phase 1.4 port decisions and done-when rows (bspline, fitspline, minimax, fitpoly, ratinterp, triginterp, numdiff, tpsolve) | giac 2.0.0.21 interpolation, spline, fitting, minimax, transport: port readings | giac, optimization.cc, bspline, fitspline, minimax, remez, fitpoly, ratinterp, floater-hormann, triginterp, numdiff, fdweights, tpsolve, modi | retriever-code | 2026-11-23 | 125 lines
R1.4-004 | phase 1.4 plan, port decisions and done-when rows for quadrature, root-finder and helper tasks | giac 2.0.0.21 numeric commands: area, romberg, evalf_int, newton, fsolve, nSolve, nDeriv, helpers (source reading) | giac, area, quadrature, romberg, gauss15, newton, fsolve, nSolve, nDeriv, interval, linspace, phase-1.4 | retriever-code | 2026-11-23 | 105 lines
R1.4-005 | phase 1.4 plan `## Interfaces` (giac numerical analysis commands) | Where a new giac command is added (phase 1.4 Interfaces) | giac, oracle, builtins, fuzz, docs, interfaces | retriever-code | 2026-11-23 | 66 lines

## Audit
- seed: median 24k, max 25k, n=10 (phase 1.4)
- previous phase 1.3: median 24k, growth 1.6%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 7793 bytes
- cookbook/INDEX.md: 9587 bytes
- rules/INDEX.md: 3242 bytes
### Carry audit — phase 1.4 (2026-11-23T19:22:34Z → open UTC, 1 sessions, 1286 requests)

| file (read) | chars | n |
|---|---|---|
| PHASE_PLAN.md | 69.5k | 6 |
| R1.4-003.md | 64.0k | 5 |
| sparse.rs | 61.9k | 4 |
| sections.md | 59.2k | 4 |
| fsolve.rs | 51.8k | 12 |
| gcd.rs | 49.7k | 14 |
| optimization.cc | 48.5k | 13 |
| splines.rs | 47.4k | 4 |
| risch.cc | 45.1k | 3 |
| integrate.rs | 43.5k | 7 |
| file (write) | chars | n |
|---|---|---|
| minimax.rs | 31.9k | 4 |
| sections.md | 28.8k | 1 |
| splines.rs | 28.0k | 2 |
| tables.rs | 27.7k | 1 |
| fits.rs | 20.8k | 1 |
| doc31.md | 20.0k | 1 |
| newton.rs | 13.7k | 2 |
| area.rs | 13.6k | 1 |
| result kind | chars | n |
|---|---|---|
| Read | 2.4M | 530 |
| bash other | 1.9M | 966 |
| Grep | 753.9k | 367 |
| tools/plan_edit.py | 259.9k | 17 |
| run.sh | 89.6k | 121 |
| Glob | 39.4k | 26 |
| Agent | 22.7k | 22 |
| Write | 18.8k | 115 |
| tools/commit_task.sh | 18.0k | 32 |
| other | 16.9k | 116 |
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,198 (n 29, prev 6,198) · retriever-digest 6,333 (n 4, prev 6,333) · retriever-web 4,600 (n 2, prev 4,600)
- outline credit: 7 outlines · 4 followed by a ranged read · 112.0k chars credited
- warm pings: 9 pings over 2 runs, 2 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.1877 + relay $0.0000 + re-arm $0.0000 = $0.1877 vs rewrites replaced $1.2042, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 730, relay 0, re-arm 18, other 34, relay cost $0.0000
- retriever re-asks: 0 of 2 retriever briefs repeat a lookup of the same run
- noise: 25 lines 1.5k chars — No such file or directory: 19 lines/1.3k, usage: : 6 lines/226
- price: read $0.20/Mtok · 1h write $7.98/Mtok · output $19.94/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 4292 chars, 1286 requests (prev 2822 chars, 2172 requests, growth 52.1%)
- flag: carry per request grew 52% vs 1.3
- flag: 2 whole-plan reads by pa-expert
- flag: 4 tool-source reads by main, pa-coder, retriever-code, retriever-digest
  /home/user/code/cas/.claude/settings.local.json  61 chars  main
  /home/user/.claude/projects/-home-user-code-cas/1b3aa32d-02f2-407f-9868-a62b20c7002c/tool-results/toolu_01TWmht7psJbdZCJHXw9nCri.txt  133 chars  retriever-digest
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  1.4k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  646 chars  retriever-code
- flag: 3 whole reads over 20.0k by main
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-core/src/builtin.rs  role: main  chars: 23.1k

## Agent runs
- T10 | pa-claude | partner | claude-opus-5-5 | medium | ctx 109k | $1.657 | saved $24.40 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 100k | $1.474 | saved $1.21 | completed | parent 1b3aa32d
- T1 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 82k | $0.030 | saved - | completed | parent ac7aecf5
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 113k | $2.439 | saved $5.23 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 46k | $0.595 | saved - | completed | parent a2e21df6
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 83k | $1.497 | saved - | completed | parent a2e21df6
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 114k | $2.334 | saved $4.80 | completed | parent 1b3aa32d
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 82k | $1.431 | saved - | completed | parent a4187229
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 127k | $2.385 | saved - | completed | parent a4187229
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 131k | $2.526 | saved $15.33 | completed | parent 1b3aa32d
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 95k | $1.864 | saved - | completed | parent aba93713
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 92k | $1.539 | saved - | completed | parent aba93713
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 88k | $1.506 | saved - | completed | parent aba93713
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 138k | $2.694 | saved $11.33 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 168k | $3.342 | saved - | completed | parent a8ed1345
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 178k | $3.820 | saved - | completed | parent a8ed1345
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 107k | $2.175 | saved $0.19 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 85k | $1.734 | saved - | completed | parent a99dd9e1
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 102k | $1.734 | saved $16.97 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 179k | $4.465 | saved - | completed | parent a0d21045
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 107k | $2.266 | saved - | completed | parent a0d21045
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 42k | $0.365 | saved - | completed | parent 1b3aa32d
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 132k | $3.525 | saved $15.46 | completed | parent 1b3aa32d
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 77k | $1.208 | saved - | completed | parent aa0bcbd8
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 99k | $1.740 | saved - | completed | parent aa0bcbd8
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 73k | $1.178 | saved - | completed | parent aa0bcbd8
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 83k | $1.408 | saved - | completed | parent aa0bcbd8
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 87k | $1.532 | saved - | completed | parent aa0bcbd8
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 10k | $0.002 | saved - | completed | parent ac548e50
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 97k | $1.876 | saved - | completed | parent aa0bcbd8
- T10 | pa-expert | expert | claude-opus-5-5 | medium | ctx 51k | $0.513 | saved - | completed | parent 1b3aa32d

## Discussions
- (none)

## Deferred
- from T10: the close may delete .run/1.4/ whole. (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-23 router: T1 next -> done
- 2026-11-23 router: T2 next -> done
- 2026-11-23 router: T3 next -> done
- 2026-11-23 router: T4 next -> done
- 2026-11-23 router: T5 next -> done
- 2026-11-23 router: T6 next -> done
- 2026-11-24 router: T7 next -> done
- 2026-11-24 router: T8 next -> done
- 2026-11-24 router: T9 next -> done
- 2026-11-24 router: T10 next -> done

## Plain-English Recap
Phase 1.4 Numerical analysis II: CAS-DOC-56 sorts 41 giac names and records CAS-517..520; tickets CAS-588..599 are done. Built: area with the rectangle, trapezoid, Simpson, Romberg and gauss15 rules (balls holding each rule's exact sum); the twelve *_solver names as option builtins, newton (Krawczyk-certified roots, fixed restart sequence instead of giac's random one), nSolve and convert's interval; linspace, linear_interpolate, interval2center, center2interval, ratinterp (faithful port), triginterp, numdiff (Lagrange-derivative weights), nDeriv; bspline, fitspline and fitpoly, exact on exact data; tpsolve (northwest corner plus MODI with giac's tie rules, matching giac on 165 problems); minimax by Remez with a full exchange step instead of giac's (which fails on abs at degree 2) and a certified error ball. Exact answers where giac prints floats, compared within giac's error. Final gate green on 1a429d0 (CI 45317637464, Fuzz 45317637578) after four fuzz findings in older code, each fixed: a harness false alarm on recursion-limit disagreement, OOM from unbounded cache-hit warnings and the trace (capped; the trace is marked truncated), a slow subresultant gcd (Brown normalisation) and slow Gaussian-integer division. Cookbook C0044-C0055. Open: R33 (coders in T3 and T7 again used local scratch binaries); a CoefficientList bug that swaps terms with ball coefficients (no ticket yet); commit_task.sh stages every changed file under logs/; docs-only commits unpushed until the next phase's first push.
