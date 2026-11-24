# PhaseEnd — Phase 1.5: Optimization II (implements Gen 1.5)
Approved: 2026-11-24 | Closed: 2026-11-24 | Planner: claude-opus-5-5/max | Tasks: 10 done, 0 superseded

## Milestone
Milestone: lpsolve and nlpsolve take giac's options (the 46 option words of ## Context: the lp entry's 31 lp_* names and nonnegint, the nlp entry's 14 nlp_* names) and simplex_reduce, find_minimum, frank_wolfe and simulated_annealing run under their giac names, answering by ## Context's Numbers rule (exact where the problem is rational, certified balls for inexact optima), each with an oracle test against giac or known optima, or a recorded departure, its doc and a fuzz entry; CAS-DOC-57 records the area; CI and Fuzz pass on the pushed master head and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-57 and one ticket per family | Done: Lific page CAS-DOC-57 "Optimization II" (draft) holds the phase design; tickets CAS-600..610 (status todo) exist, one per family, each in its module. Summary: The phase's design is now on the tracker: the 50 missing optimization names sorted by task and family, what the existing linear and nonlinear solvers already do and what each option changes in them, one section per build task with giac's algorithm and examples, and the departures and known optima the tests use. Eleven tickets carry the work, one per family. | commit: - | tasks/T1.md | logs/T1.md
- T2 — lpsolve's options: the lp_* names, maximize and assume forms, mixed integer variables, lp_method and decimal data | Done: LPSolve takes giac's lp_* options (31 option builtins, nonnegint alias), mixed integer/binary variables, every lp_method, decimal data exactly; CI green on 611c165; CAS-600, CAS-601 done. Summary: lpsolve now understands giac's option words as real names rather than loose symbols: maximize and minimize in all giac's spellings, the nonnegative, integer and binary assumptions, lists of integer or binary variables for mixed problems, and the method, verbose and presolve switches (which change no answer, since ours is always exact). Decimal coefficients are turned into the fractions giac's exact gives, so answers stay exact as giac's are. The branch-and-bound strategy and limit options still leave the call until T3 gives them meaning. | commit: 611c165 | tasks/T2.md | logs/T2.md
- T3 — lpsolve's branch and bound: giac's node and variable selection, cuts, heuristic, limits and gap | Done: LPSolve's integer search is giac's branch and bound (lpsolve.cc:1599-1931), exact, with every lp_nodeselect/lp_varselect strategy, GMI cuts, rounding heuristic, limits and gap; CI green on 6e4980d; CAS-602, CAS-603 done. Summary: Integer linear programs are now solved the way giac solves them: it picks which open branch to explore and which variable to split as giac does, tightens each subproblem with Gomory cuts, and tries rounding to find good integer points early. Every strategy option giac offers works and all give the same optimum, checked against brute force. The node, depth and iteration limits and the gap tolerance stop the search as giac's do, answering the best point found or the empty list; the time limit is accepted but bounds nothing, as we have no clock. | commit: 6e4980d | tasks/T3.md | logs/T3.md
- T4 — nlpsolve's options: the nlp_* names, initial points, limits, precision, methods, integer and binary variables | Done: NLPSolve reads the 14 nlp_* option builtins by symbol with giac's forms, sends linear programs to LPSolve (exact), and solves integer/binary nonlinear programs (enumeration, else B&B); CI green on 67bbf0c; CAS-604, CAS-605 done. Summary: nlpsolve now understands every option giac documents: starting points and regions, iteration limits, precision, tolerance, the method names, maximize and nonnegativity, and the presolve and verbose switches. Linear problems are handed to the exact linear solver, as giac does. Problems with integer or binary variables are solved by trying every integer point when the box is small, and otherwise by branch and bound on certified relaxations; exact answers come back wherever the objective is exact at integer points. Most integer problems where giac itself fails are now tested against their known optima. | commit: 67bbf0c | tasks/T4.md | logs/T4.md
- T5 — the simplex tableau: simplex_reduce, and CoefficientList's order on ball coefficients | Done: SimplexReduce (simplex_reduce, CAS-606) answers giac's sequence on exact and decimal tableaux; CoefficientList on decimal coefficients checked leading-first (CAS-607); CI green on 720e69d; both tickets done. Summary: giac's teaching-style simplex, which returns the optimum, the solution and the final tableau, is now a command, ported line by line from giac and checked against all nine of giac's probe answers first in a Python prototype. A decimal in the input leaves the other columns exact, as giac does. When the right side has a negative entry, giac returns an infeasible answer; ours leaves the call with a warning. The suspected CoefficientList term swap from phase 1.4 did not reproduce: the new test passes on CI with the leading coefficient first. | commit: 720e69d | tasks/T5.md | logs/T5.md
- T6 — univariate and polytope minimisation: find_minimum and frank_wolfe | Done: `find_minimum` answers a ball certified to hold a minimiser of f on [a,b]; `frank_wolfe` answers the iterate whose Frank–Wolfe gap is certified below tol; both registered, documented, unit and oracle tested; CI green on 4715ab0. Summary: Two giac minimisers now exist. find_minimum searches like giac (Brent over a queue of intervals) but then proves its answer with ball arithmetic, so the answer is an enclosure of the true minimiser rather than a float. frank_wolfe minimises a function over the convex hull of a list of points like giac does, but only stops once the optimality gap is proven below the tolerance, which bounds how far its value is from the minimum for convex f. | commit: 4715ab0 | tasks/T6.md | logs/T6.md
- T7 — simulated annealing: simulated_annealing | Done: `simulated_annealing` (SimulatedAnnealing) ports GSL's gsl_siman_solve as giac calls it, applying the user's cost and step functions through the evaluator with seeded draws; CI green on 77f1ddb; CAS-610 done. Summary: The last optimiser of the phase is in. It runs the classic annealing loop on whatever configuration the user's step function produces, keeps the best point it visits, and repeats exactly under a fixed seed. Its first CI run showed that feeding inexact results back into the user's step made their error bars grow until the search froze; the steps now carry exact midpoints, like the doubles GSL uses, and the help-page example lands at the known minimum under every seed tested. | commit: 77f1ddb | tasks/T7.md | logs/T7.md
- T8 — push master for the phase's final CI and Fuzz runs | Done: master pushed at a2ba5de with the local gate green; CI run 45350578782 and Fuzz run 45350578815 started on that head. Summary: The phase's work was checked locally (formatting, stable lints, the fuzz crate builds) and pushed to GitHub, which started the final CI and the 40-minute fuzzing run on the same commit. T9 reads both results. | commit: a2ba5de | tasks/T8.md | logs/T8.md
- T9 — CI and Fuzz check on the pushed head; tickets done and CAS-DOC-57 in step | Done: CI and Fuzz are green on master's pushed head a2ba5de; CAS-600..610 are done; CAS-DOC-57 §7 lists the departures the phase recorded. Summary: The phase's final push passed every CI job and every fuzz target, so nothing needed fixing. All eleven phase tickets are closed in the tracker. The design page's test section now lists the departures from giac that the tests actually record, not the ones the plan expected, and names the planned ones that turned out to match giac. | commit: - | tasks/T9.md | logs/T9.md
- T10 — phase cleanup | Done: .run/1.5/ holds nothing the repo still needs; committed files no longer point into it. Summary: The phase's scratch folder was checked item by item. Every drafted lesson is already in the cookbook, ticket and design content is in the tracker, and the rest (giac probes, gap lists, giac source, prototypes, CI logs) can be regenerated or was a one-off check. Five comments in tests and cookbook entries named scratch prototype files that the close will delete; they now name the prototype by task instead. | commit: - | tasks/T10.md | logs/T10.md

## Decisions that still bind
-

## Rules proposed
- (none)

## Cookbook entries added
C0056 | Stays tests check head and warnings; feasibility checks split equalities | test,unevaluated,oracle,equation,cas-core | 2026-11-24 | 1.5/T2 | T2.c1 CAS-600 lpsolve stays
C0057 | giac branch and bound port: GMI cuts by slack substitution and re-solve; giac's cut rhs bug | giac,port,lpsolve,integer,cuts,branch-and-bound | 2026-11-24 | 1.5/T3 | T3.c1 CAS-602
C0058 | Expr::may_contain is a bloom filter: never assert its negation | tests,expr,symbols,builtins,flaky | 2026-11-24 | 1.5/T4 | CI 45337180530 trigonometric_interpolation
C0059 | Certified KKT points under bounds: one-point ranges; simple inequalities as bounds | nlpsolve,kkt,certified,bounds,branch-and-bound | 2026-11-24 | 1.5/T4 | T4.c2 CAS-605; CI 45337180530 nlpsolve_seeded
C0061 | Certified global minimiser on an interval: float candidate, ball cover, interval Newton | minimum,certified,krawczyk,cover,brent,find_minimum | 2026-11-24 | 1.5/T6 | T6.c1 CAS-608
C0062 | Iterative float optimisers with a certified stop: test in balls, one shared work budget | frank-wolfe,certified,optimisation,limits,budget,gap | 2026-11-24 | 1.5/T6 | T6.c2 CAS-609
C0063 | A native applying user functions: Expr::apply + ev.eval, N[] for numbers | native,function,apply,eval,random,seed | 2026-11-24 | 1.5/T7 | T7.c1 CAS-610
C0064 | Iterated user steps on balls: exact midpoints, or radii compound | ball,radius,random,iteration,native,annealing | 2026-11-24 | 1.5/T7 | T7.c2 CAS-610; CI 45347780940

## Research
R1.5-001 | Optimization II port decisions and sizes (phase 1.5), read-only giac source in /home/user/code/cas/.run/1.5/src | giac 2.0.0.21 nlpsolve, find_minimum, simulated_annealing, frank_wolfe, simplex_reduce: port facts for phase 1.5 | giac, nlpsolve, find_minimum, simulated_annealing, frank_wolfe, simplex_reduce, phase-1.5, optimization | retriever-code | 2026-11-24 | 101 lines
R1.5-002 | phase 1.5 (Optimization II) port plan. Which lp_* options change answers, port size per part, done-when rows. | giac 2.0.0.21 lpsolve: options, algorithms, messages, tie-breaks | giac, lpsolve, lp, branch-and-bound, simplex, gomory, glpk, phase-1.5 | retriever-code | 2026-11-24 | 75 lines
R1.5-003 | verify and complete the earlier findings (a) to (f) for the phase 1.5 plan (Optimization II); scope excludes native/optimization/ | How a new cas native plugs in: verified facts at HEAD | cas, natives, limits, interrupt, registry, random, lambda, compile, coefficientlist | retriever-code | 2026-11-24 | 85 lines
R1.5-004 | phase 1.5 plan (Optimization II): ## Interfaces and done-whens for lp_*/nlp_* options | Optimization natives at HEAD: LPSolve, integer, simplex, NLPSolve, global, TPSolve | optimization, lpsolve, nlpsolve, simplex, integer, global, tpsolve, cas-416, cas-417, cas-418, cas-519 | retriever-code | 2026-11-24 | 102 lines

## Audit
- seed: median 24k, max 25k, n=10 (phase 1.5)
- previous phase 1.4: median 24k, growth 0.6%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 7983 bytes
- cookbook/INDEX.md: 11197 bytes
- rules/INDEX.md: 3242 bytes
### Carry audit — phase 1.5 (2026-11-24T04:40:35Z → open UTC, 1 sessions, 798 requests)

| file (read) | chars | n |
|---|---|---|
| nlpsolve.rs | 123.8k | 10 |
| lpsolve.rs | 108.6k | 21 |
| sparse.rs | 61.9k | 4 |
| gcd.rs | 49.7k | 14 |
| optimization.cc | 48.5k | 13 |
| optimization.cc | 45.9k | 7 |
| risch.cc | 45.1k | 3 |
| integrate.rs | 43.5k | 7 |
| builtin.rs | 40.9k | 8 |
| find_minimum.rs | 37.7k | 3 |
| file (write) | chars | n |
|---|---|---|
| page.md | 30.0k | 1 |
| tables.rs | 27.7k | 1 |
| proto.py | 25.1k | 1 |
| find_minimum.rs | 21.2k | 1 |
| frank_wolfe.rs | 20.7k | 1 |
| doc31.md | 20.0k | 1 |
| lpsolve.rs | 19.2k | 12 |
| edit1.py | 16.2k | 2 |
| result kind | chars | n |
|---|---|---|
| Read | 2.6M | 591 |
| bash other | 1.5M | 647 |
| Grep | 864.5k | 408 |
| tools/plan_edit.py | 241.4k | 16 |
| Glob | 40.0k | 26 |
| run.sh | 39.8k | 63 |
| Write | 16.8k | 100 |
| other | 14.6k | 88 |
| Agent | 13.3k | 13 |
| tools/commit_task.sh | 9.8k | 26 |
- spilled: 4 results, 127.0k chars on disk, 0 read whole
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,252 (n 34, prev 6,252) · retriever-digest 6,333 (n 4, prev 6,333) · retriever-web 4,600 (n 2, prev 4,600)
- outline credit: 9 outlines · 6 followed by a ranged read · 104.7k chars credited
- warm pings: 3 pings over 1 runs, 1 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.0862 + relay $0.0000 + re-arm $0.0000 = $0.0862 vs rewrites replaced $1.0845, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 396, relay 0, re-arm 8, other 34, relay cost $0.0000
- retriever re-asks: 0 of 0 retriever briefs repeat a lookup of the same run
- noise: 18 lines 1.1k chars — No such file or directory: 14 lines/1.0k, usage: : 4 lines/136
- price: read $0.20/Mtok · 1h write $8.00/Mtok · output $20.00/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 6731 chars, 798 requests (prev 5995 chars, 1286 requests, growth 12.3%)
- flag: 3 whole-plan reads by pa-expert
- flag: 4 tool-source reads by main, pa-coder, retriever-code
  /home/user/code/cas/.claude/settings.local.json  61 chars  main
  /home/user/code/cas/crates/cas-core/src/native/polytools/roots.rs  622 chars  retriever-code
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  1.4k chars  pa-coder
  /home/user/code/cas/crates/cas-core/src/native/polytools/two.rs  646 chars  retriever-code
- flag: 3 whole reads over 20.0k by main
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-poly/src/sparse.rs  role: main  chars: 20.1k
  file: /home/user/code/cas/crates/cas-core/src/builtin.rs  role: main  chars: 23.1k

## Agent runs
- T10 | pa-claude | partner | claude-opus-5-5 | medium | ctx 138k | $2.579 | saved $40.97 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 101k | $1.621 | saved - | completed | parent 1b3aa32d
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 103k | $1.825 | saved $3.16 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 167k | $4.005 | saved - | completed | parent ae3f9eec
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 124k | $2.160 | saved $5.42 | completed | parent 1b3aa32d
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 212k | $4.827 | saved - | completed | parent a7337806
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 53k | $0.582 | saved - | completed | parent a7337806
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 147k | $3.239 | saved $6.03 | completed | parent 1b3aa32d
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 193k | $4.092 | saved - | completed | parent a0051688
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 165k | $3.446 | saved - | completed | parent a0051688
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 77k | $1.217 | saved $0.19 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 84k | $1.463 | saved - | completed | parent ac8e86cc
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 111k | $2.018 | saved $3.55 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 144k | $3.176 | saved - | completed | parent a6d999c2
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 139k | $2.788 | saved - | completed | parent a6d999c2
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 87k | $1.421 | saved $0.73 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 73k | $1.184 | saved - | completed | parent a9c0410f
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 59k | $0.739 | saved - | completed | parent a9c0410f
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 45k | $0.385 | saved - | completed | parent 1b3aa32d
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 59k | $0.730 | saved - | completed | parent 1b3aa32d
- T10 | pa-expert | expert | claude-opus-5-5 | medium | ctx 57k | $0.658 | saved - | completed | parent 1b3aa32d

## Discussions
- (none)

## Deferred
- from T10: the close may delete .run/1.5/ whole. (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-24 router: T1 next -> done
- 2026-11-24 router: T2 next -> done
- 2026-11-24 router: T3 next -> done
- 2026-11-24 router: T4 next -> done
- 2026-11-24 router: T5 next -> done
- 2026-11-24 router: T6 next -> done
- 2026-11-24 router: T7 next -> done
- 2026-11-24 router: T8 next -> done
- 2026-11-24 router: T9 next -> done
- 2026-11-24 router: T10 next -> done

## Plain-English Recap
Phase 1.5 Optimization II: CAS-DOC-57 (still a draft in Lific) sorts 50 giac names; tickets CAS-600..610 are done. Built: lpsolve's 31 lp_* option builtins (nonnegint alias), maximize/assume forms, mixed integer and binary variables, every lp_method on our exact simplex with decimals made exact; giac's branch and bound ported (5 node and 4 variable strategies, Gomory cuts with the true right side where giac's rescaled cut can drop integer points, a bound-checked rounding heuristic, node/depth/iteration limits and the gap; deterministic ties; lp_timelimit accepted but no clock); nlpsolve's 14 nlp_* options, linear problems routed to LPSolve, integer NLP by enumeration or capped branch and bound (nlpinteger.rs); simplex_reduce as a line-by-line port of misc.cc with giac's mixed exact/decimal entries; find_minimum (Brent candidate, then a certified ball by interval Newton or cover); frank_wolfe with a ball-certified gap; simulated_annealing (GSL siman loop, seeded, user functions under limits, exact midpoints per step). CoefficientList is not buggy: a ball-coefficient unit test confirms leading-first order (CAS-607). Final gate green on a2ba5de with no fixes needed (CI 45350578782, Fuzz 45350578815, 56 targets). Cookbook C0056-C0064. Notes: the local stable clippy (1.99.1) is older than CI's and missed a lint; many giac integer/inequality nlpsolve probes fail, recorded as departures against known optima. Open: publish CAS-DOC-57; docs-only commits unpushed until the next phase's first push.
