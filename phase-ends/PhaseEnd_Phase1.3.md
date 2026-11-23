# PhaseEnd — Phase 1.3: Calculus, expressions and integers parity (implements Gen 1.3)
Approved: 2026-11-23 | Closed: 2026-11-23 | Planner: claude-opus-5-5/medium | Tasks: 12 done, 0 superseded

## Milestone
Milestone: The 81 giac Calculus, Expressions, Integers and list entries this plan builds (## Context, from R1.3-007) run under their giac names, each with a giac oracle test or a recorded departure, its doc chapter and, for a new builtin, a fuzz entry; Int is a recorded departure and the entries ## Context leaves out are recorded as out, with their reasons, on CAS-DOC-55; CI and Fuzz pass on the pushed master head, CAS-DOC-55 records the area and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-55 and one ticket per family | Done: CAS-DOC-55 "Calculus, expressions and integers parity" exists (## Scope, §1-§8 for T2-T9, §9 Oracle and tests); 38 tickets CAS-548..CAS-585 filed, status todo. Summary: The phase's design page is in Lific. It lists the 81 giac entries the phase builds, the ones it leaves out with their reasons, one section per build task with the algorithm, what it builds on, the oracle plan and giac's examples, and the departures. Each family of entries has a ticket in its module, so every build task's commits can name their ticket. | commit: - | tasks/T1.md | logs/T1.md
- T2 — integers: chrem, lgcd, fracmod, a certified is_prime, bernoulli_mod, icomp; Abs of a complex ball containing 0 | Done: giac's chrem, lgcd, fracmod/iratrecon, is_prime (proven, with giac's certificate), bernoulli_mod and icomp are builtins with oracle rows; ComplexBall::abs gives [0, u] around 0; CI green on ee0eab5. Summary: Six of giac's integer commands now exist and answer exactly as giac does, checked against libgiac on CI. is_prime only reports a prime after proving it (deterministically below 2^64, by the Brillhart–Lehmer–Selfridge n−1 test above); a prime it cannot prove stays unevaluated with a warning. The absolute value of a complex ball that touches 0 is now a usable interval instead of "indeterminate". All of 1.2's unpushed commits went out with this push. | commit: ee0eab5 | tasks/T2.md | logs/T2.md
- T3 — list commands: find, index, counts, deltas, sequences, folds, insert, subsop, heaps, permutations | Done: 25 giac list commands run under their giac names (CAS-555..561), with unit tests, oracle and departure rows in crates/cas-giac/tests/list_commands.rs; master pushed, CI green on 2da3db9 (run 45114421938); CAS-555..561 closed. Summary: Xcas's list commands (searching and counting, differences, builders, apply and folds, insert and subsop, heaps, next and previous permutation, the group two permutations generate) now answer as giac does. Where giac changes a variable in place (heapify, heappop) ours returns the resulting heap, and the difference is recorded. Two giac commands read undefined memory on some inputs; those inputs are tested on our side only. | commit: 2da3db9 | tasks/T3.md | logs/T3.md
- T4 — expressions: power expansion, common denominators, numerators, logarithms in a base, angles, unit vectors, algsubs, regroup | Done: 15 giac expression names built (powexpand, pow, comDenom, getNum, getDenom, simp2, regroup, algsubs, sq, logb, log2, log10, alog10, radians, degrees, normalize/unitV); CI green on d6df6d0; CAS-562..567 done. Summary: The engine now answers giac's expression commands by their giac names. Most reuse existing machinery (common denominators, norms, logarithms, Groebner elimination); algsubs follows giac's own method of eliminating variables from the substitution relation. Two answers differ in form from giac on purpose: power expansion is shown held, because our arithmetic would merge the powers back, and angle conversion answers exactly where giac gives decimals. | commit: d6df6d0 | tasks/T4.md | logs/T4.md
- T5 — expressions: variable lists, operands, evaluation, predicates, equation and list conversions, steps, Dirac, periodic | Done: 18 giac commands built (lname/indets, lvar, algvar, reorder, part, eval/evalm, evalb, epsilon2zero, isinf, isnan, compare, equal2diff, equal2list, list2exp, exp2list, linstep, simplifyDirac, periodic); CI green on 1a30dba (run 45120822285). Summary: Xcas's expression commands for listing variables and kernels, taking operands, evaluating and testing, converting equations and solution lists, rewriting steps and deltas, and making periodic functions now work. They follow giac's own source. Where giac's answer depends on its internal ordering or its print form, the oracle tests compare by value or as sets, or record the difference with giac's exact text. | commit: 1a30dba | tasks/T5.md | logs/T5.md
- T6 — calculus tools: ibpdv, ibpu, preval, function_diff, revert, sum_riemann, wz_certificate, euler_lagrange, jacobi_equation, conjugate_equation | Done: ten giac calculus names built (one alias, fonction_derivee) with unit tests, oracle and departure rows, module docs, BUILTINS 565; Integrate[f] integrates in x; master pushed. Summary: Ten Xcas calculus commands now exist in cas: step-by-step integration by parts, evaluating an antiderivative between bounds, the derivative of a function object, inverting a power series, turning a Riemann sum into its integral, Wilf–Zeilberger certificates and the Euler–Lagrange, Jacobi and conjugate-point equations of the calculus of variations. Each is a port of giac's own algorithm, compared with giac on CI. Where giac is wrong (one WZ example, coupled Euler–Lagrange equations) ours is correct and the difference is recorded. | commit: 33f610b | tasks/T6.md | logs/T6.md
- T7 — function study: domain, sign and variation tables, convexity | Done: giac's domain, convex, tabsign and tabvar run under their giac names (CAS-581..583), with unit tests, oracle and departure rows in crates/cas-giac/tests/function_study.rs; master pushed, CI green on 7b2801a (run 45136885878); CAS-581..583 closed. Summary: Xcas's function-study commands now answer as giac does: the domain of an expression as conditions, convexity as true, false or the conditions giac prints, and the sign and variation tables with giac's arrows, bars and convexity strings. All four are ports of giac's own source, checked first with Python prototypes against giac's answers. The parametric variation table and the text of giac's y' labels are recorded differences. | commit: 7b2801a | tasks/T7.md | logs/T7.md
- T8 — Kovacic's algorithm: kovacicsols | Done: kovacicsols (KovacicSols) finds Liouvillian solutions of a*y''+b*y'+c*y=0 by Kovacic's cases 1-3, ported from giac desolve.cc; CI green on 263e564 (all 5 jobs); CAS-584 done. Summary: The new kovacicsols command solves second-order linear equations with rational coefficients in closed form when such a solution exists, as giac's does. It was first written as a Python prototype that reproduced all seven of giac's example answers, then ported to Rust. Where giac's own code is wrong (case 3 with a first-derivative term, fractions in the list form), ours gives the correct answer and the difference is recorded as a departure. DSolve is unchanged. | commit: bab9f39 | tasks/T8.md | logs/T8.md
- T9 — Risch integration of exp and ln towers: risch; Int as a departure | Done: `risch(f,x)` is giac 2.0.0.21's transcendental Risch algorithm over exp/ln towers, ported via a checked sympy prototype; every antiderivative guarded by D; oracle rows, departures (Int among them), chapter; CI green on 6d00a69. Summary: The CAS now has giac's `risch` command, which decides whether a function built from exponentials and logarithms has an elementary antiderivative and finds it, writing the non-elementary exponential part with Ei. It was ported line by line from giac's source after a Python prototype matched all 87 giac probes and a random batch of 150 integrands. Every answer is checked by differentiation before it is returned, so a wrong antiderivative turns into an unevaluated call instead. The ordinary integrate command is unchanged and does not use it. | commit: 6d00a69 | tasks/T9.md | logs/T9.md
- T10 — push master for the phase's final CI and Fuzz runs | Done: master pushed at e1ceb1d after a green local gate; CI run 45171138676 and Fuzz run 45171138561 started on that head. Summary: The local checks (formatting, stable clippy, the fuzz crate's build check) passed on the committed head, and master was pushed so the phase's final CI and Fuzz runs could start. Both runs are now going on the pushed commit; the next task reads their results. | commit: e1ceb1d | tasks/T10.md | logs/T10.md
- T11 — CI and Fuzz check on the pushed head; tickets done and CAS-DOC-55 in step | Done: CI 45185597012 and Fuzz 45185596994 green on the pushed head f3fab1f; CAS-548..585 all done; CAS-DOC-55 §9 lists the departures recorded in crates/cas-giac/README.md ## Tests. Summary: The final fuzz runs found three problems, each fixed and pushed until CI and the full fuzz run passed together. Two were real engine bugs: the Fourier transform used a principal-value formula for a pole off the real line, which gave a wrong inverse, and Risch crashed on a malformed exponential with no argument. The third was a check in the fuzz harness that compared complex integrals by their real parts only. The phase's tickets are all closed and its design page lists the departures each test file records. | commit: f3fab1f | tasks/T11.md | logs/T11.md
- T12 — phase cleanup | Done: the giac lookup scripts caseval.py and aide.py live in scripts/ and the card names them; nothing else in .run/1.3/ is needed. Summary: The phase scratch folder was sorted. The two small giac lookup scripts that every phase rewrote were moved into the repo's scripts folder and listed on the card. Everything else there is regenerable, refetchable, already banked in the cookbook or documented in the code, or a one-off check, so the close can delete the folder. | commit: - | tasks/T12.md | logs/T12.md

## Decisions that still bind
- T12 — giac lookups outside tests use scripts/caseval.py and scripts/aide.py (card Paths line); planners stop rewriting them in .run/<phase>/.
-

## Rules proposed
- T9 — an expert's coder brief names logs/T<n>.c<k>.progress.md as the coder's handoff file, never TASK_PROGRESS.md.

## Cookbook entries added
C0020 | Porting giac's C integer code: truncated %, signed invmod, standalone rustc check | giac,port,prototype,integer,rustc | 2026-11-23 | 1.3/T2 | CAS-552 bernoulli_mod, CAS-548 chrem
C0021 | libgiac PARI failures abort the process after the call: test risky rows alone | giac,oracle,pari,crash,libgiac | 2026-11-23 | 1.3/T2 | CAS-551 is_prime
C0022 | Stays-unevaluated tests: write the input in canonical orderless order | test,unevaluated,orderless,cas-core | 2026-11-23 | 1.3/T2 | CAS-548 chrem stays row
C0023 | Oracle rows for a sequence-valued builtin: bracket the input | oracle,giac,sequence,makesuite | 2026-11-23 | 1.3/T3 | CAS-557
C0024 | giac list ports: no conjugation in vector products; NULL unquoted in oracle rows | giac,oracle,dot,multiply,subsop,NULL | 2026-11-23 | 1.3/T3 | CAS-558
C0025 | giac in-place heap commands: port libstdc++ stl_heap, return the stored value | giac,port,heap,in-place,crash | 2026-11-23 | 1.3/T3 | CAS-560
C0026 | giac forms the evaluator re-merges: return Hold, oracle by value after release | giac,oracle,hold,evaluator,departures | 2026-11-23 | 1.3/T4 | CAS-562 powexpand
C0027 | giac epsilon is 1e-12 by default; its boundary drifts in the shared oracle | giac,oracle,epsilon,float,epsilon2zero | 2026-11-23 | 1.3/T5 | T5.c1 CAS-569
C0028 | Fuzz-fed natives: compare decimals by to_f64 and total_cmp, never Float::to_rational | fuzz,float,panic,clippy,native | 2026-11-23 | 1.3/T5 | T5.c1 CAS-569
C0029 | Expression shapes natives meet: Range in x=a..b, Dirac derivatives, Abs of products, Times[-1,Plus] | native,cas-core,lowering,dirac,abs,range | 2026-11-23 | 1.3/T5 | T5.c2 CAS-570/572/573
C0030 | giac function answers: compare them applied, not printed (hidden ` x` param) | giac,oracle,function,lambda,clippy | 2026-11-23 | 1.3/T6 | T6.c1 CAS-577
C0031 | Gosper with symbolic parameters: per-factor resultants, B-A degree; WZ certificate via multiplier | gosper,sum,resultant,wz,performance | 2026-11-23 | 1.3/T6 | T6.c2 CAS-580
C0032 | Oracle rows on factored/signed rational answers: compare normal(ours - giac) | giac,oracle,equal,rational | 2026-11-23 | 1.3/T6 | T6.c2
C0033 | Reading giac calculus answers back: diff(f,x,2), $K_n constants, sequences, Together on radicals | giac,oracle,diff,symbols,together | 2026-11-23 | 1.3/T6 | T6.c3 CAS-575
C0034 | giac tabvar/tabsign: step_func_ reductions, columns and dedupe | giac,tabvar,tabsign,variation,port,source | 2026-11-23 | 1.3/T7 | CAS-582 T7.c2 prototype
C0035 | giac convex: symbolic Bunch-Kaufman ldl; giac's ldl command as oracle | giac,convex,ldl,bunch-kaufman,hessian,port,oracle | 2026-11-23 | 1.3/T7 | CAS-583 T7.c1
C0036 | Oracle cell compares: giac <>, Exp vs E^u, Solve on quotients, primitive conditions | giac,oracle,compare,unequal,exp,solve,format | 2026-11-23 | 1.3/T7 | CAS-581..583 T7.c3
C0037 | giac port traps: coeff highest-first, quoted '/' in lists, return false prints 0 | giac,port,coeff,quote,kovacic,prototype | 2026-11-23 | 1.3/T8 | CAS-584 T8.c1
C0038 | Kovacic port: list-form oracle rows, algebraic poles, systems by monomials | giac,oracle,dsolve,kovacic,quote,together,limits | 2026-11-23 | 1.3/T8 | CAS-584 T8.c2
C0039 | giac risch port: helper behaviours, the :693 flaw, tower ties | giac,risch,port,prototype,integrate,tower,gcd | 2026-11-23 | 1.3/T9 | CAS-585 T9.c1
C0040 | CI-only tests of a big port: one test per row group collecting every failure | tests,ci,port,oracle,risch,commit | 2026-11-23 | 1.3/T9 | CAS-585 T9.c3
C0041 | Zero-testing exp expressions: exponent atoms to fresh symbols, then Together | exp,zero-test,together,integrate,risch,derivative | 2026-11-23 | 1.3/T9 | CAS-585 T9.c6
C0042 | Ported natives indexing arguments: guard kernel arity at entry | fuzz,eval,panic,arity,risch,native | 2026-11-23 | 1.3/T11 | Fuzz 45172144522 risch Exp[]
C0043 | Fuzz numeric mismatch: check the harness side first; complex integrands by both parts | fuzz,harness,complex,fourier,nintegrate,triage | 2026-11-23 | 1.3/T11 | Fuzz 45178885251, 45136885779

## Research
R1.3-001 | phase 1.3 build tasks for Calculus (calculus of variations, Riemann sums, function derivatives, convexity, Kovacic, preval, series reversion, WZ certificates, domain). Map each giac name to the builtin or internal fn it would build on. | Calculus giac names: existing builtin hooks | calculus, dsolve, sum, limit, series, procedure, hessian, domain, phase-1.3 | retriever-code | 2026-11-23 | 130 lines
R1.3-002 | phase 1.3 plan, Expressions build tasks (rewrite, fraction, variable-list, conversion, predicate commands). Map each giac name to the closest existing builtin. | giac Expressions names: existing builtin match (en.txt, builtin.rs, native.rs) | giac, expressions, builtins, rewrite, fraction, predicates, conversion, phase-1.3 | retriever-code | 2026-11-23 | 106 lines
R1.3-003 | phase 1.3 plan, Integers build task; possible cas-num fix so Abs of a complex ball containing 0 returns [0, upper] as Arb does. | giac Integers names: existing builds, and complex Abs of a ball containing 0 | giac, integers, builtins, chinese-remainder, rational-reconstruction, primality, bits, cas-num, complex-ball, abs | retriever-code | 2026-11-23 | 72 lines
R1.3-004 | phase 1.3 plan, list-commands build task (Lific module Evaluator, CAS-DOC-28) | giac list names vs existing builtins | lists, giac, builtins, permutations, sets, index, fuzz | retriever-code | 2026-11-23 | 83 lines
R1.3-005 | phase 1.3 plan, risch (giac risch.cc port) and tabvar/tabsign; map what the code has to build on | Existing hooks for the Risch, Ei and tabvar/tabsign tasks | integrate, risch, ei, derivative, limit, sign, tabvar, tabsign, fuzz, giac-tests | retriever-code | 2026-11-23 | 107 lines
R1.3-006 | phase 1.3 plan `## Interfaces` (giac Calculus, Expressions, Integers and list commands): where a new giac command is added (names, registry, natives, fuzz, oracle, docs) | Where a new giac command is added (KNOWN verification and gaps) | giac, names, builtins, natives, fuzz, oracle, docs, phase-1.3 | retriever-code | 2026-11-23 | 60 lines
R1.3-007 | list the giac help entries still missing in the Calculus, Expressions, Integers and list areas, with kind, family, module and aide example, for the phase 1.3 "Calculus, expressions and integers parity" plan | Calculus, Expressions, Integers and list gaps for phase 1.3 (giac 2.0.0.21 aide_cas) | giac, gap, calculus, expressions, integers, lists, parity, phase-1.3 | retriever-digest | 2026-11-23 | 115 lines
R1.3-008 | map existing primitives for a giac Risch port into crates/cas-core/src/native/integrate/risch | Risch port: existing Rust primitives, file:line map | risch, integrate, cas-poly, sparse, gcd, partfrac, limits, builtins | retriever-code | 2026-11-23 | 95 lines

## Audit
- seed: median 24k, max 24k, n=12 (phase 1.3)
- previous phase 1.2: median 23k, growth 2.9%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 7793 bytes
- cookbook/INDEX.md: 7379 bytes
- rules/INDEX.md: 3096 bytes
### Carry audit — phase 1.3 (2026-11-23T03:07:47Z → open UTC, 1 sessions, 2172 requests)

| file (read) | chars | n |
|---|---|---|
| risch.cc | 126.8k | 9 |
| risch.py | 103.7k | 15 |
| desolve.cc | 83.9k | 5 |
| sparse.rs | 78.6k | 6 |
| parts.rs | 74.1k | 11 |
| misc.cc | 71.2k | 6 |
| gcd.rs | 59.1k | 16 |
| kovacic.py | 54.6k | 4 |
| tower.rs | 54.1k | 8 |
| poly.rs | 53.5k | 9 |
| file (write) | chars | n |
|---|---|---|
| tables.rs | 56.7k | 3 |
| kovacic.rs | 49.3k | 2 |
| risch.py | 48.0k | 1 |
| risch.rs | 34.6k | 9 |
| function_study.rs | 33.7k | 1 |
| doc55.md | 33.2k | 1 |
| step.rs | 30.6k | 2 |
| kovacic.py | 24.5k | 2 |
| result kind | chars | n |
|---|---|---|
| bash other | 3.3M | 1640 |
| Read | 2.5M | 507 |
| Grep | 779.4k | 374 |
| tools/plan_edit.py | 311.9k | 18 |
| run.sh | 133.9k | 162 |
| tools/card.py | 59.9k | 10 |
| tools/commit_task.sh | 36.7k | 60 |
| tools/outline.py | 33.6k | 21 |
| Glob | 32.4k | 17 |
| Agent | 30.0k | 29 |
- spilled: 1 results, 116.1k chars on disk, 0 read whole
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,252 (n 22, prev 6,252) · retriever-digest 6,321 (n 3, prev 6,321) · retriever-web 4,600 (n 2, prev 4,600)
- outline credit: 33 outlines · 15 followed by a ranged read · 419.0k chars credited
- warm pings: 8 pings over 2 runs, 2 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.1596 + relay $0.0000 + re-arm $0.0000 = $0.1596 vs rewrites replaced $1.2440, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 905, relay 0, re-arm 31, other 35, relay cost $0.0000
- retriever re-asks: 0 of 5 retriever briefs repeat a lookup of the same run
- noise: 24 lines 1.5k chars — No such file or directory: 15 lines/1.1k, usage: : 9 lines/356
- price: read $0.20/Mtok · 1h write $7.97/Mtok · output $19.92/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 3330 chars, 2172 requests (prev 3433 chars, 1358 requests, growth -3.0%)
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
- T3 | pa-claude | partner | claude-opus-5-5 | medium | ctx 132k | $3.170 | saved $44.92 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 113k | $1.974 | saved - | completed | parent 1b3aa32d
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 162k | $4.260 | saved $16.33 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 207k | $4.934 | saved - | completed | parent a7b5c05b
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 175k | $4.244 | saved - | completed | parent a7b5c05b
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 137k | $3.357 | saved $15.44 | completed | parent 1b3aa32d
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 119k | $2.283 | saved $0.44 | completed | parent acdb69c4
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 26k | $0.008 | saved - | completed | parent a2ac5c5b
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 139k | $2.544 | saved - | completed | parent acdb69c4
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 116k | $2.003 | saved - | completed | parent acdb69c4
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 129k | $2.654 | saved $14.87 | completed | parent 1b3aa32d
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 154k | $3.393 | saved - | completed | parent ac6db483
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 61k | $0.730 | saved - | completed | parent ac6db483
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 114k | $2.194 | saved $16.29 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 191k | $4.442 | saved - | completed | parent a675a6b6
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 167k | $3.422 | saved - | completed | parent a675a6b6
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 101k | $1.985 | saved $24.43 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 118k | $2.234 | saved - | completed | parent ab680641
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 204k | $6.309 | saved - | completed | parent ab680641
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 167k | $4.039 | saved - | completed | parent ab680641
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 151k | $3.119 | saved $15.85 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 97k | $1.561 | saved - | completed | parent a7482c56
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 214k | $5.100 | saved - | completed | parent a7482c56
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 190k | $4.533 | saved - | completed | parent a7482c56
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 106k | $1.710 | saved $9.41 | completed | parent 1b3aa32d
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 91k | $1.585 | saved - | completed | parent a2c3678b
- T8 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 30k | $0.011 | saved - | completed | parent a2c3678b
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 266k | $6.807 | saved - | completed | parent a2c3678b
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 209k | $6.476 | saved $44.03 | completed | parent 1b3aa32d
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 289k | $8.859 | saved $0.05 | completed | parent a76f3608
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 50k | $0.017 | saved - | completed | parent a76f3608
- T12 | pa-claude | partner | claude-opus-5-5 | medium | ctx 81k | $0.856 | saved - | completed | parent 1b3aa32d
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 306k | $7.674 | saved $0.78 | handoff | parent a76f3608
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 172k | $2.991 | saved $0.58 | completed | parent a76f3608
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 31k | $0.008 | saved - | completed | parent a87168c1
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 147k | $2.600 | saved - | completed | parent a76f3608
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 100k | $1.473 | saved - | completed | parent a76f3608
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 144k | $2.671 | saved $1.22 | completed | parent a76f3608
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 33k | $0.010 | saved - | completed | parent a4724289
- … 4 more: PY ~/.claude/pa3/pa_ledger.py report

## Discussions
- (none)

## Deferred
- from T12: the close deletes .run/1.3/; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-23 router: T1 next -> done
- 2026-11-23 router: T2 next -> done
- 2026-11-23 router: T3 next -> done
- 2026-11-23 router: T4 next -> done
- 2026-11-23 router: T5 next -> done
- 2026-11-23 router: T6 next -> done
- 2026-11-23 router: T7 next -> done
- 2026-11-23 router: T8 next -> done
- 2026-11-23 router: T9 next -> done
- 2026-11-23 router: T10 next -> done
- 2026-11-23 router: T11 next -> done
- 2026-11-23 router: T12 next -> done

## Plain-English Recap
Phase 1.3 Calculus, expressions and integers parity: CAS-DOC-55 sorts 81 giac entries; tickets CAS-548..585 are done. Built: integers (chrem, lgcd, fracmod/iratrecon, a certified is_prime with a BLS75 n-1 prover above 2^64 and giac's certificate, bernoulli_mod, icomp) and the Abs fix for a complex ball containing 0 (CAS-554, [0,u] as Arb gives); 25 list commands (find, counts, folds, subsop, heaps, permutations); 33 expression commands (powexpand returned held, comDenom, algsubs, part, eval, compare, linstep, simplifyDirac, periodic); calculus tools (ibpdv/ibpu, preval, function_diff, revert, sum_riemann, wz_certificate, euler_lagrange, domain, convex, tabsign/tabvar); kovacicsols (Kovacic cases 1-3); risch for exp/ln towers ported from giac's risch.cc, checking each answer by differentiation; Int is a departure. All were ported from giac source, most checked by a Python prototype first. Recorded departures include four places where giac is wrong (wz_certificate, euler_lagrange with coupled unknowns, Kovacic case 3 with y', linstep) and two wrong giac answers in risch. The final gate is green on f3fab1f (CI 45185597012, Fuzz 45185596994) after three fixes: Fourier rules now require a real pole (CAS-441; cost: θ(x-a) with symbolic a stays unevaluated), the risch arity guard, and the fourier_transforms harness comparing both parts. scripts/caseval.py and scripts/aide.py are kept in the repo. Cookbook C0020-C0043. Open: the developer's ruling on local throwaway example programs (R33; T9 needed 6 CI rounds); backlog CAS-586 (APRCL) and CAS-587 (risch sqrt split); the card's dangling docs/ops/INDEX.md pointer; docs-only commits unpushed until the next phase's first push.
