# PhaseEnd — Phase 1.8: Validation (implements Gen 1.8)
Approved: 2026-11-26 | Closed: 2026-11-27 | Planner: claude-opus-5-5/medium | Tasks: 12 done, 0 superseded

## Milestone
Milestone: A one-hour-per-target fuzz campaign (`gh workflow run fuzz.yml -f time_per_target=3600`) and the oracle's random-input sweeps over the generation's builtins pass on the pushed master with every finding fixed under a regression test, `bash scripts/giac_gap.sh` leaves no computational giac command of an existing module unaccounted for (built, or a departure recorded in crates/cas-giac/README.md ## Tests), and CI and Fuzz are green on the final master head (`gh run list --branch master --limit 3`).
Confirmed by the developer at the close.

## Tasks
- T1 — correct Weibull moment fit with a departure row (I2) and diff(e,x,n) as the nth derivative (CAS-636) | Done: randvar(weibull, mean=, variance=) fits the true moment equation (README departure vs giac); D[f,x,n] with n>0 is the nth derivative; an integrate_two fuzz recursion (CAS-638) fixed; CI and Fuzz green on 777f31e. Summary: The Weibull fit from a mean and a variance copied giac's equation, which forgets to square a denominator, so the law it built did not have the variance asked for; ours now solves the right equation, starting Newton from a standard approximation so it converges, and the difference from giac is recorded as a departure. diff(f,x,2) now means the second derivative, as in giac, so a test mask that hid the difference is gone. | commit: 777f31e | tasks/T1.md | logs/T1.md
- T2 — gap ledger: every missing giac entry accounted for, the 49 split into build groups, tickets filed (I2, R44) | Done: crates/cas-giac/README.md ## Tests ends with the gap ledger; scripts/gap_ledger.py re-tallies it; tickets CAS-639..648 filed. Summary: Every one of the 1004 giac entries we lack is now accounted for in the oracle crate's README: 913 out of scope with the reason, 12 already recorded departures, 29 alias spellings with what they map to, and 50 computational entries each assigned to a build task (T3-T6) or recorded as not built. A small script recounts the buckets and fails if the list changes. One tracker ticket per build family is filed. | commit: - | tasks/T2.md | logs/T2.md
- T3 — Algebra matrix leftovers (REDIM, border, semi_augment, simult, newMat, Fill, upper, lower, a2q, q2a, basis, changebase, moyal, trn) | Done: 13 giac matrix entries run under their giac names with oracle rows; moyal, Det/Inverse are README departures; charpoly fadeev and convert array fixed. CI 45751174928 green on 67580d5. Summary: The giac matrix commands we still lacked now work: reshaping, bordering, stacking, solving with a matrix right-hand side, empty and filled matrices, triangles and case changes, quadratic forms both ways, change of basis, conjugate transpose and basis. Each is checked against giac on CI. moyal is recorded as not built because giac itself does not compute it, and Det/Inverse stay inert like giac's. | commit: 67580d5 | tasks/T3.md | logs/T3.md
- T4 — Expressions, Calculus, Integers and Numerical-analysis gap leftovers (CAS-640..643) | Done: 19 giac entries run under giac's names (root atan2 has fmod modf multinomial order maple_ifactors avgRC box_constraints convert-expln additionally addtable period linabs lin shift_phase non_recursive_normal envelope), each with oracle or departure rows, a fuzz BUILTINS entry and a module doc bullet; ci.yml green on 8648565. Summary: The last expression, calculus, integer and box-constraint commands giac has and we lacked are now built, so the gap ledger has no open entries for these modules. Two of them change session state: an extra assumption (additionally) and a user-supplied transform pair (addtable). envelope is the one deliberate difference: giac draws a plot, ours returns the curve's equations. Three new cookbook entries record what tripped the coders. | commit: 8648565 | tasks/T4.md | logs/T4.md
- T5 — Language, Polynomials, Statistics and Infrastructure leftovers and kmeans/cluster options | Done: back range shift suppress at rdiv getType subs format dayofweek convert-list roots truncate divis choice built with oracle rows; kmeans/cluster output=part, string clustering, index= built; CAS-630's option list complete in README ## Tests; eval fuzz SIGSEGV (T4's Period) fixed; CI and Fuzz green. Summary: The last language, polynomial, statistics and date entries giac has that we lacked now run under their giac names. Each is checked against giac on CI, or recorded as a deliberate difference. The clustering commands gained giac's remaining options: listing the clusters themselves, clustering words by edit distance, and choosing the number of clusters by a named quality index. Along the way a crash in the fuzzer, introduced by the previous task's period command, was found and fixed, and the fuzz workflow gained a replay mode that names such crashes. | commit: 87dc694 | tasks/T5.md | logs/T5.md
- T6 — Signal processing leftovers (stft, istft, hht, boxcar, tri, sinc) and hilbert(f,t) | Done: sinc, tri, boxcar, stft, istft, hht(output=matrix) and hilbert(f,t)/hilbert(f) run under their giac names with oracle and departure rows; CI 45810495716 and Fuzz 45810495765 green on d9b6d57. Summary: The last six signal-processing entries giac has and we lacked are built, plus the transform form of hilbert beside the existing list form, each ported from giac's own source and checked against giac's answers on CI. Where our Fourier transforms cannot close a Hilbert transform giac closes (rational functions, sinc, the box), ours stays as giac's inert Hilbert(f,t) and the row is a recorded departure. A fuzz failure met on the way came from the fuzz check itself, not the engine, and the check is fixed. | commit: d9b6d57 | tasks/T6.md | logs/T6.md
- T7 — random-input oracle sweeps: polynomial tools, real roots, gcd by name, partial fractions (G1, G2) | Done: 16 proptest! sweeps (128 cases a run each) cover every G1 and G2 family; CI's oracle job runs them green (run 45820024750). Summary: Random polynomial inputs are now checked against giac on every CI run for the real-root tools, the polynomial tools, the gcd algorithms by name and the partial-fraction forms. The sweeps found one disagreement: giac's modgcd in one variable throws away the integer factor of the gcd (it says the gcd of 4x and 6x is x, while its own gcd says 2x). Ours keeps the factor like gcd does, so this is recorded as a departure and the sweep skips that input class. | commit: aeb5afd | tasks/T7.md | logs/T7.md
- T8 — random-input oracle sweeps: lists, sorts, matrices, integers, TI names, exact numerics, ASCII strings, sound codes (G3, G4, G5) | Done: 41 new proptest! sweeps (128 cases each) cover R1.8-004's G3, G4 and G5 families; CI run 45829342969 on 623e8f2 is green in every job (the oracle job runs 63 random_ sweeps, T7's included). Summary: Random inputs are now checked against giac on every CI run for the list commands, TI sorts and parts, the matrix minors, echelon forms and pivots, the integer functions, the exact numerics (linspace, interpolation, numerical derivatives, splines, transport problems), ASCII strings and sound clips. The sweeps found one bug of ours, fixed: pivot changed rows that giac leaves alone. They also found several places where giac is wrong: VAS_positive with a root at 0, linear_interpolate with negative abscissae, 8-bit mixdown, splice with a long crossfade (which crashes giac), contains("",""), and trunc with negative digits. Each of these is recorded as a departure and left out of its sweep. | commit: 623e8f2 | tasks/T8.md | logs/T8.md
- T9 — raise fuzz.yml's job cap for one-hour targets and push the campaign head | Done: fuzz.yml caps each job at 120 min (`# was 60`, header says why); fuzz/README.md notes the cap; master 76c1018 pushed, CI and Fuzz (56/56, 10 min) green. Summary: Each fuzz job may now run two hours, so the one-hour campaign targets in T10 are not cut off by the old one-hour limit. Pushing the head showed two red oracle sweeps from T8. Both were giac quirks the sweeps let through: giac's modgcd drops the integer content for any one-variable input, and giac's parser inserts a space into ".." even inside a string. The sweeps now skip both cases and pin them as departures, and the head is fully green. | commit: 76c1018 | tasks/T9.md | logs/T9.md
- T10 — one-hour-per-target fuzz campaign on T9's head; findings triaged | Done: Fuzz run 45856971637 (time_per_target=3600, head 76c1018, nothing pushed meanwhile) concluded: 55/56 targets green, eval timed out; cause found and filed as CAS-649 (Statistics). Summary: Every fuzz target ran for an hour on the head T9 pushed. All but one came through clean, and none came near the 2 GB memory limit (the highest was about 0.5 GB). The eval target found one input that runs too long: MarkovChain checks its tolerance argument by evaluating it a second time with a fresh evaluation budget, so tolerances nested inside each other double the work at every level. The fix is filed for T11. | commit: - | tasks/T10.md | logs/T10.md
- T11 — campaign findings fixed under regression tests; final gate green; tickets done | Done: CAS-649 fixed with a cas-core regression test; CI, push Fuzz and a one-hour eval re-run green on 615a65e; CAS-636..649 done; gap at 952 missing / 786 covered, matching T2's ledger. Summary: The one input the fuzz campaign found that hung, MarkovChain with tolerances nested inside each other, now returns at once: the tolerance is checked to be a number by its shape before it is evaluated again. A test nests it 30 deep and requires an answer within two seconds. The test suite, the regular fuzz run and an hour of the target that failed all pass on the pushed head, and every ticket from this phase is closed. The giac gap shrank by exactly the 52 entries built this phase, leaving only the ones recorded as not built. | commit: 615a65e | tasks/T11.md | logs/T11.md
- T12 — phase cleanup | Done: what 1.8 still needs from .run/1.8/ is in the repo (T4 prototypes, the gap ledger's list); .run/logs/gen1-gap/ is deleted (R45). Summary: The phase's scratch folder was sorted. The prototypes behind the T4 leftover commands' expected values and the exact gap list the ledger's tally script was written against now live in the repo next to the code they back; everything else there is throwaway and goes with the close. The old generation-1 gap inputs were deleted as agreed. | commit: - | tasks/T12.md | logs/T12.md

## Decisions that still bind
-

## Rules proposed
- (none)

## Cookbook entries added
C0094 | A new rewrite rule breaks held-form pins in other crates' departure lists: grep first | rules,held,departures,tests,oracle,markup | 2026-11-26 | 1.8/T1 | T1.c1 CAS-636 D/nth
C0095 | Reach-inside wrappers skip held parts of their own head: Together[Together[r]] nests to the recursion limit | fuzz,recursion-limit,together,expand,held,integrate,zero-denominator | 2026-11-26 | 1.8/T1 | Fuzz 45734501868, T1.c2 777f31e CAS-638
C0096 | Gap ledger: tally giac-gap buckets by a failing script; verify alias claims in code | giac,gap,ledger,tally,alias,validation | 2026-11-26 | 1.8/T2 | R1.8-003 hand tally
C0097 | Clippy on oracle departure tables (split + helper) and Option<Option> natives | clippy,oracle,departures,native,lint | 2026-11-26 | 1.8/T3 | T3.c1 CAS-639
C0098 | Unit-test full forms: no decimals in parse_full_form; capital user names are $A | test,full-form,parse,lowering,cas-lang,ci | 2026-11-26 | 1.8/T3 | T3.c2 CAS-639
C0099 | A new builtin captures a bare word natives match by name: grep it first | builtin,names,lowering,regression,en.txt | 2026-11-26 | 1.8/T4 | T4.c1 CAS-642
C0100 | giac oracle rows: 12-digit decimals need a tolerance; session state is process-global | giac,oracle,tolerance,session,assume,addtable | 2026-11-26 | 1.8/T4 | T4.c1 CAS-640
C0101 | Reproducing giac's ratnormal forms: distribute, never Expand; n/Abs[d] | giac,oracle,together,expand,abs,forms | 2026-11-26 | 1.8/T4 | T4.c2 CAS-640
C0102 | Two coders in one working tree: build in a detached worktree | worktree,parallel,coders,build,shared-tree | 2026-11-26 | 1.8/T5 | T5.c2 CAS-646
C0103 | giac clustering index= traps: silhouette aborts libgiac, message without period | giac,oracle,cluster,kmeans,silhouette,abort,index | 2026-11-26 | 1.8/T5 | T5.c2 CAS-646
C0104 | giac's roundoff doubles (low 5 mantissa bits cleared): mask in double-counting ports | giac,float,double,roundoff,range,port,full-form | 2026-11-26 | 1.8/T5 | T5.c1 CAS-644
C0105 | giac kmeans cluster_crit indices: numeric quirks to reproduce (WG identity, exact det, pivots) | giac,kmeans,cluster,index,port,float,determinant | 2026-11-26 | 1.8/T5 | T5.c3 CAS-646
C0106 | eval SIGSEGV, no artifact: replay cached corpus on CI; bound native recursion with ev.enter | fuzz,eval,sigsegv,stack-overflow,replay,corpus,recursion,native | 2026-11-26 | 1.8/T5 | T5.c4 CAS-640
C0107 | giac's Hilbert transform via our Fourier pair; inert giac heads need a builtin | giac,hilbert,fourier,ifourier,inert,oracle | 2026-11-26 | 1.8/T6 | T6.c1 CAS-648
C0108 | giac's hht bins: top edge in fbins-2; near-zero imaginary part flips inst_freq | giac,hht,fft,bins,ties,oracle,prototype | 2026-11-26 | 1.8/T6 | T6.c2 CAS-648
C0109 | eval fuzz idempotence: limit error on re-eval of the answer is a cache artefact | fuzz,eval,idempotence,recursion,cache,harness | 2026-11-26 | 1.8/T6 | T6.c3 Fuzz 45802274513
C0110 | proptest failure on CI: transcribe the printed cc seed into tests/<file>.proptest-regressions | proptest,regression,ci,oracle,seed | 2026-11-26 | 1.8/T7 | CI 45818464735 random_modgcd
C0111 | giac's univariate modgcd drops the gcd's integer content; gcd, heugcd keep it | giac,gcd,modgcd,oracle,departure,normalisation | 2026-11-26 | 1.8/T7 | CI 45818464735; caseval probes
C0112 | Sweeping a giac predicate (==): generate only true pairs; trunc(x,n<0) and bernoulli_mod p=2 traps | proptest,oracle,giac,equality,symbolic,trunc,sweep | 2026-11-26 | 1.8/T8 | CI 45823905195 random_type
C0113 | giac string answers with quote or backslash do not read back: compare through asc | giac,strings,oracle,sweep,quote,printing | 2026-11-26 | 1.8/T8 | T8.c3 CI 45824175835
C0114 | Sweep inputs: run a generated batch in one giac process; caseval's per-line children hide heap overruns | giac,oracle,crash,heap,caseval,sweep,splice | 2026-11-26 | 1.8/T8 | T8.c3 CI 45824175835 sound SIGABRT
C0115 | Sweeping an optimiser (tpsolve, lpsolve): compare the optimum, check our point feasible; ties are common | giac,oracle,sweep,tpsolve,lpsolve,optimization,ties,proptest,linear_interpolate | 2026-11-26 | 1.8/T8 | T8.c4 random_tpsolve
C0116 | Sweep skip guards test the drawn input, not the strategy flag; giac spaces .. in strings | proptest,oracle,sweep,giac,modgcd,strings,departure | 2026-11-26 | 1.8/T9 | CI 45831113300 random_modgcd, random_words
C0117 | A native's nested Context::eval gets fresh limits: 2^depth on nested arguments | fuzz,eval,timeout,limits,native,approximation,nesting | 2026-11-26 | 1.8/T10 | Fuzz 45856971637, T10.c1, CAS-649

## Research
R1.8-001 | phase 1.8 plan, random-input oracle sweeps over gen 1 builtins, following tests/polynomials.rs | giac oracle random sweep pattern (crates/cas-giac) for phase 1.8 | cas-giac, proptest, oracle, giac, ci, phase-1.8 | retriever-code | 2026-11-26 | 51 lines
R1.8-002 | phase 1.8 plan lookups (Weibull equation fix, CAS-636, README ## Tests gap accounting) | Weibull moment fit, diff's order argument, cas-giac README Tests | weibull, randvar, diff, CAS-636, cas-giac, README | retriever-code | 2026-11-26 | 49 lines
R1.8-003 | bucket every line of .run/1.8/gap/missing_described.txt (N, P, S, I, M, D, A, C) for the phase 1.8 gap-check plan; name the bucket of changebase, semi_augment, divis, hilbert, kmeans/cluster options | Gap buckets for the 1004 described giac entries (missing_described.txt) | giac, gen1, gap, buckets, validation-1.8 | retriever-digest | 2026-11-26 | 84 lines
R1.8-004 | Phase 1.8 (Validation) plan. List the builtin families generation 1 added in phases 1.2 to 1.7, judge each for a random-input giac sweep, and propose 3-5 sweep groups. | Generation 1 builtin families: suitability for random-input oracle sweeps against giac (phases 1.2 to 1.7) | phase-1.8, validation, giac, oracle-sweep, departures, generation-1 | retriever-digest | 2026-11-26 | 69 lines

## Audit
- seed: median 24k, max 24k, n=12 (phase 1.8)
- previous phase 1.7: median 24k, growth -1.2%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 8330 bytes
- cookbook/INDEX.md: 20882 bytes
- rules/INDEX.md: 3242 bytes
### Carry audit — phase 1.8 (2026-11-26T01:02:25Z → open UTC, 1 sessions, 2032 requests)

| file (read) | chars | n |
|---|---|---|
| PHASE_PLAN.md | 129.5k | 9 |
| print.rs | 90.0k | 19 |
| context.rs | 76.8k | 35 |
| README.md | 67.4k | 6 |
| sparse.rs | 61.9k | 4 |
| signalprocessing.cc | 59.8k | 21 |
| clustering.rs | 59.6k | 10 |
| markup.cc | 55.4k | 16 |
| clustering.rs | 49.7k | 6 |
| gcd.rs | 49.7k | 14 |
| file (write) | chars | n |
|---|---|---|
| retriever-phase-plan-rationale-t2-t10.md | 36.1k | 1 |
| random_variables.rs | 29.3k | 1 |
| tables.rs | 27.7k | 1 |
| leftovers_three.rs | 25.2k | 1 |
| leftovers_two.rs | 23.7k | 1 |
| matrices_five.rs | 20.7k | 2 |
| doc31.md | 20.0k | 1 |
| reshape.rs | 15.1k | 3 |
| result kind | chars | n |
|---|---|---|
| Read | 4.4M | 1100 |
| bash other | 2.6M | 1500 |
| Grep | 1.7M | 765 |
| tools/plan_edit.py | 186.6k | 34 |
| run.sh | 118.0k | 135 |
| Glob | 65.6k | 45 |
| Agent | 39.6k | 39 |
| tools/commit_task.sh | 23.8k | 73 |
| other | 20.5k | 169 |
| tools/cookbook_add.sh | 20.2k | 30 |
- spilled: 2 results, 84.0k chars on disk, 0 read whole
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,115 (n 73, prev 6,115) · retriever-digest 6,040 (n 13, prev 6,040) · retriever-web 4,581 (n 3, prev 4,581)
- outline credit: 21 outlines · 6 followed by a ranged read · 346.3k chars credited
- warm pings: 77 pings over 19 runs, 29 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $1.5083 + relay $0.0000 + re-arm $0.0000 = $1.5083 vs rewrites replaced $16.5429, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 930, relay 0, re-arm 53, other 34, relay cost $0.0000
- retriever re-asks: 0 of 6 retriever briefs repeat a lookup of the same run
- noise: 33 lines 2.1k chars — No such file or directory: 28 lines/1.9k, usage: : 5 lines/233
- price: read $0.20/Mtok · 1h write $7.95/Mtok · output $19.88/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 4548 chars, 2032 requests (prev 7966 chars, 1899 requests, growth -42.9%)
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
- T7 | pa-claude | partner | claude-opus-5-5 | medium | ctx 96k | $1.537 | saved $1.04 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 75k | $1.402 | saved $6.52 | completed | parent 1b3aa32d
- T1 | pa-coder | coder | claude-opus-5-5 | medium | ctx 59k | $1.005 | saved - | completed | parent a2d19854
- T1 | pa-coder | coder | claude-opus-5-5 | medium | ctx 75k | $1.403 | saved $0.64 | completed | parent a2d19854
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 86k | $1.551 | saved $1.33 | completed | parent 1b3aa32d
- T2 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 34k | $0.011 | saved - | completed | parent a9de18a0
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 102k | $1.991 | saved $12.51 | completed | parent 1b3aa32d
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 132k | $3.239 | saved $0.06 | completed | parent a3f0bde1
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 18k | $0.005 | saved - | completed | parent aad27632
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 87k | $1.803 | saved - | completed | parent a3f0bde1
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 89k | $1.602 | saved $9.13 | completed | parent 1b3aa32d
- T4 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 43k | $0.012 | saved - | completed | parent a6bb98b5
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 174k | $4.466 | saved - | completed | parent a6bb98b5
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 185k | $5.023 | saved - | completed | parent a6bb98b5
- T4 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 31k | $0.011 | saved - | completed | parent a736c9e6
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 182k | $5.186 | saved $32.64 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 226k | $5.708 | saved $0.57 | completed | parent af9d85e0
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 222k | $6.313 | saved - | completed | parent af9d85e0
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 308k | $8.837 | saved $1.77 | handoff | parent af9d85e0
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 83k | $1.568 | saved - | completed | parent af9d85e0
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 59k | $0.993 | saved - | completed | parent af9d85e0
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 142k | $3.242 | saved $9.70 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 134k | $2.612 | saved - | completed | parent ac2a9add
- T6 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 27k | $0.010 | saved - | completed | parent a1880c37
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 209k | $5.529 | saved - | completed | parent ac2a9add
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 78k | $1.345 | saved - | completed | parent ac2a9add
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 81k | $1.472 | saved $0.36 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 45k | $0.548 | saved - | completed | parent a28c8747
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 55k | $0.735 | saved - | completed | parent a28c8747
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 126k | $3.287 | saved $6.19 | completed | parent 1b3aa32d
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 109k | $2.467 | saved - | completed | parent a557951b
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 88k | $1.762 | saved - | completed | parent a557951b
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 76k | $1.263 | saved - | completed | parent a557951b
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 89k | $1.459 | saved - | completed | parent a557951b
- T8 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 10k | $0.002 | saved - | completed | parent a0c1be6e
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 63k | $1.071 | saved $0.00 | completed | parent 1b3aa32d
- T10 | pa-expert | expert | claude-opus-5-5 | medium | ctx 71k | $1.214 | saved - | completed | parent 1b3aa32d
- T10 | pa-coder | coder | claude-opus-5-5 | medium | ctx 54k | $0.734 | saved - | completed | parent a4338732
- T11 | pa-expert | expert | claude-opus-5-5 | medium | ctx 72k | $1.277 | saved - | completed | parent 1b3aa32d
- … 3 more: PY ~/.claude/pa3/pa_ledger.py report

## Discussions
- (none)

## Deferred
- from T12: the close deletes .run/1.8/ whole; nothing in it is referenced anymore. (unconsumed)

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-26 router: T1 next -> done
- 2026-11-26 router: T2 next -> done
- 2026-11-26 router: T3 next -> done
- 2026-11-26 router: T4 next -> done
- 2026-11-26 router: T5 next -> done
- 2026-11-26 router: T6 next -> done
- 2026-11-26 router: T7 next -> done
- 2026-11-26 router: T8 next -> done
- 2026-11-26 router: T9 next -> done
- 2026-11-26 router: T10 next -> done
- 2026-11-26 router: T11 next -> done
- 2026-11-26 router: T12 next -> done

## Plain-English Recap
Phase 1.8 Validation: final gate green on 615a65e (CI 45885347583, push Fuzz 45885347505, one-hour eval re-run 45891167338). CAS-636..649 are done. Built and fixed:
- Fixes: the Weibull moment fit uses the correct equation, with giac's as a departure row (CAS-637); diff(e,x,n) takes the nth derivative via rule D/nth (CAS-636); nested Together on a zero denominator ends (CAS-638); pivot leaves rows with 0 in the pivot column alone (CAS-528); period's recursion is bounded; MarkovChain refuses a non-numeric eps before re-evaluating it, which removes the 2^depth cost the campaign found (CAS-649).
- Gap ledger (end of ## Tests in crates/cas-giac/README.md, scripts/gap_ledger.py, input kept in crates/cas-giac/gap/missing_described-2026-11-26.txt): each of the 1004 missing entries sits in exactly one bucket. 50 were computational. T3-T6 built 46 of them plus 15 alias spellings, so coverage went from 734 to 786 and 952 are still missing. Not built (R112): GF (with order in a finite field), hasard, taux_accroissement, moyal (a stub in giac 1.9). The rest are out of scope: new module, plotting, syntax, I/O.
- Built under giac names: the T3 matrix commands, 19 T4 expression/calculus/integer entries, the T5 list/date/polynomial entries, choice, kmeans/cluster output=part, index= and string clustering, and the T6 sinc/tri/boxcar/stft/istft/hht and hilbert(f,t).
- 57 random oracle sweeps (T7 16, T8 41), 128 cases each, on every CI run. About a dozen giac quirks are recorded as departures (one-variable modgcd drops content; ".." inside strings gets a space; mixdown on 8-bit clips; splice crossfade overrun).
- One-hour campaign (Fuzz 45856971637 on 76c1018): 55 of 56 targets clean, peak rss 498 MB of 2048; the eval timeout became CAS-649. Every fuzz run after T4 had been cancelled before finishing, which hid a period crash; the replay mode -f replay=true now exists, and a fuzz run counts only if it completes. fuzz.yml's job cap is 120 min (was 60).
- Cleanup: T4 prototypes are in crates/cas-core/proto/leftovers/; .run/logs/gen1-gap/ is deleted (R45). Cookbook C0094-C0117.
Open: CAS-650 backlog (eight more functions re-evaluate an argument under fresh limits); for gen 2, markov/lpsolve sweeps and a residue rule for rational Hilbert transforms; Fill inside begin…end; euler_lagrange_assume.rs:45's needless rewrite; publish CAS-DOC-57/58/59; commit_task.sh still sweeps logs/; log-only commits unpushed.
