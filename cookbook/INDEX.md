<!-- Cookbook index. One line per entry, appended by tools/cookbook_add.sh. Grep this file by tag or title; never read it whole; never read an entry you did not find here. -->
<!-- C<nnnn> | <title> | <tags> | <date> | <phase>/<task> | origin: <source section or report id> -->
C0001 | Performance A/B on CI with the Workload workflow | benchmarks,ci,workload,gh | 2026-11-22 | - | gen3.2 memory benchmarks-on-ci.md
C0002 | Diagnosing slow, hanging or OOM fuzz inputs | fuzz,gdb,ci,debugging,timeout,oom | 2026-11-22 | - | gen3.2 memory ci-diagnostics-with-gdb.md
C0003 | Run CI's stable clippy locally (nightly is the default) | clippy,toolchain,ci,lint | 2026-11-22 | - | gen3.2 memory ci-on-github.md
C0004 | gh API rate limit: one slow watcher per run | gh,ci,rate-limit,watch | 2026-11-22 | - | gen3.2 memory ci-on-github.md
C0005 | Locating an eval fuzz panic and replaying it locally | fuzz,eval,panic,replay,cargo-fuzz,backtrace | 2026-11-22 | 1.1/T1 | T1 CAS-525
C0006 | Numeric oracle tests: closed form first, oracle only inside its own error | giac,oracle,numeric,tolerance,bvpsolve,departures | 2026-11-22 | 1.1/T2 | CAS-520 boundary_values
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
