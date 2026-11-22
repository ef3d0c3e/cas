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
