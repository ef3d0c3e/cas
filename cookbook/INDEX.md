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
C0044 | giac garbage answers change between runs: departure rows check ours only | giac,oracle,departures,nondeterministic,garbage | 2026-11-23 | 1.4/T2 | CAS-589 rombergt on a*x
C0045 | Oracle rows comparing our balls with giac's [lo..hi] intervals | giac,interval,ball,oracle,convert,overlap | 2026-11-23 | 1.4/T3 | T3.c2 CAS-592 convert interval
C0046 | Porting a float iteration to rounded Rationals: round every repeated factor | port,giac,newton,rational,rounding,iteration | 2026-11-23 | 1.4/T3 | T3.c3 CAS-591 newton
C0047 | Decimals in full-form unit tests as hex floats; caseval's ' | ' is a newline | tests,float,decimal,fullform,caseval,departures | 2026-11-23 | 1.4/T4 | T4.c1 CAS-593
C0048 | Recurrences over Expr on symbolic inputs: cap the point count, not only the work | fuzz,timeout,limits,interpolation,symbolic | 2026-11-23 | 1.4/T4 | T4.c3 CAS-595
C0049 | giac piecewise answers in oracle rows: branch pairs, conditions by sampling | giac,oracle,piecewise,conditions,tests | 2026-11-23 | 1.4/T5 | T5.c1 CAS-596 bspline
C0050 | giac symbolic linsolve: free parameters make comparisons skip, port the skip | giac,port,linsolve,prototype,tpsolve | 2026-11-23 | 1.4/T6 | CAS-598 tpsolve zero supply
C0051 | Remez exchange: giac's one-peak-per-zero-interval step fails; use a full exchange | remez,minimax,approximation,giac,port | 2026-11-23 | 1.4/T7 | CAS-599
C0052 | Eval fuzz plain vs traced disagree at a limit: cache order, not a bug | fuzz,eval,harness,recursion-limit,cache,stamps,triage | 2026-11-24 | 1.4/T9 | Fuzz 45285948112, T9.c1 2a70d96
C0053 | Per-cache-hit side data (warnings, trace steps) needs its own count bound; OOM triage via addr2line | fuzz,oom,cache,memo,warnings,trace,limits,addr2line | 2026-11-24 | 1.4/T9 | Fuzz 45293145351, T9.c2
C0054 | Subresultant gcd: gamma/lc-scale the last remainder before its primitive part | gcd,subresultant,cas-poly,fuzz,timeout,content,performance | 2026-11-24 | 1.4/T9 | Fuzz 45304609232, T9.c5 8a82c52
