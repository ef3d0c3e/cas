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
C0055 | Euclid over Z[i] (and other field-backed rings): round the quotient in the ring from leading bits | gcd,gaussian,euclid,cas-num,fuzz,timeout,performance,rounding | 2026-11-24 | 1.4/T9 | Fuzz 45311134759, T9.c6 8c9ca3e
C0056 | Stays tests check head and warnings; feasibility checks split equalities | test,unevaluated,oracle,equation,cas-core | 2026-11-24 | 1.5/T2 | T2.c1 CAS-600 lpsolve stays
C0057 | giac branch and bound port: GMI cuts by slack substitution and re-solve; giac's cut rhs bug | giac,port,lpsolve,integer,cuts,branch-and-bound | 2026-11-24 | 1.5/T3 | T3.c1 CAS-602
C0058 | Expr::may_contain is a bloom filter: never assert its negation | tests,expr,symbols,builtins,flaky | 2026-11-24 | 1.5/T4 | CI 45337180530 trigonometric_interpolation
C0059 | Certified KKT points under bounds: one-point ranges; simple inequalities as bounds | nlpsolve,kkt,certified,bounds,branch-and-bound | 2026-11-24 | 1.5/T4 | T4.c2 CAS-605; CI 45337180530 nlpsolve_seeded
C0060 | giac gen arithmetic on mixed exact/float data: an Exact|Float entry enum | giac,port,decimal,exact,prototype | 2026-11-24 | 1.5/T5 | CAS-606
C0061 | Certified global minimiser on an interval: float candidate, ball cover, interval Newton | minimum,certified,krawczyk,cover,brent,find_minimum | 2026-11-24 | 1.5/T6 | T6.c1 CAS-608
C0062 | Iterative float optimisers with a certified stop: test in balls, one shared work budget | frank-wolfe,certified,optimisation,limits,budget,gap | 2026-11-24 | 1.5/T6 | T6.c2 CAS-609
C0063 | A native applying user functions: Expr::apply + ev.eval, N[] for numbers | native,function,apply,eval,random,seed | 2026-11-24 | 1.5/T7 | T7.c1 CAS-610
C0064 | Iterated user steps on balls: exact midpoints, or radii compound | ball,radius,random,iteration,native,annealing | 2026-11-24 | 1.5/T7 | T7.c2 CAS-610; CI 45347780940
C0065 | giac string order is byte order: compare with str::cmp, not Expr's Ord | giac,strings,sort,order,compare | 2026-11-25 | 1.6/T2 | T2.c2 CAS-612 sort on strings
C0066 | A Context input a cached result depends on: a State field set via state_mut | context,cache,stamps,state-id,hook,settings | 2026-11-25 | 1.6/T3 | T3.c1 725956e
C0067 | giac printed answers hide internal shape (inv vs division): probe latex(quote(typed)) | giac,port,prototype,quote,latex,markup,shape | 2026-11-25 | 1.6/T4 | T4.c1 CAS-615
C0068 | Printer departure rows under R74: prototype printer on our lowered shape | giac,port,prototype,departures,printer,markup,mathml | 2026-11-25 | 1.6/T5 | T5.c2 CAS-616
C0069 | giac quoted products nest binary; mathml's aplatir_fois_plus multiplies them out | giac,port,prototype,quote,mathml,shape | 2026-11-25 | 1.6/T5 | T5.c3 CAS-616
C0070 | Settings with a bare value: name() goes through the session, not a native | settings,session,eval,head,native | 2026-11-25 | 1.6/T6 | T6.c1 settings
C0071 | Mode rerouting between natives: split the body; degree-mode trig traps | settings,native,recursion,factor,trig | 2026-11-25 | 1.6/T6 | T6.c2 complex_mode, angle_radian
C0072 | Natives reading data outside State: Eval::volatile, no stamp or memo; Context Debug fields | eval,memo,cache,native,context,clippy,session | 2026-11-25 | 1.6/T7 | T7.c1 CAS-620
C0073 | Writing a cas-core native: Ball path, View matches, Eval::context lifetime | native,cas-core,ball,view,borrow,rust | 2026-11-25 | 1.6/T8 | T8.c1 CAS-622
C0074 | giac's descending sorts are reversed stable ascending (ties reversed) | giac,sort,descending,ties,port | 2026-11-25 | 1.6/T8 | T8.c2 CAS-622
C0075 | giac uppercase HP/TI names are not always aliases (SVL=svl, SVD!=svd): probe both | giac,alias,svd,svl,hp,ti,port | 2026-11-25 | 1.6/T8 | T8.c3 CAS-623
C0076 | New builtins: ASCII-alnum internal names; LISTABLE never broadcasts a scalar; decimal rows within error | builtin,listable,registry,naming,oracle,decimal | 2026-11-25 | 1.6/T8 | T8.c4 CAS-621 CAS-622
C0077 | giac error texts from caseval.py: ' | ' is libgiac's newline on CI | giac,oracle,caseval,errors,ci,departures | 2026-11-25 | 1.7/T2 | CI 45565458930 invalid_calls_stay
C0078 | A CI test job that hangs: cancel the run, diff finished tests against the last green run | ci,tests,hang,gh,debugging,libtest | 2026-11-25 | 1.7/T2 | CI 45565458930 shared_subexpressions
C0079 | Walks over an Expr keep a seen set; a may_contain fast path only hides a tree walk | expr,dag,traversal,bloom,builtins,hang,simplify | 2026-11-25 | 1.7/T2 | T2.c4 2488b71 shared_subexpressions
C0080 | Sorting mixed exact and decimal values: one comparator kind per slice | sort,fuzz,float,total_cmp,panic,native | 2026-11-25 | 1.7/T3 | T3.c1 CAS-626 empirical cdf
C0081 | giac's undef/infinity answers are agreeing rows (infinity reads as ComplexInfinity) | giac,oracle,undef,infinity,departures | 2026-11-25 | 1.7/T3 | T3.c2 CAS-627
C0082 | giac ports and clippy: full binding names; Warning::Native needs &'static str | clippy,port,native,warning,lifetime | 2026-11-25 | 1.7/T3 | T3.c1 T3.c2
C0083 | Together leaves balls untouched: test for a ball before wrapping a numeric result | together,ball,native,decimal,regression | 2026-11-25 | 1.7/T3 | T3.c3 CAS-627
C0084 | giac convolution/FFT lists carry packed-FFT noise: compare within it | giac,fft,convolution,oracle,tolerance,kde | 2026-11-25 | 1.7/T5 | T5.c1 CAS-629 prototype
C0085 | giac's discrete choices made in doubles (bin index, trunc at ties): make them in f64 | giac,rounding,ties,exact,oracle,port,kde | 2026-11-25 | 1.7/T5 | CI 45597160314 kernel_density lists
C0086 | Certified MLE fits from giac's starts: ball Newton, Krawczyk, start fallback | mle,fit,newton,krawczyk,certified,ball,weibull,fitdistr | 2026-11-25 | 1.7/T5 | T5.c3 CAS-629 fitdistr
C0087 | giac answers from random starts: measure seed dependence over many prototype seeds | giac,prototype,random,seed,oracle,kmeans,local-optimum | 2026-11-25 | 1.7/T6 | T6.c1 CAS-630 R1.7-012
C0088 | giac's stored doubles keep 47 fraction bits: floors of decimal products drop by one | giac,double,float,prototype,floor,departures,oracle | 2026-11-25 | 1.7/T7 | T7.c1 f3689ae CAS-631
C0089 | Bit-exact giac ports: pow(x,2) folds, pow(1,NaN), bounded memcpy, last-sample flags | giac,port,prototype,f64,nan,sound | 2026-11-25 | 1.7/T8 | CAS-632
C0090 | Kaiser-windowed sinc tables for resampling: sizing H and L, measuring the response | resample,sinc,kaiser,signal,fft,table | 2026-11-25 | 1.7/T8 | CAS-633
C0091 | Porting GSL routines giac calls: tarball fetch, exact literals, ctypes reference | gsl,giac,port,wavelet,spline,prototype | 2026-11-25 | 1.7/T9 | T9.c1 CAS-634, T9.c2 CAS-635
