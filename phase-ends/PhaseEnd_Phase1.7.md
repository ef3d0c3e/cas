# PhaseEnd — Phase 1.7: Statistics III and signal processing II (implements Gen 1.7)
Approved: 2026-11-25 | Closed: 2026-11-26 | Planner: claude-opus-5-5/medium | Tasks: 12 done, 0 superseded

## Milestone
Milestone: giac's remaining computational statistics commands and its in-memory sound and signal commands (the 37 aide entries of ## Context, with hilbert's list form) run under their giac names, each with an oracle test against giac, a property test where giac draws at random or cannot answer, or a recorded departure, with its doc and a fuzz entry; CAS-DOC-59 records the area and CAS-521..524; CI and Fuzz pass on the pushed master head and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-59, one ticket per family and the diff(e,x,n) ticket | Done: CAS-DOC-59 "Statistics III and signal processing II" (draft) holds the phase's design; tickets CAS-624..635 (one per family) and CAS-636 (Calculus, diff(e,x,n)) exist, status todo. Summary: The design page lists the 37 giac entries plus hilbert's list form by task and family, marked command, option word, alias or departure, with what stays out and why. It records what the phase builds on (including the four signal tickets that never had a page), one section per build task with giac's sources, the settled and open design choices, the oracle plan and giac's own examples, and a closing section on how answers, seeded draws and departures are checked. Twelve family tickets and one Calculus ticket were filed. | commit: - | tasks/T1.md | logs/T1.md
- T2 — random variables: randvar and random_variable, sample over them, randmultinomial | Done: randvar (random_variable) answers a distribution object for every aide form. sample, randvector and randmultinomial draw from objects through the seeded generator. CAS-624 is done, master is pushed, and CI is green on 2488b71. Summary: Xcas's random variables now exist: randvar turns a law's name and parameters, moment options or weight lists into one object that sample, randvector and later tasks' cdf, icdf, mgf and fitdistr read. A named law is its density's call with only its parameters, as giac prints it. The binomial and the discrete laws are randvar's own canonical call, so every object prints and reads back with no new kind of value. Every draw is seeded, so srand repeats a run. Two latent bugs surfaced on the way and are fixed: a giac-language reading of `mean(…mean=5…)` as a product, and a Simplify walk that went exponential on shared expressions once new builtins shifted symbol bits. | commit: 2488b71 | tasks/T2.md | logs/T2.md
- T3 — tail probabilities, distribution functions and the logistic | Done: UTPN, UTPT, UTPC, UTPF, cdf, icdf, mgf, logistic and logistic_regression answer as giac, exact where closed; CAS-625..627 done; master pushed and CI green on ba84d6e. Summary: Xcas's upper-tail probabilities, its generic cdf/icdf/mgf over a law's name, a randvar object or a data list, the logistic function and its closed-form regression now exist. Tails are one minus our existing closed-form cumulative functions, so UTPC(2,6) is exactly exp(-3). The empirical cdf of a list deliberately differs from giac, whose answer is off by one value. A latent regression bug surfaced: decimal data left an unevaluated Together around ball slopes; it is fixed for every regression. | commit: ba84d6e | tasks/T3.md | logs/T3.md
- T4 — Markov chains: markov and randmarkov | Done: markov(M), markov(M, eps) and randmarkov's path, multi-start and random-matrix forms exist, seeded and limit-checked; CAS-628 is done, master is pushed and CI is green on e6fac66. Summary: Xcas's Markov chain commands now exist. markov finds the closed classes, their long-run probabilities, the passing states and the chance of ending in each class, listed in the same order as giac. It solves for these numbers exactly instead of iterating, so exact input gets exact answers (1/2 where giac prints a 13-digit fraction), and chains where giac's loop never ends still get an answer. randmarkov walks a chain or builds a random one through the seeded generator. | commit: e6fac66 | tasks/T4.md | logs/T4.md
- T5 — density estimation and fitting: kernel_density kde, bandwidth, bins, fitdistr | Done: kernel_density (kde) and fitdistr exist, with bandwidth and bins as option words; CAS-629 is done, master is pushed and CI is green on 94b71d8. Summary: Xcas's Gaussian kernel density estimate now answers as a list of densities on giac's grid, as an exact formula in x, or as giac's piecewise linear form. The bandwidth is given, giac's plug-in rule, or its rule of thumb. fitdistr fits eight laws to data and answers T2's distribution objects: exactly where the estimator is a closed form, and as certified numbers for the four iterative fits. giac's own mistake in the mean (dividing by N-1) is not copied; where it changes giac's answer the row is a recorded departure. | commit: 94b71d8 | tasks/T5.md | logs/T5.md
- T6 — clustering: kmeans and cluster | Done: kmeans (KMeans) and cluster (Cluster) port giac's Hartigan-Wong k-means and hierarchical clustering; CI 45608125956 green on baf6ecb; CAS-630 done. Summary: Two new statistics commands group points into clusters the way Xcas does. kmeans picks random starting centres from the session's seeded generator, so its answer can differ from Xcas's where the data has several equally good groupings; those cases are checked against the set of groupings Xcas's algorithm can reach. cluster is deterministic and matches Xcas's labels exactly. | commit: baf6ecb | tasks/T6.md | logs/T6.md
- T7 — sound clips: createwav, soundsec, samplerate, bit_depth, channels, duration, channel_data and about | Done: giac's sound clips exist as values (SoundClip[rate, depth, codes…]) with createwav, soundsec, samplerate, bit_depth, channels, duration, channel_data and about; CI 45614349596 green on 30e541e; CAS-631 done. Summary: Sound clips can now be built and read back the way Xcas does, with every sample exact (for example 7301/8192 where giac prints 0.891235351562). A clip is an ordinary expression that stores giac's integer codes and prints as a readable call, so the session, string() and the printers all round-trip it; giac's one-line description, its clip arithmetic and its indexing are recorded as departures. about(x) now reports the facts assume recorded for x. | commit: 2ad4efd | tasks/T7.md | logs/T7.md
- T8 — clip edits and resampling: set_channel_data, stereo2mono, mixdown, splice, resample | Done: set_channel_data, stereo2mono, mixdown, splice and resample answer clips as Xcas's (resample of our own design for the sinc qualities); CI 45621233797 green on 1d5f768; CAS-632 and CAS-633 done. Summary: Sound clips can now be edited the way Xcas edits them: channel data overwritten, stereo folded to mono, several clips mixed with offsets, gains, panning and fades, and two clips spliced with a crossfade, each matching giac's printed samples code for code. Resampling works too, which Debian's giac cannot do at all: the two simple qualities follow libsamplerate's linear and hold converters exactly, and the three high qualities use our own windowed-sinc filters whose passband and stopband errors were measured and are documented. | commit: 900f907 | tasks/T8.md | logs/T8.md
- T9 — signal transforms: dwt, idwt, hilbert's list form, instfreq, instphase, emd | Done: dwt and idwt (GSL's Daubechies, Haar, B-spline; centred, 2-D standard and image; padding), hilbert of a list, instfreq, instphase and emd exist with unit tests, oracle rows and departures; CI green on c6f3449; CAS-634 and CAS-635 closed. Summary: The wavelet transforms port GSL's filters and pyramid algorithm the way giac calls them, answering certified decimals, and exact values for the Haar wavelet on exact data. The analytic-signal family computes hilbert of a list through our FFT, the instantaneous frequency and phase of lists and of expressions, and giac's empirical mode decomposition with the same sifting decisions as giac but bounded by our limits. The first CI run caught a port slip in hilbert for odd lengths, fixed in the task. | commit: c6f3449 | tasks/T9.md | logs/T9.md
- T10 — push master for the phase's final CI and Fuzz runs | Done: master 476ba0c pushed to origin (c6f3449..476ba0c); CI run 45629972224 and Fuzz run 45629972139 started on it. Summary: The local gate (format, stable clippy, the fuzz crate's check) passed on the committed head, and master was pushed. Both the regular CI run and the phase's long Fuzz run started on that commit; T11 reads their results. | commit: 476ba0c | tasks/T10.md | logs/T10.md
- T11 — CI and Fuzz check on the pushed head; tickets done, CAS-DOC-59 in step, the eval target's memory read | Done: CI 45637752766 and Fuzz 45637752779 green on 0e22fdf (all jobs); three eval fuzz findings of the phase fixed with regression tests; CAS-624..635 done; CAS-DOC-59 §9 lists the README's departures. Summary: The phase's earlier fuzz runs had all been cancelled by later pushes, but three of them had already caught problems in the eval target, and the final run on the old head caught one more. Two were the fuzz harness comparing two evaluations of a random Markov chain, which draw different numbers; the harness now reseeds before each evaluation it compares. The third was a low-pass filter nested in its own cutoff, which took time exponential in the depth; filters now approximate only numeric parameters. After the fix push, CI and the full 47-minute fuzz run are green, and the eval target's peak memory (128 MB) is close to the previous phase's (135 MB), far below the 2 GB limit. | commit: - | tasks/T11.md | logs/T11.md
- T12 — phase cleanup | Done: the phase's cited prototypes and giac probe inputs live in crates/cas-core/proto/ with a rerun README; all repo citations point there. Summary: The tests and module docs written this phase named prototype scripts and giac probe files that lived only in the phase's scratch folder, which the close deletes. They are now in the repo beside the crate that cites them. A short guide there explains how to rerun them, and every reference points at the new place. Each prototype, run from its new place, reproduces what it produced before. | commit: ad44740 | tasks/T12.md | logs/T12.md

## Decisions that still bind
- T2 — the distribution object (cross-task contract for T3, T5): a named law is the parameter-only density call that stays unevaluated, at exactly these arities: NormalPDF[μ,σ], ExponentialPDF[λ], UniformPDF[a,b], PoissonPDF[λ], GeometricPDF[p], NegativeBinomialPDF[n,p], StudentTPDF[ν], ChiSquarePDF[ν], FisherPDF[m,n], GammaPDF[a,b], BetaPDF[a,b], WeibullPDF[k,λ], CauchyPDF[x₀,γ]. The binomial law is RandomVariable[Binomial,n,p] (Binomial[n,p] is the coefficient). A discrete law is RandomVariable[{w…}] (draws indices 0…k−1) or RandomVariable[{w…},{v…}]. randvar leaves its canonical forms as they are (a fixed point), and randvar(X) of an object is X. Read objects only through `random_variables::variable(e) -> Option<Variable>` (Law(Law,&[Expr]) | Discrete{weights, values}) and draw through `draw(ev,&Variable)`. Routed: the module doc is the product doc; its unit tests pin it.
- T7 — a sound clip is `SoundClip[rate, depth, {codes ch0}, {codes ch1}, …]` (surface `soundclip`, no native); codes signed, |k| ≤ 2^(b-1)-1, sample k/2^(b-1); well-formed check `sound::clip(e) -> Option<Clip>` (pub(in crate::native)) is what T8's natives take and give.
- T7 — the clip prints as that call everywhere (session, string(), hook, LaTeX); lang_session and lang_print hold because it is an ordinary named call.
- T12 — phase prototypes cited by tests live in crates/<crate>/proto/<area>/ (.py and .cas only; giac's .out regenerated with scripts/caseval.py, gitignored); card Paths updated.
-

## Rules proposed
- T7 — none
- T8 — none

## Cookbook entries added
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
C0092 | Natives re-approximating symbolic arguments: is_numeric before numeric::approximation | fuzz,timeout,native,approximation,numeric,memo,exponential | 2026-11-25 | 1.7/T11 | Fuzz 45608125834, T11.c1 501c211
C0093 | Decoding an old eval fuzz artifact with its run's harness; reseeding for draws | fuzz,eval,artifact,decode,replay,random,seed,harness | 2026-11-25 | 1.7/T11 | Fuzz 45608125834, 45601609896, T11.c1 951ab6a

## Research
R1.7-001 | verify registry lines and answer (b)-(e) for the clip value / about(a) plan | Phase 1.7 sound objects: registry, printer, fuzz, assumptions, departures (verified at HEAD) | registry, builtins, print, fuzz, assume, departures, phase-1.7 | retriever-code | 2026-11-25 | 78 lines
R1.7-002 | phase 1.7 clip-operations and signal-transforms tasks (giac-side behaviour, port sizes) | giac 2.0.0.21 clip ops and signal transforms: port map for phase 1.7 | phase-1.7, giac, audio-clip, splice, mixdown, resample, dwt, emd, instfreq, instphase | retriever-code | 2026-11-25 | 149 lines
R1.7-003 | phase 1.7 plan ## Interfaces and statistics task split (giac markov, randmarkov, kmeans, cluster, kernel_density, fitdistr, logistic_regression, UTP*, cdf, icdf, mgf, randvar, randmultinomial) | Giac statistics plan: existing statistics, matrix, test and doc interfaces | statistics, distributions, random, matrix, giac-tests, fuzz, docs | retriever-code | 2026-11-25 | 72 lines
R1.7-004 | phase 1.7 ## Interfaces and statistics split (kmeans, cluster, kde, fitdistr, markov, randmarkov, mgf, randvar, UTPC/N/T/F, cdf, icdf, logistic_regression, randmultinomial) | giac statistics and random natives: where the phase 1.7 plan attaches | giac, statistics, distributions, random, phase-1.7 | retriever-code | 2026-11-25 | 56 lines
R1.7-005 | phase 1.7 statistics tasks; giac-side behaviour, argument forms, algorithms and probe cross-checks for kernel_density, fitdistr, kmeans, cluster, markov, randmarkov, mgf, UTP*, cdf/icdf, randmultinomial, discreted, randvar, sample, logistic_regression. | giac 2.0.0.21 statistics commands: port spec for phase 1.7 | giac, statistics, kde, fitdistr, kmeans, cluster, markov, randmarkov, randvar, sample, phase-1.7 | retriever-code | 2026-11-25 | 243 lines
R1.7-006 | phase 1.7 plan, ## Interfaces for the three signal-processing tasks | Signal-processing code anchors for phase 1.7 natives (giac sound clips, dwt, idwt, emd, instfreq, instphase) | signal, natives, fft, windows, filters, correlation, threshold, moving-average, tests, fuzz, registry | retriever-code | 2026-11-25 | 64 lines
R1.7-007 | which R0-004 names for Statistics and Signal processing are still missing, and which other missing entries fall in these fields by description | R0-004 areas 9 and 16 vs missing_described (gap 2026-11-25), phase 1.7 scope | giac, gen1, gap, phase-1.7, statistics, signal, areas | retriever-digest | 2026-11-25 | 42 lines
R1.7-008 | phase 1.7 Rust CAS, giac resample(clip, rate, quality) mapping quality 0..4 to libsamplerate converters | libsamplerate converter algorithms, src_simple length, end handling, license | libsamplerate, resample, giac, sinc, license, phase-1.7 | retriever-web | 2026-11-25 | 62 lines
R1.7-009 | input for CAS-DOC-59 (T1 design page): one §n per build task T2-T9 and §9 Oracle and tests | Phase 1.7 PHASE_PLAN.md: Rationale, Risks and T2-T9 design facts | phase-1.7, statistics, signal-processing, giac-oracle, design-page | retriever-digest | 2026-11-25 | 116 lines
R1.7-010 | Rust port spec for T3 (giac 2.0.0.21 source, read-only) | Port spec: giac moyal.cc UTPN/UTPT/UTPC/UTPF, mgf, cdf, icdf; misc.cc logistic_regression; logistic | giac, moyal, distributions, utp, mgf, cdf, icdf, logistic, port-spec | retriever-code | 2026-11-25 | 159 lines
R1.7-011 | T5 of phase 1.7 (port kernel_density/kde and fitdistr to Rust; Python prototype must reproduce giac numbers) | giac 1.7 kernel_density / DPI bandwidth / fitdistr: exact formulas for the T5 port | giac, kernel_density, kde, dpi, bandwidth, fitdistr, weibull, cauchy, beta, gamma, source-read | retriever-code | 2026-11-25 | 80 lines
R1.7-012 | T6.c1 | giac 1.7 kmeans (AS 136) and cluster (hclust): port notes, seed dependence, expected answers for T6 | giac,kmeans,cluster,hclust,AS136,linkage,ward,seed,source-read,prototype | pa-coder | 2026-11-25 | 80 lines
R1.7-013 | T9 port of giac emd to Rust; Python prototype must reproduce giac output | giac emd: port spec for Python prototype | emd, gsl, cspline, signalprocessing, port-spec | retriever-code | 2026-11-25 | 93 lines

## Audit
- seed: median 24k, max 24k, n=12 (phase 1.7)
- previous phase 1.6: median 24k, growth -0.2%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 8209 bytes
- cookbook/INDEX.md: 16366 bytes
- rules/INDEX.md: 3242 bytes
### Carry audit — phase 1.7 (2026-11-25T06:35:05Z → open UTC, 1 sessions, 1899 requests)
- beside 1.7: session 659cf258 $0.00 (0 turns), not attributed

| file (read) | chars | n |
|---|---|---|
| PHASE_PLAN.md | 221.3k | 16 |
| print.rs | 92.4k | 21 |
| signalprocessing.cc | 83.6k | 28 |
| R1.7-005.md | 79.3k | 9 |
| bi2soamxq.txt | 76.2k | 4 |
| optimization.cc | 71.8k | 10 |
| context.rs | 71.6k | 31 |
| bz00ry02t.txt | 68.5k | 4 |
| sparse.rs | 61.9k | 4 |
| markup.cc | 55.4k | 16 |
| file (write) | chars | n |
|---|---|---|
| random_variables.rs | 58.8k | 3 |
| page.py | 53.0k | 2 |
| retriever-rationale-t2-t9-design.md | 47.0k | 1 |
| clustering.rs | 40.6k | 3 |
| retriever-phase-plan-rationale-t2-t10.md | 36.1k | 1 |
| proto.py | 32.1k | 4 |
| retriever-code-moyal-misc-port-spec.md | 29.9k | 2 |
| density.rs | 29.0k | 3 |
| result kind | chars | n |
|---|---|---|
| Read | 4.6M | 1015 |
| bash other | 3.4M | 1442 |
| Grep | 1.5M | 669 |
| tools/plan_edit.py | 130.5k | 35 |
| run.sh | 112.4k | 132 |
| Glob | 60.5k | 42 |
| tools/outline.py | 51.9k | 21 |
| Agent | 43.8k | 43 |
| Write | 25.5k | 156 |
| tools/card.py | 25.4k | 3 |
- spilled: 15 results, 546.5k chars on disk, 3 read whole
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,193 (n 65, prev 6,193) · retriever-digest 6,168 (n 11, prev 6,168) · retriever-web 4,581 (n 3, prev 4,581)
- outline credit: 24 outlines · 12 followed by a ranged read · 381.7k chars credited
- warm pings: 16 pings over 3 runs, 2 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.3533 + relay $0.0000 + re-arm $0.0000 = $0.3533 vs rewrites replaced $1.9209, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 761, relay 0, re-arm 35, other 35, relay cost $0.0000
- retriever re-asks: 0 of 8 retriever briefs repeat a lookup of the same run
- noise: 39 lines 2.4k chars — No such file or directory: 28 lines/1.9k, usage: : 11 lines/508
- price: read $0.20/Mtok · 1h write $7.90/Mtok · output $19.76/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 5287 chars, 1899 requests (prev 8040 chars, 1687 requests, growth -34.2%)
- flag: 3 spilled results read whole
  command: P=/usr/bin/python3; $P tools/plan_edit.py show --section Context; $P tools/plan_  reader: pa-expert
  command: P=/usr/bin/python3; $P tools/plan_edit.py show --section Context; $P tools/plan_  reader: pa-expert
  command: P=/usr/bin/python3; $P tools/plan_edit.py show --section Context; $P tools/plan_  reader: pa-expert
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
- T7 | pa-claude | partner | claude-opus-5-5 | medium | ctx 109k | $1.880 | saved $17.11 | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 159k | $2.845 | saved $0.83 | completed | parent 1b3aa32d
- T1 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 76k | $0.031 | saved - | completed | parent a7f2de09
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 205k | $5.502 | saved $17.40 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 98k | $1.743 | saved - | completed | parent a2995941
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 171k | $3.485 | saved - | completed | parent a2995941
- T2 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 32k | $0.010 | saved - | completed | parent a0085a75
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 47k | $0.634 | saved - | completed | parent a2995941
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 87k | $1.307 | saved - | completed | parent a2995941
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 122k | $2.240 | saved $9.93 | completed | parent 1b3aa32d
- T3 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 66k | $0.028 | saved - | completed | parent aaa29a23
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 187k | $3.983 | saved - | completed | parent aaa29a23
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 123k | $2.342 | saved - | completed | parent aaa29a23
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 53k | $0.665 | saved - | completed | parent aaa29a23
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 100k | $1.700 | saved $11.26 | completed | parent 1b3aa32d
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 206k | $5.240 | saved - | completed | parent a003f889
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 157k | $3.510 | saved $18.82 | completed | parent 1b3aa32d
- T5 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 64k | $0.029 | saved - | completed | parent a5144bbb
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 130k | $2.515 | saved $0.14 | completed | parent a5144bbb
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 220k | $4.823 | saved $0.14 | completed | parent a5144bbb
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 202k | $4.289 | saved - | completed | parent a5144bbb
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 65k | $0.899 | saved - | completed | parent a5144bbb
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 125k | $2.362 | saved $13.51 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 99k | $1.624 | saved $0.12 | completed | parent a7477801
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 213k | $4.871 | saved - | completed | parent a7477801
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 183k | $4.113 | saved $11.11 | completed | parent 1b3aa32d
- T7 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 16k | $0.005 | saved - | completed | parent a6711035
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 80k | $1.176 | saved - | completed | parent a6711035
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 161k | $3.252 | saved - | completed | parent a6711035
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 108k | $1.829 | saved - | completed | parent a6711035
- T7 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 33k | $0.011 | saved - | completed | parent a6fe1881
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 116k | $2.343 | saved $14.14 | completed | parent 1b3aa32d
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 181k | $4.070 | saved - | completed | parent aec99563
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 145k | $3.304 | saved - | completed | parent aec99563
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 135k | $3.014 | saved $15.48 | completed | parent 1b3aa32d
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 40k | $0.018 | saved - | completed | parent a5bb2ada
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 199k | $5.036 | saved - | completed | parent a5bb2ada
- T9 | pa-coder | coder | claude-opus-5-5 | medium | ctx 150k | $3.037 | saved $0.00 | completed | parent a5bb2ada
- T9 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 43k | $0.013 | saved - | completed | parent a7605bee
- … 7 more: PY ~/.claude/pa3/pa_ledger.py report

## Discussions
- I1 | -> Gen 2: frontends, decided 2026-11-25. First a language server over the CAS: completion, highlighting (an editor gramm | deferred:gen2
- I2 | -> 1.8: Weibull fit (random_variables.rs, T2) ports giac's equation, which drops a square (mean 12.5, var 1 gives a law | deferred:1.8

## Deferred
- from T12: the close may delete .run/1.7/; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)
- I1 (deferred:gen2): -> Gen 2: frontends, decided 2026-11-25. First a language server over the CAS: completion, highlighting (an editor gramm
- I2 (deferred:1.8): -> 1.8: Weibull fit (random_variables.rs, T2) ports giac's equation, which drops a square (mean 12.5, var 1 gives a law

## Changes
<!-- appended by plan_edit.py append-change; never hand-edited -->
- 2026-11-25 router: T1 next -> done
- 2026-11-25 router: T2 next -> done
- 2026-11-25 router: T3 next -> done
- 2026-11-25 router: T4 next -> done
- 2026-11-25 router: T5 next -> done
- 2026-11-25 router: T6 next -> done
- 2026-11-25 router: T7 next -> done
- 2026-11-25 router: T8 next -> done
- 2026-11-25 router: T9 next -> done
- 2026-11-25 router: T10 next -> done
- 2026-11-25 router: T11 next -> done
- 2026-11-25 router: T12 next -> done

## Plain-English Recap
Phase 1.7 Statistics III and signal processing II: CAS-DOC-59 (still a draft in Lific) covers 37 giac entries; tickets CAS-624..635 are done, and CAS-636 (diff(e,x,n) doesn't take the nth derivative) is filed and stays todo. Built:
- Random variables: randvar/random_variable. A named law is its density call left unevaluated (normald(5,3/2)); the binomial is randvar(binomial,n,p); discrete laws are given as lists. sample, randvector and randmultinomial use the seeded generator.
- Probabilities and fits: UTPN/T/C/F (exact where possible, UTPC(2,6)=exp(-3)), cdf, icdf, mgf, logistic, logistic_regression. cdf of a data list gives count(≤x)/n, a departure from giac. kernel_density/kde and fitdistr: closed forms are exact; density lists compare within giac's FFT noise; the gamma fit is a departure.
- markov and randmarkov, exact on exact input. kmeans (Hartigan-Wong) and cluster (hierarchical).
- Sound clips: SoundClip expressions holding exact integer codes; createwav, the property commands, about, set_channel_data, stereo2mono, mixdown, splice and resample. Resample qualities 2-4 are our own windowed-sinc filters, stopband 102-147 dB.
- Signal transforms: dwt/idwt (GSL filters), hilbert of a list, instfreq, instphase, emd (bit for bit with giac).
Fixes outside the plan: the reader's f(x)→f*x rule skips function names; a simplify walk that could take 2^80 steps; the fuzz harness reseeds before each compared evaluation; nested LowPass no longer grows exponentially. Final gate green on 0e22fdf (CI 45637752766, Fuzz 45637752779). Eval peak memory 128 MB of 2048, so the symbol_array memory watch is closed. The cited prototypes and probes moved to crates/cas-core/proto/. Cookbook C0077-C0093.
Open, carried to 1.8 in the inbox:
- giac's Weibull fit drops a square, and ours copies it. I recommend the correct fit plus a departure row.
- kmeans/cluster output=part, index= and string clustering, and hilbert(f,t), go to the gap check.
- Coder briefs restate R74 (broken twice this phase, in T6 and T7; nothing from those runs was used).
Also open: publish CAS-DOC-57/58/59; commit_task.sh still sweeps logs/ (R73/R96); docs-only commits stay unpushed until the next phase's first push.
