# GenerationEnd — Generation 1: Green master, then giac parity in the existing modules
Closed: 2026-11-27 | Phases: 8 (8 closed)

## Phases
- 1.1 | Green master | milestone: Master's CI (fmt, clippy, test, fuzz-check, oracle) and Fuzz (every target) pass on the developer's push of two fixes, the eval target's divide-by-zero under a regression test in its owning crate and `boundary_values` checking BVPSolve against the closed form, and CAS-517..524 are done in Lific; `gh run list --branch master --limit 2` shows both runs green. | closed | phase-ends/PhaseEnd_Phase1.1.md
- 1.2 | Algebra and polynomials parity | milestone: The computational giac commands of the Algebra and Polynomials modules that this CAS lacked run under their giac names, each with a giac oracle test or a recorded departure, its doc chapter and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.2.md
- 1.3 | Calculus, expressions and integers parity | milestone: The computational giac commands of the Calculus, Expressions and Integers modules (list commands included) that this CAS lacked run under their giac names, each with a giac oracle test or a recorded departure, its doc chapter and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.3.md
- 1.4 | Numerical analysis II | milestone: giac's quadrature rules, interpolation and approximation commands, numerical derivatives and named root finders run under their giac names, with certified balls where giac gives floats, each with an oracle test against giac or a closed form and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.4.md
- 1.5 | Optimization II | milestone: lpsolve and nlpsolve take giac's options, and the optimisers simplex_reduce, frank_wolfe, simulated_annealing and find_minimum run under their giac names, each with an oracle test against giac or known optima and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.5.md
- 1.6 | Strings, printers and session commands | milestone: giac's string commands, the latex and mathml printers and the session and settings commands this CAS lacked run under their giac names, each with an oracle test or a recorded departure and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.6.md
- 1.7 | Statistics III and signal processing II | milestone: giac's remaining computational statistics commands and its in-memory sound operations run under their giac names, each with an oracle test or a recorded departure and a fuzz entry; CI and Fuzz pass on the developer's push and the area's Lific tickets are done (`gh run list --branch master --limit 2`). | closed | phase-ends/PhaseEnd_Phase1.7.md
- 1.8 | Validation | milestone: A one-hour-per-target fuzz campaign and the oracle's random-input sweeps over the generation's builtins pass on the developer's push with every finding fixed under a regression test, and `bash scripts/giac_gap.sh` leaves no computational giac command of an existing module unaccounted for (built, or a departure recorded in crates/cas-giac/README.md). | closed | phase-ends/PhaseEnd_Phase1.8.md

## Decisions that still bind
-
- T12 — giac lookups outside tests use scripts/caseval.py and scripts/aide.py (card Paths line); planners stop rewriting them in .run/<phase>/.
- T3 — giac's random Newton restarts become a fixed start sequence start_j = (-1)^floor((j-1)/2)·j·3/4 (j = 1..20; giac's `j > 5` → undef kept); rows whose giac answer depends on random restarts compare ours with a closed form and giac's only as a root.
- T9 — `Trace` has a public `truncated: bool`; past `MAX_STEPS` (1<<20) steps the trace is a non-replayable prefix. Route: cas-core README trace row (done) and CAS-DOC trace contract at phase end.
- T2 — positions and lengths count Unicode scalar values (strings handled as Vec<char>); giac's byte answers are departure rows.
- T3 — the hook is part of the context state: set only through `set_language` (state_mut → new StateId); every T4-T7 form or reader goes through `Language` + a new `Form` variant, never a Context field outside State (C0066).
- T3 — natives print with the `terms` limit as the byte cap (print_capped); a print or read failure leaves the call with a `Warning::Native` carrying our text; no hook → Ok(None) only where a print/read is actually needed (cat of strings, expr("") still answer).
- T4 — T5 extends the same engine (markup/engine.rs MarkupBlock + Flags) with MathML fields and flags; the lowering (markup/lower.rs) is shared; a Form::MathML variant goes through hook.rs like Form::Latex.
- T5 — MathML presentation and content share T4's engine pass (one MarkupBlock, giac's design); mathml.cc's older layout is a separate cas-lang module over the same lowering.
- T6 — every setting lives in Context State (C0066); a setting's call form `name()`/`name(v)` is a session arm, never a native, when its bare name has a value (C0070).
- T7 — a native reading data kept outside State calls `Eval::volatile()` before reading (C0072).
- T2 — the distribution object (cross-task contract for T3, T5): a named law is the parameter-only density call that stays unevaluated, at exactly these arities: NormalPDF[μ,σ], ExponentialPDF[λ], UniformPDF[a,b], PoissonPDF[λ], GeometricPDF[p], NegativeBinomialPDF[n,p], StudentTPDF[ν], ChiSquarePDF[ν], FisherPDF[m,n], GammaPDF[a,b], BetaPDF[a,b], WeibullPDF[k,λ], CauchyPDF[x₀,γ]. The binomial law is RandomVariable[Binomial,n,p] (Binomial[n,p] is the coefficient). A discrete law is RandomVariable[{w…}] (draws indices 0…k−1) or RandomVariable[{w…},{v…}]. randvar leaves its canonical forms as they are (a fixed point), and randvar(X) of an object is X. Read objects only through `random_variables::variable(e) -> Option<Variable>` (Law(Law,&[Expr]) | Discrete{weights, values}) and draw through `draw(ev,&Variable)`. Routed: the module doc is the product doc; its unit tests pin it.
- T7 — a sound clip is `SoundClip[rate, depth, {codes ch0}, {codes ch1}, …]` (surface `soundclip`, no native); codes signed, |k| ≤ 2^(b-1)-1, sample k/2^(b-1); well-formed check `sound::clip(e) -> Option<Clip>` (pub(in crate::native)) is what T8's natives take and give.
- T7 — the clip prints as that call everywhere (session, string(), hook, LaTeX); lang_session and lang_print hold because it is an ordinary named call.
- T12 — phase prototypes cited by tests live in crates/<crate>/proto/<area>/ (.py and .cas only; giac's .out regenerated with scripts/caseval.py, gitignored); card Paths updated.

## Deferred
- from T8: the close may delete `.run/1.1/` whole. (unconsumed)
- from T12: the close deletes .run/1.2/; leave .run/logs/gen1-gap/ (GENERATION_PLAN.md:57). (unconsumed)
- from T12: the close deletes .run/1.3/; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)
- from T10: the close may delete .run/1.4/ whole. (unconsumed)
- from T10: the close may delete .run/1.5/ whole. (unconsumed)
- from T11: the close may delete .run/1.6/ whole; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)
- from T12: the close may delete .run/1.7/; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)
- I1 (deferred:gen2): -> Gen 2: frontends, decided 2026-11-25. First a language server over the CAS: completion, highlighting (an editor gramm
- I2 (deferred:1.8): -> 1.8: Weibull fit (random_variables.rs, T2) ports giac's equation, which drops a square (mean 12.5, var 1 gives a law
- from T12: the close deletes .run/1.8/ whole; nothing in it is referenced anymore. (unconsumed)

## Plain-English Recap
{{AUTHORED:recap}}
