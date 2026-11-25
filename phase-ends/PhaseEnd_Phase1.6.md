# PhaseEnd — Phase 1.6: Strings, printers and session commands (implements Gen 1.6)
Approved: 2026-11-25 | Closed: 2026-11-25 | Planner: claude-opus-5-5/medium | Tasks: 11 done, 0 superseded

## Milestone
Milestone: giac's string commands, the latex and MathML printers, the session and settings commands and area 18's plain Xcas names that this CAS lacked (the 72 aide entries of ## Context) run under their giac names, each with an oracle test against giac, a unit test where giac's library cannot show the behaviour, or a recorded departure, with its doc and a fuzz entry; CAS-DOC-58 records the area; CI and Fuzz pass on the pushed master head and the phase's Lific tickets are done (`gh run list --branch master --limit 2`).
Confirmed by the developer at the close.

## Tasks
- T1 — design page CAS-DOC-58 and one ticket per family | Done: Lific page CAS-DOC-58 "Strings, printers and session commands" (draft, like CAS-DOC-54..57 and CAS-DOC-7) and tickets CAS-611..623, one per ## Context family, all `todo`, in Language or Infrastructure. Summary: The phase's design page now exists in the tracker. It covers the scope (72 giac help entries, each marked as command, setting, session command, option builtin, alias or departure, plus what stays out and why), what the phase builds on, one section per build task with its design, giac sources, examples and departures, and the oracle plan. Thirteen tickets, one per family, give every later commit a ticket number. | commit: - | tasks/T1.md | logs/T1.md
- T2 — string commands: access, search, splitting, joining, trimming and distances | Done: char chr, asc, ord, inString, contains, mid, left lhs, right rhs, split, join, trim, strip, levenshtein, hamming are cas-core natives (native/strings.rs); size, concat, s[k], comparisons and sort take strings; CI green on a9598a5; CAS-611..613 done. Summary: The engine now has giac's string commands: building strings from character codes and back, finding and cutting substrings, splitting, joining, trimming, and edit distances. Strings also work in size, concat, indexing, comparisons and sort. A Python port of giac's index rules was checked against 127 giac answers before the Rust code was written. Both CI runs passed on the first push. | commit: a9598a5 | tasks/T2.md | logs/T2.md
- T3 — the printer and reader hook: string, str, expr, execute and cat | Done: cas-core's Context holds an optional `Language` hook (print a form, read one expression) that cas-lang implements and the session, oracle and fuzz eval target install; string str, expr execute, cat are natives; join and concat print non-strings; CI green on 4902032; CAS-614 done. Summary: Commands that turn expressions into text and back now work inside evaluation, so size(string(x+1)) or expr(string(e)) nest like any other call. The core engine cannot depend on the language crate, so it defines a small interface that the language crate fills in and installs. Installing it counts as a change of state, so nothing computed without it, or under another one, is ever reused. CI passed on the second push after one oracle row was fixed. | commit: 4902032 | tasks/T3.md | logs/T3.md
- T4 — the LaTeX printer: latex and TeX on a port of giac's markup engine | Done: latex(e) and TeX answer giac's export_latex text from a cas-lang port of gen2markup reached through the T3 hook (Form::Latex); CI green on 3c3a984 (run 45526771154); CAS-615 done. Summary: The engine can now write any expression as LaTeX, the way Xcas does. The printer is a port of giac's own markup code. It was first proved in a Python prototype that matched giac on 301 sample expressions, then written in Rust and fed by a small step that turns our expressions into the shapes giac's printer expects. Where giac and we evaluate the input differently (term order, rewrites giac makes before printing), the oracle tests record both answers as departures. | commit: 3c3a984 | tasks/T4.md | logs/T4.md
- T5 — the MathML printers: export_mathml and mathml | Done: export_mathml (default semantics, display, content) and mathml (giac's XHTML document) answer giac's text from ports of markup.cc and mathml.cc through the T3 hook; CI green on 6d1d5cf (run 45535388651); CAS-616 done. Summary: The engine can now write expressions as MathML in both of Xcas's styles: the modern one (presentation, content, or both linked by ids) and the older XHTML page. Both are ports of giac's own printer code, each first proved in a Python prototype against giac's answers on about three hundred sample expressions per form. Where giac and we evaluate an input differently, the oracle tests record both answers as departures. CI found two small bugs (a size cap counting text a run never writes, and one internal vector kind), both fixed. | commit: 6d1d5cf | tasks/T5.md | logs/T5.md
- T6 — settings: angle, complex and approximate modes, increasing_power, epsilon, cas_setup and giac's tuning words | Done: Xcas's 13 settings are session settings in the context state, read by name, set by `name := v`, `name(v)`, `cas_setup(list)`; angle, complex_mode, complex_variables, epsilon, approx_mode, autosimplify and increasing_power change answers or printing; CAS-617, CAS-618 done. Summary: The engine now has Xcas's configuration switches. Angles can be degrees or grads for the trigonometric functions, complex mode makes solve and factor work over the complex numbers, symbols can be declared complex, answers can be shown as decimals or passed through a simplifier, and sums of powers can be printed in increasing or decreasing order. The tuning words of giac's Gröbner engine and its evaluation depth are accepted and stored but change nothing, because ours evaluates fully and computes exactly. Where our defaults already differed from Xcas (term order, the default simplifier) we kept ours and recorded the difference in the oracle tests. | commit: 325b4df | tasks/T6.md | logs/T6.md
- T7 — session commands: restart, the variable lists, history, clocks, version and nodisp | Done: restart, rm_all_vars, rm_a_z, VARS, time, monotonic are session arms; ans, quest entry, version, nodisp natives over a bounded history in the Context; CI green on 325c274; CAS-619, CAS-620 done. Summary: The session now keeps a history of each statement's input and answer, which ans, quest and entry read, and has giac's commands for listing and removing the user's variables and for restarting with every setting back to its default. It also times an evaluation, reports the seconds since the session began, names its version and answers "Done" for nodisp. The history is read fresh every time, never from a cache, and it is bounded so long inputs cannot grow it without limit. | commit: 0e2c937 | tasks/T7.md | logs/T7.md
- T8 — area 18's names: type and the DOM_* constants, iPart, fPart, trunc, the TI sorts, sortperm, RandSeed, symbols and TI matrix names | Done: type with nine DOM_* answers; iPart, fPart frac, trunc; SortA sorta, SortD sortd, sortperm; RandSeed, SVL aliases; SVD, make_symbol, symbol_array, LSQ lsq, LQ builtins; QR/LU departures; CI green on 4d54483 (run 45547995014); CAS-621..623 done. Summary: The engine now answers giac's type queries, the TI calculator's integer and fractional parts and truncation (exact on exact input), the TI sorts and sortperm, and the TI matrix names for singular values, least squares and the LQ decomposition, plus commands that build symbols from strings. Each is checked against giac on CI, and where ours differs on purpose the difference is recorded with giac's own answer. The first push failed on two registry rules the build cannot see; one fix coder repaired them. | commit: 4d54483 | tasks/T8.md | logs/T8.md
- T9 — push master for the phase's final CI and Fuzz runs | Done: master b1c554d pushed and level with origin; CI run 45549018011 and Fuzz run 45549017969 started on it. Summary: The local gate (format, stable clippy, fuzz crate check) was green on the committed head, so master was pushed. GitHub started the phase's final CI and Fuzz runs on that commit; T10 reads their results. | commit: b1c554d | tasks/T9.md | logs/T9.md
- T10 — CI and Fuzz check on the pushed head; tickets done and CAS-DOC-58 in step | Done: CI 45549018011 and Fuzz 45549017969 green on pushed head b1c554d (every job); CAS-611..623 all `done`; CAS-DOC-58 §8's departures list matches crates/cas-giac/README.md ## Tests. Summary: The phase's final pushed commit passed every CI job and the full 43-minute fuzz run across all 56 targets, so nothing needed fixing. All thirteen phase tickets were already closed. The design page's list of places where we deliberately differ from giac had drifted from what the tests actually record, so it was rewritten from the test README; one README line that disagreed with its own test was corrected too. | commit: - | tasks/T10.md | logs/T10.md
- T11 — phase cleanup | Done: the printer prototypes and probe inputs live in crates/cas-lang/proto/ (rerun reproduces all three tables); the gap helpers are scripts/gap_describe.py; no repo file cites .run/1.6. Summary: The phase's scratch folder held the Python ports of giac's printers that the LaTeX and MathML work was checked against, and the code and tests cite them as where their expected rows came from. They now live in the repository next to the crate, with their probe inputs and a short guide to rerunning them, and a rerun gives the same tables byte for byte. The planner's small helpers for reading the gap lists became one script. Everything else in the scratch folder was one-off work or already lives elsewhere, so it is left for the close to delete. | commit: - | tasks/T11.md | logs/T11.md

## Decisions that still bind
- T2 — positions and lengths count Unicode scalar values (strings handled as Vec<char>); giac's byte answers are departure rows.
- T3 — the hook is part of the context state: set only through `set_language` (state_mut → new StateId); every T4-T7 form or reader goes through `Language` + a new `Form` variant, never a Context field outside State (C0066).
- T3 — natives print with the `terms` limit as the byte cap (print_capped); a print or read failure leaves the call with a `Warning::Native` carrying our text; no hook → Ok(None) only where a print/read is actually needed (cat of strings, expr("") still answer).
- T4 — T5 extends the same engine (markup/engine.rs MarkupBlock + Flags) with MathML fields and flags; the lowering (markup/lower.rs) is shared; a Form::MathML variant goes through hook.rs like Form::Latex.
- T5 — MathML presentation and content share T4's engine pass (one MarkupBlock, giac's design); mathml.cc's older layout is a separate cas-lang module over the same lowering.
- T6 — every setting lives in Context State (C0066); a setting's call form `name()`/`name(v)` is a session arm, never a native, when its bare name has a value (C0070).
- T7 — a native reading data kept outside State calls `Eval::volatile()` before reading (C0072).
-

## Rules proposed
- (none)

## Cookbook entries added
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

## Research
R1.6-001 | Map gen2markup (markup.cc:1233 to ~3600) and helpers for one Rust port producing LaTeX and both MathML forms; check against probe4.out. | gen2markup port map for latex and export_mathml (giac 2.0.0.21 markup.cc) | giac, markup, latex, mathml, port, phase-1.6, printers, priority, prepend_minus | retriever-code | 2026-11-24 | 118 lines
R1.6-002 | phase 1.6 settings and session-command tasks (match giac behaviour; departures where libgiac caseval cannot show it) | giac 2.0.0.21 settings and session commands: port reference | giac, settings, session, context, cas_setup, phase-1.6, port | retriever-code | 2026-11-24 | 89 lines
R1.6-003 | Phase 1.6 strings and Xcas aliases port; per-function file:line, argument forms, returns, index base, errors, quirks | giac 2.0.0.21 string commands and Xcas aliases: port semantics | giac, strings, aliases, port, phase-1.6 | retriever-code | 2026-11-24 | 65 lines
R1.6-004 | phase 1.6 latex task: map /home/user/code/cas/.run/1.6/src/tex.cc gen2tex and helpers for a faithful Rust port | giac 2.0.0.21 LaTeX printer map for Rust port (tex.cc) | giac, latex, tex.cc, printer, phase-1.6, port-size, probe4 | retriever-code | 2026-11-24 | 156 lines
R1.6-005 | phase 1.6 plan, mathml task (port size, rules a Python prototype must reproduce, done-when rows) | giac 2.0.0.21 MathML printers: port map (mathml.cc and markup.cc) | giac, mathml, export_mathml, gen2mathml, phase-1.6, port-map | retriever-code | 2026-11-24 | 106 lines
R1.6-006 | phase 1.6 plan ## Interfaces; verify (check) lines; answer Q1-Q3 | Phase 1.6 plan: Interfaces (session commands, settings, string/printer hook) | phase-1.6, interfaces, session, settings, print-hook, cas-giac | retriever-code | 2026-11-24 | 93 lines
R1.6-007 | phase 1.6 plan, ## Interfaces and task split (which giac string commands can be cas-core natives, which need cas-lang) | Strings and printing today, for giac's string commands and latex/mathml printers | cas-core, cas-lang, strings, print, session, interfaces, phase-1.6 | retriever-code | 2026-11-24 | 121 lines
R1.6-008 | extract for phase 1.6 plan ## Context the exact name lists per family, R0-004 notes, R0-002 findings, CAS-DOC-7 binding rules, and gap cross-check | Phase 1.6 context: strings, printers, session names (R0-002, R0-004, CAS-DOC-7, missing.txt) | phase-1.6, strings, printers, latex, TeX, mathml, export_mathml, session, assume, purge, call, output, DOM_, aliases, gap | retriever-digest | 2026-11-24 | 68 lines
R1.6-009 | extract Rationale and T2..T10 design decisions, done-when essentials, giac rows and departures for the design page (T1) | Phase 1.6 PHASE_PLAN: Rationale and T2-T10 design facts for CAS-DOC-58 | phase-1.6, design-page, CAS-DOC-58, giac, strings, printers, markup, mathml, settings, session, oracle | retriever-digest | 2026-11-25 | 182 lines
R1.6-010 | T4 coder brief, Expr -> giac-shaped LaTeX lowering (crates/cas-lang/src/markup) | cas-lang input printer: canonical evaluated forms and where each is printed | cas-lang, print, names, lower, forms, giac | retriever-code | 2026-11-25 | 55 lines

## Audit
- seed: median 24k, max 25k, n=11 (phase 1.6)
- previous phase 1.5: median 24k, growth -1.2%
- CLAUDE.md: 1109 bytes (-81)
- .claude/skills/pa-external-agents/SKILL.md: 11144 bytes (unchanged)
- .claude/skills/pa-planning/SKILL.md: 8419 bytes (unchanged)
- .claude/skills/project-architect/SKILL.md: 13748 bytes (unchanged)
- .claude-state/memory/MEMORY.md: 53 bytes (unchanged)
- HOW_WE_WORK.md: 8117 bytes
- cookbook/INDEX.md: 13249 bytes
- rules/INDEX.md: 3242 bytes
### Carry audit — phase 1.6 (2026-11-24T23:24:56Z → open UTC, 1 sessions, 1687 requests)

| file (read) | chars | n |
|---|---|---|
| markup.cc | 334.2k | 37 |
| print.rs | 196.1k | 34 |
| PHASE_PLAN.md | 187.8k | 13 |
| engine.rs | 160.7k | 17 |
| markup.py | 104.8k | 7 |
| context.rs | 92.5k | 32 |
| tests.rs | 90.8k | 19 |
| b09t3ptsa.txt | 66.0k | 4 |
| session.rs | 62.9k | 14 |
| mathml.cc | 62.0k | 4 |
| file (write) | chars | n |
|---|---|---|
| engine.rs | 111.4k | 31 |
| retriever-phase-plan-rationale-t2-t10.md | 72.3k | 2 |
| markup.py | 60.4k | 1 |
| page.md | 45.6k | 1 |
| mathml_proto.py | 39.3k | 4 |
| markup.py | 34.7k | 1 |
| convert.rs | 28.8k | 1 |
| tables.rs | 27.7k | 1 |
| result kind | chars | n |
|---|---|---|
| Read | 4.5M | 896 |
| bash other | 2.6M | 1248 |
| Grep | 1.2M | 547 |
| tools/plan_edit.py | 55.3k | 20 |
| run.sh | 54.0k | 93 |
| Glob | 48.7k | 39 |
| tools/outline.py | 39.4k | 20 |
| Agent | 34.4k | 34 |
| Write | 23.6k | 146 |
| tools/commit_task.sh | 19.7k | 50 |
- spilled: 13 results, 461.3k chars on disk, 1 read whole
- whole reads over 20.0k: 3
- seed floor: retriever-code 6,193 (n 51, prev 6,193) · retriever-digest 6,236 (n 8, prev 6,236) · retriever-web 4,600 (n 2, prev 4,600)
- outline credit: 24 outlines · 15 followed by a ranged read · 491.4k chars credited
- warm pings: 6 pings over 1 runs, 2 warmed waits over the TTL, rewrites across warmed waits 0, waits past the cap 0, cost: pings $0.0954 + relay $0.0000 + re-arm $0.0000 = $0.0954 vs rewrites replaced $1.1786, cheaper than one rewrite: yes, undelivered 0
- toasts by cause (waiting): question 0, review 0, replan 0, permission 0, input 0, crash 0, idle 0, stop 0, subagent-stop 0, model 0, other 0
- router turns by cause: loop 821, relay 0, re-arm 13, other 34, relay cost $0.0000
- retriever re-asks: 2 of 10 retriever briefs repeat a lookup of the same run
- noise: 23 lines 1.4k chars — No such file or directory: 15 lines/1.1k, usage: : 8 lines/347
- price: read $0.20/Mtok · 1h write $7.94/Mtok · output $19.84/Mtok (phase model mix)
- median requests after a read: 0
- carry/request: 5083 chars, 1687 requests (prev 13448 chars, 798 requests, growth -62.2%)
- flag: 1 spilled results read whole
  command: P=/usr/bin/python3; $P tools/plan_edit.py show --section Context; $P tools/plan_  reader: pa-expert
- flag: 2 whole-plan reads by pa-expert
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
- T11 | pa-claude | partner | claude-opus-5-5 | medium | ctx 85k | $1.067 | saved - | completed | parent 1b3aa32d
- T1 | pa-expert | expert | claude-opus-5-5 | medium | ctx 145k | $2.392 | saved $0.64 | completed | parent 1b3aa32d
- T1 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 57k | $0.023 | saved - | completed | parent adbcafb1
- T1 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 44k | $0.023 | saved - | completed | parent adbcafb1
- T2 | pa-expert | expert | claude-opus-5-5 | medium | ctx 132k | $2.741 | saved $2.99 | completed | parent 1b3aa32d
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 109k | $2.025 | saved - | completed | parent a4e4ea44
- T2 | pa-coder | coder | claude-opus-5-5 | medium | ctx 115k | $2.124 | saved - | completed | parent a4e4ea44
- T3 | pa-expert | expert | claude-opus-5-5 | medium | ctx 127k | $2.727 | saved $9.72 | completed | parent 1b3aa32d
- T3 | retriever-digest | retriever | claude-haiku-5-5 | medium | ctx 18k | $0.005 | saved - | completed | parent a0890283
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 77k | $1.112 | saved - | completed | parent a0890283
- T3 | pa-coder | coder | claude-opus-5-5 | medium | ctx 112k | $2.152 | saved - | completed | parent a0890283
- T4 | pa-expert | expert | claude-opus-5-5 | medium | ctx 127k | $2.795 | saved $17.64 | completed | parent 1b3aa32d
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 208k | $4.012 | saved $0.21 | completed | parent a2e536b2
- T4 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 45k | $0.015 | saved - | completed | parent a2e536b2
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 278k | $6.813 | saved - | completed | parent a2e536b2
- T4 | pa-coder | coder | claude-opus-5-5 | medium | ctx 98k | $1.820 | saved - | completed | parent a2e536b2
- T5 | pa-expert | expert | claude-opus-5-5 | medium | ctx 120k | $2.598 | saved $35.92 | completed | parent 1b3aa32d
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 306k | $6.456 | saved $0.85 | completed | parent ae285a00
- T5 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 14k | $0.005 | saved - | completed | parent a7e33a8a
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 196k | $4.628 | saved $0.77 | completed | parent ae285a00
- T5 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 13k | $0.004 | saved - | completed | parent a07aa828
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 313k | $8.107 | saved $0.27 | handoff | parent ae285a00
- T5 | pa-coder | coder | claude-opus-5-5 | medium | ctx 107k | $1.862 | saved - | completed | parent ae285a00
- T6 | pa-expert | expert | claude-opus-5-5 | medium | ctx 146k | $3.805 | saved $6.93 | completed | parent 1b3aa32d
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 133k | $3.022 | saved - | completed | parent a0016725
- T6 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 25k | $0.007 | saved - | completed | parent a9b04735
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 126k | $2.595 | saved - | completed | parent a0016725
- T6 | pa-coder | coder | claude-opus-5-5 | medium | ctx 123k | $2.663 | saved $0.10 | completed | parent a0016725
- T6 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 28k | $0.009 | saved - | completed | parent a6abe4d9
- T6 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 29k | $0.008 | saved - | completed | parent a6abe4d9
- T7 | pa-expert | expert | claude-opus-5-5 | medium | ctx 104k | $1.757 | saved $0.67 | completed | parent 1b3aa32d
- T7 | pa-coder | coder | claude-opus-5-5 | medium | ctx 140k | $2.761 | saved - | completed | parent a258ceb7
- T8 | pa-expert | expert | claude-opus-5-5 | medium | ctx 121k | $2.447 | saved $7.79 | completed | parent 1b3aa32d
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 134k | $2.788 | saved - | completed | parent a0964be1
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 77k | $1.137 | saved - | completed | parent a0964be1
- T8 | retriever-code | retriever | claude-haiku-5-5 | medium | ctx 17k | $0.005 | saved - | completed | parent a3f6df67
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 122k | $2.295 | saved - | completed | parent a0964be1
- T8 | pa-coder | coder | claude-opus-5-5 | medium | ctx 71k | $1.008 | saved $0.05 | completed | parent a0964be1
- T9 | pa-expert | expert | claude-opus-5-5 | medium | ctx 34k | $0.289 | saved $0.13 | completed | parent 1b3aa32d
- … 2 more: PY ~/.claude/pa3/pa_ledger.py report

## Discussions
- (none)

## Deferred
- from T11: the close may delete .run/1.6/ whole; .run/logs/gen1-gap/ stays until 1.8 (R45). (unconsumed)

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

## Plain-English Recap
Phase 1.6 Strings, printers and session commands: CAS-DOC-58 (still a draft in Lific) sorts 72 giac entries; tickets CAS-611..623 are done. Built:
- Strings: 14 string natives (char, asc, ord, inString, contains, mid, left/right, split, join as StringJoin, trim, strip, levenshtein, hamming). size, concat, indexing, comparisons and sort handle strings; positions count Unicode characters, so giac's byte answers are departures.
- A cas-core Language hook trait, implemented in cas-lang hook.rs and set via state_mut, so it gives a new state id (C0066). It backs string/str, expr/execute and cat.
- latex/TeX: a port of giac's gen2markup (markup.rs), 117 agreeing rows.
- export_mathml in giac's three forms, and mathml in mathml.cc's XHTML layout. The 17 printer departures come from giac's evaluation, not the printer.
- 13 Xcas settings in State. increasing_power defaults to 1 (giac 0) to keep existing output; autosimplify defaults to none.
- Session commands from a bounded history held outside State (restart, ans/quest/entry via Eval::volatile, VARS, rm_all_vars, time, version, nodisp).
- Area 18 names: type with DOM_*, iPart, fPart, trunc, the TI sorts, sortperm, RandSeed, make_symbol, symbol_array, SVD, LSQ, LQ. QR and LU are departures.
Each Python prototype was checked against giac first; the printer prototypes are now kept in crates/cas-lang/proto/, and the gap helpers became scripts/gap_describe.py. Final gate green on b1c554d with no fixes needed (CI 45549018011, Fuzz 45549017969, 56 targets). Cookbook C0065-C0076. Open: a ticket for diff(e,x,n), which doesn't take the nth derivative; watch symbol_array's memory (up to 1e6 symbols, never freed) on the next Fuzz run; publish CAS-DOC-57 and CAS-DOC-58; docs-only commits unpushed until the next phase's first push.
