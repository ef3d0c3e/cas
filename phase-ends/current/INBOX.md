# INBOX

<!-- The developer writes here whenever they want; nobody asks them to. Freeform bullets.
     The loop consumes it at every task boundary (the mod), with `tools/status.py inbox-consume`,
     which prints the contents, archives them to discussions/inbox-<ts>.md and empties
     this file. Items that change the plan become plan_edit.py edits with a Changes line;
     items that change a standing fact go to HOW_WE_WORK.md.
     Routing: each bullet becomes an I<n> row in discussions/INDEX.md, deferred:<P> or now.
     First match wins: `-> 3.12` / `-> Gen 4`; `Gen N` beyond this generation; a later
     plan phase id (`3.11 ...`); the word `later` (next open phase); else now. The current
     or an earlier phase, or an id not in GENERATION_PLAN.md, routes now. -->

- <note, correction, priority change, or question>
- -> Gen 2: frontends, decided 2026-11-25. First a language server over the CAS: completion, highlighting (an editor grammar file plus meaning-based colouring from the server), documentation from en.txt and the doc chapters, and error reporting. Evaluation, the defined-symbols view and similar features go through custom requests. It needs an error-tolerant parser in cas-lang. Then a desktop app as the language server's client: a view of the defined symbols, and graphs, which need a plotting layer (giac's ~121 plot commands were out of gen 1). For the developer at gen 2 planning: the new Lific module(s), and the app's toolkit (a retriever compares the options first). For the planner: how the app is tested on CI only.
