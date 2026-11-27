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
- Gen 2 carry-overs from 1.8 -> Gen 2: (a) CAS-650 backlog: eight functions re-evaluate an argument under fresh limits like the old markov eps (random.rs:147, ndsolve.rs:121, ...), each fix with its own test; (b) markov and lpsolve suit random oracle sweeps but no sweep group covers them; (c) a residue rule for rational functions would turn 5 Hilbert departures into agreeing rows; (d) Fill inside begin…end gives purge's misplaced-form error, only if gen 2 wants TI-style programs; euler_lagrange_assume.rs:45's ,t,2) rewrite is no longer needed.
