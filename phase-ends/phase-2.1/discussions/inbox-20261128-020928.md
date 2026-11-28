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
- giac 2.1 reference move -> Gen 2 later phase or developer's call: debian:sid moved libgiac0t64 to 2.1.0.1+dfsg2-1 (CI 46069824120, C0122); 2.1 T-pin keeps CI on 2.0.0.21 to match local caseval. Moving means updating local giac, re-probing departures, and fixing mathml.rs:377 (<mo>( for <mfenced>) and derivatives.rs:227.
- -> 2.2: 2.1 closed with three loose ends: random_n_deriv (derivatives.rs:273) fails intermittently in its near check (derivatives.rs:37-51); pin nDeriv(-2x+5x^2,x,1)@3/10 and nDeriv((4x+x^3)/(x^2+4),x,1)@3/10 (seed cc a9ff39f6) as rows and fix the check. CAS-661 (sturmseq) still open in Lific: set done. Confirm Fuzz 46082504901 on 5990fb7 completed green (R117). Integrate's DiracDelta/piecewise paths (integrate.rs:108-125) have no CAS-650 guard -> 2.6.
