# Phase history

Built by `tools/phaseend_index.py history` from the PhaseEnds and rebuilt at each close; never edited by hand. Each phase's PhaseEnd holds the rest.

## Phases

- 1.1 | 2026-11-22 | Green master | Phase 1.1 Green master: master's CI and Fuzz are green again on fa9ae22 (CI run 45047911216, all five jobs; Fuzz run 45047911210, all 55 targets).
- 1.2 | 2026-11-23 | Algebra and polynomials parity | Phase 1.2 Algebra and polynomials parity: CAS-DOC-54 sorts 89 giac entries (66 built, 1 alias, 5 options, 2 settings, 15 recorded departures); tickets CAS-526..547 done.
- 1.3 | 2026-11-23 | Calculus, expressions and integers parity | Phase 1.3 Calculus, expressions and integers parity: CAS-DOC-55 sorts 81 giac entries; tickets CAS-548..585 are done.
- 1.4 | 2026-11-24 | Numerical analysis II | Phase 1.4 Numerical analysis II: CAS-DOC-56 sorts 41 giac names and records CAS-517..520; tickets CAS-588..599 are done.
- 1.5 | 2026-11-24 | Optimization II | Phase 1.5 Optimization II: CAS-DOC-57 (still a draft in Lific) sorts 50 giac names; tickets CAS-600..610 are done.
- 1.6 | 2026-11-25 | Strings, printers and session commands | Phase 1.6 Strings, printers and session commands: CAS-DOC-58 (still a draft in Lific) sorts 72 giac entries; tickets CAS-611..623 are done.
- 1.7 | 2026-11-26 | Statistics III and signal processing II | Phase 1.7 Statistics III and signal processing II: CAS-DOC-59 (still a draft in Lific) covers 37 giac entries; tickets CAS-624..635 are done, and CAS-636 (diff(e,x,n) doesn't take the nth derivative)…
- 1.8 | 2026-11-27 | Validation | Phase 1.8 Validation: final gate green on 615a65e (CI 45885347583, push Fuzz 45885347505, one-hour eval re-run 45891167338).
- 2.1 | 2026-11-28 | Engine carry-overs (I3) | Phase 2.1 Engine carry-overs: closed on the developer's yes.

## This generation's recaps (2, newest first)

### 2.1: Engine carry-overs (I3) (2026-11-28)

Phase 2.1 Engine carry-overs: closed on the developer's yes. Final gate on 5990fb7: CI 46082504855 green, Fuzz 46082504901 completed green. CAS-650 class removed: T3 classified all 59 approximation( sites in native/ (8 CAS-650 + 51); 30 sites got an is_numeric guard (T2 8, T3 17 symbolic, T10 5 endpoint/guess/constant), T11 guards 7 sample-point sites with holds_nested (refuse own-head nesting depth >= 2; depth 1 answers); every guard has a 30-level nest test under 2 s (C0118). Sweeps random_markov_chains and random_lp_objectives (128 cases; giac's wrong 3-variable maximize is a departure row, C0119). Hilbert residue rule (CAS-654): the five rational rows agree with giac (C0120). Fixes from sweeps: sturmseq sign placement (CAS-661, C0123), random_is_irreducible skip by content in y (CAS-660, C0121). euler_lagrange_assume.rs:45 rewrite gone (CAS-659). CI's oracle giac pinned to 2.0.0.21+dfsg2-3+b1 from snapshot.debian.org with SHA-1 and dpkg-query checks (C0122, C0124), after sid's 2.1.0.1 broke two rows. Fill inside begin…end is backlog CAS-656. audit_tool_source_prefix fixed in PA3's audit.py (pa3-src 96a421a). Rules M13, M14 added. Open, in the inbox for 2.2: random_n_deriv's near check fails intermittently (two inputs recorded, unpinned); CAS-661 still to set done in Lific; Integrate's DiracDelta/piecewise paths unguarded (-> 2.6); moving to giac 2.1 is the developer's call; CAS-DOC-57/58/59 drafts; commit_task.sh still sweeps unnamed files.
