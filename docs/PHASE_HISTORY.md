# Phase history

Built by `tools/phaseend_index.py history` from the PhaseEnds and rebuilt at each close; never edited by hand. Each phase's PhaseEnd holds the rest.

## Phases

- 1.1 | 2026-11-22 | Green master | Phase 1.1 Green master: master's CI and Fuzz are green again on fa9ae22 (CI run 45047911216, all five jobs; Fuzz run 45047911210, all 55 targets).
- 1.2 | 2026-11-23 | Algebra and polynomials parity | Phase 1.2 Algebra and polynomials parity: CAS-DOC-54 sorts 89 giac entries (66 built, 1 alias, 5 options, 2 settings, 15 recorded departures); tickets CAS-526..547 done.

## This generation's recaps (1, newest first)

### 1.2: Algebra and polynomials parity (2026-11-23)

Phase 1.2 Algebra and polynomials parity: CAS-DOC-54 sorts 89 giac entries (66 built, 1 alias, 5 options, 2 settings, 15 recorded departures); tickets CAS-526..547 done. Built: named gcds (heugcd, modgcd, psrgcd, new EZ-GCD) and icontent; real root isolation (VAS), LMQ bounds, Sturm sequences, rational/complex root lists (ported from giac source, checked by Python prototype); cfsolve with certified complex roots, zeros, cZeros, surd, all_trig_solutions and with_sqrt settings; cpartfrac, propfrac, f2nd, convert options; rootof, fcoeff, homogeneize, total_degree, factor_xn, bezoutian, is_irreducible; in_ideal, syzygy, gbasis_vars, gbasis rur; 28 matrix names plus det_minor, adjoint_matrix, ref, pivot, keep_pivot, egvl, gauss, rat_jordan. Final gate green on 5680a2c (CI 45091838361, Fuzz 45091838357) after a fuzz timeout on (1+x^2/4)^46656 led to a combined size limit on polynomial powers (CAS-547). The flaky nlpsolve oracle row is fixed by seeding giac. Cookbook C0007-C0019 added. Open: the developer's ruling on local throwaway example programs; a cas-num fix for Abs of a complex ball containing 0; six log/doc commits unpushed until the next phase's first push.

### 1.1: Green master (2026-11-22)

Phase 1.1 Green master: master's CI and Fuzz are green again on fa9ae22 (CI run 45047911216, all five jobs; Fuzz run 45047911210, all 55 targets). The eval fuzz crash came from FromContinuedFraction (CAS-424), not CAS-517..524: a zero last denominator (-2 + 1/0) reached Rational::new; it now gives ComplexInfinity, under a regression test in cas-core (CAS-525). boundary_values now checks BVPSolve against closed forms on small grids and giac only at step 1/20; giac's n=4 error (0.044) is a recorded departure, and the accuracy comments were corrected (CAS-520). giac_gap.sh moved to scripts/ with an output-folder argument; the card lists Numerical analysis and Signal processing. CAS-517..525 are done in Lific. Decided during the phase: agents push master themselves when CI needs to run (card 705b703, after the developer removed the never-push fail-safe from CLAUDE.md). Open: .claude/settings.json has a trailing comma from the developer's edit (his to fix); the generation plan's standing constraints still describe a push stop, so each later phase brief says agents push; lific issue list caps at 500, so ticket checks use issue get; checkout@v4 Node 20 deprecation warning.
