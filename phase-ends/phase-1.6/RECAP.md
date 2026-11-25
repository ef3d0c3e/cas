## Recap
Phase 1.6 Strings, printers and session commands: CAS-DOC-58 (still a draft in Lific) sorts 72 giac entries; tickets CAS-611..623 are done. Built:
- Strings: 14 string natives (char, asc, ord, inString, contains, mid, left/right, split, join as StringJoin, trim, strip, levenshtein, hamming). size, concat, indexing, comparisons and sort handle strings; positions count Unicode characters, so giac's byte answers are departures.
- A cas-core Language hook trait, implemented in cas-lang hook.rs and set via state_mut, so it gives a new state id (C0066). It backs string/str, expr/execute and cat.
- latex/TeX: a port of giac's gen2markup (markup.rs), 117 agreeing rows.
- export_mathml in giac's three forms, and mathml in mathml.cc's XHTML layout. The 17 printer departures come from giac's evaluation, not the printer.
- 13 Xcas settings in State. increasing_power defaults to 1 (giac 0) to keep existing output; autosimplify defaults to none.
- Session commands from a bounded history held outside State (restart, ans/quest/entry via Eval::volatile, VARS, rm_all_vars, time, version, nodisp).
- Area 18 names: type with DOM_*, iPart, fPart, trunc, the TI sorts, sortperm, RandSeed, make_symbol, symbol_array, SVD, LSQ, LQ. QR and LU are departures.
Each Python prototype was checked against giac first; the printer prototypes are now kept in crates/cas-lang/proto/, and the gap helpers became scripts/gap_describe.py. Final gate green on b1c554d with no fixes needed (CI 45549018011, Fuzz 45549017969, 56 targets). Cookbook C0065-C0076. Open: a ticket for diff(e,x,n), which doesn't take the nth derivative; watch symbol_array's memory (up to 1e6 symbols, never freed) on the next Fuzz run; publish CAS-DOC-57 and CAS-DOC-58; docs-only commits unpushed until the next phase's first push.
