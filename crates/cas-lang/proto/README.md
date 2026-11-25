# Printer prototypes (Python ports of giac's printers)

Our answers are seen only on CI (R74), so expected printer shapes come from these ports, each checked
row for row against giac 2.0.0.21 through libgiac. Not compiled, not shipped. Cited by `src/markup*`,
`src/mathml*` and `crates/cas-giac/tests/mathml.rs`; departure rows: C0068.

- `latex/` (1.6/T4.c1, CAS-615): `proto/markup.py` ports markup.cc gen2markup's LaTeX arms; `parse.py`
  reads giac's printed form into giac's internal shapes; `NOTES.md` the port's notes; `probe.cas` 301 inputs.
- `mathml/` (1.6/T5.c1, CAS-616): `proto/markup.py` adds the MathML presentation and content arms;
  `probe.cas` export_mathml of latex's inputs in three forms (893 rows).
- `mathml/doc/` (1.6/T5.c3, CAS-616): `mathml_proto.py` ports mathml.cc (`mathml()`); `probe.cas` plus
  `extra.cas` 320 rows; `printsommet.txt` giac's printsommet/texprint registrations, from
  `printsommet.py <giac src dir>` (source: C0009).

Rerun (repo root): `/usr/bin/python3 scripts/caseval.py <dir>/probe.cas > <dir>/probe.out` (and `extra.cas`
to `extra.out` in `mathml/doc/`); latex's `probe.out` is read by all three. Then `<dir>/.../check.py`
prints `n/n ok` and writes `table.tsv` beside the probe (generated; not committed). The committed
tables in `src/markup/giac_*.tsv` were built from those tables.
