# TASK_PROGRESS — T5.c3 (CAS-616, giac `mathml`), handoff after commit 3fbf832

## Done so far (committed)
- a152922 CI fix (expert's message): markup/engine.rs `MarkupBlock::len(flags)` counts only texts the run writes (tex/pres/cont); all `check` sites pass flags; matrix arm counts LaTeX rows only if `flags.tex()`. Fixes markup/tests.rs:972 `text_past_the_cap_is_refused` (MathML run counted the list's LaTeX, made then cleared).
- 3fbf832 module `crates/cas-lang/src/mathml.rs` + `mathml/{convert,flat,print}.rs`; `pub use mathml::mathml` in lib.rs; markup `pub(crate) mod tree/ident`, `pub(crate) fn tree`, tree.rs items pub(crate), `my_isalpha` pub(crate), markup tests `pub(crate) mod tests` + `pub(crate) fn parse` (for the mathml tests to reuse). Build, fmt, stable clippy (C0003) and fuzz check green. Not wired, no tests.
- Prototype `.run/1.6/T5/c3/mathml_proto.py`: 320/320 (301 probe rows `probe.out` = giac `mathml(X)` on .run/1.6/T5/inputs.txt, + 19 `extra.out`), `python3 check.py`; `single.py` (aplatir once at top, as Rust) also 320/320.

## In flight
- Nothing uncommitted (no patch).

## Hypotheses rejected (evidence)
- Rows sqrt(x): giac evaluated tree is pow(x,1/2) (`<msqrt><mrow>`), not sqrt symbolic (`<msqrt>` no mrow) -> `giac_shape` sqrt->pow in mathml.rs; x^(-1/2) = pow(pow(x,1/2),-1) matches our lowered `sqrt(x)^-1`.
- Quoted products are nested binary in giac's parser (quote(2*inv(4)*x) -> 1/2 via sym_mult): check.py `nest()`.
- `pow` exponent ==1/2 compares numbers across types (quote(x^0.5) -> msqrt).
- _FUNC prints quoted (`<mi>'sin'</mi>`); poly1 vect prints `poly1[`; `pointminus` name `.-`.

## Current hypothesis
- Rust port = prototype (single-pass aplatir). Unverified by running (R74); CI judges.

## Next 5 steps
1. Tests `crates/cas-lang/src/mathml/tests.rs` (+ `#[cfg(test)] mod tests;` in mathml.rs): tsv `crates/cas-lang/src/markup/giac_mathml_doc.tsv` from `.run/1.6/T5/c3/table.tsv` 301 probe rows only (columns input, mode, shape, giac text with ' | '); shape via `crate::markup::tests::parse(shape, quote)`, quote rows nested (port check.py `nest`), 1e10/12345678901.0 Dbl overrides as markup tests; compare `document(&g, usize::MAX)` with '\n'->' | '. Plus TooDeep/TooLong (cap = len passes, len-1 fails), inline `[x,-1]`.
2. cas-core: `Form::MathmlDocument` (language.rs), hook.rs arm -> `crate::mathml::mathml`; native `MathML` in native/text.rs (one arg -> printed(ev,e,Form::MathmlDocument); 2+ args Ok(None)); builtin.rs `MATH_ML = "MathML"` near EXPORT_MATH_ML; native.rs registration; en.txt `MathML: mathml`; fuzz/fuzz_targets/eval.rs BUILTINS 610->611 (+"MathML" near "ExportMathML" :703); native/text/tests.rs stub Forms arm + test.
3. hook rows in markup tests pattern (SAME_SHAPE) for Form::MathmlDocument vs tsv.
4. crates/cas-giac/tests/mathml.rs: append `mathml_documents` (latex.rs inputs minus 17 departures: giac `mathml(x)` unquoted == ours exact) and `mathml_departures` (17 rows, bodies from `.run/1.6/T5/c3/deps.out`, ours = PREAMBLE+wrap(body)).
5. Docs: cas-lang README row (mathml), cas-core README native::text (MathML), cas-giac README mathml paragraph; log; gate; commit.

## Gotchas
- `gen` is a reserved keyword (edition 2024): method is `Writer::node`.
- commit_task.sh stages logs/ too.
- Prototype bug to fix: `printassto` test should be `gtype(front)==1` (_SEQ__VECT=1=_DOUBLE_), proto has a wrong expression; Rust uses DOUBLE (correct). No row covers it.
- Hook PreToolUse denies shell reads of files >20k chars (mathml_proto.py 37k): use Read with offset or python scripts.

## State to carry verbatim
- giac wrapper: PREAMBLE (mathml.cc:69) + "\n<math mode=\"display\" xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\n" + body + "\n\n</math><br/>\n" + "\n</body>\n</html>" + "\n". mathml(1/2) = 502 chars.
- Departures (ours body, prototype on our lowered shapes): `.run/1.6/T5/c3/deps.out` (17 rows, tab-separated).
