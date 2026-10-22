# CAS

A computer algebra system engine written in Rust. The programming API is first class;
user-facing syntaxes are frontends built on top of it.

## Crates

| Crate | Purpose |
|---|---|
| `cas-num` | Arbitrary-precision integers, rationals, modular arithmetic, ball arithmetic |
| `cas-core` | Hash-consed symbolic expressions, symbols and attributes; printing and reading full form; the evaluation loop |
| `cas-lang` | The surface language: lexer, parser, name table, lowering to expressions and the input-form printer (the session to come) |

## Development

`cargo check` and `cargo build` are fine locally. Tests and fuzzing run on GitHub
Actions only — never on the development machine (CAS-158):

- `.github/workflows/ci.yml` — format, clippy, the test suite and a nightly check
  of the fuzz crate, on every push and pull request.
- `.github/workflows/fuzz.yml` — every libFuzzer target on every push, one job
  per target in parallel, ten minutes each; a manual dispatch
  (`gh workflow run fuzz.yml -f time_per_target=3600 -f targets="ball eval"`)
  overrides the budget and the target set. The corpus is carried between runs
  and crash artifacts are uploaded on a finding.

Work is tracked in Lific (project `CAS`); see `AGENTS.md`. Commits reference their
ticket: subject `CAS-N: …`, footer `Closes: CAS-N`.

## License

GPL-3.0-only. See `LICENSE`.
