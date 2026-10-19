# CAS

A computer algebra system engine written in Rust. The programming API is first class;
user-facing syntaxes are frontends built on top of it.

## Crates

| Crate | Purpose |
|---|---|
| `cas-num` | Arbitrary-precision integers, rationals, modular arithmetic, ball arithmetic |
| `cas-core` | Hash-consed symbolic expressions, symbols and attributes |
| `cas-lang` | The surface language: lexer, parser and name table (lowering to come) |

## Development

```sh
cargo build
cargo test
cargo clippy --all-targets
```

Work is tracked in Lific (project `CAS`); see `AGENTS.md`. Commits reference their
ticket: subject `CAS-N: …`, footer `Closes: CAS-N`.

## License

GPL-3.0-only. See `LICENSE`.
