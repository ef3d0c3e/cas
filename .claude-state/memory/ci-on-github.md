---
name: ci-on-github
description: All CAS tests/fuzzing run on GitHub Actions (public repo ef3d0c3e/cas); gh is authenticated for logs; mind CI memory limits
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-28T18:48:08.180Z
---

All tests and fuzzing run on GitHub CI, never locally. The repo ef3d0c3e/cas is public, so CI minutes are free; `gh` is authenticated, use `gh run view --job <id> --log` to read results.

**Why:** the user wants the dev machine kept free of test/fuzz load and said to be careful with CI's memory usage limit.

**How to apply:** push to a branch (or `gh workflow run fuzz.yml -f ...`) and read logs with gh; keep fuzz jobs under libFuzzer's `-rss_limit_mb=2048` and the runner's RAM (CAS-162). See [[lific-is-remote-only]].

Clippy: CI's clippy job runs the **stable** toolchain (`dtolnay/rust-toolchain@stable`), but the local default is nightly and its bin dir is on PATH directly (no rustup proxies), so `cargo +stable` and `rustup run stable cargo clippy` still run nightly's clippy, which misses lints stable flags (needless_range_loop broke CI on 2026-10-28). To match CI: `S=~/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin; PATH=$S:$PATH CARGO_TARGET_DIR=target/stable cargo clippy -q -p <crate> --lib --tests -- -D warnings`.

Watching runs: `gh run watch` polls the API every 3 s by default; two watchers over a ~25 min CI + fuzz cycle, several times an hour, exhausted the 5000/h API quota on 2026-10-25 (HTTP 403 "API rate limit exceeded"; `gh api rate_limit` shows the reset and does not count). Use `gh run watch --interval 60`, or a background `until ... | grep -q completed; do sleep 60; done` loop, and one watcher per run.
