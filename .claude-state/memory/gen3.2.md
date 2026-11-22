## benchmarks-on-ci

---
name: benchmarks-on-ci
description: Measure performance with the Workload workflow on CI (A/B against a base commit on one runner); the runner has no perf and is ~1.3x slower than the dev machine
metadata:
  node_type: memory
  type: reference
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-25T15:46:53.307Z
---

Workloads (ignored `*_workload` tests) run through `.github/workflows/workload.yml`, never locally: `gh workflow run workload.yml -f test=<path> -f runs='ENV=a;ENV=b' [-f base=<commit>]`. It reports best-of-3 wall time in ms and max RSS; `base` builds an older commit on the same runner for a fair A/B. Read the table from the job log (`gh api repos/ef3d0c3e/cas/actions/jobs/<id>/logs`, lines starting `| `).

- The GitHub runner allows no `perf` (instruction counts show `-`); compare times only.
- The runner takes about 4/3 of the dev machine's time (old local numbers: Fateman 9.8 s local vs 12.62 s runner), which scales older local or giac figures in the docs.
- Comparing methods within one commit (e.g. `CAS_POLY_GCD=sparse:modular:12` vs `sparse:zippel:12`) needs no `base`.

See [[ci-on-github]].

## ci-diagnostics-with-gdb

---
name: ci-diagnostics-with-gdb
description: How to find where a slow, hanging or out-of-memory CAS fuzz input goes wrong without running anything locally — gdb backtraces, input printing and threshold panics on a throwaway CI branch
metadata:
  node_type: memory
  type: reference
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-11-05T13:36:30.509Z
---

When a fuzz input times out or is slow (eval target, CAS-166/CAS-169), timing instrumentation alone can miss the spot (the time was in `Trace::verify`, which the harness did not time). What found it: a throwaway branch (e.g. `cas-166-diag`) with a temporary `.github/workflows/diag.yml` triggered on push to that branch only, which builds the target (`cargo fuzz build --target x86_64-unknown-linux-gnu -O -a -s none eval`, env `CARGO_PROFILE_RELEASE_DEBUG: line-tables-only`), runs the input with `-runs=1 -timeout=0` in the background, and takes `sudo gdb -p $pid -batch -ex "thread apply all bt 60"` every few seconds (install `gdb` with apt). Decode the crash artifact's base64 into the workflow.

The push also triggers ci.yml and fuzz.yml on that branch: cancel those runs with `gh run cancel` right after pushing. Actions caches are scoped to their branch (plus the default branch), so a diag branch cannot restore the fuzz corpus (`corpus-<target>-*`) saved by the work branch's fuzz runs: a diagnostic that replays the corpus must run on the work branch itself, as a temporary commit force-pushed away afterwards. Never merge the diagnostic commits; rebuild the work branch without them (cherry-pick, no interactive rebase) and delete the diag branch afterwards. Read job logs with `gh api repos/ef3d0c3e/cas/actions/jobs/<id>/logs` (`gh run view --log` is sometimes empty).

The eval target's bytes decode only through its own `Unstructured` generator, so to see the expression add to the diag branch's `eval.rs` an `eprintln!` of the input (and of `report.trace.steps`) behind an env var (`CAS_FUZZ_SHOW`). For an out-of-memory (CAS-201), gdb samples showed the spot (`ball_op` → `from_rational`) but not why; what did: temporarily instrumenting the suspect function to print its operands and `panic!` past a threshold (ball precision > 100000 bits), which named the operand in one run. Make the diag branch with `git worktree add -b <name> <scratchpad>/diag <commit>` so the main tree's uncommitted work stays untouched; `git worktree remove --force` and `git push origin --delete <name>` afterwards.

Faster, and used for CAS-516: decode and time the input locally, no fuzzing run. Download the artifact (`gh run download <run> -D <scratch>/art`), copy `fuzz/fuzz_targets/eval.rs` into a scratch bin crate (deps `arbitrary`, `cas-core` by path), strip `#![no_main]`/`fuzz_target!` for a `main` that reads the file (`.leak()` the bytes), calls `context` then `expr(u, 5, …)`, prints them and times `ctx.eval`; add `eprintln!`s in `context` for the rules, limits and attributes. Then `gdb -p <pid> -batch -ex "thread apply all bt 60"` on the local process every few seconds names the hot function. Inserting names into the eval target's BUILTINS re-maps the whole cached corpus, so a timeout right after adding builtins is usually an old slowness newly reached, not the new code.

Before a gdb round trip, read the timeout's own trace: it is unsymbolized (`-s none` build), but its shape tells a lot. ~100 frames alternating between two return addresses meant one function calling itself from two sites per level, an exponential double visit (CAS-176: `fold_arithmetic` folded a power's exponent twice, 2^depth on towers). Decode the failing input from the log's `Output of std::fmt::Debug:` byte array (Python, bytes → utf-8) and look for its repeated structure. Also compare with master's last fuzz log: slow units (`Slowest unit: N s`) there mean the problem predates the branch.

Related: [[ci-on-github]].

## ci-on-github

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

## work-on-master

---
name: work-on-master
description: "For the CAS repo, commit and push directly to master (no feature branches); a failing build on master is acceptable while nothing is released"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-25T12:06:40.717Z
---

Work directly on master in ef3d0c3e/cas: commit there and push, without feature branches or asking before merging. Pushing a build that fails CI to master is fine for now.

**Why:** the user said so on 2026-10-25 ("it's better if you work directly into master for now. There aren't any release out, so pushing a failing build to master is fine for now"). Earlier sessions used branches (cas-154, cas-162, cas-170) and asked before each merge; that friction is no longer wanted.

**How to apply:** commit each issue on master and push; CI and fuzz run on every push ([[ci-on-github]]). Still fix what CI or the fuzzer finds right away. Revisit once the project has a release: "for now" means this can change.

## focus-on-cas-features

---
name: focus-on-cas-features
description: "Priorities — work on the CAS itself (features, builtins, correctness); low-level kernel speedups like CAS-30 (SIMD/unsafe) are optional and deferred"
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-25T16:01:31.090Z
---

On 2026-10-25 the user said of CAS-30 (SIMD/assembly multiplication kernels, which need `unsafe`): "Skip it, CAS-30 is completely optional, we'll do it later if it becomes a requirement, right now focus on the CAS itself."

**Why:** the user wants progress on what the CAS can do, not on squeezing the arithmetic kernels; the workspace keeps `unsafe_code = "deny"` (only cas-giac's loader exempt).

**How to apply:** leave CAS-30 in the backlog; when tickets run out, pick CAS features (builtins, algorithms the language exposes, giac parity) over kernel micro-optimisation. Don't introduce `unsafe` without asking. See [[benchmarks-on-ci]].

## autonomous-areas

---
name: autonomous-areas
description: Move from one CAS area to the next in my recommended order without asking; ask only when a brand new Lific module is needed
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-28T11:28:42.519Z
---

When an area (CAS-PLAN-n) is done, go straight on to the next one in my recommended order instead of asking which area comes next. Ask the user only when the work needs a brand new module (a Lific module that does not exist yet); the existing ones are Algebra, Balls, Calculus, Evaluator, Expressions, Infrastructure, Integers, Language, Polynomials, and Statistics (approved 2026-10-28).

**Why:** after Laplace transforms the user answered "Go with your recommended order. Don't notify me every time, you can keep working until a brand new module needs work."

**How to apply:** at an area's end, close its tickets and doc as usual, then write the next design doc and plan without an AskUserQuestion; give short progress notes, not area-choice questions. The order recommended on 2026-10-28: inequalities (Algebra), systems of ODEs (Calculus), Gröbner bases (Polynomials, with `eliminate`), all done; then statistics. On 2026-10-28 the user chose "Statistics (Recommended)": a new Statistics module, descriptive statistics and regressions first, then distributions (which need incomplete gamma and beta). See [[work-on-master]] and [[ci-on-github]].

## lific-is-remote-only

---
name: lific-is-remote-only
description: Lific issue tracker lives at https://lific.pundalik.org; never run the local lific server; CLI needs --backend http
metadata:
  node_type: memory
  type: project
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-24T08:49:37.734Z
---

The Lific issue tracker for this machine is the remote instance at https://lific.pundalik.org. The local Lific server must never run (the user said so on 2026-10-24; its systemd user service was uninstalled then). The local database at ~/.local/share/lific/lific.db is empty — the real data (projects CAS, SRV, AUDIO) is on the remote.

**Why:** The local instance is a stale leftover; starting it serves empty data and misleads tooling.

**How to apply:** MCP connects via `.mcp.json` (url must include the `/mcp` path: https://lific.pundalik.org/mcp). For the CLI, always pass `--backend http` — `LIFIC_URL` and `LIFIC_API_KEY` are set in `.claude/settings.local.json` env, but the backend never defaults to http, e.g. `lific issue list --project CAS --backend http --json`. Never run `lific start`, `lific service install`, or `lific init` here.
