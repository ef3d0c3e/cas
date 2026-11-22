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
