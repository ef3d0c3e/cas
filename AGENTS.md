<!-- lific:begin -->
## Issue tracking: Lific

This project tracks its work in Lific.

> Replace `APP` with this repo's project identifier. Run `lific project list` to find it.

Use the Lific MCP tools for issues, pages, and plans when they are connected. Use the `lific` CLI for work that needs this machine. The MCP server can run on another host, so it cannot see your files, working directory, or git history. The CLI also has every data command, for when MCP is not connected. Add `--json` for machine-readable output:

```bash
lific issue list --project APP --json      # browse issues
lific issue get APP-1                       # one issue with relations
lific issue update APP-1 --status done      # close an issue
lific search "auth flow" --project APP      # full-text search
```

Rules:
- When you finish an issue, set its status to `done`.
- File each issue in the right module, and keep it to one change.
- Use a plan for work that spans sessions.
- An issue assigned to a person (`for:` in lists) is theirs. Take only unassigned issues.
<!-- lific:end -->

## Tests and fuzzing run on GitHub Actions only

Never run `cargo test`, `cargo fuzz` or the benchmarks on this machine. They run
on GitHub Actions: `.github/workflows/ci.yml` (format, clippy, tests, a nightly
check of the fuzz crate — on every push) and `.github/workflows/fuzz.yml` (every
libFuzzer target — twice a week, or on demand with
`gh workflow run fuzz.yml -f time_per_target=<seconds> [-f targets="<names>"]`).
Local `cargo check` and `cargo build` are fine. A fuzz finding's crash artifact
is in the failed run's uploaded artifacts; fix it and add a regression test to
the owning crate's test suite, which CI then runs.
