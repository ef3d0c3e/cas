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
