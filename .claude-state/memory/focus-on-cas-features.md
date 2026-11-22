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
