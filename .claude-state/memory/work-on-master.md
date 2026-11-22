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
