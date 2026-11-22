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
