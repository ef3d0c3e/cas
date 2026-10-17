<!-- lific:begin -->
## Issue tracking: Lific

This project uses **Lific** for issue tracking and project management (local-first, single-binary, SQLite-backed).

> Replace `APP` with this repo's project identifier — run `lific project list` to find it.

**Preferred access — MCP.** If a Lific MCP server is configured in your client (see the tools/MCP config in this repo or your global config), use it directly: list/create/update issues, pages, and plans through the Lific tools.

**CLI fallback** (works without MCP; add `--json` for machine-readable output):

```bash
lific issue list --project APP --json      # browse issues
lific issue get APP-1                       # one issue with relations
lific issue update APP-1 --status done      # close an issue
lific search "auth flow" --project APP      # full-text search
```

**Conventions:**
- Mark an issue `done` as soon as you finish its work — don't leave it open.
- Group related work with **modules**; fit each issue into the right module.
- For multi-session work, use a **plan** so the next session can resume.
- Issues are self-contained work items; keep scope tight.
<!-- lific:end -->
