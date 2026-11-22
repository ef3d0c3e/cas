# How we work — cas

<!-- The card: standing facts for every agent. Never appended to; the cap is
     card.max_chars in .claude/pa.json (default 10000), enforced by the archive.
     Each role prints its slice with `tools/card.py slice <role>`.
     Written at install (interview) and edited in the same task as whatever changed it (H7). -->

## Project <!-- roles: claude expert coder planner review auditor curator -->
cas — <one line: what this project is -- edit me>. Constitution `PROJECT_CONTEXT.md`; roadmap `GENERATION_PLAN.md`;
rules `rules/INDEX.md`; techniques `cookbook/INDEX.md`; ops detail `docs/ops/INDEX.md`.

## Developer <!-- roles: claude planner review auditor curator -->
Who: ef3d0c3e. Experience: advanced. Domain: cas.
Preferences: recommendations, not questions. Plain-English recaps at phase end. A notification
whenever anything waits on them. Autonomy: full inside an approved plan; stop only at the two gates. Notification channel: toast.
The developer pushes; agents never do. They ratify rules at the next plan approval.

## Tools <!-- roles: expert coder planner review auditor curator -->
PY = /usr/bin/python3

| Name | Command | Purpose |
|---|---|---|
| launch | `PY tools/launch.py --state-json` | the state detector: the mod's read-only report; never typed by the developer |
| plan_edit | `PY tools/plan_edit.py` | the only writer of `PHASE_PLAN.md` (status, Changes, add/reopen) |
| status | `PY tools/status.py` | `.run/status.json` for the statusline; INBOX consume; waiting flags |
| task_log | `PY tools/task_log.py` | lint and finish a task summary |
| research_add | `PY tools/research_add.py` | allocate a report id, write the index line |
| rules_add | `PY tools/rules_add.py` | add, supersede, promote, retire a rule |
| skill_add | `PY tools/skill_add.py` | turn a `workflow:` gotcha into `.claude/skills/<name>/SKILL.md` |
| cookbook_add | `bash tools/cookbook_add.sh` | add a cookbook entry and its index line |
| phaseend_index | `PY tools/phaseend_index.py` | assemble and archive a PhaseEnd; `history` rebuilds `docs/PHASE_HISTORY.md` |
| genend_index | `PY tools/genend_index.py` | assemble a GenerationEnd |
| commit_task | `bash tools/commit_task.sh` | the only commit path; explicit paths, no trailers, never pushes |
| run | `bash tools/run.sh` | any command that may print >40 lines |
| curate | `PY tools/curate.py` | `memory --route\|--reindex`, `rules --reseed\|--demote`, `cookbook --demote`, `claude-md\|card --refresh`, `living-docs --entry`, `how-we-work --retire`; `--dry-run` first |
| <!-- TODO PROJECT_TOOL: no source found --> | `<!-- TODO PROJECT_TOOL_COMMAND: no source found -->` | <!-- TODO PROJECT_TOOL_PURPOSE: no source found --> |

## Skills <!-- roles: expert planner -->
<!-- one line per captured workflow; the SKILL.md is the canonical text -->
- <!-- TODO SKILL_NAME: no source found --> — <what it automates, when to invoke it>

## Paths <!-- roles: claude expert coder planner review auditor curator -->
- Oracles (the ground truth X3 names): <!-- TODO ORACLES: no source found -->
- Data: <!-- TODO DATA_PATHS: no source found -->
- Generated (never hand-edited, H1): <!-- TODO GENERATED_PATHS: no source found -->
- Hand-edited: <!-- TODO HAND_EDITED_PATHS: no source found -->
- Scratch: all of a phase's work in `.run/` goes under `.run/<phase>/` (gitignored; never the system temp), nowhere
  else in `.run/`. It is a temp folder only, never a home for long-term content: what the project needs lives in its
  proper place in the repo. The phase's cleanup task checks it and
  the close deletes it; a tool deletes its own temp folders when it finishes.

## Build / run / test <!-- roles: claude expert coder planner review auditor curator -->
- Build: `cargo build`
- Run: `<!-- TODO RUN_COMMAND: no source found -->`
- Test / the gate: `cargo test` — green means <!-- TODO GREEN_MEANS: no source found -->
- Modes: <!-- TODO TEST_MODES: no source found -->
Long gates go through `tools/run.sh` with a raised timeout, never a poll loop.

## Conventions & house style <!-- roles: claude expert coder planner review auditor curator -->
- <!-- TODO PROJECT_CONVENTION: no source found -->
<!-- project-specific only; the general house style is in the project-architect skill -->

## Environment <!-- roles: claude expert coder planner review auditor curator -->
- OS / shells: Linux / bash
- Python: /usr/bin/python3 (3.14.10)
- Pins: <!-- TODO PINS_SUMMARY: no source found --> — detail in `docs/ops/INDEX.md`
- Harness gotchas that bite here: <!-- TODO HARNESS_GOTCHAS: no source found -->

## Docs map <!-- roles: expert coder planner review auditor curator -->
- `docs/ops/INDEX.md` — setup, env, pins, per-topic ops notes
- `cookbook/INDEX.md` — techniques, grep by tag
- `rules/INDEX.md` — full rule texts
- `phase-ends/TASK_INDEX.md`, `phase-ends/RESEARCH_INDEX.md` — what was done and what was learned
- `docs/research-archive/` — migrated legacy reports
- `docs/retired/` — everything moved out of the load order

## Living docs <!-- roles: claude expert planner curator -->
- (none yet: `doc`: what it answers; when to update it)
