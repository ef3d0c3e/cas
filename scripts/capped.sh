#!/usr/bin/env bash
# Runs a command in a memory-capped systemd scope with a timeout.
#
# Every long test, fuzz or bench run goes through this: binary splitting,
# precision loops and Toom/NTT recursions can blow up memory when they go
# wrong, and the machine is shared with other jobs. The cap covers the whole
# process tree (cargo, rustc, linkers, test binaries) and swap is disabled
# inside the scope, so a runaway is killed instead of thrashing the machine.
#
#   scripts/capped.sh cargo test -p cas-num
#   CAP_MEM=4G CAP_TIME=600 scripts/capped.sh cargo fuzz run -O ball
#
# CAP_MEM   memory limit of the scope (default 8G)
# CAP_TIME  timeout in seconds (default 3600)
# CAP_JOBS  cargo build jobs, unless CARGO_BUILD_JOBS is set (default 6)
set -euo pipefail

if [ "$#" -eq 0 ]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi

exec systemd-run --user --scope --quiet \
    -p MemoryMax="${CAP_MEM:-8G}" -p MemorySwapMax=0 \
    --setenv=CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-${CAP_JOBS:-6}}" \
    -- timeout --kill-after=30 "${CAP_TIME:-3600}" "$@"
