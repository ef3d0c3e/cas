#!/usr/bin/env bash
# giac_gap.sh <outdir> -- giac's commands (aide_cas) with and without a name in
# crates/cas-lang/names/en.txt ([builtins] up to [types]). Run from the repo root.
# Writes giac.txt, ours.txt, missing.txt, present.txt (by name) and
# entries_covered.txt, entries_missing.txt (by aide_cas entry, aliases included).
set -eu
[ $# -eq 1 ] || { echo "usage: $0 <outdir>" >&2; exit 2; }
out=$1
mkdir -p "$out"
grep "^# " /usr/share/giac/aide_cas | cut -c3- | tr -d "\r" | tr " " "\n" | grep -v "^$" | sort -u > "$out/giac.txt"
sed -n '/^\[builtins\]/,/^\[types\]/p' crates/cas-lang/names/en.txt | grep -v '^#' | grep ':' \
  | cut -d: -f2- | tr ' ' '\n' | grep -v '^$' | sort -u > "$out/ours.txt"
comm -23 "$out/giac.txt" "$out/ours.txt" > "$out/missing.txt"
comm -12 "$out/giac.txt" "$out/ours.txt" > "$out/present.txt"
wc -l "$out"/giac.txt "$out"/ours.txt "$out"/missing.txt "$out"/present.txt
# entry level: an aide_cas entry (one `# n1 n2 …` line, aliases included) is covered when any of its names is ours
grep "^# " /usr/share/giac/aide_cas | cut -c3- | tr -d "\r" | awk -v ours="$out/ours.txt" -v out="$out" '
  BEGIN { while ((getline n < ours) > 0) have[n] = 1 }
  { hit = 0; for (i = 1; i <= NF; i++) if ($i in have) hit = 1; print > (out (hit ? "/entries_covered.txt" : "/entries_missing.txt")) }'
wc -l "$out"/entries_covered.txt "$out"/entries_missing.txt
