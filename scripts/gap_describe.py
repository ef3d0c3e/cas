"""gap_describe.py -- read giac_gap.sh's lists against aide_cas. Run from the repo root.

  missing <gapdir> [regex]       each missing aide entry: names | signature | English line
  all <gapdir> <regex>           every aide entry matching regex, marked M (missing) or C (covered)
  lines <olddir> <gapdir> <a-b,c> entries at those lines of an older entries_missing.txt (R0-004's
                                 .run/logs/gen1-gap), marked M if still missing in <gapdir>
<gapdir> is giac_gap.sh's outdir (entries_missing.txt). From the 1.6 planner's describe*.py, area_lines.py.
"""
import re
import sys

sys.path.insert(0, "scripts")
import aide  # noqa: E402


def missing_set(d):
    with open(f"{d}/entries_missing.txt", encoding="utf-8") as f:
        return {ln.strip() for ln in f}


def rows():
    for e in aide.entries():
        key = " ".join(e["names"])
        yield key, key + " | " + e["sig"] + " | " + e["en"]


def main(argv):
    if len(argv) >= 2 and argv[0] == "missing":
        miss = missing_set(argv[1])
        pat = re.compile(argv[2], re.I) if len(argv) > 2 else None
        for key, text in rows():
            if key in miss and (pat is None or pat.search(text)):
                print(text[:230])
    elif len(argv) == 3 and argv[0] == "all":
        miss, pat = missing_set(argv[1]), re.compile(argv[2], re.I)
        for key, text in rows():
            if pat.search(text):
                print(("M " if key in miss else "C ") + text[:200])
    elif len(argv) == 4 and argv[0] == "lines":
        with open(f"{argv[1]}/entries_missing.txt", encoding="utf-8") as f:
            old = f.read().splitlines()
        now, out = missing_set(argv[2]), []
        for part in argv[3].split(","):
            a, _, b = part.strip().partition("-")
            for i in range(int(a), int(b or a) + 1):
                e = old[i - 1].strip()
                out.append(("M " if e in now else "C ") + e)
        print(" | ".join(out))
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
