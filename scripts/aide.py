#!/usr/bin/env python3
"""aide.py <name>... -- giac help entries (aide_cas): names, signature, English line, examples.

Reads /usr/share/giac/aide_cas. With --cas, prints only the examples, one per line (a caseval.py input file).
"""
import re
import sys

AIDE = "/usr/share/giac/aide_cas"


def entries():
    cur = None
    with open(AIDE, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if line.startswith("# "):
                if cur:
                    yield cur
                cur = {"names": line[2:].split(), "sig": "", "en": "", "ex": []}
            elif cur is not None and line.strip():
                if line.startswith("0 "):
                    cur["sig"] = line[2:]
                elif line.startswith("2 "):
                    cur["en"] = line[2:]
                elif re.match(r"^-?\d+ ", line) or re.match(r"^-\d+$", line):
                    continue
                else:
                    cur["ex"].append(line)
    if cur:
        yield cur


def main():
    args = sys.argv[1:]
    cas = "--cas" in args
    want = [a for a in args if a != "--cas"]
    for e in entries():
        if any(n in want for n in e["names"]):
            if cas:
                for x in e["ex"]:
                    print(x)
            else:
                print("#", " ".join(e["names"]), "|", e["sig"], "|", e["en"][:160])
                for x in e["ex"]:
                    print("   ", x)


if __name__ == "__main__":
    main()
