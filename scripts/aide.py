#!/usr/bin/env python3
"""aide.py <name>... -- giac help entries (aide_cas): names, signature, English line, examples.

Reads /usr/share/giac/aide_cas. With --cas, prints only the examples, one per line (a caseval.py input file).
"""
import re
import sys

AIDE = "/usr/share/giac/aide_cas"


# A summary line: a language digit (1 fr, 2 en, 3 es, 4 el, 5 de, 8 zh), then its text.
SUMMARY = re.compile(r"^([1-9])(?: (.*))?$")
# A signature line: `0`, then the argument types (possibly empty).
SIGNATURE = re.compile(r"^0(?: (.*))?$")
# A related-names line: `-N`, then words (`-1 atan or Opt`, `-2 :=`, `-1 ` empty).
RELATED = re.compile(r"^-\d+((?:\s+\S+)*)\s*$")


def entries(path=AIDE):
    """aide_cas's entries, read by position (one reader for aide.py and doc_table.py).

    An entry runs from a `# names` line to the next one; blank lines inside it are
    skipped (UTPF and limit have some). Its header is its leading tagged lines, each tag
    once: summary digits and the `0` line, in either order (dwt puts `0` first), and a
    summary's untagged continuation lines when a tagged line follows them (chrem,
    l1norm, ldl, pmin wrap theirs). The first other line starts the examples, and
    every later line is an example, so an example starting with a digit stays an
    example (tar's late `1 fopen`). A `-N` line is a related-names line wherever it stands (posubLMQ lists them after
    its examples). Within an entry the first `0` and the first of each summary digit
    win (a repeated tag ends the header); a repeated `0` line right after it is kept in `dups` as (line number, text).
    Each entry: names, line, sig (None if absent), summaries {digit: text},
    en (summaries["2"] or None), related, ex, dups.
    """
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    out = []
    cur = None
    for i, line in enumerate(lines, 1):
        if line.startswith("# "):
            cur = {"names": line[2:].split(), "line": i, "sig": None, "summaries": {},
                   "related": [], "ex": [], "dups": [], "_raw": []}
            out.append(cur)
        elif cur is not None and line.strip():
            cur["_raw"].append((i, line))
    for e in out:
        raw = e.pop("_raw")
        seen = set()      # header tags read: summary digits and "0"
        examples = False
        cont = None       # the summary digit an untagged line would continue
        for k, (i, line) in enumerate(raw):
            m = RELATED.match(line)
            if m:
                e["related"] += m.group(1).split()
                continue
            if not examples:
                s = SUMMARY.match(line)
                g = SIGNATURE.match(line)
                if s and s.group(1) not in seen:
                    d = s.group(1)
                    e["summaries"][d] = (s.group(2) or "").strip()
                    seen.add(d)
                    cont = d
                    continue
                if g and "0" not in seen:
                    e["sig"] = (g.group(1) or "").strip()
                    seen.add("0")
                    cont = None
                    continue
                if g and SIGNATURE.match(raw[k - 1][1]):
                    e["dups"].append((i, line))  # border repeats its `0` line
                    continue
                if cont is not None and _continues(raw, k, seen):
                    e["summaries"][cont] += " " + line.strip()
                    continue
                examples = True
            e["ex"].append(line)
        e["en"] = e["summaries"].get("2")
    return out


def _continues(raw, k, seen):
    """Whether untagged raw[k] wraps the summary before it: the first tagged line after
    it is a summary or signature line with a tag not read yet (inertia wraps over
    three lines before its `0`)."""
    for _, line in raw[k + 1:]:
        s = SUMMARY.match(line) or SIGNATURE.match(line)
        if s:
            return line[0] not in seen
        if RELATED.match(line):
            return False
    return False


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
                print("#", " ".join(e["names"]), "|", e["sig"] or "", "|", (e["en"] or "")[:160],
                      "|", " ".join(e["related"]))
                for x in e["ex"]:
                    print("   ", x)


if __name__ == "__main__":
    main()
