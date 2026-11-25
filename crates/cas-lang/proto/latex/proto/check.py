"""check.py -- parse each giac evaluated form, print it with the markup.py port, compare with giac's latex.

Reads ../probe.out (scripts/caseval.py on ../probe.cas: each input X then latex(X)); writes ../table.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from markup import export_latex  # noqa: E402
from parse import parse  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
T4 = os.path.dirname(HERE)

# Rows whose caseval evaluated print is not the shape latex() sees (see NOTES.md "shape overrides"):
# a top-level `a=b` prints its rhs (caseval treats it as a statement), the newline string prints " | ".
SHAPE = {
    "x=1": "x=1",
    "a=b+1": "a=b+1",
    "[x,y]=[1,2]": "[x=1,y=2]",
    '"a\\nb"': '"a\\nb"',
    # quote() prints prod(a,inv(b)) and division(a,b) both as a/b: take the typed shape
    "quote(2*x*inv(3))": "2*x*inv(3)",
    "quote(x*inv(3)*y)": "x*inv(3)*y",
    "quote(x*inv(2)*inv(3))": "x*inv(2)*inv(3)",
    "quote(inv(x))": "inv(x)",
    "quote(x*y*inv(z))": "x*y*inv(z)",
    "quote(2*inv(4)*x)": "2*inv(4)*x",
}


def unquote(ans):
    if len(ans) >= 2 and ans[0] == '"' and ans[-1] == '"':
        return ans[1:-1].replace('""', '"')
    return None


def main():
    with open(os.path.join(T4, "probe.out"), encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f]
    rows, ok = [], 0
    for a, b in zip(lines[0::2], lines[1::2]):
        inp, ev = a.split("  =>  ", 1)
        linp, lat = b.split("  =>  ", 1)
        assert linp == f"latex({inp})", (linp, inp)
        want = unquote(lat)
        if want is None:
            print(f"skip (giac error): {inp}")
            continue
        mode = 'quote' if inp.startswith("quote(") else 'eval'
        shape = SHAPE.get(inp, ev)
        try:
            ours = export_latex(parse(shape, mode))
        except Exception as e:  # report, count as a failure
            ours = f"<{type(e).__name__}: {e}>"
        good = ours == want
        ok += good
        rows.append((inp, ev, want, ours, "ok" if good else "FAIL"))
        if not good:
            print(f"FAIL {inp}\n  eval: {ev}\n  giac: {want}\n  ours: {ours}")
    with open(os.path.join(T4, "table.tsv"), "w", encoding="utf-8") as f:
        f.write("input\tgiac_evaluated\tgiac_latex\tours\tok\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{ok}/{len(rows)} ok")
    return 0 if ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
