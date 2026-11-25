"""check.py -- print each T4 row with the MathML port, compare with giac's export_mathml (display, content, default).

Reads ../probe.out (scripts/caseval.py on ../probe.cas: export_mathml(X,display), (X,content), (X) per input)
and ../../latex/probe.out (the evaluated forms); writes ../table.tsv. Also re-checks T4's latex rows.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from markup import (export_latex, export_mathml, export_mathml_content,  # noqa: E402
                    export_mathml_presentation)
from parse import parse  # noqa: E402
from markup import Dbl  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
T5 = os.path.dirname(HERE)
T4 = os.path.join(os.path.dirname(T5), "latex")

# as T4/proto/check.py (NOTES.md "shape overrides")
SHAPE = {
    "x=1": "x=1",
    "a=b+1": "a=b+1",
    "[x,y]=[1,2]": "[x=1,y=2]",
    '"a\\nb"': '"a\\nb"',
    "quote(2*x*inv(3))": "2*x*inv(3)",
    "quote(x*inv(3)*y)": "x*inv(3)*y",
    "quote(x*inv(2)*inv(3))": "x*inv(2)*inv(3)",
    "quote(inv(x))": "inv(x)",
    "quote(x*y*inv(z))": "x*y*inv(z)",
    "quote(2*inv(4)*x)": "2*inv(4)*x",
}

# giac doubles printed without a point (parse.py reads them as integers; latex is the same)
DBL = {"1e10": Dbl("10000000000"), "12345678901.0": Dbl("12345678901")}


def unquote(ans):
    if len(ans) >= 2 and ans[0] == '"' and ans[-1] == '"':
        return ans[1:-1].replace('""', '"')
    return None


def read_pairs(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            a, b = line.rstrip("\n").split("  =>  ", 1)
            out.append((a, b))
    return out


def main():
    t4 = read_pairs(os.path.join(T4, "probe.out"))
    evals = dict(t4[0::2])
    lat = {a[len("latex("):-1]: b for a, b in t4[1::2]}
    probe = read_pairs(os.path.join(T5, "probe.out"))
    rows, ok, lat_ok, lat_n = [], 0, 0, 0
    for k in range(0, len(probe), 3):
        (pd, d), (pc, c), (pp, p) = probe[k:k + 3]
        inp = pp[len("export_mathml("):-1]
        assert pd == f"export_mathml({inp},display)" and pc == f"export_mathml({inp},content)", (pd, inp)
        ev = evals[inp]
        mode = 'quote' if inp.startswith("quote(") else 'eval'
        shape = SHAPE.get(inp, ev)
        g = DBL[inp] if inp in DBL else parse(shape, mode)
        wl = unquote(lat[inp])
        if wl is not None:
            lat_n += 1
            lat_ok += export_latex(g) == wl
        for kind, ans, fn in (("display", d, export_mathml_presentation), ("content", c, export_mathml_content),
                              ("default", p, None)):
            want = unquote(ans)
            if want is None:
                print(f"skip (giac error): {kind} {inp}")
                continue
            try:
                if fn is None:
                    # annotation: giac's print of the evaluated input, taken from giac's own text
                    a0 = want.find("<annotation encoding='giac'>") + len("<annotation encoding='giac'>")
                    printed = want[a0:want.rfind("</annotation>")]
                    ours = export_mathml(g, "\0")
                    ours = ours.replace("\0", printed)
                else:
                    ours = fn(g)
            except Exception as e:  # report, count as a failure
                ours = f"<{type(e).__name__}: {e}>"
            ours = ours.replace("\n", " | ")  # caseval shows a newline as ' | ' (C0047)
            good = ours == want
            ok += good
            rows.append((inp, mode, shape, kind, want, "ok" if good else "FAIL"))
            if not good:
                print(f"FAIL {kind} {inp}\n  shape: {shape}\n  giac: {want}\n  ours: {ours}")
    with open(os.path.join(T5, "table.tsv"), "w", encoding="utf-8") as f:
        f.write("input\tmode\tshape\tkind\tgiac\tok\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"latex regression {lat_ok}/{lat_n}")
    print(f"{ok}/{len(rows)} ok")
    return 0 if ok == len(rows) and lat_ok == lat_n else 1


if __name__ == "__main__":
    sys.exit(main())
