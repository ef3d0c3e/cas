"""check.py -- giac's mathml(X) (probe.out) vs mathml_proto on X's evaluated shape (T5 proto parse, T5 check.py's
overrides). Shape rule added for mathml.cc: an evaluated sqrt(a) is giac's pow(a,1/2) (markup.cc reads both alike)."""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROTO = os.path.join(os.path.dirname(HERE), "proto")
sys.path.insert(0, HERE)
sys.path.insert(0, PROTO)
_spec = importlib.util.spec_from_file_location("t5check", os.path.join(PROTO, "check.py"))
t5check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t5check)
from parse import parse  # noqa: E402
from markup import Vect, Symb, Frac, Int  # noqa: E402
from mathml_proto import mathml  # noqa: E402

T5 = os.path.dirname(HERE)
T4 = os.path.join(os.path.dirname(T5), "latex")


def eval_shape(g):
    if isinstance(g, Vect):
        return Vect([eval_shape(x) for x in g.items], g.sub)
    if isinstance(g, Symb):
        a = eval_shape(g.arg)
        if g.op == "sqrt":
            return Symb("pow", Vect([a, Frac(Int(1), Int(2))], 'seq'))
        return Symb(g.op, a)
    return g


def nest(g):
    """giac's parser makes a*b*c prod(prod(a,b),c): aplatir_fois_plus then multiplies (gen.cc:12521)."""
    if isinstance(g, Vect):
        return Vect([nest(x) for x in g.items], g.sub)
    if isinstance(g, Symb):
        a = nest(g.arg)
        if g.op == "prod" and isinstance(a, Vect) and len(a.items) > 2:
            acc = a.items[0]
            for x in a.items[1:]:
                acc = Symb("prod", Vect([acc, x], 'seq'))
            return acc
        return Symb(g.op, a)
    return g


# extra rows (extra.cas): shape text and mode
EXTRA = {"sin": ("'sin'", 'eval'), "exp": ("'exp'", 'eval'), "1e-300*x": ("1e-300*x", 'eval')}


def shape_of(inp, evals):
    mode = 'quote' if inp.startswith("quote(") else 'eval'
    if inp in t5check.DBL:
        return t5check.DBL[inp]
    if inp in EXTRA:
        return parse(*EXTRA[inp])
    if inp == "quote(x and y)":  # parse.py has no 'and'
        return Symb("and", Vect([parse("x"), parse("y")], 'seq'))
    if inp not in evals and mode == 'quote':
        return nest(parse(inp[len("quote("):-1], mode))
    g = parse(t5check.SHAPE.get(inp, evals[inp]), mode)
    return eval_shape(g) if mode == 'eval' else nest(g)


def main():
    t4 = t5check.read_pairs(os.path.join(T4, "probe.out"))
    evals = dict(t4[0::2])
    rows = t5check.read_pairs(os.path.join(HERE, "probe.out"))
    rows += [r for r in t5check.read_pairs(os.path.join(HERE, "extra.out")) if r[0] != "mathml(quote(not x))"]
    ok = 0
    out = []
    for a, b in rows:
        inp = a[len("mathml("):-1]
        want = t5check.unquote(b)
        try:
            ours = mathml(shape_of(inp, evals)).replace("\n", " | ")
        except Exception as e:  # noqa: BLE001
            ours = f"<{type(e).__name__}: {e}>"
        good = ours == want
        ok += good
        out.append((inp, want, "ok" if good else "FAIL"))
        if not good:
            i = 0
            while i < min(len(want), len(ours)) and want[i] == ours[i]:
                i += 1
            print(f"FAIL {inp}\n  shape: {t5check.SHAPE.get(inp, evals.get(inp))}\n"
                  f"  giac: ...{want[max(0, i - 40):i + 200]}\n  ours: ...{ours[max(0, i - 40):i + 200]}")
    with open(os.path.join(HERE, "table.tsv"), "w", encoding="utf-8") as f:
        f.write("input\tgiac\tok\n")
        for r in out:
            f.write("\t".join(r) + "\n")
    print(f"{ok}/{len(rows)} ok")
    return 0 if ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
