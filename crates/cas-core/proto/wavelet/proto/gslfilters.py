"""Parse GSL's wavelet filter arrays (read-only C sources) into {name: [decimal str]} and build filter tables. Origin: phase 1.7, T9.c1 scratch."""
import re, sys

SQRT1_2 = "0.70710678118654752440084436210"  # gsl_math.h M_SQRT1_2


def arrays(path):
    text = open(path, encoding="utf-8").read()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    out = {}
    for m in re.finditer(r"static const double (\w+)\[(\d+)\]\s*=\s*\{(.*?)\};", text, flags=re.S):
        vals = []
        for tok in m.group(3).split(","):
            t = tok.strip().replace(" ", "")
            if t == "M_SQRT1_2":
                t = SQRT1_2
            elif t == "-(M_SQRT1_2)":
                t = "-" + SQRT1_2
            if t:
                vals.append(t)
        assert len(vals) == int(m.group(2)), m.group(1)
        out[m.group(1)] = vals
    return out


BSPLINE = {  # member: (nc, h1, (g1, off), (h2, off), g2) from bspline.c bspline_init
    103: (6, "h1_103", ("g1_1", 2), ("h2_1", 2), "g2_103"),
    105: (10, "h1_105", ("g1_1", 0), ("h2_1", 0), "g2_105"),
    202: (6, "h1_202", ("g1_2", 6), ("h2_2", 6), "g2_202"),
    204: (10, "h1_204", ("g1_2", 4), ("h2_2", 4), "g2_204"),
    206: (14, "h1_206", ("g1_2", 2), ("h2_2", 2), "g2_206"),
    208: (18, "h1_208", ("g1_2", 0), ("h2_2", 0), "g2_208"),
    301: (4, "h1_301", ("g1_3", 8), ("h2_3", 8), "g2_301"),
    303: (8, "h1_303", ("g1_3", 6), ("h2_3", 6), "g2_303"),
    305: (12, "h1_305", ("g1_3", 4), ("h2_3", 4), "g2_305"),
    307: (16, "h1_307", ("g1_3", 2), ("h2_3", 2), "g2_307"),
    309: (20, "h1_309", ("g1_3", 0), ("h2_3", 0), "g2_309"),
}


def table(d):
    """{(type, k): (nc, h1, g1, h2, g2)} with each filter a list of nc decimal strings."""
    a = {}
    for f in ("daubechies.c", "haar.c", "bspline.c"):
        a.update(arrays(d + "/" + f))
    t = {("haar", 2): (2, a["ch_2"], a["cg_2"], a["ch_2"], a["cg_2"])}
    for k in range(4, 21, 2):
        t[("daubechies", k)] = (k, a[f"h_{k}"], a[f"g_{k}"], a[f"h_{k}"], a[f"g_{k}"])
    for k, (nc, h1, (g1, o1), (h2, o2), g2) in BSPLINE.items():
        t[("bspline", k)] = (nc, a[h1][:nc], a[g1][o1:o1 + nc], a[h2][o2:o2 + nc], a[g2][:nc])
        assert len(a[h1]) == nc and len(a[g2]) == nc and len(t[("bspline", k)][2]) == nc
    return t


if __name__ == "__main__":
    t = table(sys.argv[1])
    for key, v in t.items():
        print(key, v[0], v[1][:2])
