"""Prototype of dwt/idwt (GSL wavelet/dwt.c port) checked against giac rows. Run: python3 -I proto_dwt.py <gsl wavelet dir>. Origin: phase 1.7, T9.c1 scratch."""
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "proto"))
from gslfilters import table

T = table(sys.argv[1])


def wavelet(kind, k, centred, conv):
    nc, h1, g1, h2, g2 = T[(kind, k)]
    off = (nc >> 1) if centred else 0  # haar centred: 1 = 2>>1; daubechies member>>1 = nc>>1; bspline nc>>1
    return nc, off, [conv(x) for x in h1], [conv(x) for x in g1], [conv(x) for x in h2], [conv(x) for x in g2]


def step(w, a, n, fwd, zero):
    """dwt.c:dwt_step on a[0..n) (a list, in place)."""
    nc, off, h1, g1, h2, g2 = w
    s = [zero] * n
    nmod = nc * n - off
    nh = n >> 1
    if fwd:
        for ii, i in enumerate(range(0, n, 2)):
            h = zero; g = zero
            for k in range(nc):
                jf = (n - 1) & (i + nmod + k)
                h += h1[k] * a[jf]
                g += g1[k] * a[jf]
            s[ii] += h
            s[ii + nh] += g
    else:
        for ii, i in enumerate(range(0, n, 2)):
            ai = a[ii]; ai1 = a[ii + nh]
            for k in range(nc):
                jf = (n - 1) & (i + nmod + k)
                s[jf] += h2[k] * ai + g2[k] * ai1
    a[:n] = s


def transform(w, a, fwd, zero):
    """dwt.c:gsl_wavelet_transform on a whole list (length a power of 2)."""
    n = len(a)
    if n < 2:
        return a
    sizes = []
    i = n
    while i >= 2:
        sizes.append(i); i >>= 1
    for i in (sizes if fwd else reversed(sizes)):
        step(w, a, i, fwd, zero)
    return a


def col(m, j):
    return [r[j] for r in m]


def setcol(m, j, v):
    for r, x in zip(m, v):
        r[j] = x


def transform2d(w, m, fwd, zero):
    n = len(m)
    if n < 2:
        return m
    rows = lambda: [transform(w, m[i], fwd, zero) for i in range(n)]
    def cols():
        for j in range(n):
            setcol(m, j, transform(w, col(m, j), fwd, zero))
    if fwd:
        rows(); cols()
    else:
        cols(); rows()
    return m


def nstransform2d(w, m, fwd, zero):
    n = len(m)
    if n < 2:
        return m
    sizes = []
    i = n
    while i >= 2:
        sizes.append(i); i >>= 1
    for i in (sizes if fwd else reversed(sizes)):
        def rows():
            for j in range(i):
                step(w, m[j], i, fwd, zero)
        def cols():
            for j in range(i):
                c = col(m, j); step(w, c, i, fwd, zero); setcol(m, j, c)
        if fwd:
            rows(); cols()
        else:
            cols(); rows()
    return m


class Stay(Exception):
    pass


BSPL = {103, 105, 202, 204, 206, 208, 301, 303, 305, 307, 309}


def parse(args):
    """signalprocessing.cc:discrete_wavelet_transform's option loop. args: data then options
    (str = giac string, int, ('sym', name))."""
    kind, k, image, centred = "daubechies", 4, False, False
    for pos, o in enumerate(args[1:], 2):
        if isinstance(o, str):
            if o == "daubechies":
                kind = "daubechies"
            elif o == "haar":
                kind, k = "haar", 2
            elif o == "bspline":
                kind, k = "bspline", 103
            else:
                raise Stay("Unsupported wavelet type")
        elif isinstance(o, int):
            k = o
            if (kind == "daubechies" and (k % 2 or k < 4 or k > 20)) or (kind == "haar" and k != 2) \
                    or (kind == "bspline" and k not in BSPL):
                raise Stay("Wavelet not implemented")
        elif o in (("sym", "center"), ("sym", "centre")):
            centred = True
        elif o == ("sym", "image"):
            image = True
        else:
            raise Stay(f"Invalid argument {pos}")
    return kind, k, image, centred


def is_matrix(d):
    return len(d) > 0 and all(isinstance(r, list) for r in d) and len({len(r) for r in d}) == 1 and len(d[0]) > 0


def pow2(n):
    p = 1
    while p < n:
        p *= 2
    return p


def dwt(args, fwd, conv, zero, num):
    """num: data entry -> number or None (not numeric)."""
    kind, k, image, centred = parse(args)
    data = args[0]
    mat = is_matrix(data)
    n = max(len(data), len(data[0])) if mat else len(data)
    if n == 0:
        return [], None
    warn = None
    if n & (n - 1):
        n = pow2(n)
        warn = f"data size is not a power of 2, padding with zeros up to size {n}" + (f"×{n}" if mat else "")
    w = wavelet(kind, k, centred, conv)
    if not mat:
        a = [zero] * n
        for i, x in enumerate(data):
            v = num(x)
            if v is None:
                raise Stay("Data must be numeric")
            a[i] = v
        return transform(w, a, fwd, zero), warn
    m = [[zero] * n for _ in range(n)]
    for i, r in enumerate(data):
        for j, x in enumerate(r):
            v = num(x)
            if v is None:
                raise Stay("Failed to convert data to GSL matrix")
            m[i][j] = v
    return (nstransform2d if image else transform2d)(w, m, fwd, zero), warn


def num_float(x):
    return float(x) if isinstance(x, (int, float, F)) else None


def num_frac(x):
    return F(x) if isinstance(x, (int, F)) else (F(x) if isinstance(x, float) else None)


class Q2:
    """a + b*sqrt(2), a, b rational (the exact Haar field)."""
    def __init__(s, a, b=F(0)):
        s.a, s.b = F(a), F(b)
    def __add__(s, o): return Q2(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return Q2(s.a - o.a, s.b - o.b)
    def __mul__(s, o): return Q2(s.a * o.a + 2 * s.b * o.b, s.a * o.b + s.b * o.a)
    def __float__(s): return float(s.a) + float(s.b) * 2 ** 0.5
    def __repr__(s): return f"{s.a}+{s.b}r2"


def run(args, fwd=True):
    """(float result, exact-rational result (GSL decimal filters), exact Haar result or None, warning)."""
    fl, warn = dwt(args, fwd, lambda s: float(s), 0.0, num_float)
    fr, _ = dwt(args, fwd, lambda s: F(s), F(0), num_frac)
    kind = parse(args)[0]
    q = None
    exact = all(not isinstance(x, float) for x in flat(args[0]))
    if kind == "haar" and exact:
        q, _ = dwt(args, fwd, lambda s: Q2(0, F(1, 2) if s[0] != "-" else F(-1, 2)), Q2(0), lambda x: Q2(F(x)))
    return fl, fr, q, warn


def flat(d):
    for x in d:
        if isinstance(x, list):
            yield from flat(x)
        else:
            yield x


def giac_parse(s):
    import json
    return json.loads(s)


def close(ours, giac, scale):
    """giac prints 12 significant digits; its noise entries ~1e-16*scale."""
    if isinstance(giac, list):
        return len(ours) == len(giac) and all(close(o, g, scale) for o, g in zip(ours, giac))
    o = float(ours)
    return abs(o - giac) <= 1e-11 * max(abs(giac), 1e-300) + 1e-14 * scale


S = lambda n: ("sym", n)
M44 = [[15, 25, 5, 70], [21, 33, 12, 7], [64, 18, 50, 10], [5, 17, 35, 28]]
M4b = [[1, 2, 3, 4], [4, 5, 6, 7], [1, 0, 0, 1], [2, 2, 1, 1]]
V8 = [1, 2, 3, 4, 5, 6, 7, 8]
ROWS = [  # (args, fwd, giac answer)
    ([[1, 2, 3, 4]], 1, "[5.0,-1.73205080757,2.22044604925e-16,-1.41421356237]"),
    ([[1, 2, 3, 4], 4], 1, "[5.0,-1.73205080757,2.22044604925e-16,-1.41421356237]"),
    ([[1, 2, 3, 4], "daubechies", 4], 1, "[5.0,-1.73205080757,2.22044604925e-16,-1.41421356237]"),
    ([[1, 2, 3, 4], "haar"], 1, "[5.0,-2.0,-0.707106781187,-0.707106781187]"),
    ([[1, 2, 3, 4], "bspline"], 1, "[5.0,2.0,-0.707106781187,-0.707106781187]"),
    ([V8, 6], 1, "[12.7279220614,-1.47938839655,-4.40898590529,2.24670824513,0.0,0.0,-3.79376765089,0.965340526142]"),
    ([[0.0, 0.5, 0.25, 1.0, 0.75, -0.25, -0.5, 0.0]], 1, "[0.618718433538,1.05547557215,-0.121935118395,-0.207311706132,-0.385905771231,0.49163170741,-0.459279326772,-0.176776695297]"),
    ([M44], 1, "[[103.75,-7.23862066798,-2.32472581308,1.97117242248],[-6.5687476869,15.3767603888,-21.3966596673,-0.404199289747],[23.2444350265,12.3853602615,20.7733166849,39.7761931155],[11.7573506422,-31.6760032302,-35.1494327267,-13.9000770737]]"),
    ([[[1, 2], [3, 4]]], 1, "[[5.0,-1.0],[-2.0,-1.66533453694e-16]]"),
    ([[[1, 2], [3, 4]], S("image")], 1, "[[5.0,-1.0],[-2.0,-1.66533453694e-16]]"),
    ([[1, 2, 3]], 1, "[3.0,1.0,1.93185165258,-0.517638090205]"),
    ([[F(1, 2), F(1, 3), F(1, 4), F(1, 5)]], 1, "[0.641666666667,0.136821535725,-0.0268828906508,0.180089359908]"),
    ([[1, 2, 3, 4], S("center")], 1, "[5.0,1.73205080757,-1.41421356237,2.22044604925e-16]"),
    ([[1, 2, 3, 4], "haar", S("center")], 1, "[5.0,0.0,2.12132034356,-0.707106781187]"),
    ([V8, S("center")], 1, "[12.7279220614,2.96712783299,4.83012701892,-1.36602540378,-2.82842712475,2.22044604925e-16,4.4408920985e-16,4.4408920985e-16]"),
    ([V8, "bspline", S("center")], 1, "[12.7279220614,-3.53553390593,2.5,-2.5,4.94974746831,-0.707106781187,-0.707106781187,-0.707106781187]"),
    ([V8, "bspline", 309], 1, "[12.7279220614,-9.22277164993,-1.63037109375,1.63037109375,0.0,2.22044604925e-16,2.22044604925e-16,2.82842712475]"),
    ([V8, "bspline", 202, S("center")], 0, "[2.57842712475,-5.99264068712,4.49264068712,-1.99264068712,5.90685424949,-5.99264068712,6.82106781187,-2.99264068712]"),
    ([[[1, 2, 3], [4, 5, 6]]], 1, "[[5.25,2.02451905284,2.65629602229,-0.181422288142],[5.67163336987,2.26554445662,2.66264202014,0.0109915906878],[-2.06930757797,-0.812832717365,-1.00777222831,0.0322912811497],[-1.11267293737,-0.636056022068,-0.0167468245269,-0.507772228311]]"),
    ([M4b, S("image")], 1, "[[10.0,-1.54903810568,-0.0915063509461,-1.77451905284],[6.19615242271,-2.0,-0.591506350946,0.457531754731],[-2.79246824527,-2.47548094716,-0.341506350946,0.841506350946],[-1.89054445662,-0.841506350946,0.158493649054,-0.658493649054]]"),
    ([[]], 1, "[]"),
    ([[5]], 1, "[5.0]"),
    ([V8, "daubechies", 20], 1, "[12.7279220614,-3.44547055571,-1.4965724436,4.02631614739,-2.68835784988,1.67451012643,-0.916791219661,-0.897788181637]"),
    ([V8, "bspline", 105], 1, "[12.7279220614,5.65685424949,-2.0,-2.0,-0.707106781187,-0.707106781187,-0.707106781187,-0.707106781187]"),
    ([V8, "bspline", 208, S("center")], 1, "[12.7279220614,7.91904352212,3.2001953125,0.7998046875,-2.82842712475,2.22044604925e-16,0.0,-4.4408920985e-16]"),
    ([V8, "bspline", 301], 1, "[12.7279220614,-5.65685424949,-1.0,1.0,0.0,2.22044604925e-16,2.22044604925e-16,2.82842712475]"),
    ([V8, "daubechies", 8, S("center")], 0, "[5.55912833954,-1.50362598936,5.12420118238,-6.45767699054,2.05370954544,-2.47930049791,7.06195080586,-6.52995927066]"),
    ([M4b, S("image")], 0, "[[5.13725952642,0.152803983042,3.91177857926,-2.1887023679],[2.92075317547,2.64382938683,0.896234122635,-3.04575317547],[3.84479128115,-2.71979128115,2.10617061317,-3.24431033399],[1.19527222831,-6.71081668493,1.98774047358,-2.88725952642]]"),
    ([M4b, "haar", S("center")], 0, "[[0.707106781187,-3.41421356237,-0.292893218813,-1.41421356237],[-4.18198051534,8.8033008589,0.0251262658471,3.18198051534],[-3.18198051534,4.97487373415,-1.8033008589,2.18198051534],[1.41421356237,-1.70710678119,-0.585786437627,-0.707106781187]]"),
    ([M4b, "bspline", S("image"), S("center")], 1, "[[10.0,-0.5,1.0,-4.16333634234e-17],[2.0,-0.5,1.5,-0.5],[-1.0,-1.0,-2.0,1.0],[4.5,5.5,1.5,-0.5]]"),
    ([[1.5, 2, 3, 4], "haar"], 1, "[5.25,-1.75,-0.353553390593,-0.707106781187]"),
    ([[1, 2, 3, 4], S("center"), S("image")], 1, "[5.0,1.73205080757,-1.41421356237,2.22044604925e-16]"),
    ([[[1, 2], [3, 4], [5, 6]], "haar"], 1, "[[5.25,5.25,-1.06066017178,0.0],[-0.25,-0.25,-0.353553390593,0.0],[-1.41421356237,-1.41421356237,-1.11022302463e-16,0.0],[3.88908729653,3.88908729653,-0.5,0.0]]"),
    ([[1, 2, 3, 4, 5], "haar"], 1, "[5.3033008589,1.76776695297,-2.0,2.5,-0.707106781187,-0.707106781187,3.53553390593,0.0]"),
]
STAYS = [
    ([[1, 2, 3, 4], 2], "Wavelet not implemented"),
    ([[1, 2, 3, 4], S("haar")], "Invalid argument 2"),
    ([["x", 1, 2, 3]], "Data must be numeric"),
    ([[1, 2, 3, 4], "foo"], "Unsupported wavelet type"),
]


def main():
    bad = 0
    for args, fwd, g in ROWS:
        giac = giac_parse(g)
        fl, fr, q, warn = run(args, bool(fwd))
        scale = max([abs(x) for x in flat(args[0]) if not isinstance(x, str)] + [1]) * 8
        ok_fl = close(fl, giac, scale)
        ok_fr = close(fr, giac, scale)
        ok_q = q is None or close(q, giac, scale)
        # exact rational vs float: same up to float noise
        tag = "ok" if ok_fl and ok_fr and ok_q else "MISMATCH"
        bad += tag != "ok"
        print(tag, ("dwt" if fwd else "idwt"), args[1:], "warn=" + str(warn), "Q2" if q is not None else "")
        if tag != "ok":
            print("  giac", giac); print("  fl  ", fl); print("  fr  ", [float(x) if not isinstance(x, list) else [float(y) for y in x] for x in fr])
        if q is not None:
            print("  exact", q)
    for args, msg in STAYS:
        try:
            dwt(args, True, float, 0.0, num_float)
            print("MISMATCH stay", args); bad += 1
        except Stay as e:
            print("ok stay" if str(e) == msg else "MISMATCH stay", args[1:], e); bad += str(e) != msg
    # round trips (exact rational with GSL's decimal filters is not exactly orthogonal: check ~1e-25)
    for args in ([[1, 2, 3, 4], "haar"], [V8], [M44], [M4b, S("image")], [V8, "bspline", 309, S("center")], [M4b, "daubechies", 20]):
        kind, k, image, c = parse(args)
        f, _ = dwt(args, True, F, F(0), num_frac)
        b, _ = dwt([f] + args[1:], False, F, F(0), num_frac)
        err = max(abs(float(x - F(y))) for x, y in zip(flat(b), flat(args[0])))
        print("roundtrip", args[1:], "maxerr %.1e" % err)
        bad += err > 1e-20
    fq = run([[1, 2, 3, 4], "haar"])[2]
    bq, _ = dwt([fq, "haar"], False, lambda s: Q2(0, F(1, 2) if s[0] != "-" else F(-1, 2)), Q2(0), lambda x: x)
    print("haar exact roundtrip", bq)
    print("mismatches:", bad)


main()
