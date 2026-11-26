"""Prototype of giac's hht(..., output=matrix) (signalprocessing.cc:5391-5450 inst_freq,
inst_amp, hht; :5674-5791 _hht), in doubles, with proto_analytic's hilbert and
proto_emd's emd, checked against giac/p7.out (CAS-648, phase 1.8 T6.c2)."""
import cmath, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from proto_analytic import hilbert  # noqa: E402
from proto_emd import emd  # noqa: E402

MIDI_NOTE_0 = 8.17579891564


def inst_freq(z, sr):
    a = 1.0 if sr <= 0 else sr / (2 * math.pi)
    out = [0.0]
    for k in range(1, len(z)):
        p = z[k - 1].conjugate() * z[k]
        # -0.0 + 0.0 is +0.0: giac's arg of a negative real is pi.
        out.append(max(0.0, a * math.atan2(p.imag + 0.0, p.real)))
    return out


def hht(data, sr=0, fmin=None, fmax=None, fbins=100, logf=False):
    if isinstance(data[0], list):
        imf = data
    else:
        imf = emd(data)[0]
        if not imf:
            raise ValueError("No IMF found")
    n = len(imf[0])
    hs = [[0.0] * n for _ in range(fbins)]
    f, a = [], []
    for row in imf:
        z = hilbert([complex(v) for v in row])
        f.append(inst_freq(z, sr))
        a.append([abs(v) for v in z])
    if fmin is None:
        fmin = min(min(r) for r in f)
    if fmax is None:
        fmax = max(max(r) for r in f)
    if logf:
        fmin = max(MIDI_NOTE_0, fmin)
    r = math.log(fmax / fmin) if logf else fmax - fmin
    if not r > 0:
        return hs
    for j in range(n):
        for i in range(len(imf)):
            if logf:
                x = math.log(f[i][j] / fmin) if f[i][j] > 0 else -math.inf
            else:
                x = f[i][j] - fmin
            h = math.floor((fbins - 1) * x / r) if math.isfinite(x) else -1
            # giac (the build caseval loads) puts the top of the range one bin
            # low: x == r lands in fbins-2, never fbins-1 (p7.out; probes in
            # p7.cas rows 'top bin'); an exact interior integer stays.
            if x == r and fbins >= 2:
                h = fbins - 2
            if 0 <= h < fbins:
                hs[h][j] += a[i][j]
    return hs


def parse(text):
    return eval(text.replace("*i", "j"), {"i": 1j})


def close(a, b):
    return abs(a - b) <= 1e-10 * max(1.0, abs(b))


P = [0.7, 0.2, -0.5, -0.6, -0.2, 0.3, 0.5, -0.1, -0.3]
Q = [0.1, 0.4, 0.2, -0.3, -0.1, 0.2, 0.6, 0.1, -0.2]
T = [[1, 0, -1, 0, 1, 0, -1, 0], [1, 1, 0, 0, 1, 1, 0, 0]]
L16 = [0.1, 0.9, -0.4, 0.6, -0.8, 0.3, 0.7, -0.9, 0.2, 0.5, -0.6, 0.8, -0.1, 0.4, -0.7, 0.2]

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    giac = {l.split("  =>  ")[0]: l.split("  =>  ")[1].strip() for l in open(os.path.join(here, "giac", "p7.out"))}

    def seq(inp):
        # caseval prints both statements: [L, [k, L[k]] pairs].
        return parse(giac[inp])[0]

    rows = [
        ("L:=hht([[1,0,-1,0,1,0,-1,0],[1,1,0,0,1,1,0,0]],output=matrix);[seq([k,L[k]],k,0,99)]", hht(T)),
        ("L:=hht([[1,0,-1,0]],output=matrix,samplerate=8,frequencies=0..4);[seq([k,L[k]],k,0,99)]",
         hht([[1, 0, -1, 0]], sr=8, fmin=0, fmax=4)),
        ("L:=hht([[0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]],output=matrix);[seq([k,L[k]],k,0,99)]", hht([P])),
        ("L:=hht([[0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3],[0.1,0.4,0.2,-0.3,-0.1,0.2,0.6,0.1,-0.2]],output=matrix,frequencies=0.5..2.5);[seq([k,L[k]],k,0,99)]",
         hht([P, Q], fmin=0.5, fmax=2.5)),
        ("L:=hht([[0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]],output=matrix,samplerate=100,ln);[seq([k,L[k]],k,0,99)]",
         hht([P], sr=100, logf=True)),
        ("L:=hht([[0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]],output=matrix,samplerate=100,ln=true,frequencies=10..40);[seq([k,L[k]],k,0,99)]",
         hht([P], sr=100, logf=True, fmin=10, fmax=40)),
    ]
    bad = 0
    for inp, ours in rows:
        g = seq(inp)
        ok = len(g) == len(ours) and all(
            len(gr) == len(o) and all(close(x, y) for x, y in zip(o, gr)) for gr, o in zip(g, ours))
        bad += not ok
        print("ok " if ok else "BAD", inp[:90])
        if not ok:
            for k, (gr, o) in enumerate(zip(g, ours)):
                if gr != [round(v, 12) for v in o] and any(gr) or any(o):
                    if not all(close(x, y) for x, y in zip(o, gr)):
                        print("  row", k, "giac", gr, "ours", ["%.12g" % v for v in o])
    # list data: emd first; giac prints [size, size(L[0]), rows...]
    inp = "L:=hht([0.1,0.9,-0.4,0.6,-0.8,0.3,0.7,-0.9,0.2,0.5,-0.6,0.8,-0.1,0.4,-0.7,0.2],output=matrix);[size(L),size(L[0]),seq([k,L[k]],k,0,99)]"
    g, got = parse(giac[inp])
    o = hht(L16)
    ok = got[:2] == [100, 16] and len(g) == 100 and all(all(close(x, y) for x, y in zip(a, b)) for a, b in zip(o, g))
    bad += not ok
    print("ok " if ok else "BAD", inp[:90])
    print("bad", bad)
    # Values for the unit tests.
    for name, m in [("T", hht(T)), ("P sr=100 ln", hht([P], sr=100, logf=True))]:
        print(name, [(k, j, "%.17g" % v) for k, r in enumerate(m) for j, v in enumerate(r) if v])
