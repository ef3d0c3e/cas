"""Prototype of giac's dht (signalprocessing.cc:4957-4974) and the list forms of
instfreq/instphase (:5583-5646), in complex doubles, checked against giac/p3.out.
Origin: phase 1.7, T9 scratch."""
import cmath, math, os, re


def dft(x, inverse=False):
    n = len(x)
    s = 1 if inverse else -1
    out = [sum(x[j] * cmath.exp(s * 2j * math.pi * ((j * k) % n) / n) for j in range(n)) for k in range(n)]
    return [v / n for v in out] if inverse else out


def hilbert(x):
    n = len(x)
    u = dft([complex(v.real, 0) for v in x])
    v = [((1 if (k == 0 or k == (n + 1) // 2) else 2) * u[k]) for k in range(n // 2 + 1)]
    v += [0j] * (n - len(v))
    return dft(v, True)


def instfreq(z):
    return [cmath.phase(z[k + 1] * z[k].conjugate()) for k in range(len(z) - 1)]


def instphase(z):
    if not z:
        return []
    p = [cmath.phase(z[0])]
    for w in instfreq(z):
        p.append(p[-1] + w)
    return p


def close(a, b, tol=1e-10):
    return abs(a.real - b.real) <= tol and abs(a.imag - b.imag) <= tol


def parse(text):
    out = []
    for t in text.strip()[1:-1].split(","):
        t = t.replace("*i", "j").replace("i", "1j")
        out.append(complex(eval(t)))
    return out


if __name__ == "__main__":
    L = [0.7, 0.2, -0.5, -0.6, -0.2, 0.3, 0.5, -0.1, -0.3]
    giac = {l.split("  =>  ")[0]: l.split("  =>  ")[1] for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "giac", "p3.out"))}
    rows = [
        ("hilbert([1,0,0,0])", hilbert([1, 0, 0, 0])),
        ("hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3])", hilbert(L)),
        ("hilbert([1,2,3])", hilbert([1, 2, 3])),
        ("hilbert([1,i,0])", hilbert([1, 1j, 0])),
        ("instfreq(hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]))", instfreq(hilbert(L))),
        ("instphase(hilbert([0.7,0.2,-0.5,-0.6,-0.2,0.3,0.5,-0.1,-0.3]))", instphase(hilbert(L))),
        ("instfreq([1,i,-1,-i])", instfreq([1, 1j, -1, -1j])),
        ("instphase([1,i,-1,-i])", instphase([1, 1j, -1, -1j])),
        ("instfreq([1,exp(i/2),exp(i)])", instfreq([1, cmath.exp(0.5j), cmath.exp(1j)])),
        ("instphase([1,exp(i),exp(2*i),exp(3*i)])", instphase([cmath.exp(k * 1j) for k in range(4)])),
        ("instfreq([1,2])", instfreq([1, 2])),
        ("instfreq([1])", instfreq([1])),
        ("instphase([])", instphase([])),
        ("instphase([1.0,-1.0,1.0])", instphase([1.0, -1.0, 1.0])),
    ]
    bad = 0
    for inp, ours in rows:
        g = giac[inp].strip().replace("pi", str(math.pi)).replace("exp", "cmath.exp")
        g = parse(g) if g != "[]" else []
        ok = len(g) == len(ours) and all(close(a, b) for a, b in zip(ours, g))
        bad += not ok
        print("ok " if ok else "BAD", inp, ["%.12g%+.12gj" % (v.real, v.imag) for v in ours] if not ok else "")
    print("bad", bad)
