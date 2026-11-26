"""Prototype of giac's stft and istft (signalprocessing.cc:754-860), in complex
doubles, checked against giac/p6.out (CAS-648, phase 1.8 T6.c2)."""
import cmath, math, os


def dft(x, inverse=False):
    n = len(x)
    s = 1 if inverse else -1
    out = [sum(x[j] * cmath.exp(s * 2j * math.pi * ((j * k) % n) / n) for j in range(n)) for k in range(n)]
    return [v / n for v in out] if inverse else out


# giac's windows of N = m samples, M = m - 1 (signal.rs Window::weight).
def hamming(k, M):
    return 0.54 - 0.46 * math.cos(2 * math.pi * k / M)


def hann(k, M):
    return math.sin(math.pi * k / M) ** 2


def blackman(k, M, a=0.16):
    return (1 - a) / 2 - 0.5 * math.cos(2 * math.pi * k / M) + a / 2 * math.cos(4 * math.pi * k / M)


def bartlett_hann(k, M):
    return 0.62 - 0.48 * abs(k / M - 0.5) - 0.38 * math.cos(2 * math.pi * k / M)


def stft(x, m=128, wf=hamming):
    x = [complex(v).real for v in x]
    if len(x) % m:
        x += [0.0] * (m - len(x) % m)
    N = len(x) // m - 1
    w = [1j * cmath.exp(-1j * math.pi * j / m) for j in range(m)]
    a = [(1 - v) / 2 for v in w]
    b = [(1 + v) / 2 for v in w]
    res = []
    for i in range(N):
        v = [complex(x[i * m + 2 * j], x[i * m + 2 * j + 1]) * wf(j, m - 1) for j in range(m)]
        u = dft(v)
        res.append([a[j] * u[j] + b[j] * u[(m - j) % m].conjugate() for j in range(m)])
    return res


def istft(y):
    N, m = len(y), len(y[0])
    x = [0.0] * ((N + 1) * m)
    w = [1j * cmath.exp(1j * math.pi * j / m) for j in range(m)]
    a = [(1 + v) / 2 for v in w]
    b = [(1 - v) / 2 for v in w]
    for i, r in enumerate(y):
        c = [a[j] * r[j] + b[j] * complex(r[(m - j) % m]).conjugate() for j in range(m)]
        u = dft(c, True)
        for j in range(m):
            x[i * m + 2 * j] += u[j].real
            x[i * m + 2 * j + 1] += u[j].imag
    return x


def parse(text):
    return eval(text.replace("*i", "j"), {"i": 1j})


def leaves(v):
    return [t for e in v for t in leaves(e)] if isinstance(v, list) else [complex(v)]


def close(a, b):
    # giac's 12 digits and its FFT noise (C0006, C0084).
    return abs(a - b) <= 1e-10 * max(1.0, abs(b))


S = [0.3, -0.1, 0.7, 0.2, -0.5, 0.4, 0.9, -0.6, 0.1, 0.8, -0.2, 0.5, -0.7, 0.6, 0.2, -0.4]
E = [0.2, -0.8, 0.4, 0.6, 0.1, -0.2, 0.5, -0.3]
R = list(range(1, 9))

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    giac = {l.split("  =>  ")[0]: l.split("  =>  ")[1].strip() for l in open(os.path.join(here, "giac", "p6.out"))}
    rows = [
        ("stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3],2)", stft(E, 2)),
        ("stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3])", stft(E)),
        ("stft([1,2,3,4,5,6,7,8],4)", stft(R, 4)),
        ("stft([1,2,3,4,5,6,7,8],2)", stft(R, 2)),
        ("stft([1,2,3,4,5,6,7,8,9,10],4)", stft(list(range(1, 11)), 4)),
        ("stft([1,2,3,4],4)", stft([1, 2, 3, 4], 4)),
        ("stft([1,2,3,4,5,6,7,8],8)", stft(R, 8)),
        ("stft([1,2,3,4,5,6,7,8],4,hann_window)", stft(R, 4, hann)),
        ("stft([0.3,-0.1,0.7,0.2,-0.5,0.4,0.9,-0.6,0.1,0.8,-0.2,0.5,-0.7,0.6,0.2,-0.4],4,blackman_window)", stft(S, 4, blackman)),
        ("stft([0.3,-0.1,0.7,0.2,-0.5,0.4,0.9,-0.6,0.1,0.8,-0.2,0.5,-0.7,0.6,0.2,-0.4],8,bartlett_hann_window)", stft(S, 8, bartlett_hann)),
        ("stft([0.3,-0.1,0.7,0.2,-0.5,0.4,0.9,-0.6,0.1,0.8,-0.2,0.5,-0.7,0.6,0.2,-0.4],2)", stft(S, 2)),
        ("stft([1,2+i,3,4],2)", stft([1, 2 + 1j, 3, 4], 2)),
        ("stft([1/2,1/3,1/4,1/5],2)", stft([1 / 2, 1 / 3, 1 / 4, 1 / 5], 2)),
        ("stft([1,2,3,4,5,6,7],2)", stft(list(range(1, 8)), 2)),
        ("stft([1,2,3,4,5,6,7,8],2,hann_window)", stft(R, 2, hann)),
        ("istft(stft([1,2,3,4,5,6,7,8],2))", istft(stft(R, 2))),
        ("istft(stft([1,2,3,4,5,6,7,8],4))", istft(stft(R, 4))),
        ("istft(stft([0.2,-0.8,0.4,0.6,0.1,-0.2,0.5,-0.3],2))", istft(stft(E, 2))),
        ("istft([[1,2],[3,4]])", istft([[1, 2], [3, 4]])),
        ("istft([[1+i,2],[3,4-2*i]])", istft([[1 + 1j, 2], [3, 4 - 2j]])),
        ("istft([[0.5,-1,2,0.25]])", istft([[0.5, -1, 2, 0.25]])),
        ("istft([[1,2,3]])", istft([[1, 2, 3]])),
        ("istft(stft([0.3,-0.1,0.7,0.2,-0.5,0.4,0.9,-0.6,0.1,0.8,-0.2,0.5,-0.7,0.6,0.2,-0.4],4,blackman_window))", istft(stft(S, 4, blackman))),
    ]
    bad = 0
    for inp, ours in rows:
        g, o = leaves(parse(giac[inp])), leaves(ours)
        ok = len(g) == len(o) and all(close(a, b) for a, b in zip(o, g))
        bad += not ok
        print("ok " if ok else "BAD", inp, "" if ok else ["%.12g%+.12gj" % (v.real, v.imag) for v in o])
    print("bad", bad)
