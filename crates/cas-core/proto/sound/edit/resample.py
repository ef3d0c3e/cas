"""Origin: phase 1.7, T8 scratch. CAS-633 resample prototype: q0/q1 exact (integers), q2-q4 Kaiser sinc in f64.
Writes resample.log beside this file. Usage: python3 resample.py"""
import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, os.environ.get("RS_LOG", "resample.log")), "w")


def log(*a):
    print(*a, file=OUT)
    OUT.flush()


# quality -> (cutoff c, A dB, H zero crossings, L points per unit, passband fraction)
DESIGN = {
    2: (0.90, 100.0, 37, 512, 0.80),  # first set H 29, L 256 (resample.base.log)
    3: (0.95, 120.0, 86, 2048, 0.90),  # 74, 1024
    4: (0.98, 144.0, 312, 8192, 0.96),  # 232, 4096
}
if __name__ == "__main__" and len(sys.argv) > 1:  # overrides: q:H:L
    for spec in sys.argv[1:]:
        q, h, l = map(int, spec.split(":"))
        c, a, _, _, p = DESIGN[q]
        DESIGN[q] = (c, a, h, l, p)


def beta(a):
    return 0.1102 * (a - 8.7)


def i0(x):
    # the series the Rust code uses
    s, t, k, y = 1.0, 1.0, 1, (x / 2) ** 2
    while True:
        t *= y / (k * k)
        s += t
        if t < 1e-17 * s:
            return s
        k += 1


def g_exact(x, q):
    c, a, h, l, _ = DESIGN[q]
    b = beta(a)
    x = np.abs(np.asarray(x, dtype=float))
    w = np.where(x < h, np.i0(b * np.sqrt(np.clip(1 - (x / h) ** 2, 0, None))) / np.i0(b), 0.0)
    return np.sinc(x) * w


TABLES = {}


def table(q):
    if q not in TABLES:
        c, a, h, l, _ = DESIGN[q]
        x = np.arange(h * l + 1) / l
        t = g_exact(x, q)
        t[-1] = 0.0
        TABLES[q] = np.concatenate([t, [0.0]])
    return TABLES[q]


def g_lin(x, q):
    c, a, h, l, _ = DESIGN[q]
    t = table(q)
    p = np.abs(np.asarray(x, dtype=float)) * l
    i = np.floor(p).astype(np.int64)
    inside = i < h * l
    i = np.where(inside, i, 0)
    v = t[i] + (p - i) * (t[i + 1] - t[i])
    return np.where(inside, v, 0.0)


def measure(q):
    c, a, h, l, pf = DESIGN[q]
    t = table(q)[: h * l + 1]
    n = 1
    while n < 40 * h * l:  # >= 20 points per sidelobe (spacing 1/(2H))
        n *= 2
    seq = np.zeros(n)
    seq[: h * l + 1] = t
    seq[n - h * l :] = t[1:][::-1]
    cf = np.fft.rfft(seq).real  # C(f) at f = m*l/n, m = 0..n/2, f in [0, l/2]
    f = np.arange(cf.size) * (l / n)
    gf = cf / l * np.sinc(f / l) ** 2  # transform of the piecewise-linear kernel
    img = np.abs(cf) / l * np.sinc(1 - f / l) ** 2  # first image, f' = l - f
    fp, fs = pf / (2 * c), 1 / (2 * c)
    pas = gf[f <= fp]
    ripple = np.max(np.abs(pas - 1))
    stop = max(np.max(np.abs(gf[f >= fs])), np.max(img))
    # where the stopband really starts: last f below fs... the first f past which |G| < 10^(-A/20)
    lim = 10 ** (-a / 20)
    over = np.nonzero(np.abs(gf) >= lim)[0]
    edge = f[over[-1] + 1] if over.size else 0.0
    # passband edge: first f where |G-1| > lim
    bad = np.nonzero(np.abs(gf - 1) > lim)[0]
    pedge = f[bad[0]] if bad.size else f[-1]
    # table error at midpoints and quarter points against the exact g
    xs = (np.arange(h * l) + 0.5) / l
    xq = (np.arange(h * l) + 0.25) / l
    terr = max(np.max(np.abs(g_lin(xs, q) - g_exact(xs, q))), np.max(np.abs(g_lin(xq, q) - g_exact(xq, q))))
    bound = math.pi ** 2 / (24 * l * l)
    log(f"q{q}: c={c} A={a} beta={beta(a):.6f} H={h} L={l} fft n={n} df={l / n:.3e} (g units)")
    log(f"  passband f<= {fp:.5f} ({pf:.0%} of lower Nyquist): max|G-1| = {ripple:.3e}"
        f" ({20 * math.log10(1 + ripple):.3e} dB), target {lim:.1e}: {'OK' if ripple <= lim else 'MISS'}")
    log(f"  passband edge (|G-1| first > 10^(-A/20)) at f={pedge:.5f} = {pedge * 2 * c:.4f} of lower Nyquist")
    log(f"  stopband f>= {fs:.5f} (lower Nyquist) incl. first image: max|G| = {stop:.3e}"
        f" = {-20 * math.log10(stop):.1f} dB attenuation, target {a} dB: {'OK' if stop <= lim else 'MISS'}")
    log(f"  stopband edge (|G| < 10^(-A/20) beyond) at f={edge:.5f} = {edge * 2 * c:.4f} of lower Nyquist;"
        f" image max {np.max(img):.3e}")
    log(f"  DC gain G(0)-1 = {gf[0] - 1:.3e}")
    log(f"  table error max = {terr:.3e}, bound pi^2/(24L^2) = {bound:.3e}: {'OK' if terr <= bound else 'MISS'}")
    return ripple, stop, terr


# ---- the resampler models ----

def nlen(n, sr, nsr):
    return n * nsr // sr


def encode_exact(depth, v):  # v Fraction
    r = (1 << (depth - 1)) - 1
    v = min(max(v, Fraction(-1)), Fraction(1))
    return math.floor(r * v + Fraction(1, 2))


def encode_f64(depth, s):
    r = float((1 << (depth - 1)) - 1)
    if math.isnan(s):
        return 0
    s = min(max(s, -1.0), 1.0)
    return math.floor(r * s + 0.5)


def q01(codes, depth, sr, nsr, q):
    n = len(codes)
    den = 1 << (depth - 1)
    x = lambda j: Fraction(codes[min(j, n - 1)], den)
    out = []
    for k in range(nlen(n, sr, nsr)):
        t = Fraction(k * sr, nsr)
        if q == 0:
            s = t - 1
            if s <= 0:
                v = x(0)
            else:
                i = math.floor(s)
                f = s - i
                v = x(i) + f * (x(i + 1) - x(i))
        else:
            v = x(max(0, math.floor(t) - 1))
        out.append(encode_exact(depth, v))
    return out


def sinc_resample(codes, depth, sr, nsr, q, ties=None):
    c, a, h, l, _ = DESIGN[q]
    n = len(codes)
    rho = min(1.0, nsr / sr)
    s = c * rho
    w = h / s
    den = float(1 << (depth - 1))
    xs = np.array(codes, dtype=float) / den
    r = float((1 << (depth - 1)) - 1)
    out = []
    for k in range(nlen(n, sr, nsr)):
        t = (k * sr) / nsr
        lo = max(0, math.ceil(t - w))
        hi = min(n - 1, math.floor(t + w))
        j = np.arange(lo, hi + 1)
        u = (k * sr - j * nsr).astype(float) / nsr
        hv = s * g_lin(s * u, q)
        y = float(np.sum(xs[lo : hi + 1] * hv))
        if ties is not None:
            z = r * min(max(y, -1.0), 1.0) + 0.5
            ties.append(abs(z - round(z)))
        out.append(encode_f64(depth, y))
    return out


def resample(codes, depth, sr, nsr, q, ties=None):
    return q01(codes, depth, sr, nsr, q) if q <= 1 else sinc_resample(codes, depth, sr, nsr, q, ties)


def main():
    log("== design measurements (piecewise-linear table, exact transform C(f)/L*sinc^2(f/L)) ==")
    for q in (2, 3, 4):
        measure(q)

    log("\n== lengths floor(len*nsr/sr) ==")
    for n, sr, nsr in [(4, 44100, 88200), (4, 44100, 22050), (4, 44100, 30000), (10, 48000, 44100),
                       (1000, 48000, 44100), (3, 8000, 44100), (1, 44100, 48000), (4, 44100, 1)]:
        log(f"  len {n} {sr}->{nsr}: {nlen(n, sr, nsr)}")

    a16 = [29204, -29204, 14602, 0]
    log("\n== q0/q1 exact rows (A = createwav([0.5,-0.5,0.25,0]) codes) ==")
    rows = [(a16, 16, 44100, 88200, 0), (a16, 16, 44100, 88200, 1), (a16, 16, 44100, 22050, 0),
            (a16, 16, 44100, 22050, 1), (a16, 16, 44100, 30000, 0), (a16, 16, 44100, 30000, 1),
            (a16, 16, 44100, 132300, 0), ([29204], 16, 44100, 132300, 0), ([29204], 16, 44100, 132300, 1),
            ([64, -64, 32, 0], 8, 8000, 20000, 0), ([64, -64, 32, 0], 8, 8000, 20000, 1),
            ([4194304, -4194304, 2097152, 0], 24, 44100, 88200, 0),
            ([5841, 11681, 17522, 23363], 16, 44100, 88200, 0), ([5841, 11681, 17522, 23363], 16, 44100, 66150, 0),
            ([5841, 11681, 17522, 23363], 16, 44100, 66150, 1)]
    for codes, d, sr, nsr, q in rows:
        log(f"  {codes} d{d} {sr}->{nsr} q{q}: {resample(codes, d, sr, nsr, q)}")

    log("\n== sinc rows on A (departure rows; tie distance = |z - round(z)|, z = range*y + 1/2) ==")
    for nsr, q in [(44100, 2), (48000, 4), (22050, 2), (88200, 3)]:
        ties = []
        log(f"  A 44100->{nsr} q{q}: {resample(a16, 16, 44100, nsr, q, ties)} min tie distance {min(ties):.3e}")

    log("\n== DC: constant code 10000 (16 bit) and 3000000 (24 bit), 900 samples, 48000->44100 ==")
    for q in (2, 3, 4):
        for depth, code in [(16, 10000), (24, 3000000)]:
            c, a, h, l, _ = DESIGN[q]
            ties = []
            out = resample([code] * 900, depth, 48000, 44100, q, ties)
            rho = 44100 / 48000
            half = h / (c * rho)
            m = len(out)
            interior = [out[k] for k in range(m) if k * 48000 / 44100 - half >= 0 and k * 48000 / 44100 + half <= 899]
            log(f"  q{q} d{depth}: len {m}, half-width {half:.1f} input samples, interior {len(interior)}"
                f" samples, codes {sorted(set(interior))}, ends {out[0]} {out[-1]}")

    log("\n== tones at 48000 -> 44100, 24 bit, amplitude 1/2, 2048 samples ==")
    top = (1 << 23) - 1
    for freq in (440.0, 23500.0):
        src = [encode_f64(24, 0.5 * math.sin(2 * math.pi * freq * k / 48000)) for k in range(2048)]
        for q in (0, 1, 2, 3, 4):
            out = resample(src, 24, 48000, 44100, q)
            if q <= 1:
                half = 2.0
            else:
                c, a, h, l, _ = DESIGN[q]
                half = h / (c * 44100 / 48000)
            ks = [k for k in range(len(out)) if k * 48000 / 44100 - half >= 0 and k * 48000 / 44100 + half <= 2047]
            if freq < 1000:
                err = max(abs(out[k] - top * 0.5 * math.sin(2 * math.pi * freq * k / 44100)) for k in ks)
                log(f"  440 Hz q{q}: interior {ks[0]}..{ks[-1]} ({len(ks)}), max |code - top*sin/2| = {err:.2f} codes")
            else:
                peak = max(abs(out[k]) for k in ks)
                log(f"  23500 Hz q{q}: interior {ks[0]}..{ks[-1]}, max |code| = {peak} ({peak / top:.3e} of full scale,"
                    f" {20 * math.log10(max(peak, 1) / (top * 0.5)):.1f} dB re the tone)")
    log("\n== work: taps per output (2*ceil(H/s)+1) 48000->44100 and 44100->1 ==")
    for q in (2, 3, 4):
        c, a, h, l, _ = DESIGN[q]
        for sr, nsr in [(48000, 44100), (44100, 4410)]:
            s = c * min(1, nsr / sr)
            log(f"  q{q} {sr}->{nsr}: {2 * math.ceil(h / s) + 1} taps")


if __name__ == "__main__":
    main()
