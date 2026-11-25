"""Prototype of giac's emd (signalprocessing.cc:5131-5217) with GSL's spline
(interpolation/cspline.c cspline_init/cspline_eval, linalg/tridiag.c solve_tridiag,
interpolation/linear.c linear_eval), in plain doubles, for T9.c2 (phase 1.7 scratch)."""
import sys


def bsearch(xa, x, lo, hi):
    while hi > lo + 1:
        i = (hi + lo) // 2
        if xa[i] > x:
            hi = i
        else:
            lo = i
    return lo


def solve_tridiag(diag, off, b):
    n = len(diag)
    alpha = [0.0] * n
    gamma = [0.0] * n
    alpha[0] = diag[0]
    gamma[0] = off[0] / alpha[0] if n > 1 else 0.0
    for i in range(1, n - 1):
        alpha[i] = diag[i] - off[i - 1] * gamma[i - 1]
        gamma[i] = off[i] / alpha[i]
    if n > 1:
        alpha[n - 1] = diag[n - 1] - off[n - 2] * gamma[n - 2]
    z = [0.0] * n
    z[0] = b[0]
    for i in range(1, n):
        z[i] = b[i] - gamma[i - 1] * z[i - 1]
    c = [z[i] / alpha[i] for i in range(n)]
    x = [0.0] * n
    x[n - 1] = c[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = c[i] - gamma[i] * x[i + 1]
    return x


def spline(xa, ya, n):
    """Values at 0..n-1 of the GSL spline through (xa, ya)."""
    size = len(xa)
    out = []
    if size == 2:
        for t in range(n):
            t = float(t)
            dx = xa[1] - xa[0]
            out.append(ya[0] + (t - xa[0]) / dx * (ya[1] - ya[0]))
        return out
    m = size - 1
    sys_size = m - 1
    c = [0.0] * size
    off, diag, g = [], [], []
    for i in range(sys_size):
        hi = xa[i + 1] - xa[i]
        hip1 = xa[i + 2] - xa[i + 1]
        yi = ya[i + 1] - ya[i]
        yip1 = ya[i + 2] - ya[i + 1]
        gi = 1.0 / hi if hi != 0.0 else 0.0
        gip1 = 1.0 / hip1 if hip1 != 0.0 else 0.0
        off.append(hip1)
        diag.append(2.0 * (hip1 + hi))
        g.append(3.0 * (yip1 * gip1 - yi * gi))
    if sys_size == 1:
        c[1] = g[0] / diag[0]
    else:
        c[1:m] = solve_tridiag(diag, off[: sys_size - 1], g)
    for t in range(n):
        t = float(t)
        k = bsearch(xa, t, 0, size - 1)
        xlo, xhi = xa[k], xa[k + 1]
        dx = xhi - xlo
        dy = ya[k + 1] - ya[k]
        delx = t - xlo
        ci, cip1 = c[k], c[k + 1]
        b = (dy / dx) - dx * (cip1 + 2.0 * ci) / 3.0
        d = (cip1 - ci) / (3.0 * dx)
        out.append(ya[k] + delx * (b + delx * (ci + delx * d)))
    return out


def emd(data, imax=2**31 - 1, emin=2, vtol=0.1, max_sifts=10**6):
    n = len(data)
    x = [float(v) for v in data]
    imf = []
    sifts = 0
    while True:
        var = 1.0
        h = list(x)
        isres = False
        notimf = True
        while var > vtol or notimf:
            sifts += 1
            if sifts > max_sifts:
                raise RuntimeError("sift bound")
            maxes, maxy, mins, miny = [0.0], [h[0]], [0.0], [h[0]]
            zc = 0
            for i in range(1, n - 2):
                h1, h2, h3, h4 = h[i - 1], h[i], h[i + 1], h[i + 2]
                if h1 == 0 or h1 * h2 < 0:
                    zc += 1
                d, dl, dr = h3 - h2, h2 - h1, h4 - h3
                if (d == 0 and dl < 0 and dr > 0) or (dl < 0 and d > 0):
                    mins.append(float(i)); miny.append(h[i])
                elif (d == 0 and dl > 0 and dr < 0) or (dl > 0 and d < 0):
                    maxes.append(float(i)); maxy.append(h[i])
            mins.append(float(n - 1)); miny.append(h[n - 1])
            maxes.append(float(n - 1)); maxy.append(h[n - 1])
            ne = len(maxes) + len(mins) - 4
            if ne < emin:
                isres = True
                break
            notimf = abs(zc - ne) > 1
            up = spline(maxes, maxy, n)
            lo = spline(mins, miny, n)
            m = [0.5 * lo[i] + 0.5 * up[i] for i in range(n)]
            var = 0.0
            for i in range(n):
                var += (m[i] * m[i]) / (h[i] * h[i] + 1e-7)
            var /= n
            h = [h[i] - m[i] for i in range(n)]
        if isres:
            return imf, h, sifts
        imf.append(h)
        x = [x[i] - h[i] for i in range(n)]
        if imax == len(imf):
            return imf, x, sifts


def fmt(v):
    return "[" + ",".join("%.12g" % t for t in v) + "]"


if __name__ == "__main__":
    import math
    cases = [
        [1, 2, 1, 3, 1, 2, 1, 3, 1],
        [1, 2, 3, 4],
        [0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0],
        [float("%.12g" % (math.sin(k / 4) + 0.3 * math.sin(3 * k))) for k in range(41)],
    ]
    for c in cases:
        imf, r, s = emd(c)
        print("emd(%s) sifts=%d" % (fmt(c), s))
        print("  imfs=[" + ",".join(fmt(v) for v in imf) + "]")
        print("  residue=" + fmt(r))
