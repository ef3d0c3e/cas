"""Prototype of giac's lin (lin.cc) and period (signalprocessing.cc), 1.8 T4.c2 (CAS-641).
Checked against giac's answers to probe1.cas, probe2.cas (scripts/caseval.py)."""
import sympy as sp
from sympy import I, pi, exp, sin, cos, tan, cot, sec, csc, sinh, cosh, floor, Abs, sqrt

x, y, n, a, t = sp.symbols('x y n a t')


def together(e):
    return sp.together(sp.expand(e))


def compress(v):
    out = []
    for c, g in v:
        for k, (c2, g2) in enumerate(out):
            if sp.simplify(g2 - g) == 0:
                out[k] = (together(c2 + c), g2)
                break
        else:
            out.append((c, g))
    return out


def conv(v, w):
    res = []
    for c, g in v:
        for c2, g2 in w:
            res.append((c * c2, together(g + g2)))
        res = compress(res)
    return res


def lin(e):
    if e.is_Atom:
        return [(e, 0)]
    if e.is_Add:
        v = []
        for arg in e.args:
            v += lin(arg)
        return compress(v)
    if e.is_Mul:
        w = lin(e.args[0])
        for arg in e.args[1:]:
            w = conv(w, lin(arg))
        return w
    if e.func == exp:
        return [(1, unlin(lin(e.args[0])))]
    if e.is_Pow:
        b, k = e.args
        w = lin(b)
        if len(w) == 1:
            c, g = w[0]
            if g == 0 and c.is_Integer:
                return [(1, k * sp.log(c))]
            return [(c ** k, g * k)]
        if k.is_Integer and k >= 0:
            r = [(1, 0)]
            for _ in range(int(k)):
                r = conv(r, w)
            return r
        return [(unlin(w) ** k, 0)]
    f = [unlin(lin(arg)) for arg in e.args]
    if e.func in (sinh, cosh):
        return [(sp.Rational(1, 2), f[0]), (sp.Rational(1 if e.func == cosh else -1, 2), -f[0])]
    return [(e.func(*f, evaluate=False), 0)]


def unlin(v):
    return sp.Add(*[c * exp(g) for c, g in v if c != 0])


def trig2exp(e):
    if e.is_Atom:
        return e
    args = [trig2exp(arg) for arg in e.args]
    u = args[0] if args else None
    p, m = exp(I * u) if u is not None else None, exp(-I * u) if u is not None else None
    if e.func == sin:
        return (p - 1 / p) / (2 * I)
    if e.func == cos:
        return (p + 1 / p) / 2
    if e.func in (tan, cot):
        q = exp(2 * I * u)
        r = (q - 1) / (I * (q + 1))
        return r if e.func == tan else 1 / r
    if e.func == sec:
        return 2 / (p + 1 / p)
    if e.func == csc:
        return 2 * I / (p - 1 / p)
    return e.func(*args)


def linear(u, v):
    av = sp.diff(u, v)
    if av.has(v):
        return None
    return sp.simplify(av), sp.simplify(u - av * v)


class NotPeriodic(Exception):
    pass


def kernels(f, v, acc):
    """Kernels of f in v; returns whether v occurs outside them."""
    if f == v:
        return True
    if not f.has(v):
        return False
    if f.is_Add or f.is_Mul:
        return any([kernels(arg, v, acc) for arg in f.args])
    if f.is_Pow and f.args[1].is_Integer:
        return kernels(f.args[0], v, acc)
    acc.append(f)
    return False


def periodic(g, v, T):
    g = together(g)
    if not g.has(v):
        return True, T
    f = unlin(lin(trig2exp(g)))
    ks = []
    hasx = kernels(f, v, ks)
    if not ks:
        return False, T
    for k in ks:
        U = None
        if k.func == exp and (ab := linear(k.args[0], v)) and ab[0] != 0:
            c = sp.simplify(ab[0] / I)
            if c.has(I):
                return False, T
            U = sp.simplify(Abs(2 * pi / c))
        elif k.func == Abs or len(k.args) != 1 or k.is_Pow:
            return False, T
        else:
            ok, T = periodic(k.args[0], v, T)
            if not ok:
                return False, T
        if U is not None:
            if T is None:
                T = U
            else:
                r = sp.simplify(T / U)
                if not r.is_Rational:
                    return False, T
                T = T * r.q
    return not hasx, T


def period(e, v=x):
    ok, T = periodic(e, v, None)
    if not ok:
        return sp.oo
    return 0 if T is None else T


a = sp.Symbol('a', real=True)
CASES = [
    (period(sin(x * pi / 2) / (2 + cos(x * pi / 3 + pi / 4))), 12),
    (period(sin(x) + sin(sqrt(2) * x)), sp.oo),
    (period(sin(t), t), 2 * pi),
    (period(tan(3 * x)), pi / 3),
    (period(cos(x) ** 2), pi),
    (period(sin(x) ** 2 + cos(x)), 2 * pi),
    (period(sin(2 * x) * cos(3 * x)), 2 * pi),
    (period(x), sp.oo),
    (period(exp(I * x)), 2 * pi),
    (period(sin(x / 3) + tan(x / 2)), 6 * pi),
    (period(sin(a * x)), 2 * pi / Abs(a)),
    (period(sin(x) * cos(x)), pi),
    (period(Abs(sin(x))), sp.oo),
    (unlin(lin(exp(x) ** n)), exp(x * n)),
    (unlin(lin(sp.expand((exp(x) ** 3 + exp(x)) ** 2, deep=False))), exp(6 * x) + 2 * exp(4 * x) + exp(2 * x)),
    (unlin(lin(sp.Mul(exp(x), exp(y), evaluate=False))), exp(x + y)),
    (unlin(lin(sp.Mul(sp.Pow(exp(x), 2, evaluate=False), exp(y), evaluate=False))), exp(2 * x + y)),
    (unlin(lin(sin(x) * exp(x) ** 2)), sin(x) * exp(2 * x)),
    (unlin(lin(x * exp(x) ** 2 + 1)), x * exp(2 * x) + 1),
    (unlin(lin(sp.Pow(exp(x), sp.Rational(1, 2), evaluate=False))), exp(x / 2)),
]
bad = 0
for k, (got, want) in enumerate(CASES):
    if sp.simplify(got - want) != 0 and got != want:
        bad += 1
        print('MISMATCH', k, got, want)
print('cases', len(CASES), 'mismatches', bad)
