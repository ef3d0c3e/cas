"""giac reduce_Heaviside + linabs tail (signalprocessing.cc) on hand-made step forms."""
import itertools
import sympy as sp
from sympy import Heaviside as H, Abs, sin
x = sp.Symbol('x', real=True)
def size(e): return sp.count_ops(e, visual=False) + len(e.atoms())
def reduce(g):
    g = sp.expand(g); m = {}; rest = 0
    for term in sp.Add.make_args(g):
        hs = [f for f in sp.Mul.make_args(term) if isinstance(f, H)]
        if len(hs) == 1:
            m[hs[0].args[0]] = m.get(hs[0].args[0], 0) + term / hs[0]
        else: rest += term
    ret = 0
    for arg, coef in m.items():
        if sp.diff(arg, x) == 1:
            a = -arg.subs(x, 0); f = sp.together(coef); cst = 1; facs = []
            for fac in sp.Mul.make_args(f):
                (facs if fac.has(x) else [None]) and (facs.append(fac) if fac.has(x) else None)
                if not fac.has(x): cst *= fac
            best = None
            for sel in itertools.product([0, 1], repeat=len(facs)):
                r1 = sp.Mul(*[q for q, s in zip(facs, sel) if s]); un = [q for q, s in zip(facs, sel) if not s]
                if len(un) > 1: continue
                r2 = un[0] if un else sp.Integer(1)
                terms = list(sp.Add.make_args(r2)) if r2.is_Add else [r2]
                for ks in itertools.product([0, 1], repeat=len(terms)):
                    s1 = sp.Add(*[q for q, k in zip(terms, ks) if k]); s2 = sp.Add(*[q for q, k in zip(terms, ks) if not k])
                    rr = sp.expand(s1.subs(x, 2*a - x) + s2)
                    p1 = r1 * s1.subs(x, Abs(x - a) + a); p2 = -r1 * s1.subs(x, 2*a - x)
                    if rr == 0: h1 = h2 = 0
                    else:
                        f2 = rr.subs(x, 2*a - x); d = rr - f2
                        if sp.simplify(d) == 0: continue
                        mu = r1 * rr / d
                        h1 = rr.subs(x, Abs(x - a) + a) * sp.together(mu); h2 = sp.together(-f2 * mu)
                    s = size(sp.together(ret + cst*(h1 + p1) + rest + cst*(h2 + p2)))
                    if best is None or s < best[0]: best = (s, h1, h2, p1, p2)
            if best:
                _, h1, h2, p1, p2 = best; ret += cst*(h1 + p1); rest += cst*(h2 + p2); continue
        ret += coef * (Abs(arg) + arg) / (2*arg)
    return ret + sp.together(rest)
CASES = [  # step form, giac's answer
    (2*x*H(x) - x, Abs(x)),
    (2*sin(x)*H(x) - sin(x), sin(Abs(x))),
    (2*H(x) - 1, Abs(x)/x),
    (H(x), sp.Rational(1, 2) + Abs(x)/(2*x)),
    (x*H(x), x/2 + Abs(x)/2),
    (-x + 2*x*H(x), Abs(x)),
    (1 - x + 2*(x - 1)*H(x - 1), Abs(x - 1)),
    (1 - x + 2*(x - 1)*H(x - 1) + 2*H(x) - 1, Abs(x)/x + Abs(x - 1)),
]
bad = 0
for g, want in CASES:
    got = reduce(g)
    pts = [sp.Rational(k, 3) for k in (-7, -2, 1, 2, 5, 8)]
    if any(abs(sp.N(got.subs(x, p) - want.subs(x, p))) > 1e-12 for p in pts) or any(abs(sp.N(got.subs(x,p) - g.subs(x,p)))>1e-12 for p in pts):
        bad += 1; print('MISMATCH', g, got, want)
    else: print('ok', g, '->', sp.simplify(got))
print('mismatches', bad)
