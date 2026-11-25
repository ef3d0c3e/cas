"""T2 prototype (C0009): randvar's parameter fits and discrete forms, as giac
computes them (graphtheory.cc:612-866), with our exact arithmetic.
Prints `input => our object` and checks giac's probe2.out answers.
Origin: phase 1.7, T2 scratch; never shipped."""
import sympy as sp
from mpmath import mp, gamma as G, findroot

mp.dps = 30

def weibull(m, v):
    # giac :759: Gamma(1+2/k)/Gamma(1+1/k) - v/m^2 - 1 = 0 (no square on the
    # denominator: the law's variance is not v), Newton from max(1, 1/v).
    f = lambda k: G(1 + 2 / k) / G(1 + 1 / k) - mp.mpf(v) / mp.mpf(m) ** 2 - 1
    k = findroot(f, max(1, 1 / mp.mpf(v)))
    return k, mp.mpf(m) / G(1 + 1 / k)

def weibull_true(m, v):
    f = lambda k: G(1 + 2 / k) / G(1 + 1 / k) ** 2 - mp.mpf(v) / mp.mpf(m) ** 2 - 1
    k = findroot(f, max(1, 1 / mp.mpf(v)))
    return k, mp.mpf(m) / G(1 + 1 / k)

def fit(law, mean=None, sd=None):
    """The object of randvar(law, mean=, stddev=|variance=) (exact)."""
    if sd is not None and not sd > 0:
        return "stays: Standard deviation must be positive"
    if law == "normal":
        return f"normald({mean if mean is not None else 0},{sd if sd is not None else 1})"
    if law == "uniform":
        r3 = sp.sqrt(3)
        return f"uniformd({sp.simplify(mean - r3 * sd)},{sp.simplify(mean + r3 * sd)})"
    if law == "poisson":
        lam = mean if mean is not None else sd**2
        if mean is not None and sd is not None and sp.simplify(sd**2 - mean) != 0:
            return "stays: Bad Argument Value"
        return f"poisson({lam})" if lam > 0 else "stays: Bad Argument Value"
    if law == "exponential":
        lam = 1 / sp.Rational(mean) if mean is not None else 1 / sd
        return f"exponentiald({lam})"
    if law == "gamma":
        var = sd**2
        return f"gammad({mean**2 / var},{mean / var})"
    if law == "beta":
        var = sd**2
        fac = (var + mean**2 - mean) / var
        return f"betad({-mean * fac},{(mean - 1) * fac})"
    if law == "geometric":
        if mean is not None:
            return f"geometric({1 / sp.Rational(mean)})"
        var = sd**2
        return f"geometric({sp.nsimplify((sp.sqrt(4 * var + 1) - 1) / (2 * var))})"
    if law == "binomial":
        tmp = mean - sd**2
        n = mean**2 / tmp
        p = tmp / mean
        if n != sp.floor(n):
            n = sp.Integer(round(n))
            p = mean / n
        return f"randvar(binomial,{n},{p})"
    if law == "negbinomial":
        n = sp.Integer(round(mean**2 / (sd**2 - mean)))
        return f"negbinomial({n},{mean / sd**2})"
    raise ValueError(law)

def discrete(weights, values=None):
    """randvar(W), randvar(W, V), randvar([[v, w], ...]): the canonical
    object randvar(W[, V]); weights real >= 0, sum > 0."""
    if not weights or any(w < 0 for w in weights) or sum(weights) == 0:
        return "stays: Bad Argument Value"
    if values is not None and len(values) != len(weights):
        return "stays: Items and weights sizes do not match"
    w = "[" + ",".join(map(str, weights)) + "]"
    return f"randvar({w})" if values is None else f"randvar({w},[{','.join(map(str, values))}])"

def density(f, a, b):
    ks = list(range(a, b + 1))
    return discrete([f(sp.Integer(k)) for k in ks], ks)

def pick(weights, r):
    """The exact discrete draw: rational weights over a common denominator
    D, integer masses m_i; r drawn below sum(m_i); the first i whose
    cumulative mass exceeds r."""
    q = [sp.Rational(w) for w in weights]
    d = sp.ilcm(*[x.q for x in q])
    m = [int(x * d) for x in q]
    acc = 0
    for i, mi in enumerate(m):
        acc += mi
        if r < acc:
            return i
    raise AssertionError

def randmultinomial(p):
    return "ok" if sum(p) == 1 and all(x >= 0 for x in p) else "stays: Bad Argument Value"

R = sp.Rational
rows = [
    ("randvar(gammad,mean=12,variance=4)", fit("gamma", 12, sp.sqrt(4)), "gammad(36,3)"),
    ("randvar(binomial,mean=18,stddev=4)", fit("binomial", 18, 4), "randvar(binomial,162,1/9)"),
    ("randvar(uniform,mean=10,stddev=2)", fit("uniform", 10, 2), "uniformd(10 - 2*sqrt(3),2*sqrt(3) + 10)"),
    ("randvar(normal,mean=5,stddev=3/2)", fit("normal", 5, R(3, 2)), "normald(5,3/2)"),
    ("randvar(normal,mean=5,variance=2,stddev=1)", fit("normal", 5, 1), "normald(5,1)"),
    ("randvar(normal,mean=5,stddev=-2)", fit("normal", 5, -2), "stays: Standard deviation must be positive"),
    ("randvar(poisson,mean=3)", fit("poisson", 3), "poisson(3)"),
    ("randvar(exponential,mean=2)", fit("exponential", 2), "exponentiald(1/2)"),
    ("randvar(geometric,mean=4)", fit("geometric", 4), "geometric(1/4)"),
    ("randvar(geometric,variance=12)", fit("geometric", None, sp.sqrt(12)), "geometric(1/4)"),
    ("randvar(betad,mean=2/5,variance=1/25)", fit("beta", R(2, 5), R(1, 5)), "betad(2,3)"),
    ("randvar(binomial,mean=10,stddev=2)", fit("binomial", 10, 2), "randvar(binomial,17,10/17)"),
    ("randvar([3,1,2,5],[alpha,beta,gamma,delta])", discrete([3, 1, 2, 5], ["alpha", "beta", "gamma", "delta"]),
     "randvar([3,1,2,5],[alpha,beta,gamma,delta])"),
    ("randvar([[apple,1/3],[orange,1/4],[pear,1/5],[plum,13/60]])",
     discrete([R(1, 3), R(1, 4), R(1, 5), R(13, 60)], ['"apple"', '"orange"', '"pear"', '"plum"']),
     'randvar([1/3,1/4,1/5,13/60],["apple","orange","pear","plum"])'),
    ("randvar([1,2,3])", discrete([1, 2, 3]), "randvar([1,2,3])"),
    ("randvar(multinomial,[1/2,1/3,1/6],[a,b,c])", discrete([R(1, 2), R(1, 3), R(1, 6)], list("abc")),
     "randvar([1/2,1/3,1/6],[a,b,c])"),
]
bad = 0
for inp, ours, want in rows:
    ok = ours == want
    bad += not ok
    print(f"{inp}  =>  {ours}{'' if ok else '   !! want ' + want}")
d = density(lambda k: 1 - (k / 10) ** 2, -10, 10)
print("randvar(k->1-(k/10)^2,range=-10..10)  =>", d[:60] + "...")
k, lam = weibull(R(25, 2), 1)
print(f"randvar(weibull,mean=12.5,variance=1)  =>  weibulld({mp.nstr(k, 15)},{mp.nstr(lam, 15)})  giac 3.08574940721,13.9803128143")
assert abs(k - mp.mpf("3.08574940721")) < 1e-10 and abs(lam - mp.mpf("13.9803128143")) < 1e-9
kt, lt = weibull_true(R(25, 2), 1)
print(f"  (the moment fit itself, squared denominator: k={mp.nstr(kt, 12)}, lambda={mp.nstr(lt, 12)}; giac's law has variance {mp.nstr(lam**2 * (G(1 + 2 / k) - G(1 + 1 / k) ** 2), 6)})")
print("pick([1/2,1/3,1/6], r) for r in 0..5:", [pick([R(1, 2), R(1, 3), R(1, 6)], r) for r in range(6)])
print("randmultinomial([1/2,1/3,1/6]):", randmultinomial([R(1, 2), R(1, 3), R(1, 6)]),
      "; [1/2,1/2,1/2]:", randmultinomial([R(1, 2)] * 3))
print("rows off:", bad)
