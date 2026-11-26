#!/usr/bin/python3
# T6.c1 (CAS-630): Python port of giac kmeans (AS 136) and cluster (hclust), optimization.cc.
# Prototype only (R74). Python RNG stands in for giac's _rand. Origin: phase 1.7, T6.c1 scratch.
import os, random, re, sys
from fractions import Fraction

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'giac', '')
NSEEDS = 200

# ---------------- tiny parser for probe inputs ----------------
class Sym:
    def __init__(s, n): s.n = n
    def __repr__(s): return s.n
class Str:
    def __init__(s, v): s.v = v
class Eq:
    def __init__(s, l, r): s.l, s.r = l, r
class Rng:
    def __init__(s, a, b): s.a, s.b = a, b
class Num:
    def __init__(s, v, raw): s.v, s.raw = v, raw
class Lst:
    def __init__(s, items, raw): s.items, s.raw = items, raw

TOK = re.compile(r'\s*(\d+\.(?!\.)\d*|\d+|"[^"]*"|[A-Za-z_]\w*|\.\.|.)')

def tokenize(s):
    out, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: break
        if m.group(1).strip(): out.append((m.group(1), m.start(1), m.end(1)))
        pos = m.end()
    return out

class P:
    def __init__(s, txt):
        s.txt, s.t, s.i = txt, tokenize(txt), 0
    def peek(s): return s.t[s.i][0] if s.i < len(s.t) else None
    def nxt(s):
        s.i += 1; return s.t[s.i-1]
    def expr(s):
        l = s.rng()
        if s.peek() == '=':
            s.nxt(); return Eq(l, s.rng())
        return l
    def rng(s):
        a = s.atom()
        if s.peek() == '..':
            s.nxt(); return Rng(a, s.atom())
        return a
    def atom(s):
        tok, st, en = s.nxt()
        if tok == '[':
            items = []
            if s.peek() != ']':
                items.append(s.expr())
                while s.peek() == ',':
                    s.nxt(); items.append(s.expr())
            _, _, en2 = s.nxt()
            return Lst(items, s.txt[st:en2])
        if tok[0].isdigit():
            v = float(tok) if '.' in tok else Fraction(int(tok))
            if s.peek() == '/':
                s.nxt(); d, _, en = s.nxt(); v = v / Fraction(int(d))
            return Num(v, s.txt[st:en])
        if tok[0] == '"': return Str(tok[1:-1])
        if s.peek() == '(':
            s.nxt(); args = s.args(); s.nxt()
            return (tok, args)
        return Sym(tok)
    def args(s):
        a = [s.expr()]
        while s.peek() == ',':
            s.nxt(); a.append(s.expr())
        return a

def parse_call(line):
    seed = None
    m = re.match(r'srand\((\d+)\);(.*)', line)
    if m: seed, line = int(m.group(1)), m.group(2)
    p = P(line)
    name, args = p.atom()
    args = [a[1][0] if isinstance(a, tuple) and a[0] == 'matrix' else a for a in args]
    return seed, name, args

class GErr(Exception):
    def __init__(s, msg, typ=False):
        s.txt = 'GIAC_ERROR: %s. Error: Bad Argument %s' % (msg, 'Type' if typ else 'Value')

def is_int(g): return isinstance(g, Num) and isinstance(g.v, Fraction) and g.v.denominator == 1
def is_numeric(g):
    if isinstance(g, Num): return True
    if isinstance(g, Lst): return all(is_numeric(x) for x in g.items)
    return False

# ---------------- AS 136 (optimization.cc:9822-10609) ----------------
HUGE = 1.0e30

def optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx):
    for l in range(1, k+1):
        if itran[l] == 1: live[l] = m + 1
    for i in range(1, m+1):
        indx[0] += 1
        l1 = ic1[i]; l2 = ic2[i]; ll = l2
        if 1 < nc[l1]:
            if ncp[l1] != 0:
                de = 0.0
                for j in range(n):
                    df = a[i][j] - c[l1][j]; de = de + df*df
                d[i] = de * an1[l1]
            da = 0.0
            for j in range(n):
                db = a[i][j] - c[l2][j]; da = da + db*db
            r2 = da * an2[l2]
            for l in range(1, k+1):
                if (i < live[l1] or i < live[l2]) and l != l1 and l != ll:
                    rr = r2 / an2[l]
                    dc = 0.0
                    for j in range(n):
                        dd = a[i][j] - c[l][j]; dc = dc + dd*dd
                    if dc < rr:
                        r2 = dc * an2[l]; l2 = l
            if d[i] <= r2:
                ic2[i] = l2
            else:
                indx[0] = 0
                live[l1] = m + i; live[l2] = m + i
                ncp[l1] = i; ncp[l2] = i
                al1 = float(nc[l1]); alw = al1 - 1.0
                al2 = float(nc[l2]); alt = al2 + 1.0
                for j in range(n):
                    c[l1][j] = (c[l1][j]*al1 - a[i][j]) / alw
                    c[l2][j] = (c[l2][j]*al2 + a[i][j]) / alt
                nc[l1] -= 1; nc[l2] += 1
                an2[l1] = alw / al1
                an1[l1] = alw / (alw - 1.0) if 1.0 < alw else HUGE
                an1[l2] = alt / al2
                an2[l2] = alt / (alt + 1.0)
                ic1[i] = l2; ic2[i] = l1
        if indx[0] == m: return
    for l in range(1, k+1):
        itran[l] = 0; live[l] -= m

def qtran(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, indx):
    icoun = 0; istep = 0
    while True:
        for i in range(1, m+1):
            icoun += 1; istep += 1
            l1 = ic1[i]; l2 = ic2[i]
            if 1 < nc[l1]:
                if istep <= ncp[l1]:
                    da = 0.0
                    for j in range(n):
                        db = a[i][j] - c[l1][j]; da = da + db*db
                    d[i] = da * an1[l1]
                if istep < ncp[l1] or istep < ncp[l2]:
                    r2 = d[i] / an2[l2]
                    dd = 0.0
                    for j in range(n):
                        de = a[i][j] - c[l2][j]; dd = dd + de*de
                    if dd < r2:
                        icoun = 0; indx[0] = 0
                        itran[l1] = 1; itran[l2] = 1
                        ncp[l1] = istep + m; ncp[l2] = istep + m
                        al1 = float(nc[l1]); alw = al1 - 1.0
                        al2 = float(nc[l2]); alt = al2 + 1.0
                        for j in range(n):
                            c[l1][j] = (c[l1][j]*al1 - a[i][j]) / alw
                            c[l2][j] = (c[l2][j]*al2 + a[i][j]) / alt
                        nc[l1] -= 1; nc[l2] += 1
                        an2[l1] = alw / al1
                        an1[l1] = alw / (alw - 1.0) if 1.0 < alw else HUGE
                        an1[l2] = alt / al2
                        an2[l2] = alt / (alt + 1.0)
                        ic1[i] = l2; ic2[i] = l1
            if icoun == m: return

def kmns(a0, m, n, c0, k, iter_):
    """a0: m points (0-based); c0: k centers. Returns (ifault, ic1[1..m], nc, wss, c)."""
    if k <= 1 or m <= k: return 3, None, None, None, None
    a = [None] + [list(p) for p in a0]
    c = [None] + [list(p) for p in c0]
    ic1 = [0]*(m+1); ic2 = [0]*(m+1); d = [0.0]*(m+1)
    nc = [0]*(k+1); an1 = [0.0]*(k+1); an2 = [0.0]*(k+1)
    ncp = [0]*(k+1); itran = [0]*(k+1); live = [0]*(k+1)
    for i in range(1, m+1):
        ic1[i] = 1; ic2[i] = 2
        dt = [0.0, 0.0]
        for il in (1, 2):
            for j in range(n):
                da = a[i][j] - c[il][j]; dt[il-1] = dt[il-1] + da*da
        if dt[1] < dt[0]:
            ic1[i] = 2; ic2[i] = 1; dt[0], dt[1] = dt[1], dt[0]
        for l in range(3, k+1):
            db = 0.0
            for j in range(n):
                dc = a[i][j] - c[l][j]; db = db + dc*dc
            if db < dt[1]:
                if dt[0] <= db:
                    dt[1] = db; ic2[i] = l
                else:
                    dt[1] = dt[0]; ic2[i] = ic1[i]; dt[0] = db; ic1[i] = l
    for l in range(1, k+1):
        nc[l] = 0; c[l] = [0.0]*n
    for i in range(1, m+1):
        l = ic1[i]; nc[l] += 1
        for j in range(n): c[l][j] = c[l][j] + a[i][j]
    for l in range(1, k+1):
        if nc[l] == 0: return 1, None, None, None, None
    for l in range(1, k+1):
        aa = float(nc[l])
        for j in range(n): c[l][j] = c[l][j] / aa
        an2[l] = aa / (aa + 1.0)
        an1[l] = aa / (aa - 1.0) if 1.0 < aa else HUGE
        itran[l] = 1; ncp[l] = -1
    indx = [0]; ifault = 2
    for ij in range(1, iter_+1):
        optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx)
        if indx[0] == m: ifault = 0; break
        qtran(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, indx)
        if k == 2: ifault = 0; break
        for l in range(1, k+1): ncp[l] = 0
    wss = [0.0]*(k+1)
    for l in range(1, k+1): c[l] = [0.0]*n
    for i in range(1, m+1):
        ii = ic1[i]
        for j in range(n): c[ii][j] = c[ii][j] + a[i][j]
    for j in range(n):
        for l in range(1, k+1): c[l][j] = c[l][j] / float(nc[l])
        for i in range(1, m+1):
            ii = ic1[i]; da = a[i][j] - c[ii][j]; wss[ii] = wss[ii] + da*da
    return ifault, ic1[1:], nc[1:], wss[1:], c[1:]

def init_cluster_centers(x, ind, m, k, rng):
    """optimization.cc:10622-10646. _rand(n) -> rng.randrange(n); _rand(k,list) -> rng.sample."""
    part = {}
    for i in range(m): part.setdefault(1 if ind is None else ind[i], []).append(i)
    ch = []; maxsize = 0; lc = None
    for key in sorted(part):
        lst = part[key]
        if len(lst) > maxsize: maxsize = len(lst); lc = key
        pos = rng.randrange(len(lst))
        ch.append(lst[pos]); del lst[pos]
    if len(ch) < k: ch += rng.sample(part[lc], k - len(ch))
    assert len(ch) == k
    return [list(x[i]) for i in ch]

def kmeans_core(x, k, min_k, max_k, maxiter, rng):
    """optimization.cc:10654-10763, meth=_HARTIGAN_CRITERION (0) only. Returns (code, k, res)."""
    m, n = len(x), len(x[0])
    if k > 0:
        tries = 0
        while True:
            d = init_cluster_centers(x, None, m, k, rng)
            ifault, b, e, dev, _ = kmns(x, m, n, d, k, maxiter)
            if not (ifault == 1 and tries < 10): break
            tries += 1
        if ifault != 0: return ifault, k, None
        return 0, k, [v-1 for v in b]
    W_prev = 0.0; b_best = [i+1 for i in range(m)]; rough = False; b = None
    k = min_k
    while k <= max_k:
        tries = 0
        while True:
            d = init_cluster_centers(x, None if (k == min_k or tries > 0) else b, m, k, rng)
            ifault, b2, e, dev, _ = kmns(x, m, n, d, k, maxiter)
            if b2 is not None: b = b2
            if not (ifault == 1 and tries < 10): break
            tries += 1
        if ifault != 0: return ifault, k, None
        if not rough:
            W = 0.0
            for i in range(k): W += dev[i]
            if W > 0 and W_prev > 0:
                score = (W_prev/W - 1)*(m - k)
                if score <= 10:
                    rough = True; k -= 1; break
            W_prev = W
            b_best = list(b)
            if W == 0: break
        k += 1
    if k > max_k: k = max_k
    return 0, k, [v-1 for v in (b_best if rough else b)]

# ---------------- cluster_crit (optimization.cc:10784-11091), T5.c3 ----------------
# Values: Fraction (giac exact), float (giac double), None (undef), complex (a non-real ln),
# INF (giac's unsigned infinity, trace-WiB on a singular W).
import cmath
INF = float('inf')

def real(v):
    return v is not None and not isinstance(v, complex) and v == v

EPS = 1e-12  # giac's epsilon(contextptr)

def doublerref_det(a):
    """in_doublerref's determinant (vecteur.cc:9752-9871, rref_or_det_or_lu=1): the column's largest
    pivot, entries under epspivot zeroed (1e-12 first, then eps*|pivot|)."""
    N = [[float(v) for v in r] for r in a]; n = len(N)
    idet = 1.0; epspivot = max(EPS, 1e-13); l = c = 0
    while l < n and c < n:
        pivot = N[l][c]
        if abs(pivot) < epspivot: pivot = N[l][c] = 0.0
        if l == n - 1: return idet * pivot
        pl = l
        for lt in range(l + 1, n):
            t = N[lt][c]
            if abs(t) < epspivot: t = N[lt][c] = 0.0
            if abs(t) > abs(pivot): pivot = t; pl = lt
        if pivot == 0: return 0.0
        epspivot = abs(EPS * pivot)
        if pl != l: N[l], N[pl] = N[pl], N[l]; idet = -idet
        N[l] = [v / pivot for v in N[l]]
        idet = idet * pivot
        for lt in range(l + 1, n):
            f = -N[lt][c]
            for j in range(c + 1, n): N[lt][j] += f * N[l][j]
        c += 1; l += 1
    return idet

def det(a):
    """_det: exact on an exact matrix, else doublerref on its doubles (mrref, vecteur.cc:7414-7450)."""
    if all(isinstance(v, (int, Fraction)) for r in a for v in r):
        n = len(a); m = [[Fraction(v) for v in r] for r in a]; d = Fraction(1)
        for c in range(n):
            piv = next((r for r in range(c, n) if m[r][c] != 0), None)
            if piv is None: return Fraction(0)
            if piv != c: m[c], m[piv] = m[piv], m[c]; d = -d
            d *= m[c][c]
            for r in range(c + 1, n):
                f = m[r][c] / m[c][c]
                m[r] = [x - f*y for x, y in zip(m[r], m[c])]
        return d
    return doublerref_det(a)

def inverse(a):
    """A double matrix's inverse by Gauss-Jordan; None when singular (giac: infinity)."""
    if doublerref_det(a) == 0: return None
    n = len(a); m = [[float(v) for v in r] + [1.0 if i == j else 0.0 for j in range(n)] for i, r in enumerate(a)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[piv] = m[piv], m[c]
        pv = m[c][c]; m[c] = [v / pv for v in m[c]]
        for r in range(n):
            if r != c:
                f = m[r][c]; m[r] = [x - f*y for x, y in zip(m[r], m[c])]
    return [r[n:] for r in m]

def blocks(a):
    """The index sets of is_blockmatrix's recursive split (vecteur.cc:11907): the connected
    components of i ~ j when a[i][j] or a[j][i] is non-zero, each in increasing order."""
    n = len(a); seen = [False]*n; out = []
    for st in range(n):
        if seen[st]: continue
        comp = [st]; seen[st] = True; q = [st]
        while q:
            i = q.pop()
            for j in range(n):
                if not seen[j] and (a[i][j] != 0 or a[j][i] != 0): seen[j] = True; comp.append(j); q.append(j)
        out.append(sorted(comp))
    return out

def fdiv(a, b):
    """A double division as giac's: by 0.0 an infinity (giac prints infinity), 0/0 undef."""
    if b == 0: return None if a == 0 else (INF if a > 0 else -INF)
    return a / b

def gdiv(a, b):
    """giac's a/b with b a double: b == 1.0 leaves a (exact) as is."""
    return a if b == 1.0 else a / b

class ClusterCrit:
    def __init__(s, pts):
        s.pts = pts; s.N = len(pts); s.p = len(pts[0])
        s.mu0 = [sum((q[j] for q in pts), 0) / s.N for j in range(s.p)]
        s.T = [[1 if i == j else 0 for j in range(s.p)] for i in range(s.p)]  # _T=midn(p)
        s.ET = 0
        for q in pts: s.ET += float(sum(((a - b)**2 for a, b in zip(q, s.mu0)), 0))
        s.vmap = {}; s.optvals = {}
    def init(s, K, b, centers, wss):
        s.K = K; s.ci = [v - 1 for v in b]
        s.n = [0]*K
        for c in s.ci: s.n[c] += 1
        s.mu = [list(c) for c in centers]; s.wgss = list(wss)
        s.WGSS = 0
        for w in wss: s.WGSS += w
        s.have_T = s.have_WG = False; s.dcache = {}; s.dbmm = None
    def Tm(s):
        if not s.have_T:
            p, N = s.p, s.N
            for i in range(p):
                for j in range(i, p):
                    v = 0
                    for l in range(N): v += (s.pts[l][i] - s.mu0[i]) * (s.pts[l][j] - s.mu0[j])
                    s.T[i][j] = v
                    if i != j: s.T[j][i] = v
            s.have_T = True
        return s.T
    def TSS(s, j):
        v = 0
        for l in range(s.N): v += (s.pts[l][j] - s.mu0[j])**2
        return v
    def WG(s, k=-1):
        if not s.have_WG:
            p = s.p; s.WGk = []; s.WGsum = None
            for c in range(s.K):
                W = [[0]*p for _ in range(p)]
                for i in range(p):
                    for j in range(i, p):
                        v = 0
                        for l in range(s.N):
                            if s.ci[l] == c: v += (s.pts[l][i] - s.mu[c][i]) * (s.pts[l][j] - s.mu[c][j])
                        W[i][j] = v
                        if i != j: W[j][i] = s.T[i][j]  # giac's quirk: _T's entry (optimization.cc:10880)
                s.WGk.append(W)
                s.WGsum = W if s.WGsum is None else [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(s.WGsum, W)]
            s.have_WG = True
        return s.WGsum if k < 0 else s.WGk[k]
    def BGSS(s, j=-1):
        if j < 0:
            t = 0
            for jj in range(s.p): t += s.BGSS(jj)
            return t
        v = 0
        for c in range(s.K): v += s.n[c] * (s.mu[c][j] - s.mu0[j])**2
        return v
    def BG(s):
        p = s.p; B = [[0]*p for _ in range(p)]
        for i in range(p):
            for j in range(i, p):
                v = 0
                for c in range(s.K): v += s.n[c] * (s.mu[c][i] - s.mu0[i]) * (s.mu[c][j] - s.mu0[j])
                B[i][j] = v; B[j][i] = v
        return B
    def Delta(s, k, l):
        key = (min(k, l), max(k, l))
        if key not in s.dcache:
            s.dcache[key] = float(sum((a - b)**2 for a, b in zip(s.mu[k], s.mu[l])))**0.5
        return s.dcache[key]
    def DB(s, d):
        if s.dbmm is None:
            lo, hi = INF, -INF
            for k in range(s.K):
                for l in range(k+1, s.K):
                    D = s.Delta(k, l); lo = min(lo, D); hi = max(hi, D)
            s.dbmm = (lo, hi)
        return s.dbmm[0] if d < 0 else s.dbmm[1]
    # the indices
    def banfeld_raftery(s):
        r = 0
        for k in range(s.K):
            if s.wgss[k] == 0: return None
            r += s.n[k] * cmath.log(s.wgss[k] / s.n[k]).real
        return float(r)
    def calinski_harabasz(s):
        if s.K <= 1: return None
        return fdiv((s.N - s.K) * s.BGSS(), (s.K - 1) * s.WGSS)
    def davies_bouldin(s):
        M = [-INF]*s.K
        for k in range(s.K):
            for l in range(s.K):
                if k != l:
                    v = fdiv(s.wgss[k]/s.n[k] + s.wgss[l]/s.n[l], s.Delta(k, l))
                    if v is None: return None
                    M[k] = max(v, M[k])
        return sum(M) / s.K
    def pbm(s):
        v = fdiv(s.ET * s.DB(1), s.K * s.WGSS)
        return None if v is None else v**2
    def ratkowsky_lance(s):
        r = 0
        for j in range(s.p):
            t = s.TSS(j)
            if t == 0: return None
            r += s.BGSS(j) / t
        return float(r / s.p)**0.5
    def ray_turi(s):
        return fdiv(s.WGSS, s.N * s.DB(-1)**2)
    def scott_symons(s):
        r = 0
        for k in range(s.K):
            d = det(s.WG(k))
            if d == 0: return None
            x = d / s.n[k]
            if x < 0: return None
            r += s.n[k] * cmath.log(x).real
        return float(r)
    def det_ratio(s):
        den = det(s.WG())
        if den == 0: return None
        return gdiv(det(s.Tm()), den)
    def log_det_ratio(s):
        num = det(s.Tm()); den = det(s.WG())
        if num == 0 or den == 0: return None
        r = gdiv(num, den)
        if isinstance(r, Fraction): return 0 if r == 1 else complex(0, 0)  # ln of a rational: symbolic, not real
        return s.N * cmath.log(r) if r < 0 else s.N * cmath.log(r).real
    def ksq_detW(s):
        return s.K * s.K * det(s.WG())
    def log_ss_ratio(s):
        num, den = s.BGSS(), s.WGSS
        if num == 0 or den == 0: return None
        return cmath.log(num / den).real
    def trace_W(s): return s.WGSS
    def trace_WiB(s):
        """minv (vecteur.cc:11996) splits W into is_blockmatrix's blocks (:11907); a singular block's
        LU inverse has infinities: its trace undef (NaN) when B is 0 on it, else infinity (probes)."""
        W, B = s.WG(), s.BG(); t = 0.0
        for blk in blocks(W):
            Wb = [[W[i][j] for j in blk] for i in blk]
            inv = inverse(Wb)
            if inv is None:
                if all(B[i][j] == 0 for i in blk for j in blk): return None
                t = INF; continue
            t += sum(inv[a][b] * B[blk[b]][blk[a]] for a in range(len(blk)) for b in range(len(blk)))
        return float(t)
    def ball_hall(s):
        r = 0
        for k in range(s.K): r += s.wgss[k] / s.n[k]
        return r / s.K
    def update(s, pos, val, d):
        v = s.vmap.setdefault(pos, []); v.append(val)
        if val is None: return
        V = val
        elbow = d in (-2, 2)
        if elbow:
            if len(v) < 3: return
            if any(x is None for x in v[-3:]): return
            V = v[-1] - 2*v[-2] + v[-3]
        if not real(V): return
        Vd = float(V)
        o = s.optvals.get(pos)
        if o is None or (d < 0 and o[1] > Vd) or (d > 0 and o[1] < Vd):
            s.optvals[pos] = (s.K - 1 if elbow else s.K, Vd)
    def compute_indices(s, crit):
        pos = 0
        for bit, f, d in ((1, s.banfeld_raftery, -1), (2, s.davies_bouldin, -1), (4, s.ray_turi, -1),
                          (8, s.scott_symons, -1), (16, s.calinski_harabasz, 1), (32, s.pbm, 1),
                          (64, s.ratkowsky_lance, 1), (128, s.det_ratio, -2), (256, s.log_det_ratio, -2),
                          (512, s.log_ss_ratio, -2), (2048, s.ball_hall, 2), (1024, s.ksq_detW, 2),
                          (4096, s.trace_W, 2), (8192, s.trace_WiB, 2)):
            if crit & bit: s.update(pos, f(), d); pos += 1
    def ncopt(s):
        sugg = {}
        for k_, _ in s.optvals.values(): sugg[k_] = sugg.get(k_, 0) + 1
        best, num = 0, 0
        for k_ in sorted(sugg):
            if best == 0 or sugg[k_] > num: best, num = k_, sugg[k_]
        return best

KMEANS_INDEX = {'banfeld-raftery': 1, 'davies-bouldin': 2, 'ray-turi': 4, 'scott-symons': 8,
                'calinski-harabasz': 16, 'pbm': 32, 'ratkowsky-lance': 64, 'det': 128, 'log-det': 128,
                'ksq-detW': 1024, 'log-ss': 512, 'ball-hall': 2048, 'trace-W': 4096, 'trace-WiB': 8192,
                'all': 16383}

def kmeans_index(x, pts, meth, min_k, max_k, maxiter, rng):
    """kmeans' index path (optimization.cc:10682-10750, meth > 0): per-k runs seeded by the
    previous labels, Hartigan's rough optimum caps max_k at 2k, labels of cluster_crit's ncopt.
    Returns (code, k, res, kvals, values of the first index)."""
    m, n = len(x), len(x[0])
    W_prev = 0.0; rough = False; b = None
    ccr = ClusterCrit(pts); b_map = {}; kvals = []
    k = min_k
    while k <= max_k:
        tries = 0
        while True:
            d = init_cluster_centers(x, None if (k == min_k or tries > 0) else b, m, k, rng)
            ifault, b2, e, dev, c = kmns(x, m, n, d, k, maxiter)
            if b2 is not None: b = b2
            if not (ifault == 1 and tries < 10): break
            tries += 1
        if ifault != 0: return ifault, k, None, None, None
        if not rough:
            W = 0.0
            for i in range(k): W += dev[i]
            if W > 0 and W_prev > 0:
                if (W_prev/W - 1)*(m - k) <= 10:
                    rough = True; max_k = min(max_k, 2*k)
            W_prev = W
        ccr.init(k, b, c, dev); ccr.compute_indices(meth)
        kvals.append(k); b_map[k] = [v-1 for v in b]
        k += 1
    k = ccr.ncopt()
    return 0, k, b_map.get(k, []), kvals, ccr.vmap.get(0, [])

# ---------------- hclust (optimization.cc:11102-11443) ----------------
SINGLE, AVERAGE, COMPLETE, WEIGHTED, WARD = 1, 2, 3, 4, 5
DBL_MAX = sys.float_info.max

def longueur2(p, q):
    """at_longueur2 on two data points (plot.cc:6177): exact sum of squared differences, then evalf.
    Two scalars: distance2pp -> (e-f)^2 (+ im^2 = 0). Returns None when not a real number."""
    if isinstance(p, Lst) and isinstance(q, Lst):
        if not (is_numeric(p) and is_numeric(q)) or len(p.items) != len(q.items): return None
        res = 0
        for u, v in zip(p.items, q.items):
            t = u.v - v.v; res = res + t*t
        return float(res)
    if isinstance(p, Num) and isinstance(q, Num):
        t = p.v - q.v; return float(t*t + 0)
    return None

def levenshtein(p, q):
    """at_levenshtein on two strings (optimization.cc:11519-11536): std::string::at, so UTF-8 bytes."""
    if not (isinstance(p, Str) and isinstance(q, Str)): return None
    s, t = p.v.encode(), q.v.encode()
    if len(s) > len(t): s, t = t, s
    v0 = list(range(len(t) + 1))
    for i in range(len(s)):
        v1 = [i + 1] + [0]*len(t)
        for j in range(len(t)):
            v1[j+1] = min(v0[j+1] + 1, v1[j] + 1, v0[j] + (0 if s[i] == t[j] else 1))
        v0 = v1
    return float(v0[len(t)])

class Abort(Exception):
    """giac's assertion at optimization.cc:11325 (libgiac aborts)."""

class HClust:
    def __init__(s, pts, distf=longueur2):
        s.pts = pts; s.N = len(pts); s.S = list(range(s.N)); s.cache = {}; s.distf = distf
    def dist(s, i, j, cache):
        if i == j: return 0.0
        if cache:
            key = (min(i, j), max(i, j))
            if key not in s.cache:
                v = s.distf(s.pts[i], s.pts[j])
                if v is None: return -1.0
                s.cache[key] = v
            return s.cache[key]
        v = s.distf(s.pts[i], s.pts[j])
        return -1.0 if v is None else v
    def formula(s, a, b, x, na, nb, nx, dab, meth):
        dax = s.dist(a, x, True); dbx = s.dist(b, x, True)
        if dax < 0 or dbx < 0 or dab < 0 or dax + dbx < dab: return -1.0
        if meth == SINGLE: return min(dax, dbx)
        if meth == COMPLETE: return max(dax, dbx)
        if meth == AVERAGE: return (na*dax + nb*dbx)/float(na + nb)
        if meth == WEIGHTED: return (dax + dbx)/2.0
        if meth == WARD: return ((na+nx)*dax + (nb+nx)*dbx - nx*dab)/(na+nb+nx)
        raise AssertionError
    def k_low(s, dg):
        dprev = 0.0; max_step = -1.0; k = 0
        for i, (d, _) in enumerate(dg):
            if d - dprev > max_step: max_step = d - dprev; k = s.N - i
            dprev = d
        return k
    def mst_linkage(s, cache_d):
        N = s.N; S = s.S; c = S[0]; dg = []
        D = [0.0]*(2*N); prev_first = True
        for v in S:
            if v != c: D[v] = DBL_MAX
        for _ in range(1, N):
            minD = DBL_MAX; n = -1; off = 0 if prev_first else N
            S.remove(c)
            for v in S:
                Dval = min(D[v+off], s.dist(v, c, cache_d))
                if Dval < 0: return None
                if Dval < minD: n = v; minD = Dval
                D[v+N-off] = Dval
            assert n >= 0
            dg.append((D[n+N-off], (c, n)))
            c = n; prev_first = not prev_first
        return dg
    def nn_chain_linkage(s, meth):
        S = s.S; dg = []; chain = []; size = {v: 1.0 for v in S}
        a = b = c = None
        while len(S) > 1:
            if len(chain) <= 3:
                a = S[0]; chain = [a]; b = S[1]
            else:
                a = chain[-4]; b = chain[-3]; chain = chain[:-3]
            while True:
                mind = DBL_MAX
                for x in S:
                    if x != a:
                        d = s.dist(x, a, True)
                        if d < mind or (x == b and d == mind): c = x; mind = d
                b = a; a = c; chain.append(a)
                if len(chain) >= 3 and a == chain[-3]: break
            dab = s.dist(a, b, True)
            if dab < 0: return None
            dg.append((dab, (a, b)))
            S.remove(b)
            for x in S:
                if x == a: continue
                f = s.formula(a, b, x, int(size[a]), int(size[b]), int(size[x]), dab, meth)
                if f < 0: return None
                s.cache[(min(x, a), max(x, a))] = f
            size[a] += size[b]
        return dg
    def linkage(s, meth, K, min_k, max_k, depart=False):
        """optimization.cc:11262-11443. K>0 fixes k; K=-indf<0 runs the index path (silhouette only computes:
        mcclain-rao and dunn push k values and nothing else). depart: ours clamps max_k to N-1 on the index
        path, where giac aborts (Abort). Returns (k, ind, dg, kvals, vals) or (-1, None, None, None, None)."""
        dg = s.mst_linkage(K < 0) if meth == SINGLE else s.nn_chain_linkage(meth)
        if dg is None: return -1, None, None, None, None
        dg.sort()
        N = s.N
        par = list(range(N))
        def find(i):
            while par[i] != i: i = par[i]
            return i
        indf = -K if K < 0 else 0; sil = indf & 1
        if indf > 0:
            s.cache = {}  # true distances again (nn_chain cached Lance-Williams values)
            if depart: max_k = min(max_k, N - 1)
        k = N; best_k = 0
        k_lb = max(min_k, s.k_low(dg) if K <= 0 else 0)
        ind = None; kvals = []; vals = []; opt = None
        slh_a = [[0, 0.0] for _ in range(N)]; slh_b = [dict() for _ in range(N)]
        for d, (a, b) in dg:
            if k == (K if K > 0 else k_lb): break
            k -= 1
            if indf > 0 and k <= max_k:
                if k < max_k:
                    ca, cb = find(a), find(b)
                    if sil:
                        for i in range(N):
                            ci = find(i)
                            if ci == ca or ci == cb:
                                o = cb if ci == ca else ca
                                if o not in slh_b[i]: raise Abort()
                                c, sm = slh_b[i].pop(o)
                                slh_a[i][0] += c; slh_a[i][1] += sm
                    par[cb] = ca
                    if sil:
                        for i in range(N):
                            if find(i) != ca:
                                bm = slh_b[i]
                                if ca not in bm or cb not in bm: raise Abort()
                                x = bm.pop(cb); bm[ca][0] += x[0]; bm[ca][1] += x[1]
                else:
                    ra, rb = find(a), find(b)
                    if ra != rb: par[rb] = ra
                    if sil:
                        slh_a = [[0, 0.0] for _ in range(N)]
                        for i in range(N):
                            ci = find(i)
                            for j in range(i + 1, N):
                                cj = find(j); dd = s.dist(i, j, True)
                                if ci == cj:
                                    slh_a[i][0] += 1; slh_a[j][0] += 1; slh_a[i][1] += dd; slh_a[j][1] += dd
                                else:
                                    p = slh_b[i].setdefault(cj, [0, 0.0]); q = slh_b[j].setdefault(ci, [0, 0.0])
                                    p[0] += 1; q[0] += 1; p[1] += dd; q[1] += dd
                updated = False
                if sil:
                    slh_s = [[0, 0.0] for _ in range(N)]
                    for i in range(N):
                        ci = find(i); minb = DBL_MAX; ap = slh_a[i]
                        ai = 0.0 if ap[0] == 0 else ap[1]/ap[0]
                        for key in sorted(slh_b[i]):
                            bi = slh_b[i][key][1]/slh_b[i][key][0]
                            if bi < minb: minb = bi
                        den = ai if minb < ai else minb
                        si = 0.0 if ap[0] == 0 else ((minb - ai)/den if den != 0 else float('nan'))
                        slh_s[ci][0] += 1; slh_s[ci][1] += si
                    indval = sys.float_info.min  # DBL_MIN
                    for c, sm in slh_s:
                        if c == 0: continue
                        if sm/c > indval: indval = sm/c
                    vals.append(indval)
                    if opt is None or opt[1] < indval: opt = (k, indval); updated = True
                kvals.append(k)
                if updated and opt[0] != best_k:
                    ind = [find(i) for i in range(N)]; best_k = opt[0]
            else:
                ra, rb = find(a), find(b)
                if ra != rb: par[rb] = ra
        if ind is None: ind = [find(i) for i in range(N)]
        cmap = {}; lab = []
        for r in ind:
            if r not in cmap: cmap[r] = len(cmap)
            lab.append(cmap[r])
        return max(k, best_k), lab, dg, kvals, vals

def name2index(n):  # hclust::name2index (:11481-11487)
    return {'silhouette': 1, 'mcclain-rao': 2, 'dunn': 4, 'all': 7}.get(n, -1)

# ---------------- _kmeans / _cluster (optimization.cc:11571-11920) ----------------
def raw(g):
    if isinstance(g, Str): return '"%s"' % g.v
    return g.raw if isinstance(g, (Num, Lst)) else str(g)

def fmt_list(v): return '[' + ','.join(str(x) for x in v) + ']'

def fmt_val(v):
    """An index value as giac prints it: 12 digits, a/b, "not defined" (undef), infinity."""
    if v is None: return '"not defined"'
    if v in (INF, -INF): return 'infinity' if v > 0 else '-infinity'
    if isinstance(v, Fraction): return str(v)
    if isinstance(v, complex): return '(%r)' % v
    return '%.12g' % v

def clustering_result(data, p, out_what, k, res):
    if out_what == 2: return fmt_list(res)
    if out_what == 4: return str(k)
    if out_what == 3:
        # cyclic_colormap(c) | _POINT_POINT|_POINT_WIDTH_2 observed as 235405312 + c for c in 0..2
        pts = ','.join('point(%s,display=%d)' % (raw(data[i]), 235405312 + res[i]) for i in range(len(data)))
        return '[title="Number of clusters: %d",axes=0,%s]' % (k, pts)
    if out_what == 1:
        cl = [[] for _ in range(k)]
        for i, r in enumerate(res): cl[r].append(raw(data[i]))
        return 'list[' + ','.join(fmt_list(c) for c in cl) + ']'
    raise AssertionError

def parse_count(rh):
    """count= option, shared by both: returns ('k', k) or ('range', lo, hi)."""
    if isinstance(rh, Num) and is_int(rh):
        if rh.v <= 0: raise GErr('Expected a positive integer')
        return ('k', int(rh.v))
    if isinstance(rh, Rng):
        if not (is_int(rh.a) and is_int(rh.b)) or rh.a.v <= 0 or rh.b.v < rh.a.v:
            raise GErr('Expected a range of positive integers')
        return ('range', int(rh.a.v), int(rh.b.v))
    raise GErr('Expected an integer or range', True)

def is_sym(g, name): return isinstance(g, Sym) and g.n == name

def kmeans_cmd(args, rng):
    data = args[0].items if isinstance(args[0], Lst) else [args[0]]
    if not data: raise GErr('No data found')
    if not all(is_numeric(x) for x in data): raise GErr('Data must be numeric')
    NP = len(data); p = len(data[0].items) if isinstance(data[0], Lst) else 1
    out_what, min_k, max_k, k, maxiter = 2, 2, NP-1, 0, 100
    for opt in args[1:]:
        if isinstance(opt, Eq):
            lh, rh = opt.l, opt.r
            if is_sym(lh, 'count'):
                c = parse_count(rh)
                if c[0] == 'k': k = c[1]
                else: min_k, max_k = c[1], c[2]
            elif is_sym(lh, 'count_sup'):
                if not is_int(rh): raise GErr('Expected a positive integer')
                max_k = min(int(rh.v), NP-1)
            elif is_sym(lh, 'count_inf'):
                if not is_int(rh): raise GErr('Expected a positive integer')
                min_k = max(int(rh.v), 2)
            elif is_sym(lh, 'limit'):
                if not is_int(rh) or rh.v < 1: raise GErr('Expected a positive integer')
                maxiter = int(rh.v)
            elif is_sym(lh, 'index'):  # :11668-11686
                if is_sym(rh, 'true') or is_sym(rh, 'false'): k = -16 if rh.n == 'true' else 0
                elif isinstance(rh, Str):
                    cr = KMEANS_INDEX.get(rh.v, -1)
                    if cr < 0: raise GErr('Unknown index name')
                    k = -cr
                elif isinstance(rh, Lst) and rh.items and all(isinstance(t, Str) for t in rh.items):
                    agg = 0
                    for t in rh.items:
                        cr = KMEANS_INDEX.get(t.v, -1)
                        if cr < 0: raise GErr('Unknown index name')
                        agg |= cr
                    k = -agg
                else: raise GErr('Expected a string or list of strings', True)
            elif is_sym(lh, 'output'):
                if is_sym(rh, 'part'): out_what = 1
                elif is_sym(rh, 'list'): out_what = 2
                elif is_sym(rh, 'plot'): out_what = 3
                elif is_sym(rh, 'count'): out_what = 4
                elif is_sym(rh, 'index'): out_what = 5
                else: raise GErr('Invalid output specification')
            else: raise NotImplementedError('kmeans option %s' % lh)
        elif is_int(opt):
            k = int(opt.v)
            if k <= 0: raise GErr('Expected a positive integer')
        else: raise GErr('Invalid option')
    if out_what == 3 and p not in (2, 3): raise GErr('Plot output works only with 2D or 3D data')
    if out_what == 5 and k >= 0: raise GErr('No index specified for output')
    max_k = min(max_k, NP-1)
    if min_k > max_k: raise GErr('Invalid k-range')
    if min_k == max_k: k = max_k
    if p == 1 and not isinstance(data[0], Lst):
        # giac: mcols() reads _VECTptr of a scalar (vecteur.cc:3361): UB/crash. Prototype: 1-vectors.
        x = [[float(g.v)] for g in data]
    else:
        x = [[float(c.v) for c in g.items] for g in data]
    if k < 0:
        pts = [[c.v for c in g.items] for g in data] if isinstance(data[0], Lst) else [[g.v] for g in data]
        code, k, res, kvals, ind = kmeans_index(x, pts, -k, min_k, max_k, maxiter, rng)
        if code != 0: return '[]'
        if out_what == 5: return fmt_list('[%d,%s]' % (kv, fmt_val(v)) for kv, v in zip(kvals, ind))
        return clustering_result(data, p, out_what, k, res)
    code, k, res = kmeans_core(x, k, min_k, max_k, maxiter, rng)
    if code != 0: return '[]'
    if out_what == 5: return '[]'  # count=a..a with an index: k fixed, no index values (probe c3/p1)
    return clustering_result(data, p, out_what, k, res)

def cluster_cmd(args, want_dg=False, depart=False):
    data = args[0].items if isinstance(args[0], Lst) else [args[0]]
    if not data: raise GErr('No data found')
    distf = levenshtein if all(isinstance(x, Str) for x in data) else longueur2
    NP = len(data); p = len(data[0].items) if isinstance(data[0], Lst) else 1
    out_what, k, min_k, max_k, lmeth, indf = 2, 0, 0, 16, SINGLE, 0
    for opt in args[1:]:
        if isinstance(opt, Eq):
            lh, rh = opt.l, opt.r
            if is_sym(lh, 'count'):
                c = parse_count(rh)
                if c[0] == 'k': k = c[1]
                else: min_k, max_k = c[1], c[2]
            elif is_sym(lh, 'count_sup'):
                if not is_int(rh): raise GErr('Expected a positive integer')
                max_k = min(int(rh.v), NP-1)
            elif is_sym(lh, 'count_inf'):
                if not is_int(rh): raise GErr('Expected a positive integer')
                min_k = max(int(rh.v), 2)
            elif is_sym(lh, 'display'): pass  # drawing only
            elif is_sym(lh, 'type'):
                if not isinstance(rh, Str): raise GErr('Expected a string', True)
                lt = {'single': SINGLE, 'complete': COMPLETE, 'average': AVERAGE,
                      'weighted': WEIGHTED, 'ward': WARD}.get(rh.v)
                if lt is None: raise GErr('Unknown linkage type')
                lmeth = lt
            elif is_sym(lh, 'distance'):
                if is_sym(rh, 'levenshtein'): distf = levenshtein
                else: raise GErr('Invalid distance function specification', True)  # only _FUNC/program accepted
            elif is_sym(lh, 'output'):
                if is_sym(rh, 'part'): out_what = 1
                elif is_sym(rh, 'list'): out_what = 2
                elif is_sym(rh, 'plot'): out_what = 3
                elif is_sym(rh, 'count'): out_what = 4
                elif is_sym(rh, 'tree'): out_what = 5
                elif is_sym(rh, 'index'): out_what = 6
                else: raise GErr('Invalid output specification')
            elif is_sym(lh, 'index'):  # :11838-11853
                if is_sym(rh, 'true'): indf = 1
                elif is_sym(rh, 'false'): indf = 0
                elif isinstance(rh, Str):
                    indf = name2index(rh.v)
                    if indf < 0: raise GErr('Unknown index name')
                elif isinstance(rh, Lst) and rh.items and all(isinstance(x, Str) for x in rh.items):
                    indf = 0
                    for x in rh.items:
                        cr = name2index(x.v)
                        if cr < 0: raise GErr('Unknown index name')
                        indf |= cr
                else: raise GErr('Expected a string or list of strings', True)
            else: raise NotImplementedError('cluster option %s' % lh)
        elif is_int(opt):
            k = int(opt.v)
            if k <= 0: raise GErr('Expected a positive integer')
        else: raise GErr('Invalid option')
    if out_what == 3 and p not in (2, 3): raise GErr('Plot output works only with 2D or 3D data')
    if out_what == 6 and indf <= 0: raise GErr('No index specified for output')
    if max_k < min_k: raise GErr('Invalid k-range')
    if max_k == min_k: k = min_k
    hc = HClust(data, distf)
    k, ind, dg, kvals, vals = hc.linkage(lmeth, k if k > 0 else -indf, min_k, max_k, depart)
    if k <= 0: raise GErr('Bad distance function')
    if want_dg: return fmt_list(ind), k, dg
    if out_what == 6:  # mtran([kvals, vals])
        if not vals or len(vals) != len(kvals):
            e = GErr(''); e.txt = 'GIAC_ERROR: Unable to transpose Error: Bad Argument Value'; raise e
        return '[' + ','.join('[%d,%r]' % kv for kv in zip(kvals, vals)) + ']'
    if out_what == 5: return 'Graphic_object'  # dendrogram drawing not ported
    return clustering_result(data, p, out_what, k, ind)

def run(line, seed, depart=False):
    s, name, args = parse_call(line)
    rng = random.Random(seed)
    try:
        r = kmeans_cmd(args, rng) if name == 'kmeans' else cluster_cmd(args, depart=depart)
    except GErr as e:
        return e.txt
    except Abort:
        return '<abort>'
    return r if s is None else '%d,%s' % (s, r)

# ---------------- comparison ----------------
def labels_of(ans):
    if ans.startswith('[title='):
        return [int(v) - 235405312 for v in re.findall(r'display=(\d+)', ans)]
    m = re.search(r'\[([\d,]*)\]$', ans)
    return [int(v) for v in m.group(1).split(',')] if m and m.group(1) else None

def canon(lab):
    mp = {}
    return tuple(mp.setdefault(v, len(mp)) for v in lab)

def strip_labels(ans):
    return re.sub(r'display=\d+', 'display=_', re.sub(r'\[[\d,]*\]$', '[_]', ans))

def same(name, giac, ours):
    if name != 'kmeans' or giac.startswith('GIAC_ERROR'): return giac == ours
    lg, lo = labels_of(giac), labels_of(ours)
    if lg is None or lo is None: return giac == ours
    return canon(lg) == canon(lo) and strip_labels(giac) == strip_labels(ours)

def main():
    rows = []
    for f in ('probe4', 'probe8'):
        cas = open(SRC + f + '.cas').read().splitlines()
        out = open(SRC + f + '.out').read().splitlines()
        for c, o in zip(cas, out):
            pre = c + '  =>  '
            assert o.startswith(pre), o
            rows.append((f, c, o[len(pre):].rstrip()))
    nok = nfail = nub = 0
    seeddep = []
    for f, c, giac in rows:
        name = parse_call(c)[1]
        if giac == '<stderr>':
            ours = run(c, 0)
            print('%s | %s | giac: <crash: mcols on scalar data, UB> | ours: %s | UB' % (f, c, ours))
            nub += 1
        else:
            ours = run(c, 0); st = 'OK' if same(name, giac, ours) else 'FAIL'
            if st == 'FAIL' and name == 'kmeans':
                # giac's _rand stream is not reproducible: accept the first Python seed that reproduces it
                for sd in range(1, NSEEDS):
                    if same(name, giac, run(c, sd)):
                        ours = run(c, sd); st = 'OK(python seed %d)' % sd; break
            print('%s | %s | giac: %s | ours: %s | %s' % (f, c, giac, ours, st))
            st = 'FAIL' if st == 'FAIL' else 'OK'
            nok += st == 'OK'; nfail += st == 'FAIL'
        if name == 'kmeans' and labels_of(ours) is not None:
            outs = {}
            for sd in range(NSEEDS):
                r = run(c, sd)
                key = canon(labels_of(r)) if labels_of(r) is not None else r
                outs.setdefault(key, []).append(sd)
            if len(outs) > 1:
                seeddep.append((c, outs))
            print('    seeds 0..%d: %s' % (NSEEDS-1, 'seed-independent' if len(outs) == 1 else
                  'SEED-DEPENDENT ' + '; '.join('%s x%d' % (list(k_), len(v)) for k_, v in outs.items())))
    print('TOTAL rows=%d OK=%d FAIL=%d UB=%d seed-dependent kmeans rows=%d' % (len(rows), nok, nfail, nub, len(seeddep)))
    extras()
    nfail += t5()
    nfail += t5c3()
    return 0 if nfail == 0 else 1

NUM = re.compile(r'-?\d+(?:\.\d*)?(?:e[-+]?\d+)?')

def close(giac, ours):
    """The same text up to giac's 12-digit decimals (relative 1e-10)."""
    if NUM.sub('#', giac) != NUM.sub('#', ours): return False
    return all(abs(float(a) - float(b)) <= 1e-10 * max(1.0, abs(float(a)))
               for a, b in zip(NUM.findall(giac), NUM.findall(ours)))

def groups(ans):
    """A part answer as a set of groups (kmeans' labels are seed-dependent)."""
    body = ans[ans.index('list['):] if 'list[' in ans else ans
    return frozenset(re.findall(r'\[((?:\[[^\[\]]*\],?)+|[^\[\]]*)\]', body[5:-1]))

def t5():
    """T5 (CAS-646): output=part, string data, cluster index= against probe_clustering2 (1.8/T5)."""
    cas = open(SRC + 'probe_clustering2.cas').read().splitlines()
    out = open(SRC + 'probe_clustering2.out').read().splitlines()
    nfail = 0
    for c, o in zip(cas, out):
        pre = c + '  =>  '
        assert o.startswith(pre), o
        giac = o[len(pre):].rstrip()
        try:
            seed, name, _ = parse_call(c)
        except Exception:
            name = None
        if name not in ('cluster', 'kmeans'):
            print('T5 | %s | skip (not cluster/kmeans)' % c); continue
        ours = run(c, seed or 0)
        if name == 'kmeans' and 'index=' in c and not giac.startswith('GIAC'):
            g = giac.split(',', 1)[1] if seed is not None else giac
            st = 'OK' if reachable(re.sub(r'^srand\(\d+\);', '', c), g) else 'FAIL'
            print('T5 | %s | giac: %s | kmeans index=, over %d seeds (T5.c3) | %s' % (c, giac, NSEEDS, st))
            nfail += st == 'FAIL'; continue
        if seed is not None and not ours.startswith('GIAC_ERROR'): ours = '%d,%s' % (seed, ours) if ',' not in ours[:len(str(seed)) + 1] else ours
        if giac.startswith('<stderr>'):
            st = 'OK(abort)' if ours == '<abort>' else 'FAIL'
            print('T5 | %s | giac: abort | ours (departure, max_k<=N-1): %s | %s' % (c, run(c, 0, True), st))
        elif name == 'kmeans' and 'list[' in giac:
            seen = {}
            for sd in range(NSEEDS):
                seen.setdefault(groups(run(c, sd)), []).append(sd)
            st = 'OK' if groups(giac) in seen else 'FAIL'
            print('T5 | %s | giac: %s | partitions over %d seeds: %d %s | %s' % (c, giac, NSEEDS, len(seen),
                  ['x%d' % len(v) for v in seen.values()], st))
        else:
            g = giac.split(',', 1)[1] if seed is not None and not giac.startswith('GIAC') else giac
            o2 = ours.split(',', 1)[1] if seed is not None and not ours.startswith('GIAC') else ours
            st = 'OK' if close(g, o2) else 'FAIL'
            print('T5 | %s | giac: %s | ours: %s | %s' % (c, giac, ours, st))
        nfail += st == 'FAIL'
    print('T5 TOTAL FAIL=%d' % nfail)
    return nfail

def reachable(line, g):
    """Whether some Python seed of 0..NSEEDS-1 gives giac's answer g to line (giac's _rand stream is not
    ours): labels up to renaming, parts as sets of groups, index values within 12 digits, the whole
    vector or else each [k, value] by some seed. Returns None or 'seed' or 'per-k'."""
    def key(a):
        if a.startswith('list['): return ('part', groups(a))
        lab = labels_of(a)
        if lab is not None and not a.startswith('[['): return ('lab', canon(lab))
        return ('txt', a)
    runs = [run(line, sd) for sd in range(NSEEDS)]
    kg = key(g)
    if any(key(r) == kg if kg[0] != 'txt' else close(g, r) for r in runs): return 'seed'
    if g.startswith('[['):
        per = {}
        for r in runs:
            for kk, v in re.findall(r'\[(\d+),([^\[\]]*)\]', r): per.setdefault(kk, []).append(v)
        if all(any(close(v, w) for w in per.get(kk, [])) for kk, v in re.findall(r'\[(\d+),([^\[\]]*)\]', g)):
            return 'per-k'
    return None

def t5c3():
    """T5.c3 (CAS-646): kmeans index= against probe_clustering3, each row reachable over the seeds."""
    cas = open(SRC + 'probe_clustering3.cas').read().splitlines()
    out = open(SRC + 'probe_clustering3.out').read().splitlines()
    nfail = 0; nind = 0
    for c, o in zip(cas, out):
        pre = c + '  =>  '
        assert o.startswith(pre), o
        g = o[len(pre):].rstrip()
        line = re.sub(r'^srand\(\d+\);', '', c)
        if re.match(r'srand\(\d+\);', c) and not g.startswith('GIAC'): g = g.split(',', 1)[1]
        r = reachable(line, g)
        nind += r == 'per-k'; nfail += r is None
        if r is None: print('T5c3 FAIL | %s | giac: %s | seed 0: %s' % (c, g, run(line, 0)))
    print('T5c3 rows=%d per-k=%d FAIL=%d' % (len(cas), nind, nfail))
    return nfail

E6 = '[[0,0],[0,1],[5,5],[5,6],[10,0],[10,1]]'
P8 = '[[0,0],[1,0],[0,2],[4,4],[5,5],[4,6],[9,1],[10,0]]'

def extras():
    print('EXTRA kmeans (Python seeds 0..%d; canonical partition xcount):' % (NSEEDS-1))
    for c in ('kmeans([1,2,10,11,20,21],3)', 'kmeans([1,2,10,11,20,21])', 'kmeans([1,2,10,11],2)',
              'kmeans([1,2,10,11,20,21],output=count)', 'kmeans(%s,output=count)' % E6,
              'kmeans(%s,count=2..4,output=count)' % E6, 'kmeans(%s,output=count)' % P8,
              'kmeans(%s,3)' % P8, 'kmeans([[1,2],[3,4],[5,6]],5)'):
        outs = {}
        for sd in range(NSEEDS):
            r = run(c, sd)
            outs.setdefault(r if labels_of(r) is None else str(list(canon(labels_of(r)))), []).append(sd)
        print('  %s -> %s' % (c, '; '.join('%s x%d' % (k_, len(v)) for k_, v in outs.items())))
    print('EXTRA cluster (labels, k, sorted dendrogram heights):')
    cs = ['cluster(%s,type="%s",count=3)' % (E6, t) for t in ('single', 'complete', 'average', 'weighted', 'ward')]
    cs += ['cluster(%s,type="%s",count=2..4)' % (E6, t) for t in ('single', 'complete', 'average', 'weighted', 'ward')]
    cs += ['cluster(%s,type="ward")' % P8, 'cluster(%s,type="ward",count=3)' % P8, 'cluster(%s)' % P8,
           'cluster(%s,type="complete")' % P8, 'cluster(%s,type="average",count=3)' % P8]
    for c in cs:
        _, _, args = parse_call(c)
        lab, k, dg = cluster_cmd(args, want_dg=True)
        print('  %s -> %s k=%d heights=%s' % (c, lab, k, [repr(d) for d, _ in dg]))

if __name__ == '__main__':
    sys.exit(main())
