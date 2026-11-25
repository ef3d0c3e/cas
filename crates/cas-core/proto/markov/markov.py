"""Prototype of giac's markov (misc.cc:7082-7400) with our exact solves (CAS-628).

Entries: Fraction (exact) or float (decimal, C0060). Prints each probe row: ours vs giac.
Origin: phase 1.7, T4.c1 scratch. Run: /usr/bin/python3 markov.py <probe.out>...
"""
from fractions import Fraction as F
import re, sys

EPS = 1e-12  # giac epsilon(); giac keeps ~45 mantissa bits, so its boundary moves by a few 2^-45 (probe_d/e)


def is_dec(x):
    return isinstance(x, float)


def is_stochastic(m):
    n = len(m)
    if n == 0 or any(not isinstance(r, list) or len(r) != n for r in m):
        return None
    for r in m:
        for x in r:
            if not isinstance(x, (F, float)):
                return None
            if x < 0:
                return None

    def ones(rows):
        for r in rows:
            if any(is_dec(x) for x in r):
                if abs(sum(float(x) for x in r) - 1) > EPS:
                    return False
            elif sum(r) != 1:
                return False
        return True

    if ones(m):
        return m
    t = [list(c) for c in zip(*m)]
    return t if ones(t) else None


def adjacency(m):
    return [[j for j, x in enumerate(r) if float(x) != 0.0] for r in m]


def tarjan(g):
    n = len(g)
    index = [-1] * n
    low = [-1] * n
    stack, ins, sccs = [], [False] * n, []
    counter = [0]

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        ins[v] = True
        for w in g[v]:
            if index[w] == -1:
                strong(w)
                low[v] = min(low[v], low[w])
            elif ins[w]:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            scc = []
            while stack:
                x = stack.pop()
                ins[x] = False
                scc.append(x)
                if x == v:
                    break
            sccs.append(scc)

    for v in range(n):
        if index[v] == -1:
            strong(v)
    return sccs


def classify(g, sccs):
    rec, trans = [], []
    for c in sccs:
        s = set(c)
        (rec if all(w in s for v in c for w in g[v]) else trans).append(c)
    return rec, trans


def solve(a, b):
    """Gauss-Jordan on a square system a x = b, exact or float."""
    n = len(a)
    rows = [list(r) + [b[i]] for i, r in enumerate(a)]
    for j in range(n):
        p = max(range(j, n), key=lambda i: abs(float(rows[i][j])))
        rows[j], rows[p] = rows[p], rows[j]
        piv = rows[j][j]
        assert piv != 0
        rows[j] = [x / piv for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j] != 0:
                c = rows[i][j]
                rows[i] = [x - c * y for x, y in zip(rows[i], rows[j])]
    return [r[n] for r in rows]


def markov(m):
    M = is_stochastic(m)
    if M is None:
        return "Not a stochastic matrix!"
    n = len(M)
    dec = any(is_dec(x) for r in M for x in r)
    one = 1.0 if dec else F(1)
    g = adjacency(M)
    rec, trans = classify(g, tarjan(g))
    inv = []
    for c in rec:
        v = sorted(c)
        k = len(v)
        # pi (P_C - I) = 0, sum = 1: A = (P_C - I)^T with the last row ones
        a = [[M[v[j]][v[i]] - (one if i == j else 0) for j in range(k)] for i in range(k)]
        a[-1] = [one] * k
        rhs = [0] * (k - 1) + [one]
        pi = solve(a, rhs)
        w = [0] * n
        for s, p in zip(v, pi):
            w[s] = p
        inv.append(w)
    if len(rec) == 1:
        absorb = [[1] for _ in range(n)]
    else:
        t = sorted(s for c in trans for s in c)
        pos = {s: i for i, s in enumerate(t)}
        # (I - Q) B = R, R[t][c] = sum_{j in class c} P[t][j]
        a = [[(one if i == j else 0) - M[t[i]][t[j]] for j in range(len(t))] for i in range(len(t))]
        cols = []
        for c in rec:
            r = [sum((M[s][j] for j in c), 0 * one) for s in t]
            cols.append(solve(a, r) if t else [])
        absorb = []
        for s in range(n):
            ci = [i for i, c in enumerate(rec) if s in c]
            if ci:
                absorb.append([1 if i == ci[0] else 0 for i in range(len(rec))])
            else:
                absorb.append([col[pos[s]] for col in cols])
    return [rec, inv, trans, absorb]


def num(text):
    text = text.strip()
    if re.fullmatch(r"-?\d+/\d+|-?\d+", text):
        return F(text)
    if re.fullmatch(r"[-+0-9.e]+", text) and re.search(r"\d[-+]\d", text):
        return float(eval(text))
    return float(text)


def parse(text):
    """giac list text -> nested lists of Fraction/float; None for symbolic."""
    text = text.replace(" ", "")
    pos = [0]

    def item():
        if text[pos[0]] == "[":
            pos[0] += 1
            out = []
            if text[pos[0]] == "]":
                pos[0] += 1
                return out
            while True:
                out.append(item())
                c = text[pos[0]]
                pos[0] += 1
                if c == "]":
                    return out
        m = re.match(r"[^,\]]+", text[pos[0]:])
        pos[0] += m.end()
        tok = m.group(0)
        try:
            return num(tok)
        except ValueError:
            return tok

    items = [item()]
    while pos[0] < len(text) and text[pos[0]] == ",":
        pos[0] += 1
        items.append(item())
    return items


def close(a, b, tol):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, list) or isinstance(b, list):
        return False
    if is_dec(a) or is_dec(b):
        return abs(float(a) - float(b)) <= tol * max(1, abs(float(b)))
    return a == b


def main(paths):
    bad = 0
    for path in paths:
        for line in open(path):
            if "=>" not in line or "markov(" not in line or "randmarkov" in line:
                continue
            inp, out = [s.strip() for s in line.split("=>", 1)]
            m = re.fullmatch(r"markov\((.*)\)", inp)
            args = parse(m.group(1))
            eps = args[1] if len(args) > 1 else None
            if eps == "sqrt(2)":
                eps = 2 ** 0.5
            if eps is not None and not isinstance(eps, (F, float)):
                ours = "Bad Argument Value"
            else:
                mat = args[0]
                ours = markov(mat) if isinstance(mat, list) and all(isinstance(r, list) for r in mat) else "Not a stochastic matrix!"
            if out.startswith("GIAC_ERROR") or out.startswith("<timeout"):
                ok = isinstance(ours, str) and ours in out
                tag = "ok" if ok else ("hang/ours answers" if out.startswith("<timeout") else "DIFF")
                print(f"{tag:5} {inp} -> {ours!r} | giac {out}")
                bad += tag == "DIFF"
                continue
            giac = parse(out)
            # giac iterates (absorption, 1-state classes): within eps or 1e-3 for its float->rational
            tol = float(eps) if eps is not None and float(eps) > 0 else 1e-9
            ok = close(ours, giac, max(tol, 1e-9))
            if not ok:
                ok2 = close(fl(ours), fl(giac), max(tol, 1e-9))
                tag = "dep" if ok2 else "DIFF"
            else:
                tag = "ok"
            bad += tag == "DIFF"
            print(f"{tag:5} {inp} -> {fmt(ours)} | giac {out[:160]}")
    print("DIFF rows:", bad)
    return bad


def fl(x):
    if isinstance(x, list):
        return [fl(y) for y in x]
    return float(x) if not isinstance(x, str) else x


def fmt(x):
    if isinstance(x, list):
        return "[" + ",".join(fmt(y) for y in x) + "]"
    if isinstance(x, F):
        return str(x)
    return repr(round(x, 12)) if not isinstance(x, str) else x


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1:]) else 0)
