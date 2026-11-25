"""T8.c1 scratch (origin: phase 1.7, T8.c1): giac's set_channel_data, stereo2mono (to_mono), mixdown (mix), splice (similarity) in Python.

Clips and createwav come from ../clip/proto.py (47-bit giac doubles) and proto_f64.py (plain f64).
Run: /usr/bin/python3 proto.py  -> writes table.txt (row-by-row PASS for both models), prints one count line.
"""
import importlib.util, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)) + "/"
T7 = HERE + "../clip/"
GIAC = HERE + "../../giac/"


def load(name):
    spec = importlib.util.spec_from_file_location(name, T7 + name + ".py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


P47, P64 = load("proto"), load("proto_f64")
ClipError = P47.ClipError
GENERIC = "<generic giac error>"   # gentypeerr / gendimerr: giac prints the input


def err(m):
    raise ClipError(m)


class Model:
    """giac's four functions over one T7 model (P47 or P64); strict=True is ours: splice checks params always."""

    def __init__(self, P, strict=False):
        self.P, self.strict = P, strict

    # ---- clip helpers (signed codes as T7)
    def is_clip(self, c):
        return isinstance(c, dict) and "rate" in c

    def dec(self, c, k):
        return float(k) / float(1 << (c["depth"] - 1))

    def enc(self, depth, s):
        return P47.encode_f64(depth, s)

    def new(self, depth, rate, nc, n):
        return dict(rate=rate, depth=depth, channels=[[0] * n for _ in range(nc)])

    def copy(self, c):
        return dict(rate=c["rate"], depth=c["depth"], channels=[list(r) for r in c["channels"]])

    def n(self, c):
        return len(c["channels"][0])

    def isnum(self, x):
        return self.P.is_num(x)

    def gd(self, x):
        return self.P.gd(x)

    # ---- set_channel_data :2554-2602, audio_clip::set_channel_data/set_block :1515-1544
    def set_cd_one(self, c, chn, data, offset):
        if chn < 0 or chn >= len(c["channels"]):
            err("Invalid channel number")
        L, n = len(data), self.n(c)
        if offset >= n:
            return
        if offset > n - L:
            L = n - offset
        for pos in range(L):
            v = data[pos]
            if not self.isnum(v):
                err(f"Failed to set channel data, {P47.gprint(v)} is not a valid sample value")
            c["channels"][chn][offset + pos] = self.enc(c["depth"], 1.0 * self.gd(v))

    def set_channel_data(self, *args):
        if len(args) < 2 or len(args) > 5:
            err(GENERIC)
        c = args[0]
        if not self.is_clip(c):
            err("Expected audio clip or image")
        nc = len(c["channels"])
        data = args[1]
        if not isinstance(data, list):
            err("Invalid data specification")
        mat = P47.ckmatrix(data)
        if mat and len(data) != nc:
            err("Invalid data specification")
        offset, chn = 0, -1
        if len(args) > 2:
            if not isinstance(args[2], int):
                err("Expected an integer")
            offset = args[2]
            if offset < 0 or offset >= self.n(c):
                err("Invalid offset")
        if len(args) > 3:
            ch = args[3]
            if ch == "left":
                chn = 0
            elif ch == "right":
                chn = 1
            elif isinstance(ch, int):
                chn = ch
                if chn < 0 or chn >= nc:
                    err("Invalid channel number")
            else:
                err("Invalid channel specification")
            if mat:
                err("Invalid data specification")
        r = self.copy(c)
        if chn < 0:
            if mat:
                for i, row in enumerate(data):
                    self.set_cd_one(r, i, row, offset)
            else:
                for i in range(nc):
                    self.set_cd_one(r, i, data, offset)
        else:
            self.set_cd_one(r, chn, data, offset)
        return r

    # ---- stereo2mono :2339, to_mono :1223-1232
    def stereo2mono(self, *args):
        if len(args) != 1 or not self.is_clip(args[0]):
            err("Invalid sound data")
        c = args[0]
        nc = len(c["channels"])
        r = self.new(c["depth"], c["rate"], 1, self.n(c))
        for i in range(self.n(c)):
            tot = 0.0
            for ch in range(nc):
                tot += self.dec(c, c["channels"][ch][i])
            r["channels"][0][i] = self.enc(c["depth"], tot / nc)
        return r

    # ---- mixdown :1791-1843, is_audio_position :1777, mix :1164-1190
    def position(self, g, sr):
        if isinstance(g, int) and g >= 0:
            return g
        if self.isnum(g) and not isinstance(g, int) and self.gd(g) >= 0:
            return math.floor(0.5 + sr * self.gd(g))
        return None

    def mix(self, ret, o, offset, gain, pan, fi, fo):
        nc, onc = len(ret["channels"]), len(o["channels"])
        ratio = math.pow(10.0, gain / 20.0)
        theta = math.pi * (pan + 1.0) / 4.0
        pl, pr = math.cos(theta), math.sin(theta)
        L = self.n(o)
        if nc == 1:
            pl = 1.0
        clipping = anyclip = False
        for i in range(L):
            r = ratio
            if i < fi:
                r *= math.pow(math.cos((math.pi / 2) * (1.0 - float(i) / float(fi))), 2)
            if i >= L - fo:
                r *= math.pow(math.cos((math.pi / 2) * (1.0 - float(L - i - 1) / float(fo))), 2)
            for c in range(nc):
                s1 = self.dec(ret, ret["channels"][c][i + offset])
                s2 = (pl if c == 0 else pr) * r * self.dec(o, o["channels"][min(onc - 1, c)][i])
                s = s1 + s2
                clipping = s > 1.0 or s < -1.0
                anyclip |= clipping
                ret["channels"][c][i + offset] = self.enc(ret["depth"], s)
        return not clipping, anyclip

    def mixdown(self, *args):
        s = len(args)
        if s <= 0 or s % 2:
            err(GENERIC)
        k = s // 2
        sr = bd = nc = L = 0
        parms = []
        for i in range(k):
            c = args[2 * i]
            if not self.is_clip(c):
                err(f"Invalid argument {2 * i + 1}")
            if len(c["channels"]) > 2:
                err("Only mono and stereo clips are allowed")
            if sr == 0:
                sr = c["rate"]
            elif sr != c["rate"]:
                err("All clips must have the same sample rate")
            bd, nc = max(bd, c["depth"]), max(nc, len(c["channels"]))
            p = args[2 * i + 1]
            if isinstance(p, list) and len(p) > 5:
                err(GENERIC)
            parm = p if isinstance(p, list) else [p]
            off, gain, pan, fi, fo = 0, 0.0, 0.0, 0, 0
            if len(parm) > 0:
                off = self.position(parm[0], sr)
                if off is None:
                    err("Invalid offset specification")
            if len(parm) > 1:
                if not self.isnum(parm[1]) or self.gd(parm[1]) > 0:
                    err("Invalid gain specification")
                gain = self.gd(parm[1])
            if len(parm) > 2:
                if not self.isnum(parm[2]) or not (-1.0 <= self.gd(parm[2]) <= 1.0):
                    err("Invalid pan specification")
                pan = self.gd(parm[2])
            if len(parm) > 3:
                fi = self.position(parm[3], sr)
                if fi is None:
                    err("Invalid fade-in specification")
            if len(parm) > 4:
                fo = self.position(parm[4], sr)
                if fo is None:
                    err("Invalid fade-out specification")
            if fi + fo > self.n(c):
                err("Fade-in/out is too long")
            L = max(L, off + self.n(c))
            parms.append((c, off, gain, pan, fi, fo))
        ret = self.new(bd, sr, nc, L)
        clip_giac = clip_ours = False
        for c, off, gain, pan, fi, fo in parms:
            ok, anyc = self.mix(ret, c, off, gain, pan, fi, fo)
            clip_giac |= not ok
            clip_ours |= anyc
        self.warned = clip_ours if self.strict else clip_giac
        return ret

    # ---- splice :1738-1773, audio_clip::splice :1194-1221, similarity :1130-1149
    def same(self, a, b):
        return a["depth"] == b["depth"] and a["rate"] == b["rate"] and len(a["channels"]) == len(b["channels"])

    def similarity(self, a, b, offset):
        if offset < 0 or offset >= self.n(a):
            err("Invalid chunk position and/or length specification")
        if not self.same(a, b):
            err("Clips should have the same bit depths, sample rates and channel counts")
        l = min(self.n(a) - offset, self.n(b))
        p = 1.0
        for c in range(len(a["channels"])):
            pp = na = nb = 0.0
            for i in range(l):
                x = self.dec(a, a["channels"][c][i + offset])
                y = self.dec(b, b["channels"][c][i])
                pp += x * y
                na += x * x
                nb += y * y
            d = math.sqrt(na * nb)
            p *= (pp / d) if d != 0 else (math.nan if pp == 0 else math.copysign(math.inf, pp))
        try:
            return math.pow(p, 1.0 / len(a["channels"]))
        except ValueError:      # C pow(negative, non-integer) = NaN
            return math.nan

    def splice(self, *args):
        if len(args) < 2 or len(args) > 4:
            err(GENERIC)
        a, b = args[0], args[1]
        if not self.is_clip(a) or not self.is_clip(b):
            err("First two arguments should be audio clips")
        pos, sr, L, ex = -1, a["rate"], self.n(a), -1.0
        if len(args) > 2:
            x = args[2]
            if not self.isnum(x):
                err("Crossfade length should be a real number")
            if isinstance(x, float):
                pos = L - math.floor(0.5 + sr * self.gd(x))
            elif isinstance(x, int):
                pos = L - x
            else:
                err(GENERIC)
        if len(args) > 3:
            if not self.isnum(args[3]):
                err("Crossfade exponent should be a real number")
            ex = self.gd(args[3])
        if pos < 0:
            da, db = L / float(sr), self.n(b) / float(b["rate"])
            pos = L - math.floor(0.5 + sr * min(min(da, db) / 2.0, 0.5))
        if pos < 0 or pos >= L:
            err("Invalid splice position")
        if self.strict and not self.same(a, b):
            err("Clips should have the same bit depths, sample rates and channel counts")
        e = 0.5 + 0.5 * abs(self.similarity(a, b, pos)) if ex < 0 else ex
        nc, Lb = len(a["channels"]), self.n(b)
        ret = self.new(a["depth"], a["rate"], nc, pos + Lb)
        RL = pos + Lb
        cflen = min(L - pos, Lb)
        pos = L - cflen
        for c in range(nc):
            for i in range(min(pos, RL)):
                ret["channels"][c][i] = a["channels"][c][i]
            for j in range(Lb - cflen):
                if L + j < RL:
                    ret["channels"][c][L + j] = b["channels"][c][cflen + j]
        clipping = anyclip = False
        for c in range(nc):
            for i in range(cflen):
                x = self.dec(a, a["channels"][c][i + pos])
                y = self.dec(b, b["channels"][c][i])
                if e == 0:
                    v = x + y
                else:
                    cf = math.pow(float(i + 1) / float(cflen), e)
                    v = x * (1.0 - cf) + y * cf
                clipping = abs(v) > 1.0
                anyclip |= clipping
                if 0 <= i + pos < RL:
                    ret["channels"][c][i + pos] = self.enc(a["depth"], v)
        self.warned = anyclip if self.strict else clipping
        return ret

    # ---- channel_data / channels on results
    def cd(self, c, *opt):
        rows = [[self.dec(c, k) for k in ch] for ch in c["channels"]]
        return rows[0] if len(rows) == 1 else rows

    def channels(self, c):
        return len(c["channels"])

    def bit_depth(self, c):
        return c["depth"]

    def W(self, *a):
        return self.P.createwav(*a)


# ---------------------------------------------------------------- rows
A = [0.5, -0.5, 0.25, 0]
B = [0.1, 0.2, 0.3, 0.4]
S = [[0.5, -0.5, 0.25, 0], [0.1, 0.2, 0.3, 0.4]]
Z = [0, 0, 0, 0]
O = P47.Opt
F = P47.F


def sqrt2(P):
    import mpmath
    return P.Irr(mpmath.sqrt(2), P.trunc5(math.sqrt(2)))


def R(m):
    """rows as (file, line, fn(model)); the giac line's answer is the oracle."""
    W, O = m.W, m.P.Opt
    return [
        ("probe5", 38, lambda: m.cd(m.stereo2mono(W(S)))),
        ("probe5", 39, lambda: m.channels(m.stereo2mono(W(S)))),
        ("probe5", 40, lambda: m.cd(m.set_channel_data(W(A), [0.9, 0.8]))),
        ("probe5", 41, lambda: m.cd(m.set_channel_data(W(A), [0.9, 0.8], 1))),
        ("probe5", 42, lambda: m.cd(m.set_channel_data(W(S), [0.9], 2, "right"))),
        ("probe5", 43, lambda: m.cd(m.set_channel_data(W(A), [0.9, 0.8, 0.7, 0.6, 0.5], 2))),
        ("probe5", 44, lambda: m.cd(m.splice(W(A), W(B)))),
        ("probe5", 45, lambda: m.cd(m.splice(W(A), W(B), 2, 1))),
        ("probe5", 46, lambda: m.cd(m.splice(W(A), W(B), 2, 0.5))),
        ("probe5", 47, lambda: m.cd(m.splice(W(A), W(B), 2, 2))),
        ("probe5", 48, lambda: m.cd(m.mixdown(W(A), W(B)))),
        ("probe5", 49, lambda: m.cd(m.mixdown(W(A), 1, W(B), 1))),
        ("probe5", 50, lambda: m.cd(m.mixdown(W(A), 0.5, W(B), 0.5))),
        ("probe5", 51, lambda: m.cd(m.mixdown(W(A), [0.5, -1], W(B), [0.5, 1]))),
        ("probe5", 52, lambda: m.cd(m.mixdown(W(A), -6, W(B), 0))),
        ("probe9", 1, lambda: m.cd(m.mixdown(W(A), [1, 0], W(B), [1, 0]))),
        ("probe9", 2, lambda: m.cd(m.mixdown(W(A), [1, 0, 0], W(B), [1, 0, 0]))),
        ("probe9", 3, lambda: m.cd(m.mixdown(W(A), 0, W(B), 0))),
        ("probe9", 4, lambda: m.cd(m.mixdown(W(A), 2, W(B), 0))),
        ("probe9", 5, lambda: m.cd(m.mixdown(W(A), [0, 0], W(B), [0, 2]))),
        ("probe9", 6, lambda: m.cd(m.mixdown(W(A), [0, -6], W(B), [0, -6]))),
        ("probe9", 7, lambda: m.cd(m.mixdown(W(A), [0, 0, -1], W(B), [0, 0, 1]))),
        ("probe9", 8, lambda: m.channels(m.mixdown(W(A), [0, 0, -1], W(B), [0, 0, 1]))),
        ("probe9", 9, lambda: len(m.cd(m.mixdown(W(A), 0.5, W(B), 0.5)))),
        ("probe9", 10, lambda: len(m.cd(m.mixdown(W(A), 0.0001, W(B), 0)))),
        ("probe9", 11, lambda: m.cd(m.splice(W(A), W(B), 0))),
        ("probe9", 12, lambda: m.cd(m.splice(W(A), W(B), 1))),
        ("probe9", 13, lambda: m.cd(m.splice(W(A), W(B), 3))),
        ("probe9", 14, lambda: m.cd(m.splice(W(A), W(B), 4))),
        ("probe9", 15, lambda: m.cd(m.splice(W(A), W(B), 5))),
        ("probe9", 16, lambda: m.cd(m.splice(W(A), W(B, O("samplerate", 22050))))),
        ("probe9", 17, lambda: m.cd(m.splice(W([[0.5, -0.5, 0.25, 0], Z]), W(B), 2))),
        ("probe9", 18, lambda: m.cd(m.splice(W(A), W(B), 2, -1))),
        ("probe9", 45, lambda: m.cd(m.stereo2mono(W([0.5, 0.25])))),
        ("probe9", 46, lambda: m.cd(m.stereo2mono([0.5, 0.25]))),
        # sound/edit/probe1 (this task's probes)
        ("../sound/edit/probe1", 1, lambda: len(m.cd(m.mixdown(W(A), 0.0, W(B), 0)))),
        ("../sound/edit/probe1", 2, lambda: len(m.cd(m.mixdown(W(A), F(1, 2), W(B), 0)))),
        ("../sound/edit/probe1", 3, lambda: len(m.cd(m.mixdown(W(A), sqrt2(m.P), W(B), 0)))),
        ("../sound/edit/probe1", 4, lambda: len(m.cd(m.mixdown(W(A), 0, W(B), 1e-5)))),
        ("../sound/edit/probe1", 5, lambda: P47.text(m.mixdown(W(A), 0))),
        ("../sound/edit/probe1", 6, lambda: m.mixdown(W(A))),
        ("../sound/edit/probe1", 7, lambda: m.mixdown(W(A), 0, W(B))),
        ("../sound/edit/probe1", 8, lambda: m.mixdown(W(A), 0, [1, 2], 0)),
        ("../sound/edit/probe1", 9, lambda: m.cd(m.mixdown(W(A), 0, W(S), 0))),
        ("../sound/edit/probe1", 10, lambda: m.cd(m.mixdown(W(A), [0, 0, 0.5], W(S), 0))),
        ("../sound/edit/probe1", 11, lambda: m.cd(m.mixdown(W(A), [0, 0, -0.5], W(B), [1, -3, 0.5]))),
        ("../sound/edit/probe1", 12, lambda: m.cd(m.mixdown(W(S), [0, 0, -0.5], W(B), [1, -3, 0.5]))),
        ("../sound/edit/probe1", 13, lambda: m.cd(m.mixdown(W(A), [0, 0, 0, 2, 2], W(B), 1))),
        ("../sound/edit/probe1", 14, lambda: m.cd(m.mixdown(W(A), [0, 0, 0, 1, 0], W(B), [0, 0, 0, 0, 3]))),
        ("../sound/edit/probe1", 15, lambda: m.mixdown(W(A), [0, 0, 0, 3, 2], W(B), 1)),
        ("../sound/edit/probe1", 16, lambda: m.mixdown(W(A), [0, 0, 0, -1], W(B), 1)),
        ("../sound/edit/probe1", 17, lambda: m.mixdown(W(A), [0, 0, 0, 0, "x"], W(B), 1)),
        ("../sound/edit/probe1", 18, lambda: m.mixdown(W(A), [0, 0, 2], W(B), 1)),
        ("../sound/edit/probe1", 19, lambda: m.mixdown(W(A), [0, 0, 0, 0, 0, 0], W(B), 1)),
        ("../sound/edit/probe1", 20, lambda: m.mixdown(W(A), 0, W([0.1, 0.2], O("samplerate", 22050)), 0)),
        ("../sound/edit/probe1", 21, lambda: m.mixdown(W([0.1, 0.2], O("channels", 3)), 0, W(B), 0)),
        ("../sound/edit/probe1", 22, lambda: m.cd(m.mixdown(W([0.5, -0.5], O("bit_depth", 8)), 0, W(B), 1))),
        ("../sound/edit/probe1", 23, lambda: m.bit_depth(m.mixdown(W([0.5, -0.5], O("bit_depth", 8)), 0, W(B), 1))),
        ("../sound/edit/probe1", 24, lambda: m.cd(m.mixdown(W(A), 0, W(A), 0, W(A), 0))),
        ("../sound/edit/probe1", 25, lambda: m.cd(m.mixdown(W([0.9, -0.9, 0.9]), 0, W([0.9, -0.9, 0.1]), 0))),
        ("../sound/edit/probe1", 26, lambda: m.cd(m.mixdown(W(A), "x", W(B), 0))),
        ("../sound/edit/probe1", 27, lambda: m.cd(m.mixdown(W(A), [], W(B), 0))),
        ("../sound/edit/probe1", 28, lambda: m.cd(m.mixdown(W(A), [0, "x"], W(B), 0))),
        ("../sound/edit/probe1", 29, lambda: m.cd(m.mixdown(W(A), 0, W(B), 0, W(A), 1))),
        ("../sound/edit/probe1", 30, lambda: m.cd(m.mixdown(W([0.5, 0.25], O("bit_depth", 24)), 0, W([0.1, 0.3], O("bit_depth", 24)), 0))),
        ("../sound/edit/probe1", 31, lambda: m.cd(m.stereo2mono(W([0.5, -0.25], O("bit_depth", 8))))),
        ("../sound/edit/probe1", 32, lambda: m.cd(m.stereo2mono(W([[0.5, 0.25], [0.1, 0.2], [0.3, -0.4]])))),
        ("../sound/edit/probe1", 33, lambda: m.cd(m.stereo2mono(W([[0.7, 0.25], [0.3, -0.2]], O("bit_depth", 24))))),
        ("../sound/edit/probe1", 34, lambda: m.stereo2mono(W(A), W(B))),
        ("../sound/edit/probe1", 35, lambda: m.stereo2mono("x")),
        ("../sound/edit/probe1", 36, lambda: m.cd(m.splice(W(Z), W(B)))),
        ("../sound/edit/probe1", 37, lambda: m.cd(m.splice(W(Z), W(B), 2, 1))),
        ("../sound/edit/probe1", 38, lambda: m.cd(m.splice(W(A), W(Z)))),
        ("../sound/edit/probe1", 39, lambda: m.cd(m.splice(W(S), W([B, [-0.5, -0.5, 0.5, 0.5]])))),
        ("../sound/edit/probe1", 40, lambda: m.cd(m.splice(W(S), W([B, [0.5, 0.5, 0.5, 0.5]])))),
        ("../sound/edit/probe1", 41, lambda: m.cd(m.splice(W(A), W(B, O("samplerate", 22050)), 2, 1))),
        ("../sound/edit/probe1", 42, lambda: m.cd(m.splice(W(A), W(B, O("bit_depth", 8)), 2, 1))),
        ("../sound/edit/probe1", 43, lambda: m.cd(m.splice(W(A), W(B, O("samplerate", 22050)), 2, 0))),
        ("../sound/edit/probe1", 44, lambda: m.cd(m.splice(W(A), W(B), 1.5))),
        ("../sound/edit/probe1", 45, lambda: m.cd(m.splice(W(A), W(B), 0.00004))),
        ("../sound/edit/probe1", 46, lambda: m.cd(m.splice(W(A), W(B), F(3, 2)))),
        ("../sound/edit/probe1", 47, lambda: m.cd(m.splice(W(A), W(B), "x"))),
        ("../sound/edit/probe1", 48, lambda: m.cd(m.splice(W(A), W(B), 2, "x"))),
        ("../sound/edit/probe1", 49, lambda: m.cd(m.splice(W(A), W(B), 2, F(1, 2)))),
        ("../sound/edit/probe1", 50, lambda: m.cd(m.splice(W(A), W(B), 2, 0))),
        ("../sound/edit/probe1", 51, lambda: m.cd(m.splice(W(A), W(B), 2, 0, 1))),
        ("../sound/edit/probe1", 52, lambda: m.cd(m.splice(W(A)))),
        ("../sound/edit/probe1", 53, lambda: m.cd(m.splice(W(A), [1, 2]))),
        ("../sound/edit/probe1", 54, lambda: m.cd(m.splice(W(A), W(B), -3))),
        ("../sound/edit/probe1", 55, lambda: m.cd(m.splice(W(A), W([0.9])))),
        ("../sound/edit/probe1", 56, lambda: m.cd(m.splice(W(A), W([0.9]), 3, 1))),
        ("../sound/edit/probe1", 57, lambda: m.cd(m.splice(W([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]), W([0.9, -0.9]), 6, 1))),
        ("../sound/edit/probe1", 58, lambda: m.cd(m.splice(W([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]), W([0.9, -0.9]), 7, 0))),
        ("../sound/edit/probe1", 59, lambda: m.cd(m.splice(W([0.9] * 4), W([0.9] * 4), 2, 0))),
        ("../sound/edit/probe1", 60, lambda: m.cd(m.splice(W(A, O("bit_depth", 24)), W(B, O("bit_depth", 24)), 2, 1))),
        ("../sound/edit/probe1", 61, lambda: m.cd(m.splice(W(A, O("bit_depth", 8)), W(B, O("bit_depth", 8))))),
        ("../sound/edit/probe1", 62, lambda: m.cd(m.set_channel_data(W(A), [0.9, 0.8], 3))),
        ("../sound/edit/probe1", 63, lambda: m.set_channel_data(W(A), [0.9, 0.8], 4)),
        ("../sound/edit/probe1", 64, lambda: m.set_channel_data(W(A), [0.9, 0.8], -1)),
        ("../sound/edit/probe1", 65, lambda: m.set_channel_data(W(A), [0.9, 0.8], 1.0)),
        ("../sound/edit/probe1", 66, lambda: m.set_channel_data(W(A), [0.9, 0.8], 0, 1)),
        ("../sound/edit/probe1", 67, lambda: m.set_channel_data(W(A), [0.9, 0.8], 0, "right")),
        ("../sound/edit/probe1", 68, lambda: m.cd(m.set_channel_data(W(A), [0.9, 0.8], 0, "left"))),
        ("../sound/edit/probe1", 69, lambda: m.set_channel_data(W(A), [0.9, 0.8], 0, "x")),
        ("../sound/edit/probe1", 70, lambda: m.set_channel_data(W(A), 0.9)),
        ("../sound/edit/probe1", 71, lambda: P47.text(m.set_channel_data(W(A), [[0.9, 0.8]]))),
        ("../sound/edit/probe1", 72, lambda: m.cd(m.set_channel_data(W(A), [[0.9, 0.8]], 1))),
        ("../sound/edit/probe1", 73, lambda: m.set_channel_data(W(S), [[0.9, 0.8], [0.1, 0.2]], 0, 0)),
        ("../sound/edit/probe1", 74, lambda: m.cd(m.set_channel_data(W(S), [[0.9, 0.8], [0.1, 0.2]], 1))),
        ("../sound/edit/probe1", 75, lambda: m.cd(m.set_channel_data(W(S), [0.9, 0.8]))),
        ("../sound/edit/probe1", 76, lambda: m.set_channel_data(W(A), [0.9, "x"])),
        ("../sound/edit/probe1", 77, lambda: m.cd(m.set_channel_data(W(A), [0.9, "x"], 3))),
        ("../sound/edit/probe1", 78, lambda: m.set_channel_data(W(A), [[0.9], [0.8, 0.7]])),
        ("../sound/edit/probe1", 79, lambda: m.cd(m.set_channel_data(W(A), [sqrt_half(m.P), F(1, 3), -2, pi_(m.P)]))),
        ("../sound/edit/probe1", 80, lambda: m.cd(m.set_channel_data(W(A), []))),
        ("../sound/edit/probe1", 81, lambda: m.cd(m.set_channel_data(W(A), [0.9], 0, 0, 7))),
        ("../sound/edit/probe1", 82, lambda: m.set_channel_data(W(A))),
        ("../sound/edit/probe1", 83, lambda: m.set_channel_data([1, 2], [0.9])),
        ("../sound/edit/probe1", 84, lambda: m.cd(m.set_channel_data(W([0.5, -0.5], O("bit_depth", 8)), [0.9, -1]))),
        ("../sound/edit/probe1", 85, lambda: m.cd(m.set_channel_data(W([0.5, -0.5], O("bit_depth", 24)), [0.9, -1]))),
        ("../sound/edit/probe1", 86, lambda: m.cd(m.set_channel_data(W(A), [1, 1]))),
    ]


def sqrt_half(P):
    import mpmath
    return P.Irr(mpmath.sqrt(2) / 2, P.trunc5(math.sqrt(2) / 2))


def pi_(P):
    import mpmath
    return P.Irr(mpmath.pi, P.trunc5(math.pi))


NUM = re.compile(r"-?\d+(?:\.\d*)?(?:e[-+]?\d+)?")


def flatnums(x):
    if isinstance(x, list):
        for y in x:
            yield from flatnums(y)
    else:
        yield x


def render(x):
    if isinstance(x, list):
        return "[" + ",".join(render(y) for y in x) + "]"
    if isinstance(x, float):
        return P47.gprint(x)
    return str(x)


def giac_row(f, ln):
    line = open(GIAC + f + ".out").read().split("\n")[ln - 1]
    return line.split("  =>  ", 1)


def judge(fn, gout):
    try:
        ours = fn()
    except (P47.ClipError, P64.ClipError) as e:
        m = str(e)
        if m == GENERIC:
            return "ERROR: generic", gout.startswith("GIAC_ERROR: ") and "Error: " in gout and ". Error" not in gout.split("|")[0]
        return "ERROR: " + m, gout.startswith("GIAC_ERROR: " + m + ".")
    if isinstance(ours, str):
        return ours, ours == gout
    r = render(ours)
    gn = [float(t) for t in NUM.findall(gout)]
    on = list(flatnums(ours)) if isinstance(ours, list) else [ours]
    good = NUM.sub("#", r) == NUM.sub("#", gout) and len(gn) == len(on) and all(
        abs(float(o) - g) <= 1e-11 * max(1, abs(g)) for o, g in zip(on, gn))
    return r, good


def main():
    models = [("47", Model(P47)), ("f64", Model(P64)), ("ours", Model(P64, strict=True))]
    lines, count = [], {k: 0 for k, _ in models}
    rows = {k: R(m) for k, m in models}
    n = len(rows["47"])
    for j in range(n):
        f, ln, _ = rows["47"][j]
        inp, gout = giac_row(f, ln)
        res = {}
        for k, m in models:
            m.warned = False
            r, good = judge(rows[k][j][2], gout)
            res[k] = (r, good, m.warned)
            count[k] += good
        flags = " ".join(f"{k}={'PASS' if g else 'FAIL'}{'+warn' if w else ''}" for k, (r, g, w) in res.items())
        lines.append(f"- {f}:{ln} `{inp[:160]}`\n    {flags}\n    giac: {gout[:200]}\n    47: {res['47'][0][:200]}")
        if res["ours"][0] != res["47"][0]:
            lines.append(f"    ours: {res['ours'][0][:200]}")
    with open(HERE + "table.txt", "w") as t:
        t.write("\n".join(lines) + "\n")
        t.write(f"PASS 47-bit {count['47']}/{n}, f64 {count['f64']}/{n}, ours {count['ours']}/{n}\n")
    print(f"PASS 47-bit {count['47']}/{n}, f64 {count['f64']}/{n}, ours {count['ours']}/{n} (table.txt)")


if __name__ == "__main__":
    main()
