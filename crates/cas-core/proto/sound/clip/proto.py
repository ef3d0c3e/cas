"""T7.c1 scratch (origin: phase 1.7, T7.c1): Python prototype of giac's audio clip (createwav, channel_data, soundsec, ...).

Two paths per clip:
  giac path : f64 like giac, with giac's gen doubles (low 5 mantissa bits zeroed, see trunc5)
  exact path: Fractions / mpmath, no f64; codes compared with the giac path, differences flagged.
"""
import math, os, re, struct, sys
from fractions import Fraction as F
import mpmath
mpmath.mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
GIAC = HERE + "/../../giac/"

class ClipError(Exception):
    pass

def trunc5(x):
    """giac gen double: low 5 mantissa bits cleared (probe: (1.0+2.0^-48)-1.0 => 0.0, floor(0.001*8000) => 7)."""
    if not math.isfinite(x):
        return x
    b = struct.unpack("<Q", struct.pack("<d", x))[0] & ~0x1F
    return struct.unpack("<d", struct.pack("<Q", b))[0]

class Irr:
    """irrational datum: exact value (mpmath) and giac's own float evaluation."""
    def __init__(self, exact, giac_f):
        self.exact, self.giac_f = exact, giac_f

class Opt:
    def __init__(self, name, value):
        self.name, self.value = name, value

def is_num(x):
    return isinstance(x, (int, float, F, Irr)) and not isinstance(x, bool)

def gd(x):  # evalf_double into a gen
    return x.giac_f if isinstance(x, Irr) else trunc5(float(x))

def ex(x):  # exact value
    if isinstance(x, Irr):
        return x.exact
    return F(x) if isinstance(x, (int, float)) else x

def gprint(x):
    if isinstance(x, list):
        return "[" + ",".join(gprint(y) for y in x) + "]"
    if isinstance(x, float):
        s = format(x, ".12g")
        return s if ("." in s or "e" in s or "inf" in s) else s + ".0"
    return str(x)

def ckmatrix(d):
    return len(d) > 0 and all(isinstance(r, list) for r in d) and len({len(r) for r in d}) == 1 and len(d[0]) > 0

def flat(d):
    for y in d:
        if isinstance(y, list):
            yield from flat(y)
        else:
            yield y

def absmax(vals):
    m = None
    for y in vals:
        if not is_num(y):
            return None
        a = abs(mpmath.mpf(ex(y))) if isinstance(y, Irr) else abs(ex(y))
        m = a if m is None or a > m else m
    return m

def encode_f64(depth, s):
    rng = (1 << (depth - 1)) - 1
    if math.isnan(s):   # giac: ratio inf (peak 0), 0*inf = NaN, (int)floor(NaN) -> INT_MIN, low bytes 0
        return 0
    v = max(-1.0, min(1.0, s)) if (s < -1.0 or s > 1.0) else s
    return int(math.floor(rng * v + 0.5))

def encode_exact(depth, s):
    rng = (1 << (depth - 1)) - 1
    if s is None:
        return 0
    v = max(-1, min(1, s))
    return int(mpmath.floor(rng * v + mpmath.mpf(1) / 2)) if isinstance(v, mpmath.mpf) else math.floor(rng * v + F(1, 2))

def createwav(*args):
    if len(args) == 1 and isinstance(args[0], int):
        if args[0] < 1:
            raise ClipError("Expected a length in samples")
        args, data, length = (), None, args[0]
    elif len(args) == 1 and isinstance(args[0], list):
        args, data, length = (), args[0], 0
    else:
        data, length = None, 0
    bd, sr, nc, norm, secs = 16, 44100, 0, 1.0, -1.0
    norm_x, secs_x = F(1), None
    for i, it in enumerate(args, 1):
        if isinstance(it, Opt):
            n, rh = it.name, it.value
            isint = isinstance(rh, int)
            if n == "channels":
                if not isint or rh < 1:
                    raise ClipError("Number of channels should be a positive integer")
                nc = rh
            elif n == "bit_depth":
                if not isint or rh not in (8, 16, 24):
                    raise ClipError("Bit depth should be either 8, 16 or 24")
                bd = rh
            elif n == "samplerate":
                if not isint or rh < 1:
                    raise ClipError("Sample rate should be a positive integer")
                sr = rh
            elif n == "length":
                if not isint or rh < 1:
                    raise ClipError("Length should be a positive integer")
                length = rh
            elif n == "normalize":
                if not is_num(rh):
                    raise ClipError("Normalization level should be a real number")
                norm, norm_x = -gd(rh), -ex(rh)
            elif n == "duration":
                if not is_num(rh):
                    raise ClipError("Duration should be a real number")
                secs, secs_x = gd(rh), ex(rh)
            elif n == "channel_data":
                if not isinstance(rh, list):
                    raise ClipError("Channel data should be either a matrix, vector, or string")
                data = rh
            else:
                raise ClipError(f"Invalid argument {i}")
        elif isinstance(it, list) and data is None:
            data = it
        else:
            raise ClipError(f"Invalid argument {i}")
    len_exact = length
    if secs >= 0 and length == 0:
        length = 1 + math.floor(secs * sr)
        len_exact = 1 + math.floor(secs_x * sr)
    is_mat = False
    if data is not None:
        is_mat = ckmatrix(data)
        if length == 0:
            length = len_exact = len(data[0]) if is_mat else len(data)
        if nc == 0:
            nc = len(data) if is_mat else 1
        elif is_mat and len(data) != nc:
            raise ClipError("Number of channels does not match the number of data rows")
    if nc == 0:
        raise ClipError("Number of channels not set")
    if length == 0:
        raise ClipError("Audio length not set")
    clip = dict(rate=sr, depth=bd, channels=[[0] * length for _ in range(nc)],
                exact=[[0] * len_exact for _ in range(nc)], len_exact=len_exact)
    if data is not None:
        ratio, ratio_x = 1.0, F(1)
        if norm > 0:
            m = absmax(list(flat(data)))
            if m is None:
                raise ClipError("Sample data should be numeric")
            md = trunc5(float(m))
            ratio = 10.0 ** (-norm / 20.0) / md if md != 0 else math.inf
            ratio_x = None if m == 0 else mpmath.power(10, mpmath.mpf(-norm_x) / 20) / mpmath.mpf(m)
        for c in range(nc):
            row = data[c] if is_mat else data
            for key, L in (("channels", length), ("exact", len_exact)):
                for p, y in enumerate(row[:L]):
                    if not is_num(y):
                        raise ClipError(f"Failed to set channel data, {gprint(y)} is not a valid sample value")
                    if key == "channels":
                        clip[key][c][p] = encode_f64(bd, ratio * gd(y))
                    else:
                        if ratio_x is None:
                            s = None
                        elif isinstance(ratio_x, F) and not isinstance(y, Irr):
                            s = ratio_x * ex(y)
                        else:
                            s = mpmath.mpf(ratio_x) * mpmath.mpf(ex(y))
                        clip[key][c][p] = encode_exact(bd, s)
    if clip["exact"] != clip["channels"]:
        DIFFS.append((clip["channels"], clip["exact"]))
    return clip

DIFFS = []

def is_clip(c):
    return isinstance(c, dict) and "rate" in c

def need_clip(c):
    if not is_clip(c):
        raise ClipError("Invalid sound data")

def decode(clip, k):
    return F(k, 1 << (clip["depth"] - 1))

def samplerate(c): need_clip(c); return c["rate"]
def bit_depth(c):
    if not is_clip(c): raise ClipError("Expected audio clip or image")
    return c["depth"]
def channels(c):
    if not is_clip(c): raise ClipError("Expected audio clip or image")
    return len(c["channels"])
def length(c): return len(c["channels"][0])
def duration(c): need_clip(c); return F(length(c), c["rate"])
def text(c):
    return "a sound clip with %d samples at %d Hz (%d bit, %s)" % (
        length(c), c["rate"], c["depth"], "mono" if len(c["channels"]) == 1 else "stereo")
def peak(c):
    return max([F(0)] + [decode(c, k) for ch in c["channels"] for k in ch])
def peak_dbfs(c):
    p = peak(c)
    return -math.inf if p == 0 else 20 * math.log10(p)
def about(c):
    nc = len(c["channels"])
    chd = {1: " (mono)", 2: " (stereo)"}.get(nc, "")
    log = [f"Bit depth: {c['depth']}", f"Channels: {nc}{chd}", f"Sample rate: {c['rate']}",
           f"Length: {length(c)} samples", "Duration: %g seconds" % (length(c) / c["rate"]),
           "Level: %g dB" % peak_dbfs(c)]
    return text(c), log

def index(c, chn):  # clip[chn] (operator[], one channel, whole length)
    if chn >= len(c["channels"]) or chn < 0:
        raise ClipError("Invalid channel index")
    return [decode(c, k) for k in c["channels"][chn]]

class Seq(list): pass        # giac sequence answer (a,b)

def channel_data(c, *opts):
    if not is_clip(c):
        raise ClipError("Expected audio clip or image")
    L, sr = length(c), c["rate"]
    chan, start, slen, asm = -1, 0, L, False
    for i, o in enumerate(opts, 2):
        if isinstance(o, int): chan = o
        elif o == "left": chan = 0
        elif o == "right": chan = 1
        elif o == "matrix": asm = True
        elif isinstance(o, Opt):
            if o.name == "range":
                r = o.value
                if isinstance(r, list):
                    if len(r) != 2: raise ClipError("Invalid range")
                    if not all(isinstance(y, int) for y in r): raise ClipError("Invalid range")
                    start, slen = r
                elif isinstance(r, tuple):     # a..b
                    a, b = r
                    start = math.floor(0.5 + sr * gd(a))
                    slen = math.floor(0.5 + sr * trunc5(float(ex(b) - ex(a))))
                    assert start == math.floor(F(1, 2) + sr * ex(a)) and slen == math.floor(F(1, 2) + sr * (ex(b) - ex(a)))
                else: raise ClipError("Invalid range")
        else:
            raise ClipError(f"Invalid argument {i}")
    slen = min(slen, L - start)
    if start < 0 or start >= L or slen <= 0:
        raise ClipError("Invalid slice bounds")
    get = lambda ch: [decode(c, k) for k in c["channels"][ch][start:start + slen]]
    if chan < 0:
        rows = [get(ch) for ch in range(len(c["channels"]))]
        if len(rows) == 1: return rows[0]
        return rows if asm else Seq(rows)
    if chan >= len(c["channels"]):
        raise ClipError("Invalid channel number")
    return get(chan)

def soundsec(n, N=44100):
    if not is_num(n) or gd(n) <= 0 or not isinstance(N, int) or N < 1:
        raise ClipError("Invalid sound parameters")
    nr = float(N) * gd(n)
    cnt_giac = sum(1 for _ in range(math.ceil(nr) + 1) if _ < nr)
    cnt = math.ceil(ex(n) * N)          # #k with k < n*N
    assert cnt == cnt_giac, (n, N, cnt, cnt_giac)
    return [F(k, N) for k in range(cnt)]

# ---------------------------------------------------------------- rows
def sin_data():
    out = []
    t880pi = trunc5(880 * math.pi)
    for k in range(8):
        t = trunc5(k / 8000.0)
        out.append(Irr(mpmath.sin(mpmath.mpf(11) * k / 100 * mpmath.pi), trunc5(math.sin(trunc5(t880pi * t)))))
    return out

W = createwav
A = [0.5, -0.5, 0.25, 0]
M = [[0.5, -0.5, 0.25, 0], [0.1, 0.2, 0.3, 0.4]]
S6 = [0.5, -0.5, 0.25, 0, 1, 0.75]
O = Opt
ROWS = [
 ("probe5", 1, lambda: text(W(A))),
 ("probe5", 3, lambda: '"' + text(W(A)) + '"'),
 ("probe5", 4, lambda: length(W(A))),
 ("probe5", 5, lambda: index(W(A), 0)),
 ("probe5", 6, lambda: about(W(A))[0]),
 ("probe5", 7, lambda: channels(W(A))),
 ("probe5", 8, lambda: bit_depth(W(A))),
 ("probe5", 9, lambda: samplerate(W(A))),
 ("probe5", 10, lambda: duration(W(A))),
 ("probe5", 11, lambda: channel_data(W(A))),
 ("probe5", 12, lambda: channel_data(W(M))),
 ("probe5", 13, lambda: channel_data(W(M), "matrix")),
 ("probe5", 14, lambda: channel_data(W(M), 2)),
 ("probe5", 15, lambda: channel_data(W(M), "right")),
 ("probe5", 16, lambda: channel_data(W(M), "left")),
 ("probe5", 17, lambda: channel_data(W(S6, O("samplerate", 2)), O("range", (1, 2)))),
 ("probe5", 18, lambda: channel_data(W(S6, O("samplerate", 2)), O("range", (1.0, 2.5)))),
 ("probe5", 19, lambda: channel_data(W(S6), O("range", [1, 3]))),
 ("probe5", 20, lambda: channels(W(M))),
 ("probe5", 21, lambda: text(W([F(1, 2), F(-1, 2)]))),
 ("probe5", 22, lambda: channel_data(W([F(1, 2), F(-1, 2), F(1, 3)]))),
 ("probe5", 23, lambda: channel_data(W(A, O("bit_depth", 8)))),
 ("probe5", 24, lambda: channel_data(W(A, O("bit_depth", 24)))),
 ("probe5", 25, lambda: channel_data(W(A, O("bit_depth", 32)))),
 ("probe5", 26, lambda: channel_data(W([0.123456789, -0.5]))),
 ("probe5", 27, lambda: channel_data(W([2, -1, 0.5]))),
 ("probe5", 28, lambda: channel_data(W([2, -1, 0.5], O("normalize", -3)))),
 ("probe5", 29, lambda: channel_data(W([0.2, -0.1, 0.05], O("normalize", 0)))),
 ("probe5", 30, lambda: channel_data(W(O("duration", 0.0001)))),
 ("probe5", 31, lambda: channel_data(W(O("duration", 0.0001), O("samplerate", 20000)))),
 ("probe5", 32, lambda: channel_data(W(O("channels", 2), O("size", 3)))),
 ("probe5", 33, lambda: channel_data(W(O("size", 3)))),
 ("probe5", 34, lambda: samplerate(W([0.1], O("samplerate", 48000)))),
 ("probe5", 35, lambda: duration(W([0.1, 0.2, 0.3], O("samplerate", 2)))),
 ("probe5", 36, lambda: text(W([0.1, 0.2], O("channels", 2)))),
 ("probe5", 37, lambda: channel_data(W([0.1, 0.2], O("channels", 2)))),
 ("probe5", 62, lambda: len(soundsec(1))),
 ("probe5", 63, lambda: soundsec(1)[0:4]),
 ("probe5", 64, lambda: len(soundsec(1, 22100))),
 ("probe5", 65, lambda: soundsec(1, 22100)[0:4]),
 ("probe5", 66, lambda: soundsec(F(1, 1000), 8000)),
 ("probe5", 67, lambda: soundsec(0.5, 4)),
 ("probe5", 68, lambda: soundsec(2, 3)),
 ("probe9", 25, lambda: channel_data(W([0.5, -0.5], O("normalize", -1)))),
 ("probe9", 26, lambda: channel_data(W([0.5, -0.5], O("normalize", 1)))),
 ("probe9", 27, lambda: channel_data(W([0.5, -0.5], O("normalize", -100)))),
 ("probe9", 28, lambda: channel_data(W([0.2, -0.1], O("normalize", 0)))),
 ("probe9", 29, lambda: channel_data(W([1.5, -0.1], O("normalize", 0)))),
 ("probe9", 30, lambda: channel_data(W([0, 0, 0]))),
 ("probe9", 31, lambda: channel_data(W(["x", 1]))),
 ("probe9", 32, lambda: channel_data(W([0.5, 0.25], O("channels", 3)))),
 ("probe9", 33, lambda: text(W([0.5, 0.25], O("channels", 3)))),
 ("probe9", 34, lambda: text(W([[0.5, 0.25], [0.1, 0.2], [0.3, 0.4]]))),
 ("probe9", 35, lambda: text(W([[0.5, 0.25], [0.1]]))),
 ("probe9", 36, lambda: channel_data(W([[0.5, 0.25], [0.1]]), "matrix")),
 ("probe9", 37, lambda: text(W([0.5, 0.25], O("samplerate", 0)))),
 ("probe9", 38, lambda: text(W([0.5, 0.25], O("samplerate", 1.5)))),
 ("probe9", 39, lambda: text(W([], O("samplerate", 8000)))),
 ("probe9", 40, lambda: text(W(O("duration", 0.001), O("channels", 1), O("samplerate", 8000)))),
 ("probe9", 41, lambda: channel_data(W(O("duration", 0.001), O("channels", 1), O("samplerate", 8000)))),
 ("probe9", 42, lambda: channel_data(W(O("size", 3), O("channels", 1)))),
 ("probe9", 43, lambda: channel_data(W([0.5, 0.25], O("duration", 0.0001)))),
 ("probe9", 44, lambda: channel_data(W(sin_data()))),
 ("probe9", 47, lambda: samplerate([1, 2])),
 ("probe9", 48, lambda: channel_data(W([0.5, 0.25]), 0)),
 ("probe9", 49, lambda: channel_data(W([0.5, 0.25]), 1)),
 ("probe9", 50, lambda: channel_data(W([[0.5, 0.25], [0.1, 0.2]]), 1)),
 ("probe9", 51, lambda: channel_data(W([[0.5, 0.25], [0.1, 0.2]]), 0)),
 ("probe9", 52, lambda: channel_data(W(M), "right", O("range", (1, 2)))),
 ("probe9", 53, lambda: channel_data(W(S6), O("range", (4, 10)))),
 ("probe9", 54, lambda: channel_data(W(S6), O("range", [4, 10]))),
 ("probe9", 55, lambda: channel_data(W(S6, O("samplerate", 2)), O("range", (0.5, 1.0)))),
 ("probe10", 21, lambda: length(W(M))),
 ("probe10", 22, lambda: index(W(M), 1)),
 ("probe10", 23, lambda: index(W(M), 0)[2]),
 ("probe10", 27, lambda: text(W(W([0.5, -0.5])))),
 ("probe10", 28, lambda: channel_data(W(W([0.5, -0.5]), O("bit_depth", 8)))),
 ("probe11", 34, lambda: about(W([0.5, 0.25], O("channels", 2)))[0]),
 # sound/clip/dur.cas (caseval, this task): the duration length rule under giac's 47-bit gen doubles
 ("../sound/clip/dur", 4, lambda: len(channel_data(W(O("duration", 0.002), O("channels", 1), O("samplerate", 8000))))),
 ("../sound/clip/dur", 5, lambda: len(channel_data(W(O("duration", F(1, 1000)), O("channels", 1), O("samplerate", 8000))))),
 ("../sound/clip/dur", 6, lambda: len(channel_data(W(O("duration", 0.5), O("channels", 1), O("samplerate", 4))))),
 ("../sound/clip/dur", 7, lambda: len(channel_data(W(O("duration", 1), O("channels", 1), O("samplerate", 10))))),
 ("../sound/clip/dur", 8, lambda: len(channel_data(W(O("duration", 0.25), O("channels", 1), O("samplerate", 10))))),
 ("../sound/clip/dur", 9, lambda: len(channel_data(W(O("length", 3), O("duration", 1), O("channels", 1), O("samplerate", 10))))),
 ("../sound/clip/dur", 10, lambda: len(channel_data(W([1, 2], O("duration", 1), O("samplerate", 10))))),
]
# rows of the grep not run here (other natives/operators, later T7 coders): probe5 2 (type), 38-39 (stereo2mono),
# 40-43 (set_channel_data), 44-47 (splice), 48-52 (mixdown), 53-57 (resample), 58-61 (+, unary -, *);
# probe9 1-24 (mixdown, splice, +, ==, -), 45-46 (stereo2mono); probe10 20, 24-26 (concat, -, abs, evalf).

NUM = re.compile(r"-?\d+(?:\.\d*)?(?:e[-+]?\d+)?")

def render(x, exact=True):
    if isinstance(x, Seq):
        return ",".join(render(r) for r in x)
    if isinstance(x, list):
        return "[" + ",".join(render(r) for r in x) + "]"
    return str(x)

def flatnums(x):
    if isinstance(x, list):
        for y in x: yield from flatnums(y)
    else:
        yield x

def giac_row(f, ln):
    line = open(GIAC + f + ".out").read().split("\n")[ln - 1]
    inp, out = line.split("  =>  ", 1)
    return inp, out

def main():
    ok = bad = 0
    table = []
    for f, ln, fn in ROWS:
        inp, gout = giac_row(f, ln)
        DIFFS.clear()
        try:
            ours = fn()
            r = render(ours) if not isinstance(ours, str) else ours
            if isinstance(ours, str):
                good = ours == gout
            else:
                shape_o = NUM.sub("#", re.sub(r"(-?\d+)/(\d+)", "#", r))
                shape_g = NUM.sub("#", gout)
                gn = [float(t) for t in NUM.findall(gout)]
                on = list(flatnums(ours)) if isinstance(ours, list) else [ours]
                good = shape_o == shape_g and len(gn) == len(on) and all(
                    abs(float(o) - g) <= 1e-11 * max(1, abs(g)) for o, g in zip(on, gn))
        except ClipError as e:
            r = "ERROR: " + str(e)
            good = gout.startswith("GIAC_ERROR: " + str(e) + ".")
        ok += good; bad += not good
        if DIFFS:
            r += "  [EXACT PATH DIFFERS: " + "; ".join(f"giac {a} exact {b}" for a, b in DIFFS)[:400] + "]"
        table.append((f, ln, inp, r, good))
        print(("OK  " if good else "BAD ") + f"{f}:{ln}  {inp}\n      ours: {r[:300]}\n      giac: {gout[:300]}")
    print(f"matched {ok}/{ok + bad}")
    with open(HERE + "/table.txt", "w") as t:
        for f, ln, inp, r, good in table:
            t.write(f"- {f}:{ln} `{inp}` => `{r}`{'' if good else '  MISMATCH'}\n")

if __name__ == "__main__":
    main()
