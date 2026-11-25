"""Python port of giac markup.cc gen2markup: LaTeX, MathML presentation and content arms
(giac 2.0.0.21, markup.cc (source: C0009)). Extends 1.6/T4's latex-only port (../../latex/proto/markup.py).

Each function cites the markup.cc lines it ports. Scheme and unit arms are dropped.
Gen tree: Int, Bool, Frac, Dbl, Cplx, Idnt, Str, Vect, Symb (giac-shaped, see parse.py).
"""
import re
from dataclasses import dataclass

# ---- tree -------------------------------------------------------------------------------------


@dataclass
class Int:
    v: int


@dataclass
class Bool:
    v: bool


@dataclass
class Frac:
    num: object
    den: object


@dataclass
class Dbl:
    text: str
    neg: bool = False

    @property
    def val(self):
        return -float(self.text) if self.neg else float(self.text)


@dataclass
class Cplx:
    re: object
    im: object


@dataclass
class Idnt:
    name: str


@dataclass
class Str:
    s: str


@dataclass
class Func:            # giac _FUNC: a bare builtin ('Gamma')
    name: str


@dataclass
class Vect:
    items: list
    sub: str = 'list'


@dataclass
class Symb:
    op: str
    arg: object


I = Cplx(Int(0), Int(1))

# ---- constants :35-117 ------------------------------------------------------------------------

MML_ITIMES = "<mo>&it;</mo>"          # :35
MML_CDOT = "<mo>&middot;</mo>"        # :36
MML_DOT = "<mo>.</mo>"
MML_TIMES = "<mo>&times;</mo>"
MML_PLUS = "<mo>+</mo>"
MML_MINUS = "<mo>&minus;</mo>"        # :40
MML_APPLY = "<mo>&af;</mo>"           # :41
MML_ICOMMA = "<mo>&ic;</mo>"
MML_I = "<mi mathvariant='normal'>i</mi>"
MML_E = "<mi mathvariant='normal'>e</mi>"
MML_D = "<mi mathvariant='normal'>d</mi>"

TEX_ITIMES = " "
TEX_MATHI = "\\mathrm{i}"
TEX_MATHE = "\\mathrm{e}"
TEX_MATHD = "\\mathrm{d}"

P_APPLY, P_UNARY, P_EXP, P_MUL, P_ADD, P_SET, P_INEQ, P_EQ, P_COMP = 1, 2, 3, 4, 5, 6, 7, 8, 9
P_NOT, P_AND, P_XOR, P_OR, P_COND, P_ASSGN, P_OTHER = 10, 11, 12, 13, 14, 15, 100

T_GENERAL, T_NUM_EXACT, T_NUM_APPROX, T_SUBTYPE_IDNT, T_FRACTION = 0, 1, 2, 4, 8
T_ROOT, T_NEG, T_SUMPROD, T_FACTORIAL, T_MATRIX, T_ELEMAPP = 16, 32, 64, 128, 256, 512
T_FUNC, T_POWER, T_DERIVATIVE, T_HAS_SUBSCRIPT, T_LEADING_DIGIT, T_IDNT_NAME = (
    1024, 2048, 4096, 8192, 16384, 32768)

F_LATEX, F_TOPLEVEL, F_PRES, F_CONTENT, F_ELEMPOW, F_FACTOR = 1, 2, 8, 16, 32, 64

HEADER = "mode='display' xmlns='http://www.w3.org/1998/Math/MathML'"  # :3745


@dataclass
class MB:              # MarkupBlock :179-192 (no scheme)
    priority: int = 0
    type: int = T_GENERAL
    neg: bool = False
    latex: str = ""
    markup: str = ""
    content: str = ""
    split_pos: int = -1
    split_pos_tex: int = 0
    appl: bool = False

    def ctype(self, t):
        return (self.type & t) != 0


class Idc:             # giac's `int &idc`, one counter per export
    def __init__(self):
        self.v = 0

    def inc(self):
        self.v += 1
        return self.v

# ---- small giac helpers -----------------------------------------------------------------------

_ATOF = re.compile(r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[+-]?(?:inf(?:inity)?|nan))",
                   re.IGNORECASE)


def atof(s):
    m = _ATOF.match(s)
    return float(m.group(1)) if m else 0.0


def my_isalpha(c):
    return c.isascii() and c.isalpha() or ord(c) > 127


def is_greek_letter(s):  # :141
    return s in {"mu", "nu", "pi", "Pi", "xi", "Xi", "chi", "phi", "Phi", "eta", "rho", "tau", "psi",
                 "Psi", "beta", "zeta", "alpha", "delta", "Delta", "gamma", "Gamma", "kappa", "theta",
                 "Theta", "sigma", "Sigma", "Omega", "omega", "lambda", "Lambda", "epsilon",
                 "varepsilon"}


def is_double_letter(s):  # :167
    return len(s) == 2 and my_isalpha(s[0]) and s[0] == s[1]


def is_digits(s):  # :171
    return all(c.isdigit() for c in s)


def is_integer(g):
    return isinstance(g, Int)


def numval(g):
    if isinstance(g, Int):
        return g.v
    if isinstance(g, Dbl):
        return g.val
    if isinstance(g, Frac) and isinstance(g.num, Int) and isinstance(g.den, Int):
        from fractions import Fraction
        return Fraction(g.num.v, g.den.v)
    return None


def is_one(g):
    return numval(g) == 1


def is_minus_one(g):
    return numval(g) == -1


def is_zero(g):
    return numval(g) == 0


def is_positive(g):
    v = numval(g)
    return v is None or v >= 0


def is_sym(g, *ops):
    return isinstance(g, Symb) and g.op in ops


def is_inf(g):
    return isinstance(g, Idnt) and g.name in ("infinity", "unsigned_inf")

# ---- MathML tags :241-302 ---------------------------------------------------------------------

CONTENT_ELEMENTS = {"cn", "ci", "cs", "csymbol", "cerror", "share", "bind", "bvar", "semantics", "apply",
                    "domainofapplication", "interval", "condition", "lowlimit", "uplimit", "degree",
                    "momentabout", "logbase", "inverse", "lambda", "piecewise", "piece", "otherwise",
                    "set", "list", "vector", "matrix", "matrixrow", "declare", "reln", "fn"}


def mml_tag(tag, s, idc=0, *attrs):  # :248
    r = "id" if (tag[0] == 'c' or tag in CONTENT_ELEMENTS) else "xref"
    a = ""
    for k, v in zip(attrs[0::2], attrs[1::2]):
        if k and v:
            a += " " + k + "='" + v + "'"
    if idc == 0:
        return "<" + tag + a + ">" + s + "</" + tag + ">"
    return "<" + tag + " " + r + "='id" + str(idc) + "'" + a + ">" + s + "</" + tag + ">"


def mml_csymbol(s, cd):  # :272
    return mml_tag("csymbol", s, 0, "cd" if cd else "", cd)


def extract_id(ml, content):  # :277
    s = ml.content if content else ml.markup
    attr = " id='" if content else " xref='"
    gt = s.find(">")
    st = s.find(attr)
    if st < 0 or st > gt:
        return ""
    end = s.find("'", st + len(attr))
    ret = s[st + len(attr):end]
    s = s[:st] + s[end + 1:]
    if content:
        ml.content = s
    else:
        ml.markup = s
    return ret


def insert_id(s, idc, content):  # :293
    if idc == 0:
        return s
    pos = s.find(">")
    if pos <= 0:
        return s
    if s[pos - 1] == '/':
        pos -= 1
    return s[:pos] + (" id='id" if content else " xref='id") + str(idc) + "'" + s[pos:]

# ---- strings :304-468 -------------------------------------------------------------------------


def trim_string(s):  # :304
    indent, i = 0, 0
    while i < len(s):
        c = s[i]
        if c == ' ':
            indent += 1
        elif c not in '\n\t':
            break
        i += 1
    j = len(s) - 1
    while j >= 0 and s[j] in ' \n\t':
        j -= 1
    return s[i:j + 1], indent


def str_to_tex(s_orig, quote):  # :333
    s, _ = trim_string(s_orig)
    ret = "".join("\\" + c if c in "\\%_${}^" else c for c in s)
    return "``" + ret + "''" if quote else ret


def str_to_mml(s_orig, ind=False, tag="", idc=0):  # :382
    s, indent = trim_string(s_orig)
    ret = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    if tag:
        ret = mml_tag(tag, ret, idc)
    if ind and indent > 0:
        ret = mml_tag("mphantom", mml_tag("mtext", "x" * indent)) + ret
    return ret


def string2markup(ml, s_orig, flags, idc):  # :408 (quote set, not code, not error)
    lines = []
    last = 0
    while True:
        pos = s_orig.find("\n", last)
        if pos < 0:
            break
        s = s_orig[last:pos]
        if s.strip(' ') != "":
            lines.append(s)
        last = pos + 1
    lines.append(s_orig[last:])
    if flags & F_LATEX:
        ml.latex = "\\\\".join("\\text{" + str_to_tex(l, True) + "}" for l in lines)
        if len(lines) > 1:
            ml.latex = "\\begin{array}{l}" + ml.latex + "\\end{array}"
    if flags & F_CONTENT:
        ml.content = str_to_mml(s_orig, False, "cs", idc.inc())
    if flags & F_PRES:
        if len(lines) == 1:
            ml.markup = mml_tag("ms", str_to_mml(lines[0]), idc.v)
        else:
            ml.markup = mml_tag("mtable", "".join("<mtr><mtd>" + str_to_mml(l, False, "ms") + "</mtd></mtr>"
                                                  for l in lines), idc.v, "columnalign", "left")

# ---- identifiers :621-934 (non-unit) ----------------------------------------------------------


def idnt2markup(s_orig, tex, idc=0):
    """:621 -> (text, isname); tex=False is the MathML arm (typ 2/3)."""
    if s_orig.endswith("_") and is_greek_letter(s_orig[:-1]):  # :845
        s = s_orig[:-1]
        if s in ("phi", "epsilon"):
            s = "var" + s
        return ("\\" + s + "\\_" if tex else "&" + s + ";_"), False
    i = len(s_orig) - 1
    while i >= 0 and s_orig[i].isdigit():
        i -= 1
    s = s_orig[:i + 1]
    ssub = ""
    if i < len(s_orig) - 1 and (s == "log" or len(s) == 1 or is_greek_letter(s) or is_double_letter(s)):
        ssub = s_orig[i + 1:]
    if ssub == "":
        pos = s_orig.find("_")
        if pos >= 0:
            s = s_orig[:pos]
            au = s_orig[pos + 1:]
            if len(s) == 1 or is_greek_letter(s) or is_double_letter(s) or is_digits(au):
                ssub = au
    if ssub == "":
        s = s_orig
    ln = len(s)
    mdf = ""
    if is_greek_letter(s):
        if tex and s in ("phi", "epsilon"):
            s = "var" + s
        s = "\\" + s if tex else "&" + s + ";"
        ln = 1
    len_sub = len(ssub)
    if ssub and is_greek_letter(ssub):
        if tex and ssub in ("phi", "epsilon"):
            ssub = "var" + ssub
        ssub = "\\" + ssub if tex else "&" + ssub + ";"
        len_sub = 1
    if is_double_letter(s):
        s = s[0]
        if tex:
            s = "\\mathbf{" + s + "}"
            ln = 1
        else:
            mdf = "mathvariant"
    if ln > 3 and s[0] == '`' and s[1] == ' ' and s[-1] == '`':
        s = s[2:2 + ln - 3]
    sid = idc if ssub == "" else 0
    if tex:
        ret = s if ln == 1 else "\\mathrm{" + s + "}"
    else:
        ret = mml_tag("mi", s, sid) if ln == 1 else mml_tag("mi", s, sid, mdf, "bold")
    if ssub:
        tokens = ssub.split("_")
        multi = all(t != "" and (t == "0" or atof(t) != 0) for t in tokens)
        if multi:
            ssub = ("," if tex else "<mo>,<mo>").join(tokens)
        if tex:
            ret += "_{" + (ssub if (len_sub == 1 or atof(ssub) != 0) else "\\mathrm{" + ssub + "}") + "}"
        else:
            ret = mml_tag("msub", ret + ("<mrow>" + ssub + "</mrow>" if multi else
                                         ("<mn>" + ssub + "</mn>" if (ssub == "0" or atof(ssub) != 0)
                                          else "<mi>" + ssub + "</mi>")), idc)
    elif tex:
        ret = ret.replace("_", "\\_")
    return ret, ln > 1

# ---- named functions :1106-1226 ---------------------------------------------------------------

ELEMFUNCS = {"ln", "sin", "cos", "sec", "csc", "tan", "cot", "asec", "acsc", "asin", "acos", "atan", "acot",
             "sinh", "cosh", "tanh", "asinh", "acosh", "atanh"}
ELEMFUNCS |= {f.upper() for f in ELEMFUNCS}

# name -> (tex, content, presentation)
FUNC_ELEM = {"ln": ("\\ln ", "ln", "ln"), "sin": ("\\sin ", "sin", "sin"), "cos": ("\\cos ", "cos", "cos"),
             "tan": ("\\tan ", "tan", "tan"), "cot": ("\\cot ", "cot", "cot"),
             "sinh": ("\\sinh ", "sinh", "sinh"), "cosh": ("\\cosh ", "cosh", "cosh"),
             "tanh": ("\\tanh ", "tanh", "tanh"), "asin": ("\\arcsin ", "arcsin", "arcsin"),
             "acos": ("\\arccos ", "arccos", "arccos"), "atan": ("\\arctan ", "arctan", "arctan"),
             "acot": ("\\operatorname{arccot}", "arccot", "arccot"), "sec": ("\\sec ", "sec", "sec"),
             "csc": ("\\csc ", "csc", "csc"), "asec": ("\\operatorname{arcsec}", "arcsec", "arcsec"),
             "acsc": ("\\operatorname{arccsc}", "arccsc", "arccsc"),
             "asinh": ("\\operatorname{arsinh}", "arcsinh", "arsinh"),
             "acosh": ("\\operatorname{arcosh}", "arccosh", "arcosh"),
             "atanh": ("\\operatorname{artanh}", "arctanh", "artanh")}
FUNC_CONTENT_ONLY = {"id": "<ident/>", "gcd": "<gcd/>", "lcm": "<lcm/>", "irem": "<rem/>", "rem": "<rem/>",
                     "iquo": "<quotient/>", "quo": "<quotient/>", "arg": "<arg/>", "max": "<max/>",
                     "min": "<min/>", "det": "<determinant/>"}
# name -> (tex, content or None (generic ci), presentation)
FUNC_SPECIAL = {"Dirac": ("\\operatorname{\\delta}", None, "<mi>&delta;</mi>"),
                "Heaviside": ("\\operatorname{\\theta}", None, "<mi>&theta;</mi>"),
                "Gamma": ("\\Gamma", mml_csymbol("gamma", "hypergeo0"), "<mi mathvariant='normal'>&Gamma;</mi>"),
                "Beta": ("\\Beta", mml_csymbol("beta", "hypergeo0"), "<mi mathvariant='normal'>&Beta;</mi>"),
                "euler": ("\\phi", mml_csymbol("euler", "integer2"), "<mi mathvariant='normal'>&phi;</mi>"),
                "Airy_Ai": ("\\operatorname{Ai}", mml_csymbol("Ai", "airy"), "<mi>Ai</mi>"),
                "Airy_Bi": ("\\operatorname{Bi}", mml_csymbol("Bi", "airy"), "<mi>Bi</mi>"),
                "Psi": ("\\Psi", None, "<mi mathvariant='normal'>&Psi;</mi>"),
                "Zeta": ("\\zeta", None, "<mi mathvariant='normal'>&zeta;</mi>")}


def is_elemfunc(g):
    return isinstance(g, Symb) and g.op in ELEMFUNCS


def func2markup(name, typ, idc=0):  # :1128 ; typ 0 tex, 2 content, 3 presentation
    low = name.lower() if name.upper() == name and name.lower() in FUNC_ELEM else name
    has_id = False
    if low in FUNC_ELEM:
        t, c, p = FUNC_ELEM[low]
        ret = t if typ == 0 else ("<" + c + "/>" if typ == 2 else "<mi>" + p + "</mi>")
    else:
        if typ == 2 and name in FUNC_CONTENT_ONLY:
            ret = FUNC_CONTENT_ONLY[name]
        elif name in FUNC_SPECIAL and not (typ == 2 and FUNC_SPECIAL[name][1] is None):
            t, c, p = FUNC_SPECIAL[name]
            ret = t if typ == 0 else (c if typ == 2 else p)
        else:
            if typ == 2:
                ret = mml_tag("ci", name, idc, "type", "function")
            else:
                ret, _ = idnt2markup(name, typ == 0, idc)
            has_id = True
        if typ == 0 and len(ret) > 1 and my_isalpha(ret[0]):  # :1219
            ret = "\\operatorname{" + ret + "}"
    if typ != 0 and not has_id:
        ret = insert_id(ret, idc, typ == 2)
    return ret

# ---- block helpers :958-1051 ------------------------------------------------------------------


def parenthesize(ml, flags):  # :958
    if flags & F_LATEX:
        ml.latex = "\\left(" + ml.latex + "\\right)"
    if flags & F_PRES:
        ml.markup = "<mfenced>" + ml.markup + "</mfenced>"
    ml.priority = 0


def prepend_minus(ml, flags, circled=False, mrow=True):  # :972
    if not ml.neg:
        return
    if ml.priority > P_MUL:
        parenthesize(ml, flags)
    if flags & F_LATEX:
        ml.latex = ("\\ominus " if circled else "-") + ml.latex
    if flags & F_CONTENT:
        if mrow:
            i = extract_id(ml, True)
            ml.content = mml_tag("apply", "<minus/>" + ml.content, 0, "id" if i else "", i)
        else:
            ml.content = "<apply><minus/>" + ml.content + "</apply>"
    if flags & F_PRES:
        m = "<mo>&ominus;</mo>" if circled else MML_MINUS
        if mrow:
            i = extract_id(ml, False)
            ml.markup = mml_tag("mrow", m + ml.markup, 0, "xref" if i else "", i)
        else:
            ml.markup = m + ml.markup
    ml.neg = False
    ml.priority = P_MUL
    ml.type = T_NEG


def flatten_operands(g):  # :936
    if not isinstance(g.arg, Vect):
        return [g.arg]
    ops = []
    for op in g.arg.items:
        if isinstance(op, Symb) and op.op == g.op:
            ops.extend(flatten_operands(op))
        else:
            ops.append(op)
    return ops


def get_leftright(args, ml, flags, idc):  # :1039
    left = gen2markup(args[0], flags, idc)
    prepend_minus(left, flags)
    if ml is not None and left.priority >= ml.priority:
        parenthesize(left, flags)
    right = gen2markup(args[-1], flags, idc)
    prepend_minus(right, flags)
    if ml is not None and right.priority >= ml.priority:
        parenthesize(right, flags)
    return left, right


def is_partialdiff(g):  # :1053
    df, deg = g, []
    while is_sym(df, "derive") and isinstance(df.arg, Vect):
        a = df.arg.items
        if len(a) != 2 or not is_integer(a[-1]) or not is_positive(a[-1]):
            return None
        deg.append(Int(a[-1].v + 1))
        df = a[0]
    if not isinstance(df, Idnt) or not deg:
        return None
    return df, deg[::-1]


def get_derive_vars(args, vars_):  # :1069
    n = len(args)
    if n == 3 and isinstance(args[1], Idnt) and is_integer(args[-1]) and args[-1].v > 0:
        k = args[-1].v
        if vars_ and vars_[-1][0] == args[1]:
            vars_[-1][1] += k
        else:
            vars_.append([args[1], k])
        return True
    for g in args[1:]:
        if is_sym(g, "dollar"):
            h, m = g.arg.items[0], g.arg.items[-1]
            if not isinstance(h, Idnt) or not is_integer(m) or m.v <= 0:
                return False
            if vars_ and vars_[-1][0] == h:
                vars_[-1][1] += m.v
            vars_.append([h, m.v])
        elif isinstance(g, Idnt):
            if vars_ and vars_[-1][0] == g:
                vars_[-1][1] += 1
            else:
                vars_.append([g, 1])
        elif isinstance(g, Vect):
            if not get_derive_vars(g.items, vars_):
                return False
        else:
            return False
    return True


def ckmatrix(items):
    if not items or not all(isinstance(r, Vect) for r in items):
        return False
    n = len(items[0].items)
    return n > 0 and all(len(r.items) == n and not any(isinstance(x, Vect) for x in r.items)
                         for r in items)

# ---- gen2markup :1233-3627 --------------------------------------------------------------------


def gen2markup(g, flags_orig, idc):
    flags = flags_orig
    toplevel = (flags & F_TOPLEVEL) != 0
    isfactor = (flags & F_FACTOR) != 0
    tex = flags & F_LATEX
    pres = flags & F_PRES
    cont = flags & F_CONTENT
    flags &= ~(F_TOPLEVEL | F_FACTOR)
    ml = MB()
    if isinstance(g, Bool):  # :1308-1321
        name = "true" if g.v else "false"
        if cont:
            ml.content = insert_id("<" + name + ">", idc.inc(), True)
        if pres:
            ml.markup, _ = idnt2markup(name, False, idc.v)
        if tex:
            ml.latex, _ = idnt2markup(name, True)
        ml.type = T_SUBTYPE_IDNT | T_LEADING_DIGIT
        return ml
    if isinstance(g, Int):  # :1323-1337
        ml.type = T_NUM_EXACT | T_LEADING_DIGIT
        ml.neg = g.v < 0
        s = str(abs(g.v))
        if cont:
            ml.content = mml_tag("cn", s, idc.inc(), "type", "integer")
        if pres:
            ml.markup = mml_tag("mn", s, idc.v)
        if tex:
            ml.latex = s
        return ml
    if isinstance(g, Dbl):  # :1362-1407
        ml.type = T_NUM_APPROX | T_LEADING_DIGIT
        ml.neg = g.neg
        s = g.text
        if cont:
            ml.content = mml_tag("cn", s, idc.inc(), "type", "double")
        cpos = s.find('e')
        if cpos >= 0 and atof(s[:cpos]) != 0 and atof(s[cpos + 2:]) != 0:
            ml.priority = P_MUL
            ex = s[cpos + 2:].lstrip('0')
            mant = s[:cpos]
            if '.' not in mant:
                mant += ".0"
            if pres:
                ml.markup = mml_tag("mrow", "<mn>" + mant + "</mn><mo>&times;</mo><msup><mn>10</mn><mrow>"
                                    + (MML_MINUS if s[cpos + 1] == '-' else "") + "<mn>" + ex
                                    + "</mn></mrow></msup>", idc.v)
            if tex:
                ml.latex = mant + "\\times10^{" + ("-" if s[cpos + 1] == '-' else "") + ex + "}"
        else:
            if pres:
                ml.markup = mml_tag("mn", s, idc.v)
            if tex:
                ml.latex = s
        return ml
    if isinstance(g, Cplx):  # :1408-1462
        if is_zero(g.im):
            return gen2markup(g.re, flags, idc)
        isone = abs(numval(g.im)) == 1
        if is_zero(g.re):
            if isone:
                if cont:
                    ml.content = insert_id("<imaginaryi/>", idc.v, True)  # no ++idc (:1415)
                if pres:
                    ml.markup = insert_id(MML_I, idc.v, False)
                if tex:
                    ml.latex = TEX_MATHI
                ml.neg = is_minus_one(g.im)
            else:
                tmp = gen2markup(g.im, flags | F_FACTOR, idc)
                ml.neg = tmp.neg
                ml.priority = P_MUL
                if tmp.priority >= ml.priority:
                    parenthesize(tmp, flags)
                if cont:  # giac puts tmp.markup in the content (:1430)
                    ml.content = mml_tag("apply", "<times/><imaginaryi/>" + tmp.markup, idc.inc())
                if pres:
                    ml.markup = mml_tag("mrow", tmp.markup + MML_ITIMES + MML_I, idc.v)
                if tex:
                    ml.latex = tmp.latex + TEX_ITIMES + TEX_MATHI
            return ml
        ml.priority = P_ADD
        left = gen2markup(g.re, flags, idc)
        prepend_minus(left, flags)
        if left.priority >= ml.priority:
            parenthesize(left, flags)
        right = gen2markup(g.im, flags, idc)
        if right.priority >= P_MUL:
            parenthesize(right, flags)
        if cont:
            ml.content = mml_tag("apply", ("<minus/>" if right.neg else "<plus/>") + left.content
                                 + ("<imaginaryi/>" if isone else
                                    mml_tag("apply", "<times/><imaginaryi/>" + right.content)), idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", left.markup + (MML_MINUS if right.neg else MML_PLUS)
                                + ("" if isone else right.markup + MML_ITIMES) + MML_I, idc.v)
        if tex:
            ml.latex = (left.latex + ("-" if right.neg else "+") + ("" if isone else right.latex + TEX_ITIMES)
                        + TEX_MATHI)
        return ml
    if isinstance(g, Frac):  # :1463-1480
        ml.type = T_FRACTION
        ml.priority = P_MUL
        if isinstance(g.num, Cplx) and not is_zero(g.num.im):
            return gen2markup(Symb("division", Vect([g.num, g.den], 'seq')), flags, idc)
        ml.neg = (g.num.v < 0) != (g.den.v < 0)
        ld, rd = str(abs(g.num.v)), str(abs(g.den.v))
        if cont:
            ml.content = mml_tag("cn", ld + "<sep/>" + rd, idc.inc(), "type", "rational")
        if pres:
            ml.markup = mml_tag("mfrac", "<mn>" + ld + "</mn><mn>" + rd + "</mn>", idc.v)
        if tex:
            ml.latex = "\\frac{" + ld + "}{" + rd + "}"
        return ml
    if isinstance(g, Str):  # :1481-1485
        string2markup(ml, g.s, flags, idc)
        return ml
    if isinstance(g, Idnt):  # :1486-1542
        if is_inf(g):
            c, p, t = "<infinity/>", "&infin;", "\\infty "
        elif g.name == "pi":
            c, p, t = "<pi/>", "&pi;", "\\pi "
        elif g.name == "euler_gamma":
            c, p, t = "<eulergamma/>", "&gamma;", "\\gamma "
        else:
            if cont:
                ml.content = mml_tag("ci", g.name, idc.inc())
            has_sub = False
            isname = False
            if pres:
                ml.markup, isname = idnt2markup(g.name, False, idc.v)
                has_sub = "<msub>" in ml.markup
            if tex:
                ml.latex, isname = idnt2markup(g.name, True)
                has_sub = "_{" in ml.latex
            if has_sub:
                ml.type |= T_HAS_SUBSCRIPT
            if isname:
                ml.type |= T_IDNT_NAME
            return ml
        if cont:
            ml.content = insert_id(c, idc.inc(), True)
        if pres:
            ml.markup = mml_tag("mi", p, idc.v)
        if tex:
            ml.latex = t
        return ml
    if isinstance(g, Func):  # :1647-1657
        if cont:
            ml.content = func2markup(g.name, 2, idc.inc())
        if pres:
            ml.markup = func2markup(g.name, 3, idc.v)
        if tex:
            ml.latex = func2markup(g.name, 0)
        ml.type = T_FUNC
        return ml
    if isinstance(g, Vect):  # :1658-1773
        st = g.sub
        if st in ('list', 'matrix') and ckmatrix(g.items):
            ml.type = T_MATRIX
            rows = []
            for r in g.items:
                cells = []
                if cont:
                    ml.content += "<matrixrow>"
                if pres:
                    ml.markup += "<mtr>"
                for x in r.items:
                    tmp = gen2markup(x, flags, idc)
                    prepend_minus(tmp, flags)
                    if cont:
                        ml.content += tmp.content
                    if pres:
                        ml.markup += "<mtd>" + tmp.markup + "</mtd>"
                    cells.append(tmp.latex)
                if cont:
                    ml.content += "</matrixrow>"
                if pres:
                    ml.markup += "</mtr>"
                rows.append("&".join(cells))
            if cont:
                ml.content = mml_tag("matrix", ml.content, idc.inc())
            ld, rd = ("[", "]") if st == 'list' else ("(", ")")
            if pres:
                ml.markup = mml_tag("mfenced", "<mtable>" + ml.markup + "</mtable>", idc.v, "open", ld,
                                    "close", rd)
            if tex:
                ml.latex = ("\\left" + ld + "\\begin{array}{" + "c" * len(g.items[0].items) + "}"
                            + "\\\\".join(rows) + "\\end{array}\\right" + rd)
            return ml
        parts = []
        for x in g.items:
            tmp = gen2markup(x, flags, idc)
            prepend_minus(tmp, flags)
            if cont:
                ml.content += tmp.content
            if pres:
                ml.markup += tmp.markup
            parts.append(tmp.latex)
        ml.latex = ",".join(parts)
        if st == 'seq':
            if cont:
                ml.content = mml_tag("list", ml.content, idc.inc())
            ld = ld_tex = "" if toplevel else "("
            rd = rd_tex = "" if toplevel else ")"
        elif st == 'set':
            if not g.items:
                if cont:
                    ml.content = insert_id("<emptyset/>", idc.inc(), True)
                if pres:
                    ml.markup = mml_tag("mi", "&empty;", idc.v)
                ml.latex = "\\emptyset "
                return ml
            if cont:
                ml.content = mml_tag("set", ml.content, idc.inc())
            ld, rd, ld_tex, rd_tex = "{", "}", "\\{", "\\}"
        else:  # list, poly1, matrix-subtype non-matrix
            if cont:
                ml.content = mml_tag("vector", ml.content, idc.inc())
            ld = ld_tex = "["
            rd = rd_tex = "]"
        if pres:
            ml.markup = mml_tag("mfenced", ml.markup, idc.v, "open", ld, "close", rd)
        ml.latex = ("\\left" + ld_tex if ld_tex else "") + ml.latex + ("\\right" + rd_tex if rd_tex else "")
        return ml
    if not isinstance(g, Symb):
        raise TypeError(g)
    return symb2markup(g, flags, isfactor, ml, idc)


def binop(ml, flags, idc, left, right, csymb, op_mml, op_tex):
    """The common tail of the comparison arms (:2760-2767, :2797-2804)."""
    if flags & F_CONTENT:
        ml.content = mml_tag("apply", csymb + left.content + right.content, idc.inc())
    if flags & F_PRES:
        ml.markup = mml_tag("mrow", left.markup + "<mo>" + op_mml + "</mo>" + right.markup, idc.v)
    if flags & F_LATEX:
        ml.latex = left.latex + op_tex + right.latex
    return ml


def symb2markup(g, flags, isfactor, ml, idc):  # :1774-3611
    op = g.op
    vectarg = isinstance(g.arg, Vect)
    isbinary = vectarg and len(g.arg.items) == 2
    tex = flags & F_LATEX
    pres = flags & F_PRES
    cont = flags & F_CONTENT
    if op == "neg":  # :1777
        ml = gen2markup(g.arg, flags, idc)
        ml.neg = not ml.neg
        return ml
    if op == "epsilon":  # :1782
        ml.latex = "\\epsilon " if tex else ""
        if pres:
            ml.markup = "&epsilon;"
        if cont:
            ml.content = mml_tag("ci", "epsilon", idc.inc())
        return ml
    if op in ("plus", "pointplus"):  # :1793-1847
        if vectarg:
            ml.priority = P_ADD
            for k, x in enumerate(flatten_operands(g)):
                tmp = gen2markup(x, flags, idc)
                if not tmp.neg and k > 0:
                    if pres:
                        ml.markup += MML_PLUS if op == "plus" else "<mo>&oplus;</mo>"
                    if tex:
                        ml.latex += "+" if op == "plus" else "\\oplus "
                prepend_minus(tmp, flags, op == "pointplus", False)
                if tmp.priority >= ml.priority:
                    parenthesize(tmp, flags)
                ml.content += tmp.content
                ml.markup += tmp.markup
                ml.latex += tmp.latex
            if cont:
                ml.content = mml_tag("apply", ("<plus/>" if op == "plus" else "<ci>.+</ci>") + ml.content,
                                     idc.inc())
            if pres:
                ml.markup = mml_tag("mrow", ml.markup, idc.v)
        elif op == "plus":
            ml = gen2markup(g.arg, flags, idc)
            if not ml.neg and ml.priority == 0:
                if cont:
                    i = extract_id(ml, True)
                    ml.content = mml_tag("apply", "<plus/>" + ml.content, 0, "id", i)
                if pres:
                    i = extract_id(ml, False)
                    ml.markup = mml_tag("mrow", MML_PLUS + ml.markup, 0, "xref", i)
                if tex:
                    ml.latex = "+" + ml.latex
                ml.priority = P_ADD
        return ml
    if op == "pointminus" and isbinary:  # :1848
        ml.priority = P_ADD
        left, right = get_leftright(g.arg.items, None, flags, idc)
        prepend_minus(left, flags)
        if left.priority >= ml.priority:
            parenthesize(left, flags)
        if cont:
            ml.content = mml_tag("apply", "<ci>" + (".+" if right.neg else ".-") + "</ci>" + left.content
                                 + right.content, idc.inc())
        if right.priority >= ml.priority:
            parenthesize(right, flags)
        if pres:
            ml.markup = mml_tag("mrow", left.markup + ("<mo>&oplus;</mo>" if right.neg else "<mo>&ominus;</mo>")
                                + right.markup, idc.v)
        if tex:
            ml.latex = left.latex + ("\\oplus " if right.neg else "\\ominus ") + right.latex
        return ml
    if op in ("prod", "ampersand_times") and vectarg:  # :1869-2032
        return prod2markup(g, flags, ml, idc)
    if op == "pointprod" and vectarg:  # :2033
        ml.priority = P_MUL
        for k, x in enumerate(flatten_operands(g)):
            tmp = gen2markup(x, flags, idc)
            prepend_minus(tmp, flags)
            if tmp.priority >= ml.priority:
                parenthesize(tmp, flags)
            ml.content += tmp.content
            ml.markup += ("" if k == 0 else "<mo>&odot;</mo>") + tmp.markup
            ml.latex += ("" if k == 0 else "\\odot ") + tmp.latex
        if cont:
            ml.content = mml_tag("apply", "<ci>.*</ci>" + ml.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", ml.markup, idc.v)
        return ml
    if op == "inv":  # :2082
        tmp = gen2markup(g.arg, flags, idc)
        ml.neg = tmp.neg
        ml.type = T_FRACTION
        ml.priority = P_MUL
        if cont:
            ml.content = mml_tag("apply", "<divide/><cn type='integer'>1</cn>" + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mfrac", "<mn>1</mn>" + tmp.markup, idc.v)
        if tex:
            ml.latex = "\\frac1{" + tmp.latex + "}"
        return ml
    if op in ("division", "rdiv") and isbinary:  # :2097
        a = g.arg.items
        ml.type = T_FRACTION
        ml.priority = P_MUL
        num = gen2markup(a[0], flags, idc)
        den = gen2markup(a[-1], flags, idc)
        ml.neg = num.neg != den.neg
        if cont:
            ml.content = mml_tag("apply", "<divide/>" + num.content + den.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mfrac", num.markup + den.markup, idc.v)
        if tex:
            ml.latex = "\\frac{" + num.latex + "}{" + den.latex + "}"
        return ml
    if op == "conj":  # :2130
        tmp = gen2markup(g.arg, flags, idc)
        if isfactor and tmp.priority > P_MUL:
            parenthesize(tmp, flags)
        if tmp.ctype(T_LEADING_DIGIT):
            ml.type = T_LEADING_DIGIT
        ml.priority = P_UNARY
        ml.neg = tmp.neg
        if cont:
            ml.content = mml_tag("apply", "<conjugate/>" + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mover", tmp.markup + "<mo>&#xaf;</mo>", idc.v)
        if tex:
            ml.latex = "\\overline{" + tmp.latex + "}"
        return ml
    if op in ("exp", "EXP"):  # :2148
        if is_one(g.arg):
            ml.priority = 0
            if cont:
                ml.content = insert_id("<exponentiale/>", idc.inc(), True)
            if pres:
                ml.markup = insert_id(MML_E, idc.v, False)
            if tex:
                ml.latex = TEX_MATHE
        else:
            tmp = gen2markup(g.arg, flags, idc)
            prepend_minus(tmp, flags)
            ml.priority = P_EXP
            ml.type = T_POWER
            if cont:
                ml.content = mml_tag("apply", "<exp/>" + tmp.content, idc.inc())
            if pres:
                ml.markup = mml_tag("msup", MML_E + tmp.markup, idc.v)
            if tex:
                ml.latex = TEX_MATHE + "^{" + tmp.latex + "}"
        return ml
    if op == "sq":  # :2175
        ml.priority = P_EXP
        ml.type = T_POWER
        tmp = gen2markup(g.arg, flags, idc)
        if tmp.priority >= ml.priority:
            parenthesize(tmp, flags)
        if cont:
            ml.content = mml_tag("apply", "<power/>" + tmp.content + "<cn type='integer'>2</cn>", idc.inc())
        if pres:
            ml.markup = mml_tag("msup", tmp.markup + "<mn>2</mn>", idc.v)
        if tex:
            ml.latex = tmp.latex + "^{2}"
        return ml
    if op == "sqrt":  # :2191
        ml.priority = P_UNARY
        ml.type = T_ROOT
        tmp = gen2markup(g.arg, flags, idc)
        prepend_minus(tmp, flags)
        if cont:
            ml.content = mml_tag("apply", "<root/>" + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("msqrt", tmp.markup, idc.v)
        if tex:
            ml.latex = "\\sqrt{" + tmp.latex + "}"
        return ml
    if op in ("pow", "matpow", "pointpow") and isbinary:  # :2224-2280
        ml.priority = P_EXP
        a = g.arg.items
        left = gen2markup(a[0], flags, idc)
        prepend_minus(left, flags)
        e = numval(a[-1])
        if op == "pow" and e is not None and 2 * e == 1:  # :2230
            ml.type = T_ROOT
            ml.priority = P_UNARY
            if cont:
                ml.content = mml_tag("apply", "<root/>" + left.content, idc.inc())
            if pres:
                ml.markup = mml_tag("msqrt", left.markup, idc.v)
            if tex:
                ml.latex = "\\sqrt{" + left.latex + "}"
            return ml
        right = gen2markup(a[-1], flags, idc)
        prepend_minus(right, flags)
        if cont:
            ml.content = mml_tag("apply", ("<ci>.^</ci>" if op == "pointpow" else "<power/>") + left.content
                                 + right.content, idc.inc())
        if ((flags & F_ELEMPOW) and op == "pow" and right.priority == 0
                and not (e is not None and abs(e) == 1)
                and left.ctype(T_ELEMAPP) and left.split_pos >= 0):  # :2248
            ml.type = T_ELEMAPP
            ml.appl = left.appl
            if pres:
                pos = left.markup.find(">")
                sp = left.split_pos
                ml.markup = (left.markup[:pos + 1] + "<msup>" + left.markup[pos + 1:pos + 1 + sp] + right.markup
                             + "</msup>" + left.markup[sp + pos + 1:])
            if tex:
                k = left.split_pos_tex
                ml.latex = left.latex[:k] + "^{" + right.latex + "}" + left.latex[k:]
        else:
            ml.type = T_POWER
            if left.ctype(T_LEADING_DIGIT):
                ml.type |= T_LEADING_DIGIT
            if left.appl or left.priority >= ml.priority:
                parenthesize(left, flags)
            if pres:
                ml.markup = mml_tag("msup", left.markup + right.markup, idc.v)
            if tex:
                ml.latex = left.latex + "^{" + right.latex + "}"
        return ml
    if op in ("abs", "floor", "ceil"):  # :2281
        tmp = gen2markup(g.arg, flags, idc)
        if op != "abs":
            prepend_minus(tmp, flags)
        ld, rd = {"abs": ("|", "|"), "floor": ("&lfloor;", "&rfloor;"), "ceil": ("&lceil;", "&rceil;")}[op]
        ldt, rdt = {"abs": ("|", "|"), "floor": ("\\lfloor ", "\\rfloor "),
                    "ceil": ("\\lceil ", "\\rceil ")}[op]
        if cont:
            ml.content = mml_tag("apply", {"abs": "<abs/>", "floor": "<floor/>", "ceil": "<ceiling/>"}[op]
                                 + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mfenced", tmp.markup, idc.v, "open", ld, "close", rd)
        if tex:
            ml.latex = "\\left" + ldt + tmp.latex + "\\right" + rdt
        return ml
    if op in ("integrate", "int", "Int", "sum", "somme", "add", "product", "mul"):  # :2375-2487
        return bigop2markup(g, flags, ml, idc)
    if op in ("limit", "limite") and vectarg and len(g.arg.items) > 1:  # :2488-2548
        a = g.arg.items
        ml.priority = P_MUL
        ml.type = T_SUMPROD
        dirn = 0
        e = gen2markup(a[0], flags, idc)
        ml.neg = e.neg
        if e.priority > ml.priority:
            parenthesize(e, flags)
        elif pres:
            e.markup = MML_APPLY + e.markup
        var, dest = MB(), MB()
        if is_sym(a[1], "equal"):
            var = gen2markup(a[1].arg.items[0], flags, idc)
            dest = gen2markup(a[1].arg.items[-1], flags, idc)
            if len(a) == 3 and is_integer(a[2]):
                dirn = a[2].v
        elif len(a) >= 3:
            var = gen2markup(a[1], flags, idc)
            dest = gen2markup(a[2], flags, idc)
            if len(a) == 4 and is_integer(a[3]):
                dirn = a[3].v
        prepend_minus(var, flags)
        prepend_minus(dest, flags)
        if cont:
            s = ""
            if var.content:
                s += mml_tag("bvar", var.content)
                if dest.content:
                    s += mml_tag("condition", mml_tag("tendsto", var.content + dest.content, 0,
                                                      "" if dirn == 0 else "type",
                                                      "below" if dirn < 0 else "above"))
            s += e.content
            ml.content = mml_tag("apply", "<limit/>" + s, idc.inc())
        if pres:
            m = "<munder><mi>lim</mi><mrow>" + var.markup + "<mo>&rarr;</mo>"
            if dirn == 0:
                m += dest.markup
            else:
                m += "<msup>" + dest.markup + (MML_MINUS if dirn < 0 else MML_PLUS) + "</msup>"
            m += "</mrow></munder>" + e.markup
            ml.markup = mml_tag("mrow", m, idc.v)
        if tex:
            ml.latex = ("\\lim_{" + var.latex + "\\to " + dest.latex
                        + ("" if dirn == 0 else ("^{-}" if dirn < 0 else "^{+}")) + "}" + e.latex)
        return ml
    if op in ("derive", "deriver"):  # :2560-2727
        r = derive2markup(g, flags, ml, vectarg, idc)
        if r is not None:
            return r
    if op in ("same", "different", "sto", "array_sto", "equal") and isbinary:  # :2728
        ml.priority = (P_ASSGN if op in ("sto", "array_sto") else (P_EQ if op == "equal" else P_COMP))
        if op in ("sto", "array_sto"):
            right, left = get_leftright(g.arg.items, ml, flags, idc)
        else:
            left, right = get_leftright(g.arg.items, ml, flags, idc)
        opm = {"same": "&#x2a75;", "different": "&#x2260;", "equal": "=", "sto": "&#x2254;",
               "array_sto": "&larr;"}[op]
        opt = {"same": "==", "different": "\\neq ", "equal": "=", "sto": ":=", "array_sto": "\\leftarrow"}[op]
        cs = {"same": "==", "different": "<neq/>", "equal": "<eq/>", "sto": ":=", "array_sto": "=&lt;"}[op]
        if cs[0] != '<':
            cs = mml_tag("ci", cs)
        return binop(ml, flags, idc, left, right, cs, opm, opt)
    if op in ("inferieur_strict", "inferieur_egal", "superieur_strict", "superieur_egal") \
            and isbinary:  # :2770
        ml.priority = P_INEQ
        left, right = get_leftright(g.arg.items, ml, flags, idc)
        opm, opt, cs = {"inferieur_strict": ("&lt;", "<", "<lt/>"),
                        "inferieur_egal": ("&leq;", "\\leq ", "<leq/>"),
                        "superieur_strict": ("&gt;", ">", "<gt/>"),
                        "superieur_egal": ("&geq;", "\\geq ", "<geq/>")}[op]
        return binop(ml, flags, idc, left, right, cs, opm, opt)
    if op == "interval" and isbinary:  # :2942
        ml.priority = P_SET
        left, right = get_leftright(g.arg.items, ml, flags, idc)
        if cont:
            ml.content = mml_tag("interval", left.content + right.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", left.markup + "<mo>&#x2025;</mo>" + right.markup, idc.v)
        if tex:
            ml.latex = left.latex + "\\ldots " + right.latex
        return ml
    if op == "factorial":  # :3317
        ml.priority = P_UNARY
        ml.type = T_FACTORIAL
        tmp = gen2markup(g.arg, flags, idc)
        prepend_minus(tmp, flags)
        if tmp.priority >= ml.priority:
            parenthesize(tmp, flags)
        if cont:
            ml.content = mml_tag("apply", "<factorial/>" + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", tmp.markup + "<mo>!</mo>", idc.v)
        if tex:
            ml.latex = tmp.latex + "\\mathclose{!}"
        return ml
    if op in ("binomial", "nCr", "comb") and isbinary:  # :3334
        left = gen2markup(g.arg.items[0], flags, idc)
        prepend_minus(left, flags)
        right = gen2markup(g.arg.items[-1], flags, idc)
        prepend_minus(right, flags)
        if cont:
            ml.content = mml_tag("apply", mml_csymbol("binomial", "combinat1") + left.content + right.content,
                                 idc.inc())
        if pres:
            ml.markup = mml_tag("mfenced", "<mfrac linethickness='0'>" + left.markup + right.markup + "</mfrac>",
                                idc.v)
        if tex:
            ml.latex = "\\binom{" + left.latex + "}{" + right.latex + "}"
        return ml
    if op == "at" and isbinary:  # :3434
        a = g.arg.items
        ml.priority = P_APPLY
        ml.type = T_HAS_SUBSCRIPT
        tmp = gen2markup(a[0], flags, idc)
        prepend_minus(tmp, flags)
        if tmp.priority > ml.priority or tmp.ctype(T_HAS_SUBSCRIPT):
            parenthesize(tmp, flags)
        if isinstance(a[-1], Vect) and a[-1].sub == 'seq':
            idx = a[-1].items
            sep, sept = MML_ICOMMA, " "
            for x in idx:
                if (not (is_integer(x) and 1 <= x.v <= 9)) and (
                        not isinstance(x, Idnt) or ((len(x.name) != 1 or not my_isalpha(x.name[0]))
                                                    and not is_greek_letter(x.name))):
                    sep, sept = "<mo>,</mo>", ","
                    break
            subs = []
            for k, x in enumerate(idx):
                sub = gen2markup(x, flags, idc)
                prepend_minus(sub, flags)
                ml.content += sub.content
                ml.markup += ("" if k == 0 else sep) + sub.markup
                subs.append(sub.latex)
            if cont:
                ml.content = mml_tag("apply", "<selector/>" + tmp.content + ml.content, idc.inc())
            if pres:
                ml.markup = mml_tag("msub", tmp.markup + mml_tag("mrow", ml.markup), idc.v)
            if tex:
                ml.latex = tmp.latex + "_{" + sept.join(subs) + "}"
        else:
            sub = gen2markup(a[-1], flags, idc)
            prepend_minus(sub, flags)
            if cont:
                ml.content = mml_tag("apply", "<selector/>" + tmp.content + sub.content, idc.inc())
            if pres:
                ml.markup = mml_tag("msub", tmp.markup + sub.markup, idc.v)
            if tex:
                ml.latex = tmp.latex + "_{" + sub.latex + "}"
        return ml
    # general SYMB :3532-3611
    ml.priority = P_APPLY
    if op == "of":
        a = g.arg.items
        tmp = gen2markup(a[0], flags, idc)
        prepend_minus(tmp, flags)
        if not tmp.ctype(T_DERIVATIVE) and tmp.priority >= ml.priority:
            parenthesize(tmp, flags)
        ml.content = tmp.content if cont else ""
        ml.markup = tmp.markup + MML_APPLY if pres else ""
        ml.latex = tmp.latex if tex else ""
        gp = a[-1]
    else:
        gp = g.arg
        if cont:
            ml.content = func2markup(op, 2, idc.inc())
        if pres:
            ml.markup = func2markup(op, 3, idc.v)
        if tex:
            ml.latex = func2markup(op, 0)
    if isinstance(gp, Vect) and gp.sub == 'seq':
        # args are NOT prepend_minus'ed here: f(-x,y) prints f(x,y) (giac bug)
        vs = vsc = ""
        vst = []
        for x in gp.items:
            tmp = gen2markup(x, flags, idc)
            vsc += tmp.content
            vs += tmp.markup
            vst.append(tmp.latex)
        if cont:
            ml.content = mml_tag("apply", ml.content + vsc, idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", ml.markup + mml_tag("mfenced", vs), idc.v)
        if tex:
            ml.latex += "\\left(" + ",".join(vst) + "\\right)"
    else:
        tmp = gen2markup(gp, flags, idc)
        prepend_minus(tmp, flags)
        if is_elemfunc(g):
            ml.type = T_ELEMAPP
            ml.split_pos = len(ml.markup)
            ml.split_pos_tex = len(ml.latex)
            ml.appl = (tmp.priority == 0 and not tmp.ctype(T_NUM_EXACT) and not tmp.ctype(T_NUM_APPROX)
                       and not tmp.ctype(T_SUBTYPE_IDNT) and not tmp.ctype(T_IDNT_NAME))
        if not ml.appl:
            parenthesize(tmp, flags)
        if cont:
            ml.content = mml_tag("apply", ml.content + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("mrow", ml.markup + MML_APPLY + tmp.markup, idc.v)
        if tex:
            ml.latex += tmp.latex
    return ml


def prod2markup(g, flags, ml, idc):  # :1869-2032
    tex = flags & F_LATEX
    pres = flags & F_PRES
    cont = flags & F_CONTENT
    args = flatten_operands(g)
    ncnst, cnsts = [], []
    for i in range(len(args) - 1, -1, -1):
        x = args[i]
        if is_integer(x) or isinstance(x, Frac):
            ncnst.append(x)
            del args[i]
        elif isinstance(x, Idnt) and x.name in ("pi", "euler_gamma"):
            cnsts.append(x)
            del args[i]
        elif isinstance(x, Cplx) and is_zero(x.re) and not is_one(x.im):
            if is_minus_one(x.im):
                ml.neg = not ml.neg
            else:
                cnsts.append(x.im)
            args[i] = I
    args = ncnst + cnsts + args
    ml.priority = P_MUL
    num_count = den_count = nc = dc = neg_count = 0
    hasleadingfrac = True
    np_ = dp = -1
    for k, x in enumerate(args):
        if is_sym(x, "inv"):
            den_count += 1
            if not hasleadingfrac:
                continue
            if is_integer(x.arg) and dp < 0:
                dp = k
            else:
                dp = -1
                hasleadingfrac = False
        else:
            num_count += 1
            if hasleadingfrac and np_ < 0 and is_integer(x):
                np_ = k
    if num_count == 1 or dp < 0:
        hasleadingfrac = False
    if hasleadingfrac:
        di = args[dp].arg
        del args[dp]
        if np_ > dp:
            np_ -= 1
        if np_ < 0:
            args.insert(0, Frac(Int(1), di))
        else:
            ni = args[np_]
            del args[np_]
            args.insert(0, Frac(ni, di))
        num_count = sum(1 for x in args if not is_sym(x, "inv"))
        den_count = len(args) - num_count
    pden, pnum = MB(), MB()
    num = numc = numt = den = denc = dent = ""
    for x in args:
        isinv = is_sym(x, "inv")
        flg = flags
        if (isinv and den_count > 1) or (not isinv and num_count > 1):
            flg |= F_FACTOR
        fc = x.arg if isinv else x
        tmp = gen2markup(fc, flg, idc)
        if tmp.neg:
            neg_count += 1
        if isinv:
            dc += 1
            prev = pden
        else:
            nc += 1
            prev = pnum
        if (((not isinv and num_count > 1) or (isinv and den_count > 1))
                and (tmp.priority > ml.priority or (tmp.priority == ml.priority and not isinstance(fc, Cplx)))
                and (not tmp.ctype(T_FRACTION)
                     and (not tmp.ctype(T_SUMPROD) or (not isinv and nc < num_count)
                          or (isinv and dc < den_count)))):
            parenthesize(tmp, flags)
        is_cdot = (tmp.ctype(T_LEADING_DIGIT) or tmp.ctype(T_SUBTYPE_IDNT) or tmp.ctype(T_FACTORIAL)
                   or tmp.ctype(T_IDNT_NAME) or prev.ctype(T_IDNT_NAME) or prev.ctype(T_SUBTYPE_IDNT)
                   or (tmp.ctype(T_FRACTION) and prev.ctype(T_FRACTION))
                   or (prev.ctype(T_ELEMAPP) and prev.appl))
        sign = MML_CDOT if is_cdot else MML_ITIMES
        sign_tex = "\\cdot " if is_cdot else TEX_ITIMES
        if isinv:
            pden = tmp
            denc += tmp.content
            den += ("" if den == "" else sign) + tmp.markup
            dent += ("" if dent == "" else sign_tex) + tmp.latex
        else:
            pnum = tmp
            numc += tmp.content
            num += ("" if num == "" else sign) + tmp.markup
            numt += ("" if numt == "" else sign_tex) + tmp.latex
    if neg_count % 2:
        ml.neg = not ml.neg
    if numc == "":
        numc = "<cn type='integer'>1</cn>"
    if num == "":
        num = "<mn>1</mn>"
    if numt == "":
        numt = "1"
    if cont:
        if denc == "":
            ml.content = mml_tag("apply", "<times/>" + numc, idc.inc())
        else:
            ml.type = T_FRACTION
            ml.content = mml_tag("apply", "<divide/>" + numc + denc, idc.inc())
    if pres:
        if den == "":
            ml.markup = mml_tag("mrow", num, idc.v)
        else:
            ml.type = T_FRACTION
            ml.markup = mml_tag("mfrac", (mml_tag("mrow", num) if nc > 1 else num)
                                + (mml_tag("mrow", den) if dc > 1 else den), idc.v)
    if tex:
        if dent == "":
            ml.latex = numt
        else:
            ml.type = T_FRACTION
            ml.latex = "\\frac{" + numt + "}{" + dent + "}"
    return ml


def bigop2markup(g, flags, ml, idc):  # :2375-2487
    tex = flags & F_LATEX
    pres = flags & F_PRES
    cont = flags & F_CONTENT
    op = g.op
    isint = op in ("integrate", "int")
    issum = op in ("sum", "somme", "add")
    isprod = op in ("product", "mul")
    big = "<mo>" + ("&int;" if isint else ("&sum;" if issum else "&prod;")) + "</mo>"
    big_tex = "\\int " if isint else ("\\sum " if issum else "\\prod ")
    has_lb = has_ub = False
    lb, ub, var = MB(), MB(), MB()
    ml.priority = P_MUL
    ml.type = T_SUMPROD
    if isinstance(g.arg, Vect):
        a = g.arg.items
        e = gen2markup(a[0], flags, idc)
        if len(a) > 1:
            if is_sym(a[1], "equal") and is_sym(a[1].arg.items[-1], "interval"):
                lh = a[1].arg.items[0]
                low, up = a[1].arg.items[-1].arg.items[0], a[1].arg.items[-1].arg.items[-1]
                lb = gen2markup(low, flags, idc)
                has_lb = True
                ub = gen2markup(up, flags, idc)
                has_ub = True
                var = gen2markup(lh, flags, idc)
            else:
                var = gen2markup(a[1], flags, idc)
            prepend_minus(var, flags)
            if isint:
                if var.priority != 0:
                    parenthesize(var, flags)
                if pres:
                    var.markup = MML_ITIMES + MML_D + var.markup
                if tex:
                    var.latex = "\\," + TEX_MATHD + var.latex
        if not has_lb and len(a) > 2:
            lb = gen2markup(a[2], flags, idc)
            has_lb = True
        if not has_ub and len(a) > 3:
            ub = gen2markup(a[3], flags, idc)
            has_ub = True
        if has_lb:
            prepend_minus(lb, flags)
        if has_ub:
            prepend_minus(ub, flags)
        if has_lb and not isint:
            if pres:
                lb.markup = mml_tag("mrow", var.markup + "<mo>=</mo>" + lb.markup)
            if tex:
                lb.latex = var.latex + "=" + lb.latex
        if not has_lb and not isint:
            lb = MB(**var.__dict__)
            has_lb = True
        if has_lb and not has_ub:
            big = mml_tag("munder", big + lb.markup)
            big_tex += "_{" + lb.latex + "}"
        elif has_lb and has_ub:
            big = mml_tag("munderover", big + lb.markup + ub.markup)
            big_tex += "_{" + lb.latex + "}^{" + ub.latex + "}"
    else:
        e = gen2markup(g.arg, flags, idc)
    if isprod:
        prepend_minus(e, flags)
    else:
        ml.neg = e.neg
    if ((isprod and not e.ctype(T_FRACTION) and e.priority >= ml.priority)
            or (not isprod and (e.ctype(T_NEG) or e.priority > ml.priority))):
        parenthesize(e, flags)
    if cont:
        s = ""
        if var.content:
            s += mml_tag("bvar", var.content)
        if has_lb and not has_ub:
            s += mml_tag("domainofapplication" if isint else "condition", lb.content)
        elif has_lb and has_ub:
            s += mml_tag("lowlimit", lb.content) + mml_tag("uplimit", ub.content)
        s += e.content
        ml.content = mml_tag("apply", ("<int/>" if isint else ("<sum/>" if issum else "<product/>")) + s,
                             idc.inc())
    if pres:
        ml.markup = mml_tag("mrow", big + e.markup + (var.markup if isint else ""), idc.v)
    if tex:
        ml.latex = big_tex + e.latex + (var.latex if isint else "")
    return ml


def derive2markup(g, flags, ml, vectarg, idc):  # :2560-2727 ; None = fall through
    tex = flags & F_LATEX
    pres = flags & F_PRES
    cont = flags & F_CONTENT
    pd = is_partialdiff(g)
    if pd is not None:
        f, deg = pd
        ml.priority = P_EXP
        ml.type = T_POWER | T_DERIVATIVE
        left = gen2markup(f, flags, idc)
        right = gen2markup(Vect(deg, 'seq'), flags, idc)
        if tex:
            ml.latex = left.latex + "^{" + right.latex + "}"
        if cont:
            ml.content = mml_tag("apply", "<partialdiff/>" + right.content + left.content, idc.inc())
        if pres:
            ml.markup = mml_tag("msup", left.markup + right.markup, idc.v)
        return ml
    arg = g.arg
    n = 0
    while is_sym(arg, "derive", "deriver"):
        arg = arg.arg
        n += 1
    if not (n == 0 or not vectarg):
        return None
    ml.type = T_DERIVATIVE
    tmp = gen2markup(arg.items[0] if vectarg else arg, flags, idc)
    simp = tmp.priority == 0
    if tmp.priority > P_APPLY:
        parenthesize(tmp, flags)
    ml.neg = tmp.neg
    if n > 0:
        ml.priority = P_EXP
        ml.type |= T_POWER
        n += 1
        N = str(n)
        if cont:
            ml.content = mml_tag("apply", mml_csymbol("nthdiff", "calculus1") + "<cn>" + N + "</cn>" + tmp.content,
                                 idc.inc())
        if pres:
            pm = {2: "<mo>&Prime;</mo>", 3: "<mo>&tprime;</mo>"}.get(n) or mml_tag("mfenced", "<mn>" + N + "</mn>")
            ml.markup = mml_tag("msup", tmp.markup + pm, idc.v)
        if tex:
            ml.latex = tmp.latex + ({2: "''", 3: "'''"}.get(n) or "^{(" + N + ")}")
        return ml
    if n == 0 and not vectarg:
        ml.priority = P_EXP
        if cont:
            ml.content = mml_tag("apply", "<diff/>" + tmp.content, idc.inc())
        if pres:
            ml.markup = mml_tag("msup", tmp.markup + "<mo>&prime;</mo>", idc.v)
        if tex:
            ml.latex = tmp.latex + "'"
        return ml
    if vectarg and len(arg.items) > 1:
        vars_ = []
        if get_derive_vars(arg.items, vars_):
            nv = len(vars_)
            ml.priority = P_MUL
            ml.type |= T_FRACTION
            if nv == 1 and vars_[0][1] == 1:
                var = gen2markup(vars_[0][0], flags, idc)
                if cont:
                    ml.content = mml_tag("apply", "<diff/><bvar>" + var.content + "</bvar>" + tmp.content,
                                         idc.inc())
                if pres:
                    ml.markup = (mml_tag("mfrac", "<mrow>" + MML_D + tmp.markup + "</mrow><mrow>" + MML_D
                                         + var.markup + "</mrow>", idc.v) if simp else
                                 mml_tag("mrow", mml_tag("mfrac", MML_D + mml_tag("mrow", var.markup))
                                         + MML_APPLY + tmp.markup, idc.v))
                if tex:
                    ml.latex = ("\\frac{" + TEX_MATHD + tmp.latex + "}{" + TEX_MATHD + var.latex + "}" if simp
                                else "\\frac{" + TEX_MATHD + "}{" + TEX_MATHD + var.latex + "}" + tmp.latex)
                return ml
            ds = "<mo>&part;</mo>" if nv > 1 else MML_D
            ds_tex = "\\partial " if nv > 1 else TEX_MATHD
            parts = []
            for c, (v, k) in enumerate(vars_):
                var = gen2markup(v, flags, idc)
                s = str(k)
                if cont:
                    ml.content += mml_tag("bvar", var.content + ("" if k == 1 else
                                                                 "<degree><cn>" + s + "</cn></degree>"), idc.v)
                if pres:
                    ml.markup += (ds + var.markup if k == 1 else
                                  mml_tag("msup", "<mrow>" + ds + var.markup + "</mrow><mn>" + s + "</mn>"))
                    if c + 1 < nv:
                        ml.markup += MML_ITIMES
                parts.append(ds_tex + var.latex + ("" if k == 1 else "^{" + s + "}"))
            body = TEX_ITIMES.join(parts)
            ns = str(sum(k for _, k in vars_))
            if cont:
                ml.content = mml_tag("apply", ("<partialdiff/>" if nv > 1 else "<diff/>") + ml.content
                                     + "<degree><cn type='integer'>" + ns + "</cn></degree>" + tmp.content,
                                     idc.inc())
            if pres:
                ml.markup = (mml_tag("mfrac", "<mrow><msup>" + ds + "<mn>" + ns + "</mn></msup>" + tmp.markup
                                     + "</mrow><mrow>" + ml.markup + "</mrow>", idc.v) if simp else
                             mml_tag("mrow", mml_tag("mfrac", "<msup>" + ds + "<mn>" + ns + "</mn></msup>"
                                                     + "<mrow>" + ml.markup + "</mrow>") + MML_APPLY + tmp.markup,
                                     idc.v))
            if tex:
                ml.latex = ("\\frac{" + ds_tex + "^{" + ns + "}" + tmp.latex + "}{" + body + "}" if simp
                            else "\\frac{" + ds_tex + "^{" + ns + "}}{" + body + "}" + tmp.latex)
            return ml
    return None

# ---- entries :3629-3776 -----------------------------------------------------------------------


def export_latex(g):  # :3629
    ml = gen2markup(g, F_TOPLEVEL | F_LATEX | F_ELEMPOW, Idc())
    prepend_minus(ml, F_LATEX)
    return ml.latex


def export_mathml_content(g):  # :3747
    flags = F_TOPLEVEL | F_CONTENT
    ml = gen2markup(g, flags, Idc())
    prepend_minus(ml, flags)
    return "<math " + HEADER + ">" + ml.content + "</math>"


def export_mathml_presentation(g):  # :3755
    flags = F_TOPLEVEL | F_PRES | F_ELEMPOW
    ml = gen2markup(g, flags, Idc())
    prepend_minus(ml, flags)
    return "<math " + HEADER + ">" + ml.markup + "</math>"


def export_mathml(g, printed):  # :3765 ; printed = the input's print (annotation)
    flags = F_TOPLEVEL | F_PRES | F_CONTENT | F_ELEMPOW
    ml = gen2markup(g, flags, Idc())
    prepend_minus(ml, flags)
    return ("<math " + HEADER + "><semantics>" + ml.markup
            + "<annotation-xml encoding='MathML-Content'>" + ml.content + "</annotation-xml>"
            + "<annotation encoding='giac'>" + str_to_mml(printed) + "</annotation></semantics></math>")
