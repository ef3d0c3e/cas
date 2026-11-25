"""Python port of giac markup.cc gen2markup, LaTeX arms only (giac 2.0.0.21, markup.cc (source: C0009)).

Each function cites the markup.cc lines it ports. MathML, scheme and unit arms are dropped.
Gen tree: Int, Bool, Frac, Dbl, Cplx, Idnt, Str, Vect, Symb (giac-shaped, see parse.py).
"""
import re
from dataclasses import dataclass, field

# ---- tree -------------------------------------------------------------------------------------

@dataclass
class Int:
    v: int

@dataclass
class Bool:
    v: bool

@dataclass
class Frac:            # giac _FRAC; num may be Cplx (then printed as division, :1466)
    num: object
    den: object

@dataclass
class Dbl:             # giac _DOUBLE_; text = giac print of |value| (12 significant digits)
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
class Vect:            # sub: 'list' (subtype 0), 'seq', 'set', 'matrix', 'poly1'
    items: list
    sub: str = 'list'

@dataclass
class Symb:            # giac symbolic: sommet name, feuille (a Vect 'seq' when several args)
    op: str
    arg: object

I = Cplx(Int(0), Int(1))

# ---- constants :49-115 ------------------------------------------------------------------------

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

F_TOPLEVEL, F_ELEMPOW, F_FACTOR = 2, 32, 64


@dataclass
class MB:              # MarkupBlock :179-192, latex fields only
    priority: int = 0
    type: int = T_GENERAL
    neg: bool = False
    latex: str = ""
    split_pos: int = -1
    split_pos_tex: int = 0
    appl: bool = False

    def ctype(self, t):
        return (self.type & t) != 0

# ---- small giac helpers -----------------------------------------------------------------------

_ATOF = re.compile(r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[+-]?(?:inf(?:inity)?|nan))",
                   re.IGNORECASE)


def atof(s):
    """C atof: longest float prefix, 0 when none."""
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


def is_digits(s):  # :171 (true for "")
    return all(c.isdigit() for c in s)


def is_integer(g):
    return isinstance(g, Int)


def numval(g):
    """Python number for Int/Frac(int)/Dbl, else None."""
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


def is_positive(g):     # giac is_positive: >= 0 for numbers
    v = numval(g)
    return v is None or v >= 0


def is_sym(g, *ops):
    return isinstance(g, Symb) and g.op in ops


def vectarg_of(g):
    return isinstance(g.arg, Vect)


def args_of(g):
    return g.arg.items if isinstance(g.arg, Vect) else [g.arg]


def is_inf(g):
    return isinstance(g, Idnt) and g.name in ("infinity", "unsigned_inf")

# ---- strings :304-468 -------------------------------------------------------------------------


def trim_string(s):  # :304 -> (trimmed, indent)
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


def str_to_tex(s_orig, quote):  # :333 (ind=false for QUOTE strings)
    s, _ = trim_string(s_orig)
    ret = "".join("\\" + c if c in "\\%_${}^" else c for c in s)
    return "``" + ret + "''" if quote else ret


def string2markup(ml, s_orig, quote=True):  # :408, tex arm
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
    ml.latex = "\\\\".join("\\text{" + str_to_tex(l, quote) + "}" for l in lines)
    if len(lines) > 1:
        ml.latex = "\\begin{array}{l}" + ml.latex + "\\end{array}"

# ---- identifiers :621-934 (tex, non-unit) -----------------------------------------------------


def idnt2markup(s_orig):
    """-> (latex, isname)."""
    if s_orig.endswith("_") and is_greek_letter(s_orig[:-1]):  # :845
        s = s_orig[:-1]
        if s in ("phi", "epsilon"):
            s = "var" + s
        return "\\" + s + "\\_", False
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
    if is_greek_letter(s):
        if s in ("phi", "epsilon"):
            s = "var" + s
        s = "\\" + s
        ln = 1
    len_sub = len(ssub)
    if ssub and is_greek_letter(ssub):
        if ssub in ("phi", "epsilon"):
            ssub = "var" + ssub
        ssub = "\\" + ssub
        len_sub = 1
    if is_double_letter(s):
        s = "\\mathbf{" + s[0] + "}"
        ln = 1
    if ln > 3 and s[0] == '`' and s[1] == ' ' and s[-1] == '`':
        s = s[2:2 + ln - 3]
    ret = s if ln == 1 else "\\mathrm{" + s + "}"
    if ssub:
        tokens = ssub.split("_")
        if all(t != "" and (t == "0" or atof(t) != 0) for t in tokens):
            ssub = ",".join(tokens)
        ret += "_{" + (ssub if (len_sub == 1 or atof(ssub) != 0) else "\\mathrm{" + ssub + "}") + "}"
    else:
        ret = ret.replace("_", "\\_")
    return ret, ln > 1

# ---- named functions :1106-1226 ---------------------------------------------------------------

ELEMFUNCS = {"ln", "sin", "cos", "sec", "csc", "tan", "cot", "asec", "acsc", "asin", "acos", "atan", "acot",
             "sinh", "cosh", "tanh", "asinh", "acosh", "atanh"}
ELEMFUNCS |= {f.upper() for f in ELEMFUNCS}

FUNC_TEX = {"ln": "\\ln ", "sin": "\\sin ", "cos": "\\cos ", "tan": "\\tan ", "cot": "\\cot ",
            "sinh": "\\sinh ", "cosh": "\\cosh ", "tanh": "\\tanh ", "asin": "\\arcsin ",
            "acos": "\\arccos ", "atan": "\\arctan ", "acot": "\\operatorname{arccot}",
            "sec": "\\sec ", "csc": "\\csc ", "asec": "\\operatorname{arcsec}",
            "acsc": "\\operatorname{arccsc}", "asinh": "\\operatorname{arsinh}",
            "acosh": "\\operatorname{arcosh}", "atanh": "\\operatorname{artanh}"}
FUNC_SPECIAL = {"Dirac": "\\operatorname{\\delta}", "Heaviside": "\\operatorname{\\theta}",
                "Gamma": "\\Gamma", "Beta": "\\Beta", "euler": "\\phi", "Airy_Ai": "\\operatorname{Ai}",
                "Airy_Bi": "\\operatorname{Bi}", "Psi": "\\Psi", "Zeta": "\\zeta"}


def is_elemfunc(g):
    return isinstance(g, Symb) and g.op in ELEMFUNCS


def func2markup(name):  # :1128, tex
    low = name.lower() if name.upper() == name and name.lower() in FUNC_TEX else name
    if low in FUNC_TEX:
        return FUNC_TEX[low]
    if name in FUNC_SPECIAL:
        ret = FUNC_SPECIAL[name]
    else:
        ret, _ = idnt2markup(name)
    if len(ret) > 1 and my_isalpha(ret[0]):  # :1219
        ret = "\\operatorname{" + ret + "}"
    return ret

# ---- block helpers :958-1051 ------------------------------------------------------------------


def parenthesize(ml):  # :958
    ml.latex = "\\left(" + ml.latex + "\\right)"
    ml.priority = 0


def prepend_minus(ml, circled=False):  # :972
    if not ml.neg:
        return
    if ml.priority > P_MUL:
        parenthesize(ml)
    ml.latex = ("\\ominus " if circled else "-") + ml.latex
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


def get_leftright(args, ml, flags):  # :1039 ; ml None = no wrapping
    left = gen2markup(args[0], flags)
    prepend_minus(left)
    if ml is not None and left.priority >= ml.priority:
        parenthesize(left)
    right = gen2markup(args[-1], flags)
    prepend_minus(right)
    if ml is not None and right.priority >= ml.priority:
        parenthesize(right)
    return left, right


def is_partialdiff(g):  # :1053 -> (f, deg) or None
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


def get_derive_vars(args, vars_):  # :1069 ; vars_ = list of [idnt, count]
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
            vars_.append([h, m.v])  # giac pushes even after merging (:1086-1087)
        elif isinstance(g, Idnt):
            if vars_ and vars_[-1][0] == g:
                vars_[-1][1] += 1
            else:
                vars_.append([g, 1])
        elif isinstance(g, Vect):
            if not get_derive_vars(g.items, vars_):  # recursion skips g.items[0] (giac quirk)
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


def gen2markup(g, flags_orig):
    flags = flags_orig
    toplevel = (flags & F_TOPLEVEL) != 0
    isfactor = (flags & F_FACTOR) != 0
    flags &= ~(F_TOPLEVEL | F_FACTOR)
    ml = MB()
    if isinstance(g, Bool):  # :1308-1321 (boolean _INT_ subtype)
        ml.latex, _ = idnt2markup("true" if g.v else "false")
        ml.type = T_SUBTYPE_IDNT | T_LEADING_DIGIT
        return ml
    if isinstance(g, Int):  # :1323-1337
        ml.type = T_NUM_EXACT | T_LEADING_DIGIT
        ml.neg = g.v < 0
        ml.latex = str(abs(g.v))
        return ml
    if isinstance(g, Dbl):  # :1362-1407
        ml.type = T_NUM_APPROX | T_LEADING_DIGIT
        ml.neg = g.neg
        s = g.text
        cpos = s.find('e')
        if cpos >= 0 and atof(s[:cpos]) != 0 and atof(s[cpos + 2:]) != 0:
            ml.priority = P_MUL
            ex = s[cpos + 2:].lstrip('0')
            mant = s[:cpos]
            if '.' not in mant:
                mant += ".0"
            ml.latex = mant + "\\times10^{" + ("-" if s[cpos + 1] == '-' else "") + ex + "}"
        else:
            ml.latex = s
        return ml
    if isinstance(g, Cplx):  # :1408-1462
        if is_zero(g.im):
            return gen2markup(g.re, flags)
        isone = abs(numval(g.im)) == 1
        if is_zero(g.re):
            if isone:
                ml.latex = TEX_MATHI
                ml.neg = is_minus_one(g.im)
            else:
                tmp = gen2markup(g.im, flags | F_FACTOR)
                ml.neg = tmp.neg
                ml.priority = P_MUL
                if tmp.priority >= ml.priority:
                    parenthesize(tmp)
                ml.latex = tmp.latex + TEX_ITIMES + TEX_MATHI
            return ml
        ml.priority = P_ADD
        left = gen2markup(g.re, flags)
        prepend_minus(left)
        if left.priority >= ml.priority:
            parenthesize(left)
        right = gen2markup(g.im, flags)
        if right.priority >= P_MUL:
            parenthesize(right)
        ml.latex = (left.latex + ("-" if right.neg else "+") + ("" if isone else right.latex + TEX_ITIMES)
                    + TEX_MATHI)
        return ml
    if isinstance(g, Frac):  # :1463-1480
        ml.type = T_FRACTION
        ml.priority = P_MUL
        if isinstance(g.num, Cplx) and not is_zero(g.num.im):
            return gen2markup(Symb("division", Vect([g.num, g.den], 'seq')), flags)
        ml.neg = (g.num.v < 0) != (g.den.v < 0)
        ml.latex = "\\frac{" + str(abs(g.num.v)) + "}{" + str(abs(g.den.v)) + "}"
        return ml
    if isinstance(g, Str):  # :1481-1485
        string2markup(ml, g.s, quote=True)
        return ml
    if isinstance(g, Idnt):  # :1486-1542
        if is_inf(g):
            ml.latex = "\\infty "
        elif g.name == "pi":
            ml.latex = "\\pi "
        elif g.name == "euler_gamma":
            ml.latex = "\\gamma "
        else:
            ml.latex, isname = idnt2markup(g.name)
            if "_{" in ml.latex:
                ml.type |= T_HAS_SUBSCRIPT
            if isname:
                ml.type |= T_IDNT_NAME
        return ml
    if isinstance(g, Vect):  # :1658-1773
        st = g.sub
        if st in ('list', 'matrix') and ckmatrix(g.items):
            ml.type = T_MATRIX
            rows = []
            for r in g.items:
                cells = []
                for x in r.items:
                    tmp = gen2markup(x, flags)
                    prepend_minus(tmp)
                    cells.append(tmp.latex)
                rows.append("&".join(cells))
            ld, rd = ("[", "]") if st == 'list' else ("(", ")")
            ml.latex = ("\\left" + ld + "\\begin{array}{" + "c" * len(g.items[0].items) + "}"
                        + "\\\\".join(rows) + "\\end{array}\\right" + rd)
            return ml
        parts = []
        for x in g.items:
            tmp = gen2markup(x, flags)
            prepend_minus(tmp)
            parts.append(tmp.latex)
        ml.latex = ",".join(parts)
        if st == 'seq':
            ld, rd = ("", "") if toplevel else ("(", ")")
        elif st == 'set':
            if not g.items:
                ml.latex = "\\emptyset "
                return ml
            ld, rd = "\\{", "\\}"
        else:  # list, poly1, matrix-subtype non-matrix
            ld, rd = "[", "]"
        ml.latex = ("\\left" + ld if ld else "") + ml.latex + ("\\right" + rd if rd else "")
        return ml
    if not isinstance(g, Symb):
        raise TypeError(g)
    return symb2markup(g, flags, isfactor, ml)


def symb2markup(g, flags, isfactor, ml):  # :1774-3611
    op = g.op
    vectarg = isinstance(g.arg, Vect)
    isbinary = vectarg and len(g.arg.items) == 2
    if op == "neg":  # :1777
        ml = gen2markup(g.arg, flags)
        ml.neg = not ml.neg
        return ml
    if op == "epsilon":  # :1782
        ml.latex = "\\epsilon "
        return ml
    if op in ("plus", "pointplus"):  # :1793-1847
        if vectarg:
            ml.priority = P_ADD
            for k, x in enumerate(flatten_operands(g)):
                tmp = gen2markup(x, flags)
                if not tmp.neg and k > 0:
                    ml.latex += "+" if op == "plus" else "\\oplus "
                prepend_minus(tmp, op == "pointplus")
                if tmp.priority >= ml.priority:
                    parenthesize(tmp)
                ml.latex += tmp.latex
        elif op == "plus":
            ml = gen2markup(g.arg, flags)
            if not ml.neg and ml.priority == 0:
                ml.latex = "+" + ml.latex
                ml.priority = P_ADD
        return ml
    if op == "pointminus" and isbinary:  # :1848
        ml.priority = P_ADD
        left, right = get_leftright(g.arg.items, None, flags)
        prepend_minus(left)
        if left.priority >= ml.priority:
            parenthesize(left)
        if right.priority >= ml.priority:
            parenthesize(right)
        ml.latex = left.latex + ("\\oplus " if right.neg else "\\ominus ") + right.latex
        return ml
    if op in ("prod", "ampersand_times") and vectarg:  # :1869-2032
        return prod2markup(g, flags, ml)
    if op == "pointprod" and vectarg:  # :2033
        ml.priority = P_MUL
        for k, x in enumerate(flatten_operands(g)):
            tmp = gen2markup(x, flags)
            prepend_minus(tmp)
            if tmp.priority >= ml.priority:
                parenthesize(tmp)
            ml.latex += ("" if k == 0 else "\\odot ") + tmp.latex
        return ml
    if op == "inv":  # :2082
        tmp = gen2markup(g.arg, flags)
        ml.neg = tmp.neg
        ml.type = T_FRACTION
        ml.priority = P_MUL
        ml.latex = "\\frac1{" + tmp.latex + "}"
        return ml
    if op in ("division", "rdiv") and isbinary:  # :2097
        a = g.arg.items
        ml.type = T_FRACTION
        ml.priority = P_MUL
        num = gen2markup(a[0], flags)
        den = gen2markup(a[-1], flags)
        ml.neg = num.neg != den.neg
        ml.latex = "\\frac{" + num.latex + "}{" + den.latex + "}"
        return ml
    if op == "conj":  # :2130
        tmp = gen2markup(g.arg, flags)
        if isfactor and tmp.priority > P_MUL:
            parenthesize(tmp)
        if tmp.ctype(T_LEADING_DIGIT):
            ml.type = T_LEADING_DIGIT
        ml.priority = P_UNARY
        ml.neg = tmp.neg
        ml.latex = "\\overline{" + tmp.latex + "}"
        return ml
    if op in ("exp", "EXP"):  # :2148
        if is_one(g.arg):
            ml.priority = 0
            ml.latex = TEX_MATHE
        else:
            tmp = gen2markup(g.arg, flags)
            prepend_minus(tmp)
            ml.priority = P_EXP
            ml.type = T_POWER
            ml.latex = TEX_MATHE + "^{" + tmp.latex + "}"
        return ml
    if op == "sq":  # :2175
        ml.priority = P_EXP
        ml.type = T_POWER
        tmp = gen2markup(g.arg, flags)
        if tmp.priority >= ml.priority:
            parenthesize(tmp)
        ml.latex = tmp.latex + "^{2}"
        return ml
    if op == "sqrt":  # :2191
        ml.priority = P_UNARY
        ml.type = T_ROOT
        tmp = gen2markup(g.arg, flags)
        prepend_minus(tmp)
        ml.latex = "\\sqrt{" + tmp.latex + "}"
        return ml
    if op in ("pow", "matpow", "pointpow") and isbinary:  # :2224-2280
        ml.priority = P_EXP
        a = g.arg.items
        left = gen2markup(a[0], flags)
        prepend_minus(left)
        e = numval(a[-1])
        if op == "pow" and e is not None and 2 * e == 1:  # is_one(2*exponent) :2230
            ml.type = T_ROOT
            ml.priority = P_UNARY
            ml.latex = "\\sqrt{" + left.latex + "}"
            return ml
        right = gen2markup(a[-1], flags)
        prepend_minus(right)
        if ((flags & F_ELEMPOW) and op == "pow" and right.priority == 0
                and not (e is not None and abs(e) == 1)
                and left.ctype(T_ELEMAPP) and left.split_pos >= 0):  # :2248
            ml.type = T_ELEMAPP
            ml.appl = left.appl
            k = left.split_pos_tex
            ml.latex = left.latex[:k] + "^{" + right.latex + "}" + left.latex[k:]
        else:
            ml.type = T_POWER
            if left.ctype(T_LEADING_DIGIT):
                ml.type |= T_LEADING_DIGIT
            if left.appl or left.priority >= ml.priority:
                parenthesize(left)
            ml.latex = left.latex + "^{" + right.latex + "}"
        return ml
    if op in ("abs", "floor", "ceil"):  # :2281
        tmp = gen2markup(g.arg, flags)
        if op != "abs":
            prepend_minus(tmp)  # abs drops the child's pending minus (giac quirk)
        ld, rd = {"abs": ("|", "|"), "floor": ("\\lfloor ", "\\rfloor "),
                  "ceil": ("\\lceil ", "\\rceil ")}[op]
        ml.latex = "\\left" + ld + tmp.latex + "\\right" + rd
        return ml
    if op in ("integrate", "int", "Int", "sum", "somme", "add", "product", "mul"):  # :2375-2487
        return bigop2markup(g, flags, ml)
    if op in ("limit", "limite") and vectarg and len(g.arg.items) > 1:  # :2488-2548
        a = g.arg.items
        ml.priority = P_MUL
        ml.type = T_SUMPROD
        dirn = 0
        e = gen2markup(a[0], flags)
        ml.neg = e.neg
        if e.priority > ml.priority:
            parenthesize(e)
        var = dest = MB()
        if is_sym(a[1], "equal"):
            var = gen2markup(a[1].arg.items[0], flags)
            dest = gen2markup(a[1].arg.items[-1], flags)
            if len(a) == 3 and is_integer(a[2]):
                dirn = a[2].v
        elif len(a) >= 3:
            var = gen2markup(a[1], flags)
            dest = gen2markup(a[2], flags)
            if len(a) == 4 and is_integer(a[3]):
                dirn = a[3].v
        prepend_minus(var)
        prepend_minus(dest)
        ml.latex = ("\\lim_{" + var.latex + "\\to " + dest.latex
                    + ("" if dirn == 0 else ("^{-}" if dirn < 0 else "^{+}")) + "}" + e.latex)
        return ml
    if op in ("derive", "deriver"):  # :2560-2727
        r = derive2markup(g, flags, ml, vectarg)
        if r is not None:
            return r  # None: giac falls through with ml's derive-set fields kept
    if op in ("same", "different", "sto", "array_sto", "equal") and isbinary:  # :2728
        ml.priority = (P_ASSGN if op in ("sto", "array_sto") else (P_EQ if op == "equal" else P_COMP))
        if op in ("sto", "array_sto"):
            right, left = get_leftright(g.arg.items, ml, flags)
        else:
            left, right = get_leftright(g.arg.items, ml, flags)
        opt = {"same": "==", "different": "\\neq ", "equal": "=", "sto": ":=",
               "array_sto": "\\leftarrow"}[op]
        ml.latex = left.latex + opt + right.latex
        return ml
    if op in ("inferieur_strict", "inferieur_egal", "superieur_strict", "superieur_egal") \
            and isbinary:  # :2770
        ml.priority = P_INEQ
        left, right = get_leftright(g.arg.items, ml, flags)
        opt = {"inferieur_strict": "<", "inferieur_egal": "\\leq ", "superieur_strict": ">",
               "superieur_egal": "\\geq "}[op]
        ml.latex = left.latex + opt + right.latex
        return ml
    if op == "interval" and isbinary:  # :2942
        ml.priority = P_SET
        left, right = get_leftright(g.arg.items, ml, flags)
        ml.latex = left.latex + "\\ldots " + right.latex
        return ml
    if op == "factorial":  # :3317
        ml.priority = P_UNARY
        ml.type = T_FACTORIAL
        tmp = gen2markup(g.arg, flags)
        prepend_minus(tmp)
        if tmp.priority >= ml.priority:
            parenthesize(tmp)
        ml.latex = tmp.latex + "\\mathclose{!}"
        return ml
    if op in ("binomial", "nCr", "comb") and isbinary:  # :3334
        left = gen2markup(g.arg.items[0], flags)
        prepend_minus(left)
        right = gen2markup(g.arg.items[-1], flags)
        prepend_minus(right)
        ml.latex = "\\binom{" + left.latex + "}{" + right.latex + "}"
        return ml
    if op == "at" and isbinary:  # :3434
        a = g.arg.items
        ml.priority = P_APPLY
        ml.type = T_HAS_SUBSCRIPT
        tmp = gen2markup(a[0], flags)
        prepend_minus(tmp)
        if tmp.priority > ml.priority or tmp.ctype(T_HAS_SUBSCRIPT):
            parenthesize(tmp)
        if isinstance(a[-1], Vect) and a[-1].sub == 'seq':
            idx = a[-1].items
            sept = " "
            for x in idx:
                if (not (is_integer(x) and 1 <= x.v <= 9)) and (
                        not isinstance(x, Idnt) or ((len(x.name) != 1 or not my_isalpha(x.name[0]))
                                                    and not is_greek_letter(x.name))):
                    sept = ","
                    break
            subs = []
            for x in idx:
                sub = gen2markup(x, flags)
                prepend_minus(sub)
                subs.append(sub.latex)
            ml.latex = tmp.latex + "_{" + sept.join(subs) + "}"
        else:
            sub = gen2markup(a[-1], flags)
            prepend_minus(sub)
            ml.latex = tmp.latex + "_{" + sub.latex + "}"
        return ml
    # general SYMB :3532-3611
    ml.priority = P_APPLY
    if op == "of":
        a = g.arg.items
        tmp = gen2markup(a[0], flags)
        prepend_minus(tmp)
        if not tmp.ctype(T_DERIVATIVE) and tmp.priority >= ml.priority:
            parenthesize(tmp)
        ml.latex = tmp.latex
        gp = a[-1]
    else:
        gp = g.arg
        ml.latex = func2markup(op)
    if isinstance(gp, Vect) and gp.sub == 'seq':
        # args are NOT prepend_minus'ed here: f(-x,y) prints f\left(x,y\right) (giac bug)
        ml.latex += "\\left(" + ",".join(gen2markup(x, flags).latex for x in gp.items) + "\\right)"
    else:
        tmp = gen2markup(gp, flags)
        prepend_minus(tmp)
        if is_elemfunc(g):
            ml.type = T_ELEMAPP
            ml.split_pos = 0  # markup length: 0 in a latex-only run
            ml.split_pos_tex = len(ml.latex)
            ml.appl = (tmp.priority == 0 and not tmp.ctype(T_NUM_EXACT) and not tmp.ctype(T_NUM_APPROX)
                       and not tmp.ctype(T_SUBTYPE_IDNT) and not tmp.ctype(T_IDNT_NAME))
        if not ml.appl:
            parenthesize(tmp)
        ml.latex += tmp.latex
    return ml


def prod2markup(g, flags, ml):  # :1869-2032
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
            args.insert(0, Frac(Int(1), di))  # giac fraction(1,di): not reduced
        else:
            ni = args[np_]
            del args[np_]
            args.insert(0, Frac(ni, di))
        num_count = sum(1 for x in args if not is_sym(x, "inv"))
        den_count = len(args) - num_count
    pden, pnum = MB(), MB()
    numt = dent = ""
    for x in args:
        isinv = is_sym(x, "inv")
        flg = flags
        if (isinv and den_count > 1) or (not isinv and num_count > 1):
            flg |= F_FACTOR
        fc = x.arg if isinv else x
        tmp = gen2markup(fc, flg)
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
            parenthesize(tmp)
        is_cdot = (tmp.ctype(T_LEADING_DIGIT) or tmp.ctype(T_SUBTYPE_IDNT) or tmp.ctype(T_FACTORIAL)
                   or tmp.ctype(T_IDNT_NAME) or prev.ctype(T_IDNT_NAME) or prev.ctype(T_SUBTYPE_IDNT)
                   or (tmp.ctype(T_FRACTION) and prev.ctype(T_FRACTION))
                   or (prev.ctype(T_ELEMAPP) and prev.appl))
        sign = "\\cdot " if is_cdot else TEX_ITIMES
        if isinv:
            pden = tmp
            dent += ("" if dent == "" else sign) + tmp.latex
        else:
            pnum = tmp
            numt += ("" if numt == "" else sign) + tmp.latex
    if neg_count % 2:
        ml.neg = not ml.neg
    if numt == "":
        numt = "1"
    if dent == "":
        ml.latex = numt
    else:
        ml.type = T_FRACTION
        ml.latex = "\\frac{" + numt + "}{" + dent + "}"
    return ml


def bigop2markup(g, flags, ml):  # :2375-2487
    op = g.op
    isint = op in ("integrate", "int")
    issum = op in ("sum", "somme", "add")
    isprod = op in ("product", "mul")
    big_tex = "\\int " if isint else ("\\sum " if issum else "\\prod ")
    has_lb = has_ub = False
    lb = ub = var = MB()
    ml.priority = P_MUL
    ml.type = T_SUMPROD
    if isinstance(g.arg, Vect):
        a = g.arg.items
        e = gen2markup(a[0], flags)
        if len(a) > 1:
            if is_sym(a[1], "equal") and is_sym(a[1].arg.items[-1], "interval"):
                lh = a[1].arg.items[0]
                low, up = a[1].arg.items[-1].arg.items[0], a[1].arg.items[-1].arg.items[-1]
                lb = gen2markup(low, flags)
                has_lb = True
                ub = gen2markup(up, flags)
                has_ub = True
                var = gen2markup(lh, flags)
            else:
                var = gen2markup(a[1], flags)
            prepend_minus(var)
            if isint:
                if var.priority != 0:
                    parenthesize(var)
                var.latex = "\\," + TEX_MATHD + var.latex
        if not has_lb and len(a) > 2:
            lb = gen2markup(a[2], flags)
            has_lb = True
        if not has_ub and len(a) > 3:
            ub = gen2markup(a[3], flags)
            has_ub = True
        if has_lb:
            prepend_minus(lb)
        if has_ub:
            prepend_minus(ub)
        if has_lb and not isint:
            lb.latex = var.latex + "=" + lb.latex
        if not has_lb and not isint:
            lb = var
            has_lb = True
        if has_lb and not has_ub:
            big_tex += "_{" + lb.latex + "}"
        elif has_lb and has_ub:
            big_tex += "_{" + lb.latex + "}^{" + ub.latex + "}"
    else:
        e = gen2markup(g.arg, flags)
    if isprod:
        prepend_minus(e)
    else:
        ml.neg = e.neg
    if ((isprod and not e.ctype(T_FRACTION) and e.priority >= ml.priority)
            or (not isprod and (e.ctype(T_NEG) or e.priority > ml.priority))):
        parenthesize(e)
    ml.latex = big_tex + e.latex + (var.latex if isint else "")
    return ml


def derive2markup(g, flags, ml, vectarg):  # :2560-2727 ; None = fall through
    pd = is_partialdiff(g)
    if pd is not None:
        f, deg = pd
        ml.priority = P_EXP
        ml.type = T_POWER | T_DERIVATIVE
        left = gen2markup(f, flags)
        right = gen2markup(Vect(deg, 'seq'), flags)
        ml.latex = left.latex + "^{" + right.latex + "}"
        return ml
    arg = g.arg
    n = 0
    while is_sym(arg, "derive", "deriver"):
        arg = arg.arg
        n += 1
    if not (n == 0 or not vectarg):
        return None
    ml.type = T_DERIVATIVE
    tmp = gen2markup(arg.items[0] if vectarg else arg, flags)
    simp = tmp.priority == 0
    if tmp.priority > P_APPLY:
        parenthesize(tmp)
    ml.neg = tmp.neg
    if n > 0:
        ml.priority = P_EXP
        ml.type |= T_POWER
        n += 1
        ml.latex = tmp.latex + ({2: "''", 3: "'''"}.get(n) or "^{(" + str(n) + ")}")
        return ml
    if n == 0 and not vectarg:
        ml.priority = P_EXP
        ml.latex = tmp.latex + "'"
        return ml
    if vectarg and len(arg.items) > 1:
        vars_ = []
        if get_derive_vars(arg.items, vars_):
            nv = len(vars_)
            ml.priority = P_MUL
            ml.type |= T_FRACTION
            if nv == 1 and vars_[0][1] == 1:
                var = gen2markup(vars_[0][0], flags)
                ml.latex = ("\\frac{" + TEX_MATHD + tmp.latex + "}{" + TEX_MATHD + var.latex + "}" if simp
                            else "\\frac{" + TEX_MATHD + "}{" + TEX_MATHD + var.latex + "}" + tmp.latex)
                return ml
            ds_tex = "\\partial " if nv > 1 else TEX_MATHD
            parts = []
            for v, k in vars_:
                var = gen2markup(v, flags)
                parts.append(ds_tex + var.latex + ("" if k == 1 else "^{" + str(k) + "}"))
            body = TEX_ITIMES.join(parts)
            ns = str(sum(k for _, k in vars_))
            ml.latex = ("\\frac{" + ds_tex + "^{" + ns + "}" + tmp.latex + "}{" + body + "}" if simp
                        else "\\frac{" + ds_tex + "^{" + ns + "}}{" + body + "}" + tmp.latex)
            return ml
    return None


def export_latex(g):  # :3629
    ml = gen2markup(g, F_TOPLEVEL | F_ELEMPOW)
    prepend_minus(ml)
    return ml.latex
