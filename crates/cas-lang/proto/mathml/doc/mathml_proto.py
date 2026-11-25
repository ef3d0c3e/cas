"""mathml_proto.py -- Python port of giac's `mathml` (mathml.cc, live Linux #else branch :63-2032).

Tree: T5's markup.py classes (Int, Bool, Frac, Dbl, Cplx, Idnt, Str, Func, Vect, Symb).
`gprint` ports giac's print of a gen (symbolic.cc add_print*, unary.cc printsommetasoperator, usual.cc
printas*), used by the <mtext> fallbacks, integrate bounds and number texts.
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "proto"))
from markup import Int, Bool, Frac, Dbl, Cplx, Idnt, Str, Func, Vect, Symb, my_isalpha  # noqa: E402

# ---- op tables ---------------------------------------------------------------------------------

# sommet print names where giac's differ from its internal name (the `_xxx_s` strings)
PNAME = {"plus": "+", "prod": "*", "pow": "^", "neg": "-", "equal": "=", "same": "==", "different": "!=",
         "interval": "..", "inferieur_strict": "<", "inferieur_egal": "<=", "superieur_strict": ">",
         "superieur_egal": ">=", "division": "/", "derive": "diff", "ou": "or", "et": "and",
         "minus": "minus", "pointminus": ".-", "factorial": "factorial", "equal2": "%="}

# functions with a texprint (grep of define_unary_function_eval* in giac's sources, printsommet.py)
TEXPRINT = {"ampersand_times", "and", "at", "binary_minus", "complement", "composepow", "compose", "conj",
            "cross", "decrement", "derive", "deuxpoints", "divcrement", "division", "div", "dollar", "equal2",
            "equal", "et", "exp", "im", "increment", "inferieur_egal", "inferieur_strict", "integrate",
            "intersect", "interval", "leftopen_interval", "leftrightopen_interval", "limit", "lim", "minus",
            "MOD", "moyal", "multcrement", "normalmod", "NTHROOT", "of", "oufr", "ou", "pointdivision",
            "pointpow", "pointprod", "polar_complex", "program", "quote_pow", "re", "rightopen_interval",
            "same", "sqrt", "sto", "struct_dot", "superieur_egal", "superieur_strict", "symmetric_difference",
            "tilocal", "ti_semi", "union", "xor"}


def _printsommet_names():
    # printsommet.txt: printsommet.py's output on giac 2.0.0.21's source (C0009), frozen at 1.6/T11
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "printsommet.txt"), encoding="utf-8") as f:
        return {ln.split()[0] for ln in f if ln.split()[1] != "0"}


PRINTSOMMET = _printsommet_names()

ANDOR = [("and", "&and;"), ("et", "&and;"), ("ou", "&or;"), ("xor", "&veebar;"), ("union", "&cup;"),
         ("intersect", "&cap;"), ("symmetric_difference", "&Delta;")]  # :1592-1593

# ---- numbers -----------------------------------------------------------------------------------


def numval(g):
    from fractions import Fraction
    if isinstance(g, Int):
        return g.v
    if isinstance(g, Bool):
        return int(g.v)
    if isinstance(g, Dbl):
        return g.val
    if isinstance(g, Frac) and isinstance(g.num, Int) and isinstance(g.den, Int):
        return Fraction(g.num.v, g.den.v)
    return None


def is_zero(g):
    if isinstance(g, Cplx):
        return is_zero(g.re) and is_zero(g.im)
    return numval(g) == 0


def is_one(g):
    return numval(g) == 1


def is_positive(g):  # giac is_positive: symbolic -> true (unknown sign), number >= 0
    v = numval(g)
    return v is None or v >= 0


def neg_num(g):
    if isinstance(g, Int):
        return Int(-g.v)
    if isinstance(g, Bool):
        return Int(-int(g.v))
    if isinstance(g, Dbl):
        return Dbl(g.text, not g.neg)
    if isinstance(g, Frac):
        return Frac(neg_num(g.num), g.den)
    if isinstance(g, Cplx):
        return Cplx(neg_num(g.re), neg_num(g.im))
    return Symb("neg", g)


def abs_num(g):
    return neg_num(g) if not is_positive(g) else g


def gtype(g):  # giac type order: _INT_ 0 _DOUBLE_ 1 _ZINT 2 _REAL 3 _CPLX 4 _POLY 5 _IDNT 6 _VECT 7 _SYMB 8
    if isinstance(g, (Int, Bool)):
        return 0 if -2**31 <= numval(g) < 2**31 else 2
    if isinstance(g, Dbl):
        return 1
    if isinstance(g, Cplx):
        return 4
    if isinstance(g, Idnt):
        return 6
    if isinstance(g, Vect):
        return 7
    if isinstance(g, Symb):
        return 8
    if isinstance(g, Frac):
        return 10
    if isinstance(g, Str):
        return 12
    if isinstance(g, Func):
        return 13
    raise TypeError(g)


def is_int(g):
    return isinstance(g, (Int, Bool))


def need_par_func(op):  # usual.cc:8233-8242 (the _FUNC branch, also reached by a sommet)
    if op in ("pow", "division", "prod"):
        return False
    if op in ("neg", "inv", "minus", "and", "et", "ou", "oufr", "xor", "same", "equal", "equal2",
              "superieur_egal", "superieur_strict", "inferieur_egal", "inferieur_strict"):
        return True
    return op in PRINTSOMMET


def need_parenthesis(g):  # usual.cc:8218
    if is_int(g):
        return numval(g) < 0
    if isinstance(g, Cplx):
        rg, ig = neg_num(g.re), neg_num(g.im)
        if is_zero(rg):
            return numval(ig) > 0
        if is_zero(ig):
            return numval(rg) > 0
        return True
    if isinstance(g, Frac):
        return True
    if isinstance(g, Symb):
        return need_par_func(g.op)
    if isinstance(g, Func):
        return need_par_func(g.name)
    return False


def args_of(g):
    return g.arg.items if isinstance(g.arg, Vect) else [g.arg]


# ---- giac print --------------------------------------------------------------------------------


def pname(op):
    return PNAME.get(op, op)


def print_atom(g):  # gen::print gen.cc:14569
    if is_int(g):
        return "true" if g == Bool(True) else "false" if g == Bool(False) else str(g.v)
    if isinstance(g, Dbl):
        return ("-" if g.neg else "") + g.text
    if isinstance(g, Frac):
        return gprint(g.num) + "/" + gprint(g.den)
    if isinstance(g, Cplx):  # gen.cc:14606-14638
        re_, im_ = g.re, g.im
        if is_zero(im_):
            return gprint(re_)
        if is_zero(re_):
            if is_one(im_):
                return "i"
            if numval(im_) == -1:
                return "-i"
            return gprint(im_) + "*i"
        if is_one(im_):
            return gprint(re_) + "+i"
        if numval(im_) == -1:
            return gprint(re_) + "-i"
        if numval(im_) is not None and numval(im_) < 0:
            return gprint(re_) + "-" + gprint(neg_num(im_)) + "*i"
        return gprint(re_) + "+" + gprint(im_) + "*i"
    if isinstance(g, Idnt):
        return g.name
    if isinstance(g, Str):
        return '"' + g.s + '"'
    if isinstance(g, Func):
        return "'" + g.name + "'"
    if isinstance(g, Vect):  # print_VECT gen.cc:13422
        if not g.items:
            return {"seq": "seq[]", "set": "set[ ]"}.get(g.sub, "[]")
        inner = ",".join(gprint(x) for x in g.items)
        return {"seq": inner, "set": "set[" + inner + "]"}.get(g.sub, "[" + inner + "]")
    raise TypeError(g)


def gprint(g):
    s = []
    add_print(s, g)
    return "".join(s)


def _last(s):
    return s[-1][-1] if s and s[-1] else (s[-2][-1] if len(s) > 1 and s[-2] else "")


def _setlast(s, ch):
    t = "".join(s)
    s.clear()
    if t:
        s.append(t[:-1] + ch)


def add_print(s, g):  # symbolic.cc:601
    if isinstance(g, Symb):
        add_print_symbolic(s, g)
        return
    if isinstance(g, Int) and gtype(g) == 0:  # add_print_int :512 (an _INT_; a _ZINT prints as below)
        if g.v < 0:
            t = "".join(s)
            if t.endswith("+"):
                _setlast(s, "-")
            elif t.endswith("-"):
                if len(t) == 1:
                    s.clear()
                else:
                    _setlast(s, "+")
            else:
                s.append("-")
        s.append(str(abs(g.v)))
        return
    if isinstance(g, Vect) and g.sub == 'list':
        s.append("[" + ",".join(gprint(x) for x in g.items) + "]")
        return
    if isinstance(g, Frac) and is_int(g.num) and is_int(g.den):
        add_print(s, g.num)
        s.append("/")
        add_print(s, g.den)
        return
    tmp = print_atom(g)
    t = "".join(s)
    if t.endswith("+") and tmp.startswith("-"):
        s.clear()
        s.append(t[:-1] + tmp)
        return
    if t.endswith("-") and tmp.startswith("-"):
        s.clear()
        s.append(t[:-1] + "+" + tmp[1:])
        return
    s.append(tmp)


def add_print_symbolic(s, g):  # symbolic.cc:542
    op, f = g.op, g.arg
    if op in PRINTSOMMET and op not in ("plus", "prod", "pow"):
        s.append(printsommet(op, f))
        return
    if isinstance(f, Vect) and not f.items:
        s.append(pname(op) + {"list": "([])", "set": "(set[ ])"}.get(f.sub, "(NULL)"))
        return
    if op == "prod":
        return add_print_prod(s, g)
    if op == "plus":
        return add_print_plus(s, g)
    if op == "pow":
        return add_print_pow(s, g)
    if op == "neg":
        return add_print_neg(s, f)
    if op == "inv":
        return add_print_inv(s, f)
    if not isinstance(f, Vect) or isinstance(f.items[0], Vect):
        s.append(pname(op) + "(")
        add_print(s, f)
        s.append(")")
        return
    s.append(pname(op) + "(")
    if f.sub != 'seq':
        s.append({"set": "set["}.get(f.sub, "["))
    for k, x in enumerate(f.items):
        if k:
            s.append(",")
        add_print(s, x)
    if f.sub != 'seq':
        s.append("]")
    s.append(")")


def add_print_plus(s, g):  # symbolic.cc:175
    f = g.arg
    if not isinstance(f, Vect):
        s.append("+")
        return add_print(s, f)
    v = f.items
    if not v:
        return
    if len(v) == 1:
        s.append("+")
    add_print(s, v[0])
    for x in v[1:]:
        s.append("+")
        par = isinstance(x, Symb) and x.op in ("superieur_egal", "superieur_strict", "inferieur_strict",
                                                "inferieur_egal")
        if par:
            s.append("(")
        add_print(s, x)
        if par:
            s.append(")")


def add_print_neg(s, f):  # symbolic.cc:209
    t = "".join(s)
    if t.endswith("+"):
        _setlast(s, "-")
    elif t.endswith("-"):
        _setlast(s, "+")
    else:
        s.append("-")
    if not isinstance(f, Cplx):
        if not isinstance(f, Symb) or f.op in ("inv", "prod") or not need_parenthesis(f):
            return add_print(s, f)
    s.append("(")
    add_print(s, f)
    s.append(")")


def add_print_inv(s, f):  # symbolic.cc:231
    isneg = False
    if isinstance(f, Symb) and f.op == "neg" and not isinstance(f.arg, Vect):
        f, isneg = f.arg, True
    if gtype(f) < 4 and numval(f) is not None and numval(f) < 0 or (gtype(f) < 4 and numval(f) == 0):
        f, isneg = neg_num(f), not isneg
    t = "".join(s)
    if isneg and t.endswith("-"):
        if len(t) == 1:
            s.clear()
        else:
            _setlast(s, "+")
        isneg = False
    t = "".join(s)
    if isneg and t.endswith("+"):
        _setlast(s, "-")
        isneg = False
    bt = isinstance(f, Symb) and (f.op in ("plus", "inv", "prod") or need_par_func(f.op))
    if not bt and not isinstance(f, Cplx):
        s.append("-1/" if isneg else "1/")
        return add_print(s, f)
    s.append("-1" if isneg else "1/(")  # giac's "-1" (no "/(") kept
    add_print(s, f)
    s.append(")")


def is_inequation(g):
    return isinstance(g, Symb) and g.op in ("superieur_egal", "superieur_strict", "inferieur_strict",
                                            "inferieur_egal", "different", "same", "equal")


def add_print_prod(s, g):  # symbolic.cc:268
    f = g.arg
    if not isinstance(f, Vect):
        return add_print(s, f)
    for i, e in enumerate(f.items):
        if not isinstance(e, Symb):
            if i:
                s.append("*")
            if isinstance(e, Cplx) and need_parenthesis(e):
                s.append("(" + gprint(e) + ")")
            else:
                add_print(s, e)
        elif e.op == "inv":
            ff = e.arg
            if i:
                if isinstance(ff, Cplx) or (isinstance(ff, Symb) and (
                        ff.op in ("plus", "prod", "inv") or need_par_func(ff.op)
                        or (ff.op == "pow" and isinstance(args_of(ff)[0], Symb) and args_of(ff)[0].op == "inv"))):
                    s.append("/(")
                    add_print(s, ff)
                    s.append(")")
                else:
                    s.append("/")
                    add_print(s, ff)
            else:
                add_print(s, e)
        else:
            if i:
                s.append("*")
            if e.op in ("plus", "neg") or is_inequation(e):
                s.append("(")
                add_print(s, e)
                s.append(")")
            else:
                add_print(s, e)


def add_print_pow(s, g):  # symbolic.cc:344
    f = g.arg
    if not isinstance(f, Vect) or len(f.items) != 2:
        s.append("pow(")
        add_print(s, f)
        s.append(")")
        return
    arg, pui = f.items
    argpar = ((gtype(arg) > 4) or not is_positive(arg)) and not isinstance(arg, Idnt)
    if gtype(pui) > 3 and pui == Frac(Int(1), Int(2)):
        s.append("sqrt(")
        add_print(s, arg)
        s.append(")")
        return
    if gtype(pui) > 3 and pui == Frac(Int(-1), Int(2)):
        s.append("1/sqrt(")
        add_print(s, arg)
        s.append(")")
        return
    if isinstance(arg, Idnt) or (isinstance(arg, Symb) and arg.op != "neg" and arg.op not in PRINTSOMMET):
        argpar = isinstance(arg, Symb) and arg.op == "inv"
        if argpar:
            s.append("(")
        if isinstance(pui, (Symb, Frac, Cplx)):
            add_print(s, arg)
            if argpar:
                s.append(")")
            s.append("^(")
            add_print(s, pui)
            s.append(")")
            return
        add_print(s, arg)
        if argpar:
            s.append(")")
        s.append("^")
        return add_print(s, pui)
    if argpar:
        s.append("(")
    add_print(s, arg)
    if argpar:
        s.append(")")
    s.append("^")
    puipar = isinstance(pui, (Symb, Frac, Cplx)) or (isinstance(pui, Vect) and pui.sub == 'seq')
    if puipar:
        s.append("(")
    add_print(s, pui)
    if puipar:
        s.append(")")


def printsommetasoperator(f, sop):  # unary.cc:599
    if not isinstance(f, Vect):
        return gprint(f)
    if (my_isalpha(sop[0]) and ord(sop[0]) > 32) or sop[0] in "%.":
        sop = " " + sop + " "
    v = f.items
    if not v:
        return ""
    a = v[0]
    if isinstance(a, Frac) and sop != "=":
        s = "(" + gprint(a) + ")"
    elif sop == "=" or isinstance(a, Idnt) or (gtype(a) <= 4 and is_positive(a)):
        s = gprint(a)
    else:
        s = "(" + gprint(a) + ")"
    for b in v[1:]:
        if (isinstance(b, (Symb, Frac, Cplx)) or (isinstance(b, Vect) and b.sub == 'seq')) and sop != " . ":
            s += sop + "(" + gprint(b) + ")"
        else:
            s += sop + gprint(b)
    return s


def printsommet(op, f):
    if op == "same":  # usual.cc:6507
        return "(" + printsommetasoperator(f, "==") + ")"
    if op == "different":  # usual.cc:5711
        return printsommetasoperator(f, "<>")
    if op in ("equal", "interval", "inferieur_strict", "inferieur_egal", "superieur_strict", "superieur_egal",
              "minus", "union", "intersect", "xor", "equal2"):
        return printsommetasoperator(f, pname(op))
    if op == "of":  # printasof_ usual.cc
        a = args_of(Symb(op, f))
        return gprint(a[0]) + "(" + gprint(a[1]) + ")"
    if op == "at":
        a = args_of(Symb(op, f))
        return gprint(a[0]) + "[" + gprint(a[1]) + "]"
    if op == "factorial":  # usual.cc:7721
        if isinstance(f, Idnt) or (gtype(f) <= 1 and is_positive(f)):
            return gprint(f) + "!"
        return "(" + gprint(f) + ")!"
    if op == "sqrt":  # usual.cc:1396
        return "sqrt(" + gprint(f) + ")"
    return pname(op) + "(" + gprint(f) + ")"  # other printsommets: not ported


# ---- aplatir_fois_plus gen.cc:12505 ------------------------------------------------------------


def aplatir_plus(g):  # gen.cc:12474
    if not (isinstance(g, Symb) and g.op == "plus" and isinstance(g.arg, Vect)):
        return g
    items = g.arg.items
    if not items:
        return Int(0)
    v = []
    it, itend = 0, len(items)
    while True:
        for k in range(itend - 1, it, -1):
            v.append(aplatir(items[k]))
        first = items[it]
        if not (isinstance(first, Symb) and first.op == "plus"):
            v.append(aplatir(first))
            break
        if not isinstance(first.arg, Vect):
            v.append(first)
            break
        items, it, itend = first.arg.items, 0, len(first.arg.items)
    v.reverse()
    return Symb("plus", Vect(v, 'seq'))


def rat(q):
    q = Fraction(q)
    return Int(q.numerator) if q.denominator == 1 else Frac(Int(q.numerator), Int(q.denominator))


def frac2symb(f):  # gen.cc:1952 _FRAC2_SYMB
    if is_one(f.num):
        return Symb("inv", f.den)
    if numval(f.num) == -1:
        return Symb("neg", Symb("inv", f.den))
    return Symb("prod", Vect([f.num, Symb("inv", f.den)], 'seq'))


def sym_mult(a, b):  # gen.cc:7228, the branches a printed tree reaches
    if isinstance(a, Symb) and a.op == "inv" and a.arg == b or isinstance(b, Symb) and b.op == "inv" and b.arg == a:
        return Int(1)
    for x in (a, b):
        if isinstance(x, (Int, Bool, Dbl)) and is_zero(x):
            return x
    if is_one(a):
        return b
    if is_one(b):
        return a
    if isinstance(a, Symb) and a.op == "neg":
        if isinstance(b, Symb) and b.op == "neg":
            return sym_mult(a.arg, b.arg)
        return Symb("neg", sym_mult(a.arg, b))
    if isinstance(b, Symb) and b.op == "neg":
        return Symb("neg", sym_mult(a, b.arg))
    exact = lambda x: isinstance(x, (Int, Bool, Frac))  # noqa: E731
    inv_int = lambda x: isinstance(x, Symb) and x.op == "inv" and isinstance(x.arg, Int)  # noqa: E731
    if isinstance(a, Frac):  # :7462-7469
        if inv_int(b):
            return rat(numval(a) / b.arg.v)
        if exact(b):
            return rat(numval(a) * numval(b))
        if not isinstance(b, (Symb, Idnt)):
            return Symb("prod", Vect([a, b], 'seq'))  # numeric product with a double: not ported
        return sym_mult(frac2symb(a), b)
    if isinstance(b, Frac):  # :7470-7478
        if inv_int(a):
            return rat(numval(b) / a.arg.v)
        if exact(a):
            return rat(numval(a) * numval(b))
        if not isinstance(a, (Symb, Idnt)):
            return Symb("prod", Vect([a, b], 'seq'))
        return sym_mult(a, frac2symb(b))
    if isinstance(a, Int) and inv_int(b):  # :7479-7486 (gaussian integers not ported)
        return rat(Fraction(a.v, b.arg.v))
    if isinstance(b, Int) and inv_int(a):
        return rat(Fraction(b.v, a.arg.v))
    if gtype(a) <= 3 and numval(a) < 0:  # :7487-7490 (the _INT_ < 0 branches after are not reached)
        return Symb("neg", sym_mult(neg_num(a), b))
    if gtype(b) <= 3 and numval(b) < 0:
        return Symb("neg", sym_mult(a, neg_num(b)))
    if isinstance(a, Symb) and isinstance(b, Symb):
        if a.op == "prod" and b.op == "prod":
            return Symb("prod", Vect(args_of(a) + args_of(b), 'seq'))
        if a.op == "prod":
            return Symb("prod", Vect(args_of(a) + [b], 'seq'))
        if b.op == "prod":
            return Symb("prod", Vect([a] + args_of(b), 'seq'))
        return Symb("prod", Vect([a, b], 'seq'))
    if isinstance(b, Symb):
        if b.op == "prod":
            return Symb("prod", Vect([a] + args_of(b), 'seq'))
        return Symb("prod", Vect([a, b], 'seq'))
    if isinstance(a, Symb):
        if a.op == "prod":
            return Symb("prod", Vect(args_of(a) + [b], 'seq'))
        return Symb("prod", Vect([b, a] if gtype(b) == 0 else [a, b], 'seq'))
    return Symb("prod", Vect([a, b], 'seq'))


def aplatir(g):  # gen.cc:12505
    if isinstance(g, Vect):
        return Vect([aplatir(x) for x in g.items], g.sub)
    if not isinstance(g, Symb):
        return g
    if g.op == "pow" and isinstance(g.arg, Vect) and g.arg.sub != 'seq':
        return Symb("pow", Vect(g.arg.items, 'seq'))
    if g.op == "plus":
        return aplatir_plus(g)
    f = g.arg
    if g.op == "prod" and isinstance(f, Vect) and len(f.items) == 2:
        return sym_mult(aplatir(f.items[0]), aplatir(f.items[1]))
    return Symb(g.op, aplatir(f))


# ---- mathml.cc ---------------------------------------------------------------------------------


def string2mathml(m):  # :111
    for a, b in (("&", "&amp;"), ("<", "&lt;"), (">", "&gt;"), ("\n", "</mi></mtd></mtr><mtr><mtd><mi>")):
        m = m.replace(a, b)
    return '<mtable columnalign="left"><mtr><mtd><mi>' + m + "</mi></mtd></mtr></mtable>"


def ckmatrix(items):  # vecteur.cc ckmatrix(m, false)
    if not items or not all(isinstance(r, Vect) for r in items):
        return False
    n = len(items[0].items)
    return n > 0 and all(len(r.items) == n and not any(isinstance(x, Vect) for x in r.items) for r in items)


def matrix2mathml(rows):  # :200
    s = "<mrow><mo>(</mo><mtable>"
    for r in rows:
        s += "<mtr>" + "".join("<mtd>" + gen2mathml(x) + "</mtd>" for x in r.items) + "</mtr>"
    return s + "</mtable><mo>)</mo></mrow>"


def vect2mathml(v, sub):  # :221
    begin, end = {"seq": ("", ""), "set": ("set[", "]"), "poly1": ("poly1[", "]")}.get(sub, ("[", "]"))
    return begin + "<mo>,</mo>".join(gen2mathml(x) for x in v) + end


def needs_times(g):  # :256 (eval taken as the identity)
    return gtype(g) <= 3


def is_numtype(g):
    return gtype(g) <= 3


def pow_times(it):  # :292-299 / :338-345 test
    if isinstance(it, Symb) and it.op == "pow":
        a = args_of(it)
        if numval(a[-1]) != Fraction(1, 2) and is_int(a[0]) and is_positive(a[0]):
            return True
    return False


def prod_vect2mathml(v):  # :274
    if not v:
        return "<mn>1</mn>"
    s = ""
    prec_int = False
    for k, it in enumerate(v):
        if isinstance(it, Symb) and it.op in ("plus", "neg"):
            s += " <mo>(</mo>" + gen2mathml(it) + " <mo>)</mo>"
        elif isinstance(it, Cplx):
            if not is_zero(it.re) and not is_zero(it.im):
                s += "<mo>(</mo>" + gen2mathml(it) + "<mo>)</mo>"
            else:
                s += gen2mathml(it)
        else:
            if prec_int and pow_times(it):
                if not s.endswith("<mo>&times;</mo>"):
                    s += "<mo>&times;</mo>"
            s += gen2mathml(it)
        prec_int = is_numtype(it)
        if k + 1 == len(v):
            return s
        if needs_times(v[k + 1]) and prec_int:
            s += "<mo>&times;</mo>"


def prod_vect2mathml_no_bra(v):  # :321
    if not v:
        return "<mn>1</mn>"
    if len(v) == 1:
        return gen2mathml(v[0])
    s = ""
    prec_int = False
    for k, it in enumerate(v):
        if isinstance(it, Symb) and it.op in ("plus", "neg"):
            s += " <mo>(</mo>" + gen2mathml(it) + " <mo>)</mo>"
        elif isinstance(it, Cplx):
            if not is_zero(it.re) and not is_zero(it.im):
                s += "<mo>(</mo>" + gen2mathml(it) + "<mo>)</mo>"
            else:
                s += gen2mathml(it)
        elif prec_int and isinstance(it, Symb) and it.op == "pow":
            if pow_times(it):
                s += "<mo>&times;</mo>"
            s += gen2mathml(it)
        else:
            s += gen2mathml(it)
        prec_int = is_numtype(it)
        if k + 1 == len(v):
            return s
        if needs_times(v[k + 1]) and prec_int:
            s += "<mo>&times;</mo>"


def printsommetasoperator_mml(f, sop):  # :369
    if not isinstance(f, Vect):
        return "<mi>" + gprint(f) + "</mi>"
    v = f.items
    if not v:
        return ""
    a = v[0]
    eq = sop in ("<mo>=</mo>", "<mo>≈</mo>")
    if isinstance(a, Frac):
        s = gen2mathml(a)
    elif not need_parenthesis(a) or (isinstance(a, Symb) and a.op == "of") or eq:
        s = gen2mathml(a)
    else:
        s = "<mo>(</mo>" + gen2mathml(a) + "<mo>)</mo>"
    for b in v[1:]:
        if not isinstance(b, Symb) or eq:
            s += sop + gen2mathml(b)
        else:
            s += sop + "<mo>(</mo>" + gen2mathml(b) + "<mo>)</mo>"
    return s


def is_equal(g):
    return isinstance(g, Symb) and g.op == "equal"


def equal2arrow(g):  # :448
    if not is_equal(g):
        return gen2mathml(g)
    v = args_of(g)
    return "<mrow>" + gen2mathml(v[0]) + "<mo>&rarr;</mo>" + gen2mathml(v[1]) + "</mrow>"


def printaslimit(f):  # :457
    s = "<mo>lim</mo>"
    if not isinstance(f, Vect):
        return s + gen2mathml(f)
    v = f.items
    if not v:
        return s
    if len(v) == 1:
        return s + gen2mathml(v[0])
    if len(v) == 2:
        return "<munder>" + s + equal2arrow(v[1]) + "</munder>" + gen2mathml(v[0])
    if len(v) == 3:
        if is_one(v[2]):
            return ("<munder>" + s + "<mrow><msup>" + equal2arrow(v[1]) + "<mo>+</mo></msup></mrow></munder>"
                    + gen2mathml(v[0]))
        if numval(v[2]) == -1:
            return ("<munder>" + s + "<mrow><msup>" + equal2arrow(v[1]) + "<mo>-</mo></msup></mrow></munder>"
                    + gen2mathml(v[0]))
        return "<munder>" + s + equal2arrow(v[1]) + "</munder>" + gen2mathml(v[0])
    return s


def printasintegrate(f):  # :485
    s = "<mo>&int;</mo>"
    if not isinstance(f, Vect):
        return s + gen2mathml(f)
    v = f.items
    if not v:
        return s
    if len(v) == 1:
        return s + gen2mathml(v[0])
    if len(v) == 2:
        return s + gen2mathml(v[0]) + "<mi>d</mi>" + gen2mathml(v[1])
    if len(v) == 4:
        return ("<msubsup>" + s + "<mi>" + gprint(v[2]) + "</mi><mi>" + gprint(v[3]) + "</mi></msubsup>"
                + gen2mathml(v[0]) + "<mi>d</mi><mi>" + gprint(v[1]) + "</mi>")
    return s


def printasderive(f):  # :509
    if not isinstance(f, Vect) or len(f.items) < 2:
        if isinstance(f, Idnt) or gtype(f) <= 4 or (isinstance(f, Symb) and f.op == "derive"):
            return gen2mathml(f) + "'"
        return "<msup><mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow><mi>'</mi></msup>"
    v = f.items
    par = isinstance(v[0], Symb)
    lp, rp = ("<mo>(</mo>", "<mo>)</mo>") if par else ("", "")
    if len(v) > 2:
        return ("<mfrac><mrow><msup><mo>&part;</mo><mrow>" + gen2mathml(v[2]) + "</mrow></msup>" + lp
                + gen2mathml(v[0]) + rp + "</mrow><mrow><mo>&part;</mo><msup>" + gen2mathml(v[1]) + "<mrow>"
                + gen2mathml(v[2]) + "</mrow></msup></mrow></mfrac>")
    return ("<mfrac><mrow><mo>&part;</mo>" + lp + gen2mathml(v[0]) + rp + "</mrow><mrow><mo>&part;</mo>"
            + gen2mathml(v[1]) + "</mrow></mfrac>")


def printassto(f, sop):  # :433
    if not isinstance(f, Vect) or len(f.items) != 2:
        return sop + "<mo>(</mo>" + gen2mathml(f) + "<mo>)</mo>"
    s = gen2mathml(f.items[1]) + "<mi>:=</mi>"
    if gtype(f.items[0]) == 1 and not isinstance(f.items[0], Dbl) or gtype(f.items[0]) == 2:  # type==_SEQ__VECT(=1)
        return s + "<mo>(</mo>" + gen2mathml(f.items[0]) + "<mo>)</mo>"
    return s + gen2mathml(f.items[0])


def mathml_print(g):  # :526
    op, f = g.op, g.arg
    if op in ("equal", "equal2"):
        chk = args_of(g)[1] if isinstance(f, Vect) and len(f.items) > 1 else Int(0)
        eq = not isinstance(chk, Dbl) and not (isinstance(chk, Cplx) and (isinstance(chk.re, Dbl) or isinstance(chk.im, Dbl)))
        return printsommetasoperator_mml(f, "<mo>=</mo>" if eq else "<mo>≈</mo>")
    tab = {"different": "<mo>≠</mo>", "inferieur_egal": "<mo>&le;</mo>", "superieur_egal": "<mo>&ge;</mo>",
           "inferieur_strict": "<mo>&lt;</mo>", "superieur_strict": "<mo>&gt;</mo>"}
    if op in tab:
        return printsommetasoperator_mml(f, tab[op])
    if op == "re":
        return "<mi>Re</mi><mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow>"
    if op == "im":
        return "<mi>Im</mi><mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow>"
    if op == "exp":
        return "<msup><mi>e</mi><mrow>" + gen2mathml(f) + "</mrow></msup>"
    if op == "sqrt":
        return "<msqrt>" + gen2mathml(f) + "</msqrt>"
    if op == "integrate":
        return printasintegrate(f)
    if op == "derive":
        return printasderive(f)
    if op == "sto":
        return printassto(f, "<mi>:=</mi>")
    if op == "limit":
        return printaslimit(f)
    s = gprint(g).replace("<", "&lt;").replace(">", "&gt;").replace("\n", " ")
    return "<mtext>" + s + "</mtext>"


def printassum(g):  # :600
    v = args_of(g)
    if len(v) != 4:
        return "<mtext>" + gprint(g) + "</mtext>"
    return ("<munderover><mi>&Sigma;</mi><mrow>" + gen2mathml(v[1]) + "<mo>=</mo>" + gen2mathml(v[2])
            + "</mrow><mrow>" + gen2mathml(v[3]) + "</mrow></munderover>" + gen2mathml(v[0]))


def idnt2mathml_(s0):  # :1808
    if not s0:
        return s0
    n = len(s0)
    j = n - 1
    while j >= 2:
        if s0[j].isascii() and s0[j].isalpha():  # giac's signed char: a non-ASCII byte is not > 32
            break
        j -= 1
    s, sadd = s0[:j + 1], (s0[j + 1:] if j < n - 1 else "")
    two = {"mu", "nu", "pi", "xi", "Xi"}
    if len(s) == 2:
        if s in two:
            return "&" + s + ";" + sadd
        special = {"im": "&Im;", "re": "&Re;", "RR": "ℝ", "QQ": "ℚ", "ZZ": "ℤ", "CC": "ℂ", "NN": "ℕ", "HH": "ℍ",
                   "PP": "ℙ"}
        if s in special:
            return special[s] + sadd
    greek = {3: {"chi", "phi", "Phi", "eta", "rho", "tau", "psi", "Psi"}, 4: {"beta", "zeta"},
             5: {"alpha", "delta", "Delta", "gamma", "Gamma", "kappa", "theta", "Theta", "sigma", "Sigma",
                 "Omega", "omega", "aleph", "nabla"}, 6: {"lambda", "Lambda", "approx"}, 7: {"epsilon", "product"}}
    if s in greek.get(len(s), ()):
        return "&" + s + ";" + sadd
    return s0


def idnt2mathml(s0):  # :1870
    n = len(s0)
    for j in range(n - 1, 0, -1):
        if j < n - 1 and s0[j] == "_":
            return "<msub><mi>" + idnt2mathml_(s0[:j]) + "</mi><mn>" + s0[j + 1:] + "</mn></msub>"
    return "<mi>" + idnt2mathml_(s0) + "</mi>"


PER_NODE = True  # giac flattens at each symbolic (:1929); False: once at the top, as the Rust port does


def gen2mathml(e):  # :1879
    if isinstance(e, Symb) and e.op == "sum":
        return printassum(e)
    if isinstance(e, Symb) and e.op == "abs":
        return "<mo>&VerticalBar;</mo>" + gen2mathml(e.arg) + "<mo>&VerticalBar;</mo>"
    if isinstance(e, (Int, Bool)):
        return "<mn>" + print_atom(e) + "</mn>"
    if isinstance(e, Dbl):
        return "<mn>" + print_atom(e) + "</mn>"
    if isinstance(e, Cplx):  # :1899
        re_, im_ = e.re, e.im
        part_re, part_im = "", "<mi>i</mi>"
        if not is_zero(re_):
            part_re = "<mn>" + gprint(re_) + "</mn>"
        if not is_zero(im_):
            if is_positive(im_):
                if not is_zero(re_):
                    part_re += "<mo>+</mo>"
            else:
                part_re += "<mo>-</mo>"
        if not is_one(neg_num(im_)) and not is_one(im_):
            part_im = "<mn>" + gprint(abs_num(im_)) + "</mn>" + part_im
        if is_zero(im_):
            part_im = ""
            if is_zero(re_):
                part_re = "<mn>0.0</mn>"
        return part_re + part_im
    if isinstance(e, Idnt):
        if e.name == "infinity":
            return "<mn>&infin;</mn>"
        if e.name == "pi":
            return "<mi>&pi;</mi>"
        if e.name == "undef":
            return "<mi>undef</mi>"
        return idnt2mathml(e.name)
    if isinstance(e, Symb):
        t = aplatir(e) if PER_NODE else e
        if not isinstance(t, Symb):
            return gen2mathml(t)
        return symbolic2mathml(t)
    if isinstance(e, Vect):
        if e.sub != 'seq' and ckmatrix(e.items):
            return matrix2mathml(e.items)
        return vect2mathml(e.items, e.sub)
    if isinstance(e, Frac):
        return "<mfrac><mrow>" + gen2mathml(e.num) + "</mrow><mrow>" + gen2mathml(e.den) + "</mrow></mfrac>"
    if isinstance(e, Str):
        return string2mathml(e.s)
    if isinstance(e, Func):
        return "<mi>" + print_atom(e) + "</mi>"
    raise TypeError(e)


def symbolic2mathml(g):  # :1571
    op, f = g.op, g.arg
    opstring = pname(op)
    if op in ("of", "at"):
        a = args_of(g)
        if op == "of":
            return gen2mathml(a[0]) + "<mo>(</mo>" + gen2mathml(a[1]) + "<mo>)</mo>"
        return "<msub><mrow>" + gen2mathml(a[0]) + "</mrow><mrow>" + gen2mathml(a[1]) + "</mrow></msub>"
    if op == "program" and isinstance(f, Vect) and len(f.items) >= 3:
        if not (isinstance(f.items[2], Symb) and f.items[2].op == "local"):
            return gen2mathml(f.items[0]) + "<mo>→</mo>" + gen2mathml(f.items[2])
    if op == "complement":
        return gen2mathml(f) + "<mo>&#x2201;</mo>"
    andor = next((m for o, m in ANDOR if o == op), None)
    if andor is None and opstring != "/" and (op in TEXPRINT or op == "different"):
        return mathml_print(g)
    if op == "member" and isinstance(f, Vect) and len(f.items) == 2:
        return gen2mathml(f.items[0]) + "<mtext>∈</mtext>" + gen2mathml(f.items[1])
    if op == "vector":
        s = '<mover><mrow accent="true">'
        s += gen2mathml(f) if not isinstance(f, Vect) else "".join(gen2mathml(x) for x in f.items)
        return s + '</mrow><mo mathsize="100%">&rarr;</mo></mover>'
    # rootof with a value (:1625) not ported: needs giac's evaluation
    if isinstance(f, Vect) and not f.items:
        return "<mi>" + opstring + "()</mi>"
    if not isinstance(f, Vect) or isinstance(f.items[0], Vect):
        if op == "factorial":
            if isinstance(f, Symb):
                return "<mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow><mtext>!</mtext>"
            return gen2mathml(f) + "<mtext>!</mtext>"
        if op in ("neg", "plus"):
            if not isinstance(f, Symb) or (f.op in ("at", "of") and isinstance(args_of(f)[0], Idnt)):
                return "<mo>" + opstring + "</mo>" + gen2mathml(f)
            if f.op in ("inv", "pow"):
                return "<mo>" + opstring + "</mo>" + gen2mathml(f)
            return "<mo>" + opstring + "</mo><mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow>"
        if op == "inv":
            return "<mfrac><mrow><mn>1</mn></mrow><mrow>" + gen2mathml(f) + "</mrow></mfrac>"
        if op == "pow":
            return ("<msup><mrow>" + gen2mathml(f.items[0]) + "</mrow><mrow>" + gen2mathml(f.items[1])
                    + "</mrow></msup>")
        return "<mi>" + opstring + "</mi><mrow><mo>(</mo>" + gen2mathml(f) + "<mo>)</mo></mrow>"
    v = f.items
    l = len(v)
    if andor is not None and l:
        s = ""
        for i, x in enumerate(v):
            par = isinstance(x, Symb) and x.op != "complement" and x.op != op
            s += "<mrow>" + ("(" if par else "") + gen2mathml(x) + (")" if par else "") + "</mrow>"
            if i != l - 1:
                s += "<mo> " + andor + " </mo>"
        return s
    if op == "plus":
        s = ""
        for i, x in enumerate(v):
            if isinstance(x, Symb) and x.op == "neg":
                if isinstance(x.arg, Symb) and x.arg.op == "plus":
                    s += "<mo>-</mo><mrow><mo>(</mo>" + gen2mathml(x.arg) + "<mo>)</mo></mrow>"
                else:
                    s += "<mo>-</mo><mrow>" + gen2mathml(x.arg) + "</mrow>"
            elif is_int(x) and not is_positive(x):
                s += "<mo>-</mo><mn>" + gprint(neg_num(x)) + "</mn>"
            else:
                adds = gen2mathml(x)
                if i and not adds.startswith("<mo>-"):
                    s += "<mo>" + opstring + "</mo>"
                s += adds
        return s
    if op == "prod":
        s = ""
        num, den = [], []
        for x in v:
            if isinstance(x, Symb) and x.op == "inv":
                den.append(x.arg)
            else:
                if den:
                    s += ("<mfrac><mrow>" + prod_vect2mathml_no_bra(num) + "</mrow><mrow>"
                          + prod_vect2mathml_no_bra(den) + "</mrow></mfrac>")
                    s += "<mo>&times;</mo>"
                    num, den = [], []
                num.append(x)
        if not den:
            return s + prod_vect2mathml(num)
        return (s + "<mfrac><mrow>" + prod_vect2mathml_no_bra(num) + "</mrow><mrow>"
                + prod_vect2mathml_no_bra(den) + "</mrow></mfrac>")
    if op == "pow":
        base, ex = v[0], v[-1]
        if numval(ex) == Fraction(1, 2):  # gen == compares numbers across types
            return "<msqrt><mrow>" + gen2mathml(base) + "</mrow></msqrt>"
        if numval(ex) == Fraction(-1, 2):
            return "<mfrac><mn>1</mn><msqrt>" + gen2mathml(base) + "</msqrt></mfrac>"
        s0, s1 = gen2mathml(v[0]), gen2mathml(v[1])
        s_bra = '<msup><mfenced open="(" close=")"><mrow>(' + s0 + ")</mrow></mfenced><mrow>" + s1 + "</mrow></msup>"
        s_no_bra = "<msup><mrow> " + s0 + "</mrow><mrow>" + s1 + "</mrow></msup>"
        if isinstance(base, Symb):
            mf = base.arg
            if isinstance(mf, Vect) and not mf.items:
                return s_bra
            if not isinstance(mf, Vect) or isinstance(mf.items[0], Vect):
                return s_bra
            return s_bra if len(mf.items) > 1 else s_no_bra
        if isinstance(base, Frac):
            return s_bra
        if isinstance(base, Cplx):
            return s_no_bra if is_zero(base.im) else s_bra
        if is_int(base):
            return s_no_bra if is_positive(base) else s_bra
        return s_no_bra
    if opstring == "/":
        return "<mfrac><mrow>" + gen2mathml(v[0]) + "</mrow><mrow>" + gen2mathml(v[1]) + "</mrow></mfrac>"
    return "<mi>" + opstring + "</mi><mrow><mo>(</mo>" + ",".join(gen2mathml(x) for x in v) + "<mo>)</mo></mrow>"


PREAMBLE = ('<?xml version="1.0" encoding="iso-8859-1"?>\n<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1 plus '
            'MathML 2.0//EN"\n"http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd" [\n<!ENTITY mathml '
            '"http://www.w3.org/1998/Math/MathML">\n]>\n<html xmlns="http://www.w3.org/1999/xhtml">\n<body>\n')  # :69
END = "\n</body>\n</html>"  # :70


def ingen2mathml(g, html5):  # :1978
    if html5:
        return ('\n<math display="inline" xmlns="http://www.w3.org/1998/Math/MathML">\n\n' + gen2mathml(g)
                + "\n\n</math>\n")
    return '\n<math mode="display" xmlns="http://www.w3.org/1998/Math/MathML">\n\n' + gen2mathml(g) + "\n\n</math><br/>\n"


def mathml(g):  # _mathml :1987 (one argument)
    if not PER_NODE:
        g = aplatir(g)
    if isinstance(g, Vect) and len(g.items) == 2 and isinstance(g.items[1], Int) and gtype(g.items[1]) == 0:
        return ingen2mathml(g.items[0], g.items[1].v != 0)
    if isinstance(g, Vect) and g.sub == 'seq' and len(g.items) > 2 and isinstance(g.items[-1], Int):
        return ingen2mathml(Vect(g.items[:-1], g.sub), g.items[-1].v != 0)
    return PREAMBLE + ingen2mathml(g, False) + END + "\n"
