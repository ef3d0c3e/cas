"""Parser of giac's printed form into the markup.py tree, with giac's internal shapes.

mode 'eval' (an evaluated answer): n/m ints -> Frac; 1/b -> inv(b); a/b -> prod(a,inv(b)); numeric
  complex folding (2*i -> Cplx(0,2), 1+Cplx(0,2) -> Cplx(1,2)); Cplx/int stays prod(Cplx,inv(int))
  (x*i/2 needs it, see NOTES.md).
mode 'quote' (a quote(...) answer): a/b -> division(a,b) (giac's parser shape before eval).
Both: '/' binds tighter than '*', unary minus tighter than '*' and looser than '^' (giac prints
quote(-x*y) as (-x)*y and quote(2*x/3) latex is 2 \\frac{x}{3}); a-b -> plus(a,neg(b)); f(x) for an
unknown f -> of(f,x); x[1] -> at(x,1); set[..], poly1[..] vect subtypes.
"""
import re
from markup import Int, Bool, Frac, Dbl, Cplx, Idnt, Str, Vect, Symb, I, numval

BUILTINS = {"sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan", "acot", "asec", "acsc",
            "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "ln", "exp", "sqrt", "abs", "floor",
            "ceil", "erf", "Gamma", "Psi", "Dirac", "Heaviside", "Zeta", "Beta", "rootof", "integrate",
            "int", "sum", "product", "limit", "diff", "derive", "factorial", "binomial", "sto", "inv"}
RENAME = {"diff": "derive", "log": "ln"}

TOK = re.compile(r"""\s*(?:
    (?P<num>\d+\.\d*(?:[eE][+-]?\d+)?|\d+(?:[eE][+-]?\d+)?(?!\.\.)|\.\d+(?:[eE][+-]?\d+)?)
  | (?P<str>"(?:\\.|[^"\\])*")
  | (?P<id>[A-Za-z_][A-Za-z0-9_]*)
  | (?P<op>:=|==|<>|<=|>=|\.\.|\.-|[-+*/^()\[\],=<>!'])
  )""", re.VERBOSE)


def tokenize(s):
    out, pos = [], 0
    s = s.rstrip()
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m or m.end() == pos:
            raise SyntaxError(f"bad char at {pos}: {s[pos:pos+10]!r}")
        pos = m.end()
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
    return out


def is_num(g):
    return isinstance(g, (Int, Dbl, Frac))


def negate(g):
    """unary minus with numeric folding (giac's parser/eval make -3 an integer)."""
    if isinstance(g, Int):
        return Int(-g.v)
    if isinstance(g, Dbl):
        return Dbl(g.text, not g.neg)
    if isinstance(g, Frac) and isinstance(g.num, Int):
        return Frac(Int(-g.num.v), g.den)
    if isinstance(g, Cplx):
        return Cplx(negate(g.re), negate(g.im))
    return Symb("neg", g)


class Parser:
    def __init__(self, text, mode):
        self.t = tokenize(text)
        self.i = 0
        self.mode = mode

    def peek(self, k=0):
        j = self.i + k
        return self.t[j] if j < len(self.t) else (None, None)

    def take(self, val=None):
        tok = self.peek()
        if val is not None and tok[1] != val:
            raise SyntaxError(f"expected {val!r} got {tok!r} at {self.i}")
        self.i += 1
        return tok

    def at(self, *vals):
        k, v = self.peek()
        return k == 'op' and v in vals

    def parse(self):
        g = self.expr()
        if self.i != len(self.t):
            raise SyntaxError(f"trailing {self.t[self.i:]}")
        return g

    def expr(self):  # assignment
        lhs = self.equality()
        if self.at(":="):
            self.take()
            rhs = self.expr()
            return Symb("sto", Vect([rhs, lhs], 'seq'))
        return lhs

    def equality(self):
        g = self.comparison()
        while self.at("=", "==", "<>"):
            op = {"=": "equal", "==": "same", "<>": "different"}[self.take()[1]]
            g = Symb(op, Vect([g, self.comparison()], 'seq'))
        return g

    def comparison(self):
        g = self.interval()
        while self.at("<", ">", "<=", ">="):
            op = {"<": "inferieur_strict", "<=": "inferieur_egal", ">": "superieur_strict",
                  ">=": "superieur_egal"}[self.take()[1]]
            g = Symb(op, Vect([g, self.interval()], 'seq'))
        return g

    def interval(self):
        g = self.additive()
        if self.at(".."):
            self.take()
            g = Symb("interval", Vect([g, self.additive()], 'seq'))
        return g

    def additive(self):
        terms = [self.multiplicative()]
        while self.at("+", "-", ".-"):
            op = self.take()[1]
            rhs = self.multiplicative()
            if op == ".-":
                terms = [Symb("pointminus", Vect([self.mk_plus(terms), rhs], 'seq'))]
            else:
                terms.append(rhs if op == "+" else negate(rhs))
        return self.mk_plus(terms)

    def mk_plus(self, terms):
        if len(terms) == 1:
            return terms[0]
        if self.mode == 'eval' and len(terms) == 2 and is_num(terms[0]) and isinstance(terms[1], Cplx) \
                and numval(terms[1].re) == 0:
            return Cplx(terms[0], terms[1].im)
        return Symb("plus", Vect(terms, 'seq'))

    def multiplicative(self):
        facs = [self.division()]
        while self.at("*"):
            self.take()
            facs.append(self.division())
        if len(facs) == 1:
            return facs[0]
        if self.mode == 'eval' and len(facs) == 2 and is_num(facs[0]) and facs[1] == I:
            return Cplx(Int(0), facs[0])
        return Symb("prod", Vect(facs, 'seq'))

    def division(self):
        g = self.unary()
        while self.at("/"):
            self.take()
            d = self.unary()
            if self.mode == 'quote':
                g = Symb("division", Vect([g, d], 'seq'))
            elif isinstance(g, Int) and isinstance(d, Int):
                g = Frac(g, d)
            elif isinstance(g, Int) and g.v == 1:
                g = Symb("inv", d)
            else:
                g = Symb("prod", Vect([g, Symb("inv", d)], 'seq'))
        return g

    def unary(self):
        if self.at("-"):
            self.take()
            return negate(self.unary())
        if self.at("+"):
            self.take()
            return Symb("plus", self.unary())
        return self.power()

    def power(self):
        base = self.postfix()
        if self.at("^"):
            self.take()
            return Symb("pow", Vect([base, self.unary()], 'seq'))
        return base

    def postfix(self):
        g = self.primary()
        while True:
            if self.at("("):
                args = self.arglist("(", ")")
                arg = args[0] if len(args) == 1 else Vect(args, 'seq')
                if isinstance(g, Idnt) and g.name in BUILTINS | set(RENAME):
                    g = Symb(RENAME.get(g.name, g.name), arg)
                else:
                    g = Symb("of", Vect([g, arg], 'seq'))
            elif self.at("["):
                args = self.arglist("[", "]")
                g = Symb("at", Vect([g, args[0] if len(args) == 1 else Vect(args, 'seq')], 'seq'))
            elif self.at("!"):
                self.take()
                g = Symb("factorial", g)
            elif self.at("'"):
                self.take()
                g = Symb("derive", g)
            else:
                return g

    def arglist(self, lo, hi):
        self.take(lo)
        args = []
        if not self.at(hi):
            args.append(self.expr())
            while self.at(","):
                self.take()
                args.append(self.expr())
        self.take(hi)
        return args

    def primary(self):
        kind, val = self.take()
        if kind == 'num':
            if re.fullmatch(r"\d+", val):
                return Int(int(val))
            return Dbl(val)
        if kind == 'str':
            body = val[1:-1]
            return Str(re.sub(r'\\(.)', lambda m: {"n": "\n", "t": "\t"}.get(m.group(1), m.group(1)), body))
        if kind == 'id':
            if val in ("set", "poly1") and self.at("["):
                return Vect(self.arglist("[", "]"), {"set": "set", "poly1": "poly1"}[val])
            if val == "i":
                return I
            if val in ("true", "false"):
                return Bool(val == "true")
            return Idnt(val)
        if val == "'":  # quoted identifier 'Gamma'
            k, name = self.take()
            self.take("'")
            return Idnt(name)
        if val == "(":
            args = []
            args.append(self.expr())
            while self.at(","):
                self.take()
                args.append(self.expr())
            self.take(")")
            return args[0] if len(args) == 1 else Vect(args, 'seq')
        if val == "[":
            self.i -= 1
            return Vect(self.arglist("[", "]"), 'list')
        raise SyntaxError(f"unexpected {val!r}")


def parse(text, mode='eval'):
    return Parser(text, mode).parse()
