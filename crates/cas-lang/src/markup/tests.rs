//! Tests of the LaTeX printer. The engine on giac's own answers: the 301
//! rows of `giac_latex.tsv` (1.6/T4.c1), each giac's printed answer read
//! back into giac's shape by a port of the prototype's `parse.py`, and
//! printed as giac's `latex()` did. Lowering: inputs read and evaluated
//! here, printed through the hook as giac prints its own answer, except
//! where our answer has another shape or order ("departure:"). The limits.

use cas_core::{Context, Expr, Form, Integer, Symbol};

use super::engine::Engine;
use super::latex;
use super::tree::{Gen, Sub};
use crate::hook::Hook;
use crate::names::Names;
use crate::print::{Balls, MAX_PRINT_DEPTH, PrintErrorKind};

// ---- giac's printed answers as trees (the prototype's parse.py) --------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Tok {
    Num,
    Str,
    Id,
    Op,
}

fn tokenize(s: &str) -> Vec<(Tok, String)> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    let digits = |i: &mut usize| {
        let start = *i;
        while b.get(*i).is_some_and(u8::is_ascii_digit) {
            *i += 1;
        }
        *i > start
    };
    let exponent = |i: &mut usize| {
        if matches!(b.get(*i), Some(b'e' | b'E')) {
            let mut j = *i + 1;
            if matches!(b.get(j), Some(b'+' | b'-')) {
                j += 1;
            }
            let start = j;
            while b.get(j).is_some_and(u8::is_ascii_digit) {
                j += 1;
            }
            if j > start {
                *i = j;
            }
        }
    };
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        let start = i;
        let kind =
            if c.is_ascii_digit() || (c == b'.' && b.get(i + 1).is_some_and(u8::is_ascii_digit)) {
                digits(&mut i);
                if b.get(i) == Some(&b'.') {
                    i += 1;
                    digits(&mut i);
                }
                exponent(&mut i);
                Tok::Num
            } else if c == b'"' {
                i += 1;
                while b[i] != b'"' {
                    i += if b[i] == b'\\' { 2 } else { 1 };
                }
                i += 1;
                Tok::Str
            } else if c.is_ascii_alphabetic() || c == b'_' {
                while b
                    .get(i)
                    .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
                {
                    i += 1;
                }
                Tok::Id
            } else {
                let two = s.get(i..i + 2).unwrap_or("");
                i += if [":=", "==", "<>", "<=", ">=", "..", ".-"].contains(&two) {
                    2
                } else {
                    assert!("-+*/^()[],=<>!'".contains(char::from(c)), "bad char in {s}");
                    1
                };
                Tok::Op
            };
        out.push((kind, s[start..i].to_owned()));
    }
    out
}

fn is_num(g: &Gen) -> bool {
    matches!(g, Gen::Int(_) | Gen::Dbl { .. } | Gen::Frac(..))
}

/// Unary minus with giac's numeric folding.
fn negate(g: Gen) -> Gen {
    match g {
        Gen::Int(n) => Gen::Int(-n),
        Gen::Dbl { text, neg } => Gen::Dbl { text, neg: !neg },
        Gen::Frac(n, d) if matches!(*n, Gen::Int(_)) => Gen::Frac(Box::new(negate(*n)), d),
        Gen::Cplx(re, im) => Gen::Cplx(Box::new(negate(*re)), Box::new(negate(*im))),
        g => Gen::symb("neg", g),
    }
}

/// giac's printed answer as giac's tree: `eval` shapes (`n/m` a fraction,
/// `1/b` `inv(b)`, `a/b` `prod(a, inv(b))`, numeric complexes folded) or,
/// for `quote`, `a/b` as `division(a, b)`.
struct Parser {
    t: Vec<(Tok, String)>,
    i: usize,
    quote: bool,
}

impl Parser {
    fn peek(&self) -> Option<&(Tok, String)> {
        self.t.get(self.i)
    }

    fn take(&mut self) -> (Tok, String) {
        self.i += 1;
        self.t[self.i - 1].clone()
    }

    fn expect(&mut self, v: &str) {
        assert_eq!(self.take().1, v);
    }

    fn at(&self, vals: &[&str]) -> bool {
        matches!(self.peek(), Some((Tok::Op, v)) if vals.contains(&v.as_str()))
    }

    fn expr(&mut self) -> Gen {
        let lhs = self.equality();
        if self.at(&[":="]) {
            self.take();
            let rhs = self.expr();
            return Gen::symb_n("sto", vec![rhs, lhs]);
        }
        lhs
    }

    fn equality(&mut self) -> Gen {
        let mut g = self.comparison();
        while self.at(&["=", "==", "<>"]) {
            let op = match self.take().1.as_str() {
                "=" => "equal",
                "==" => "same",
                _ => "different",
            };
            g = Gen::symb_n(op, vec![g, self.comparison()]);
        }
        g
    }

    fn comparison(&mut self) -> Gen {
        let mut g = self.interval();
        while self.at(&["<", ">", "<=", ">="]) {
            let op = match self.take().1.as_str() {
                "<" => "inferieur_strict",
                "<=" => "inferieur_egal",
                ">" => "superieur_strict",
                _ => "superieur_egal",
            };
            g = Gen::symb_n(op, vec![g, self.interval()]);
        }
        g
    }

    fn interval(&mut self) -> Gen {
        let g = self.additive();
        if self.at(&[".."]) {
            self.take();
            return Gen::symb_n("interval", vec![g, self.additive()]);
        }
        g
    }

    fn additive(&mut self) -> Gen {
        let mut terms = vec![self.multiplicative()];
        while self.at(&["+", "-", ".-"]) {
            let op = self.take().1;
            let rhs = self.multiplicative();
            match op.as_str() {
                ".-" => {
                    let lhs = self.mk_plus(std::mem::take(&mut terms));
                    terms = vec![Gen::symb_n("pointminus", vec![lhs, rhs])];
                }
                "+" => terms.push(rhs),
                _ => terms.push(negate(rhs)),
            }
        }
        self.mk_plus(terms)
    }

    fn mk_plus(&self, mut terms: Vec<Gen>) -> Gen {
        if terms.len() == 1 {
            return terms.remove(0);
        }
        if !self.quote
            && let [a, Gen::Cplx(re, im)] = &terms[..]
            && is_num(a)
            && re.is_zero()
        {
            return Gen::Cplx(Box::new(a.clone()), im.clone());
        }
        Gen::symb_n("plus", terms)
    }

    fn multiplicative(&mut self) -> Gen {
        let mut facs = vec![self.division()];
        while self.at(&["*"]) {
            self.take();
            facs.push(self.division());
        }
        if facs.len() == 1 {
            return facs.remove(0);
        }
        if !self.quote && facs.len() == 2 && is_num(&facs[0]) && facs[1] == Gen::i() {
            return Gen::Cplx(Box::new(Gen::int(0)), Box::new(facs.remove(0)));
        }
        Gen::symb_n("prod", facs)
    }

    fn division(&mut self) -> Gen {
        let mut g = self.unary();
        while self.at(&["/"]) {
            self.take();
            let d = self.unary();
            g = if self.quote {
                Gen::symb_n("division", vec![g, d])
            } else if matches!((&g, &d), (Gen::Int(_), Gen::Int(_))) {
                Gen::Frac(Box::new(g), Box::new(d))
            } else if g.is_one() && matches!(g, Gen::Int(_)) {
                Gen::symb("inv", d)
            } else {
                Gen::symb_n("prod", vec![g, Gen::symb("inv", d)])
            };
        }
        g
    }

    fn unary(&mut self) -> Gen {
        if self.at(&["-"]) {
            self.take();
            return negate(self.unary());
        }
        if self.at(&["+"]) {
            self.take();
            return Gen::symb("plus", self.unary());
        }
        self.power()
    }

    fn power(&mut self) -> Gen {
        let base = self.postfix();
        if self.at(&["^"]) {
            self.take();
            return Gen::symb_n("pow", vec![base, self.unary()]);
        }
        base
    }

    fn postfix(&mut self) -> Gen {
        const BUILTINS: &[&str] = &[
            "sin",
            "cos",
            "tan",
            "cot",
            "sec",
            "csc",
            "asin",
            "acos",
            "atan",
            "acot",
            "asec",
            "acsc",
            "sinh",
            "cosh",
            "tanh",
            "asinh",
            "acosh",
            "atanh",
            "ln",
            "exp",
            "sqrt",
            "abs",
            "floor",
            "ceil",
            "erf",
            "Gamma",
            "Psi",
            "Dirac",
            "Heaviside",
            "Zeta",
            "Beta",
            "rootof",
            "integrate",
            "int",
            "sum",
            "product",
            "limit",
            "diff",
            "derive",
            "factorial",
            "binomial",
            "sto",
            "inv",
            "log",
        ];
        let mut g = self.primary();
        loop {
            if self.at(&["("]) {
                let mut args = self.arglist("(", ")");
                let arg = if args.len() == 1 {
                    args.remove(0)
                } else {
                    Gen::Vect(args, Sub::Seq)
                };
                g = match g {
                    Gen::Idnt(name) if BUILTINS.contains(&name.as_str()) => {
                        let op = match name.as_str() {
                            "diff" => "derive",
                            "log" => "ln",
                            n => n,
                        };
                        Gen::symb(op, arg)
                    }
                    g => Gen::symb_n("of", vec![g, arg]),
                };
            } else if self.at(&["["]) {
                let mut args = self.arglist("[", "]");
                let arg = if args.len() == 1 {
                    args.remove(0)
                } else {
                    Gen::Vect(args, Sub::Seq)
                };
                g = Gen::symb_n("at", vec![g, arg]);
            } else if self.at(&["!"]) {
                self.take();
                g = Gen::symb("factorial", g);
            } else if self.at(&["'"]) {
                self.take();
                g = Gen::symb("derive", g);
            } else {
                return g;
            }
        }
    }

    fn arglist(&mut self, lo: &str, hi: &str) -> Vec<Gen> {
        self.expect(lo);
        let mut args = Vec::new();
        if !self.at(&[hi]) {
            args.push(self.expr());
            while self.at(&[","]) {
                self.take();
                args.push(self.expr());
            }
        }
        self.expect(hi);
        args
    }

    fn primary(&mut self) -> Gen {
        let (kind, val) = self.take();
        match kind {
            Tok::Num if val.bytes().all(|c| c.is_ascii_digit()) => {
                Gen::Int(val.parse::<Integer>().expect("digits"))
            }
            Tok::Num => Gen::Dbl {
                text: val,
                neg: false,
            },
            Tok::Str => {
                let body = &val[1..val.len() - 1];
                let mut s = String::new();
                let mut chars = body.chars();
                while let Some(c) = chars.next() {
                    if c == '\\' {
                        match chars.next() {
                            Some('n') => s.push('\n'),
                            Some('t') => s.push('\t'),
                            Some(c) => s.push(c),
                            None => {}
                        }
                    } else {
                        s.push(c);
                    }
                }
                Gen::Str(s)
            }
            Tok::Id => match val.as_str() {
                // poly1[…] prints as a list.
                "set" | "poly1" if self.at(&["["]) => {
                    let sub = if val == "set" { Sub::Set } else { Sub::List };
                    Gen::Vect(self.arglist("[", "]"), sub)
                }
                "i" => Gen::i(),
                "true" => Gen::Bool(true),
                "false" => Gen::Bool(false),
                _ => Gen::Idnt(val),
            },
            Tok::Op if val == "'" => {
                let (_, name) = self.take();
                self.expect("'");
                Gen::Idnt(name)
            }
            Tok::Op if val == "(" => {
                let mut args = vec![self.expr()];
                while self.at(&[","]) {
                    self.take();
                    args.push(self.expr());
                }
                self.expect(")");
                if args.len() == 1 {
                    args.remove(0)
                } else {
                    Gen::Vect(args, Sub::Seq)
                }
            }
            Tok::Op if val == "[" => {
                self.i -= 1;
                Gen::Vect(self.arglist("[", "]"), Sub::List)
            }
            Tok::Op => panic!("unexpected {val}"),
        }
    }
}

fn parse(text: &str, quote: bool) -> Gen {
    let mut p = Parser {
        t: tokenize(text),
        i: 0,
        quote,
    };
    let g = p.expr();
    assert_eq!(p.i, p.t.len(), "trailing tokens in {text}");
    g
}

// ---- the engine on giac's answers --------------------------------------

#[test]
fn the_engine_writes_giacs_latex_of_giacs_answers() {
    let engine = Engine {
        max_len: usize::MAX,
    };
    let mut failures = Vec::new();
    let mut rows = 0;
    for line in include_str!("giac_latex.tsv").lines() {
        if line.starts_with('#') {
            continue;
        }
        let [input, mode, shape, want] = line.split('\t').collect::<Vec<_>>()[..] else {
            panic!("a row of four columns: {line}");
        };
        rows += 1;
        let got = engine
            .export_latex(&parse(shape, mode == "quote"))
            .expect("no cap");
        if got != want {
            failures.push(format!("{input} ({shape}): {got:?}, giac {want:?}"));
        }
    }
    assert_eq!(rows, 301);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn the_test_parser_keeps_giacs_shapes() {
    let x = || Gen::idnt("x");
    let inv = |g| Gen::symb("inv", g);
    assert_eq!(
        parse("x/2", false),
        Gen::symb_n("prod", vec![x(), inv(Gen::int(2))])
    );
    assert_eq!(parse("1/x", false), inv(x()));
    assert_eq!(
        parse("2/3", false),
        Gen::Frac(Box::new(Gen::int(2)), Box::new(Gen::int(3)))
    );
    assert_eq!(
        parse("x/3", true),
        Gen::symb_n("division", vec![x(), Gen::int(3)])
    );
    assert_eq!(
        parse("1+2*i", false),
        Gen::Cplx(Box::new(Gen::int(1)), Box::new(Gen::int(2)))
    );
    assert_eq!(
        parse("-x*y", false),
        Gen::symb_n("prod", vec![Gen::symb("neg", x()), Gen::idnt("y")])
    );
}

// ---- lowering: our answers in giac's shapes -----------------------------

fn context() -> Context {
    let mut context = Context::new();
    context.set_language(Some(Hook::english()));
    context
}

/// `text` read, evaluated, and printed in LaTeX through the hook.
fn latex_of(context: &Context, text: &str) -> String {
    let hook = context.language().expect("installed");
    let e = context
        .eval(&hook.read(text, context).expect("one expression"))
        .expect("evaluates");
    hook.print(&e, Form::Latex, context, usize::MAX)
        .unwrap_or_else(|error| panic!("{text}: {error}"))
}

#[test]
fn answers_print_as_giac_prints_its_own() {
    let context = context();
    // Expected texts are giac's latex() of the same input (table.tsv),
    // except "departure:" rows, where our answer differs in shape or order.
    let rows = [
        ("1/2", "\\frac{1}{2}"),
        ("-1/2", "-\\frac{1}{2}"),
        ("x^2", "x^{2}"),
        ("sqrt(x)", "\\sqrt{x}"),
        // departure: our sums put constants first (giac: \sqrt{x+1}).
        ("sqrt(x+1)", "\\sqrt{1+x}"),
        ("sqrt(2)/2", "\\frac{\\sqrt{2}}{2}"),
        // departure: constants first, ascending degree (giac:
        // \frac{x^{4}-1}{x^{2}+3}).
        ("(x^4-1)/(x^2+3)", "\\frac{-1+x^{4}}{3+x^{2}}"),
        ("x/2", "\\frac{x}{2}"),
        ("2/x", "\\frac{2}{x}"),
        ("x*y", "x y"),
        ("2*x", "2 x"),
        ("-x", "-x"),
        ("x-y", "x-y"),
        ("-x-y", "-x-y"),
        // departure: ascending degree (giac: 2 x^{2}-3 x+1).
        ("2*x^2-3*x+1", "1-3 x+2 x^{2}"),
        ("1/x", "\\frac1{x}"),
        ("x^(-2)", "\\frac1{x^{2}}"),
        // departure: constants first (giac: \frac1{x+1}).
        ("1/(x+1)", "\\frac1{1+x}"),
        // departure: constants first (giac: \left(x+1\right)^{2}).
        ("(x+1)^2", "\\left(1+x\\right)^{2}"),
        // departure: constants first, factors in our order (giac:
        // \left(x+1\right) \left(x-1\right)).
        ("(x+1)*(x-1)", "\\left(-1+x\\right) \\left(1+x\\right)"),
        ("sin(x)", "\\sin x"),
        ("cos(x)^2", "\\cos ^{2}x"),
        ("sin(x^2)", "\\sin \\left(x^{2}\\right)"),
        // departure: our products order cos before sin (giac:
        // \sin x\cdot \cos x).
        ("sin(x)*cos(x)", "\\cos x\\cdot \\sin x"),
        ("exp(x)", "\\mathrm{e}^{x}"),
        ("exp(2*x)", "\\mathrm{e}^{2 x}"),
        ("ln(x)", "\\ln x"),
        ("abs(x)", "\\left|x\\right|"),
        ("floor(x)", "\\left\\lfloor x\\right\\rfloor "),
        ("asin(x)", "\\arcsin x"),
        ("pi", "\\pi "),
        ("infinity", "\\infty "),
        ("i", "\\mathrm{i}"),
        ("1+2*i", "1+2 \\mathrm{i}"),
        ("3.14", "3.14"),
        ("1.5e-7", "1.5\\times10^{-7}"),
        ("-2.5", "-2.5"),
        ("3.0", "3.0"),
        ("0.001", "0.001"),
        ("1e-5", "1.0\\times10^{-5}"),
        ("1.23e20", "1.23\\times10^{20}"),
        ("1e100", "1.0\\times10^{100}"),
        ("12345678901.0", "12345678901"),
        ("1.0/3", "0.333333333333"),
        ("[1,2,3]", "\\left[1,2,3\\right]"),
        (
            "[[1,2],[3,4]]",
            "\\left[\\begin{array}{cc}1&2\\\\3&4\\end{array}\\right]",
        ),
        ("x=1", "x=1"),
        ("x^2+y^2=1", "x^{2}+y^{2}=1"),
        // departure: giac's evaluator turns x<y into y>x.
        ("x<y", "x<y"),
        ("f(x)", "f\\left(x\\right)"),
        ("f(x,y)", "f\\left(x,y\\right)"),
        (
            "diff(f(x),x)",
            "\\frac{\\mathrm{d}}{\\mathrm{d}x}f\\left(x\\right)",
        ),
        ("int(f(x),x)", "\\int f\\left(x\\right)\\,\\mathrm{d}x"),
        ("\"abc\"", "\\text{``abc''}"),
        ("alpha", "\\alpha"),
        ("abc", "\\mathrm{abc}"),
        ("x_1", "x_{1}"),
        ("a_n", "a_{n}"),
        ("n!", "n\\mathclose{!}"),
        ("x^y^z", "x^{y^{z}}"),
        ("a/(b*c)", "\\frac{a}{b c}"),
        ("2/3*x", "\\frac{2}{3} x"),
        // departure: ours keeps Times[-1, Plus[1, x]] (giac: -x-1).
        ("-(x+1)", "-\\left(1+x\\right)"),
        ("x*(-y)", "-x y"),
        ("1..2", "1\\ldots 2"),
        ("true", "\\mathrm{true}"),
        ("undef", "\\mathrm{undef}"),
    ];
    let mut failures = Vec::new();
    for (text, want) in rows {
        let got = latex_of(&context, text);
        if got != want {
            failures.push(format!("{text}: {got:?}, expected {want:?}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// ---- limits ---------------------------------------------------------------

#[test]
fn too_deep_is_refused_before_any_work() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let list = |e: Expr| Expr::apply(&Symbol::LIST.into(), &[e]);
    let mut e = Expr::from(1);
    // An atom is 1 deep.
    for _ in 1..MAX_PRINT_DEPTH {
        e = list(e);
    }
    assert!(latex(&e, names, balls, usize::MAX).is_ok());
    let error = latex(&list(e), names, balls, usize::MAX).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooDeep);
}

#[test]
fn text_past_the_cap_is_refused() {
    let names = Names::english();
    let balls = Balls::Shown {
        working_precision: 53,
    };
    let e = Expr::apply(
        &Symbol::LIST.into(),
        &(0..100).map(Expr::from).collect::<Vec<_>>(),
    );
    let text = latex(&e, names, balls, usize::MAX).unwrap();
    assert_eq!(latex(&e, names, balls, text.len()).unwrap(), text);
    let error = latex(&e, names, balls, text.len() - 1).unwrap_err();
    assert_eq!(
        (error.kind, error.part),
        (PrintErrorKind::TooLong, e.clone())
    );
    // A long atom deep inside is refused as its block is made.
    let long = Expr::string(&"a".repeat(1000));
    let e = Expr::apply(&Symbol::PLUS.into(), &[long, Expr::from(Symbol::new("x"))]);
    let error = latex(&e, names, balls, 100).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooLong);
    // Through the hook, as text.
    let context = context();
    let hook = context.language().unwrap();
    assert!(hook.print(&e, Form::Latex, &context, 100).is_err());
}
