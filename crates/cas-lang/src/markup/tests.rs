//! Tests of the LaTeX and `MathML` printers. The engine on giac's own
//! answers: the 301 rows of `giac_latex.tsv` (1.6/T4.c1) and the 893 of
//! `giac_mathml.tsv` (1.6/T5.c1), each giac's printed answer read back into
//! giac's shape by a port of the prototype's `parse.py`, and printed as
//! giac's `latex()` and `export_mathml` did. Lowering: inputs read and
//! evaluated here, printed through the hook as giac prints its own answer,
//! except where our answer has another shape or order ("departure:"), and
//! in giac's `mathml` document (`giac_mathml_doc.tsv`, [`crate::mathml`]).
//! The limits.

use cas_core::{Context, Expr, Form, Integer, Symbol};

use super::engine::Engine;
use super::tree::{Gen, Sub};
use super::{MathmlKind, export_mathml, latex, mml};
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
                // poly1[…]: a list to gen2markup, poly1[ to mathml.cc.
                "set" | "poly1" if self.at(&["["]) => {
                    let sub = if val == "set" { Sub::Set } else { Sub::Poly1 };
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
                // A bare builtin evaluates to giac's _FUNC (:1647).
                Gen::Func(name)
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

pub(crate) fn parse(text: &str, quote: bool) -> Gen {
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
    let engine = Engine::new(usize::MAX);
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

/// giac's text of the annotation in a default-form row: its input form.
fn annotation(text: &str) -> &str {
    let open = "<annotation encoding='giac'>";
    let from = text.find(open).expect("an annotation") + open.len();
    &text[from..text.rfind("</annotation>").expect("closed")]
}

#[test]
fn the_engine_writes_giacs_mathml_of_giacs_answers() {
    let engine = Engine::new(usize::MAX);
    let mut failures = Vec::new();
    let mut rows = 0;
    for line in include_str!("giac_mathml.tsv").lines() {
        if line.starts_with('#') {
            continue;
        }
        let [input, mode, shape, kind, want] = line.split('\t').collect::<Vec<_>>()[..] else {
            panic!("a row of five columns: {line}");
        };
        rows += 1;
        // giac's doubles printed without a point (the prototype's check.py).
        let g = match input {
            "1e10" => Gen::Dbl {
                text: "10000000000".to_owned(),
                neg: false,
            },
            "12345678901.0" => Gen::Dbl {
                text: "12345678901".to_owned(),
                neg: false,
            },
            _ => parse(shape, mode == "quote"),
        };
        let got = match kind {
            "display" => engine.export_mathml_presentation(&g),
            "content" => engine.export_mathml_content(&g),
            // The annotation is the printed input, here giac's own.
            _ => engine
                .export_mathml(&g, "\0")
                .map(|t| t.replace('\0', annotation(want))),
        }
        .expect("no cap")
        // giac's table shows a newline as ' | ' (C0047).
        .replace('\n', " | ");
        if got != want {
            failures.push(format!("{kind} {input} ({shape}): {got:?}, giac {want:?}"));
        }
    }
    assert_eq!(rows, 893);
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
    assert_eq!(parse("'Gamma'", true), Gen::Func("Gamma".to_owned()));
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

/// The inputs of [`answers_print_as_giac_prints_its_own`] whose answers
/// have giac's shape: their `MathML` is giac's (`giac_mathml.tsv`).
const SAME_SHAPE: [&str; 58] = [
    "1/2",
    "-1/2",
    "x^2",
    "sqrt(x)",
    "sqrt(2)/2",
    "x/2",
    "2/x",
    "x*y",
    "2*x",
    "-x",
    "x-y",
    "-x-y",
    "1/x",
    "x^(-2)",
    "sin(x)",
    "cos(x)^2",
    "sin(x^2)",
    "exp(x)",
    "exp(2*x)",
    "ln(x)",
    "abs(x)",
    "floor(x)",
    "asin(x)",
    "pi",
    "infinity",
    "i",
    "1+2*i",
    "3.14",
    "1.5e-7",
    "-2.5",
    "3.0",
    "0.001",
    "1e-5",
    "1.23e20",
    "1e100",
    "12345678901.0",
    "1.0/3",
    "[1,2,3]",
    "[[1,2],[3,4]]",
    "x=1",
    "x^2+y^2=1",
    "f(x)",
    "f(x,y)",
    "diff(f(x),x)",
    "int(f(x),x)",
    "\"abc\"",
    "alpha",
    "abc",
    "x_1",
    "a_n",
    "n!",
    "x^y^z",
    "a/(b*c)",
    "2/3*x",
    "x*(-y)",
    "1..2",
    "true",
    "undef",
];

/// `text` read, evaluated, and printed in `form` through the hook, with
/// its input form.
fn printed(context: &Context, text: &str, form: Form) -> (String, String) {
    let hook = context.language().expect("installed");
    let e = context
        .eval(&hook.read(text, context).expect("one expression"))
        .expect("evaluates");
    let print = |form| {
        hook.print(&e, form, context, usize::MAX)
            .unwrap_or_else(|error| panic!("{text}: {error}"))
    };
    (print(form), print(Form::Input))
}

/// giac's `MathML` row of `input` in `kind`.
fn giac_mathml(input: &str, kind: &str) -> &'static str {
    include_str!("giac_mathml.tsv")
        .lines()
        .find_map(|line| match line.split('\t').collect::<Vec<_>>()[..] {
            [i, _, _, k, want] if i == input && k == kind => Some(want),
            _ => None,
        })
        .unwrap_or_else(|| panic!("no {kind} row of {input}"))
}

/// giac's `MathML` of the same inputs; in the default form the annotation is
/// our input form, escaped, and giac's is put in its place to compare.
#[test]
fn answers_in_mathml_as_giac_writes_its_own() {
    let context = context();
    let mut failures = Vec::new();
    for text in SAME_SHAPE {
        for (form, kind) in [
            (Form::MathmlDisplay, "display"),
            (Form::MathmlContent, "content"),
            (Form::Mathml, "default"),
        ] {
            let want = giac_mathml(text, kind);
            let (got, input) = printed(&context, text, form);
            let got = got.replace('\n', " | ");
            let mut ok = true;
            let got = if kind == "default" {
                let ours = annotation(&got).to_owned();
                ok = ours == mml::str_to_mml(&input).replace('\n', " | ");
                got.replace(
                    &format!("<annotation encoding='giac'>{ours}</annotation>"),
                    &format!(
                        "<annotation encoding='giac'>{}</annotation>",
                        annotation(want)
                    ),
                )
            } else {
                got
            };
            if !ok || got != want {
                failures.push(format!("{kind} {text} ({input}): {got:?}, giac {want:?}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// giac's `mathml` document of `input` (`giac_mathml_doc.tsv`).
fn giac_mathml_doc(input: &str) -> &'static str {
    include_str!("giac_mathml_doc.tsv")
        .lines()
        .find_map(|line| match line.split('\t').collect::<Vec<_>>()[..] {
            [i, _, _, want] if i == input => Some(want),
            _ => None,
        })
        .unwrap_or_else(|| panic!("no mathml row of {input}"))
}

/// giac's `mathml` (the document, `crate::mathml`) of the same inputs.
#[test]
fn answers_in_a_mathml_document_as_giac_writes_its_own() {
    let context = context();
    let mut failures = Vec::new();
    for text in SAME_SHAPE {
        let want = giac_mathml_doc(text);
        let (got, input) = printed(&context, text, Form::MathmlDocument);
        let got = got.replace('\n', " | ");
        if got != want {
            failures.push(format!("{text} ({input}): {got:?}, giac {want:?}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Departures: our answers' shapes or orders (as in
/// [`answers_print_as_giac_prints_its_own`] and cas-giac's `latex.rs`) in
/// presentation and content, the texts of the prototype
/// (`proto/mathml/proto/markup.py`) on our shapes; giac's are in
/// `giac_mathml.tsv`.
#[test]
#[allow(clippy::too_many_lines)]
fn mathml_departures() {
    let context = context();
    let rows: &[(&str, &str, &str)] = &[
        (
            "sqrt(x+1)",
            r"<msqrt><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></msqrt>",
            r"<apply id='id4'><root/><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        (
            "(x^4-1)/(x^2+3)",
            r"<mfrac><mrow><mo>&minus;</mo><mn>1</mn><mo>+</mo><msup><mi>x</mi><mn>4</mn></msup></mrow><mrow><mn>3</mn><mo>+</mo><msup><mi>x</mi><mn>2</mn></msup></mrow></mfrac>",
            r"<apply id='id11'><divide/><apply id='id5'><plus/><apply><minus/><cn id='id1' type='integer'>1</cn></apply><apply id='id4'><power/><ci id='id2'>x</ci><cn id='id3' type='integer'>4</cn></apply></apply><apply id='id10'><plus/><cn id='id6' type='integer'>3</cn><apply id='id9'><power/><ci id='id7'>x</ci><cn id='id8' type='integer'>2</cn></apply></apply></apply>",
        ),
        (
            "2*x^2-3*x+1",
            r"<mrow><mn>1</mn><mo>&minus;</mo><mrow><mn>3</mn><mo>&it;</mo><mi>x</mi></mrow><mo>+</mo><mrow><mn>2</mn><mo>&it;</mo><msup><mi>x</mi><mn>2</mn></msup></mrow></mrow>",
            r"<apply id='id10'><plus/><cn id='id1' type='integer'>1</cn><apply><minus/><apply id='id4'><times/><cn id='id2' type='integer'>3</cn><ci id='id3'>x</ci></apply></apply><apply id='id9'><times/><cn id='id5' type='integer'>2</cn><apply id='id8'><power/><ci id='id6'>x</ci><cn id='id7' type='integer'>2</cn></apply></apply></apply>",
        ),
        (
            "1/(x+1)",
            r"<mfrac><mn>1</mn><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfrac>",
            r"<apply id='id4'><divide/><cn type='integer'>1</cn><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        (
            "(x+1)^2",
            r"<msup><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced><mn>2</mn></msup>",
            r"<apply id='id5'><power/><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply><cn id='id4' type='integer'>2</cn></apply>",
        ),
        (
            "(x+1)*(x-1)",
            r"<mrow><mfenced><mrow><mo>&minus;</mo><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced><mo>&it;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced></mrow>",
            r"<apply id='id7'><times/><apply id='id3'><plus/><apply><minus/><cn id='id1' type='integer'>1</cn></apply><ci id='id2'>x</ci></apply><apply id='id6'><plus/><cn id='id4' type='integer'>1</cn><ci id='id5'>x</ci></apply></apply>",
        ),
        (
            "sin(x)*cos(x)",
            r"<mrow><mrow><mi>cos</mi><mo>&af;</mo><mi>x</mi></mrow><mo>&middot;</mo><mrow><mi>sin</mi><mo>&af;</mo><mi>x</mi></mrow></mrow>",
            r"<apply id='id7'><times/><apply id='id3'><cos id='id1'/><ci id='id2'>x</ci></apply><apply id='id6'><sin id='id4'/><ci id='id5'>x</ci></apply></apply>",
        ),
        (
            "-(x+1)",
            r"<mrow><mo>&minus;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced></mrow>",
            r"<apply id='id3'><minus/><apply><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        (
            "x<y",
            r"<mrow><mi>x</mi><mo>&lt;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><lt/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        (
            "x<=y",
            r"<mrow><mi>x</mi><mo>&leq;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><leq/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        (
            "x!=y",
            r"<mrow><mi>x</mi><mo>&#x2260;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><neq/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        (
            "ceil(x)",
            r"<mfenced open='&lceil;' close='&rceil;'><mi>x</mi></mfenced>",
            r"<apply id='id2'><ceiling/><ci id='id1'>x</ci></apply>",
        ),
        (
            "binomial(n,k)",
            r"<mfenced><mfrac linethickness='0'><mi>n</mi><mi>k</mi></mfrac></mfenced>",
            r"<apply id='id3'><csymbol cd='combinat1'>binomial</csymbol><ci id='id1'>n</ci><ci id='id2'>k</ci></apply>",
        ),
        (
            "x^(2/3)",
            r"<msup><mi>x</mi><mfrac><mn>2</mn><mn>3</mn></mfrac></msup>",
            r"<apply id='id3'><power/><ci id='id1'>x</ci><cn id='id2' type='rational'>2<sep/>3</cn></apply>",
        ),
        (
            "x*sqrt(2)",
            r"<mrow><msqrt><mn>2</mn></msqrt><mo>&it;</mo><mi>x</mi></mrow>",
            r"<apply id='id4'><times/><apply id='id2'><root/><cn id='id1' type='integer'>2</cn></apply><ci id='id3'>x</ci></apply>",
        ),
        (
            "sum(k,k,1,n)",
            r"<mrow><mfrac><mn>1</mn><mn>2</mn></mfrac><mo>&it;</mo><mi>n</mi><mo>&it;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>n</mi></mrow></mfenced></mrow>",
            r"<apply id='id6'><times/><cn id='id1' type='rational'>1<sep/>2</cn><ci id='id2'>n</ci><apply id='id5'><plus/><cn id='id3' type='integer'>1</cn><ci id='id4'>n</ci></apply></apply>",
        ),
        (
            "diff(f(x),x,2)",
            r"<mrow><mi>derive</mi><mfenced><mrow><mfrac><mi mathvariant='normal'>d</mi><mrow><mi>x</mi></mrow></mfrac><mo>&af;</mo><mrow><mi>f</mi><mo>&af;</mo><mo>&af;</mo><mfenced><mi>x</mi></mfenced></mrow></mrow><mn>2</mn></mfenced></mrow>",
            r"<apply id='id13'><ci id='id6' type='function'>derive</ci><apply id='id11'><diff/><bvar><ci id='id10'>x</ci></bvar><apply id='id9'><ci id='id7'>f</ci><ci id='id8'>x</ci></apply></apply><cn id='id12' type='integer'>2</cn></apply>",
        ),
    ];
    let math = |body: &str| {
        format!("<math mode='display' xmlns='http://www.w3.org/1998/Math/MathML'>{body}</math>")
    };
    let mut failures = Vec::new();
    for (text, display, content) in rows {
        for (form, want) in [
            (Form::MathmlDisplay, display),
            (Form::MathmlContent, content),
        ] {
            let (got, _) = printed(&context, text, form);
            if got != math(want) {
                failures.push(format!(
                    "{form:?} {text}: {got:?}, expected {:?}",
                    math(want)
                ));
            }
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
    let error = latex(&list(e.clone()), names, balls, usize::MAX).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooDeep);
    for kind in [
        MathmlKind::Default,
        MathmlKind::Display,
        MathmlKind::Content,
    ] {
        assert!(export_mathml(&e, names, balls, usize::MAX, kind).is_ok());
        let error = export_mathml(&list(e.clone()), names, balls, usize::MAX, kind).unwrap_err();
        assert_eq!(error.kind, PrintErrorKind::TooDeep);
    }
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
    for kind in [
        MathmlKind::Default,
        MathmlKind::Display,
        MathmlKind::Content,
    ] {
        let text = export_mathml(&e, names, balls, usize::MAX, kind).unwrap();
        assert_eq!(
            export_mathml(&e, names, balls, text.len(), kind).unwrap(),
            text
        );
        let error = export_mathml(&e, names, balls, text.len() - 1, kind).unwrap_err();
        assert_eq!(
            (error.kind, error.part),
            (PrintErrorKind::TooLong, e.clone())
        );
    }
    // A long atom deep inside is refused as its block is made.
    let long = Expr::string(&"a".repeat(1000));
    let e = Expr::apply(&Symbol::PLUS.into(), &[long, Expr::from(Symbol::new("x"))]);
    let error = latex(&e, names, balls, 100).unwrap_err();
    assert_eq!(error.kind, PrintErrorKind::TooLong);
    // Through the hook, as text.
    let context = context();
    let hook = context.language().unwrap();
    assert!(hook.print(&e, Form::Latex, &context, 100).is_err());
    for kind in [
        MathmlKind::Default,
        MathmlKind::Display,
        MathmlKind::Content,
    ] {
        let error = export_mathml(&e, names, balls, 100, kind).unwrap_err();
        assert_eq!(error.kind, PrintErrorKind::TooLong);
    }
    assert!(hook.print(&e, Form::Mathml, &context, 100).is_err());
}

/// `increasing_power` off (CAS-618): giac's `latex` and `export_mathml` of
/// the same answers with the flag off (1.6/T6.c3 probes); with it on, ours
/// as before.
#[test]
fn sums_of_one_symbol_by_decreasing_degree() {
    let mut context = context();
    let cubic = "x^3+3*x^2+3*x+1";
    assert_eq!(latex_of(&context, cubic), "1+3 x+3 x^{2}+x^{3}");
    context.set_increasing_power(false);
    for (text, want) in [
        (cubic, "x^{3}+3 x^{2}+3 x+1"),
        ("1-x^2", "-x^{2}+1"),
        ("2*x^2-x+5", "2 x^{2}-x+5"),
        ("(x+1)^2", "\\left(x+1\\right)^{2}"),
        ("sin(x^2+2*x+1)", "\\sin \\left(x^{2}+2 x+1\\right)"),
    ] {
        assert_eq!(latex_of(&context, text), want, "{text}");
    }
    // Two symbols: the order stays.
    let two = latex_of(&context, "x^2+y+1");
    context.set_increasing_power(true);
    assert_eq!(latex_of(&context, "x^2+y+1"), two);
    context.set_increasing_power(false);
    let (display, _) = printed(&context, cubic, Form::MathmlDisplay);
    assert_eq!(
        display,
        "<math mode='display' xmlns='http://www.w3.org/1998/Math/MathML'><mrow><msup><mi>x</mi><mn>3</mn></msup><mo>+</mo><mrow><mn>3</mn><mo>&it;</mo><msup><mi>x</mi><mn>2</mn></msup></mrow><mo>+</mo><mrow><mn>3</mn><mo>&it;</mo><mi>x</mi></mrow><mo>+</mo><mn>1</mn></mrow></math>"
    );
    let (semantic, _) = printed(&context, cubic, Form::MathmlContent);
    assert_eq!(
        semantic,
        "<math mode='display' xmlns='http://www.w3.org/1998/Math/MathML'><apply id='id13'><plus/><apply id='id3'><power/><ci id='id1'>x</ci><cn id='id2' type='integer'>3</cn></apply><apply id='id8'><times/><cn id='id4' type='integer'>3</cn><apply id='id7'><power/><ci id='id5'>x</ci><cn id='id6' type='integer'>2</cn></apply></apply><apply id='id11'><times/><cn id='id9' type='integer'>3</cn><ci id='id10'>x</ci></apply><cn id='id12' type='integer'>1</cn></apply></math>"
    );
    let (default, input) = printed(&context, cubic, Form::Mathml);
    assert_eq!(input, "x^3 + 3*x^2 + 3*x + 1");
    assert!(default.contains("<annotation-xml encoding='MathML-Content'><apply id='id13'><plus/><apply id='id3'><power/>"), "{default}");
}
