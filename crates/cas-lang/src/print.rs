//! The input form: expressions as text of the language (CAS-43), see
//! [`print()`].

use std::fmt::{self, Write as _};

use cas_core::{Complex, Expr, Integer, Rational, Symbol, View};

use crate::lexer::lex;
use crate::literal::name;
use crate::lower::is_value;
use crate::names::{Names, SyntaxForm, Target, escape_name, unescape_name};
use crate::parser::MAX_DEPTH;
use crate::syntax::SyntaxKind;

/// The deepest expression [`print()`] writes, by [`Expr::depth`].
///
/// The text of a deeper one could nest beyond [`MAX_DEPTH`], which the
/// parser refuses: each level of an expression takes at most two levels of
/// the parser (`a + (…)`, `if … then … end`), an atom at most two more
/// (`x*(2 + 3*i/4)`), and the statement one.
pub const MAX_PRINT_DEPTH: usize = (MAX_DEPTH - 1) / 2;

/// Why [`print()`] cannot write an expression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrintError {
    /// The part the language cannot write; the whole expression when it is
    /// too deep.
    pub part: Expr,
    /// Why.
    pub kind: PrintErrorKind,
}

/// What a [`PrintError`] is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PrintErrorKind {
    /// A symbol without a name in the language: an internal builtin on its
    /// own (`Plus`, `Null`), or a symbol whose name would read back as
    /// another one: a reserved name (`sin` as a user's symbol), a name that
    /// is not an identifier, or a capitalised one without its `$` (§5).
    NoName,
    /// An application without syntax: an operator with another number of
    /// arguments (`Plus[x]`), a pattern, a rule or a statement.
    NoSyntax,
    /// A parameter of a lambda that is not a user's name.
    Parameter,
    /// Nested deeper than [`MAX_PRINT_DEPTH`].
    TooDeep,
}

impl fmt::Display for PrintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind {
            PrintErrorKind::NoName => write!(f, "the language has no name for {:?}", self.part),
            PrintErrorKind::NoSyntax => match self.part.head() {
                Some(head) => {
                    let n = self.part.args().len();
                    let s = if n == 1 { "" } else { "s" };
                    write!(
                        f,
                        "the language has no syntax for {head:?} with {n} argument{s}"
                    )
                }
                None => write!(f, "the language has no syntax for {:?}", self.part),
            },
            PrintErrorKind::Parameter => {
                write!(f, "{:?} cannot be the parameter of a lambda", self.part)
            }
            PrintErrorKind::TooDeep => write!(f, "nested more than {MAX_PRINT_DEPTH} deep"),
        }
    }
}

impl std::error::Error for PrintError {}

/// Writes `e` as text that [`parse`](crate::parse) and
/// [`lower`](crate::lower) read back as `e`, with the names of `names`: the
/// input form, so output can be pasted as input.
///
/// Lowering builds what the text says and does no arithmetic. So the text
/// gives back the expression itself where the language writes its structure
/// (`f(x)`, `x^2`, `[a, b]`, `0 < x <= 1`), and otherwise the arithmetic
/// that evaluates to it:
///
/// - a rational is a quotient, `3/4`, and a complex number a sum or a
///   product, `2 + 3*i`, `-i/2`;
/// - a product has its sign and coefficient in front and its factors with
///   negative integer exponents in a denominator: `Times[-3/4, x,
///   Power[y, -2]]` is `-3*x/(4*y^2)`, and `Power[x, -1]` is `1/x`;
/// - a term with a negative coefficient is subtracted: `x - 2*y`.
///
/// Evaluation makes the same expression of both (CAS-DOC-4): products are
/// flat and orderless, and numbers combine.
///
/// The other conventions:
///
/// - Parentheses only where the precedence table (CAS-DOC-7 §3) needs them,
///   with two more for reading: an operand of an arithmetic operator never
///   starts with a sign (`2*(-x)`, `x^(-1/2)`), and `if … end` is
///   parenthesised as an operand.
/// - `*` in every product. Implicit products bind tighter than `*` and `/`
///   (`1/2x` is `1/(2x)`), so the printer leaves them to people.
/// - Names as the table spells them in its locale, user names without the
///   `$` of capitalised ones (§5). A function's name as the base of a power
///   is parenthesised, `(sin)^2`, since `sin^2` is an error (§3.2).
/// - `call(f, …)` for an application that the call-or-product rule would
///   read as a product (a head that occurs in its single argument, or a
///   constant head) and for one whose head is not a name, an application, an
///   index or a lambda (§3.1).
/// - `=` for `Equation` and `==` for `Equal`; an `And` of comparisons that
///   share their operands is a chain, `0 < x <= 1`.
/// - `+infinity` for `Infinity` and `-infinity` for `DirectedInfinity[-1]`,
///   whole forms (§5); so the negated complex infinity is `-(infinity)`.
///
/// What the language cannot write is an error ([`PrintError`]): internal
/// builtins on their own (`Plus`), operators with another number of
/// arguments (`Plus[x]`), patterns, rules and statements, symbols whose names
/// would read back as other symbols, and expressions nested deeper than
/// [`MAX_PRINT_DEPTH`].
///
/// ```
/// use cas_core::{Expr, Rational, Symbol};
/// use cas_lang::{Names, print};
///
/// let names = Names::english();
/// let [x, y] = ["x", "y"].map(|n| Expr::from(Symbol::new(n)));
/// let apply = |head: Symbol, args: &[Expr]| Expr::apply(&head.into(), args);
///
/// let coefficient = Expr::from(Rational::new((-3).into(), 4.into()));
/// let e = apply(Symbol::TIMES, &[coefficient, x.clone(), y.pow(-2)]);
/// assert_eq!(print(&e, names).unwrap(), "-3*x/(4*y^2)");
///
/// // `x(x + 1)` would read as a product.
/// let e = Expr::apply(&x, &[apply(Symbol::PLUS, &[x.clone(), Expr::from(1)])]);
/// assert_eq!(print(&e, names).unwrap(), "call(x, x + 1)");
/// ```
///
/// # Errors
///
/// A [`PrintError`] for a part the language cannot write.
pub fn print(e: &Expr, names: &Names) -> Result<String, PrintError> {
    if usize::from(e.depth()) > MAX_PRINT_DEPTH {
        return Err(PrintError {
            part: e.clone(),
            kind: PrintErrorKind::TooDeep,
        });
    }
    let mut printer = Printer {
        names,
        out: String::new(),
    };
    printer.expr(e)?;
    Ok(printer.out)
}

/// How tightly a printed form binds, weakest first: the levels of
/// CAS-DOC-7 §3 that the printer writes.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Prec {
    /// Lambdas, and `if … end` as an operand.
    Lowest,
    Or,
    And,
    Not,
    /// `=` and the comparisons.
    Compare,
    Range,
    Sum,
    /// Products and quotients, rationals among them.
    Product,
    /// A sign in front: `-x`, `-3`, `+infinity`.
    Prefix,
    Power,
    /// Names, numbers, strings, brackets, calls, indexes, factorials.
    Postfix,
}

/// How an expression prints: how tightly it binds, and whether it starts
/// with a sign.
#[derive(Clone, Copy, Debug)]
struct Shape {
    prec: Prec,
    sign: bool,
}

/// How an application prints.
enum Form<'a> {
    Sum(&'a [Expr]),
    /// `Times[…]`, or a lone `Power` with a negative integer exponent.
    Product(&'a [Expr]),
    /// `b^n` or `b.^n`.
    Power(&'a Expr, &'static str, &'a Expr),
    /// `a.*b.*c` or `a./b`.
    Elementwise(&'a [Expr], &'static str),
    Factorial(&'a Expr),
    /// `a = b` and the comparisons.
    Relation(&'a Expr, &'static str, &'a Expr),
    /// `And` of comparisons that share their operands: `a < b <= c`.
    Chain(&'a [Expr]),
    /// `and` or `or`, binding at the first level, with operands at the
    /// second.
    Logical(&'a [Expr], &'static str, Prec, Prec),
    Not(&'a Expr),
    Range(&'a Expr, &'a Expr),
    /// A list `[…]` or a set `{…}`.
    Bracketed(&'a [Expr], char, char),
    Index(&'a Expr, &'a [Expr]),
    /// `x -> body`, with the parameters.
    Lambda(&'a [Expr], &'a Expr),
    If(&'a [Expr]),
    MinusInfinity,
    /// Anything else: `head(args…)` or `call(head, args…)`.
    Call(&'a Expr, &'a [Expr]),
}

/// How the application `e` prints.
fn form(e: &Expr) -> Form<'_> {
    let View::Apply { head, args } = e.view() else {
        unreachable!("only applications have a form");
    };
    let Some(symbol) = head.as_symbol() else {
        return Form::Call(head, args);
    };
    match (symbol, args) {
        (Symbol::PLUS, [_, _, ..]) => Form::Sum(args),
        (Symbol::TIMES, [_, _, ..]) => Form::Product(args),
        (Symbol::POWER, [_, n]) if is_negative_integer(n) => Form::Product(std::slice::from_ref(e)),
        (Symbol::POWER, [b, n]) => Form::Power(b, "^", n),
        (Symbol::ELEMENTWISE_POWER, [b, n]) => Form::Power(b, ".^", n),
        (Symbol::ELEMENTWISE_TIMES, [_, _, ..]) => Form::Elementwise(args, ".*"),
        (Symbol::ELEMENTWISE_DIVIDE, [_, _]) => Form::Elementwise(args, "./"),
        (Symbol::FACTORIAL, [n]) => Form::Factorial(n),
        (Symbol::AND, [_, _, ..]) if is_chain(args) => Form::Chain(args),
        (Symbol::AND, [_, _, ..]) => Form::Logical(args, " and ", Prec::And, Prec::Not),
        (Symbol::OR, [_, _, ..]) => Form::Logical(args, " or ", Prec::Or, Prec::And),
        (Symbol::NOT, [p]) => Form::Not(p),
        (Symbol::RANGE, [a, b]) => Form::Range(a, b),
        (Symbol::LIST, _) => Form::Bracketed(args, '[', ']'),
        (Symbol::FINITE_SET, _) => Form::Bracketed(args, '{', '}'),
        (Symbol::AT, [base, indexes @ ..]) => Form::Index(base, indexes),
        (Symbol::FUNCTION, [params, body]) if has_head(params, Symbol::LIST) => {
            Form::Lambda(params.args(), body)
        }
        (Symbol::IF, [_, _] | [_, _, _]) => Form::If(args),
        (Symbol::DIRECTED_INFINITY, [d]) if d.as_i64() == Some(-1) => Form::MinusInfinity,
        (_, [a, b]) => match relation(symbol) {
            Some(op) => Form::Relation(a, op, b),
            None => Form::Call(head, args),
        },
        _ => Form::Call(head, args),
    }
}

/// The operator of `Equation` and of the comparisons.
fn relation(s: Symbol) -> Option<&'static str> {
    Some(match s {
        Symbol::EQUATION => "=",
        Symbol::EQUAL => "==",
        Symbol::UNEQUAL => "!=",
        Symbol::LESS => "<",
        Symbol::LESS_EQUAL => "<=",
        Symbol::GREATER => ">",
        Symbol::GREATER_EQUAL => ">=",
        Symbol::ELEMENT => "in",
        _ => return None,
    })
}

/// Whether the arguments of an `And` are comparisons that share operands, as
/// lowering makes `0 < x <= 1`: `And[Less[0, x], LessEqual[x, 1]]`.
fn is_chain(tests: &[Expr]) -> bool {
    let comparison = |t: &Expr| {
        t.args().len() == 2
            && t.head()
                .and_then(Expr::as_symbol)
                .is_some_and(|s| s != Symbol::EQUATION && relation(s).is_some())
    };
    tests.iter().all(comparison) && tests.windows(2).all(|w| w[0].args()[1] == w[1].args()[0])
}

fn has_head(e: &Expr, s: Symbol) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(s)
}

fn is_negative_integer(e: &Expr) -> bool {
    match e.view() {
        View::Int(n) => n < 0,
        View::Integer(n) => n.is_negative(),
        _ => false,
    }
}

/// How `e` prints, from its structure. A part the language cannot write has
/// some shape: printing reports it.
fn shape(e: &Expr) -> Shape {
    let (prec, sign) = match e.view() {
        View::Symbol(Symbol::INFINITY) => (Prec::Prefix, true),
        View::Int(n) if n < 0 => (Prec::Prefix, true),
        View::Integer(n) if n.is_negative() => (Prec::Prefix, true),
        View::Rational(q) => (Prec::Product, q.is_negative()),
        View::Complex(z) => complex_shape(z),
        View::Apply { .. } => match form(e) {
            Form::Sum(_) => (Prec::Sum, false),
            Form::Product(factors) => {
                let p = product(factors);
                match (p.numerator.len(), p.denominator.is_empty()) {
                    (1, true) if p.negative && !p.writes_one() => (Prec::Prefix, true),
                    _ => (Prec::Product, p.negative),
                }
            }
            Form::Power(..) => (Prec::Power, false),
            Form::Elementwise(args, _) => {
                let first = shape(&args[0]);
                (Prec::Product, first.prec >= Prec::Prefix && first.sign)
            }
            Form::Relation(..) | Form::Chain(_) => (Prec::Compare, false),
            Form::Range(..) => (Prec::Range, false),
            Form::Not(_) => (Prec::Not, false),
            Form::Logical(_, _, prec, _) => (prec, false),
            Form::Lambda(..) | Form::If(_) => (Prec::Lowest, false),
            Form::MinusInfinity => (Prec::Prefix, true),
            Form::Factorial(_) | Form::Bracketed(..) | Form::Index(..) | Form::Call(..) => {
                (Prec::Postfix, false)
            }
        },
        _ => (Prec::Postfix, false),
    };
    Shape { prec, sign }
}

/// `i`, `-i`, `3*i`, `i/2`, `2 + 3*i`.
fn complex_shape(z: &Complex) -> (Prec, bool) {
    if !z.re().is_zero() {
        return (Prec::Sum, z.re().is_negative());
    }
    let negative = z.im().is_negative();
    match (z.im().abs().is_one(), negative) {
        (true, false) => (Prec::Postfix, false),
        (true, true) => (Prec::Prefix, true),
        (false, _) => (Prec::Product, negative),
    }
}

/// A product as printed: a sign, a numerator and a denominator.
struct Product<'a> {
    negative: bool,
    /// Never empty: `1` when nothing else is left.
    numerator: Vec<Item<'a>>,
    denominator: Vec<Item<'a>>,
}

impl Product<'_> {
    /// Whether the sign is written `-1*`: before an integer or `i`, a bare
    /// `-` would make one number of them (`-3`, `-i`), which reads back as
    /// that number rather than as the product `Times[-1, 3]`.
    fn writes_one(&self) -> bool {
        let absorbs_sign = |e: &Expr| match e.view() {
            View::Int(_) | View::Integer(_) => true,
            View::Complex(z) => z.re().is_zero() && z.im().abs().is_one(),
            _ => false,
        };
        self.negative && matches!(self.numerator.first(), Some(Item::Factor(f)) if absorbs_sign(f))
    }
}

/// A factor of a [`Product`].
enum Item<'a> {
    /// A positive integer, from the coefficient.
    Number(Integer),
    /// The imaginary unit, from the coefficient.
    Imaginary,
    Factor(&'a Expr),
    /// `b^n` for `n >= 2`, in a denominator.
    Power(&'a Expr, Integer),
}

/// The parts of the product of `factors`. A first factor that is a real or
/// an imaginary number is the coefficient: its sign goes in front, its
/// numerator and `i` to the numerator, its denominator to the denominator.
/// A factor `b^n` with an integer `n < 0` goes to the denominator as
/// `b^(-n)`.
fn product(factors: &[Expr]) -> Product<'_> {
    let mut p = Product {
        negative: false,
        numerator: Vec::new(),
        denominator: Vec::new(),
    };
    let mut rest = factors;
    if let Some((first, tail)) = factors.split_first() {
        let coefficient = match first.view() {
            View::Complex(z) if z.re().is_zero() => Some((z.im().clone(), true)),
            _ => first.to_rational().map(|q| (q, false)),
        };
        if let Some((q, imaginary)) = coefficient {
            rest = tail;
            p.negative = q.is_negative();
            let (numerator, denominator) = q.into_parts();
            let magnitude = numerator.abs();
            // `1*x` for `Times[1, x]`: without the 1, it would print as `x`.
            let one = !p.negative && !imaginary && denominator.is_one();
            if one || !magnitude.is_one() {
                p.numerator.push(Item::Number(magnitude));
            }
            if imaginary {
                p.numerator.push(Item::Imaginary);
            }
            if !denominator.is_one() {
                p.denominator.push(Item::Number(denominator));
            }
        }
    }
    for f in rest {
        match inverse_power(f) {
            Some((b, n)) if n.is_one() => p.denominator.push(Item::Factor(b)),
            Some((b, n)) => p.denominator.push(Item::Power(b, n)),
            None => p.numerator.push(Item::Factor(f)),
        }
    }
    if p.numerator.is_empty() {
        p.numerator.push(Item::Number(Integer::from(1)));
    }
    p
}

/// `(b, -n)` for `Power[b, n]` with an integer `n < 0`.
fn inverse_power(e: &Expr) -> Option<(&Expr, Integer)> {
    match e.view() {
        View::Apply { head, args: [b, n] }
            if head.as_symbol() == Some(Symbol::POWER) && is_negative_integer(n) =>
        {
            Some((b, -n.to_integer()?))
        }
        _ => None,
    }
}

/// What a subtracted term stands for: `x - 3`, `x - 3*i`, `x - 2*y`.
enum Magnitude<'a> {
    /// A positive rational.
    Number(Rational),
    /// A positive multiple of `i`.
    Imaginary(Rational),
    /// A product, without its sign.
    Product(Product<'a>),
}

/// The magnitude of a term that is written subtracted: a negative number, a
/// negative multiple of `i`, a product with a negative coefficient.
/// Lowering reads `a - t` as `a + (-1)*t`, which evaluates to the term.
fn negation(term: &Expr) -> Option<Magnitude<'_>> {
    match term.view() {
        View::Complex(z) if z.re().is_zero() && z.im().is_negative() => {
            Some(Magnitude::Imaginary(-z.im()))
        }
        View::Apply { .. } => match form(term) {
            Form::Product(factors) if has_head(term, Symbol::TIMES) => {
                let mut p = product(factors);
                p.negative.then(|| {
                    p.negative = false;
                    Magnitude::Product(p)
                })
            }
            _ => None,
        },
        _ => {
            let q = term.to_rational()?;
            q.is_negative().then(|| Magnitude::Number(-q))
        }
    }
}

/// Whether `s` occurs in `e` other than as the head of an application:
/// where it would make `s(e)` read as a product (CAS-DOC-7 §3.1). A head
/// written as `call(s, …)` would count too, but `s` is written so only when
/// it occurs in that call's argument already, or is a constant, which this
/// is never asked about.
fn occurs(s: Symbol, e: &Expr) -> bool {
    if !e.may_contain(s) {
        return false;
    }
    match e.view() {
        View::Symbol(t) => t == s,
        View::Apply { head, args } => {
            (head.as_symbol().is_none() && occurs(s, head)) || args.iter().any(|a| occurs(s, a))
        }
        _ => false,
    }
}

/// Whether `text` lexes as one name and is in the form [`name`] gives.
fn is_name(text: &str) -> bool {
    let lexed = lex(text);
    lexed.len() == 1
        && lexed.kind(0) == SyntaxKind::Ident
        && lexed.errors().is_empty()
        && name(text) == text
}

/// Whether an argument of `subst` would read as rules: an equation or a list
/// (CAS-DOC-7 §6).
fn reads_as_rules(e: &Expr) -> bool {
    has_head(e, Symbol::LIST) || (has_head(e, Symbol::EQUATION) && e.args().len() == 2)
}

/// The spelling of a target that every table names: the imaginary unit,
/// `infinity` and `call`.
fn spelled<'n>(names: &'n Names, target: &Target) -> &'n str {
    names
        .spelling(target)
        .expect("the English table names it, and locales only change spellings")
}

struct Printer<'a> {
    names: &'a Names,
    out: String,
}

impl<'a> Printer<'a> {
    fn write(&mut self, x: impl fmt::Display) {
        write!(self.out, "{x}").expect("writing to a String");
    }

    /// The name `s` is written with: a builtin's spelling, or a user's name
    /// without its `$`. `None` if the language has none, or if that name
    /// would read back as another symbol.
    fn symbol_name(&self, s: Symbol) -> Option<&'a str> {
        if s.is_builtin() {
            return self.names.spelling(&Target::Builtin(s));
        }
        let user = unescape_name(s.name());
        let reads_back =
            escape_name(user) == s.name() && is_name(user) && !self.names.is_reserved(user);
        reads_back.then_some(user)
    }

    /// Writes `e` where operands bind at least as tightly as `min`,
    /// parenthesised if it binds weaker, or if `no_sign` and it would start
    /// with a sign.
    fn operand(&mut self, e: &Expr, min: Prec, no_sign: bool) -> Result<(), PrintError> {
        let shape = shape(e);
        if shape.prec < min || (no_sign && shape.sign) {
            self.out.push('(');
            self.expr(e)?;
            self.out.push(')');
            Ok(())
        } else {
            self.expr(e)
        }
    }

    /// Writes `e` without parentheses around it.
    fn expr(&mut self, e: &Expr) -> Result<(), PrintError> {
        match e.view() {
            View::Symbol(Symbol::INFINITY) => {
                self.out.push('+');
                let infinity = Target::Builtin(Symbol::COMPLEX_INFINITY);
                self.out.push_str(spelled(self.names, &infinity));
            }
            View::Symbol(s) => match self.symbol_name(s) {
                Some(n) => self.out.push_str(n),
                None => {
                    return Err(PrintError {
                        part: e.clone(),
                        kind: PrintErrorKind::NoName,
                    });
                }
            },
            View::Int(n) => self.write(n),
            View::Integer(n) => self.write(n),
            View::Rational(q) => self.write(q),
            View::Complex(z) => self.complex(z),
            View::String(s) => self.string(s),
            View::Apply { .. } => return self.application(e),
            _ => {
                return Err(PrintError {
                    part: e.clone(),
                    kind: PrintErrorKind::NoSyntax,
                });
            }
        }
        Ok(())
    }

    /// `2 + 3*i`, `1/2 - i`, `-3*i/4`.
    fn complex(&mut self, z: &Complex) {
        let (re, im) = (z.re(), z.im());
        if re.is_zero() {
            if im.is_negative() {
                self.out.push('-');
            }
        } else {
            self.write(re);
            self.out
                .push_str(if im.is_negative() { " - " } else { " + " });
        }
        self.imaginary(&im.abs());
    }

    /// `m*i` for a positive rational `m`: `i`, `3*i`, `i/2`, `3*i/4`.
    fn imaginary(&mut self, m: &Rational) {
        if !m.numer().is_one() {
            self.write(m.numer());
            self.out.push('*');
        }
        self.out.push_str(spelled(self.names, &Target::Imaginary));
        if !m.denom().is_one() {
            self.out.push('/');
            self.write(m.denom());
        }
    }

    /// `"…"`, with the escapes of the lexer: `\"`, `\\`, `\n`, `\t`, and
    /// `\u{…}` for the other control characters.
    fn string(&mut self, s: &str) {
        self.out.push('"');
        for c in s.chars() {
            match c {
                '"' => self.out.push_str("\\\""),
                '\\' => self.out.push_str("\\\\"),
                '\n' => self.out.push_str("\\n"),
                '\t' => self.out.push_str("\\t"),
                c if c.is_control() => self.write(format_args!("\\u{{{:x}}}", u32::from(c))),
                c => self.out.push(c),
            }
        }
        self.out.push('"');
    }

    fn application(&mut self, e: &Expr) -> Result<(), PrintError> {
        match form(e) {
            Form::Sum(terms) => self.sum(terms),
            Form::Product(factors) => self.product(&product(factors)),
            Form::Power(base, op, exponent) => {
                self.base(base)?;
                self.out.push_str(op);
                self.operand(exponent, Prec::Power, true)
            }
            Form::Elementwise(operands, op) => {
                for (i, a) in operands.iter().enumerate() {
                    if i > 0 {
                        self.out.push_str(op);
                    }
                    self.operand(a, Prec::Prefix, i > 0)?;
                }
                Ok(())
            }
            Form::Factorial(n) => {
                self.operand(n, Prec::Postfix, false)?;
                self.out.push('!');
                Ok(())
            }
            Form::Relation(a, op, b) => {
                self.operand(a, Prec::Range, false)?;
                self.relation(op, b)
            }
            Form::Chain(tests) => {
                self.operand(&tests[0].args()[0], Prec::Range, false)?;
                for t in tests {
                    let op = t.head().and_then(Expr::as_symbol).and_then(relation);
                    self.relation(op.expect("a comparison"), &t.args()[1])?;
                }
                Ok(())
            }
            Form::Range(a, b) => {
                // `a + 1 .. b`: spaced when an end is a sum.
                let spaced = shape(a).prec == Prec::Sum || shape(b).prec == Prec::Sum;
                self.operand(a, Prec::Sum, false)?;
                self.out.push_str(if spaced { " .. " } else { ".." });
                self.operand(b, Prec::Sum, false)
            }
            Form::Not(p) => {
                self.out.push_str("not ");
                self.operand(p, Prec::Not, false)
            }
            Form::Logical(operands, op, _, min) => {
                for (i, a) in operands.iter().enumerate() {
                    if i > 0 {
                        self.out.push_str(op);
                    }
                    self.operand(a, min, false)?;
                }
                Ok(())
            }
            Form::Bracketed(items, open, close) => {
                self.out.push(open);
                self.comma_separated(items)?;
                self.out.push(close);
                Ok(())
            }
            Form::Index(base, indexes) => {
                self.operand(base, Prec::Postfix, false)?;
                self.out.push('[');
                self.comma_separated(indexes)?;
                self.out.push(']');
                Ok(())
            }
            Form::Lambda(params, body) => self.lambda(params, body),
            Form::If(args) => self.if_expr(args),
            Form::MinusInfinity => {
                self.out.push('-');
                let infinity = Target::Builtin(Symbol::COMPLEX_INFINITY);
                self.out.push_str(spelled(self.names, &infinity));
                Ok(())
            }
            Form::Call(head, args) => self.call(e, head, args),
        }
    }

    /// ` op b`, the right-hand side of a relation.
    fn relation(&mut self, op: &str, b: &Expr) -> Result<(), PrintError> {
        self.out.push(' ');
        self.out.push_str(op);
        self.out.push(' ');
        self.operand(b, Prec::Range, false)
    }

    fn comma_separated(&mut self, items: &[Expr]) -> Result<(), PrintError> {
        for (i, item) in items.iter().enumerate() {
            if i > 0 {
                self.out.push_str(", ");
            }
            self.expr(item)?;
        }
        Ok(())
    }

    /// Terms joined by `+`, or by `-` for the negative ones.
    fn sum(&mut self, terms: &[Expr]) -> Result<(), PrintError> {
        self.operand(&terms[0], Prec::Product, false)?;
        for t in &terms[1..] {
            if let Some(magnitude) = negation(t) {
                self.out.push_str(" - ");
                self.magnitude(&magnitude)?;
            } else {
                self.out.push_str(" + ");
                self.operand(t, Prec::Product, true)?;
            }
        }
        Ok(())
    }

    /// What follows the `-` of a subtracted term.
    fn magnitude(&mut self, m: &Magnitude<'_>) -> Result<(), PrintError> {
        match m {
            Magnitude::Number(q) => self.write(q),
            Magnitude::Imaginary(q) => self.imaginary(q),
            Magnitude::Product(p) => {
                return match (p.numerator.as_slice(), p.denominator.as_slice()) {
                    // The right operand of a binary `-` binds at least like a
                    // product.
                    ([item], []) => self.item(item, Prec::Product),
                    _ => self.product(p),
                };
            }
        }
        Ok(())
    }

    /// `-2*x*y/(3*z)`: the sign, the numerator, then the denominator,
    /// parenthesised when it has more than one factor.
    fn product(&mut self, p: &Product<'_>) -> Result<(), PrintError> {
        if p.negative {
            self.out.push('-');
            if p.writes_one() {
                self.out.push_str("1*");
            }
        }
        for (i, item) in p.numerator.iter().enumerate() {
            if i > 0 {
                self.out.push('*');
            }
            if i == 0 && p.negative && !p.writes_one() {
                // After the sign comes the operand of a prefix `-`; and
                // `-infinity` is `DirectedInfinity[-1]` (CAS-DOC-7 §5), so
                // the negated complex infinity is `-(infinity)`.
                if matches!(item, Item::Factor(f) if f.as_symbol() == Some(Symbol::COMPLEX_INFINITY))
                {
                    self.out.push('(');
                    self.item(item, Prec::Lowest)?;
                    self.out.push(')');
                } else {
                    self.item(item, Prec::Power)?;
                }
            } else {
                self.item(item, Prec::Prefix)?;
            }
        }
        if let [d] = p.denominator.as_slice() {
            self.out.push('/');
            self.item(d, Prec::Prefix)?;
        } else if !p.denominator.is_empty() {
            self.out.push_str("/(");
            for (i, d) in p.denominator.iter().enumerate() {
                if i > 0 {
                    self.out.push('*');
                }
                self.item(d, Prec::Prefix)?;
            }
            self.out.push(')');
        }
        Ok(())
    }

    fn item(&mut self, item: &Item<'_>, min: Prec) -> Result<(), PrintError> {
        match item {
            Item::Number(n) => self.write(n),
            Item::Imaginary => self.out.push_str(spelled(self.names, &Target::Imaginary)),
            Item::Factor(f) => return self.operand(f, min, true),
            Item::Power(b, n) => {
                self.base(b)?;
                self.out.push('^');
                self.write(n);
            }
        }
        Ok(())
    }

    /// The base of a power: parenthesised unless it binds like a postfix
    /// form, and always when it is a function's name, since `sin^2` is an
    /// error (CAS-DOC-7 §3.2).
    fn base(&mut self, b: &Expr) -> Result<(), PrintError> {
        let function = b
            .as_symbol()
            .is_some_and(|s| s.is_builtin() && !is_value(s));
        if function {
            self.out.push('(');
            self.expr(b)?;
            self.out.push(')');
            Ok(())
        } else {
            self.operand(b, Prec::Postfix, false)
        }
    }

    /// `f(args…)`, `f(a)(b)`, `v[0](x)`, `(x -> x^2)(3)`, or
    /// `call(head, args…)` when the head cannot be written in front of the
    /// parentheses (CAS-DOC-7 §3.1).
    fn call(&mut self, e: &Expr, head: &Expr, args: &[Expr]) -> Result<(), PrintError> {
        match head.view() {
            View::Symbol(s) if s != Symbol::INFINITY => {
                let Some(name) = self.symbol_name(s) else {
                    // An internal builtin applied without its syntax.
                    let (part, kind) = if s.is_builtin() {
                        (e, PrintErrorKind::NoSyntax)
                    } else {
                        (head, PrintErrorKind::NoName)
                    };
                    return Err(PrintError {
                        part: part.clone(),
                        kind,
                    });
                };
                let product =
                    args.len() == 1 && ((s.is_builtin() && is_value(s)) || occurs(s, &args[0]));
                if !product {
                    self.out.push_str(name);
                    return self.arguments(args, s == Symbol::REPLACE_ALL);
                }
            }
            View::Apply { .. } => match form(head) {
                Form::Call(..) | Form::Index(..) => {
                    self.expr(head)?;
                    return self.arguments(args, false);
                }
                Form::Lambda(..) => {
                    self.out.push('(');
                    self.expr(head)?;
                    self.out.push(')');
                    return self.arguments(args, false);
                }
                _ => {}
            },
            _ => {}
        }
        let call = Target::Syntax(SyntaxForm::Call);
        self.out.push_str(spelled(self.names, &call));
        self.out.push('(');
        self.expr(head)?;
        for a in args {
            self.out.push_str(", ");
            self.expr(a)?;
        }
        self.out.push(')');
        Ok(())
    }

    /// `(args…)`. After the first argument of `subst`, an equation or a list
    /// reads as rules; parentheses keep it a value.
    fn arguments(&mut self, args: &[Expr], rules: bool) -> Result<(), PrintError> {
        self.out.push('(');
        for (i, a) in args.iter().enumerate() {
            if i > 0 {
                self.out.push_str(", ");
            }
            if rules && i > 0 && reads_as_rules(a) {
                self.out.push('(');
                self.expr(a)?;
                self.out.push(')');
            } else {
                self.expr(a)?;
            }
        }
        self.out.push(')');
        Ok(())
    }

    /// `x -> body`, `(x, y) -> body`, `() -> body`.
    fn lambda(&mut self, params: &[Expr], body: &Expr) -> Result<(), PrintError> {
        let mut names = Vec::with_capacity(params.len());
        for p in params {
            let user = p.as_symbol().filter(|s| !s.is_builtin());
            match user.and_then(|s| self.symbol_name(s)) {
                Some(n) => names.push(n),
                None => {
                    return Err(PrintError {
                        part: p.clone(),
                        kind: PrintErrorKind::Parameter,
                    });
                }
            }
        }
        if let [n] = names.as_slice() {
            self.out.push_str(n);
        } else {
            self.out.push('(');
            self.out.push_str(&names.join(", "));
            self.out.push(')');
        }
        self.out.push_str(" -> ");
        self.expr(body)
    }

    /// `if c then a elif d then b else e end`: an `If` in the last branch is
    /// an `elif`.
    fn if_expr(&mut self, mut args: &[Expr]) -> Result<(), PrintError> {
        self.out.push_str("if ");
        loop {
            self.expr(&args[0])?;
            self.out.push_str(" then");
            self.branch(&args[1])?;
            match &args[2..] {
                [] => break,
                [rest] if !rest.is_atom() && matches!(form(rest), Form::If(_)) => {
                    self.out.push_str(" elif ");
                    args = rest.args();
                }
                [rest] => {
                    self.out.push_str(" else");
                    self.branch(rest)?;
                    break;
                }
                _ => unreachable!("If has two or three arguments"),
            }
        }
        self.out.push_str(" end");
        Ok(())
    }

    /// A branch: its value as a statement, nothing for `Null`.
    fn branch(&mut self, value: &Expr) -> Result<(), PrintError> {
        if value.as_symbol() != Some(Symbol::NULL) {
            self.out.push(' ');
            self.expr(value)?;
        }
        Ok(())
    }
}

/// The arithmetic that the text of [`print()`] writes out, folded. Products
/// are flattened, their numbers multiplied into one coefficient in front and
/// their other factors sorted. Reciprocals are distributed over the factors
/// of products and powers as written, before their numbers combine, and
/// taken of nonzero numbers: `1/(2*y^2)` gives `Times[1/2, Power[y, -2]]`.
/// Sums of numbers are added: `2 + 3*i`.
///
/// [`print()`] promises `fold_arithmetic(&read) == fold_arithmetic(e)` for
/// the expression `read` its text lowers to, and evaluation takes the same
/// steps. This is test support for that promise, for the unit tests and the
/// `lang_print` fuzz target, until evaluation exists to check it with.
#[doc(hidden)]
#[must_use]
pub fn fold_arithmetic(e: &Expr) -> Expr {
    let View::Apply { head, args } = e.view() else {
        return e.clone();
    };
    if let (Some(Symbol::POWER), [b, n]) = (head.as_symbol(), args) {
        if fold_arithmetic(n).as_i64() == Some(-1) {
            return reciprocal(b);
        }
    }
    let head = fold_arithmetic(head);
    let args: Vec<Expr> = args.iter().map(fold_arithmetic).collect();
    match (head.as_symbol(), args.as_slice()) {
        (Some(Symbol::TIMES), _) => fold_times(args),
        (Some(Symbol::PLUS), [_, _, ..]) if args.iter().all(Expr::is_number) => {
            let sum = args
                .iter()
                .filter_map(Expr::to_complex)
                .fold(Complex::from(0), |s, z| &s + &z);
            Expr::from(sum)
        }
        _ => Expr::apply(&head, &args),
    }
}

/// `Times[factors…]` of folded factors, folded: see [`fold_arithmetic`].
fn fold_times(factors: Vec<Expr>) -> Expr {
    let mut coefficient = Complex::from(1);
    let mut rest = Vec::new();
    for f in factors {
        let inner = if has_head(&f, Symbol::TIMES) {
            f.args().to_vec()
        } else {
            vec![f]
        };
        for g in inner {
            match g.to_complex() {
                Some(z) => coefficient = &coefficient * &z,
                None => rest.push(g),
            }
        }
    }
    rest.sort();
    let coefficient = Expr::from(coefficient);
    if rest.is_empty() {
        return coefficient;
    }
    if coefficient.as_i64() != Some(1) {
        rest.insert(0, coefficient);
    } else if rest.len() == 1 {
        return rest.pop().expect("one factor");
    }
    Expr::apply(&Symbol::TIMES.into(), &rest)
}

/// `Power[base, -1]`, folded: see [`fold_arithmetic`]. Distributing over
/// `base` as written keeps the factors a zero would absorb: `1/(2*(0*y))`
/// keeps its `1/2`.
fn reciprocal(base: &Expr) -> Expr {
    if let View::Apply { head, args } = base.view() {
        match (head.as_symbol(), args) {
            (Some(Symbol::TIMES), _) => return fold_times(args.iter().map(reciprocal).collect()),
            (Some(Symbol::POWER), [b, n]) => {
                let exponent = fold_times(vec![fold_arithmetic(n), Expr::from(-1)]);
                return match exponent.as_i64() {
                    Some(-1) => reciprocal(b),
                    Some(1) => fold_arithmetic(b),
                    _ => Expr::apply(&Symbol::POWER.into(), &[fold_arithmetic(b), exponent]),
                };
            }
            _ => {}
        }
    }
    let base = fold_arithmetic(base);
    match base.to_complex() {
        Some(z) if !z.is_zero() => Expr::from(z.recip()),
        _ => Expr::apply(&Symbol::POWER.into(), &[base, Expr::from(-1)]),
    }
}

#[cfg(test)]
mod tests;
