#![allow(clippy::cast_precision_loss, reason = "small test values")]

use proptest::prelude::*;

use crate::build::symbols;
use crate::{Expr, Integer, Rational, Symbol};

fn apply(head: Symbol, args: &[Expr]) -> Expr {
    Expr::apply(&head.into(), args)
}

fn plus(args: &[Expr]) -> Expr {
    apply(Symbol::PLUS, args)
}

fn times(args: &[Expr]) -> Expr {
    apply(Symbol::TIMES, args)
}

fn pow(b: &Expr, n: impl Into<Expr>) -> Expr {
    apply(Symbol::POWER, &[b.clone(), n.into()])
}

fn q(n: i64, d: i64) -> Expr {
    Expr::from(Rational::new(n.into(), d.into()))
}

fn int(n: i64) -> Expr {
    n.into()
}

/// Infix snapshots: `(expression, expected)`.
fn check(cases: &[(Expr, &str)]) {
    for (e, expected) in cases {
        assert_eq!(&e.to_string(), expected, "{e:?}");
    }
}

#[test]
fn sums_and_negative_terms() {
    let [x, y] = symbols(["x", "y"]);
    check(&[
        (&x + 2 * y.pow(3), "x + 2*y^3"),
        (plus(&[int(-3), x.clone()]), "-3 + x"),
        (plus(&[x.clone(), int(-3)]), "x - 3"),
        (times(&[int(-1), x.clone()]), "-x"),
        (times(&[int(-2), x.clone()]), "-2*x"),
        (plus(&[x.clone(), times(&[int(-1), y.clone()])]), "x - y"),
        (plus(&[x.clone(), times(&[int(-2), y.clone()])]), "x - 2*y"),
        (plus(&[times(&[int(-1), x.clone()]), y.clone()]), "-x + y"),
        (plus(&[x.clone(), int(i64::MIN)]), "x - 9223372036854775808"),
        (times(&[int(i64::MIN), x.clone()]), "-9223372036854775808*x"),
    ]);
}

#[test]
fn no_operand_starts_with_a_sign() {
    let [x, y] = symbols(["x", "y"]);
    let neg = |e: &Expr| times(&[int(-1), e.clone()]);
    check(&[
        (times(&[int(2), int(-3)]), "2*(-3)"),
        (times(&[x.clone(), neg(&y)]), "x*(-y)"),
        (pow(&int(-3), 2), "(-3)^2"),
        (pow(&neg(&x), 2), "(-x)^2"),
        (pow(&x, -1), "x^(-1)"),
        (pow(&x, neg(&y)), "x^(-y)"),
        (pow(&plus(&[neg(&x), y.clone()]), 2), "(-x + y)^2"),
        // Nested unary minus.
        (neg(&neg(&x)), "-(-x)"),
        (neg(&int(-3)), "-(-3)"),
        // Not `-3`, which is the integer.
        (neg(&int(3)), "-1*3"),
        (times(&[int(-1), int(0), x.clone()]), "-1*0*x"),
        (plus(&[x.clone(), neg(&int(3))]), "x - 1*3"),
        (times(&[int(-2), int(-3)]), "-2*(-3)"),
        (plus(&[x.clone(), neg(&neg(&y))]), "x - (-y)"),
        (plus(&[x.clone(), neg(&int(-3))]), "x - (-3)"),
        (neg(&plus(&[x.clone(), y.clone()])), "-(x + y)"),
        (times(&[int(2), plus(&[neg(&x), y.clone()])]), "2*(-x + y)"),
        (
            plus(&[x.clone(), plus(&[neg(&y), x.clone()])]),
            "x + (-y + x)",
        ),
    ]);
}

#[test]
fn powers() {
    let [x, y, z] = symbols(["x", "y", "z"]);
    check(&[
        (pow(&pow(&x, &y), &z), "(x^y)^z"),
        (pow(&x, pow(&y, &z)), "x^y^z"),
        (pow(&times(&[int(2), x.clone()]), 3), "(2*x)^3"),
        (pow(&plus(&[x.clone(), int(1)]), 2), "(x + 1)^2"),
        (pow(&Expr::call(Symbol::SIN, [&x]), 2), "Sin[x]^2"),
        (pow(&x, q(1, 2)), "x^(1/2)"),
        (pow(&x, q(-3, 2)), "x^(-3/2)"),
        (pow(&q(1, 2), &x), "(1/2)^x"),
        (pow(&q(-1, 2), &x), "(-1/2)^x"),
        // Not a rational: an integer over an integer stays a power.
        (pow(&x, times(&[int(1), pow(&int(2), -1)])), "x^(1*2^(-1))"),
    ]);
}

#[test]
fn numbers() {
    let [x, y] = symbols(["x", "y"]);
    let big = Expr::from(Integer::from(1i128 << 70));
    let neg_big = Expr::from(-Integer::from(1i128 << 70));
    check(&[
        (big.clone(), "1180591620717411303424"),
        (neg_big.clone(), "-1180591620717411303424"),
        (
            plus(&[x.clone(), neg_big.clone()]),
            "x - 1180591620717411303424",
        ),
        (
            times(&[int(2), neg_big.clone()]),
            "2*(-1180591620717411303424)",
        ),
        (times(&[neg_big, x.clone()]), "-1180591620717411303424*x"),
        (q(1, 2), "1/2"),
        (q(-3, 2), "-3/2"),
        (plus(&[x.clone(), q(1, 2)]), "x + 1/2"),
        (plus(&[x.clone(), q(-1, 2)]), "x - 1/2"),
        (plus(&[q(-1, 2), x.clone()]), "-1/2 + x"),
        (times(&[q(1, 2), x.clone()]), "1/2*x"),
        (times(&[q(-1, 2), x.clone()]), "-1/2*x"),
        (
            plus(&[y.clone(), times(&[q(-1, 2), x.clone()])]),
            "y - 1/2*x",
        ),
        (times(&[x.clone(), q(-1, 2)]), "x*(-1/2)"),
        (times(&[int(2), q(3, 4)]), "2*(3/4)"),
        (times(&[int(-1), q(1, 2)]), "-1*(1/2)"),
        (
            plus(&[y.clone(), times(&[int(-1), q(1, 2)])]),
            "y - 1*(1/2)",
        ),
        (times(&[q(1, 2), pow(&x, -1)]), "1/2/x"),
        (
            Expr::apply(&Expr::symbol("f"), &[q(1, 2), big]),
            "f[1/2, 1180591620717411303424]",
        ),
        // An integer over an integer would read as a rational.
        (times(&[int(3), pow(&int(2), -1)]), "3*2^(-1)"),
        (times(&[int(-3), pow(&int(2), -1)]), "-3*2^(-1)"),
        (times(&[int(-1), pow(&int(2), -1)]), "-2^(-1)"),
        (
            plus(&[x.clone(), times(&[int(-3), pow(&int(2), -1)])]),
            "x - 3*2^(-1)",
        ),
        (times(&[x.clone(), pow(&int(2), -1)]), "x/2"),
        (times(&[int(3), pow(&int(-2), -1)]), "3/(-2)"),
        (times(&[int(3), pow(&int(2), -2)]), "3/2^2"),
    ]);
    assert_eq!(format!("{:?}", plus(&[x, q(-3, 2)])), "Plus[x, -3/2]");
}

#[test]
fn quotients() {
    let [x, y, z] = symbols(["x", "y", "z"]);
    check(&[
        (times(&[x.clone(), pow(&y, -1)]), "x/y"),
        (times(&[pow(&y, -1), x.clone()]), "x/y"),
        // The sign is the first numerator factor, wherever it is.
        (times(&[pow(&y, -1), int(-2)]), "-2/y"),
        (
            plus(&[x.clone(), times(&[pow(&y, -1), int(-2)])]),
            "x - 2/y",
        ),
        (times(&[int(2), times(&[pow(&y, -1), int(-1)])]), "2*(-1/y)"),
        // A denominator group is printed like the product it reads as.
        (times(&[pow(&int(-1), -1), pow(&y, -1)]), "1/(-y)"),
        (
            times(&[x.clone(), pow(&int(-2), -1), pow(&y, -1)]),
            "x/(-2*y)",
        ),
        (
            times(&[x.clone(), pow(&y, -1), pow(&pow(&z, -1), -1)]),
            "x/(y/z)",
        ),
        // Only negative exponents: `y^0` is a factor.
        (times(&[x.clone(), pow(&y, 0)]), "x*y^0"),
        (times(&[x.clone(), pow(&y, -2)]), "x/y^2"),
        (times(&[x.clone(), pow(&y, -1), pow(&z, -3)]), "x/(y*z^3)"),
        (times(&[int(2), pow(&x, -1)]), "2/x"),
        (times(&[pow(&x, -1), pow(&y, -1)]), "1/(x*y)"),
        (times(&[int(-1), pow(&x, -1)]), "-1/x"),
        (times(&[int(-3), x.clone(), pow(&y, -1)]), "-3*x/y"),
        (
            times(&[x.clone(), pow(&plus(&[y.clone(), z.clone()]), -1)]),
            "x/(y + z)",
        ),
        (
            times(&[x.clone(), pow(&times(&[int(2), y.clone()]), -1)]),
            "x/(2*y)",
        ),
        (times(&[x.clone(), pow(&int(-2), -1)]), "x/(-2)"),
        (times(&[x.clone(), pow(&pow(&y, &z), -2)]), "x/(y^z)^2"),
        (
            plus(&[x.clone(), times(&[int(-1), pow(&y, -1)])]),
            "x - 1/y",
        ),
        (pow(&times(&[x.clone(), pow(&y, -1)]), 2), "(x/y)^2"),
        // Only integer exponents make a denominator.
        (times(&[x.clone(), pow(&y, &z)]), "x*y^z"),
        (&x / &y, "x/y"),
        (&x / (2 * &y), "x/(2*y)"),
    ]);
}

#[test]
fn full_form_where_infix_does_not_apply() {
    let [x, y, z] = symbols(["x", "y", "z"]);
    let f = Symbol::new("f");
    check(&[
        (Expr::call(f, [x.clone(), &y + &z]), "f[x, y + z]"),
        (f.call([-&x]), "f[-x]"),
        (
            Expr::apply(&f.call([&x]), std::slice::from_ref(&y)),
            "f[x][y]",
        ),
        (
            Expr::apply(&(&x + &y), std::slice::from_ref(&z)),
            "(x + y)[z]",
        ),
        (Expr::apply(&x.pow(2), &[]), "(x^2)[]"),
        (plus(std::slice::from_ref(&x)), "Plus[x]"),
        (plus(&[]), "Plus[]"),
        (times(&[int(-1)]), "Times[-1]"),
        (apply(Symbol::POWER, std::slice::from_ref(&x)), "Power[x]"),
        (
            apply(Symbol::POWER, &[x.clone(), y.clone(), z.clone()]),
            "Power[x, y, z]",
        ),
        (crate::build::list([int(1), -&x]), "{1, -x}"),
        (crate::build::list(Vec::<Expr>::new()), "{}"),
        // Unflattened operators keep their structure.
        (
            plus(&[x.clone(), plus(&[y.clone(), z.clone()])]),
            "x + (y + z)",
        ),
        (
            plus(&[plus(&[x.clone(), y.clone()]), z.clone()]),
            "(x + y) + z",
        ),
        (
            times(&[x.clone(), times(&[y.clone(), z.clone()])]),
            "x*(y*z)",
        ),
    ]);
}

#[test]
fn full_form() {
    let [x, y] = symbols(["x", "y"]);
    let e = Expr::apply(&(&x + 1), &[x.pow(-1), Expr::string("s"), y]);
    let expected = "Plus[1, x][Power[x, -1], \"s\", y]";
    assert_eq!(e.full_form().to_string(), expected);
    assert_eq!(format!("{e:?}"), expected);
    assert_eq!(format!("{:?}", apply(Symbol::LIST, &[])), "List[]");
}

#[test]
fn strings_are_escaped() {
    let s = Expr::string("a\"b\\c\nd\te\r\u{7}é");
    let expected = r#""a\"b\\c\nd\te\r\u{7}é""#;
    assert_eq!(s.to_string(), expected);
    assert_eq!(format!("{s:?}"), expected);
}

#[test]
fn deep_expressions_use_constant_stack() {
    let [x] = symbols(["x"]);
    let f = Symbol::new("f");
    let n = 100_000;
    let (mut args, mut heads, mut powers) = (x.clone(), x.clone(), x.clone());
    for _ in 0..n {
        args = f.call([args]);
        heads = Expr::apply(&heads, &[]);
        powers = pow(&powers, 2);
    }
    assert_eq!(args.to_string().len(), 1 + 3 * n);
    assert_eq!(format!("{args:?}").len(), 1 + 3 * n);
    assert_eq!(heads.to_string().len(), 1 + 2 * n);
    assert_eq!(format!("{heads:?}").len(), 1 + 2 * n);
    // `(…((x^2)^2)…)^2`
    assert_eq!(powers.to_string().len(), 1 + 4 * n - 2);
    assert_eq!(format!("{powers:?}").len(), 1 + 10 * n);
}

/// Parses the infix output back, strictly: an operand that starts with a
/// sign after an operator is a syntax error. Builds with [`Expr::apply`], so
/// nothing is normalised.
struct Parser<'a> {
    tokens: Vec<&'a str>,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn parse(s: &'a str) -> Result<Expr, String> {
        let mut tokens = Vec::new();
        let mut rest = s;
        while let Some(c) = rest.chars().next() {
            let len = if c.is_ascii_alphanumeric() {
                rest.find(|c: char| !c.is_ascii_alphanumeric())
                    .unwrap_or(rest.len())
            } else {
                1
            };
            if c != ' ' {
                tokens.push(&rest[..len]);
            }
            rest = &rest[len..];
        }
        let mut p = Parser { tokens, pos: 0 };
        let e = p.sum()?;
        match p.peek() {
            None => Ok(e),
            Some(t) => Err(format!("trailing {t:?} in {s:?}")),
        }
    }

    fn peek(&self) -> Option<&'a str> {
        self.tokens.get(self.pos).copied()
    }

    fn eat(&mut self, t: &str) -> bool {
        let found = self.peek() == Some(t);
        self.pos += usize::from(found);
        found
    }

    fn expect(&mut self, t: &str) -> Result<(), String> {
        if self.eat(t) {
            Ok(())
        } else {
            Err(format!("expected {t:?} at {:?}", self.peek()))
        }
    }

    fn sum(&mut self) -> Result<Expr, String> {
        let negated = self.eat("-");
        let mut terms = vec![self.product(negated)?];
        loop {
            if self.eat("+") {
                terms.push(self.product(false)?);
            } else if self.eat("-") {
                terms.push(self.product(true)?);
            } else {
                break;
            }
        }
        Ok(if terms.len() == 1 {
            terms.pop().unwrap()
        } else {
            plus(&terms)
        })
    }

    /// A product, after a minus sign if `negated`: the inverse of how the
    /// printer writes a negative coefficient.
    fn product(&mut self, negated: bool) -> Result<Expr, String> {
        let literal = self.at_literal();
        let mut factors = vec![self.power()?];
        // A literal over a literal is a rational.
        if literal && self.peek() == Some("/") {
            self.pos += 1;
            if self.at_literal() {
                let d = self.power()?.to_integer().unwrap();
                let n = factors[0].to_integer().unwrap();
                factors[0] = Expr::from(Rational::new(n, d));
            } else {
                factors.push(pow(&self.power()?, -1));
            }
        }
        loop {
            if self.eat("*") {
                factors.push(self.power()?);
            } else if self.eat("/") {
                factors.push(pow(&self.power()?, -1));
            } else {
                break;
            }
        }
        if negated {
            match factors[0].to_rational() {
                Some(q) if q.is_positive() => factors[0] = Expr::from(-q),
                _ => factors.insert(0, int(-1)),
            }
        }
        Ok(if factors.len() == 1 {
            factors.pop().unwrap()
        } else {
            times(&factors)
        })
    }

    /// Whether the next tokens are a bare integer literal, not the base of a
    /// power or the head of an application.
    fn at_literal(&self) -> bool {
        let digits = self
            .peek()
            .is_some_and(|t| t.starts_with(|c: char| c.is_ascii_digit()));
        let next = self.tokens.get(self.pos + 1).copied();
        digits && next != Some("^") && next != Some("[")
    }

    fn power(&mut self) -> Result<Expr, String> {
        let base = self.postfix()?;
        if self.eat("^") {
            Ok(pow(&base, self.power()?))
        } else {
            Ok(base)
        }
    }

    fn postfix(&mut self) -> Result<Expr, String> {
        let mut e = self.primary()?;
        while self.eat("[") {
            e = Expr::apply(&e, &self.items("]")?);
        }
        Ok(e)
    }

    fn primary(&mut self) -> Result<Expr, String> {
        let t = self.peek().ok_or("unexpected end")?;
        self.pos += 1;
        match t {
            "(" => {
                let e = self.sum()?;
                self.expect(")")?;
                Ok(e)
            }
            "{" => Ok(apply(Symbol::LIST, &self.items("}")?)),
            _ if t.starts_with(|c: char| c.is_ascii_digit()) => t
                .parse::<Integer>()
                .map(Expr::from)
                .map_err(|e| e.to_string()),
            _ if t.starts_with(|c: char| c.is_ascii_alphabetic()) => Ok(Expr::symbol(t)),
            _ => Err(format!("unexpected {t:?}")),
        }
    }

    fn items(&mut self, close: &str) -> Result<Vec<Expr>, String> {
        let mut items = Vec::new();
        if self.eat(close) {
            return Ok(items);
        }
        loop {
            items.push(self.sum()?);
            if self.eat(close) {
                return Ok(items);
            }
            self.expect(",")?;
        }
    }
}

/// A value for `e`: arithmetic for `Plus`, `Times` and two-argument `Power`;
/// any other application is an order-sensitive function of its head and
/// arguments.
fn value(e: &Expr) -> f64 {
    let sym = |s: Symbol| match s.name() {
        "a" => 1.5,
        "b" => 0.625,
        "c" => 2.75,
        name => name.len() as f64,
    };
    if let Some(q) = e.to_rational() {
        let f = |n: &Integer| n.to_string().parse::<f64>().unwrap();
        return f(q.numer()) / f(q.denom());
    }
    if let Some(s) = e.as_symbol() {
        return sym(s);
    }
    let (head, args) = (e.head().unwrap(), e.args());
    let values = args.iter().map(value);
    match (head.as_symbol(), args.len()) {
        (Some(Symbol::PLUS), _) => values.sum(),
        (Some(Symbol::TIMES), _) => values.product(),
        (Some(Symbol::POWER), 2) => value(&args[0]).powf(value(&args[1])),
        _ => (1..)
            .zip(values)
            .fold(value(head), |acc, (i, v)| acc + f64::from(i) * v),
    }
}

fn same_value(a: f64, b: f64) -> bool {
    if a.is_finite() && b.is_finite() {
        (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
    } else {
        // Division by zero or overflow on both sides.
        !a.is_finite() && !b.is_finite()
    }
}

/// Trees over `Plus`, `Times`, `Power` (often to a negative integer), `List`
/// and user functions, with compound heads.
fn tree() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "b", "c"]).prop_map(Expr::symbol),
        (-3i64..4).prop_map(Expr::from),
        (-3i64..4, 1i64..4).prop_map(|(n, d)| q(n, d)),
        prop::sample::select(vec![-(1i128 << 64), 1i128 << 63])
            .prop_map(|n| Expr::from(Integer::from(n))),
    ];
    leaf.prop_recursive(5, 40, 4, |inner| {
        let head = prop_oneof![
            4 => prop::sample::select(vec![
                Symbol::PLUS,
                Symbol::TIMES,
                Symbol::POWER,
                Symbol::LIST,
                Symbol::new("f"),
            ])
            .prop_map(Expr::from),
            1 => inner.clone(),
        ];
        prop_oneof![
            3 => (head, prop::collection::vec(inner.clone(), 0..4))
                .prop_map(|(h, a)| Expr::apply(&h, &a)),
            1 => (inner.clone(), -3i64..0).prop_map(|(b, n)| pow(&b, n)),
            1 => (prop::collection::vec(inner, 1..4), -3i64..0).prop_map(|(mut f, c)| {
                f.insert(0, int(c));
                times(&f)
            }),
        ]
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(3000))]

    #[test]
    fn infix_reparses_to_the_same_text_and_value(e in tree()) {
        let s = e.to_string();
        let parsed = Parser::parse(&s).map_err(TestCaseError::fail)?;
        prop_assert_eq!(parsed.to_string(), s.clone(), "{:?}", e);
        let (v, w) = (value(&e), value(&parsed));
        prop_assert!(same_value(v, w), "{} for {:?}: {} vs {}", s, e, v, w);
    }

    #[test]
    fn no_sign_follows_an_operator(e in tree()) {
        let s = e.to_string();
        for bad in ["--", "+ -", "- -", "*-", "/-", "^-"] {
            prop_assert!(!s.contains(bad), "{:?} in {}", bad, s);
        }
    }

    #[test]
    fn full_form_matches_a_recursive_model(e in tree()) {
        fn model(e: &Expr) -> String {
            match (e.head(), e.to_rational(), e.as_symbol()) {
                (Some(h), ..) => {
                    let args: Vec<String> = e.args().iter().map(model).collect();
                    format!("{}[{}]", model(h), args.join(", "))
                }
                (_, Some(q), _) if q.is_integer() => q.numer().to_string(),
                (_, Some(q), _) => format!("{}/{}", q.numer(), q.denom()),
                (_, _, Some(s)) => s.name().to_owned(),
                _ => unreachable!(),
            }
        }
        prop_assert_eq!(e.full_form().to_string(), model(&e));
    }
}
