//! Our expressions as giac-shaped trees ([`Gen`]), in the shapes giac's
//! evaluator leaves (cookbook C0067): what `gen2markup` reads decides the
//! text, so `x/2` must be `prod(x, inv(2))` as giac evaluates it, not a
//! product with the factor 1/2.
//!
//! - `Times[-1, y]` is `neg(y)`; an integer coefficient a factor, its sign
//!   in it (`-2*x`).
//! - Factors `b^-n` go to a denominator: one factor `inv(den)` at the end,
//!   `den` a product when there are several (`a/(b*c)`); a lone `x^-2` is
//!   `inv(x^2)`.
//! - A rational coefficient `±1/q` puts `q` in the denominator (`x/2` is
//!   `\frac{x}{2}`); `p/q` with `p ≠ ±1` stays a fraction factor (`2/3*x` is
//!   `\frac{2}{3} x`) unless there is a denominator already, which then
//!   takes `q` (`2/(3*x)`).
//! - `x^(1/2)` is `sqrt(x)`, `x^(-1/2)` `sqrt(x)^-1`; `E^x` and `Exp[x]`
//!   are `exp(x)`, `E` `exp(1)`.
//! - Complex numbers are `_CPLX`; balls and decimals giac's doubles, with
//!   giac's text of the decimal they show (12 significant digits).
//! - Lists are vectors (equal-length lists of lists print as matrices),
//!   `Range` is `interval`, the relations giac's (`=` `equal`, `==`
//!   `same`), `D` `derive`, `Derivative[n][f][x]` at a symbol
//!   `derive(f(x), x, n)`, at another argument `f'(a)`, a partial one
//!   `derive(derive(f, 0), 1)(x, y)`.
//! - `True`, `Indeterminate`, `ComplexInfinity`, `Infinity` are `true`,
//!   `undef`, `infinity`, `+infinity`.
//! - Other builtins are giac's operator of that name (a table where ours
//!   differs, else the name's spelling), user functions `of(f, args)`.

use cas_core::{Expr, Integer, Rational, Symbol, View, ball_decimal};
use cas_num::ball::Ball;

use super::tree::{Gen, Sub};
use crate::names::{Names, Target, unescape_name};
use crate::print::{Balls, MAX_PRINT_DEPTH, PrintError, PrintErrorKind};

/// Lowers expressions with the names of `names`, balls shown as `balls`.
pub(super) struct Lowering<'a> {
    pub(super) names: &'a Names,
    pub(super) balls: Balls,
}

/// giac's operator for the builtins whose giac name is not our spelling,
/// or that markup.cc recognises by name (so a locale's spellings do not
/// change the text).
fn giac_op(s: Symbol) -> Option<&'static str> {
    Some(match s {
        Symbol::EXP => "exp",
        Symbol::LOG => "ln",
        Symbol::SIN => "sin",
        Symbol::COS => "cos",
        Symbol::TAN => "tan",
        Symbol::SEC => "sec",
        Symbol::CSC => "csc",
        Symbol::COT => "cot",
        Symbol::ARC_SIN => "asin",
        Symbol::ARC_COS => "acos",
        Symbol::ARC_TAN => "atan",
        Symbol::ARC_SEC => "asec",
        Symbol::ARC_CSC => "acsc",
        Symbol::ARC_COT => "acot",
        Symbol::SINH => "sinh",
        Symbol::COSH => "cosh",
        Symbol::TANH => "tanh",
        Symbol::ARC_SINH => "asinh",
        Symbol::ARC_COSH => "acosh",
        Symbol::ARC_TANH => "atanh",
        Symbol::ABS => "abs",
        Symbol::FLOOR => "floor",
        Symbol::CEILING => "ceil",
        Symbol::SQRT => "sqrt",
        Symbol::CONJUGATE => "conj",
        Symbol::FACTORIAL => "factorial",
        Symbol::BINOMIAL => "binomial",
        Symbol::SUM => "sum",
        Symbol::PRODUCT => "product",
        Symbol::INTEGRATE => "integrate",
        Symbol::LIMIT => "limit",
        Symbol::D => "derive",
        Symbol::GAMMA => "Gamma",
        Symbol::BETA => "Beta",
        Symbol::ZETA => "Zeta",
        Symbol::POLY_GAMMA => "Psi",
        Symbol::DIRAC_DELTA => "Dirac",
        Symbol::HEAVISIDE_THETA => "Heaviside",
        Symbol::AIRY_AI => "Airy_Ai",
        Symbol::AIRY_BI => "Airy_Bi",
        Symbol::EULER_PHI => "euler",
        Symbol::ERF => "erf",
        Symbol::ROOT_OF => "rootof",
        _ => return None,
    })
}

/// giac's operator of a relation.
fn relation(s: Symbol) -> Option<&'static str> {
    Some(match s {
        Symbol::EQUATION => "equal",
        Symbol::EQUAL => "same",
        Symbol::UNEQUAL => "different",
        Symbol::LESS => "inferieur_strict",
        Symbol::LESS_EQUAL => "inferieur_egal",
        Symbol::GREATER => "superieur_strict",
        Symbol::GREATER_EQUAL => "superieur_egal",
        _ => return None,
    })
}

/// A rational as giac's integer or fraction.
fn rational(q: &Rational) -> Gen {
    if q.denom().is_one() {
        Gen::Int(q.numer().clone())
    } else {
        Gen::Frac(
            Box::new(Gen::Int(q.numer().clone())),
            Box::new(Gen::Int(q.denom().clone())),
        )
    }
}

/// `(b, n)` for `Power[b, -n]` with an integer `n > 0`.
fn inverse_power(e: &Expr) -> Option<(&Expr, Integer)> {
    let View::Apply { head, args: [b, n] } = e.view() else {
        return None;
    };
    if head.as_symbol() != Some(Symbol::POWER) {
        return None;
    }
    let n = n.to_integer().filter(Integer::is_negative)?;
    Some((b, -n))
}

/// The digits and exponent of a decimal's text, `-` aside: `1.5e-7` is
/// `("15", -8)`.
fn parse_decimal(s: &str) -> Option<(String, i64)> {
    let (mantissa, exponent) = match s.find(['e', 'E']) {
        Some(i) => (
            &s[..i],
            s[i + 1..].trim_start_matches('+').parse::<i64>().ok()?,
        ),
        None => (s, 0),
    };
    let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = format!("{int}{frac}");
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let frac_len = i64::try_from(frac.len()).ok()?;
    Some((digits, exponent.checked_sub(frac_len)?))
}

/// giac's text of a double (`print_DOUBLE_`, gen.cc): C's `%.12g` of
/// `digits · 10^exponent`, and `.0` after an integer below 2^30: `3.0`,
/// `0.333333333333`, `1e-05`, `1.23e+20`, `12345678901`.
fn double_text(digits: &str, exponent: i64) -> String {
    const PRECISION: usize = 12;
    let digits = digits.trim_start_matches('0');
    if digits.is_empty() {
        return "0.0".to_owned();
    }
    let mut sig: Vec<u8> = digits.bytes().map(|c| c - b'0').collect();
    // The decimal exponent of the first digit.
    let mut x = exponent.saturating_add(i64::try_from(sig.len()).unwrap_or(i64::MAX) - 1);
    if sig.len() > PRECISION {
        let up = sig[PRECISION] >= 5;
        sig.truncate(PRECISION);
        if up {
            if let Some(k) = sig.iter().rposition(|&d| d != 9) {
                sig[k] += 1;
                sig.truncate(k + 1);
            } else {
                sig = vec![1];
                x = x.saturating_add(1);
            }
        }
    }
    while sig.last() == Some(&0) {
        sig.pop();
    }
    let s: String = sig.iter().map(|&d| char::from(b'0' + d)).collect();
    let text = if !(-4..12).contains(&x) {
        let (first, rest) = s.split_at(1);
        let point = if rest.is_empty() { "" } else { "." };
        let sign = if x < 0 { '-' } else { '+' };
        format!("{first}{point}{rest}e{sign}{:02}", x.unsigned_abs())
    } else if x < 0 {
        let zeros = usize::try_from(-x - 1).unwrap_or(0);
        format!("0.{}{s}", "0".repeat(zeros))
    } else {
        let int_len = usize::try_from(x + 1).unwrap_or(0);
        if s.len() > int_len {
            format!("{}.{}", &s[..int_len], &s[int_len..])
        } else {
            format!("{s}{}", "0".repeat(int_len - s.len()))
        }
    };
    if !text.contains(['.', 'e']) && text.parse::<u64>().is_ok_and(|v| v < 1 << 30) {
        text + ".0"
    } else {
        text
    }
}

impl Lowering<'_> {
    /// `e` as giac's tree.
    pub(super) fn tree(&self, e: &Expr) -> Result<Gen, PrintError> {
        let ball_error = || PrintError {
            part: e.clone(),
            kind: PrintErrorKind::Ball,
        };
        Ok(match e.view() {
            View::Symbol(s) => self.symbol(s),
            View::Int(n) => Gen::int(n),
            View::Integer(n) => Gen::Int(n.clone()),
            View::Rational(q) => rational(q),
            View::Complex(z) => Gen::Cplx(Box::new(rational(z.re())), Box::new(rational(z.im()))),
            View::Ball(b) => self.ball(b).ok_or_else(ball_error)?,
            View::ComplexBall(z) => Gen::Cplx(
                Box::new(self.ball(z.re()).ok_or_else(ball_error)?),
                Box::new(self.ball(z.im()).ok_or_else(ball_error)?),
            ),
            View::String(s) => Gen::Str(s.to_owned()),
            View::Apply { head, args } => self.apply(head, args)?,
            _ => {
                return Err(PrintError {
                    part: e.clone(),
                    kind: PrintErrorKind::NoSyntax,
                });
            }
        })
    }

    fn all(&self, args: &[Expr]) -> Result<Vec<Gen>, PrintError> {
        args.iter().map(|a| self.tree(a)).collect()
    }

    /// The argument of an operator: the one, or a sequence of several.
    fn arg(&self, args: &[Expr]) -> Result<Gen, PrintError> {
        match args {
            [a] => self.tree(a),
            _ => Ok(Gen::Vect(self.all(args)?, Sub::Seq)),
        }
    }

    /// giac's name for a builtin: [`giac_op`], else its spelling, else its
    /// internal name.
    fn op_name(&self, s: Symbol) -> String {
        giac_op(s)
            .or_else(|| self.names.spelling(&Target::Builtin(s)))
            .unwrap_or_else(|| s.name())
            .to_owned()
    }

    fn symbol(&self, s: Symbol) -> Gen {
        match s {
            Symbol::TRUE => Gen::Bool(true),
            Symbol::FALSE => Gen::Bool(false),
            Symbol::PI => Gen::idnt("pi"),
            Symbol::E => Gen::symb("exp", Gen::int(1)),
            Symbol::EULER_GAMMA => Gen::idnt("euler_gamma"),
            Symbol::I => Gen::i(),
            Symbol::INDETERMINATE => Gen::idnt("undef"),
            Symbol::COMPLEX_INFINITY => Gen::idnt("infinity"),
            Symbol::INFINITY => Gen::symb("plus", Gen::idnt("infinity")),
            s if s.is_builtin() => Gen::Idnt(self.op_name(s)),
            s => Gen::idnt(unescape_name(s.name())),
        }
    }

    /// A ball as a double: giac's text of the decimal it shows, and its
    /// sign. A ball that shows no digit is `0.0`. `None` for a decimal
    /// beyond `i64` exponents.
    fn ball(&self, b: &Ball) -> Option<Gen> {
        let neg = b.mid().is_negative();
        let magnitude = if neg { b.neg() } else { b.clone() };
        let shown = match self.balls {
            Balls::Exact { working_precision } => magnitude.to_literal(working_precision)?,
            Balls::Shown { working_precision } => ball_decimal(&magnitude, working_precision),
        };
        let text = if shown.starts_with('[') {
            "0.0".to_owned()
        } else {
            let (digits, exponent) = parse_decimal(&shown)?;
            double_text(&digits, exponent)
        };
        Some(Gen::Dbl { text, neg })
    }

    fn apply(&self, head: &Expr, args: &[Expr]) -> Result<Gen, PrintError> {
        let Some(s) = head.as_symbol() else {
            return self.applied(head, args);
        };
        if !s.is_builtin() {
            let f = Gen::idnt(unescape_name(s.name()));
            return Ok(Gen::symb_n("of", vec![f, self.arg(args)?]));
        }
        Ok(match (s, args) {
            (Symbol::PLUS, [a]) => Gen::symb("plus", self.tree(a)?),
            (Symbol::PLUS, _) => Gen::symb_n("plus", self.all(args)?),
            (Symbol::TIMES, _) => self.product(args)?,
            (Symbol::POWER, [b, n]) => self.power(b, n)?,
            (Symbol::LIST, _) => Gen::Vect(self.all(args)?, Sub::List),
            (Symbol::FINITE_SET, _) => Gen::Vect(self.all(args)?, Sub::Set),
            (Symbol::SEQUENCE, _) => Gen::Vect(self.all(args)?, Sub::Seq),
            (Symbol::RANGE, [a, b]) => Gen::symb_n("interval", vec![self.tree(a)?, self.tree(b)?]),
            (Symbol::AT, [b, idx @ ..]) if !idx.is_empty() => {
                Gen::symb_n("at", vec![self.tree(b)?, self.arg(idx)?])
            }
            (Symbol::DIRECTED_INFINITY, [d]) if d.as_i64() == Some(-1) => {
                Gen::symb("neg", Gen::idnt("infinity"))
            }
            (Symbol::DIRECTED_INFINITY, [d]) if d.as_i64() == Some(1) => {
                Gen::symb("plus", Gen::idnt("infinity"))
            }
            (Symbol::DECIMAL, [n, exponent])
                if n.to_integer().is_some() && exponent.as_i64().is_some() =>
            {
                let n = n.to_integer().expect("checked");
                let exponent = exponent.as_i64().expect("checked");
                Gen::Dbl {
                    text: double_text(&n.abs().to_string(), exponent),
                    neg: n.is_negative(),
                }
            }
            // `diff(e, [x, n])` is giac's `diff(e, x$n)`.
            (Symbol::D, [f, vars @ ..]) => {
                let mut items = vec![self.tree(f)?];
                for v in vars {
                    items.push(match v.view() {
                        View::Apply { head, args: [x, n] }
                            if head.as_symbol() == Some(Symbol::LIST) && n.as_i64().is_some() =>
                        {
                            Gen::symb_n("dollar", vec![self.tree(x)?, self.tree(n)?])
                        }
                        _ => self.tree(v)?,
                    });
                }
                Gen::symb_n("derive", items)
            }
            // `int(f, [x, a, b])` is giac's `int(f, x, a, b)`.
            (Symbol::INTEGRATE | Symbol::SUM | Symbol::PRODUCT, [f, spec, rest @ ..]) => {
                let mut items = vec![self.tree(f)?];
                match spec.view() {
                    View::Apply { head, args: bounds }
                        if head.as_symbol() == Some(Symbol::LIST) =>
                    {
                        items.extend(self.all(bounds)?);
                    }
                    _ => items.push(self.tree(spec)?),
                }
                items.extend(self.all(rest)?);
                Gen::symb_n(&self.op_name(s), items)
            }
            (s, [a, b]) if relation(s).is_some() => {
                let op = relation(s).expect("checked");
                Gen::symb_n(op, vec![self.tree(a)?, self.tree(b)?])
            }
            _ => Gen::Symb(self.op_name(s), Box::new(self.arg(args)?)),
        })
    }

    /// An application whose head is not a symbol: derivatives of a
    /// function, else `of(head, args)`.
    fn applied(&self, head: &Expr, args: &[Expr]) -> Result<Gen, PrintError> {
        if let Some(g) = self.derivative(head, args)? {
            return Ok(g);
        }
        Ok(Gen::symb_n("of", vec![self.tree(head)?, self.arg(args)?]))
    }

    /// `Derivative[n][f]` (`f'`, `f''`, …), `Derivative[n][f][x]` at a
    /// symbol `derive(f(x), x, n)` as giac evaluates `diff(f(x), x, n)`, at
    /// another argument `f'(a)`, and `Derivative[n1, n2][f][x, y]`
    /// `derive(derive(f, 0), 1)(x, y)`. `None` for another head, and for
    /// orders whose nesting would pass [`MAX_PRINT_DEPTH`].
    fn derivative(&self, head: &Expr, args: &[Expr]) -> Result<Option<Gen>, PrintError> {
        let orders = |d: &Expr| -> Option<Vec<u64>> {
            let View::Apply { head, args: ns } = d.view() else {
                return None;
            };
            if head.as_symbol() != Some(Symbol::DERIVATIVE) {
                return None;
            }
            ns.iter()
                .map(|n| n.as_i64().and_then(|n| u64::try_from(n).ok()))
                .collect()
        };
        let max = u64::try_from(MAX_PRINT_DEPTH).unwrap_or(u64::MAX);
        let nested = |f: Gen, n: u64| (0..n).fold(f, |g, _| Gen::symb("derive", g));
        // `Derivative[n][f]` on its own.
        if let (Some(ns), [f]) = (orders(head), args)
            && let [n] = ns[..]
            && n <= max
        {
            return Ok(Some(nested(self.tree(f)?, n)));
        }
        let View::Apply { head: d, args: [f] } = head.view() else {
            return Ok(None);
        };
        let Some(ns) = orders(d).filter(|ns| ns.len() == args.len()) else {
            return Ok(None);
        };
        if let ([n], [a]) = (&ns[..], args) {
            if a.as_symbol().is_some_and(|s| !s.is_builtin()) {
                let call = self.tree(&Expr::apply(f, args))?;
                let mut items = vec![call, self.tree(a)?];
                if *n >= 2 {
                    items.push(Gen::Int(Integer::from(*n)));
                }
                return Ok(Some(if *n == 0 {
                    items.swap_remove(0)
                } else {
                    Gen::symb_n("derive", items)
                }));
            }
            if *n > max {
                return Ok(None);
            }
            let g = nested(self.tree(f)?, *n);
            return Ok(Some(Gen::symb_n("of", vec![g, self.tree(a)?])));
        }
        if ns.iter().sum::<u64>() > max {
            return Ok(None);
        }
        let mut g = self.tree(f)?;
        for (k, &n) in ns.iter().enumerate() {
            for _ in 0..n {
                let k = Gen::Int(Integer::from(i64::try_from(k).unwrap_or(i64::MAX)));
                g = Gen::symb_n("derive", vec![g, k]);
            }
        }
        Ok(Some(Gen::symb_n("of", vec![g, self.arg(args)?])))
    }

    /// `b^n`: `sqrt(b)` for `n = 1/2`, `exp(n)` for `b = E`, `inv(b^-n)`
    /// for a negative integer `n`, `sqrt(b)^-1` for `n = -1/2`.
    fn power(&self, b: &Expr, n: &Expr) -> Result<Gen, PrintError> {
        let half = Rational::new(1.into(), 2.into());
        let q = n.to_rational();
        if q.as_ref() == Some(&half) {
            return Ok(Gen::symb("sqrt", self.tree(b)?));
        }
        if b.as_symbol() == Some(Symbol::E) {
            return Ok(Gen::symb("exp", self.tree(n)?));
        }
        if q.as_ref() == Some(&-&half) {
            let root = Gen::symb("sqrt", self.tree(b)?);
            return Ok(Gen::symb_n("pow", vec![root, Gen::int(-1)]));
        }
        if let Some(n) = n.to_integer().filter(Integer::is_negative) {
            return Ok(Gen::symb("inv", self.denominator(b, -n)?));
        }
        Ok(Gen::symb_n("pow", vec![self.tree(b)?, self.tree(n)?]))
    }

    /// `b^n` in a denominator: `b` for `n = 1`.
    fn denominator(&self, b: &Expr, n: Integer) -> Result<Gen, PrintError> {
        let b = self.tree(b)?;
        Ok(if n.is_one() {
            b
        } else {
            Gen::symb_n("pow", vec![b, Gen::Int(n)])
        })
    }

    /// `Times[…]` as giac evaluates the product (module doc).
    fn product(&self, factors: &[Expr]) -> Result<Gen, PrintError> {
        let (coefficient, rest) = match factors.split_first() {
            Some((c, rest)) if c.to_rational().is_some() => (c.to_rational(), rest),
            _ => (None, factors),
        };
        let (mut num, mut den) = (Vec::new(), Vec::new());
        for f in rest {
            match inverse_power(f) {
                Some((b, n)) => den.push(self.denominator(b, n)?),
                None => num.push(self.tree(f)?),
            }
        }
        let mut minus = false;
        if let Some(q) = coefficient {
            let (p, d) = q.into_parts();
            let negative = p.is_negative();
            if den.is_empty() && d.is_one() {
                if p.abs().is_one() {
                    minus = negative;
                } else {
                    num.insert(0, Gen::Int(p));
                }
            } else if den.is_empty() && p.abs().is_one() {
                den.push(Gen::Int(d));
                minus = negative;
            } else if den.is_empty() {
                num.insert(0, Gen::Frac(Box::new(Gen::Int(p)), Box::new(Gen::Int(d))));
            } else {
                if !d.is_one() {
                    den.insert(0, Gen::Int(d));
                }
                if !p.abs().is_one() || (negative && num.is_empty()) {
                    // giac's `-1/x` is `prod(-1, inv(x))`.
                    num.insert(0, Gen::Int(p));
                } else {
                    minus = negative;
                }
            }
        }
        match den.len() {
            0 => {}
            1 => num.push(Gen::symb("inv", den.remove(0))),
            _ => num.push(Gen::symb("inv", Gen::symb_n("prod", den))),
        }
        let body = match num.len() {
            0 => Gen::int(1),
            1 => num.remove(0),
            _ => Gen::symb_n("prod", num),
        };
        Ok(if minus { Gen::symb("neg", body) } else { body })
    }
}
