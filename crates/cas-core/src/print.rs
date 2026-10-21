//! Printing: full form (`Debug`, [`Expr::full_form`]) and infix (`Display`).
//!
//! Both printers run on an explicit stack, so expressions of any depth print
//! in constant native stack.

use std::fmt::{self, Write as _};

use cas_num::ball::{Ball, ComplexBall, DecimalStyle};

use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// An expression in full form, `Plus[x, Times[2, Power[y, 3]]]`: every
/// application as `head[args…]`, and a complex number as `Complex[re, im]`;
/// unambiguous. See [`Expr::full_form`].
#[derive(Clone, Copy, Debug)]
pub struct FullForm<'a>(&'a Expr);

impl Expr {
    /// The full form, for `Display`: `format!("{}", e.full_form())`. The same
    /// as `Debug`.
    ///
    /// ```
    /// use cas_core::build::symbols;
    ///
    /// let [x, y] = symbols(["x", "y"]);
    /// let e = &x + 2 * y.pow(3);
    /// assert_eq!(e.full_form().to_string(), "Plus[x, Times[2, Power[y, 3]]]");
    /// assert_eq!(format!("{e:?}"), "Plus[x, Times[2, Power[y, 3]]]");
    /// ```
    #[must_use]
    pub fn full_form(&self) -> FullForm<'_> {
        FullForm(self)
    }
}

impl fmt::Display for FullForm<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut stack = vec![Piece::Expr(self.0)];
        while let Some(piece) = stack.pop() {
            match piece {
                Piece::Text(text) => f.write_str(text)?,
                Piece::Expr(e) => match e.view() {
                    View::Apply { head, args } => {
                        stack.push(Piece::Text("]"));
                        for (i, a) in args.iter().enumerate().rev() {
                            stack.push(Piece::Expr(a));
                            if i > 0 {
                                stack.push(Piece::Text(", "));
                            }
                        }
                        stack.extend([Piece::Text("["), Piece::Expr(head)]);
                    }
                    View::Complex(z) => write!(f, "Complex[{}, {}]", z.re(), z.im())?,
                    // Exact: hexadecimal floats (CAS-DOC-5 §6.2).
                    View::Ball(b) => full_form_ball(b, f)?,
                    View::ComplexBall(z) => {
                        f.write_str("Complex[")?;
                        full_form_ball(z.re(), f)?;
                        f.write_str(", ")?;
                        full_form_ball(z.im(), f)?;
                        f.write_str("]")?;
                    }
                    _ => atom(e, f)?,
                },
            }
        }
        Ok(())
    }
}

enum Piece<'a> {
    Text(&'static str),
    Expr(&'a Expr),
}

/// `Ball[mid, rad]` with both in hexadecimal, as the reader reads them back.
fn full_form_ball(b: &Ball, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "Ball[{}, {}]", b.mid().to_hex(), b.rad().to_hex())
}

/// The precision balls print at without a context: 53 bits, the default
/// working precision.
const DISPLAY_PRECISION: u64 = 53;

/// A ball in decimal (CAS-DOC-5 §6.1, revised 2026-10-20). Known are the
/// digits correct to ±1 in the last place, at most those of the larger of
/// `prec` and the ball's accuracy, in bits. It shows the shortest decimal
/// inside it when that has no more digits than are known, and the known
/// digits otherwise: `0.1 + 0.2` shows as `0.3` and `1/3.0` as
/// `0.333333333333333`. Always with a point: `2.0`.
///
/// ```
/// use cas_core::ball_decimal;
/// use cas_num::ball::Ball;
///
/// let tenth = Ball::from_decimal("0.1", 53).unwrap();
/// assert_eq!(ball_decimal(&tenth, 53), "0.1");
/// let third = Ball::from_decimal("0.3", 53).unwrap();
/// let sum = tenth.add(&Ball::from_decimal("0.2", 53).unwrap(), 69).round_to_accuracy(53);
/// assert_eq!(ball_decimal(&sum, 53), "0.3");
/// assert!(sum.overlaps(&third));
/// ```
#[must_use]
pub fn ball_decimal(b: &Ball, prec: u64) -> String {
    let accuracy = u64::try_from(b.rel_accuracy_bits()).unwrap_or(0);
    let bits = if b.is_exact() {
        prec
    } else {
        prec.max(accuracy)
    };
    let digits = (bits.saturating_mul(30_103) / 100_000).max(1);
    let known = b.to_decimal(digits, DecimalStyle::Digits);
    let count = significant_digits(&known);
    if count > 0 {
        if let Some(shortest) = b.to_shortest_decimal(count) {
            return shortest;
        }
    }
    with_point(known)
}

/// The significant digits of a decimal `to_decimal` printed: leading zeros
/// aside, trailing ones counted; 0 for `[± R]`.
fn significant_digits(s: &str) -> u64 {
    if s.starts_with('[') {
        return 0;
    }
    let mantissa = s.split(['e', 'E']).next().unwrap_or(s);
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    u64::try_from(digits.trim_start_matches('0').len()).unwrap_or(u64::MAX)
}

/// A decimal with a point in its mantissa: `2` is `2.0`, `1e+20` is
/// `1.0e+20`.
fn with_point(s: String) -> String {
    if s.starts_with('[') {
        return s;
    }
    let (mantissa, exponent) = match s.find(['e', 'E']) {
        Some(i) => s.split_at(i),
        None => (s.as_str(), ""),
    };
    if mantissa.contains('.') {
        return s;
    }
    format!("{mantissa}.0{exponent}")
}

/// A complex ball as `re + im*I`, its parts as [`ball_decimal`] writes
/// them, the real part left out when it is exactly 0.
fn complex_ball_decimal(z: &ComplexBall, prec: u64) -> String {
    let im = ball_decimal(z.im(), prec);
    if z.re().is_exact() && z.re().mid().is_zero() {
        return format!("{im}*I");
    }
    let re = ball_decimal(z.re(), prec);
    match im.strip_prefix('-') {
        Some(magnitude) => format!("{re} - {magnitude}*I"),
        None => format!("{re} + {im}*I"),
    }
}

/// Symbols by name, numbers in decimal (`-3/2` for a rational, `1 - 3/4*I`
/// for a complex number), strings quoted and escaped.
fn atom(e: &Expr, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match e.view() {
        View::Symbol(s) => f.write_str(s.name()),
        View::Int(x) => write!(f, "{x}"),
        View::Integer(x) => write!(f, "{x}"),
        View::Rational(x) => write!(f, "{x}"),
        View::Complex(z) => write!(f, "{z}"),
        View::Ball(b) => f.write_str(&ball_decimal(b, DISPLAY_PRECISION)),
        View::ComplexBall(z) => f.write_str(&complex_ball_decimal(z, DISPLAY_PRECISION)),
        View::String(s) => string(s, f),
        View::Apply { .. } => unreachable!("not an atom"),
    }
}

/// `"…"`, escaping `\`, `"` and control characters.
fn string(s: &str, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    f.write_char('"')?;
    for c in s.chars() {
        match c {
            '"' => f.write_str("\\\"")?,
            '\\' => f.write_str("\\\\")?,
            '\n' => f.write_str("\\n")?,
            '\t' => f.write_str("\\t")?,
            '\r' => f.write_str("\\r")?,
            c if c.is_control() => write!(f, "\\u{{{:x}}}", u32::from(c))?,
            c => f.write_char(c)?,
        }
    }
    f.write_char('"')
}

/// Full form.
impl fmt::Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.full_form(), f)
    }
}

/// Infix: `x + 2*y^3`.
///
/// `Plus`, `Times` and `Power` print as operators, with parentheses only
/// where precedence needs them; `List` prints as `{a, b}`; everything else
/// in full form, `f[x, y]`.
///
/// - Numbers print in decimal, a rational as `3/2`, which binds like a
///   quotient: `x^(1/2)`, `2*(3/4)`, and `1/2*x` as a leading coefficient.
/// - A complex number prints as `2 + 3*I`, `3*I` or `I`, and binds like the
///   sum, product or atom it reads as: `x*(2 + 3*I)`, `(3*I)^2`, `I^2`. The
///   symbol `I` prints as `I` too; full form tells them apart.
/// - A term with a negative coefficient is subtracted: `x - 2*y`, and
///   `Times[-1, x]` is `-x`.
/// - Factors with a negative integer exponent are a denominator:
///   `x/(y*z^2)`, except that an integer over an integer would read as a
///   rational, so `Times[3, Power[2, -1]]` is `3*2^(-1)`. Likewise
///   `Times[-1, 3]` is `-1*3`, apart from the integer `-3`.
/// - An operand never starts with a minus sign after an operator:
///   `2*(-3)`, `(-3)^2`, `x^(-1)`, `-(-x)`.
/// - Nested applications of the same operator, which the builder flattens,
///   are parenthesised: `Plus[x, Plus[y, z]]` is `x + (y + z)`.
/// - `Plus` and `Times` with fewer than two arguments, and `Power` with
///   other than two, print in full form.
///
/// ```
/// use cas_core::build::symbols;
///
/// let [x, y, z] = symbols(["x", "y", "z"]);
/// assert_eq!((&x + 2 * y.pow(3)).to_string(), "x + 2*y^3");
/// assert_eq!((&x - &y / z.pow(2)).to_string(), "x - y/z^2");
/// assert_eq!((-(&x + &y)).pow(-1).to_string(), "(-(x + y))^(-1)");
/// ```
impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut stack = vec![Task::Expr {
            e: self,
            min: Prec::Lowest,
            no_minus: false,
        }];
        let mut expansion = Vec::new();
        while let Some(task) = stack.pop() {
            match task {
                Task::Text(s) => f.write_str(s)?,
                Task::Abs(n) => write!(f, "{}", n.unsigned_abs())?,
                Task::Magnitude(e) => match e.view() {
                    View::Int(n) => write!(f, "{}", n.unsigned_abs())?,
                    View::Integer(n) => write!(f, "{}", n.unsigned_abs())?,
                    View::Rational(q) => write!(f, "{}", q.abs())?,
                    // Only a negative imaginary number has a sign to drop.
                    View::Complex(z) => write!(f, "{}", -z)?,
                    View::Ball(b) => f.write_str(&ball_decimal(&b.neg(), DISPLAY_PRECISION))?,
                    View::ComplexBall(z) => {
                        f.write_str(&complex_ball_decimal(&z.neg(), DISPLAY_PRECISION))?;
                    }
                    _ => unreachable!("only numbers have a magnitude"),
                },
                Task::Expr { e, min, no_minus } => {
                    if prec(e) < min || (no_minus && leading_minus(e)) {
                        expansion.extend([
                            Task::Text("("),
                            Task::Expr {
                                e,
                                min: Prec::Lowest,
                                no_minus: false,
                            },
                            Task::Text(")"),
                        ]);
                    } else if e.is_atom() {
                        atom(e, f)?;
                    } else {
                        expand(e, &mut expansion);
                    }
                }
                Task::Negated(e) => match form(e) {
                    Form::Number => expansion.push(Task::Magnitude(e)),
                    Form::Times(args) => times(args, true, &mut expansion),
                    _ => unreachable!("only integers and products are negated"),
                },
                Task::Product { factors, negated } => product(&factors, negated, &mut expansion),
            }
            stack.extend(expansion.drain(..).rev());
        }
        Ok(())
    }
}

/// Binding strength, weakest first.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Prec {
    Lowest,
    Sum,
    /// Products and quotients.
    Product,
    Power,
    /// Atoms, lists and applications in full form.
    Atom,
}

enum Task<'a> {
    Text(&'static str),
    /// The absolute value of an integer.
    Abs(i64),
    /// The absolute value of a number atom.
    Magnitude(&'a Expr),
    /// `e`, parenthesised if it binds weaker than `min`, or if `no_minus`
    /// and it would start with a minus sign.
    Expr {
        e: &'a Expr,
        min: Prec,
        no_minus: bool,
    },
    /// `e` without its leading minus sign, which the caller has written. Only
    /// for a number or a product with [`leading_minus`].
    Negated(&'a Expr),
    /// The product of `factors`, see [`product`].
    Product {
        factors: Vec<Factor<'a>>,
        negated: bool,
    },
}

/// What the infix printer sees in an expression.
enum Form<'a> {
    /// A number atom.
    Number,
    Plus(&'a [Expr]),
    Times(&'a [Expr]),
    Power(&'a Expr, &'a Expr),
    List(&'a [Expr]),
    /// Other atoms, and applications in full form.
    Other,
}

fn form(e: &Expr) -> Form<'_> {
    match e.view() {
        View::Int(_)
        | View::Integer(_)
        | View::Rational(_)
        | View::Complex(_)
        | View::Ball(_)
        | View::ComplexBall(_) => Form::Number,
        View::Apply { head, args } => match (head.as_symbol(), args) {
            (Some(Symbol::PLUS), [_, _, ..]) => Form::Plus(args),
            (Some(Symbol::TIMES), [_, _, ..]) => Form::Times(args),
            (Some(Symbol::POWER), [b, n]) => Form::Power(b, n),
            (Some(Symbol::LIST), _) => Form::List(args),
            _ => Form::Other,
        },
        _ => Form::Other,
    }
}

fn prec(e: &Expr) -> Prec {
    match form(e) {
        Form::Plus(_) => Prec::Sum,
        Form::Times(_) => Prec::Product,
        Form::Power(..) => Prec::Power,
        Form::Number => number_prec(e),
        Form::List(_) | Form::Other => Prec::Atom,
    }
}

/// A rational reads as a quotient, and a complex number as the sum
/// (`2 + 3*I`), product (`3*I`) or atom (`I`) it prints as.
fn number_prec(e: &Expr) -> Prec {
    match e.view() {
        View::Rational(_) => Prec::Product,
        View::Complex(z) if !z.re().is_zero() => Prec::Sum,
        View::Complex(z) if !z.im().abs().is_one() => Prec::Product,
        View::ComplexBall(z) if z.re().is_exact() && z.re().mid().is_zero() => Prec::Product,
        View::ComplexBall(_) => Prec::Sum,
        _ => Prec::Atom,
    }
}

/// Whether `e`, printed without parentheses, starts with a minus sign. Only
/// numbers and products can: a sum with a negative first term is an operand
/// only where its precedence already parenthesises it. Precedence ignores the
/// sign; operands that must not start with one say so with `no_minus`.
fn leading_minus(e: &Expr) -> bool {
    match form(e) {
        Form::Number => is_negative(e),
        Form::Times(args) => args
            .iter()
            .find(|f| inverse_power(f).is_none())
            .is_some_and(is_negative),
        _ => false,
    }
}

/// Whether `e` is a negative number, or a negative imaginary one (`-3*I`),
/// which also prints with a sign that a product or a sum can take out. A
/// complex number with a real part binds like a sum, so precedence already
/// parenthesises it where a sign would matter.
fn is_negative(e: &Expr) -> bool {
    match e.view() {
        View::Int(n) => n < 0,
        View::Integer(n) => n.is_negative(),
        View::Rational(q) => q.is_negative(),
        View::Complex(z) => z.re().is_zero() && z.im().is_negative(),
        View::Ball(b) => ball_decimal(b, DISPLAY_PRECISION).starts_with('-'),
        View::ComplexBall(z) => {
            z.re().is_exact()
                && z.re().mid().is_zero()
                && ball_decimal(z.im(), DISPLAY_PRECISION).starts_with('-')
        }
        _ => false,
    }
}

fn is_rational(e: &Expr) -> bool {
    matches!(e.view(), View::Rational(_))
}

fn is_integer(e: &Expr) -> bool {
    matches!(e.view(), View::Int(_) | View::Integer(_))
}

/// The tasks that print the application `e`, in order.
fn expand<'a>(e: &'a Expr, out: &mut Vec<Task<'a>>) {
    match form(e) {
        Form::Plus(args) => {
            out.push(operand(&args[0], Prec::Product, false));
            for t in &args[1..] {
                if leading_minus(t) {
                    out.extend([Task::Text(" - "), Task::Negated(t)]);
                } else {
                    out.extend([Task::Text(" + "), operand(t, Prec::Product, true)]);
                }
            }
        }
        Form::Times(args) => times(args, false, out),
        Form::Power(b, n) => out.extend([
            operand(b, Prec::Atom, true),
            Task::Text("^"),
            operand(n, Prec::Power, true),
        ]),
        Form::List(items) => {
            out.push(Task::Text("{"));
            comma_separated(items, out);
            out.push(Task::Text("}"));
        }
        Form::Number | Form::Other => {
            let (head, args) = (e.head().expect("an application"), e.args());
            out.extend([operand(head, Prec::Atom, true), Task::Text("[")]);
            comma_separated(args, out);
            out.push(Task::Text("]"));
        }
    }
}

fn operand(e: &Expr, min: Prec, no_minus: bool) -> Task<'_> {
    Task::Expr { e, min, no_minus }
}

/// Arguments or list items: no parentheses needed between commas.
fn comma_separated<'a>(items: &'a [Expr], out: &mut Vec<Task<'a>>) {
    for (i, e) in items.iter().enumerate() {
        if i > 0 {
            out.push(Task::Text(", "));
        }
        out.push(operand(e, Prec::Lowest, false));
    }
}

/// A product of at least two factors (see [`product`]).
fn times<'a>(args: &'a [Expr], negated: bool, out: &mut Vec<Task<'a>>) {
    let factors: Vec<Factor<'a>> = args.iter().map(Factor::Plain).collect();
    product(&factors, negated, out);
}

/// Factors with a negative integer exponent are a denominator, and a
/// negative integer first in the numerator a sign (`-x`, `-2*x`). With
/// `negated`, the sign is left out.
fn product<'a>(factors: &[Factor<'a>], negated: bool, out: &mut Vec<Task<'a>>) {
    let mut numerator = Vec::new();
    // Each with the factor it came from.
    let mut denominator = Vec::new();
    for &f in factors {
        match f {
            Factor::Plain(e) => match inverse_power(e) {
                Some((b, -1)) => denominator.push((Factor::Plain(b), f)),
                Some((b, n)) => denominator.push((Factor::Power(b, n), f)),
                None => numerator.push(f),
            },
            Factor::Power(..) => numerator.push(f),
        }
    }
    // `n/d` for integers `n` and `d >= 0` is how a rational prints.
    if let ([Factor::Plain(n)], [(Factor::Plain(d), original)]) = (&numerator[..], &denominator[..])
    {
        if is_integer(n) && is_integer(d) && !is_negative(d) {
            numerator.push(*original);
            denominator.clear();
        }
    }
    let denominator: Vec<Factor<'a>> = denominator.into_iter().map(|(f, _)| f).collect();
    signed_product(&numerator, negated, out);
    match denominator[..] {
        [] => {}
        [d] => {
            out.push(Task::Text("/"));
            d.emit(out);
        }
        _ => {
            // Printed like the product it reads as, which may have a
            // denominator of its own: a task, so nesting uses no stack.
            out.extend([
                Task::Text("/("),
                Task::Product {
                    factors: denominator,
                    negated: false,
                },
                Task::Text(")"),
            ]);
        }
    }
}

/// A factor as printed in a product.
#[derive(Clone, Copy)]
enum Factor<'a> {
    Plain(&'a Expr),
    /// `b^-n`, for `Power[b, n]` with `n < -1` in a denominator.
    Power(&'a Expr, i64),
}

impl<'a> Factor<'a> {
    /// The number that the factor prints as, if it is a bare number.
    fn number(self) -> Option<&'a Expr> {
        match self {
            Factor::Plain(e) if e.is_number() => Some(e),
            _ => None,
        }
    }

    fn emit(self, out: &mut Vec<Task<'a>>) {
        match self {
            Factor::Plain(e) => out.push(operand(e, Prec::Power, true)),
            Factor::Power(b, n) => {
                out.extend([operand(b, Prec::Atom, true), Task::Text("^"), Task::Abs(n)]);
            }
        }
    }
}

/// `factors` joined by `*`, a negative number first as the sign, and `1`
/// for no factors. With `negated`, the sign is left out.
fn signed_product<'a>(factors: &[Factor<'a>], negated: bool, out: &mut Vec<Task<'a>>) {
    let mut items = Vec::new();
    let mut rest = factors;
    if let Some((&first, tail)) = factors.split_first() {
        if let Some(c) = first.number().filter(|c| is_negative(c)) {
            rest = tail;
            // `Times[-1, 3]` is `-1*3`: `-3` is the integer.
            let literal_next = tail
                .first()
                .and_then(|f| f.number())
                .is_some_and(|m| !is_negative(m));
            if c.as_i64() == Some(-1) && !literal_next {
                if !negated {
                    out.push(Task::Text("-"));
                }
            } else if negated {
                items.push(Task::Magnitude(c));
            } else {
                items.push(operand(c, Prec::Lowest, false));
            }
        }
    }
    let mut first = items.is_empty();
    if first && rest.is_empty() {
        out.push(Task::Text("1"));
    }
    out.extend(items);
    for &f in rest {
        match f {
            // A leading rational needs no parentheses: `1/2*x`.
            Factor::Plain(e) if first && is_rational(e) => {
                out.push(operand(e, Prec::Product, true));
            }
            _ => {
                if !first {
                    out.push(Task::Text("*"));
                }
                f.emit(out);
            }
        }
        first = false;
    }
}

/// `(b, n)` for `Power[b, n]` with an integer `n < 0`.
fn inverse_power(e: &Expr) -> Option<(&Expr, i64)> {
    match form(e) {
        Form::Power(b, n) => n.as_i64().filter(|&n| n < 0).map(|n| (b, n)),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
