//! Printing: full form (`Debug`, [`Expr::full_form`]) and infix (`Display`).
//!
//! Both printers run on an explicit stack, so expressions of any depth print
//! in constant native stack.

use std::fmt::{self, Write as _};

use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// An expression in full form, `Plus[x, Times[2, Power[y, 3]]]`: every
/// application as `head[args…]`, unambiguous. See [`Expr::full_form`].
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

/// Symbols by name, integers in decimal, strings quoted and escaped.
fn atom(e: &Expr, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match e.view() {
        View::Symbol(s) => f.write_str(s.name()),
        View::Int(x) => write!(f, "{x}"),
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
/// - A term with a negative coefficient is subtracted: `x - 2*y`, and
///   `Times[-1, x]` is `-x`.
/// - Factors with a negative integer exponent are a denominator:
///   `x/(y*z^2)`.
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
                    Form::Int(n) => expansion.push(Task::Abs(n)),
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
    /// `e`, parenthesised if it binds weaker than `min`, or if `no_minus`
    /// and it would start with a minus sign.
    Expr {
        e: &'a Expr,
        min: Prec,
        no_minus: bool,
    },
    /// `e` without its leading minus sign, which the caller has written. Only
    /// for an integer or a product with [`leading_minus`].
    Negated(&'a Expr),
    /// The product of `factors`, see [`product`].
    Product {
        factors: Vec<Factor<'a>>,
        negated: bool,
    },
}

/// What the infix printer sees in an expression.
enum Form<'a> {
    Int(i64),
    Plus(&'a [Expr]),
    Times(&'a [Expr]),
    Power(&'a Expr, &'a Expr),
    List(&'a [Expr]),
    /// Other atoms, and applications in full form.
    Other,
}

fn form(e: &Expr) -> Form<'_> {
    match e.view() {
        View::Int(n) => Form::Int(n),
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
        Form::Int(_) | Form::List(_) | Form::Other => Prec::Atom,
    }
}

/// Whether `e`, printed without parentheses, starts with a minus sign. Only
/// integers and products can: a sum with a negative first term is an operand
/// only where its precedence already parenthesises it. Precedence ignores the
/// sign; operands that must not start with one say so with `no_minus`.
fn leading_minus(e: &Expr) -> bool {
    match form(e) {
        Form::Int(n) => n < 0,
        Form::Times(args) => args
            .iter()
            .find(|f| inverse_power(f).is_none())
            .is_some_and(is_negative),
        _ => false,
    }
}

fn is_negative(e: &Expr) -> bool {
    e.as_i64().is_some_and(|n| n < 0)
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
        Form::Int(_) | Form::Other => {
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
    let mut denominator = Vec::new();
    for &f in factors {
        match f {
            Factor::Plain(e) => match inverse_power(e) {
                Some((b, -1)) => denominator.push(Factor::Plain(b)),
                Some((b, n)) => denominator.push(Factor::Power(b, n)),
                None => numerator.push(f),
            },
            Factor::Power(..) => numerator.push(f),
        }
    }
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
    /// The integer that the factor prints as, if it is a bare integer.
    fn int(self) -> Option<(&'a Expr, i64)> {
        match self {
            Factor::Plain(e) => e.as_i64().map(|n| (e, n)),
            Factor::Power(..) => None,
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

/// `factors` joined by `*`, a negative integer first as the sign, and `1`
/// for no factors. With `negated`, the sign is left out.
fn signed_product<'a>(factors: &[Factor<'a>], negated: bool, out: &mut Vec<Task<'a>>) {
    let mut items = Vec::new();
    let mut rest = factors;
    if let Some((&first, tail)) = factors.split_first() {
        if let Some((c, n)) = first.int().filter(|&(_, n)| n < 0) {
            rest = tail;
            // `Times[-1, 3]` is `-1*3`: `-3` is the integer.
            let literal_next = tail
                .first()
                .and_then(|f| f.int())
                .is_some_and(|(_, m)| m >= 0);
            if n == -1 && !literal_next {
                if !negated {
                    out.push(Task::Text("-"));
                }
            } else if negated {
                items.push(Task::Abs(n));
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
        if !first {
            out.push(Task::Text("*"));
        }
        first = false;
        f.emit(out);
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
