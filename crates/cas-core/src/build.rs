//! Building expressions from Rust: operators, calls and builtin functions.
//!
//! Everything here is context-free (CAS-DOC-3, D1 and D3). Builtin heads get
//! the cheap structural normalisation their fixed attributes call for: `Flat`
//! heads are flattened and `Orderless` heads sorted into the canonical order,
//! the evaluator's, which keeps the lists of a `Times` in their order.
//! Nothing is evaluated: `x + x` stays `Plus[x, x]` and `1 + 2` stays
//! `Plus[1, 2]`.
//!
//! ```
//! use cas_core::build::{cos, sin, symbols};
//!
//! let [x, y] = symbols(["x", "y"]);
//! let e = &x + 2 * y.pow(3);
//! assert_eq!(format!("{e:?}"), "Plus[x, Times[2, Power[y, 3]]]");
//!
//! // Flattened and sorted: `+` builds one `Plus`.
//! let e = sin(&x).pow(2) + (&y + cos(&x).pow(2)) + 1;
//! assert_eq!(format!("{e:?}"), "Plus[1, y, Power[Cos[x], 2], Power[Sin[x], 2]]");
//!
//! // `-` and `/` in terms of `Plus`, `Times` and `Power`.
//! assert_eq!(format!("{:?}", &x - &y), "Plus[x, Times[-1, y]]");
//! assert_eq!(format!("{:?}", &x / &y), "Times[x, Power[y, -1]]");
//! assert_eq!(format!("{:?}", x - 3), "Plus[-3, x]");
//! ```

use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use cas_num::{Complex, Integer};

use crate::attributes::Attributes;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

impl Expr {
    /// `head[args…]`, normalised if `head` is a builtin: the arguments of a
    /// `Flat` head are flattened (`Plus[a, Plus[b, c]]` → `Plus[a, b, c]`)
    /// and those of an `Orderless` head sorted into the canonical order.
    /// Other heads are applied as given, like [`Expr::apply`].
    ///
    /// ```
    /// use cas_core::{Expr, Symbol};
    ///
    /// let (f, x, y) = (Symbol::new("f"), Symbol::new("x"), Symbol::new("y"));
    /// assert_eq!(format!("{:?}", Expr::call(f, [y, x])), "f[y, x]");
    /// assert_eq!(format!("{:?}", Expr::call(Symbol::TIMES, [y, x])), "Times[x, y]");
    /// ```
    #[must_use]
    pub fn call<A: Into<Expr>>(head: impl Into<Expr>, args: impl IntoIterator<Item = A>) -> Expr {
        normalized(&head.into(), args.into_iter().map(Into::into).collect())
    }

    /// The symbol named `name`, as an expression.
    #[must_use]
    pub fn symbol(name: &str) -> Expr {
        Symbol::new(name).into()
    }

    /// `Power[self, exponent]`. Not simplified: `x.pow(1)` is `Power[x, 1]`.
    #[must_use]
    pub fn pow(&self, exponent: impl Into<Expr>) -> Expr {
        Expr::apply(&Symbol::POWER.into(), &[self.clone(), exponent.into()])
    }
}

impl Symbol {
    /// `self[args…]`, normalised like [`Expr::call`].
    #[must_use]
    pub fn call<A: Into<Expr>>(self, args: impl IntoIterator<Item = A>) -> Expr {
        Expr::call(self, args)
    }

    /// `Power[self, exponent]`.
    #[must_use]
    pub fn pow(self, exponent: impl Into<Expr>) -> Expr {
        Expr::from(self).pow(exponent)
    }
}

/// The normalisation of [`Expr::call`].
fn normalized(head: &Expr, mut args: Vec<Expr>) -> Expr {
    let attributes = head
        .as_symbol()
        .and_then(Symbol::builtin_attributes)
        .unwrap_or(Attributes::NONE);
    if attributes.contains(Attributes::FLAT) && args.iter().any(|a| a.head() == Some(head)) {
        args = flatten(head, &args);
    }
    // The evaluator's canonical order: in `Times`, lists after the other
    // factors in their own order, as `*` does not commute on them (R13); a
    // plain sort made `Q R` into `R Q` (CAS-344).
    if attributes.contains(Attributes::ORDERLESS)
        && let Some(sorted) = crate::eval::sort_arguments(head, &args)
    {
        args = sorted;
    }
    Expr::apply(head, &args)
}

/// Splices nested applications of `head`, at any depth, into one argument
/// list. Iterative, so deep nesting uses constant stack.
pub(crate) fn flatten(head: &Expr, args: &[Expr]) -> Vec<Expr> {
    let mut out = Vec::with_capacity(args.len());
    let mut stack = vec![args.iter()];
    while let Some(level) = stack.last_mut() {
        match level.next() {
            Some(a) if a.head() == Some(head) => stack.push(a.args().iter()),
            Some(a) => out.push(a.clone()),
            None => {
                stack.pop();
            }
        }
    }
    out
}

fn plus(a: Expr, b: Expr) -> Expr {
    normalized(&Symbol::PLUS.into(), vec![a, b])
}

fn times(a: Expr, b: Expr) -> Expr {
    normalized(&Symbol::TIMES.into(), vec![a, b])
}

/// `Times[-1, e]`; the negated literal for a number atom.
fn negate(e: Expr) -> Expr {
    match e.view() {
        View::Int(n) => match n.checked_neg() {
            Some(n) => n.into(),
            None => Expr::from(-Integer::from(n)),
        },
        View::Integer(n) => Expr::from(-n),
        View::Rational(q) => Expr::from(-q),
        View::Complex(z) => Expr::from(-z),
        // Exact for balls too.
        View::Ball(b) => Expr::from(b.neg()),
        View::ComplexBall(z) => Expr::from(z.neg()),
        _ => times((-1).into(), e),
    }
}

impl From<&Expr> for Expr {
    fn from(e: &Expr) -> Expr {
        e.clone()
    }
}

/// The operators for `Expr`, `&Expr` and `Symbol` on the left, and `i64` on
/// the left of those three.
macro_rules! operators {
    ($($op:ident $method:ident $assign:ident $assign_method:ident => $build:expr;)*) => {$(
        impl<T: Into<Expr>> $op<T> for Expr {
            type Output = Expr;

            fn $method(self, rhs: T) -> Expr {
                $build(self, rhs.into())
            }
        }

        impl<T: Into<Expr>> $op<T> for &Expr {
            type Output = Expr;

            fn $method(self, rhs: T) -> Expr {
                $build(self.clone(), rhs.into())
            }
        }

        impl<T: Into<Expr>> $op<T> for Symbol {
            type Output = Expr;

            fn $method(self, rhs: T) -> Expr {
                $build(self.into(), rhs.into())
            }
        }

        impl $op<Expr> for i64 {
            type Output = Expr;

            fn $method(self, rhs: Expr) -> Expr {
                $build(self.into(), rhs)
            }
        }

        impl $op<&Expr> for i64 {
            type Output = Expr;

            fn $method(self, rhs: &Expr) -> Expr {
                $build(self.into(), rhs.clone())
            }
        }

        impl $op<Symbol> for i64 {
            type Output = Expr;

            fn $method(self, rhs: Symbol) -> Expr {
                $build(self.into(), rhs.into())
            }
        }

        impl<T: Into<Expr>> $assign<T> for Expr {
            fn $assign_method(&mut self, rhs: T) {
                *self = $build(std::mem::replace(self, Expr::from(0)), rhs.into());
            }
        }
    )*};
}

operators! {
    Add add AddAssign add_assign => plus;
    Sub sub SubAssign sub_assign => |a, b| plus(a, negate(b));
    Mul mul MulAssign mul_assign => times;
    Div div DivAssign div_assign => |a, b: Expr| times(a, b.pow(-1));
}

/// `Times[-1, e]`, or a negative literal.
impl Neg for Expr {
    type Output = Expr;

    fn neg(self) -> Expr {
        negate(self)
    }
}

impl Neg for &Expr {
    type Output = Expr;

    fn neg(self) -> Expr {
        negate(self.clone())
    }
}

impl Neg for Symbol {
    type Output = Expr;

    fn neg(self) -> Expr {
        negate(self.into())
    }
}

/// One `Plus` of all the terms; `0` if there are none, the term itself if
/// there is one.
impl Sum for Expr {
    fn sum<I: Iterator<Item = Expr>>(iter: I) -> Expr {
        fold(Symbol::PLUS, 0, iter.collect())
    }
}

impl<'a> Sum<&'a Expr> for Expr {
    fn sum<I: Iterator<Item = &'a Expr>>(iter: I) -> Expr {
        iter.cloned().sum()
    }
}

/// One `Times` of all the factors; `1` if there are none, the factor itself
/// if there is one.
impl Product for Expr {
    fn product<I: Iterator<Item = Expr>>(iter: I) -> Expr {
        fold(Symbol::TIMES, 1, iter.collect())
    }
}

impl<'a> Product<&'a Expr> for Expr {
    fn product<I: Iterator<Item = &'a Expr>>(iter: I) -> Expr {
        iter.cloned().product()
    }
}

fn fold(head: Symbol, empty: i64, mut args: Vec<Expr>) -> Expr {
    match args.len() {
        0 => empty.into(),
        1 => args.pop().expect("one argument"),
        _ => normalized(&head.into(), args),
    }
}

/// `n` symbols, by name: `let [x, y] = symbols(["x", "y"]);`.
#[must_use]
pub fn symbols<const N: usize>(names: [&str; N]) -> [Expr; N] {
    names.map(Expr::symbol)
}

/// The imaginary unit: the complex number `0 + 1·I`, not the symbol
/// [`Symbol::I`], which evaluates to it.
#[must_use]
pub fn i() -> Expr {
    Expr::from(Complex::I)
}

/// `List[items…]`.
#[must_use]
pub fn list<A: Into<Expr>>(items: impl IntoIterator<Item = A>) -> Expr {
    Expr::call(Symbol::LIST, items)
}

macro_rules! functions {
    ($($name:ident $symbol:ident;)*) => {$(
        #[doc = concat!("[`Symbol::", stringify!($symbol), "`]`[x]`.")]
        #[must_use]
        pub fn $name(x: impl Into<Expr>) -> Expr {
            Expr::call(Symbol::$symbol, [x])
        }
    )*};
}

functions! {
    exp EXP;
    log LOG;
    sin SIN;
    cos COS;
    tan TAN;
}

#[cfg(test)]
mod tests;
