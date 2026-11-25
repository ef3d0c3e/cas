//! The tree the LaTeX and `MathML` printers work on: giac's `gen`, in the shapes giac's
//! evaluator leaves (markup.cc reads `g.type`, `g._SYMBptr->sommet` and
//! `feuille`). Built by [`lower`](super::lower) from our expressions, and by
//! the tests from giac's printed answers.

use cas_core::Integer;

/// A giac `gen`, as far as the printers look at it.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Gen {
    /// `_INT_`/`_ZINT`.
    Int(Integer),
    /// A boolean (`_INT_` of subtype `_INT_BOOLEAN`).
    Bool(bool),
    /// `_FRAC`: an integer, or a gaussian integer (`Cplx`), over an integer.
    Frac(Box<Gen>, Box<Gen>),
    /// `_DOUBLE_`: `text` is giac's print of the magnitude (12 significant
    /// digits, `1e-05`, `3.0`), `neg` its sign.
    Dbl { text: String, neg: bool },
    /// `_CPLX`: real and imaginary parts, numbers.
    Cplx(Box<Gen>, Box<Gen>),
    /// `_IDNT`, by name.
    Idnt(String),
    /// `_STRNG`.
    Str(String),
    /// `_FUNC`: a builtin by its giac name, not applied (`Gamma`).
    Func(String),
    /// `_VECT` with its subtype.
    Vect(Vec<Gen>, Sub),
    /// `_SYMB`: the operator's name (giac's `sommet`) and its argument, a
    /// `Vect` of subtype `Seq` when it has several.
    Symb(String, Box<Gen>),
}

/// The subtypes of a `_VECT` the printer tells apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Sub {
    /// A list, `[…]` (subtype 0).
    List,
    /// A sequence, `(…)`, bare at the top level (`_SEQ__VECT`).
    Seq,
    /// A set, `{…}` (`_SET__VECT`).
    Set,
    /// A dense polynomial's coefficients, `poly1[…]` (`_POLY1__VECT`): a
    /// list to gen2markup, `poly1[` to mathml.cc (`begin_VECT_string`).
    /// Only giac's shapes in the tests make one; our lowering never does.
    #[cfg_attr(not(test), allow(dead_code))]
    Poly1,
}

impl Gen {
    pub(crate) fn int(n: i64) -> Gen {
        Gen::Int(Integer::from(n))
    }

    /// The imaginary unit.
    pub(crate) fn i() -> Gen {
        Gen::Cplx(Box::new(Gen::int(0)), Box::new(Gen::int(1)))
    }

    pub(crate) fn idnt(name: &str) -> Gen {
        Gen::Idnt(name.to_owned())
    }

    /// `op(arg)`.
    pub(crate) fn symb(op: &str, arg: Gen) -> Gen {
        Gen::Symb(op.to_owned(), Box::new(arg))
    }

    /// `op(args…)`: the argument a sequence.
    pub(crate) fn symb_n(op: &str, args: Vec<Gen>) -> Gen {
        Gen::symb(op, Gen::Vect(args, Sub::Seq))
    }

    /// The arguments of a `Symb`: its sequence's items, or the one argument.
    pub(crate) fn args(&self) -> &[Gen] {
        match self {
            Gen::Symb(_, arg) => match &**arg {
                Gen::Vect(items, _) => items,
                a => std::slice::from_ref(a),
            },
            _ => &[],
        }
    }

    /// Whether this is a `Symb` of one of `ops`.
    pub(crate) fn is_sym(&self, ops: &[&str]) -> bool {
        matches!(self, Gen::Symb(op, _) if ops.contains(&op.as_str()))
    }

    /// Whether this is a number equal to `p/q` (giac's `is_one`,
    /// `is_zero`, `is_minus_one` and `2*e == 1`).
    pub(crate) fn equals(&self, p: i64, q: i64) -> bool {
        match self {
            Gen::Int(n) => q == 1 && *n == Integer::from(p),
            Gen::Frac(a, b) => match (&**a, &**b) {
                (Gen::Int(a), Gen::Int(b)) => a * &Integer::from(q) == &Integer::from(p) * b,
                _ => false,
            },
            #[allow(clippy::cast_precision_loss, clippy::float_cmp)]
            Gen::Dbl { text, neg } => text.parse::<f64>().is_ok_and(|v| {
                let v = if *neg { -v } else { v };
                v == p as f64 / q as f64
            }),
            _ => false,
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.equals(0, 1)
    }

    pub(crate) fn is_one(&self) -> bool {
        self.equals(1, 1)
    }

    pub(crate) fn is_minus_one(&self) -> bool {
        self.equals(-1, 1)
    }
}
