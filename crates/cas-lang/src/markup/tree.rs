//! The tree the LaTeX printer works on: giac's `gen`, in the shapes giac's
//! evaluator leaves (markup.cc reads `g.type`, `g._SYMBptr->sommet` and
//! `feuille`). Built by [`lower`](super::lower) from our expressions, and by
//! the tests from giac's printed answers.

use cas_core::Integer;

/// A giac `gen`, as far as the LaTeX printer looks at it.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Gen {
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
    /// `_VECT` with its subtype.
    Vect(Vec<Gen>, Sub),
    /// `_SYMB`: the operator's name (giac's `sommet`) and its argument, a
    /// `Vect` of subtype `Seq` when it has several.
    Symb(String, Box<Gen>),
}

/// The subtypes of a `_VECT` the printer tells apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Sub {
    /// A list, `[…]` (subtype 0).
    List,
    /// A sequence, `(…)`, bare at the top level (`_SEQ__VECT`).
    Seq,
    /// A set, `{…}` (`_SET__VECT`).
    Set,
}

impl Gen {
    pub(super) fn int(n: i64) -> Gen {
        Gen::Int(Integer::from(n))
    }

    /// The imaginary unit.
    pub(super) fn i() -> Gen {
        Gen::Cplx(Box::new(Gen::int(0)), Box::new(Gen::int(1)))
    }

    pub(super) fn idnt(name: &str) -> Gen {
        Gen::Idnt(name.to_owned())
    }

    /// `op(arg)`.
    pub(super) fn symb(op: &str, arg: Gen) -> Gen {
        Gen::Symb(op.to_owned(), Box::new(arg))
    }

    /// `op(args…)`: the argument a sequence.
    pub(super) fn symb_n(op: &str, args: Vec<Gen>) -> Gen {
        Gen::symb(op, Gen::Vect(args, Sub::Seq))
    }

    /// The arguments of a `Symb`: its sequence's items, or the one argument.
    pub(super) fn args(&self) -> &[Gen] {
        match self {
            Gen::Symb(_, arg) => match &**arg {
                Gen::Vect(items, _) => items,
                a => std::slice::from_ref(a),
            },
            _ => &[],
        }
    }

    /// Whether this is a `Symb` of one of `ops`.
    pub(super) fn is_sym(&self, ops: &[&str]) -> bool {
        matches!(self, Gen::Symb(op, _) if ops.contains(&op.as_str()))
    }

    /// Whether this is a number equal to `p/q` (giac's `is_one`,
    /// `is_zero`, `is_minus_one` and `2*e == 1`).
    pub(super) fn equals(&self, p: i64, q: i64) -> bool {
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

    pub(super) fn is_zero(&self) -> bool {
        self.equals(0, 1)
    }

    pub(super) fn is_one(&self) -> bool {
        self.equals(1, 1)
    }

    pub(super) fn is_minus_one(&self) -> bool {
        self.equals(-1, 1)
    }
}
