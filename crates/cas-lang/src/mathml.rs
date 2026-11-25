//! `MathML` as giac's `mathml` writes it (CAS-616): a port of giac's own
//! `MathML` printer (giac 2.0.0.21, `src/mathml.cc`, the live Linux branch
//! :63-2032), older than `export_mathml`'s (`markup.cc`, [`crate::markup`])
//! and written apart from it: presentation only, no ids, and an XHTML
//! document around the formula.
//!
//! - [`mathml`]: `_mathml` (:1987) with one argument. The text is
//!   `gen2mathmlfull` (:1983): the XHTML preamble (:69), the formula in
//!   `<math mode="display" …>` with a `<br/>` (`ingen2mathml` :1980), the
//!   end (:70) and a newline. A list `[e, n]`, `n` a machine integer, is
//!   the formula alone (:1999-2000), `<math display="inline" …>` when `n` is
//!   not 0 (:1979).
//! - `convert`: `gen2mathml` (:1879-1975) over the tree giac's evaluator
//!   would leave ([`Gen`], from [`crate::markup`]'s lowering): numbers
//!   `<mn>`, complex numbers by parts (:1899), identifiers (`idnt2mathml`
//!   :1870, greek letters :1808), strings in a one-cell `<mtable>` (:111),
//!   vectors and matrices (:200, :221), and `symbolic2mathml` (:1571-1773):
//!   `of`/`at` (:1576), relations, `re`, `im`, `exp`, `sqrt`, integrals,
//!   derivatives, `:=` and limits (`mathml_print` :526-598, :403-523),
//!   `and`/`or` (:1592), one-argument heads, `neg` (:1637), `inv` (:1645),
//!   `plus` (:1674), `prod` with fractions and `&times;` (:1698, :274-364),
//!   `pow` with `<msqrt>` and the bracketed base (:1723), and the general
//!   `<mi>f</mi><mrow><mo>(</mo>…<mo>)</mo></mrow>` (:1766); `sum` (:600),
//!   `abs` (:607) first.
//! - `flat`: giac's `aplatir_fois_plus` (gen.cc :12505), run first, as
//!   giac runs it on each symbolic (:1929): nested sums flattened, a
//!   product of two multiplied out (`sym_mult` gen.cc :7228), so `2*inv(4)`
//!   is `1/2`.
//! - `print`: giac's print of a `gen` (`symbolic.cc` `add_print`, unary.cc
//!   `printsommetasoperator`) for the texts giac takes printed: the
//!   `<mtext>` of a head without a `MathML` form (`1 .. 2`), integral
//!   bounds, `-n` in a sum.
//!
//! The text keeps giac's quirks: a space in `<msup><mrow> <mi>x</mi>`, the
//! bracketed base written twice (`<mfenced …>(…)</mfenced>`), raw commas
//! between the general application's arguments, unescaped `&` in
//! `<mtext>`. Our `x^(1/2)` (lowered as `sqrt(x)`) is giac's evaluated
//! `pow(x, 1/2)` here, since this printer writes the two apart. The
//! prototype that fixed each step against giac's `mathml()` on the 301
//! inputs of `markup`'s tables is 1.6/T5.c3's (`.run/1.6/T5/c3/`); its
//! table is the test data (`markup/giac_mathml_doc.tsv`).
//!
//! Not ported: giac's evaluation inside the printer (`eval` of a limit's,
//! sum's or integral's body, :468, :499, :604, and `needs_times` :256 are
//! the identity here), `rootof` values (:1625), modular vectors, spreadsheets
//! and graphics (`pnt`, SVG), `sym_mult` with doubles, and giac's size
//! guard (`max_prettyprint_equation`, :1996): ours are the cap and depth
//! checks [`crate::markup`] uses.

mod convert;
mod flat;
mod print;

use cas_core::Expr;

use crate::markup::tree::{Gen, Sub};
use crate::names::Names;
use crate::print::{Balls, PrintError, PrintErrorKind};

/// The XHTML around the formula (:69).
const PREAMBLE: &str = "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1 plus MathML 2.0//EN\"\n\
\"http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd\" [\n\
<!ENTITY mathml \"http://www.w3.org/1998/Math/MathML\">\n\
]>\n<html xmlns=\"http://www.w3.org/1999/xhtml\">\n<body>\n";

/// After it (:70), and `gen2mathmlfull`'s newline (:1984).
const END: &str = "\n</body>\n</html>\n";

/// `e` as giac's `mathml(e)` writes the same answer: an XHTML document
/// with the formula in presentation `MathML`. Names and balls as
/// [`latex`](crate::latex) has them.
///
/// # Errors
///
/// [`PrintErrorKind::TooDeep`] for `e` nested deeper than
/// [`MAX_PRINT_DEPTH`](crate::print::MAX_PRINT_DEPTH), refused before any work;
/// [`PrintErrorKind::TooLong`] as soon as a part's text is longer than
/// `max_len` bytes; [`PrintErrorKind::Ball`] for a decimal whose exponent
/// is beyond `i64`.
pub fn mathml(e: &Expr, names: &Names, balls: Balls, max_len: usize) -> Result<String, PrintError> {
    let g = crate::markup::tree(e, names, balls)?;
    document(&g, max_len).map_err(|convert::TooLong| PrintError {
        part: e.clone(),
        kind: PrintErrorKind::TooLong,
    })
}

/// `_mathml` (:1987-2006) on giac's tree.
fn document(g: &Gen, max_len: usize) -> Result<String, convert::TooLong> {
    let g = flat::aplatir(&giac_shape(g));
    let writer = convert::Writer::new(max_len);
    let text = match &g {
        // :1999-2005: `[e, n]` and `(a, b, …, n)`, the formula alone.
        Gen::Vect(items, sub) => match (items.as_slice(), machine_int(items.last())) {
            ([front, _], Some(n)) => inline(&writer.node(front)?, n),
            ([.., _], Some(n)) if *sub == Sub::Seq && items.len() > 2 => {
                let head = Gen::Vect(items[..items.len() - 1].to_vec(), Sub::Seq);
                inline(&writer.node(&head)?, n)
            }
            _ => full(&writer.node(&g)?),
        },
        _ => full(&writer.node(&g)?),
    };
    writer.checked(text)
}

/// A `_INT_`'s value: an integer giac holds in 32 bits, or a boolean.
fn machine_int(g: Option<&Gen>) -> Option<i64> {
    match g? {
        Gen::Int(n) => n.to_i64().filter(|v| i32::try_from(*v).is_ok()),
        Gen::Bool(b) => Some(i64::from(*b)),
        _ => None,
    }
}

/// `ingen2mathml(g, n)` (:1978-1981).
fn inline(formula: &str, n: i64) -> String {
    if n == 0 {
        format!(
            "\n<math mode=\"display\" xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\n{formula}\n\n</math><br/>\n"
        )
    } else {
        format!(
            "\n<math display=\"inline\" xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\n{formula}\n\n</math>\n"
        )
    }
}

/// `gen2mathmlfull` (:1983).
fn full(formula: &str) -> String {
    [PREAMBLE, &inline(formula, 0), END].concat()
}

/// Our lowering writes `x^(1/2)` as `sqrt(x)`, as `markup.cc` takes both;
/// giac's evaluator leaves `pow(x, 1/2)`, which this printer writes apart
/// from a `sqrt` (:1725, :407).
fn giac_shape(g: &Gen) -> Gen {
    match g {
        Gen::Vect(items, sub) => Gen::Vect(items.iter().map(giac_shape).collect(), *sub),
        Gen::Symb(op, arg) if op == "sqrt" => Gen::symb_n(
            "pow",
            vec![
                giac_shape(arg),
                Gen::Frac(Box::new(Gen::int(1)), Box::new(Gen::int(2))),
            ],
        ),
        Gen::Symb(op, arg) => Gen::symb(op, giac_shape(arg)),
        _ => g.clone(),
    }
}
