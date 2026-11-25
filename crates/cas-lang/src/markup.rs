//! LaTeX as giac writes it (CAS-615): a port of giac's `gen2markup`
//! (giac 2.0.0.21, `src/markup.cc`), the converter behind Xcas's `latex()`,
//! LaTeX arms only, so `latex(e)` gives giac's text for giac's answer.
//!
//! - [`lower`]: our expression as the tree giac's evaluator would leave
//!   ([`Gen`](tree::Gen)): `x/2` is `prod(x, inv(2))`, `x^(1/2)` `sqrt(x)`,
//!   `Times[-1, y]` `neg(y)`. The shape decides the text, so this is where
//!   our forms meet giac's.
//! - [`engine`]: `gen2markup` over that tree, one function per giac branch,
//!   each citing its markup.cc lines (`:1869` is line 1869): atoms
//!   (:1308-1542), vectors (:1658-1773), `neg` (:1777), `plus` (:1793),
//!   `prod` (:1869-2032), `inv` (:2082), `division` (:2097), `exp` (:2148),
//!   `sqrt` (:2191), `pow` (:2224-2280), `abs`/`floor`/`ceil` (:2281),
//!   integrals, sums and products (:2375-2487), `limit` (:2488-2548),
//!   derivatives (:2560-2727), relations (:2728, :2770), `interval`
//!   (:2942), `factorial` (:3317), `binomial` (:3334), `at` (:3434), and
//!   the general application (:3532-3611). A [`MarkupBlock`](engine::MarkupBlock)
//!   carries a text with its priority, kind bits and a pending minus;
//!   [`Flags`](engine::Flags) select LaTeX, the top level, powers of
//!   elementary functions on their names (`ELEMPOW`) and factors. `MathML`
//!   (T5) adds its fields and flags to the same engine.
//! - [`ident`]: identifiers (`x_{1}`, `\alpha`, `\mathrm{abc}`), function
//!   names (`\sin `, `\operatorname{arsinh}`) and strings.
//!
//! The text keeps giac's quirks: trailing spaces (`\pi `, `\sin `,
//! `\cdot `), `\frac1{…}` from `inv` only, `abs` dropping its argument's
//! minus. The prototype that fixed each step against giac's `latex()` on 301
//! inputs is 1.6/T4.c1's (`.run/1.6/T4/proto/`); its table is the engine's
//! test data (`markup/giac_latex.tsv`).

mod engine;
mod ident;
mod lower;
mod mml;
mod tree;

use cas_core::Expr;

use crate::names::Names;
use crate::print::{Balls, MAX_PRINT_DEPTH, PrintError, PrintErrorKind, print_capped};

/// `e` in LaTeX as giac's `latex()` writes the same answer (`export_latex`,
/// markup.cc :3629): `x/2` is `\frac{x}{2}`, `cos(x)^2` `\cos ^{2}x`. Names
/// of builtins without a giac form are spelled from `names`, balls written
/// as `balls` shows them, in giac's 12 digits.
///
/// # Errors
///
/// [`PrintErrorKind::TooDeep`] for `e` nested deeper than
/// [`MAX_PRINT_DEPTH`], refused before any work;
/// [`PrintErrorKind::TooLong`] as soon as a part's text is longer than
/// `max_len` bytes; [`PrintErrorKind::Ball`] for a decimal whose exponent
/// is beyond `i64`.
pub fn latex(e: &Expr, names: &Names, balls: Balls, max_len: usize) -> Result<String, PrintError> {
    let g = tree(e, names, balls)?;
    engine::Engine::new(max_len)
        .export_latex(&g)
        .map_err(|engine::TooLong| too_long(e))
}

/// The forms of giac's `export_mathml` (markup.cc :3747-3776).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MathmlKind {
    /// `export_mathml(e)` (:3765): presentation in `<semantics>`, its tags
    /// referring (`xref`) to the content's ids, the content in
    /// `<annotation-xml encoding='MathML-Content'>`, and the input form in
    /// `<annotation encoding='giac'>`.
    Default,
    /// `export_mathml(e, display)` (:3755): presentation only.
    Display,
    /// `export_mathml(e, content)` (:3747): content only.
    Content,
}

/// `e` in MathML as giac's `export_mathml` writes the same answer, in the
/// form `kind`; names and balls as [`latex`] has them. The annotation of
/// [`MathmlKind::Default`] is our input form of `e` (giac's is its own),
/// escaped (`str_to_mml`).
///
/// # Errors
///
/// As [`latex`]'s, and the input form's for the annotation.
pub fn export_mathml(
    e: &Expr,
    names: &Names,
    balls: Balls,
    max_len: usize,
    kind: MathmlKind,
) -> Result<String, PrintError> {
    let g = tree(e, names, balls)?;
    let engine = engine::Engine::new(max_len);
    match kind {
        MathmlKind::Default => {
            let printed = print_capped(e, names, balls, max_len)?;
            engine.export_mathml(&g, &printed)
        }
        MathmlKind::Display => engine.export_mathml_presentation(&g),
        MathmlKind::Content => engine.export_mathml_content(&g),
    }
    .map_err(|engine::TooLong| too_long(e))
}

/// `e` as giac's tree, refused when deeper than [`MAX_PRINT_DEPTH`].
fn tree(e: &Expr, names: &Names, balls: Balls) -> Result<tree::Gen, PrintError> {
    if usize::from(e.depth()) > MAX_PRINT_DEPTH {
        return Err(PrintError {
            part: e.clone(),
            kind: PrintErrorKind::TooDeep,
        });
    }
    lower::Lowering { names, balls }.tree(e)
}

fn too_long(e: &Expr) -> PrintError {
    PrintError {
        part: e.clone(),
        kind: PrintErrorKind::TooLong,
    }
}

#[cfg(test)]
mod tests;
