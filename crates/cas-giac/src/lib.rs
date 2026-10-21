//! libgiac, the engine of Xcas, as a test oracle (CAS-DOC-8 §12, CAS-138).
//!
//! Xcas is the reference for this CAS's defaults, so tests compare with its
//! answers. The library is loaded at run time, so this crate builds without
//! it, and checks skip when it is missing ([`Giac::shared`] is `None`).
//!
//! - [`Giac`]: the library's two entry points. [`Giac::eval`] is
//!   `caseval`, about 10 ms a call whatever it computes. [`Giac::eval_fast`]
//!   is `nws_caseval`, the NumWorks calculator's, from about 10 µs, with its
//!   calculator formatting undone. Calls take a global lock: giac keeps
//!   global state, and answers in a static buffer.
//! - [`ours`]: the value of a text in this CAS: parsed, lowered with the
//!   English names, its decimals made balls, evaluated in a new context.
//!   giac's answers read back this way, since the names follow Xcas
//!   (CAS-DOC-7 §5).
//! - [`Giac::value`] reads giac's answer back so; [`Giac::agrees`] compares
//!   it with this CAS's value of the same input, and [`Giac::answers`] with
//!   the value of an expected text. On a difference, the error shows
//!   giac's text and both values in the language.
//!
//! giac's state persists between calls: an assignment, even `x=1/2` (which
//! giac reads as one), changes the later answers. [`Giac::restart`]
//! forgets it.
//!
//! ```
//! use cas_giac::Giac;
//!
//! let Some(giac) = Giac::shared() else {
//!     return; // No libgiac: nothing to check.
//! };
//! assert_eq!(giac.eval("factor(x^2-1)").unwrap(), "(x-1)*(x+1)");
//! giac.agrees("diff(sin(x^2), x)").unwrap();
//! ```

mod format;
#[allow(
    unsafe_code,
    reason = "calling a C library loaded at run time: the one way to reach giac"
)]
mod library;

use std::fmt;
use std::sync::OnceLock;

use cas_core::{Context, Expr};
use cas_lang::{Names, decimal_balls, lower, parse, print};

use library::Library;

/// The environment variable naming the library to load instead of
/// `libgiac.so.0`.
pub const LIBRARY_VARIABLE: &str = "CAS_GIAC_LIBRARY";

/// The library loaded by default, from the usual paths.
pub const DEFAULT_LIBRARY: &str = "libgiac.so.0";

/// What went wrong with an oracle check.
#[derive(Clone, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum OracleError {
    /// The input has a NUL, or giac gave no answer.
    NoAnswer(String),
    /// A text this CAS cannot read or evaluate: the text, and why.
    Unreadable(String, String),
    /// The two sides differ.
    Differs {
        /// What giac was asked.
        input: String,
        /// giac's answer, as it printed it.
        giac: String,
        /// giac's answer read back and evaluated, in the language.
        giac_value: String,
        /// This CAS's value, in the language.
        ours: String,
    },
}

impl fmt::Display for OracleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OracleError::NoAnswer(input) => write!(f, "giac gave no answer to {input:?}"),
            OracleError::Unreadable(text, why) => write!(f, "cannot read {text:?}: {why}"),
            OracleError::Differs {
                input,
                giac,
                giac_value,
                ours,
            } => write!(
                f,
                "{input}: giac gives {giac} ({giac_value} here), this CAS {ours}"
            ),
        }
    }
}

impl std::error::Error for OracleError {}

/// A loaded libgiac.
#[derive(Clone, Copy, Debug)]
pub struct Giac {
    library: Library,
}

impl Giac {
    /// The library at `path`; a file name searches the usual paths.
    ///
    /// # Errors
    ///
    /// The loader's message when the library or its entry points are
    /// missing.
    pub fn load(path: &str) -> Result<Giac, String> {
        Library::load(path).map(|library| Giac { library })
    }

    /// The library named by [`LIBRARY_VARIABLE`], or [`DEFAULT_LIBRARY`],
    /// loaded once for the process; `None` when it cannot be loaded, so that
    /// checks skip.
    #[must_use]
    pub fn shared() -> Option<&'static Giac> {
        static SHARED: OnceLock<Option<Giac>> = OnceLock::new();
        SHARED
            .get_or_init(|| {
                let path =
                    std::env::var(LIBRARY_VARIABLE).unwrap_or_else(|_| DEFAULT_LIBRARY.to_owned());
                Giac::load(&path).ok()
            })
            .as_ref()
    }

    /// giac's answer to `input`, through `caseval`, as it prints it.
    ///
    /// # Errors
    ///
    /// [`OracleError::NoAnswer`] for an input with a NUL, or no answer.
    pub fn eval(&self, input: &str) -> Result<String, OracleError> {
        self.library
            .call(input, false)
            .ok_or_else(|| OracleError::NoAnswer(input.to_owned()))
    }

    /// giac's answer to `input`, through `nws_caseval`, its calculator
    /// formatting undone: `√(…)` is `sqrt(…)`, the appended approximation
    /// goes, matrix rows get their commas. A one-row matrix reads back as a
    /// vector, and an equation with a decimal right side loses it.
    ///
    /// # Errors
    ///
    /// As [`Giac::eval`].
    pub fn eval_fast(&self, input: &str) -> Result<String, OracleError> {
        self.library
            .call(input, true)
            .map(|answer| format::undo_nws(&answer))
            .ok_or_else(|| OracleError::NoAnswer(input.to_owned()))
    }

    /// giac's `restart`: forgets assignments, so later answers are those of
    /// a new session.
    pub fn restart(&self) {
        let _ = self.library.call("restart", false);
    }

    /// giac's answer to `input`, read back in the language and evaluated
    /// ([`ours`]).
    ///
    /// # Errors
    ///
    /// [`OracleError::NoAnswer`], or [`OracleError::Unreadable`] when this
    /// CAS cannot read the answer.
    pub fn value(&self, input: &str) -> Result<Expr, OracleError> {
        ours(&format::for_language(&self.eval(input)?))
    }

    /// Whether giac's answer to `input` is, once read back and evaluated,
    /// this CAS's value of the same `input`: for builtins with Xcas's names.
    ///
    /// # Errors
    ///
    /// [`OracleError::Differs`] showing both, or the errors of reading
    /// either side.
    pub fn agrees(&self, input: &str) -> Result<(), OracleError> {
        let expected = ours(input)?;
        self.compare(input, &expected)
    }

    /// Whether giac's answer to `input` is, once read back and evaluated,
    /// the value of `expected` in this CAS.
    ///
    /// # Errors
    ///
    /// As [`Giac::agrees`].
    pub fn answers(&self, input: &str, expected: &str) -> Result<(), OracleError> {
        let expected = ours(expected)?;
        self.compare(input, &expected)
    }

    /// giac's answer to `input` against `expected`: the same expression, or
    /// equal by this CAS's `==`, which takes accurate overlapping balls as
    /// equal (CAS-DOC-5 §9.1, D3), as decimals from the two sides are.
    fn compare(&self, input: &str, expected: &Expr) -> Result<(), OracleError> {
        let answer = self.eval(input)?;
        let value = ours(&format::for_language(&answer))?;
        let equal = Expr::apply(
            &cas_core::Symbol::EQUAL.into(),
            &[value.clone(), expected.clone()],
        );
        if value == *expected
            || Context::new().eval(&equal).ok().and_then(|e| e.as_symbol())
                == Some(cas_core::Symbol::TRUE)
        {
            return Ok(());
        }
        Err(OracleError::Differs {
            input: input.to_owned(),
            giac: answer,
            giac_value: written(&value),
            ours: written(expected),
        })
    }
}

/// `e` in the language, or in the evaluator's notation where the language
/// cannot write it.
fn written(e: &Expr) -> String {
    print(e, Names::english()).unwrap_or_else(|_| e.to_string())
}

/// The value of `text` in this CAS: one statement, parsed and lowered with
/// the English names, its decimals balls at 53 bits, evaluated in a new
/// context.
///
/// # Errors
///
/// [`OracleError::Unreadable`] for a text with syntax or lowering errors,
/// not one statement, or whose evaluation fails.
pub fn ours(text: &str) -> Result<Expr, OracleError> {
    let unreadable = |why: String| OracleError::Unreadable(text.to_owned(), why);
    let parsed = parse(text);
    if let Some(e) = parsed.errors().first() {
        return Err(unreadable(format!("{e:?}")));
    }
    let lowered = lower(&parsed, Names::english());
    if lowered.has_errors() {
        return Err(unreadable(format!("{:?}", lowered.diagnostics)));
    }
    let [statement] = &lowered.statements[..] else {
        return Err(unreadable("not one statement".to_owned()));
    };
    let e = decimal_balls(&statement.expr, cas_core::DEFAULT_WORKING_PRECISION);
    Context::new()
        .eval(&e)
        .map_err(|error| unreadable(error.to_string()))
}
