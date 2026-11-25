//! The language as cas-core's [`Language`] hook (CAS-614): natives that
//! turn expressions into text and back, `string`, `expr`, `cat` and `join`,
//! `latex` (CAS-615, [`Form::Latex`]), `export_mathml` and `mathml`
//! (CAS-616), call it through
//! the context, since cas-core cannot depend on this crate. The session and
//! the oracle install it ([`Context::set_language`]); it is part of the
//! context's state, so a context with another hook has another id and
//! shares no cached values.

use std::sync::{Arc, LazyLock};

use cas_core::{Context, Expr, Form, Language, Symbol};

use crate::literal::decimal_balls;
use crate::lower::{Severity, lower};
use crate::markup::{MathmlKind, export_mathml, latex};
use crate::mathml::mathml;
use crate::names::Names;
use crate::parser::parse;
use crate::print::{Balls, print_capped};

/// The language with the names of one locale: the input form and LaTeX
/// printed with balls as they show at the context's working precision, and
/// text read as one expression.
#[derive(Clone, Debug)]
pub struct Hook {
    names: Names,
}

impl Hook {
    /// The hook with `names`.
    #[must_use]
    pub fn new(names: Names) -> Hook {
        Hook { names }
    }

    /// The hook with the English names, one instance shared by all callers.
    #[must_use]
    pub fn english() -> Arc<dyn Language> {
        static ENGLISH: LazyLock<Arc<dyn Language>> =
            LazyLock::new(|| Arc::new(Hook::new(Names::english().clone())));
        Arc::clone(&ENGLISH)
    }
}

impl Language for Hook {
    fn print(
        &self,
        e: &Expr,
        form: Form,
        context: &Context,
        max_len: usize,
    ) -> Result<String, String> {
        let balls = Balls::Shown {
            working_precision: context.working_precision(),
        };
        match form {
            Form::Input => {
                print_capped(e, &self.names, balls, max_len).map_err(|error| error.to_string())
            }
            Form::Latex => latex(e, &self.names, balls, max_len).map_err(|error| error.to_string()),
            Form::Mathml | Form::MathmlDisplay | Form::MathmlContent => {
                let kind = match form {
                    Form::Mathml => MathmlKind::Default,
                    Form::MathmlDisplay => MathmlKind::Display,
                    _ => MathmlKind::Content,
                };
                export_mathml(e, &self.names, balls, max_len, kind)
                    .map_err(|error| error.to_string())
            }
            Form::MathmlDocument => {
                mathml(e, &self.names, balls, max_len).map_err(|error| error.to_string())
            }
            _ => Err(format!("no printer for the form {form:?}")),
        }
    }

    /// One expression: syntax and lowering errors, no statement or several,
    /// and statements (`:=`, `assume`, `purge`) are refused.
    fn read(&self, text: &str, context: &Context) -> Result<Expr, String> {
        let parsed = parse(text);
        if let Some(error) = parsed.errors().first() {
            return Err(error.to_string());
        }
        let lowered = lower(&parsed, &self.names);
        if let Some(error) = lowered
            .diagnostics
            .iter()
            .find(|d| d.severity == Severity::Error)
        {
            return Err(error.to_string());
        }
        let [statement] = &lowered.statements[..] else {
            return Err(format!(
                "one expression expected, not {} statements",
                lowered.statements.len()
            ));
        };
        let head = statement.expr.head().and_then(Expr::as_symbol);
        if matches!(
            head,
            Some(
                Symbol::SET
                    | Symbol::SET_DELAYED
                    | Symbol::UP_SET_DELAYED
                    | Symbol::ASSUME
                    | Symbol::CLEAR
            )
        ) {
            return Err("a statement, not an expression".to_owned());
        }
        Ok(decimal_balls(&statement.expr, context.working_precision()))
    }
}

#[cfg(test)]
mod tests;
