//! Random variates of the distributions (CAS-DOC-33 §4, CAS-372): Xcas's
//! `randnorm`, `randexp`, `randbinomial`, `randpoisson`, `randgeometric`,
//! `randchisquare`, `randstudent` and `randfisher`, and the laws of
//! `randvector`, `randmatrix` and `randpoly`, named by their densities.
//!
//! A variate is the inverse cumulative function of its law (CAS-DOC-23) at
//! a uniform `u` in `(0, 1)`, an exact ball, made one number by `N` where
//! it keeps an exact factor, as the normal's `√2`: the inverse checks the
//! parameters, and an invalid one leaves the call with its warning
//! (`randbinomial(10, 2)`, where giac errs, and `randexp(-1)`, which giac
//! makes negative); a symbol leaves it too. `u` is drawn again where it is
//! 0, and, up to 8 times, where a discrete law's inverse stays though its
//! parameters are numbers, balls among them: its comparisons with `u`
//! left undecided by their radii. With exact parameters they decide.

use super::unit;
use crate::context::EvalError;
use crate::eval::Eval;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, EvalError>;

/// The draws of one variate at most.
const TRIES: usize = 8;

/// A law of CAS-DOC-23.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(in crate::native) enum Law {
    Normal,
    Exponential,
    Uniform,
    Binomial,
    Poisson,
    Geometric,
    NegativeBinomial,
    StudentT,
    ChiSquare,
    Fisher,
}

impl Law {
    /// The law whose density, or mass, is `head`: giac's names of the laws
    /// in `randvector` (`normald`, `binomial`, …).
    pub(in crate::native) fn of_density(head: Symbol) -> Option<Law> {
        Some(match head {
            Symbol::NORMAL_PDF => Law::Normal,
            Symbol::EXPONENTIAL_PDF => Law::Exponential,
            Symbol::UNIFORM_PDF => Law::Uniform,
            Symbol::BINOMIAL => Law::Binomial,
            Symbol::POISSON_PDF => Law::Poisson,
            Symbol::GEOMETRIC_PDF => Law::Geometric,
            Symbol::NEGATIVE_BINOMIAL_PDF => Law::NegativeBinomial,
            Symbol::STUDENT_T_PDF => Law::StudentT,
            Symbol::CHI_SQUARE_PDF => Law::ChiSquare,
            Symbol::FISHER_PDF => Law::Fisher,
            _ => return None,
        })
    }

    /// Its inverse cumulative function.
    fn inverse(self) -> Symbol {
        match self {
            Law::Normal => Symbol::NORMAL_ICDF,
            Law::Exponential => Symbol::EXPONENTIAL_ICDF,
            Law::Uniform => Symbol::UNIFORM_ICDF,
            Law::Binomial => Symbol::BINOMIAL_ICDF,
            Law::Poisson => Symbol::POISSON_ICDF,
            Law::Geometric => Symbol::GEOMETRIC_ICDF,
            Law::NegativeBinomial => Symbol::NEGATIVE_BINOMIAL_ICDF,
            Law::StudentT => Symbol::STUDENT_T_ICDF,
            Law::ChiSquare => Symbol::CHI_SQUARE_ICDF,
            Law::Fisher => Symbol::FISHER_ICDF,
        }
    }

    /// Whether it takes `n` parameters: the normal 0, standard, or 2.
    fn takes(self, n: usize) -> bool {
        match self {
            Law::Normal => n == 0 || n == 2,
            Law::Uniform | Law::Binomial | Law::NegativeBinomial | Law::Fisher => n == 2,
            Law::Exponential | Law::Poisson | Law::Geometric | Law::StudentT | Law::ChiSquare => {
                n == 1
            }
        }
    }

    fn discrete(self) -> bool {
        matches!(
            self,
            Law::Binomial | Law::Poisson | Law::Geometric | Law::NegativeBinomial
        )
    }
}

/// A variate of `law` with `parameters`; `None` when the inverse stays.
pub(in crate::native) fn variate(ev: &mut Eval<'_>, law: Law, parameters: &[Expr]) -> Native {
    if !law.takes(parameters.len()) {
        return Ok(None);
    }
    let balls = parameters.iter().all(Expr::is_number)
        && parameters
            .iter()
            .any(|p| matches!(p.view(), View::Ball(_) | View::ComplexBall(_)));
    let retry = law.discrete() && balls;
    let inverse = Expr::from(law.inverse());
    for _ in 0..TRIES {
        let u = unit(ev);
        if u.mid().is_zero() {
            continue;
        }
        let mut args = parameters.to_vec();
        args.push(Expr::from(u));
        let warned = ev.warnings_given();
        let value = ev.eval(&Expr::apply(&inverse, &args))?;
        if !value.contains(&inverse) {
            // The normal's inverse keeps √2 exact beside its ball: one
            // number by `N`.
            if !value.is_number() && value.is_numeric() {
                return ev.eval(&Expr::call(Symbol::N, [value])).map(Some);
            }
            return Ok(Some(value));
        }
        // Invalid parameters warn, and symbols stay: no draw tells more.
        if !retry || ev.warnings_given() != warned {
            return Ok(None);
        }
    }
    Ok(None)
}

/// `RandomNormal[]`, `RandomNormal[μ, σ]`.
pub(in crate::native) fn random_normal(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Normal, e.args())
}

/// `RandomExponential[λ]`.
pub(in crate::native) fn random_exponential(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Exponential, e.args())
}

/// `RandomBinomial[n, p]`.
pub(in crate::native) fn random_binomial(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Binomial, e.args())
}

/// `RandomPoisson[λ]`.
pub(in crate::native) fn random_poisson(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Poisson, e.args())
}

/// `RandomGeometric[p]`.
pub(in crate::native) fn random_geometric(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Geometric, e.args())
}

/// `RandomChiSquare[ν]`.
pub(in crate::native) fn random_chi_square(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::ChiSquare, e.args())
}

/// `RandomStudentT[ν]`.
pub(in crate::native) fn random_student_t(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::StudentT, e.args())
}

/// `RandomFisher[ν₁, ν₂]`.
pub(in crate::native) fn random_fisher(ev: &mut Eval<'_>, e: &Expr) -> Native {
    variate(ev, Law::Fisher, e.args())
}
