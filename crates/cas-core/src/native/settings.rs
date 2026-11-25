//! Xcas's settings (CAS-545, CAS-617, CAS-618): their bare names read the
//! context, and `CasSetup[]`, giac's `cas_setup()`, lists them. The session
//! sets them; a native only reads the context.
//!
//! `cas_setup()` is giac's list (prog.cc `cas_setup` builder, giac
//! 2.0.0.21): `[approx_mode, complex_variables, complex_mode, angle, 0,
//! [epsilon, 1e-15], Digits, [1, recursion, 0, eval_level],
//! increasing_power, with_sqrt, all_trig_solutions, 1]`; giac's 7th element
//! is its threads, ours 1, and its recursion limit ours.

use cas_num::ball::{Ball, Float};

use crate::context::Context;
use crate::eval::Eval;
use crate::expr::Expr;
use crate::symbol::Symbol;

type Native = Result<Option<Expr>, crate::context::EvalError>;

/// The exact decimal of a finite `x`, 0 otherwise.
fn decimal(x: f64) -> Expr {
    Float::from_f64(x).map_or_else(|| Expr::from(0), |f| Expr::from(Ball::exact(f)))
}

fn flag(on: bool) -> Expr {
    Expr::from(i64::from(on))
}

/// The value of the setting `s`, read by its bare name.
pub(super) fn value(context: &Context, s: Symbol) -> Option<Expr> {
    Some(match s {
        // CAS-151.
        Symbol::DIGITS => {
            let digits = crate::print::shown_digits(context.working_precision());
            Expr::from(i64::try_from(digits).unwrap_or(i64::MAX))
        }
        // CAS-545.
        Symbol::ALL_TRIG_SOLUTIONS => flag(context.all_trig_solutions()),
        Symbol::WITH_SQRT => flag(context.with_sqrt()),
        // CAS-617, CAS-618.
        Symbol::APPROX_MODE => flag(context.approx_mode()),
        Symbol::COMPLEX_MODE => flag(context.complex_mode()),
        Symbol::COMPLEX_VARIABLES => flag(context.complex_variables()),
        Symbol::ANGLE_RADIAN => Expr::from(context.angle() as i64),
        Symbol::KEEP_ALGEXT => flag(context.keep_algext()),
        Symbol::INCREASING_POWER => flag(context.increasing_power()),
        Symbol::EPSILON => decimal(context.epsilon()),
        Symbol::EVAL_LEVEL => Expr::from(i64::try_from(context.eval_level()).unwrap_or(i64::MAX)),
        Symbol::AUTO_SIMPLIFY => context
            .autosimplify()
            .cloned()
            .unwrap_or_else(|| Expr::from(0)),
        Symbol::RUR_CERTIFY => Expr::from(context.rur_certify().saturating_add(1)),
        Symbol::GBASIS_MAX_PAIRS => Expr::from(context.gbasis_max_pairs()),
        Symbol::GBASIS_SIMULT_PRIMES => Expr::from(context.gbasis_simult_primes()),
        _ => return None,
    })
}

/// `CasSetup[]`: the settings as giac's `cas_setup()` lists them; with
/// arguments it stays (the session sets them).
pub(super) fn cas_setup(ev: &mut Eval<'_>, e: &Expr) -> Native {
    if !e.args().is_empty() {
        return Ok(None);
    }
    let context = ev.context();
    let read = |s| value(context, s).unwrap_or_else(|| Expr::from(0));
    let list = |items: Vec<Expr>| Expr::call(Symbol::LIST, items);
    Ok(Some(list(vec![
        read(Symbol::APPROX_MODE),
        read(Symbol::COMPLEX_VARIABLES),
        read(Symbol::COMPLEX_MODE),
        read(Symbol::ANGLE_RADIAN),
        Expr::from(0),
        list(vec![read(Symbol::EPSILON), decimal(1e-15)]),
        read(Symbol::DIGITS),
        list(vec![
            Expr::from(1),
            Expr::from(i64::from(context.limits().recursion)),
            Expr::from(0),
            read(Symbol::EVAL_LEVEL),
        ]),
        read(Symbol::INCREASING_POWER),
        read(Symbol::WITH_SQRT),
        read(Symbol::ALL_TRIG_SOLUTIONS),
        Expr::from(1),
    ])))
}
