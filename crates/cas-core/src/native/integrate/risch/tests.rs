//! Tests of `Risch` (CAS-585) on the probe rows of T9.c1
//! (.run/1.3/T9/probe.out, giac 2.0.0.21): each answer equals giac's as a
//! value and its derivative is the integrand; partial answers split into
//! the primitive and `Risch[rem, x]`, both giac's; what giac leaves stays.
//! Each test collects every failing row and panics once with the list; a
//! row that panics is one failure, with its message, and the hook prints
//! the panic's backtrace (CAS-585 T9.c4).
//!
//! Departures, not tested (the port cannot reach giac's answer by design):
// departure: 1/(x^2-2) — the resultant factor 8t^2-1 needs a sqrt split (giac's with_sqrt), not ported
// departure: 1/(x^3+1) — a degree-2 resultant factor over Q(i) needs a sqrt(3) split, not ported
// departure: 1/(x^4+1) — degree-2 resultant factors need a sqrt(2) split, not ported
// departure: x^3/(x^2+x+1) — a degree-2 resultant factor needs a sqrt(3) split, not ported
// departure: exp(x)/(x^2+1) — giac writes re/im of Ei(x ± i) in its re(halftan) form with tan(1/2)
// departure: 1/(exp(x)+exp(-x)) — giac's -pi+atan(exp(x)) carries the constant of its re(halftan) step
// departure: exp(x)/(exp(2*x)+1) — the same: giac's -pi+atan(exp(x)), a constant apart at best
// departure: tower-size ties — kernels of one size keep our canonical order (Exp before Log, as giac's
//   printout for exp/ln; ln(2x-1)/ln(3x-1) and exp(2x)/ln(x-2) ties of T9.c1's batch differ), not probe rows

use std::panic::{self, AssertUnwindSafe};
use std::sync::Once;

use crate::context::Context;
use crate::expr::{Expr, View};
use crate::symbol::Symbol;

/// A panic hook printing each panic with a forced backtrace, installed once
/// per test binary, so a CI log names the panicking frame.
fn backtrace_hook() {
    static HOOK: Once = Once::new();
    HOOK.call_once(|| {
        panic::set_hook(Box::new(|info| {
            eprintln!("{info}\n{}", std::backtrace::Backtrace::force_capture());
        }));
    });
}

/// `check` on one row under `catch_unwind`: its failure, or the row's panic
/// as one; the context is made anew after a panic.
fn row(
    ctx: &mut Context,
    input: &str,
    check: impl FnOnce(&Context) -> Option<String>,
) -> Option<String> {
    match panic::catch_unwind(AssertUnwindSafe(|| check(&*ctx))) {
        Ok(failure) => failure,
        Err(payload) => {
            *ctx = Context::new();
            let message = payload
                .downcast_ref::<&str>()
                .map(|m| (*m).to_owned())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "a panic without a message".to_owned());
            Some(format!("{input}: panicked: {message}"))
        }
    }
}

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err:?}"))
}

/// `e` with each `E^u` written `Exp[u]` (C0036).
fn exps(e: &Expr) -> Expr {
    e.replace(|z| match z.view() {
        View::Apply { head, args: [b, u] }
            if head.as_symbol() == Some(Symbol::POWER) && b.as_symbol() == Some(Symbol::E) =>
        {
            Some(Expr::call(Symbol::EXP, [u.clone()]))
        }
        _ => None,
    })
}

/// Whether `d` is 0: by `Together`, or after `TrigToExp` and
/// `TranscendentalSimplify`.
fn vanishes(ctx: &Context, d: &Expr) -> bool {
    let d = exps(&eval(ctx, d));
    if eval(ctx, &Expr::call(Symbol::TOGETHER, [d.clone()])) == Expr::from(0) {
        return true;
    }
    let t = exps(&eval(ctx, &Expr::call(Symbol::TRIG_TO_EXP, [d])));
    let t = eval(ctx, &Expr::call(Symbol::TRANSCENDENTAL_SIMPLIFY, [t]));
    eval(ctx, &Expr::call(Symbol::TOGETHER, [t])) == Expr::from(0)
}

/// Whether `D[prim, x] + rem − f` is 0.
fn derivative_holds(ctx: &Context, prim: &Expr, rem: &Expr, f: &Expr) -> bool {
    let x = read("x");
    let d = Expr::call(Symbol::D, [prim.clone(), x]);
    vanishes(ctx, &(d + rem.clone() - f.clone()))
}

/// `Risch[f, x]`, evaluated.
fn risch(ctx: &Context, f: &Expr) -> Expr {
    eval(ctx, &Expr::call(Symbol::RISCH, [f.clone(), read("x")]))
}

fn is_risch(e: &Expr) -> bool {
    e.head() == Some(&Expr::from(Symbol::RISCH))
}

/// `(prim, rem)` of `prim + Risch[rem, x]`: the terms without `Risch`, the
/// sum of the integrands of those with it.
fn split(e: &Expr) -> (Expr, Expr) {
    let terms: Vec<Expr> = if e.head() == Some(&Expr::from(Symbol::PLUS)) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    };
    let (rems, prims): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(is_risch);
    let rem: Vec<Expr> = rems.iter().map(|r| r.args()[0].clone()).collect();
    (
        Expr::call(Symbol::PLUS, prims),
        Expr::call(Symbol::PLUS, rem),
    )
}

/// Each `(f, giac's answer)` integrates to giac's value, with `f` as its
/// derivative.
fn integrated(rows: &[(&str, &str)]) {
    backtrace_hook();
    let mut ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in rows {
        wrong.extend(row(&mut ctx, input, |ctx| {
            let f = eval(ctx, &read(input));
            let found = risch(ctx, &f);
            let expected = eval(ctx, &read(expected));
            if found.distinct().any(is_risch) {
                return Some(format!("{input}: left {found:?}, expected {expected:?}"));
            }
            let equal = vanishes(ctx, &(found.clone() - expected.clone()));
            let derivative = derivative_holds(ctx, &found, &Expr::from(0), &f);
            (!equal || !derivative).then(|| {
                format!(
                    "{input}: {found:?}, expected {expected:?} (equal {equal}, derivative {derivative})"
                )
            })
        }));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// risch(ln(x),x) = -x+x*ln(x) and the other required rows (`PHASE_PLAN` T9).
#[test]
fn required_rows() {
    integrated(REQUIRED);
}

/// The other probe rows giac integrates.
#[test]
fn other_probe_rows() {
    integrated(OTHERS);
}

/// exp(x^2), 1/ln(x), sqrt(x), x^2*exp(-x^2), (x^2+1)/(x^3+x+1) and the
/// other rows giac leaves stay `Risch[f, x]`.
#[test]
fn unevaluated_rows() {
    backtrace_hook();
    let mut ctx = Context::new();
    let mut wrong = Vec::new();
    for input in UNEVALUATED {
        wrong.extend(row(&mut ctx, input, |ctx| {
            let f = eval(ctx, &read(input));
            let found = risch(ctx, &f);
            (!is_risch(&found) || found.args() != [f.clone(), read("x")])
                .then(|| format!("{input}: {found:?}, expected it to stay"))
        }));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Partial answers: giac's primitive and remainder, as values, and
/// `D[prim] + rem = f`.
#[test]
fn partial_rows() {
    backtrace_hook();
    let mut ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in PARTIAL {
        wrong.extend(row(&mut ctx, input, |ctx| {
            let f = eval(ctx, &read(input));
            let found = risch(ctx, &f);
            let (prim, rem) = split(&found);
            let (gprim, grem) = split(&read(expected));
            let (gprim, grem) = (eval(ctx, &gprim), eval(ctx, &grem));
            let prim_equal = vanishes(ctx, &(prim.clone() - gprim.clone()));
            let rem_equal = vanishes(ctx, &(rem.clone() - grem.clone()));
            let derivative = derivative_holds(ctx, &prim, &rem, &f);
            (rem == Expr::from(0) || !prim_equal || !rem_equal || !derivative).then(|| {
                format!(
                    "{input}: {found:?}, expected {gprim:?} + Risch[{grem:?}, x] \
                     (prim {prim_equal}, rem {rem_equal}, derivative {derivative})"
                )
            })
        }));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The required rows (`PHASE_PLAN` T9): `(f, giac's answer)`.
const REQUIRED: &[(&str, &str)] = &[
    // probe row 1: risch(ln(x),x) => -x+x*ln(x)
    ("Log[x]", "Plus[Times[-1, x], Times[x, Log[x]]]"),
    // probe row 2: risch(x*exp(x),x) => x*exp(x)-exp(x)
    (
        "Times[x, Exp[x]]",
        "Plus[Times[x, Exp[x]], Times[-1, Exp[x]]]",
    ),
    // probe row 3: risch(exp(x)*sin(x),x) => -1/2*cos(x)*exp(x)+1/2*exp(x)*sin(x)
    (
        "Times[Exp[x], Sin[x]]",
        "Plus[Times[Times[Times[-1, Power[2, -1]], Cos[x]], Exp[x]], Times[Times[1/2, Exp[x]], Sin[x]]]",
    ),
    // probe row 4: risch(exp(2*x)/(1+exp(x)),x) => exp(x)-ln(exp(x)+1)
    (
        "Times[Exp[Times[2, x]], Power[Plus[1, Exp[x]], -1]]",
        "Plus[Exp[x], Times[-1, Log[Plus[Exp[x], 1]]]]",
    ),
    // probe row 5: risch(ln(ln(x))/x,x) => -ln(x)+ln(x)*ln(ln(x))
    (
        "Times[Log[Log[x]], Power[x, -1]]",
        "Plus[Times[-1, Log[x]], Times[Log[x], Log[Log[x]]]]",
    ),
    // probe row 6: risch(exp(x)*ln(x)+exp(x)/x,x) => exp(x)*ln(x)
    (
        "Plus[Times[Exp[x], Log[x]], Times[Exp[x], Power[x, -1]]]",
        "Times[Exp[x], Log[x]]",
    ),
    // probe row 7: risch(exp(x)/x,x) => Ei(x)
    ("Times[Exp[x], Power[x, -1]]", "ExpIntegralEi[x]"),
    // probe row 8: risch(exp(x)/x^2,x) => (x*Ei(x)-exp(x))/x
    (
        "Times[Exp[x], Power[Power[x, 2], -1]]",
        "Times[Plus[Times[x, ExpIntegralEi[x]], Times[-1, Exp[x]]], Power[x, -1]]",
    ),
];

/// The other probe rows giac integrates.
const OTHERS: &[(&str, &str)] = &[
    // probe row 14: risch(x*exp(x^2),x) => 1/2*exp(x^2)
    ("Times[x, Exp[Power[x, 2]]]", "Times[1/2, Exp[Power[x, 2]]]"),
    // probe row 15: risch(ln(x)^2,x) => x*ln(x)^2-2*x*ln(x)+2*x
    (
        "Power[Log[x], 2]",
        "Plus[Plus[Times[x, Power[Log[x], 2]], Times[-1, Times[Times[2, x], Log[x]]]], Times[2, x]]",
    ),
    // probe row 16: risch(1/(x*ln(x)),x) => ln(ln(x))
    ("Times[1, Power[Times[x, Log[x]], -1]]", "Log[Log[x]]"),
    // probe row 17: risch(exp(exp(x))*exp(x),x) => exp(exp(x))
    ("Times[Exp[Exp[x]], Exp[x]]", "Exp[Exp[x]]"),
    // probe row 18: risch(1/(1+exp(x)),x) => x-ln(exp(x)+1)
    (
        "Times[1, Power[Plus[1, Exp[x]], -1]]",
        "Plus[x, Times[-1, Log[Plus[Exp[x], 1]]]]",
    ),
    // probe row 19: risch(exp(x)/(1+exp(x))^2,x) => -1/(exp(x)+1)
    (
        "Times[Exp[x], Power[Power[Plus[1, Exp[x]], 2], -1]]",
        "Times[-1, Power[Plus[Exp[x], 1], -1]]",
    ),
    // probe row 20: risch(ln(x)/x,x) => 1/2*ln(x)^2
    (
        "Times[Log[x], Power[x, -1]]",
        "Times[1/2, Power[Log[x], 2]]",
    ),
    // probe row 21: risch(x^2*exp(x),x) => x^2*exp(x)-2*x*exp(x)+2*exp(x)
    (
        "Times[Power[x, 2], Exp[x]]",
        "Plus[Plus[Times[Power[x, 2], Exp[x]], Times[-1, Times[Times[2, x], Exp[x]]]], Times[2, Exp[x]]]",
    ),
    // probe row 22: risch(exp(-x)*x^3,x) => -x^3*exp(-x)-3*x^2*exp(-x)-6*x*exp(-x)-6*exp(-x)
    (
        "Times[Exp[Times[-1, x]], Power[x, 3]]",
        "Plus[Plus[Plus[Times[Times[-1, Power[x, 3]], Exp[Times[-1, x]]], Times[-1, Times[Times[3, Power[x, 2]], Exp[Times[-1, x]]]]], Times[-1, Times[Times[6, x], Exp[Times[-1, x]]]]], Times[-1, Times[6, Exp[Times[-1, x]]]]]",
    ),
    // probe row 23: risch(ln(x)^3,x) => x*ln(x)^3-3*x*ln(x)^2+6*x*ln(x)-6*x
    (
        "Power[Log[x], 3]",
        "Plus[Plus[Plus[Times[x, Power[Log[x], 3]], Times[-1, Times[Times[3, x], Power[Log[x], 2]]]], Times[Times[6, x], Log[x]]], Times[-1, Times[6, x]]]",
    ),
    // probe row 24: risch(x*ln(x)^2,x) => 1/2*x^2*ln(x)^2-1/2*x^2*ln(x)+1/4*x^2
    (
        "Times[x, Power[Log[x], 2]]",
        "Plus[Plus[Times[Times[1/2, Power[x, 2]], Power[Log[x], 2]], Times[-1, Times[Times[1/2, Power[x, 2]], Log[x]]]], Times[1/4, Power[x, 2]]]",
    ),
    // probe row 25: risch(ln(x)^2/x,x) => 1/3*ln(x)^3
    (
        "Times[Power[Log[x], 2], Power[x, -1]]",
        "Times[1/3, Power[Log[x], 3]]",
    ),
    // probe row 27: risch(ln(x^2+1),x) => (-i)*ln((-i)*x-1)+i*ln(i*x-1)-2*x+x*ln(x^2+1)
    (
        "Log[Plus[Power[x, 2], 1]]",
        "Plus[Plus[Plus[Times[Times[-1, I], Log[Plus[Times[Times[-1, I], x], Times[-1, 1]]]], Times[I, Log[Plus[Times[I, x], Times[-1, 1]]]]], Times[-1, Times[2, x]]], Times[x, Log[Plus[Power[x, 2], 1]]]]",
    ),
    // probe row 29: risch(1/(x*ln(x)^2),x) => -1/ln(x)
    (
        "Times[1, Power[Times[x, Power[Log[x], 2]], -1]]",
        "Times[-1, Power[Log[x], -1]]",
    ),
    // probe row 30: risch(ln(x)/(x*(ln(x)+1)),x) => -ln(ln(x)+1)+ln(x)
    (
        "Times[Log[x], Power[Times[x, Plus[Log[x], 1]], -1]]",
        "Plus[Times[-1, Log[Plus[Log[x], 1]]], Log[x]]",
    ),
    // probe row 31: risch(1/(x*ln(x)*ln(ln(x))),x) => ln(ln(ln(x)))
    (
        "Times[1, Power[Times[Times[x, Log[x]], Log[Log[x]]], -1]]",
        "Log[Log[Log[x]]]",
    ),
    // probe row 32: risch(ln(x)/(x+1)^2,x) => -ln(x+1)+ln(x)+1/(-x-1)*ln(x)
    (
        "Times[Log[x], Power[Power[Plus[x, 1], 2], -1]]",
        "Plus[Plus[Times[-1, Log[Plus[x, 1]]], Log[x]], Times[Times[1, Power[Plus[Times[-1, x], Times[-1, 1]], -1]], Log[x]]]",
    ),
    // probe row 33: risch(ln(x)/x^2,x) => -1/x-1/x*ln(x)
    (
        "Times[Log[x], Power[Power[x, 2], -1]]",
        "Plus[Times[-1, Power[x, -1]], Times[-1, Times[Times[1, Power[x, -1]], Log[x]]]]",
    ),
    // probe row 34: risch(exp(-x),x) => -exp(-x)
    ("Exp[Times[-1, x]]", "Times[-1, Exp[Times[-1, x]]]"),
    // probe row 35: risch(exp(x)+exp(-x),x) => (exp(x)^2-1)/exp(x)
    (
        "Plus[Exp[x], Exp[Times[-1, x]]]",
        "Times[Plus[Power[Exp[x], 2], Times[-1, 1]], Power[Exp[x], -1]]",
    ),
    // probe row 36: risch(exp(2*x)/exp(x),x) => exp(x)
    ("Times[Exp[Times[2, x]], Power[Exp[x], -1]]", "Exp[x]"),
    // probe row 37: risch(exp(x)*exp(2*x),x) => 1/3*exp(x)^3
    (
        "Times[Exp[x], Exp[Times[2, x]]]",
        "Times[1/3, Power[Exp[x], 3]]",
    ),
    // probe row 40: risch(1/(exp(x)-1),x) => -x+ln(exp(x)-1)
    (
        "Times[1, Power[Plus[Exp[x], Times[-1, 1]], -1]]",
        "Plus[Times[-1, x], Log[Plus[Exp[x], Times[-1, 1]]]]",
    ),
    // probe row 41: risch(exp(3*x)/(exp(x)+1),x) => 1/2*exp(x)^2-exp(x)+ln(exp(x)+1)
    (
        "Times[Exp[Times[3, x]], Power[Plus[Exp[x], 1], -1]]",
        "Plus[Plus[Times[1/2, Power[Exp[x], 2]], Times[-1, Exp[x]]], Log[Plus[Exp[x], 1]]]",
    ),
    // probe row 42: risch(exp(-2*x)+x*exp(-x),x) => -x*exp(-x)-1/2*exp(-x)^2-exp(-x)
    (
        "Plus[Exp[Times[-2, x]], Times[x, Exp[Times[-1, x]]]]",
        "Plus[Plus[Times[Times[-1, x], Exp[Times[-1, x]]], Times[-1, Times[1/2, Power[Exp[Times[-1, x]], 2]]]], Times[-1, Exp[Times[-1, x]]]]",
    ),
    // probe row 44: risch(exp(x)*exp(exp(x))*exp(exp(exp(x))),x) => exp(exp(exp(x)))
    (
        "Times[Times[Exp[x], Exp[Exp[x]]], Exp[Exp[Exp[x]]]]",
        "Exp[Exp[Exp[x]]]",
    ),
    // probe row 45: risch((2*x+1)*exp(x^2+x),x) => exp(x^2+x)
    (
        "Times[Plus[Times[2, x], 1], Exp[Plus[Power[x, 2], x]]]",
        "Exp[Plus[Power[x, 2], x]]",
    ),
    // probe row 46: risch(x^3*exp(x^2),x) => 1/2*x^2*exp(x^2)-1/2*exp(x^2)
    (
        "Times[Power[x, 3], Exp[Power[x, 2]]]",
        "Plus[Times[Times[1/2, Power[x, 2]], Exp[Power[x, 2]]], Times[-1, Times[1/2, Exp[Power[x, 2]]]]]",
    ),
    // probe row 47: risch(exp(x+1),x) => exp(1)*exp(x)
    ("Exp[Plus[x, 1]]", "Times[Exp[1], Exp[x]]"),
    // probe row 48: risch(x*exp(2*x+1),x) => 1/2*x*exp(1)*exp(2*x)-1/4*exp(1)*exp(2*x)
    (
        "Times[x, Exp[Plus[Times[2, x], 1]]]",
        "Plus[Times[Times[Times[1/2, x], Exp[1]], Exp[Times[2, x]]], Times[-1, Times[Times[1/4, Exp[1]], Exp[Times[2, x]]]]]",
    ),
    // probe row 49: risch(exp(x)/(x-1),x) => Ei(x-1)*exp(1)
    (
        "Times[Exp[x], Power[Plus[x, Times[-1, 1]], -1]]",
        "Times[ExpIntegralEi[Plus[x, Times[-1, 1]]], Exp[1]]",
    ),
    // probe row 50: risch(exp(2*x+1)/x,x) => Ei(2*x)*exp(1)
    (
        "Times[Exp[Plus[Times[2, x], 1]], Power[x, -1]]",
        "Times[ExpIntegralEi[Times[2, x]], Exp[1]]",
    ),
    // probe row 51: risch(exp(x)/x^3,x) => (x^2*Ei(x)-x*exp(x)-exp(x))/(2*x^2)
    (
        "Times[Exp[x], Power[Power[x, 3], -1]]",
        "Times[Plus[Plus[Times[Power[x, 2], ExpIntegralEi[x]], Times[-1, Times[x, Exp[x]]]], Times[-1, Exp[x]]], Power[Times[2, Power[x, 2]], -1]]",
    ),
    // probe row 52: risch(exp(-x)/x,x) => Ei(-x)
    (
        "Times[Exp[Times[-1, x]], Power[x, -1]]",
        "ExpIntegralEi[Times[-1, x]]",
    ),
    // probe row 53: risch(exp(x)/(x^2-1),x) => (-Ei(x+1)+Ei(x-1)*exp(1)^2)/(2*exp(1))
    (
        "Times[Exp[x], Power[Plus[Power[x, 2], Times[-1, 1]], -1]]",
        "Times[Plus[Times[-1, ExpIntegralEi[Plus[x, 1]]], Times[ExpIntegralEi[Plus[x, Times[-1, 1]]], Power[Exp[1], 2]]], Power[Times[2, Exp[1]], -1]]",
    ),
    // probe row 55: risch(x*exp(x)/(x+1)^2,x) => exp(x)/(x+1)
    (
        "Times[Times[x, Exp[x]], Power[Power[Plus[x, 1], 2], -1]]",
        "Times[Exp[x], Power[Plus[x, 1], -1]]",
    ),
    // probe row 56: risch(exp(i*x)/x,x) => Ci(x)+i*Si(x)
    (
        "Times[Exp[Times[I, x]], Power[x, -1]]",
        "Plus[CosIntegral[x], Times[I, SinIntegral[x]]]",
    ),
    // probe row 57: risch(sin(x)/x,x) => Si(x)
    ("Times[Sin[x], Power[x, -1]]", "SinIntegral[x]"),
    // probe row 58: risch(cos(x)/x,x) => Ci(x)
    ("Times[Cos[x], Power[x, -1]]", "CosIntegral[x]"),
    // probe row 59: risch(cos(x)*exp(2*x),x) => 2/5*cos(x)*exp(2*x)+1/5*exp(2*x)*sin(x)
    (
        "Times[Cos[x], Exp[Times[2, x]]]",
        "Plus[Times[Times[2/5, Cos[x]], Exp[Times[2, x]]], Times[Times[1/5, Exp[Times[2, x]]], Sin[x]]]",
    ),
    // probe row 60: risch(sin(x)^2,x) => (x*tan(x/2)^4+2*x*tan(x/2)^2+x+2*tan(x/2)^3-2*tan(x/2))/(2*tan(x/2)^4+4*tan(x/2)^2+2)
    (
        "Power[Sin[x], 2]",
        "Times[Plus[Plus[Plus[Plus[Times[x, Power[Tan[Times[x, Power[2, -1]]], 4]], Times[Times[2, x], Power[Tan[Times[x, Power[2, -1]]], 2]]], x], Times[2, Power[Tan[Times[x, Power[2, -1]]], 3]]], Times[-1, Times[2, Tan[Times[x, Power[2, -1]]]]]], Power[Plus[Plus[Times[2, Power[Tan[Times[x, Power[2, -1]]], 4]], Times[4, Power[Tan[Times[x, Power[2, -1]]], 2]]], 2], -1]]",
    ),
    // probe row 61: risch(x*cos(x),x) => x*sin(x)+cos(x)
    ("Times[x, Cos[x]]", "Plus[Times[x, Sin[x]], Cos[x]]"),
    // probe row 62: risch(a*exp(a*x),x) => exp(a*x)
    ("Times[a, Exp[Times[a, x]]]", "Exp[Times[a, x]]"),
    // probe row 63: risch(exp(a*x),x) => exp(a*x)/a
    ("Exp[Times[a, x]]", "Times[Exp[Times[a, x]], Power[a, -1]]"),
    // probe row 64: risch(x*exp(a*x),x) => (a*x*exp(a*x)-exp(a*x))/a^2
    (
        "Times[x, Exp[Times[a, x]]]",
        "Times[Plus[Times[Times[a, x], Exp[Times[a, x]]], Times[-1, Exp[Times[a, x]]]], Power[Power[a, 2], -1]]",
    ),
    // probe row 65: risch(ln(a*x),x) => -x+x*ln(a*x)
    (
        "Log[Times[a, x]]",
        "Plus[Times[-1, x], Times[x, Log[Times[a, x]]]]",
    ),
    // probe row 66: risch(1/(x+a),x) => ln(a+x)
    ("Times[1, Power[Plus[x, a], -1]]", "Log[Plus[a, x]]"),
    // probe row 67: risch(atan(x),x) => (2*x*atan(x)-ln(x^2+1))/2
    (
        "ArcTan[x]",
        "Times[Plus[Times[Times[2, x], ArcTan[x]], Times[-1, Log[Plus[Power[x, 2], 1]]]], Power[2, -1]]",
    ),
    // probe row 68: risch(x*atan(x),x) => (2*x^2*atan(x)-2*x+(-i)*ln((-i)*x-1)+i*ln(i*x-1))/4
    (
        "Times[x, ArcTan[x]]",
        "Times[Plus[Plus[Plus[Times[Times[2, Power[x, 2]], ArcTan[x]], Times[-1, Times[2, x]]], Times[Times[-1, I], Log[Plus[Times[Times[-1, I], x], Times[-1, 1]]]]], Times[I, Log[Plus[Times[I, x], Times[-1, 1]]]]], Power[4, -1]]",
    ),
    // probe row 71: risch(1/(x^2-1),x) => 1/2*ln(x-1)-1/2*ln(x+1)
    (
        "Times[1, Power[Plus[Power[x, 2], Times[-1, 1]], -1]]",
        "Plus[Times[1/2, Log[Plus[x, Times[-1, 1]]]], Times[-1, Times[1/2, Log[Plus[x, 1]]]]]",
    ),
    // probe row 72: risch(1/(x^2+1),x) => (-i)/2*ln((-i)*x-1)+i/2*ln(i*x-1)
    (
        "Times[1, Power[Plus[Power[x, 2], 1], -1]]",
        "Plus[Times[Times[Times[-1, I], Power[2, -1]], Log[Plus[Times[Times[-1, I], x], Times[-1, 1]]]], Times[Times[I, Power[2, -1]], Log[Plus[Times[I, x], Times[-1, 1]]]]]",
    ),
    // probe row 74: risch(1/(x^3-x),x) => 1/2*ln(x^2-1)-ln(x)
    (
        "Times[1, Power[Plus[Power[x, 3], Times[-1, x]], -1]]",
        "Plus[Times[1/2, Log[Plus[Power[x, 2], Times[-1, 1]]]], Times[-1, Log[x]]]",
    ),
    // probe row 75: risch((3*x^2+1)/(x^3+x+1),x) => ln(x^3+x+1)
    (
        "Times[Plus[Times[3, Power[x, 2]], 1], Power[Plus[Plus[Power[x, 3], x], 1], -1]]",
        "Log[Plus[Plus[Power[x, 3], x], 1]]",
    ),
    // probe row 76: risch(x/(x^2+1)^2,x) => 1/(-2*x^2-2)
    (
        "Times[x, Power[Power[Plus[Power[x, 2], 1], 2], -1]]",
        "Times[1, Power[Plus[Times[-2, Power[x, 2]], Times[-1, 2]], -1]]",
    ),
    // probe row 77: risch(1/(x^2+1)^2,x) => (-i)/4*ln((-i)*x-1)+i/4*ln(i*x-1)-x/(-2*x^2-2)
    (
        "Times[1, Power[Power[Plus[Power[x, 2], 1], 2], -1]]",
        "Plus[Plus[Times[Times[Times[-1, I], Power[4, -1]], Log[Plus[Times[Times[-1, I], x], Times[-1, 1]]]], Times[Times[I, Power[4, -1]], Log[Plus[Times[I, x], Times[-1, 1]]]]], Times[-1, Times[x, Power[Plus[Times[-2, Power[x, 2]], Times[-1, 2]], -1]]]]",
    ),
    // probe row 82: risch(x*exp(x)*ln(x),x) => x*exp(x)*ln(x)+Ei(x)-exp(x)*ln(x)-exp(x)
    (
        "Times[Times[x, Exp[x]], Log[x]]",
        "Plus[Plus[Plus[Times[Times[x, Exp[x]], Log[x]], ExpIntegralEi[x]], Times[-1, Times[Exp[x], Log[x]]]], Times[-1, Exp[x]]]",
    ),
    // probe row 83: risch(exp(x)*(1/x-1/x^2),x) => exp(x)/x
    (
        "Times[Exp[x], Plus[Times[1, Power[x, -1]], Times[-1, Times[1, Power[Power[x, 2], -1]]]]]",
        "Times[Exp[x], Power[x, -1]]",
    ),
    // probe row 84: risch((exp(x)+1)/(exp(x)+x),x) => ln(x+exp(x))
    (
        "Times[Plus[Exp[x], 1], Power[Plus[Exp[x], x], -1]]",
        "Log[Plus[x, Exp[x]]]",
    ),
    // probe row 85: risch(exp(x)/(exp(x)+1)+ln(x),x) => x*ln(x)-x+ln(exp(x)+1)
    (
        "Plus[Times[Exp[x], Power[Plus[Exp[x], 1], -1]], Log[x]]",
        "Plus[Plus[Times[x, Log[x]], Times[-1, x]], Log[Plus[Exp[x], 1]]]",
    ),
    // probe row 86: risch(2^x,x) => exp(x*ln(2))/ln(2)
    (
        "Power[2, x]",
        "Times[Exp[Times[x, Log[2]]], Power[Log[2], -1]]",
    ),
];

/// Rows giac integrates in part: `(f, prim + Risch[rem, x])`.
const PARTIAL: &[(&str, &str)] = &[
    // probe row 26: risch(ln(x+1)*ln(x),x) => (-x+x*ln(x))*ln(x+1)+integrate(-x/(x+1)*ln(x)+x/(x+1),x)
    (
        "Times[Log[Plus[x, 1]], Log[x]]",
        "Plus[Times[Plus[Times[-1, x], Times[x, Log[x]]], Log[Plus[x, 1]]], Risch[Plus[Times[Times[Times[-1, x], Power[Plus[x, 1], -1]], Log[x]], Times[x, Power[Plus[x, 1], -1]]], x]]",
    ),
    // probe row 28: risch(ln(ln(x)),x) => x*ln(ln(x))+integrate(-1/ln(x),x)
    (
        "Log[Log[x]]",
        "Plus[Times[x, Log[Log[x]]], Risch[Times[-1, Power[Log[x], -1]], x]]",
    ),
    // probe row 81: risch(exp(x)/(x*ln(x))-exp(x)/(x*ln(x)^2)+exp(x)/x,x) => (Ei(x)*ln(x)+exp(x))/ln(x)+integrate((-exp(x)*x+exp(x))/(ln(x)*x),x)
    (
        "Plus[Plus[Times[Exp[x], Power[Times[x, Log[x]], -1]], Times[-1, Times[Exp[x], Power[Times[x, Power[Log[x], 2]], -1]]]], Times[Exp[x], Power[x, -1]]]",
        "Plus[Times[Plus[Times[ExpIntegralEi[x], Log[x]], Exp[x]], Power[Log[x], -1]], Risch[Times[Plus[Times[Times[-1, Exp[x]], x], Exp[x]], Power[Times[Log[x], x], -1]], x]]",
    ),
];

/// Rows giac leaves: the call stays.
const UNEVALUATED: &[&str] = &[
    // probe row 9: risch(exp(x^2),x)
    "Exp[Power[x, 2]]",
    // probe row 10: risch(1/ln(x),x)
    "Times[1, Power[Log[x], -1]]",
    // probe row 11: risch(sqrt(x),x)
    "Sqrt[x]",
    // probe row 12: risch(x^2*exp(-x^2),x)
    "Times[Power[x, 2], Exp[Times[-1, Power[x, 2]]]]",
    // probe row 13: risch((x^2+1)/(x^3+x+1),x)
    "Times[Plus[Power[x, 2], 1], Power[Plus[Plus[Power[x, 3], x], 1], -1]]",
    // probe row 43: risch(exp(exp(x)),x)
    "Exp[Exp[x]]",
    // probe row 69: risch(x^(1/2)*exp(x),x)
    "Times[Power[x, 1/2], Exp[x]]",
    // probe row 70: risch(exp(x)*sqrt(x),x)
    "Times[Exp[x], Sqrt[x]]",
    // probe row 87: risch(x^x,x)
    "Power[x, x]",
];
