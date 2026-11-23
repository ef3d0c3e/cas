//! Xcas's `risch` against giac (CAS-585): giac's answers read back and
//! compared with ours (`Giac::agrees`), or, where the forms differ, by
//! value: the difference made 0 by `Together` (after `TrigToExp` and
//! `TranscendentalSimplify` for the trigonometric rows), else the same
//! change between two points (an antiderivative is a constant apart at
//! most: giac's `ln(i x − 1)` is ours `ln(x + i)`). giac prints what it
//! cannot integrate as `integrate(rem,x)`, which our `Integrate` would
//! evaluate on reading back: its partial answers are read with
//! `integrate(` as the inert `remains(`, and compared by their parts, the
//! primitive and the remainder. Then the departures, each with giac's
//! answer as libgiac prints it. The rows are T9.c1's probes
//! (.run/1.3/T9/probe.out) and its seeded batch (batch.cas).

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
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

/// `e` evaluated, or `None` when the evaluation fails.
fn eval(ctx: &Context, e: &Expr) -> Option<Expr> {
    ctx.eval(e).ok()
}

/// Whether `d` is 0: by `Together`, or after `TrigToExp` and
/// `TranscendentalSimplify`.
fn vanishes(ctx: &Context, d: &Expr) -> bool {
    let Some(d) = eval(ctx, d).map(|d| exps(&d)) else {
        return false;
    };
    if eval(ctx, &Expr::call(Symbol::TOGETHER, [d.clone()])) == Some(Expr::from(0)) {
        return true;
    }
    let Some(t) = eval(ctx, &Expr::call(Symbol::TRIG_TO_EXP, [d])).map(|t| exps(&t)) else {
        return false;
    };
    let Some(t) = eval(ctx, &Expr::call(Symbol::TRANSCENDENTAL_SIMPLIFY, [t])) else {
        return false;
    };
    eval(ctx, &Expr::call(Symbol::TOGETHER, [t])) == Some(Expr::from(0))
}

/// Whether `d` changes by less than 10⁻⁹ from `x = 7/2` to `x = 9/2`,
/// with `a = 2/3`, its real and imaginary parts apart (C0013).
fn constant(ctx: &Context, d: &Expr) -> bool {
    let (x, a) = (Expr::symbol("x"), Expr::symbol("a"));
    let point = |v: i64| {
        d.substitute(&[
            (x.clone(), Expr::from(v) / Expr::from(2)),
            (a.clone(), Expr::from(2) / Expr::from(3)),
        ])
    };
    let gap = Expr::call(Symbol::N, [point(9) - point(7)]);
    let bound = Expr::call(Symbol::N, [Expr::from(1) / Expr::from(1_000_000_000_i64)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let size = Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]);
        eval(ctx, &Expr::call(Symbol::LESS, [size, bound.clone()])).and_then(|e| e.as_symbol())
            == Some(Symbol::TRUE)
    })
}

/// Whether `e` is a remainder: ours `Risch[rem, x]`, giac's
/// `integrate(rem,x)` read as `remains(rem,x)`.
fn is_remainder(e: &Expr) -> bool {
    e.head()
        .and_then(Expr::as_symbol)
        .is_some_and(|h| h == Symbol::RISCH || h == Symbol::new("remains"))
}

/// `(prim, rem)` of `prim + Risch[rem, x]`: the terms without a remainder,
/// the sum of the integrands of those with one.
fn split(e: &Expr) -> (Expr, Expr) {
    let terms: Vec<Expr> = if e.head() == Some(&Expr::from(Symbol::PLUS)) {
        e.args().to_vec()
    } else {
        vec![e.clone()]
    };
    let (rems, prims): (Vec<Expr>, Vec<Expr>) = terms.into_iter().partition(is_remainder);
    let rem: Vec<Expr> = rems.iter().map(|r| r.args()[0].clone()).collect();
    (
        Expr::call(Symbol::PLUS, prims),
        Expr::call(Symbol::PLUS, rem),
    )
}

/// Whether `ours` and `theirs` have the same remainder and primitives
/// equal, or a constant apart.
fn same(ctx: &Context, ours: &Expr, theirs: &Expr) -> bool {
    let (prim, rem) = split(ours);
    let (gprim, grem) = split(theirs);
    let d = prim - gprim;
    vanishes(ctx, &(rem - grem)) && (vanishes(ctx, &d) || constant(ctx, &d))
}

/// giac's answer to `input` in this CAS, `integrate(` read as `remains(`.
fn theirs(giac: &Giac, input: &str) -> Result<Expr, String> {
    let answer = giac.eval(input).map_err(|e| e.to_string())?;
    ours(&answer.replace("integrate(", "remains(")).map_err(|e| e.to_string())
}

/// The integrand of `risch(f,x)`.
fn integrand(input: &str) -> Result<Expr, String> {
    let f = input
        .strip_prefix("risch(")
        .and_then(|s| s.strip_suffix(",x)"))
        .ok_or_else(|| format!("{input}: not risch(f,x)"))?;
    ours(f).map_err(|e| e.to_string())
}

/// Whether `D[prim, x] + rem − f` is 0 for `answer = prim + Risch[rem, x]`.
fn derivative_holds(ctx: &Context, answer: &Expr, f: &Expr) -> bool {
    let (prim, rem) = split(answer);
    let d = Expr::call(Symbol::D, [prim, Expr::symbol("x")]);
    vanishes(ctx, &(d + rem - f.clone()))
}

/// Whether `e` holds a `Risch` call.
fn has_risch(e: &Expr) -> bool {
    e.distinct()
        .any(|z| z.head() == Some(&Expr::from(Symbol::RISCH)))
}

/// `input` against giac: `agrees`, or the same by [`same`]; `partial`
/// rows must keep a remainder, the others none.
fn agreement(giac: &Giac, input: &str, partial: bool) -> Result<(), String> {
    if !partial && giac.agrees(input).is_ok() {
        return Ok(());
    }
    let mine = ours(input).map_err(|e| e.to_string())?;
    let theirs = theirs(giac, input)?;
    if has_risch(&mine) != partial {
        return Err(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
    }
    if same(&Context::new(), &mine, &theirs) {
        Ok(())
    } else {
        Err(format!("{input}: ours {mine:?}, giac's {theirs:?}"))
    }
}

/// Each of `inputs` against giac, every failure listed (C0040).
fn check(inputs: &[&str], partial: bool) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| agreement(giac, input, partial).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The probe rows giac integrates (probe.out row numbers): logarithms,
/// exponentials, `Ei`, `Si` and `Ci`, trigonometric functions through
/// exponentials, parameters, `atan` and rational functions whose residues
/// are in ℚ(i).
#[test]
fn agreements() {
    check(
        &[
            "risch(ln(x),x)",                               // 1
            "risch(x*exp(x),x)",                            // 2
            "risch(exp(x)*sin(x),x)",                       // 3
            "risch(exp(2*x)/(1+exp(x)),x)",                 // 4
            "risch(ln(ln(x))/x,x)",                         // 5
            "risch(exp(x)*ln(x)+exp(x)/x,x)",               // 6
            "risch(exp(x)/x,x)",                            // 7
            "risch(exp(x)/x^2,x)",                          // 8
            "risch(x*exp(x^2),x)",                          // 14
            "risch(ln(x)^2,x)",                             // 15
            "risch(1/(x*ln(x)),x)",                         // 16
            "risch(exp(exp(x))*exp(x),x)",                  // 17
            "risch(1/(1+exp(x)),x)",                        // 18
            "risch(exp(x)/(1+exp(x))^2,x)",                 // 19
            "risch(ln(x)/x,x)",                             // 20
            "risch(x^2*exp(x),x)",                          // 21
            "risch(exp(-x)*x^3,x)",                         // 22
            "risch(ln(x)^3,x)",                             // 23
            "risch(x*ln(x)^2,x)",                           // 24
            "risch(ln(x)^2/x,x)",                           // 25
            "risch(ln(x^2+1),x)",                           // 27
            "risch(1/(x*ln(x)^2),x)",                       // 29
            "risch(ln(x)/(x*(ln(x)+1)),x)",                 // 30
            "risch(1/(x*ln(x)*ln(ln(x))),x)",               // 31
            "risch(ln(x)/(x+1)^2,x)",                       // 32
            "risch(ln(x)/x^2,x)",                           // 33
            "risch(exp(-x),x)",                             // 34
            "risch(exp(x)+exp(-x),x)",                      // 35
            "risch(exp(2*x)/exp(x),x)",                     // 36
            "risch(exp(x)*exp(2*x),x)",                     // 37
            "risch(1/(exp(x)-1),x)",                        // 40
            "risch(exp(3*x)/(exp(x)+1),x)",                 // 41
            "risch(exp(-2*x)+x*exp(-x),x)",                 // 42
            "risch(exp(x)*exp(exp(x))*exp(exp(exp(x))),x)", // 44
            "risch((2*x+1)*exp(x^2+x),x)",                  // 45
            "risch(x^3*exp(x^2),x)",                        // 46
            "risch(exp(x+1),x)",                            // 47
            "risch(x*exp(2*x+1),x)",                        // 48
            "risch(exp(x)/(x-1),x)",                        // 49
            "risch(exp(2*x+1)/x,x)",                        // 50
            "risch(exp(x)/x^3,x)",                          // 51
            "risch(exp(-x)/x,x)",                           // 52
            "risch(exp(x)/(x^2-1),x)",                      // 53
            "risch(x*exp(x)/(x+1)^2,x)",                    // 55
            "risch(exp(i*x)/x,x)",                          // 56
            "risch(sin(x)/x,x)",                            // 57
            "risch(cos(x)/x,x)",                            // 58
            "risch(cos(x)*exp(2*x),x)",                     // 59
            "risch(sin(x)^2,x)",                            // 60
            "risch(x*cos(x),x)",                            // 61
            "risch(a*exp(a*x),x)",                          // 62
            "risch(exp(a*x),x)",                            // 63
            "risch(x*exp(a*x),x)",                          // 64
            "risch(ln(a*x),x)",                             // 65
            "risch(1/(x+a),x)",                             // 66
            "risch(atan(x),x)",                             // 67
            "risch(x*atan(x),x)",                           // 68
            "risch(1/(x^2-1),x)",                           // 71
            "risch(1/(x^2+1),x)",                           // 72
            "risch(1/(x^3-x),x)",                           // 74
            "risch((3*x^2+1)/(x^3+x+1),x)",                 // 75
            "risch(x/(x^2+1)^2,x)",                         // 76
            "risch(1/(x^2+1)^2,x)",                         // 77
            "risch(x*exp(x)*ln(x),x)",                      // 82
            "risch(exp(x)*(1/x-1/x^2),x)",                  // 83
            "risch((exp(x)+1)/(exp(x)+x),x)",               // 84
            "risch(exp(x)/(exp(x)+1)+ln(x),x)",             // 85
            "risch(2^x,x)",                                 // 86
        ],
        false,
    );
}

/// The probe rows giac integrates in part: the primitive and the
/// remainder, each giac's.
#[test]
fn partial_agreements() {
    check(
        &[
            "risch(ln(x+1)*ln(x),x)",                                // 26
            "risch(ln(ln(x)),x)",                                    // 28
            "risch(exp(x)/(x*ln(x))-exp(x)/(x*ln(x)^2)+exp(x)/x,x)", // 81
        ],
        true,
    );
}

/// Ours, as decided, on a departure row.
enum Decided {
    /// The call stays: `Risch` (or the user's `Int`) unevaluated.
    Stays,
    /// An answer without a remainder whose derivative is the integrand.
    Integrates,
    /// An answer, maybe partial or the call itself, with
    /// `D[prim] + rem = f`: not giac's split.
    Holds,
}

/// The departure rows: `(input, giac's text, ours as decided)`.
const DEPARTURES: [(&str, &str, Decided); 22] = [
    // giac leaves these (probe rows 9-13, 43, 69, 70, 87) and prints
    // `integrate(f,x)`, which our `Integrate` evaluates on reading back;
    // ours stays `risch(f,x)`.
    ("risch(exp(x^2),x)", "integrate(exp(x^2),x)", Decided::Stays),
    ("risch(1/ln(x),x)", "integrate(1/ln(x),x)", Decided::Stays),
    ("risch(sqrt(x),x)", "integrate(sqrt(x),x)", Decided::Stays),
    (
        "risch(x^2*exp(-x^2),x)",
        "integrate(x^2*exp(-x^2),x)",
        Decided::Stays,
    ),
    (
        "risch((x^2+1)/(x^3+x+1),x)",
        "integrate((x^2+1)/(x^3+x+1),x)",
        Decided::Stays,
    ),
    (
        "risch(exp(exp(x)),x)",
        "integrate(exp(exp(x)),x)",
        Decided::Stays,
    ),
    (
        "risch(x^(1/2)*exp(x),x)",
        "integrate(sqrt(x)*exp(x),x)",
        Decided::Stays,
    ),
    (
        "risch(exp(x)*sqrt(x),x)",
        "integrate(exp(x)*sqrt(x),x)",
        Decided::Stays,
    ),
    ("risch(x^x,x)", "integrate(x^x,x)", Decided::Stays),
    // A Rothstein–Trager factor of degree 2 split by a square root
    // (giac's with_sqrt factoring), not ported: the call stays, or
    // for x^3/(x^2+x+1) its polynomial part integrates and the rest
    // stays (probe rows 73, 78-80).
    (
        "risch(1/(x^2-2),x)",
        "-sqrt(2)/4*ln(x+sqrt(2))+2*sqrt(2)/8*ln(x-sqrt(2))",
        Decided::Stays,
    ),
    (
        "risch(1/(x^3+1),x)",
        "-(i*3*sqrt(3)+3)/18*ln(2*x+(-i)*sqrt(3)-1)+1/3*ln(x+1)-((-i)*3*sqrt(3)+3)/18*\
         ln(2*x+i*sqrt(3)-1)",
        Decided::Stays,
    ),
    (
        "risch(1/(x^4+1),x)",
        "-sqrt(2)*(1+i)/8*ln((1+i)*x+(-i)*sqrt(2))-sqrt(2)*(1-i)/8*ln((1+i)*x-sqrt(2))+\
         2*sqrt(2)*(1-i)/16*ln((1+i)*x+sqrt(2))+2*sqrt(2)*(1+i)/16*ln((1+i)*x+i*sqrt(2))",
        Decided::Stays,
    ),
    (
        "risch(x^3/(x^2+x+1),x)",
        "-i*sqrt(3)/3*ln(2*x+(-i)*sqrt(3)+1)-(-i)*sqrt(3)/3*ln(2*x+i*sqrt(3)+1)+1/2*x^2-x",
        Decided::Holds,
    ),
    // giac's last step re(halftan(…)) writes re/im of Ei(x ± i) with
    // tan(1/2), and carries -pi before atan(exp(x)): ours integrates
    // without that step (probe rows 54, 38, 39).
    (
        "risch(exp(x)/(x^2+1),x)",
        "(im(Ei(x+i))*tan(1/2)^2-im(Ei(x+i))-im(Ei(x-i))*tan(1/2)^2+im(Ei(x-i))+\
         2*re(Ei(x+i))*tan(1/2)+2*re(Ei(x-i))*tan(1/2))/(2*tan(1/2)^2+2)",
        Decided::Integrates,
    ),
    (
        "risch(1/(exp(x)+exp(-x)),x)",
        "-pi+atan(exp(x))",
        Decided::Integrates,
    ),
    (
        "risch(exp(x)/(exp(2*x)+1),x)",
        "-pi+atan(exp(x))",
        Decided::Integrates,
    ),
    // Tower ties (T9.c1 batch rows 19 and 89): kernels of one size
    // keep our canonical order, giac the order of its own printout, so
    // the split differs; ours passes the D check or stays.
    (
        "risch((2/(1*x+0))*ln(3*x+-1)^2+(-1/(1*x+1))*exp(2*x+-1)*ln(2*x+-1)+\
         (2*x^2)*1/(exp(1*x+0)+1),x)",
        "-Ei(2*x+2)*ln(2*x-1)/exp(1)^3+integrate((4*exp(x)*x*exp(1)*exp(2)-2*exp(x)*\
         exp(1)*exp(2)+4*x*exp(1)*exp(2)-2*exp(1)*exp(2))/(2*exp(x)*x^2*exp(1)*exp(2)-\
         exp(x)*x*exp(1)*exp(2)+2*x^2*exp(1)*exp(2)-x*exp(1)*exp(2))*ln(3*x-1)^2+\
         (2*exp(x)*x*Ei(2*(x+1))+4*x^4*exp(1)*exp(2)-2*x^3*exp(1)*exp(2)+2*x*Ei(2*(x+1)))/\
         (2*exp(x)*x^2*exp(1)*exp(2)-exp(x)*x*exp(1)*exp(2)+2*x^2*exp(1)*exp(2)-x*exp(1)*\
         exp(2)),x)",
        Decided::Holds,
    ),
    (
        "risch((-3/(3*x+-1))*exp(1*x+0)/(exp(2*x+1)+2)+(-1/(2*x+2))*ln(x^2+-1*x+-2)+\
         ((-2*x+-2)/(3*x+2))*exp(2*x+2)*ln(1*x+-2),x)",
        "-1/2*ln(x+1)*ln(x+1)+integrate((-12*exp(x)^4*x^3*exp(1)^3-20*exp(x)^4*x^2*\
         exp(1)^3-4*exp(x)^4*x*exp(1)^3+4*exp(x)^4*exp(1)^3-24*exp(x)^2*x^3*exp(1)^2-\
         40*exp(x)^2*x^2*exp(1)^2-9*exp(x)^2*x^2*exp(1)-8*exp(x)^2*x*exp(1)^2-3*exp(x)^2*\
         x*exp(1)+8*exp(x)^2*exp(1)^2+2*exp(x)^2*exp(1)-18*x^2-6*x+4)/(18*exp(x)^2*x^3*\
         exp(1)+24*exp(x)^2*x^2*exp(1)+2*exp(x)^2*x*exp(1)-4*exp(x)^2*exp(1)+36*x^3+\
         48*x^2+4*x-8)*ln(x-2)+(9*exp(x)^2*x^2*exp(1)*ln(x+1)+3*exp(x)^2*x*exp(1)*\
         ln(x+1)-2*exp(x)^2*exp(1)*ln(x+1)-18*exp(x)*x^2-30*exp(x)*x-12*exp(x)+18*x^2*\
         ln(x+1)+6*x*ln(x+1)-4*ln(x+1))/(18*exp(x)^2*x^3*exp(1)+24*exp(x)^2*x^2*exp(1)+\
         2*exp(x)^2*x*exp(1)-4*exp(x)^2*exp(1)+36*x^3+48*x^2+4*x-8),x)",
        Decided::Holds,
    ),
    // giac's wrong answers (batch rows 30 and 65): its primitive fails
    // differentiation (risch.cc:693 drops the one allowed logarithm;
    // -x-x), ours fails the D check and stays.
    (
        "risch((-1*x^1)*ln(1*x+1)+ln(3*x+-2)^2/(1*x^1)+(-3)*exp(x^2+0*x+-2),x)",
        "1/4*x^2-1/2*x+1/2*ln(x+1)-1/2*x^2*ln(x+1)+ln(x)*ln(3*x-2)^2+\
         integrate(-3*exp(x^2-2),x)",
        Decided::Stays,
    ),
    (
        "risch((-1)*1/(exp(3*x+-1)+1)+1/(exp(2*x+1)+2)/(1/(2*x+1)),x)",
        "-x-x+integrate((2*exp(x)^3*x+exp(x)^3-exp(x)^2*exp(1)^2+2*x*exp(1)-exp(1))/\
         (exp(x)^5*exp(1)+2*exp(x)^3+exp(x)^2*exp(1)^2+2*exp(1)),x)",
        Decided::Stays,
    ),
    // giac's Int is integrate unevaluated (probe5.out:22, probe2.out:9);
    // a capitalised name here is the user's call (R1.2-001), which stays.
    ("Int(sin(x),x)", "integrate(sin(x),x)", Decided::Stays),
    ("Int(sin(t),t)", "integrate(sin(t),t)", Decided::Stays),
];

/// Where giac's answer is not ours: giac's text, and ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, answer, decided) in DEPARTURES {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
        let mine = match ours(input) {
            Ok(mine) => mine,
            Err(e) => {
                wrong.push(format!("{input}: {e}"));
                continue;
            }
        };
        let holds = |f: Result<Expr, String>| match f {
            Ok(f) => derivative_holds(&ctx, &mine, &f),
            Err(_) => false,
        };
        let fine = match decided {
            Decided::Stays => {
                let head = if input.starts_with("Int(") {
                    Symbol::new("$Int")
                } else {
                    Symbol::RISCH
                };
                mine.head().and_then(Expr::as_symbol) == Some(head)
            }
            Decided::Integrates => !has_risch(&mine) && holds(integrand(input)),
            Decided::Holds => holds(integrand(input)),
        };
        if !fine {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
